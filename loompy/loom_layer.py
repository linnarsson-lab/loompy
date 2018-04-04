import numpy as np
from typing import *
import scipy
from loompy import timestamp


class MemoryLoomLayer():
	def __init__(self, name: str, matrix: np.ndarray) -> None:
		self.name = name
		self.shape = matrix.shape
		self.values = matrix

	def __getitem__(self, slice: Tuple[Union[int, slice], Union[int, slice]]) -> np.ndarray:
		return self.values[slice]

	def __setitem__(self, slice: Tuple[Union[int, slice], Union[int, slice]], data: np.ndarray) -> None:
		self.values[slice] = data

	def sparse(self, rows: np.ndarray, cols: np.ndarray) -> scipy.sparse.coo_matrix:
		return scipy.sparse.coo_matrix(self.values[rows, :][:, cols])

	def permute(self, ordering: np.ndarray, *, axis: int) -> None:
		if axis == 0:
			self.values = self.values[ordering, :]
		elif axis == 1:
			self.values = self.values[:, ordering]
		else:
			raise ValueError("axis must be 0 or 1")


class LoomLayer():
	def __init__(self, name: str, ds: Any) -> None:
		self.ds = ds
		self.name = name
		self.shape = ds.shape

	def last_modified(self) -> str:
		"""
		Return a compact ISO8601 timestamp (UTC timezone) indicating when the file was last modified

		If a layer has no timestamps and the loom file is connected in read-write mode, it will be initialised to current time.

		If a layer has no timestamps and the loom file is connected in read-only mode, "19700101T000000Z" (start of Unix Time) is returned.
		"""
		if self.name == "":
			if "last_modified" in self.ds._file["/matrix"].attrs:
				return self.ds._file["/matrix"].attrs["last_modified"]
			elif self.ds.mode == 'r+':
				self.ds._file["/matrix"].attrs["last_modified"] = timestamp()
				self.ds._file.flush()
				return self.ds._file["/matrix"].attrs["last_modified"]

		if self.name != "":
			if "last_modified" in self.ds._file["/layers/" + self.name].attrs:
				return self.ds._file["/layers/" + self.name].attrs["last_modified"]
			elif self.ds.mode == 'r+':
				self.ds._file["/layers/" + self.name].attrs["last_modified"] = timestamp()
				self.ds._file.flush()
				return self.ds._file["/layers/" + self.name].attrs["last_modified"]

		return "19700101T000000Z"

	def __getitem__(self, slice: Tuple[Union[int, slice], Union[int, slice]]) -> np.ndarray:
		if self.name == "":
			return self.ds._file['/matrix'].__getitem__(slice)
		return self.ds._file['/layers/' + self.name].__getitem__(slice)

	def __setitem__(self, slice: Tuple[Union[int, slice], Union[int, slice]], data: np.ndarray) -> None:
		if self.ds is None:
			raise IOError("Layer not modified: LoomConnection is None")
		elif self.ds.closed:
			raise IOError("Layer not modified: cannot modify closed LoomConnection")
		elif self.ds.mode != 'r+':
			raise IOError("Layer not modified: cannot modify loom file when connected in read-only mode")
		else:
			if self.name == "":
				self.ds._file['/matrix'][slice] = data
				self.ds._file["/matrix"].attrs["last_modified"] = timestamp()
				self.ds._file.attrs["last_modified"] = timestamp()
				self.ds._file.flush()
			else:
				self.ds._file['/layers/' + self.name][slice] = data
				self.ds._file["/layers/" + self.name].attrs["last_modified"] = timestamp()
				self.ds._file.attrs["last_modified"] = timestamp()
				self.ds._file.flush()

	def sparse(self, rows: np.ndarray, cols: np.ndarray) -> scipy.sparse.coo_matrix:
		n_genes = self.ds.shape[0] if rows is None else rows.shape[0]
		n_cells = self.ds.shape[1] if cols is None else cols.shape[0]

		data: List[np.ndarray] = []
		row: List[np.ndarray] = []
		col: List[np.ndarray] = []
		for (ix, selection, view) in self.ds.scan(items=cols, axis=1, layers=[self.name]):
			if rows is not None:
				vals = view.layers[self.name][rows, :]
			else:
				vals = view.layers[self.name][:, :]
			nonzeros = np.where(vals > 0)
			data.append(vals[nonzeros])
			row.append(nonzeros[0])
			col.append(nonzeros[1] + ix)
		return scipy.sparse.coo_matrix((np.concatenate(data), (np.concatenate(row), np.concatenate(col))), shape=(n_genes, n_cells))

	def resize(self, size: Tuple[int, int], axis: int = None) -> None:
		"""Resize the dataset, or the specified axis.

		The dataset must be stored in chunked format; it can be resized up to the "maximum shape" (keyword maxshape) specified at creation time.
		The rank of the dataset cannot be changed.
		"Size" should be a shape tuple, or if an axis is specified, an integer.

		BEWARE: This functions differently than the NumPy resize() method!
		The data is not "reshuffled" to fit in the new shape; each axis is grown or shrunk independently.
		The coordinates of existing data are fixed.
		"""
		if self.ds is None:
			raise IOError("Dataset not resized: LoomConnection is None")
		elif self.ds.mode != 'r+':
			raise IOError("Dataset not resized: cannot modify loom file when connected in read-only mode")
		else:
			if self.name == "":
				self.ds._file['/matrix'].resize(size, axis)
			else:
				self.ds._file['/layers/' + self.name].resize(size, axis)

	def map(self, f_list: List[Callable[[np.ndarray], int]], axis: int = 0, chunksize: int = 1000, selection: np.ndarray = None) -> List[np.ndarray]:
		"""
		Apply a function along an axis without loading the entire dataset in memory.

		Args:
			f_list (list of func):		Function(s) that takes a numpy ndarray as argument

			axis (int):		Axis along which to apply the function (0 = rows, 1 = columns)

			chunksize (int): Number of rows (columns) to load per chunk

			selection (array of bool): Columns (rows) to include

		Returns:
			numpy.ndarray result of function application

			If you supply a list of functions, the result will be a list of numpy arrays. This is more
			efficient than repeatedly calling map() one function at a time.
		"""
		if hasattr(f_list, '__call__'):
			raise ValueError("f_list must be a list of functions, not a function itself")

		result = []
		if axis == 0:
			rows_per_chunk = chunksize
			for i in range(len(f_list)):
				result.append(np.zeros(self.shape[0]))
			ix = 0
			while ix < self.shape[0]:
				rows_per_chunk = min(self.shape[0] - ix, rows_per_chunk)
				if selection is not None:
					chunk = self[ix:ix + rows_per_chunk, :][:, selection]
				else:
					chunk = self[ix:ix + rows_per_chunk, :]
				for i in range(len(f_list)):
					result[i][ix:ix + rows_per_chunk] = np.apply_along_axis(f_list[i], 1, chunk)
				ix = ix + rows_per_chunk
		elif axis == 1:
			cols_per_chunk = chunksize
			for i in range(len(f_list)):
				result.append(np.zeros(self.shape[1]))
			ix = 0
			while ix < self.shape[1]:
				cols_per_chunk = min(self.shape[1] - ix, cols_per_chunk)
				if selection is not None:
					chunk = self[:, ix:ix + cols_per_chunk][selection, :]
				else:
					chunk = self[:, ix:ix + cols_per_chunk]
				for i in range(len(f_list)):
					result[i][ix:ix + cols_per_chunk] = np.apply_along_axis(f_list[i], 0, chunk)
				ix = ix + cols_per_chunk
		return result

	def permute(self, ordering: np.ndarray, *, axis: int) -> None:
		if self.ds is None:
			raise IOError("Layer not permuted: LoomConnection is None")
		elif self.ds.closed:
			raise IOError("Layer not permuted: cannot modify closed LoomConnection")
		elif self.ds.mode != 'r+':
			raise IOError("Layer not permuted: cannot modify loom file when connected in read-only mode")
		else:
			if self.name == "":
				obj = self.ds._file['/matrix']
			else:
				obj = self.ds._file['/layers/' + self.name]
			if axis == 0:
				chunksize = 5000
				start = 0
				while start < self.shape[1]:
					submatrix = obj[:, start:start + chunksize]
					obj[:, start:start + chunksize] = submatrix[ordering, :]
					start = start + chunksize
			elif axis == 1:
				chunksize = 100000000 // self.shape[1]
				start = 0
				while start < self.shape[0]:
					submatrix = obj[start:start + chunksize, :]
					obj[start:start + chunksize, :] = submatrix[:, ordering]
					start = start + chunksize
			else:
				raise ValueError("axis must be 0 or 1")