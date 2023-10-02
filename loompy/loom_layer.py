import numpy as np
from typing import *
import scipy
from loompy import timestamp


class MemoryLoomLayer():
	"""
	A layer residing in memory (without a corresponding layer on disk), typically
	as part of a :class:`loompy.LoomView`. MemoryLoomLayer supports a subset of 
	the operations suported for regular layers.
	"""
	def __init__(self, name: str, matrix: np.ndarray) -> None:
		self.name = name  #: Name of the layer
		self.shape = matrix.shape  #: Shape of the layer
		self.values = matrix

	def __getitem__(self, slice: Tuple[Union[int, slice], Union[int, slice]]) -> np.ndarray:
		return self.values[slice]

	def __setitem__(self, slice: Tuple[Union[int, slice], Union[int, slice]], data: np.ndarray) -> None:
		self.values[slice] = data

	def sparse(self, rows: np.ndarray, cols: np.ndarray) -> scipy.sparse.coo_matrix:
		"""
		Return the layer as :class:`scipy.sparse.coo_matrix`
		"""
		return scipy.sparse.coo_matrix(self.values[rows, :][:, cols])

	def permute(self, ordering: np.ndarray, *, axis: int) -> None:
		"""
		Permute the layer along an axis

		Args:
			axis: The axis to permute (0, permute the rows; 1, permute the columns)
			ordering: The permutation vector
		"""
		if axis == 0:
			self.values = self.values[ordering, :]
		elif axis == 1:
			self.values = self.values[:, ordering]
		else:
			raise ValueError("axis must be 0 or 1")


class LoomLayer():
	"""
	Represents a layer (matrix) of values in the loom file, which can be accessed by slicing.
	"""
	
	def __init__(self, name: str, ds: Any) -> None:
		self.ds = ds  #: The :class:`.LoomConnection` object this layer belongs to
		self.name = name  #: Name of the layer (str)
		self.shape = ds.shape  #: Shape of the layer, tuple of (n_rows, n_cols)
		self.dtype = ""  #: Datatype of the layer (str)
		if name == "":
			self.dtype = self.ds._file["/matrix"].dtype
		else:
			self.dtype = self.ds._file["/layers/" + self.name].dtype

	def last_modified(self) -> str:
		"""
		Return a compact ISO8601 timestamp (UTC timezone) indicating when the file was last modified

		Note: if the layer does not contain a timestamp, and the mode is 'r+', a new timestamp will be set and returned.
		Otherwise, the current time in UTC will be returned.
		"""
		if self.name == "":
			if "last_modified" in self.ds._file["/matrix"].attrs:
				return self.ds._file["/matrix"].attrs["last_modified"]
			elif self.ds._file.mode == 'r+':
				self.ds._file["/matrix"].attrs["last_modified"] = timestamp()
				self.ds._file.flush()
				return self.ds._file["/matrix"].attrs["last_modified"]

		if self.name != "":
			if "last_modified" in self.ds._file["/layers/" + self.name].attrs:
				return self.ds._file["/layers/" + self.name].attrs["last_modified"]
			elif self.ds._file.mode == 'r+':
				self.ds._file["/layers/" + self.name].attrs["last_modified"] = timestamp()
				self.ds._file.flush()
				return self.ds._file["/layers/" + self.name].attrs["last_modified"]

		return timestamp()

	def __getitem__(self, slice: Tuple[Union[int, slice], Union[int, slice]]) -> np.ndarray:
		if self.name == "":
			return self.ds._file['/matrix'].__getitem__(slice)
		return self.ds._file['/layers/' + self.name].__getitem__(slice)

	def __setitem__(self, slice: Tuple[Union[int, slice], Union[int, slice]], data: np.ndarray) -> None:
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

	def sparse(self, rows: np.ndarray = None, cols: np.ndarray = None, dtype = None) -> scipy.sparse.coo_matrix:
		if rows is not None:
			if np.issubdtype(rows.dtype, np.bool_):
				rows = np.where(rows)[0]
		if cols is not None:
			if np.issubdtype(cols.dtype, np.bool_):
				cols = np.where(cols)[0]
				
		n_genes = self.ds.shape[0] if rows is None else rows.shape[0]
		n_cells = self.ds.shape[1] if cols is None else cols.shape[0]
		# Calculate sparse data length to be able to reserve proper sized arrays beforehand
		nnonzero = 0
		for (ix, selection, view) in self.ds.scan(items=cols, axis=1, layers=[self.name], what=["layers"], batch_size=4096):
			if rows is not None:
				vals = view.layers[self.name][rows, :]
			else:
				vals = view.layers[self.name][:, :]
			nnonzero += np.count_nonzero(vals)
		data = np.empty((nnonzero,), dtype=dtype) #(data: List[np.ndarray] = []
		row = np.empty((nnonzero,), dtype=('uint16' if self.ds.shape[0] < 2**16 else 'uint32')) #row : List[np.ndarray] = []
		col = np.empty((nnonzero,), dtype=('uint32' if self.ds.shape[1] < 2**32 else 'uint64')) #col: List[np.ndarray] = []
		i = 0
		ci = 0
		for (ix, selection, view) in self.ds.scan(items=cols, axis=1, layers=[self.name], what=["layers"], batch_size=4096):
			if rows is not None:
				vals = view.layers[self.name][rows, :]
			else:
				vals = view.layers[self.name][:, :]
			if dtype:
				vals = vals.astype(dtype)
			nonzeros = np.where(vals != 0)
			n = len(nonzeros[0])
			data[ci:ci+n] = vals[nonzeros] #data.append(vals[nonzeros])
			row[ci:ci+n] = nonzeros[0] #row.append(nonzeros[0])
			col[ci:ci+n] = (nonzeros[1]+i) #col.append(nonzeros[1] + i)
			ci += n
			i += selection.shape[0]
		return scipy.sparse.coo_matrix((data, (row, col)), shape=(n_genes, n_cells), dtype=dtype)
		#return scipy.sparse.coo_matrix((np.concatenate(data, dtype=dtype), (np.concatenate(row), np.concatenate(col))), shape=(n_genes, n_cells), dtype=dtype)

	def _resize(self, size: Tuple[int, int], axis: int = None) -> None:
		"""Resize the dataset, or the specified axis.

		The dataset must be stored in chunked format; it can be resized up to the "maximum shape" (keyword maxshape) specified at creation time.
		The rank of the dataset cannot be changed.
		"Size" should be a shape tuple, or if an axis is specified, an integer.

		BEWARE: This functions differently than the NumPy resize() method!
		The data is not "reshuffled" to fit in the new shape; each axis is grown or shrunk independently.
		The coordinates of existing data are fixed.
		"""
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

	def _permute(self, ordering: np.ndarray, *, axis: int) -> None:
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
