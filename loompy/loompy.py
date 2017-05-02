# Copyright (c) 2015 Sten Linnarsson
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np
from typing import *
import h5py
import os.path
from scipy.io import mmread
from shutil import copyfile
import logging
import time


def strip(s: str) -> str:
	if s[0:2] == "b'" and s[-1] == "'":
		return s[2:-1]
	return s


class LoomAttributeManager(object):
	__slots__ = ('f',)

	def __init__(self, f: h5py.File) -> None:
		self.f = f

	def __contains__(self, name: str) -> bool:
		return self.f.attrs.__contains__(name)

	def __setitem__(self, name: str, value: str) -> None:
		if type(value) == bytes:
			self.f.attrs[name] = value
		else:
			self.f.attrs[name] = value.encode('utf-8')
		self.f.flush()

	def __getitem__(self, name: str) -> str:
		val = self.f.attrs[name]
		if type(val) == bytes:
			val = val.decode('utf-8')
		else:
			val = str(val)

		# Fix cosmetic bugs accidentally introduced by Python 2/3 bug
		if val[0:2] == "b'" and val[-1] == "'":
			val = val[2:-1]

		return val

	def __iter__(self) -> Iterator[str]:
		for val in self.f.attrs:
			yield val

	def __len__(self) -> int:
		return len(self.f.attrs)

	def get(self, name: str, default: str = None) -> str:
		if self.__contains__(name):
			return self[name]
		else:
			return default


class LoomConnection:
	def __init__(self, filename: str, mode: str = 'r+') -> None:
		"""
		Establish a connection to a .loom file.

		Args:
			filename:			Name of the .loom file to open
			mode:				read/write mode, accepts 'r+' (read/write) or
								'r' (read-only), defaults to 'r+' without arguments,
								and to 'r' with incorrect arguments

		Returns:
			Nothing.

		Row and column attributes are loaded into memory for fast access.
		"""

		# make sure a valid mode was passed, if not default to read-only
		# because you probably are doing something that you don't want to
		if mode != 'r+' and mode != 'r':
			logging.warn("Wrong mode passed to LoomConnection, using read-only to not destroy data")
			mode = 'r'
		self.mode = mode
		self.filename = filename
		self._file = h5py.File(filename, mode)
		self.shape = [0, 0]  # The correct shape gets assigned when the layers are loaded

		if self._file.__contains__("/matrix"):
			self.layer = {
				"@DEFAULT": LoomLayer(self, "@DEFAULT", self._file["/matrix"].dtype)
			}
			self.shape = self._file["/matrix"].shape
			if self._file.__contains__("/layers"):
				for key in self._file["/layers"].keys():
					self.layer[key] = LoomLayer(self, key, self._file["/layers/" + key].dtype)
		else:
			self.layer = {}

		self.row_attrs = {}  # type: Dict[str, np.ndarray]
		for key in self._file['row_attrs'].keys():
			self._load_attr(key, axis=0)
			v = self.row_attrs[key]
			if type(v[0]) is np.str_ and len(v[0]) >= 3 and v[0][:2] == "b'" and v[0][-1] == "'":
				logging.warn("Unicode bug detected in row %s" % key)
				if mode == 'r+':
					logging.warn("Fixing unicode bug by re-setting row attribute '" + key + "'")
					self._save_attr(key, np.array([x[2:-1] for x in v]), axis=0)
					self._load_attr(key, axis=0)

		self.col_attrs = {}  # type: Dict[str, np.ndarray]
		for key in self._file['col_attrs'].keys():
			self._load_attr(key, axis=1)
			v = self.col_attrs[key]
			if type(v[0]) is np.str_ and len(v[0]) >= 3 and v[0][:2] == "b'" and v[0][-1] == "'":
				logging.warn("Unicode bug detected in column %s" % key)
				if mode == 'r+':
					logging.warn("Fixing unicode bug by re-setting column attribute '" + key + "'")
					self._save_attr(key, np.array([x[2:-1] for x in v]), axis=1)
					self._load_attr(key, axis=1)

		self.attrs = LoomAttributeManager(self._file)

	def _save_attr(self, name: str, values: np.ndarray, axis: int) -> None:
		"""
		Save an attribute to the file, nothing else

		Remarks:
			Handles unicode to ascii conversion (lossy, but HDF5 supports only ascii)
			Does not update the attribute cache (use _load_attr for this)
		"""
		if self.mode != "r+":
			raise IOError("Cannot save attributes when connected in read-only mode")
		if values.dtype.type is np.str_:
			values = np.array([x.encode('ascii', 'ignore') for x in values])

		a = ["/row_attrs/", "/col_attrs/"][axis]
		if self.shape[axis] != 0 and len(values) != self.shape[axis]:
			raise ValueError("Attribute must have exactly %d values" % self.shape[axis])
		if self._file[a].__contains__(name):
			del self._file[a + name]
		self._file[a + name] = values
		self._file.flush()

	def _load_attr(self, name: str, axis: int) -> None:
		"""
		Load an attribute from the file, nothing else

		Remarks:
			Handles ascii to unicode conversion
			Updates the attribute cache as well as the class attributes
		"""
		a = ["/row_attrs/", "/col_attrs/"][axis]

		if self._file[a][name].dtype.kind == 'S':
			vals = np.array([x.decode('utf8') for x in self._file[a][name][:]])
		else:
			vals = self._file[a][name][:]

		if axis == 0:
			self.row_attrs[name] = vals
			if not hasattr(LoomConnection, name):
				setattr(self, name, self.row_attrs[name])
		else:
			self.col_attrs[name] = vals
			if not hasattr(LoomConnection, name):
				setattr(self, name, self.col_attrs[name])

	def _repr_html_(self) -> str:
		"""
		Return an HTML representation of the loom file, showing the upper-left 10x10 corner.
		"""
		rm = min(10, self.shape[0])
		cm = min(10, self.shape[1])
		html = "<p>"
		if self.attrs.__contains__("title"):
			html += "<strong>" + self.attrs["title"] + "</strong> "
		html += "(" + str(self.shape[0]) + " genes, " + str(self.shape[1]) + " cells, " + str(len(self.layer)) + " layers)<br/>"
		html += self._file.filename + "<br/>"
		if self.attrs.__contains__("description"):
			html += "<em>" + self.attrs["description"] + "</em><br/>"
		html += "<table>"
		# Emit column attributes
		for ca in self.col_attrs.keys():
			html += "<tr>"
			for ra in self.row_attrs.keys():
				html += "<td>&nbsp;</td>"  # Space for row attrs
			html += "<td><strong>" + ca + "</strong></td>"  # Col attr name
			for v in self.col_attrs[ca][:cm]:
				html += "<td>" + str(v) + "</td>"
			if self.shape[1] > cm:
				html += "<td>...</td>"
			html += "</tr>"

		# Emit row attribute names
		html += "<tr>"
		for ra in self.row_attrs.keys():
			html += "<td><strong>" + ra + "</strong></td>"  # Row attr name
		html += "<td>&nbsp;</td>"  # Space for col attrs
		for v in range(cm):
			html += "<td>&nbsp;</td>"
		if self.shape[1] > cm:
			html += "<td>...</td>"
		html += "</tr>"

		# Emit row attr values and matrix values
		for row in range(rm):
			html += "<tr>"
			for ra in self.row_attrs.keys():
				html += "<td>" + str(self.row_attrs[ra][row]) + "</td>"
			html += "<td>&nbsp;</td>"  # Space for col attrs

			for v in self[row, :cm]:
				html += "<td>" + str(v) + "</td>"
			if self.shape[1] > cm:
				html += "<td>...</td>"
			html += "</tr>"
		# Emit ellipses
		if self.shape[0] > rm:
			html += "<tr>"
			for v in range(rm + 1 + len(self.row_attrs.keys())):
				html += "<td>...</td>"
			if self.shape[1] > cm:
				html += "<td>...</td>"
			html += "</tr>"
		html += "</table>"
		return html

	def __getitem__(self, slice: Tuple[Union[int, slice], Union[int, slice]]) -> np.ndarray:
		"""
		Get a slice of the main matrix.

		Args:
			slice:		A slice object (see http://docs.h5py.org/en/latest/high/dataset.html)

		Returns:
			A numpy matrix
		"""
		return self.layer["@DEFAULT"][slice]

	def __setitem__(self, slice: Tuple[Union[int, slice], Union[int, slice]], data: np.ndarray) -> None:
		"""
		Assign a slice of the main matrix.

		Args:
			slice:		A slice object (see http://docs.h5py.org/en/latest/high/dataset.html)

		Returns:
			Nothing.
		"""
		self.layer["@DEFAULT"][slice] = data

	def close(self) -> None:
		"""
		Close the connection. After this, the connection object becomes invalid.
		"""
		self._file.close()
		self._file = None
		self.row_attrs = {}
		self.col_attrs = {}
		self.shape = [0, 0]

	def set_layer(self, name: str, matrix: np.ndarray, chunks: Tuple[int, int] = (64, 64), chunk_cache: int = 512, dtype: str = "float32", compression_opts: int = 2) -> None:
		if self.mode != "r+":
			raise IOError("Cannot save attributes when connected in read-only mode")
		if not np.isfinite(matrix).all():
			raise ValueError("INF and NaN not allowed in loom matrix")

		if not self._file.__contains__("/layers"):
			self._file.create_group("/layers")

		# make sure chunk size is not bigger than actual matrix size
		chunks = (min(chunks[0], matrix.shape[0]), min(chunks[1], matrix.shape[1]))
		path = "/layers/" + name
		if name == "@DEFAULT":
			path = "/matrix"
		if self._file.__contains__(path):
			del self._file[path]
		
		# Save the main matrix
		if compression_opts is None:
			self._file.create_dataset(
				path,
				data=matrix.astype(dtype),
				maxshape=(matrix.shape[0], None),
				chunks=chunks,
				fletcher32=False
			)
		else:
			self._file.create_dataset(
				path,
				data=matrix.astype(dtype),
				maxshape=(matrix.shape[0], None),
				chunks=chunks,
				fletcher32=False,
				compression="gzip",
				shuffle=False,
				compression_opts=compression_opts
			)
		
		self.layer[name] = LoomLayer(self, name, dtype)
		if name == "@DEFAULT":
			self.shape = matrix.shape
		self._file.flush()

	def add_columns(self, submatrix: np.ndarray, col_attrs: Dict[str, np.ndarray], fill_values: Dict[str, np.ndarray] = None) -> None:
		"""
		Add columns of data and attribute values to the dataset.

		Args:
			submatrix (dict or numpy.ndarray):
				Either:
				1) A N-by-M matrix of float32s (N rows, M columns) in this case columns are added at the @DEFAULT layer
				2) A dict {layer_name : matrix} specified so that the matrix (N, M) will be added to layer `layer_name`

			col_attrs (dict):
				Column attributes, where keys are attribute names and values are numpy arrays (float or string) of length M

		Returns:
			Nothing.

		Notes
		-----
		- This will modify the underlying HDF5 file, which will interfere with any concurrent readers.
		- Column attributes in the file that are NOT provided, will be deleted.
		- Array with Nan should not be provided

		"""
		if self.mode != "r+":
			raise IOError("Cannot add columns when connected in read-only mode")

		if not type(submatrix) == dict:
			submatrix_dict = dict()
			submatrix_dict["@DEFAULT"] = submatrix
		else:
			submatrix_dict = cast(dict, submatrix)  # equivalent to submatrix_dict = submatrix # only avoids problems with type checker
			submatrix = submatrix_dict["@DEFAULT"]
		
		# for k, v in submatrix_dict.items():
		# 	if not np.isfinite(v).all():
		# 		raise ValueError("INF and NaN not allowed in loom matrix")

		if submatrix.shape[0] != self.shape[0]:
			raise ValueError("New submatrix must have same number of rows as existing matrix")

		did_remove = False
		todel = []  # type: List[str]
		for key, vals in col_attrs.items():
			if key not in self.col_attrs:
				if fill_values is not None:
					if fill_values == "auto":
						fill_with = np.zeros(1, dtype=col_attrs[key].dtype)[0]
					else:
						fill_with = fill_values[key]
					self.set_attr(key, np.array([fill_with] * self.shape[1]), axis=1)
				else:
					did_remove = True
					todel.append(key)
			if len(vals) != submatrix.shape[1]:
				raise ValueError("Each column attribute must have exactly %s values" % submatrix.shape[1])
		for key in todel:
			del col_attrs[key]
		if did_remove:
			logging.warn("Some column attributes were removed: " + ",".join(todel))

		todel = []
		did_remove = False
		for key in self.col_attrs.keys():
			if key not in col_attrs:
				if fill_values is not None:
					if fill_values == "auto":
						fill_with = np.zeros(1, dtype=self.col_attrs[key].dtype)[0]
					else:
						fill_with = fill_values[key]
					col_attrs[key] = np.array([fill_with] * submatrix.shape[1])
				else:
					did_remove = True
					todel.append(key)
		for key in todel:
			self.delete_attr(key, axis=1)
		if did_remove:
			logging.warn("Some column attributes were removed: " + ",".join(todel))

		n_cols = submatrix.shape[1] + self.shape[1]
		for key, vals in col_attrs.items():
			vals = np.array(vals)
			if vals.dtype.type is np.str_:
				vals = np.array([x.encode('ascii', 'ignore') for x in vals])
			temp = self._file['/col_attrs/' + key][:]
			casting_rule_dtype = np.result_type(temp, vals)
			if vals.dtype != casting_rule_dtype:
				vals = vals.astype(casting_rule_dtype)
			if temp.dtype != casting_rule_dtype:
				temp = temp.astype(casting_rule_dtype)
			temp.resize((n_cols,))
			temp[self.shape[1]:] = vals
			del self._file['/col_attrs/' + key]
			self._file['/col_attrs/' + key] = temp
			self.col_attrs[key] = self._file['/col_attrs/' + key]
		
		# Add the columns layerwise
		for key in self.layer.keys():
			self.layer[key].resize(n_cols, axis=1)
			self.layer[key][:, self.shape[1]:n_cols] = submatrix_dict[key].astype(self.layer[key].dtype)
			self._file.flush()

		self.shape = [self.shape[0], n_cols]

	def add_loom(self, other_file: str, key: str = None, fill_values: Dict[str, np.ndarray] = None) -> None:
		"""
		Add the content of another loom file

		Args:
			other_file (str):	filename of the loom file to append
			fill_values (dict): default values to use for missing attributes (or None to drop missing attrs, or 'auto' to fill with sensible defaults)

		Returns:
			Nothing, but adds the loom file. Note that the other loom file must have exactly the same
			number of rows, in the same order, and must have exactly the same column attributes.
			The all the contents including layers but ignores layers in `other_file` that are not already persent in self
		"""
		if self.mode != "r+":
			raise IOError("Cannot add data when connected in read-only mode")
		# Connect to the loom files
		other = connect(other_file)
		# Verify that the row keys are identical
		if key is not None:
			pk1 = other.row_attrs[key]
			pk2 = self.row_attrs[key]
			for ix, val in enumerate(pk1):
				if pk2[ix] != val:
					raise ValueError("Primary keys are not identical")
		diff_layers = set(self.layer.keys()) - set(other.layer.keys())
		if len(diff_layers) > 0:
			raise ValueError("%s is missing a layer, cannot merge with current file. layers missing:%s" % (other_file, diff_layers))

		for (ix, selection, vals) in other.batch_scan_layers(axis=1, layers=self.layer.keys()):
			ca = {key: v[selection] for key, v in other.col_attrs.items()}
			self.add_columns(vals, ca, fill_values)
		other.close()

	def delete_attr(self, name: str, axis: int = 0, raise_on_missing: bool = True) -> None:
		"""
		Permanently delete an existing attribute and all its values

		Args:

			name (str): 	Name of the attribute to remove
			axis (int):		Axis of the attribute (0 = rows, 1 = columns)

		Returns:
			Nothing.
		"""
		if self.mode != "r+":
			raise IOError("Cannot delete attributes when connected in read-only mode")
		if axis == 0:
			if name not in self.row_attrs:
				if raise_on_missing:
					raise KeyError("Row attribute " + name + " does not exist")
				else:
					return
			del self.row_attrs[name]
			del self._file['/row_attrs/' + name]
			if hasattr(self, name):
				delattr(self, name)

		elif axis == 1:
			if name not in self.col_attrs:
				if raise_on_missing:
					raise KeyError("Column attribute " + name + " does not exist")
				else:
					return
			del self.col_attrs[name]
			del self._file['/col_attrs/' + name]
			if hasattr(self, name):
				delattr(self, name)
		else:
			raise ValueError("Axis must be 0 or 1")

		self._file.flush()

	def set_attr(self, name: str, values: np.ndarray, axis: int = 0, dtype: str = None) -> None:
		"""
		Create or modify an attribute.

		Args:
			name (str): 			Name of the attribute
			values (numpy.ndarray):	Array of values of length equal to the axis length
			axis (int):				Axis of the attribute (0 = rows, 1 = columns)

		Returns:
			Nothing.

		This will overwrite any existing attribute of the same name.
		"""
		if self.mode != "r+":
			raise IOError("Cannot save attributes when connected in read-only mode")
		if dtype is not None:
			raise DeprecationWarning("Data type should no longer be provided")

		self.delete_attr(name, axis, raise_on_missing=False)
		self._save_attr(name, values, axis)
		self._load_attr(name, axis)

	def get_edges(self, name: str, axis: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
		if axis == 0:
			return (self._file["/row_edges/" + name + "/a"][:], self._file["/row_edges/" + name + "/b"][:], self._file["/row_edges/" + name + "/w"][:])
		if axis == 1:
			return (self._file["/col_edges/" + name + "/a"][:], self._file["/col_edges/" + name + "/b"][:], self._file["/col_edges/" + name + "/w"][:])
		raise ValueError("Axis must be 0 or 1")

	def set_edges(self, name: str, a: np.ndarray, b: np.ndarray, w: np.ndarray, axis: int = 0) -> None:
		if self.mode != "r+":
			raise IOError("Cannot save edges when connected in read-only mode")
		if not a.dtype.kind == 'i':
			raise ValueError("Nodes must be integers")
		if not b.dtype.kind == 'i':
			raise ValueError("Nodes must be integers")
		if axis == 1:
			if a.max() > self.shape[1] or a.min() < 0:
				raise ValueError("Nodes out of range")
			if b.max() > self.shape[1] or b.min() < 0:
				raise ValueError("Nodes out of range")
			if self._file.__contains__("/col_edges/" + name):
				del self._file["/col_edges/" + name + "/a"]
				del self._file["/col_edges/" + name + "/b"]
				del self._file["/col_edges/" + name + "/w"]
			self._file["/col_edges/" + name + "/a"] = a
			self._file["/col_edges/" + name + "/b"] = b
			self._file["/col_edges/" + name + "/w"] = w
		elif axis == 0:
			if a.max() > self.shape[0] or a.min() < 0:
				raise ValueError("Nodes out of range")
			if b.max() > self.shape[0] or b.min() < 0:
				raise ValueError("Nodes out of range")
			if self._file.__contains__("/row_edges/" + name):
				del self._file["/row_edges/" + name + "/a"]
				del self._file["/row_edges/" + name + "/b"]
				del self._file["/row_edges/" + name + "/w"]
			self._file["/row_edges/" + name + "/a"] = a
			self._file["/row_edges/" + name + "/b"] = b
			self._file["/row_edges/" + name + "/w"] = w
		else:
			raise ValueError("Axis must be 0 or 1")

	def batch_scan(self, cells: np.ndarray = None, genes: np.ndarray = None, axis: int = 0, batch_size: int = 1000) -> Iterable[Tuple[int, np.ndarray, np.ndarray]]:
		"""Performs a batch scan of the loom file

		Args
		----
		cells: np.ndarray
			the indexes [1,2,3,..,1000] of the cells to select
		genes: np.ndarray
			the indexes [1,2,3,..,1000] of the genes to select
		axis: int
			0:rows or 1:cols
		batch_size: int
			the chuncks returned at every element of the iterator

		Returns
		------
		Iterable that yields triplets
		(ix, indexes, vals)

		ix: int
			first position / how many rows/cols have been yielded alredy
		indexes: np.ndarray[int]
			the indexes with the same numbering of the input args cells / genes (i.e. np.arange(len(ds.shape[axis])))
			this is ix + selection
		vals: np.ndarray
			the matrix corresponding to the chunk
		"""
		if cells is None:
			cells = np.fromiter(range(self.shape[1]), dtype='int')
		if genes is None:
			genes = np.fromiter(range(self.shape[0]), dtype='int')
		if axis == 1:
			cols_per_chunk = batch_size
			ix = 0
			while ix < self.shape[1]:
				cols_per_chunk = min(self.shape[1] - ix, cols_per_chunk)

				selection = cells - ix
				# Pick out the cells that are in this batch
				selection = selection[np.where(np.logical_and(selection >= 0, selection < cols_per_chunk))[0]]
				if selection.shape[0] == 0:
					ix += cols_per_chunk
					continue

				# Load the whole chunk from the file, then extract genes and cells using fancy indexing
				vals = self[:, ix:ix + cols_per_chunk]
				vals = vals[genes, :]
				vals = vals[:, selection]

				yield (ix, ix + selection, vals)
				ix += cols_per_chunk

		if axis == 0:
			rows_per_chunk = batch_size
			ix = 0
			while ix < self.shape[0]:
				rows_per_chunk = min(self.shape[0] - ix, rows_per_chunk)

				selection = genes - ix
				# Pick out the genes that are in this batch
				selection = selection[np.where(np.logical_and(selection >= 0, selection < rows_per_chunk))[0]]
				if selection.shape[0] == 0:
					ix += rows_per_chunk
					continue

				# Load the whole chunk from the file, then extract genes and cells using fancy indexing
				vals = self[ix:ix + rows_per_chunk, :]
				vals = vals[selection, :]
				vals = vals[:, cells]
				yield (ix, ix + selection, vals)
				ix += rows_per_chunk

	def batch_scan_layers(self, cells: np.ndarray = None, genes: np.ndarray = None, axis: int = 0, batch_size: int = 1000, layers: Iterable = None) -> Iterable[Tuple[int, np.ndarray, Dict]]:
		"""Performs a batch scan of the loom file dealing with multiple layer files

		Args
		----
		cells: np.ndarray
			the indexes [1,2,3,..,1000] of the cells to select
		genes: np.ndarray
			the indexes [1,2,3,..,1000] of the genes to select
		axis: int
			0:rows or 1:cols
		batch_size: int
			the chuncks returned at every element of the iterator
		layers: iterable
			if specified it will batch scan only accross some of the layers of the loom file 
			i.g. if layers = ["@DEFAULT"] batch_scan_layers is equivalent to batch_scan

		Returns
		------
		Iterable that yields triplets
		(ix, indexes, vals)

		ix: int
			first position / how many rows/cols have been yielded alredy
		indexes: np.ndarray[int]
			the indexes with the same numbering of the input args cells / genes (i.e. np.arange(len(ds.shape[axis])))
			this is ix + selection
		vals: Dict[layername, np.ndarray]
			a dictionary of the matrixes corresponding to the chunks of different layers
		"""
		if cells is None:
			cells = np.fromiter(range(self.shape[1]), dtype='int')
		if genes is None:
			genes = np.fromiter(range(self.shape[0]), dtype='int')
		if layers is None:
			layers = self.layer.keys()
		if axis == 1:
			cols_per_chunk = batch_size
			ix = 0
			while ix < self.shape[1]:
				cols_per_chunk = min(self.shape[1] - ix, cols_per_chunk)

				selection = cells - ix
				# Pick out the cells that are in this batch
				selection = selection[np.where(np.logical_and(selection >= 0, selection < cols_per_chunk))[0]]
				if selection.shape[0] == 0:
					ix += cols_per_chunk
					continue

				# Load the whole chunk from the file, then extract genes and cells using fancy indexing
				vals = dict()
				for key in layers:
					vals[key] = self.layer[key][:, ix:ix + cols_per_chunk]
					vals[key] = vals[key][genes, :]
					vals[key] = vals[key][:, selection]

				yield (ix, ix + selection, vals)
				ix += cols_per_chunk

		if axis == 0:
			rows_per_chunk = batch_size
			ix = 0
			while ix < self.shape[0]:
				rows_per_chunk = min(self.shape[0] - ix, rows_per_chunk)

				selection = genes - ix
				# Pick out the genes that are in this batch
				selection = selection[np.where(np.logical_and(selection >= 0, selection < rows_per_chunk))[0]]
				if selection.shape[0] == 0:
					ix += rows_per_chunk
					continue

				# Load the whole chunk from the file, then extract genes and cells using fancy indexing
				vals = dict()
				for key in layers:
					vals[key] = self.layer[key][ix:ix + rows_per_chunk, :]
					vals[key] = vals[key][selection, :]
					vals[key] = vals[key][:, cells]
				yield (ix, ix + selection, vals)
				ix += rows_per_chunk

	def map(self, f_list: List[Callable[[np.ndarray], int]], axis: int = 0, chunksize: int = 1000, selection: np.ndarray = None) -> List[np.ndarray]:
		"""
		Apply a function along an axis without loading the entire dataset in memory.

		Args:
			f (list of func):		Function(s) that takes a numpy ndarray as argument

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

	def permute(self, ordering: np.ndarray, axis: int) -> None:
		"""
		Permute the dataset along the indicated axis.

		Args:
			ordering (list of int): 	The desired order along the axis

			axis (int):					The axis along which to permute

		Returns:
			Nothing.
		"""
		if self._file.__contains__("tiles"):
			del self._file['tiles']

		ordering = list(np.array(ordering).flatten())  # Flatten the ordering, in case we got a column vector
		if axis == 0:
			for layer in self.layer:
				if layer == "@DEFAULT":
					obj = self._file['/matrix']
				else:
					obj = self._file['/layers/' + layer]
				chunksize = 5000
				start = 0
				while start < self.shape[1]:
					submatrix = obj[:, start:start + chunksize]
					obj[:, start:start + chunksize] = submatrix[ordering, :]
					start = start + chunksize
			for key in list(self.row_attrs.keys()):
				self.set_attr(key, self.row_attrs[key][ordering], axis=0)
			self._file.flush()
		if axis == 1:
			for layer in self.layer:
				if layer == "@DEFAULT":
					obj = self._file['/matrix']
				else:
					obj = self._file['/layers/' + layer]
				chunksize = 100000000 // self.shape[1]
				start = 0
				while start < self.shape[0]:
					submatrix = obj[start:start + chunksize, :]
					obj[start:start + chunksize, :] = submatrix[:, ordering]
					start = start + chunksize
			for key in list(self.col_attrs.keys()):
				self.set_attr(key, self.col_attrs[key][ordering], axis=1)
			self._file.flush()


class LoomLayer():
	def __init__(self, ds: LoomConnection, name: str, dtype: str) -> None:
		self.ds = ds
		self.name = name
		self.dtype = dtype
		self.shape = ds.shape

	def __getitem__(self, slice: Tuple[Union[int, slice], Union[int, slice]]) -> np.ndarray:
		if self.name == "@DEFAULT":
			return self.ds._file['/matrix'].__getitem__(slice)
		return self.ds._file['/layers/' + self.name].__getitem__(slice)

	def __setitem__(self, slice: Tuple[Union[int, slice], Union[int, slice]], data: np.ndarray) -> None:
		if self.name == "@DEFAULT":
			self.ds._file['/matrix'].__setitem__(slice, data.astype(self.dtype))
		else:
			self.ds._file['/layers/' + self.name].__setitem__(slice, data.astype(self.dtype))

	def resize(self, size: Tuple[int, int], axis: int = None) -> None:
		"""Resize the dataset, or the specified axis.
		
		The dataset must be stored in chunked format; it can be resized up to the "maximum shape" (keyword maxshape) specified at creation time.
		The rank of the dataset cannot be changed.
		"Size" should be a shape tuple, or if an axis is specified, an integer.
		
		BEWARE: This functions differently than the NumPy resize() method!
		The data is not "reshuffled" to fit in the new shape; each axis is grown or shrunk independently.
		The coordinates of existing data are fixed.
		"""
		if self.name == "@DEFAULT":
			self.ds._file['/matrix'].resize(size, axis)
		else:
			self.ds._file['/layers/' + self.name].resize(size, axis)
		

def create(filename: str, matrix: np.ndarray, row_attrs: Dict[str, np.ndarray], col_attrs: Dict[str, np.ndarray], file_attrs: Dict[str, str] = None, chunks: Tuple[int, int] = (64, 64), chunk_cache: int = 512, dtype: str = "float32", compression_opts: int = 2) -> LoomConnection:
	"""
	Create a new .loom file from the given data.

	Args:
		filename (str):         The filename (typically using a `.loom` file extension)
		matrix (numpy.ndarray): Two-dimensional (N-by-M) numpy ndarray of float values
		row_attrs (dict):       Row attributes, where keys are attribute names and values
								are numpy arrays (float or string) of length N
		col_attrs (dict):       Column attributes, where keys are attribute names and
								values are numpy arrays (float or string) of length M
		chunks (tuple):         The chunking of the matrix. Small chunks are slow
								when loading a large batch of rows/columns in sequence,
								but fast for single column/row retrieval.
								Defaults to (64,64).
		chunk_cache (int):      Sets the chunk cache used by the HDF5 format inside
								the loom file, in MB. If the cache is too small to
								contain all chunks of a row/column in memory, then
								sequential row/column access will be a lot slower.
								Defaults to 512.
		matrix_dtype (str):     Dtype of the matrix. Default float32 (uint16, float16 could be used)
		compression_opts (int): Strenght of the gzip compression. Default None.
	Returns:
		LoomConnection to created loom file.
	"""
	if file_attrs is None:
		file_attrs = {}

	# Create the file (empty).
	f = h5py.File(name=filename, mode='w')
	f.create_group('/layers')
	f.create_group('/row_attrs')
	f.create_group('/col_attrs')
	f.flush()
	f.close()

	ds = connect(filename)
	ds.set_layer("@DEFAULT", matrix, chunks, chunk_cache, dtype, compression_opts)

	for key, vals in row_attrs.items():
		ds.set_attr(key, vals, axis=0)

	for key, vals in col_attrs.items():
		ds.set_attr(key, vals, axis=1)

	for vals in file_attrs:
		ds.attrs[vals] = file_attrs[vals]

	# store creation date
	currentTime = time.localtime(time.time())
	ds.attrs['creation_date'] = time.strftime('%Y/%m/%d %H:%M:%S', currentTime)
	ds.attrs['chunks'] = str(chunks)
	return ds


def create_from_cellranger(folder: str, loom_file: str, cell_id_prefix: str = '', sample_annotation: Dict[str, np.ndarray] = None, genome: str = 'mm10') -> LoomConnection:
	"""
	Create a .loom file from 10X Genomics cellranger output

	Args:
		folder (str):				path to the cellranger output folder (usually called `outs`)
		loom_file (str):			full path of the resulting loom file
		cell_id_prefix (str):		prefix to add to cell IDs (e.g. the sample id for this sample)
		sample_annotation (dict): 	dict of additional sample attributes
		genome (str):				genome build to load (e.g. 'mm10')

	Returns:
		Nothing, but creates loom_file
	"""
	if sample_annotation is None:
		sample_annotation = {}
	matrix_folder = os.path.join(folder, 'filtered_gene_bc_matrices', genome)
	matrix = mmread(os.path.join(matrix_folder, "matrix.mtx")).astype("float32").todense()

	barcodes = np.loadtxt(os.path.join(matrix_folder, "barcodes.tsv"), delimiter="\t", dtype="unicode")
	col_attrs = {"CellID": np.array([(cell_id_prefix + strip(bc)[2:-1]) for bc in barcodes])}

	temp = np.loadtxt(os.path.join(matrix_folder, "genes.tsv"), delimiter="\t", dtype="unicode")
	row_attrs = {"Accession": temp[:, 0], "Gene": temp[:, 1]}

	for key in sample_annotation.keys():
		col_attrs[key] = np.array([sample_annotation[key]] * matrix.shape[1])

	tsne_file = os.path.join(folder, "analysis", "tsne", "projection.csv")
	# In cellranger V2 the file moved one level deeper
	if not os.path.exists(tsne_file):
		tsne_file = os.path.join(folder, "analysis", "tsne", "2_components", "projection.csv")
	if os.path.exists(tsne_file):
		tsne = np.loadtxt(tsne_file, usecols=(1, 2), delimiter=',', skiprows=1)
		col_attrs["_tSNE1"] = tsne[:, 0].astype('float64')
		col_attrs["_tSNE2"] = tsne[:, 1].astype('float64')

	pca_file = os.path.join(folder, "analysis", "pca", "projection.csv")
	if os.path.exists(pca_file):
		pca = np.loadtxt(pca_file, usecols=(1, 2), delimiter=',', skiprows=1)
		col_attrs["_PC1"] = pca[:, 0].astype('float64')
		col_attrs["_PC2"] = pca[:, 1].astype('float64')

	kmeans = np.loadtxt(os.path.join(folder, "analysis", "kmeans", "10_clusters", "clusters.csv"), usecols=(1, ), delimiter=',', skiprows=1)
	col_attrs["_KMeans_10"] = kmeans.astype('float64')

	return create(loom_file, matrix, row_attrs, col_attrs)


def combine(files: List[str], output_file: str, key: str = None, file_attrs: Dict[str, str] = None) -> None:
	"""
	Combine two or more loom files and save as a new loom file

	Args:
		files (list of str):	the list of input files (full paths)
		output_file (str):		full path of the output loom file
		key (string):			Row attribute to use to verify row ordering
		file_attrs (dict):		file attributes (title, description, url, etc.)

	Returns:
		Nothing, but creates a new loom file combining the input files.

	The input files must (1) have exactly the same number of rows and in the same order, (2) have
	exactly the same sets of row and column attributes.
	"""
	if file_attrs is None:
		file_attrs = {}

	if len(files) == 0:
		raise ValueError("The input file list was empty")

	copyfile(files[0], output_file)

	ds = connect(output_file)
	for a in file_attrs:
		ds.attrs[a] = file_attrs[a]

	if len(files) >= 2:
		for f in files[1:]:
			ds.add_loom(f, key)
	ds.close()


def connect(filename: str, mode: str = 'r+') -> LoomConnection:
	"""
	Establish a connection to a .loom file.

	Args:
		filename (str):     Name of the .loom file to open
		mode (str):         read/write mode, accepts 'r+' (read/write) or
							'r' (read-only), defaults to 'r+'

	Returns:
		A LoomConnection instance.
	"""
	return LoomConnection(filename, mode)
