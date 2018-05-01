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
import scipy.sparse
from shutil import copyfile
import logging
import time
import loompy
from loompy import deprecated, timestamp


class LoomConnection:
	@property
	def mode(self):
		return self._file.mode

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
		"""

		# make sure a valid mode was passed, if not default to read-only
		# because you probably are doing something that you don't want to
		if mode != 'r+' and mode != 'r':
			raise ValueError("Mode must be either 'r' or 'r+'")
		self.filename = filename
		self._file = h5py.File(filename, mode)
		self._closed = False
		try:
			if "matrix" in self._file:
				self.shape = self._file["/matrix"].shape
			else:
				self.shape = (0, 0)
			self.layers = loompy.LayerManager(self)
			self.view = loompy.ViewManager(self)
			self.ra = loompy.AttributeManager(self, axis=0)
			self.ca = loompy.AttributeManager(self, axis=1)
			self.attrs = loompy.FileAttributeManager(self._file)
			self.row_graphs = loompy.GraphManager(self, axis=0)
			self.col_graphs = loompy.GraphManager(self, axis=1)

			# Compatibility
			self.layer = self.layers
			self.row_attrs = self.ra
			self.col_attrs = self.ca

		except Exception as e:
			logging.warn("initialising LoomConnection to %s failed, closing file connection", filename)
			self.close()
			raise e

	def last_modified(self) -> str:
		"""
		Return an ISO8601 timestamp when the file was last modified

		Note: if the file has no timestamp, and mode is 'r+', a new timestamp is created and returned.
		Otherwise, the current time in UTC is returned
		"""
		if "last_modified" in self._file.attrs:
			return self._file.attrs["last_modified"]
		elif self._file.mode == "r+":
			# Make sure the file has modification timestamps
			self._file.attrs["last_modified"] = timestamp()
			self._file.flush()
			return self._file.attrs["last_modified"]
		return timestamp()

	def get_changes_since(self, timestamp: str) -> Dict[str, List]:
		rg = []
		cg = []
		ra = []
		ca = []
		layers = []

		if self.last_modified() > timestamp:
			if self.row_graphs.last_modified() > timestamp:
				for name in self.row_graphs.keys():
					if self.row_graphs.last_modified(name) > timestamp:
						rg.append(name)
			if self.col_graphs.last_modified() > timestamp:
				for name in self.col_graphs.keys():
					if self.col_graphs.last_modified(name) > timestamp:
						cg.append(name)
			if self.ra.last_modified() > timestamp:
				for name in self.ra.keys():
					if self.ra.last_modified(name) > timestamp:
						ra.append(name)
			if self.ca.last_modified() > timestamp:
				for name in self.ca.keys():
					if self.ca.last_modified(name) > timestamp:
						ca.append(name)
			if self.layers.last_modified() > timestamp:
				for name in self.layers.keys():
					if self.layers.last_modified(name) > timestamp:
						layers.append(name)
		return {"row_graphs": rg, "col_graphs": cg, "row_attrs": ra, "col_attrs": ca, "layers": layers}

	def __enter__(self) -> Any:
		"""
		Context manager, to support "with loompy.connect(..)" construct
		"""
		return self

	def __exit__(self, type: Any, value: Any, traceback: Any) -> None:
		"""
		Context manager, to support "with loompy.connect(..)" construct
		"""
		if not self.closed:
			self.close(True)

	def _repr_html_(self) -> str:
		"""
		Return an HTML representation of the loom file, showing the upper-left 10x10 corner.
		"""
		if not self.closed:
			return loompy.to_html(self)
		else:
			return "This LoomConnection has been closed"

	def __getitem__(self, slice_: Any) -> np.ndarray:
		"""
		Get a slice of the main matrix.

		Args:
			slice:		A slice object (see http://docs.h5py.org/en/latest/high/dataset.html), or np.ndarray, or int

		Returns:
			A numpy ndarray matrix
		"""
		if type(slice_) is str:
			return self.layers[slice_]
		if type(slice_) is not tuple:
			raise ValueError("Slice must be a 2-tuple")
		return self.layers[""][slice_]

	def __setitem__(self, slice_: Any, data: np.ndarray) -> None:
		"""
		Assign a slice of the main matrix.

		Args:
			slice_:		A slice object (see http://docs.h5py.org/en/latest/high/dataset.html), or np.ndarray, or int
			data:		A matrix corresponding to the slice, of the same datatype as the main matrix

		Returns:
			Nothing.
		"""
		if type(slice_) is str:
			self.layers[slice_] = data
		else:
			self.layers[""][slice_] = data

	def sparse(self, rows: np.ndarray = None, cols: np.ndarray = None, layer: str = None) -> scipy.sparse.coo_matrix:
		"""
		Return the main matrix as a scipy.sparse.coo_matrix, without loading dense matrix in RAM

		Args:
			rows:		Rows to include, or None to include all
			cols:		Columns to include, or None to include all
			layer:		Layer to return, or None to return the default layer

		Returns:
			scipy.sparse.coo_matrix
		"""
		if layer is None:
			return self.layers[""].sparse(rows=rows, cols=cols)
		else:
			return self.layers[layer].sparse(rows=rows, cols=cols)

	def close(self, suppress_warning: bool = False) -> None:
		"""
		Close the connection. After this, the connection object becomes invalid. Warns user if called after closing.

		Args:
			suppress_warning:		Suppresses warning message if True (defaults to false)
		"""
		if self._file is None:
			if not suppress_warning:
				# Warn user that they're being paranoid
				# and should clean up their code
				logging.warn("Connection to %s is already closed", self.filename)
		else:
			self._file.close()
			self._file = None
		self.layers = None
		self.ra = None
		self.row_attrs = None
		self.ca = None
		self.col_attrs = None
		self.row_graphs = None
		self.col_graphs = None
		self.shape = (0, 0)
		self._closed = True

	@property
	def closed(self) -> bool:
		return self._closed

	def set_layer(self, name: str, matrix: np.ndarray, chunks: Tuple[int, int] = (64, 64), chunk_cache: int = 512, dtype: str = "float32", compression_opts: int = 2) -> None:
		"""
		**DEPRECATED** - Use `ds.layer.Name = matrix` or `ds.layer[`Name`] = matrix` instead
		"""
		deprecated("'set_layer' is deprecated. Use 'ds.layer.Name = matrix' or 'ds.layer['Name'] = matrix' instead")
		self.layers[name] = matrix

	def add_columns(self, layers: Union[np.ndarray, Dict[str, np.ndarray], loompy.LayerManager], col_attrs: Dict[str, np.ndarray], *, fill_values: Dict[str, np.ndarray] = None) -> None:
		"""
		Add columns of data and attribute values to the dataset.

		Args:
			layers (dict or numpy.ndarray or LayerManager):
				Either:
				1) A N-by-M matrix of float32s (N rows, M columns) in this case columns are added at the default layer
				2) A dict {layer_name : matrix} specified so that the matrix (N, M) will be added to layer `layer_name`
				3) A LayerManager object (such as what is returned by view.layers)
			col_attrs (dict):
				Column attributes, where keys are attribute names and values are numpy arrays (float or string) of length M
			fill_values: dictionary of values to use if a column attribute is missing, or "auto" to fill with zeros or empty strings

		Returns:
			Nothing.

		Notes
		-----
		- This will modify the underlying HDF5 file, which will interfere with any concurrent readers.
		- Column attributes in the file that are NOT provided, will be deleted (unless fill value provided).
		- Array with Nan should not be provided

		"""
		if self._file.mode != "r+":
			raise IOError("Cannot add columns when connected in read-only mode")

		layers_dict: Dict[str, np.ndarray] = None
		if isinstance(layers, np.ndarray):
			layers_dict = {"": layers}
		elif isinstance(layers, loompy.LayerManager):
			layers_dict = {k: v[:, :] for k, v in layers.items()}
		elif isinstance(layers, dict):
			layers_dict = layers
		else:
			raise ValueError("Invalid type for layers argument")
		n_cols = 0
		for layer, matrix in layers_dict.items():
			if layer not in self.layers.keys():
				raise ValueError(f"Layer {layer} does not exist in the target loom file")
			if matrix.shape[0] != self.shape[0]:
				raise ValueError(f"Layer {layer} has {matrix.shape[0]} rows but file has {self.shape[0]}")
			if n_cols == 0:
				n_cols = matrix.shape[1]
			elif matrix.shape[1] != n_cols:
				raise ValueError(f"Layer {layer} has {matrix.shape[1]} columns but the first layer had {n_cols}")
		for layer in self.layers.keys():
			if layer not in layers_dict.keys():
				raise ValueError(f"Layer {layer} does not exist in the target loom file")

		did_remove = False
		todel = []  # type: List[str]
		for key, vals in col_attrs.items():
			if key not in self.col_attrs:
				if fill_values is not None:
					if fill_values == "auto":
						fill_with = np.zeros(1, dtype=col_attrs[key].dtype)[0]
					else:
						fill_with = fill_values[key]
					self.ca[key] = np.array([fill_with] * self.shape[1])
				else:
					did_remove = True
					todel.append(key)
			if len(vals) != n_cols:
				raise ValueError(f"Each column attribute must have exactly {n_cols} values, but {key} had {len(vals)}")
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
					col_attrs[key] = np.array([fill_with] * n_cols)
				else:
					did_remove = True
					todel.append(key)
		for key in todel:
			del self.ca[key]  # delete_attr(key, axis=1)
		if did_remove:
			logging.warn("Some column attributes were removed: " + ",".join(todel))

		n_cols = n_cols + self.shape[1]
		old_n_cols = self.shape[1]
		# Must set new shape here, otherwise the attribute manager will complain
		self.shape = (self.shape[0], n_cols)
		for key, vals in col_attrs.items():
			self.ca[key] = np.concatenate([self.ca[key], vals])

		# Add the columns layerwise
		for key in self.layers.keys():
			self.layers[key].resize(n_cols, axis=1)
			self.layers[key][:, old_n_cols:n_cols] = layers_dict[key]
		self._file.flush()

	def add_loom(self, other_file: str, key: str = None, fill_values: Dict[str, np.ndarray]=None, batch_size: int=1000, convert_attrs: bool=False) -> None:
		"""
		Add the content of another loom file

		Args:
			other_file (str):       filename of the loom file to append
			key:                    Primary key to use to align rows in the other file with this file
			fill_values (dict):     default values to use for missing attributes (or None to drop missing attrs, or 'auto' to fill with sensible defaults)
			batch_size (int):       the batch size used by batchscan (limits the number of rows/columns read in memory)
			convert_attrs (bool):   convert file attributes that differ between files into column attributes

		Returns:
			Nothing, but adds the loom file. Note that the other loom file must have exactly the same
			number of rows, and must have exactly the same column attributes.
			The all the contents including layers but ignores layers in `other_file` that are not already persent in self
		"""
		if self._file.mode != "r+":
			raise IOError("Cannot add data when connected in read-only mode")
		# Connect to the loom files
		other = connect(other_file)
		# Verify that the row keys can be aligned
		ordering = None
		if key is not None:
			# This was original Sten's version but it creates a 400M entries array in memory
			# ordering = np.where(other.row_attrs[key][None, :] == self.row_attrs[key][:, None])[1]

			def ixs_thatsort_a2b(a: np.ndarray, b: np.ndarray, check_content: bool=True) -> np.ndarray:
				"This is super duper magic sauce to make the order of one list to be like another"
				if check_content:
					assert len(np.intersect1d(a, b)) == len(a), "The two arrays are not matching"
				return np.argsort(a)[np.argsort(np.argsort(b))]

			ordering = ixs_thatsort_a2b(a=other.row_attrs[key], b=self.row_attrs[key])
			pk1 = sorted(other.row_attrs[key])
			pk2 = sorted(self.row_attrs[key])
			for ix, val in enumerate(pk1):
				if pk2[ix] != val:
					raise ValueError("Primary keys are not 1-to-1 alignable!")
		diff_layers = set(self.layers.keys()) - set(other.layers.keys())
		if len(diff_layers) > 0:
			raise ValueError("%s is missing a layer, cannot merge with current file. layers missing:%s" % (other_file, diff_layers))

		if convert_attrs:
			# Prepare to replace any global attribute that differ between looms or is missing in either loom with a column attribute.
			globalkeys = set(self.attrs)
			globalkeys.update(other.attrs)
			for globalkey in globalkeys:
				if globalkey in self.attrs and globalkey in other.attrs and self.attrs[globalkey] == other.attrs[globalkey]:
					continue
				if globalkey not in self.col_attrs:
					self_value = self.attrs[globalkey] if globalkey in self.attrs else np.zeros(1, dtype=other.attrs[globalkey].dtype)[0]
					self.col_attrs[globalkey] = np.array([self_value] * self.shape[1])
				if globalkey not in other.col_attrs:
					other_value = other.attrs[globalkey] if globalkey in other.attrs else np.zeros(1, dtype=self.attrs[globalkey].dtype)[0]
					other.col_attrs[globalkey] = np.array([other_value] * other.shape[1])
				if globalkey in self.attrs:
					delattr(self.attrs, globalkey)

		for (ix, selection, vals) in other.batch_scan_layers(axis=1, layers=self.layers.keys(), batch_size=batch_size):
			ca = {key: v[selection] for key, v in other.col_attrs.items()}
			if ordering is not None:
				vals = {key: val[ordering, :] for key, val in vals.items()}
			self.add_columns(vals, ca, fill_values=fill_values)
		other.close()

	def delete_attr(self, name: str, axis: int = 0) -> None:
		"""
		**DEPRECATED** - Use `del ds.ra.key` or `del ds.ca.key` instead, where `key` is replaced with the attribute name
		"""
		deprecated("'delete_attr' is deprecated. Use 'del ds.ra.key' or 'del ds.ca.key' instead")
		if axis == 0:
			del self.ra[name]
		else:
			del self.ca[name]

	def set_attr(self, name: str, values: np.ndarray, axis: int = 0, dtype: str = None) -> None:
		"""
		**DEPRECATED** - Use `ds.ra.key = values` or `ds.ca.key = values` instead
		"""
		deprecated("'set_attr' is deprecated. Use 'ds.ra.key = values' or 'ds.ca.key = values' instead")
		if axis == 0:
			self.ra[name] = values
		else:
			self.ca[name] = values

	def list_edges(self, *, axis: int) -> List[str]:
		"""
		**DEPRECATED** - Use `ds.row_graphs.keys()` or `ds.col_graphs.keys()` instead
		"""
		deprecated("'list_edges' is deprecated. Use 'ds.row_graphs.keys()' or 'ds.col_graphs.keys()' instead")
		if axis == 0:
			return self.row_graphs.keys()
		elif axis == 1:
			return self.col_graphs.keys()
		else:
			return []

	def get_edges(self, name: str, *, axis: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
		"""
		**DEPRECATED** - Use `ds.row_graphs[name]` or `ds.col_graphs[name]` instead
		"""
		deprecated("'get_edges' is deprecated. Use 'ds.row_graphs[name]' or 'ds.col_graphs[name]' instead")
		if axis == 0:
			g = self.row_graphs[name]
			return (g.row, g.col, g.data)
		if axis == 1:
			g = self.col_graphs[name]
			return (g.row, g.col, g.data)
		raise ValueError("Axis must be 0 or 1")

	def set_edges(self, name: str, a: np.ndarray, b: np.ndarray, w: np.ndarray, *, axis: int) -> None:
		"""
		**DEPRECATED** - Use `ds.row_graphs[name] = g` or `ds.col_graphs[name] = g` instead
		"""
		deprecated("'set_edges' is deprecated. Use 'ds.row_graphs[name] = g' or 'ds.col_graphs[name] = g' instead")
		try:
			g = scipy.sparse.coo_matrix((w, (a, b)), (self.shape[axis], self.shape[axis]))
		except Exception:
			raise ValueError("Input arrays could not be converted to a sparse matrix")
		if axis == 0:
			self.row_graphs[name] = g
		elif axis == 1:
			self.col_graphs[name] = g
		else:
			raise ValueError("axis must be 0 (rows) or 1 (columns)")

	def scan(self, *, items: np.ndarray = None, axis: int = None, layers: Iterable = None, key: str = None, batch_size: int = 8 * 64) -> Iterable[Tuple[int, np.ndarray, loompy.LoomView]]:
		"""
		Scan across one axis and return batches of rows (columns) as LoomView objects

		Args
		----
		items: np.ndarray
			the indexes [0, 2, 13, ... ,973] of the rows/cols to include along the axis
			OR: boolean mask array giving the rows/cols to include
		axis: int
			0:rows or 1:cols
		batch_size: int
			the chuncks returned at every element of the iterator
		layers: iterable
			if specified it will batch scan only across some of the layers of the loom file
			if layers == None, all layers will be scanned
			if layers == [""] or "", only the default layer will be scanned
		key:
			Name of primary key attribute. If specified, return the values sorted by the key

		Returns
		------
		Iterable that yields triplets
		(ix, indexes, view)

		ix: int
			first position / how many rows/cols have been yielded alredy
		indexes: np.ndarray[int]
			the indexes with the same numbering of the input args cells / genes (i.e. np.arange(len(ds.shape[axis])))
			this is ix + selection
		view: LoomView
			a view corresponding to the current chunk
		"""
		if axis is None:
			raise ValueError("Axis must be given (0 = rows, 1 = cols)")
		if layers is None:
			layers = self.layers.keys()
		if layers == "":
			layers = [""]

		if (items is not None) and (np.issubdtype(items.dtype, np.bool_)):
			items = np.where(items)[0]

		ordering: np.ndarray = None
		vals: Dict[str, loompy.MemoryLoomLayer] = {}
		if axis == 1:
			if key is not None:
				ordering = np.argsort(self.ra[key])
			else:
				ordering = np.arange(self.shape[0])
			if items is None:
				items = np.fromiter(range(self.shape[1]), dtype='int')
			cols_per_chunk = batch_size
			ix = 0
			while ix < self.shape[1]:
				cols_per_chunk = min(self.shape[1] - ix, cols_per_chunk)
				selection = items - ix
				# Pick out the cells that are in this batch
				selection = selection[np.where(np.logical_and(selection >= 0, selection < cols_per_chunk))[0]]
				if selection.shape[0] == 0:
					ix += cols_per_chunk
					continue

				# Load the whole chunk from the file, then extract genes and cells using fancy indexing
				for layer in layers:
					temp = self.layers[layer][:, ix:ix + cols_per_chunk]
					temp = temp[ordering, :]
					temp = temp[:, selection]
					vals[layer] = loompy.MemoryLoomLayer(layer, temp)
				lm = loompy.LayerManager(None)
				for key, layer in vals.items():
					lm[key] = loompy.MemoryLoomLayer(key, layer)
				view = loompy.LoomView(lm, self.ra[ordering], self.ca[ix + selection], self.row_graphs[ordering], self.col_graphs[ix + selection], filename=self.filename, file_attrs=self.attrs)
				yield (ix, ix + selection, view)
				ix += cols_per_chunk
		elif axis == 0:
			if key is not None:
				ordering = np.argsort(self.ca[key])
			else:
				ordering = np.arange(self.shape[1])
			if items is None:
				items = np.fromiter(range(self.shape[0]), dtype='int')
			rows_per_chunk = batch_size
			ix = 0
			while ix < self.shape[0]:
				rows_per_chunk = min(self.shape[0] - ix, rows_per_chunk)
				selection = items - ix
				# Pick out the genes that are in this batch
				selection = selection[np.where(np.logical_and(selection >= 0, selection < rows_per_chunk))[0]]
				if selection.shape[0] == 0:
					ix += rows_per_chunk
					continue

				# Load the whole chunk from the file, then extract genes and cells using fancy indexing
				for layer in layers:
					temp = self.layers[layer][ix:ix + rows_per_chunk, :]
					temp = temp[:, ordering]
					temp = temp[selection, :]
					vals[layer] = loompy.MemoryLoomLayer(layer, temp)
				lm = loompy.LayerManager(None)
				for key, layer in vals.items():
					lm[key] = loompy.MemoryLoomLayer(key, layer)
				view = loompy.LoomView(lm, self.ra[ix + selection], self.ca[ordering], self.row_graphs[ix + selection], self.col_graphs[ordering], filename=self.filename, file_attrs=self.attrs)
				yield (ix, ix + selection, view)
				ix += rows_per_chunk
		else:
			raise ValueError("axis must be 0 or 1")

	def batch_scan(self, cells: np.ndarray = None, genes: np.ndarray = None, axis: int = 0, batch_size: int = 1000, layer: str = None) -> Iterable[Tuple[int, np.ndarray, np.ndarray]]:
		"""
		**DEPRECATED** - Use `scan` instead
		"""
		deprecated("'batch_scan' is deprecated. Use 'scan' instead")
		if cells is None:
			cells = np.fromiter(range(self.shape[1]), dtype='int')
		if genes is None:
			genes = np.fromiter(range(self.shape[0]), dtype='int')
		if layer is None:
			layer = ""
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
				vals = self.layers[layer][:, ix:ix + cols_per_chunk]
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
				vals = self.layers[layer][ix:ix + rows_per_chunk, :]
				vals = vals[selection, :]
				vals = vals[:, cells]
				yield (ix, ix + selection, vals)
				ix += rows_per_chunk

	def batch_scan_layers(self, cells: np.ndarray = None, genes: np.ndarray = None, axis: int = 0, batch_size: int = 1000, layers: Iterable = None) -> Iterable[Tuple[int, np.ndarray, Dict]]:
		"""
		**DEPRECATED** - Use `scan` instead
		"""
		deprecated("'batch_scan_layers' is deprecated. Use 'scan' instead")
		if cells is None:
			cells = np.fromiter(range(self.shape[1]), dtype='int')
		if genes is None:
			genes = np.fromiter(range(self.shape[0]), dtype='int')
		if layers is None:
			layers = self.layers.keys()
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
					vals[key] = self.layers[key][:, ix:ix + cols_per_chunk]
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
					vals[key] = self.layers[key][ix:ix + rows_per_chunk, :]
					vals[key] = vals[key][selection, :]
					vals[key] = vals[key][:, cells]
				yield (ix, ix + selection, vals)
				ix += rows_per_chunk

	def map(self, f_list: List[Callable[[np.ndarray], int]], *, axis: int = 0, chunksize: int = 1000, selection: np.ndarray = None) -> List[np.ndarray]:
		"""
		Apply a function along an axis without loading the entire dataset in memory.

		Args:
			f (list of func):		Function(s) that takes a numpy ndarray as argument

			axis (int):		Axis along which to apply the function (0 = rows, 1 = columns)

			chunksize (int): Number of rows (columns) to load per chunk

			selection (array of bool): Columns (rows) to include

		Returns:
			numpy.ndarray result of function application
			The result is a list of numpy arrays, one per supplied function in f_list.
			This is more efficient than repeatedly calling map() one function at a time.
		"""
		return self.layers[""].map(f_list, axis, chunksize, selection)

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
		self.layers.permute(ordering, axis=axis)
		if axis == 0:
			self.row_attrs.permute(ordering)
			self.row_graphs.permute(ordering)
		if axis == 1:
			self.col_attrs.permute(ordering)
			self.col_graphs.permute(ordering)

	def export(self, out_file: str, layer: str = None, format: str = "tab") -> None:
		"""
		Export the specified layer and row/col attributes as tab-delimited file
		"""
		if format != "tab":
			raise NotImplementedError("Only 'tab' is supported")

		with open(out_file, "w") as f:
			# Emit column attributes
			for ca in self.col_attrs.keys():
				for ra in self.row_attrs.keys():
					f.write("\t")
				f.write(ca + "\t")
				for v in self.col_attrs[ca]:
					f.write(str(v) + "\t")
				f.write("\n")

			# Emit row attribute names
			for ra in self.row_attrs.keys():
				f.write(ra + "\t")
			f.write("\t")
			for v in range(self.shape[1]):
				f.write("\t")
			f.write("\n")

			# Emit row attr values and matrix values
			for row in range(self.shape[0]):
				for ra in self.row_attrs.keys():
					f.write(str(self.row_attrs[ra][row]) + "\t")
				f.write("\t")

				if layer is None:
					for v in self[row, :]:
						f.write(str(v) + "\t")
				else:
					for v in self.layers[layer][row, :]:
						f.write(str(v) + "\t")
				f.write("\n")


def _create_sparse(filename: str, matrix: np.ndarray, row_attrs: Dict[str, np.ndarray], col_attrs: Dict[str, np.ndarray], *, layers: Dict[str, np.ndarray] = None, file_attrs: Dict[str, str] = None) -> None:
	"""
	Create a new loom file from a sparse matrix input
	"""
	if os.path.exists(filename):
		raise FileExistsError("Cannot overwrite existing file " + filename)
	logging.info("Converting to csc format")
	matrix = matrix.tocsc()
	window = 6400
	ix = 0
	while ix < matrix.shape[1]:
		window = min(window, matrix.shape[1] - ix)
		if window == 0:
			break
		ca = {key: val[ix:ix + window] for (key, val) in col_attrs.items()}
		create_append(filename, matrix[:, ix:ix + window].toarray(), row_attrs, ca, file_attrs=file_attrs)
		ix += window


def create_append(filename: str, layers: Union[np.ndarray, Dict[str, np.ndarray], loompy.LayerManager], row_attrs: Dict[str, np.ndarray], col_attrs: Dict[str, np.ndarray], *, file_attrs: Dict[str, str] = None, fill_values: Dict[str, np.ndarray] = None) -> None:
	"""
	Append columns to a loom file, or create a new loom file if it doesn't exist

	Args:
		filename (str):         The filename (typically using a `.loom` file extension)
		layers (np.ndarray or Dict[str, np.ndarray] or LayerManager):
								Two-dimensional (N-by-M) numpy ndarray of float values
								Or dictionary of named layers, each an N-by-M ndarray
								or LayerManager, each layer an N-by-M ndarray
		row_attrs (dict):       Row attributes, where keys are attribute names and values
								are numpy arrays (float or string) of length N
		col_attrs (dict):       Column attributes, where keys are attribute names and
								values are numpy arrays (float or string) of length M
		file_attrs (dict):      Global attributes, where keys are attribute names and
								values are strings
	Returns:
		Nothing
	"""
	if os.path.exists(filename):
		with connect(filename) as ds:
			ds.add_columns(layers, col_attrs, fill_values=fill_values)
	else:
		create(filename, layers, row_attrs, col_attrs, file_attrs=file_attrs)


def create(filename: str, layers: Union[np.ndarray, Dict[str, np.ndarray], loompy.LayerManager], row_attrs: Dict[str, np.ndarray], col_attrs: Dict[str, np.ndarray], *, file_attrs: Dict[str, str] = None) -> None:
	"""
	Create a new .loom file from the given data.

	Args:
		filename (str):         The filename (typically using a `.loom` file extension)
		layers (np.ndarray or scipy.sparse or Dict[str, np.ndarray] or LayerManager):
								Two-dimensional (N-by-M) numpy ndarray of float values
								Or sparse matrix
								Or dictionary of named layers, each an N-by-M ndarray
								or LayerManager, each layer an N-by-M ndarray
		row_attrs (dict):       Row attributes, where keys are attribute names and values
								are numpy arrays (float or string) of length N
		col_attrs (dict):       Column attributes, where keys are attribute names and
								values are numpy arrays (float or string) of length M
		file_attrs (dict):      Global attributes, where keys are attribute names and
								values are strings
	Returns:
		Nothing

	Remarks:
		If the file exists, it will be overwritten. See create_append for a function that will append to existing files.
	"""
	if filename.startswith("~/"):
		filename = os.path.expanduser(filename)
	if file_attrs is None:
		file_attrs = {}

	if isinstance(layers, np.ndarray):
		layers = {"": layers}
	elif scipy.sparse.issparse(layers):
		_create_sparse(filename, layers, row_attrs, col_attrs, file_attrs=file_attrs)
		return
	elif isinstance(layers, loompy.LayerManager):
		layers = {k: v[:, :] for k, v in layers.items()}
	if "" not in layers:
		raise ValueError("Data for default layer must be provided")

	# Create the file (empty).
	# Yes, this might cause an exception, which we prefer to send to the caller
	f = h5py.File(name=filename, mode='w')
	f.create_group('/layers')
	f.create_group('/row_attrs')
	f.create_group('/col_attrs')
	f.create_group('/row_graphs')
	f.create_group('/col_graphs')
	f.flush()
	f.close()

	try:
		with connect(filename) as ds:
			for key, vals in layers.items():
				ds.layer[key] = vals

			for key, vals in row_attrs.items():
				ds.ra[key] = vals

			for key, vals in col_attrs.items():
				ds.ca[key] = vals

			for vals in file_attrs:
				ds.attrs[vals] = file_attrs[vals]

			# store creation date
			currentTime = time.localtime(time.time())
			ds.attrs['CreationDate'] = timestamp()
			ds.attrs["LOOM_SPEC_VERSION"] = loompy.loom_spec_version

	except ValueError as ve:
		ds.close(suppress_warning=True)
		os.remove(filename)
		raise ve


def create_from_cellranger(indir: str, outdir: str = None, genome: str = None) -> None:
	"""
	Create a .loom file from 10X Genomics cellranger output

	Args:
		indir (str):	path to the cellranger output folder (the one that contains 'outs')
		outdir (str):	output folder wher the new loom file should be saved (default to indir)
		genome (str):	genome build to load (e.g. 'mm10'; if None, determine species from outs folder)

	Returns:
		LoomConnection to created loom file.
	"""
	if outdir is None:
		outdir = indir
	sampleid = os.path.split(os.path.abspath(indir))[-1]
	matrix_folder = os.path.join(indir, 'outs', 'filtered_gene_bc_matrices')
	if genome is None:
		genome = [f for f in os.listdir(matrix_folder) if not f.startswith(".")][0]
	matrix_folder = os.path.join(matrix_folder, genome)
	matrix = mmread(os.path.join(matrix_folder, "matrix.mtx")).astype("float32").todense()

	with open(os.path.join(matrix_folder, "genes.tsv"), "r") as f:
		lines = f.readlines()
	accession = np.array([x.split("\t")[0] for x in lines]).astype("str")
	gene = np.array([x.split("\t")[1].strip() for x in lines]).astype("str")
	with open(os.path.join(matrix_folder, "barcodes.tsv"), "r") as f:
		lines = f.readlines()
	cellids = np.array([sampleid + ":" + x.strip() for x in lines]).astype("str")

	col_attrs = {"CellID": cellids}
	row_attrs = {"Accession": accession, "Gene": gene}

	tsne_file = os.path.join(indir, "outs", "analysis", "tsne", "projection.csv")
	# In cellranger V2 the file moved one level deeper
	if not os.path.exists(tsne_file):
		tsne_file = os.path.join(indir, "outs", "analysis", "tsne", "2_components", "projection.csv")
	if os.path.exists(tsne_file):
		tsne = np.loadtxt(tsne_file, usecols=(1, 2), delimiter=',', skiprows=1)
		col_attrs["X"] = tsne[:, 0].astype('float32')
		col_attrs["Y"] = tsne[:, 1].astype('float32')

	clusters_file = os.path.join(indir, "outs", "analysis", "clustering", "graphclust", "clusters.csv")
	if os.path.exists(clusters_file):
		labels = np.loadtxt(clusters_file, usecols=(1, ), delimiter=',', skiprows=1)
		col_attrs["ClusterID"] = labels.astype('int') - 1

	create(os.path.join(outdir, sampleid + ".loom"), matrix, row_attrs, col_attrs, file_attrs={"Genome": genome})


def combine(files: List[str], output_file: str, key: str = None, file_attrs: Dict[str, str]=None, batch_size: int=1000, convert_attrs: bool=False) -> None:
	"""
	Combine two or more loom files and save as a new loom file

	Args:
		files (list of str):    the list of input files (full paths)
		output_file (str):      full path of the output loom file
		key (string):           Row attribute to use to verify row ordering
		file_attrs (dict):      file attributes (title, description, url, etc.)
		batch_size (int):       limits the batch or cols/rows read in memory (default: 1000)
		convert_attrs (bool):   convert file attributes that differ between files into column attributes

	Returns:
		Nothing, but creates a new loom file combining the input files.

	The input files must (1) have exactly the same number of rows, (2) have
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
			ds.add_loom(f, key, batch_size=batch_size, convert_attrs=convert_attrs)
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

	Remarks:
		This function should typically be called as a context manager:

			with loompy.connect(filename) as ds:
				...do something...

		This ensures that the file will be closed automatically when the context block ends
	"""
	return LoomConnection(filename, mode)
