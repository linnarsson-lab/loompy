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

import gzip
import logging
import os.path
import re
from shutil import copyfile
from typing import Tuple, Union, Any, Dict, List, Iterable, Callable, Optional

import h5py
import numpy as np
import numpy_groupies.aggregate_numpy as npg
import scipy.sparse
from scipy.io import mmread
import scipy.sparse as sparse

import loompy
from loompy import deprecated, timestamp
from loompy.metadata_loaders import make_row_attrs_from_gene_metadata, load_sample_metadata
from .cell_calling import call_cells


class LoomConnection:
	'''
	A connection to a Loom file on disk. Typically LoomConnection objects are created using one of the
	functions on the loompy module, such as :func:`loompy.connect` or :func:`loompy.new`. LoomConnection
	objects are context managers and should normally be
	wrapped in a ``with`` block:

	.. highlight:: python
	.. code-block:: python

		import loompy
		with loompy.connect("mydata.loom") as ds:
			print(ds.ca.keys())

	Inside the ``with`` block, you can access the dataset (here using the variable ``ds``). When execution
	leaves the ``with`` block, the connection is automatically closed, freeing up resources.
	'''
	def __init__(self, filename: str, mode: str = 'r+', *, validate: bool = True) -> None:
		"""
		Establish a connection to a Loom file.

		Args:
			filename:			Name of the .loom file to open
			mode:				read/write mode, accepts 'r+' (read/write) or
								'r' (read-only), defaults to 'r+' without arguments,
								and to 'r' with incorrect arguments
			validate:			Validate that the file conforms with the Loom specification
		Returns:
			Nothing.
		"""
		if not os.path.exists(filename):
			raise IOError(f"File '{filename}' not found")
		# make sure a valid mode was passed
		if mode != 'r+' and mode != 'r':
			raise ValueError("Mode must be either 'r' or 'r+'")
		self.filename = filename  #: Path to the file (as given when the LoomConnection was created)

		# Validate the file
		if validate:
			lv = loompy.LoomValidator()
			if not lv.validate(filename):
				raise ValueError("\n".join(lv.errors) + f"\n{filename} does not appear to be a valid Loom file according to Loom spec version '{lv.version}'")

		self._file = h5py.File(filename, mode)
		self._closed = False
		if "matrix" in self._file:
			self.shape = self._file["/matrix"].shape  #: Shape of the dataset (n_rows, n_cols)
		else:
			self.shape = (0, 0)
		self.layers = loompy.LayerManager(self)
		self.view = loompy.ViewManager(self)  #: Create a view of the file by slicing this attribute, like ``ds.view[:100, :100]``
		self.ra = loompy.AttributeManager(self, axis=0)  #: Row attributes, dict-like with np.ndarray values
		self.ca = loompy.AttributeManager(self, axis=1)  #: Column attributes, dict-like with np.ndarray values
		self.attrs = loompy.GlobalAttributeManager(self._file)  #: Global attributes
		self.row_graphs = loompy.GraphManager(self, axis=0)  #: Row graphs, dict-like with values that are :class:`scipy.sparse.coo_matrix` objects
		self.col_graphs = loompy.GraphManager(self, axis=1)  #: Column graphs, dict-like with values that are :class:`scipy.sparse.coo_matrix` objects

		# Compatibility
		self.layer = self.layers
		self.row_attrs = self.ra
		self.col_attrs = self.ca

	@property
	def mode(self) -> str:
		"""
		The access mode of the connection ('r' or 'r+')
		"""
		return self._file.mode

	def last_modified(self) -> str:
		"""
		Return an ISO8601 timestamp indicating when the file was last modified

		Returns:
			An ISO8601 timestamp indicating when the file was last modified

		Remarks:
			If the file has no timestamp, and mode is 'r+', a new timestamp is created and returned.
			Otherwise, the current time in UTC is returned
		"""
		if "last_modified" in self.attrs:
			return self.attrs["last_modified"]
		elif self.mode == "r+":
			# Make sure the file has modification timestamps
			self.attrs["last_modified"] = timestamp()
			return self.attrs["last_modified"]
		return timestamp()

	def get_changes_since(self, timestamp: str) -> Dict[str, List]:
		"""
		Get a summary of the parts of the file that changed since the given time

		Args:
			timestamp:	ISO8601 timestamp

		Return:
			dict:	Dictionary like ``{"row_graphs": rg, "col_graphs": cg, "row_attrs": ra, "col_attrs": ca, "layers": layers}`` listing the names of objects that were modified since the given time
		"""
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
		if self.shape[0] == 0 or self.shape[1] == 0:
			raise ValueError("Newly created loom file must be filled with data before leaving the 'with' statement")
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
		Return the main matrix or specified layer as a scipy.sparse.coo_matrix, without loading dense matrix in RAM

		Args:
			rows:		Rows to include, or None to include all
			cols:		Columns to include, or None to include all
			layer:		Layer to return, or None to return the default layer

		Returns:
			Sparse matrix (:class:`scipy.sparse.coo_matrix`)
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
		self.layers = None  # type: ignore
		self.ra = None  # type: ignore
		self.row_attrs = None  # type: ignore
		self.ca = None  # type: ignore
		self.col_attrs = None  # type: ignore
		self.row_graphs = None  # type: ignore
		self.col_graphs = None  # type: ignore
		self.shape = (0, 0)
		self._closed = True

	@property
	def closed(self) -> bool:
		"""
		True if the connection is closed.
		"""
		return self._closed

	def set_layer(self, name: str, matrix: np.ndarray, chunks: Tuple[int, int] = (64, 64), chunk_cache: int = 512, dtype: str = "float32", compression_opts: int = 2) -> None:
		"""
		**DEPRECATED** - Use `ds.layer.Name = matrix` or `ds.layer[`Name`] = matrix` instead
		"""
		deprecated("'set_layer' is deprecated. Use 'ds.layer.Name = matrix' or 'ds.layer['Name'] = matrix' instead")
		self.layers[name] = matrix

	def add_columns(self, layers: Union[np.ndarray, Dict[str, np.ndarray], loompy.LayerManager], col_attrs: Dict[str, np.ndarray], *, row_attrs: Dict[str, np.ndarray] = None, fill_values: Dict[str, np.ndarray] = None) -> None:
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
			row_attrs (dict):
				Optional row attributes, where keys are attribute names and values are numpy arrays (float or string) of length M
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

		# If this is an empty loom file, just assign the provided row and column attributes, and set the shape
		is_new = self.shape == (0, 0)
		if is_new:
			if row_attrs is None:
				raise ValueError("row_attrs must be provided when adding to an empty (new) Loom file")
			for k, v in row_attrs.items():
				self.ra[k] = v
				self.shape = (self.ra[k].shape[0], self.shape[1])
			if len(self.ca) == 0:
				for k, v in col_attrs.items():
					v = np.array(v)
					self.ca[k] = np.zeros(0, v.dtype)

		layers_dict: Dict[str, np.ndarray] = {}
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
			if not is_new and layer not in self.layers.keys():
				raise ValueError(f"Layer {layer} does not exist in the target loom file")
			if matrix.shape[0] != self.shape[0]:
				raise ValueError(f"Layer {layer} has {matrix.shape[0]} rows but file has {self.shape[0]}")
			if n_cols == 0:
				n_cols = matrix.shape[1]
			elif matrix.shape[1] != n_cols:
				raise ValueError(f"Layer {layer} has {matrix.shape[1]} columns but the first layer had {n_cols}")

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
			logging.debug("Some column attributes were removed: " + ",".join(todel))

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
			logging.debug("Some column attributes were removed: " + ",".join(todel))

		if is_new:
			for k, v in layers_dict.items():
				self.layers[k] = v
			for k, v in col_attrs.items():
				self.ca[k] = v
		else:
			n_cols = n_cols + self.shape[1]
			old_n_cols = self.shape[1]
			# Must set new shape here, otherwise the attribute manager will complain
			self.shape = (self.shape[0], n_cols)
			todel = []
			for key, vals in col_attrs.items():
				vals = np.array(vals)
				if vals.shape[1:] != self.col_attrs[key].shape[1:]:
					logging.debug(f"Removing attribute {key} because shape {vals.shape} did not match existing shape {self.col_attrs[key].shape} beyond first dimension")
					todel.append(key)
				else:
					self.ca[key] = np.concatenate([self.ca[key], vals])
			for key in todel:
				del self.ca[key]

			# Add the columns layerwise
			for key in self.layers.keys():
				self.layers[key]._resize(n_cols, axis=1)
				self.layers[key][:, old_n_cols:n_cols] = layers_dict[key]
		self._file.flush()

	def add_loom(self, other_file: str, key: str = None, fill_values: Dict[str, np.ndarray] = None, batch_size: int = 1000, convert_attrs: bool = False, include_graphs: bool = False) -> None:
		"""
		Add the content of another loom file

		Args:
			other_file:       filename of the loom file to append
			key:                    Primary key to use to align rows in the other file with this file
			fill_values:     default values to use for missing attributes (or None to drop missing attrs, or 'auto' to fill with sensible defaults)
			batch_size:       the batch size used by batchscan (limits the number of rows/columns read in memory)
			convert_attrs:   convert file attributes that differ between files into column attributes
			include_graphs:  if true, include all the column graphs from other_file that are also present in this file

		Returns:
			Nothing, but adds the loom file. Note that the other loom file must have exactly the same
			number of rows, and must have exactly the same column attributes.
			Adds all the contents including layers but ignores layers in `other_file` that are not already present in self
			Note that graphs are normally not added, unless include_graphs == True, in which case column graphs are added
		"""
		if self._file.mode != "r+":
			raise IOError("Cannot add data when connected in read-only mode")
		# Connect to the loom files
		with loompy.connect(other_file) as other:
			# Verify that the row keys can be aligned
			ordering = None
			if key is not None:
				# This was original Sten's version but it creates a 400M entries array in memory
				# ordering = np.where(other.row_attrs[key][None, :] == self.row_attrs[key][:, None])[1]

				def ixs_thatsort_a2b(a: np.ndarray, b: np.ndarray, check_content: bool = True) -> np.ndarray:
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
				self.add_columns(vals, ca, fill_values=fill_values)  # type: ignore
			
			if include_graphs:
				for gname in self.col_graphs.keys():
					if gname in other.col_graphs:
						g1 = self.col_graphs[gname]
						g2 = other.col_graphs[gname]
						n = self.shape[1]
						m = other.shape[1]
						g3 = scipy.sparse.coo_matrix((np.concatenate([g1.data, g2.data]), (np.concatenate([g1.row, g2.row + n]), np.concatenate([g1.col, g2.col + n]))), shape=(n + m, n + m))
						self.col_graphs[gname] = g3

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

	def scan(self, *, items: np.ndarray = None, axis: int = None, layers: Iterable = None, key: str = None, batch_size: int = 8 * 64, what: List[str] = ["col_attrs", "row_attrs", "layers", "col_graphs", "row_graphs"]) -> Iterable[Tuple[int, np.ndarray, loompy.LoomView]]:
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
		Iterable that yields triplets of (ix, indexes, view) where

		ix: int
			first position / how many rows/cols have been yielded alredy
		indexes: np.ndarray[int]
			the indexes with the same numbering of the input args cells / genes (i.e. ``np.arange(len(ds.shape[axis]))``)
			this is ``ix + selection``
		view: LoomView
			a view corresponding to the current chunk
		"""
		if "layers" not in what:
			raise ValueError("Layers must be included in 'what' parameter (but you can select specific layers using 'layers')")
		if axis is None:
			raise ValueError("Axis must be given (0 = rows, 1 = cols)")
		if layers is None:
			layers = self.layers.keys()
		if layers == "":
			layers = [""]

		if (items is not None) and (type(items) != int) and (np.issubdtype(items.dtype, np.bool_)):
			items = np.where(items)[0]

		ordering: Union[np.ndarray, slice] = None
		vals: Dict[str, loompy.MemoryLoomLayer] = {}
		if axis == 1:
			if key is not None:
				ordering = np.argsort(self.ra[key])
			else:
				# keep everything in original order
				ordering = slice(None)
			if items is None:
				items = np.fromiter(range(self.shape[1]), dtype='int')
			elif type(items) == int:
				items = np.fromiter(range(items, self.shape[1]), dtype='int')
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

				if selection.shape[0] == cols_per_chunk:
					selection = None  # Meaning, select all columns

				# Load the whole chunk from the file, then extract genes and cells using fancy indexing
				for layer in layers:
					temp = self.layers[layer][:, ix:ix + cols_per_chunk]
					temp = temp[ordering, :]
					if selection is not None:
						temp = temp[:, selection]
					vals[layer] = loompy.MemoryLoomLayer(layer, temp)
				lm = loompy.LayerManager(None)
				for key, layer in vals.items():
					lm[key] = loompy.MemoryLoomLayer(key, layer)
				ra = self.ra[ordering] if "row_attrs" in what else {}
				if "col_attrs" in what:
					if selection is not None:
						ca = self.ca[ix + selection]
					else:
						ca = self.ca[ix: ix + cols_per_chunk]
				else:
					ca = {}
				rg = self.row_graphs[ordering] if "row_graphs" in what else None
				if "col_graphs" in what:
					if selection is not None:
						cg = self.col_graphs[ix + selection]
					else:
						cg = self.col_graphs[ix: ix + cols_per_chunk]
				else:
					cg = None
				view = loompy.LoomView(lm, ra, ca, rg, cg, filename=self.filename, file_attrs=self.attrs)
				if selection is not None:
					yield (ix, ix + selection, view)
				else:
					yield (ix, ix + np.arange(cols_per_chunk), view)
				ix += cols_per_chunk
		elif axis == 0:
			if key is not None:
				ordering = np.argsort(self.ca[key])
			else:
				# keep everything in original order
				ordering = slice(None)
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

				if selection.shape[0] == rows_per_chunk:
					selection = None  # Meaning, select all rows

				# Load the whole chunk from the file, then extract genes and cells using fancy indexing
				for layer in layers:
					temp = self.layers[layer][ix:ix + rows_per_chunk, :]
					temp = temp[:, ordering]
					if selection is not None:
						temp = temp[selection, :]
					vals[layer] = loompy.MemoryLoomLayer(layer, temp)
				lm = loompy.LayerManager(None)
				for key, layer in vals.items():
					lm[key] = loompy.MemoryLoomLayer(key, layer)
				
				if "row_attrs" in what:
					if selection is not None:
						ra = self.ra[ix + selection]
					else:
						ra = self.ra[ix: ix + rows_per_chunk]
				else:
					ra = {}
				ca = self.ca[ordering] if "col_attrs" in what else {}
				if "row_graphs" in what:
					if selection is not None:
						rg = self.row_graphs[ix + selection]
					else:
						rg = self.row_graphs[ix: ix + rows_per_chunk]
				else:
					rg = None
				cg = self.col_graphs[ordering] if "col_graphs" in what else None
				view = loompy.LoomView(lm, ra, ca, rg, cg, filename=self.filename, file_attrs=self.attrs)
				if selection is not None:
					yield (ix, ix + selection, view)
				else:
					yield (ix, ix + np.arange(rows_per_chunk), view)
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
			f:		Function(s) that takes a numpy ndarray as argument

			axis:		Axis along which to apply the function (0 = rows, 1 = columns)

			chunksize: Number of rows (columns) to load per chunk

			selection: Columns (rows) to include

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
		self.layers._permute(ordering, axis=axis)
		if axis == 0:
			self.row_attrs._permute(ordering)
			self.row_graphs._permute(ordering)
		if axis == 1:
			self.col_attrs._permute(ordering)
			self.col_graphs._permute(ordering)

	def aggregate(self, out_file: str = None, select: np.ndarray = None, group_by: Union[str, np.ndarray] = "Clusters", \
	                    aggr_by: str = "mean", aggr_ca_by: Dict[str, str] = None, layer: str = "", aggr_ca_if_equal: bool = False) -> np.ndarray:
		"""
		Aggregate the Loom file by applying aggregation functions to the main matrix as well as to the column attributes

		Args:
			out_file	The name of the output Loom file (will be appended to if it exists)
			select		Bool array giving the columns to include (or None, to include all)
			group_by	The column attribute to group by, or an np.ndarray of integer group labels
			aggr_by 	The aggregation function for the main matrix
			aggr_ca_by	A dictionary of aggregation functions for the column attributes (or None to skip)
			layer		The name of the layer to aggregate. Defaults to main layer
			aggr_ca_if_equal If True, scalar column attributes not in aggr_ca_by will be transferred when they have the same value within each aggregate  group

		Returns:
			m			Aggregated main matrix

		Remarks:
			aggr_by gives the aggregation function for the main matrix
			aggr_ca_by is a dictionary with column attributes as keys and aggregation functionas as values
			
			Aggregation functions can be any valid aggregation function from here: https://github.com/ml31415/numpy-groupies

			In addition, you can specify:
				"tally" to count the number of occurences of each value of a categorical attribute
		
		"""
		ca = {}  # type: Dict[str, np.ndarray]
		if select is not None:
			raise ValueError("The 'select' argument is deprecated")
		if isinstance(group_by, np.ndarray):
			labels = group_by
		else:
			labels = (self.ca[group_by]).astype('int')
		_, zero_strt_sort_noholes_lbls = np.unique(labels, return_inverse=True)
		n_groups = len(set(labels))
		if aggr_ca_if_equal:
			for key in self.ca.keys():
				if (aggr_ca_by is None or key not in aggr_ca_by) and np.isscalar(self.ca[key][0]):
					value_by_pos = self.ca[key]
					nvalues_per_lbl = npg.aggregate(zero_strt_sort_noholes_lbls, value_by_pos, lambda v: len(set(v)))
					if np.all(nvalues_per_lbl == 1):
						ca[key] = npg.aggregate(zero_strt_sort_noholes_lbls, self.ca[key], func='first', fill_value=self.ca[key][0])
		if aggr_ca_by is not None:
			for key in self.ca.keys():
				if key not in aggr_ca_by:
					continue
				func = aggr_ca_by[key]
				if func == "tally":
					for val in set(self.ca[key]):
						if np.issubdtype(type(val), np.str_):
							valnew = val.replace("/", "-")  # Slashes are not allowed in attribute names
							valnew = valnew.replace(".", "_")  # Nor are periods
						else:
							valnew = val
						ca[key + "_" + str(valnew)] = npg.aggregate(zero_strt_sort_noholes_lbls, (self.ca[key] == val).astype('int'), func="sum", fill_value=0)
				elif func == "mode":
					def mode(x):  # type: ignore
						return scipy.stats.mode(x)[0][0]
					ca[key] = npg.aggregate(zero_strt_sort_noholes_lbls, self.ca[key], func=mode, fill_value=0).astype('str')
				elif func == "mean":
					ca[key] = npg.aggregate(zero_strt_sort_noholes_lbls, self.ca[key], func=func, fill_value=0)
				elif func == "first":
					ca[key] = npg.aggregate(zero_strt_sort_noholes_lbls, self.ca[key], func=func, fill_value=self.ca[key][0])

		m = np.empty((self.shape[0], n_groups))
		for (_, selection, view) in self.scan(axis=0, layers=[layer]):
			vals_aggr = npg.aggregate(zero_strt_sort_noholes_lbls, view[layer][:, :], func=aggr_by, axis=1, fill_value=0)
			m[selection, :] = vals_aggr

		if out_file is not None:
			loompy.create(out_file, m, self.ra, ca)

		return m

	def export(self, out_file: str, layer: str = None, format: str = "tab") -> None:
		"""
		Export the specified layer and row/col attributes as tab-delimited file.

		Args:
			out_file:	Path to the output file
			layer:	Name of the layer to export, or None to export the main matrix
			format: Desired file format (only 'tab' is supported)
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


def create_append(filename: str, layers: Union[np.ndarray, Dict[str, np.ndarray], loompy.LayerManager], row_attrs: Dict[str, np.ndarray], col_attrs: Dict[str, np.ndarray], *, file_attrs: Dict[str, str] = None, fill_values: Dict[str, np.ndarray] = None) -> None:
	"""
	**DEPRECATED** - Use `new` instead; see https://github.com/linnarsson-lab/loompy/issues/42
	"""
	deprecated("'create_append' is deprecated. See https://github.com/linnarsson-lab/loompy/issues/42")
	if os.path.exists(filename):
		with connect(filename) as ds:
			ds.add_columns(layers, col_attrs, fill_values=fill_values)
	else:
		create(filename, layers, row_attrs, col_attrs, file_attrs=file_attrs)


def new(filename: str, *, file_attrs: Optional[Dict[str, str]] = None) -> LoomConnection:
	"""
	Create an empty Loom file, and return it as a context manager.
	"""
	if filename.startswith("~/"):
		filename = os.path.expanduser(filename)
	if file_attrs is None:
		file_attrs = {}

	# Create the file (empty).
	# Yes, this might cause an exception, which we prefer to send to the caller
	f = h5py.File(name=filename, mode='w')
	f.create_group('/attrs')  # v3.0.0
	f.create_group('/layers')
	f.create_group('/row_attrs')
	f.create_group('/col_attrs')
	f.create_group('/row_graphs')
	f.create_group('/col_graphs')
	f.flush()
	f.close()

	ds = connect(filename, validate=False)
	for vals in file_attrs:
		if file_attrs[vals] is None:
			ds.attrs[vals] = "None"
		else:
			ds.attrs[vals] = file_attrs[vals]
	# store creation date
	ds.attrs['CreationDate'] = timestamp()
	ds.attrs["LOOM_SPEC_VERSION"] = loompy.loom_spec_version
	return ds


def create(filename: str, layers: Union[np.ndarray, Dict[str, np.ndarray], loompy.LayerManager], row_attrs: Union[loompy.AttributeManager, Dict[str, np.ndarray]], col_attrs: Union[loompy.AttributeManager, Dict[str, np.ndarray]], *, file_attrs: Dict[str, str] = None) -> None:
	"""
	Create a new Loom file from the given data.

	Args:
		filename (str):         The filename (typically using a ``.loom`` file extension)
		layers:					One of the following:

								* Two-dimensional (N-by-M) numpy ndarray of float values
								* Sparse matrix (e.g. :class:`scipy.sparse.csr_matrix`)
								* Dictionary of named layers, each an N-by-M ndarray or sparse matrix
								* A :class:`.LayerManager`, with each layer an N-by-M ndarray
		row_attrs (dict):       Row attributes, where keys are attribute names and values
								are numpy arrays (float or string) of length N
		col_attrs (dict):       Column attributes, where keys are attribute names and
								values are numpy arrays (float or string) of length M
		file_attrs (dict):      Global attributes, where keys are attribute names and
								values are strings
	Returns:
		Nothing

	Remarks:
		If the file exists, it will be overwritten.
	"""

	if isinstance(row_attrs, loompy.AttributeManager):
		row_attrs = {k: v[:] for k, v in row_attrs.items()}
	if isinstance(col_attrs, loompy.AttributeManager):
		col_attrs = {k: v[:] for k, v in col_attrs.items()}

	if isinstance(layers, np.ndarray) or scipy.sparse.issparse(layers):
		layers = {"": layers}
	elif isinstance(layers, loompy.LayerManager):
		layers = {k: v[:, :] for k, v in layers.items()}
	if "" not in layers:
		raise ValueError("Data for default layer must be provided")

	# Sanity checks
	shape = layers[""].shape  # type: ignore
	if shape[0] == 0 or shape[1] == 0:
		raise ValueError("Main matrix cannot be empty")
	for name, layer in layers.items():
		if layer.shape != shape:  # type: ignore
			raise ValueError(f"Layer '{name}' is not the same shape as the main matrix")
	for name, ra in row_attrs.items():
		if len(ra) != shape[0]:
			raise ValueError(f"Row attribute '{name}' is not the same length ({len(ra)}) as number of rows in main matrix ({shape[0]})")
	for name, ca in col_attrs.items():
		if len(ca) != shape[1]:
			raise ValueError(f"Column attribute '{name}' is not the same length ({len(ca)}) as number of columns in main matrix ({shape[1]})")

	try:
		with new(filename, file_attrs=file_attrs) as ds:
			ds.layer[""] = layers[""]
			for key, vals in layers.items():
				if key == "":
					continue
				ds.layer[key] = vals

			for key, vals in row_attrs.items():
				ds.ra[key] = vals

			for key, vals in col_attrs.items():
				ds.ca[key] = vals

	except ValueError as ve:
		if os.path.exists(filename):
			os.remove(filename)
		raise ve


def create_from_cellranger(indir: str, outdir: str = None, genome: str = None, file_attrs: Dict[str, str] = None, selected_cellids: List[str] = None) -> str:
	"""
	Create a .loom file from 10X Genomics cellranger output

	Args:
		indir (str):	path to the cellranger output folder (the one that contains 'outs')
		outdir (str):	output folder where the new loom file should be saved (default to indir)
		genome (str):	genome build to load (e.g. 'mm10'; if None, determine species from outs folder)
		file_attrs:	dict of global file attributes, or None

	Returns:
		path (str):		Full path to the created loom file.

	Remarks:
		The resulting file will be named ``{sampleID}.loom``, where the sampleID is the one given by cellranger.
	"""
	if outdir is None:
		outdir = indir
	sampleid = os.path.split(os.path.abspath(indir))[-1]
	matrix_folder = os.path.join(indir, 'outs', 'filtered_gene_bc_matrices')
	if os.path.exists(matrix_folder):
		if genome is None:
			genome = [f for f in os.listdir(matrix_folder) if not f.startswith(".")][0]
		matrix_folder = os.path.join(matrix_folder, genome)
		matrix = mmread(os.path.join(matrix_folder, "matrix.mtx")).todense()
		genelines = open(os.path.join(matrix_folder, "genes.tsv"), "r").readlines()
		bclines = open(os.path.join(matrix_folder, "barcodes.tsv"), "r").readlines()
	else:  # cellranger V3 file locations
		if genome is None:
			genome = ""  # Genome is not visible from V3 folder
		matrix_folder = os.path.join(indir, 'outs', 'filtered_feature_bc_matrix')
		matrix = mmread(os.path.join(matrix_folder, "matrix.mtx.gz")).todense()
		genelines = [l.decode() for l in gzip.open(os.path.join(matrix_folder, "features.tsv.gz"), "r").readlines()]
		bclines = [l.decode() for l in gzip.open(os.path.join(matrix_folder, "barcodes.tsv.gz"), "r").readlines()]

	accession = np.array([x.split("\t")[0] for x in genelines]).astype("str")
	gene = np.array([x.split("\t")[1].strip() for x in genelines]).astype("str")
	cellids = np.array([sampleid + ":" + x.strip() for x in bclines]).astype("str")

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

	path = os.path.join(outdir, sampleid + ".loom")
	if file_attrs == None:
		file_attrs = {}
	if not "Genome" in file_attrs:
		if genome is None:
			invocationfile = os.path.join(indir, "_invocation")
			for line in open(invocationfile):
				m = re.search('reference_path.+"([^"]+)"', line)
				if m:
					genome = m.group(1)
					break
		if genome:
			file_attrs["Genome"] = genome
	cmdlinefile = os.path.join(indir, "_cmdline")
	if not "Cmdline" in file_attrs and os.path.exists(cmdlinefile):
		cmdline = open(cmdlinefile).read().strip()
		file_attrs["Cmdline"] = cmdline
	versionsfile = os.path.join(indir, "_versions")
	if not "Versions" in file_attrs and os.path.exists(versionsfile):
		versions = open(versionsfile).read().strip().replace("\n", " ")
		file_attrs["Versions"] = versions
	create(path, matrix, row_attrs, col_attrs, file_attrs=file_attrs)
	return path


def create_from_tsv(out_file: str, tsv_file: str, row_metadata_loomfile: str = None, row_metadata_attr: str = "Accession", delim: str = "\t", \
                    dtype: str = "float32", sample_id: str = "", file_attrs: Dict[str, str] = None, \
                    col_metadata_tsv: str = None, metadata_delim: str = '\t') -> None:
	"""
	Create a .loom file from .tsv file

	Args:
		out_file:               path to the newly created .loom file (will be overwritten if it exists)
		tsv_file:		input tab separated data matrix file. Header line should contain cell IDs.
		row_metadata_loomfile:  path to loomfile that will supply gene data and their order. Use to make the new loomfile conform with existing loomfile(s).
		row_metadata_attr:      row attribute of row_metadata_loomfile to use to match with gene IDs of tsv_file
                delim:                  delimiter of expression matrix file
		dtype:                  requested type of loomfile data matrix
		sample_id:              string to use as prefix for cell IDs, or nothing if header fields already include sample IDs
		file_attrs:             dict of global loomfile attributes
                col_metadata_tsv:       metadata for cells. Header line shoud be names of attributes. First column should be CellIDs. Order has to match data matrix file.
                metadata_delim:         delimiter of tsv metadata file
	"""
	id2rowidx = None
	row_attrs = {}
	if row_metadata_loomfile:
		with loompy.connect(row_metadata_loomfile, "r") as ds:
			for attr in ("Accession", "Gene", "Chromosome", "Start", "End"):
				if attr in ds.ra:
					row_attrs[attr] = ds.ra[attr][:]
			nrows = ds.shape[0]
		id2rowidx = { n : i for i, n in enumerate(row_attrs[row_metadata_attr]) }
	with (gzip.open(tsv_file, "rt") if tsv_file.endswith(".gz") else open(tsv_file, "r")) as fd:
		headerrow = fd.readline().rstrip().split(delim)
		datarow1 = fd.readline().rstrip().split(delim)
		headerfirstcellid = 1 if len(datarow1)==len(headerrow) else 0
		if not row_metadata_loomfile:
			nrows = 1
			for line in fd:
				nrows += 1
	geneids = []
	with (gzip.open(tsv_file, "rt") if tsv_file.endswith(".gz") else open(tsv_file, "r")) as fd:
		headerrow = fd.readline().rstrip().split(delim)
		headerrow = [re.sub(r'^"(.+)"$', r'\1', f) for f in headerrow]
		cellids = np.array([ sample_id + cellid for cellid in headerrow[headerfirstcellid:] ]).astype('str')
		matrix = np.zeros([nrows, len(cellids)], dtype=dtype)
		nlines = nnomatch = 0
		for inrowidx, line in enumerate(fd):
			nlines += 1
			row = line.rstrip().split(delim)
			geneid = re.sub(r'^"(.+)"$', r'\1', row[0])
			if id2rowidx:
				if not geneid in id2rowidx:
					nnomatch += 1
					continue
				rowidx = id2rowidx[geneid]
			else:
				rowidx = inrowidx
			datarow = [ float(v) for v in row[1:] ]
			matrix[rowidx, :] = np.array(datarow).astype(dtype)
			geneids.append(geneid)
	if len(row_attrs) == 0:
		row_attrs['Gene'] = geneids
	col_attrs = {"CellID": cellids}
	if col_metadata_tsv:
		with (gzip.open(col_metadata_tsv, "rt") if col_metadata_tsv.endswith(".gz") else open(col_metadata_tsv, "r")) as fd:
			cm_attrs = fd.readline().rstrip().split(metadata_delim)
			cm_attrs = [re.sub(r'^"(.+)"$', r'\1', a) for a in cm_attrs]
			for cm_attr in cm_attrs:
				col_attrs[cm_attr] = []
			cmrowidx = 0
			line = fd.readline()
			while line:
				cm_values = line.rstrip().split(metadata_delim)
				cm_values = [re.sub(r'^"(.+)"$', r'\1', v) for v in cm_values]
				cmcellid = cm_values[0]
				if cmcellid != cellids[cmrowidx]:
					raise ValueError("CellID %s at row %s in %s does not match with corresponding header column %s of %s" % \
					                 (cmcellid, cmrowidx, col_metadata_tsv, cellids[cmrowidx], tsv_file))
				cm_idx0 = 0 if len(cm_values) == len(cm_attrs) else 1
				for colidx, cm_attr in enumerate(cm_attrs):
					col_attrs[cm_attr].append(cm_values[cm_idx0 + colidx])
				cmrowidx += 1
				line = fd.readline()
	create(out_file, matrix, row_attrs, col_attrs, file_attrs=file_attrs)
	print("No match for %s / %s genes in input file. Size of output loomfile: %s" % (nnomatch, nlines, matrix.shape) )

def create_from_matrix_market(out_file: str, sample_id: str, layer_paths: Dict[str, str], row_metadata_path: str, column_metadata_path: str, delim: str = "\t", skip_row_headers: bool = False, skip_colums_headers: bool = False, file_attrs: Dict[str, str] = None, matrix_transposed: bool = False) -> None:
	"""
	Create a .loom file from .mtx matrix market format

	Args:
		out_file:                    path to the newly created .loom file (will be overwritten if it exists)
		sample_id:                   string to use as prefix for cell IDs
		layer_paths:                 dict mapping layer names to paths to the corresponding matrix file (usually with .mtx extension)
		row_metadata_path:           path to the row (usually genes) metadata file
		column_metadata_path:        path to the column (usually cells) metadata file
		delim:                       delimiter used for metadata (default: "\t")
		skip_row_headers:            if true, skip first line in rows metadata file
		skip_column_headers:         if true, skip first line in columns metadata file
		file_attrs:                  dict of global loomfile attributes, or None
		matrix_transposed:           if true, the main matrix is transposed
	
	Remarks:
		layer_paths should typically map the empty string to a matrix market file: {"": "path/to/filename.mtx"}.
		To create a multilayer loom file, map multiple named layers {"": "path/to/layer1.mtx", "layer2": "path/to/layer2.mtx"}
		Note: the created file MUST have a main layer named "". If no such layer is given, BUT all given layers are the same
		datatype, then a main layer will be created as the sum of the other layers. For example, {"spliced": "spliced.mtx", "unspliced": "unspliced.mtx"}
		will create three layers, "", "spliced", and "unspliced", where "" is the sum of the other two.
	"""
	layers: Dict[str, Union[np.ndarray, scipy.sparse.coo_matrix]] = {}

	for name, path in layer_paths.items():
		matrix = mmread(path)
		if matrix_transposed:
			matrix = matrix.T
		layers[name] = matrix
	if "" not in layers:
		main_matrix = None
		for name, matrix in layers.items():
			if main_matrix is None:
				main_matrix = matrix.copy()
			else:
				main_matrix = main_matrix + matrix
		layers[""] = main_matrix

	genelines = open(row_metadata_path, "r").readlines()
	bclines = open(column_metadata_path, "r").readlines()

	accession = np.array([x.split("\t")[0] for x in genelines]).astype("str")
	if(len(genelines[0].split("\t")) > 1):
		gene = np.array([x.split("\t")[1].strip() for x in genelines]).astype("str")
		row_attrs = {"Accession": accession, "Gene": gene}
	else:
		row_attrs = {"Accession": accession}

	cellids = np.array([sample_id + ":" + x.strip() for x in bclines]).astype("str")
	col_attrs = {"CellID": cellids}

	create(out_file, layers[""], row_attrs, col_attrs, file_attrs=file_attrs)

	if len(layers) > 1:
		with loompy.connect(out_file) as ds:
			for name, layer in layers.items():
				if name == "":
					continue
				ds[name] = layer


def create_from_star(indir : str, outfile : str, sample_id : str, \
                     cell_filter : str = "star", expected_n_cells : int = 0, min_total_umis : int = 0, ambient_pthreshold : float = 1.0, \
                     sample_metadata_file : str = None, gtf_file : str = None, main_layer : str = "velosum", extra_layers = None, file_attrs = None):
	"""
		Create a .loom file from STARsolo output
		Args:
			  indir (str):	              path to STARsolo output folder (the one that contains 'Solo.out')
			  outfile (str):              path and name of the new loom file
			  sample_id (str):            sample_id (typically 10Xxxx_x)
			  cell_filter (str):          'emptydrops', 'star', 'none', 'combine', 'Gene', 'GeneFull', 'GeneFull_ExonOverIntron', 'GeneFull_Ex50pAS'
			                              emptydrops: with parameters expected_n_cells, min_total_umis, and ambient_pthreshold
			                              star: uses the barcodes from 'Solo.out/Velocyto/filtered' subdir of indir (Same as 'Gene')
			                              combine: combines cells from all extra_layers and fills layer's non-valid cells with zeroes
			                              Gene... : uses the cells from that Solo.out filtered subdir. Will be added to extra_layers.
			                              none: all cells are kept - Lots of cells!
			  sample_metadata_file (str): file of sample metadata or path to sqlite3 database. Used for gene annotation.
			  gtf_file (str):             path to gtf file used by STARsolo. Used for cell annotations.
			  main_layer(str):            The STAR output matrix to use for the main loom layer (e.g. 'Gene', or 'GeneFull').
			                              'velosum' adds up spliced, unspliced, and ambiguous.
			  extra_layers(list):         Names of additional layers to add in data (e.g. 'GeneFull').
			  file_attrs:                 dict of global file attributes, or None
		Returns:
		           nothing

	"""
	gtypes = extra_layers if extra_layers else []
	if cell_filter in ("Gene", "GeneFull", "GeneFull_ExonOverIntron", "GeneFull_Ex50pAS") and cell_filter not in gtypes:
		gtypes.append(cell_filter)
	if main_layer in ("Gene", "GeneFull", "GeneFull_ExonOverIntron", "GeneFull_Ex50pAS") and main_layer not in gtypes:
		gtypes.append(main_layer)
	subdir = "raw" if cell_filter == "none" else "filtered"
	velodir = os.path.join(indir, "Solo.out", "Velocyto", subdir) if not indir.endswith("Solo.out") else os.path.join(indir, "Velocyto", subdir) 
	accessions = [ l.split('\t')[0] for l in open(os.path.join(velodir, "features.tsv")).readlines() ]
	genes = [ l.split('\t')[1] for l in open(os.path.join(velodir, "features.tsv")).readlines() ]
	gtype2bcs = {}
	n_genes = len(accessions)
	layers = {}
	for gtype in gtypes:
		mtxpath = os.path.join(indir, "Solo.out", gtype, subdir, "matrix.mtx")
		bcspath = os.path.join(indir, "Solo.out", gtype, subdir, "barcodes.tsv")
		if os.path.exists(mtxpath) and os.path.exists(bcspath):
			bcs = [ l.rstrip() for l in open(bcspath).readlines() ]
			gtype2bcs[gtype] = bcs
			mtx = np.loadtxt(mtxpath, skiprows=3, delimiter=' ')
			layers[gtype] = sparse.csc_matrix((mtx[:,2], (mtx[:,0]-1, mtx[:,1]-1)), shape = (n_genes, len(bcs))).todense()
	velobcspath = os.path.join(velodir, "barcodes.tsv")
	velobcs = [ l.rstrip() for l in open(velobcspath).readlines() ]
	gtype2bcs['spliced'] = gtype2bcs['unspliced'] = gtype2bcs['ambiguous'] = velobcs
	velomtxshape = (n_genes, len(velobcs))
	common_mtx = os.path.join(velodir, "matrix.mtx")
	if os.path.exists(common_mtx):
		mtx = np.loadtxt(common_mtx, skiprows=3, delimiter=' ')
		layers['spliced'] = allspliced = sparse.csc_matrix((mtx[:,2], (mtx[:,0]-1, mtx[:,1]-1)), shape = velomtxshape).todense()
		layers['unspliced'] = sparse.csc_matrix((mtx[:,3], (mtx[:,0]-1, mtx[:,1]-1)), shape = velomtxshape).todense()
		layers['ambiguous'] = sparse.csc_matrix((mtx[:,4], (mtx[:,0]-1, mtx[:,1]-1)), shape = velomtxshape).todense()
	else: # STAR >= 2.7.9
		for vtype in ('spliced', 'unspliced', 'ambiguous'):
			mtx = np.loadtxt(os.path.join(velodir, vtype + ".mtx"), skiprows=3, delimiter=' ')
			layers[vtype] = sparse.csc_matrix((mtx[:,2], (mtx[:,0]-1, mtx[:,1]-1)), shape = velomtxshape).todense()
	if cell_filter == "emptydrops":
		allspliced = layers['spliced']
		spliced_total_umis = np.array(allspliced.sum(axis=0))[0]
		ambient_umis, ambient_pvalue = call_cells(allspliced.tocsc(), expected_n_cells)
		valid_cell_idxs = (ambient_pvalue <= ambient_pthreshold) | (spliced_total_umis >= min_total_umis)
		valid_pvalue = ambient_pvalue[valid_cell_idxs]
		valid_bcs = gtype2bcs['spliced'][valid_cell_idxs]
	elif cell_filter == "combine":
		all_bcs = set()
		for bcs in gtype2bcs.values():
			all_bcs.update(bcs)
		valid_bcs = all_bcs
	elif cell_filter in gtypes:
		valid_bcs = gtype2bcs[cell_filter]
	elif cell_filter !="none":
		if cell_filter != "star":
			print (f"Unknown cell_filter type: {cell_filter} - using cells from Solo.out/Velocyto/filtered/")
		valid_bcs = velobcs # gtype2bcs["Gene"]

	valid_bcs = list(valid_bcs)
	n_valid_cells = len(valid_bcs)
	new_shape = (n_genes, n_valid_cells)
	cellids = np.array([f"{sample_id}:{v_bc}" for v_bc in valid_bcs])
	ca = { "CellID": cellids }
	for gtype in layers:
		valid_cells_layer = np.full(new_shape, 0, layers[gtype].dtype)
		ca_valid = np.full(n_valid_cells, 0)
		#print (valid_cells_layer.shape, layers[gtype].shape)
		for fromidx, bc in enumerate(gtype2bcs[gtype]):
			try:
				toidx = valid_bcs.index(bc)
				#print (fromidx, toidx, bc, layers[gtype].shape, valid_cells_layer.shape)
				valid_cells_layer[:,toidx] = layers[gtype][:,fromidx].reshape(n_genes,)
				ca_valid[toidx] = 1
			except ValueError:
				pass #print ("Err:", fromidx, bc, len(valid_bcs))
		layers[gtype] = valid_cells_layer
		#print (gtype, valid_cells_layer.shape, valid_cells_layer.sum())
		ca["Valid_" + gtype] = ca_valid
		#print ("ca", gtype, ca_valid.shape)
	if main_layer == "velosum":
		layers[''] = layers['spliced'] + layers['unspliced'] + layers['ambiguous']
	else:
		layers[''] = layers[main_layer]
		del layers[main_layer]
	total_umis = layers[''].sum(axis=0)
	ca["BarcodeTotalUMIs"] = total_umis
	if cell_filter == "emptydrops":
		ca["AmbientPValue"] = valid_pvalue
		ca["AmbientUMIs"] = np.full(n_valid_cells, ambient_umis)
	if sample_metadata_file:
		try:
			sample_metadata = load_sample_metadata(sample_metadata_file, sample_id)
			for key, value in sample_metadata.items():
				ca[key] = np.full(n_valid_cells, value)
		except:
			print("No sample metadata for %s in %s. Skipping annotation." % (sample_id, sample_metadata_file))
	if gtf_file:
		ra = make_row_attrs_from_gene_metadata(gtf_file, accessions)
	else:
		ra = { "Gene": np.array(genes), "Accession": np.array(accessions) }
	create(filename=outfile, layers=layers, row_attrs=ra, col_attrs=ca, file_attrs=file_attrs)

def create_from_kallistobus(out_file: str, in_dir: str, tr2g_file: str, whitelist_file: str, file_attrs: Dict[str, str] = None, layers: Dict[str, str] = None):
	"""
	Create a loom file from a kallisto-bus output folder.

	Args:
		out_file				Full path to the loom file to be created
		in_dir					Full path to the kallisto-bus directory (containing output.bus, matrix.ec and transcripts.txt)
		whitelist_file			Full path to the barcode whitelist file (e.g. 10xv2_whitelist.txt)
		file_attrs				Optional dictionary of global attributes
		layers					Dict of {layer_name: capture_file_path} to define extra layers
	"""

	pass

	# def import(file: str, key: str)
	# def demote()
		
def combine(files: List[str], output_file: str, key: str = None, file_attrs: Dict[str, str] = None, batch_size: int = 1000, convert_attrs: bool = False) -> None:
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
		Note that the loom files must have exactly the same
		number of rows, and must have exactly the same column attributes.
		Named layers not present in the first file are discarded.

		.. warning::
			If you don't give a ``key`` argument, the files will be combined without changing
			the ordering of rows or columns. Row attributes will be taken from the first file.
			Hence, if rows are not in the same order in all files, the result may be meaningless.

		To guard against this issue, you are strongly advised to provide a ``key`` argument,
		which is used to sort files while merging. The ``key`` should be the name of a row
		attribute that contains a unique value for each row. For example, to order rows by
		the attribute ``Accession``:

		.. highlight:: python
		.. code-block:: python

			import loompy
			loompy.combine(files, key="Accession")

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


def combine_faster(files: List[str], output_file: str, file_attrs: Dict[str, str] = None, selections: List[np.ndarray] = None, key: str = None, skip_attrs: List[str] = None) -> None:
	"""
	Combine loom files and save as a new loom file

	Args:
		files (list of str):    the list of input files (full paths)
		output_file (str):      full path of the output loom file
		file_attrs (dict):      file attributes (title, description, url, etc.)
		selections:				list of indicator arrays (one per file; or None to include all cells for the file) determining which columns to include from each file, or None to include all cells from all files
		key:					row attribute to use as key for ordering the rows, or None to skip ordering
		skip_attrs:				list of column attributes that should not be included in the output

	Returns:
		Nothing, but creates a new loom file combining the input files.

		Note that the loom files must have exactly the same
		number of rows, in exactly the same order, and must have exactly the same column attributes.
		Values in layers missing from one or more files will be replaced by zeros

		.. warning::
			The files will be combined without changing
			the ordering of rows or columns. Row attributes will be taken from the first file.
			Hence, if rows are not in the same order in all files, the result may be meaningless.
	
	Remarks:
		This version assumes that the individual files will fit in memory. If you run out of memory, try the standard combine() method.
	"""
	if file_attrs is None:
		file_attrs = {}
	if skip_attrs is None:
		skip_attrs = []

	if len(files) == 0:
		raise ValueError("The input file list was empty")

	if selections is None:
		selections = [None for _ in files]  # None means take all cells from the file

	n_cells = 0
	n_genes = 0
	for f, s in zip(files, selections):
		with loompy.connect(f, "r") as ds:
			if n_genes == 0:
				n_genes = ds.shape[0]
			elif n_genes != ds.shape[0]:
				raise ValueError(f"All files must have exactly the same number of rows, but {f} had {ds.shape[0]} rows while previous files had {n_genes}")
			if s is None:
				n_cells += ds.shape[1]
			else:
				n_cells += s.sum()

	col_attrs: Dict[str, np.ndarray] = {}
	ix = 0
	with loompy.new(output_file, file_attrs=file_attrs) as dsout:
		for f, s in zip(files, selections):
			with loompy.connect(f, "r") as ds:
				if key is not None:
					ordering = np.argsort(ds.ra[key])
				dsout.shape = (ds.shape[0], n_cells)  # Not really necessary to set this for each file, but no harm either; needed in order to make sure the first layer added will be the right shape
				n_selected = s.sum() if s is not None else ds.shape[1]
				j = 0
				batch_size = 500_000_000 // ds.shape[0] // 4
				for (_, _, view) in ds.scan(items=s, axis=1, key=key, what=["layers"], batch_size=batch_size):
					logging.debug(j)
					for layer in ds.layers.keys():
						# Create the layer if it doesn't exist
						if layer not in dsout.layers:
							dsout[layer] = ds[layer].dtype.name
						# Make sure the dtype didn't change between files
						if dsout.layers[layer].dtype != ds[layer].dtype:
							raise ValueError(f"Each layer must be same datatype in all files, but {layer} of type {ds[layer].dtype} in {f} differs from previous files where it was {dsout[layer].dtype}")
						dsout[layer][:, ix + j: ix + j + view.shape[1]] = view[layer][:, :]
					j += view.shape[1]
				for attr, vals in ds.ca.items():
					if attr in skip_attrs:
						continue
					if attr in col_attrs:
						if col_attrs[attr].dtype != vals.dtype:
							raise ValueError(f"Each column attribute must be same datatype in all files, but {attr} is {vals.dtype} in {f} but was {col_attrs[attr].dtype} in previous files")
					else:
						shape = list(vals.shape)
						shape[0] = n_cells
						col_attrs[attr] = np.zeros(shape, dtype=vals.dtype)
					col_attrs[attr][ix: ix + n_selected] = vals[s]
				for attr, vals in ds.ra.items():
					if attr not in dsout.ra:
						if key is None:
							dsout.ra[attr] = vals
						else:
							dsout.ra[attr] = vals[ordering]
			ix = ix + n_selected
		for attr, vals in col_attrs.items():
			dsout.ca[attr] = vals


def connect(filename: str, mode: str = 'r+', *, validate: bool = True, spec_version: str = "3.0.0") -> LoomConnection:
	"""
	Establish a connection to a .loom file.

	Args:
		filename:		Path to the Loom file to open
		mode:			Read/write mode, 'r+' (read/write) or 'r' (read-only), defaults to 'r+'
		validate:		Validate the file structure against the Loom file format specification
		spec_version:	The loom file spec version to validate against (e.g. "2.0.1" or "old")
	Returns:
		A LoomConnection instance.

	Remarks:
		This function should typically be used as a context manager (i.e. inside a ``with``-block):

		.. highlight:: python
		.. code-block:: python

			import loompy
			with loompy.connect("mydata.loom") as ds:
				print(ds.ca.keys())

		This ensures that the file will be closed automatically when the context block ends

		Note: if validation is requested, an exception is raised if validation fails.
	"""
	return LoomConnection(filename, mode, validate=validate)

