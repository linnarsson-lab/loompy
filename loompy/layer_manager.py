from typing import *
import numpy as np
import loompy
import logging


class LayerManager:
	"""
	Manage a set of layers with a backing HDF5 file store
	"""
	def __init__(self, ds: Any) -> None:  # Note: can't give type for ds because it will be circular and mypy doesn't support it
		setattr(self, "!ds", ds)
		storage: Dict[str, np.ndarray] = {}
		setattr(self, "!storage", storage)

		if ds is not None:
			self.__dict__["storage"][""] = None
			if "layers" in ds._file:
				for key in self.ds._file["layers"].keys():
					self.__dict__["storage"][key] = None
			elif self.ds.mode == 'r+':
				ds._file.create_group('/layers')
			else:
				logging.warn("Missing layers group not added: cannot modify loom file when connected in read-only mode")

	def last_modified(self, name: str = None) -> str:
		"""
		Return a compact ISO8601 timestamp (UTC timezone) indicating when the layer was last modified

		Note: if name is None, the modification time of the most recently modified layer is returned

		If a layer has no timestamps and the loom file is connected in read-write mode, it will be initialised to current time.

		If a layer has no timestamps and the loom file is connected in read-only mode, "19700101T000000Z" (start of Unix Time) is returned.
		"""
		if name is not None:
			return self[name].last_modified()
		ts = None
		for name in self.keys():
			if ts is None:
				ts = self[name].last_modified()
			else:
				_ts = self[name].last_modified()
				if _ts > ts:
					ts = _ts
		return ts

	def keys(self) -> List[str]:
		return list(self.__dict__["storage"].keys())

	def items(self) -> Iterable[Tuple[str, np.ndarray]]:
		for key in self.keys():
			yield (key, self[key])

	def __len__(self) -> int:
		return len(self.keys())

	def __contains__(self, name: str) -> bool:
		return name in self.keys()

	def __iter__(self) -> Iterator[str]:
		for key in self.keys():
			yield key

	def __getitem__(self, thing: Any) -> np.ndarray:
		"""
		Access a layer by name, or slice through all the layers

		Args:
			thing:		if string, return the specified layer ("" is the default layer)
						if slice 2-tuple, return a new LayerManager with all layers sliced
		"""
		if type(thing) is str:
			return self.__getattr__(thing)
		else:
			# Assume some kind of slice
			lm = LayerManager(None)
			for key, layer in self.items():
				lm[key] = loompy.MemoryLoomLayer(key, layer[thing])
			return lm

	def __getattr__(self, name: str) -> np.ndarray:
		try:
			vals = self.__dict__["storage"][name]
			if vals is None:
				# Read values from the HDF5 file
				return loompy.LoomLayer(name, self.ds)
			return vals
		except KeyError:
			raise AttributeError(f"'{type(self)}' object has no attribute '{name}'")

	def __setitem__(self, name: str, val: np.ndarray) -> None:
		return self.__setattr__(name, val)

	def __setattr__(self, name: str, val: np.ndarray) -> None:
		if name.startswith("!"):
			super(LayerManager, self).__setattr__(name[1:], val)
		elif self.ds is None:
			raise IOError("Layer not set: LoomConnection is None")
		elif self.ds.closed:
			raise IOError("Layer not set: cannot modify closed LoomConnection")
		elif self.ds.mode != 'r+':
			raise IOError("Layer not set: cannot modify loom file when connected in read-only mode")
		else:
			if self.ds is not None:
				matrix = val
				if not np.isfinite(matrix).all():
					raise ValueError("INF and NaN not allowed in loom matrix")
				if not (np.issubdtype(matrix.dtype, np.integer) or np.issubdtype(matrix.dtype, np.floating)):
					raise ValueError("Matrix elements must be integer or float")
				if not self.ds._file.__contains__("/layers"):
					self.ds._file.create_group("/layers")

				# make sure chunk size is not bigger than actual matrix size
				chunks = (min(64, matrix.shape[0]), min(64, matrix.shape[1]))
				path = "/layers/" + name
				if name == "":
					path = "/matrix"
				if self.ds._file.__contains__(path):
					del self.ds._file[path]

				# Save the matrix
				self.ds._file.create_dataset(
					path,
					data=matrix,
					maxshape=(matrix.shape[0], None),
					chunks=chunks,
					fletcher32=False,
					compression="gzip",
					shuffle=False,
					compression_opts=2
				)
				if name == "":
					self.ds.shape = matrix.shape
				self.ds._file.flush()
				self.__dict__["storage"][name] = None

	def __delitem__(self, name: str) -> None:
		return self.__delattr__(name)

	def __delattr__(self, name: str) -> None:
		if self.ds is None:
			if name in self.__dict__["storage"]:
				del self.__dict__["storage"][name]
		elif self.ds.closed:
			raise IOError("Layer not deleted: cannot modify closed LoomConnection")
		elif self.ds.mode != 'r+':
			raise IOError("Layer not deleted: cannot modify loom file when connected in read-only mode")
		if name == "":
			raise ValueError("Layer not deleted: cannot delete default layer")
		else:
			path = "/layers/" + name
			if self.ds._file.__contains__(path):
				del self.ds._file[path]
			self.ds._file.flush()

	def permute(self, ordering: np.ndarray, *, axis: int) -> None:
		if self.ds is None:
			raise IOError("Layer not permuted: LoomConnection is None")
		elif self.ds.closed:
			raise IOError("Layer not permuted: cannot modify closed LoomConnection")
		elif self.ds.mode != 'r+':
			raise IOError("Layer not permuted: cannot modify loom file when connected in read-only mode")
		for key in self.keys():
			self[key].permute(ordering, axis=axis)
