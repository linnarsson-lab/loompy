from typing import *
import numpy as np
import loompy
from loompy import timestamp


class AttributeManager:
	"""
	Manage a set of attributes (either row or column) with a backing HDF5 file store
	"""
	def __init__(self, ds: Any, *, axis: int) -> None:  # Note: can't give type for ds because it will be circular and mypy doesn't support it
		setattr(self, "!axis", axis)
		setattr(self, "!ds", ds)
		storage: Dict[str, np.ndarray] = {}
		setattr(self, "!storage", storage)

		if ds is not None:
			a = ["/row_attrs/", "/col_attrs/"][self.axis]
			for key in ds._file[a].keys():
				self.__dict__["storage"][key] = None

	def keys(self) -> List[str]:
		"Return the attribute names"
		return list(self.__dict__["storage"].keys())

	def items(self) -> Iterable[Tuple[str, np.ndarray]]:
		"Return an iterator over attribute (name, value) tuples"
		for key in self.keys():
			yield (key, self[key])

	def __len__(self) -> int:
		"Return the number of attributes"
		return len(self.keys())

	def __contains__(self, name: str) -> bool:
		"Return True if attribute exists"
		return name in self.keys()

	def __iter__(self) -> Iterator[str]:
		for key in self.keys():
			yield key

	def last_modified(self, name: str = None) -> str:
		"""
		Return a compact ISO8601 timestamp (UTC timezone) indicating when an attribute was last modified

		Note: if no attribute name is given (the default), the modification time of the most recently modified attribute will be returned

		If an attribute has no timestamp and the loom file is connected in read-write mode, it will be initialised to current time.

		If an attribute has no timestamp and the loom file is connected in read-only mode, "19700101T000000Z" (start of Unix Time) is returned.
		"""
		a = ["/row_attrs/", "/col_attrs/"][self.axis]
		ds = self.ds
		if ds is not None:
			if name is None:
				if "last_modified" in ds._file[a].attrs:
					return ds._file[a].attrs["last_modified"]
				elif ds.mode == "r+":
					ds._file[a].attrs["last_modified"] = timestamp()
					ds._file.flush()
					return ds._file[a].attrs["last_modified"]
			else:
				if "last_modified" in ds._file[a + name].attrs:
					return ds._file[a + name].attrs["last_modified"]
				elif ds.mode == "r+":
					ds._file[a + name].attrs["last_modified"] = timestamp()
					ds._file.flush()
					return ds._file[a + name].attrs["last_modified"]
		# Fallback: start of Unix Time
		return "19700101T000000Z"

	def __getitem__(self, thing: Any) -> np.ndarray:
		"""
		Return a named attribute, or a slice through all the attributes

		Args:
			thing:		if string, return the named attribute
						if slice, np.ndarray or int, return a slice through all the attributes
		"""
		if type(thing) is slice or type(thing) is np.ndarray or type(thing) is int:
			am = AttributeManager(None, axis=self.axis)
			for key, val in self.items():
				am[key] = val[thing]
			return am
		else:
			return self.__getattr__(thing)

	def __getattr__(self, name: str) -> np.ndarray:
		"""
		Return the named attribute

		Args:
			name (str) 	Name of the attribute

		Remarks:
			The values will be loaded from file, and properly HTML unescaped
		"""
		try:
			vals = self.__dict__["storage"][name]
			if vals is None:
				# Read values from the HDF5 file
				a = ["/row_attrs/", "/col_attrs/"][self.axis]
				vals = loompy.materialize_attr_values(self.ds._file[a][name][:])
				self.__dict__["storage"][name] = vals
			return vals
		except KeyError:
			raise AttributeError(f"'{type(self)}' object has no attribute '{name}'")

	def __setitem__(self, name: str, val: np.ndarray) -> None:
		"""
		Set the value of a named attribute
		"""
		return self.__setattr__(name, val)

	def __setattr__(self, name: str, val: np.ndarray) -> None:
		"""
		Set the value of a named attribute

		Args:
			name (str) 			Name of the attribute
			val (np.ndarray)	Value of the attribute

		Remarks:
			Length must match the corresponding matrix dimension
			The values are automatically HMTL escaped and converted to ASCII for storage
		"""
		if name.startswith("!"):
			super(AttributeManager, self).__setattr__(name[1:], val)
		elif self.ds is None:
			raise IOError("Attribute not set: LoomConnection is None")
		elif self.ds.closed:
			raise IOError("Attribute not set: cannot modify closed LoomConnection")
		elif self.ds.mode != 'r+':
			raise IOError("Attribute not set: cannot modify loom file when connected in read-only mode")
		else:
			ds = self.ds
			values = loompy.normalize_attr_values(val)
			a = ["/row_attrs/", "/col_attrs/"][self.axis]
			if ds.shape[self.axis] != 0 and values.shape[0] != ds.shape[self.axis]:
				raise ValueError(f"Attribute must have exactly {ds.shape[self.axis]} values but {len(values)} were given")
			if ds._file[a].__contains__(name):
				del ds._file[a + name]
			ds._file[a + name] = values  # TODO: for 2D arrays, use block compression along columns/rows
			ds._file[a + name].attrs["last_modified"] = timestamp()
			ds._file[a].attrs["last_modified"] = timestamp()
			ds._file.attrs["last_modified"] = timestamp()
			ds._file.flush()
			self.__dict__["storage"][name] = loompy.materialize_attr_values(ds._file[a][name][:])

	def __delitem__(self, name: str) -> None:
		"""
		Remove a named attribute
		"""
		return self.__delattr__(name)

	def __delattr__(self, name: str) -> None:
		"""
		Remove a named attribute
		"""
		if self.ds is None:
			raise IOError("Attribute not deleted: LoomConnection is None")
		elif self.ds.closed:
			raise IOError("Attribute not deleted: cannot modify closed LoomConnection")
		elif self.ds.mode != 'r+':
			raise IOError("Attribute not deleted: cannot modify loom file when connected in read-only mode")
		else:
			a = ["/row_attrs/", "/col_attrs/"][self.axis]
			if self.ds._file[a].__contains__(name):
				del self.ds._file[a + name]
				self.ds._file.flush()
			if name in self.__dict__["storage"]:
				del self.__dict__["storage"][name]

	def permute(self, ordering: np.ndarray) -> None:
		"""
		Permute all the attributes in the collection

		Remarks:
			This permutes the order of the values for each attribute in the file
		"""
		if self.ds is None:
			raise IOError("Attribute order not permuted: LoomConnection is None")
		elif self.ds.closed:
			raise IOError("Attribute order not permuted: cannot modify closed LoomConnection")
		elif self.ds.mode != 'r+':
			raise IOError("Attribute order not permuted: cannot modify loom file when connected in read-only mode")
		else:
			for key in self.keys():
				self[key] = self[key][ordering]
