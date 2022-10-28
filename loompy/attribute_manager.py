from multiprocessing.sharedctypes import Value
from typing import *
import numpy as np
import h5py
import loompy
from loompy import timestamp
from .utils import compare_loom_spec_version


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
			for key in self.ds._file[a].keys():
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
		Note: if the attributes do not contain a timestamp, and the mode is 'r+', a new timestamp is created and returned.
		Otherwise, the current time in UTC will be returned.
		"""
		a = ["/row_attrs/", "/col_attrs/"][self.axis]

		if self.ds is not None:
			if name is None:
				if "last_modified" in self.ds._file[a].attrs:
					return self.ds._file[a].attrs["last_modified"]
				elif self.ds._file.mode == 'r+':
					self.ds._file[a].attrs["last_modified"] = timestamp()
					self.ds._file.flush()
					return self.ds._file[a].attrs["last_modified"]
			if name is not None:
				if "last_modified" in self.ds._file[a + name].attrs:
					return self.ds._file[a + name].attrs["last_modified"]
				elif self.ds._file.mode == 'r+':
					self.ds._file[a + name].attrs["last_modified"] = timestamp()
					self.ds._file.flush()
					return self.ds._file[a + name].attrs["last_modified"]
		return timestamp()

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
				if val is None:
					raise ValueError(key)
				am[key] = val[thing]
			return am
		elif type(thing) is tuple:
			# A tuple of strings giving alternative names for attributes
			result: np.ndarray = None
			for t in thing:
				if t in self.__dict__["storage"]:
					if result is None:
						result = self.__getattr__(t)
					else:
						vals = self.__getattr__(t)
						if vals.dtype != result.dtype:
							raise AttributeError(f"Cannot stack attributes of different types ({vals.dtype} and {result.dtype})")
						result = np.vstack((result, vals)).transpose()
			if result is None:
				raise AttributeError(f"'{type(self)}' object has no attribute {thing}")
			else:
				return result
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
		elif "/" in name:
			raise KeyError("Attribute name cannot contain slash (/)")
		else:
			if self.ds is not None:
				values = loompy.normalize_attr_values(val, compare_loom_spec_version(self.ds._file, "3.0.0") >= 0)
				a = ["/row_attrs/", "/col_attrs/"][self.axis]
				if self.ds.shape[self.axis] != 0 and values.shape[0] != self.ds.shape[self.axis]:
					raise ValueError(f"Attribute '{name}' must have exactly {self.ds.shape[self.axis]} values but {len(values)} were given")
				if self.ds._file[a].__contains__(name):
					del self.ds._file[a + name]

				self.ds._file.create_dataset(
					a + name,
					data=values,
					dtype=h5py.special_dtype(vlen=str) if values.dtype == np.object_ else values.dtype,
					maxshape=(values.shape[0], ) if len(values.shape) == 1 else (values.shape[0], None),
					fletcher32=False,
					compression="gzip",
					shuffle=False,
					compression_opts=2
				)
				self.ds._file[a + name].attrs["last_modified"] = timestamp()
				self.ds._file[a].attrs["last_modified"] = timestamp()
				self.ds._file.attrs["last_modified"] = timestamp()
				self.ds._file.flush()
				self.__dict__["storage"][name] = loompy.materialize_attr_values(self.ds._file[a][name][:])
			else:
				self.__dict__["storage"][name] = val

	def __delitem__(self, name: str) -> None:
		"""
		Remove a named attribute
		"""
		return self.__delattr__(name)

	def __delattr__(self, name: str) -> None:
		"""
		Remove a named attribute
		"""
		if self.ds is not None:
			a = ["/row_attrs/", "/col_attrs/"][self.axis]
			if self.ds._file[a].__contains__(name):
				del self.ds._file[a + name]
				self.ds._file.flush()
		if name in self.__dict__["storage"]:
			del self.__dict__["storage"][name]

	def _permute(self, ordering: np.ndarray) -> None:
		"""
		Permute all the attributes in the collection

		Remarks:
			This permutes the order of the values for each attribute in the file
		"""
		for key in self.keys():
			self[key] = self[key][ordering]

	def get(self, name: str, default: np.ndarray) -> np.ndarray:
		"""
		Return the value for a named attribute if it exists, else default.
		Default has to be a numpy array of correct size.
		"""

		if name in self:
			return self[name]
		else:
			if not isinstance(default, np.ndarray):
				raise ValueError(f"Default must be an np.ndarray with exactly {self.ds.shape[self.axis]} values")

			if default.shape[0] != self.ds.shape[self.axis]:
				raise ValueError(f"Default must be an np.ndarray with exactly {self.ds.shape[self.axis]} values but {len(default)} were given")

			return default
