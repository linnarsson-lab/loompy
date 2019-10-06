from typing import Dict, List, Iterable, Tuple, Any, Iterator
import scipy.sparse as sparse
import numpy as np
import loompy
import warnings
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	import h5py
from .utils import compare_loom_spec_version


class GlobalAttributeManager(object):
	def __init__(self, f: h5py.File) -> None:
		setattr(self, "!f", f)
		storage: Dict[str, str] = {}
		setattr(self, "!storage", storage)
		if "attrs" not in self.f:
			for key, val in f.attrs.items():
				materialized = loompy.materialize_attr_values(val)
				self.__dict__["storage"][key] = materialized
		else:
			for key, val in f["attrs"].items():
				materialized = loompy.materialize_attr_values(val[()])
				self.__dict__["storage"][key] = materialized

	def keys(self) -> List[str]:
		return list(self.__dict__["storage"].keys())

	def items(self) -> Iterable[Tuple[str, sparse.coo_matrix]]:
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
		return self.__getattr__(thing)

	def __getattr__(self, name: str) -> np.ndarray:
		try:
			return self.__dict__["storage"][name]
		except KeyError:
			if self.f is not None:
				if loompy.compare_loom_spec_version(self.f, "3.0.0") < 0:
					if name in self.f.attrs:
						val = self.f.attrs[name]
					else:
						raise AttributeError(f"File has no global attribute '{name}'")
				else:
					if name in self.f["attrs"]:
						val = self.f["attrs"][name]
					else:
						raise AttributeError(f"File has no global attribute '{name}'")
				materialized = loompy.materialize_attr_values(val)
				self.__dict__["storage"][name] = materialized
				return materialized

	def __setitem__(self, name: str, val: Any) -> None:
		return self.__setattr__(name, val)

	def __setattr__(self, name: str, val: Any) -> None:
		if name.startswith("!"):
			super(GlobalAttributeManager, self).__setattr__(name[1:], val)
		elif "/" in name:
			raise KeyError("Attribute name cannot contain slash (/)")
		else:
			if self.f is not None:
				if loompy.compare_loom_spec_version(self.f, "3.0.0") < 0 and "attrs" not in self.f["/"]:
					normalized = loompy.normalize_attr_values(val, False)
					self.f.attrs[name] = normalized
					self.f.flush()
					val = self.f.attrs[name]
					# Read it back in to ensure it's synced and normalized
					normalized = loompy.materialize_attr_values(val)
					self.__dict__["storage"][name] = normalized
				else:
					normalized = loompy.normalize_attr_values(val, True)
					if name in self.f["attrs"]:
						del self.f["attrs"][name]
					if not np.isscalar(normalized) and normalized.dtype == np.object_:
						self.ds._file.create_dataset("attrs/" + name, data=normalized, dtype=h5py.special_dtype(vlen=str))
					else:
						self.f["attrs"][name] = normalized
					self.f.flush()
					val = self.f["attrs"][name][()]
					# Read it back in to ensure it's synced and normalized
					normalized = loompy.materialize_attr_values(val)
					self.__dict__["storage"][name] = normalized

	def __delattr__(self, name: str) -> None:
		if name.startswith("!"):
			super(GlobalAttributeManager, self).__delattr__(name[1:])
		else:
			if self.f is not None:
				if loompy.compare_loom_spec_version(self.f, "3.0.0") < 0:
					if name in self.f.attrs:
						del self.f.attrs[name]
				else:
					if name in self.f["attrs"]:
						del self.f["attrs"][name]
			del self.__dict__["storage"][name]

	def get(self, name: str, default: Any = None) -> np.ndarray:
		"""
		Return the value for a named attribute if it exists, else default.
		If default is not given, it defaults to None, so that this method never raises a KeyError.
		"""
		if name in self:
			return self[name]
		else:
			return default
