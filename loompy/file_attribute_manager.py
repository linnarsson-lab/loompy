from typing import *
import h5py
import scipy.sparse as sparse
import numpy as np
import loompy


class FileAttributeManager(object):
	def __init__(self, f: h5py.File) -> None:
		setattr(self, "!f", f)
		storage: Dict[str, str] = {}
		setattr(self, "!storage", storage)
		for key, val in f.attrs.items():
			materialized = loompy.materialize_attr_values(val)
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
				if name in self.f.attrs:
					val = self.f.attrs[name]
					materialized = loompy.materialize_attr_values(val)
					self.__dict__["storage"][name] = materialized
					return materialized
			raise AttributeError(f"'{type(self)}' object has no attribute '{name}'")

	def __setitem__(self, name: str, val: Any) -> None:
		return self.__setattr__(name, val)

	def __setattr__(self, name: str, val: Any) -> None:
		if name.startswith("!"):
			super(FileAttributeManager, self).__setattr__(name[1:], val)
		else:
			if self.f is not None:
				normalized = loompy.normalize_attr_values(val)
				self.f.attrs[name] = normalized
				self.f.flush()
				val = self.f.attrs[name]
				# Read it back in to ensure it's synced and normalized
				normalized = loompy.materialize_attr_values(val)
				self.__dict__["storage"][name] = normalized
