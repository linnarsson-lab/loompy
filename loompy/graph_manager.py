import scipy.sparse as sparse
import numpy as np
from typing import *
from loompy import timestamp


def _renumber(a: np.ndarray, keys: np.ndarray, values: np.ndarray) -> np.ndarray:
	"""
	Renumber 'a' by replacing any occurrence of 'keys' by the corresponding 'values'
	"""
	ordering = np.argsort(keys)
	keys = keys[ordering]
	values = keys[ordering]
	index = np.digitize(a.ravel(), keys, right=True)
	return(values[index].reshape(a.shape))


class GraphManager:
	"""
	Manage a set of graphs (either for rows or columns) with a backing HDF5 file store
	"""
	def __init__(self, ds: Any, *, axis: int) -> None:
		setattr(self, "!axis", axis)
		setattr(self, "!ds", ds)
		storage: Dict[str, np.ndarray] = {}
		setattr(self, "!storage", storage)

		if ds is not None:
			# Patch old files that use the old naming convention
			if ds._file.mode == "r+":
				if "row_graphs" not in ds._file:
					ds._file.create_group('/row_graphs')
				if "col_graphs" not in ds._file:
					ds._file.create_group('/col_graphs')
				if "row_edges" in ds._file:
					for key in ds._file["row_edges"]:
						ds._file["row_graphs"][key] = ds._file["row_edges"][key]
					del ds._file["row_edges"]
				if "col_edges" in ds._file:
					for key in ds._file["col_edges"]:
						ds._file["col_graphs"][key] = ds._file["col_edges"][key]
					del ds._file["col_edges"]

			a = ["row_graphs", "col_graphs"][self.axis]
			if a in ds._file:
				for key in ds._file[a]:
					self.__dict__["storage"][key] = None
			else:
				if ds.mode == "r+":
					ds._file.create_group(a)

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

	def last_modified(self, name: str = None) -> str:
		"""
		Return a compact ISO8601 timestamp (UTC timezone) indicating when a graph was last modified

		Note: if no graph name is given (the default), the modification time of the most recently modified graph will be returned
		Note: if the graphs do not contain a timestamp, and the mode is 'r+', a new timestamp is created and returned.
		Otherwise, the current time in UTC will be returned.
		"""
		a = ["row_graphs", "col_graphs"][self.axis]

		if name is None:
			if "last_modified" in self.ds._file[a].attrs:
				return self.ds._file[a].attrs["last_modified"]
			elif self.ds._file.mode == 'r+':
				self.ds._file[a].attrs["last_modified"] = timestamp()
				self.ds._file.flush()
				return self.ds._file[a].attrs["last_modified"]
		if name is not None:
			if "last_modified" in self.ds._file[a + name].attrs:
				return self.ds._file[a][name].attrs["last_modified"]
			elif self.ds._file.mode == 'r+':
				self.ds._file[a][name].attrs["last_modified"] = timestamp()
				self.ds._file.flush()
				return self.ds._file[a][name].attrs["last_modified"]
		return timestamp()

	def __getitem__(self, thing: Any) -> sparse.coo_matrix:
		if type(thing) is slice or type(thing) is np.ndarray or type(thing) is int:
			gm = GraphManager(None, axis=self.axis)
			for key, g in self.items():
				# Slice the graph matrix properly without making it dense
				(a, b, w) = (g.row, g.col, g.data)
				indices = np.arange(g.shape[0])[thing]
				mask = np.logical_and(np.in1d(a, indices), np.in1d(b, indices))
				a = a[mask]
				b = b[mask]
				w = w[mask]
				d = dict(zip(np.sort(indices), np.arange(indices.shape[0])))
				a = np.array([d[x] for x in a])
				b = np.array([d[x] for x in b])
				gm[key] = sparse.coo_matrix((w, (a, b)), shape=(len(indices), len(indices)))
			return gm
		elif type(thing) is tuple:
			# A tuple of strings giving alternative names for graphs
			for t in thing:
				if t in self.__dict__["storage"]:
					return self.__getattr__(t)
			raise AttributeError(f"'{type(self)}' object has no attribute {thing}")
		else:
			return self.__getattr__(thing)

	def __getattr__(self, name: str) -> sparse.coo_matrix:
		try:
			g = self.__dict__["storage"][name]
			if g is None:
				# Read values from the HDF5 file
				a = ["row_graphs", "col_graphs"][self.axis]
				r = self.ds._file[a][name]["a"]
				c = self.ds._file[a][name]["b"]
				w = self.ds._file[a][name]["w"]
				g = sparse.coo_matrix((w, (r, c)), shape=(self.ds.shape[self.axis], self.ds.shape[self.axis]))
				self.__dict__["storage"][name] = g
			return g
		except KeyError:
			raise AttributeError(f"'{type(self)}' object has no graph '{name}' on axis {self.axis}")

	def __setitem__(self, name: str, g: sparse.coo_matrix) -> None:
		return self.__setattr__(name, g)

	def __setattr__(self, name: str, g: sparse.coo_matrix) -> None:
		if name.startswith("!"):
			super(GraphManager, self).__setattr__(name[1:], g)
		elif "/" in name:
			raise KeyError("Graph name cannot contain slash (/)")
		else:
			g = sparse.coo_matrix(g)
			if self.ds is not None:
				a = ["row_graphs", "col_graphs"][self.axis]
				if g.shape[0] != self.ds.shape[self.axis] or g.shape[1] != self.ds.shape[self.axis]:
					raise ValueError(f"Adjacency matrix shape for axis {self.axis} must be ({self.ds.shape[self.axis]},{self.ds.shape[self.axis]}) but shape was {g.shape}")
				if name in self.ds._file[a]:
					del self.ds._file[a][name]["a"]
					del self.ds._file[a][name]["b"]
					del self.ds._file[a][name]["w"]
					del self.ds._file[a][name]
				self.ds._file[a].create_group(name)
				self.ds._file[a][name]["a"] = g.row
				self.ds._file[a][name]["b"] = g.col
				self.ds._file[a][name]["w"] = g.data
				self.ds._file[a][name].attrs["last_modified"] = timestamp()
				self.ds._file[a].attrs["last_modified"] = timestamp()
				self.ds._file.attrs["last_modified"] = timestamp()
				self.ds._file.flush()
				self.__dict__["storage"][name] = g
			else:
				self.__dict__["storage"][name] = g

	def __delitem__(self, name: str) -> None:
		return self.__delattr__(name)

	def __delattr__(self, name: str) -> None:
		if self.ds is not None:
			a = ["row_graphs", "col_graphs"][self.axis]
			if self.ds._file[a].__contains__(name):
				del self.ds._file[a][name]["a"]
				del self.ds._file[a][name]["b"]
				del self.ds._file[a][name]["w"]
				del self.ds._file[a][name]
				self.ds._file.flush()
		if name in self.__dict__["storage"]:
			del self.__dict__["storage"][name]
	
	def _permute(self, ordering: np.ndarray) -> None:
		for name in self.keys():
			g = self[name]
			(a, b, w) = (g.row, g.col, g.data)
			a = _renumber(a, np.array(ordering), np.arange(g.shape[1]))
			b = _renumber(b, np.array(ordering), np.arange(g.shape[1]))
			g = sparse.coo_matrix((w, (a, b)), g.shape)
			self[name] = g
