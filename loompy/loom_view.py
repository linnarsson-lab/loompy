from typing import *
import loompy
import numpy as np


class LoomView:
	"""
	An in-memory loom dataset
	"""
	def __init__(self, layers: loompy.LayerManager, row_attrs: loompy.AttributeManager, col_attrs: loompy.AttributeManager, row_graphs: loompy.GraphManager, col_graphs: loompy.GraphManager, *, filename: str, file_attrs: loompy.FileAttributeManager) -> None:
		self.filename = filename
		self.view = loompy.ViewManager(self)
		self.layers = layers
		self.shape = [layer.shape for (name, layer) in layers.items()][0]
		self.ra = row_attrs
		self.ca = col_attrs
		self.row_graphs = row_graphs
		self.col_graphs = col_graphs
		self.attrs = file_attrs

		# Compatibility with loompy v1.x
		self.layer = layers
		self.row_attrs = row_attrs
		self.col_attrs = col_attrs

	def __getitem__(self, slice: Tuple[Union[int, np.ndarray, slice], Union[int, np.ndarray, slice]]) -> np.ndarray:
		"""
		Get a slice of the main matrix.
		Args:
			slice:		A 2D slice object (see http://docs.h5py.org/en/latest/high/dataset.html) or np.ndarrays or ints
		Returns:
			A numpy matrix
		"""
		return self.layers[""][slice]

	def _repr_html_(self) -> str:
		"""
		Return an HTML representation of the loom view, showing the upper-left 10x10 corner.
		"""
		return loompy.to_html(self)

	def permute(self, ordering: np.ndarray, *, axis: int) -> None:
		"""
		Permute the view, by permuting its layers, attributes and graphs

		Args:
			ordering (np.ndarray):	The desired ordering along the axis
			axis (int):				0, permute rows; 1, permute columns
		"""
		for layer in self.layers:
			layer.permute(ordering, axis=axis)
		for g in (self.row_graphs, self.col_graphs)[axis]:
			g.permute(ordering)
		for a in (self.row_attrs, self.col_attrs)[axis]:
			a.permute(ordering)
