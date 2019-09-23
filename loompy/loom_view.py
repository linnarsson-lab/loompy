from typing import *
import loompy
import numpy as np


class LoomView:
	"""
	An in-memory loom dataset
	"""
	def __init__(self, layers: loompy.LayerManager, row_attrs: loompy.AttributeManager, col_attrs: loompy.AttributeManager, row_graphs: loompy.GraphManager, col_graphs: loompy.GraphManager, *, filename: str = "", file_attrs: loompy.GlobalAttributeManager = None) -> None:
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

	def __getitem__(self, slice_: Union[str, Tuple[Union[int, np.ndarray, slice], Union[int, np.ndarray, slice]]]) -> np.ndarray:
		"""
		Get a slice of the main matrix.
		Args:
			slice:		A 2D slice object (see http://docs.h5py.org/en/latest/high/dataset.html) or np.ndarrays or ints
		Returns:
			A numpy matrix
		"""
		if type(slice_) is str:
			return self.layers[slice_]
		else:
			return self.layers[""][slice_]

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
		if axis not in (0, 1):
			raise ValueError("Axis must be 0 (rows) or 1 (columns)")
		for layer in self.layers.values():
			layer._permute(ordering, axis=axis)
		if axis == 0:
			if self.row_graphs is not None:
				for g in self.row_graphs.values():
					g._permute(ordering)
			for a in self.row_attrs.values():
				a._permute(ordering)
		elif axis == 1:
			if self.col_graphs is not None:
				for g in self.col_graphs.values():
					g._permute(ordering)
			for a in self.col_attrs.values():
				a._permute(ordering)
