import numpy as np
import loompy
from typing import *


class ViewManager:
	"""
	Create views by slicing an underlying LoomConnection or LoomView
	"""
	def __init__(self, ds: Any) -> None:
		self.ds = ds

	def __getitem__(self, slice_: Tuple[Union[slice, np.ndarray, int], Union[slice, np.ndarray, int]]) -> loompy.LoomView:
		"""
		Create a new view by slicing through the loom file or view

		Args:
			slice_ (2-tuple of slice, int or np.ndarray): 	How to slice the file or view

		Returns:
			A LoomView object, an in-memory representation of the sliced file
		"""
		if type(slice_) is not tuple or len(slice_) != 2:
			raise ValueError("Views require slices along two dimensions")

		rows = slice_[0]
		cols = slice_[1]

		ra = self.ds.ra[rows]
		row_graphs = self.ds.row_graphs[rows]
		ca = self.ds.ca[cols]
		col_graphs = self.ds.col_graphs[cols]
		layers = self.ds.layer[rows, cols]

		return loompy.LoomView(layers, ra, ca, row_graphs, col_graphs, filename=self.ds.filename, file_attrs=self.ds.attrs)
