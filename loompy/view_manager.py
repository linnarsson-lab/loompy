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
		if type(slice_) is not tuple or len(slice_) is not 2:
			raise ValueError("Views require slices along two dimensions")

		rows = slice_[0]
		cols = slice_[1]

		ra = self.ds.ra[rows]
		row_graphs = self.ds.row_graphs[rows]
		ca = self.ds.ca[cols]
		col_graphs = self.ds.col_graphs[cols]
		layers = self.ds.layer[rows, cols]

		return loompy.LoomView(layers, ra, ca, row_graphs, col_graphs, filename=self.ds.filename, file_attrs=self.ds.attrs)
