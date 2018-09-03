import h5py
from typing import *
import logging
import numpy as np
import loompy


class LoomValidator:
	def __init__(self, version: str = "2.0.1") -> None:
		"""
		Args:
			version: 		The Loom file format version to validate against
		"""
		self.version = version
		if version != "2.0.1":
			raise ValueError("This validator can only validate against Loom spec 2.0.1")
		self.errors: List[str] = []
		self.summary: List[str] = []

	def _check(self, condition: bool, message: str) -> bool:
		if not condition:
			self.errors.append(message)
		return condition
	
	def validate(self, path: str, strictness: str = "speconly") -> bool:
		valid1 = True
		with h5py.File(path) as f:
			valid1 = self.validate_spec(f)
			if not valid1:
				self.errors.append("For help, see http://linnarssonlab.org/loompy/format/")

		valid2 = True
		if strictness == "conventions":
			with loompy.connect(path) as ds:
				valid2 = self.validate_conventions(ds)
				if not valid2:
					self.errors.append("For help, see http://linnarssonlab.org/loompy/conventions/")

		return valid1 and valid2

	def validate_conventions(self, ds: loompy.LoomConnection) -> bool:
		"""
		Validate the LoomConnection object against the attribute name/dtype conventions.
		"""
		(n_genes, n_cells) = ds.shape
		if self._check("ClusterID" in ds.ca, "Column attribute 'ClusterID' is missing"):
			self._check(np.issubdtype(ds.ca.ClusterID.dtype, np.int_), "Column attribute 'ClusterID' must be integer dtype")
			self._check(len(np.unique(ds.ca.ClusterID)) == np.max(ds.ca.ClusterID) and np.min(ds.ca.ClusterID) == 0, "Column attribute 'ClusterID' must be integers 0, 1, 2, ... with no missing values")
			self._check(ds.ca.ClusterID.shape == (n_cells,), f"Column attribute 'ClusterID' must be 1-dimensional array of {n_cells} elements")

		if "ClusterName" in ds.ca:
			self._check(ds.ca.ClusterName.dtype == object and np.issubdtype(ds.ca.ClusterName[0].dtype, np.str_), "Column attribute 'ClusterName' must be an array of strings")
			self._check(ds.ca.ClusterName.shape == (n_cells,), f"Column attribute 'ClusterName' must be 1-dimensional array of {n_cells} elements")
			one_to_one = True
			for cid in np.unique(ds.ca.ClusterID):
				if len(np.unique(ds.ca.ClusterName[ds.ca.ClusterID == cid])) != 1:
					one_to_one = False
					break
			for cn in np.unique(ds.ca.ClusterName):
				if len(np.unique(ds.ca.ClusterID[ds.ca.ClusterName == cn])) != 1:
					one_to_one = False
					break
			if not one_to_one:
				self._check(False, "ClusterName must correspond 1:1 with ClusterID")
				
		if self._check("CellID" in ds.ca, "Column attribute 'CellID' is missing"):
			self._check(ds.ca.CellID.dtype == object and np.issubdtype(ds.ca.CellID[0].dtype, np.str_), f"Column attribute 'CellID' must be an array of strings, not '{ds.ca.CellID[0].dtype}'")
			self._check(ds.ca.CellID.shape == (n_cells,), f"Column attribute 'CellID' must be 1-dimensional array of {n_cells} elements")
			self._check(len(np.unique(ds.ca.CellID)) == n_cells, "Column attribute 'CellID' cannot contain duplicate values")

		if "Valid" in ds.ca:
			self._check(np.issubdtype(ds.ca.Valid.dtype, np.int_), f"Column attribute 'Valid' must be integer dtype, not '{ds.ca.Valid.dtype}'")
			valids = np.unique(ds.ca.Valid)
			self._check(np.all(np.isin(ds.ca.Valid, [0, 1])), "Column attribute 'Valid' must be integers 0 or 1 only")
			self._check(ds.ca.Valid.shape == (n_cells,), f"Column attribute 'Valid' must be 1-dimensional array of {n_cells} elements")

		if "Outliers" in ds.ca:
			self._check(np.issubdtype(ds.ca.Outliers.dtype, np.int_), f"Column attribute 'Outliers' must be integer dtype, not '{ds.ca.Outliers.dtype}'")
			self._check(np.all(np.isin(ds.ca.Outliers, [0, 1])), "Column attribute 'Outliers' must be integers 0 or 1 only")
			self._check(ds.ca.Outliers.shape == (n_cells,), f"Column attribute 'Outliers' must be 1-dimensional array of {n_cells} elements")
			
		if self._check("Accession" in ds.ra, "Row attribute 'Accession' is missing"):
			self._check(ds.ra.Accession.dtype == object and np.issubdtype(ds.ra.Accession[0].dtype, np.str_), f"Row attribute 'Accession' must be an array of strings, not '{ds.ra.Accession[0].dtype}'")
			self._check(ds.ra.Accession.shape == (n_genes,), f"Row attribute 'Accession' must be 1-dimensional array of {n_genes} elements")
			self._check(len(np.unique(ds.ra.Accession)) == n_genes, "Row attribute 'Accession' cannot contain duplicate values")

		if self._check("Gene" in ds.ra, "Row attribute 'Gene' is missing"):
			self._check(ds.ra.Gene.dtype == object and np.issubdtype(ds.ra.Gene[0].dtype, np.str_), f"Row attribute 'Gene' must be an array of strings, not '{ds.ra.Gene[0].dtype}'")
			self._check(ds.ra.Gene.shape == (n_genes,), f"Row attribute 'Gene' must be 1-dimensional array of {n_genes} elements")

		if "Valid" in ds.ra:
			self._check(np.issubdtype(ds.ra.Valid.dtype, np.int_), f"Row attribute 'Valid' must be integer dtype, not '{ds.ra.Valid.dtype}'")
			valids = np.unique(ds.ra.Valid)
			self._check(np.all(np.isin(ds.ra.Valid, [0, 1])), "Row attribute 'Valid' must be integers 0 or 1 only")
			self._check(ds.ra.Valid.shape == (n_cells,), f"Row attribute 'Valid' must be 1-dimensional array of {n_cells} elements")

		if "Selected" in ds.ra:
			self._check(np.issubdtype(ds.ra.Selected.dtype, np.int_), f"Row attribute 'Selected' must be integer dtype, not '{ds.ra.Selected.dtype}'")
			valids = np.unique(ds.ra.Selected)
			self._check(np.all(np.isin(ds.ra.Selected, [0, 1])), "Row attribute 'Selected' must be integers 0 or 1 only")
			self._check(ds.ra.Selected.shape == (n_cells,), f"Row attribute 'Selected' must be 1-dimensional array of {n_cells} elements")

		return len(self.errors) == 0
		
	def validate_spec(self, file: h5py.File) -> bool:
		matrix_types = ["float16", "float32", "float64", "int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64"]
		vertex_types = ["int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64"]
		weight_types = ["float16", "float32", "float64"]

		def delay_print(text: str) -> None:
			self.summary.append(text)

		def dt(t: str) -> str:
			if str(t).startswith("|S"):
				return f"string"
			return str(t)

		width_ra = max([len(x) for x in (file["row_attrs"].keys())])
		width_ca = max([len(x) for x in (file["col_attrs"].keys())])
		width_globals = max([len(x) for x in file.attrs.keys()])
		width_layers = 0
		if "layers" in file and len(file["layers"]) > 0:
			width_layers = max([len(x) for x in file["layers"].keys()])
		width_layers = max(width_layers, len("Main matrix"))
		width = max(width_ca, width_ra, width_globals)

		delay_print("Global attributes:")
		for key, value in file.attrs.items():
			delay_print(f"{key: >{width}}: {value} {str(type(value))}")

		if self._check("matrix" in file, "Main matrix missing"):
			self._check(file["matrix"].dtype in matrix_types, f"Main matrix dtype={file['matrix'].dtype} is not allowed")
			shape = file["matrix"].shape
			delay_print(f"Layers shape={shape}:")
			delay_print(f"{'Main matrix': >{width}} {file['matrix'].dtype}")

		if "layers" in file:
			for layer in file["layers"]:
				self._check(file["layers"][layer].shape == shape, f"Layer '{layer}' shape {file['layers'][layer].shape} does not match main matrix shape {shape}")
				self._check(file["layers"][layer].dtype in matrix_types, f"Layer '{layer}' dtype={file['layers'][layer].dtype} is not allowed")
				delay_print(f"{layer: >{width}} {file['layers'][layer].dtype}")

		delay_print("Row attributes:")
		if self._check("row_attrs" in file, "'row_attrs' group is missing"):
			for ra in file["row_attrs"]:
				self._check(file["row_attrs"][ra].shape[0] == shape[0], f"Row attribute '{ra}' shape {file['row_attrs'][ra].shape[0]} first dimension does not match row dimension {shape}")
				self._check(file["row_attrs"][ra].dtype in matrix_types or np.issubdtype(file['row_attrs'][ra].dtype, np.string_), f"Row attribute '{ra}' dtype {file['row_attrs'][ra].dtype} is not allowed")
				ra_shape = file['row_attrs'][ra].shape
				delay_print(f"{ra: >{width}} {dt(file['row_attrs'][ra].dtype)} {ra_shape if len(ra_shape) > 1 else ''}")
			if len(file["row_attrs"]) == 0:
				delay_print("    (none)")

		delay_print("Column attributes:")
		if self._check("col_attrs" in file, "'col_attrs' group is missing"):
			for ca in file["col_attrs"]:
				self._check(file["col_attrs"][ca].shape[0] == shape[1], f"Column attribute '{ca}' shape {file['col_attrs'][ca].shape[0]} first dimension does not match column dimension {shape}")
				self._check(file["col_attrs"][ca].dtype in matrix_types or np.issubdtype(file["col_attrs"][ca].dtype, np.string_), f"Column attribute '{ca}' dtype {file['col_attrs'][ca].dtype} is not allowed")
				ca_shape = file['col_attrs'][ca].shape
				delay_print(f"{ca: >{width}} {dt(file['col_attrs'][ca].dtype)} {ca_shape if len(ca_shape) > 1 else ''}")
			if len(file["col_attrs"]) == 0:
				delay_print("    (none)")

		delay_print("Row graphs:")
		if self._check("row_graphs" in file, "'row_graphs' group is missing"):
			for g in file["row_graphs"]:
				self._check("a" in file["row_graphs"][g], f"Row graph '{g}' is missing vector 'a', denoting start vertices")
				self._check(file["row_graphs"][g]['a'].dtype in vertex_types, f"/row_graphs/{g}/a.dtype {file['row_graphs'][g]['a'].dtype} must be integer")
				self._check("b" in file["row_graphs"][g], f"Row graph '{g}' is missing vector 'b', denoting end vertices")
				self._check(file["row_graphs"][g]['b'].dtype in vertex_types, f"/row_graphs/{g}/b.dtype {file['row_graphs'][g]['b'].dtype} must be integer")
				self._check("w" in file["row_graphs"][g], f"Row graph '{g}' is missing vector 'w', denoting vertex weights")
				self._check(file["row_graphs"][g]['w'].dtype in weight_types, f"/row_graphs/{g}/w.dtype {file['row_graphs'][g]['w'].dtype} must be float")
				self._check(file['row_graphs'][g]['a'].shape[0] == file['row_graphs'][g]['b'].shape[0] and file['row_graphs'][g]['a'].shape[0] == file['row_graphs'][g]['w'].shape[0], f"Row graph '{g}' sparse vectors a, b and w must have equal length")
				delay_print(f"    '{g}' with {file['row_graphs'][g]['a'].shape[0]} edges")
			if len(file["row_graphs"]) == 0:
				delay_print("    (none)")

		delay_print("Column graphs:")
		if self._check("col_graphs" in file, "'col_graphs' group is missing"):
			for g in file["col_graphs"]:
				self._check("a" in file["col_graphs"][g], f"Column graph '{g}' is missing vector 'a', denoting start vertices")
				self._check(file["col_graphs"][g]['a'].dtype in vertex_types, f"/col_graphs/{g}/a.dtype {file['col_graphs'][g]['a'].dtype} must be integer")
				self._check("b" in file["col_graphs"][g], f"Column graph '{g}' is missing vector 'b', denoting end vertices")
				self._check(file["col_graphs"][g]['b'].dtype in vertex_types, f"/col_graphs/{g}/b.dtype {file['col_graphs'][g]['b'].dtype} must be integer")
				self._check("w" in file["col_graphs"][g], f"Column graph '{g}' is missing vector 'w', denoting vertex weights")
				self._check(file["col_graphs"][g]['w'].dtype in weight_types, f"/col_graphs/{g}/w.dtype {file['col_graphs'][g]['w'].dtype} must be float")
				self._check(file['col_graphs'][g]['a'].shape[0] == file['col_graphs'][g]['b'].shape[0] and file['col_graphs'][g]['a'].shape[0] == file['col_graphs'][g]['w'].shape[0], f"Column graph '{g}' sparse vectors a, b and w must have equal length")
				delay_print(f"    '{g}' with {file['col_graphs'][g]['a'].shape[0]} edges")
			if len(file["col_graphs"]) == 0:
				delay_print("    (none)")

		return len(self.errors) == 0
