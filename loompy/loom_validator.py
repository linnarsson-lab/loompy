import h5py
from typing import *
import logging
import numpy as np


class LoomValidator:
	def __init__(self, strictness: str = "speconly", version: str = "2.0.1") -> None:
		"""
		Args:
			strictness:		"speconly" (basic conformance to the Loom file format specification)
							"conventions" (basic conformance, plus adherence to attribute conventions)
			version: 		The Loom file format version to validate against
		"""
		self.strictness = strictness
		self.version = version
		if version != "2.0.1":
			raise ValueError("This validator can only validate against Loom spec 2.0.1")

	def validate_spec(self, file: h5py.File, verbose: bool = True) -> None:
		matrix_types = ["float16", "float32", "float64", "int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64"]
		vertex_types = ["int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64"]
		weight_types = ["float16", "float32", "float64"]
		has_errors = False
		to_print = []

		def check(condition: bool, message: str) -> None:
			if not condition:
				logging.error(message)
				has_errors = True

		def delay_print(text: str) -> None:
			to_print.append(text)

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

		check("matrix" in file, "Main matrix missing")
		check(file["matrix"].dtype in matrix_types, f"Main matrix dtype={file['matrix'].dtype} is not allowed")
		shape = file["matrix"].shape
		delay_print(f"Layers shape={shape}:")
		delay_print(f"{'Main matrix': >{width}} {file['matrix'].dtype}")

		if "layers" in file:
			for layer in file["layers"]:
				check(file["layers"][layer].shape == shape, f"Layer '{layer}' shape {file['layers'][layer].shape} does not match main matrix shape {shape}")
				check(file["layers"][layer].dtype in matrix_types, f"Layer '{layer}' dtype={file['layers'][layer].dtype} is not allowed")
				if verbose:
					delay_print(f"{layer: >{width}} {file['layers'][layer].dtype}")

		delay_print("Row attributes:")
		check("row_attrs" in file, "'row_attrs' group is missing")
		for ra in file["row_attrs"]:
			check(file["row_attrs"][ra].shape[0] == shape[0], f"Row attribute '{ra}' shape {file['row_attrs'][ra].shape[0]} first dimension does not match row dimension {shape}")
			check(file["row_attrs"][ra].dtype in matrix_types or np.issubdtype(file['row_attrs'][ra].dtype, np.string_), f"Row attribute '{ra}' dtype {file['row_attrs'][ra].dtype} is not allowed")
			if verbose:
				ra_shape = file['row_attrs'][ra].shape
				delay_print(f"{ra: >{width}} {dt(file['row_attrs'][ra].dtype)} {ra_shape if len(ra_shape) > 1 else ''}")
		if len(file["row_attrs"]) == 0:
			delay_print("    (none)")

		delay_print("Column attributes:")
		check("col_attrs" in file, "'col_attrs' group is missing")
		for ca in file["col_attrs"]:
			check(file["col_attrs"][ca].shape[0] == shape[1], f"Column attribute '{ca}' shape {file['col_attrs'][ca].shape[0]} first dimension does not match column dimension {shape}")
			check(file["col_attrs"][ca].dtype in matrix_types or np.issubdtype(file["col_attrs"][ca].dtype, np.string_), f"Column attribute '{ca}' dtype {file['col_attrs'][ca].dtype} is not allowed")
			ca_shape = file['col_attrs'][ca].shape
			delay_print(f"{ca: >{width}} {dt(file['col_attrs'][ca].dtype)} {ca_shape if len(ca_shape) > 1 else ''}")
		if len(file["col_attrs"]) == 0:
			delay_print("    (none)")

		delay_print("Row graphs:")
		check("row_graphs" in file, "'row_graphs' group is missing")
		for g in file["row_graphs"]:
			check("a" in file["row_graphs"][g], f"Row graph '{g}' is missing vector 'a', denoting start vertices")
			check(file["row_graphs"][g]['a'].dtype in vertex_types, f"/row_graphs/{g}/a.dtype {file['row_graphs'][g]['a'].dtype} must be integer")
			check("b" in file["row_graphs"][g], f"Row graph '{g}' is missing vector 'b', denoting end vertices")
			check(file["row_graphs"][g]['b'].dtype in vertex_types, f"/row_graphs/{g}/b.dtype {file['row_graphs'][g]['b'].dtype} must be integer")
			check("w" in file["row_graphs"][g], f"Row graph '{g}' is missing vector 'w', denoting vertex weights")
			check(file["row_graphs"][g]['w'].dtype in weight_types, f"/row_graphs/{g}/w.dtype {file['row_graphs'][g]['w'].dtype} must be float")
			check(file['row_graphs'][g]['a'].shape[0] == file['row_graphs'][g]['b'].shape[0] and file['row_graphs'][g]['a'].shape[0] == file['row_graphs'][g]['w'].shape[0], f"Row graph '{g}' sparse vectors a, b and w must have equal length")
			delay_print(f"    '{g}' with {file['row_graphs'][g]['a'].shape[0]} edges")
		if len(file["row_graphs"]) == 0:
			delay_print("    (none)")

		delay_print("Column graphs:")
		check("col_graphs" in file, "'col_graphs' group is missing")
		for g in file["col_graphs"]:
			check("a" in file["col_graphs"][g], f"Column graph '{g}' is missing vector 'a', denoting start vertices")
			check(file["col_graphs"][g]['a'].dtype in vertex_types, f"/col_graphs/{g}/a.dtype {file['col_graphs'][g]['a'].dtype} must be integer")
			check("b" in file["col_graphs"][g], f"Column graph '{g}' is missing vector 'b', denoting end vertices")
			check(file["col_graphs"][g]['b'].dtype in vertex_types, f"/col_graphs/{g}/b.dtype {file['col_graphs'][g]['b'].dtype} must be integer")
			check("w" in file["col_graphs"][g], f"Column graph '{g}' is missing vector 'w', denoting vertex weights")
			check(file["col_graphs"][g]['w'].dtype in weight_types, f"/col_graphs/{g}/w.dtype {file['col_graphs'][g]['w'].dtype} must be float")
			check(file['col_graphs'][g]['a'].shape[0] == file['col_graphs'][g]['b'].shape[0] and file['col_graphs'][g]['a'].shape[0] == file['col_graphs'][g]['w'].shape[0], f"Column graph '{g}' sparse vectors a, b and w must have equal length")
			delay_print(f"    '{g}' with {file['col_graphs'][g]['a'].shape[0]} edges")
		if len(file["col_graphs"]) == 0:
			delay_print("    (none)")

		if verbose:
			for line in to_print:
				print(line)
		assert not has_errors

