import loompy
import numpy as np
from typing import *
import scipy.sparse as sparse


def test_new() -> None:
	with loompy.new("test.loom") as ds:
		m = np.zeros((20, 100))
		ra = {"Gene": [x for x in "ABCDEFGHILJKLMNOPQRS"]}
		ca = {"Cell": np.arange(100)}
		ds.add_columns(m, ca, row_attrs=ra)
		ds.add_columns(m, ca, row_attrs=ra)
	with loompy.connect("test.loom") as ds:
		assert(ds.shape == (20, 200))


def test_sparse() -> None:
	G = 1000
	C = 100
	S = sparse.eye(G, C)

	loompy.create('test.loom', S, {'g_id': np.arange(G)}, {'c_id': np.arange(C)})
	with loompy.connect("test.loom") as ds:
		ds["layer"] = S
		assert(np.all(ds[:, :] == S.toarray()))
		assert(np.all(ds.sparse().data == S.tocoo().data))
		assert(np.all(ds.layers["layer"][:, :] == S.toarray()))
		assert(np.all(ds.layers["layer"].sparse().data == S.tocoo().data))
