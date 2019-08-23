import loompy
import numpy as np
import scipy.sparse as sparse

print(loompy.call_cells(sparse.csc_matrix(np.random.randint(0, 10, size=(500, 100000), dtype=np.uint16)), expected_n_cells=5000))
