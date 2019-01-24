import os
from tempfile import NamedTemporaryFile
from unittest import TestCase

import numpy as np
import loompy


class LoomConnectionTests(TestCase):
    def setUp(self) -> None:
        self.file = NamedTemporaryFile(suffix=".loom")
        self.file.close()
        loompy.create(
            self.file.name,
            np.random.random((5, 5)),
            row_attrs={
                "key": np.fromiter(range(5), dtype=np.int)
            },
            col_attrs={
                "key": np.fromiter(range(5), dtype=np.int)
            })

    def tearDown(self) -> None:
        os.remove(self.file.name)

    def test_scan_with_default_ordering(self) -> None:
        with loompy.connect(self.file.name) as ds:
            for axis in [0, 1]:
                _, _, view = next(iter(ds.scan(axis=axis)))
                no_ordering_data = view[:, :]

                _, _, view = next(iter(ds.scan(axis=axis, key="key")))
                original_ordering_data = view[:, :]

        np.testing.assert_almost_equal(no_ordering_data, original_ordering_data,
                                       err_msg="Default ordering should same as in file")
