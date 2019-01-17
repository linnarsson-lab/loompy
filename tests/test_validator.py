import os
from tempfile import NamedTemporaryFile
from unittest import TestCase

import numpy as np
import loompy
from loompy import LoomValidator


class LoomValidatorTests(TestCase):
    def test_file_with_empty_col_attrs_is_valid(self) -> None:
        f = NamedTemporaryFile(suffix=".loom")
        f.close()
        loompy.create(f.name, np.zeros((5, 5)), {}, {})
        try:
            self.assertTrue(
                LoomValidator().validate(f.name),
                "File with empty col_attrs or row_attrs should be valid"
            )
        finally:
            os.remove(f.name)
