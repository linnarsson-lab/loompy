import os
from tempfile import NamedTemporaryFile
from unittest import TestCase

import h5py
import numpy as np

from loompy import GlobalAttributeManager


class GlobalAttributeManagerTests(TestCase):
	VALUE_IN_FILE = np.arange(3)

	def setUp(self):
		f = NamedTemporaryFile(suffix="loom")
		f.close()
		self.filename = f.name
		self.file = h5py.File(f.name)
		self.file.attrs["arr"] = self.VALUE_IN_FILE

	def tearDown(self):
		self.file.close()
		os.remove(self.filename)

	def test_get(self):
		m = GlobalAttributeManager(self.file)
		default = np.arange(2, 5)

		val = m.get("arr")
		np.testing.assert_array_equal(val, self.VALUE_IN_FILE)

		val = m.get("arr", default)
		np.testing.assert_array_equal(val, self.VALUE_IN_FILE)

		val = m.get("missing")
		self.assertIsNone(val)

		val = m.get("missing", default)
		np.testing.assert_array_equal(val, default)
