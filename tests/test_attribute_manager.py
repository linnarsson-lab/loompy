from unittest import TestCase
import numpy as np
from loompy import AttributeManager


class AttributeManagerTests(TestCase):
	def setUp(self) -> None:
		self.ds = self.DummyLoomConnection()

	def test_get(self) -> None:
		m = AttributeManager(self.ds, axis=0)
		default = np.array([2., 3, 4, 5, 6, 7, 8])

		# should return value for attribute in file
		val = m.get("arr", default)
		np.testing.assert_array_equal(
			val, np.array([1., 2, 3, 4, 5, 6, 7])
		)

		# should return default for missing attribute
		val = m.get("missing", default)
		np.testing.assert_array_equal(val, default)

	def test_get_raises_on_invalid_default_value(self) -> None:
		m = AttributeManager(self.ds, axis=0)

		# wrong type
		with self.assertRaises(ValueError):
			m.get("missing", None)

		# wrong size
		with self.assertRaises(ValueError):
			m.get("missing", np.array([1., 2, 3, 4, 5]))

		m = AttributeManager(self.ds, axis=1)
		with self.assertRaises(ValueError):
			m.get("missing", np.array([1., 2, 3, 4, 5, 6, 7]))

	class DummyLoomConnection:
		shape = (7, 5)
		_file = {
			"/row_attrs/": {
				"arr": np.array([1., 2, 3, 4, 5, 6, 7])
			},
			"/col_attrs/": {
				"arr": np.array([1., 2, 3, 4, 5])
			}
		}
