import unittest

from main import overlapBool
from main import overlapLength
from main import overlapRange


class OverlapLengthTest(unittest.TestCase):
    def test_no_overlap(self) -> None:
        self.assertEquals(overlapLength([1, 50], [51, 300]), 0)

    def test_partial_overlap(self) -> None:
        self.assertEquals(overlapLength([1, 10], [5, 20]), 6)

    def test_partial_overlap_by_one(self) -> None:
        self.assertEqual(overlapLength([1, 5], [5, 20]), 1)

    def test_nested_overlap(self) -> None:
        self.assertEqual(overlapLength([1, 20], [5, 10]), 6)

    def test_complete_overlap(self) -> None:
        self.assertEqual(overlapLength([1, 5], [1, 5]), 5)

    def test_one_overlap(self) -> None:
        self.assertEqual(overlapLength([1, 1], [1, 20]), 1)


class OverlapBoolTest(unittest.TestCase):
    def test_no_overlap(self) -> None:
        self.assertFalse(overlapBool([1, 10], [11, 15]))

    def test_complete_overlap(self) -> None:
        self.assertTrue(overlapBool([1, 10], [1, 10]))

    def test_partial_overlap(self) -> None:
        self.assertTrue(overlapBool([1, 10], [5, 20]))


class OverlapRangeTest(unittest.TestCase):
    def test_no_overlap(self) -> None:
        self.assertEqual(overlapRange([1, 50], [51, 300]), [])

    def test_partial_overlap(self) -> None:
        self.assertEqual(overlapRange([1, 10], [5, 20]), [5, 10])

    def test_partial_overlap_by_one(self) -> None:
        self.assertEqual(overlapRange([1, 5], [5, 20]), [5, 5])

    def test_nested_overlap(self) -> None:
        self.assertEqual(overlapRange([1, 20], [5, 10]), [5, 10])

    def test_complete_overlap(self) -> None:
        self.assertEqual(overlapRange([1, 5], [1, 5]), [1, 5])

    def test_one_overlap(self) -> None:
        self.assertEqual(overlapRange([1, 1], [1, 20]), [1, 1])


if __name__ == '__main__':
    import pytest

    pytest.main(['test.py'])
    unittest.main()
