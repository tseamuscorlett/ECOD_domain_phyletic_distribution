import unittest

from main import overlapBool
from main import overlapLength
from main import overlapRange
import portion as P


class OverlapLengthTest(unittest.TestCase):
    def test_no_overlap(self) -> None:
        self.assertEquals(overlapLength(P.closed(1, 50), P.closed(51, 300)), 0)

    def test_partial_overlap(self) -> None:
        self.assertEquals(overlapLength(P.closed(1, 10), P.closed(5, 20)), 6)

    def test_partial_overlap_by_one(self) -> None:
        self.assertEqual(overlapLength(P.closed(1, 5), P.closed(5, 20)), 1)

    def test_nested_overlap(self) -> None:
        self.assertEqual(overlapLength(P.closed(1, 20), P.closed(5, 10)), 6)

    def test_complete_overlap(self) -> None:
        self.assertEqual(overlapLength(P.closed(1, 5), P.closed(1, 5)), 5)

    def test_one_overlap(self) -> None:
        self.assertEqual(overlapLength(P.closed(1, 1), P.closed(1, 20)), 1)


class OverlapBoolTest(unittest.TestCase):
    def test_no_overlap(self) -> None:
        self.assertFalse(overlapBool(P.closed(1, 10), P.closed(11, 15)))

    def test_complete_overlap(self) -> None:
        self.assertTrue(overlapBool(P.closed(1, 10), P.closed(1, 10)))

    def test_partial_overlap(self) -> None:
        self.assertTrue(overlapBool(P.closed(1, 10), P.closed(5, 20)))


class OverlapRangeTest(unittest.TestCase):
    def test_no_overlap(self) -> None:
        self.assertTrue(overlapRange(P.closed(1, 50), P.closed(51, 300)).empty)

    def test_partial_overlap(self) -> None:
        test_range = overlapRange(P.closed(1, 10), P.closed(5, 20))
        self.assertEqual(test_range.lower, 5)
        self.assertEqual(test_range.upper, 10)

    def test_partial_overlap_by_one(self) -> None:
        test_range = overlapRange(P.closed(1, 5), P.closed(5, 20))
        self.assertEqual(test_range.lower, 5)
        self.assertEqual(test_range.upper, 5)

    def test_nested_overlap(self) -> None:
        test_range = overlapRange(P.closed(1, 20), P.closed(5, 10))
        self.assertEqual(test_range.lower, 5)
        self.assertEqual(test_range.upper, 10)

    def test_complete_overlap(self) -> None:
        test_range = overlapRange(P.closed(1, 5), P.closed(1, 5))
        self.assertEqual(test_range.lower, 1)
        self.assertEqual(test_range.upper, 5)

    def test_one_overlap(self) -> None:
        test_range = overlapRange(P.closed(1, 1), P.closed(1, 20))
        self.assertEqual(test_range.lower, 1)
        self.assertEqual(test_range.upper, 1)


if __name__ == '__main__':
    file_path = 'data/GB_GCA_000008085.1_protein.txt'
    import pytest

    pytest.main(['test.py'])
    unittest.main()
