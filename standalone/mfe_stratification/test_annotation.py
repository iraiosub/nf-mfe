"""Boundary and classification tests; synthetic coordinates are BED-style."""
import unittest
import pandas as pd
from stratify_mfe import IntervalIndex, annotate


class AnnotationTests(unittest.TestCase):
    def setUp(self):
        self.index = {('exon', '1', '+'): IntervalIndex([(10, 20), (18, 30), (50, 60)]),
                      ('gene', '1', '+'): IntervalIndex([(0, 100)])}

    def pair(self, ll, lr, rl, rr, **kwargs):
        row = dict(lchr='chr1', ll=ll, lr=lr, lstrand='+', rchr='1', rl=rl, rr=rr, rstrand='+')
        row.update(kwargs)
        return annotate(pd.DataFrame([row]), self.index).iloc[0]

    def test_full_exon_boundaries(self):
        row = self.pair(10, 20, 50, 60)
        self.assertEqual(row.annotation_class, 'exonic')
        self.assertEqual(row.loop_size_nt, 30)

    def test_union_of_exons_is_not_single_exon(self):
        self.assertEqual(self.pair(10, 30, 50, 60).annotation_class, 'mixed_or_partial')

    def test_intronic_and_partial(self):
        self.assertEqual(self.pair(30, 40, 60, 70).annotation_class, 'intronic')
        self.assertEqual(self.pair(19, 35, 60, 70).annotation_class, 'mixed_or_partial')
        self.assertEqual(self.pair(110, 120, 130, 140).annotation_class, 'unannotated')

    def test_strand_and_special_loops(self):
        self.assertEqual(self.pair(10, 20, 50, 60, lstrand='-', rstrand='-').annotation_class, 'unannotated')
        self.assertEqual(self.pair(50, 60, 10, 20).loop_size_nt, 30)
        self.assertEqual(self.pair(50, 60, 10, 20, lstrand='-', rstrand='-').loop_size_nt, 30)
        self.assertEqual(self.pair(10, 20, 50, 60, lstrand='-', rstrand='-').loop_size_nt, 30)
        self.assertEqual(self.pair(10, 20, 15, 25).loop_size_bin, 'overlapping')
        self.assertEqual(self.pair(10, 20, 20, 30).loop_size_nt, 0)
        self.assertEqual(self.pair(10, 20, 50, 60, rstrand='-').loop_size_bin, 'opposite_strand')
        self.assertEqual(self.pair(10, 20, 50, 60, rchr='2').loop_size_bin, 'trans_chromosomal')


if __name__ == '__main__':
    unittest.main()
