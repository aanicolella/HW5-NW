# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    # create NWA object and run alignment
    nwa = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', -10, -1)
    nwa.align(seq1, seq2)
    
    # create ground truth matrices for nwa._align_matrix, nwa._gapA_matrix, and nwa._gapB_matrix
    gt_align_matrix = np.array([[0., -np.inf, -np.inf, -np.inf],
                                [-np.inf, 5., -11., -13.],
                                [-np.inf, -12., 4., -8.],
                                [-np.inf, -12., -1., 5.],
                                [-np.inf, -14., -6., 4.]])
    
    gt_gapA_matrix = np.array([[-10., -np.inf, -np.inf, -np.inf],
                                [-11., -12., -6., -7.],
                                [-12., -13., -14., -7.],
                                [-13., -14., -15., -12.],
                                [-14., -15., -16., -17.]])
    
    gt_gapB_matrix = np.array([[-10., -11., -12., -13.],
                                [-np.inf, -12., -13., -14.],
                                [-np.inf, -6., -14., -15.],
                                [-np.inf, -7., -7., -16.],
                                [-np.inf, -8., -8., -6.]])
    # assert correct values for all three
    np.testing.assert_array_equal(nwa._align_matrix, gt_align_matrix)
    np.testing.assert_array_equal(nwa._gapA_matrix, gt_gapA_matrix)
    np.testing.assert_array_equal(nwa._gapB_matrix, gt_gapB_matrix)

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    # create NWA object and run alignment
    nwa = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', -10, -1)
    score, seq3_align, seq4_align = nwa.align(seq3, seq4)

    # assert values == ground truth for score and alignment sequences
    assert score == 17
    assert seq3_align == 'MAVHQLIRRP'
    assert seq4_align == 'M---QLIRHP'





