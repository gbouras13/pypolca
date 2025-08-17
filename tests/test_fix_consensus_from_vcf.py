"""
Unit tests for dnaapler's fix_consensus_from_vcf.py file
"""

import pypolca.utils.fix_consensus_from_vcf


def test_edit_distance():
    """
    The edit_distance function returns a tuple of (subs, indels).
    """
    ed = pypolca.utils.fix_consensus_from_vcf.edit_distance
    assert ed("", "") == (0, 0)
    assert ed("", "ACGT") == (0, 4)
    assert ed("ACGT", "") == (0, 4)
    assert ed("ACAC", "GTGT") == (4, 0)
    assert ed("AGCTACG", "AGCTACG") == (0, 0)
    assert ed("AGCTACG", "AGCCACG") == (1, 0)
    assert ed("AGCTACG", "AGCTAACG") == (0, 1)
    assert ed("AGCTACG", "AGCTAAACG") == (0, 2)
    assert ed("AGCTACG", "AGCTAAAACG") == (0, 3)
    assert ed("AGCTACG", "ATCTAAAACG") == (1, 3)
    assert ed("GTTTTTTTTTTTAGTGC", "GTTTTTTTTTTAGTGC") == (0, 1)


def test_run_length_encoding():
    rle = pypolca.utils.fix_consensus_from_vcf.run_length_encoding
    assert rle("") == []
    assert rle("A") == [("A", 1)]
    assert rle("AA") == [("A", 2)]
    assert rle("AAA") == [("A", 3)]
    assert rle("AAC") == [("A", 2), ("C", 1)]
    assert rle("AACCC") == [("A", 2), ("C", 3)]
    assert rle("AACCCGG") == [("A", 2), ("C", 3), ("G", 2)]
    assert rle("AACCCGGGT") == [("A", 2), ("C", 3), ("G", 3), ("T", 1)]


def test_is_homopolymer_change_1():
    """
    True homopolymer changes
    """
    ihc = pypolca.utils.fix_consensus_from_vcf.is_homopolymer_change
    assert ihc("CGGGGGGGACT", "CGGGGGGGGACT", 4)
    assert ihc("GTTTTTTTTTTTAGTGC", "GTTTTTTTTTTAGTGC", 4)
    assert ihc("AGGGGGGGGGTGCTGT", "AGGGGGGGGTGCTGT", 4)
    assert ihc("AAAAAA", "AAAA", 4)
    assert ihc("GAAAAAAAAAA", "GAAAAAAAAA", 4)


def test_is_homopolymer_change_2():
    """
    Testing the length threhold
    """
    ihc = pypolca.utils.fix_consensus_from_vcf.is_homopolymer_change
    assert ihc("CGGGGA", "CGGGGGA", 3)
    assert ihc("CGGGGA", "CGGGGGA", 4)
    assert not ihc("CGGGGA", "CGGGGGA", 5)
    assert not ihc("CGGGGA", "CGGGGGA", 6)
    assert ihc("AGGGGGTG", "AGGGGTG", 3)
    assert ihc("AGGGGGTG", "AGGGGTG", 4)
    assert not ihc("AGGGGGTG", "AGGGGTG", 5)
    assert not ihc("AGGGGGTG", "AGGGGTG", 6)
    assert ihc("ACGT", "ACCGT", 1)
    assert ihc("ACGGT", "ACGT", 1)
    assert not ihc("ACGT", "ACCGT", 2)
    assert not ihc("ACGGT", "ACGT", 2)


def test_is_homopolymer_change_3():
    """
    Non-homopolymer changes
    """
    ihc = pypolca.utils.fix_consensus_from_vcf.is_homopolymer_change
    assert not ihc("A", "C", 4)
    assert not ihc("ACGT", "ACCT", 4)
    assert not ihc("AGGATCAGCGAGT", "ACGATCAGCGACT", 4)
    assert not ihc("TACGA", "TAGA", 4)
    assert not ihc("TGA", "TGCA", 4)
    assert not ihc("ACAGCACGACACG", "ACGCACACTACG", 4)


def test_is_homopolymer_change_4():
    """
    Too short to be a homopolymer change
    """
    ihc = pypolca.utils.fix_consensus_from_vcf.is_homopolymer_change
    assert not ihc("CAAAAAAAACCG", "CAAAAAAAAACCG", 9)
    assert not ihc("TCCCCCCCCCTCGCA", "TCCCCCCCCTCGCA", 9)
    assert not ihc("TGGGGGGGGAGGTTT", "TGGGGGGGGGAGGTTT", 9)
    assert not ihc("AAAAAAAA", "AAAAAAAAA", 9)
    assert not ihc("ATTTTTTTTTAAAATTT", "ATTTTTTTTTTAAAATTT", 10)
    assert not ihc("TCCCCCCC", "TCCCCCCCC", 8)
    assert not ihc("ATGACTACGTTTTTTTTTGGA", "ATGACTACGTTTTTTTTTTGGA", 10)


def test_is_homopolymer_change_5():
    """
    A homopolymer change plus another change doesn't count
    """
    ihc = pypolca.utils.fix_consensus_from_vcf.is_homopolymer_change
    assert not ihc("GCCCCCCCTGACTAT", "GCCCCCCCCTGACGAT", 4)
    assert not ihc("CTACGCAAAAAAAACCG", "CTAAGCAAAAAAAAACCG", 4)
    assert not ihc("CCCCCCCTGACT", "CCCCCCCCTGACG", 4)
    assert not ihc("CGCAAAAAAAA", "AGCAAAAAAAAA", 4)


def test_is_homopolymer_change_6():
    """
    Two homopolymer changes together count, as long as both meet the length threshold
    """
    ihc = pypolca.utils.fix_consensus_from_vcf.is_homopolymer_change
    assert ihc("GAAATATTTTGA", "GAAAATATTTTTGA", 3)
    assert not ihc("GAAATATTTTGA", "GAAAATATTTTTGA", 4)
    assert ihc("CGGGGGACTGTTTTTTTTAGTGC", "CGGGGGGACTGTTTTTTTAGTGC", 4)
    assert ihc("CGGGGGACTGTTTTTTTTAGTGC", "CGGGGGGACTGTTTTTTTAGTGC", 5)
    assert not ihc("CGGGGGACTGTTTTTTTTAGTGC", "CGGGGGGACTGTTTTTTTAGTGC", 6)


def test_is_homopolymer_change_7():
    """
    No change doesn't count (shouldn't occur)
    """
    ihc = pypolca.utils.fix_consensus_from_vcf.is_homopolymer_change
    assert not ihc("", "", 4)
    assert not ihc("ACGCG", "ACGCG", 4)
    assert not ihc("TCGAAAAAAAAAAATG", "TCGAAAAAAAAAAATG", 4)
