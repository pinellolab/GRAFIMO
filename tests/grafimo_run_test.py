from grafimo.constructVG import build_vg, indexVG
from grafimo.motif_ops import build_motif_MEME
from grafimo.score_sequences import compute_results
from grafimo.extract_regions import get_seqs

import multiprocessing as mp
import numpy as np
import pandas as pd
import pytest
import os


def test_vg_construct():

    test_ref = "test_data/input/test.fa"
    test_vcf = "test_data/input/test.vcf.gz"

    expected_vg = "test_data/expected_results/expected.vg"
    test_vg = "test.vg"

    # create vg and check results
    done = build_vg(test_vg, test_ref, test_vcf, "x", 1)
    print(os.stat(test_vg).st_size)
    print(os.stat(expected_vg).st_size)

    assert ((done == 0) and (os.stat(test_vg).st_size == os.stat(expected_vg).st_size))

# end of test_vg_construction()


def test_vg_index():

    test_vg = "test.vg"
    test_xg = "test.xg"
    test_gbwt = "test.gbwt"
    test_vcf = "test_data/input/test.vcf.gz"

    expected_xg = "test_data/expected_results/expected.xg"
    expected_gbwt = "test_data/expected_results/expected.gbwt"

    # index vg and check results
    done = indexVG(test_vg, test_vcf, 1, True)

    assert ((done == 0) and 
            (os.stat(test_xg).st_size == os.stat(expected_xg).st_size) and
            (os.stat(test_gbwt).st_size == os.stat(expected_gbwt).st_size))

# end of test_vg_index()


def test_sequence_extraction():

    vg = "test.xg"
    region = "x:0-20"
    width = 19

    seqs_extracted = "seqs_extracted.tsv"
    expected_seqs = "test_data/expected_results/expected_seqs.tsv"

    # extract sequences and check results
    query = "vg find -x {0} -E -p {1} -K {2} > {3}".format(vg, region, width, seqs_extracted)
    get_seqs(query)

    result = pd.read_csv(seqs_extracted, sep='\t', header=None).sort_values([1,2,3])
    result.index = range(len(result))  # with sort_values the indexes are not ordered
    expected = pd.read_csv(expected_seqs, sep='\t', header=None).sort_values([1,2,3])
    expected.index = range(len(expected))

    assert (result.equals(expected))

# end of test_sequence_extraction()


def test_motif_processing():

    er_motif_file = "test_data/expected_results/motif_processing_test.txt"
    infile_meme = "test_data/input/MA0139.1.meme"  # CTCF motif in MEME format

    # read the expected processed motif
    er_motif = np.loadtxt(er_motif_file).astype(int)

    # process the motif in MEME format
    proc_motif_meme = build_motif_MEME(infile_meme, "UNIF", 0.1, False,
                                       mp.cpu_count(), False)[0].getMotif_scoreMatrix()

    # check correctness
    assert (proc_motif_meme == er_motif).all()

# end of test_motif_processing_()


def test_scoring():

    infile_meme = "test_data/input/MA0139.1.meme"  # CTCF motif in MEME format
    motif = build_motif_MEME(infile_meme, "UNIF", 0.1, False, mp.cpu_count(), False)[0]

    input_seqs = "test_data/input/"

    # read the expected results
    expected_results = pd.read_csv("test_data/expected_results/scoring_results.tsv", sep="\t", index_col=0)

    # test scoring procedure
    results = compute_results(motif, input_seqs, None, True)

    results.to_csv("scoring_test.tsv", sep="\t")
    results = pd.read_csv("scoring_test.tsv", sep='\t', index_col=0)

    # check correctness
    assert (results.equals(expected_results))

# end of test_scoring()

