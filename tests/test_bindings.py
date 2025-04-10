import strobealign


def test_logger():
    logger = strobealign.Logger.get()
    # Enable debug logging
    logger.set_level(strobealign.LOG_LEVELS.LOG_DEBUG)
    # Only log on error
    logger.set_level(strobealign.LOG_LEVELS.LOG_ERROR)


def test_references():
    refs = strobealign.References.from_fasta("tests/phix.fasta")
    assert len(refs) == 1
    assert refs[0].name == "NC_001422.1"
    assert refs[0].sequence.startswith("GAGTTTTATC")


def test_index_parameters():
    params = strobealign.IndexParameters.from_read_length(100)
    assert isinstance(params.syncmer.k, int)
    assert isinstance(params.syncmer.s, int)
    assert isinstance(params.syncmer.t, int)
    assert isinstance(params.randstrobe.w_min, int)
    assert isinstance(params.randstrobe.w_max, int)


def test_reverse_complement():
    assert strobealign.reverse_complement("") == ""
    assert strobealign.reverse_complement("A") == "T"
    assert strobealign.reverse_complement("AAAACCCGGT") == "ACCGGGTTTT"


def test_indexing_and_match_finding():
    refs = strobealign.References.from_fasta("tests/phix.fasta")
    index_parameters = strobealign.IndexParameters.from_read_length(100)
    index = strobealign.StrobemerIndex(refs, index_parameters)
    index.populate()

    # Find NAMs for a single query sequence (matches in forward orientation on
    # the reference)
    query = "TGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCG"
    randstrobes = strobealign.randstrobes_query(query, index_parameters)
    # For this test, we ignore the randstrobes for the reverse-complemented query
    randstrobes = randstrobes[0]
    hits = strobealign.find_hits(randstrobes, index, use_mcs=False)
    for hit in hits:
        reference_index = index.reference_index(hit.position)
        ref = refs[reference_index].sequence
        reference_start = index.get_strobe1_position(hit.position)
        assert not hit.is_partial
        reference_end = reference_start + index.strobe2_offset(hit.position) + index.k
        ref_aligned = ref[reference_start:reference_end]
        query_aligned = query[hit.query_start:hit.query_end]
        assert ref_aligned == query_aligned

    matches = strobealign.hits_to_matches(hits, index)
    assert matches
    nams = strobealign.merge_matches_into_nams(matches, index.k, sort=False, is_revcomp=False)
    assert nams
    print(nams)


def test_index_find():
    refs = strobealign.References.from_fasta("tests/phix.fasta")
    index_parameters = strobealign.IndexParameters.from_read_length(100)
    index = strobealign.StrobemerIndex(refs, index_parameters)
    index.populate()

    query = "TGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCG"
    query_randstrobes = strobealign.randstrobes_query(query, index_parameters)[0]
    assert query_randstrobes
    # First randstrobe must be found
    assert index.find(query_randstrobes[0].hash)

    n = 0
    for qr in query_randstrobes:
        for rs in index.find(qr.hash):
            n += 1
            assert rs.hash == qr.hash
    # Ensure the for loop did test something
    assert n > 1
