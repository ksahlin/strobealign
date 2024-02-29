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


def test_indexing_and_nams_finding():
    refs = strobealign.References.from_fasta("tests/phix.fasta")
    index_parameters = strobealign.IndexParameters.from_read_length(100)
    index = strobealign.StrobemerIndex(refs, index_parameters)
    index.populate()

    # Find NAMs for a single query sequence
    query = "TGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCG"
    randstrobes = strobealign.randstrobes_query(query, index_parameters)
    nams = strobealign.find_nams(randstrobes, index)
    assert nams
    for nam in nams:
        ref = refs[nam.reference_index].sequence
        ref_aligned = ref[nam.ref_start:nam.ref_end]
        query_aligned = query[nam.query_start:nam.query_end]
        score = nam.score


def test_index_find():
    refs = strobealign.References.from_fasta("tests/phix.fasta")
    index_parameters = strobealign.IndexParameters.from_read_length(100)
    index = strobealign.StrobemerIndex(refs, index_parameters)
    index.populate()

    query = "TGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCG"
    query_randstrobes = strobealign.randstrobes_query(query, index_parameters)
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
