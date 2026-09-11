from tests.official_data import fetch_official_data


def test_fetches_biopython_fixture_with_verified_checksum(tmp_path):
    path = fetch_official_data("biopython_dups.fasta", tmp_path)
    assert path.read_text().startswith(">")
    assert path.stat().st_size == 129
