"""Release regression checks reject unapproved citation or grant changes."""
from pathlib import Path
import copy
import json
import pytest
import verify_metadata as metadata

ROOT=Path(__file__).resolve().parents[1]

def test_complete_approved_release_metadata():
    result=metadata.metadata_check(ROOT)
    assert result['status']=='PASS' and result['publication_authorized'] is False

@pytest.mark.parametrize('change',['wrong_author','invented_version','blanket_mit'])
def test_citation_rejects_unapproved_metadata(change):
    citation=json.loads((ROOT/'CITATION.cff').read_text(encoding='utf-8'))
    citation=copy.deepcopy(citation)
    if change=='wrong_author':citation['authors'][0]['family-names']='Unapproved'
    elif change=='invented_version':citation['version']='1.0.0'
    else:citation['license']=['MIT']
    with pytest.raises(AssertionError):metadata.citation_check(citation)

def test_mit_legal_wording_tamper_is_rejected(tmp_path):
    target=tmp_path/'LICENSES/MIT.txt';target.parent.mkdir()
    original=(ROOT/'LICENSES/MIT.txt').read_bytes();target.write_bytes(original)
    metadata.mit_check(tmp_path)
    target.write_bytes(original.replace(b'free of charge',b'for a fee'))
    with pytest.raises(AssertionError):metadata.mit_check(tmp_path)

def test_independent_empirical_validation_substitution_is_rejected():
    altered=metadata.AI.replace('do not constitute independent empirical validation','constitute independent empirical validation')
    with pytest.raises(AssertionError):metadata.disclosure_check(altered)
