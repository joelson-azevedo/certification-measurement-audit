"""Portable release checks; original numerical fixture tests are separate."""
from pathlib import Path
import json
import pytest
import verify_release as release

def test_manifest_complete_and_hash_verified():assert release.manifest_check()>50

def test_claim_cells_and_contexts():
    unique,occurrences=release.claims_check()
    assert unique>=50 and occurrences>unique

def test_approved_pdf_figure_and_link_integrity():
    result=release.documentation_check()
    assert result['pdf_copies_verified']==2 and result['curated_figures_verified']==3

def test_data_reconstruction_and_classical_results():
    data,ctt=release.data_and_classical_check()
    assert data['raw_score_mismatches']==data['pass_fail_mismatches']==0
    assert sum(r['item_correlations_checked'] for r in ctt)==120

def test_selected_model_calculations_reproduce_frozen_results():
    result=release.model_calculation_check()
    assert result['parameter_refitting'] is False and result['bootstrap_intervals_recomputed'] is False

def test_manifest_rejects_tampering_and_unlisted_file(tmp_path):
    file=tmp_path/'evidence.csv';file.write_text('fixed\n')
    manifest=tmp_path/'RELEASE_MANIFEST.sha256'
    manifest.write_text(release.sha(file)+'  evidence.csv\n')
    assert release.manifest_check(tmp_path)==1
    file.write_text('changed\n')
    with pytest.raises(AssertionError):release.manifest_check(tmp_path)
    file.write_text('fixed\n');(tmp_path/'unlisted.txt').write_text('unexpected')
    with pytest.raises(AssertionError):release.manifest_check(tmp_path)

def test_manifest_rejects_parent_path(tmp_path):
    (tmp_path/'RELEASE_MANIFEST.sha256').write_text('0'*64+'  ../outside.txt\n')
    with pytest.raises(AssertionError):release.manifest_check(tmp_path)
