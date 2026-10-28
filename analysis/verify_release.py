"""Portable, network-free Tier 1 verification of this local review candidate."""
from pathlib import Path
import argparse
import hashlib
import importlib
import importlib.metadata
import json
import re
import os
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
from verify_metadata import metadata_check, install_network_guard
install_network_guard()
import numpy as np
import pandas as pd
from pypdf import PdfReader
from PIL import Image
from src.data_integrity import Findings,verify_manifest,validate_structure,validate_responses,reconstruct_histories
from src.scoring import make_long,reconstruct_scores,compare_scores
from src.ctt import form_matrices,corrected_item_total,cronbach_alpha,kr20
from src.linking_transformations import transform_parameters,probabilities
from src.score_distribution import normal_quadrature,marginal_distribution
from src.irt_observed_score_equating import equate_observed
from src.irt_true_score_equating import score_roots
from src.classification_probabilities import passing_curves
from src.decision_consistency import integrated_models,joint_classification

RUNTIME={'verification','.test_runs','__pycache__','.pytest_cache','.venv'}
DISCLOSURE='Independent blind analysis of a synthetic professional certification dataset.'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(path):return pd.read_csv(ROOT/path,float_precision='round_trip')
def static_files(root=ROOT):
    return {p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and not set(p.relative_to(root).parts)&RUNTIME and p.name!='RELEASE_MANIFEST.sha256'}

def manifest_check(root=ROOT):
    expected={}
    for line in (root/'RELEASE_MANIFEST.sha256').read_text().splitlines():
        digest,name=line.split('  ',1);p=(root/name).resolve()
        assert re.fullmatch('[0-9a-f]{64}',digest) and p.is_relative_to(root.resolve()) and p.is_file(),name
        assert name not in expected and sha(p)==digest,name
        expected[name]=digest
    assert set(expected)==static_files(root),(set(expected)^static_files(root))
    return len(expected)

def environment_check():
    versions={}
    for line in (ROOT/'environment/requirements-quick.txt').read_text().splitlines():
        if not line.strip() or line.startswith('#'):continue
        name,wanted=line.split('==');actual=importlib.metadata.version(name)
        assert actual==wanted,(name,wanted,actual)
        versions[name]=actual
    return versions

def source_value(source):
    p=ROOT/source['packaged_artifact'];field=source.get('packaged_field',source['field'])
    if p.suffix=='.csv':return pd.read_csv(p,float_precision='round_trip').iloc[int(source['row'])-1][field]
    if field.startswith('lines:'):
        a,b=map(int,field[6:].split('-'));return '\n'.join(p.read_text(encoding='utf-8').splitlines()[a-1:b])
    data=json.loads(p.read_text(encoding='utf-8'))
    if field=='$':return data
    for part in field.split('.'):data=data[int(part)] if isinstance(data,list) else data[part]
    return data

def claims_check():
    registry=json.loads((ROOT/'docs/claim_registry.json').read_text(encoding='utf-8'))
    index=read('docs/PUBLIC_CLAIM_EVIDENCE_MAP.csv')
    for ident,c in registry.items():
        for source in c['sources']:
            actual=source_value(source)
            assert actual==source['value'],(ident,actual,source['value'])
            if 'format' in source:
                assert format(float(actual)*source['multiplier'],source['format'])==c['display'],ident
    for row in index.itertuples():
        assert row.claim_id in registry
        file,lineno=row.published_location.rsplit(':',1)
        line=(ROOT/file).read_text(encoding='utf-8').splitlines()[int(lineno)-1]
        assert line==row.value_or_statement and '['+row.claim_id+']' in line
    # Check same numbers' different meanings using source row roles.
    def row_for(key):
        source=registry[key]['sources'][0]
        return read(source['packaged_artifact']).iloc[int(source['row'])-1]
    for key,form in [('Q09','A'),('Q10','B')]:
        row=row_for(key);assert row['form']==form and row.population=='first attempts'
    for key,form in [('Q07','A'),('Q08','B')]:
        row=row_for(key);assert row['form']==form and row.item_id=='ANC_18' and row.specification=='primary_rest_eap'
        assert registry[key]['sources'][0]['field']=='signed_reference_mean_difference' and registry[key]['value']>0
    assert registry['Q13']['sources'][0]['field']=='observed_inverse_B_continuous'
    for key in ['Q13','Q14','Q15']:
        row=row_for(key);assert row.specification=='S18_SL' and row.reference_distribution=='A_standard_normal' and row.reference_A_score==42
    row=row_for('Q16');assert row.scenario=='O' and row.pair_type=='cross_form' and row.cut_first==row.cut_second==42
    row=row_for('Q19');assert row.scenario=='O' and row.A_expected_score_location==42 and row.cut_A==row.cut_B==42
    return len(registry),len(index)

def data_and_classical_check():
    data=ROOT/'data/proposed';findings=Findings()
    assert verify_manifest(data,findings).match.all() and not findings.critical
    admin=pd.read_csv(data/'input/exam_administrations.csv',dtype=str,keep_default_na=False)
    meta=pd.read_csv(data/'input/item_metadata.csv',dtype=str,keep_default_na=False)
    key=pd.read_csv(data/'input/answer_key.csv',dtype=str,keep_default_na=False)
    checks=validate_structure(admin,meta,key,findings)
    validate_responses(admin,meta,findings);history=reconstruct_histories(admin,findings)
    assert checks.passed.all() and not findings.critical
    long=make_long(admin,meta,key,history);scored=reconstruct_scores(admin,long,history)
    compare=compare_scores(scored)
    expected=json.loads((ROOT/'results/stage1/audit_summary.json').read_text())
    got=dict(administrations=len(scored),first_attempts=int(scored.is_first_attempt.sum()),retakes=int((~scored.is_first_attempt).sum()),
             raw_score_mismatches=int(compare.raw_mismatch.sum()),pass_fail_mismatches=int(compare.pass_fail_mismatch.sum()),
             percent_score_mismatches=int(compare.percent_mismatch.sum()),omissions=int(long.omitted.sum()))
    assert all(v==expected[k] for k,v in got.items())
    assert admin.administration_id.str.fullmatch(r'ADM_\d{6}').all() and admin.candidate_id.str.fullmatch(r'CAN_\d{5}').all()
    assert set(admin.region)=={'North','South'}
    assert admin.candidate_id.nunique()==expected['unique_candidates']
    first=scored.loc[scored.is_first_attempt];long_first=long.loc[long.is_first_attempt]
    rates=read('results/stage1/operational_pass_rates.csv');stats=read('results/stage2/test_statistics_by_form.csv')
    items=read('results/stage2/ctt_item_statistics.csv');reproduced=[]
    for form in ['A','B']:
        sit,members,x,_=form_matrices(first,long_first,meta,form)
        r,status=corrected_item_total(x);old=items.loc[items.form.eq(form)].set_index('item_id').loc[members]
        np.testing.assert_allclose(r,old.corrected_item_total_r,atol=2e-14,rtol=0,equal_nan=True)
        original=stats.loc[stats.form.eq(form)].iloc[0]
        alpha,kr=cronbach_alpha(x),kr20(x)
        np.testing.assert_allclose([alpha,kr],[original.cronbach_alpha,original.kr20],atol=2e-14,rtol=0)
        old_rate=rates.loc[rates.population.eq('first attempts')&rates.form.eq(form)].iloc[0]
        rate=float(sit.audited_pass_fail.eq('PASS').mean())
        assert len(sit)==old_rate.administrations and abs(rate-old_rate.pass_rate)<1e-15
        reproduced.append(dict(form=form,first_attempts=len(sit),pass_rate=rate,alpha=alpha,kr20=kr,item_correlations_checked=len(r)))
    return got,reproduced

def model_calculation_check():
    parameters=read('results/stage4/irt_item_parameters.csv');links=read('results/stage7/linking_constants_by_method.csv')
    models={}
    for form in ['A','B']:
        p=parameters.loc[parameters.model_id.eq(form+'_full_signed_60_missing_ignored_q81')]
        assert len(p)==60 and p.item_id.is_unique and p.model_numerically_admissible.all()
        models[form]=(p.intercept_d.to_numpy(),p.discrimination_a.to_numpy())
    link=links.loc[links.anchor_specification.eq('S18')&links.method.eq('Stocking-Lord')&links.grid_name.eq('primary')].iloc[0]
    da,aa=models['A'];db0,ab0=models['B'];db,ab=transform_parameters(db0,ab0,link.scale_A,link.shift_B)
    grid=np.linspace(-3,3,23)
    np.testing.assert_allclose(probabilities(grid,db,ab),probabilities((grid-link.shift_B)/link.scale_A,db0,ab0),atol=2e-15,rtol=0)
    theta,w=normal_quadrature(321)
    distA=marginal_distribution(da,aa,theta,w);distB=marginal_distribution(db,ab,theta,w)
    inverse=float(equate_observed(np.array([42.]),distA['pmf'],distB['pmf'])['equated'][0])
    old=read('results/stage8/provisional_cut_correspondence.csv')
    old=old.loc[old.specification.eq('S18_SL')&old.reference_distribution.eq('A_standard_normal')&old.reference_A_score.eq(42)].iloc[0]
    assert abs(inverse-old.observed_inverse_B_continuous)<1e-9
    roots=score_roots(42,da,aa);assert len(roots)==1
    pa=passing_curves(np.array(roots),da,aa,[42])[0][42]
    pb=passing_curves(np.array(roots),db,ab,[42])[0][42]
    conditional=float(joint_classification(pa,pb)['agreement'][0])
    near=read('results/stage9/near_cut_classification_uncertainty.csv')
    near=near.loc[near.scenario.eq('O')&near.A_expected_score_location.eq(42)].iloc[0]
    assert abs(conditional-near.agreement)<1e-10
    frame,_,logs=integrated_models(da,aa,db,ab,link.scale_A,link.shift_B,'A_standard_normal',
                                  [dict(scenario='O',cut_B=42)],[321,641])
    global_value=float(frame.loc[frame.pair_type.eq('cross_form')].agreement.iloc[0])
    old=read('results/stage9/cross_form_decision_consistency.csv')
    old=old.loc[old.linking_method.eq('Stocking-Lord')&old.reference_population.eq('A_standard_normal')&old.scenario.eq('O')].iloc[0]
    assert abs(global_value-old.agreement)<1e-10
    return dict(observed_inverse_B_at_A42=inverse,conditional_agreement_at_A_expected42=conditional,
                integrated_cross_form_agreement_original_rules=global_value,parameter_refitting=False,
                bootstrap_intervals_recomputed=False,model_recomputed_from='frozen q81 item parameters and provisional S18 link')

def documentation_check():
    provenance=read('docs/FILE_PROVENANCE.csv')
    for r in provenance.itertuples():
        assert sha(ROOT/r.packaged_file)==r.packaged_sha256
        if r.change=='none':assert r.original_sha256==r.packaged_sha256
    for name,count in [('EXECUTIVE_MEASUREMENT_BRIEF',5),('EXECUTIVE_ONE_PAGE_SUMMARY',1)]:
        p=ROOT/'reports'/(name+'.pdf');reader=PdfReader(p)
        assert len(reader.pages)==count and reader.metadata.subject==DISCLOSURE
        text=' '.join('\n'.join(page.extract_text() for page in reader.pages).split())
        assert DISCLOSURE in text and 'not an authorized conversion or recommended new passing cut' in text
        assert not re.search(r'[A-Z]:[\\/]|[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}',str(reader.metadata)+text)
    for p in (ROOT/'figures').glob('*.png'):
        with Image.open(p) as im:im.verify()
    for p in ROOT.rglob('*.md'):
        if set(p.relative_to(ROOT).parts)&RUNTIME:continue
        text=p.read_text(encoding='utf-8')
        links=re.findall(r'\]\(([^)]+)\)',text)+re.findall(r'^\[[^]]+\]:\s+(\S+)',text,re.M)
        for target in links:
            if re.match(r'https?://',target):continue
            name=target.split('#')[0]
            if not name:continue
            dest=(p.parent/name).resolve()
            assert dest.is_relative_to(ROOT) and dest.exists(),(str(p.relative_to(ROOT)),target)
        assert not re.search(r'[A-Z]:[\\/]|\{\{|\bTODO\b|\bTBD\b',text)
    for n in ['README.md','CASE_STUDY.md','METHODOLOGY.md']:
        assert 'synthetic' in (ROOT/n).read_text(encoding='utf-8').lower()
    case=(ROOT/'CASE_STUDY.md').read_text(encoding='utf-8')
    words=len(re.findall(r"\b[\w'-]+\b",re.sub(r'^\[[A-Z]+\d+\]:.*$', '',case,flags=re.M)))
    assert 3000<=words<=4500,words
    return dict(case_study_words=words,code_provenance_rows=len(provenance),pdf_copies_verified=2,curated_figures_verified=3)

def verify(run_tests=False):
    install_network_guard()
    started=time.perf_counter();files=manifest_check();versions=environment_check();release_metadata=metadata_check(ROOT)
    registry_count,occurrences=claims_check();data,ctt=data_and_classical_check();model=model_calculation_check();docs=documentation_check()
    modules=[]
    for p in sorted((ROOT/'analysis/src').glob('*.py')):
        if p.stem=='__init__':continue
        module=importlib.import_module('src.'+p.stem)
        assert Path(module.__file__).resolve().is_relative_to(ROOT/'analysis')
        modules.append(p.stem)
    result=dict(status='PASS',tier='Tier 1 only',release_metadata=release_metadata,python_network_operations_disabled=True,manifest_files=files,claim_definitions=registry_count,
                claim_occurrences=occurrences,data_reconstruction=data,classical_reproduction=ctt,
                selected_model_calculations=model,documentation=docs,dependencies=versions,imported_modules=modules,
                full_analytical_replay_available=False,elapsed_seconds=time.perf_counter()-started)
    target=ROOT/'verification';target.mkdir(exist_ok=True)
    if run_tests:
        os.environ['PYTEST_DISABLE_PLUGIN_AUTOLOAD']='1'
        import pytest
        previous=Path.cwd()
        try:
            os.chdir(ROOT)
            returncode=pytest.main(['-q','-p','no:cacheprovider','--junitxml='+str(target/'test_results.xml')])
        finally:os.chdir(previous)
        assert returncode==0,'Packaged tests failed'
        result['tests_returncode']=int(returncode)
        result['elapsed_seconds']=time.perf_counter()-started
    (target/'tier1_result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run-tests',action='store_true');args=parser.parse_args()
    print(json.dumps(verify(args.run_tests),indent=2))
