"""Offline release metadata and mixed-license checks; no statistical estimation."""
from pathlib import Path
import csv
import hashlib
import json
import sys
import tomllib

MIT_SHA256='ec66f5fb009c0ec29db46645a4c33960371a665cb96b8f9778de74e276b3d0dc'
AUTHOR='Joelson Azevêdo'
AI='''This independent measurement science project was conducted under the direction of Joelson Azevêdo, with AI assistance in software development, statistical implementation, automated testing, technical documentation, and editorial production.

Research objectives, methodological decisions, interpretation of findings, and publication decisions were reviewed and directed by the author.

The project includes numerical verification, automated testing, and reproducibility checks. These procedures do not constitute independent empirical validation of the certification decisions.'''
_GUARD_INSTALLED=False

def install_network_guard():
    """Deny Python socket creation, connect and DNS operations in this process."""
    global _GUARD_INSTALLED
    if _GUARD_INSTALLED:return
    def deny(event,args):
        if event.startswith('socket.'):
            raise RuntimeError('Network operation blocked by Tier 1 verification policy')
    sys.addaudithook(deny);_GUARD_INSTALLED=True

def citation_check(citation):
    assert citation['cff-version']=='1.2.0'
    assert citation['type']=='software'
    assert citation['title']=='Can You Defend Your Pass/Fail Decisions? An Independent Measurement Audit of a Professional Certification Examination'
    assert isinstance(citation['message'],str) and 'MIT' in citation['message'] and 'CC BY 4.0' in citation['message']
    assert citation['authors']==[{'given-names':'Joelson','family-names':'Azevêdo','affiliation':'Independent Researcher'}]
    assert citation['license']==['MIT','CC-BY-4.0']
    assert set(citation)=={'cff-version','message','type','title','authors','license','abstract','keywords'}
    assert isinstance(citation['abstract'],str) and 'synthetic' in citation['abstract']
    assert isinstance(citation['keywords'],list) and all(isinstance(k,str) and k for k in citation['keywords'])
    return True

def mit_check(root):
    p=Path(root)/'LICENSES/MIT.txt'
    assert hashlib.sha256(p.read_bytes()).hexdigest()==MIT_SHA256,'MIT text/notice differs from the locally verified exact grant'
    assert 'Copyright (c) 2026 '+AUTHOR in p.read_text(encoding='utf-8')

def disclosure_check(text):
    assert AI in text,'Approved author-directed AI disclosure missing or modified'

def license_for(path):
    p=Path(path)
    if p.suffix=='.py' or path in {'.gitignore','pyproject.toml','environment/requirements-quick.txt','LICENSES/MIT.txt'} or path.startswith('analysis/configs/'):return 'MIT'
    return 'CC-BY-4.0'

def metadata_check(root):
    root=Path(root)
    # CFF is encoded with JSON flow syntax: a valid YAML 1.2 subset, parsed with stdlib JSON.
    citation=json.loads((root/'CITATION.cff').read_text(encoding='utf-8'))
    citation_check(citation);mit_check(root)
    disclosure_check((root/'README.md').read_text(encoding='utf-8'))
    disclosure_check((root/'AI_ASSISTED_WORKFLOW_DISCLOSURE.md').read_text(encoding='utf-8'))
    with (root/'docs/LICENSE_SCOPE.csv').open(encoding='utf-8',newline='') as f:scope=list(csv.DictReader(f))
    names=[r['file'] for r in scope]
    static={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and not set(p.relative_to(root).parts)&{'verification','.test_runs','__pycache__','.pytest_cache','.venv'}}
    assert names==sorted(set(names)) and set(names)==static,'License inventory does not match exact static candidate'
    for r in scope:
        assert r['license']==license_for(r['file']) and r['copyright_holder']==AUTHOR
        assert r['approval'].startswith('Author Stage 15')
    policy=(root/'LICENSE_POLICY.md').read_text(encoding='utf-8')
    assert 'https://creativecommons.org/licenses/by/4.0/legalcode' in policy
    assert 'not relicensed' in policy and 'AWAITING SEPARATE EXPLICIT AUTHORIZATION' in policy
    cc=(root/'LICENSES/CC-BY-4.0-NOTICE.md').read_text(encoding='utf-8')
    assert 'not the canonical legal text' in cc
    metadata=tomllib.loads((root/'pyproject.toml').read_text(encoding='utf-8'))['tool']['audit_release']
    assert metadata['author']==AUTHOR and metadata['affiliation']=='Independent Researcher'
    assert metadata['intended_repository']=='https://github.com/joelson-azevedo/certification-measurement-audit'
    assert metadata['intended_visibility']=='public' and metadata['publication_authorization']=='awaiting separate explicit authorization'
    return dict(status='PASS',author=AUTHOR,affiliation='Independent Researcher',
        citation_parse='JSON-compatible YAML 1.2 flow syntax',cff_schema_oriented_checks='PASS; full formal CFF schema validator unavailable offline',
        licenses=['MIT','CC-BY-4.0'],license_scope_files=len(scope),canonical_mit_text_verified=True,
        cc_legal_text_bundled=False,approved_ai_disclosure_exact=True,publication_authorized=False)
