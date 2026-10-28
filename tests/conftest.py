"""Keep synthetic test fixtures inside this independent analyst workspace."""

from pathlib import Path
from uuid import uuid4


def pytest_configure(config):
    if config.option.basetemp is None:
        workspace = Path(__file__).resolve().parents[1]
        target = (workspace / ".test_runs" / uuid4().hex).resolve()
        # pytest may remove an existing basetemp recursively: select a fresh
        # nonexistent leaf and verify containment before handing it to pytest.
        if not target.is_relative_to(workspace) or target.exists():
            raise RuntimeError("Unsafe test temporary directory")
        target.parent.mkdir(exist_ok=True)
        config.option.basetemp = str(target)
