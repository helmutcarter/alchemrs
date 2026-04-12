from __future__ import annotations

import shutil
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
PYTHON_ROOT = REPO_ROOT / "python"
PACKAGE_ROOT = PYTHON_ROOT / "alchemrs"
BUILT_EXTENSION = REPO_ROOT / "target" / "debug" / "_alchemrs.dll"
STAGED_EXTENSION = PACKAGE_ROOT / "_alchemrs.pyd"


def pytest_sessionstart(session) -> None:
    if str(PYTHON_ROOT) not in sys.path:
        sys.path.insert(0, str(PYTHON_ROOT))

    if BUILT_EXTENSION.exists():
        if not STAGED_EXTENSION.exists() or BUILT_EXTENSION.stat().st_mtime > STAGED_EXTENSION.stat().st_mtime:
            shutil.copy2(BUILT_EXTENSION, STAGED_EXTENSION)
