from __future__ import annotations

import importlib.machinery
import shutil
import sys
import sysconfig
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
PYTHON_ROOT = REPO_ROOT / "python"
PACKAGE_ROOT = PYTHON_ROOT / "alchemrs"


def pytest_sessionstart(session) -> None:
    if str(PYTHON_ROOT) not in sys.path:
        sys.path.insert(0, str(PYTHON_ROOT))

    built_extension = find_built_extension()
    if built_extension is not None:
        staged_extension = PACKAGE_ROOT / f"_alchemrs{extension_suffix()}"
        if (
            not staged_extension.exists()
            or built_extension.stat().st_mtime > staged_extension.stat().st_mtime
        ):
            shutil.copy2(built_extension, staged_extension)


def find_built_extension() -> Path | None:
    debug_dir = REPO_ROOT / "target" / "debug"
    suffixes = tuple(set(importlib.machinery.EXTENSION_SUFFIXES) | {".dll"})
    candidates = []
    for pattern in ("_alchemrs*", "lib_alchemrs*"):
        candidates.extend(debug_dir.glob(pattern))

    candidates = [
        path for path in candidates if path.is_file() and path.name.endswith(suffixes)
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda path: path.stat().st_mtime)


def extension_suffix() -> str:
    return sysconfig.get_config_var("EXT_SUFFIX") or importlib.machinery.EXTENSION_SUFFIXES[0]
