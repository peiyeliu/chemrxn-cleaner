"""Keep CITATION.cff in sync with the package version."""

from __future__ import annotations

import pathlib
import re
from typing import Tuple

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
PACKAGE_INIT = REPO_ROOT / "chemrxn_cleaner" / "__init__.py"
CITATION_FILE = REPO_ROOT / "CITATION.cff"


def read_package_version() -> str:
    """Read __version__ from chemrxn_cleaner/__init__.py without importing."""
    init_text = PACKAGE_INIT.read_text(encoding="utf-8")
    match = re.search(r'__version__\s*=\s*[\'"]([^\'"]+)[\'"]', init_text)
    if not match:
        raise RuntimeError("Unable to find __version__ in chemrxn_cleaner/__init__.py")
    return match.group(1)


def update_citation_version(new_version: str) -> Tuple[bool, str]:
    """
    Replace the version field in CITATION.cff.

    Returns (updated, message).
    """
    citation_text = CITATION_FILE.read_text(encoding="utf-8")
    updated_text, count = re.subn(
        r"^version:\s*\".*\"$",
        f'version: "{new_version}"',
        citation_text,
        flags=re.MULTILINE,
    )
    if count == 0:
        raise RuntimeError("No version field found in CITATION.cff")
    if citation_text == updated_text:
        return False, f"CITATION.cff already at version {new_version}"

    CITATION_FILE.write_text(updated_text, encoding="utf-8")
    return True, f"Updated CITATION.cff to version {new_version}"


def main() -> None:
    new_version = read_package_version()
    updated, message = update_citation_version(new_version)
    print(message)
    if not updated:
        return


if __name__ == "__main__":
    main()
