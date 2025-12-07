import pathlib
import re

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]


def read_package_version() -> str:
    init_text = (REPO_ROOT / "chemrxn_cleaner" / "__init__.py").read_text(
        encoding="utf-8"
    )
    match = re.search(r'__version__\s*=\s*[\'"]([^\'"]+)[\'"]', init_text)
    assert match, "Unable to find __version__ in chemrxn_cleaner/__init__.py"
    return match.group(1)


def read_citation_version() -> str:
    citation_text = (REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r'^version:\s*"([^"]+)"$', citation_text, flags=re.MULTILINE)
    assert match, "Unable to find version field in CITATION.cff"
    return match.group(1)


def test_citation_version_matches_package() -> None:
    assert read_citation_version() == read_package_version()
