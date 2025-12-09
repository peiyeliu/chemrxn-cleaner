# chemrxn_cleaner/io/loader_registry.py
from __future__ import annotations

import logging
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Literal,
    Optional,
    Protocol,
    Sequence,
    overload,
    runtime_checkable,
)

from ..types import ReactionRecord

logger = logging.getLogger(__name__)


@runtime_checkable
class InputLoader(Protocol):
    def __call__(self, source: Any, /, **kwargs: Any) -> Iterable[ReactionRecord]:
        """Load reactions from a source into an iterable of ReactionRecords.

        Args:
            source: Data source such as a path string, file object, DataFrame,
                or list. Paths are preferred for CLI usability.
            **kwargs: Additional loader-specific options (e.g., delimiter,
                column mapping).

        Returns:
            Iterable of ``ReactionRecord`` instances.
        """
        ...


_INPUT_FORMAT_REGISTRY: Dict[str, InputLoader] = {}


class InputFormatError(ValueError):
    """Error raised when registering or retrieving an input format fails."""

    pass


def register_input_format(
    name: str,
    loader: InputLoader,
    *,
    overwrite: bool = False,
) -> None:
    """Register a new input format.

    Args:
        name: Input format name (e.g., ``"uspto"``, ``"ord"``, ``"csv"``).
            Prefer lowercase with underscores.
        loader: Callable with signature ``loader(source, **kwargs)`` returning
            an iterable of ``ReactionRecord`` instances.
        overwrite: When True, replace any existing loader registered under the
            same name.

    Raises:
        InputFormatError: If the name already exists and ``overwrite`` is False
            or the loader does not satisfy the protocol.
    """
    key = name.strip().lower()
    if not key:
        raise InputFormatError("Input format name must be a non-empty string.")

    if key in _INPUT_FORMAT_REGISTRY and not overwrite:
        raise InputFormatError(
            f"Input format '{key}' already exists. Use overwrite=True to replace it."
        )

    if not isinstance(loader, InputLoader):
        # Type check is intentionally loose but can help catch issues early
        raise InputFormatError(
            f"Loader for format '{key}' does not match InputLoader protocol. "
            f"Expected callable(source, **kwargs) -> Iterable[ReactionRecord]."
        )

    _INPUT_FORMAT_REGISTRY[key] = loader
    logger.debug("Registered input format '%s' (overwrite=%s)", key, overwrite)


def get_input_format(name: str) -> InputLoader:
    """Retrieve the loader registered for the given name.

    Args:
        name: Format name used during registration (e.g., ``"uspto"``).

    Returns:
        The loader callable for the specified format.

    Raises:
        InputFormatError: If the format is not registered.
    """
    key = name.strip().lower()
    try:
        return _INPUT_FORMAT_REGISTRY[key]
    except KeyError:
        available = ", ".join(sorted(_INPUT_FORMAT_REGISTRY.keys())) or "<none>"
        raise InputFormatError(
            f"Unknown input format '{key}'. Available formats: {available}"
        ) from None


@overload
def load_reactions(
    source: str, *, fmt: Literal["uspto"], keep_meta: bool = False
) -> List[ReactionRecord]: ...


@overload
def load_reactions(
    source: str,
    *,
    fmt: Literal["ord"],
    generate_if_missing: bool = True,
    allow_incomplete: bool = True,
    canonical: bool = True,
    meta_extractor: Optional[Callable[[Any], Dict[str, Any]]] = None,
) -> List[ReactionRecord]: ...


@overload
def load_reactions(
    source: str,
    *,
    fmt: Literal["csv"],
    reactant_columns: Sequence[str] = (),
    product_columns: Sequence[str] = (),
    reagent_columns: Optional[Sequence[str]] = None,
    reaction_smiles_column: Optional[str] = None,
    delimiter: str = ",",
    skip_lines: int = 0,
    mapper: Optional[
        Callable[[ReactionRecord, Dict[str, Any]], Optional[ReactionRecord]]
    ] = None,
) -> List[ReactionRecord]: ...


@overload
def load_reactions(
    source: str,
    *,
    fmt: Literal["json"],
    mapper: Callable[[Any], Optional[ReactionRecord]],
) -> List[ReactionRecord]: ...


@overload
def load_reactions(
    source: Any,
    *,
    fmt: str,
    **kwargs: Any,
) -> List[ReactionRecord]: ...


def load_reactions(
    source: Any,
    *,
    fmt: str,
    **kwargs: Any,
) -> List[ReactionRecord]:
    """Parse an external data source using a registered format loader.

    Args:
        source: External data source such as a file path, DataFrame, or list.
        fmt: Registered input format name (e.g., ``"uspto"``, ``"ord"``, ``"csv"``).
        **kwargs: Extra arguments forwarded to the selected loader.

    Returns:
        List of ``ReactionRecord`` instances produced by the loader.
    """
    logger.info("Loading reactions using format '%s'", fmt)
    loader = get_input_format(fmt)
    return list(loader(source, **kwargs))


def _register_builtin_formats() -> None:
    """Register built-in input formats on import."""
    from .loader import (
        load_csv,
        load_json,
        load_ord,
        load_uspto,
    )

    register_input_format("uspto", load_uspto, overwrite=True)
    register_input_format("ord", load_ord, overwrite=True)
    register_input_format("csv", load_csv, overwrite=True)
    register_input_format("json", load_json, overwrite=True)


# Optional: auto-register built-in formats when the module is imported
_register_builtin_formats()
