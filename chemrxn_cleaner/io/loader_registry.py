# chemrxn_cleaner/io/loader_registry.py
from __future__ import annotations

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
    Tuple,
    Union,
    overload,
    runtime_checkable,
)

from ..types import ReactionRecord


@runtime_checkable
class InputLoader(Protocol):
    def __call__(
        self, source: Any, /, **kwargs: Any
    ) -> Iterable[Union[ReactionRecord, Tuple[str, Dict[str, Any]]]]:
        """
        Parameters
        ----------
        source:
            Any agreed-upon data source. Can be a path string, file object, DataFrame,
            list, etc.
            Recommendation: most loaders should use a str path for easy CLI use.
        **kwargs:
            Extra configuration options (e.g., delimiter, column mapping).

        Returns
        -------
        Iterable[Union[ReactionRecord, Tuple[str, Dict[str, Any]]]]
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
    """
    Register a new input format.

    Parameters
    ----------
    name:
        Input format name, e.g., "uspto", "ord", "csv", "my_lab_format".
        Suggested style: lowercase letters and underscores.
    loader:
        A callable with signature similar to:
            loader(source, **kwargs) -> Iterable[Tuple[str, Dict[str, Any]]]
    overwrite:
        If False (default), raise when name already exists.
        If True, overwrite an existing loader with the same name.

    Raises
    ------
    InputFormatError
        When name already exists and overwrite=False
        or loader does not match the protocol.
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
            f"Expected callable(source, **kwargs) -> "
            f"Iterable[ReactionRecord | Tuple[str, Dict[str, Any]]]."
        )

    _INPUT_FORMAT_REGISTRY[key] = loader


def get_input_format(name: str) -> InputLoader:
    """
    Retrieve the loader registered for the given name.

    Parameters
    ----------
    name:
        The format name used during registration, e.g., "uspto".

    Returns
    -------
    InputLoader

    Raises
    ------
    InputFormatError
        When the format is not registered.
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
) -> List[Union[ReactionRecord, Tuple[str, Dict[str, Any]]]]: ...


def load_reactions(
    source: Any,
    *,
    fmt: str,
    **kwargs: Any,
) -> List[Union[ReactionRecord, Tuple[str, Dict[str, Any]]]]:
    """
    Parse an external data source into ReactionRecord objects
    or (reaction_smiles, metadata) tuples using the given format.

    Parameters
    ----------
    source:
        External data source, e.g.:
        - Path to a text file of USPTO reaction SMILES
        - Path to an ORD JSON file
        - Path to a CSV file
        - An in-memory DataFrame, list, etc.
    fmt:
        Input format name, e.g., "uspto", "ord", "csv".
    **kwargs:
        Extra arguments forwarded to the loader.

    Returns
    -------
    List[Union[ReactionRecord, Tuple[str, Dict[str, Any]]]]
    """
    loader = get_input_format(fmt)
    return list(loader(source, **kwargs))


def _register_builtin_formats() -> None:
    """
    Register built-in input formats when the package is imported.
    Can be called in __init__.py: loader_registry._register_builtin_formats()
    """
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
