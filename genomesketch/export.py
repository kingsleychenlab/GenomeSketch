"""CSV / JSON export helpers.

These use the standard library plus (optionally) pandas so that tabular data
can be written without any network access.
"""

from __future__ import annotations

import csv
import io
import json
from typing import Any, Iterable, List, Mapping


def to_json(data: Any, path: str | None = None, *, indent: int = 2) -> str:
    """Serialise ``data`` to JSON text, optionally writing it to ``path``."""
    text = json.dumps(data, indent=indent, default=_json_default)
    if path is not None:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
    return text


def _json_default(obj: Any) -> Any:
    # Support dataclasses and objects exposing as_dict().
    if hasattr(obj, "as_dict"):
        return obj.as_dict()
    if hasattr(obj, "__dict__"):
        return {k: v for k, v in vars(obj).items() if not k.startswith("_")}
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serialisable")


def to_csv(
    rows: Iterable[Mapping[str, Any]],
    path: str | None = None,
    *,
    fieldnames: List[str] | None = None,
) -> str:
    """Write an iterable of dict rows as CSV text, optionally to ``path``."""
    rows = list(rows)
    if fieldnames is None:
        fieldnames = []
        for row in rows:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)

    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        writer.writerow({k: row.get(k, "") for k in fieldnames})
    text = buffer.getvalue()

    if path is not None:
        with open(path, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
    return text
