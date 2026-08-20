"""Portable canonical command-log handling."""

from __future__ import annotations

import re
from pathlib import Path


def portable_command_log(text: str, repository_root: Path | str) -> str:
    """Replace the repository root in tool output and reject other local paths."""

    root = str(Path(repository_root).resolve())
    normalized = text.replace(root, "${REPOSITORY_ROOT}")
    workstation_path = re.search(r"(?:^|[\s:=('\"])(/(?:home|Users|tmp)/[^\s]*)", normalized)
    if workstation_path is not None:
        raise ValueError(f"Command log contains an absolute workstation path: {workstation_path.group(1)}")
    return normalized
