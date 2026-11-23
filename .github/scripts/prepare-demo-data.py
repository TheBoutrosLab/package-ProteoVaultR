#!/usr/bin/env python3
"""Prepare checksum-verified public demo databases for documentation builds."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import time
from typing import TypedDict
import urllib.error
import urllib.request


class Dataset(TypedDict):
    """One pinned public demo database recorded in the manifest."""

    id: str
    filename: str
    url: str
    sha256: str
    bytes: int


class Manifest(TypedDict):
    """The demo datasets required for documentation builds."""

    datasets: list[Dataset]


def valid_file(path: Path, entry: Dataset) -> bool:
    """Check both the expected size and SHA-256 of a cached database."""
    if not path.is_file() or path.stat().st_size != entry["bytes"]:
        return False
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return digest == entry["sha256"]


def prepare_file(entry: Dataset, directory: Path) -> None:
    """Reuse verified data or atomically replace it with a verified download."""
    destination = directory / entry["filename"]
    if valid_file(destination, entry):
        print(f"Verified cached file: {destination.name}", flush=True)
        return
    for attempt in range(3):
        temporary: Path | None = None
        try:
            print(f"Downloading {entry['url']} (attempt {attempt + 1})", flush=True)
            with tempfile.NamedTemporaryFile(dir=directory, suffix=".part", delete=False) as output:
                temporary = Path(output.name)
                with urllib.request.urlopen(entry["url"], timeout=60) as response:
                    shutil.copyfileobj(response, output, length=1024 * 1024)
            if not valid_file(temporary, entry):
                raise ValueError(f"Size or SHA-256 mismatch: {entry['filename']}")
            os.replace(temporary, destination)
            print(f"Verified downloaded file: {destination.name}", flush=True)
            return
        except (OSError, urllib.error.URLError, ValueError):
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)


def main() -> None:
    """Prepare every required input listed in the dataset manifest."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(".github/demo-data.json")
    )
    parser.add_argument(
        "--directory",
        type=Path,
        default=Path("vignettes/demo-data")
    )
    args = parser.parse_args()
    manifest: Manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    entries = manifest["datasets"]

    for entry in entries:
        if Path(entry["filename"]).name != entry["filename"]:
            raise ValueError("Dataset filenames must not contain a directory")

    args.directory.mkdir(parents=True, exist_ok=True)
    for entry in entries:
        prepare_file(entry, args.directory)


if __name__ == "__main__":
    main()
