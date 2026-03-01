"""
Support for external tutorial sources (git, archive) in systemtests.
"""

import hashlib
import logging
import os
import shutil
import subprocess
import tarfile
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

# Cache directory for fetched tutorials. Can be overridden via PRECICE_EXTERNAL_CACHE_DIR env.
_DEFAULT_CACHE = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "precice-tutorials"
PRECICE_EXTERNAL_CACHE_DIR = Path(os.environ.get("PRECICE_EXTERNAL_CACHE_DIR", _DEFAULT_CACHE))


@dataclass
class TutorialSource:
    """Describes where a tutorial is sourced from."""

    type: str  # "local" | "git" | "archive"
    url: Optional[str] = None
    ref: Optional[str] = None
    subdir: Optional[str] = None

    @classmethod
    def local(cls) -> "TutorialSource":
        return cls(type="local")

    @classmethod
    def from_dict(cls, data: dict) -> "TutorialSource":
        if data is None or data.get("type") == "local":
            return cls.local()
        return cls(
            type=data["type"],
            url=data.get("url"),
            ref=data.get("ref"),
            subdir=data.get("subdir"),
        )


def _cache_key(prefix: str, url: str, ref: Optional[str] = None, subdir: Optional[str] = None) -> str:
    """Generate a short content-addressable cache key."""
    parts = [url]
    if ref:
        parts.append(ref)
    if subdir:
        parts.append(subdir)
    raw = f"{prefix}:{':'.join(parts)}"
    return hashlib.sha256(raw.encode()).hexdigest()[:16]


def fetch_git_repo(url: str, ref: str, cache_dir: Path, subdir: Optional[str] = None) -> Path:
    """
    Clone or update a git repository and return the path to the checkout.
    If subdir is given, returns the path to that subdirectory within the repo.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = _cache_key("git", url, ref, subdir)
    checkout = cache_dir / key

    if checkout.exists():
        try:
            subprocess.run(
                ["git", "-C", str(checkout), "fetch", "origin", ref, "--depth", "1"],
                check=True,
                capture_output=True,
                timeout=120,
            )
            subprocess.run(
                ["git", "-C", str(checkout), "checkout", "FETCH_HEAD"],
                check=True,
                capture_output=True,
                timeout=60,
            )
        except subprocess.CalledProcessError as e:
            logging.warning(f"Git fetch/checkout failed for {url}, recloning: {e}")
            shutil.rmtree(checkout, ignore_errors=True)

    if not checkout.exists():
        result = subprocess.run(
            ["git", "clone", "--depth", "1", "--branch", ref, url, str(checkout)],
            capture_output=True,
            text=True,
            timeout=300,
        )
        if result.returncode != 0:
            # Fallback: branch may not exist (e.g. repo uses develop/master instead of main)
            # Clone without branch, then fetch and checkout ref (with common aliases)
            shutil.rmtree(checkout, ignore_errors=True)
            logging.debug(
                f"git clone --branch {ref} failed ({result.stderr}), trying clone + fetch"
            )
            subprocess.run(
                ["git", "clone", "--depth", "1", url, str(checkout)],
                check=True,
                capture_output=True,
                timeout=300,
            )
            refs_to_try = [ref]
            if ref == "main":
                refs_to_try.extend(["develop", "master"])
            elif ref == "master":
                refs_to_try.extend(["main", "develop"])
            last_err = None
            for r in refs_to_try:
                res = subprocess.run(
                    ["git", "-C", str(checkout), "fetch", "origin", r, "--depth", "1"],
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
                if res.returncode == 0:
                    subprocess.run(
                        ["git", "-C", str(checkout), "checkout", "FETCH_HEAD"],
                        check=True,
                        capture_output=True,
                        timeout=60,
                    )
                    break
                last_err = res.stderr
            else:
                raise RuntimeError(
                    f"Could not fetch ref '{ref}' (tried {refs_to_try}): {last_err}"
                )

    if subdir:
        subpath = checkout / subdir
        if not subpath.is_dir():
            raise FileNotFoundError(f"Subdirectory {subdir} not found in {url} (ref {ref})")
        return subpath
    return checkout


def fetch_archive(url: str, cache_dir: Path, subdir: Optional[str] = None) -> Path:
    """
    Download and extract an archive (tar.gz, tar, zip) and return the path.
    """
    import urllib.request

    cache_dir.mkdir(parents=True, exist_ok=True)
    key = _cache_key("archive", url, subdir=subdir)
    extract_dir = cache_dir / key

    if extract_dir.exists():
        return extract_dir / subdir if subdir else extract_dir

    with tempfile.NamedTemporaryFile(delete=False, suffix=".tar.gz") as tmp:
        tmp_path = Path(tmp.name)
    try:
        logging.info(f"Downloading {url}")
        urllib.request.urlretrieve(url, tmp_path)

        extract_dir.mkdir(parents=True, exist_ok=True)
        if url.endswith(".tar.gz") or url.endswith(".tgz") or url.endswith(".tar"):
            with tarfile.open(tmp_path, "r:*") as tf:
                tf.extractall(extract_dir)
        else:
            import zipfile

            with zipfile.ZipFile(tmp_path, "r") as zf:
                zf.extractall(extract_dir)
    finally:
        tmp_path.unlink(missing_ok=True)

    if subdir:
        subpath = extract_dir / subdir
        if not subpath.is_dir():
            raise FileNotFoundError(f"Subdirectory {subdir} not found in {url}")
        return subpath
    return extract_dir


def resolve_tutorial_root(
    path: Path,
    source: TutorialSource,
    cache_dir: Path,
) -> Path:
    """
    Resolve the filesystem path to the tutorial root.

    For local sources, returns path as-is (already under PRECICE_TUTORIAL_DIR).
    For git/archive sources, fetches the repository/archive and returns the path
    to the tutorial directory. The tutorial name (path.name) is used as the
    subdirectory within the fetched content.
    """
    if source.type == "local":
        return path

    if source.type == "git":
        if not source.url or not source.ref:
            raise ValueError("git source requires 'url' and 'ref'")
        root = fetch_git_repo(source.url, source.ref, cache_dir, source.subdir)
        return root / path.name

    if source.type == "archive":
        if not source.url:
            raise ValueError("archive source requires 'url'")
        root = fetch_archive(source.url, cache_dir, source.subdir)
        return root / path.name

    raise ValueError(f"Unknown source type: {source.type}")
