from importlib import metadata as importlib_metadata

from .model import FreeSurferSeg


def get_version() -> str:
    try:
        return importlib_metadata.version(__name__)
    except importlib_metadata.PackageNotFoundError:  # pragma: no cover
        return "unknown"


__version__: str = get_version()

__all__ = ["FreeSurferSeg"]
