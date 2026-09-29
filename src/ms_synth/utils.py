"""
utils.py
========
Utility helpers: logging setup, reproducibility seeding, YAML/JSON
config loading, and common I/O conveniences.

All other modules import ``get_logger`` and ``load_config`` from here.
"""

from __future__ import annotations

import json
import logging
import os
import random
import sys
from pathlib import Path
from typing import Any

import numpy as np
import yaml


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

_LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def setup_logging(
    level: str = "INFO",
    log_file: str | Path | None = None,
) -> None:
    """Configure the root logger for the framework.

    Parameters
    ----------
    level : str
        One of "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL".
    log_file : str or Path, optional
        If provided, a FileHandler is also added (in addition to stdout).

    Notes
    -----
    Call this once at the start of ``main_pipeline.py`` before importing
    any other framework module so that every module's logger inherits the
    handler.
    """
    numeric_level = getattr(logging, level.upper(), logging.INFO)
    handlers: list[logging.Handler] = [
        logging.StreamHandler(sys.stdout),
    ]
    if log_file is not None:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file, mode="a"))

    logging.basicConfig(
        level=numeric_level,
        format=_LOG_FORMAT,
        datefmt=_DATE_FORMAT,
        handlers=handlers,
        force=True,
    )


def get_logger(name: str) -> logging.Logger:
    """Return a named logger (child of the root logger).

    Parameters
    ----------
    name : str
        Typically ``__name__`` of the calling module.

    Returns
    -------
    logging.Logger
    """
    return logging.getLogger(name)


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

def set_seed(seed: int = 42) -> None:
    """Set random seeds for Python, NumPy, and (if available) PyTorch.

    Parameters
    ----------
    seed : int
        Master random seed.  Stored in ``os.environ["MS_PERC_SEED"]`` so
        subprocesses can read it.
    """
    os.environ["MS_PERC_SEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    # Optional: torch
    try:
        import torch  # type: ignore
        torch.manual_seed(seed)
    except ImportError:
        pass
    logger = get_logger(__name__)
    logger.debug("Global seed set to %d", seed)


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_config(path: str | Path) -> dict[str, Any]:
    """Load a YAML or JSON configuration file.

    Parameters
    ----------
    path : str or Path
        Path to the config file.  Extension determines the parser
        (``.yaml`` / ``.yml`` → PyYAML; ``.json`` → json).

    Returns
    -------
    dict
        Parsed configuration dictionary.

    Raises
    ------
    FileNotFoundError
        If the config file does not exist.
    ValueError
        If the file extension is not recognised.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    suffix = path.suffix.lower()
    if suffix in {".yaml", ".yml"}:
        with open(path, "r", encoding="utf-8") as fh:
            cfg = yaml.safe_load(fh)
    elif suffix == ".json":
        with open(path, "r", encoding="utf-8") as fh:
            cfg = json.load(fh)
    else:
        raise ValueError(f"Unsupported config format: '{suffix}'.  Use .yaml or .json")

    logger = get_logger(__name__)
    logger.info("Loaded config from %s", path)
    return cfg


def ensure_dir(path: str | Path) -> Path:
    """Create a directory (and parents) if it does not exist.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    Path
        The (possibly newly created) directory.
    """
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def save_dict_json(data: dict[str, Any], path: str | Path) -> None:
    """Serialise a dictionary to a JSON file (numpy-safe).

    Parameters
    ----------
    data : dict
    path : str or Path
    """
    path = Path(path)
    ensure_dir(path.parent)

    def _default(obj: Any) -> Any:
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        raise TypeError(f"Object of type {type(obj)} is not JSON serialisable")

    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, default=_default)
