"""Tiny YAML config helper shared by the command-line scripts."""
from __future__ import annotations

import os

import yaml


def load_config(path: str) -> dict:
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    cfg.setdefault("paths", {})
    return cfg


def p(cfg: dict, key: str, *parts: str) -> str:
    """Resolve cfg['paths'][key] joined with ``parts``; creates nothing."""
    return os.path.join(cfg["paths"][key], *parts)


def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path
