from pathlib import Path

import scanpy as sc
import yaml


def get_root_dir():
    """walk up the directories to find the root directory"""
    d = Path.cwd()
    while not (d / ".git").is_dir():
        d = d.parent
    return d


def load(path=None):
    if not path:
        root = get_root_dir()
        path = root / "workflow" / "config.yml"

    with path.open("r") as f:
        return yaml.safe_load(f)


def setup_scanpy(outs):
    sc.settings.figdir = outs
    sc.settings.verbosity = 0
