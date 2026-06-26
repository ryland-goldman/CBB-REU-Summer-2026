"""Build a per-eval LINACSIM_OUT_DIR sandbox: an isolated config/ copy + empty logs/,
with the shared fieldmaps/ symlinked back to the repo. See docs/per-eval-isolation-plan.md.
"""

import os
import shutil

from sim.helpers.tools import REPO_ROOT


def _link(target, link):                           # idempotent symlink
    if os.path.islink(link) or os.path.exists(link):
        return
    os.symlink(target, link)


def make_out_dir(out_dir, src_root=REPO_ROOT):
    """Populate a LINACSIM_OUT_DIR sandbox: own config/ copy + empty logs/, fieldmaps shared.
    NO-OP when out_dir == src_root (the unset / REPO_ROOT case -- never sandbox the repo itself)."""
    if os.path.abspath(out_dir) == os.path.abspath(src_root):
        return                                     # plain run: leave repo config/logs/fieldmaps as-is
    os.makedirs(out_dir, exist_ok=True)
    shutil.copytree(f"{src_root}/config", f"{out_dir}/config", dirs_exist_ok=True)  # isolated copy
    os.makedirs(f"{out_dir}/logs", exist_ok=True)
    _link(f"{src_root}/fieldmaps", f"{out_dir}/fieldmaps")     # shared maps, NOT redirected
