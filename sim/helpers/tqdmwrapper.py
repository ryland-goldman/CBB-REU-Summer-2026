"""Progress bars for the long engine runs.

Each stage runs in its own subprocess whose stdout is captured to the pipeline log by
main.py while stderr stays on the terminal -- so a tqdm bar written to stderr shows live
while the engine's verbose stdout lands in the log. WarpX stages get their bar from
lume-warpx's `w.run(progress=...)`; Impact-T (linac4-8) drives the bar below from a
background poll of its `fort.18` longitudinal-position output.
"""

import contextlib
import sys
import threading
import time


@contextlib.contextmanager
def progress_bar(total=None, desc="", unit="it"):
    """A tqdm bar on the terminal (stderr). Fills to `total` on exit so it never hangs short."""
    from tqdm import tqdm
    bar = tqdm(total=total, desc=desc, unit=unit, ncols=88, leave=True,
               file=sys.stderr, disable=not sys.stderr.isatty())
    try:
        yield bar
    finally:
        if bar.total and bar.n < bar.total:
            bar.n = bar.total
            bar.refresh()
        bar.close()


@contextlib.contextmanager
def impact_progress(fort18_path, total_length_m, desc="linac4-8"):
    """Drive a tqdm bar (in metres of beam travel) from Impact-T's fort.18 while it runs.

    Impact-T appends the bunch centroid z to fort.18 each step; a daemon thread polls the
    last value and advances the bar. Wrap I.run() in this context.
    """
    import os
    with progress_bar(total=round(total_length_m, 3), desc=desc, unit="m") as bar:
        stop = threading.Event()

        def _poll():
            last = 0.0
            while not stop.is_set():
                try:
                    if os.path.isfile(fort18_path):
                        with open(fort18_path, "rb") as fh:
                            fh.seek(0, os.SEEK_END)
                            size = fh.tell()
                            fh.seek(max(0, size - 4096))
                            tail = fh.read().decode("ascii", "replace").splitlines()
                        for line in reversed(tail):
                            parts = line.split()
                            if parts:
                                z = float(parts[0])      # col 1 = reference z [m]
                                if z > last:
                                    bar.update(z - last)
                                    last = z
                                break
                except Exception:
                    pass
                time.sleep(0.3)

        t = threading.Thread(target=_poll, daemon=True)
        t.start()
        try:
            yield bar
        finally:
            stop.set()
            t.join(timeout=1.0)
