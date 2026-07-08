"""Progress bars for the long engine runs.

Each stage's stdout is captured to the pipeline log while stderr (this bar) stays live on
the terminal.
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
def impact_progress(fort18_path, total_length_m, desc="linac5-8"):
    """Drive a tqdm bar (in metres of beam travel) from Impact-T's fort.18 while it runs."""
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
                            if len(parts) >= 2:
                                z = float(parts[1])      # col 2 = reference z [m] (col 1 is time [s])
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


def g4bl_progress(line_iter, total, desc="converter"):
    """Drive a tqdm bar from g4bl's `Event N Completed` stdout lines, re-yielding each line."""
    import re
    pat = re.compile(r"Event\s+(\d+)\s+Completed")
    with progress_bar(total=total, desc=desc, unit="ev") as bar:
        last = 0
        for line in line_iter:
            m = pat.search(line)
            if m:
                n = int(m.group(1))
                if n > last:
                    bar.update(n - last)
                    last = n
            yield line
