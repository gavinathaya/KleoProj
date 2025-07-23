"""
Miscellaneous Functions (:mod:`kleopy.miscellaneous`)
==========================================================
Miscellaneous functions used in KleoPy.

Note: Functions in this module does NOT depend on other modules in KleoPy.

Functions
---------

progress_bar(step, total_steps, fill_char, width) -> None      : Shows progress bar of a running process.
"""
import sys

def progress_bar(step, total_steps, *, fill_char='━', width=40):
    """
    Print a progress bar to the console.

    Parameters
    ----------
    """
    # Filler characters: ━ ╍ █
    # Names in order: "U+2501" "U+2505" "U+2588"
    # Edge characters: ╢ ╟
    # Names in order: "U+2562" "U+2560"
    progress = f"{step+1}/{total_steps} {'\033[92m'}╢{fill_char*((step+1)*width//total_steps):<{width}}╟{'\033[0m'}"
    
    #Line overwrite logic
    if step + 1 == total_steps:
        sys.stdout.write('\r' + progress + '\n')
    else:
        sys.stdout.write('\r' + progress)
    sys.stdout.flush()