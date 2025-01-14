# ___Version___ #
from pkg_resources import get_distribution

try:
    __version__ = get_distribution("make_prg").version
except:
    __version__ = "0.3.0"

__all__ = ["from_msa", "subcommands", "io_utils", "seq_utils"]
