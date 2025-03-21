
import nsforest as ns
from . import preprocessing as pp
from . import nsforesting
from . import evaluating as ev
from . import plotting as pl

__version__ = "4.1"
NSFOREST_VERSION = __version__

__all__ = ["ns", "pp", "nsforesting", "ev", "pl", NSFOREST_VERSION]
