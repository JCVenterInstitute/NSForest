
from . import preprocessing as pp
from . import nsforesting
from . import evaluating as ev
from . import plotting as pl
import nsforest as ns

__all__ = ["ns", "pp", "nsforesting", "ev", "pl"]

__version__ = "4.1"

# import logging

# logging.basicConfig(filename="nsforest.log",
#                     level=logging.INFO,
#                     format="%(asctime)s - %(levelname)s - %(message)s",
#                     datefmt="%Y-%m-%d %H:%M:%S"  # Excludes milliseconds
# )

#     logger = logging.getLogger("nsforest")
#     logger.info("testing info")
#     logger.debug("Harmless debug Message")
#     logger.info("Just an information")
#     logger.warning("Its a Warning")
#     logger.error("Did you try to divide by zero")
#     logger.critical("Internet is down")