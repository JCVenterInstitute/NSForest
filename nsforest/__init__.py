
import nsforest as ns
from . import preprocessing as pp
from . import nsforesting
from . import evaluating as ev
from . import plotting as pl

__all__ = ["ns", "pp", "nsforesting", "ev", "pl"]

__version__ = "4.0.0"

# # Figuring out ev
# from nsforest import utils
# from nsforest import evaluating as ev
# from nsforest.evaluating import _do_calc
# from nsforest.preprocessing._add_ann import preprocessing_medians
# results = pd.read_csv("test/test_results.csv")
# results["NSForest_markers"] = utils.str_to_list(results["NSForest_markers"])
# markers_dict = dict(zip(results["clusterName"], results["NSForest_markers"]))
# fractions_df = ns.ev.on_target_fraction(adata, markers_dict, "cluster", "medians_cluster", "test/", "test")
# results = results.merge(fractions_df, on = "clusterName", how = "left")
