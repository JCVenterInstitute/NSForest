
import glob
import pandas as pd

# Defining paths
dir_path = "/expanse/projects/jcl125/nsforest/NSForest/"
h5ad = "gtex_granular"
output_folder = dir_path + "NSForest_outputs_" + h5ad + "/"

# Combining *_results.csv files
files = glob.glob(output_folder + "*/*_results.csv")
df = pd.DataFrame()
for file in files: 
    data = pd.read_csv(file)
    df = pd.concat([df, data])
print(dir_path + "results_" + h5ad + ".csv")
df.to_csv(dir_path + "results_" + h5ad + ".csv", index = False)

# Combining output files
files = glob.glob(output_folder + "outputs/*.out")
df = pd.DataFrame()
for file in files: 
    dictionary = {"filename": [file.split("outputs/")[-1]]}
    inputs = False
    prepare = False
    nsforest = False
    with open(file) as f:
        for line in f:
            if "Input values" in line: 
                inputs = True
            elif inputs: 
                # dictionary["arguments.csv"] = [line]
                dictionary["cluster"] = [line.split("'cluster_list': ['")[-1].split("'")[0]]
                inputs = False
            elif "Preparing data..." in line: 
                prepare = True
            elif prepare: 
                dictionary["time_prepare"] = [float(line.replace("-", "").replace(" ", "").split(".")[0])/60]
                prepare = False
            elif "running NSForest on all clusters" in line: 
                nsforest = True
            elif nsforest and "---" in line: 
                dictionary["time_nsforest"] = [float(line.replace("-", "").replace(" ", "").split(".")[0])/60]
        df = pd.concat([df, pd.DataFrame(dictionary)])
df = df.sort_values("time_nsforest", ascending = False)
df["time_prepare"] = [str(round(val, 1)) + " minutes" for val in df["time_prepare"]]
df["time_nsforest"] = [str(round(val, 1)) + " minutes" for val in df["time_nsforest"]]
print(dir_path + "output_" + h5ad + ".csv")
df.to_csv(dir_path + "output_" + h5ad + ".csv", index = False)
