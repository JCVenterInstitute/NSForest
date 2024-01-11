
import argparse
def _parse_args(args, **kwargs): 
  parser = argparse.ArgumentParser(**kwargs,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("-a", metavar = "arguments", type=str, 
                      help="arguments csv, if you want to specify more nsforest arguments than just filename and cluster")
  
  parser.add_argument("-f", metavar = "filename", type=str, 
                      help="anndata object")
  
  parser.add_argument("-c", metavar = "cluster", type=str, 
                      help="specify the header name of clusters in adata.obs")
  
  return parser.parse_args()
