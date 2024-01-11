
import argparse

def _parse_args(args, **kwargs): 

  parser = argparse.ArgumentParser(**kwargs,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)

  parser.add_argument("filename",
                      type=str,
                      help="")
  
  parser.add_argument("cluster",
                      type=str,
                      help="")
  # parser.add_argument("-c", 
  #                     type = int, nargs = '?', 
  #                     help = "test")
  # parser.add_argument('-f', 
  #                     type = str, nargs='?',
  #                     help = "MedPC file to be parsed")
  # parser.add_argument('-d', '--directory')
  # parser.add_argument('-o', '--stdout')

  return parser.parse_args()