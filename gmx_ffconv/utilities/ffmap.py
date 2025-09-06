from .utilities_ffmap.mapping_standard import run_ffmap_standard
from .utilities_ffmap.mapping_naming import run_ffmap_naming



def run_ffmap(args):
    if args.consistent_naming:
        run_ffmap_naming(args)
    else:
        run_ffmap_standard(args)