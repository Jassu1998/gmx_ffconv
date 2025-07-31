
import argparse
from pathlib import Path
from gmx_ffconv.utilities.ffmap import run_ffmap
from gmx_ffconv.utilities.groconv import run_groconv

def main():
    parser = argparse.ArgumentParser(description="Converts an all-atom system between 2 force fields")
    subparsers = parser.add_subparsers(dest="command", required=True)
    # --- ffmap subcommand ---
    ffmap_parser = subparsers.add_parser("ffmap", help="Finds mapping between two force fields for each molecule type")
    ffmap_parser.add_argument("-itp1", required=True, help="First ITP file (path), corresponding to force field used in .gro file",type=Path)
    ffmap_parser.add_argument("-itp2", required=True, help="Second ITP file (path)", type=Path)
    ffmap_parser.add_argument("-name", required=True, help="Name of the molecule, does not need to match itp files")
    ffmap_parser.add_argument("--duplicate", action="store_true",
                              help="Skip graph matching, create a mapping where everything is kept in same order. Useful for proteins, solvent,ions ")
    ffmap_parser.add_argument("--all_mappings", action="store_true",
                              help="Obtain all mappings, not recommended")
    ffmap_parser.add_argument("--validate", action="store_true",
                              help="Carry out conversion in both directions")
    ffmap_parser.set_defaults(func=run_ffmap)


    # --- groconv subcommand ---
    groconv_parser =subparsers.add_parser("groconv", help="Converts gro to match new topology")
    groconv_parser.add_argument("-name", nargs="+",required=True, help="Molecule names separated by spaces")
    groconv_parser.add_argument("-nmol",nargs="+", required=True, help="Molecule counts separated by spaces")
    groconv_parser.add_argument("-coordfile", required=True, help="Input .gro file")
    groconv_parser.add_argument("-mapping_dir", default=".", help="Directory containing mapping CSV files")
    groconv_parser.add_argument("-output", required=True, help="Output .gro file name")
    groconv_parser.add_argument("--validate", action="store_true", help="Generate back-converted structure")
    groconv_parser.add_argument("--norename",action="store_true", help="Do not rename the reordered gro")


    args = parser.parse_args()
    if args.command is None:
        parser._print_help()
        parser.exit(1)

    if args.command == "ffmap":
        run_ffmap(args)
    elif args.command == "groconv":
        run_groconv(args)

if __name__ == "__main__":
    main()