# gmx_ffconv
A semi-automated force field converter for GROMACS

Usage:

gmx_ffconv ffmap [-h] -itp1 ITP1 -itp2 ITP2 -name NAME [--duplicate]
                        [--all_mappings]

options:
  -h, --help      show this help message and exit  
  -itp1 ITP1      First ITP file (path), corresponding to force field used in
                  .gro file  
  -itp2 ITP2      Second ITP file (path)  
  -name NAME      Name of the molecule  
  --duplicate     Skip graph matching, create a mapping where everything is
                  kept in same order. Useful for proteins, solvent and ions  
  --all_mappings  Obtain all mappings, not recommended  


gmx_ffconv groconv [-h] -name NAME [NAME ...] -nmol NMOL [NMOL ...]
                          -coordfile COORDFILE [-mapping_dir MAPPING_DIR]
                          -output OUTPUT  

options:
  -h, --help            show this help message and exit  
  -name NAME [NAME ...]
                        Molecule names separated by spaces, matching those used with ffmap  
  -nmol NMOL [NMOL ...]
                        Molecule counts separated by spaces  
  -coordfile COORDFILE  Input .gro file  
  -mapping_dir MAPPING_DIR
                        Directory containing mapping CSV files  
  -output OUTPUT        Output .gro file name

# *Example:*

System containing a membrane with 12 DPPC, 12 CHL, 12 DOPE, 2 K+, 2 Cl-, 1332 OPC. The initial force field is Amber, being converted to CHARMM.

First, obtain a mapping for DPPC, CHL, DOPE,K+, Cl-, OPC as follows:  

gmx_ffconv -itp1 PPPC_AMBER.itp -itp2 DPPC_CHARMM.itp -name DPPC  
gmx_ffconv -itp1 CHL_AMBER.itp -itp2 CHL1_CHARMM.itp -name CHL  
gmx_ffconv -itp1 OOPE_AMBER.itp -itp2 DOPE_CHARMM.itp -name DOPE  
gmx_ffconv -itp1 K+.itp -itp2 K+.itp --duplicate -name POT  
gmx_ffconv -itp1 Cl-.itp -itp2 Cl-.itp --duplicate -name CL  
gmx_ffconv -itp1 OPC.itp -itp2 OPC.itp --duplicate -name SOL  

Now, you can convert it to CHARMM to as follows:

gmx_ffconv groconv -coordfile MEMB_AMBER.gro -nmol 12 12 12 2 2 1332 -name DPPC CHL DOPE POT CL SOL -output MEMB_CHARMM.gro

I would recommend moving the mapping files into another directory, perhaps a user created database for their most commonly used mappings.   
As groconv does accept a mapping directory (e.g. AMBER_TO_CHARMM) as input, this enables users to create a reusable database for their most common mappings instead of having to always run ffmap.

Now, you should run ffmap again, with this time the order of itp files reversed to create a mapping from CHARMM to AMBER (ignoring the ones with --duplicate)  
After, you should run:

gmx_ffconv groconv -coordfile MEMB_CHARMM.gro -nmol 12 12 12 2 2 1332 -name DPPC CHL DOPE POT CL SOL -output backconv_MEMB_AMBER.gro

Validation:

To validate the mapping obtained, users should carry out a single-point energy calculation on the original coordinate file and the backconverted structure.


To do:
- perform automatically the backconversion, instead of requiring the user to first do a forward conversion and then a backwards one.

If you find my tool useful, please cite:
