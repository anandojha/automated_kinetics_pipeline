from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from openeye import oechem
from rdkit import Chem
import argparse

def protein_ligand_pdbs(protein_ligand_file, ligand_resname, ligand_file="ligand.pdb", protein_file="protein.pdb"):
    """
    Splits a PDB file containing a protein-ligand complex into two separate PDB files: 
    one for the ligand and one for the protein.

    Parameters:
    protein_ligand_file (str): Path to the input PDB file containing both the protein and ligand.
    ligand_resname (str): The residue name of the ligand in the PDB file.
    ligand_file (str): Path to the output PDB file where the ligand atoms will be saved.
    protein_file (str): Path to the output PDB file where the protein atoms will be saved.
    
    Functionality:
    - The function reads the input `protein_ligand_file` line by line.
    - Only lines that start with "ATOM" or "HETATM" are processed, ensuring that only atom data is included in the output.
    - The function checks the residue name (located in columns 18-20 of the PDB format).
    - If the residue name matches `ligand_resname`, the line is written to the `ligand_file`.
    - If the residue name does not match `ligand_resname`, it is assumed to be part of the protein, and the line is written to the `protein_file`.
    
    """
    with open(protein_ligand_file, 'r') as protein_ligand_file:
        with open(ligand_file, 'w') as ligand_file, open(protein_file, 'w') as protein_file:
            for line in protein_ligand_file:
                # Process only lines that start with "ATOM" or "HETATM"
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    residue_name = line[17:20].strip()
                    if residue_name == ligand_resname:
                        ligand_file.write(line)
                    else:
                        protein_file.write(line)

def pdb_to_sdf(ligand_file="ligand.pdb", sdf_file="ligand.sdf"):
    """
    Converts a PDB file to an SDF file using OpenEye's OEChem toolkit.
    
    Parameters:
    - ligand_file: str, path to the input PDB file
    - sdf_file: str, path to the output SDF file
    """
    # Create input and output molecule streams
    ifs = oechem.oemolistream()
    ofs = oechem.oemolostream()
    # Open input PDB file and output SDF file
    if not ifs.open(ligand_file):
        oechem.OEThrow.Fatal(f"Unable to open the input PDB file: {ligand_file}")
    if not ofs.open(sdf_file):
        oechem.OEThrow.Fatal(f"Unable to create the output SDF file: {sdf_file}")
    # Convert PDB to SDF using a generator for reading molecules
    for mol in ifs.GetOEGraphMols():
        oechem.OEWriteMolecule(ofs, mol)
    # Close the streams
    ifs.close()
    ofs.close()
    print(f"Successfully converted {ligand_file} to {sdf_file}.")

def main():
    parser = argparse.ArgumentParser(description="Process PDB files and generate ligand images.")
    parser.add_argument("--protein_ligand_file", default="protein_ligand.pdb", help="Path to the input PDB protein_ligand file.")
    parser.add_argument("--ligand_resname", default="BEN", help="Residue name of the ligand.")
    parser.add_argument("--ligand_file", default="ligand.pdb", help="Ligand PDB file.")
    parser.add_argument("--protein_file", default="protein.pdb", help="Protein PDB file.")
    parser.add_argument("--sdf_file", default="ligand.sdf", help="Ligand SDF file.")
    parser.add_argument("--dpi", type=int, nargs=2, default=(1000, 1000), help="DPI for the image.")
    args = parser.parse_args()
    # Split protein and ligand
    protein_ligand_pdbs(protein_ligand_file=args.protein_ligand_file, ligand_resname=args.ligand_resname, 
                        ligand_file=args.ligand_file, protein_file=args.protein_file)
    # Convert PDB to SDF
    pdb_to_sdf(ligand_file=args.ligand_file, sdf_file=args.sdf_file)

if __name__ == "__main__":
    main()
