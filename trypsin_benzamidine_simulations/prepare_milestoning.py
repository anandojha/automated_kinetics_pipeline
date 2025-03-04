from parmed import openmm as parmed_openmm
from simtk.openmm import XmlSerializer
from simtk.openmm.app import PDBFile
import subprocess
import argparse
import shutil
import glob
import ast
import os

def get_pqr(xml_file, pdb_file, pqr_file_path):
    """
    Load system from an XML file, coordinates from a PDB file, filter out unwanted residues,
    and save the filtered structure as a PQR file.
    
    Args:
        xml_file (str): Path to the XML file containing system data.
        pdb_file (str): Path to the PDB file containing coordinates.
        pqr_file_path (str): Path where the output PQR file will be saved.
    """
    residues = ["HOH", "WAT", "NA", "CL"]
    # Load the system from the XML file
    with open(xml_file, 'r') as file:
        system = XmlSerializer.deserialize(file.read())
    # Load the coordinates from the PDB file
    pdb = PDBFile(pdb_file)
    positions = pdb.getPositions(asNumpy=True)
    # Load the structure using ParmEd
    structure = parmed_openmm.load_topology(pdb.topology, system)
    structure.positions = positions
    # Function to filter out unwanted residues from the structure
    def filter_structure(structure, residues):
        filtered_indices = []
        for residue in structure.residues:
            if residue.name not in residues:
                filtered_indices.extend(atom.idx for atom in residue.atoms)
        # Create a new structure with filtered atoms
        filtered_structure = structure[filtered_indices]
        return filtered_structure
    # Filter out unwanted residues
    filtered_structure = filter_structure(structure, residues)
    # Save the filtered structure as a PQR file
    if os.path.exists(pqr_file_path):
        os.remove(pqr_file_path)
    filtered_structure.save(pqr_file_path, format='pqr')
    print(f"PQR file saved as '{pqr_file_path}'")
    
def split_pqr_file(pqr_file_path, hetatm_pqr_path, atom_pqr_path):
    """
    Splits a PQR file into separate HETATM and ATOM entries.

    Parameters:
    pqr_file_path (str): Path to the input PQR file.
    hetatm_pqr_path (str): Path to save the HETATM entries.
    atom_pqr_path (str): Path to save the ATOM entries.
    """
    # Read the contents of the PQR file
    with open(pqr_file_path, 'r') as pqr_file:
        pqr_contents = pqr_file.readlines()
    # Split the contents into HETATM and ATOM entries
    hetatm_lines = [line for line in pqr_contents if line.startswith("HETATM")]
    atom_lines = [line for line in pqr_contents if line.startswith("ATOM")]
    # Save HETATM entries to a separate PQR file
    with open(hetatm_pqr_path, 'w') as hetatm_file:
        hetatm_file.writelines(hetatm_lines)
    # Save ATOM entries to another PQR file
    with open(atom_pqr_path, 'w') as atom_file:
        atom_file.writelines(atom_lines)
    print(f"HETATM entries saved to: {hetatm_pqr_path}")
    print(f"ATOM entries saved to: {atom_pqr_path}")
    
def modify_ligand_pqr(hetatm_pqr_path):
    """
    Reads a PQR file, modifies its content, and overwrites the original file.
    
    Args:
    hetatm_pqr_path (str): Path to the PQR file.
    """
    try:
        # Read the original content
        with open(hetatm_pqr_path, 'r') as file:
            pqr_lines = file.readlines()
        modified_lines = []
        for index, line in enumerate(pqr_lines):
            parts = line.split()
            # Update atom number and residue number
            parts[1] = str(index + 1)  # Atom number
            parts[5] = str(index + 1)  # Residue number
            # Remove chain ID and ensure correct formatting
            parts[4] = ''  # Chain ID is blank
            # Reformat the line with correct spacing
            modified_line = f"{parts[0]:<6}{parts[1]:>5}  {parts[2]:<4} {parts[3]:<4}{parts[4]:>1}{parts[5]:>4}  {parts[6]:>9}  {parts[7]:>6}  {parts[8]:>6} {parts[9]:>8} {parts[10]:>8}\n"
            modified_lines.append(modified_line)
        # Overwrite the original file with modified content
        with open(hetatm_pqr_path, 'w') as file:
            file.writelines(modified_lines)
    except Exception as e:
        print(f"An error occurred: {e}")

def create_pdb_protein_ligand(input_pdb, output_pdb):
    """
    Filters out specified entries from a PDB file and saves the filtered content to a new file.

    Args:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to the output PDB file.
    """
    entries = ["HOH", "WAT", "CRYST", "CONECT", "NA", "CL", "REMARK", "TER", "END"]
    # Read the input PDB file
    with open(input_pdb, 'r') as file:
        lines = file.readlines()
    # Filter lines to exclude certain entries
    filtered_lines = []
    for line in lines:
        if line.startswith(tuple(entries)):
            continue
        if line.startswith("ATOM") or line.startswith("HETATM"):
            residue_name = line[17:20].strip()
            if residue_name in entries:
                continue
        filtered_lines.append(line)
    # Write the filtered lines to the output PDB file
    with open(output_pdb, 'w') as file:
        file.writelines(filtered_lines)
    print(f"PDB file saved as '{output_pdb}'")

def extract_atoms_from_pdb(pdb_file):
    """
    Reads a PDB file and extracts atom types for protein and ligand atoms, counts each entry, and prints unique atom types.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        tuple: Two lists containing atom types for protein and ligand, respectively.
    """
    protein_atoms = []
    ligand_atoms = []
    # Read the PDB file
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                atom_type = line[12:16].strip()
                protein_atoms.append(atom_type[0] if atom_type[0].isalpha() else atom_type)
            elif line.startswith("HETATM"):
                atom_type = line[12:16].strip()
                ligand_atoms.append(atom_type[0] if atom_type[0].isalpha() else atom_type)
    # Create sets for unique atom types
    unique_protein_atoms = set(protein_atoms)
    unique_ligand_atoms = set(ligand_atoms)
    #print("Protein atom types:", protein_atoms)
    #print("Ligand atom types:", ligand_atoms)
    print("Total protein atoms:", len(protein_atoms))
    print("Total ligand atoms:", len(ligand_atoms))
    print("Unique protein atom types:", unique_protein_atoms)
    print("Unique ligand atom types:", unique_ligand_atoms)
    return protein_atoms, ligand_atoms

def map_atom_radii(atom_list, atom_radii):
    """
    Maps atom types to their corresponding radii.

    Args:
        atom_list (list): List of atom types.
        atom_radii (dict): Dictionary of atom types to their corresponding radii.

    Returns:
        list: List of radii corresponding to the atom types in the input list.
    """
    radii_list = []
    for atom in atom_list:
        if atom in atom_radii:
            radii_list.append(atom_radii[atom])
        else:
            radii_list.append(None)  
    return radii_list

def update_pqr_radii_inplace(pqr_path, new_radii_list):
    """
    Updates the radii in a PQR file with new radii values, overwriting the original file while preserving original formatting.

    Args:
        pqr_path (str): Path to the PQR file to be updated.
        new_radii_list (list): List of new radii values.
    """
    updated_lines = []
    with open(pqr_path, 'r') as file:
        lines = file.readlines()
    for line, new_radii in zip(lines, new_radii_list):
        # Locate the last whitespace before the radius and slice the string there
        last_space_index = line.rfind(' ') + 1
        # Preserve the formatting by reconstructing the line with the new radius
        updated_line = line[:last_space_index] + f"{new_radii:.4f}\n"
        updated_lines.append(updated_line)
    # Overwrite the original file with the updated radii
    with open(pqr_path, 'w') as file:
        file.writelines(updated_lines)
    print(f"Updated PQR file saved over the original file at {pqr_path}")

def move_intermediates(dest_folder, patterns):
    """
    Move files that match specific patterns to a destination folder. 
    If the folder doesn't exist, it will be created. The folder will not be deleted if it already exists.
    
    Parameters:
    - dest_folder: The name of the folder to be created.
    - patterns: A list of filename patterns to search for and move.
    """
    # Check if the destination folder exists; if not, create it
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)
        print(f"Created new folder: {dest_folder}")
    else:
        print(f"Folder already exists: {dest_folder}")
    # Find files that match the patterns and move them to the folder
    for pattern in patterns:
        for file_name in glob.glob(pattern):
            if os.path.exists(file_name):
                shutil.move(file_name, os.path.join(dest_folder, os.path.basename(file_name)))
                print(f"Moved file {file_name} to {dest_folder}")
            else:
                print(f"No files matching pattern {pattern} were found.")

def run_seekr2_prepare(input_xml, seekr2_path):
    """
    Runs the SEEKR2 prepare command with the specified input XML file and SEEKR2 installation directory.

    Args:
        input_xml (str): Path to the input XML file.
        seekr2_path (str): Path to the directory where SEEKR2 is installed.
    
    Returns:
        None
    """
    # Construct the full path to the prepare.py script
    prepare_script = os.path.join(seekr2_path, "seekr2", "prepare.py")
    # Construct the command to be executed
    command = ["python", prepare_script, input_xml]
    try:
        # Execute the command using subprocess
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Print the output of the command
        print(result.stdout.decode())
    except subprocess.CalledProcessError as e:
        print(f"Error running SEEKR2 prepare: {e.stderr.decode()}")
        raise

def main():
    parser = argparse.ArgumentParser(description="Process PQR and PDB files.")
    parser.add_argument('--xml_file', type=str, default='receptor_ligand.xml', help="XML file for system.")
    parser.add_argument('--pdb_file', type=str, default='receptor_ligand.pdb', help="PDB file for structure.")
    parser.add_argument('--pqr_file', type=str, default='complex.pqr', help="Output PQR file.")
    parser.add_argument('--hetatm_pqr', type=str, default='ligand.pqr', help="PQR file for HETATM records.")
    parser.add_argument('--atom_pqr', type=str, default='receptor.pqr', help="PQR file for ATOM records.")
    parser.add_argument('--output_pdb', type=str, default='receptor_ligand_BD.pdb', help="Output PDB file after filtering.")
    parser.add_argument('--move_intermediates', type=str, default='yes', choices=['yes', 'no'], help="Specify whether to move intermediate files (yes or no). Default is yes.")
    parser.add_argument('--intermediates_folder', type=str, default='intermediates', help="Folder where intermediate files will be moved.")
    parser.add_argument('--intermediates_patterns', nargs='+', default=["complex.pqr", "ligand.pqr", "receptor.pqr", "receptor_ligand_BD.pdb", "input.xml", "receptor_ligand_min.pdb", "*.txt"], help="Patterns for files to be moved as intermediates.")
    parser.add_argument('--atom_radii', type=str, default="{'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80, 'F': 1.47, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98, 'Fe': 1.80, 'Cu': 1.40, 'Zn': 1.39, 'Ca': 1.00, 'Mg': 1.18}", help="Dictionary of atom radii, specified as a string.")
    parser.add_argument('--input_xml', type=str, default="input.xml", help="Path to the input XML file.")
    parser.add_argument('--seekr2_path', type=str, default="/home/USERNAME/seekr2", help="Path to the directory where SEEKR2 is installed.")
    args = parser.parse_args()
    atom_radii = ast.literal_eval(args.atom_radii)
    get_pqr(xml_file=args.xml_file, pdb_file=args.pdb_file, pqr_file_path=args.pqr_file)
    split_pqr_file(pqr_file_path=args.pqr_file, hetatm_pqr_path=args.hetatm_pqr, atom_pqr_path=args.atom_pqr)
    modify_ligand_pqr(hetatm_pqr_path=args.hetatm_pqr)
    create_pdb_protein_ligand(input_pdb=args.pdb_file, output_pdb=args.output_pdb)
    protein_atoms, ligand_atoms = extract_atoms_from_pdb(pdb_file=args.output_pdb)
    protein_atom_radii = map_atom_radii(atom_list=protein_atoms, atom_radii=atom_radii)
    ligand_atom_radii = map_atom_radii(atom_list=ligand_atoms, atom_radii=atom_radii)
    update_pqr_radii_inplace(pqr_path=args.atom_pqr, new_radii_list=protein_atom_radii)
    update_pqr_radii_inplace(pqr_path=args.hetatm_pqr, new_radii_list=ligand_atom_radii)
    run_seekr2_prepare(input_xml=args.input_xml, seekr2_path=args.seekr2_path)
    if args.move_intermediates == 'yes':
        move_intermediates(dest_folder=args.intermediates_folder, patterns=args.intermediates_patterns)

if __name__ == "__main__":
    main()
