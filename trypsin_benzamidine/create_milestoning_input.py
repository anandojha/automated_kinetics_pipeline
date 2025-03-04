import mdtraj as md
import numpy as np
import argparse
import os

def process_pdb(input_pdb, processed_pdb):
    # Keywords to filter out
    keywords = ["REMARK", "CRYST", "TER", "HOH", "WAT", "CONECT", "END", "Na", "Cl"]
    resnames = ["XX0", "XX1"] 
    with open(input_pdb, 'r') as infile, open(processed_pdb, 'w') as outfile:
        for line in infile:
            # Check if the line contains any of the keywords
            if any(keyword in line for keyword in keywords):
                # Skip the line if it does not belong to the resnames
                if not any(resname in line for resname in resnames):
                    continue
            # Write the line if it passes the filter
            outfile.write(line)

def get_full_path(filename):
    # Get the current working directory
    current_directory = os.getcwd()
    # Join the current directory with the filename
    full_path = os.path.join(current_directory, filename)
    return full_path

def parse_pdb(processed_pdb):
    """Parse the PDB file to separate protein and ligand atoms."""
    with open(processed_pdb, 'r') as file:
        lines = file.readlines()
    protein_atoms = []
    ligand_atoms = []
    for line in lines:
        if line.startswith("HETATM") and " CA " in line:  
            atom_index = int(line[6:11].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            protein_atoms.append((atom_index, np.array([x, y, z])))
        elif line.startswith("HETATM") and "XX1" in line:
            atom_index = int(line[6:11].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            ligand_atoms.append((atom_index, np.array([x, y, z])))
    return protein_atoms, ligand_atoms

def extract_alpha_carbon_and_ligand_indices(protein_atoms, ligand_atoms, threshold):
    """Extract indices of alpha carbons close to the ligand and indices of all ligand atoms."""
    close_c_alpha_atom_indices = []
    ligand_atom_indices = [atom_index for atom_index, _ in ligand_atoms]
    for atom_index, ca_pos in protein_atoms:
        for _, lig_pos in ligand_atoms:
            distance = np.linalg.norm(ca_pos - lig_pos)
            if distance <= threshold:
                close_c_alpha_atom_indices.append(atom_index)
                break  
    return close_c_alpha_atom_indices, ligand_atom_indices

def adjust_alpha_carbon_indices(alpha_carbon_indices):
    """Adjust alpha carbon indices by subtracting one and calculate the length."""
    adjusted_indices = [index - 1 for index in alpha_carbon_indices]
    length_of_indices = len(adjusted_indices)
    return adjusted_indices, length_of_indices

def adjust_ligand_indices(ligand_indices):
    """Adjust ligand indices by subtracting one and calculate the length."""
    adjusted_indices = [index - 2 for index in ligand_indices]
    length_of_indices = len(adjusted_indices)
    return adjusted_indices, length_of_indices

def create_sequential_ligand_list(length):
    """Create a sequential list of integers from 0 to length-1."""
    return list(range(length))

def com_com_distance(processed_pdb, receptor_alpha_indices, ligand_atom_indices, com_com_file):
    """
    Calculate the COM-COM distance between receptor alpha carbons and ligand atoms for a single PDB frame.
    
    Parameters:
        processed_pdb (str): Path to the PDB file.
        receptor_alpha_indices (list): List of alpha carbon indices (0-based).
        ligand_atom_indices (list): List of ligand atom indices (0-based).
        
    Returns:
        float: COM-COM distance between receptor and ligand atoms in nanometers.
    """
    # Load PDB 
    traj = md.load(processed_pdb)
    topology = traj.topology
    # Lists to store atom details for verification
    receptor_atom_details = []
    ligand_atom_details = []
    # Collect XYZ coordinates and atom details for receptor alpha carbons
    receptor_positions = traj.xyz[0, receptor_alpha_indices, :]  # Shape: (n_receptor_atoms, 3)
    for idx, position in zip(receptor_alpha_indices, receptor_positions):
        atom = topology.atom(idx)
        receptor_atom_details.append([atom.index, atom.name, position])
    # Collect XYZ coordinates and atom details for ligand atoms
    ligand_positions = traj.xyz[0, ligand_atom_indices, :]  # Shape: (n_ligand_atoms, 3)
    for idx, position in zip(ligand_atom_indices, ligand_positions):
        atom = topology.atom(idx)
        ligand_atom_details.append([atom.index, atom.name, position])
    # Calculate COM of receptor alpha carbons
    receptor_com = receptor_positions.mean(axis=0)  # Shape: (3,)
    # Calculate COM of ligand atoms
    ligand_com = ligand_positions.mean(axis=0)  # Shape: (3,)
    # Calculate COM-COM distance in nanometers (nm)
    com_com_distance = np.linalg.norm(receptor_com - ligand_com)
    # Save the details and distance to a file
    with open(com_com_file, "w") as f:
        f.write("Receptor atom details (index, name, coordinates):\n")
        for detail in receptor_atom_details:
            f.write(f"{detail}\n")
        f.write("\nLigand atom details (index, name, coordinates):\n")
        for detail in ligand_atom_details:
            f.write(f"{detail}\n")
        f.write(f"\nCOM-COM distance in nanometers: {com_com_distance}\n")
    return com_com_distance

def calculate_com_com_distance(processed_pdb, receptor_alpha_indices, ligand_atom_indices):
    """
    Calculate the COM-COM distance between receptor alpha carbons and ligand atoms for a single PDB frame.
    
    Parameters:
        processed_pdb (str): Path to the PDB file.
        receptor_alpha_indices (list): List of alpha carbon indices (0-based).
        ligand_atom_indices (list): List of ligand atom indices (0-based).
        
    Returns:
        float: COM-COM distance between receptor and ligand atoms in nanometers.
    """
    # Load PDB 
    traj = md.load(processed_pdb)
    topology = traj.topology
    # Collect XYZ coordinates for receptor alpha carbons
    receptor_positions = traj.xyz[0, receptor_alpha_indices, :]  # Shape: (n_receptor_atoms, 3)
    # Collect XYZ coordinates for ligand atoms
    ligand_positions = traj.xyz[0, ligand_atom_indices, :]  # Shape: (n_ligand_atoms, 3)
    # Calculate COM of receptor alpha carbons
    receptor_com = receptor_positions.mean(axis=0)  # Shape: (3,)
    # Calculate COM of ligand atoms
    ligand_com = ligand_positions.mean(axis=0)  # Shape: (3,)
    # Calculate and return COM-COM distance in nanometers (nm)
    com_com_distance = np.linalg.norm(receptor_com - ligand_com)
    return com_com_distance

def write_input_xml(
    filename,
    calculation_type,
    md_output_interval,
    md_steps_per_anchor,
    temperature,
    pressure,
    ensemble,
    root_directory,
    md_program,
    constraints,
    rigidWater,
    hydrogenMass,
    integrator,
    timestep,
    nonbonded_cutoff,
    receptor_indices,
    ligand_indices_openMM,
    radii,
    system_filename,
    receptor_pqr_filename,
    ligand_pqr_filename,
    ligand_indices_BD,
    num_b_surface_trajectories,
    n_threads,
    comment_browndye_settings):  
    # Convert the lists to strings for XML formatting
    receptor_indices_str = ", ".join(map(str, receptor_indices))
    ligand_indices_openMM_str = ", ".join(map(str, ligand_indices_openMM))
    ligand_indices_BD_str = ", ".join(map(str, ligand_indices_BD))

    # Initialize the input_anchors content
    input_anchors_content = ""
    for i, radius in enumerate(radii):
        # Set bound_state and bulk_anchor based on the position in the list
        bound_state = "False"
        bulk_anchor = "False"
        
        if i == 0:
            bound_state = "True"
        elif i == len(radii) - 1:
            bulk_anchor = "True"
            system_filename = ""  # Assuming the last one does not have a system_filename
        # Append the input_anchor content for each radius
        input_anchors_content += f"""
                <input_anchor class="Spherical_cv_anchor">
                    <radius>{radius}</radius>
                    <lower_milestone_radius/>
                    <upper_milestone_radius/>
                    <starting_forcefield_params class="Forcefield_params">
                        <system_filename>{system_filename}</system_filename>
                        <box_vectors/>
                        <pdb_coordinates_filename></pdb_coordinates_filename>
                    </starting_forcefield_params>
                    <bound_state>{bound_state}</bound_state>
                    <bulk_anchor>{bulk_anchor}</bulk_anchor>
                </input_anchor>"""
    # Handle the hydrogenMass tag based on the value
    hydrogenMass_tag = "<hydrogenMass/>" if hydrogenMass == 1 else f"<hydrogenMass>{hydrogenMass}</hydrogenMass>"
    # Construct ions content
    ion_details = [
        {"radius": 1.2, "charge": -1.0, "conc": 0.01},
        {"radius": 0.9, "charge": 1.0, "conc": 0.01}
    ]
    ions_content = ""
    for ion in ion_details:
        ions_content += f"""
            <ion class="Ion">
                <radius>{ion['radius']}</radius>
                <charge>{ion['charge']}</charge>
                <conc>{ion['conc']}</conc>
            </ion>"""
    # Define the browndye_settings_input content
    browndye_settings_input_content = f"""
    <browndye_settings_input class="Browndye_settings_input">
        <binary_directory></binary_directory>
        <receptor_pqr_filename>{receptor_pqr_filename}</receptor_pqr_filename>
        <ligand_pqr_filename>{ligand_pqr_filename}</ligand_pqr_filename>
        <apbs_grid_spacing>0.5</apbs_grid_spacing>
        <receptor_indices>[{receptor_indices_str}]</receptor_indices>
        <ligand_indices>[{ligand_indices_BD_str}]</ligand_indices>
        <ions>{ions_content}
        </ions>
        <num_b_surface_trajectories>{num_b_surface_trajectories}</num_b_surface_trajectories>
        <n_threads>{n_threads}</n_threads>
    </browndye_settings_input>"""
    # Conditionally comment out the browndye_settings_input section
    if comment_browndye_settings:
        browndye_settings_input_content = f"<!--{browndye_settings_input_content}-->"
    # Combine all parts into the final XML content
    xml_content = f"""<?xml version="1.0" ?>
<model_input class='Model_input'>
    <calculation_type>{calculation_type}</calculation_type>
    <calculation_settings class="MMVT_input_settings">
        <md_output_interval>{md_output_interval}</md_output_interval>
        <md_steps_per_anchor>{md_steps_per_anchor}</md_steps_per_anchor>
    </calculation_settings>
    <temperature>{temperature}</temperature>
    <pressure>{pressure}</pressure>
    <ensemble>{ensemble}</ensemble>
    <root_directory>{root_directory}</root_directory>
    <md_program>{md_program}</md_program> 
    <constraints>{constraints}</constraints>
    <rigidWater>{rigidWater}</rigidWater>
    {hydrogenMass_tag}
    <integrator_type>{integrator}</integrator_type>
    <timestep>{timestep}</timestep>
    <nonbonded_cutoff>{nonbonded_cutoff}</nonbonded_cutoff>
    <cv_inputs>
        <cv_input class="Spherical_cv_input">
            <group1>[{receptor_indices_str}]</group1>
            <group2>[{ligand_indices_openMM_str}]</group2>
            <input_anchors>{input_anchors_content}
            </input_anchors>
        </cv_input>
    </cv_inputs>
    {browndye_settings_input_content}
</model_input>"""
    # Write the final content to the input.xml file
    with open(filename, "w") as f:
        f.write(xml_content)

def get_arguments():
    parser = argparse.ArgumentParser(description="Generate input.xml for SEEKR.")
    # Add all necessary arguments for the script with default values
    parser.add_argument('--input_pdb', type=str, default="receptor_ligand_min.pdb", help="Input PDB file path.")  
    parser.add_argument('--processed_pdb', type=str, default="receptor_ligand_BD.pdb", help="Output processed PDB file path.")
    parser.add_argument('--threshold', type=float, default=6.0, help="Threshold for alpha carbon and ligand distance.")
    parser.add_argument('--com_com_file', type=str, default="com_com.txt", help="COM-COM distance between the alphs-C atoms of the bound state and the ligand atoms.")
    parser.add_argument('--xml_output', type=str, default="input.xml", help="Output XML file path.")
    parser.add_argument('--system_filename', type=str, default="receptor_ligand.xml", help="Serialized system file for ESPALOMA.")
    parser.add_argument('--root_directory', type=str, default="SEEKR_SIMULATION", help="Root directory for the simulation.")
    parser.add_argument('--receptor_pqr', type=str, default="receptor.pqr", help="Receptor PQR filename.")
    parser.add_argument('--ligand_pqr', type=str, default="ligand.pqr", help="Ligand PQR filename.")
    parser.add_argument('--radii', nargs='+', type=float, default=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,1.3], help="List of radii.")
    parser.add_argument('--num_b_surface_trajectories', type=int, default=100000, help="Number of B-surface trajectories.")
    parser.add_argument('--n_threads', type=int, default=12, help="Number of threads for BrownDye.")
    parser.add_argument('--calculation_type', type=str, default='mmvt', help="Calculation type (e.g., mmvt).")
    parser.add_argument('--md_output_interval', type=int, default=10000, help="MD output interval.")
    parser.add_argument('--md_steps_per_anchor', type=int, default=200000, help="MD steps per anchor.")
    parser.add_argument('--temperature', type=float, default=300.0, help="Temperature in Kelvin.")
    parser.add_argument('--pressure', type=float, default=1.0, help="Pressure in atmospheres.")
    parser.add_argument('--ensemble', type=str, default='nvt', help="Ensemble type (e.g., nvt, npt).")
    parser.add_argument('--md_program', type=str, default='openmm', help="MD program to use.")
    parser.add_argument('--constraints', type=str, default='HBonds', help="Constraints for MD.")
    parser.add_argument('--rigidWater', type=str, default='True', help="Use rigid water or not.")
    parser.add_argument('--hydrogenMass', type=float, default=1.0, help="Hydrogen mass.")
    parser.add_argument('--integrator', type=str, default='langevin', help="Integrator type.")
    parser.add_argument('--timestep', type=float, default=0.002, help="Timestep in picoseconds.")
    parser.add_argument('--nonbonded_cutoff', type=float, default=1.0, help="Non-bonded cutoff distance.")
    parser.add_argument('--comment_browndye_settings', type=str, choices=['yes', 'no'], default='yes', help="Comment out browndye settings (yes or no). Default is yes.")
    return parser.parse_args()

def main():
    args = get_arguments()
    process_pdb(input_pdb=args.input_pdb, processed_pdb=args.processed_pdb)
    protein_atoms, ligand_atoms = parse_pdb(processed_pdb=args.processed_pdb)
    receptor_alpha_indices, ligand_indices = extract_alpha_carbon_and_ligand_indices(protein_atoms=protein_atoms, ligand_atoms=ligand_atoms, threshold=args.threshold)
    print("Receptor Alpha Carbon Indices for the given PDB file:", receptor_alpha_indices)
    print("Ligand Atom Indices for the given PDB file:", ligand_indices)
    adjusted_receptor_alpha_indices, receptor_alpha_carbon_length = adjust_alpha_carbon_indices(receptor_alpha_indices)
    print("Receptor Alpha Carbon Indices for model.xml SEEKR input file (0-based):", adjusted_receptor_alpha_indices)
    adjusted_ligand_indices, ligand_length = adjust_ligand_indices(ligand_indices)
    browndye_ligand_list = create_sequential_ligand_list(ligand_length)
    print("Ligand Atom Indices for model.xml SEEKR input file (0-based):", adjusted_ligand_indices)
    com_com_distance(processed_pdb=args.processed_pdb, receptor_alpha_indices=adjusted_receptor_alpha_indices, ligand_atom_indices=adjusted_ligand_indices, com_com_file=args.com_com_file)
    com_com_distance_nm = calculate_com_com_distance(processed_pdb=args.processed_pdb, receptor_alpha_indices=adjusted_receptor_alpha_indices, ligand_atom_indices=adjusted_ligand_indices)
    print("Sanity check: Ensure that the first radius value in the radii list is close to the initial COM-COM distance between the receptor alpha-carbon atoms and ligand atoms, which is:", com_com_distance_nm, "nm. For more details, please check the", args.com_com_file, "file.")
    comment_browndye = True if args.comment_browndye_settings == 'yes' else False
    write_input_xml(
        filename=get_full_path(args.xml_output),
        calculation_type=args.calculation_type,
        md_output_interval=args.md_output_interval,
        md_steps_per_anchor=args.md_steps_per_anchor,
        temperature=args.temperature,
        pressure=args.pressure,
        ensemble=args.ensemble,
        root_directory=get_full_path(args.root_directory), 
        md_program=args.md_program,
        constraints=args.constraints,
        rigidWater=args.rigidWater,
        hydrogenMass=args.hydrogenMass,
        integrator=args.integrator,
        timestep=args.timestep,
        nonbonded_cutoff=args.nonbonded_cutoff,
        receptor_indices=adjusted_receptor_alpha_indices,
        ligand_indices_openMM=adjusted_ligand_indices,
        radii=args.radii,
        system_filename=get_full_path(args.system_filename),  
        receptor_pqr_filename=get_full_path(args.receptor_pqr),  
        ligand_pqr_filename=get_full_path(args.ligand_pqr),  
        ligand_indices_BD=browndye_ligand_list,
        num_b_surface_trajectories=args.num_b_surface_trajectories,
        n_threads=args.n_threads,
        comment_browndye_settings=comment_browndye )

if __name__ == "__main__":
    main()
