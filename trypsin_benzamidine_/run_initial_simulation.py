import subprocess
import argparse

def run_initial_simulations(seekrtools_path, input_XML_file, pdb_file, pulling_scheme, pulling_velocity):
    """
    Run the seekrtools HIDR command with specified parameters.
    
    Parameters:
        seekrtools_path (str): Path to the seekrtools directory.
        input_XML_file (str): Path to the model.xml file.
        pdb_file (str): Path to the PDB file.
        pulling_scheme (str): Pulling dynamics method to implement (e.g., SMD, metaD).
        pulling_velocity (float): Pulling velocity value to use.
    """
    command = [
        "python", f"{seekrtools_path}/seekrtools/hidr/hidr.py", 
        "any", input_XML_file, 
        "-M", pulling_scheme, 
        "-p", pdb_file, 
        "-v", str(pulling_velocity)
    ]
    try:
        # Execute the command
        result = subprocess.run(command, check=True, capture_output=False, text=True)
        print("Command output:", result.stdout)
        print("Command completed successfully.")
    except subprocess.CalledProcessError as e:
        print("Error occurred while executing the command:", e.stderr)

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Run the seekrtools HIDR command with specified parameters. HIDR uses enhanced sampling to populate all accessible anchors in SEEKR with initial structures.")
    parser.add_argument("--seekrtools_path", type=str, default="/home/USERNAME/seekrtools", help="Path to the seekrtools directory. Default: /home/USERNAME/seekrtools")
    parser.add_argument("--input_XML_file", type=str, default="SEEKR_SIMULATION/model.xml", help="Path to the model.xml file. Default: SEEKR_SIMULATION/model.xml")
    parser.add_argument("--pdb_file", type=str, default="receptor_ligand.pdb", help="Path to the PDB file containing the structure to be placed in the correct anchors. Default: receptor_ligand.pdb")
    parser.add_argument("--pulling_scheme", type=str, default="SMD", help="Pulling dynamics method to implement. Options include 'SMD', 'RAMD', 'MetaD', 'MetaDxyz'. Default: SMD")
    parser.add_argument("--pulling_velocity", type=float, default=0.5, help="Translation velocity in nm/ns to pull the system along an SMD trajectory for distance-based CVs. Default: 0.01 nm/ns")
    args = parser.parse_args()

    run_initial_simulations(args.seekrtools_path, args.input_XML_file, args.pdb_file, args.pulling_scheme, args.pulling_velocity)

if __name__ == "__main__":
    main()

