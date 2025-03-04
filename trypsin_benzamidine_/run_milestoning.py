import subprocess
import argparse

def run_milestoning_sim(seekr2_path, input_XML_file):
    """
    Run the milestoning simulations for receptor-ligand complexes.
    
    Parameters:
        seekr2_path (str): Path to the SEEKR2 directory.
        input_XML_file (str): Path to the model.xml file.
    """
    seekr2_run_script = f"{seekr2_path}/seekr2/run.py"
    command = ["python", seekr2_run_script, "any", input_XML_file]
    try:
        # Execute the command
        result = subprocess.run(command, check=True, capture_output=False, text=True)
        print("Command output:", result.stdout)
        print("Milestoning simulation completed successfully.")
    except subprocess.CalledProcessError as e:
        print("Error occurred while executing the command:", e.stderr)

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description="Run a Milestoning simulation using SEEKR2 framework.")
    parser.add_argument("--seekr2_path", type=str, default="/home/USERNAME/seekr2", help="Path to the SEEKR2 directory. Default: /home/USERNAME/seekr2")
    parser.add_argument("--input_XML_file", type=str, default="SEEKR_SIMULATION/model.xml", help="Path to the model.xml file. Default: SEEKR_SIMULATION/model.xml")
    args = parser.parse_args()
    run_milestoning_sim(args.seekr2_path, args.input_XML_file)

if __name__ == "__main__":
    main()
