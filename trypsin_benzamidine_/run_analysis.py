import argparse
import subprocess

def run_milestoning_analysis(seekr2_path, xml_file, output_file):
    """ Run the milestoning analysis for receptor-ligand complexes. """
    analyze_script = f"{seekr2_path}/seekr2/analyze.py"
    command = ["python", analyze_script, xml_file]
    try:
        with open(output_file, "w") as out:
            result = subprocess.run(command, check=True, stdout=out, stderr=subprocess.PIPE, text=True)
        print("Milestoning analysis completed successfully.")
    except subprocess.CalledProcessError as e:
        print("Error occurred while executing the command:", e.stderr)

def main():
    parser = argparse.ArgumentParser(description="Run milestoning analysis with SEEKR2 framework.")
    parser.add_argument("--seekr2_path", type=str, default="/mnt/home/USERNAME/seekr2", help="Path to the SEEKR2 directory. Default: /mnt/home/USERNAME/seekr2")
    parser.add_argument("--xml_file", type=str, default="SEEKR_SIMULATION/model.xml", help="Path to the model.xml file. Default: SEEKR_SIMULATION/model.xml")
    parser.add_argument("--output_file", type=str, default="SEEKR_SIMULATION/analyze.out", help="Path to the output file. Default: SEEKR_SIMULATION/analyze.out")
    args = parser.parse_args()
    run_milestoning_analysis(args.seekr2_path, args.xml_file, args.output_file)

if __name__ == "__main__":
    main()
