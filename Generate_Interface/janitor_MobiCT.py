import subprocess
import argparse

parser = argparse.ArgumentParser(description='This script is used to clean and organize MobiCT cfDNA outputs')
parser.add_argument('-i', '--input', type=str, help='path to the directory containing the output files of MobiCT')
parser.add_argument('-o', '--output', type=str, help='path to the desired output directory')
args = parser.parse_args()

# Creating the arborescence
command = "mkdir -p "+args.output
subprocess.run(command, shell=True)

command = "mkdir -p "+args.output+"/stats/"
subprocess.run(command, shell=True)

command = "mkdir -p "+args.output+"/variants/"
subprocess.run(command, shell=True)

command = "mkdir -p "+args.output+"/alignment/"
subprocess.run(command, shell=True)

# Moving files
command = "rsync -rvL "+args.input+"/*_consensusMerge* "+args.output+"/alignment/"
print(command)
subprocess.run(command, shell=True)

command = "rsync -rvL "+args.input+"/*_vep.vcf "+args.output+"/variants/"
subprocess.run(command, shell=True)

command = "rsync -rvL "+args.input+"/*metrics* "+args.output+"/stats/"
subprocess.run(command, shell=True)

command = "rsync -rvL "+args.input+"/*reporteAfter "+args.output+"/stats/"
subprocess.run(command, shell=True)

