#!/usr/bin/python
 
import getopt
import os
import pickle
import subprocess
from subprocess import PIPE
import sys
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
from Bio.Alphabet import generic_dna


def run(in_folder, command_file, out_folder):

    reads = {}
    results = {}

    for filename in [os.path.splitext(f)[0] for f in os.listdir(in_folder) if f.endswith("msf")]:
        in_file = os.path.join(in_folder, filename + ".tfa")
        out_file = "temp_out_" + os.path.basename(command_file)
        out_file_msf = out_file + ".msf"
        ref_file = os.path.join(in_folder, filename + ".msf")
        
        print in_file
        print out_file

        commands = ["/bin/bash", command_file, 
                    in_file,
                    out_file]
        subprocess.call(subprocess.list2cmdline(commands), \
            shell=True)

        reads[filename] = {}

        dn = open(os.devnull, 'w')
        commands = ["seqret", "-sequence", out_file, "-outseq", out_file_msf, 
            "-osformat2", "msf"]

        subprocess.call(subprocess.list2cmdline(commands), \
            shell=True, stdout=dn, stderr=dn)

        commands = ["bali_score", ref_file, out_file_msf]

        p = subprocess.Popen(commands, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate()
        try:
            sp = float((output.split("\n")[5]).split('=')[1])
            tc = float((output.split("\n")[7]).split('=')[1])
        except IndexError:
            print subprocess.list2cmdline(commands)
            print output
            print err
            continue

        #os.unlink(out_file)

        result = {
            "scores": [sp, tc]
        }

        results[filename] = result

    out_path = os.path.join(out_folder, os.path.basename(command_file))
    out = open(out_path, "wb")
    pickle.dump((reads, results), out, -1)


def main(argv):

    in_folder = ""
    ref_folder = ""
    command_file = ""
    out_folder = ""

    try:
        opts, args = getopt.getopt(argv,"i:c:o:",[
            "in-folder=",
            "ref-folder=",
            "command-file=",
            "out-folder="
        ])
    except getopt.GetoptError:
        print '[ERROR] Syntax: run.py -i in-folder -c command-file -o out_folder'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-c", "--command-file"):
            command_file = arg
        elif opt in ("-i", "--in-folder"):
            in_folder = arg
        elif opt in ("-o", "--out-folder"):
            out_folder = arg

    run(in_folder, command_file, out_folder)



if __name__ == "__main__":
    main(sys.argv[1:])
