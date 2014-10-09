#!/usr/bin/python
 
import getopt
import os
import pickle
import subprocess
from subprocess import PIPE
import sys
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio.Alphabet import generic_dna

def run(in_folder, ref_folder, command_file, out_folder):

    reads = {}
    results = {}

    for filename in os.listdir(in_folder):
        in_file = os.path.join(in_folder, filename)
        out_file = "temp_out_" + os.path.basename(command_file)
        ref_filename = filename.replace("test-only", "ref")
        ref_file = os.path.join(ref_folder, ref_filename)

        test = list(SeqIO.parse(open(in_file), "fasta"))
        ref = list(SeqIO.parse(open(ref_file), "fasta"))
        clean_ref = [SeqRecord(n.seq.ungap("-"),
            id = n.id, 
            name= n.name, 
            description = n.description) 
        for n in ref]
        
        test.extend(clean_ref)

        full_in_file = "temp_in_" + os.path.basename(command_file)
        SeqIO.write(test, open(full_in_file, "w"), "fasta")

        commands = ["/bin/bash", command_file, 
                    full_in_file,
                    out_file]
        subprocess.call(subprocess.list2cmdline(commands), \
            shell=True)

        names = [n.name for n in ref]
        out = AlignIO.read(open(out_file), "fasta")
        in_ref = [n for n in out if n.name in names]
        
        SeqIO.write(in_ref, open(out_file, "w"), "fasta")
        summary = AlignInfo.SummaryInfo(ref)
        

        length = len(ref)

        read_info = {"length" : length,
        }

        reads[filename] = read_info

        commands = ["qscore", "-test", out_file, "-ref", 
            ref_file, "-cline", "-modeler", "-seqdiffwarn", "-truncname", "-ignoretestcase", "-ignorerefcase"]

        p = subprocess.Popen(commands, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate()
        print output
        try:
            q = float(output.split(";")[2].split('=')[1])
            tc = float(output.split(";")[3].split('=')[1])
            cline = float(output.split(";")[4].split('=')[1])
            modeler = float(output.split(";")[5].split('=')[1])
        except IndexError:
            print subprocess.list2cmdline(commands)
            print output
            print err
            continue

        #os.unlink(out_file)

        result = {
            "scores": [q, tc, cline, modeler]
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
        opts, args = getopt.getopt(argv,"i:r:c:o:",[
            "in-folder=",
            "ref-folder=",
            "command-file=",
            "out-folder="
        ])
    except getopt.GetoptError:
        print '[ERROR] Syntax: run.py -i in-folder -r ref-folder -c command-file -o out_folder'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-c", "--command-file"):
            command_file = arg
        elif opt in ("-i", "--in-folder"):
            in_folder = arg
        elif opt in ("-r", "--ref-folder"):
            ref_folder = arg
        elif opt in ("-o", "--out-folder"):
            out_folder = arg

    run(in_folder, ref_folder, command_file, out_folder)



if __name__ == "__main__":
    main(sys.argv[1:])