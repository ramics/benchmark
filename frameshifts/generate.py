import random
import subprocess
import sys
from Bio import Seq, SeqIO, SeqRecord, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna


def run_art(sequence):

    in_file = "temp_art.fasta"
    out_file = "temp_aln"
    coverage = "200"

    ungapped_record = SeqRecord.SeqRecord(sequence.seq.ungap("-"),
        id=sequence.id,
        name=sequence.name,
        description="")

    SeqIO.write([ungapped_record], open(in_file, "w"), "fasta")

    commands = ['art_454', '-a', '-A', in_file, out_file, coverage]
    subprocess.check_call(subprocess.list2cmdline(commands), shell=True)

    aln_filename = out_file + '.aln'
    fastq_filename = out_file + '.fq'
    
    return (aln_filename, fastq_filename)

def gap_reads(intermediate_alns):

    global_pos = 0

    gapped_alignment = []
    alignment_positions = []
    names = []

    for f in intermediate_alns:
        for read in f["reads"]:
            gapped_alignment.append("")
            alignment_positions.append(0)
            names.append(read.name)

    while True:
        in_insert = False
        finished = True

        idx = -1
        for f in intermediate_alns:
            for i in range(len(f["reads"])):
                idx += 1

                ref = f["alignments"][i][0]

                if alignment_positions[idx] >= len(ref):
                    continue 

                finished = False

                if (
                    global_pos < f["start_positions"][i]
                    or f["sequence"][global_pos] == "-"
                ):
                    continue

                if ref[alignment_positions[idx]] == "-":
                    in_insert = True

        if finished:
            break

        idx = -1
        for f in intermediate_alns:
            for i in range(len(f["reads"])):
                idx += 1

                ref = f["alignments"][i][0]
                query = f["alignments"][i][1]

                if (
                    global_pos < f["start_positions"][i]
                    or f["sequence"][global_pos] == "-"
                    or alignment_positions[idx] >= len(ref)
                ):
                    gapped_alignment[idx] += "-"
                    continue

                if in_insert:
                    if ref[alignment_positions[idx]] == "-":
                        gapped_alignment[idx] += query[alignment_positions[idx]]
                        alignment_positions[idx] += 1
                    else:
                        gapped_alignment[idx] += "-"  
                else:
                    gapped_alignment[idx] += query[alignment_positions[idx]]
                    alignment_positions[idx] += 1

        if not in_insert:
            global_pos += 1

    final_gapped_alignment = []
    for i, string_alignment in enumerate(gapped_alignment):
        name = names[i]
        final_gapped_alignment.append(
            SeqRecord.SeqRecord(Seq.Seq(string_alignment), id=name, name=name,
                description="")
        )
    return final_gapped_alignment

def aln_parse(aln_filename, fastq_filename):
    """ Convert an ART alignment output to a BioPython one """

    reads = []
    aligns = []
    start_positions = []
    fastqs = SeqIO.parse(open(fastq_filename), "fastq")
    
    with open(aln_filename, 'r') as fp:
        # Skip the first two lines - no info we want there
        line = fp.readline()
        line = fp.readline()
        line = fp.readline()
        # Deal with headers and get sequence length
        header_array = str.split(line, '\t')
        if str.rstrip(header_array[0]) != '@SQ':
            print "ART file header corrupted: ", header_array[0]
            sys.exit(2)
        header_array = str.split(line, '\t')
        sequence_length = int(header_array[2])
        while line:

            # Find the next record
            while not line.startswith(">"):
                line = str.rstrip(fp.readline())

            name_array = str.split(line, '\t')

            strand = str.rstrip(name_array[3])
            
            if strand == '-':
                continue
            
            name = name_array[1]

            start_pos = int(name_array[2])

            read = fastqs.next()
      
            clean_seq = Seq.Seq(str.rstrip(fp.readline()), generic_dna)
            dirty_seq = Seq.Seq(str.rstrip(fp.readline()), generic_dna)
             
            align = [   SeqRecord.SeqRecord(clean_seq, id=read.name, name=read.name,
                            description=""),
                        SeqRecord.SeqRecord(dirty_seq, id=read.name, name=read.name,
                            description="")]

            reads.append(read)
            aligns.append(align)
            start_positions.append(start_pos)

            line = fp.readline()

    return reads, aligns, start_positions

def generate_read_alignment(alignment_file, in_file, out_file):
    alignment = AlignIO.read(open(alignment_file), "fasta")

    final_aln = []
    final_reads = []

    intermediate_alns = []
    for seq in alignment:
        aln_file, fq_file = run_art(seq)
        reads, aligns, start_positions = aln_parse(aln_file, fq_file)
        intermediate_alns.append({"sequence": seq, 
            "reads": reads, 
            "alignments": aligns, 
            "start_positions": start_positions
        })
        final_reads.extend(reads)

    gapped_alignment = gap_reads(intermediate_alns)
    
    random.shuffle(final_reads)

    for seq in final_aln:
        print seq.id, len(seq)
    final_msa = MultipleSeqAlignment(gapped_alignment)
    SeqIO.write(final_reads, open(in_file, "w"), "fasta")
    AlignIO.write(final_msa, open(out_file, "w"), "fasta")

def main(alignment_file):

    in_file = alignment_file + ".in.fasta"
    out_file = alignment_file + ".out.fasta"
    generate_read_alignment(alignment_file, in_file, out_file)

if __name__ == "__main__":
  main(sys.argv[1])
