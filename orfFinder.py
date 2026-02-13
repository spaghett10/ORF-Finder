from Bio import SeqIO
from Bio.Seq import Seq
import re

import os # to verify file extensions and not have BioPython returning warnings

def main():
    fasta_file = input("FASTA sequence file: ")  # take input of file of genetic data
    output_file = input("Desired output file name (no extension): ") # take output file name desired

    sequence, seq_length = parse(fasta_file) # will raise various errors if fasta file and/or content is invalid

    GContent = GC_calc(sequence) # calculate GC content

    output(fasta_file, output_file, GContent, seq_length, sequence) # generate output file


def parse(file):

    # Verify file type
    input_tup = os.path.splitext(file)

    if input_tup[1] != ".fna" and input_tup[1] != ".fasta":
        raise TypeError("Incorrect file type input")

    # Parse verified FASTA file, not accepting files with more than 10 FASTA sequences or no sequences
    else:
        records = list(SeqIO.parse(file, "fasta"))

        if len(records) > 10:
            raise ValueError("Over 10 separate FASTA sequences detected")
        elif len(records) < 1:
            raise ValueError("No FASTA sequences detected")

        # Create Sequence object with total sequence from all FASTA sequences
        totalsequence = Seq("")
        for seq_record in SeqIO.parse(file, "fasta"):
            totalsequence += seq_record.seq

        return totalsequence, len(totalsequence)


# calcluation of GC content for total sequence;
# seq argument is a Seq object, not a string sequence
def GC_calc(seq):
    upper_seq = seq.upper() #ensure uniform casing
    GC_content = ((upper_seq.count("G") + upper_seq.count("C")) / len(upper_seq)) * 100
    return round(GC_content, 4)


# Finding all possible protein sequences within DNA sequence
# seq argument is a Seq object, not a string sequence
def orf_finder(seq):

    start_codon = r"ATG" # regex pattern for start codon
    stop_codon = r"TAG|TAA|TGA" # regex pattern for stop codon

    # Start/stop codon match iterables put into lists, case independent
    start_list = list(re.finditer(start_codon, str(seq), re.IGNORECASE))
    stop_list = list(re.finditer(stop_codon, str(seq), re.IGNORECASE))

    # If no start or stop codons are found, end function
    if not start_list or not stop_list:
        raise ValueError

    else:

        orfs = []

        # Outer loop for iterating through start codons
        for orf1 in start_list:

            # Inner loop for iterating through stop codons and calculating nt difference between start and stop codons
            for orf2 in stop_list:
                diff = orf2.end() - orf1.start()

                # Filtering for sequences where stop codon comes after start codon and nt difference is divisible by 3
                if diff >= 30 and diff % 3 == 0:
                    orf_seq = seq[orf1.start():orf2.end()]

                    # Add translated sequence and it's start/stop indices to list that is nested within master list
                    orf = (orf_seq, orf1.start()+1, orf2.end()+1)
                    orfs.append(orf)

                    # end inner loop after valid first sequence approached to move onto next start codon
                    break

                else:
                    continue

    # return nested list of orfs if any found
    if not orfs:
        raise ValueError
    else:
        return orfs


# output file code
def output(in_file, out_file, GC, seq_len, seq):

    # Create output file with stats: title line, sequence nt length, GC content
    with open(f"{out_file}.txt", "a") as file:
        file.write(f"Summary Report for {in_file} File\n\n")
        file.write(f"Sequence Length: {seq_len} nucleotides\n")
        file.write(f"Total GC Content: {GC}%\n\n")

        try:
            sense_orfs = orf_finder(seq)
        except ValueError:
            file.write(f"No Potential Protein Sequences Found in Forward Strand\n\n")
        else:
            # Print out forward strand protein sequence(s) and start/stop locations from nested orf list
            for sq, stt, stp in sense_orfs:
                p = sq.translate()
                file.write(f"Strand: +\nLength: {len(str(p))} amino acids\nSequence: {str(p)}\nStart: {stt}\nEnd: {stp}\n\n")

        try:
            antisense_orfs = orf_finder(seq.reverse_complement())
        except ValueError:
            file.write(f"No Potential Protein Sequences Found in Reverse Strand")
        else:
            # Print out reverse strand protein sequence(s) and start/stop locations
            for sq, stt, stp in antisense_orfs:
                p = sq.translate()
                stt = seq_len - (stt-1) # Orient start indice to leading strand nt
                stp = seq_len - (stp-1) # Orient stop indice to leading strand nt
                file.write(f"Strand: -\nLength: {len(str(p))} amino acids\nSequence: {str(p)}\nStart: {stt}\nEnd: {stp}\n\n")

if __name__ == "__main__":
    main()
