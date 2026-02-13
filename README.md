# INITIAL ANALYSIS OF A DNA SEQUENCE

# ORF-Finder
The goal of this program is to take a DNA sequence and output a new file listing the length of the sequence, GC% content, and all potential ORFs and their locations in the sequence, thus identifying possible protein sequences.

## Inputting a DNA sequence
This program reads the DNA sequence from a FASTA file format only, a common file type for sequencing files. The parse() function within the program translates the sequence within a FASTA file into a readable Sequence object through the imported Biopython library. The function returns both this Sequence object and it's length, representing the total nucleotide count which is later included as one of the sequence stats in the output file.

>**IMPORTANT NOTE:**\
_The program only allows a maximum of 10 separate FASTA sequences within the input file. A ValueError will be returned if the file contains more than 10 sequences. This was done so as to not overwhelm an average CPU with too long of a sequence._

## Calculating GC% Content
Guanine and cytosine are two of the four nucleotides present within a DNA sequence. Calculating the percentage proportion that guanine and cytosine make up of the entire sequence can tell a biologist a number of factors about a sequence of unknown origin. GC% content can confer insights into structural stability, expression patterns, organismal origin, and experimental limitations (particularly for PCR primer design). Most importantly, GC% content can be used to discern between prokaryotes and eukaryotes as well as between different bacterial species which can be helpful while studying microbiomes, for example.\

The GC_calc() function calculates GC% content simply by tabulating the number of cytosine and guanine nucleotides encountered while iterating through the entire sequence and then dividing this by the total length of the sequence. This can be represented by the formula below.
```
GC_content = ((seq.count("G") + seq.count("C")) / len(seq)) * 100
```

## Finding all possible open reading frames (ORFs)
Since DNA serves as the template that ultimately writes proteins, it is essential in biology to be able to identify all possible protein sequences. Not all DNA is ultimately transcribed and then translated into a protein, but because we know that specific codons (a triplet of nucleotides) code for the start and end of a protein sequence we can deduce where a potential protein may start and end.

>**IMPORTANT NOTE:**\
 _This program utilizes the standard codon translation table which may be the incorrect table for certain sequences depending on the organism from where the sequence originates. This is something that could be improved with the program in the future by allowing input from the user on the desired translation table._

The orf_finder() function works by finding all possible start ("ATG") and stop ("TAG|TAA|TGA") codons using regular expressions. From there, it uses a nested loop to first iterate through each start codon and subsequently iterate through all the stop codons. If the sequence between a start codon and a stop codon is greater than 3 nucleotides and is divisible by 3, this sequence and it's location within the greater sequence are saved to lists. Once a sequence is found with a specific start codon, the inner loop is broken so the outer loop can continue onto the next potential start codon since a sequence cannot have multiple stop codons in it. The function concludes once all start codons have been iterated through and returns three lists: a list of the potential protein sequences, a list of all these sequences' start locations, and a list of all these sequences' stop locations. If no protein sequences are found, the function will return a ValueError stating "No protein sequences to be found" which will be incorporated into the output file upon generation.

## Output of new file
Upon running the program, the user can input the desired name of the output file. Typically, this is a .txt file name. The program will output this new file once the program has concluded. If the program is run again and the same output file name is used, the program will append the results to that file rather than replace the contents of the initial file. It should look like the below:
```
Summary Report for test.fasta File

Sequence Length: 96 nucleotides
Total GC Content: 47.9167%

Potential Protein Sequences:

Length: 14 amino acids
Sequence: MGMTPRLGLESLLE
Start: 25
End: 70

Length: 12 amino acids
Sequence: MTPRLGLESLLE
Start: 31
End: 76
```
If no protein sequences are found, GC% content will still be calculated and the output file will look like the below:
```
Summary Report for test.fasta File

Sequence Length: 96 nucleotides
Total GC Content: 47.9167%

No protein sequences to be found
```
