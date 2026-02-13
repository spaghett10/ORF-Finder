import pytest
from project import parse, GC_calc, orf_finder, output
from Bio.Seq import Seq
from pathlib import Path

@pytest.fixture
def test_data_dir():
    return Path(__file__).parent / "test_data"

# Testing parse() function
class TestParse:
    def test_incorrect_file_type(self, test_data_dir):
        #test for non-fasta files
        with pytest.raises(TypeError, match = "Incorrect file type input"):
            parse(str(test_data_dir / "test_incorrect_filetype.txt"))

    def test_empty_fasta(self, test_data_dir):
        # test with empty fasta file
        with pytest.raises(ValueError):
            parse(str(test_data_dir / "test_no_fasta_sequences.fna"))

    def test_parse_file_not_found(self, test_data_dir):
        # test parsing non-existent file
        with pytest.raises(FileNotFoundError):
            parse(str(test_data_dir / "nonexistent.fasta"))

    def test_parse_invalid_content(self, test_data_dir):
        # test for files that don't contain required number of sequences
        with pytest.raises(ValueError, match = "Over 10 separate FASTA sequences detected"):
            parse(str(test_data_dir / "test_too_many_fasta_sequences.fna"))
        with pytest.raises(ValueError, match = "No FASTA sequences detected"):
            parse(str(test_data_dir / "test_no_fasta_sequences.fna"))

    def test_parse_simple(self, test_data_dir):
        # test with single line fasta
        sequence, length = parse(str(test_data_dir / "test_simple.fasta"))
        assert isinstance(sequence, Seq)
        assert length == len(sequence)

    def test_parse_multiline(self, test_data_dir):
        # test with multi-line fasta file
        sequence, length = parse(str(test_data_dir / "test_multiple.fasta"))
        assert isinstance(sequence, Seq)
        assert length == len(sequence)
        assert sequence == Seq("ATGGCATAGATGTCCTAA")

    def test_real_file(self, test_data_dir):
        # test with real fasta file
        sequence, length = parse(str(test_data_dir / "test_real_seq.fna"))
        assert isinstance(sequence, Seq)
        assert length == len(sequence)

# Testing GC_Calc function
@pytest.mark.parametrize("seq, GC_content", [
    (Seq("GCGCCCCCGGGGCGGCGGCCCGGGGGCGCGCGCCCGCGC"), 100.0), #100% GC Content
    (Seq("ATTTTATTATATATAATTATAATATATTATAAAAA"), 0.0), #0% GC Content
    (Seq("GCGCGGCATTATATTGGCCGCAATAT"), 50.0), #50% GC Content
    (Seq("gcgcgGCattatATTggccgcaaTAt"), 50.0), #Handles mixed casing
    (Seq("ATGGCAGGCATGACACCTCGCCTGGGCCTGGAATCCCTGCTGGAATAACTAGTGA"), 56.3636) #Mixed% GC content
])
def test_GCcontent(seq, GC_content):
    assert GC_calc(seq) == GC_content

# Testing orf_finder function
class TestOrfFinder:
    def test_one_orf(self):
        # finding one orf
        result = orf_finder(Seq("ATGGATTACCACGGGATCTAG"))
        assert isinstance(result, list)
        assert len(result) == 1
        assert len(result[0]) == 3

    def test_multiple_orfs(self):
        # finding multiple orfs
        result = orf_finder(Seq("ATGGCAGGCATGACACCTCGCCTGGGCCTGGAATCCCTGCTGGAATAACTAGTGA"))
        assert isinstance(result, list)
        assert len(result) == 2
        assert len(result[0]) == 3 and len(result[1]) == 3

    def test_multiple_stop_codons(self):
        # finds one orf given sequence with multiple stop codons (earlier stop codon is where sequence should end)
        result = orf_finder(Seq("ATGGCATAGTAA"))
        assert result[0][0] == Seq("ATGGCATAG")

    def test_multiple_start_codons(self):
        # finds two orfs given sequence with multiple start codons and one stop codon
        result = orf_finder(Seq("ATGGCATTTATGCGGCCGTAA"))
        assert len(result) == 2
        assert result[0][0] == Seq("ATGGCATTTATGCGGCCGTAA")
        assert result[1][0] == Seq("ATGCGGCCGTAA")

    def test_no_start_codon(self):
        # finding no orfs when only stop codons present
        with pytest.raises(ValueError):
            orf_finder(Seq("GGGTAACCATAG"))

    def test_no_stop_codon(self):
        # finding no orfs when only start codons present
        with pytest.raises(ValueError):
            orf_finder(Seq("GGGATGCCACTCTCCATGTATA"))

    def test_no_codons(self):
        # finding no start or stop codons
        with pytest.raises(ValueError):
            orf_finder(Seq("GGGGCCCCAAAA"))

    def test_not_divisible_by3(self):
        # finding no orfs when both start and stop codons present, but not divisible by 3 (unable to be translated to protein)
        with pytest.raises(ValueError):
            orf_finder(Seq("ATGTACCCTAG"))

# Testing output function
class TestOutput:
    def test_output_creates_file(self, tmp_path, test_data_dir):
        # test that output file is created
        output_file = tmp_path / "test_output"
        seq, seq_len = parse(str(test_data_dir / "test_simple.fasta"))
        gc_content = GC_calc(seq)
        output("test_simple.fasta", str(output_file), gc_content, seq_len, seq)

        assert output_file.with_suffix('.txt').exists()

    def test_output_no_orfs_message(self, tmp_path, test_data_dir):
         # test that output file contains no protein sequences found messaging
        output_file = tmp_path / "test_output"
        seq, seq_len = parse(str(test_data_dir / "test_no_orf.fasta"))
        gc_content = GC_calc(seq)
        output("test_no_orf.fasta", str(output_file), gc_content, seq_len, seq)

        content = f"No Potential Protein Sequences Found in Forward Strand\n\nNo Potential Protein Sequences Found in Reverse Strand"
        assert content in output_file.with_suffix('.txt').read_text()

    # test for various content in output file
    @pytest.mark.parametrize("content", [
    (f"Summary Report for test_simple.fasta File"), #title line
    (f"Total GC Content: 56.3636%\n\n"), #GC content line
    (f"Sequence Length: 55 nucleotides\n"), #sequence length line
    ])

    def test_output_file_content(self, content, tmp_path, test_data_dir):
        output_file = tmp_path / "test_output"
        seq, seq_len = parse(str(test_data_dir / "test_simple.fasta"))
        gc_content = GC_calc(seq)
        output("test_simple.fasta", str(output_file), gc_content, seq_len, seq)

        assert content in output_file.with_suffix('.txt').read_text()


# Testing entire program
class TestIntegration:
    def test_full_workflow_with_orfs(self, test_data_dir):
        # test real fasta file with orfs
        seq, seq_len = parse(str(test_data_dir / "test_real_seq.fna"))
        assert seq_len > 0

        gc_content = GC_calc(seq)
        assert 0 <= gc_content <= 100

        orfs = orf_finder(seq)
        assert len(orfs) == 13

    def test_reverse_complement_orfs(self, test_data_dir):
        # test that orfs identifes in lagging strand
        seq, _ = parse(str(test_data_dir / "test_real_seq.fna"))

        rev_orfs = orf_finder(seq.reverse_complement())
        assert len(rev_orfs) == 10

    def test_full_workflow_no_orfs(self, tmp_path, test_data_dir):
        # test valid fasta file with no orfs and correct messaging displayed
        output_file = tmp_path / "test_output"

        seq, seq_len = parse(str(test_data_dir / "test_no_orf.fasta"))
        assert seq_len > 0

        gc_content = GC_calc(seq)
        assert 0 <= gc_content <= 100

        with pytest.raises(ValueError):
            orf_finder(seq)

        with pytest.raises(ValueError):
            orf_finder(seq.reverse_complement())

        output("test_no_orf.fasta", str(output_file), gc_content, seq_len, seq)

        content1 = "No Potential Protein Sequences Found in Forward Strand"
        content2 = "No Potential Protein Sequences Found in Reverse Strand"

        assert content1 in output_file.with_suffix('.txt').read_text()
        assert content2 in output_file.with_suffix('.txt').read_text()
