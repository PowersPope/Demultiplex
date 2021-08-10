
# from _typeshed import FileDescriptor, ReadableBuffer


def add_indexes_to_reads(file1: str, file2: str, file3: str, file4: str) -> str:
    """Open all of the files and iter through each of them line by line. Adding the indexes
    from R2 and R3 to the header of File1 and File4 lines."""

    return read1_with_indexes, read2_with_indexes

# Input: 4 fastq Files

# Output: 2 reads with indexes attached to the headers of each read

def bucket_addition(header1: str, header2: str) -> str:
    """Take the header for each read after the addition from add_indexes_to_reads and check
    to make sure that the indexes are correctly dual unique matching. If not put into bad or index swapped.
    If it passes the correct parameters then put the indexes into its corresponding buckets and readfile orientation.
    Example R1/R4"""

    return addition_of_records_to_correct_bucket

# Input: headers of records

# Ouptut: puts them into buckets corresponding to the header indexes



def correct_error_singlebase_call(sequence1: str, sequence2: str) -> str:
    """Test to see if there is an N in the either sequence. If so then at the position of the N in the string,
    check the same reverse position in the opposite string. Ex. N is at pos 25 in sequence1 the pos to look at
    in sequence2 is -[8 - 2] which gives pos -6. Then you take that nucleotide and run it against your revdict to fill in."""

    rev_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

    return string_of_index_with_correct_sequence_nomore_N


# Input: NATGC, GCATG

# Output: CATGC, GCATG

def convert_phred(letter: str) -> int:
    '''Convert a phred score to an integer value'''
    return ord(letter) - 33

def reverse_compliment(sequence: str) -> str:
    """Reverse compliment a string"""
    matches = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    reverse = ''
    for num in range(1, len(sequence)+1):
        reverse = reverse + matches[sequence[-num]]
    return reverse