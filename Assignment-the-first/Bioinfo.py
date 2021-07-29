
def add_indexes_to_reads(file1: str, file2: str, file3: str, file4: str) -> str:
    """Open all of the files and iter through each of them line by line. Adding the indexes
    from R2 and R3 to the header of File1 and File4 lines."""

    return "None"

def bucket_addition(header1: str, header2: str) -> str:
    """Take the header for each read after the addition from add_indexes_to_reads and check
    to make sure that the indexes are correctly dual unique matching. If not put into bad or index swapped.
    If it passes the correct parameters then put the indexes into its corresponding buckets and readfile orientation.
    Example R1/R4"""

    return "None"

def correct_error_singlebase_call(sequence1: str, sequence2: str) -> str:
    """Test to see if there is an N in the either sequence. If so then at the position of the N in the string,
    check the same reverse position in the opposite string. Ex. N is at pos 25 in sequence1 the pos to look at
    in sequence2 is -[8 - 2] which gives pos -6. Then you take that nucleotide and run it against your revdict to fill in."""

    rev_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

    return "None"



