#!/usr/bin/env python


files = R1, R2, R3, R4

index_dict = {48 keys and values} #24 indexes normal and 24 indexes rev comp as keys and their corresponding values 

open R1, R2, R3, R4 as read1, index1, index2, read2: #R1-4 is the illumina read files

    open files 1-52 as output1 to output52: #I will have to write out the 52 output files that are needed to write to.

        read 4 lines: #Read 4 lines in all of the files at once and store that into a dict

        dict = {read1: [4 lines], index1: [4 lines], index2: [4 lines], read2: [4 lines} #store all lines in a dict with a list as a value of the 4 lines

        check quality score mean for all three store as variables in the dictionary

            if index1 in index_dict: #Check to see if index one is in the dictionary of all possible indexes

                if index2 == index_dict[index1]: #check to see if the index2 matches the values of the dict
                    add index1 + "-" + index2 to header in read1 and read2 #Add the indexes to the header of both reads
                    
                    if quality score is at or above cut off for all things:
                        write read1 to read1file and read2 to read2file in bucket of index1_match #Write this to the appropriate files in the correct index bucket

                    else:
                        write to junk

                else: #If they don't match then it indexed hopped
                    add index1 + "-" + index2 to header in read1 and read2 #Add the headers to the 2 reads
                    if quality score is at or above cut off for all things:
                        write these into index hopping #Write these to the appropriate files in the index hopping bucket

                    else:
                        write to junk

            else: #If the index is not in the index_dict then it is trash
                add index1 + "-" + index2 to header in read1 and read2 #Add the indexes to the two reads

                write to trash_bucket #Write these to the appropriate files in the trash bucket
