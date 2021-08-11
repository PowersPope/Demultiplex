#!/bin/env python
import argparse
import Bioinfo
import numpy as np
import os 
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Takes the 4 sequencing files from an Illumina run and demultiplexes them.")
    parser.add_argument("-f1", "--file1", help="Read1 fastq file.", required=True, type=str)
    parser.add_argument("-f2", "--file2", help="Index1 fastq file.", required=True, type=str)
    parser.add_argument("-f3", "--file3", help="Index2 fastq file.", required=True, type=str)
    parser.add_argument("-f4", "--file4", help="Read2 fastq file.", required=True, type=str)
    parser.add_argument("-i", "--indexes", help='A file where each line is an index, no header or footer', required=True, type=str)
    parser.add_argument("-c", "--cutoff", help="Add quality score cutoff")
    return parser.parse_args()

# Pass the command line specified files to the arg parser
args = get_args()


####################################### Create the file dict and index-pair dicts that we will access #########################################################################

# Instantiate the dictionary that will hold my records
holding_dict = {'read1': [], 'index1': [], 'index2': [], 'read2': []}


# Create dict with indexes as keys and a tuple with (handle, open(outputfile, 'w')) as values
file_index_dict ={} #Empty dict to be filled
index_set = set() # Create an empty set that will hold the unique indexes
with open(args.indexes, 'r') as index_file:
    header = index_file.readline() #Get rid of the header line so that columns are in the dictionary
    for line in index_file:
        # Pull out the handle (Ex. B1, A11)
        index_handle = (line.strip()).split('\t')[3]
        # Pull out the index 
        index_string = (line.strip()).split('\t')[4]
        # This is to cehck and see if the directory exists already. If not, then create it. If not then do nothing.
        file_string = 'index{}/'.format(index_handle)
        os.makedirs(file_string, exist_ok=True)

        # Create the key value pair described at the top of this. The first one is the forward index, the second one is the reverse comp of it
        file_index_dict[index_string] = (index_handle, open(file_string + 'read1.fastq', 'w'))
        file_index_dict[Bioinfo.reverse_compliment(index_string)] = (index_handle, open(file_string + 'read2.fastq', 'w'))
        
        # Add all of the unique indexes to a set
        index_set.add((line.strip()).split('\t')[4])


# All possible index pairs (dual and swapped) dictionary that holds the counts for them.
new_bin_dict = {} # New empty dict for the counts
access_index = list(index_set) # Make an iterable list of the set
storage_index = set() # Set to keep the combined dual indexes together
# Create all possible pairs 24 * 24 (576)
for foo in range(0, len(index_set)): # First variable in the pair
    for red_foo in range(0, len(index_set)): # Second variable in the pair
        combine_index = access_index[foo] + "-" + Bioinfo.reverse_compliment(access_index[red_foo]) # Combine them together
        new_bin_dict[combine_index] = 0 # Add to dict with a 0 value
        if access_index[foo] == access_index[red_foo]: # Create a set of all the unique dual indexes
            storage_index.add(combine_index) # Add to a set if it is completely unique dual index (This is for later stats)

# Add in the bad and index_hopping file paths
os.makedirs('bad', exist_ok=True)
file_index_dict['bad1'] = ('bad', open("bad/read1.fastq", 'w'))
file_index_dict['bad2'] = ('bad', open("bad/read2.fastq", 'w'))
os.makedirs('index_hopping', exist_ok=True)
file_index_dict['index_hopping1'] = ('index_hopping', open("index_hopping/read1.fastq", 'w'))
file_index_dict['index_hopping2'] = ('index_hopping', open("index_hopping/read2.fastq", 'w'))


# Add the bad and index_hopping counts to the count dict
new_bin_dict['bad'] = 0
new_bin_dict['index_hopping'] = 0

#######################  Sort records into their buckets/files loop ##############################################################################

# Open all of the command line passed files.
with gzip.open(args.file1, 'rt') as read1, gzip.open(args.file2, 'rt') as index1, gzip.open(args.file3, 'rt') as index2, gzip.open(args.file4, 'rt') as read2:
    
    # Instantiate the line count
    count = 0
    
    #Zip all the lines of the 4 files together so that they will be iterated through together.
    for r1, i1, i2, r2 in zip(read1, index1, index2, read2):
        #Increment the line number
        count += 1

        #If the count is % 4 then you have 4 lines in your dictionary and it can now be sorted into groups.
        if count % 4 == 0:

            #Add in the quality score and then the dictionary is read to be processed.
            holding_dict['read1'].append(r1.strip())
            holding_dict['index1'].append(i1.strip())
            holding_dict['read2'].append(r2.strip())
            holding_dict['index2'].append(i2.strip())

            
            # Make the header index addition that can be compared later
            # reverse_index = Bioinfo.reverse_compliment(holding_dict['index2'][1])
            index_addition = holding_dict['index1'][1] + '-' + holding_dict['index2'][1]

            
            # Add the index to the qname 
            holding_dict['read1'][0] = holding_dict['read1'][0] + ' ' + index_addition
            holding_dict['read2'][0] = holding_dict['read2'][0] + ' ' + index_addition
            
            # Check to make sure that the indexes and the reads pass the quality score that we have set 
            if ('N' not in holding_dict['index1'][1] and 'N' not in holding_dict['index2'][1]) and index_addition in new_bin_dict:
                
                # Check to make sure that the index isn't hopping
                if file_index_dict[holding_dict['index1'][1]][0] == file_index_dict[holding_dict['index2'][1]][0]:
                    

                    # Convert quality scores to qScore then take the mean and run check variables
                    qual_index1 = np.mean([Bioinfo.convert_phred(num) for num in holding_dict['index1'][3]])
                    qual_index2 = np.mean([Bioinfo.convert_phred(num) for num in holding_dict['index2'][3]])
                    # Check against cutoff
                    mean1_test = qual_index1 >= int(args.cutoff)
                    mean2_test = qual_index2 >= int(args.cutoff)


                    # Check to see if the mean qScore is over 30 (I could add a modular part to this instead of it being hardcoded)
                    if mean1_test or mean2_test:

                        ####### Add to the appropriate index dual pair file
                        for line_f1 in holding_dict['read1']:
                            file_index_dict[holding_dict['index1'][1]][1].write(line_f1 + '\n')
                        for line_f2 in holding_dict['read2']:
                            file_index_dict[holding_dict['index2'][1]][1].write(line_f2 + '\n')
                        # Increment the dual index specific count
                        new_bin_dict[index_addition] += 1
                    else:

                        ###### Add to the bad files if they fail the mean qScore check
                        for line_f1 in holding_dict['read1']:
                            file_index_dict['bad1'][1].write(line_f1+'\n')
                        for line_f2 in holding_dict['read2']:
                            file_index_dict['bad2'][1].write(line_f2+'\n')
                        # Increment the bad specific counter
                        new_bin_dict['bad'] += 1
                else:

                    ##### Add to index hopping files if they don't pass the reverse comp/dual index test
                    for line_f1 in holding_dict['read1']:
                        file_index_dict['index_hopping1'][1].write(line_f1+'\n')
                    for line_f2 in holding_dict['read2']:
                        file_index_dict['index_hopping2'][1].write(line_f2+'\n')
                    # Increment the index_hopping total and the specific index hopped pair total
                    new_bin_dict['index_hopping'] += 1
                    new_bin_dict[index_addition] += 1
            else:

                ##### Add to bad files if they have Ns in index
                for line_f1 in holding_dict['read1']:
                    file_index_dict['bad1'][1].write(line_f1+'\n')
                for line_f2 in holding_dict['read2']:
                    file_index_dict['bad2'][1].write(line_f2+'\n')
                # Increment bad file count
                new_bin_dict['bad'] += 1

            # Now that all of the processing is done. This will make the dictionary values blank again. So we can restart the process. 
            # Much better computationally to .clear() then to just make a new list. The old list would have stayed in livrary and slowed everything down.
            holding_dict['read1'].clear()
            holding_dict['index1'].clear()
            holding_dict['read2'].clear()
            holding_dict['index2'].clear()

        else:
            #We don't have enough lines into the dictionary values, add line to list value.    
            holding_dict['read1'].append(r1.strip())
            holding_dict['index1'].append(i1.strip())
            holding_dict['read2'].append(r2.strip())
            holding_dict['index2'].append(i2.strip())
            
# Close all of the opened files
for key in file_index_dict.keys():
    file_index_dict[key][1].close()


#Stats file
count = count/4 # Get accurate record numbers
total_index = 0 
for number_of_indexpair in list(storage_index): 
    total_index += new_bin_dict[number_of_indexpair] # Sum all dual index pairs

# Write to the stats file.
with open('total_stats.txt', 'w') as stats:
    stats.write('Stats Analysis\n\n')
    stats.write('Count Total: ' + str(count) + '\n')
    stats.write('Cutoff: ' + args.cutoff + '\n')
    stats.write('Index\tCount\tPercent\n')
    for key in new_bin_dict.keys():
        stats.write(key + "\t" + str(new_bin_dict[key]) + "\t" + str(round(new_bin_dict[key]/count * 100, 2)) + '%\n')
    stats.write('Total Dual Indexes\t' + str(total_index) + '\t' + str(round(total_index/count * 100,2)) + '%\n' )
    
# Yay! We Finished!
print("All closed!")