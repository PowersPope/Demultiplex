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
    return parser.parse_args()

# Pass the command line specified files to the arg parser
args = get_args()


####################################### Create the file dict and index-pair dicts that we will access #########################################################################

# Instantiate the dictionary that will hold my records
holding_dict = {'read1': [], 'index1': [], 'index2': [], 'read2': []}

# Create a list of the indexes with their reverse compliment next to them.
list_of_indexes = []
index_set = set()
# count_bin_dict = {} # A dictionary to store counts 
with open(args.indexes, 'r') as full_indexes:
    for line in full_indexes:
        list_of_indexes.append(line.strip())
        # count_bin_dict[line.strip()] = 0
        list_of_indexes.append(Bioinfo.reverse_compliment(line.strip()))
        index_set.add(line.strip())

# All possible pair dictionary
new_bin_dict = {}
access_index = list(index_set)
storage_index = set() # Set to keep the combined dual indexes together
count = 0
for foo in range(0, len(index_set)):
    count += 1
    for red_foo in range(0, len(index_set)):
        combine_index = access_index[foo] + "-" + access_index[red_foo]
        new_bin_dict[combine_index] = 0
        if access_index[foo] == access_index[red_foo]:
            storage_index.add(combine_index)



# Create a dictionary with file paths for the indexes
file_index_dict ={} #Empty dict to be filled
index_range = list(range(1, len(list_of_indexes)+1)) # list of n numbers of indexes
file_range = [1,2] * int(len(list_of_indexes)/2) # list of 1,2 * n/2 so that we get read1, and read2 designations
count = 0 
for new_index in list_of_indexes:
    count += 1 # increment count
    file_string = "index{}/read{}.fastq".format(index_range[0], file_range[0])
    os.makedirs(os.path.dirname(file_string), exist_ok=True)
    file_index_dict[new_index] =  open(file_string, 'w') # create the dictionary entry
    file_range.pop(0) #remove the file_range first item after every addition 
    if count % 2 == 0: # if the count is divisible by 2 then remove the first item from the index range list
        index_range.pop(0)


# Add in the bad and index_hopping file paths
file_index_dict['bad1'] = open("bad/read1.fastq", 'w')
file_index_dict['bad2'] = open("bad/read2.fastq", 'w')
file_index_dict['index_hopping1'] = open("index_hopping/read1.fastq", 'w')
file_index_dict['index_hopping2'] = open("index_hopping/read2.fastq", 'w')


# Add the bad and index_hopping counts to the count dict
new_bin_dict['bad'] = 0
new_bin_dict['index_hopping'] = 0

#######################  Sort records into their buckets/files ##############################################################################

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
            
            
            #Initalize an empty list for the qual scores to go into
            # qual_read1 = []
            # qual_read2 = []
            qual_index1 = []
            qual_index2 = []

            
            # Make the header index addition that can be compared later
            index_addition = holding_dict['index1'][1] + '-' + Bioinfo.reverse_compliment(holding_dict['index2'][1])
            
            # Add the index to the qname 
            holding_dict['read1'][0] = holding_dict['read1'][0] + ' ' + index_addition
            holding_dict['read2'][0] = holding_dict['read2'][0] + ' ' + index_addition
            
            # Check to make sure that the indexes and the reads pass the quality score that we have set 
            if ('N' not in holding_dict['index1'][1] and 'N' not in holding_dict['index2'][1]) and index_addition in new_bin_dict:

                if holding_dict['index1'][1] == Bioinfo.reverse_compliment(holding_dict['index2'][1]):
                    

                    #Check Quality score of the reads and indexes
                    for num in holding_dict['index1'][3]:
                        qual_index1.append(Bioinfo.convert_phred(num))
                    for num in holding_dict['read2'][3]:
                        qual_index2.append(Bioinfo.convert_phred(num))
                    # for num in holding_dict['read1'][3]:
                    #     qual_read1.append(Bioinfo.convert_phred(num))
                    # for num in holding_dict['index2'][3]:
                    #     qual_read2.append(Bioinfo.convert_phred(num))
                    


                    # Calculate the mean for the indexes
                    mean_i1 = np.array(qual_index1).mean()
                    mean_i2 = np.array(qual_index2).mean()
                    # mean1 = np.array(qual_read1).mean()
                    # mean2 = np.array(qual_read2).mean()

                    # Check to see if the mean qScore is over 30 (I could add a modular part to this instead of it being hardcoded)
                    if mean_i1 >= 30 or mean_i2 >= 30:

                        ####### Add to the appropriate index dual pair file
                        for line_f1 in holding_dict['read1']:
                            file_index_dict[holding_dict['index1'][1]].write(line_f1 + '\n')
                        for line_f2 in holding_dict['read2']:
                            file_index_dict[holding_dict['index2'][1]].write(line_f2 + '\n')
                        # Increment the dual index specific count
                        new_bin_dict[index_addition] += 1
                    else:

                        ###### Add to the bad files if they fail the mean qScore check
                        for line_f1 in holding_dict['read1']:
                            file_index_dict['bad1'].write(line_f1+'\n')
                        for line_f2 in holding_dict['read2']:
                            file_index_dict['bad2'].write(line_f2+'\n')
                        # Increment the bad specific counter
                        new_bin_dict['bad'] += 1
                else:

                    ##### Add to index hopping files if they don't pass the reverse comp/dual index test
                    for line_f1 in holding_dict['read1']:
                        file_index_dict['index_hopping1'].write(line_f1+'\n')
                    for line_f2 in holding_dict['read2']:
                        file_index_dict['index_hopping2'].write(line_f2+'\n')
                    # Increment the index_hopping total and the specific index hopped pair total
                    new_bin_dict['index_hopping'] += 1
                    new_bin_dict[index_addition] += 1
            else:

                ##### Add to bad files if they have Ns in index
                for line_f1 in holding_dict['read1']:
                    file_index_dict['bad1'].write(line_f1+'\n')
                for line_f2 in holding_dict['read2']:
                    file_index_dict['bad2'].write(line_f2+'\n')
                # Increment bad file count
                new_bin_dict['bad'] += 1

            #Now that all of the processing is done. This will make the dictionary values blank again. So we can restart the process.
            holding_dict['read1'] = []
            holding_dict['index1'] = []
            holding_dict['read2'] = []
            holding_dict['index2'] = []

        else:
            #We don't have enough lines into the dictionary values, add line to list value.    
            holding_dict['read1'].append(r1.strip())
            holding_dict['index1'].append(i1.strip())
            holding_dict['read2'].append(r2.strip())
            holding_dict['index2'].append(i2.strip())
            
# Close all of the opened files
for key in file_index_dict.keys():
    file_index_dict[key].close()


#Stats file
count = count/4 # Get accurate record numbers
total_index = 0 
for number_of_indexpair in list(storage_index): 
    total_index += new_bin_dict[number_of_indexpair] # Sum all dual index pairs

# Write to the stats file.
with open('total_stats.txt', 'w') as stats:
    stats.write('Stats Analysis\n\n')
    stats.write('Count Total: ' + str(count) + '\n')
    stats.write('Index\tCount\tPercent\n')
    for key in new_bin_dict.keys():
        stats.write(key + "\t" + str(new_bin_dict[key]) + "\t" + str(new_bin_dict[key]/count * 100) + '%\n')
    stats.write('Total Dual Indexes\t' + str(total_index) + '\t' + str(total_index/count * 100) + '%\n' )
    




print("All closed!")