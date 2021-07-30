#!/usr/bin/env python

from numpy import number
import Bioinfo
import gzip
import argparse
import matplotlib.pyplot as plt

def get_args():
    parsers = argparse.ArgumentParser(description='File that needs to be added to compute stats.')
    parsers.add_argument('-f', '--file', help='reads.fastq.gz', type=str, required=True)
    parsers.add_argument('-n', '--number', help='Put the number of nucleotides of the reads. Default = 101bp', default=101, type=int)
    return parsers.parse_args()

args = get_args()

file = args.file
file_ouptut = file.split('/')[-1]

print(file)
print(file_ouptut)
def init_list(lst: list=[], read_length=args.number, value=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''

    for i in range(read_length):
        if value == list():
            value = list()
        lst.append(value)
    return lst
    

def populate_list(load_file):
    """Update with your own docstring"""
    initial = init_list([])
    with gzip.open(load_file, 'rt') as open_file:
        file_pos = 0 
        for line in open_file:
            if file_pos % 500000 == 0:
                print(f'{file_pos} lines done')
            file_pos +=1        
            if file_pos % 4 == 0:
                for i in range(0, len(initial)):
                    initial[i] += Bioinfo.convert_phred(line.strip('\n')[i])
    return initial, file_pos

my_list, num_lines = populate_list(args.file)

mean = [(my_list[i]/(num_lines/4)) for i in range(len(my_list))]

# def variance_and_stdev(load_file, mean_list, record_numbers):
#     """Update with your own docstring"""
#     var = init_list([])
#     stdev = init_list([])
#     record_numbers = record_numbers/4
#     with gzip.open(load_file, 'rt') as open_file:
#         file_pos = 0 
#         for line in open_file:
#             # if file_pos % 500000 == 0:
#             #     print(f'{file_pos} lines done')
#             file_pos +=1        
#             if file_pos % 4 == 0:
#                 for i in range(0, len(var)):
#                     var[i] += ((Bioinfo.convert_phred(line.strip('\n')[i]) - mean_list[i])**2)
#     new_var = [x/record_numbers for x in var]
#     stdev = [number**0.5 for number in new_var]

#     return new_var, stdev

# var, stdev = variance_and_stdev(args.file, mean, num_lines)

with open('{}-mean.txt'.format(file_ouptut), 'w') as f:
    
    f.write('# Base Pair\tMean\n')
    for i in range(len(my_list)):
        f.write(str(i) + '\t' + str(mean[i])+'\n')

x = [x for x in range(args.number)]
y = [(my_list[i]/(num_lines/4)) for i in range(len(my_list))]

plt.figure(figsize=(20,10))
plt.bar(x, mean, align='center', alpha=0.5, ecolor='black')
plt.title('Mean qScore for each bp position in {}'.format(file_ouptut))
plt.xlabel('BP position')
plt.ylabel('Mean qScore')
plt.savefig('{}-stats.png'.format(file_ouptut))