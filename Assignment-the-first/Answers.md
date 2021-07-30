# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 |
| 1294_S1_L008_R4_001.fastq.gz | Read 2 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    !['R1 Graph'](./distribution/1294_S1_L008_R1_001.fastq.gz-stats.png)
    !['R2 Graph'](./distribution/1294_S1_L008_R2_001.fastq.gz-stats.png)
    !['R3 Graph'](./distribution/1294_S1_L008_R3_001.fastq.gz-stats.png)
    !['R4 Graph'](./distribution/1294_S1_L008_R4_001.fastq.gz-stats.png)

    2. ```20 qScores are a favorable cutoff for our purposes. As long as the whole read has a 20 qScore overall, then it will be of high quality enough to proceed forward with. A qScore of 20 will have a 1 in 100 chance of having an incorrect read associated with it. Which is good for the enitre read, since the there are 101 bps within the read as a whole. If only 1 nucleotide is wrong, that will be corrected for by the other reads which should have the correct nucleotide.```
    
    3. 
    ```
    zcat 1294_S1_L008_R[2,3]_001.fastq.gz |sed -n '2~4p' | awk '$0 ~ /N/ {sum+=1} END{print sum}' 
    7304664
    ```
    
## Part 2
1. Define the problem

We have 4 reads from the illumina sequencer flow cell run. R1, R2, R3, R4. You need to combine the index and the read, so that way we can see where the read came from. As the index were assigned before amplification and sequencing. They should be able to give us more information. 

2. Describe output

We have to put these reads into their index associated bins. There will be a forward read file, and a reverse comp read file. If there is ine swapping we will put those reads into their own bucket/files. If the read qualitys do not meet our quality score cutoff then they will be put into a trash bucket. 

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).

4. Pseudocode

Look at pseudocode.py which is included in the same github repo.

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement


Look at the Bioinfo.py to check high level functions