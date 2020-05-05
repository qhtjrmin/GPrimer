# GPrimer
GPrimer: a fast GPU-based pipeline for primerdesign for qPCR experiments

# 1. Running
## 1.1. Prepare input sequence data
 The mRNA sequence DBs for six species (human, mouse, rat, zebrafish, cow, and pig) from the NCBI Reference Sequence (RefSeq) database (http://www.ncbi.nlm.nih.gov/refseq) were used as input in our experiments. They have NM as the prefix of GenBank accession number. The DBs that we used are in “input_sequence” folder (human and mouse sequence were excluded due to GitHub capacity issues). If users want to take their own sequence as the input, they can use it. The format of input should be ‘sid(\t)sequence’. 
 
The example of input:
```
1	CCAGGCATTGGGACAGTCGTATGTATCTGGCAAGAGCAAAGCTGCA…
2	ATCGGCGCTGGGTAGGCGCTCTGGTGCTCGCCGAGGACACTCCGCT…
3	GCGGGCACAGTCCGGGAGCCGCTGTCGCCCGGCCAGCCTGAGGCGG…
```

## 1.2. Get GPrimer execution file and run
The query for running GPrimer:
```
./gprimer -i <input_sequence> -o <final_output_path> -d <working_directory> -t <num_of_threads(option)> -g <num_of_gpus(option)> -w <is_write_intermidiate_files(option)> -p1 <is_change_parameters_single(option))> -p2 <is_change_parameters_pair(option)>
```
- Parameters
  - input_sequence: path of the input sequence file (necessary)
  - final_output_path: path of the output file (necessary)
  - working_directory: path of working directory (necessary)
  - num_of_threads: the number of CPU threads exploited (in default 20)
  - num_of_gpus: the number of GPUs exploited (in default 1)
  - is_write_intermidate_files: whether it writes the intermidiate files (0: no, 1: yes, in default 0)
  - is_change_parameters_single: whether it changes the paramerters for single filtering (0: no, 1: yes, in default 0)
  - is_change_parameters_pair: whether it changes the parameters for pair filtering (0: no, 1: yes, in defaul 0)
