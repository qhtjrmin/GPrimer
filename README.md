# GPrimer
GPrimer: a fast GPU-based pipeline for primerdesign for qPCR experiments

==================================

# 1. Running
## 1.1. Prepare input sequence data
 The mRNA sequence DBs for six species (human, mouse, rat, zebrafish, cow, and pig) from the NCBI Reference Sequence (RefSeq) database (http://www.ncbi.nlm.nih.gov/refseq) were used as input in our experiments. They have NM as the prefix of GenBank accession number. The DBs that we used are in “input_sequence” folder (human and mouse sequence were excluded due to GitHub capacity issues). If users want to take their own sequence as the input, they can use it. The format of input should be ‘sid(\t)sequence’. 
 
- The example of input:
```
1	CCAGGCATTGGGACAGTCGTATGTATCTGGCAAGAGCAAAGCTGCA…
2	ATCGGCGCTGGGTAGGCGCTCTGGTGCTCGCCGAGGACACTCCGCT…
3	GCGGGCACAGTCCGGGAGCCGCTGTCGCCCGGCCAGCCTGAGGCGG…
```

## 1.2. Get GPrimer execution file and run
- The query for running GPrimer:
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
  
 - The example of query:
```
./gprimer –i ./input_sequence/s_scrofa_refseq_181107.txt –o ./inter/output.txt –d ./inter/ 
```

## 1.3. The parameters for filter
 The number of CPU threads and GPUs that users want to exploit can be set. They will be limited by the performance of the machine that used in your experiments. In our experiment, the number of CPU threads exploited was fixed as 20 and the number of GPUs exploited is changed for performance comparison as varying the number of GPUs.
- The example command for using 4 GPUs:
```
./gprimer –i ./input_sequence/s_scrofa_refseq_181107.txt –o ./inter/output.txt –d ./inter/ -g 4
```
- The example command for using 8 GPUs:
```
./gprimer –i ./input_sequence/s_scrofa_refseq_181107.txt –o ./inter/output.txt –d ./inter/ -g 8
```

# 2. The parameters for filtering in GPrimer
## 2.1. The parameters for single filtering
The default parameters for single filtering are shown in the below table. It is the same with the online default parameters in MRPrimerW2 (http://mrprimerw2.com)

|Parameters|Values|
|:---:|:---:|
|primer length (bp)|19-23|
|melting temperature (TM,℃)|58-62|
|GC content(%)|40-60|
|self-complementarity|<5-mer|
|3' self-complementarity|<4-mer|
|contiguous residue|<6-mer|
|Gibbs free energy (∆G, kcal/mol)|>=-9|
|hairpin|<3-mer|

If users wat to change the parameters, set the option -p1 as 1. Then, they can enter the paramter values that they want.

The example command for setting your own single filtering parameters:
```
./gprimer -i ./input_sequence/s_scrofa_refseq_181107.txt -o ./inter/output.txt -d ./inter/ -p1 1
```

## 2.2. The parameters about pair filtering
 The default parameters for pair filtering is shown in the below table. It is the same with the online default parameters in MRPrimerW2 (http://mrprimerw2.com)
 
|Parameters|Values|
|:---:|:---:|
|length difference|<=5-mer|
|TM difference (℃)|<=3|
|product size (bp)|100-250|
|pair-complementarity|<5-mer|
|3' pair-complementarity|<4-mer|

If users want to change the parameters for pair filtering, set the option –p2 as 1. Then, they can enter the parameter values that they want.

The example command for setting your own single filtering parameters:
```
./gprimer -i ./input_sequence/s_scrofa_refseq_181107.txt -o ./inter/output.txt -d ./inter/ -p2 1
```

## 3. The format of output file
 The final output file is stored in the output path that you set. In GPrimer, the output file is written by CPU threads in Step 5 and thyen it is sorted after being gathered as one file. The format of output file is shown in the below.
 
 The example of output file:
 ```
 1   ATTGACGATGCTTGGCGAGA+GTAAGGCCAAGTCAGTCCACT+1+94+238 8.458069801330566
1   ACAGTGAGGTACAGAGAGATGGA+GTAAGGCCAAGTCAGTCCACT+1+56+238  8.513699531555176
1   TGAGGTACAGAGAGATGGAAGGA+GTAAGGCCAAGTCAGTCCACT+1+60+238  8.522356033325195
…				…				…
1011    TCTGAAGAGCCGACGAGTC+CAAGAGCAGCAACACCAGGG+1011+14+116    12.908910751342773
1011    TCTGAAGAGCCGACGAGTC+CAGAGCATAAGGAGCCCGGA+1011+14+194    12.910165786743164
1011    TCTGAAGAGCCGACGAGTC+GGCAAGAGCAGCAACACCA+1011+14+118 13.148364067077637
1011-960    CAATCCCAACACAACCAACGC+ACACGCTGACATTCACCTCCT+1011+906+1054   7.154574394226074
1011-960    CAATCCCAACACAACCAACGC+ACACGCTGACATTCACCTCCT+960+874+1022    7.154574394226074
1011-960    GGAGCTGTCTCGTCTTATGCA+TGGTCATCCCCACAAAAGCT+1011+205+389 7.270344257354736
```

The second columns represents ‘f.P+r.P+sid+f.pos+r.pos’. The f.P is forward primer and r.P is reverse primer. The sid means the target sequence id of f.P and r.P. The f.pos and r.pos are the target position of f.P and r.P in the sid. The first column is the common target sidset of the f.P and r.P. The last column represents penalty score of the primer pair (f.P and r.P) and primer pairs with low scores have high rank for the corresponding target sequence. Thus, users can identify top-1 or top-n primer pairs for the target sequences.
