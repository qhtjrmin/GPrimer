# GPrimer
GPrimer: a fast GPU-based pipeline for primerdesign for qPCR experiments
Paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04133-4

# 1. Get GPrimer
To get and compline GPrimer, follow the quires that listed below
```
$ git clone https://github.com/qhtjrmin/GPrimer.git
```
You can directly use the binary file in 'execution_file' folder for execution of GPrimer without compilation.

## Environment
GPrimer needs the following software to run on the system:

- CUDA toolkit version 8 or higher.
- Nvidia driver (v384 or higher)
- gcc/g++ 4.8.x or later


# 2. Running
## 2.1. Prepare input sequence data
 The input DBs that can be used are in “input_sequence” folder. In the case of mouse and human RefSeq data, they should be decompressed as the following command to use for input sequence DB.
 
 ```
 $ gunzip m_musculus_refseq_181107.txt.gz 
 $ gunzip h_sapiens_refseq_181107.txt.gz
 ```
 
 If you want to take your own sequence as the input, then use it. The format of input should be ‘sid(\t)sequence’. 
 
- The example of input:
```
1	CCAGGCATTGGGACAGTCGTATGTATCTGGCAAGAGCAAAGCT…
2	ATCGGCGCTGGGTAGGCGCTCTGGTGCTCGCCGAGGACACTCC…
3	GCGGGCACAGTCCGGGAGCCGCTGTCGCCCGGCCAGCCTGAGG…
```

## 2.2. Get GPrimer execution file and run
- The query for running GPrimer:
```
$ ./gprimer -i <input_sequence> -o <final_output_path> -d <working_directory> -t <num_of_threads(option)> -g <num_of_gpus(option)> -w <is_write_intermidiate_files(option)> -b <buffer_size(option)> -p1 <is_change_parameters_single(option))> -p2 <is_change_parameters_pair(option)>
```
- Parameters
  - input_sequence: path of the input sequence file (necessary)
  - final_output_path: path of the output file (necessary)
  - working_directory: path of working directory (necessary)
  - num_of_threads: the number of CPU threads exploited (in default 20)
  - num_of_gpus: the number of GPUs exploited (in default 1)
  - is_write_intermidate_files: whether it writes the intermidiate files (0: no, 1: yes, in default 0)
  - buffer_size: the buffer utilization percetange(%) of main memory for linux sort (in default 30, do not write '%' together)
  - is_change_parameters_single: whether it changes the paramerters for single filtering (0: no, 1: yes, in default 0)
  - is_change_parameters_pair: whether it changes the parameters for pair filtering (0: no, 1: yes, in defaul 0)
  
- The example of query:
```
$ ./gprimer –i ./input_sequence/s_scrofa_refseq_181107.txt –o ./inter/output.txt –d ./inter/ 
```

## 2.3. The parameters for filter
 The number of CPU threads and GPUs that you want to exploit can be set. They will be limited by the performance of the machine that used in the experiments. 
 
- The example command for using 4 GPUs:
```
$ ./gprimer –i ../input_sequence/s_scrofa_refseq_181107.txt –o ./inter/output.txt –d ./inter/ -g 4
```
- The example command for using 8 GPUs:
```
$ ./gprimer –i ../input_sequence/s_scrofa_refseq_181107.txt –o ./inter/output.txt –d ./inter/ -g 8
```

# 3. The parameters for filtering in GPrimer
## 3.1. The parameters for single filtering
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

If you want to change the parameters, set the option -p1 as 1. Then, you can enter the paramter values that you want.

The example command for setting your own single filtering parameters:
```
$ ./gprimer -i ../input_sequence/s_scrofa_refseq_181107.txt -o ./inter/output.txt -d ./inter/ -p1 1
```

## 3.2. The parameters about pair filtering
 The default parameters for pair filtering is shown in the below table. It is the same with the online default parameters in MRPrimerW2 (http://mrprimerw2.com)
 
|Parameters|Values|
|:---:|:---:|
|length difference|<=5-mer|
|TM difference (℃)|<=3|
|product size (bp)|100-250|
|pair-complementarity|<5-mer|
|3' pair-complementarity|<4-mer|

If you want to change the parameters for pair filtering, set the option –p2 as 1. Then, you can enter the parameter values that you want.

The example command for setting your own single filtering parameters:
```
$ ./gprimer -i ../input_sequence/s_scrofa_refseq_181107.txt -o ./inter/output.txt -d ./inter/ -p2 1
```

## 4. The format of output file
 The final output file is stored in the output path that you set. In GPrimer, the output file is written by CPU threads in Step 5 and then it is sorted after being gathered as one file. The format of output file is shown in the below.
 
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

The first column represents common target sidset of the primer pair (forward primer and reverser primer) in second column. The second column represents ‘f.P+r.P+sid+f.pos+r.pos’. The f.P is forward primer and r.P is reverse primer. The sid means the target sequence id of f.P and r.P. The f.pos and r.pos are the target position of f.P and r.P in the sid. The last column represents penalty score of the primer pair (f.P and r.P). Here, primer pairs with low scores have high rank for the corresponding target sequence. Thus, you can identify top-1 or top-n primer pairs for the target sequences.
