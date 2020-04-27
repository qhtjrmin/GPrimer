# GPrimer
GPrimer: a fast GPU-based pipeline for primerdesign for qPCR experiments

# Running
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
