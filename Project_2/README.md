## Project 2 - Designed and examined DNA strings (probes) with fuzzy matching and simulations

#### *Languages and Tools: Part I*
* Python
  * csv
  * argparse
  * itertools
  * fuzzywuzzy (python-Levenshtein)
  * numba
  * NumPy
  * typing
* Bash
  * Slurm
---------------
### Design Contiguous Strings and Get a Cartesian Product (Equivalent to an SQL Cross-Join)
```python3
# Get contiguous substrings from given strings, create Cartesian product, get fuzzy scores
def get_substring(string, len_k):
    '''Function to get all substrings of a specific length (k) -- ie string of nucleotides (probes) and a specified string length (len_k)'''
    substr = [string[x: y] for x in range(len(string)) for y in range(x + 1, len(string) + 1) if len(string[x:y]) == len_k]
    return substr

# add flags for running .py script on the command line
parser = argparse.ArgumentParser(description="Program to calculate cartesian product and filter by Levenshtein scores.")

parser.add_argument("-p", "--strings", help="Provide a file containing [redacted] sequences (do not include a header row) that you wish to turn into strings.", required=True, type=str)
parser.add_argument("-k", "--kmer", help="Specify the length of the thing.", required=True, type=int)
parser.add_argument("-c", "--cdr3s", help="Provide a file containing some things (do not include a header row) that you wish to test for Levenshtein scores.", required=True, type=str)
parser.add_argument("-t", "--threshold", help="Specify a Levenshtein threshold (between 0 and 100). The threshold is not inclusive.", required=True, type=int)
parser.add_argument("-o", "--out", help="Name output file.", required=False, type=str)

args = parser.parse_args()

with open(STRINGS, newline='') as f1, open(THINGS, newline='') as f2:
    reader1 = csv.reader(f1)
    reader2 = csv.reader(f2)

    thing1 = list(reader1) # Convert to lists
    thing2 = list(reader2)

#######################

# get list of sequences
sequences = [] # list for the test cdr3s, contains all 129 enhanced TCRs and 100 others
for i in thing1:
    for j in i:
        sequences.append(j)

seq_dict = {}
# Use get_substring function to get each contiguous substring of each seq
for i in CDR3s:
    if i not in CDR3_dict:
        seq_dict[i] = get_substring(i, KMER)
    else:
        break

# Add the potential strings to a list
thing_list = []
for k, v in seq_dict.items():
    for i in v:
        thing_list.append(i)

########################

f3 = open(OUT, 'w') # Open a file that will hold k-v pairs that meet the threshold

string_array = np.array(thing_list).flatten()
thing_array = np.array(thing2).flatten()

for cartesian_tuples in chunked_cartesian_product(*[string_array, thing_array], chunk_size=5):
    for string, seq in cartesian_tuples:
        if fuzz.partial_ratio(string, seq) > THRESH:
            f3.write(str(string).strip('[]') + ',' + str(seq).strip('[]') + ',' + str(fuzz.partial_ratio(string, seq)) + '\n')

f3.close()
```
### Get Fuzzy Scores by Running .py Script on Cluster with Slurm (Workload Manager)
```bash
#!/bin/bash

#SBATCH --job-name=some_job
#SBATCH --output=slurm-%j-%x.out
#SBATCH --nodes=1                ### Node count required for the job
#SBATCH --ntasks=1		### Cores
#SBATCH --time=15-00:00:00       ### Days-HH:MM:SS
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=100
#SBATCH --mail-user=name@email.com
#SBATCH --mail-type=ALL

source deactivate
source activate my_root

cd /home/directory

STRINGS=/home/STRINGS.csv
THINGS=/home/THINGS.csv
OUT=/home/RESULTS.csv

/usr/bin/time -v ./GET_SEQS_PASS_THRESH.py -p $STRINGS -k 30 -c $THINGS -t 90 -o $OUT
```
---------------
#### *Languages and Tools: Part II*
* Python
  * pandas
  * scikit-learn
  * argparse
  * NumPy
  * SciPy
  * itertools
  * csv
  * multiprocessing
  * fuzzywuzzy (python-Levenshtein)
* Bash
  * Slurm
---------------
### Logistic Regression with Best Strings from Cartesian Product
```python3
parser = argparse.ArgumentParser(description="Program to produce a model for things of differing lengths and Levenshtein distances.")

parser.add_argument("-p", "--thing", help="Specify thing length. It corresponds with the input thing file.", required=True, type=int)
parser.add_argument("-i", "--input", help="Input file containing all things with a Levenshtein distance of at least 90.", required=True, type=str)
parser.add_argument("-s", "--seqs", help="Provide the path to the sequences.", required=True, type=str)
parser.add_argument("-e", "--enhanced", help="Proved the path to the dataset containing all sequences.", required=True, type=str)
parser.add_argument("-ct", "--counts", help="Specify how many sequences one thing can target.", required=True, type=int)

args = parser.parse_args()

# read in enhanced seqs
df = pd.read_csv(SEQS, names=['seq_thing'], low_memory=False)
enhanced = list(enhanced_df['seq_thing'])

df = pd.read_csv(ENHANCED, low_memory=False, usecols=['seq_thing','sample','status'])
encode = {True:1, False:0}
df['status'] = df['status'].map(encode)
CMV_df = pd.merge(enhanced_df,df)

# read in filtered df
kmer = pd.read_csv(INPUT, low_memory=False, names=['things','seq_thing','edit_dist'])
kmer = kmer[kmer['seq_thing'].map(len) >= 18] # drop any seqs smaller than the smallest thing

# reorder df by thing and split into individual dfs
kmer_levs = [v for k, v in kmer.groupby('things')]

things = kmer['things'].unique()
things = pd.DataFrame(things, columns=['things'])

# merge the seqs with the input seqs (output from Levenshtein algorithm)
dfs = []
for i in range(len(kmer_levs)):
    dfs.append(df)
    kmer_levs[i] = kmer_levs[i].merge(dfs[i], on='seq_thing')

kmer_levs2 = [v for k, v in kmer_df1.groupby('sample')]

dfs2 = []
for i in range(len(kmer_levs2)):
    # dfs2.append(things)
    kmer_levs2[i] = kmer_levs[i].merge(things)

kmer_df2 = pd.concat(kmer_levs2)

grouped = kmer_df1.groupby(['sample'], as_index=False).agg({'seq_thing':'nunique',
                                                    'things':'nunique',
                                                    'edit_dist':'nunique',
                                                    'status':'first'})

df_sub = df[df['seq_thing'].str.contains('|'.join(enhanced), na=False)]

# list of dataframes by sample
dfs_sub = [v for k, v in df_sub.groupby('sample')]

thing_ct = []
status = []
unique_seqs = []

# get sample names
samples = list(df_sub['sample'].unique())

# get dataframe of all seqs, but only select the sample names currently in the dataset
full_samples =  df[df['sample'].isin(samples)]
# split the above dataframe into a list of dataframes by sample name
full_sample_dfs = [v for k, v in full_samples.groupby('sample')]

# get unique counts
for i in range(len(full_sample_dfs)):
    thing_ct.append(len(pd.merge(enhanced_df, full_sample_dfs[i])))
    status.append(full_sample_dfs[i]['status'].iloc[0])
    unique_seqs.append(len(full_sample_dfs[i]['seq_thing'].unique())) # CHANGE TO FULL_SAMPLE_DFS LATER

features = pd.DataFrame(
{'thing_Count':thing_ct,
 'status':status,
 'sample':samples,
 'Unique_seqs':unique_seqs})
```
### Get Strings of Specified Lengths Using the Above .py Script
```bash
#!/bin/bash

#SBATCH --job-name=feats
#SBATCH --output=slurm-%j-%x.out
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --ntasks=1		### Cores
#SBATCH --time=15-00:00:00       ### Days-HH:MM:SS
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=100
#SBATCH --mail-user=name@email.com
#SBATCH --mail-type=ALL

source deactivate
source activate my_root

cd /home/model

INPUT=/home/model/input.csv
SEQS=/home/model/seqs.csv
THING=/home/model/thing.csv

/usr/bin/time -v ./KMER_MODEL.py -p 30 -i $INPUT -s $SEQS -e $THING
```
