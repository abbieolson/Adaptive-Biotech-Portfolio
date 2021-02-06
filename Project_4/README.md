## Project 4 - Tool to identify outputs between similar runs of a clinical test

#### *Languages and Tools*
* Python
  * pandas
  * NumPy
  * scipy
  * seaborn
  * Plotly
  * os
  * glob
  * collections
  * gzip
  * argparse
* Bash
---------------
### Helper Functions
```python3
def extend_pad(l, n, pad=0):
    '''Function for padding lists with 0'''
    if len(l) >= n:
        del l[n:]
    else:
        l.extend([pad] * (n - len(l)))

def get_df(dir):
    '''Function to get dataframe'''

    all_files = glob.glob(dir + "/*thing*_thing*")

    li = []
    for fi in all_files:
        df = pd.read_csv(fi,  sep='\t', skiprows=21, low_memory=False)
        # add the file name as a new column
        df['thing'] = os.path.basename(fi)
        # extract just the thing name
        df['thing'] = df['thing'].str.split('_').str[5]
        li.append(df)

    # concatenate into one large dataframe
    df = pd.concat(li, axis=0, ignore_index=True)
    return df

def get_calculated_df(df):
    '''Function to get dataframe for coefficient of varation plot'''
    # split the dataframe into individual dfs by thinge (thing sequence)
    dfs = [x for _, x in df.groupby('thing')]

    thing_list = []
    input_list = []

    for i in range(len(dfs)):
        # add the input template counts of the technical things where the thinge is identified to a list
        input_list.append(dfs[i]['inputTemplateEstimate'].to_list())

    # if a thing doesn't contain a thinge, just add a count of 0
    for i in input_list:
        extend_pad(i, 4)

    CV = []
    mean = []
    for i in input_list:
        CV.append(variation(i))
        mean.append(np.mean(i))

    calc_df = pd.DataFrame(
    {'CV':CV,
     'mean':np.log10(mean)})
    return calc_df

def df_crossjoin(df1, df2, **kwargs):
    '''Cross join (cartesian product) between two dataframes.'''
    df1['_tmpkey'] = 1
    df2['_tmpkey'] = 1

    res = pd.merge(df1, df2, on='_tmpkey', **kwargs).drop('_tmpkey', axis=1)
    res.index = pd.MultiIndex.from_product((df1.index, df2.index))

    df1.drop('_tmpkey', axis=1, inplace=True)
    df2.drop('_tmpkey', axis=1, inplace=True)

    return res

def using_npsort(df):
    '''Drop any duplicate rows that are column agnostic'''
    df1 = pd.DataFrame(np.sort(df.values, axis=1), index=df.index).drop_duplicates()
    df2 = df.loc[df1.index]
    return df2
```
### Main Script
```python3
#!/usr/bin/env python3

# import packages
from helpers import *
import argparse

### MAIN PROGRAM
#!/usr/bin/env python3

# import packages
from helpers import *
import argparse

### MAIN PROGRAM
def main():
    '''Program to get the proper input files for the differential abundance algorithm.'''
    parser = argparse.ArgumentParser(description="Program to merge lines from MIRA output files with the MIRA Excel sheet based on sequence and experiment ID and DNA sequence (enter '-h' after the script name to see all flags).")
    parser.add_argument("-fc", "--thing", help="Path to the file containing a list of thingcell IDs (include file in path).", required=True, type=str)
    parser.add_argument("-d", "--donor", help="Path to the donor directories.", required=True, type=str)
    parser.add_argument("-op", "--out_path", help="Path to eventual file that will go into the differential abundance algorithm.", required=True, type=str)

    # arguments that contain defaults and aren't required
    parser.add_argument("-fr", "--frac", help="Partial string for extracting names from files (likely 'Frac').", default="Frac", required=False, type=str)
    parser.add_argument("-fe", "--fi_ext", help="Specify the file extension of the input files (likely '.gz').", default=".gz",required=False, type=str)
    parser.add_argument("-fg", "--figs", help="Path to the directory for output figures.", required=False, type=str)
    parser.add_argument("-n", "--nrows", help="Number of header rows that comprise the summary statistics (likely 21).", default=21, required=False, type=int)

    args = parser.parse_args()

    # get summary statistics for every thing
    filenames = []
    lines = []
    splitted = []
    dfs = []

    thing_list = [line.rstrip('\n') for line in open(args.thing)]

    # grab every file in each donor directory
    for name in glob.glob(args.donor + '/' + 'donor' + '*/*/*'):
        if FRAC in name:
            filenames.append(name)

    # add the eventual column header to the filename 'thing__'
    for name in filenames:
        with gzip.open(name, 'rt', newline='') as f_in:
            lines.append('thing__' + name)
            for i in range(N):
                line = next(f_in).strip()
                line = line[1:].replace('=', '__')
                lines.append(line)

    # separate each summary stat into an individual list
    for item in lines:
        if file_ext in item:
            splitted.append([])
        splitted[-1].append(item)

    # put the sum stats in one pandas df columns, will eventual split into two
    for i in splitted:
        df = pd.DataFrame({'col':i})
        dfs.append(df)

    for i in range(len(dfs)):
        # split the columns in two and transpose so the column headers are in the correct location
        dfs[i] = dfs[i].col.str.split('__',expand=True,)
        dfs[i] = pd.concat([dfs[i][0], dfs[i][1].str.split(',', expand=True)], axis=1)
        dfs[i] = dfs[i].T

        # modify the headers
        new_head = dfs[i].iloc[0]
        dfs[i] = dfs[i][1:]
        dfs[i].columns = new_head

        # fill any None types with NaN and eventually with 0
        dfs[i] = dfs[i].fillna(value=np.nan)
        dfs[i]['thing'] = dfs[i]['thing'].str.extract(r'^(?:[^_]+_){5}([^. ]+)', expand=True)
        dfs[i] = dfs[i].replace('', 0).replace(np.nan, 0).replace('N/A', 0)

    # concatenate all individual dataframes
    merged = pd.concat(dfs)

    # identify specific donors in the concatenated dataframe
    thing_df = merged[merged['thing'].str.contains('Donor', na=False)]

    # create a specific donor column
    thing_df['Donor'] = thing_df['thing'].str.rsplit('_', 1, expand=True).drop(0, 1)

    # split the concatenated dataframe into dataframes by donor
    thing_dfs = [x for _, x in thing_df.groupby('Donor')]

    # perform an SQL style crossjoin on each dataframe
    # prepare the dataframes for the differential abundance algorithm
    count = 0
    for i in range(len(thing_dfs)):
        count += 1
        thing_dfs[i]['thingcell'] = thing_list[i]
        thing_dfs[i] = thing_dfs[i][['thing', 'thingcell']]
        thing_dfs[i] = df_crossjoin(thing_dfs[i], thing_dfs[i].copy())
        thing_dfs[i] = thing_dfs[i][thing_dfs[i]['thing_x'] != thing_dfs[i]['thing_y']]
        thing_dfs[i].insert(0, 'comparison', thing_dfs[i]['thing_x'] + '__' + thing_dfs[i]['thing_y'])
        thing_dfs[i].rename(columns={'old_col':'new_col',
                           'thing_x':'sample 1',
                           'thingcell_x':'thingcell 1',
                           'thing_y':'sample 2',
                           'thingcell_y':'thingcell 2'}, inplace=True)
        thing_dfs[i].reset_index(drop=True, inplace=True)
        thing_dfs[i] = thing_dfs[i][~thing_dfs[i].index.isin(using_npsort(thing_dfs[i][['sample 1', 'sample 2']]).index)]
        thing_dfs[i].reset_index(drop=True, inplace=True)

        # write each dataframe to a diffab text file
        thing_dfs[i].to_csv(args.out_path + '/' + f'donor{count}_sumstats.txt', index=None, sep='\t')

if __name__ == "__main__":
    main()
#### END
```
### Run the Above Script on the Cluster with Slurm
```bash
#!/bin/bash

#SBATCH --job-name=some_job
#SBATCH --output=slurm-%j-%x.out
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --ntasks=1             
#SBATCH --time=25-00:00:00      ### Days-HH:MM:SS
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=100
#SBATCH --mail-user=email@email.com
#SBATCH --mail-type=ALL

source deactivate
source activate my_root

FLOWCELL=/home/dna.corp.adaptivebiotech.com/aolson/CARTRIDGE_PROJ/INPUT_FI/thingcell_IDs.txt
DONOR=/home/dna.corp.adaptivebiotech.com/clinkem/cartridge
OUT_PATH=/home/dna.corp.adaptivebiotech.com/aolson/CARTRIDGE_PROJ/OUT

/usr/bin/time -v ./get_diffab.py -fc $FLOWCELL -d $DONOR -op $OUT_PATH
```
### Various plot scripts

```python3
alt = thing_df[thing_df['thing'].str.contains('AltReagent')]
alt = alt[~alt['thing'].str.contains('FINA')]

pcr = thing_df[thing_df['thing'].str.contains('thing*_|QFPCR')]
pcr = pcr[~pcr['thing'].str.contains('FINA')]

thing = thing_df[thing_df['thing'].str.contains('thing*_|thing')]
thing = thing[~thing['thing'].str.contains('FINA')]

thing = thing_df[thing_df['thing'].str.contains('FINA')]

thing = pd.concat([alt.assign(dataset='AltReagent'),
           pcr.assign(dataset='QFPCR'),
           thing.assign(dataset='thing'),
           thing.assign(dataset='FINA')])

thing['thing'] = pd.to_numeric(thing['thing'])

thing_dfs_plot = [x for _, x in thing.groupby('Donor')]

################################################

# box plots for each set of things by donor
fig1 = px.box(thing_dfs_plot[0], x="dataset", y="thing", color="dataset", points="all", labels={"dataset": "thing"},
             title="Donor 1: Clonlity from Replicates")
fig1.write_image(FIGS + '/' + 'donor1_thing.png')
fig1.show()

fig2 = px.box(thing_dfs_plot[1], x="dataset", y="thing", color="dataset", points="all", labels={"dataset": "thing"},
             title="Donor 2: Clonlity from Replicates")
fig2.write_image(FIGS + '/' + 'donor2_thing.png')
fig2.show()

fig3 = px.box(thing_dfs_plot[2], x="dataset", y="thing", color="dataset", points="all", labels={"dataset": "thing"},
             title="Donor 3: Clonlity from Replicates")
fig3.write_image(FIGS + '/' + 'donor3_thing.png')
fig3.show()

fig4 = px.box(thing_dfs_plot[3], x="dataset", y="thing", color="dataset", points="all", labels={"dataset": "thing"},
             title="Donor 4: Clonlity from Replicates")
fig4.write_image(FIGS + '/' + 'donor4_thing.png')
fig4.show()

###############################################

# Plots for Deliverable 1 - CLonality per Donor
don1_thing = thing_df[thing_df['thing'].str.contains('Donor1', na=False)]
don2_thing = thing_df[thing_df['thing'].str.contains('Donor2', na=False)]
don3_thing = thing_df[thing_df['thing'].str.contains('Donor3', na=False)]
don4_thing = thing_df[thing_df['thing'].str.contains('Donor4', na=False)]

fig1 = px.bar(don1_thing, x="thing", y="thing", color="thing", title="Clonality per Replicate (Donor 1)", width=1400, height=700)
fig1.write_image('/Users/aolsen/IMMUNOSEQ/figs/donor1_thing.png')
fig1.show()

fig2 = px.bar(don2_thing, x="thing", y="thing", color="thing", title="Clonality per Replicate (Donor 2)", width=1400, height=700)
fig2.write_image('/Users/aolsen/IMMUNOSEQ/figs/donor2_thing.png')
fig2.show()

fig3 = px.bar(don3_thing, x="thing", y="thing", color="thing", title="Clonality per Replicate (Donor 3)", width=1400, height=700)
fig3.write_image('/Users/aolsen/IMMUNOSEQ/figs/donor3_thing.png')
fig3.show()

fig4 = px.bar(don4_thing, x="thing", y="thing", color="thing", title="Clonality per Replicate (Donor 4)", width=1400, height=700)
fig4.write_image('/Users/aolsen/IMMUNOSEQ/figs/donor4_thing.png')
fig4.show()

#################################################

# Precision plot from full data per donor
PATH = '/Users/aolsen/IMMUNOSEQ/cartridge/'

donor1 = PATH + 'donor1/'
donor2 = PATH + 'donor2/'
donor3 = PATH + 'donor3/'
donor4 = PATH + 'donor4/'

all_files = glob.glob(PATH + "/*thing*_thing*")

donor1_df = get_df(donor1)
donor1_df['thingcell'] = '8a7a94db75d3dca90176300240ce3779'
donor1_calc = get_calculated_df(donor1_df)
donor1_calc.to_csv(PATH + 'donor1_thing97A-D.csv')

donor2_df = get_df(donor2)
donor2_df['thingcell'] = '8a7a94db75d3dca901763e90cac4495d'
donor2_calc = get_calculated_df(donor2_df)
donor2_calc.to_csv(PATH + 'donor2_thing97A-D.csv')

donor3_df = get_df(donor3)
donor3_df['thingcell'] = '8a7a94db75d3dca901763e90c96248ed'
donor3_calc = get_calculated_df(donor3_df)
donor3_calc.to_csv(PATH + 'donor3_thing97A-D.csv')

donor4_df = get_df(donor4)
donor4_df['thingcell'] = '8a7a94db75d3dca9017643bc4d2b379c'
donor4_calc = get_calculated_df(donor4_df)
donor4_calc.to_csv(PATH + 'donor4_thing97A-D.csv')

concatenated = pd.concat([donor1_df.assign(dataset='donor1'),
                          donor2_df.assign(dataset='donor2'),
                          donor3_df.assign(dataset='donor3'),
                          donor4_df.assign(dataset='donor4')])

concatenated.to_csv(PATH + 'combined_donors.csv')

fig = plt.gcf()

# Change seaborn plot size
fig.set_size_inches(24, 16)

sns.scatterplot(x='mean', y='CV', data=concatenated,
                hue='dataset', style='dataset')
fig.savefig(PATH + 'out_plot.png')
```
