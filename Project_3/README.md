## Project 3 - Tool for laboratory researchers to easily merge and migrate inconsistent Excel data

#### *Languages and Tools: Part I*
* Python
  * pandas
  * NumPy
  * os
  * csv
  * re
  * glob
  * date
  * argparse
* Bash
  * Slurm
---------------
### Helper Functions
* Python
* NumPy
* pandas
* glob
* date
* os
```python3
def find(sample):
    '''Function to extract IDs that were entered into the spreadsheet inconsistently or improperly.'''
    thing= re.findall("e[A-Z]+[0-9]+", str(sample))

    for i in sample:
        if re.findall("e[A-Z]+[0-9]+_", str(i)):
            ID.append(i.split('_')[0])

        elif re.findall("redacted", str(i)):
            ID.append(i.split()[-1])

        elif re.findall("redacted", str(i)):
            ID.append(i.split()[-1])

        elif re.findall("e[A-Z]+[0-9]+ [0-9]+ [0-9]+", str(i)):
            ID.extend(('redacted'))

    thing= np.unique(ID)
    return ID

def excel_to_df(excel_path):
    '''Converts manually entered Excel sheet to a dictionary.'''
    thing_excel = pd.read_excel(excel_path, sheet_name='all') # read full excel file in
    thing_excel = thing_excel.sort_index(ascending=False) # reverse order so that the newest things are first
    thing_excel['exp'] = thing_excel['thing exp']
    thing_excel = thing_excel.drop('thing exp', axis = 1)
    return thing_excel

def glob_files(thing_path, file_ext, thing_files, fields):
    '''Function to identify all headers from the thing output files in the current directory.'''
    # path to thing file directory
    os.chdir(thing_path)
    # use glob to get all thing output files with a .txt extension
    for file in glob.glob(f'*{file_ext}'):
        # save file names with .txt extension in a list
        thing_files.append(file)
    # iterate over each file in the thing_files list
    for name in thing_files:
        try:
            # open each file
            with open(name, 'r+', newline='') as f_in:
                try:
                    # get headers from each file
                    reader = csv.reader(f_in, delimiter='\t', quotechar='"')
                    headers = next(reader)
                    # if the column header isn't already in the list of headers, add it - for getting all fields
                    for col in headers:
                        if col not in fields:
                            fields.append(col)
                # this exception occurs when an iterator's next() method signals that there are no further values
                except StopIteration:
                    break
                    raise Exception("Reached end of iterator. Continuing.")
        # pass over any files without the proper permissions
        except PermissionError:
            break
            raise Exception("Looks like you don't have access to at least one file in this directory. Double-check your permissions settings.")

def file_date(out_file):
    '''Function to add the date to the output of the merge.'''
    date_of_file = out_file + date.today().strftime('_%Y_%m_%d') + '.csv'
    return date_of_file

def write_excel(out_file, out_path, keep_df, extra_df):
    '''Function for writing the merged data to an Excel file.'''
    date_file = out_file + date.today().strftime('_%Y_%m_%d') + '.xlsx'
    writer = pd.ExcelWriter(out_path + '/' + date_file, engine='xlsxwriter')
    keep_df.to_excel(writer, sheet_name='thing_COLS')
    extra_df.to_excel(writer, sheet_name='EXTRA_COLS')
    writer.save()
```

### Main Tool
* Python
* pandas
* NumPy
* csv
* re
* argparse
```python3
#!/usr/bin/env python3

# import packages
from helpers import *
import argparse

### MAIN PROGRAM
def main():
    '''Program to merge a thing excel sheet with thing out files.'''
    parser = argparse.ArgumentParser(description="Program to merge lines from thing output files with the thing Excel sheet based on sequence and other thing and string sequence (enter '-h' after the script name to see all flags).")
    parser.add_argument("-e", "--extension", help="File extension of the desired files (likely .txt).", required=True, type=str)
    parser.add_argument("-mp", "--thing_path", help="Path to the root directory containing the output files from the thing assay (provide just the path).", required=True, type=str)
    parser.add_argument("-ep", "--excel_path", help="Path to the excel sheet (provide both the path and the file).", required=True, type=str)
    parser.add_argument("-op", "--out_path", help="Path to eventual merged file (provide just the path).", required=True, type=str)
    parser.add_argument("-of", "--out_file", help="Name the output merged file (do not include a file extension).", required=True, type=str)
    parser.add_argument("-c", "--col_path", help="Path to .txt file containing relevant column headers (probably bleepblorp.txt).", required=True, type=str)

    # arguments that contain defaults and aren't required (chance they may change in the future, likely due to the Excel file)
    parser.add_argument("-m_id", "--thing_id", help="Specify the bleep blorp column header from the thing files (likely 'thingID').", default="thingID",required=False, type=str)
    parser.add_argument("-m_seq", "--thing_seq", help="Specify the string sequence column header from the thing files.", default="bleep", required=False, type=str)
    parser.add_argument("-e_id", "--excel_id", help="Specify the bleep blorp column header from the Excel file.", default="bloop", required=False, type=str)
    parser.add_argument("-e_seq", "--excel_seq", help="Specify the string sequence column header from the Excel file.", default="blorp", required=False, type=str)
    parser.add_argument("-id_c", "--id_col", help="Number of the thing file column where sample IDs are located.", default=0, required=False, type=int)
    args = parser.parse_args()

    # convert the Excel sheet to a dictionary
    excel_df = excel_to_df(args.excel_path)
    # covert excel dataframe to dictionary
    excel_dict = excel_df.to_dict('records')
    # get list of unique bleep blorps
    exps_excel = list(excel_df['exp'].unique())

    # processed list of bleep blorps
    IDs = find_IDs(exps_excel)

    # join all bleep blorps for re search later
    pattern = '|'.join(IDs)
    # list for all fields from globbed files
    fields = []
    # list to hold thing file names
    thing_files = []

    # glob the files in the current directory and identify all unique fields
    glob_files(args.thing_path, args.extension, thing_files, fields)

    # add file path field
    fields.append('file_name')

    # add date to the output file
    out_with_date = file_date(args.out_file)

    # open matched file
      with open(args.out_path + '/' + out_with_date, 'a+', newline='') as f_matched:
          # initialize a writer for the matched rows between the excel sheet and thing files
          writer_matched = csv.DictWriter(f_matched, fieldnames=fields)
          # write headers from the fields list for both files. files without certian fields will be written out as empty strings
          writer_matched.writeheader()

          # iterate over each file in the globbed directory
          for name in thing_files:
              # open each file and search for IDs
              with open(name, 'r+', newline='') as f_in1:
                  # get just the ID columns
                  ID_col = list(zip(*csv.reader(f_in1, delimiter= '\t')))[args.id_col]
                  # sarch the ID column for IDs from the Excel sheet
                  if re.search(pattern, str(ID_col)):
                      # open all fields
                      with open(name, 'r+', newline='') as f_in2:
                          # convert each line to a dictionary
                          for line1 in csv.DictReader(f_in2, delimiter='\t'):
                              # if an ID is recognized in the line, iterate over it
                              if re.search(pattern, str(line1)):
                                  # add a field for the original thing file
                                  line1.update({'file_name':name})
                                  thing_file_seq = line1[args.thing_seq]
                                  # dict lookup thingID key in the globbed thing files
                                  thing_file_thingID = line1[args.thing_id]

                                  # for each line in the excel dictionary, each line of excel file is in records
                                  for line2 in excel_dict:
                                      # check if the seq and thing columnsin the excel sheet are in the seq and thingID columns from file
                                      if (line2[args.excel_seq] == thing_file_seq and line2[args.excel_id] in thing_file_thingID):
                                          # merge the two dictionaries
                                          line3 = {**line1, **line2}
                                          # find any field headers that aren't already in the field list
                                          for k in line3.keys():
                                              if k not in fields:
                                                  fields.append(k)
                                          # write line to file if a match is found
                                          writer_matched.writerow(line3)
                                      else:
                                          # pass if no match
                                          pass

                  else:
                      pass


    # get headers for the "keep" sheet of the Excel file
    with open(args.col_path, 'r') as f:
        header_list = [line.rstrip('\n') for line in f]
    header_list.append('file_name')

    # get all headers from the merge, add them to a list
    thing_headers = []
    with open(args.out_path + '/' + out_with_date, newline='') as f:
      reader = csv.reader(f)
      row1 = next(reader)
    for col in row1:
        thing_headers.append(col)

    # read in df
    thing_data = pd.read_csv(args.out_path + '/' + out_with_date, low_memory=False, names=thing_headers)
    thing_data = thing_data.drop(0)

    # add columns to the dataframe if they don't exist
    for col in header_list:
        if col in thing_data.columns:
            pass
        else: # add columns containing NaNs if they don't already exists in the file
            thing_data[col] = np.nan

    thing_copy = thing_data.copy()
    # subset of dataframe to only inclue the columns2parse
    keep = thing_copy[header_list]
    # remaining columns from dataframe that weren't included in the Columns2parse file
    extra = thing_copy.drop(header_list, axis=1)

    # write final Excel sheet
    write_excel(args.out_file, args.out_path, keep, extra)

if __name__ == "__main__":
    main()
##### END
```
### Run the Above Script on the Cluster with Slurm
* Bash
* Slurm
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

THING_PATH=/home/thing_dir
EXCEL_PATH=/home/input_dir
OUT_PATH=/home/out_path
COLUMNS=/home/columns.txt
OUT_FILE="OUT"

/usr/bin/time -v ./THING.py -e .txt -c $COLUMNS -mp $THING_PATH -ep $EXCEL_PATH -op $OUT_PATH -of $OUT_FILE
```
