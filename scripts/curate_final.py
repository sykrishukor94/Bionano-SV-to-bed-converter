# load packages
import sys
import pandas as pd
#from pandas import DataFrame, Series
import os 
import pathlib
working_dir = pathlib.Path().absolute()
os.chdir(working_dir)

# helper function to print string
# def print_string(string):
#     string += "processed"
#     return string.strip()

# color function which adds 'color' column to df
def color_sv_by_type(type_str):
    type_str = type_str.strip()
    if "del" in type_str:
        return "255,0,0"
    elif "dup" in type_str:
         # ignore duplications to stay in concordance with ctrl_df['Type']
        return "0,0,255"
    elif "ins" in type_str:
        return "0,0,255"
    elif "inv" in type_str:
        return "0,255,255"
    elif "trans" in type_str:
        return "255,0,255"
    elif "mono" in type_str:
        return "255,0,0"
    elif "tri" in type_str:
        return "0,0,255"
    elif "rear" in type_str:
        return "128,0,128"
    else:
        return(type_str)

# main function
def main(input_filepath, output_filepath):

    print("Argument List: ", str(sys.argv))

    # how to call the function: python "executable.py filepath" "input filepath" "output filepath"
    # example: python "C:\Users\mmiesen\OneDrive - Bionano Genomics\Desktop\Bionano\myCuration\executable.py" "\C:\Users\mmiesen\OneDrive - Bionano Genomics\Desktop\Bionano\myCuration\Postnatal_Raw_all_cluster_variant.txt" "C:\Users\mmiesen\OneDrive - Bionano Genomics\Desktop\Bionano\myCuration\output\postnatal_sv_knowledge_base.bed"

    # LOAD DATA: filepaths and reading csv as a dataframe
    df = pd.read_table(input_filepath, skiprows=1) #\t for tab separated , for comma separated
    
    print("Originial imported df columns: ", df.columns)   

    # declaring column types as strings and lower case
    df.columns = df.columns.str.lower()
    df = df.astype({'size':str,'type':str})

    # call color column function
    df['color'] = df.apply(lambda x: color_sv_by_type(x['type']), axis=1)
    # making name column
    df['name'] = df['size'].astype(str) + "_" + df["replicate"] + "_" + df["type"] + "_" + df["classification"] #TO ADD: caseid, sex, clinical indication, (size without decimal), orientation of translocation, VAF, ethnicity

    #other columns
    df['score'] = ""
    df['strand'] = ""
    df['blockCount'] = ""
    df['blockSizes'] = ""
    df['type'].unique()
    
    #sex chromosomes: chr23 is chrX and chr24 is chrY
    df['chrom1'] = df['chrom1'].replace(23, "X")
    df['chrom2'] = df['chrom2'].replace(23, "X")
    df['chrom1'] = df['chrom1'].replace(24, "Y")
    df['chrom2'] = df['chrom2'].replace(24, "Y")

    # subsets inverse of inter-translocations
    # non_trans_df = df.loc[~df['type'].str.contains('interchr')] 
    non_trans_df = df[~df['type'].isin(['trans_intrachr_common','trans_interchr_common','translocation_interchr','translocation_intrachr','trans_interchr_segdupe','trans_intrachr_segdupe'])]
    # subsets inter-translocations
    # interchr_trans_df = df.loc[df['type'].str.contains('interchr')]
    trans_df = df[df['type'].isin(['trans_intrachr_common','trans_interchr_common','translocation_interchr','translocation_intrachr','trans_interchr_segdupe','trans_intrachr_segdupe'])]


    # making 1st bed data frame (non-translocations)
    if non_trans_df.shape[0] > 0:
        bed_df1 = non_trans_df[['chrom1', 'position1', 'position2', 'name', 'score', 'strand', 'blockCount', 'blockSizes', 'color']]
        bed_df1["chrom1"] = "chr" + bed_df1["chrom1"].astype(str)
        bed_df1 = bed_df1.fillna(0)
        bed_df1 = bed_df1.astype({'position1':'int64'})
        bed_df1 = bed_df1.astype({'position2':'int64'})
        # bed_df1.to_csv(output_filepath, index=False)
        print("Shape of non-translocations (bed_df1): ", bed_df1.shape)
    else:
        print("non_trans_df is empty")
     
    # making 2nd bed data frame (inter-chromosomal translocations) where one row is split into two
        # original row: position1 = x, position2 = y
        # row1: position1 = x-6667, position2 = x
        # row2: position1 = y-6667, position2 = y
    
    if trans_df.shape[0] > 0:
        #split the original row containing "interchr" into two separate rows
        new_rows = []
        for _, row_to_split in trans_df.iterrows():
            # first row that's split
            new_row_1 = row_to_split.copy()
            new_row_1['position1'] = row_to_split['position1'] - 6667
            new_row_1['position2'] = row_to_split['position1']
            new_row_1['chrom1'] = row_to_split['chrom1']
            # second row that's split
            new_row_2 = row_to_split.copy()
            new_row_2['position1'] = row_to_split['position2'] - 6667
            new_row_2['position2'] = row_to_split['position2']
            new_row_2['chrom1'] = row_to_split['chrom2']
            # putting together
            new_rows.append(new_row_1.to_frame().T)
            new_rows.append(new_row_2.to_frame().T)
        bed_df2 = pd.concat(new_rows)
        print("bed_df2 first 5 rows: ", bed_df2.head())
        bed_df2 = bed_df2[['chrom1', 'position1', 'position2', 'name', 'score', 'strand', 'blockCount', 'blockSizes', 'color']]
        # bed_df2 = interchr_trans_df[['chrom1', 'position1', 'position2', 'replicate', 'name', 'color']]
        bed_df2["chrom1"] = "chr" + bed_df2["chrom1"].astype(str)
        bed_df2 = bed_df2.fillna(0)
        bed_df2 = bed_df2.astype({'position1':'int64'})
        bed_df2 = bed_df2.astype({'position2':'int64'})
        print("Shape of interchr & intrachr translocations (bed_df2): ", bed_df2.shape)

        # merge non-translocations together
        final_df = pd.concat([bed_df1, bed_df2])
        final_df = final_df.rename(columns={'chrom1': 'chrom', 'position1': 'chromStart', 'position2': 'chromEnd'})
        final_df.to_csv(output_filepath, sep='\t', index=False)
    else:
        print("interchr_trans_df is empty")


# call main function
if __name__ == "__main__":
    
    # print("Enter the input file path: ")
    # input_filepath = input().strip()

    # print("Enter the output file path: ")
    # output_filepath = input().strip()
    
    # alert
    if len(sys.argv) != 3:
        print("Use valid commands as such: python filepathfor_executable.py input_filepath output_filepath")
        sys.exit(1)

    input_filepath = sys.argv[1]
    output_filepath = sys.argv[2] 

    main(input_filepath, output_filepath)