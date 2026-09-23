import pandas as pd 
import numpy as np
import re
import argparse


parser = argparse.ArgumentParser(
                    prog='format ICTV viral family genome sizes',
					usage='python3 generate_family_sizes.py -i ICTV_Virus_Properties.csv',
					epilog='')
# add options
parser.add_argument("-i", "--input")

args = parser.parse_args()
#print(args.input)
if (args.input == None):
        print(parser.usage)
        exit(0)


virus_properties = args.input

def convert_range_to_avg(char_var):
    if char_var == "uncertain":
        return "NA"
    #print(char_var)
    char_var=char_var.replace("(+)", "").replace("(-)","")
    char_list=re.split(r'[,-]+', char_var)
    #print(char_list)
    #char_list=char_var.split("-")
    num_list = list(map(float, char_list))
    #print(num_list)
    avg=np.mean(num_list)*1000
    #print(avg)
    return avg



ictv_df = pd.read_csv(virus_properties)

#print(ictv_df.columns)
ictv_df["Genome Size"] = ictv_df["Genome size (kb/kbp)"].apply(convert_range_to_avg)


sub_ictv_df = ictv_df[["Family","Genome Size"]]
#print(sub_ictv_df.head())
sub_ictv_df.to_csv("new_viralFamily_genomeSize.txt", sep="\t", index=False)
print("output written to: new_viralFamily_genomeSize.txt")
