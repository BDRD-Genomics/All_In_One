import pandas as pd 
import argparse


parser = argparse.ArgumentParser(
                    prog='format ICTV viral family genome sizes',
					usage='python3 merge_viralFamily_genomeSizes.py -n newFile -o oldFile',
					epilog='')
# add options
parser.add_argument("-n", "--newfile")
parser.add_argument("-o", "--oldfile")

args = parser.parse_args()
#print(args.input)
if (args.newfile == None or args.oldfile == None):
        print(parser.usage)
        exit(0)


primary_df = pd.read_csv(args.newfile, sep="\t")
old_df = pd.read_csv(args.oldfile, sep="\t")

print(primary_df.columns)
print(old_df.columns)


merged_df = pd.concat([primary_df, old_df])
merged_df = merged_df.drop_duplicates(subset=["Family"], keep="first")
merged_df.to_csv("viralFamily_genomeSize.txt", sep="\t", index=False)