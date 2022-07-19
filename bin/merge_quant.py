import pandas as pd
import glob

mfiles = glob.glob("*quant.sf")
mprefix = [e.split('-')[0] for e in mfiles]

final_df = pd.read_csv(mfiles[0], sep="\t")
final_df = final_df.add_suffix(mprefix[0])

final_df.rename(columns = {'Name'+mprefix[0]:'Name'}, inplace = True)

for i in range(1, len(mfiles)):
    mfile = mfiles[i]
    df_temp = pd.read_csv(mfile, sep="\t")
    df_temp = df_temp.add_suffix(mprefix[i])
    df_temp.rename(columns = {'Name'+mprefix[i]:'Name'}, inplace = True)
    final_df = final_df.merge(df_temp, on='Name', how='left')

final_df["sum_counts"] = final_df.filter(regex='NumReads*').sum(1)
final_df.to_csv("merged_quant.sf")
