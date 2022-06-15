# importing pandas
import pandas as pd

# reading pickle
df = pd.read_pickle("simulations_df.pkl")

# set the max columns to none
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)

print(df)