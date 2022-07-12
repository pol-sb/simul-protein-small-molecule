# importing pandas
import pandas as pd

# reading pickle
df = pd.read_pickle("simulations_df.pkl")

# set the max columns to none
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)
pd.set_option("display.width", None)

# conc_query = df.loc[((df['conc'] == 20.0) & (df['temp'] == 320.0))]


print(df)

# print(
#     conc_query[
#         [
#             "protein",
#             "small_molec",
#             "conc",
#             "lambd",
#             "sigma",
#             "temp",
#             "idp_plat_avg",
#             "drg_plat_avg",
#         ]
#     ]
# )
