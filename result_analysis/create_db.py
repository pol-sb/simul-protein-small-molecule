import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from pathlib import Path

file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(root))
import modules.utils as ut

np.set_printoptions(precision=3, threshold=1)
pd.set_option("display.max_colwidth", 15)
pd.set_option("display.max_columns", 15)


df = ut.generate_db()

print(df[df.temp == 323.0][df.lambd == 0.350][df.sigma == 0.450])
