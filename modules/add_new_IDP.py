import pandas as pd


SEQ_Q5_8_20 = (
    "MSKGPGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQ"
    "PYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPY"
    "QGRGDQPYQGRGDQPYQGY"
)

SEQ_Q5_8_24 = (
    "MSKGPGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQ"
    "PYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPY"
    "QGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGY"
)

SEQ_Q5_8_30 = (
    "MSKGPGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQ"
    "PYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPY"
    "QGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQGRGDQPYQG"
    "RGDQPYQGRGDQPYQGRGDQPYQGY"
)

prot_fil = pd.read_pickle("proteins_old.pkl")
print(prot_fil)

prot_fil.loc[len(prot_fil.index)] = [0.2, 6.5, 0.1, list(SEQ_Q5_8_20)]
prot_fil.loc[len(prot_fil.index)] = [0.2, 6.5, 0.1, list(SEQ_Q5_8_24)]
prot_fil.loc[len(prot_fil.index)] = [0.2, 6.5, 0.1, list(SEQ_Q5_8_30)]

# Be careful, hardcoded values in this section!!
prot_fil = prot_fil.rename(index={36: "Q5-8_20"})
prot_fil = prot_fil.rename(index={37: "Q5-8_24"})
prot_fil = prot_fil.rename(index={38: "Q5-8_30"})


prot_fil.to_pickle("./proteins.pkl")

prot_new = pd.read_pickle("proteins.pkl")
print(prot_new)
