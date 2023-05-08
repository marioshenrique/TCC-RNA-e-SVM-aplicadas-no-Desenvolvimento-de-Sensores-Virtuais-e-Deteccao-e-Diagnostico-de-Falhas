import pandas as pd

# Lendo o arquivo e armazenando no dataframe "df"
#df = pd.read_csv("Simulacoes_12s_tratada.csv", nrows=4000000)
df = pd.read_csv("Simulacoes_12s_acuracia_tratada.csv")

# Inicializando o dataframe "base_dados"
colunas = ["CB(k-3)", "CB(k-2)", "CB(k-1)", "Tr(k-3)", "Tr(k-2)", "Tr(k-1)", "Tc(k-3)", "Tc(k-2)", "Tc(k-1)", "Tr0(k-3)", "Tr0(k-2)", "Tr0(k-1)", "CA0(k-3)", "CA0(k-2)", "CA0(k-1)", "Qc(k-3)", "Qc(k-2)", "Qc(k-1)", "qr(k-3)", "qr(k-2)", "qr(k-1)", "CB(K)"]
base_dados = pd.DataFrame(columns=colunas)

# Construindo o dataframe "base_dados"
linhas = []
for i in range(3, len(df)):
    print(len(df)-i)
    linha = {
        "CB(k-3)": df["CB (kmol/m3)"][i - 3],
        "CB(k-2)": df["CB (kmol/m3)"][i - 2],
        "CB(k-1)": df["CB (kmol/m3)"][i - 1],
        "Tr(k-3)": df["Tr (K)"][i - 3],
        "Tr(k-2)": df["Tr (K)"][i - 2],
        "Tr(k-1)": df["Tr (K)"][i - 1],
        "Tc(k-3)": df["Tc (K)"][i - 3],
        "Tc(k-2)": df["Tc (K)"][i - 2],
        "Tc(k-1)": df["Tc (K)"][i - 1],
        "Tr0(k-3)": df["Tr0 (K)"][i - 3],
        "Tr0(k-2)": df["Tr0 (K)"][i - 2],
        "Tr0(k-1)": df["Tr0 (K)"][i - 1],
        "CA0(k-3)": df["CA0 (kmol/m3)"][i - 3],
        "CA0(k-2)": df["CA0 (kmol/m3)"][i - 2],
        "CA0(k-1)": df["CA0 (kmol/m3)"][i - 1],
        "Qc(k-3)": df["Qc (kJ/min)"][i - 3],
        "Qc(k-2)": df["Qc (kJ/min)"][i - 2],
        "Qc(k-1)": df["Qc (kJ/min)"][i - 1],
        "qr(k-3)": df["qr (m3/min)"][i - 3],
        "qr(k-2)": df["qr (m3/min)"][i - 2],
        "qr(k-1)": df["qr (m3/min)"][i - 1],
        "CB(K)": df["CB (kmol/m3)"][i],
    }
    linhas.append(pd.DataFrame(linha, index=[0]))

base_dados = pd.concat(linhas, ignore_index=True)
base_dados.to_csv('Simulacoes_12s_acuracia_tratada_3_atraso.csv', index=False)