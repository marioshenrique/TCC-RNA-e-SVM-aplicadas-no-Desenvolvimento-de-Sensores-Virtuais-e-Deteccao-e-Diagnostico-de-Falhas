import pandas as pd

# Lendo o arquivo e armazenando no dataframe "df"
df = pd.read_csv("Simulacoes_3s_acuracia_tratada.csv")

# Inicializando o dataframe "base_dados"
colunas = ["CB(k-1)", "Tr(k-1)", "Tc(k-1)", "Tr0(k-1)", "CA0(k-1)", "Qc(k-1)", "qr(k-1)", "CB(K)"]
base_dados = pd.DataFrame(columns=colunas)

# Construindo o dataframe "base_dados"
linhas = []
for i in range(1, len(df)):
    print(len(df) - i)
    linha = {
        "CB(k-1)": df["CB (kmol/m3)"][i - 1],
        "Tr(k-1)": df["Tr (K)"][i - 1],
        "Tc(k-1)": df["Tc (K)"][i - 1],
        "Tr0(k-1)": df["Tr0 (K)"][i - 1],
        "CA0(k-1)": df["CA0 (kmol/m3)"][i - 1],
        "Qc(k-1)": df["Qc (kJ/min)"][i - 1],
        "qr(k-1)": df["qr (m3/min)"][i - 1],
        "CB(K)": df["CB (kmol/m3)"][i],
    }
    linhas.append(pd.DataFrame(linha, index=[0]))

base_dados = pd.concat(linhas, ignore_index=True)
base_dados.to_csv('Simulacoes_3s_acuracia_tratada_1_atraso.csv', index=False)