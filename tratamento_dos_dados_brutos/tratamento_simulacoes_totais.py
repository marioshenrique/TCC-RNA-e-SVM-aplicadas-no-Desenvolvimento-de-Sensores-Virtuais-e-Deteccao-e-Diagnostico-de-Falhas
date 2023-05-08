import pandas as pd

# Ler o arquivo CSV e armazenar em um dataframe "df"
df = pd.read_csv("Simulacoes_12s_acuracia.csv")

# Alterar a ordem das colunas, movendo "t (min)" para o final do dataframe
cols = list(df.columns)
cols.remove("t (min)")
cols.append("t (min)")
df = df[cols]

# Criar a coluna "t (s)" com a conversão de minutos para segundos
df["t (s)"] = df["t (min)"] * 60

# Exportar o dataframe "df" como "Simulacoes_softsensor_tratada.csv"
df.to_csv("Simulacoes_12s_acuracia_tratada.csv", index=False)