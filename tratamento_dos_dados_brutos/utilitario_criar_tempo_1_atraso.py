import pandas as pd

"Código responsável por pegar os dados da simulação de soft sensor e gerar uma nova planilha com 1 tempo de atrasado"

# Lendo o arquivo "tempo_amostragem_20s_modificado.csv" e armazenando no dataframe "df"
df = pd.read_csv("Simulacoes_115s_acuracia_tratada.csv")

# Criando um novo dataframe "base_dados" com as colunas especificadas
base_dados = pd.DataFrame(columns=["CB(k-1)", "Tr(k-1)", "Tc(k-1)", "Tr0(k-1)", "CA0(k-1)", "Qc(k-1)", "qr(k-1)", "CB(K)"])

# Construindo o dataframe "base_dados" conforme especificado
for i in range(len(df) - 1):
    print(i)
    base_dados.loc[i] = [
        df["CB (kmol/m3)"][i],
        df["Tr (K)"][i],
        df["Tc (K)"][i],
        df["Tr0 (K)"][i],
        df["CA0 (kmol/m3)"][i],
        df["Qc (kJ/min)"][i],
        df["qr (m3/min)"][i],
        df["CB (kmol/m3)"][i + 1],
    ]

# Exportando o dataframe "base_dados" como "dados_20s_1_atraso.csv"
base_dados.to_csv("Simulacoes_115s_acuracia_tratada_1_atraso.csv",index = False)