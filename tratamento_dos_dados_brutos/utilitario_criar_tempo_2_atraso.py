import pandas as pd

# Ler o arquivo e armazenar no dataframe "df"
df = pd.read_csv('tempo_amostragem_20s_perfiltemporal.csv')

# Criar um novo dataframe "base_dados"
base_dados = pd.DataFrame(columns=['CB(k-2)', 'CB(k-1)', 'Tr(k-2)', 'Tr(k-1)', 'Tc(k-2)', 'Tc(k-1)', 'Tr0(k-2)', 'Tr0(k-1)', 'CA0(k-2)', 'CA0(k-1)', 'Qc(k-2)', 'Qc(k-1)', 'qr(k-2)', 'qr(k-1)', 'CB(K)'])

# Construir o dataframe "base_dados" seguindo a lógica especificada
for i in range(2, len(df)):
    base_dados.loc[i-2] = [
        df.loc[i-2, 'CB (kmol/m3)'],
        df.loc[i-1, 'CB (kmol/m3)'],
        df.loc[i-2, 'Tr (K)'],
        df.loc[i-1, 'Tr (K)'],
        df.loc[i-2, 'Tc (K)'],
        df.loc[i-1, 'Tc (K)'],
        df.loc[i-2, 'Tr0 (K)'],
        df.loc[i-1, 'Tr0 (K)'],
        df.loc[i-2, 'CA0 (kmol/m3)'],
        df.loc[i-1, 'CA0 (kmol/m3)'],
        df.loc[i-2, 'Qc (kJ/min)'],
        df.loc[i-1, 'Qc (kJ/min)'],
        df.loc[i-2, 'qr (m3/min)'],
        df.loc[i-1, 'qr (m3/min)'],
        df.loc[i, 'CB (kmol/m3)']
    ]

# Exportar o dataframe "base_dados" como "dados_20s_2_atraso.csv"
base_dados.to_csv('dados_20s_2_atraso_perfiltemporal.csv', index=False)