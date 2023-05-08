# Importar a biblioteca pandas
import pandas as pd

"Esse código é responsável pela simulação dos tempos de amostragens definidos"

# Ler o arquivo "Simulacoes_softsensor_tratada.csv" e armazenar os dados em um DataFrame chamado "df"
df = pd.read_csv("Simulacoes_softsensor_tratada.csv")

# Criar uma nova coluna chamada "contador" com valores inteiros em ordem crescente, começando com 1
df["contador"] = range(1, len(df) + 1)

# Selecionar as linhas onde o valor da coluna "contador" é múltiplo de 334
tempo_amostragem_20s = df[df["contador"] % 334 == 0] #define o tempo de amostragem

# Exportar o DataFrame "tempo_amostragem_20s" para um arquivo csv chamado "tempo_amostragem_20s.csv"
tempo_amostragem_20s.to_csv("tempo_amostragem_20s_perfiltemporal.csv", index=False)