import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tensorflow.keras.models import load_model
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error

# Carregando o modelo treinado
model = load_model("modelo_treinado_9s_3_atraso.h5")

# Carregando os dados de treinamento
dados_treinamento = pd.read_csv("Simulacoes_9s_tratada_3_atraso.csv")

# Carregando os dados de validação
dados_validacao = pd.read_csv("Simulacoes_9s_acuracia_tratada_3_atraso.csv", nrows=9000)

# Preparando os dados para a predição
X = dados_validacao[['CB(k-3)', 'CB(k-2)', 'CB(k-1)', 'Tr(k-3)', 'Tr(k-2)', 'Tr(k-1)', 'Tc(k-3)', 'Tc(k-2)', 'Tc(k-1)', 'Tr0(k-3)', 'Tr0(k-2)', 'Tr0(k-1)', 'CA0(k-3)', 'CA0(k-2)', 'CA0(k-1)', 'Qc(k-3)', 'Qc(k-2)', 'Qc(k-1)', 'qr(k-3)', 'qr(k-2)', 'qr(k-1)']]

# Preparando os dados de treinamento
X_train = dados_treinamento[['CB(k-3)', 'CB(k-2)', 'CB(k-1)', 'Tr(k-3)', 'Tr(k-2)', 'Tr(k-1)', 'Tc(k-3)', 'Tc(k-2)', 'Tc(k-1)', 'Tr0(k-3)', 'Tr0(k-2)', 'Tr0(k-1)', 'CA0(k-3)', 'CA0(k-2)', 'CA0(k-1)', 'Qc(k-3)', 'Qc(k-2)', 'Qc(k-1)', 'qr(k-3)', 'qr(k-2)', 'qr(k-1)']]

# Ajustando e normalizando os dados de entrada com o MinMaxScaler
scaler = MinMaxScaler()
scaler.fit(X_train)
X_scaled = scaler.transform(X)

# Realizando a predição
y_pred = model.predict(X_scaled)
y_true = dados_validacao["CB(K)"]

# Criando a coluna "t(s)"
dados_validacao["t(s)"] = np.arange(27, 27 + 9 * len(dados_validacao), 9)

# Salvando os dados validação em um arquivo csv
dados_validacao.to_csv("simulacao_soft_real.csv", index=False)

# Plotando o gráfico
plt.figure(figsize=(12, 5), dpi=200)
plt.plot(dados_validacao["t(s)"], y_true, linewidth=5, label="Valor Real")
plt.plot(dados_validacao["t(s)"], y_pred, linestyle='none', marker='o', markersize=3, label="Valor Predito")
plt.xlabel("tempo (s)")
plt.ylabel("Concentração de B (kmol/m³)")
plt.xticks(np.arange(0, max(dados_validacao["t(s)"]) + 1, 5000))
plt.title('RNA Tempo de Amostragem 9s e 3 Tempos de Atraso - Validação')
plt.legend()
plt.savefig("grafico_soft_real.png", dpi=300)
plt.show()

mse = mean_squared_error(y_true, y_pred)
r2 = r2_score(y_true, y_pred)
mae = mean_absolute_error(y_true, y_pred)

print(f'MSE (validação): {mse}')
print(f'R2 (validação): {r2}')
print(f'MAE (validação): {mae}')