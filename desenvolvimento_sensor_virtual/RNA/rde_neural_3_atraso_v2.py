import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense
from tensorflow.keras.callbacks import EarlyStopping

# 1. Carregar os dados de treinamento
dados = pd.read_csv("Simulacoes_9s_tratada_3_atraso.csv")
X_train = dados[['CB(k-3)','CB(k-2)','CB(k-1)','Tr(k-3)','Tr(k-2)','Tr(k-1)','Tc(k-3)','Tc(k-2)','Tc(k-1)','Tr0(k-3)','Tr0(k-2)','Tr0(k-1)','CA0(k-3)','CA0(k-2)','CA0(k-1)','Qc(k-3)','Qc(k-2)','Qc(k-1)','qr(k-3)','qr(k-2)','qr(k-1)']]
y_train = dados['CB(K)']

# 2. Carregar os dados de validação
dados_validacao = pd.read_csv("Simulacoes_9s_acuracia_tratada_3_atraso.csv")
X_validacao = dados_validacao[['CB(k-3)','CB(k-2)','CB(k-1)','Tr(k-3)','Tr(k-2)','Tr(k-1)','Tc(k-3)','Tc(k-2)','Tc(k-1)','Tr0(k-3)','Tr0(k-2)','Tr0(k-1)','CA0(k-3)','CA0(k-2)','CA0(k-1)','Qc(k-3)','Qc(k-2)','Qc(k-1)','qr(k-3)','qr(k-2)','qr(k-1)']]
y_validacao = dados_validacao['CB(K)']

# 3. Normalizar os dados
scaler = MinMaxScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_validacao_scaled = scaler.transform(X_validacao)

# 4. Criar o modelo de rede neural
model = Sequential()
model.add(Dense(16, activation='relu', input_shape=(X_train_scaled.shape[1],)))
model.add(Dense(8, activation='relu'))
model.add(Dense(1))

model.compile(optimizer='adam', loss='mse')

# 5. Treinar o modelo com os dados de treinamento e validação
early_stop = EarlyStopping(monitor='val_loss', patience=10)
model.fit(X_train_scaled, y_train, epochs=100, batch_size=32, verbose=0, validation_data=(X_validacao_scaled, y_validacao), callbacks=[early_stop])

# 6. Salvar o modelo treinado
model.save('modelo_treinado_9s_3_atraso.h5')

# 7. Carregar o modelo treinado (quando necessário)
# model = load_model('modelo_treinado.h5')

# 8. Avaliar o desempenho do modelo no conjunto de validação
y_pred_validacao = model.predict(X_validacao_scaled)
mse_validacao = mean_squared_error(y_validacao, y_pred_validacao)
r2_validacao = r2_score(y_validacao, y_pred_validacao)
mae_validacao = mean_absolute_error(y_validacao, y_pred_validacao)

print(f'MSE (validação): {mse_validacao}')
print(f'R2 (validação): {r2_validacao}')
print(f'MAE (validação): {mae_validacao}')

# 9. Salvar os valores reais e estimativas em arquivos CSV
y_validacao.to_csv('Valores_reais_validacao.csv', index=False, header=True)
pd.DataFrame(y_pred_validacao, columns=['Estimativas']).to_csv('valores_estimativas_validacao.csv', index=False)

# 10. Gerar e salvar o gráfico de dispersão em alta qualidade
plt.scatter(y_validacao, y_pred_validacao, s=10, c='blue', alpha=0.5)
plt.plot([y_validacao.min(), y_validacao.max()], [y_validacao.min(), y_validacao.max()], 'orange')
plt.xlabel('Valores Reais')
plt.ylabel('Estimativas')
plt.title('RNA Tempo de Amostragem 9s e 3 Tempos de Atraso - Validação')
plt.savefig('grafico_dispersao_9s_3_atraso_validacao.png', dpi=300)
plt.show()