import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn import svm
from sklearn.model_selection import train_test_split

# 1. Carregar os dados
dados = pd.read_csv("Simulacoes_3s_tratada_1_atraso.csv")
X = dados[['CB(k-1)', 'Tr(k-1)', 'Tc(k-1)', 'Tr0(k-1)', 'CA0(k-1)', 'Qc(k-1)', 'qr(k-1)']]
y = dados['CB(K)']

# 2. Dividir os dados em treinamento e teste
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 3. Normalizar os dados
scaler = MinMaxScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# 4. Criar o modelo de Máquinas de Vetores de Suporte
model = svm.SVR(kernel='rbf', C=100, gamma='scale')

# 5. Treinar o modelo com os dados de treinamento
model.fit(X_train_scaled, y_train)

# 8. Avaliar o desempenho do modelo no conjunto de teste
y_pred_test = model.predict(X_test_scaled)
mse_test = mean_squared_error(y_test, y_pred_test)
r2_test = r2_score(y_test, y_pred_test)
mae_test = mean_absolute_error(y_test, y_pred_test)

print(f'MSE (teste): {mse_test}')
print(f'R2 (teste): {r2_test}')
print(f'MAE (teste): {mae_test}')

# 9. Salvar os valores reais e estimativas em arquivos CSV
y_test.to_csv('Valores_reais_teste.csv', index=False, header=True)
pd.DataFrame(y_pred_test, columns=['Estimativas']).to_csv('valores_estimativas_teste.csv', index=False)

# 10. Gerar e salvar o gráfico de dispersão em alta qualidade
plt.scatter(y_test, y_pred_test, s=10, c='blue', alpha=0.5)
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'orange')
plt.xlabel('Valores Reais')
plt.ylabel('Estimativas')
plt.title('SVM Tempo de Amostragem 3s e 1 Tempo de Atraso')
plt.savefig('grafico_dispersao_3s_1_atraso_SVM.png', dpi=300)
plt.show()