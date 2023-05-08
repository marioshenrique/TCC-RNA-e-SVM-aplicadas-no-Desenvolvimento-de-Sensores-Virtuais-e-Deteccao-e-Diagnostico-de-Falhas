import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn import svm
from sklearn.model_selection import train_test_split, GridSearchCV

# 1. Carregar os dados
dados = pd.read_csv("Simulacoes_12s_tratada_3_atraso.csv")
X = dados[['CB(k-3)','CB(k-2)','CB(k-1)','Tr(k-3)','Tr(k-2)','Tr(k-1)','Tc(k-3)','Tc(k-2)','Tc(k-1)','Tr0(k-3)','Tr0(k-2)','Tr0(k-1)','CA0(k-3)','CA0(k-2)','CA0(k-1)','Qc(k-3)','Qc(k-2)','Qc(k-1)','qr(k-3)','qr(k-2)','qr(k-1)','CB(K)']]
y = dados['CB(K)']

# 2. Dividir os dados em treinamento e teste
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 3. Normalizar os dados
scaler = MinMaxScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# 4. Criar o modelo de Máquinas de Vetores de Suporte
model = svm.SVR()

# 5. Definir os hiperparâmetros a serem testados
#param_grid = {'kernel': ['linear', 'rbf', 'poly'],'C': [0.5, 1, 10, 100],'gamma': ['scale', 'auto', 0.1, 1, 10]}

param_grid = {'kernel': ['rbf'],'C': [0.5],'gamma': [1]}

# 6. Aplicar a técnica de validação cruzada
grid_search = GridSearchCV(model, param_grid, cv=5, scoring='neg_mean_squared_error', verbose=1)
grid_search.fit(X_train_scaled, y_train)

# 7. Imprimir os melhores hiperparâmetros encontrados
print("Melhores hiperparâmetros encontrados:")
print(grid_search.best_params_)

# 8. Avaliar o desempenho do modelo com os melhores hiperparâmetros no conjunto de teste
best_model = grid_search.best_estimator_
y_pred_test = best_model.predict(X_test_scaled)
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
plt.title('SVM Tempo de Amostragem 12s e 3 Tempo de Atraso')
plt.savefig('grafico_dispersao_12s_3_atraso.png', dpi=300)
plt.show()