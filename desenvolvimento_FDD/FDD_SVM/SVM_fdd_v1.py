import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, jaccard_score, confusion_matrix, classification_report
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

# Configurando a qualidade das figuras
plt.rcParams['figure.dpi'] = 300

# Lendo o arquivo e armazenando as colunas desejadas em um dataframe
file_path = "Simulacoes_9s_falhas.csv"
columns = ['Tr (K)', 'Tc (K)', 'Tr0 (K)', 'CA0 (kmol/m3)', 'Qc (kJ/min)', 'qr (m3/min)', 'Classes_Falhas']
df = pd.read_csv(file_path, usecols=columns)

# Reduzindo o tamanho do conjunto de dados
df = df.sample(frac=0.3, random_state=42)  # Reduz o conjunto de dados para 30% do tamanho original

# Dividindo os dados em atributos previsores e classes
X = df.iloc[:, :-1].values
y = df.iloc[:, -1].values

# Dividindo os dados em conjuntos de treino e teste
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Normalizando os dados
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# Treinando o classificador LinearSVC
classifier = SVC(random_state=42)
classifier.fit(X_train, y_train)

# Realizando as previsões
y_pred = classifier.predict(X_test)

# Calculando as métricas de desempenho
accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred, average=None)
recall = recall_score(y_test, y_pred, average='weighted')
f1 = f1_score(y_test, y_pred, average='weighted')
jaccard = jaccard_score(y_test, y_pred, average='weighted')
error_rate = 1 - accuracy
conf_matrix = confusion_matrix(y_test, y_pred)

# Exibindo a matriz de confusão como uma figura
class_labels = classifier.classes_
plt.figure(figsize=(10, 7))
sns.heatmap(conf_matrix, annot=True, cmap="YlGnBu", fmt='g', cbar=False, xticklabels=class_labels, yticklabels=class_labels)
plt.xlabel('Classe prevista')
plt.ylabel('Classe verdadeira')
plt.title('Matriz de Confusão')
plt.show()

# Exibindo as métricas de desempenho como uma figura
metrics_data = {'Métrica': ['Acurácia', 'Recall', 'F1-Score', 'Índice de Jaccard', 'Taxa de Erro'],
                'Valor': [accuracy, recall, f1, jaccard, error_rate]}

metrics_df = pd.DataFrame(metrics_data)

plt.figure(figsize=(10, 6))
sns.barplot(x='Métrica', y='Valor', data=metrics_df)
plt.title('Métricas de Desempenho')
plt.ylim(0, 1)
plt.show()

print("Acurácia: {:.2f}".format(accuracy))
print("Recall: {:.2f}".format(recall))
print("F1-Score: {:.2f}".format(f1))
print("Índice de Jaccard: {:.2f}".format(jaccard))
print("Taxa de Erro: {:.2f}".format(error_rate))

# Imprimindo a precisão por classe
print("\nPrecisão por classe:")
for i, label in enumerate(class_labels):
    print("Classe {}: {:.2f}".format(label, precision[i]))

# Imprimindo o relatório de classificação completo
print("\nRelatório de classificação:")
print(classification_report(y_test, y_pred))