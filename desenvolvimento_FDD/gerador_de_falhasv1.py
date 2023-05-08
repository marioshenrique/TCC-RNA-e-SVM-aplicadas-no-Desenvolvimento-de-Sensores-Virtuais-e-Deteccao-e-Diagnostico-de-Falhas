import pandas as pd
import numpy as np

# Ler o arquivo e armazenar em um dataframe chamado "df"
df = pd.read_csv("simulacoes_9s.csv")

# Criar um novo dataframe chamado "Simulacoes_9s_falhas"
Simulacoes_9s_falhas = pd.DataFrame(columns=['CA (kmol/m3)', 'Tr (K)', 'Tc (K)', 'Tr0 (K)', 'CA0 (kmol/m3)', 'Qc (kJ/min)', 'qr (m3/min)', 'Classes_Falhas'])

num_normal = 0
num_falha1 = 0
num_falha2 = 0
num_falha3 = 0
num_falha4 = 0
num_falha5 = 0
num_falha6 = 0

# Função para realizar a perturbação
def perturbar(valor):
    return valor * (1 + np.random.uniform(-0.9, 0.9))

for index, row in df.iterrows():
    row = pd.Series(row)
    # Gerar um número aleatório (0 ou 1) para decidir se a linha sofrerá perturbação
    perturbacao = np.random.randint(0, 2)

    # Se não houver perturbação
    if perturbacao == 0:
        row['Classes_Falhas'] = 'normal'
        num_normal = num_normal + 1
        Simulacoes_9s_falhas = pd.concat([Simulacoes_9s_falhas, row.to_frame().T], ignore_index=True)
    else:
        # Escolher uma coluna aleatoriamente para perturbar
        coluna_perturbada = np.random.choice(['CA0 (kmol/m3)', 'Tr0 (K)', 'Qc (kJ/min)', 'qr (m3/min)', 'Tr (K)', 'Tc (K)'])
        row[coluna_perturbada] = perturbar(row[coluna_perturbada])

        # Atribuir a classe de falha correspondente
        if coluna_perturbada == 'CA0 (kmol/m3)':
            row['Classes_Falhas'] = 'Falha_1'
            num_falha1 = num_falha1 + 1
        elif coluna_perturbada == 'Tr0 (K)':
            row['Classes_Falhas'] = 'Falha_2'
            num_falha2 = num_falha2 + 1
        elif coluna_perturbada == 'Qc (kJ/min)':
            row['Classes_Falhas'] = 'Falha_3'
            num_falha3 = num_falha3 + 1
        elif coluna_perturbada == 'qr (m3/min)':
            row['Classes_Falhas'] = 'Falha_4'
            num_falha4 = num_falha4 + 1
        elif coluna_perturbada == 'Tr (K)':
            row['Classes_Falhas'] = 'Falha_5'
            num_falha5 = num_falha5 + 1
        elif coluna_perturbada == 'Tc (K)':
            row['Classes_Falhas'] = 'Falha_6'
            num_falha6 = num_falha6 + 1
        
        print('Dados normais: '+str(num_normal))
        print('Dados Falha 1: '+str(num_falha1))
        print('Dados Falha 2: '+str(num_falha2))
        print('Dados Falha 3: '+str(num_falha3))
        print('Dados Falha 4: '+str(num_falha4))
        print('Dados Falha 5: '+str(num_falha5))
        print('Dados Falha 6: '+str(num_falha6))
    
        Simulacoes_9s_falhas = pd.concat([Simulacoes_9s_falhas, row.to_frame().T], ignore_index=True)
        
Simulacoes_9s_falhas.to_csv("Simulacoes_9s_falhas.csv", index=False)

print('Dados normais: '+str(num_normal))
print('Dados Falha 1: '+str(num_falha1))
print('Dados Falha 2: '+str(num_falha2))
print('Dados Falha 3: '+str(num_falha3))
print('Dados Falha 4: '+str(num_falha4))
print('Dados Falha 5: '+str(num_falha5))
print('Dados Falha 6: '+str(num_falha6))