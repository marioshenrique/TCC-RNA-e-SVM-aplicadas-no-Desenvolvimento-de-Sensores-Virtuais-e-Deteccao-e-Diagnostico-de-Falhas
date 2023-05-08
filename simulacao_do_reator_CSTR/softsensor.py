import psycopg2

from reatorcstr import BalancoEnergiaCSTR

"Inicialização dos parâmetros fixos do reator CSTR"

Vr = 0.01 #volume do reator (m³)
ro_r = 934.2 #densidade da corrente de reagente (kg/m³)
cp_r = 3.01 #capacidade calorífica do reagente (kJ/kgK)
mc = 5 #massa de fluido refrigerante (kg)
cp_c = 2 #capacidade calorífica do fluido refrigerante (kJ/kgK)
Ar = 0.215 #área da superfície de transferência de calor (m²)
U = 67.2 #coeficiente global de transferência de calor (kJ/min m²K)
k01 = 2.145e+10 #fator pré-exponencial da reação 1 (1/min)
k02 = 2.145e+10 #fator pré-exponencial da reação 2 (1/min)
k03 = 1.5072e+8 #fator pré-exponencial da reação 3 (1/min kmol)
E1_R = 9758.3 #razão entre a energia de ativação da reação 1 e R (K)
E2_R = 9758.3 #razão entre a energia de ativação da reação 2 e R (K)
E3_R = 8560 #razão entre a energia de atiação da reação 3 e R (K)
h1 = -4200 #Entalpia da reação 1 (kJ/kmol)
h2 = 11000 #Entalpia da reação 1 (kJ/kmol)
h3 = 41850 #Entalpia da reação 1 (kJ/kmol)

reator1 = BalancoEnergiaCSTR(Vr,ro_r,cp_r,mc,cp_c,Ar,U,k01,k02,k03,E1_R,E2_R,E3_R,h1,h2,h3)

"Parâmetros para execução da simulação de start-up de um reator CSTR"

Ca_estacionario = 2.194305281355689 #concentração do componente A no interior do reator(kmol/m³)
Cb_estacionario = 1.0950409066440403 #concentração inicial de B no interior do reator (kmol/m³)
Tr_estacionario = 391.3895074432161 #temperatura inicial no interior do reator (K)
Tc_estacionario = 401.39053689010575 #temperatura inicial do fluido refrigerante na saída (K)
Ca0_estacionario = 5.1 #concentração de A na corrente de alimentação do reator (kmol/m³)
Tr0_estacionario = 387.05 #temperatura da corrente de alimentação de reagentes (K)
Qc_estacionario = 144.5 #calor retirado ou fornecido pelo fluido refrigerante (kJ/min)
criterio_erro = 0.0001 
H = 0.001 #passo do método iterativo (min)
qr_estacionario = 0.0032 #Vazão volumétrica da corrente de alimentação de reagentes (m³/min)

nome_tabela = 'simulacoes_totais' #Nome da tabela que deverá armazenar os dados no banco de dados SQL

nome_tabela_20s = 'simulacoes_20s'
nome_tabela_75s = 'simulacoes_75s'
nome_tabela_95s = 'simulacoes_95s'
nome_tabela_115s = 'simulacoes_115s'
nome_tabela = 'simulacoes_teste' #Nome da tabela que deverá armazenar os dados no banco de dados SQL

qntd_simulacoes = 5000

"Simulação do comportamento do sistema após uma pertubação no sistema ideal de operação"

reator1.gera_dados_softsensor(Ca_estacionario,Cb_estacionario,Tr_estacionario,Tc_estacionario,
                              Ca0_estacionario,Tr0_estacionario,Qc_estacionario,criterio_erro,H,qr_estacionario,
                              qntd_simulacoes,nome_tabela_20s,nome_tabela_75s,nome_tabela_95s,nome_tabela_115s,
                              nome_tabela)