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

Ca_inicial = 2.194305281355689 #concentração do componente A no interior do reator(kmol/m³)
Cb_inicial = 1.0950409066440403 #concentração inicial de B no interior do reator (kmol/m³)
Tr_inicial = 391.3895074432161 #temperatura inicial no interior do reator (K)
Tc_inicial = 401.39053689010575 #temperatura inicial do fluido refrigerante na saída (K)
Ca0 = 5.1 #concentração de A na corrente de alimentação do reator (kmol/m³)
Tr0 = 387.05 #temperatura da corrente de alimentação de reagentes (K)
Qc = 144.5 #calor retirado ou fornecido pelo fluido refrigerante (kJ/min)
criterio_erro = 0.0001 
H = 0.001 #passo do método iterativo (min)
qr = 0.0032 #Vazão volumétrica da corrente de alimentação de reagentes (m³/min)

var_pert = 2 #variável a ser perturbada (0 = qr, 1 = Qc, 2 = Ca0, 3 = Tr0)
porc_pert = -2.58

arquivo_vs_t = 'perturbacao_Tr0_-258.csv'

"Simulação do comportamento do sistema após uma pertubação no sistema ideal de operação"

reator1.Simu_perturbacoes_unitario(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,H,qr,arquivo_vs_t,var_pert,porc_pert)