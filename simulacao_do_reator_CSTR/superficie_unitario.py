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

"Abaixo temos a simulação da superfície dentro de um intervalo especificado de qr"
qr = 0.0029 #limite superior de qr

while (qr >= 0.0010): #critério com o limite inferior de qr
    
    Ca_inicial = 5.1 #concentração do componente A no interior do reator(kmol/m³)
    Cb_inicial = 0 #concentração inicial de B no interior do reator (kmol/m³)
    Tr_inicial = 370 #temperatura inicial no interior do reator (K)
    Tc_inicial = 350 #temperatura inicial do fluido refrigerante na saída (K)
    Ca0 = 5.1 #concentração de A na corrente de alimentação do reator (kmol/m³)
    Tr0 = 387.05 #temperatura da corrente de alimentação de reagentes (K)
    criterio_erro = 0.000001 
    H = 0.001 #passo do método iterativo (min)
    #qr = 0.0036  #vazão volumétrica da corrente de alimentação de reagentes (m³/min)
    liminf_Qc = -500 #limite inferior do range de energia térmica removida do fluido refrigerante (kJ/min)
    limsup_Qc = 500#limite superior do range de energia térmica removida do fluido refrigerante (kJ/min)
    H_Qc = 0.5 #passo de varredura do intervalo de energia térmica removida do fluido refrigerante (kJ/min)
    arquivo = 'superficie_qr_'+str(qr)+'_Qcinf_'+str(liminf_Qc)+'_Qcsup_'+str(limsup_Qc)+'.csv'
    reator1.Simu_superficieCB_unitario(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,criterio_erro,H,qr,liminf_Qc,limsup_Qc,H_Qc,arquivo)
    
    qr = qr - 0.0001
    qr = round(qr,4)
