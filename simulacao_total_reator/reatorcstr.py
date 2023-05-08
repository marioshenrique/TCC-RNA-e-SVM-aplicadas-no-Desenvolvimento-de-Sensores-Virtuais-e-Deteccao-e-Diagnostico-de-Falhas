import numpy as np
import matplotlib.pyplot as plt
import random as random
import pandas as pd

class BalancoEnergiaCSTR():
    
    def __init__(self,Vr,ro_r,cp_r,mc,cp_c,Ar,U,k01,k02,k03,E1_R,E2_R,E3_R,h1,h2,h3):
        """Inicialização dos parâmetros fixos do reator CSTR"""
        self.Vr = Vr #volume do reator (m³)
        self.ro_r = ro_r #densidade da corrente de reagente (kg/m³)
        self.cp_r = cp_r #capacidade calorífica do reagente (kJ/kgK)
        self.mc = mc #massa de fluido refrigerante (kg)
        self.cp_c = cp_c #capacidade calorífica do fluido refrigerante (kJ/kgK)
        self.Ar = Ar #área da superfície de transferência de calor (m²)
        self.U = U #coeficiente global de transferência de calor (kJ/min m²K)
        self.k01 = k01 #fator pré-exponencial da reação 1 (1/min)
        self.k02 = k02 #fator pré-exponencial da reação 2 (1/min)
        self.k03 = k03 #fator pré-exponencial da reação 3 (1/min kmol)
        self.E1_R = E1_R #razão entre a energia de atiação da reação 1 e R (K)
        self.E2_R = E2_R #razão entre a energia de atiação da reação 2 e R (K)
        self.E3_R = E3_R #razão entre a energia de atiação da reação 3 e R (K)
        self.h1 = h1 #Entalpia da reação 1 (kJ/kmol)
        self.h2 = h2 #Entalpia da reação 2 (kJ/kmol)
        self.h3 = h3 #Entalpia da reação 3 (kJ/kmol)
        
        """Vetores para plotagem dos gráficos"""
        #gráficos de perfis de concentração com o tempo
        self.Cavst_plot = []
        self.Cbvst_plot = []
        self.Trvst_plot = []
        self.Tcvst_plot = []
        self.tempo_plot = []
        
        "Vetores auxiliares (Simu_estacionarioCSTR)"
        #graficos de valores de estado estacionário
        
        #gráficos vs qr
        self.Cavsqr_max = []
        self.qrCa_max = []
        self.Cbvsqr_max = []
        self.qrCb_max = []
        self.Trvsqr_max = []
        self.qrTr_max = []
        self.Tcvsqr_max = []
        self.qrTc_max = []
        #armazena o valor de Qc nos pontos máximos de Cb
        self.Qc_Cbmax_vs_qr = []
        
        self.Cavsqr_plot = []
        self.qrCa_plot = []
        self.Cbvsqr_plot = []
        self.qrCb_plot = []
        self.Trvsqr_plot = []
        self.qrTr_plot = []
        self.Tcvsqr_plot = []
        self.qrTc_plot = []
        self.Qc_Cbmax_vs_qr_plot = []
        
        #gráficos vs Qc
        self.Cavsqc_plot = []
        self.Cbvsqc_plot = []
        self.Trvsqr_plot = []
        self.Tcvsqc_plot = []
        self.qc_plot = []
        
        "Vetores auxiliares (Simu_superficieCSTR)"
        #Cb em função de qr e Qc
        self.Cb_sup_plot = []
        self.Ca_sup_plot = []
        self.Tr_sup_plot = []
        self.Tc_sup_plot = []
        self.qr_sup_plot = []
        self.Qc_sup_plot = []
        
        "vetores auxiliares (Simu_temporalCSTR)"
        self.Cb_estac_i = []
        self.tempo_estac = []
        
        "variáveis auxiliares (__RK_gerar_dados)"
        self.Ca_t = 0
        self.Cb_t = 0
        self.Tr_t = 0
        self.Tc_t = 0
        
        "vetores auxiliares (__RK_tempo_amostr)"
        self.Cavst_tamostr = []
        self.Cbvst_tamostr = []
        self.Trvst_tamostr = []
        self.Tcvst_tamostr = []
        self.tempo_tamostr = []
        
    
    @staticmethod
    def __calc_constvel(Tr_n,a,b):
        #a = fator pré-exponencial
        #b = E/R
        k = a*np.exp(-b*((1)/(Tr_n)))
        return k
    
    def __tempo_estac_CB(self,erro):
        #erro = erro entre um determinado valor e o valor do estado estacionário final
        cb_estac = self.Cbvst_plot[-1]
        n = 0
        tam = len(self.Cbvst_plot)
        while (n < tam):
            cb_estac_i = self.Cbvst_plot[n]
            erro_i = (np.abs(cb_estac_i - cb_estac)/cb_estac)*100 #erro percentual entre valor estacionário final e valor de Cb no tempo t =  n
            if erro_i <= erro:
                self.Cb_estac_i.clear()
                self.tempo_estac.clear()
                self.Cb_estac_i.append(cb_estac_i)
                self.tempo_estac.append(self.tempo_plot[n])
                return self.tempo_plot[n]
            n = n + 1
    
    def __tempo_estac_CB_tamostr(self,erro):
        #erro = erro entre um determinado valor e o valor do estado estacionário final
        cb_estac = self.Cbvst_tamostr[-1]
        n = 0
        tam = len(self.Cbvst_tamostr)
        while (n < tam):
            cb_estac_i = self.Cbvst_tamostr[n]
            erro_i = (np.abs(cb_estac_i - cb_estac)/cb_estac)*100 #erro percentual entre valor estacionário final e valor de Cb no tempo t =  n
            if erro_i <= erro:
                self.Cb_estac_i.clear()
                self.tempo_estac.clear()
                self.Cb_estac_i.append(cb_estac_i)
                self.tempo_estac.append(self.tempo_tamostr[n])
                return self.tempo_tamostr[n]
            n = n + 1
         
    def __RK__estacionario(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr):
        """Processo iterativo para construção dos perfis de concentração com o tempo"""
        "Esse bloco é responsável pela geração de dados de estado estacionário. Os dados são armazenados"
        "após a satisfação de um critério de parada."
        
        #Ca0: concentração do componente A na corrente de entrada (kmol/m³)
        #Tr0: Temperatura da corrente de alimentação de reagentes (K)
        #Qc: Calor transferido do/para fluido refrigerante (kJ/min)
        #Ca_inicial: concentração inicial de A no interior do reator
        #Cb_inicial: concentração inicial de B no interior do reator
        #Tc_inicial: temperatura inicial do fluido refrigerante na saída
        #Tr_inicial: temperatura inicial no interior do reator
        #criterio_erro: critério de parada (estado estacionário)
        #H: passo do método iterativo (min)
        #qr: vazão volumétrica da alimentação de reagentes
        #Tr = temperatura do reator (K)
        
        t_inicial = 0
        t = 0
        self.Cavst_plot.clear()
        self.Cbvst_plot.clear()
        self.Tcvst_plot.clear()
        self.Trvst_plot.clear()
        self.tempo_plot.clear()
        self.Cavst_plot.append(Ca_inicial)
        self.Cbvst_plot.append(Cb_inicial)
        self.Tcvst_plot.append(Tc_inicial)
        self.Trvst_plot.append(Tr_inicial)
        self.tempo_plot.append(t_inicial)
        Ca_n = Ca_inicial
        Cb_n = Cb_inicial
        Tr_n = Tr_inicial
        Tc_n = Tc_inicial
        "Cálculo das velocidades específicas de reação"
        k1_n = self.__calc_constvel(Tr_n,a = self.k01,b = self.E1_R)
        k2_n = self.__calc_constvel(Tr_n,a = self.k02,b = self.E2_R)
        k3_n = self.__calc_constvel(Tr_n,a = self.k03,b = self.E3_R)
        "Cálculo do calor de reação"
        hr_n = (self.h1*k1_n*Ca_n) + (self.h2*k2_n*Cb_n) +(self.h3*k3_n*((Ca_n)**2))
        "Cálculo dos valores vetoriais usando RK de primeira ordem"
        Ca_n1 = Ca_n + H*((((qr)/(self.Vr))*(Ca0 - Ca_n))-(k1_n*Ca_n)-(k3_n*((Ca_n)**2)))
        Cb_n1 = Cb_n + H*(-((qr)/(self.Vr))*Cb_n+(k1_n*Ca_n)-(k2_n*Cb_n))
        Tr_n1 = Tr_n + H*((((qr)/(self.Vr))*(Tr0 - Tr_n))-((hr_n)/(self.ro_r*self.cp_r))+(((self.Ar*self.U)/(self.Vr*self.ro_r*self.cp_r))*(Tc_n - Tr_n)))
        Tc_n1 = Tc_n + H*(((1)/(self.mc*self.cp_c))*(Qc + (self.Ar*self.U*(Tr_n - Tc_n))))
        "construção dos vetores"
        t = t + H
        self.Cavst_plot.append(Ca_n1)
        self.Cbvst_plot.append(Cb_n1)
        self.Tcvst_plot.append(Tc_n1)
        self.Trvst_plot.append(Tr_n1)
        self.tempo_plot.append(t)
        "Cálculo dos erros"
        erro_ca = np.abs(Ca_n1 - Ca_n)
        erro_cb = np.abs(Cb_n1 - Cb_n)
        erro_tr = np.abs(Tr_n1 - Tr_n)
        erro_tc = np.abs(Tc_n1 - Tc_n)
        
        soma_erro = erro_ca + erro_cb + erro_tr + erro_tc
        
        #print(erro_total)
        "atualização dos parâmetros"
        Ca_n = Ca_n1
        Cb_n = Cb_n1
        Tr_n = Tr_n1
        Tc_n = Tc_n1
        
        while (soma_erro > criterio_erro):
            "Cálculo das velocidades específicas de reação"
            k1_n = self.__calc_constvel(Tr_n,a = self.k01,b = self.E1_R)
            k2_n = self.__calc_constvel(Tr_n,a = self.k02,b = self.E2_R)
            k3_n = self.__calc_constvel(Tr_n,a = self.k03,b = self.E3_R)
            "Cálculo do calor de reação"
            hr_n = (self.h1*k1_n*Ca_n) + (self.h2*k2_n*Cb_n) +(self.h3*k3_n*((Ca_n)**2))
            "Cálculo dos valores vetoriais usando RK de primeira ordem"
            Ca_n1 = Ca_n + H*((((qr)/(self.Vr))*(Ca0 - Ca_n))-(k1_n*Ca_n)-(k3_n*((Ca_n)**2)))
            Cb_n1 = Cb_n + H*(-((qr)/(self.Vr))*Cb_n+(k1_n*Ca_n)-(k2_n*Cb_n))
            Tr_n1 = Tr_n + H*((((qr)/(self.Vr))*(Tr0 - Tr_n))-((hr_n)/(self.ro_r*self.cp_r))+(((self.Ar*self.U)/(self.Vr*self.ro_r*self.cp_r))*(Tc_n - Tr_n)))
            Tc_n1 = Tc_n + H*(((1)/(self.mc*self.cp_c))*(Qc + (self.Ar*self.U*(Tr_n - Tc_n))))
            "construção dos vetores"
            t = t + H
            self.Cavst_plot.append(Ca_n1)
            self.Cbvst_plot.append(Cb_n1)
            self.Tcvst_plot.append(Tc_n1)
            self.Trvst_plot.append(Tr_n1)
            self.tempo_plot.append(t)
            "Cálculo dos erros"
            erro_ca = np.abs(Ca_n1 - Ca_n)
            erro_cb = np.abs(Cb_n1 - Cb_n)

            erro_tr = np.abs(Tr_n1 - Tr_n)

            erro_tc = np.abs(Tc_n1 - Tc_n)
            
            soma_erro = erro_ca + erro_cb + erro_tr + erro_tc
            "atualização dos parâmetros"
            Ca_n = Ca_n1
            Cb_n = Cb_n1
            Tr_n = Tr_n1
            Tc_n = Tc_n1
            
            
    def __RK_tempo_amostr(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,t_final,H,qr):
        """Processo iterativo para geração de dados de simulação"""
        
        #Ca0: concentração do componente A na corrente de entrada (kmol/m³)
        #Tr0: Temperatura da corrente de alimentação de reagentes (K)
        #Qc: Calor transferido do/para fluido refrigerante (kJ/min)
        #Ca_inicial: concentração inicial de A no interior do reator
        #Cb_inicial: concentração inicial de B no interior do reator
        #Tc_inicial: temperatura inicial do fluido refrigerante na saída
        #Tr_inicial: temperatura inicial no interior do reator
        #Ca0: concentração de A na corrente de alimentação do reator
        #Tr0: temperatura da corrente de alimentação de reagentes
        #Qc: calor retirado ou fornecido pelo fluido refrigerante
        #t_final: tempo final
        #H: passo do método iterativo (min)
        #qr: vazão volumétrica da alimentação de reagentes
        #Tr = temperatura do reator (K)
        
        t = 0
        t_inicial = 0
        self.Cavst_tamostr.clear()
        self.Cbvst_tamostr.clear()
        self.Tcvst_tamostr.clear()
        self.Trvst_tamostr.clear()
        self.tempo_tamostr.clear()
        self.Cavst_tamostr.append(Ca_inicial)
        self.Cbvst_tamostr.append(Cb_inicial)
        self.Tcvst_tamostr.append(Tc_inicial)
        self.Trvst_tamostr.append(Tr_inicial)
        self.tempo_tamostr.append(t_inicial)
        Ca_n = Ca_inicial
        Cb_n = Cb_inicial
        Tr_n = Tr_inicial
        Tc_n = Tc_inicial
        while (t < t_final):
            "Cálculo das velocidades específicas de reação"
            k1_n = self.__calc_constvel(Tr_n,a = self.k01,b = self.E1_R)
            k2_n = self.__calc_constvel(Tr_n,a = self.k02,b = self.E2_R)
            k3_n = self.__calc_constvel(Tr_n,a = self.k03,b = self.E3_R)
            "Cálculo do calor de reação"
            hr_n = (self.h1*k1_n*Ca_n) + (self.h2*k2_n*Cb_n) +(self.h3*k3_n*((Ca_n)**2))
            "Cálculo dos valores vetoriais usando RK de primeira ordem"
            Ca_n1 = Ca_n + H*((((qr)/(self.Vr))*(Ca0 - Ca_n))-(k1_n*Ca_n)-(k3_n*((Ca_n)**2)))
            Cb_n1 = Cb_n + H*(-((qr)/(self.Vr))*Cb_n+(k1_n*Ca_n)-(k2_n*Cb_n))
            Tr_n1 = Tr_n + H*((((qr)/(self.Vr))*(Tr0 - Tr_n))-((hr_n)/(self.ro_r*self.cp_r))+(((self.Ar*self.U)/(self.Vr*self.ro_r*self.cp_r))*(Tc_n - Tr_n)))
            Tc_n1 = Tc_n + H*(((1)/(self.mc*self.cp_c))*(Qc + (self.Ar*self.U*(Tr_n - Tc_n))))
            "construção dos vetores"
            t = t + H
            self.Cavst_tamostr.append(Ca_n1)
            self.Cbvst_tamostr.append(Cb_n1)
            self.Tcvst_tamostr.append(Tc_n1)
            self.Trvst_tamostr.append(Tr_n1)
            self.tempo_tamostr.append(t)
            "atualização dos parâmetros"
            Ca_n = Ca_n1
            Cb_n = Cb_n1
            Tr_n = Tr_n1
            Tc_n = Tc_n1
    
    def __RK_gerar_dados(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,t_final,H,qr):
        
        #Ca0: concentração do componente A na corrente de entrada (kmol/m³)
        #Tr0: Temperatura da corrente de alimentação de reagentes (K)
        #Qc: Calor transferido do/para fluido refrigerante (kJ/min)
        #Ca_inicial: concentração inicial de A no interior do reator
        #Cb_inicial: concentração inicial de B no interior do reator
        #Tc_inicial: temperatura inicial do fluido refrigerante na saída
        #Tr_inicial: temperatura inicial no interior do reator
        #Ca0: concentração de A na corrente de alimentação do reator
        #Tr0: temperatura da corrente de alimentação de reagentes
        #Qc: calor retirado ou fornecido pelo fluido refrigerante
        #t_final: tempo final
        #H: passo do método iterativo (min)
        #qr: vazão volumétrica da alimentação de reagentes (m³/min)
        #Tr = temperatura do reator (K)
        
        """Processo iterativo para geração de valores em determinados intervalos de tempo"""
    
        t = 0
        Ca_n = Ca_inicial
        Cb_n = Cb_inicial
        Tr_n = Tr_inicial
        Tc_n = Tc_inicial
        
        flag = True
        while flag:
            if (t < t_final):
                "Cálculo das velocidades específicas de reação"
                k1_n = self.__calc_constvel(Tr_n,a = self.k01,b = self.E1_R)
                k2_n = self.__calc_constvel(Tr_n,a = self.k02,b = self.E2_R)
                k3_n = self.__calc_constvel(Tr_n,a = self.k03,b = self.E3_R)
                "Cálculo do calor de reação"
                hr_n = (self.h1*k1_n*Ca_n) + (self.h2*k2_n*Cb_n) +(self.h3*k3_n*((Ca_n)**2))
                "Cálculo dos valores vetoriais usando RK de primeira ordem"
                Ca_n1 = Ca_n + H*((((qr)/(self.Vr))*(Ca0 - Ca_n))-(k1_n*Ca_n)-(k3_n*((Ca_n)**2)))
                Cb_n1 = Cb_n + H*(-((qr)/(self.Vr))*Cb_n+(k1_n*Ca_n)-(k2_n*Cb_n))
                Tr_n1 = Tr_n + H*((((qr)/(self.Vr))*(Tr0 - Tr_n))-((hr_n)/(self.ro_r*self.cp_r))+(((self.Ar*self.U)/(self.Vr*self.ro_r*self.cp_r))*(Tc_n - Tr_n)))
                Tc_n1 = Tc_n + H*(((1)/(self.mc*self.cp_c))*(Qc + (self.Ar*self.U*(Tr_n - Tc_n))))
                #atualização dos parâmetros
                Ca_n = Ca_n1
                Cb_n = Cb_n1
                Tr_n = Tr_n1
                Tc_n = Tc_n1
                t = t + H
            elif (t >= t_final):
                "Cálculo das velocidades específicas de reação"
                k1_n = self.__calc_constvel(Tr_n,a = self.k01,b = self.E1_R)
                k2_n = self.__calc_constvel(Tr_n,a = self.k02,b = self.E2_R)
                k3_n = self.__calc_constvel(Tr_n,a = self.k03,b = self.E3_R)
                "Cálculo do calor de reação"
                hr_n = (self.h1*k1_n*Ca_n) + (self.h2*k2_n*Cb_n) +(self.h3*k3_n*((Ca_n)**2))
                "Cálculo dos valores vetoriais usando RK de primeira ordem"
                Ca_n1 = Ca_n + H*((((qr)/(self.Vr))*(Ca0 - Ca_n))-(k1_n*Ca_n)-(k3_n*((Ca_n)**2)))
                Cb_n1 = Cb_n + H*(-((qr)/(self.Vr))*Cb_n+(k1_n*Ca_n)-(k2_n*Cb_n))
                Tr_n1 = Tr_n + H*((((qr)/(self.Vr))*(Tr0 - Tr_n))-((hr_n)/(self.ro_r*self.cp_r))+(((self.Ar*self.U)/(self.Vr*self.ro_r*self.cp_r))*(Tc_n - Tr_n)))
                Tc_n1 = Tc_n + H*(((1)/(self.mc*self.cp_c))*(Qc + (self.Ar*self.U*(Tr_n - Tc_n))))
                
                "valores das variáveis no tempo final"
                
                self.Ca_t = Ca_n1
                self.Cb_t = Cb_n1
                self.Tr_t = Tr_n1
                self.Tc_t = Tc_n1
                flag = False
            else:
                pass
    
    """Bloco que faz tudo"""
    
    "O bloco Simu_temporalCSTR realizada a simulação de um estado estacionário a partir de condições iniciais definidas"
    
    def Simu_temporalCSTR(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr,e_estac_fin,
                          arquivo_vs_t):
        
        
        #executar o processo iterativo e construir os vetores para plotagem dos gráficos
        self.__RK__estacionario(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr)
        
        
        print('-------------------------------------------------------')
        print('Dados da iteração')
        print('Cb (concentração de b no reator (kmol/m³)): {}'.format(self.Cbvst_plot[-1]))
        print('Ca (concentração de a no reator (kmol/m³)): {}'.format(self.Cavst_plot[-1]))
        print('Tr (temperatura do reator (K)): {}'.format(self.Trvst_plot[-1]))
        print('Tc (temperatura do fluido refrigerante (K)): {}'.format(self.Tcvst_plot[-1]))
        print('Tr0 (temperatura da corrente de alimentação (K)): {}'.format(Tr0))
        print('Ca0 (concentração de A na corrente de alimentação (kmol/m³)): {}'.format(Ca0))
        print('Qc (calor trocado com o fluido refrigerante antes da entrada no reator (kJ/min)): {}'.format(Qc))
        print('qr (vazão volumétrica da corrente de alimentação de reagentes (m³/min)): {}'.format(qr))
        print('-------------------------------------------------------')
        
        #plotagem do gráfico
        "Plots dos gráficos de concentração"
        fig, graf = plt.subplots()
        #subplot da concentração de A
        graf.plot(self.tempo_plot,self.Cavst_plot,label='CA (kmol/m³)')
        #subplot da concentração de B
        graf.plot(self.tempo_plot,self.Cbvst_plot,label='CB (kmol/m³)')
        graf.set_xlabel('tempo (min)')
        graf.set_ylabel('concentração molar (kmol/m³)')
        graf.set_title('concentração em função do tempo')
        graf.legend()
        
        "Plots dos gráficos de temperatura"
        fig, graf1 = plt.subplots()
        #subplot de valores de temperatura do reator Tr
        graf1.plot(self.tempo_plot,self.Trvst_plot,label='Tr (K)')
        #subplot de valores de temperatura de saída do fluido refrigerante Tc
        graf1.plot(self.tempo_plot,self.Tcvst_plot,label='Tc (K)')
        graf1.set_xlabel('tempo (min)')
        graf1.set_ylabel('Temperatura (K)')
        graf1.set_title('Temperatura em função do tempo')
        graf1.legend()
        
        "Salvar os dados de simulação em planilha"
        #dataframe com dados das variaveis em função do tempo
        planilha_variaveis_vs_t = pd.DataFrame(data = {'CA (kmol/m3)':self.Cavst_plot,
                                             'CB (kmol/m3)':self.Cbvst_plot,
                                             'Tr (K)':self.Trvst_plot,
                                             'Tc (K)':self.Tcvst_plot,
                                             't (min)':self.tempo_plot})
        planilha_variaveis_vs_t.to_csv(arquivo_vs_t,index = False)

        
        "Determinação do tempo para atingir o estado estacionário"
        aux1 = self.__tempo_estac_CB(e_estac_fin)
        a = (self.Cb_estac_i[-1]/self.Cbvst_plot[-1])*100
        
        
        "Exibição do relatório da simulação"
        print('--------------------------------------------------')
        print('Simulação do start up do reator CSTR')
        print('Cb atinge '+str(a)+'% do valor de estado estacionário final em '+str(aux1)+'min.')
        
    def Simu_estacionarioCBvsqr(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,criterio_erro,H,
                              liminf_qr,limsup_qr,liminf_Qc,limsup_Qc,H_qr,H_Qc,arquivo):
        "Este bloco realiza a simulação do estado estacionário para diferentes cenários de vazão volumétrica de alimentação, qr,"
        "e diferentes cenários de calor removido da corrente de fluido refrigerante, Qc."
        
        #Ca_inicial = concentração inicial de A no reator
        #Cb_inicial = concentração inicial de B no reator
        #Tr_inicial = Temperatura inicial do reator
        #Tc_inicial = Temperatura inicial do fluido refrigerante
        #Ca0 = concentração de A na corrente de alimentação
        #Tr0 = temperatura da corrente de alimentação
        #Qc = energia térmica removida da corrente de fluido refrigerante
        #criterio_erro = critério adotado para o estabelecimento da condição de estado estacionário
        #H = passo do método de resolução das equações diferenciais
        #liminf_qr = limite inferior da vazão volumétrica da corrente de alimentação
        #limsup_qr = limite superior da vazão volumétrica da corrente de alimentação
        #liminf_Qc = limite inferior da energia térmica removida do fluido refrigerante
        #limsup_Qc = limite superior da energia térmica removida do fluido refrigerante
        #H_qr = passo de varredura do intervalo de vazões volumétricas da corrente de alimentação
        #H_Qc = passo de varredura do intervalo de energia térmica removida do fluido refrigerante
        
        #Análise para diferentes valores de qr (para cada valor de qr é simulado todas as condições de Qc e determinado o
        # valor máximo de Cb que pode ser obtido no qr analisado), estamos tentando encontrar o ponto ótimo de produção de B
        
        #limpeza dos vetores auxiliares
        
        self.Cavsqr_plot.clear()
        self.Cbvsqr_plot.clear()
        self.Qc_Cbmax_vs_qr_plot.clear()
        self.Trvsqr_plot.clear()
        self.Tcvsqr_plot.clear()
        
        iteracao = 1
        
        qr = liminf_qr
        
        while (qr <= limsup_qr):
            
            Qc = liminf_Qc
            
            self.__RK__estacionario(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr)
            "--------------------------------"
            print('==============  iteracao: '+str(iteracao)+'  =================')
            print('Cb (kmol/m³): {}'.format(self.Cbvst_plot[-1]))
            print('Ca (kmol/m³): {}'.format(self.Cavst_plot[-1]))
            print('Tr (K): {}'.format(self.Trvst_plot[-1]))
            print('Tc (K): {}'.format(self.Tcvst_plot[-1]))
            print('Tr0 (K): {}'.format(Tr0))
            print('Ca0 (kmol/m³): {}'.format(Ca0))
            print('Qc (kJ/min): {}'.format(Qc))
            print('qr (m³/min): {}'.format(qr))
            "---------------------------------------"
            
            self.Cavsqr_max.append(self.Cavst_plot[-1])
            self.Cbvsqr_max.append(self.Cbvst_plot[-1])
            self.Trvsqr_max.append(self.Trvst_plot[-1])
            self.Tcvsqr_max.append(self.Tcvst_plot[-1])
            
            
            Qc = Qc + H_Qc
            
            iteracao = iteracao + 1
            
            
            while (Qc <= limsup_Qc):
        
                #executar o processo iterativo e construir os vetores para plotagem
                self.__RK__estacionario(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr)
                
                "--------------------------------"
                print('==============  iteracao: '+str(iteracao)+'  =================')
                print('Cb (kmol/m³): {}'.format(self.Cbvst_plot[-1]))
                print('Ca (kmol/m³): {}'.format(self.Cavst_plot[-1]))
                print('Tr (K): {}'.format(self.Trvst_plot[-1]))
                print('Tc (K): {}'.format(self.Tcvst_plot[-1]))
                print('Tr0 (K): {}'.format(Tr0))
                print('Ca0 (kmol/m³): {}'.format(Ca0))
                print('Qc (kJ/min): {}'.format(Qc))
                print('qr (m³/min): {}'.format(qr))
                "---------------------------------------"
                

                if self.Cbvst_plot[-1] >= self.Cbvsqr_max[-1]:
                    self.Cbvsqr_max.clear()
                    self.Cbvsqr_max.append(self.Cbvst_plot[-1])
                    
                    self.Qc_Cbmax_vs_qr.clear()
                    self.Qc_Cbmax_vs_qr.append(Qc)
            
                    self.Cavsqr_max.clear()
                    self.Cavsqr_max.append(self.Cavst_plot[-1])
                    
                    self.Trvsqr_max.clear()
                    self.Trvsqr_max.append(self.Trvst_plot[-1])
                
                    self.Tcvsqr_max.clear()
                    self.Tcvsqr_max.append(self.Tcvst_plot[-1])

                    
                Qc = Qc + H_Qc
                    
                iteracao = iteracao + 1

                
            # inserindo os dados para plotagem do gráfico de ca vs qr
            self.Cavsqr_plot.append(self.Cavsqr_max[-1])
            self.qrCa_plot.append(qr)
            
            #inserindo os dados para plotagem do gráfico de cb vs qr
            self.Cbvsqr_plot.append(self.Cbvsqr_max[-1])
            self.qrCb_plot.append(qr)
            
            #inserindo os dados para plotagem do gráfico de Qc vs qr
            self.Qc_Cbmax_vs_qr_plot.append(self.Qc_Cbmax_vs_qr[-1])
            #inserindo os dados para plotagem do gráfico de Tr vs qr
            self.Trvsqr_plot.append(self.Trvsqr_max[-1])
            self.qrTr_plot.append(qr)
            
            #inserindo os dados para plotagem do gráfico de Tc vs qr
            self.Tcvsqr_plot.append(self.Tcvsqr_max[-1])
            self.qrTc_plot.append(qr)
                
            qr = qr + H_qr    
        
        print('-------------------------------------------------------')
        
        #plotagem do gráfico
        "Plots dos gráficos de concentração"
        fig, graf = plt.subplots()
        #subplot da concentração de A_s
        graf.plot(self.qrCa_plot,self.Cavsqr_plot,label='CA_s (kmol/m³)')
        graf.set_xlabel('qr (m³/min)')
        graf.set_ylabel('concentração molar (kmol/m³)')
        graf.set_title('concentração de A em função de qr')
        graf.legend()
        
        fig, graf2 = plt.subplots()
        #subplot da concentração de A_s
        graf2.plot(self.qrCb_plot,self.Cbvsqr_plot,label='CB_s (kmol/m³)')
        graf2.set_xlabel('qr (m³/min)')
        graf2.set_ylabel('concentração molar (kmol/m³)')
        graf2.set_title('concentração máxima de B em função de qr')
        graf2.legend()
        
        "Plots dos gráficos de temperatura"
        fig, graf1 = plt.subplots()
        #subplot de valores de temperatura do reator Tr_s
        graf1.plot(self.qrTr_plot,self.Trvsqr_plot,label='Tr_s (K)')
        #subplot de valores de temperatura de saída do fluido refrigerante Tc_s
        graf1.plot(self.qrTc_plot,self.Tcvsqr_plot,label='Tc_s (K)')
        graf1.set_xlabel('qr (m³/min)')
        graf1.set_ylabel('Temperatura (K)')
        graf1.set_title('Temperaturas em função de qr')
        graf1.legend()
        
        "Plot do gráfico de Qc"
        fig, graf3 = plt.subplots()
        #subplot de valores de Qc_s
        graf3.plot(self.qrCb_plot,self.Qc_Cbmax_vs_qr_plot,label='Qc_s (kJ/min)')
        graf3.set_xlabel('qr (m³/min)')
        graf3.set_ylabel('Qc (kJ/min)')
        graf3.set_title('Qc em função de qr')
        graf3.legend()
        
        "Salvando os dados em planilhas"
        #planilha de dados de concentração em função de qr
        planilha_variaveis_vs_qr = pd.DataFrame({'CA_s(kmol/m3)':self.Cavsqr_plot,
                                                 'CB_s(kmol/m3)':self.Cbvsqr_plot,
                                                 'Tr_s(K)':self.Trvsqr_plot,
                                                 'Tc_s(K)':self.Tcvsqr_plot,
                                                 'qr_s(m3/min)':self.qrCa_plot,
                                                 'Qc_s(kJ/min)':self.Qc_Cbmax_vs_qr_plot})
        planilha_variaveis_vs_qr.to_csv(arquivo,index = False)
        
    def Simu_estacionarioCBvsQc(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,criterio_erro,H,
                              liminf_qr,limsup_qr,liminf_Qc,limsup_Qc,H_qr,H_Qc,arquivo):
        "Este bloco realiza a simulação do estado estacionário para diferentes cenários de vazão volumétrica de alimentação, qr,"
        "e diferentes cenários de calor removido da corrente de fluido refrigerante, Qc."
        
        #Ca_inicial = concentração inicial de A no reator
        #Cb_inicial = concentração inicial de B no reator
        #Tr_inicial = Temperatura inicial do reator
        #Tc_inicial = Temperatura inicial do fluido refrigerante
        #Ca0 = concentração de A na corrente de alimentação
        #Tr0 = temperatura da corrente de alimentação
        #Qc = energia térmica removida da corrente de fluido refrigerante
        #criterio_erro = critério adotado para o estabelecimento da condição de estado estacionário
        #H = passo do método de resolução das equações diferenciais
        #liminf_qr = limite inferior da vazão volumétrica da corrente de alimentação
        #limsup_qr = limite superior da vazão volumétrica da corrente de alimentação
        #liminf_Qc = limite inferior da energia térmica removida do fluido refrigerante
        #limsup_Qc = limite superior da energia térmica removida do fluido refrigerante
        #H_qr = passo de varredura do intervalo de vazões volumétricas da corrente de alimentação
        #H_Qc = passo de varredura do intervalo de energia térmica removida do fluido refrigerante
        
        #Análise para diferentes valores de qr (para cada valor de qr é simulado todas as condições de Qc e determinado o
        # valor máximo de Cb que pode ser obtido no qr analisado), estamos tentando encontrar o ponto ótimo de produção de B
        iteracao = 1
        
        Qc = liminf_Qc
        
        while (Qc <= limsup_Qc):
            
            qr = liminf_qr
            
            self.__RK__estacionario(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr)
            "--------------------------------"
            print('==============  iteracao: '+str(iteracao)+'  =================')
            print('Cb (kmol/m³): {}'.format(self.Cbvst_plot[-1]))
            print('Ca (kmol/m³): {}'.format(self.Cavst_plot[-1]))
            print('Tr (K): {}'.format(self.Trvst_plot[-1]))
            print('Tc (K): {}'.format(self.Tcvst_plot[-1]))
            print('Tr0 (K): {}'.format(Tr0))
            print('Ca0 (kmol/m³): {}'.format(Ca0))
            print('Qc (kJ/min): {}'.format(Qc))
            print('qr (m³/min): {}'.format(qr))
            "---------------------------------------"
            
            self.Cavsqr_max.append(self.Cavst_plot[-1])
            self.Cbvsqr_max.append(self.Cbvst_plot[-1])
            self.Trvsqr_max.append(self.Trvst_plot[-1])
            self.Tcvsqr_max.append(self.Tcvst_plot[-1])
            
            qr = qr + H_qr
            
            iteracao = iteracao + 1
            
            
            while (qr <= limsup_qr):
        
                #executar o processo iterativo e construir os vetores para plotagem
                self.__RK__estacionario(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr)
                
                "--------------------------------"
                print('==============  iteracao: '+str(iteracao)+'  =================')
                print('Cb (kmol/m³): {}'.format(self.Cbvst_plot[-1]))
                print('Ca (kmol/m³): {}'.format(self.Cavst_plot[-1]))
                print('Tr (K): {}'.format(self.Trvst_plot[-1]))
                print('Tc (K): {}'.format(self.Tcvst_plot[-1]))
                print('Tr0 (K): {}'.format(Tr0))
                print('Ca0 (kmol/m³): {}'.format(Ca0))
                print('Qc (kJ/min): {}'.format(Qc))
                print('qr (m³/min): {}'.format(qr))
                "---------------------------------------"
                

                if self.Cbvst_plot[-1] >= self.Cbvsqr_max[-1]:
                    self.Cbvsqr_max.clear()
                    self.Cbvsqr_max.append(self.Cbvst_plot[-1])
                    
                if self.Cavst_plot[-1] >= self.Cavsqr_max[-1]:
                    self.Cavsqr_max.clear()
                    self.Cavsqr_max.append(self.Cavst_plot[-1])
                    
                if self.Trvst_plot[-1] >= self.Trvsqr_max[-1]: 
                    self.Trvsqr_max.clear()
                    self.Trvsqr_max.append(self.Trvst_plot[-1])
                
                if self.Tcvst_plot[-1] >= self.Tcvsqr_max[-1]:
                    self.Tcvsqr_max.clear()
                    self.Tcvsqr_max.append(self.Tcvst_plot[-1])

                    
                qr = qr + H_qr
                    
                iteracao = iteracao + 1

                
            # inserindo os dados para plotagem do gráfico de ca vs qr
            self.Cavsqr_plot.append(self.Cavsqr_max[-1])
            self.qrCa_plot.append(Qc)
            
            #inserindo os dados para plotagem do gráfico de cb vs qr
            self.Cbvsqr_plot.append(self.Cbvsqr_max[-1])
            self.qrCb_plot.append(Qc)
            
            #inserindo os dados para plotagem do gráfico de Tr vs qr
            self.Trvsqr_plot.append(self.Trvsqr_max[-1])
            self.qrTr_plot.append(Qc)
            
            #inserindo os dados para plotagem do gráfico de Tc vs qr
            self.Tcvsqr_plot.append(self.Tcvsqr_max[-1])
            self.qrTc_plot.append(Qc)
                
            Qc = Qc + H_Qc  
        
        print('-------------------------------------------------------')
        
        #plotagem do gráfico
        "Plots dos gráficos de concentração"
        fig, graf = plt.subplots()
        #subplot da concentração de A_s
        graf.plot(self.qrCa_plot,self.Cavsqr_plot,label='CA_s (kmol/m³)')
        graf.set_xlabel('Qc (kJ/min)')
        graf.set_ylabel('concentração molar (kmol/m³)')
        graf.set_title('concentração de A em função de Qc')
        graf.legend()
        
        fig, graf2 = plt.subplots()
        #subplot da concentração de A_s
        graf2.plot(self.qrCb_plot,self.Cbvsqr_plot,label='CB_s (kmol/m³)')
        graf2.set_xlabel('Qc (kJ/min)')
        graf2.set_ylabel('concentração molar (kmol/m³)')
        graf2.set_title('concentração de B em função de Qc')
        graf2.legend()
        
        "Plots dos gráficos de temperatura"
        fig, graf1 = plt.subplots()
        #subplot de valores de temperatura do reator Tr_s
        graf1.plot(self.qrTr_plot,self.Trvsqr_plot,label='Tr_s (K)')
        #subplot de valores de temperatura de saída do fluido refrigerante Tc_s
        graf1.plot(self.qrTc_plot,self.Tcvsqr_plot,label='Tc_s (K)')
        graf1.set_xlabel('Qc (kJ/min)')
        graf1.set_ylabel('Temperatura (K)')
        graf1.set_title('Temperatura em função de Qc')
        graf1.legend()
        
        "Salvando os dados em planilhas"
        #planilha de dados de concentração em função de qr
        planilha_variaveis_vs_qr = pd.DataFrame({'CA (kmol/m3':self.Cavsqr_plot,
                                                 'CB (kmol/m3':self.Cbvsqr_plot,
                                                 'Tr (K)':self.Trvsqr_plot,
                                                 'Tc (K)':self.Tcvsqr_plot,
                                                 'Qc (kJ/min)':self.qrCa_plot})
        planilha_variaveis_vs_qr.to_csv(arquivo,index = False)
        
    def Simu_superficieCB(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,criterio_erro,H,
                              liminf_qr,limsup_qr,liminf_Qc,limsup_Qc,H_qr,H_Qc,arquivo):
        "Este bloco realiza a simulação do estado estacionário para diferentes cenários de vazão volumétrica de alimentação, qr,"
        "e diferentes cenários de calor removido da corrente de fluido refrigerante, Qc."
        
        #Ca_inicial = concentração inicial de A no reator
        #Cb_inicial = concentração inicial de B no reator
        #Tr_inicial = Temperatura inicial do reator
        #Tc_inicial = Temperatura inicial do fluido refrigerante
        #Ca0 = concentração de A na corrente de alimentação
        #Tr0 = temperatura da corrente de alimentação
        #Qc = energia térmica removida da corrente de fluido refrigerante
        #criterio_erro = critério adotado para o estabelecimento da condição de estado estacionário
        #H = passo do método de resolução das equações diferenciais
        #liminf_qr = limite inferior da vazão volumétrica da corrente de alimentação
        #limsup_qr = limite superior da vazão volumétrica da corrente de alimentação
        #liminf_Qc = limite inferior da energia térmica removida do fluido refrigerante
        #limsup_Qc = limite superior da energia térmica removida do fluido refrigerante
        #H_qr = passo de varredura do intervalo de vazões volumétricas da corrente de alimentação
        #H_Qc = passo de varredura do intervalo de energia térmica removida do fluido refrigerante
        
        #Análise para diferentes valores de qr (para cada valor de qr é simulado todas as condições de Qc e determinado o
        # valor máximo de Cb que pode ser obtido no qr analisado), estamos tentando encontrar o ponto ótimo de produção de B
        
        #esvaziamento dos vetores auxiliares
        
        self.Cb_sup_plot.clear()
        self.Ca_sup_plot.clear()
        self.Tr_sup_plot.clear()
        self.Tc_sup_plot.clear()
        self.qr_sup_plot.clear()
        self.Qc_sup_plot.clear()
        
        iteracao = 1
        
        qr = liminf_qr
        
        while (qr <= limsup_qr):
            
            Qc = liminf_Qc
            
            self.__RK__estacionario(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr)
            "--------------------------------"
            print('==============  iteracao: '+str(iteracao)+'  =================')
            print('Cb (kmol/m³): {}'.format(self.Cbvst_plot[-1]))
            print('Ca (kmol/m³): {}'.format(self.Cavst_plot[-1]))
            print('Tr (K): {}'.format(self.Trvst_plot[-1]))
            print('Tc (K): {}'.format(self.Tcvst_plot[-1]))
            print('Tr0 (K): {}'.format(Tr0))
            print('Ca0 (kmol/m³): {}'.format(Ca0))
            print('Qc (kJ/min): {}'.format(Qc))
            print('qr (m³/min): {}'.format(qr))
            "---------------------------------------"
            
            self.Cb_sup_plot.append(self.Cbvst_plot[-1])
            self.Ca_sup_plot.append(self.Cavst_plot[-1])
            self.Tr_sup_plot.append(self.Trvst_plot[-1])
            self.Tc_sup_plot.append(self.Tcvst_plot[-1])
            self.qr_sup_plot.append(qr)
            self.Qc_sup_plot.append(Qc)
            
            
            Qc = Qc + H_Qc
            
            iteracao = iteracao + 1
            
            
            while (Qc <= limsup_Qc):
        
                #executar o processo iterativo e construir os vetores para plotagem
                self.__RK__estacionario(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr)
                
                "--------------------------------"
                print('==============  iteracao: '+str(iteracao)+'  =================')
                print('Cb (kmol/m³): {}'.format(self.Cbvst_plot[-1]))
                print('Ca (kmol/m³): {}'.format(self.Cavst_plot[-1]))
                print('Tr (K): {}'.format(self.Trvst_plot[-1]))
                print('Tc (K): {}'.format(self.Tcvst_plot[-1]))
                print('Tr0 (K): {}'.format(Tr0))
                print('Ca0 (kmol/m³): {}'.format(Ca0))
                print('Qc (kJ/min): {}'.format(Qc))
                print('qr (m³/min): {}'.format(qr))
                "---------------------------------------"
                
                self.Cb_sup_plot.append(self.Cbvst_plot[-1])
                self.Ca_sup_plot.append(self.Cavst_plot[-1])
                self.Tr_sup_plot.append(self.Trvst_plot[-1])
                self.Tc_sup_plot.append(self.Tcvst_plot[-1])
                self.qr_sup_plot.append(qr)
                self.Qc_sup_plot.append(Qc)
                
                Qc = Qc + H_Qc
                
                iteracao = iteracao + 1

            "Plots dos gráficos de concentração vs Qc"
            fig, graf = plt.subplots()
            #subplot da concentração de A_s
            graf.plot(self.Qc_sup_plot,self.Ca_sup_plot,label='CA_s (kmol/m³)')
            graf.set_xlabel('Qc (kJ/min)')
            graf.set_ylabel('concentração molar (kmol/m³)')
            graf.set_title('CA em função de Qc '+'(qr = '+str(qr)+')')
            graf.legend()
            
            fig, graf2 = plt.subplots()
            #subplot da concentração de A_s
            graf2.plot(self.Qc_sup_plot,self.Cb_sup_plot,label='CB_s (kmol/m³)')
            graf2.set_xlabel('Qc (kJ/min)')
            graf2.set_ylabel('concentração molar (kmol/m³)')
            graf2.set_title('CB em função de Qc '+'(qr = '+str(qr)+')')
            graf2.legend()
            
            
            "Plots dos gráficos de temperatura"
            fig, graf3 = plt.subplots()
            #subplot de valores de temperatura do reator Tr_s
            graf3.plot(self.Qc_sup_plot,self.Tr_sup_plot,label='Tr_s (K)')
            #subplot de valores de temperatura de saída do fluido refrigerante Tc_s
            graf3.plot(self.Qc_sup_plot,self.Tc_sup_plot,label='Tc_s (K)')
            graf3.set_xlabel('Qc (kJ/min)')
            graf3.set_ylabel('Temperatura (K)')
            graf3.set_title('T em função de Qc '+'(qr = '+str(qr)+')')
            graf3.legend()                     
                
            qr = qr + H_qr    
        
        print('-------------------------------------------------------')
        
        
        "Salvando os dados em planilhas"
        #planilha de dados de concentração em função de qr
        planilha_superficie_cb = pd.DataFrame({'CB_s(kmol/m3)':self.Cb_sup_plot,
                                               'CA_s(kmol/m3)':self.Ca_sup_plot,
                                               'Tr_s(K)':self.Tr_sup_plot,
                                               'Tc_s(K)':self.Tc_sup_plot,
                                               'qr_s(m3/min)':self.qr_sup_plot,
                                               'Qc_s(kJ/min)':self.Qc_sup_plot})
        planilha_superficie_cb.to_csv(arquivo,index = False)
    "----------------------------------------------------------------------------------------------------"
    
    def Simu_superficieCB_unitario(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,criterio_erro,H,
                              qr,liminf_Qc,limsup_Qc,H_Qc,arquivo):
        "Este bloco realiza a simulação do estado estacionário para diferentes cenários de vazão volumétrica de alimentação, qr,"
        "e diferentes cenários de calor removido da corrente de fluido refrigerante, Qc."
        
        #Ca_inicial = concentração inicial de A no reator
        #Cb_inicial = concentração inicial de B no reator
        #Tr_inicial = Temperatura inicial do reator
        #Tc_inicial = Temperatura inicial do fluido refrigerante
        #Ca0 = concentração de A na corrente de alimentação
        #Tr0 = temperatura da corrente de alimentação
        #Qc = energia térmica removida da corrente de fluido refrigerante
        #criterio_erro = critério adotado para o estabelecimento da condição de estado estacionário
        #H = passo do método de resolução das equações diferenciais
        #liminf_Qc = limite inferior da energia térmica removida do fluido refrigerante
        #limsup_Qc = limite superior da energia térmica removida do fluido refrigerante
        #H_Qc = passo de varredura do intervalo de energia térmica removida do fluido refrigerante
        
        #Análise para diferentes valores de qr (para cada valor de qr é simulado todas as condições de Qc e determinado o
        # valor máximo de Cb que pode ser obtido no qr analisado), estamos tentando encontrar o ponto ótimo de produção de B
        
        #esvaziamento dos vetores auxiliares
        
        self.Cb_sup_plot.clear()
        self.Ca_sup_plot.clear()
        self.Tr_sup_plot.clear()
        self.Tc_sup_plot.clear()
        self.qr_sup_plot.clear()
        self.Qc_sup_plot.clear()
        
        iteracao = 1
        
            
        Qc = liminf_Qc
            
        self.__RK__estacionario(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr)
        "--------------------------------"
        print('==============  iteracao: '+str(iteracao)+'  =================')
        print('Cb (kmol/m³): {}'.format(self.Cbvst_plot[-1]))
        print('Ca (kmol/m³): {}'.format(self.Cavst_plot[-1]))
        print('Tr (K): {}'.format(self.Trvst_plot[-1]))
        print('Tc (K): {}'.format(self.Tcvst_plot[-1]))
        print('Tr0 (K): {}'.format(Tr0))
        print('Ca0 (kmol/m³): {}'.format(Ca0))
        print('Qc (kJ/min): {}'.format(Qc))
        print('qr (m³/min): {}'.format(qr))
        "---------------------------------------"
        
        self.Cb_sup_plot.append(self.Cbvst_plot[-1])
        self.Ca_sup_plot.append(self.Cavst_plot[-1])
        self.Tr_sup_plot.append(self.Trvst_plot[-1])
        self.Tc_sup_plot.append(self.Tcvst_plot[-1])
        self.qr_sup_plot.append(qr)
        self.Qc_sup_plot.append(Qc)
        
        
        Qc = Qc + H_Qc
        
        iteracao = iteracao + 1
        
        
        while (Qc <= limsup_Qc):
    
            #executar o processo iterativo e construir os vetores para plotagem
            self.__RK__estacionario(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr)
            
            "--------------------------------"
            print('==============  iteracao: '+str(iteracao)+'  =================')
            print('Cb (kmol/m³): {}'.format(self.Cbvst_plot[-1]))
            print('Ca (kmol/m³): {}'.format(self.Cavst_plot[-1]))
            print('Tr (K): {}'.format(self.Trvst_plot[-1]))
            print('Tc (K): {}'.format(self.Tcvst_plot[-1]))
            print('Tr0 (K): {}'.format(Tr0))
            print('Ca0 (kmol/m³): {}'.format(Ca0))
            print('Qc (kJ/min): {}'.format(Qc))
            print('qr (m³/min): {}'.format(qr))
            "---------------------------------------"
            
            self.Cb_sup_plot.append(self.Cbvst_plot[-1])
            self.Ca_sup_plot.append(self.Cavst_plot[-1])
            self.Tr_sup_plot.append(self.Trvst_plot[-1])
            self.Tc_sup_plot.append(self.Tcvst_plot[-1])
            self.qr_sup_plot.append(qr)
            self.Qc_sup_plot.append(Qc)
            
            Qc = Qc + H_Qc
            
            iteracao = iteracao + 1
            
        "Plots dos gráficos de concentração vs Qc"
        fig, graf = plt.subplots()
        #subplot da concentração de A_s
        graf.plot(self.Qc_sup_plot,self.Ca_sup_plot,label='CA_s (kmol/m³)')
        graf.set_xlabel('Qc (kJ/min)')
        graf.set_ylabel('concentração molar (kmol/m³)')
        graf.set_title('CA em função de Qc '+'(qr = '+str(qr)+')')
        graf.legend()
        fig.savefig('qr_'+str(qr)+'_CA_.png', format='png')
        plt.close(fig)
        
        fig, graf2 = plt.subplots()
        #subplot da concentração de A_s
        graf2.plot(self.Qc_sup_plot,self.Cb_sup_plot,label='CB_s (kmol/m³)')
        graf2.set_xlabel('Qc (kJ/min)')
        graf2.set_ylabel('concentração molar (kmol/m³)')
        graf2.set_title('CB em função de Qc '+'(qr = '+str(qr)+')')
        graf2.legend()
        fig.savefig('qr_'+str(qr)+'_CB_.png', format='png')
        plt.close(fig)
            
            
        "Plots dos gráficos de temperatura"
        fig, graf3 = plt.subplots()
        #subplot de valores de temperatura do reator Tr_s
        graf3.plot(self.Qc_sup_plot,self.Tr_sup_plot,label='Tr_s (K)')
        #subplot de valores de temperatura de saída do fluido refrigerante Tc_s
        graf3.plot(self.Qc_sup_plot,self.Tc_sup_plot,label='Tc_s (K)')
        graf3.set_xlabel('Qc (kJ/min)')
        graf3.set_ylabel('Temperatura (K)')
        graf3.set_title('T em função de Qc '+'(qr = '+str(qr)+')')
        graf3.legend() 
        fig.savefig('qr_'+str(qr)+'_Tr_Tc_.png', format='png')
        plt.close(fig)
         
        
        print('-------------------------------------------------------')
        
        
        "Salvando os dados em planilhas"
        #planilha de dados de concentração em função de qr
        planilha_superficie_cb = pd.DataFrame({'CB_s(kmol/m3)':self.Cb_sup_plot,
                                               'CA_s(kmol/m3)':self.Ca_sup_plot,
                                               'Tr_s(K)':self.Tr_sup_plot,
                                               'Tc_s(K)':self.Tc_sup_plot,
                                               'qr_s(m3/min)':self.qr_sup_plot,
                                               'Qc_s(kJ/min)':self.Qc_sup_plot})
        planilha_superficie_cb.to_csv(arquivo,index = False)
    "----------------------------------------------------------------------------------------------------"
    
    def tempo_estacionario(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,t_final,H,qr,e_estac_fin):
        
        """A simulação temporal fornece os plots de temperatura do fluido refrigerante em contato com o reator (Tc) e 
        a temperatura do reator (Tr). Além disso, fornece o tempo que o sistema leva para atingir um determinado valor em relação
        ao estado estacionário final."""
        
        #executar o processo iterativo e construir os vetores para plotagem
        self.__RK_tempo_amostr(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,t_final,H,qr)
        
        #plotagem do gráfico
        "Plots dos gráficos de concentração"
        fig, graf = plt.subplots()
        #subplot da concentração de A
        graf.plot(self.tempo_tamostr,self.Cavst_tamostr,label='CA (kmol/m³)')
        #subplot da concentração de B
        graf.plot(self.tempo_tamostr,self.Cbvst_tamostr,label='CB (kmol/m³)')
        graf.set_xlabel('tempo (min)')
        graf.set_ylabel('concentração molar (kmol/m³)')
        graf.set_title('concentração em função do tempo')
        graf.legend()
        
    "----------------------------------------------------------------------------------------------------"
    
    def gerar_dados(self,n_iter,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,t_final,H,qr,e_estac_fin,
                       nome_arq,t_perturbacao,porc_inf,porc_sup):
        "Bloco para geração de dados para treinamento das máquinas de aprendizado"
        "Esse bloco investiga como a concentração CB é influenciada por perturbações em Qr"
        
        #os argumentos dessa função possuem o seguinte significado:
        #n_iter = quantidade de dados que devem ser gerados nas perturbações
        #t_final = tempo final da simulação de start-up do reator (min)
        #H = passo do método iterativo de RK (min)
        #Ca_inicial = concentração de A no interior do reator no start-up (kmol/m³)
        #Cb_inicial = concentração inicial de B no interior do reator no start-up (kmol/m³)
        #Tr_inicial = temperatura inicial no interior do reator no start-up (K)
        #Tc_inicial = temperatura inicial do fluido refrigerante na saída no start-up (K)
        #Ca0 = concentração de A na corrente de alimentação do reator no start-up (kmol/m³)
        #Tr0 = temperatura da corrente de alimentação de reagentes no start-up (K)
        #Qc = calor retirado ou fornecido pelo fluido refrigerante no start-up (kJ/min)
        #porc_inf = porcentagem de perturbação inferior da variável em uso (%)
        #porc_sup = porcentagem de perturbação superior da variável em uso (%)
        #qr: vazão volumétrica da alimentação de reagentes (m³/min)
        #e_estac_fin = porcentagem utilizada para determinação do tempo que a variável avaliada leva pra apresentar um erro de x% do valor do estado estacionário final
        #tipo_pert = tipo de perturbação: somente qr ou qr e Qc, simultaneamente, ou apenas Qc (qr, qrQc, Qc)
        #count_f = número de simulações para determinação do tempo ótimo de amostragem
        #t_final_amostr = tempo final de amostragem para determinação do tempo ótimo de amostragem (min)
        #t_perturbacao = tempo de cada perturbação (min)
        
        "Análise do estado estacionário inicial(start-up do reator) e geração de dados estacionários"
        self.Simu_temporalCSTR(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,t_final,H,qr,e_estac_fin)
        
        Cb_atual = self.Cbvst_plot[-1] #valor de estado estacionário da concentração de B (kmol/m³)
        Ca_atual = self.Cavst_plot[-1] #valor de estado estacionário da concentração de A (kmol/m³)
        Tr_atual = self.Trvst_plot[-1] #valor de estado estacionário da temperatura do reator (K)
        Tc_atual = self.Tcvst_plot[-1] #valor de estado estacionário da temperatura do fluido refrigerante em contato com o reator (K)
        qr_atual = qr #Vazão volumétrica da corrente de alimentação de reagentes (m³/min)
        Tr0_atual = Tr0
        Ca0_atual = Ca0
        Qc_atual = Qc
        
        "Dados da influência de perturbações em qr sobre Cb"
        
        "Vetores contendo os resultados da simulação dinâmica"
        
        dados_Tr0_k2 = [] #temperatura da corrente de alimentação no instante t-2
        dados_Tr0_k1 = [] #temperatura da corrente de alimentação no instante t-1
        dados_Ca0_k2 = [] #concentração de A na alimentação no instante t-2
        dados_Ca0_k1 = [] #concentração de A na alimentação no instante t-1
        dados_Qc_k2 = [] #calor retirado ou fornecido pelo fluido refrigerante no instante t-2
        dados_Qc_k1 = [] #calor retirado ou fornecido pelo fluido refrigerante no instante t-2
        dados_qr_k2 = [] #vazão volumétrica da alimentação de reagentes no instante t-2
        dados_qr_k1 = [] #vazão volumétrica da alimentação de reagentes no instante t-1
        dados_Cb_k2 = [] #concentração de B no reator no instante t-2
        dados_Cb_k1 = [] #concentração de B no reator no instante t-1
        dados_Cb_k = [] #concentração de B no reator no instante t
        dados_t = [] #dados de tempo para plotagem do perfil temporal
        
        #construção do vetor da primeira iteração
        t = 0
        dados_Cb_k2.append(Cb_atual)
        dados_Tr0_k2.append(Tr0_atual)
        dados_Ca0_k2.append(Ca0_atual)
        dados_qr_k2.append(qr_atual)
        dados_Qc_k2.append(Qc_atual)
        
        #construção do vetor da segunda iteração
        
        #bloco para seleção de qual variável será perturbada
        
        porc_perturbacao_Tr0 = 0
        porc_perturbacao_Ca0 = 0
        porc_perturbacao_Qc = 0
        porc_perturbacao_qr = 0
        
        boolean_Tr0 = False
        boolean_Ca0 = False
        boolean_Qc = False
        boolean_qr = False
        
        ativacao_Tr0 = random.randrange(0,2)         
        ativacao_Ca0 = random.randrange(0,2)
        ativacao_Qc = random.randrange(0,2)
        ativacao_qr = random.randrange(0,2)
        
        if ativacao_Tr0 == 1:
            boolean_Tr0 = True
        if ativacao_Ca0 == 1:
            boolean_Ca0 = True
        if ativacao_Qc == 1:
            boolean_Qc = True
        if ativacao_qr == 1:
            boolean_qr = True
        
        if boolean_Tr0:
            porc_perturbacao_Tr0 = random.uniform(porc_inf,porc_sup)
            Tr0_atual = Tr0_atual + Tr0_atual*((porc_perturbacao_Tr0)/(100))
        if boolean_Ca0:
            porc_perturbacao_Ca0 = random.uniform(porc_inf,porc_sup)
            Ca0_atual = Ca0_atual + Ca0_atual*((porc_perturbacao_Ca0)/(100))
        if boolean_Qc:
            porc_perturbacao_Qc = random.uniform(porc_inf,porc_sup)
            Qc_atual = Qc_atual + Qc_atual*((porc_perturbacao_Qc)/(100))
        if boolean_qr:
            porc_perturbacao_qr = random.uniform(porc_inf,porc_sup)
            qr_atual = qr_atual + qr_atual*((porc_perturbacao_qr)/(100))
        
        dados_Tr0_k1.append(Tr0_atual)
        dados_Ca0_k1.append(Ca0_atual)
        dados_qr_k1.append(qr_atual)
        dados_Qc_k1.append(Qc_atual)
        dados_Cb_k1.append(Cb_atual)
        
        self.__RK_gerar_dados(Ca_atual,Cb_atual,Tr_atual,Tc_atual,Ca0_atual,Tr0_atual,
                              Qc_atual,t_perturbacao,H,qr_atual)
        
        
        Ca_atual = self.Ca_t
        Cb_atual = self.Cb_t
        Tr_atual = self.Tr_t
        Tc_atual = self.Tc_t
        dados_Cb_k.append(Cb_atual)
        t = t + t_perturbacao
        dados_t.append(t)
        
        
        "Bloco para geração de dados para perturbação"
        n = 1
        while (n <= (n_iter - 2)):
            #bloco para seleção de qual variável será perturbada
            
            dados_Cb_k2.append(dados_Cb_k1[-1])
            dados_Tr0_k2.append(dados_Tr0_k1[-1])
            dados_Ca0_k2.append(dados_Ca0_k1[-1])
            dados_qr_k2.append(dados_qr_k1[-1])
            dados_Qc_k2.append(dados_Qc_k2[-1])
            
            porc_perturbacao_Tr0 = 0
            porc_perturbacao_Ca0 = 0
            porc_perturbacao_Qc = 0
            porc_perturbacao_qr = 0
            
            boolean_Tr0 = False
            boolean_Ca0 = False
            boolean_Qc = False
            boolean_qr = False
            
            ativacao_Tr0 = random.randrange(0,2)         
            ativacao_Ca0 = random.randrange(0,2)
            ativacao_Qc = random.randrange(0,2)
            ativacao_qr = random.randrange(0,2)
            
            if ativacao_Tr0 == 1:
                boolean_Tr0 = True
            if ativacao_Ca0 == 1:
                boolean_Ca0 = True
            if ativacao_Qc == 1:
                boolean_Qc = True
            if ativacao_qr == 1:
                boolean_qr = True
            
            if boolean_Tr0:
                porc_perturbacao_Tr0 = random.uniform(porc_inf,porc_sup)
                Tr0_atual = Tr0_atual + Tr0_atual*((porc_perturbacao_Tr0)/(100))
            if boolean_Ca0:
                porc_perturbacao_Ca0 = random.uniform(porc_inf,porc_sup)
                Ca0_atual = Ca0_atual + Ca0_atual*((porc_perturbacao_Ca0)/(100))
            if boolean_Qc:
                porc_perturbacao_Qc = random.uniform(porc_inf,porc_sup)
                Qc_atual = Qc_atual + Qc_atual*((porc_perturbacao_Qc)/(100))
            if boolean_qr:
                porc_perturbacao_qr = random.uniform(porc_inf,porc_sup)
                qr_atual = qr_atual + qr_atual*((porc_perturbacao_qr)/(100))
            
            dados_Tr0_k1.append(Tr0_atual)
            dados_Ca0_k1.append(Ca0_atual)
            dados_qr_k1.append(qr_atual)
            dados_Qc_k1.append(Qc_atual)
            dados_Cb_k1.append(dados_Cb_k[-1])
                
            self.__RK_gerar_dados(Ca_atual,Cb_atual,Tr_atual,Tc_atual,Ca0_atual,Tr0_atual,
                                  Qc_atual,t_perturbacao,H,qr_atual)
            
            Ca_atual = self.Ca_t
            Cb_atual = self.Cb_t
            Tr_atual = self.Tr_t
            Tc_atual = self.Tc_t
            
            dados_Cb_k.append(Cb_atual)
            t = t + t_perturbacao
            dados_t.append(t)
            
            n = n + 1
            
        "Salvar os dados em planilhas de Excel"
        
        d = {'qr (k-2)':dados_qr_k2,'qr (k-1)':dados_qr_k1,'Qc (k-2)':dados_Qc_k2,'Qc (k-1)':dados_Qc_k1,'Ca0 (k-2)':dados_Ca0_k2,
             'Ca0 (k-1)':dados_Ca0_k1,'Tr0 (k-2)':dados_Tr0_k2,'Tr0 (k-1)':dados_Tr0_k1,'Cb (k-2)':dados_Cb_k2,'Cb (k-1)':dados_Cb_k1,
             'Cb (k)':dados_Cb_k}
        dados = pd.DataFrame(data = d)
        dados.to_csv(nome_arq,index=False)
        
        "Gráficos dos perfis temporais da simulação dinâmica"
        
        fig, graf2 = plt.subplots()
        #subplot de valores temporais de Cb (k-2)
        graf2.plot(dados_t,dados_Cb_k,label='Cb (kmol/m³)')
        graf2.set_xlabel('tempo (min)')
        graf2.set_ylabel('Cb (k) kmol/m³')
        graf2.set_title('Cb (k) em função do tempo')
        graf2.legend()
        
        fig, graf3 = plt.subplots()
        #subplot de valores temporais de qr (k-1)
        graf3.plot(dados_t,dados_qr_k1,label='qr (m³/min)')
        graf3.set_xlabel('tempo (min)')
        graf3.set_ylabel('qr (k-1) m³/min')
        graf3.set_title('qr (k-1) em função do tempo')
        graf3.legend()
        
        fig, graf4 = plt.subplots()
        #subplot de valores temporais de Qc (k-1)
        graf4.plot(dados_t,dados_Qc_k1,label='Qc (kJ/min)')
        graf4.set_xlabel('tempo (min)')
        graf4.set_ylabel('Qc (k-1) Qj/min')
        graf4.set_title('Qc (k-1) em função do tempo')
        graf4.legend()
        
        fig, graf5 = plt.subplots()
        #subplot de valores temporais de Tr0 (k-1)
        graf5.plot(dados_t,dados_Tr0_k1,label='Tr0 (Kelvin)')
        graf5.set_xlabel('tempo (min)')
        graf5.set_ylabel('Tr0 (k-1) Kelvin')
        graf5.set_title('Tr0 (k-1) em função do tempo')
        graf5.legend()
        
        fig, graf6 = plt.subplots()
        #subplot de valores temporais de Ca0 (k-1)
        graf6.plot(dados_t,dados_Ca0_k1,label='Ca0 (kmol/m³)')
        graf6.set_xlabel('tempo (min)')
        graf6.set_ylabel('Ca0 (k-1) kmol/m³')
        graf6.set_title('Ca0 (k-1) em função do tempo')
        graf6.legend()
        
        #"bloco para cálculo do tempo ótimo de amostragem"
        #nesse bloco são realizadas perturbações e a determinação do tempo que o sistema leva para
        #atingir um novo estado estacionário a cada perturbação. É construída uma lista e o valor ótimo
        #é definido como sendo a média dos tempos de cada perturbação.
        #count_i = 0 
        #count_f = 100 #número de simulações
        #while (count_i <= count_f):
            #porc_perturbacao_qr = random.uniform(porc_inf_qr,porc_sup_qr)
            #qr_1 = qr + qr*((porc_perturbacao_qr)/(100))
            
            #geração dos dados de simulação
            #self.__RK_tempo_amostr(Ca_estac,Cb_estac,Tr_estac,Tc_estac,Ca0,Tr0,Qc,t_final_amostr,H,qr_1)
            
            #Determinação do tempo com erro de e_estac_fin% do valor de estado estacionário final
            #aux1 = self.__tempo_estac_CB_tamostr(e_estac_fin) #fator 60 converte resultado em minutos para segundos
            #t_amostra.append(aux1)
            #count_i = count_i + 1
        
            #gerando dados em segundos para visualização
            #t_amostra_segundos = []
            #for i in t_amostra:
                #t_amostra_segundos.append(60*i)