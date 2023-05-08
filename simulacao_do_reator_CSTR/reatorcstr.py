import numpy as np
import matplotlib.pyplot as plt
import random as random
import pandas as pd
import psycopg2
from decimal import Decimal

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
        self.Tr0vst_plot = []
        self.Ca0vst_plot = []
        self.Qcvst_plot = []
        self.qrvst_plot = []
        self.tempo_plot = []
        self.contador = []
        
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
        
        "vetores auxiliares (Simu_perturbacoes_unitario)"
        self.Cavst_simuperturbacoesunitario = []
        self.Cbvst_simuperturbacoesunitario = []
        self.Trvst_simuperturbacoesunitario = []
        self.Tcvst_simuperturbacoesunitario = []
        self.tempo_simuperturbacoesunitario = []
        
        "vetores auxiliares da geração de dados de soft sensor"
        #20 segundos
        self.Cavst_20s = []
        self.Cbvst_20s = []
        self.Tcvst_20s = []
        self.Trvst_20s = []
        self.Tr0vst_20s = []
        self.Ca0vst_20s = []
        self.Qcvst_20s = []
        self.qrvst_20s = []
        self.tempo_20s = []
        
        #75 segundos
        self.Cavst_75s = []
        self.Cbvst_75s = []
        self.Tcvst_75s = []
        self.Trvst_75s = []
        self.Tr0vst_75s = []
        self.Ca0vst_75s = []
        self.Qcvst_75s = []
        self.qrvst_75s = []
        self.tempo_75s = []
        
        #95 segundos
        self.Cavst_95s = []
        self.Cbvst_95s = []
        self.Tcvst_95s = []
        self.Trvst_95s = []
        self.Tr0vst_95s = []
        self.Ca0vst_95s = []
        self.Qcvst_95s = []
        self.qrvst_95s = []
        self.tempo_95s = []
        
        #115 segundos
        self.Cavst_115s = []
        self.Cbvst_115s = []
        self.Tcvst_115s = []
        self.Trvst_115s = []
        self.Tr0vst_115s = []
        self.Ca0vst_115s = []
        self.Qcvst_115s = []
        self.qrvst_115s = []
        self.tempo_115s = []
             
    
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
    
    def __RK__gera_dados_softsensor(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr,
                                    t_inicial,contador):
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
        
        t = t_inicial
        #t = round(t,3)
        
        count = contador
        
        self.Cavst_plot.clear()
        self.Cbvst_plot.clear()
        self.Tcvst_plot.clear()
        self.Trvst_plot.clear()
        self.Tr0vst_plot.clear()
        self.Ca0vst_plot.clear()
        self.Qcvst_plot.clear()
        self.qrvst_plot.clear()
        self.tempo_plot.clear()
        self.contador.clear()
        '-------------------start-------------'
        #20 segundos
        self.Cavst_20s.clear()
        self.Cbvst_20s.clear()
        self.Tcvst_20s.clear()
        self.Trvst_20s.clear()
        self.Tr0vst_20s.clear()
        self.Ca0vst_20s.clear()
        self.Qcvst_20s.clear()
        self.qrvst_20s.clear()
        self.tempo_20s.clear()
        
        #75 segundos
        self.Cavst_75s.clear()
        self.Cbvst_75s.clear()
        self.Tcvst_75s.clear()
        self.Trvst_75s.clear()
        self.Tr0vst_75s.clear()
        self.Ca0vst_75s.clear()
        self.Qcvst_75s.clear()
        self.qrvst_75s.clear()
        self.tempo_75s.clear()
        
        #95 segundos
        self.Cavst_95s.clear()
        self.Cbvst_95s.clear()
        self.Tcvst_95s.clear()
        self.Trvst_95s.clear()
        self.Tr0vst_95s.clear()
        self.Ca0vst_95s.clear()
        self.Qcvst_95s.clear()
        self.qrvst_95s.clear()
        self.tempo_95s.clear()
        
        #115 segundos
        self.Cavst_115s.clear()
        self.Cbvst_115s.clear()
        self.Tcvst_115s.clear()
        self.Trvst_115s.clear()
        self.Tr0vst_115s.clear()
        self.Ca0vst_115s.clear()
        self.Qcvst_115s.clear()
        self.qrvst_115s.clear()
        self.tempo_115s.clear()
        '-----------------------end------------'
        
        self.Cavst_plot.append(Ca_inicial)
        self.Cbvst_plot.append(Cb_inicial)
        self.Tcvst_plot.append(Tc_inicial)
        self.Trvst_plot.append(Tr_inicial)
        self.Tr0vst_plot.append(Tr0)
        self.Ca0vst_plot.append(Ca0)
        self.Qcvst_plot.append(Qc)
        self.qrvst_plot.append(qr)
        self.tempo_plot.append(t_inicial)
        
        '--------------------start--------------'
        #Coleta de dados a cada 0,050 minutos (3 segundos)
        if (count % 50 == 0):
            self.Cavst_20s.append(Ca_inicial)
            self.Cbvst_20s.append(Cb_inicial)
            self.Tcvst_20s.append(Tc_inicial)
            self.Trvst_20s.append(Tr_inicial)
            self.Tr0vst_20s.append(Tr0)
            self.Ca0vst_20s.append(Ca0)
            self.Qcvst_20s.append(Qc)
            self.qrvst_20s.append(qr)
            self.tempo_20s.append(t_inicial)
        
        #Coleta de dados a cada 1,25 minutos (75 segundos)
        if (count % 100 == 0):
            self.Cavst_75s.append(Ca_inicial)
            self.Cbvst_75s.append(Cb_inicial)
            self.Tcvst_75s.append(Tc_inicial)
            self.Trvst_75s.append(Tr_inicial)
            self.Tr0vst_75s.append(Tr0)
            self.Ca0vst_75s.append(Ca0)
            self.Qcvst_75s.append(Qc)
            self.qrvst_75s.append(qr)
            self.tempo_75s.append(t_inicial)
        
        #Coleta de dados a cada 1,583 minutos (94,98 segundos)
        if (count % 150 == 0):     
            self.Cavst_95s.append(Ca_inicial)
            self.Cbvst_95s.append(Cb_inicial)
            self.Tcvst_95s.append(Tc_inicial)
            self.Trvst_95s.append(Tr_inicial)
            self.Tr0vst_95s.append(Tr0)
            self.Ca0vst_95s.append(Ca0)
            self.Qcvst_95s.append(Qc)
            self.qrvst_95s.append(qr)
            self.tempo_95s.append(t_inicial)
            
        #Coleta de dados a cada 1,916 minutos (114,96 segundos)
        if (count % 200 == 0):    
            
            self.Cavst_115s.append(Ca_inicial)
            self.Cbvst_115s.append(Cb_inicial)
            self.Tcvst_115s.append(Tc_inicial)
            self.Trvst_115s.append(Tr_inicial)
            self.Tr0vst_115s.append(Tr0)
            self.Ca0vst_115s.append(Ca0)
            self.Qcvst_115s.append(Qc)
            self.qrvst_115s.append(qr)
            self.tempo_115s.append(t_inicial)
        '------------------------end----------'
        
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
        t = round(t,3)
        
        count = count + 1
        
        self.Cavst_plot.append(Ca_n1)
        self.Cbvst_plot.append(Cb_n1)
        self.Tcvst_plot.append(Tc_n1)
        self.Trvst_plot.append(Tr_n1)
        self.Tr0vst_plot.append(Tr0)
        self.Ca0vst_plot.append(Ca0)
        self.Qcvst_plot.append(Qc)
        self.qrvst_plot.append(qr)
        self.tempo_plot.append(t)
        self.contador.append(count)
        
        '--------------------------------start-----------------------'
        #Coleta de dados a cada 0,333 minutos (19,98 segundos)
        if (count % 50 == 0):
            self.Cavst_20s.append(Ca_n1)
            self.Cbvst_20s.append(Cb_n1)
            self.Tcvst_20s.append(Tc_n1)
            self.Trvst_20s.append(Tr_n1)
            self.Tr0vst_20s.append(Tr0)
            self.Ca0vst_20s.append(Ca0)
            self.Qcvst_20s.append(Qc)
            self.qrvst_20s.append(qr)
            self.tempo_20s.append(t)
        
        #Coleta de dados a cada 1,25 minutos (75 segundos)
        if (count % 100 == 0):
            
            self.Cavst_75s.append(Ca_n1)
            self.Cbvst_75s.append(Cb_n1)
            self.Tcvst_75s.append(Tc_n1)
            self.Trvst_75s.append(Tr_n1)
            self.Tr0vst_75s.append(Tr0)
            self.Ca0vst_75s.append(Ca0)
            self.Qcvst_75s.append(Qc)
            self.qrvst_75s.append(qr)
            self.tempo_75s.append(t)
        
        #Coleta de dados a cada 1,583 minutos (94,98 segundos)
        if (count % 150 == 0):
            
            self.Cavst_95s.append(Ca_n1)
            self.Cbvst_95s.append(Cb_n1)
            self.Tcvst_95s.append(Tc_n1)
            self.Trvst_95s.append(Tr_n1)
            self.Tr0vst_95s.append(Tr0)
            self.Ca0vst_95s.append(Ca0)
            self.Qcvst_95s.append(Qc)
            self.qrvst_95s.append(qr)
            self.tempo_95s.append(t)
            
        #Coleta de dados a cada 1,916 minutos (114,96 segundos)
        if (count % 200 == 0):
            
            self.Cavst_115s.append(Ca_n1)
            self.Cbvst_115s.append(Cb_n1)
            self.Tcvst_115s.append(Tc_n1)
            self.Trvst_115s.append(Tr_n1)
            self.Tr0vst_115s.append(Tr0)
            self.Ca0vst_115s.append(Ca0)
            self.Qcvst_115s.append(Qc)
            self.qrvst_115s.append(qr)
            self.tempo_115s.append(t)
        '-------------------------------------end---------------------'
        
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
            t = round(t,3)
            
            count = count + 1
            
            self.Cavst_plot.append(Ca_n1)
            self.Cbvst_plot.append(Cb_n1)
            self.Tcvst_plot.append(Tc_n1)
            self.Trvst_plot.append(Tr_n1)
            self.Tr0vst_plot.append(Tr0)
            self.Ca0vst_plot.append(Ca0)
            self.Qcvst_plot.append(Qc)
            self.qrvst_plot.append(qr)
            self.tempo_plot.append(t)
            self.contador.append(count)
            
            '--------------------------------start-----------------------'
            #Coleta de dados a cada 0,333 minutos (19,98 segundos)
            #Coleta de dados a cada 0,333 minutos (19,98 segundos)
            if (count % 50 == 0):
                self.Cavst_20s.append(Ca_n1)
                self.Cbvst_20s.append(Cb_n1)
                self.Tcvst_20s.append(Tc_n1)
                self.Trvst_20s.append(Tr_n1)
                self.Tr0vst_20s.append(Tr0)
                self.Ca0vst_20s.append(Ca0)
                self.Qcvst_20s.append(Qc)
                self.qrvst_20s.append(qr)
                self.tempo_20s.append(t)
            
            #Coleta de dados a cada 1,25 minutos (75 segundos)
            if (count % 100 == 0):
                
                self.Cavst_75s.append(Ca_n1)
                self.Cbvst_75s.append(Cb_n1)
                self.Tcvst_75s.append(Tc_n1)
                self.Trvst_75s.append(Tr_n1)
                self.Tr0vst_75s.append(Tr0)
                self.Ca0vst_75s.append(Ca0)
                self.Qcvst_75s.append(Qc)
                self.qrvst_75s.append(qr)
                self.tempo_75s.append(t)
            
            #Coleta de dados a cada 1,583 minutos (94,98 segundos)
            if (count % 150 == 0):
                
                self.Cavst_95s.append(Ca_n1)
                self.Cbvst_95s.append(Cb_n1)
                self.Tcvst_95s.append(Tc_n1)
                self.Trvst_95s.append(Tr_n1)
                self.Tr0vst_95s.append(Tr0)
                self.Ca0vst_95s.append(Ca0)
                self.Qcvst_95s.append(Qc)
                self.qrvst_95s.append(qr)
                self.tempo_95s.append(t)
                
            #Coleta de dados a cada 1,916 minutos (114,96 segundos)
            if (count % 200 == 0):
                
                self.Cavst_115s.append(Ca_n1)
                self.Cbvst_115s.append(Cb_n1)
                self.Tcvst_115s.append(Tc_n1)
                self.Trvst_115s.append(Tr_n1)
                self.Tr0vst_115s.append(Tr0)
                self.Ca0vst_115s.append(Ca0)
                self.Qcvst_115s.append(Qc)
                self.qrvst_115s.append(qr)
                self.tempo_115s.append(t)
            '-------------------------------------end---------------------'
            
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
    
    def __RK__perturbacoes_unitario(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,t_inicial,t_final,H,qr):
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
        #t_inicial = tempo inicial de simulação (min)
        #t_final = tempo final de simulação (min)
        #H: passo do método iterativo (min)
        #qr: vazão volumétrica da alimentação de reagentes
        #Tr = temperatura do reator (K)
        
        t = t_inicial
        
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
        
        "atualização dos parâmetros"
        Ca_n = Ca_n1
        Cb_n = Cb_n1
        Tr_n = Tr_n1
        Tc_n = Tc_n1
        
        while (t <= t_final):
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
            
    "Bloco para armazenamento dos dados gerados nas simulações totais do reator CSTR para construção do SoftSensor"
    
    #inserir dados de simulações totais 
    def inserir_dados(self,nome_tabela):
        try:
            connection = psycopg2.connect(
                dbname="TCC",
                user="postgres",
                password="pesquisa00",
                host="localhost",
                port="5432"
            )
            cursor = connection.cursor()

            insert_query = """
            INSERT INTO """+str(nome_tabela)+""" ("CA (kmol/m3)", "CB (kmol/m3)", "Tr (K)", "Tc (K)", "t (min)","Tr0 (K)","CA0 (kmol/m3)","Qc (kJ/min)","qr (m3/min)")
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s);
            """

            for i in range(len(self.Cavst_plot)):
                
                ca = self.Cavst_plot[i]
                cb = self.Cbvst_plot[i]
                tr = self.Trvst_plot[i]
                tc = self.Tcvst_plot[i]
                tr0 = self.Tr0vst_plot[i]
                ca0 = self.Ca0vst_plot[i]
                qc = self.Qcvst_plot[i]
                qr = self.qrvst_plot[i]
                t = self.tempo_plot[i]

                cursor.execute(insert_query, (ca, cb, tr, tc, t, tr0, ca0, qc, qr))

            connection.commit()
            cursor.close()
        except (Exception, psycopg2.Error) as error:
            print(f"Erro ao inserir dados: {error}")
        finally:
            if connection:
                connection.close()
    
    #Inserir dados de tempo de amostragem de 20 segundos
    def inserir_dados_20s(self,nome_tabela):
        try:
            connection = psycopg2.connect(
                dbname="TCC",
                user="postgres",
                password="pesquisa00",
                host="localhost",
                port="5432"
            )
            cursor = connection.cursor()

            insert_query = """
            INSERT INTO """+str(nome_tabela)+""" ("CA (kmol/m3)", "CB (kmol/m3)", "Tr (K)", "Tc (K)", "t (min)","Tr0 (K)","CA0 (kmol/m3)","Qc (kJ/min)","qr (m3/min)")
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s);
            """

            for i in range(len(self.Cavst_20s)):
                
                ca = self.Cavst_20s[i]
                cb = self.Cbvst_20s[i]
                tr = self.Trvst_20s[i]
                tc = self.Tcvst_20s[i]
                tr0 = self.Tr0vst_20s[i]
                ca0 = self.Ca0vst_20s[i]
                qc = self.Qcvst_20s[i]
                qr = self.qrvst_20s[i]
                t = self.tempo_20s[i]

                cursor.execute(insert_query, (ca, cb, tr, tc, t, tr0, ca0, qc, qr))

            connection.commit()
            cursor.close()
        except (Exception, psycopg2.Error) as error:
            print(f"Erro ao inserir dados: {error}")
        finally:
            if connection:
                connection.close()
    #Inserir dados de tempo de amostragem de 75 segundos
    def inserir_dados_75s(self,nome_tabela):
        try:
            connection = psycopg2.connect(
                dbname="TCC",
                user="postgres",
                password="pesquisa00",
                host="localhost",
                port="5432"
            )
            cursor = connection.cursor()

            insert_query = """
            INSERT INTO """+str(nome_tabela)+""" ("CA (kmol/m3)", "CB (kmol/m3)", "Tr (K)", "Tc (K)", "t (min)","Tr0 (K)","CA0 (kmol/m3)","Qc (kJ/min)","qr (m3/min)")
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s);
            """

            for i in range(len(self.Cavst_75s)):
                
                ca = self.Cavst_75s[i]
                cb = self.Cbvst_75s[i]
                tr = self.Trvst_75s[i]
                tc = self.Tcvst_75s[i]
                tr0 = self.Tr0vst_75s[i]
                ca0 = self.Ca0vst_75s[i]
                qc = self.Qcvst_75s[i]
                qr = self.qrvst_75s[i]
                t = self.tempo_75s[i]

                cursor.execute(insert_query, (ca, cb, tr, tc, t, tr0, ca0, qc, qr))

            connection.commit()
            cursor.close()
        except (Exception, psycopg2.Error) as error:
            print(f"Erro ao inserir dados: {error}")
        finally:
            if connection:
                connection.close()
    #inserir dados de 95 segundos
    def inserir_dados_95s(self,nome_tabela):
        try:
            connection = psycopg2.connect(
                dbname="TCC",
                user="postgres",
                password="pesquisa00",
                host="localhost",
                port="5432"
            )
            cursor = connection.cursor()

            insert_query = """
            INSERT INTO """+str(nome_tabela)+""" ("CA (kmol/m3)", "CB (kmol/m3)", "Tr (K)", "Tc (K)", "t (min)","Tr0 (K)","CA0 (kmol/m3)","Qc (kJ/min)","qr (m3/min)")
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s);
            """

            for i in range(len(self.Cavst_95s)):
                
                ca = self.Cavst_95s[i]
                cb = self.Cbvst_95s[i]
                tr = self.Trvst_95s[i]
                tc = self.Tcvst_95s[i]
                tr0 = self.Tr0vst_95s[i]
                ca0 = self.Ca0vst_95s[i]
                qc = self.Qcvst_95s[i]
                qr = self.qrvst_95s[i]
                t = self.tempo_95s[i]

                cursor.execute(insert_query, (ca, cb, tr, tc, t, tr0, ca0, qc, qr))

            connection.commit()
            cursor.close()
        except (Exception, psycopg2.Error) as error:
            print(f"Erro ao inserir dados: {error}")
        finally:
            if connection:
                connection.close()
    #inserir dados de 115 segundos
    def inserir_dados_115s(self,nome_tabela):
        try:
            connection = psycopg2.connect(
                dbname="TCC",
                user="postgres",
                password="pesquisa00",
                host="localhost",
                port="5432"
            )
            cursor = connection.cursor()

            insert_query = """
            INSERT INTO """+str(nome_tabela)+""" ("CA (kmol/m3)", "CB (kmol/m3)", "Tr (K)", "Tc (K)", "t (min)","Tr0 (K)","CA0 (kmol/m3)","Qc (kJ/min)","qr (m3/min)")
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s);
            """

            for i in range(len(self.Cavst_115s)):
                
                ca = self.Cavst_115s[i]
                cb = self.Cbvst_115s[i]
                tr = self.Trvst_115s[i]
                tc = self.Tcvst_115s[i]
                tr0 = self.Tr0vst_115s[i]
                ca0 = self.Ca0vst_115s[i]
                qc = self.Qcvst_115s[i]
                qr = self.qrvst_115s[i]
                t = self.tempo_115s[i]

                cursor.execute(insert_query, (ca, cb, tr, tc, t, tr0, ca0, qc, qr))

            connection.commit()
            cursor.close()
        except (Exception, psycopg2.Error) as error:
            print(f"Erro ao inserir dados: {error}")
        finally:
            if connection:
                connection.close()
    
    """Bloco que faz tudo"""
    
    "O bloco Simu_temporalCSTR realizada a simulação de um estado estacionário a partir de condições iniciais definidas"
    
    def Simu_temporalCSTR(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr,
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
        print('Tempo de simulação: '+str(self.tempo_plot[-1])+' min.')
        print('-------------------------------------------------------')
        
        #plotagem do gráfico
        "Plots dos gráficos de concentração"
        fig, graf = plt.subplots()
        plt.figure(figsize=(18,18))
        #subplot da concentração de A
        graf.plot(self.tempo_plot,self.Cavst_plot,label='CA (kmol/m³)')
        #subplot da concentração de B
        graf.plot(self.tempo_plot,self.Cbvst_plot,label='CB (kmol/m³)')
        graf.set_xlabel('tempo (min)')
        graf.set_ylabel('concentração molar (kmol/m³)')
        graf.set_title('concentração em função do tempo (-2,58% Tr0)')
        graf.legend()
        
        "Plots dos gráficos de temperatura"
        fig, graf1 = plt.subplots()
        plt.figure(figsize=(18,18))
        #subplot de valores de temperatura do reator Tr
        graf1.plot(self.tempo_plot,self.Trvst_plot,label='Tr (K)')
        #subplot de valores de temperatura de saída do fluido refrigerante Tc
        graf1.plot(self.tempo_plot,self.Tcvst_plot,label='Tc (K)')
        graf1.set_xlabel('tempo (min)')
        graf1.set_ylabel('Temperatura (K)')
        graf1.set_title('Temperatura em função do tempo (-2,58% Tr0)')
        graf1.legend()
        
        "Salvar os dados de simulação em planilha"
        #dataframe com dados das variaveis em função do tempo
        planilha_variaveis_vs_t = pd.DataFrame(data = {'CA (kmol/m3)':self.Cavst_plot,
                                             'CB (kmol/m3)':self.Cbvst_plot,
                                             'Tr (K)':self.Trvst_plot,
                                             'Tc (K)':self.Tcvst_plot,
                                             't (min)':self.tempo_plot})
        planilha_variaveis_vs_t.to_csv(arquivo_vs_t,index = False)
        
    
    def Simu_superficieCB_unitario(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,criterio_erro,H,
                              qr,liminf_Qc,limsup_Qc,H_Qc,arquivo):
        "Este bloco é responsável pela obtenção de todos os estados estacionários dentro dos intervalos de operação de"
        "qr e Qc. A partir dele são construídas as superfícies de Concentração e temperatura em função de qr e Qc"
        
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
    
    def Simu_perturbacoes_unitario(self,Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,H,qr,
                          arquivo_vs_t,var_pert,porc_pert):
        
        #var_pert = variável a ser perturbada
        #porc_pert = porcentagem de perturbação
        
        porc_pert = porc_pert/100
        
        "Perturbação na variável escolhida (-10% e +10%)"
        
        self.Cavst_simuperturbacoesunitario.clear()
        self.Cbvst_simuperturbacoesunitario.clear()
        self.Trvst_simuperturbacoesunitario.clear()
        self.Tcvst_simuperturbacoesunitario.clear()
        self.tempo_simuperturbacoesunitario.clear()
        
        #Realização da simulação nas condições ideais de operação
        t_inicial = 0
        t_final = 5
        self.__RK__perturbacoes_unitario(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,t_inicial,t_final,H,qr)
        
        for i in self.Cavst_plot:
            self.Cavst_simuperturbacoesunitario.append(i)    
        for i in self.Cbvst_plot:
            self.Cbvst_simuperturbacoesunitario.append(i)
        for i in self.Trvst_plot:
            self.Trvst_simuperturbacoesunitario.append(i)
        for i in self.Tcvst_plot:
            self.Tcvst_simuperturbacoesunitario.append(i)
        for i in self.tempo_plot:
            self.tempo_simuperturbacoesunitario.append(i)
        
        Ca_inicial = self.Cavst_simuperturbacoesunitario[-1]
        Cb_inicial = self.Cbvst_simuperturbacoesunitario[-1]
        Tr_inicial = self.Trvst_simuperturbacoesunitario[-1]
        Tc_inicial = self.Tcvst_simuperturbacoesunitario[-1]
        

        
        #Realização da perturbação nas condições ideais de operação (-10%)
        
        if var_pert == 0:
            #pertubação em qr
            qr = qr + (qr*porc_pert)
        elif var_pert == 1:
            #perturbação em Qc
            Qc = Qc + (Qc*porc_pert)
        elif var_pert == 2:
            #perturbação em Ca0
            Ca0 = Ca0 + (Ca0*porc_pert)
        elif var_pert == 3:
            #perturbação em Tr0
            Tr0 = Tr0 + (Tr0*porc_pert)
        
        t_inicial = t_final
        t_final = 60
        
        self.__RK__perturbacoes_unitario(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,t_inicial,t_final,H,qr)
        
        for i in self.Cavst_plot:
            self.Cavst_simuperturbacoesunitario.append(i)    
        for i in self.Cbvst_plot:
            self.Cbvst_simuperturbacoesunitario.append(i)
        for i in self.Trvst_plot:
            self.Trvst_simuperturbacoesunitario.append(i)
        for i in self.Tcvst_plot:
            self.Tcvst_simuperturbacoesunitario.append(i)
        for i in self.tempo_plot:
            self.tempo_simuperturbacoesunitario.append(i)
        
        
        print('-------------------------------------------------------')
        print('Dados da iteração após perturbação (-10%)')
        print('Cb (concentração de b no reator (kmol/m³)): {}'.format(self.Cbvst_simuperturbacoesunitario[-1]))
        print('Ca (concentração de a no reator (kmol/m³)): {}'.format(self.Cavst_simuperturbacoesunitario[-1]))
        print('Tr (temperatura do reator (K)): {}'.format(self.Trvst_simuperturbacoesunitario[-1]))
        print('Tc (temperatura do fluido refrigerante (K)): {}'.format(self.Tcvst_simuperturbacoesunitario[-1]))
        print('Tr0 (temperatura da corrente de alimentação (K)): {}'.format(Tr0))
        print('Ca0 (concentração de A na corrente de alimentação (kmol/m³)): {}'.format(Ca0))
        print('Qc (calor trocado com o fluido refrigerante antes da entrada no reator (kJ/min)): {}'.format(Qc))
        print('qr (vazão volumétrica da corrente de alimentação de reagentes (m³/min)): {}'.format(qr))
        print('Tempo de simulação: '+str(self.tempo_simuperturbacoesunitario[-1])+' min.')
        print('-------------------------------------------------------')
        
        #plotagem do gráfico
        "Plots dos gráficos de concentração"
        fig, graf = plt.subplots()
        plt.figure(figsize=(18,5))
        #subplot da concentração de A
        graf.plot(self.tempo_simuperturbacoesunitario,self.Cavst_simuperturbacoesunitario,label='CA (kmol/m³)')
        #subplot da concentração de B
        graf.plot(self.tempo_simuperturbacoesunitario,self.Cbvst_simuperturbacoesunitario,label='CB (kmol/m³)')
        graf.set_xlabel('tempo (min)')
        graf.set_ylabel('concentração molar (kmol/m³)')
        graf.set_title('concentração em função do tempo (-2,58% Tr0)')
        graf.legend()
        
        "Plots dos gráficos de temperatura"
        fig, graf1 = plt.subplots()
        plt.figure(figsize=(18,5))
        #subplot de valores de temperatura do reator Tr
        graf1.plot(self.tempo_simuperturbacoesunitario,self.Trvst_simuperturbacoesunitario,label='Tr (K)')
        #subplot de valores de temperatura de saída do fluido refrigerante Tc
        graf1.plot(self.tempo_simuperturbacoesunitario,self.Tcvst_simuperturbacoesunitario,label='Tc (K)')
        graf1.set_xlabel('tempo (min)')
        graf1.set_ylabel('Temperatura (K)')
        graf1.set_title('Temperatura em função do tempo (-2,58% Tr0)')
        graf1.legend()
        
        "Salvar os dados de simulação em planilha"
        #dataframe com dados das variaveis em função do tempo
        planilha_variaveis_vs_t = pd.DataFrame(data = {'CA (kmol/m3)':self.Cavst_simuperturbacoesunitario,
                                             'CB (kmol/m3)':self.Cbvst_simuperturbacoesunitario,
                                             'Tr (K)':self.Trvst_simuperturbacoesunitario,
                                             'Tc (K)':self.Tcvst_simuperturbacoesunitario,
                                             't (min)':self.tempo_simuperturbacoesunitario})
        planilha_variaveis_vs_t.to_csv(arquivo_vs_t,index = False)
    
    def gera_dados_softsensor(self,Ca_estacionario,Cb_estacionario,Tr_estacionario,Tc_estacionario,
                              Ca0_estacionario,Tr0_estacionario,Qc_estacionario,criterio_erro,H,qr_estacionario,
                              qntd_simulacoes,nome_tabela_20s,nome_tabela_75s,nome_tabela_95s,nome_tabela_115s,
                              nome_tabela):
        
        #var_pert = variável a ser perturbada
        #porc_pert = porcentagem de perturbação
        #qntd_simulacoes = quantidade de simulações de perturbações que devem se geradas
        #nome_tabela = o nome da tabela armazenada no database "TCC" que armazenará os dados da simulação
        #view_graf = visualizar os gráficos (1 = Não, 0 = Sim)
        #planilha_dados = nome da planilha contendo os dados salvos em CSV
        
        #Definindo os contadores de perturbações
        
        t_inicial = 0
        
        contador = 0
        
        cont_pert_Tr0 = 0
        cont_pert_Ca0 = 0
        cont_pert_Qc = 0
        cont_pert_qr = 0
        
        self.Cavst_simuperturbacoesunitario.clear()
        self.Cbvst_simuperturbacoesunitario.clear()
        self.Trvst_simuperturbacoesunitario.clear()
        self.Tcvst_simuperturbacoesunitario.clear()
        self.tempo_simuperturbacoesunitario.clear()
        
        cont = 0
        
        Ca_inicial = Ca_estacionario
        Cb_inicial = Cb_estacionario
        Tr_inicial = Tr_estacionario
        Tc_inicial = Tc_estacionario
        
        "Realização da simulação em torno das condições ideais de operação"
        
        while (cont <= qntd_simulacoes):
            
            "Definindo o ponto inicial da simulação como sendo o estado estacionário ideal"
            
            Ca0 = Ca0_estacionario
            Tr0 = Tr0_estacionario
            Qc = Qc_estacionario
            qr = qr_estacionario
            
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
                cont_pert_Tr0 = cont_pert_Tr0 + 1
            if ativacao_Ca0 == 1:
                boolean_Ca0 = True
                cont_pert_Ca0 = cont_pert_Ca0 + 1
            if ativacao_Qc == 1:
                boolean_Qc = True
                cont_pert_Qc = cont_pert_Qc + 1
            if ativacao_qr == 1:
                boolean_qr = True
                cont_pert_qr = cont_pert_qr + 1
            
            if boolean_Tr0:
                porc_perturbacao_Tr0 = random.uniform(-2.58,2.58)
                Tr0 = Tr0 + Tr0*((porc_perturbacao_Tr0)/(100))
                Tr0 = round(Tr0,4)
            if boolean_Ca0:
                porc_perturbacao_Ca0 = random.uniform(-10,10)
                Ca0 = Ca0 + Ca0*((porc_perturbacao_Ca0)/(100))
                Ca0 = round(Ca0,4)
            if boolean_Qc:
                porc_perturbacao_Qc = random.uniform(-10,10)
                Qc = Qc + Qc*((porc_perturbacao_Qc)/(100))
                Qc = round(Qc,4)
            if boolean_qr:
                porc_perturbacao_qr = random.uniform(-10,10)
                qr = qr + qr*((porc_perturbacao_qr)/(100))
                qr = round(qr,6)
            
            self.__RK__gera_dados_softsensor(Ca_inicial,Cb_inicial,Tr_inicial,Tc_inicial,Ca0,Tr0,Qc,criterio_erro,H,qr,
                                             t_inicial,contador)
            
            contador = self.contador[-1]
            
            
            for i in self.Cavst_plot:
                self.Cavst_simuperturbacoesunitario.append(i)
            for i in self.Cbvst_plot:
                self.Cbvst_simuperturbacoesunitario.append(i)
            for i in self.Trvst_plot:
                self.Trvst_simuperturbacoesunitario.append(i)
            for i in self.Tcvst_plot:
                self.Tcvst_simuperturbacoesunitario.append(i)
            for i in self.tempo_plot:
                self.tempo_simuperturbacoesunitario.append(i)
                
            Ca_inicial = self.Cavst_simuperturbacoesunitario[-1]
            Cb_inicial = self.Cbvst_simuperturbacoesunitario[-1]
            Tr_inicial = self.Trvst_simuperturbacoesunitario[-1]
            Tc_inicial = self.Tcvst_simuperturbacoesunitario[-1]
            
            
            print('-------------------------------------------------------')
            print('DADOS DA ITERAÇÃO: '+str(cont))
            print('Cb (concentração de b no reator (kmol/m³)): {}'.format(self.Cbvst_simuperturbacoesunitario[-1]))
            print('Ca (concentração de a no reator (kmol/m³)): {}'.format(self.Cavst_simuperturbacoesunitario[-1]))
            print('Tr (temperatura do reator (K)): {}'.format(self.Trvst_simuperturbacoesunitario[-1]))
            print('Tc (temperatura do fluido refrigerante (K)): {}'.format(self.Tcvst_simuperturbacoesunitario[-1]))
            print('Tr0 (temperatura da corrente de alimentação (K)): {}'.format(Tr0))
            print('Ca0 (concentração de A na corrente de alimentação (kmol/m³)): {}'.format(Ca0))
            print('Qc (calor trocado com o fluido refrigerante antes da entrada no reator (kJ/min)): {}'.format(Qc))
            print('qr (vazão volumétrica da corrente de alimentação de reagentes (m³/min)): {}'.format(qr))
            print('Tempo de simulação: '+str(self.tempo_simuperturbacoesunitario[-1])+' min.')
            print('DADOS DE PERTURBAÇÕES')
            print('Nº de perturbações em Tr0: '+str(cont_pert_Tr0))
            print('Nº de perturbações em Ca0: '+str(cont_pert_Ca0))
            print('Nº de perturbações em Qc: '+str(cont_pert_Qc))
            print('Nº de perturbações em qr: '+str(cont_pert_qr))
            print('-------------------------------------------------------')
            
            "inserir os dados na base de dados"
            
            #inserir dados de simulações totais
            #self.inserir_dados(nome_tabela)
            
            #Dados de tempo de amostragem de 20 segundos
            self.inserir_dados_20s(nome_tabela_20s) #Insere os dados na tabela "nome_tabela"
            #Dados de tempo de amostragem de 75 segundos
            self.inserir_dados_75s(nome_tabela_75s)
            #Dados de tempo de amostragem de 95 segundos
            self.inserir_dados_95s(nome_tabela_95s)
            #Dados de tempo de amostragem de 115 segundos
            self.inserir_dados_115s(nome_tabela_115s)
            
            self.Cavst_simuperturbacoesunitario.clear()
            self.Cbvst_simuperturbacoesunitario.clear()
            self.Trvst_simuperturbacoesunitario.clear()
            self.Tcvst_simuperturbacoesunitario.clear()
            self.tempo_simuperturbacoesunitario.clear()
            
            t_inicial = self.tempo_plot[-1] + H
            
            cont = cont + 1      
        
        #Conectando com o banco de dados PostegreSQL    
        connection = psycopg2.connect(
            dbname="TCC",
            user="postgres",
            password="pesquisa00",
            host="localhost",
            port="5432"
        )
        
        "DADOS DE TEMPO DE AMOSTRAGEM DE 20 SEGUNDOS"
        
        # Executar a consulta para selecionar e ordenar todas as colunas da tabela "simulacoes_totais"
        query = "SELECT * FROM simulacoes_20s ORDER BY \"t (min)\" ASC;"
        
        # Ler os dados usando pandas
        df = pd.read_sql_query(query, connection)
        
        # Fechar a conexão com o banco de dados
        connection.close()
        
        # Salvar os dados em um arquivo CSV
        df.to_csv('Simulacoes_3s_acuracia.csv', index=False)      
          
        """
        Ca_plot = df['CA (kmol/m3)'] #armazena dados de concentração de A
        Cb_plot = df['CB (kmol/m3)'] #armazena dados de concentração de B
        Tr_plot = df['Tr (K)'] #armazena dados de temperatura do reator 
        Tc_plot = df['Tc (K)'] #armazena dados de temperatura da jaqueta térmica
        Tr0_plot = df['Tr0 (K)'] #armazena dados de temperatura da corrente de alimentação
        Ca0_plot = df['CA0 (kmol/m3)'] #armazena dados da concentração do reagente A
        Qc_plot = df['Qc (kJ/min)'] #armazena dados de Qc
        qr_plot = df['qr (m3/min)'] #armazena dados de vazão de alimentação de reagentes
        t_plot = df['t (min)'] #armazena dados de tempo
        
        #tr0 = self.Tr0vst_plot[i]
        #ca0 = self.Ca0vst_plot[i]
        #qc = self.Qcvst_plot[i]
        #qr = self.qrvst_plot[i]
        #t = self.tempo_plot[i]
        
        #plotagem do gráfico
        "Plots dos gráficos de concentração"
        dpi = 300  # Resolução em pontos por polegada (DPI)
        fig_size = (18, 10)  # Tamanho da figura em polegadas (largura, altura)
       
        #Variável controlada
        fig, graf = plt.subplots(figsize = fig_size, dpi=dpi)
        graf.plot(t_plot,Cb_plot)
        graf.set_xlabel('tempo (min)')
        graf.set_ylabel('concentração molar (kmol/m³)')
        graf.set_title('VARIÁVEL CONTROLADA - 20 segundos')
        graf.legend()
        plt.savefig("concentracao_b_20s.png", dpi=dpi, bbox_inches="tight")
        
        #Variável distúrbio
        fig, graf2 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf2.plot(t_plot,Ca0_plot)
        graf2.set_xlabel('tempo (min)')
        graf2.set_ylabel('concentração molar (kmol/m³)')
        graf2.set_title('VARIÁVEL DE DISTÚRBIO - 20 segundos')
        graf2.legend()
        plt.savefig("concentracao_ca0_20s.png", dpi=dpi, bbox_inches="tight")
        
        fig, graf3 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf3.plot(t_plot,Tr0_plot)
        graf3.set_xlabel('tempo (min)')
        graf3.set_ylabel('Temperatura (K)')
        graf3.set_title('VARIÁVEL DE DISTÚRBIO - 20 segundos')
        graf3.legend()
        plt.savefig("temperatura_Tr0_20s.png", dpi=dpi, bbox_inches="tight")
        
        #Variável manipulada
        fig, graf4 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf4.plot(t_plot,qr_plot)
        graf4.set_xlabel('tempo (min)')
        graf4.set_ylabel('Vazão (m³/min)')
        graf4.set_title('VARIÁVEL MANIPULADA - 20 segundos')
        graf4.legend()
        plt.savefig("vazao_qr_20s.png", dpi=dpi, bbox_inches="tight")
        
        fig, graf5 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf5.plot(t_plot,Qc_plot)
        graf5.set_xlabel('tempo (min)')
        graf5.set_ylabel('Taxa de troca térmica (kJ/min)')
        graf5.set_title('VARIÁVEL MANIPULADA - 20 segundos')
        graf5.legend()
        plt.savefig("trocatermica_Qc_20s.png", dpi=dpi, bbox_inches="tight")
        """
        
        #Conectando com o banco de dados PostegreSQL    
        connection = psycopg2.connect(
            dbname="TCC",
            user="postgres",
            password="pesquisa00",
            host="localhost",
            port="5432"
        )
        
        "DADOS DE TEMPO DE AMOSTRAGEM DE 75 SEGUNDOS"
        
        # Executar a consulta para selecionar e ordenar todas as colunas da tabela "simulacoes_totais"
        query = "SELECT * FROM simulacoes_75s ORDER BY \"t (min)\" ASC;"
        
        # Ler os dados usando pandas
        df1 = pd.read_sql_query(query, connection)
        
        # Fechar a conexão com o banco de dados
        connection.close()
        
        # Salvar os dados em um arquivo CSV
        df1.to_csv('Simulacoes_6s_acuracia.csv', index=False)   
        
        """
        
        Ca_plot = df1['CA (kmol/m3)'] #armazena dados de concentração de A
        Cb_plot = df1['CB (kmol/m3)'] #armazena dados de concentração de B
        Tr_plot = df1['Tr (K)'] #armazena dados de temperatura do reator 
        Tc_plot = df1['Tc (K)'] #armazena dados de temperatura da jaqueta térmica
        Tr0_plot = df1['Tr0 (K)'] #armazena dados de temperatura da corrente de alimentação
        Ca0_plot = df1['CA0 (kmol/m3)'] #armazena dados da concentração do reagente A
        Qc_plot = df1['Qc (kJ/min)'] #armazena dados de Qc
        qr_plot = df1['qr (m3/min)'] #armazena dados de vazão de alimentação de reagentes
        t_plot = df1['t (min)'] #armazena dados de tempo
        
        #tr0 = self.Tr0vst_plot[i]
        #ca0 = self.Ca0vst_plot[i]
        #qc = self.Qcvst_plot[i]
        #qr = self.qrvst_plot[i]
        #t = self.tempo_plot[i]
        
        #plotagem do gráfico
        "Plots dos gráficos de concentração"
        dpi = 300  # Resolução em pontos por polegada (DPI)
        fig_size = (18, 10)  # Tamanho da figura em polegadas (largura, altura)
       
        #Variável controlada
        fig, graf = plt.subplots(figsize = fig_size, dpi=dpi)
        graf.plot(t_plot,Cb_plot)
        graf.set_xlabel('tempo (min)')
        graf.set_ylabel('concentração molar (kmol/m³)')
        graf.set_title('VARIÁVEL CONTROLADA - 75 segundos')
        graf.legend()
        plt.savefig("concentracao_b_75s.png", dpi=dpi, bbox_inches="tight")
        
        #Variável distúrbio
        fig, graf2 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf2.plot(t_plot,Ca0_plot)
        graf2.set_xlabel('tempo (min)')
        graf2.set_ylabel('concentração molar (kmol/m³)')
        graf2.set_title('VARIÁVEL DE DISTÚRBIO - 75s')
        graf2.legend()
        plt.savefig("concentracao_ca0_75s.png", dpi=dpi, bbox_inches="tight")
        
        fig, graf3 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf3.plot(t_plot,Tr0_plot)
        graf3.set_xlabel('tempo (min)')
        graf3.set_ylabel('Temperatura (K)')
        graf3.set_title('VARIÁVEL DE DISTÚRBIO')
        graf3.legend()
        plt.savefig("temperatura_Tr0_75s.png", dpi=dpi, bbox_inches="tight")
        
        #Variável manipulada
        fig, graf4 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf4.plot(t_plot,qr_plot)
        graf4.set_xlabel('tempo (min)')
        graf4.set_ylabel('Vazão (m³/min)')
        graf4.set_title('VARIÁVEL MANIPULADA - 75 segundos')
        graf4.legend()
        plt.savefig("vazao_qr_75s.png", dpi=dpi, bbox_inches="tight")
        
        fig, graf5 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf5.plot(t_plot,Qc_plot)
        graf5.set_xlabel('tempo (min)')
        graf5.set_ylabel('Taxa de troca térmica (kJ/min)')
        graf5.set_title('VARIÁVEL MANIPULADA - 75 segundos')
        graf5.legend()
        plt.savefig("trocatermica_Qc_75s.png", dpi=dpi, bbox_inches="tight")
        """
        
        #Conectando com o banco de dados PostegreSQL    
        connection = psycopg2.connect(
            dbname="TCC",
            user="postgres",
            password="pesquisa00",
            host="localhost",
            port="5432"
        )
        
        "DADOS DE TEMPO DE AMOSTRAGEM DE 95 SEGUNDOS"
        
        # Executar a consulta para selecionar e ordenar todas as colunas da tabela "simulacoes_totais"
        query = "SELECT * FROM simulacoes_95s ORDER BY \"t (min)\" ASC;"
        
        # Ler os dados usando pandas
        df = pd.read_sql_query(query, connection)
        
        # Fechar a conexão com o banco de dados
        connection.close()
        
        # Salvar os dados em um arquivo CSV
        df.to_csv('Simulacoes_9s_acuracia.csv', index=False)    
        
        """
        
        Ca_plot = df['CA (kmol/m3)'] #armazena dados de concentração de A
        Cb_plot = df['CB (kmol/m3)'] #armazena dados de concentração de B
        Tr_plot = df['Tr (K)'] #armazena dados de temperatura do reator 
        Tc_plot = df['Tc (K)'] #armazena dados de temperatura da jaqueta térmica
        Tr0_plot = df['Tr0 (K)'] #armazena dados de temperatura da corrente de alimentação
        Ca0_plot = df['CA0 (kmol/m3)'] #armazena dados da concentração do reagente A
        Qc_plot = df['Qc (kJ/min)'] #armazena dados de Qc
        qr_plot = df['qr (m3/min)'] #armazena dados de vazão de alimentação de reagentes
        t_plot = df['t (min)'] #armazena dados de tempo
        
        #tr0 = self.Tr0vst_plot[i]
        #ca0 = self.Ca0vst_plot[i]
        #qc = self.Qcvst_plot[i]
        #qr = self.qrvst_plot[i]
        #t = self.tempo_plot[i]
        
        #plotagem do gráfico
        "Plots dos gráficos de concentração"
        dpi = 300  # Resolução em pontos por polegada (DPI)
        fig_size = (18, 10)  # Tamanho da figura em polegadas (largura, altura)
       
        #Variável controlada
        fig, graf = plt.subplots(figsize = fig_size, dpi=dpi)
        graf.plot(t_plot,Cb_plot)
        graf.set_xlabel('tempo (min)')
        graf.set_ylabel('concentração molar (kmol/m³)')
        graf.set_title('VARIÁVEL CONTROLADA - 95 segundos')
        graf.legend()
        plt.savefig("concentracao_b_95s.png", dpi=dpi, bbox_inches="tight")
        
        #Variável distúrbio
        fig, graf2 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf2.plot(t_plot,Ca0_plot)
        graf2.set_xlabel('tempo (min)')
        graf2.set_ylabel('concentração molar (kmol/m³)')
        graf2.set_title('VARIÁVEL DE DISTÚRBIO - 95 segundos')
        graf2.legend()
        plt.savefig("concentracao_ca0_95s.png", dpi=dpi, bbox_inches="tight")
        
        fig, graf3 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf3.plot(t_plot,Tr0_plot)
        graf3.set_xlabel('tempo (min)')
        graf3.set_ylabel('Temperatura (K)')
        graf3.set_title('VARIÁVEL DE DISTÚRBIO - 95 segundos')
        graf3.legend()
        plt.savefig("temperatura_Tr0_95s.png", dpi=dpi, bbox_inches="tight")
        
        #Variável manipulada
        fig, graf4 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf4.plot(t_plot,qr_plot)
        graf4.set_xlabel('tempo (min)')
        graf4.set_ylabel('Vazão (m³/min)')
        graf4.set_title('VARIÁVEL MANIPULADA - 95 segundos')
        graf4.legend()
        plt.savefig("vazao_qr_95s.png", dpi=dpi, bbox_inches="tight")
        
        fig, graf5 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf5.plot(t_plot,Qc_plot)
        graf5.set_xlabel('tempo (min)')
        graf5.set_ylabel('Taxa de troca térmica (kJ/min)')
        graf5.set_title('VARIÁVEL MANIPULADA - 95 segundos')
        graf5.legend()
        plt.savefig("trocatermica_Qc_95s.png", dpi=dpi, bbox_inches="tight")
        """
        #Conectando com o banco de dados PostegreSQL    
        connection = psycopg2.connect(
            dbname="TCC",
            user="postgres",
            password="pesquisa00",
            host="localhost",
            port="5432"
        )
        
        "DADOS DE TEMPO DE AMOSTRAGEM DE 115 SEGUNDOS"
        
        # Executar a consulta para selecionar e ordenar todas as colunas da tabela "simulacoes_totais"
        query = "SELECT * FROM simulacoes_115s ORDER BY \"t (min)\" ASC;"
        
        # Ler os dados usando pandas
        df = pd.read_sql_query(query, connection)
        
        # Fechar a conexão com o banco de dados
        connection.close()
        
        # Salvar os dados em um arquivo CSV
        df.to_csv('Simulacoes_12s_acuracia.csv', index=False)      
        """
        Ca_plot = df['CA (kmol/m3)'] #armazena dados de concentração de A
        Cb_plot = df['CB (kmol/m3)'] #armazena dados de concentração de B
        Tr_plot = df['Tr (K)'] #armazena dados de temperatura do reator 
        Tc_plot = df['Tc (K)'] #armazena dados de temperatura da jaqueta térmica
        Tr0_plot = df['Tr0 (K)'] #armazena dados de temperatura da corrente de alimentação
        Ca0_plot = df['CA0 (kmol/m3)'] #armazena dados da concentração do reagente A
        Qc_plot = df['Qc (kJ/min)'] #armazena dados de Qc
        qr_plot = df['qr (m3/min)'] #armazena dados de vazão de alimentação de reagentes
        t_plot = df['t (min)'] #armazena dados de tempo
        
        #tr0 = self.Tr0vst_plot[i]
        #ca0 = self.Ca0vst_plot[i]
        #qc = self.Qcvst_plot[i]
        #qr = self.qrvst_plot[i]
        #t = self.tempo_plot[i]
        
        #plotagem do gráfico
        "Plots dos gráficos de concentração"
        dpi = 300  # Resolução em pontos por polegada (DPI)
        fig_size = (18, 10)  # Tamanho da figura em polegadas (largura, altura)
       
        #Variável controlada
        fig, graf = plt.subplots(figsize = fig_size, dpi=dpi)
        graf.plot(t_plot,Cb_plot)
        graf.set_xlabel('tempo (min)')
        graf.set_ylabel('concentração molar (kmol/m³)')
        graf.set_title('VARIÁVEL CONTROLADA - 115 segundos')
        graf.legend()
        plt.savefig("concentracao_b_115s.png", dpi=dpi, bbox_inches="tight")
        
        #Variável distúrbio
        fig, graf2 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf2.plot(t_plot,Ca0_plot)
        graf2.set_xlabel('tempo (min)')
        graf2.set_ylabel('concentração molar (kmol/m³)')
        graf2.set_title('VARIÁVEL DE DISTÚRBIO - 115 segundos')
        graf2.legend()
        plt.savefig("concentracao_ca0_115s.png", dpi=dpi, bbox_inches="tight")
        
        fig, graf3 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf3.plot(t_plot,Tr0_plot)
        graf3.set_xlabel('tempo (min)')
        graf3.set_ylabel('Temperatura (K)')
        graf3.set_title('VARIÁVEL DE DISTÚRBIO - 115 segundos')
        graf3.legend()
        plt.savefig("temperatura_Tr0_115s.png", dpi=dpi, bbox_inches="tight")
        
        #Variável manipulada
        fig, graf4 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf4.plot(t_plot,qr_plot)
        graf4.set_xlabel('tempo (min)')
        graf4.set_ylabel('Vazão (m³/min)')
        graf4.set_title('VARIÁVEL MANIPULADA - 115 segundos')
        graf4.legend()
        plt.savefig("vazao_qr_115s.png", dpi=dpi, bbox_inches="tight")
        
        fig, graf5 = plt.subplots(figsize = fig_size, dpi=dpi)
        graf5.plot(t_plot,Qc_plot)
        graf5.set_xlabel('tempo (min)')
        graf5.set_ylabel('Taxa de troca térmica (kJ/min)')
        graf5.set_title('VARIÁVEL MANIPULADA - 115 segundos')
        graf5.legend()
        plt.savefig("trocatermica_Qc_115s.png", dpi=dpi, bbox_inches="tight")
        """
    "----------------------------------------------------------------------------------------------------"
    
