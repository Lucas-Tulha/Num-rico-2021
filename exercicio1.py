#Módulos Necessários
import numpy as np
from math import exp, cos, sin
from metodos import exercicio1_1
#---------------------------------------------------------------------------------------------------

#Funções de X (solução explícita)
def x1(t):
    return exp(-t) * sin(t) + exp(-3*t) * cos(3*t)

def x2(t):
    return exp(-t) * cos(t) + exp(-3*t) * sin(3*t)

def x3(t):
    return -1 * exp(-t) * sin(t) + exp(-3*t) * cos(3*t)

def x4(t):
    return -1* exp(-t) * cos(t) + exp(-3*t) * sin(3*t)
#---------------------------------------------------------------------------------------------------

#Constantes
T0 = 0
Tf = 2
lista_n = [20, 40, 80, 160, 320, 640]

A = np.array([[-2, -1, -1, -2], [1, -2, 2, 1], [-1, -2, -2, -1], [2, -1, 1, -2]])
x0 = np.array([1, 1, 1, -1])

lista_x_exp = [x1, x2, x3, x4]
#---------------------------------------------------------------------------------------------------


#Calculo de Ri



#Executar Função
exercicio1_1(T0, Tf, lista_n, x0, A, lista_x_exp)
