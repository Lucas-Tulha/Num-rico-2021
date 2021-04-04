#Módulos Necessários
import numpy as np
from math import exp, cos, sin
from metodos import exercicio1_1, exercicio1_2
#---------------------------------------------------------------------------------------------------
######################################### EXERCÍCIO 1 PARTE 1 ######################################

#Funções de X (solução explícita)
def x1(t):
    return exp(-t) * sin(t) + exp(-3*t) * cos(3*t)

def x2(t):
    return exp(-t) * cos(t) + exp(-3*t) * sin(3*t)

def x3(t):
    return -1 * exp(-t) * sin(t) + exp(-3*t) * cos(3*t)

def x4(t):
    return -1 * exp(-t) * cos(t) + exp(-3*t) * sin(3*t)

#Constantes
T0_1_1 = 0
Tf_1_1 = 2
lista_n = [20, 40, 80, 160, 320, 640]

A = np.array([[-2, -1, -1, -2], [1, -2, 2, -1], [-1, -2, -2, -1], [2, -1, 1, -2]])
x0_1_1= np.array([1, 1, 1, -1])

lista_x_exp = [x1, x2, x3, x4]
#---------------------------------------------------------------------------------------------------
######################################### EXERCÍCIO 1 PARTE 2 ######################################

#Derivada de X
def x_derivada(t, x):
     return 2*t + (x-t**2)**2

#Função X Explícita
def x_explicita(t):
    return t**2 + 1/(1 - t)

#Jacobiano
def jacobinano_1_2(x, t, h):
    return -2*h*(x - t**2)

#G(X, t, h)
#Monta polinômio que representa G(x, t, h)
def montar_G(h, t, x0):
    polinomio = []

    termo_constante = -2*h*t -h*(t**4) -x0
    polinomio.append(termo_constante)

    termo_vetorial = 1 + 2*h*(t**2)
    polinomio.append(termo_vetorial)

    termo_quadratico = -h
    polinomio.append(termo_quadratico)

    return polinomio

#Constantes
T0_1_2 = 1.1
Tf_1_2 = 3.0
n = 5000
x0_1_2 = -8.79
#---------------------------------------------------------------------------------------------------
#Execução dos Exerícios
exercicio1_1(T0_1_1, Tf_1_1, lista_n, x0_1_1, A, lista_x_exp)

exercicio1_2(T0_1_2, Tf_1_2, n, x0_1_2, x_explicita, montar_G)
