#Módulos Necessários
import numpy as np
from math import exp, cos, sin
from metodos import exercicio3
#---------------------------------------------------------------------------------------------------
######################################### EXERCÍCIO 3 PARTE 1 ######################################

#Taxas de Crescimento (Derivadas)
def derivada_x(x, y, z):
    return x*(1 - 0.001 * x - 0.001 * y - 0.015 * z)

def derivada_y(x, y, z):
    return y*(1 - 0.0015 * x - 0.001 * y - 0.001 * z)

def derivada_z(x, y, z, alpha):
    return z*(-1 + alpha * x + 0.0005 * y)


#Constantes
lista_alpha = [0.001, 0.002, 0.0033, 0.0036, 0.005, 0.0055]
x0 = 500
y0 = 500
z0 = 10

n = 7000
T0 = 0
lista_Tf = [100, 100, 500, 500, 2000, 2000]
#---------------------------------------------------------------------------------------------------

#Execução dos Exerícios
exercicio3(lista_alpha, x0, y0, z0, T0, lista_Tf, derivada_x, derivada_y, derivada_z, n)
