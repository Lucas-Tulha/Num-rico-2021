#Módulos Necessários
import numpy as np
from math import exp, cos, sin
from metodos import exercicio2_1, exercicio2_2, exercicio2_3, exercicio2_4
#---------------------------------------------------------------------------------------------------
######################################### EXERCÍCIO 2 PARTE 1 ######################################

#Taxas de Crescimento (Derivadas)
def derivada_x(x, y, alpha, lambida):
    return lambida * x - alpha * x * y

def derivada_y(x, y, beta, gama):
    return beta * x * y - gama * y

#Constantes
lambida = 2/3
alpha = 4/3
beta = 1.0
gama = 1.0
x0 = 1.5
y0 = 1.5
T0 = 0.0
Tf = 10.0
n = 7000
#---------------------------------------------------------------------------------------------------
######################################### EXERCÍCIO 2 PARTE 2 ######################################

#Monta G(U)
def montar_G_duplo(X, Y, X0, Y0, h):
    G = []
    G0 = X - h*((2/3)*X -(4/3)*X*Y) - X0
    G1 = Y - h*(X*Y -Y) - Y0
    G.append([G0])
    G.append([G1])

    return np.array(G)

#Derivadas das derivadas de x e de y
def derivada_f1_x(x, y):
    return (2/3) - (4/3) * y

def derivada_f1_y(x, y):
    return (-4/3) * x

def derivada_f2_x(x, y):
    return y

def derivada_f2_y(x, y):
    return x -1

#Monta Jacobiano
def montar_J(x, y, h):
    J = []

    J00 = 1 - h*(derivada_f1_x(x, y))
    J01 = -h*(derivada_f1_y(x, y))
    J.append([J00, J01])

    J10 = -h*(derivada_f2_x(x, y))
    J11 = 1 - h*(derivada_f2_y(x, y))
    J.append([J10, J11])

    return np.array(J)
#---------------------------------------------------------------------------------------------------
######################################### EXERCÍCIO 2 PARTE 3 ######################################

#Constantes
lista_n = [250, 500, 1000, 2000, 4000]
#---------------------------------------------------------------------------------------------------
#Execução dos Exerícios
#exercicio2_1(lambida, alpha, beta, gama, x0, y0, T0, Tf, n, derivada_x, derivada_y, plotar=True)

#exercicio2_2(x0, y0, T0, Tf, n, montar_G_duplo, montar_J, plotar=True)

#exercicio2_3(lambida, alpha, beta, gama, x0, y0, T0, Tf, lista_n, derivada_x, derivada_y, montar_G_duplo, montar_J)

exercicio2_4(T0, Tf, n, derivada_x, derivada_y, x0, y0, alpha, lambida, beta, gama)
