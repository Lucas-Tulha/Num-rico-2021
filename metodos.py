#Módulos Necessários
import numpy as np
from matplotlib import pyplot as plt
from math import exp, cos, sin
#---------------------------------------------------------------------------------------------------

#Funções

#Calcula o passo para este n
def calcular_passo(T0, Tf, n):
    return (Tf - T0) / n

#Aproximar valores
def aproximar(valor):
    return int(valor *100000000000) / 100000000000

#Calcular valores de T que serão abordados neste intervalo
def calcular_lista_t(T0, Tf, h):
    lista_t = []
    t = T0

    while t <= Tf:
        lista_t.append(t)
        t = aproximar(t + h)

    return lista_t
#---------------------------------------------------------------------------------------------------
######################## EXERCICO 1 PARTE 1 ##########################################################

#RUNGE-KUTTA, nessa parte vao ficar todos os metodos pedidos como de Runge-Kutta

#Calcular taxa de crescimento de X num dado ponto
#(específico para o exercício 1.1)
def calcular_k1_1(x, A):
    return np.matmul(A, x)

#Calcular próximo ponto com o desvio do método de Runge-Kutta de Ordem 4
def calcular_prox_ponto(x, k, h):
    return (h * k) + x

#Calcula uma lista de valores para x usando o Método de Runge-Kutta de Ordem 4
def RK4(x0, h, lista_t, A):
    #Inicia a lista de valores de vetores x
    lista_x = [x0]

    xi = x0
    for t in lista_t[1:]:
        #k1 (taxa de crescimento deste ponto)
        k1 = calcular_k1_1(xi, A)
        x_prov = calcular_prox_ponto(xi, k1, h / 2) #x provisório
        #k2
        k2 = calcular_k1_1(x_prov, A)
        x_prov = calcular_prox_ponto(xi, k2, h / 2)
        #k3
        k3 = calcular_k1_1(x_prov, A)
        x_prov = calcular_prox_ponto(xi, k3, h)
        #k4
        k4 = calcular_k1_1(x_prov, A)

        #média ponderada dos k
        kp = (k1 + 2*k2 + 2*k3 + k4) / 6

        #valor de x estimado para este t:
        x = calcular_prox_ponto(xi, kp, h)
        lista_x.append(x)

        #Atualiza o valor de xi para x
        xi = x

    return lista_x

#Calcula os erros de cada aproximação dos vetores x em cada t, tomando como erro
#a maior diferença entre cada um dos 4 valores de x explícito e aproximado
#(específico para o exercício 1.1)
def calcular_erros1_1(lista_x_exp, lista_x, lista_t):
    lista_erro = []

    for i in range(len(lista_t)):
        maior_erro = 0
        for j in range(0, 4):
            erro = abs(lista_x_exp[j](lista_t[i]) - lista_x[i][j])

            if erro >= maior_erro:
                maior_erro = erro

        lista_erro.append(maior_erro)

    return lista_erro, max(lista_erro)

#Resolve o exercício 1, iterando o cálculo do Método de Runge-Kutta e dos erros
#obtidos para cada n e plotando um gráfico
def exercicio1_1(T0, Tf, lista_n, x0, A, lista_x_exp):
    lista_erro_max = []
    for n in lista_n:
        #Parâmetros do intervalo
        h = calcular_passo(T0, Tf, n)
        lista_t = calcular_lista_t(T0, Tf, h)

        #Runge-Kutta de Ordem 4 para calcular os valores dos vetores x
        lista_x = RK4(x0, h, lista_t, A)

        #Cálculo dos erros, usando os valores explícitos de x, incluindo lista de maiores erros
        lista_erro, erro_max = calcular_erros1_1(lista_x_exp, lista_x, lista_t)
        lista_erro_max.append(erro_max)

        #Plotagem do Gráfico
        plt.plot(lista_t, lista_erro)
        plt.xlabel('t')
        plt.ylabel('Erro(t,n=' + str(n) + ')')
        plt.title('Erro dos valores aproximados de x para cada t num passo ' + str(h))
        plt.show()

    #Cálculo da Razão dos Erros
    lista_R = []
    for i in range (0,5):
        #print(lista_erro_max[i])
        lista_R.append(lista_erro_max[i]/lista_erro_max[i+1])
        print("O valor de R"+str(i+1)+" eh igual a " + str(lista_R[i]))
    #plt.show() #Essa parte é necessária?
#---------------------------------------------------------------------------------------------------
######################## EXERCICO 1 PARTE 2 ##########################################################

#Valor de um polinômio num dado ponto
def valor_polinomio(polinomio, x):
    #Polinomio de Grau 0
    if len(polinomio) == 0:
        return 0

    valor = polinomio[0]
    for i in range(1, len(polinomio)):
        valor += pow(x, i) * polinomio[i]

    return valor

def derivada_polinomio(polinomio):
    #Polinomio de Grau 0
    if len(polinomio) == 0:
        raise Exception("Nao esqueca de definir o polinomio!\n")
        return None

    d_polinomio = []
    for i in range(1, len(polinomio)):
        coef = i * polinomio[i]
        d_polinomio.append(coef)

    return d_polinomio

#Método de Newton
def metodo_newton(Xi, polinomio, max_iteracoes, n_iteracoes=0, d_polinomio = None):

    #Estorou o limite de iteracoes
    if n_iteracoes > max_iteracoes:
        #print("Numero maximo de iteracoes atingido")
        return Xi

    n_iteracoes += 1

    numerador = 0
    denominador = 0

    #Polinômio
    numerador += valor_polinomio(polinomio, Xi)

    if d_polinomio == None:
        d_polinomio = derivada_polinomio(polinomio)

    denominador += valor_polinomio(d_polinomio, Xi)

    #Xi já é uma raiz
    if numerador == 0:
        return Xi

    #Ponto estacionário (Exception Handling)
    if denominador == 0:
        raise ZeroDivisionError("A derivada da funcao vale 0 no ponto " + str(Xi) +'...\n')
        return None

    #Valor desta Iteração
    X = Xi - (numerador / denominador)

    return metodo_newton(Xi=X, polinomio=polinomio, max_iteracoes=max_iteracoes, n_iteracoes=n_iteracoes, d_polinomio=d_polinomio)


#Método de Euler Implícito
def euler_imp(lista_t, x0, montar_G, h):
    lista_x = [x0]
    xi = x0
    for t in lista_t[1:]:
        polinomio_G = montar_G(h, t, xi)
        x = metodo_newton(Xi=xi, polinomio=polinomio_G, max_iteracoes=7)
        lista_x.append(x)
        xi = x

    return lista_x

#Calcula os Erros de cada aproximação em relação à fórmula explícita
def calcular_erros1_2(lista_x, lista_x_exp):
    lista_erros = []
    for i in range(len(lista_x)):
        erro = abs(lista_x[i] - lista_x_exp[i])
        lista_erros.append(erro)

    return lista_erros

#Resolve o Exerício 1.2
def exercicio1_2(T0, Tf, n, x0, x_explicita, montar_G):
    #Parâmetros do Intervalo
    h = calcular_passo(T0, Tf, n)
    lista_t = calcular_lista_t(T0, Tf, h)

    #Método de Euler Implícito
    lista_x = euler_imp(lista_t, x0, montar_G, h)

    #Lista dos valores corretos de X
    lista_x_exp = []
    for t in lista_t:
        x_exp = x_explicita(t)
        lista_x_exp.append(x_exp)

    lista_erros = calcular_erros1_2(lista_x, lista_x_exp)

    plt.plot(lista_t, lista_x)
    plt.xlabel('t')
    plt.ylabel('x(t)')
    plt.title("Solucao aproximada")
    plt.show()

    plt.plot(lista_t, lista_x_exp)
    plt.xlabel('t')
    plt.ylabel('x*(t)')
    plt.title("Solucao explicita")
    plt.show()

    plt.plot(lista_t, lista_erros)
    plt.xlabel('t')
    plt.ylabel('Erro(t)')
    plt.title("Erro")
    plt.show()
#---------------------------------------------------------------------------------------------------
######################## EXERCICO 2 PARTE 1 ##########################################################

#Euler Explícito para duas funções
def euler_expl_duplo(h, lambida, alpha, beta, gama, x0, y0, derivada_x, derivada_y, lista_t):
    lista_x = [x0]
    lista_y = [y0]

    xi = x0
    yi = y0
    for i in range(len(lista_t) -1):
        x = xi + h * derivada_x(xi, yi, alpha, lambida)
        lista_x.append(x)
        y = yi + h * derivada_y(xi, yi, beta, gama)
        lista_y.append(y)

        xi = x
        yi = y

    return lista_x, lista_y


#Resolve o Exercício 2.1
def exercicio2_1(lambida, alpha, beta, gama, x0, y0, T0, Tf, n, derivada_x, derivada_y, plotar):
    #Parâmetros do intervalo
    h = calcular_passo(T0, Tf, n)
    lista_t = calcular_lista_t(T0, Tf, h)

    #Resolução Euler Explícito de duas Funções
    lista_x, lista_y = euler_expl_duplo(h, lambida, alpha, beta, gama, x0, y0, derivada_x, derivada_y, lista_t)

    if plotar:
        #Plotagem dos Gráficos
        plt.plot(lista_t, lista_x, label="Coelhos", color = "red")
        plt.plot(lista_t, lista_y, label="Raposas", color = "blue")
        plt.title("Populações de Raposas e Coelhos em função do tempo - Euler Explícito")
        plt.show()

        plt.plot(lista_x, lista_y)
        plt.xlabel("Coelhos")
        plt.ylabel("Raposas")
        plt.title("População de Raposas em função da População de Coelhos - Euler Explícito")
        plt.show()

    return lista_x, lista_y
#---------------------------------------------------------------------------------------------------
######################## EXERCICO 2 PARTE 2 ##########################################################

def metodo_newton_duplo(Ui, montar_G_duplo, montar_J, h, max_iteracoes=7, n_iteracoes=0, Ui0=None, Ui0_iniciado=False):
    #Estorou o limite de iteracoes
    if n_iteracoes > max_iteracoes:
        #print("Numero maximo de iteracoes atingido")
        return Ui

    n_iteracoes += 1

    #Primeira Iteração
    #(Precisamos do booleano Ui0_iniciado pois não é possível fazer Ui0 == None
    #diretamente depois que Ui0 já está iniciado)
    if not Ui0_iniciado:
        Ui0 = Ui
        Ui0_iniciado = True

    G = montar_G_duplo(Ui[0][0], Ui[1][0], Ui0[0][0], Ui0[1][0], h)

    J = montar_J(Ui[0][0], Ui[1][0], h)
    J_inv = np.linalg.inv(J)

    U = Ui - np.matmul(J_inv, G)

    return metodo_newton_duplo(Ui=U, montar_G_duplo=montar_G_duplo, montar_J=montar_J, h=h,
    max_iteracoes=max_iteracoes, n_iteracoes=n_iteracoes, Ui0=Ui0, Ui0_iniciado=Ui0_iniciado)

#Resolve o sistema de duas EDOs com o método de Euler Implícito
def euler_impl_duplo(x0, y0, montar_G_duplo, h, montar_J, lista_t):
    lista_x = [x0]
    lista_y = [y0]

    xi = x0
    yi = y0

    for i in range(len(lista_t) -1):
        Ui = np.array([[xi], [yi]])

        U = metodo_newton_duplo(Ui, montar_G_duplo, montar_J, h, max_iteracoes=7)
        lista_x.append(U[0][0])
        lista_y.append(U[1][0])

        xi = U[0][0]
        yi = U[1][0]

    return lista_x, lista_y

#Resolve o Exercício 2.2
def exercicio2_2(x0, y0, T0, Tf, n, montar_G_duplo, montar_J, plotar):
    #Parâmetros do intervalo
    h = calcular_passo(T0, Tf, n)
    lista_t = calcular_lista_t(T0, Tf, h)

    lista_x, lista_y = euler_impl_duplo(x0, y0, montar_G_duplo, h, montar_J, lista_t)

    if plotar:
        #Plotagem dos Gráficos
        plt.plot(lista_t, lista_x, label="Coelhos", color = "red")
        plt.plot(lista_t, lista_y, label="Raposas", color = "blue")
        plt.title("Populações de Raposas e Coelhos em função do tempo - Euler Implícito")
        plt.show()

        plt.plot(lista_x, lista_y)
        plt.xlabel("Coelhos")
        plt.ylabel("Raposas")
        plt.title("População de Raposas em função da População de Coelhos - Euler Implícito")
        plt.show()

    return lista_x, lista_y
#---------------------------------------------------------------------------------------------------
######################## EXERCICO 2 PARTE 3 ##########################################################

#Resolve o Exerício 2.3
def exercicio2_3(lambida, alpha, beta, gama, x0, y0, T0, Tf, lista_n, derivada_x, derivada_y, montar_G_duplo, montar_J):
    for n in lista_n:
        #Parâmetros do intervalo
        h = calcular_passo(T0, Tf, n)
        lista_t = calcular_lista_t(T0, Tf, h)

        lista_x_2_1, lista_y_2_1 = exercicio2_1(lambida, alpha, beta, gama, x0, y0, T0, Tf, n, derivada_x, derivada_y, plotar=False)

        lista_x_2_2, lista_y_2_2 = exercicio2_2(x0, y0, T0, Tf, n, montar_G_duplo, montar_J, plotar=False)

        lista_erro_x = []
        lista_erro_y = []
        for i in range((len(lista_x_2_1))):
            #Ximp - Xexp
            erro_x = lista_x_2_2[i] - lista_x_2_1[i]
            lista_erro_x.append(erro_x)

            #Yimp - Yexp
            erro_y = lista_y_2_2[i] - lista_y_2_1[i]
            lista_erro_y.append(erro_y)

        #Plotaagem do Gráfico
        plt.plot(lista_t, lista_erro_x, label="Diferença em X", color = "blue")
        plt.plot(lista_t, lista_erro_x, label="Diferença em Y", color = "red")
        plt.title("Diferença entre Euler Implícito e Explícito")
        plt.show()
#---------------------------------------------------------------------------------------------------
######################## EXERCICO 2 PARTE 4 ##########################################################

#Aplicação do Método de Runge Kutta para um sistema de duas equações
def RK4_duplo(lista_t, h, x0, y0, derivada_x, derivada_y, alpha, lambida, beta, gama):
    lista_x = [x0]
    lista_y = [y0]

    xi = x0
    yi = y0
    for i in range(len(lista_t) - 1):

        #k1
        k1_x = derivada_x(xi, yi, alpha, lambida)
        k1_y = derivada_y(xi, yi, beta, gama)
        x_prov = calcular_prox_ponto(xi, k1_x, h / 2) #x provisório
        y_prov = calcular_prox_ponto(yi, k1_y, h / 2) #y provisório

        #k2
        k2_x = derivada_x(x_prov, y_prov, alpha, lambida)
        k2_y = derivada_y(x_prov, y_prov, beta, gama)
        x_prov = calcular_prox_ponto(xi, k2_x, h / 2)
        y_prov = calcular_prox_ponto(yi, k2_y, h / 2)

        #k3
        k3_x = derivada_x(x_prov, y_prov, alpha, lambida)
        k3_y = derivada_y(x_prov, y_prov, beta, gama)
        x_prov = calcular_prox_ponto(xi, k3_x, h)
        y_prov = calcular_prox_ponto(yi, k3_y, h)

        #k4
        k4_x = derivada_x(x_prov, y_prov, alpha, lambida)
        k4_y = derivada_y(x_prov, y_prov, beta, gama)

        #média ponderada dos k
        kp_x = (k1_x + 2*k2_x + 2*k3_x + k4_x) / 6
        kp_y = (k1_y + 2*k2_y + 2*k3_y + k4_y) / 6

        #valor de x estimado para este t:
        x = calcular_prox_ponto(xi, kp_x, h)
        y = calcular_prox_ponto(yi, kp_y, h)

        lista_x.append(x)
        lista_y.append(y)

        #Atualiza o valor de xi e yi
        xi = x
        yi = y

    return lista_x, lista_y

#Resolve o Exercício 2.4
def exercicio2_4(T0, Tf, n, derivada_x, derivada_y, x0, y0, alpha, lambida, beta, gama):
    #Parâmetros do Intervalo
    h = calcular_passo(T0, Tf, n)
    lista_t = calcular_lista_t(T0, Tf, h)

    lista_x, lista_y = RK4_duplo(lista_t, h, x0, y0, derivada_x, derivada_y, alpha, lambida, beta, gama)

    #Plotagem dos Gráficos
    plt.plot(lista_t, lista_x, label="Coelhos", color = "red")
    plt.plot(lista_t, lista_y, label="Raposas", color = "blue")
    plt.title("Populações de Raposas e Coelhos em função do tempo - Runge Kutta 4")
    plt.show()

    plt.plot(lista_x, lista_y)
    plt.xlabel("Coelhos")
    plt.ylabel("Raposas")
    plt.title("População de Raposas em função da População de Coelhos - Runge Kutta 4")
    plt.show()
#---------------------------------------------------------------------------------------------------
######################## EXERCICO 3 PARTE 1 ##########################################################

def euler_expl_triplo(h, lista_t, x0, derivada_x, y0, derivada_y, z0, derivada_z, alpha):
    lista_x = [x0]
    lista_y = [y0]
    lista_z = [z0]

    xi = x0
    yi = y0
    zi = z0
    for i in range(len(lista_t) -1):
        x = xi + h * derivada_x(xi, yi, zi)
        lista_x.append(x)
        y = yi + h * derivada_y(xi, yi, zi)
        lista_y.append(y)
        z = zi + h * derivada_z(xi, yi, zi, alpha)
        lista_z.append(z)

        xi = x
        yi = y
        zi = z

    return lista_x, lista_y, lista_z


def RK4_triplo(h, lista_t, x0, derivada_x, y0, derivada_y, z0, derivada_z, alpha):
        lista_x = [x0]
        lista_y = [y0]
        lista_z = [z0]

        xi = x0
        yi = y0
        zi = z0
        for i in range(len(lista_t) - 1):

            #k1
            k1_x = derivada_x(xi, yi, zi)
            k1_y = derivada_y(xi, yi, zi)
            k1_z = derivada_z(xi, yi, zi, alpha)
            x_prov = calcular_prox_ponto(xi, k1_x, h / 2) #x provisório
            y_prov = calcular_prox_ponto(yi, k1_y, h / 2) #y provisório
            z_prov = calcular_prox_ponto(zi, k1_z, h / 2) #z provisório

            #k2
            k2_x = derivada_x(x_prov, y_prov, z_prov)
            k2_y = derivada_y(x_prov, y_prov, z_prov)
            k2_z = derivada_z(x_prov, y_prov, z_prov, alpha)
            x_prov = calcular_prox_ponto(xi, k2_x, h / 2)
            y_prov = calcular_prox_ponto(yi, k2_y, h / 2)
            z_prov = calcular_prox_ponto(zi, k2_z, h / 2)

            #k3
            k3_x = derivada_x(x_prov, y_prov, z_prov)
            k3_y = derivada_y(x_prov, y_prov, z_prov)
            k3_z = derivada_z(x_prov, y_prov, z_prov, alpha)
            x_prov = calcular_prox_ponto(xi, k3_x, h)
            y_prov = calcular_prox_ponto(yi, k3_y, h)
            z_prov = calcular_prox_ponto(zi, k3_z, h)

            #k4
            k4_x = derivada_x(x_prov, y_prov, z_prov)
            k4_y = derivada_y(x_prov, y_prov, z_prov)
            k4_z = derivada_z(x_prov, y_prov, z_prov, alpha)

            #média ponderada dos k
            kp_x = (k1_x + 2*k2_x + 2*k3_x + k4_x) / 6
            kp_y = (k1_y + 2*k2_y + 2*k3_y + k4_y) / 6
            kp_z = (k1_z + 2*k2_z + 2*k3_z + k4_z) / 6

            #valor de x estimado para este t:
            x = calcular_prox_ponto(xi, kp_x, h)
            y = calcular_prox_ponto(yi, kp_y, h)
            z = calcular_prox_ponto(zi, kp_z, h)

            lista_x.append(x)
            lista_y.append(y)
            lista_z.append(z)

            #Atualiza o valor de xi, yi e zi
            xi = x
            yi = y
            zi = z

        return lista_x, lista_y, lista_z

#Resolve o Exercício 3.1
def exercicio3_1(lista_alpha, x0, y0, z0, T0, lista_Tf, derivada_x, derivada_y, derivada_z, n):
    for i in range(len(lista_alpha)):
        alpha = lista_alpha[i]
        Tf = lista_Tf[i]

        #Parâmetros do Intervalo
        h = calcular_passo(T0, Tf, n)
        lista_t =  calcular_lista_t(T0, Tf, h)

        #Euler Explícito
        lista_x, lista_y, lista_z = euler_expl_triplo(h, lista_t, x0, derivada_x, y0, derivada_y, z0, derivada_z, alpha)

        #Plotagem dos Gráficos
        axes = plt.axes(projection = '3d')
        axes.plot3D(lista_x, lista_y, lista_z, 'purple')
        plt.title("Modelo Presa-Predador - Euler Explícito - com alpha = " + str(alpha) + " e intervalo [" + str(T0) + ',' + str(Tf) + ']')
        axes.set_xlabel("Coelhos")
        axes.set_ylabel("Lebres")
        axes.set_zlabel("Raposas")
        plt.show()

        #Runge Kutta 4
        lista_x, lista_y, lista_z = RK4_triplo(h, lista_t, x0, derivada_x, y0, derivada_y, z0, derivada_z, alpha)

        #Plotagem do Gráfico
        axes = plt.axes(projection = '3d')
        axes.plot3D(lista_x, lista_y, lista_z, 'purple')
        plt.title("Modelo Presa-Predador - Runge Kutta 4 - com alpha = " + str(alpha) + " e intervalo [" + str(T0) + ',' + str(Tf) + ']')
        axes.set_xlabel("Coelhos")
        axes.set_ylabel("Lebres")
        axes.set_zlabel("Raposas")
        plt.show()
