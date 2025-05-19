from sympy import symbols, Eq, solve, diff
import numpy as np

# Definir las coordenadas naturales chi, eta
chi, eta = symbols('chi eta')
X = 0.8
Y = 0.3

# Funciones de interpolación (ya definidas)
N1 = (1 - chi) * (1 - eta) / 4
N2 = (1 + chi) * (1 - eta) / 4
N3 = (1 + chi) * (1 + eta) / 4
N4 = (1 - chi) * (1 + eta) / 4

# Asumimos desplazamientos para los nodos (ejemplo simple)
# Vector de desplazamientos (en este caso, 8 valores)

coords = np.array([
    [0.0, 0.0],   # Nodo 1
    [1.0, 1.0],   # Nodo 2
    [1.0, 0.0],   # Nodo 3
    [0.0, 1.0]    # Nodo 4
])

# Formulamos las ecuaciones basadas en las funciones de interpolación
# Suponiendo que el vector de desplazamientos se ha formulado correctamente, lo que queremos es un sistema de ecuaciones.
eq1 = Eq(N1*coords[0][0] + N2*coords[1][0]+ N3*coords[2][0]+N4*coords[3][0],X)  # Primera ecuación basada en el primer nodo
eq2 = Eq(N1*coords[0][1] + N2*coords[1][1]+ N3*coords[2][1]+N4*coords[3][1],Y)  # Segunda ecuación basada en el segundo nodo

# Resolver el sistema de ecuaciones
soluciones = solve((eq1, eq2), (chi, eta))

# Mostrar las soluciones
if soluciones:
    # Si hay soluciones, se muestran
    for sol in soluciones:
        chi_valor = sol[0]  # Accede a chi (primer valor en la tupla)
        eta_valor = sol[1]  # Accede a eta (segundo valor en la tupla)
    print(f"Solución para chi: {chi_valor}")
    print(f"Solución para eta: {eta_valor}")
else:
    print("No se encontraron soluciones.")


chi = chi_valor
eta = eta_valor

print(f"X: {X}")
print(f"Y: {Y}")

#Pregunta D   

# Definir las coordenadas naturales chi y eta
dN_dchi = 0.25 * np.array([
    -(1 - eta), (1 - eta), (1 + eta), -(1 + eta)
])

dN_deta = 0.25 * np.array([
    -(1 - chi), -(1 + chi), (1 + chi), (1 - chi)
])

# Jacobiano
J = np.zeros((2, 2))
for i in range(4):
    J[0, 0] += dN_dchi[i] * coords[i, 0]
    J[0, 1] += dN_deta[i] * coords[i, 0]
    J[1, 0] += dN_dchi[i] * coords[i, 1]
    J[1, 1] += dN_deta[i] * coords[i, 1]

# Inversa del Jacobiano
J_inv = np.linalg.inv(J)
print("Jacobiano:")
print(J)
print('Determinante del Jacobiano:')
print(np.linalg.det(J))
