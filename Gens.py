import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import math
import scipy.special

# **********************************
# 1 Generation of reads
# **********************************

# %% volumes
# T1:
V11 = 0.5  # random.uniform(0, 1)
V12 = 1 - V11

# T2:
V21 = V11 * 0.6  # пока насильно уменьшаем объем первой ткани
V22 = 1 - V21

# %% Number of reads
N = 100

# %% Probabilities
# Variants:
L = 20
Theta = np.zeros((2, L))
np.random.seed(0)
result = np.random.uniform(size=L, low=0, high=1)
array_len = []
for i in range(1, 101):
    array_len.append(2 / 99 * i / 100)
array_len.sort(reverse=True)
print('Отрезки вероятностей:', array_len)
m = 0
for k in result:
    i = 1
    g = array_len[0]
    while g < k:
        i += 1
        g += array_len[i - 1]
    if m < int(L / 2):
        Theta[0, m] = i / 100
    else:
        Theta[1, m] = i / 100
    m += 1
print("Theta", Theta)

# Theta = np.array([[0.4, 0, 0.6, 0.2, 0.3, 0, 0.5, 0.3, 0, 0.5, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.7, 0.0, 0.5],
#                   [0.0, 0.4, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.3, 0.0, 0.8, 0.1, 0.6, 0, 0.3, 0.1, 0.5, 0, 0.7, 0]])


# Tissues
p11 = V11 / (V11 + V12)  # 1 index - time;  2 index - tissue
p21 = V21 / (V21 + V22)
print('P1/2', p11, p21)
# %% Number of reads 1 and 2 types
Reads = np.zeros((2, L))  # [0] N11, [1] N21
for variants in range(L):
    Reads[0, variants] = sum(np.random.binomial(1, p11 * Theta[0, variants] + (1 - p11) * Theta[1, variants], N))
    Reads[1, variants] = sum(np.random.binomial(1, p21 * Theta[0, variants] + (1 - p21) * Theta[1, variants], N))
A = (Reads[0, :] + Reads[1, :]) * 0.5
# print("Result 1", p11 * Theta[0, :] + (1 - p11) * Theta[1, :])
# print("Result 2", p21 * Theta[0, :] + (1 - p21) * Theta[1, :])
print('Метод макс. правдоподобия PiK', A)

# **********************************
# 2 Volume ratios
# **********************************
TissueA = []
TissueB = []
target = 3  # 1 Tissue
Theta2 = np.zeros((2, L))
p1 = V11 / (V11 + V12)
p2 = V21 / (V21 + V22)
for variants in range(L):
    # Пусть рид из первой ткани
    PiK = p1 * Theta2[0, variants]
    D = (p2 * Reads[0, variants]) ** 2 + (p1 * N) ** 2 + (p1 * Reads[1, variants]) ** 2 + (p2 * N) ** 2 - 6 * p1 * p2 * \
        Reads[0, variants] * N + 2 * p1 * p2 * Reads[0, variants] * Reads[1, variants] + 2 * (p2 ** 2) * Reads[
            0, variants] * N + 2 * (p1 ** 2) * Reads[1, variants] * N + 2 * p1 * p2 * (N ** 2) - 6 * p1 * p2 * Reads[
            1, variants] * N
    ThetaK1 = (p2 * Reads[0, variants] + p1 * Reads[1, variants] + p1 * N + p2 * N - math.sqrt(D)) / (4 * p1 * p2 * N)

    alf1 = scipy.special.comb(N, Reads[0, variants]) * ((p1 * ThetaK1) ** Reads[0, variants]) * (
                (1 - p1 * ThetaK1) ** (N - Reads[0, variants])) * scipy.special.comb(N, Reads[1, variants]) * (
                      (p2 * ThetaK1) ** Reads[1, variants]) * (
                      (1 - p2 * ThetaK1) ** (N - Reads[1, variants]))


    # Пусть рид из второй ткани
    PiK = (1 - p1) * Theta2[1, variants]
    D1 = ((1 - p2) * Reads[0, variants]) ** 2 + ((1 - p1) * N) ** 2 + ((1 - p1) * Reads[1, variants]) ** 2 + (
            (1 - p2) * N) ** 2 - 6 * (1 - p1) * (1 - p2) * Reads[0, variants] * N + 2 * (1 - p1) * (1 - p2) * Reads[
             0, variants] * Reads[1, variants] + 2 * ((1 - p2) ** 2) * Reads[0, variants] * N + 2 * ((1 - p1) ** 2) * \
         Reads[1, variants] * N + 2 * (1 - p1) * (1 - p2) * (N ** 2) - 6 * (1 - p1) * (1 - p2) * Reads[1, variants] * N
    ThetaK2 = ((1 - p2) * Reads[0, variants] + (1 - p1) * Reads[1, variants] + (1 - p1) * N + (1 - p2) * N - math.sqrt(
        D1)) / (4 * (1 - p1) * (1 - p2) * N)

    alf2 = scipy.special.comb(N, Reads[0, variants]) * (((1 - p1) * ThetaK2) ** Reads[0, variants]) * ((
            1 - (1 - p1) * ThetaK2) ** (N - Reads[0, variants])) * \
          scipy.special.comb(N, Reads[1, variants]) * (((1 - p2) * ThetaK2) ** Reads[1, variants]) * ((
            1 - (1 - p2) * ThetaK2) ** (N - Reads[1, variants]))

    if alf1 > alf2:
        Theta2[0, variants] = ThetaK1
    else:
        Theta2[1, variants] = ThetaK2

    if (Reads[0, variants] - Reads[1, variants]) > 0:
        TissueA.append(variants)
    elif (Reads[0, variants] - Reads[1, variants]) < 0:
        TissueB.append(variants)
print("Theta2", Theta2)

# %%
# ************************************************
# 3 Visualisation: change in observed frequencies
# ************************************************

figure(figsize=(8, 6))
x = np.array(range(0, 2))
y = np.zeros((L, 2))

y[:, 0] = (Reads[0, :]).transpose()
y[:, 1] = (Reads[1, :]).transpose()

print("Наблюдаемые частоты в t1 и t2: \n", y)
plt.title("Observed frequencies plot")
plt.xlabel("X - time")
plt.ylabel("Y - frequencies")

for i, array in enumerate(y):
    if Theta2[0, i] > 0:
        plt.plot(x, array, color='#6667AB', marker="o", label=f"#{i}")
    else:
        plt.plot(x, array, color='red', marker="o", label=f"#{i}")


plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.show()


