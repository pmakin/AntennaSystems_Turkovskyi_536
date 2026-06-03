import numpy
import matplotlib.pyplot as plt
import math

variant = 21
lambd = 0.034
a = 0.14
b = 0.16

F1e = [1]
FCe = [1]
FE = [1]
steps_E = [0]
SGP_e = 0
hp_angle_e = 0

max_x_FE = []
max_y_FE = []
min_x_FE = []
min_y_FE = []

for teta in numpy.arange(0.00001, numpy.pi / 2, 0.00001):
    mn1 = abs((1 + numpy.cos(teta)) / 2)
    u_e = (numpy.pi * b * math.sin(teta)) / lambd
    mn2 = abs(math.sin(u_e) / u_e)
    mn3 = mn1 * mn2

    F1e.append(mn1)
    FCe.append(mn2)
    FE.append(mn3)

    if 0.707 < mn3 < 0.708 and SGP_e == 0:
        hp_angle_e = math.degrees(teta)
        SGP_e = 2 * hp_angle_e

    steps_E.append(math.degrees(teta))

for i in range(1, len(FE) - 1):
    if FE[i] > FE[i - 1] and FE[i] > FE[i + 1]:
        max_x_FE.append(steps_E[i])
        max_y_FE.append(FE[i])

for i in range(1, len(FE) - 1):
    if FE[i] < FE[i - 1] and FE[i] < FE[i + 1]:
        min_x_FE.append(steps_E[i])
        min_y_FE.append(0)

F1h = [1]
FCh = [1]
FH = [1]
steps_H = [0]
SGP_h = 0
hp_angle_h = 0

max_x_FH = []
max_y_FH = []
min_x_FH = []
min_y_FH = []

for teta in numpy.arange(0.00001, numpy.pi / 2, 0.00001):
    mn1 = abs((1 + numpy.cos(teta)) / 2)
    u_h = (numpy.pi * a * math.sin(teta)) / lambd
    v_h = (2 * a * math.sin(teta)) / lambd
    denominator = 1 - v_h ** 2

    if abs(denominator) < 1e-6:
        mn2 = numpy.pi / 4
    else:
        mn2 = abs(math.cos(u_h) / denominator)

    mn3 = mn1 * mn2

    F1h.append(mn1)
    FCh.append(mn2)
    FH.append(mn3)

    if 0.707 < mn3 < 0.708 and SGP_h == 0:
        hp_angle_h = math.degrees(teta)
        SGP_h = 2 * hp_angle_h

    steps_H.append(math.degrees(teta))

for i in range(1, len(FH) - 1):
    if FH[i] > FH[i - 1] and FH[i] > FH[i + 1]:
        max_x_FH.append(steps_H[i])
        max_y_FH.append(FH[i])

for i in range(1, len(FH) - 1):
    if FH[i] < FH[i - 1] and FH[i] < FH[i + 1]:
        min_x_FH.append(steps_H[i])
        min_y_FH.append(0)

print(f"Варіант = {variant}")
print(f"Довжина хвилі = {lambd} (м)")
print(f"Розмір розкриву рупора a*b = {a}*{b} (м)\n")

print("--- Площина E ---")
print("Табл. 1 - Значення нульових кутів")
print("| № |   θ   | FE(θ) |")
for i, (x, y) in enumerate(zip(min_x_FE, min_y_FE), 1):
    print(f"| {i} | {x:5.2f} |   {y}   |")

print("\nТабл. 2 - Значення максимальних кутів")
print("| № |   θ   | FE(θ) |")
for i, (x, y) in enumerate(zip(max_x_FE, max_y_FE), 1):
    print(f"| {i} | {x:5.2f} | {y:5.2f} |")

print(f"\nШирина головної пелюстки в площині E = {SGP_e:.2f}°\n")

print("--- Площина H ---")
print("Табл. 1 - Значення нульових кутів")
print("| № |   θ   | FH(θ) |")
for i, (x, y) in enumerate(zip(min_x_FH, min_y_FH), 1):
    print(f"| {i} | {x:5.2f} |   {y}   |")

print("\nТабл. 2 - Значення максимальних кутів")
print("| № |   θ   | FH(θ) |")
for i, (x, y) in enumerate(zip(max_x_FH, max_y_FH), 1):
    print(f"| {i} | {x:5.2f} | {y:5.2f} |")

print(f"\nШирина головної пелюстки в площині H = {SGP_h:.2f}°")

# E
plt.figure('Площина E', figsize=(10, 6))
plt.plot(steps_E, F1e, label=r'$F_{1e}(\theta)$', linewidth=1)
plt.plot(steps_E, FCe, label=r'$F_C(\theta)$', linewidth=1)
plt.plot(steps_E, FE, label=r'$F_E(\theta)$', linewidth=1)

plt.plot([0, hp_angle_e], [0.707, 0.707], color='red', linestyle='--', linewidth=0.8, alpha=0.7)
plt.plot([hp_angle_e, hp_angle_e], [0, 0.707], color='red', linestyle='--', linewidth=0.8, alpha=0.7)

plt.scatter(hp_angle_e, 0.707, color='red', zorder=5, label='Рівень половинної потужності в площині E')
plt.annotate(f'(0.707, {hp_angle_e:.2f}°)', xy=(hp_angle_e, 0.707), xytext=(hp_angle_e + 4, 0.74),
             arrowprops=dict(arrowstyle='->', lw=1.5))

plt.scatter(max_x_FE, max_y_FE, color='black', zorder=5, label=r'$\theta_{max}$ E')
plt.scatter(min_x_FE, min_y_FE, color='blue', zorder=5, label=r'$\theta_{min}$ E')

plt.title(fr'ДС пірамідального рупору в площині E з параметрами: $\lambda$ = {lambd}, Bp = {b}')
plt.xlabel(r'$\Theta^{\circ}$')
plt.ylabel(r'$|FE(\theta^{\circ})|, |F1e(\theta^{\circ})|, |Fc(\theta^{\circ})|$')
plt.xlim(0, 90)
plt.ylim(-0.01, 1.01)
plt.xticks(numpy.arange(0, 91, 2))
plt.grid(True, linestyle='-', alpha=0.6)
plt.legend(loc='upper right')
plt.savefig("Лабораторна2_E.jpg", dpi=600)

# H
plt.figure('Площина H', figsize=(10, 6))
plt.plot(steps_H, F1h, label=r'$F_{1h}(\theta)$', linewidth=1)
plt.plot(steps_H, FCh, label=r'$F_C(\theta)$', linewidth=1)
plt.plot(steps_H, FH, label=r'$F_H(\theta)$', linewidth=1)

plt.plot([0, hp_angle_h], [0.707, 0.707], color='red', linestyle='--', linewidth=0.8, alpha=0.7)
plt.plot([hp_angle_h, hp_angle_h], [0, 0.707], color='red', linestyle='--', linewidth=0.8, alpha=0.7)

plt.scatter(hp_angle_h, 0.707, color='red', zorder=5, label='Рівень половинної потужності в площині H')
plt.annotate(f'(0.707, {hp_angle_h:.2f}°)', xy=(hp_angle_h, 0.707), xytext=(hp_angle_h + 4, 0.74),
             arrowprops=dict(arrowstyle='->', lw=1.5))

plt.scatter(max_x_FH, max_y_FH, color='black', zorder=5, label=r'$\theta_{max}$ H')
plt.scatter(min_x_FH, min_y_FH, color='blue', zorder=5, label=r'$\theta_{min}$ H')

plt.title(fr'ДС пірамідального рупору в площині H з параметрами: $\lambda$ = {lambd}, Ap = {a}')
plt.xlabel(r'$\Theta^{\circ}$')
plt.ylabel(r'$|Fh(\theta^{\circ})|, |F1h(\theta^{\circ})|, |Fc(\theta^{\circ})|$')
plt.xlim(0, 90)
plt.ylim(-0.01, 1.01)
plt.xticks(numpy.arange(0, 91, 2))
plt.grid(True, linestyle='-', alpha=0.6)
plt.legend(loc='upper right')

plt.savefig("Лабораторна2_H.jpg", dpi=600)
plt.show()