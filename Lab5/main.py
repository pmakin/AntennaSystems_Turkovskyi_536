import numpy as np
import matplotlib.pyplot as plt

lam = 3.4
N = 6
d = 3.2
slot_length = 1.6
v = 1
lambda_v = 4.0

k = 2 * np.pi / lam

theta_deg = np.linspace(-90, 90, 5000)
theta = np.radians(theta_deg)

x1 = (np.pi * slot_length / lam) * np.sin(theta)
F1H = np.ones_like(theta)
mask1 = np.abs(x1) > 1e-12
F1H[mask1] = np.sin(x1[mask1]) / x1[mask1]
F1H = np.abs(F1H)

psi = k * d * np.sin(theta) - (2 * np.pi / lambda_v) * d + v * np.pi
num = np.sin(N * psi / 2)
den = N * np.sin(psi / 2)

FHC = np.ones_like(theta)
mask2 = np.abs(den) > 1e-12
FHC[mask2] = num[mask2] / den[mask2]
FHC = np.abs(FHC)

FH = F1H * FHC
FH = FH / np.max(FH)

peaks = []
for i in range(1, len(FH) - 1):
    if FH[i] > FH[i - 1] and FH[i] > FH[i + 1]:
        peaks.append(i)

main_maxima_indices = []
side_lobes_indices = []

for idx in peaks:
    if FHC[idx] > 0.95:
        main_maxima_indices.append(idx)
    else:
        side_lobes_indices.append(idx)

for idx in main_maxima_indices:
    print(f"Кут θ = {theta_deg[idx]:.2f}°, Результуюча амплітуда F_H = {FH[idx]:.3f}")
print("==================================================")

abs_max_idx = main_maxima_indices[1]  # це наш глобальний максимум
half_power = 0.707

left_idx = abs_max_idx
while left_idx > 0 and FH[left_idx] >= half_power:
    left_idx -= 1

right_idx = abs_max_idx
while right_idx < len(FH) - 1 and FH[right_idx] >= half_power:
    right_idx += 1

beamwidth = theta_deg[right_idx] - theta_deg[left_idx]
print(f"Ширина головної пелюстки (0.707): {beamwidth:.2f}°")

highest_side_lobe_idx = side_lobes_indices[np.argmax(FH[side_lobes_indices])]
side_lobe_val = FH[highest_side_lobe_idx]

print(f"Максимальний рівень бокових пелюсток:")
print(f"Амплітуда: {side_lobe_val:.3f} (на куті {theta_deg[highest_side_lobe_idx]:.2f}°)")
print(f"Дб: {20 * np.log10(side_lobe_val):.2f} дБ")
print("==================================================")

print("\nТаблиця 1: Нулі ДС:")
for m in range(-15, 16):
    if m % N == 0:
        continue
    arg = (lam * m) / (d * N) + lam / lambda_v - (v * lam) / (2 * d)
    if -1 <= arg <= 1:
        theta0 = np.degrees(np.arcsin(arg))
        print(f"m = {m:3d} -> θ0 = {theta0:.2f}°")

print("\nТаблиця 2: Максимуми бокових пелюсток:")
for n in range(-15, 16):
    arg = ((2 * n + 1) * lam) / (2 * d * N) + lam / lambda_v - (v * lam) / (2 * d)
    if -1 <= arg <= 1:
        theta_max = np.degrees(np.arcsin(arg))
        theta_max_rad = np.arcsin(arg)

        x1_max = (np.pi * slot_length / lam) * np.sin(theta_max_rad)
        f1_val = 1.0 if abs(x1_max) < 1e-12 else abs(np.sin(x1_max) / x1_max)
        fhc_val = abs(1 / (N * np.sin(((2 * n + 1) * np.pi) / (2 * N))))

        total_amplitude = f1_val * fhc_val
        print(f"n = {n:3d} -> θmax = {theta_max:.2f}°, Теоретична амплітуда F = {total_amplitude:.3f}")

plt.figure(figsize=(12, 7))
plt.plot(theta_deg, F1H, '--', label='F1H(θ) (ДС однієї щілини)', alpha=0.7)
plt.plot(theta_deg, FHC, ':', label='FHC(θ) (Множник решітки)', alpha=0.7)
plt.plot(theta_deg, FH, linewidth=2.5, color='green', label='F_H(θ) = F1H * FHC (Результуюча)')

plt.axhline(y=half_power, color='red', linestyle='--', label='Рівень 0.707')
plt.axvline(theta_deg[left_idx], color='red', linestyle=':')
plt.axvline(theta_deg[right_idx], color='red', linestyle=':')

for idx in main_maxima_indices:
    plt.plot(theta_deg[idx], FH[idx], 'ro')
    plt.text(theta_deg[idx], FH[idx] + 0.03, f"{theta_deg[idx]:.1f}°\n(F={FH[idx]:.2f})", ha='center', fontsize=9)

plt.title('Нормована діаграма спрямованості ХЩА', fontsize=12)
plt.xlabel('Кут θ (градуси)')
plt.ylabel('Нормована амплітуда')
plt.xlim(-90, 90)
plt.ylim(0, 1.1)
plt.xticks(np.arange(-90, 91, 15))
plt.grid(True)
plt.legend(loc='upper right')
plt.savefig("ДС_ХЩА.jpg", dpi=600, bbox_inches='tight')
plt.show()