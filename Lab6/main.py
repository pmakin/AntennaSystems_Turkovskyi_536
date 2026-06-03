import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv

lam = 0.034
D = 1.4
f = 0.4

k = 2 * np.pi / lam
R0 = D / 2
p = 2 * f
v = 3.5 * (R0 / p)

theta_deg = np.linspace(-20, 20, 8000)
theta = np.radians(theta_deg)

u = (k * R0) * np.sin(np.abs(theta))

def safe_div(a, b):
    eps = 1e-12
    denom = np.where(np.abs(b) < eps, eps, b)
    return a / denom

J0_u = jv(0, u)
J1_u = jv(1, u)
J2_u = jv(2, u)

J0_v = jv(0, v)
J1_v = jv(1, v)
J1_15v = jv(1, v * 1.5)
J2_v = jv(2, v)
J2_15v = jv(2, v * 1.5)

term1 = 0.74 * safe_div(v * J1_v * J0_u - u * J1_u * J0_v, (v**2 - u**2))
term2 = 0.26 * safe_div(J1_u, u)
term3 = 0.25 * safe_div(u * J1_u * J2_15v - (v * 1.5) * J1_15v * J2_u, ((v * 1.5)**2 - u**2))

norm_denom = 0.74 * (J1_v / v) + 0.13
norm = 1.0 / norm_denom

FE = (np.cos(theta / 2) ** 2) * (term1 + term2 + term3) * norm
FH = (np.cos(theta / 2) ** 2) * (term1 + term2 - term3) * norm

FE = np.abs(FE)
FH = np.abs(FH)
FE /= np.max(FE)
FH /= np.max(FH)

level = 0.707
def get_beamwidth_24(theta_deg, F):
    pos_idx = np.where((theta_deg >= 0) & (F >= level))[0]
    right_edge = theta_deg[pos_idx[-1]]
    return 2 * right_edge

bw_24_E = get_beamwidth_24(theta_deg, FE)
bw_24_H = get_beamwidth_24(theta_deg, FH)

print(f"Двобічна ширина пелюстки 2θ_0.5 (E): {bw_24_E:.2f}°")
print(f"Двобічна ширина пелюстки 2θ_0.5 (H): {bw_24_H:.2f}°")

plt.figure(figsize=(10, 5))
plt.plot(theta_deg, FH, label='H-площина')
plt.plot(theta_deg, FE, label='E-площина', linestyle='--')
plt.axhline(level, color='r', linestyle='--', label='Рівень 0.707')
plt.xlim(-5, 5)
plt.xlabel('θ, градуси')
plt.ylabel('Нормована ДС')
plt.title('іаграма спрямованості ДзА')
plt.grid(True)
plt.legend()
plt.savefig("ДС_ДзА.jpg", dpi=600, bbox_inches='tight')
plt.show()