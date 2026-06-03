import numpy
import matplotlib.pyplot as plt
import math

lambd = 0.034
l = 0.237
ksi = 1 + lambd / (2 * l)
h = 0.24

steps_s = [0]
Fb_s = [1]
FH_s = [1]
FE_s = [1]
steps_d = [0]
Fb_d = [1]
FC_d = [1]
FH_d = [1]
FE_d = [1]

for teta in numpy.arange(0.01, numpy.pi / 2, 0.0001):
    arg_den = numpy.pi * l / lambd * (ksi - 1)
    arg_num = numpy.pi * l / lambd * (ksi - numpy.cos(teta))

    if abs(ksi - numpy.cos(teta)) < 1e-10:
        mn_b = abs((ksi - 1) / numpy.sin(arg_den) * numpy.pi * l / lambd)
    else:
        mn_b = abs((ksi - 1) / numpy.sin(arg_den) * numpy.sin(arg_num) / (ksi - numpy.cos(teta)))

    mn_c = abs(numpy.cos(numpy.pi * h / lambd * numpy.sin(teta)))
    cos_t = abs(numpy.cos(teta))

    Fb_s.append(mn_b)
    FH_s.append(mn_b)
    FE_s.append(mn_b * cos_t)
    steps_s.append(math.degrees(teta))

    Fb_d.append(mn_b)
    FC_d.append(mn_c)
    FH_d.append(mn_b * mn_c)
    FE_d.append(mn_b * mn_c * cos_t)
    steps_d.append(math.degrees(teta))

max_FH_s = max(FH_s); FH_s = [v/max_FH_s for v in FH_s]
max_FE_s = max(FE_s); FE_s = [v/max_FE_s for v in FE_s]
max_Fb   = max(Fb_s);  Fb_s = [v/max_Fb  for v in Fb_s]; Fb_d = [v/max_Fb for v in Fb_d]
max_FH_d = max(FH_d); FH_d = [v/max_FH_d for v in FH_d]; FC_d = [v/max_FH_d for v in FC_d]
max_FE_d = max(FE_d); FE_d = [v/max_FE_d for v in FE_d]

def half_power(steps, F):
    idx = [i for i, v in enumerate(F) if v >= 0.707]
    if idx:
        return 2 * steps[idx[-1]], F[idx[-1]]
    return 0, 0

def find_extrema(steps, F):
    mx, my, mnx, mny = [], [], [], []
    for i in range(1, len(F) - 1):
        if F[i] > F[i - 1] and F[i] > F[i + 1]:
            mx.append(steps[i]); my.append(F[i])
        if F[i] < F[i - 1] and F[i] < F[i + 1]:
            mnx.append(steps[i]); mny.append(0)
    return mx, my, mnx, mny

SGP_H_s, fS_H_s = half_power(steps_s, FH_s)
SGP_E_s, fS_E_s = half_power(steps_s, FE_s)
SGP_H_d, fS_H_d = half_power(steps_d, FH_d)
SGP_E_d, fS_E_d = half_power(steps_d, FE_d)

mx_Hs, my_Hs, mnx_Hs, mny_Hs = find_extrema(steps_s, FH_s)
mx_Es, my_Es, mnx_Es, mny_Es = find_extrema(steps_s, FE_s)
mx_Hd, my_Hd, mnx_Hd, mny_Hd = find_extrema(steps_d, FH_d)
mx_Ed, my_Ed, mnx_Ed, mny_Ed = find_extrema(steps_d, FE_d)

print("Варіант = 21")
print(f"Довжина хвилі, lambda = {lambd} м")
print(f"Довжина стрижня, l = {l} м")
print(f"Коефіцієнт уповільнення, ksi = {round(ksi, 4)}")
print(f"Відстань між стрижнями, h = {h} м")

print("\nОдиночний стрижень:")
print(f"в площині H = {SGP_H_s:.2f}°")
print(f"в площині E = {SGP_E_s:.2f}°")
print(f"Перша бічна пелюстка H = {my_Hs[0]:.3f}, {mx_Hs[0]:.2f}°" if my_Hs else "")
print(f"Перша бічна пелюстка E = {my_Es[0]:.3f}, {mx_Es[0]:.2f}°" if my_Es else "")

print(f"\nДвострижнева антена:")
print(f"в площині H = {SGP_H_d:.2f}°")
print(f"в площині E = {SGP_E_d:.2f}°")
print(f"Перша бічна пелюстка H = {my_Hd[0]:.3f}, {mx_Hd[0]:.2f}°" if my_Hd else "")
print(f"Перша бічна пелюстка E = {my_Ed[0]:.3f}, {mx_Ed[0]:.2f}°" if my_Ed else "")

valid_mnx_Es = [x for x in mnx_Es if x > 0.5]
valid_mx_Es  = [(x, y) for x, y in zip(mx_Es, my_Es) if y > 0.01]

print("\nТабл. 1 - Нульові кути (площина E, одиночний)")
print("----------------")
print("| № |  θ°  |FE(θ)|")
print("----------------")
for i, x in enumerate(valid_mnx_Es[:3]):
    print(f"| {i+1} |{x:6.2f}|  0  |")
print("----------------")

print("\nТабл. 2 - Максимальні кути (площина E, одиночний)")
print("----------------")
print("| № |  θ°  |FE(θ)|")
print("----------------")
for i, (x, y) in enumerate(valid_mx_Es[:3]):
    print(f"| {i+1} |{x:6.2f}|{y:5.3f}|")
print("----------------")

def annotate_hp(ax, sgp, fs, side='right'):
    if not sgp:
        return
    dx = 3 if side == 'right' else -15
    ax.plot(sgp/2, fs, 'ro', markersize=4, label="Рівень половинної потужності")
    ax.annotate(f'({fs:.3f}, {sgp/2:.2f}°)',
                xy=(sgp / 2, fs),
                xytext=(sgp / 2 + dx, fs + 0.05),
                arrowprops=dict(arrowstyle='->', color='black'), fontsize=6)
    ax.plot([0, sgp / 2], [fs, fs], 'r--', linewidth=0.5)
    ax.plot([sgp / 2, sgp / 2], [0, fs], 'r--', linewidth=0.5)

def style_ax(ax, ylabel):
    ax.set_xlabel('θ°', fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_xticks(numpy.arange(0, 100, 2))
    ax.set_yticks(numpy.arange(0, 1.2, 0.1))
    ax.tick_params(labelsize=7)
    ax.set_ylim(-0.01, 1.01)
    ax.set_xlim(0, 90.5)
    ax.legend(loc="upper right", fontsize=7)
    ax.grid(which='both', linestyle='-', linewidth=0.2, color='gray')

fig, ax = plt.subplots(figsize=(20/2.54, 12/2.54))
ax.plot(steps_s, Fb_s, linewidth=0.7, label='$F_{б}(θ)$')
ax.plot(steps_s, FH_s, linewidth=1.2, label='$F_{H}(θ)$')
annotate_hp(ax, SGP_H_s, fS_H_s, 'right')
ax.plot(mx_Hs, my_Hs, 'o', markersize=4, color='black', label='$θ_{max}$')
ax.plot(mnx_Hs, mny_Hs, 'o', markersize=4, color='blue', label='$θ_{min}$')
style_ax(ax, '|FH(θ)|')
fig.tight_layout()
fig.savefig("ДС_одиночний_H.jpg", dpi=600, bbox_inches='tight')
plt.show()

fig, ax = plt.subplots(figsize=(20/2.54, 12/2.54))
ax.plot(steps_s, Fb_s, linewidth=0.7, label='$F_{б}(θ)$')
ax.plot(steps_s, [abs(numpy.cos(numpy.deg2rad(t))) for t in steps_s],
        linewidth=0.7, linestyle='-', label='$cos(θ)$')
ax.plot(steps_s, FE_s, linewidth=1.2, label='$F_{E}(θ)$')
annotate_hp(ax, SGP_E_s, fS_E_s, 'left')
ax.plot(mx_Es, my_Es, 'o', markersize=4, color='black', label='$θ_{max}$')
ax.plot(mnx_Es, mny_Es, 'o', markersize=4, color='blue', label='$θ_{min}$')
style_ax(ax, '|FE(θ)|')
fig.tight_layout()
fig.savefig("ДС_одиночний_E.jpg", dpi=600, bbox_inches='tight')
plt.show()

fig, ax = plt.subplots(figsize=(20/2.54, 12/2.54))
ax.plot(steps_d, Fb_d, linewidth=0.7, label='$F_{б}(θ)$')
ax.plot(steps_d, FC_d, linewidth=0.7, linestyle='-', label='$F_{C}(θ)$')
ax.plot(steps_d, FH_d, linewidth=1.2, label='$F_{H}(θ)$')
annotate_hp(ax, SGP_H_d, fS_H_d, 'right')
ax.plot(mx_Hd, my_Hd, 'o', markersize=4, color='black', label='$θ_{max}$')
ax.plot(mnx_Hd, mny_Hd, 'o', markersize=4, color='blue', label='$θ_{min}$')
style_ax(ax, '|FH(θ)|')
fig.tight_layout()
fig.savefig("ДС_двострижнева_H.jpg", dpi=600, bbox_inches='tight')
plt.show()

fig, ax = plt.subplots(figsize=(20/2.54, 12/2.54))
ax.plot(steps_d, Fb_d, linewidth=0.7, label='$F_{б}(θ)$')
ax.plot(steps_d, FC_d, linewidth=0.7, linestyle='-', label='$F_{C}(θ)$')
ax.plot(steps_d, FE_d, linewidth=1.2, label='$F_{E}(θ)$')
annotate_hp(ax, SGP_E_d, fS_E_d, 'left')
ax.plot(mx_Ed, my_Ed, 'o', markersize=4, color='black', label='$θ_{max}$')
ax.plot(mnx_Ed, mny_Ed, 'o', markersize=4, color='blue', label='$θ_{min}$')
style_ax(ax, '|FE(θ)|')
fig.tight_layout()
fig.savefig("ДС_двострижнева_E.jpg", dpi=600, bbox_inches='tight')
plt.show()