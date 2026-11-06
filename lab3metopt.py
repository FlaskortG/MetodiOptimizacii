# 1. Определение нерегулярных точек
def g(x, y):
    return y**2 - x**3 + x**4

def dg_dx(x, y):
    return -3*x**2 + 4*x**3

def dg_dy(x, y):
    return 2*y

def is_zero(v, eps=1e-9):
    return abs(v) < eps

def find_irregular_points():
    pts = []
    for xv in [0.0, 0.75]:
        yv = 0.0
        if is_zero(dg_dx(xv, yv)) and is_zero(dg_dy(xv, yv)) and is_zero(g(xv, yv)):
            pts.append((xv, yv))
    return pts

irregular = find_irregular_points()
print("1. Нерегулярные точки:", irregular)


# 2. Целевая функция и ограничение
def f(x, y):
    return x

def constraint(x, y):
    return g(x, y)


# 4. Частные производные Лагранжа
def dL_dx(x, y, lam):
    return 1 + lam * (-3*x**2 + 4*x**3)

def dL_dy(x, y, lam):
    return lam * 2*y

def dL_dlam(x, y, lam):
    return y**2 - x**3 + x**4


# 5. Система уравнений
def solve_system():
    solutions = []
    for xv in [0.0, 1.0]:
        yv = 0.0
        if not is_zero(dL_dlam(xv, yv, 0.0)):
            continue
        coef = -3*xv**2 + 4*xv**3
        if is_zero(coef):
            continue
        lamv = -1.0 / coef
        if is_zero(dL_dy(xv, yv, lamv)):
            solutions.append((xv, yv, lamv))
    return solutions

stationary = solve_system()
print("6. Стационарные точки (x, y, lambda):", stationary)


# 7. Проверка ограничения
candidates = []
for (xv, yv, lamv) in stationary:
    if is_zero(constraint(xv, yv)):
        candidates.append((xv, yv))
for p in irregular:
    candidates.append(p)

unique = []
for p in candidates:
    if not any(is_zero(p[0]-q[0]) and is_zero(p[1]-q[1]) for q in unique):
        unique.append(p)
candidates = unique

print("7. Допустимые кандидаты:", candidates)


# 8. Характер точек: сравнение значений f
values = [(p, f(p[0], p[1])) for p in candidates]
print("8. Значения f в кандидатах:")
for p, v in values:
    print("   ", p, "->", v)

fmin = min(values, key=lambda t: t[1])
fmax = max(values, key=lambda t: t[1])
print("   Минимум:  точка", fmin[0], " f =", fmin[1])
print("   Максимум: точка", fmax[0], " f =", fmax[1])


# 9. Оптимальные значения целевой функции
print("9. min f =", fmin[1], "в точке", fmin[0])
print("   max f =", fmax[1], "в точке", fmax[0])


# 10. Ответ
print("10. Ответ:")
print("    min f =", fmin[1], "в точке", fmin[0])
print("    max f =", fmax[1], "в точке", fmax[0])