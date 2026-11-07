import math

EPS = 1e-8
H = 0.1
GRID_X = (-3.0, 3.0)
GRID_Y = (-3.0, 3.0)

def f(x, y):
    return x


def g(x, y):
    return y**2 - x**3 + x**4


def dg_dx(x, y):
    return -3*x**2 + 4*x**3


def dg_dy(x, y):
    return 2*y


def d2g_dx2(x, y):
    return -6*x + 12*x**2


def d2g_dxdy(x, y):
    return 0.0


def d2g_dy2(x, y):
    return 2.0



def is_zero(v, eps=EPS):
    return abs(v) < eps


def frange(a, b, step):
    out = []
    v = a
    while v <= b + 1e-12:
        out.append(round(v, 10))
        v += step
    if a <= 0.0 <= b:
        out.append(0.0)
    out = sorted(set(out))
    return out


def solve_linear(A, b):
    n = len(A)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-14:
            return None
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        for j in range(col, n + 1):
            M[col][j] /= pv
        for r in range(n):
            if r != col and abs(M[r][col]) > 0:
                fac = M[r][col]
                for j in range(col, n + 1):
                    M[r][j] -= fac * M[col][j]
    return [M[i][n] for i in range(n)]


# 1. НЕРЕГУЛЯРНЫЕ ТОЧКИ
def newton_irregular(x0, y0):
    x, y = x0, y0
    for _ in range(100):
        F1 = dg_dx(x, y)
        F2 = dg_dy(x, y)
        J = [[d2g_dx2(x, y), d2g_dxdy(x, y)],
             [d2g_dxdy(x, y), d2g_dy2(x, y)]]
        d = solve_linear(J, [-F1, -F2])
        if d is None:
            return None
        x += d[0]
        y += d[1]
        if abs(d[0]) + abs(d[1]) < 1e-12:
            break
    return (x, y)


def find_irregular_points():
    found = []

    def add_if_new(x, y):
        for (px, py) in found:
            if is_zero(x - px, 1e-5) and is_zero(y - py, 1e-5):
                return
        found.append((x, y))

    xs = frange(*GRID_X, H)
    ys = frange(*GRID_Y, H)

    for x0 in xs:
        for y0 in ys:
            grad_norm2 = dg_dx(x0, y0)**2 + dg_dy(x0, y0)**2
            if grad_norm2 < 1.0:
                sol = newton_irregular(x0, y0)
                if sol is None:
                    continue
                cx, cy = sol
                if (is_zero(dg_dx(cx, cy), 1e-6) and is_zero(dg_dy(cx, cy), 1e-6) and is_zero(g(cx, cy), 1e-6)):
                    add_if_new(cx, cy)
    return found


# 4–6. СИСТЕМА ЛАГРАНЖА
def dL_dx(x, y, lam):
    return 1.0 + lam * dg_dx(x, y)

def dL_dy(x, y, lam):
    return lam * dg_dy(x, y)

def dL_dlam(x, y, lam):
    return g(x, y)


def jacobian_L(x, y, lam):
    J11 = lam * d2g_dx2(x, y)
    J12 = lam * d2g_dxdy(x, y)
    J13 = dg_dx(x, y)
    J21 = lam * d2g_dxdy(x, y)
    J22 = lam * d2g_dy2(x, y)
    J23 = dg_dy(x, y)
    J31 = dg_dx(x, y)
    J32 = dg_dy(x, y)
    J33 = 0.0
    return [[J11, J12, J13],
            [J21, J22, J23],
            [J31, J32, J33]]


def newton_lagrange(x0, y0, lam0):
    x, y, lam = x0, y0, lam0
    for _ in range(200):
        F1 = dL_dx(x, y, lam)
        F2 = dL_dy(x, y, lam)
        F3 = dL_dlam(x, y, lam)
        J = jacobian_L(x, y, lam)
        d = solve_linear(J, [-F1, -F2, -F3])
        if d is None:
            return None
        x += d[0]; y += d[1]; lam += d[2]
        if abs(d[0]) + abs(d[1]) + abs(d[2]) < 1e-12:
            break
    return (x, y, lam)


def solve_lagrange_system():
    found = []

    def add_if_new(x, y, lam):
        for (px, py, pl) in found:
            if is_zero(x - px, 1e-5) and is_zero(y - py, 1e-5) \
               and is_zero(lam - pl, 1e-5):
                return
        found.append((x, y, lam))

    xs = frange(*GRID_X, H)
    ys = frange(*GRID_Y, H)

    for x0 in xs:
        for y0 in ys:
            for lam0 in (-1.0, 0.0, 1.0):
                sol = newton_lagrange(x0, y0, lam0)
                if sol is None:
                    continue
                x, y, lam = sol
                if (is_zero(dL_dx(x, y, lam), 1e-5) and is_zero(dL_dy(x, y, lam), 1e-5) and is_zero(dL_dlam(x, y, lam), 1e-5)):
                    add_if_new(x, y, lam)
    return found


# 7–10. ПРОВЕРКА, СРАВНЕНИЕ, ОТВЕТ
def main():
    print("=== Вариант 18: f(x,y)=x,  y^2 = x^3(1-x) ===\n")

    irregular = find_irregular_points()
    print("1. Нерегулярные точки:", irregular)

    stationary = solve_lagrange_system()

    regular = []
    for (x, y, lam) in stationary:
        if not (is_zero(dg_dx(x, y), 1e-5) and is_zero(dg_dy(x, y), 1e-5)):
            regular.append((x, y, lam))
    print("6. Регулярные стационарные точки (x, y, lam):", regular)

    candidates = []
    for (x, y, lam) in regular:
        if is_zero(g(x, y), 1e-5):
            candidates.append((x, y))
    for (x, y) in irregular:
        candidates.append((x, y))

    unique = []
    for p in candidates:
        if not any(is_zero(p[0]-q[0], 1e-5) and is_zero(p[1]-q[1], 1e-5)
                   for q in unique):
            unique.append(p)
    candidates = unique

    print("7. Допустимые кандидаты:", candidates)

    values = [(p, f(p[0], p[1])) for p in candidates]
    print("8. Значения f в кандидатах:")
    for p, v in values:
        print("   ", p, "->", v)

    if not values:
        print("Нет кандидатов.")
        return

    fmin = min(values, key=lambda t: t[1])
    fmax = max(values, key=lambda t: t[1])

    print("9. min f =", fmin[1], "в точке", fmin[0])
    print("   max f =", fmax[1], "в точке", fmax[0])
    print("10. Ответ:")
    print("    min f =", fmin[1], "в точке", fmin[0])
    print("    max f =", fmax[1], "в точке", fmax[0])


if __name__ == "__main__":
    main()
