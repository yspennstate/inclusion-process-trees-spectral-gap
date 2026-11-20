"""Check of the S(2,2,2) separation in Lemma 9.

SIP_k on a weighted graph: a particle at x jumps to y at rate w_xy * eta_x * (a_y + eta_y). Reversible measure
pi(eta) = prod_x Gamma(a_x + eta_x) / (Gamma(a_x) eta_x!). We symmetrise -L_k with pi^{1/2} and compute all
eigenvalues at 50 significant digits with mpmath, for k = 1 and k = 2, and print gap_1, gap_2 and their ratio.
Also prints the exact count of eigenvalues below t = 1/200000 by an LDL^T inertia count in rational arithmetic
(sympy Rational).

usage: python3 check_spider.py
"""
import itertools

import mpmath as mp
from sympy import Rational as Q

mp.mp.dps = 50

# spider S(2,2,2): centre 0, arms 0-1-2, 0-3-4, 0-5-6
EDGES = [(0, 1), (1, 2), (0, 3), (3, 4), (0, 5), (5, 6)]
A = [Q(1, 10000), Q(1, 10000), Q(1, 40), Q(1, 10000), Q(1, 40), Q(1, 10000), Q(1, 40)]
W = {(0, 1): 16, (0, 3): 16, (0, 5): 16, (1, 2): 1, (3, 4): 1, (5, 6): 1}


def configs(n, k):
    out = []
    for c in itertools.combinations_with_replacement(range(n), k):
        e = [0] * n
        for x in c:
            e[x] += 1
        out.append(tuple(e))
    return out


def pi_weight(e):
    p = Q(1)
    for x, m in enumerate(e):
        for j in range(m):
            p *= (A[x] + j) / (j + 1)
    return p


def generator(k):
    n = len(A)
    S = configs(n, k)
    ix = {s: i for i, s in enumerate(S)}
    d = len(S)
    G = [[Q(0)] * d for _ in range(d)]          # G = -L_k as a rate matrix acting on functions
    for s in S:
        i = ix[s]
        for (x, y), w in W.items():
            for u, v in ((x, y), (y, x)):
                if s[u]:
                    t = list(s); t[u] -= 1; t[v] += 1
                    r = w * s[u] * (A[v] + s[v])
                    G[i][i] += r
                    G[i][ix[tuple(t)]] -= r
    return S, G


def spectrum(k):
    S, G = generator(k)
    P = [pi_weight(s) for s in S]
    d = len(S)
    M = mp.matrix(d, d)
    for i in range(d):
        for j in range(d):
            if G[i][j] != 0:
                M[i, j] = mp.sqrt(mp.mpf(P[i].p) / P[i].q) * mp.mpf(G[i][j].p) / G[i][j].q \
                    / mp.sqrt(mp.mpf(P[j].p) / P[j].q)
    asym = max(abs(M[i, j] - M[j, i]) for i in range(d) for j in range(d))
    ev = sorted(mp.eigsy(M, eigvals_only=True))
    return ev, asym


def inertia_below(k, t):
    S, G = generator(k)
    P = [pi_weight(s) for s in S]
    d = len(S)
    B = [[P[i] * G[i][j] for j in range(d)] for i in range(d)]
    for i in range(d):
        for j in range(d):
            assert B[i][j] == B[j][i]
        B[i][i] -= t * P[i]
    neg = 0
    for c in range(d):
        p = B[c][c]
        assert p != 0, 'zero pivot'
        if p < 0:
            neg += 1
        for i in range(c + 1, d):
            if B[i][c] != 0:
                f = B[i][c] / p
                for j in range(c + 1, d):
                    B[i][j] -= f * B[c][j]
    return neg


def dist_table():
    n = len(A)
    adj = {x: set() for x in range(n)}
    for x, y in EDGES:
        adj[x].add(y); adj[y].add(x)
    D = {}
    for s in range(n):
        seen = {s: 0}; fr = [s]
        while fr:
            nf = []
            for u in fr:
                for v in adj[u]:
                    if v not in seen:
                        seen[v] = seen[u] + 1; nf.append(v)
            fr = nf
        for t, d in seen.items():
            D[s, t] = d
    return adj, D


def lam22(ahat):
    """Bottom of -M for two independent RW(ahat) walkers (rate c_xz ahat_z) on unordered pairs at distance >= 2,
    killed on any move to distance <= 1; symmetric in L2(ahat x ahat). Also returns gap_RW(ahat)."""
    n = len(A)
    adj, D = dist_table()
    c = {}
    for (x, y), w in W.items():
        c[x, y] = c[y, x] = w
    pairs = [(x, y) for x in range(n) for y in range(x + 1, n) if D[x, y] >= 2]
    ix = {p: i for i, p in enumerate(pairs)}
    m = len(pairs)
    M = mp.matrix(m, m)
    for (x, y) in pairs:
        i = ix[x, y]
        for mover, other in ((x, y), (y, x)):
            for z in adj[mover]:
                r = mp.mpf(c[mover, z]) * ahat[z]
                M[i, i] += r
                if D[z, other] >= 2:
                    q = (min(z, other), max(z, other))
                    # symmetrise with sqrt of ahat x ahat
                    M[i, ix[q]] -= r * mp.sqrt(ahat[x] * ahat[y]) / mp.sqrt(ahat[q[0]] * ahat[q[1]])
    ev2 = sorted(mp.eigsy(M, eigvals_only=True))
    R = mp.matrix(n, n)
    for x in range(n):
        for z in adj[x]:
            r = mp.mpf(c[x, z]) * ahat[z]
            R[x, x] += r
            R[x, z] -= r * mp.sqrt(ahat[x]) / mp.sqrt(ahat[z])
    ev1 = sorted(mp.eigsy(R, eigvals_only=True))
    return ev2[0], ev1[1], len(pairs)


if __name__ == '__main__':
    ahat = [mp.mpf(a.p) / a.q * 10000 for a in A]
    l22, grw, npairs = lam22(ahat)
    print('small-weight limit, ahat = a/1e-4: pairs', npairs, 'lambda_22 =', mp.nstr(l22, 15), 'gap_RW =',
          mp.nstr(grw, 15), 'ratio', mp.nstr(l22 / grw, 15))
    e1, s1 = spectrum(1)
    e2, s2 = spectrum(2)
    print('k=1 dim', len(e1), 'symmetry defect', mp.nstr(s1, 3), 'lowest', [mp.nstr(v, 12) for v in e1[:3]])
    print('k=2 dim', len(e2), 'symmetry defect', mp.nstr(s2, 3), 'lowest', [mp.nstr(v, 12) for v in e2[:4]])
    g1, g2 = e1[1], e2[1]
    print('gap_1 =', mp.nstr(g1, 20))
    print('gap_2 =', mp.nstr(g2, 20))
    print('gap_2/gap_1 =', mp.nstr(g2 / g1, 20))
    t = Q(1, 200000)
    print('exact counts below t=1/200000: k=1', inertia_below(1, t), ' k=2', inertia_below(2, t))
