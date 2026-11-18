#!/usr/bin/env python3
"""Generates dev/ref/*.csv: random networks and, for every node, the change statistic computed BY DEFINITION: the statistic z(y)
is counted over the whole outcome vector from its definition (a count of configurations), and the change statistic of node i is
z(y with y_i = 1) - z(y with y_i = 0). This is independent of the plugin's formulas (which are the local, incremental
formulas of the ALAAM literature). Used once to make the fixtures that cscripts/test_nwalaam_changestats.do compares the
plugin with. Usage: python3 dev/gen_ref_definitions.py dev/ref
"""
import sys, random, math, itertools

out = sys.argv[1]
random.seed(20261118)


def z_und(name, par, A, y, c, k, b):
    n = len(y)
    N = [[j for j in range(n) if A[i][j]] for i in range(n)]
    d = [len(x) for x in N]
    ones = [i for i in range(n) if y[i]]
    dplus = [sum(y[v] for v in N[j]) for j in range(n)]
    tri = [t for t in itertools.combinations(range(n), 3) if A[t[0]][t[1]] and A[t[0]][t[2]] and A[t[1]][t[2]]]
    tp = lambda i, v: sum(1 for u in range(n) if A[i][u] and A[v][u])
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if A[i][j]]
    if name == "density": return sum(y)
    if name == "activity": return sum(d[i] for i in ones)
    if name == "twostar": return sum(d[i] * (d[i] - 1) / 2 for i in ones)
    if name == "threestar": return sum(d[i] * (d[i] - 1) * (d[i] - 2) / 6 for i in ones)
    if name == "gwactivity": return sum(math.exp(-par * d[i]) for i in ones)
    if name == "contagion": return sum(1 for i, j in pairs if y[i] and y[j])
    if name == "alter2star1": return sum(sum(d[v] - 1 for v in N[i]) for i in ones)
    if name == "alter2star2": return sum(tp(i, v) for i in ones for v in ones if i < v)
    if name == "partneractivity": return sum(d[i] + d[j] - 2 for i, j in pairs if y[i] and y[j])
    if name == "partnerpartner":
        s = 0
        for cc in range(n):
            nb = [v for v in N[cc] if True]
            for a, b in itertools.combinations(N[cc], 2):
                if y[a] and y[cc] and y[b]: s += 2
        return s
    if name == "triangle1": return sum(y[a] + y[b] + y[c_] for a, b, c_ in tri)
    if name == "triangle2": return sum(math.comb(y[a] + y[b] + y[c_], 2) for a, b, c_ in tri)
    if name == "triangle3": return sum(1 for a, b, c_ in tri if y[a] and y[b] and y[c_])
    if name == "gwcontagion": return sum(math.exp(-par * dplus[j]) for j in ones)
    if name == "logcontagion": return sum(math.log(dplus[j] + 1) for j in ones)
    if name == "powercontagion": return sum(dplus[j] ** (1.0 / par) for j in ones)
    if name == "ego_c": return sum(c[i] for i in ones)
    if name == "alter_c": return sum(c[u] for i in ones for u in N[i])
    if name == "same_k": return sum(1 for i in ones for u in N[i] if k[u] == k[i])
    if name == "diff_k": return sum(1 for i in ones for u in N[i] if k[u] != k[i])
    if name == "samecontagion_k": return sum(1 for i, j in pairs if y[i] and y[j] and k[i] == k[j])
    # two-paths i-u-v (v != i) between the outcome node i and any v: with the same / a different category of i and v
    paths = [(i, u, v) for i in range(n) for u in N[i] for v in N[u] if v != i]
    if name == "samepartneractivity_k": return sum(1 for i, u, v in paths if y[i] and k[i] == k[v])
    if name == "diffpartneractivity_k": return sum(1 for i, u, v in paths if y[i] and k[i] != k[v])
    if name == "sameindirect_k": return sum(1 for i, u, v in paths if y[i] and y[v] and i < v and k[i] == k[v])
    if name == "diffindirect_k": return sum(1 for i, u, v in paths if y[i] and y[v] and i < v and k[i] != k[v])
    if name == "alterbinarytwostar1_b": return sum(1 for i, u, v in paths if y[i] and b[u])
    if name == "alterbinarytwostar2_b": return sum(1 for i, u, v in paths if y[i] and b[u] and y[v] and i < v)
    raise KeyError(name)


def z_dir(name, par, A, y, c, k, b):
    n = len(y)
    out_ = [[j for j in range(n) if A[i][j]] for i in range(n)]
    in_ = [[j for j in range(n) if A[j][i]] for i in range(n)]
    o = [len(x) for x in out_]; ii = [len(x) for x in in_]
    ones = [i for i in range(n) if y[i]]
    rec = [sum(1 for u in out_[i] if A[u][i]) for i in range(n)]
    # transitive triples a->b, b->c, a->c and directed 3-cycles a->b->c->a
    trans = [(a, b, cc) for a in range(n) for b in out_[a] for cc in out_[b] if cc != a and A[a][cc]]
    cyc = [(a, b, cc) for a in range(n) for b in out_[a] for cc in out_[b] if cc != a and A[cc][a] and a < b and a < cc]  # each cycle once (a = smallest node)
    if name == "density": return sum(y)
    if name == "sender": return sum(o[i] for i in ones)
    if name == "receiver": return sum(ii[i] for i in ones)
    if name == "reciprocity": return sum(rec[i] for i in ones)
    if name == "instar2": return sum(ii[i] * (ii[i] - 1) / 2 for i in ones)
    if name == "instar3": return sum(ii[i] * (ii[i] - 1) * (ii[i] - 2) / 6 for i in ones)
    if name == "outstar2": return sum(o[i] * (o[i] - 1) / 2 for i in ones)
    if name == "outstar3": return sum(o[i] * (o[i] - 1) * (o[i] - 2) / 6 for i in ones)
    if name == "mixedstar": return sum(ii[i] * o[i] - rec[i] for i in ones)
    if name == "mixedsource": return sum(sum(o[u] - (1 if A[u][i] else 0) for u in out_[i]) for i in ones)
    if name == "mixedsink": return sum(sum(ii[u] - (1 if A[i][u] else 0) for u in in_[i]) for i in ones)
    if name == "contagion": return sum(1 for i in range(n) for u in out_[i] if y[i] and y[u])
    if name == "contagionrecip": return sum(1 for i in range(n) for u in out_[i] if i < u and A[u][i] and y[i] and y[u])
    if name == "transitive1": return sum(y[b] for a, b, cc in trans)
    if name == "transitived1": return sum(y[a] for a, b, cc in trans)
    if name == "transitiveu1": return sum(y[cc] for a, b, cc in trans)
    if name == "transitive3": return sum(1 for a, b, cc in trans if y[a] and y[b] and y[cc])
    if name == "cyclic1": return sum(y[a] + y[b] + y[cc] for a, b, cc in cyc)
    if name == "cyclic3": return sum(1 for a, b, cc in cyc if y[a] and y[b] and y[cc])
    if name == "alterinstar2": return sum(math.comb(sum(y[v] for v in out_[u]), 2) for u in range(n))
    if name == "alteroutstar2": return sum(math.comb(sum(y[v] for v in in_[u]), 2) for u in range(n))
    if name == "gwsender": return sum(math.exp(-par * o[i]) for i in ones)
    if name == "gwreceiver": return sum(math.exp(-par * ii[i]) for i in ones)
    if name == "ego_c": return sum(c[i] for i in ones)
    if name == "alter_c": return sum(sum(c[u] for u in out_[i]) + sum(c[u] for u in in_[i]) for i in ones)
    # n_i: the number of ties of i to nodes with the outcome (a mutual tie counts twice, as in contagion)
    nn = [sum(y[u] for u in out_[i]) + sum(y[u] for u in in_[i]) for i in range(n)]
    if name == "gwcontagion": return sum(math.exp(-par * nn[i]) for i in ones)
    if name == "logcontagion": return sum(math.log(nn[i] + 1) for i in ones)
    if name == "powercontagion": return sum(nn[i] ** (1.0 / par) for i in ones)
    if name == "same_k": return sum(1 for i in ones for u in out_[i] + in_[i] if k[u] == k[i])
    if name == "diff_k": return sum(1 for i in ones for u in out_[i] + in_[i] if k[u] != k[i])
    if name == "samecontagion_k": return sum(1 for i in range(n) for u in out_[i] if y[i] and y[u] and k[i] == k[u])
    if name == "sendermatch_k": return sum(1 for i in ones for u in out_[i] if k[u] == k[i])
    if name == "sendermismatch_k": return sum(1 for i in ones for u in out_[i] if k[u] != k[i])
    if name == "receivermatch_k": return sum(1 for i in ones for u in in_[i] if k[u] == k[i])
    if name == "receivermismatch_k": return sum(1 for i in ones for u in in_[i] if k[u] != k[i])
    if name == "reciprocitymatch_k": return sum(1 for i in ones for u in out_[i] if A[u][i] and k[u] == k[i])
    if name == "reciprocitymismatch_k": return sum(1 for i in ones for u in out_[i] if A[u][i] and k[u] != k[i])
    raise KeyError(name)


def write_case(name, n, directed, p, effects):
    A = [[0] * n for _ in range(n)]
    edges = []
    for i in range(n):
        for j in (range(n) if directed else range(i + 1, n)):
            if i != j and random.random() < p:
                A[i][j] = 1
                if not directed: A[j][i] = 1
                edges.append((i + 1, j + 1))
    y = [1 if random.random() < 0.4 else 0 for _ in range(n)]
    c = [round(random.gauss(0, 1), 4) for _ in range(n)]
    k = [random.randint(1, 3) for _ in range(n)]
    b = [random.randint(0, 1) for _ in range(n)]
    z = z_dir if directed else z_und
    with open("%s/%s_edges.csv" % (out, name), "w") as f:
        f.write("from,to\n" + "".join("%d,%d\n" % e for e in edges))
    with open("%s/%s_nodes.csv" % (out, name), "w") as f:
        f.write("y,c,k,b\n" + "".join("%d,%s,%d,%d\n" % (y[i], c[i], k[i], b[i]) for i in range(n)))
    with open("%s/%s_delta.csv" % (out, name), "w") as f:
        f.write(",".join(e[0] if e[1] is None else "%s_%s" % (e[0], e[1]) for e in effects) + "\n")
        for i in range(n):
            y1 = list(y); y1[i] = 1
            y0 = list(y); y0[i] = 0
            row = []
            for nm, par in effects:
                base = nm
                row.append(repr(float(z(base, par, A, y1, c, k, b) - z(base, par, A, y0, c, k, b))))
            f.write(",".join(row) + "\n")
    print(name, n, len(edges), "ties,", sum(y), "ones")


und = [("density", None), ("activity", None), ("twostar", None), ("threestar", None), ("gwactivity", 0.5), ("contagion", None),
       ("alter2star1", None), ("alter2star2", None), ("partneractivity", None), ("partnerpartner", None), ("triangle1", None),
       ("triangle2", None), ("triangle3", None), ("gwcontagion", 0.5), ("logcontagion", None), ("powercontagion", 2.0),
       ("ego_c", None), ("alter_c", None), ("same_k", None), ("diff_k", None), ("samecontagion_k", None),
       ("samepartneractivity_k", None), ("diffpartneractivity_k", None), ("sameindirect_k", None), ("diffindirect_k", None),
       ("alterbinarytwostar1_b", None), ("alterbinarytwostar2_b", None)]
dire = [("density", None), ("sender", None), ("receiver", None), ("reciprocity", None), ("instar2", None), ("instar3", None),
        ("outstar2", None), ("outstar3", None), ("mixedstar", None), ("mixedsource", None), ("mixedsink", None), ("contagion", None),
        ("contagionrecip", None), ("transitive1", None), ("transitive3", None), ("transitived1", None), ("transitiveu1", None),
        ("cyclic1", None), ("cyclic3", None), ("alterinstar2", None), ("alteroutstar2", None), ("gwsender", 0.5), ("gwreceiver", 0.5),
        ("ego_c", None), ("alter_c", None),
        ("gwcontagion", 0.5), ("logcontagion", None), ("powercontagion", 2.0), ("same_k", None), ("diff_k", None), ("samecontagion_k", None),
        ("sendermatch_k", None), ("sendermismatch_k", None), ("receivermatch_k", None), ("receivermismatch_k", None),
        ("reciprocitymatch_k", None), ("reciprocitymismatch_k", None)]
write_case("und_a", 26, False, 0.2, und)
write_case("und_b", 18, False, 0.35, und)
write_case("dir_a", 22, True, 0.15, dire)
write_case("dir_b", 16, True, 0.3, dire)
