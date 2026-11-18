#!/usr/bin/env python3
"""Generates dev/ref_alaamee/*: random networks and, for every node, the change statistic of every effect as computed by the
Python ALAAMEE software (Stivala; MIT licence). ALAAMEE is NOT part of this repository; give the path of its python/ directory:
    python3 -I dev/gen_ref_alaamee.py <path to ALAAMEE>/python dev/ref_alaamee
Only the numbers it computes are recorded here, no ALAAMEE code. Files per case X:
  X_edges.csv (from,to), X_nodes.csv (y,c,k,b; missing = empty), X_set.csv (setting network), X_dist.csv (distance matrix, n rows),
  X_spec.csv (name,code,param,cov: the plugin code of each column of X_delta.csv, its parameter and covariate kind
  0 none 1 c 2 k 3 b 4 setting 5 distance), X_delta.csv (header = names), X_meta.csv (n,directed,n1)
"""
import sys, os, random, math
from functools import partial

src, out = sys.argv[1], sys.argv[2]
sys.path.insert(0, src)
os.makedirs(out, exist_ok=True)
import numpy as np
from Graph import Graph
from Digraph import Digraph
from BipartiteGraph import BipartiteGraph, MODE_A, MODE_B
import changeStatisticsALAAM as U
import changeStatisticsALAAMdirected as D
import changeStatisticsALAAMbipartite as B

random.seed(20261118)
tmp = os.path.join(out, "_tmp")
os.makedirs(tmp, exist_ok=True)


def na(p, v, nav):
    return nav if random.random() < p else v


def make_case(name, kind, n, p, n1=None, mode=None):
    """kind: und, dir, bip"""
    edges = []
    if kind == "bip":
        for i in range(n1):
            for j in range(n1, n):
                if random.random() < p: edges.append((i + 1, j + 1))
    else:
        for i in range(n):
            for j in (range(n) if kind == "dir" else range(i + 1, n)):
                if i != j and random.random() < p: edges.append((i + 1, j + 1))
    sedges = [(i + 1, j + 1) for i in range(n) for j in range(i + 1, n) if random.random() < 0.12]
    b = [na(0.08, random.randint(0, 1), -1) for _ in range(n)]
    k = [na(0.08, random.randint(1, 3), -1) for _ in range(n)]
    c = [na(0.08, round(random.gauss(0, 1), 4), float("nan")) for _ in range(n)]
    y = [1 if random.random() < 0.4 else 0 for _ in range(n)]
    if kind == "bip":
        for i in range(n):
            inmode = (i < n1) if mode == "A" else (i >= n1)
            if not inmode: y[i] = -1
    else:
        for i in range(n):
            if random.random() < 0.07: y[i] = -1
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            dist[i][j] = dist[j][i] = round(random.uniform(0.1, 2.0), 4)
    f = lambda s: os.path.join(tmp, name + s)
    with open(f(".net"), "w") as fh:
        if kind == "bip": fh.write("*Vertices %d %d\n" % (n, n1))
        else: fh.write("*Vertices %d\n" % n)
        for i in range(n): fh.write('%d "%d"\n' % (i + 1, i + 1))
        fh.write("*Arcs\n" if kind == "dir" else "*Edges\n")
        for e in edges: fh.write("%d %d\n" % e)
    with open(f("_set.net"), "w") as fh:
        fh.write("*Vertices %d\n" % n)
        for i in range(n): fh.write('%d "%d"\n' % (i + 1, i + 1))
        fh.write("*Edges\n")
        for e in sedges: fh.write("%d %d\n" % e)
    with open(f("_b.txt"), "w") as fh: fh.write("b\n" + "".join("%d\n" % v for v in b))
    with open(f("_k.txt"), "w") as fh: fh.write("k\n" + "".join("%d\n" % v for v in k))
    with open(f("_c.txt"), "w") as fh: fh.write("c\n" + "".join(("NA\n" if v != v else "%r\n" % v) for v in c))
    args = dict(binattr_filename=f("_b.txt"), contattr_filename=f("_c.txt"), catattr_filename=f("_k.txt"))
    if kind == "dir": G = Digraph(f(".net"), **args)
    elif kind == "bip": G = BipartiteGraph(f(".net"), **args)
    else: G = Graph(f(".net"), **args)
    GS = Graph(f("_set.net"))
    dm = np.array(dist)
    # effect list: (name, plugin code, param, cov kind, function (G, A, i))
    E = []
    def add(nm, code, par, cov, fn): E.append((nm, code, par, cov, fn))
    if kind == "und":
        add("density", 1, 0, 0, U.changeDensity); add("activity", 2, 0, 0, U.changeActivity)
        add("twostar", 3, 0, 0, U.changeTwoStar); add("threestar", 4, 0, 0, U.changeThreeStar)
        add("gwactivity", 5, 0.5, 0, partial(U.changeGWActivity, 0.5)); add("contagion", 6, 0, 0, U.changeContagion)
        add("alter2star1", 7, 0, 0, U.changePartnerActivityTwoPath); add("alter2star2", 8, 0, 0, U.changeIndirectPartnerAttribute)
        add("partneractivity", 9, 0, 0, U.changePartnerAttributeActivity); add("partnerpartner", 10, 0, 0, U.changePartnerPartnerAttribute)
        add("triangle1", 11, 0, 0, U.changeTriangleT1); add("triangle2", 12, 0, 0, U.changeTriangleT2); add("triangle3", 13, 0, 0, U.changeTriangleT3)
        add("gwcontagion", 14, 0.5, 0, partial(U.changeGWContagion, 0.5)); add("logcontagion", 15, 0, 0, U.changeLogContagion)
        add("powercontagion", 16, 2.0, 0, partial(U.changePowerContagion, 2.0))
        add("ego_c", 17, 0, 1, partial(U.changeoOc, "c")); add("ego_b", 17, 0, 3, partial(U.changeoOb, "b"))
        add("alter_c", 18, 0, 1, partial(U.changeo_Oc, "c")); add("alter_b", 18, 0, 3, partial(U.changeo_Ob, "b"))
        add("same_k", 19, 0, 2, partial(U.changeoO_Osame, "k")); add("diff_k", 20, 0, 2, partial(U.changeoO_Odiff, "k"))
        add("samecontagion_k", 21, 0, 2, partial(U.changeoO_OsameContagion, "k"))
        add("samepartneractivity_k", 22, 0, 2, partial(U.changeSamePartnerActivityTwoPath, "k"))
        add("diffpartneractivity_k", 23, 0, 2, partial(U.changeDiffPartnerActivityTwoPath, "k"))
        add("sameindirect_k", 24, 0, 2, partial(U.changeSameIndirectPartnerAttribute, "k"))
        add("diffindirect_k", 25, 0, 2, partial(U.changeDiffIndirectPartnerAttribute, "k"))
        add("alterbinarytwostar1_b", 26, 0, 3, partial(U.changeAlterBinaryTwoStar1, "b"))
        add("alterbinarytwostar2_b", 27, 0, 3, partial(U.changeAlterBinaryTwoStar2, "b"))
        add("settinghomophily", 40, 0, 4, partial(U.changeSettingHomophily, GS))
        add("geographichomophily", 41, 0, 5, partial(U.changeGeographicHomophily, dm))
        add("contagiondist", 42, 0, 5, partial(U.changeContagionDist, dm))
    elif kind == "dir":
        add("density", 1, 0, 0, U.changeDensity)
        add("ego_c", 17, 0, 1, partial(U.changeoOc, "c")); add("ego_b", 17, 0, 3, partial(U.changeoOb, "b"))
        for nm, code, fn in [("sender", 102, D.changeSender), ("receiver", 103, D.changeReceiver), ("reciprocity", 104, D.changeReciprocity),
                ("instar2", 105, D.changeEgoInTwoStar), ("instar3", 106, D.changeEgoInThreeStar), ("outstar2", 107, D.changeEgoOutTwoStar),
                ("outstar3", 108, D.changeEgoOutThreeStar), ("mixedstar", 109, D.changeMixedTwoStar), ("mixedsource", 110, D.changeMixedTwoStarSource),
                ("mixedsink", 111, D.changeMixedTwoStarSink), ("contagion", 112, D.changeContagion), ("contagionrecip", 113, D.changeContagionReciprocity),
                ("transitive1", 114, D.changeTransitiveTriangleT1), ("transitive3", 115, D.changeTransitiveTriangleT3),
                ("transitived1", 116, D.changeTransitiveTriangleD1), ("transitiveu1", 117, D.changeTransitiveTriangleU1),
                ("cyclic1", 118, D.changeCyclicTriangleC1), ("cyclic3", 119, D.changeCyclicTriangleC3),
                ("alterinstar2", 120, D.changeAlterInTwoStar2), ("alteroutstar2", 121, D.changeAlterOutTwoStar2),
                ("logcontagion", 126, D.changeLogContagion)]:
            add(nm, code, 0, 0, fn)
        add("gwsender", 122, 0.5, 0, partial(D.changeGWSender, 0.5)); add("gwreceiver", 123, 0.5, 0, partial(D.changeGWReceiver, 0.5))
        add("gwcontagion", 125, 0.5, 0, partial(D.changeGWContagion, 0.5)); add("powercontagion", 127, 2.0, 0, partial(D.changePowerContagion, 2.0))
        for nm, code, fn in [("sendermatch_k", 131, D.changeSenderMatch), ("sendermismatch_k", 132, D.changeSenderMismatch),
                ("receivermatch_k", 133, D.changeReceiverMatch), ("receivermismatch_k", 134, D.changeReceiverMismatch),
                ("reciprocitymatch_k", 135, D.changeReciprocityMatch), ("reciprocitymismatch_k", 136, D.changeReciprocityMismatch)]:
            add(nm, code, 0, 2, partial(fn, "k"))
        add("settinghomophily", 40, 0, 4, partial(U.changeSettingHomophily, GS))
        add("geographichomophily", 41, 0, 5, partial(U.changeGeographicHomophily, dm))
    else:
        m = MODE_A if mode == "A" else MODE_B
        add("density", 1, 0, 0, partial(B.changeBipartiteDensity, m)); add("activity", 2, 0, 0, partial(B.changeBipartiteActivity, m))
        add("twostar", 3, 0, 0, partial(B.changeBipartiteEgoTwoStar, m)); add("threestar", 4, 0, 0, partial(B.changeBipartiteEgoThreeStar, m))
        add("gwactivity", 5, 0.5, 0, partial(B.changeBipartiteGWActivity, m, 0.5))
        add("alter2star1", 7, 0, 0, partial(B.changeBipartiteAlterTwoStar1, m)); add("alter2star2", 8, 0, 0, partial(B.changeBipartiteAlterTwoStar2, m))
        add("fourcycle1", 28, 0, 0, partial(B.changeBipartiteFourCycle1, m)); add("fourcycle2", 29, 0, 0, partial(B.changeBipartiteFourCycle2, m))
        add("ego_c", 17, 0, 1, partial(U.changeoOc, "c")); add("ego_b", 17, 0, 3, partial(U.changeoOb, "b"))
        add("alter_c", 18, 0, 1, partial(U.changeo_Oc, "c")); add("alter_b", 18, 0, 3, partial(U.changeo_Ob, "b"))
        add("samepartneractivity_k", 22, 0, 2, partial(B.changeBpAlterSameTwoStar1, m, "k"))
        add("diffpartneractivity_k", 23, 0, 2, partial(B.changeBpAlterDiffTwoStar1, m, "k"))
        add("sameindirect_k", 24, 0, 2, partial(B.changeBpAlterSameTwoStar2, m, "k"))
        add("diffindirect_k", 25, 0, 2, partial(B.changeBpAlterDiffTwoStar2, m, "k"))
        add("alterbinarytwostar1_b", 26, 0, 3, partial(B.changeBpAlterBinaryTwoStar1, m, "b"))
        add("alterbinarytwostar2_b", 27, 0, 3, partial(B.changeBpAlterBinaryTwoStar2, m, "b"))
        add("settinghomophily", 40, 0, 4, partial(U.changeSettingHomophily, GS))
        add("geographichomophily", 41, 0, 5, partial(U.changeGeographicHomophily, dm))
    # change statistics: A with A_i = 0, others as observed (missing counts as not 1)
    rows = []
    for i in range(n):
        if y[i] < 0 and kind == "bip": rows.append(None); continue
        A = np.array(y, dtype=int)
        A[i] = 0
        rows.append([float(fn(G, A, i)) for (_, _, _, _, fn) in E])
    w = lambda s: open(os.path.join(out, name + s), "w")
    with w("_edges.csv") as fh: fh.write("from,to\n" + "".join("%d,%d\n" % e for e in edges))
    with w("_nodes.csv") as fh:
        fh.write("y,c,k,b\n")
        for i in range(n):
            fh.write("%s,%s,%s,%s\n" % ("" if y[i] < 0 else y[i], "" if c[i] != c[i] else repr(c[i]), "" if k[i] < 0 else k[i], "" if b[i] < 0 else b[i]))
    with w("_set.csv") as fh: fh.write("from,to\n" + "".join("%d,%d\n" % e for e in sedges))
    with w("_dist.csv") as fh: fh.write("\n".join(",".join(repr(v) for v in r) for r in dist) + "\n")
    with w("_spec.csv") as fh:
        fh.write("name,code,param,cov\n" + "".join("%s,%d,%r,%d\n" % (nm, code, par, cov) for (nm, code, par, cov, _) in E))
    with w("_delta.csv") as fh:
        fh.write(",".join(e[0] for e in E) + "\n")
        for r in rows: fh.write(",".join("." for _ in E) if r is None else ",".join(repr(v) for v in r)); fh.write("\n")
    with w("_meta.csv") as fh:
        fh.write("n,directed,n1\n%d,%d,%d\n" % (n, int(kind == "dir"), n1 or 0))
    print(name, kind, n, len(edges), "ties,", len(E), "effects")


make_case("x_und_a", "und", 28, 0.2)
make_case("x_und_b", "und", 17, 0.35)
make_case("x_dir_a", "dir", 22, 0.15)
make_case("x_dir_b", "dir", 15, 0.3)
make_case("x_bip_b", "bip", 22, 0.25, n1=13, mode="B")
make_case("x_bip_a", "bip", 20, 0.3, n1=8, mode="A")
