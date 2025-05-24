from eznf import modeler
import argparse
import itertools
import math


def encode(n, g, forced_sym=False):
    """
    Encode the existence of N points in the plane (general position) without g-gons.
    Args:
        N (int): number of points
        g (int): number of sides forbidden gon
        forced_sym (bool, int): forces k-fold symmetry when provided a positive integer.
    Returns:
        modeler.Modeler: encoded model"""
    enc = modeler.Modeler()
    for (p, q, r) in itertools.combinations(range(n), 3):
        enc.add_var(f"cc_{p, q, r}")

    # linear order on the points
    for (p, q) in itertools.permutations(range(n), 2):
        enc.add_var(f"<_{p, q}")

    # non-degeneracy of the order
    for (p, q) in itertools.combinations(range(n), 2):
        enc.add_clause([f"<_{p, q}", f"<_{q, p}"])
        enc.add_clause([f"-<_{p, q}", f"-<_{q, p}"])

    # transitivity <
    for (p, q, r) in itertools.permutations(range(n), 3):
        enc.add_clause([f"-<_{p, q}", f"-<_{q, r}", f"<_{p, r}"])
        enc.add_clause([f"<_{p, q}", f"<_{q, r}", f"-<_{p, r}"])

    # cyclic symmetry and anti-symmetry
    def cc(p, q, r):
        res = tuple(sorted((p, q, r)))
        sgn = (-1) ** (len([(u, v) for u, v in itertools.combinations((p, q, r), 2) if u > v]))
        return sgn*enc.v(f"cc_{res}")
        
    # Enforce symmetry
    sym_map = {}
    if forced_sym:
        layers = [list(range(forced_sym*i, forced_sym*(i+1))) for i in range(n // forced_sym)]
        for i in range(n % forced_sym):
            layers.append([n - (n % forced_sym) + i])
        
        for layer in layers:
            for j, el in enumerate(layer):
                sym_map[el] = layer[(j+1) % len(layer)]

        print(sym_map)

        # for tri in itertools.combinations(range(n), 3):
        #     sym_tri = (sym_map[tri[0]], sym_map[tri[1]], sym_map[tri[2]])
        #     enc.add_var_equivalence(f"cc_{tri}", f"cc_{sym_tri}") TODO: add back once allsat-script supports this.

        for tri in itertools.combinations(range(n), 3):
            sym_tri = (sym_map[tri[0]], sym_map[tri[1]], sym_map[tri[2]])
            enc.add_clause([-cc(*tri), cc(*sym_tri)])
            enc.add_clause([-cc(*sym_tri), cc(*tri)])

    print("Clauses for ordered signotope axioms...")
    cnt = 0
    for (p, q, r, s) in itertools.permutations(range(n), 4):
        cnt += 1
        if cnt % 100 == 0:
            print(f"at {cnt}/{math.perm(n, 4)}.... ({cnt/math.perm(n, 4)*100:.2f}%)", end="\r")
        if q < s and r < s:
            enc.add_clause([f"-<_{p, q}", f"-<_{p, r}", f"-<_{p, s}",  cc(p,q,r), -cc(p,q,s),  cc(p,r,s)])
        enc.add_clause([f"-<_{p, r}", f"-<_{q, r}", f"-<_{r, s}",  cc(p,q,r), -cc(p,r,s),  cc(q,r,s)])
    print('\n')

    # convex quadrilaterals
    cnt = 0
    for (p, q, r, s) in itertools.combinations(range(n), 4):
        cnt += 1
        if cnt % 100 == 0:
            print(f"at {cnt}/{math.comb(n, 4)}.... ({cnt/math.comb(n, 4)*100:.2f}%)", end="\r")
        enc.add_var(f"conv_{p,q,r,s}")
        enc.constraint(f"conv_{p,q,r,s} <-> ((cc_{p, q, r} <-> cc_{p, r, s}) <-> (cc_{p, q, s} <-> cc_{q, r, s}))")
    print('\n')


    def lex_smallest_rot(seq):
        def nxt(seq):
            return tuple(sym_map[i] for i in seq)
        rots = [seq]
        while True:
            nx = nxt(rots[-1])
            if nx == rots[0]:
                return min(rots)
            rots.append(nx)

    # no g-gons
    cnt = 0
    for comb in itertools.combinations(range(n), g):
        cnt += 1
        if cnt % 100 == 0:
            print(f"at {cnt}/{math.comb(n, g)}.... ({cnt/math.comb(n, g)*100:.2f}%)", end="\r")
        if forced_sym and lex_smallest_rot(comb) != comb:
            continue
        enc.add_clause([f"-conv_{qd}" for qd in itertools.combinations(comb, 4)])
    print('\n')

    # Convex hull structure
    if forced_sym:
        print(f"layers: {layers}")
        for i in range(len(layers)-1):
            if len(layers[i]) <= 2:
                continue
            for j in layers[i+1]:
                for k in range(len(layers[i])):
                    enc.add_clause([cc(layers[i][k], layers[i][(k+1)%(len(layers[i]))], j)])
                    
        for layer in layers:
            if len(layer) <= 2:
                continue
            enc.add_clause([cc(layer[j], layer[(j+1)%(len(layer))], layer[(j+2)%(len(layer))]) for j in range(len(layer))])

    # same quadrant constraint
    if forced_sym == 4:
        for i in range(n):
            if i and i % 4 == 0:
                enc.add_clause([-cc(0, 2, i)])
                enc.add_clause([cc(1, 3, i)])

    if forced_sym and n % forced_sym == 1:
        for i in range(n-1):
            if i and i % forced_sym == 0:
                enc.add_clause([ cc(0, i, n-1)])
                enc.add_clause([-cc(1, i, n-1)])
    return enc

def main():
    argparser = argparse.ArgumentParser(description="Encode the Erdos-Szekeres g-gon problem")
    argparser.add_argument("-n", "--n", type=int, default=16, help="Number of vertices")
    argparser.add_argument("-g", "--g", type=int, default=6, help="Number of sides")
    argparser.add_argument("-s", "--sym", type=int, default=0, help="Force k-fold symmetry")

    n = argparser.parse_args().n
    g = argparser.parse_args().g
    sym = argparser.parse_args().sym

    encoding = encode(n, g, forced_sym=sym)
    encoding.serialize(f"formulas/es_{n}_{g}_{sym}sym.cnf")
    print(f"Serialized to: formulas/es_{n}_{g}_{sym}sym.cnf")

if __name__ == "__main__":
    main()
