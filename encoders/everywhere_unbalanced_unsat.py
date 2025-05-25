from eznf import modeler
import itertools
import argparse

def encode(n_points, d):
    enc = modeler.Modeler()
    points = list(range(n_points))
    point_triples = list(itertools.combinations(points, 3))

    # introduce basic A, B, C variables for orientations
    for p1, p2, p3 in point_triples:
        enc.add_var(f"A_{p1, p2, p3}", f"{p1} is above the line from {p2} to {p3}. Alternatively, p1-p2-p3 is a counterclockwise turn.",
        )
        enc.add_var(f"B_{p1, p2, p3}", f"{p1} is below the line from {p2} and {p3}.")
        enc.add_var(f"C_{p1, p2, p3}", f"{p1} is collinear to the the line through {p2} and {p3}.")
        enc.exactly_one([f"A_{p1, p2, p3}", f"B_{p1, p2, p3}", f"C_{p1, p2, p3}"])


    point_quadruples = list(itertools.combinations(points, 4))
    for p1, p2, p3, p4 in point_quadruples:
        # C_{p1p2p3} ^ C_{p2p3p4} => C_{p1p2p4}
        enc.add_clause([f"-C_{p1, p2, p3}", f"-C_{p2, p3, p4}", f"C_{p1, p2, p4}"])
        # C_{p1 p2 p3} ^ C_{p2 p3 p4} => C_{p1 p3 p4}
        enc.add_clause([f"-C_{p1, p2, p3}", f"-C_{p2, p3, p4}", f"C_{p1, p3, p4}"])
        # C_{p1, p2, p3} ^ C_{p1, p2, p4} => C_{p1, p3, p4}
        enc.add_clause([f"-C_{p1, p2, p3}", f"-C_{p1, p2, p4}", f"C_{p1, p3, p4}"])
        # C_{p1, p2, p3} ^ C_{p1, p3, p4} => C_{p2, p3, p4}
        enc.add_clause([f"-C_{p1, p2, p3}", f"-C_{p1, p3, p4}", f"C_{p2, p3, p4}"])

        # Signotope axioms:
        # (Eq 1, left): A_{p1 p2 p3} v B_{p1 p2 p4} v A_{p1 p3 p4}
        # (Eq 1, right): B_{p1 p2 p3} v A_{p1 p2 p4} v B_{p1 p3 p4}
        enc.add_clause([f"C_{p1, p2, p3}", f"A_{p1, p2, p3}", f"B_{p1, p2, p4}", f"A_{p1, p3, p4}"])
        enc.add_clause([f"C_{p1, p2, p3}", f"B_{p1, p2, p3}", f"A_{p1, p2, p4}", f"B_{p1, p3, p4}"])

        # (Eq 2, left): A_{p1 p2 p3} v B_{p1 p3 p4} v A_{p2 p3 p4}
        # (Eq 2, right): B_{p1 p2 p3} v A_{p1 p3 p4} v B_{p2 p3 p4}
        enc.add_clause([f"C_{p1, p2, p3}", f"A_{p1, p2, p3}", f"B_{p1, p3, p4}", f"A_{p2, p3, p4}"])
        enc.add_clause([f"C_{p1, p2, p3}", f"B_{p1, p2, p3}", f"A_{p1, p3, p4}", f"B_{p2, p3, p4}"])

        # Relations between A, B, and C.
        # A_{p1 p2 p3} ^ C_{p2 p3 p4} => A_{p1 p2 p4}
        # A_{p1 p2 p3} ^ C_{p2 p3 p4} => A_{p1 p3 p4}
        # B_{p1 p2 p3} ^ C_{p2 p3 p4} => B_{p1 p2 p4}
        # B_{p1 p2 p3} ^ C_{p2 p3 p4} => B_{p1 p3 p4}
        enc.add_clause([f"-A_{p1, p2, p3}", f"-C_{p2, p3, p4}", f"A_{p1, p2, p4}"])
        enc.add_clause([f"-A_{p1, p2, p3}", f"-C_{p2, p3, p4}", f"A_{p1, p3, p4}"])
        enc.add_clause([f"-B_{p1, p2, p3}", f"-C_{p2, p3, p4}", f"B_{p1, p2, p4}"])
        enc.add_clause([f"-B_{p1, p2, p3}", f"-C_{p2, p3, p4}", f"B_{p1, p3, p4}"])
        # A_{p1 p2 p4} ^ C_{p1 p2 p3} => A_{p1 p3 p4}
        # A_{p1 p2 p4} ^ C_{p1 p2 p3} => A_{p2 p3 p4}
        # B_{p1 p2 p4} ^ C_{p1 p2 p3} => B_{p1 p3 p4}
        # B_{p1 p2 p4} ^ C_{p1 p2 p3} => B_{p2 p3 p4}
        enc.add_clause([f"-A_{p1, p2, p4}", f"-C_{p1, p2, p3}", f"A_{p1, p3, p4}"])
        enc.add_clause([f"-A_{p1, p2, p4}", f"-C_{p1, p2, p3}", f"A_{p2, p3, p4}"])
        enc.add_clause([f"-B_{p1, p2, p4}", f"-C_{p1, p2, p3}", f"B_{p1, p3, p4}"])
        enc.add_clause([f"-B_{p1, p2, p4}", f"-C_{p1, p2, p3}", f"B_{p2, p3, p4}"])

    # Issue:
    # We're only building variables for p1 < p2 < p3.
    # but then we want to count the number of p3's that are above a pair {p1, p2}.
    # So we want to count A_{p1, p2, p3} for every p3 > p2
    # But also, if there's a p3 between p1 and p2, then we want to count B_{p1, p3, p2}
    # if p3 comes before p1, then consider that
    # A_{p3, p2, p1} = B_{p1, p2, p3}

    # Degeneracy cardinality constraints
    point_pairs = list(itertools.combinations(points, 2))
    for idx, (p1, p2) in enumerate(point_pairs):
        variables_for_above = [enc.v(f"A_{p3, p1, p2}") for p3 in points[:p1]]
        variables_for_above += [enc.v(f"B_{p1, p3, p2}") for p3 in points[p1 + 1 : p2]]
        variables_for_above += [enc.v(f"A_{p1, p2, p3}") for p3 in points[p2 + 1 :]]

        variables_for_below = [enc.v(f"B_{p3, p1, p2}") for p3 in points[:p1]]
        variables_for_below += [enc.v(f"A_{p1, p3, p2}") for p3 in points[p1 + 1 : p2]]
        variables_for_below += [enc.v(f"B_{p1, p2, p3}") for p3 in points[p2 + 1 :]]


        enc.add_svar(f"cA({p1}, {p2})", "COUNTING_VARS", variables=variables_for_above)
        enc.add_svar(f"cB({p1}, {p2})", "COUNTING_VARS", variables=variables_for_below)
        for k in range(0, n_points - 1):
            for i in range(-d+1, d):
                if k + i >= 0 and k + i < n_points - 1:  # and (2*k+i < n_points-1):
                    enc.add_clause(
                        [-enc.v(f"cA({p1}, {p2})_{k}"), -enc.v(f"cB({p1}, {p2})_{i+k}")]
                    )

       

    # symmetry breaking
    for p2 in range(1, n_points):
        for p3 in range(p2 + 1, n_points):
            enc.add_clause([-enc.v(f"B_{0, p2, p3}")])

    return enc

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Encode the everywhere unbalanced points problem (unsat case)")
    parser.add_argument("-n", type=int, help="Number of points")
    parser.add_argument("-k", type=int, default=2, help="imbalance (default: 2)")
    args = parser.parse_args()
    n_points = args.n
    imbalance = args.k
    encoding = encode(n_points, imbalance)
    encoding.serialize(f"everywhere_unbalanced_unsat_{n_points}_{imbalance}.cnf")
