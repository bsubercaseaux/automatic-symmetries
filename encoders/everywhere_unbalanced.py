from itertools import combinations

from pysat.card import *
import itertools


class CC_system:
    rot = 1
    encode_id = 1
    orients = {}
    n = 0
    point_to_id = {}
    id_to_point = {}
    label_to_class = {}
    class_to_label = {}
    equiv_classes = []

    def __init__(self, rot=1, id=1):
        self.rot = rot
        self.encode_id = id

    def add_point(self, label):
        ##Populate new points
        new_points = [self.n + i for i in range(self.rot)]
        self.label_to_class[label] = len(self.equiv_classes)
        self.class_to_label[len(self.equiv_classes)] = label
        for i in range(self.rot):
            self.point_to_id[(self.label_to_class[label], i)] = new_points[i]
            self.id_to_point[new_points[i]] = (self.label_to_class[label], i)

        ##Populate orients
        for idx in range(self.rot):
            p3 = (self.label_to_class[label], idx)
            for i1, i2 in itertools.combinations(list(range(self.n)) + new_points[:idx], 2):
                p1 = self.id_to_point[i1]
                p2 = self.id_to_point[i2]
                triple, inv = self.normalize((p1, p2, p3))
                # print(f"Triple: {(p1, p2, p3)} Normalized: {triple} inv: {inv}")
                if not triple in self.orients:
                    # print("Unpopulated triple")
                    self.orients[triple] = [self.encode_id, self.encode_id + 1, self.encode_id + 2]
                    self.encode_id += 3
                self.orients[(p1, p2, p3)] = [self.orients[triple][2 * inv], self.orients[triple][1],
                                              self.orients[triple][2 * (1 - inv)]]
        ##Update parameters
        self.equiv_classes.append([(self.n, i) for i in range(self.rot)])
        self.n += self.rot

    def normalize(self, triple):
        ##Returns a normalized triple, along with a flag indicating if orientation has been flipped
        ##Points are (class, rot) pairs
        inv = len([(u, v) for u, v in itertools.combinations(triple, 2) if u[0] > v[0]])
        p1, p2, p3 = sorted(triple, key=lambda x: x[0])
        np1 = (p1[0], 0)
        np2 = (p2[0], (p2[1] - p1[1]) % self.rot)
        np3 = (p3[0], (p3[1] - p1[1]) % self.rot)
        if np2 < np3:
            return (np1, np2, np3), inv % 2
        return (np1, np3, np2), (inv + 1) % 2

    def get_orients(self, tri):
        i1, i2, i3 = tri
        p1 = self.id_to_point[i1]
        p2 = self.id_to_point[i2]
        p3 = self.id_to_point[i3]
        triple, inv = self.normalize((p1, p2, p3))
        if inv == 0:
            return self.orients[triple]
        else:
            x, y, z = self.orients[triple]
            return [z, y, x]

    def get_orient_lit(self, tri, ori):
        return self.get_orients(tri)[ori]

    def encode_uniquess(self, cnf, id):
        ##Encodes that every triple has exactly one orientation
        for triple in self.orients:
            if triple == self.normalize(triple)[0]:
                id = append_eq(cnf, self.orients[triple], 1, id)
        return id

    def is_canonical(self, triple):
        i1, i2, i3 = triple
        p1 = self.id_to_point[i1]
        p2 = self.id_to_point[i2]
        p3 = self.id_to_point[i3]
        return (p1, p2, p3) == self.normalize((p1, p2, p3))[0]

def append_eq(cnf, lits, k, id):
    ##Encodes sum x \in lits x = k
    ##Returns the next available identifier.
    if k == 1:
        eq = CardEnc.equals(lits, k, id - 1, encoding=EncType.pairwise)
    else:
        eq = CardEnc.equals(lits, k, id - 1, encoding=EncType.kmtotalizer)
    for e in eq:
        for l in e:
            id = max(id, abs(l))
        cnf.append(e)
    return id + 1


def append_am(cnf, lits, k, id):
    ##Encodes sum x \in lits x <= k
    ##Returns the next available identifier.
    if k == 1:
        eq = CardEnc.atmost(lits, k, id - 1, encoding=EncType.pairwise)
    else:
        eq = CardEnc.atmost(lits, k, id - 1, encoding=EncType.kmtotalizer)
    for e in eq:
        for l in e:
            id = max(id, abs(l))
        cnf.append(e)
    return id + 1


def cpt_amo(lits, cnf, id):
    ##Compact encoding of x_1 + x_2 + ... + x_n <= 1
    ##cnf: CNF object
    ##id: first available identifier
    if len(lits) <= 4:
        id = append_am(cnf, lits, 1, id)
        return id
    else:
        pref = lits[0:3] + [-id]
        suff = [id] + lits[3:]
        id = cpt_amo(pref, cnf, id + 1)
        id = cpt_amo(suff, cnf, id)
        return id


def sort_triple(tri):
    ##returns sorted triple and parity of inversions
    return sorted(tri), len([(u, v) for (u, v) in itertools.combinations(tri, 2) if u > v]) % 2


def get_orient_lit(tri, ori, orients):
    ##Given a triple tri, returns the literal i such that i <-> (a, b, c) has orientation ori
    cur, inv = sort_triple(tri)
    if ori == 1 or inv == 0:
        return orients[tuple(cur)][ori]
    else:
        return orients[tuple(cur)][(ori + 2) % 4]


def cyclic_permutations(it):
    ##Generates all permutations up to cyclic rotations
    pivot = it[0]
    syms = []
    for perm in itertools.permutations(it):
        if perm[0] == pivot:
            syms.append(perm)
    return syms


def append_diff(lits1, lits2, k, cnf, id):
    ##Encodes sum lits1 ~ k + sum lits2 for ~ \in {>=, <=}
    ##requires len(lits1) == len(lits2)
    ##id: first available id
    ##returns; first available id
    s = len(lits1)
    assert len(lits1) == len(lits2), "Invalid discrepancy encoding"

    xcomp = [[j + (i + 2) * (i + 1) // 2 + id - 1 for j in range(i + 2)] for i in range(s)]
    id += (s + 2) * (s + 1) // 2 - 1
    ycomp = [[j + (i + 2) * (i + 1) // 2 + id - 1 for j in range(i + 2)] for i in range(s)]
    id += (s + 2) * (s + 1) // 2 - 1

    lits = [lits1, lits2]
    comps = [xcomp, ycomp]
    # print(comps[0])
    # print(comps[1])
    for idx in range(2):
        lit = lits[idx]
        comp = comps[idx]

        ##Base case
        cnf.append([-lit[0], comp[0][1]])
        cnf.append([-lit[0], -comp[0][0]])
        cnf.append([lit[0], -comp[0][1]])
        cnf.append([lit[0], comp[0][0]])

        for i in range(1, s):
            for p in range(i + 2):
                if p == 0:
                    cnf.append([-comp[i][p], -lit[i]])
                    cnf.append([-comp[i][p], comp[i - 1][p]])
                    cnf.append([comp[i][p], -comp[i - 1][p], lit[i]])
                else:
                    cnf.append([-comp[i][p], -lit[i], comp[i - 1][p - 1]])
                    cnf.append([comp[i][p], -lit[i], -comp[i - 1][p - 1]])
                    if p <= i:
                        cnf.append([comp[i][p], lit[i], -comp[i - 1][p]])
                        cnf.append([-comp[i][p], lit[i], comp[i - 1][p]])
                    else:
                        cnf.append([lit[i], -comp[i][p]])
    for i in range(s + 1):
        for j in range(s + 1):
            if abs(i - j) < k:
                cnf.append([-comps[0][-1][i], -comps[1][-1][j]])

    return id


def encode_balance(a, b, point_config: CC_system, k, cnf, id):
    ##Encodes that the line joining points (a, b) is k-unbalanced
    pos = []
    neg = []
    for p in range(point_config.n):
        if p in [a, b]:
            continue
        pos.append(point_config.get_orient_lit((a, b, p), 2))
        neg.append(point_config.get_orient_lit((a, b, p), 0))

    ##Want to encode pos >= k + neg or neg >= k + pos, surely there's a better way?

    # print(f"Encoding discrepancy at: {a, b}")
    # print(f"Positive ids: {pos}")
    # print(f"Negative ids: {neg}")

    id = append_diff(pos, neg, k, cnf, id)

    return id


def encode_rot_sym(points, orients, rot, cnf):
    ##Encodes that the point configuration has a rotational symmetry of cycle = rot
    ##Currently only asserts equivalence of equivalence classes.
    n = len(points)
    assert n % rot == 0, "Total number of points is not divisible by cycle of rotation."
    equiv_size = n // rot
    base = [i for i in range(equiv_size)]
    print(base)
    for triple in itertools.combinations(base, 3):
        a, b, c = triple
        for i in range(1, 3):
            offset = equiv_size * i
            for k in range(3):
                cnf.append([-orients[triple][k], orients[((a + offset) % n, (b + offset) % n, (c + offset) % n)][k]])
                cnf.append([orients[triple][k], -orients[((a + offset) % n, (b + offset) % n, (c + offset) % n)][k]])


def get_signotope_axiom(lits, point_config: CC_system):
    ##Returns a list of signotope axioms
    assert len(lits) == 4
    clauses = []
    i, j, k, l = lits
    s1 = (i, j, k)  ##abc
    s2 = (i, k, l)  ##acd
    s3 = (i, j, l)  ##abd
    s4 = (j, k, l)  ##bcd
    configs = [[s1, s2, s3], [s1, s4, s2]]
    for config in configs:
        t1, t2, t3 = config  ##Implement t1 ~ 0 \land t2 ~ 0 \implies t3 ~ 0 where ~ = {>, <}
        for sign in [0, 2]:
            l1_pos = point_config.get_orient_lit(t1, sign)
            l2_pos = point_config.get_orient_lit(t2, sign)
            l3_pos = point_config.get_orient_lit(t3, sign)
            clauses.append([-l1_pos, -l2_pos, l3_pos])
    return clauses


def encode_n4_axioms(n, orients, cnf):
    ##Uses O(n^4) clauses, assumes points are in increasing x-coordinate
    ##Forbidding non-realizable patterns: Section 3.3 in https://arxiv.org/pdf/2403.00737
    for i, j, k, l in itertools.combinations(range(n), 4):
        s1 = (i, j, k)  ##abc
        s2 = (i, k, l)  ##acd
        s3 = (i, j, l)  ##abd
        s4 = (j, k, l)  ##bcd
        configs = [[s1, s2, s3], [s1, s4, s2]]
        for config in configs:
            t1, t2, t3 = config  ##Implement t1 ~ 0 \land t2 ~ 0 \implies t3 ~ 0 where ~ = {>, <, =}
            for sign in [0, 2]:
                cnf.append([-orients[t1][sign], -orients[t2][sign], orients[t3][sign]])
                cnf.append([-orients[t1][sign], -orients[t2][1], orients[t3][sign]])
                cnf.append([-orients[t1][1], -orients[t2][sign], orients[t3][sign]])
            cnf.append([-orients[t1][1], -orients[t2][1], orients[t3][1]])


def encode_matroid_axioms(n, orients, cnf):
    ##Uses O(n^5) clauses, does not assume position of points
    ##Encodes general axioms for an oriented matroid

    ##3-term relations
    for b1, b2, a1, a2, a3 in itertools.permutations(range(n), 5):
        t1 = (b1, a2, a3)
        t2 = (a1, b2, a3)
        t3 = (b2, a2, a3)
        t4 = (b1, a1, a3)
        t5 = (a1, a2, a3)
        t6 = (b1, b2, a3)
        for s1 in [-1, 0, 1]:
            for s2 in [-1, 0, 1]:
                if s1 * s2 < 0:
                    continue
                for s3 in [-1, 0, 1]:
                    for s4 in [-1, 0, 1]:
                        if s3 * s4 < 0:
                            continue
                        for s5 in [-1, 0, 1]:
                            if s5 == 0:
                                continue
                            l1 = get_orient_lit(t1, s1 + 1, orients)
                            l2 = get_orient_lit(t2, s2 + 1, orients)
                            l3 = get_orient_lit(t3, s3 + 1, orients)
                            l4 = get_orient_lit(t4, s4 + 1, orients)
                            l5 = get_orient_lit(t5, s5 + 1, orients)
                            l6 = get_orient_lit(t6, (s5 + 3) % 4, orients)
                            cnf.append([-l1, -l2, -l3, -l4,
                                        -l5, -l6])


def order_collinear_axiom(lits, point_config: CC_system, order, cnf):
    assert len(lits) == 4
    a, b, c, d = lits
    assert a < b < c, "Collinearity only handled for monotone triples"

    col = point_config.get_orient_lit((a, b, c), 1)
    pairs = list(itertools.combinations((a, b, c), 2))
    p1 = pairs[0]

    ## Transitivity of collinearity
    for i in range(1, 3):
        pair = pairs[i]
        s1 = point_config.get_orient_lit((p1[0], p1[1], d), 1)
        s2 = point_config.get_orient_lit((pair[0], pair[1], d), 1)
        cnf.append([-col, -s1, s2])
        cnf.append([-col, s1, -s2])

    for k in [-1, 1]:
        ##Note that we may assume WLOG that none of these are collinear by clauses above
        abd = point_config.get_orient_lit((a, b, d), 2)
        acd = point_config.get_orient_lit((a, c, d), 2)
        bcd = point_config.get_orient_lit((b, c, d), 2)
        abd_col = point_config.get_orient_lit((a, b, d), 1)

        ##Case 1: a -> b -> c, in this case all orientations should be equal
        case_1 = [-col, k * order[(a, b)], k * order[(b, c)]]
        cnf.append(case_1 + [-abd, acd])
        cnf.append(case_1 + [abd, -acd])
        cnf.append(case_1 + [-abd, bcd])
        cnf.append(case_1 + [abd, -bcd])

        ##Case 2: a -> c <- b, in this case abd = acd = -bcd
        case_2 = [-col, abd_col, k * order[(a, b)], k * order[(a, c)], -k * order[(b, c)]]
        cnf.append(case_2 + [-abd, acd])
        cnf.append(case_2 + [abd, -acd])
        cnf.append(case_2 + [-abd, -bcd])
        cnf.append(case_2 + [abd, bcd])
        #
        ##Case 3: c <- a -> b, in this case abd = -acd = -bcd
        case_3 = [-col, abd_col, k * order[(a, b)], -k * order[(a, c)]]
        cnf.append(case_3 + [-abd, -acd])
        cnf.append(case_3 + [abd, acd])
        cnf.append(case_3 + [-abd, -bcd])
        cnf.append(case_3 + [abd, bcd])


def collinear_axiom(lits, orients, cohere, cnf):
    assert len(lits) == 4
    a, b, c, d = lits

    ##Encodes that if (a, b, c) are collinear, then their relationship to d is the same
    ##This really also asserts that if (a, b, c) is collinear, then a < b < c is necessarily monotone
    ##It is important that this is preserved under rotation

    assert a < b < c, "Collinearity only handled for monotone triples"
    col = orients[(a, b, c)][1]
    pairs = list(itertools.combinations((a, b, c), 2))
    p1 = pairs[0]
    ## Transitivity of collinearity
    for i in range(1, 3):
        pair = pairs[i]
        s1 = get_orient_lit((p1[0], p1[1], d), 1, orients)
        s2 = get_orient_lit((pair[0], pair[1], d), 1, orients)
        cnf.append([-col, -s1, s2])
        cnf.append([-col, s1, -s2])

    ## Enforcing coherence
    for i in range(1, 3):
        pair = pairs[i]
        mult = cohere[(a, b, c)][i - 1]
        for dir in [0, 2]:
            s1 = get_orient_lit((p1[0], p1[1], d), dir, orients)
            s2_pos = get_orient_lit((pair[0], pair[1], d), dir, orients)
            s2_neg = get_orient_lit((pair[0], pair[1], d), (dir + 2) % 4, orients)
            cnf.append([-col, -s1, -mult, s2_pos])
            cnf.append([-col, -s1, mult, s2_neg])

    ##Necessary conditions
    for dir in [0, 2]:
        abd = get_orient_lit((a, b, d), dir, orients)

        ##Case works, either c is to the left of a, between a and b, or to the right of b.

        ##Cae 1: To the right of b
        bcd = get_orient_lit((b, c, d), dir, orients)
        acd = get_orient_lit((a, c, d), dir, orients)
        cnf.append([-col, -abd, -bcd, acd])

        ##Case 2: To the left of a
        cad = get_orient_lit((c, a, d), dir, orients)
        cbd = get_orient_lit((c, b, d), dir, orients)
        cnf.append([-col, -abd, -cad, cbd])

        ##Case 3: Between a and b
        cnf.append([-col, -abd, -cad, -bcd])


def interior_axiom(lits, orients, cnf):
    assert len(lits) == 4
    a, b, c, d = lits
    s1 = get_orient_lit((b, c, d), 2, orients)
    s2 = get_orient_lit((c, a, d), 2, orients)
    s3 = get_orient_lit((a, b, d), 2, orients)
    s4 = get_orient_lit((a, b, c), 2, orients)
    cnf.append([-s1, -s2, -s3, s4])


def transitive_axiom(lits, orients, cnf):
    assert len(lits) == 5
    p, q, r, s, t = lits
    for dir in [2]:
        s1 = get_orient_lit((t, s, p), dir, orients)
        s2 = get_orient_lit((t, s, q), dir, orients)
        s3 = get_orient_lit((t, s, r), dir, orients)
        s4 = get_orient_lit((t, p, q), dir, orients)
        s5 = get_orient_lit((t, q, r), dir, orients)
        s6 = get_orient_lit((t, p, r), dir, orients)
    cnf.append([-s1, -s2, -s3, -s4, -s5, s6])


def not_inside(lits, point_config: CC_system, cnf):
    ##For lits = (a, b, c, d), encodes that a is not in the triangle formed by b, c, d
    a, b, c, d = lits
    col = point_config.get_orient_lit((a, b, c), 1)
    for dir in [0, 2]:
        s1 = point_config.get_orient_lit((a, b, c), dir)
        s2 = point_config.get_orient_lit((a, c, d), dir)
        s3 = point_config.get_orient_lit((a, b, d), dir)
        cnf.append([col, -s1, -s2, s3])

def get_order_lit(p, q, ordering):
    ##Returns literal encoding p \prec q
    assert p != q, "cannot order equivalent points"
    return ordering[(p, q)] if p < q else -ordering[(q, p)]


def encode_signotope_axioms(point_config: CC_system, cnf, id):
    ## Encodes signotope axioms on unsorted points by encoding
    ## an ordering on the points that respects the x-coordinates
    ## should only use O(n^4) clauses

    n = point_config.n
    ordering = {}
    for pair in itertools.combinations(range(n), 2):
        ordering[pair] = id
        id += 1

    ## Transitivity of ordering
    for lits in itertools.permutations(range(n), 3):
        monotone = []
        for i in range(2):
            monotone.append(
                -ordering[(lits[i], lits[i + 1])] if lits[i] < lits[i + 1] else ordering[(lits[i + 1], lits[i])])
        cnf.append(monotone + [ordering[(lits[1], lits[2])] if lits[1] < lits[2] else -ordering[(lits[2], lits[1])]])

    ##Collinearity: Enforced for every set of 4 points (a, b, c, d) where a < b < c
    for triple in itertools.combinations(range(n), 3):
        # if not point_config.is_canonical(triple):
        #     continue
        for p in range(n):
            if p in triple:
                continue
            order_collinear_axiom((triple[0], triple[1], triple[2], p), point_config, ordering, cnf)

    ## Standard signotope axioms
    for p, q, r, s in itertools.permutations(range(n), 4):
        signotope_axioms = get_signotope_axiom((p, q, r, s), point_config)
        pq = get_order_lit(p, q, ordering)
        qr = get_order_lit(q, r, ordering)
        rs = get_order_lit(r, s, ordering)
        mont = [-pq, -qr, -rs]
        for c in signotope_axioms:
            cnf.append(mont + c)

    ##Enforce convex hull ordering
    for lits in itertools.combinations(range(n), 4):
        not_inside(lits, point_config, cnf)
    return id


def encode_cc_axioms(n, orients, cnf, id):
    ##Axioms 1-3 are implicitly encoded

    ##Collinearity: auxiliary variables encoding coherence of collinear triples
    cohere = {}
    for triple in itertools.combinations(range(n), 3):
        cohere[triple] = [id, id + 1]
        ##Cannot be (-1, +1)
        cnf.append([-(id + 1), id])

        id += 2
    ##Collinearity: Enforced for every set of 4 points (a, b, c, d) where a < b < c
    for triple in itertools.combinations(range(n), 3):
        for p in range(n):
            if p in triple:
                continue
            collinear_axiom((triple[0], triple[1], triple[2], p), orients, cohere, cnf)

    ##Strict interior axiom
    for cur in itertools.combinations(range(n), 4):
        a, b, c, d = cur
        interior_axiom((a, b, c, d), orients, cnf)
        interior_axiom((a, b, d, c), orients, cnf)

    ##Strict transitivity axiom
    for pair in itertools.permutations(range(n), 2):
        s, t = pair
        for triple in itertools.combinations(range(n), 3):
            if len(set(triple) & set(pair)) > 0:
                continue
            p, q, r = triple
            transitive_axiom((p, q, r, s, t), orients, cnf)
            transitive_axiom((p, r, q, s, t), orients, cnf)
    return id


def encode_rot(n, orients, cycle, cnf):
    assert n % cycle == 0, "Number of points not divisible by cycle of rotation."
    class_size = n // cycle
    for triple in itertools.combinations(range(n), 3):
        a, b, c = triple
        for i in range(1, cycle):
            offset = i * class_size
            an = (a + offset) % n
            bn = (b + offset) % n
            cn = (c + offset) % n
            for s in range(3):
                lit = get_orient_lit((an, bn, cn), s, orients)
                cnf.append([-orients[triple][s], lit])
                cnf.append([orients[triple][s], -lit])


def encode(n, balance, rot=1):
    assert n % rot == 0, "Incompatible rotation"
    cnf = CNF()
    point_config = CC_system(rot, 1)
    for i in range(n // rot):
        point_config.add_point(i)
    print(f"Unique orientation variables: {point_config.encode_id - 1} Total points: {point_config.n}")
    # print(point_config.orients)

    ##Fixes
    cnf.append([point_config.get_orient_lit((0, 1, 2), 2)])

    id = point_config.encode_id
    id = point_config.encode_uniquess(cnf, id)
    id = encode_signotope_axioms(point_config, cnf, id)

    for a, b in itertools.combinations(range(n), 2):
        if point_config.id_to_point[a][1] == 0:
            id = encode_balance(a, b, point_config, balance, cnf, id)
    cnf.to_file(f"formulas/everywhere_unbalanced_{n}_{balance}_{rot}.cnf")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Encoder')
    parser.add_argument('-n', '--npoints', type=int, help='Number of points', required=True)
    parser.add_argument('-k', '--unbalanced_size', type=int, help='Size of inequality', required=True)
    parser.add_argument('-s', '--symmetry', type=int, help='Cycle of rotation symmetry', default=1)

    args = parser.parse_args()
    n = args.npoints
    k = args.unbalanced_size
    s = args.symmetry
    assert n >= 3, "Number of points smaller than 3"
    encode(n, k, s)
