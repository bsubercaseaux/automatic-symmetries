from eznf import modeler
import itertools


def encode(N):
    enc = modeler.Modeler()
    for (p, q, r) in itertools.permutations(range(N), 3):
        enc.add_var(f"cc_{p, q, r}")
        
    # basic axioms
    for (p, q, r) in itertools.permutations(range(N), 3):
        enc.add_clause([f"-cc_{p, q, r}", f"cc_{q, r, p}"])
        enc.add_clause([f"cc_{p, q, r}", f"cc_{p, r, q}"])
        enc.add_clause([f"-cc_{p, q, r}", f"-cc_{p, r, q}"])
                    
    for (p, q) in itertools.permutations(range(N), 2):
        enc.add_var(f"<_{p, q}")
        
    for (p, q) in itertools.permutations(range(N), 2):
        enc.add_clause([f"<_{p, q}", f"<_{q, p}"])
        enc.add_clause([f"-<_{p, q}", f"-<_{q, p}"])
        
    # transitivity <
    for (p, q, r) in itertools.permutations(range(N), 3):
        enc.add_clause([f"-<_{p, q}", f"-<_{q, r}", f"<_{p, r}"])
        enc.add_clause([f"<_{p, q}", f"<_{q, r}", f"-<_{p, r}"])
    

    def cc(p, q, r):
        return enc.v(f"cc_{p, q, r}")
    
    # ordered signotope axioms
    for (p, q, r, s) in itertools.permutations(range(N), 4):
        enc.add_clause([f"-<_{p, q}", f"-<_{p, r}", f"-<_{p, s}",  cc(p,q,r), -cc(p,q,s),  cc(p,r,s)])
        enc.add_clause([f"-<_{p, r}", f"-<_{q, r}", f"-<_{r, s}",  cc(p,q,r), -cc(p,r,s),  cc(q,r,s)])


    # the clauses of the transitivity axiom that will be negated using Tseitin
    trans_clauses = []
    for (p, q, r, t, s) in itertools.permutations(range(N), 5):
        trans_clauses.append([-cc(t, s, p), -cc(t, s, q), -cc(t, s, r), -cc(t, p, q), -cc(t, q, r), cc(t, p, r)])

        
    for cls in trans_clauses:
        enc.add_var(f"y_{cls}", description=f"repreents that {cls} is falsified")
        for lit in cls:
            enc.add_clause([f"-y_{cls}", -lit])
        
    enc.add_clause([f"y_{cls}" for cls in trans_clauses])

    return enc
    
encoding = encode(5)
encoding.serialize("axiom_proof_5.cnf")    