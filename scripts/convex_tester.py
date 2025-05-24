import itertools
import argparse
import functools
import re


def parse_orientations(orientations_file):
    # 1-pass: parse each line, record triples + orientation and track max index
    triples = []
    max_i = 0
    pat = re.compile(r'([AB])_\((\d+), (\d+), (\d+)\)')
    with open(orientations_file, 'r') as f:
        for L in f:
            m = pat.match(L.strip())
            if not m:
                print(f"Invalid line: {L.strip()}")
                continue
            label, ps, qs, rs = m.groups()
            p, q, r = int(ps), int(qs), int(rs)
            triples.append((p, q, r, label == 'A'))
            max_i = max(max_i, p, q, r)

    # now build 3D matrix directly
    n = max_i
    cc = [[[False]*(n+1) for _ in range(n+1)] for _ in range(n+1)]
    for p, q, r, val in triples:
        cc[p][q][r] = val
    return cc, n
    
    
def count_polygons(orientations_file):
    cc, n = parse_orientations(orientations_file)
    
    @functools.lru_cache(maxsize=None)
    def convex(p, q, r, s):
        return (cc[p][q][r] == cc[p][r][s]) == (cc[p][q][s] == cc[q][r][s])
        
    convex4 = [quad for quad in itertools.combinations(range(1,n+1),4)
           if convex(*quad)]
           
    conv_polys = {4: convex4}
    counts = {4: len(convex4), 5: 0, 6: 0, 7:0}
    for k in (5,6,7):
        next_polys = []
        for poly in conv_polys[k-1]:
            last = poly[-1]
            for p in range(last+1, n+1):
                # only need to check each triple from poly plus p
                if all(convex(a,b,c,p) for a,b,c in itertools.combinations(poly,3)):
                    next_polys.append(poly + (p,))
        conv_polys[k] = next_polys
        counts[k] = len(next_polys)
        
    return counts
        # print(f"Number of convex {k}-gons: {cnt}")
                
                
if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Count convex polygons in an orientation file")
    argparser.add_argument("-f", "--file", type=str, required=True, help="Path to the input orientation file")
    file = argparser.parse_args().file
    print(count_polygons(file))