import itertools
import argparse
from run_realizer import test_and_output
from convex_tester import count_polygons

# Takes a file with potentially multiple solutions (v-lines) and decodes it.

argparser = argparse.ArgumentParser()
argparser.add_argument(
    "-f", "--file", type=str, required=True, help="Path to the input file"
)
argparser.add_argument("-n", type=int, required=True, help="Number of points")
argparser.add_argument("-o", "--output", type=str, required=True, help="Output folder")
argparser.add_argument("-r", "--realizer", type=str, help="path to realizer")
argparser.add_argument(
    "--counts", action="store_true", help="Count the number of induced convex polygons"
)

argparser.add_argument("--fix", type=str, help="File to fix points for realization")
argparser.add_argument("--sym", type=str, help="File with symmetry orbits")

argparser.add_argument(
    "--allsat", action="store_true", help="Multiple solutions in the input file"
)

file = argparser.parse_args().file
N = argparser.parse_args().n
output_folder = argparser.parse_args().output
allsat = argparser.parse_args().allsat


M = {}
for p, q, r in itertools.combinations(range(1, N + 1), 3):
    M[len(M) + 1] = (p, q, r)

IM = {}
for p, q, r in itertools.combinations(range(N), 3):
    IM[(p, q, r)] = len(IM) + 1

with open(file, "r", encoding="utf-8") as f:
    if allsat:
        cnt = 0
        realizable_cnt = 0
        lines = [line for line in f if line.startswith("v ")]
        printed_realizability_msg = False
        print(f"Found {len(lines)} solutions in {file}")
        # layers = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15]]
        # cross = set()
        # for a, b, c in itertools.combinations(range(N), 3):
        #     if len({(a) // 4, (b) // 4, (c) // 4}) == 3:
        #         cross.add(IM[(a, b, c)])

        for line in lines:

            tokens = line[:-1].split(" ")
            tokens = [int(t) for t in tokens[1:-1]]  # exclude the 0 at the end
            cnt += 1
            content = []
            with open(f"{output_folder}/.tmp.or", "w", encoding="utf-8") as f2:
                for t in tokens:
                    if abs(t) > len(M):
                        continue
                    if t > 0:
                        f2.write(f"A_{M[t]}\n")
                        content.append(f"A_{M[t]}\n")
                    else:
                        f2.write(f"B_{M[-t]}\n")
                        content.append(f"B_{M[-t]}\n")

            output_name = f"N_{N}_sol_{cnt}.or"
            if argparser.parse_args().counts:
                counts = count_polygons(f"{output_folder}/.tmp.or")
                count_content = []
                for c, count in counts.items():
                    count_content.append(f"{count} {c}-gons")
                counts_str = f"{counts[4]}_{counts[5]}_{counts[6]}"
                output_name = f"N_{N}_sol_{counts_str}_{cnt}.or"
                print(
                    f"Solution {cnt} has {sum(counts.values())} convex polygons: {', '.join(count_content)}"
                )

            with open(f"{output_folder}/{output_name}", "w", encoding="utf-8") as f2:
                f2.write("".join(content))

            if argparser.parse_args().realizer is None:
                if not printed_realizability_msg:
                    print("No realizer path provided, skipping realizability test.")
                    printed_realizability_msg = True
                continue
            res = test_and_output(
                f"{output_folder}/{output_name}",
                timeout=2 if N < 20 else 6,
                points_output_file=f"realizations/{output_name}.txt",
                realizer_path=argparser.parse_args().realizer,
                fix_file=argparser.parse_args().fix,
                sym_file=argparser.parse_args().sym,
            )
            if res:
                realizable_cnt += 1
                print(f"Realizable solution {cnt} found!")
            else:
                print(f"Solution {cnt} is (seemingly) not realizable.")
        if argparser.parse_args().realizer is None:
            print(f"Decoded {cnt} solutions! didn't test realizability.")
        else:
            print(f"Decoded {cnt} solutions! {realizable_cnt} were realizable.")
    else:  # TODO: avoid code duplication
        lines = [line for line in f if line.startswith("v ")]
        itkns = []
        for line in lines:
            tokens = line[:-1].split(" ")
            itkns.extend([int(t) for t in tokens[1:]])
        with open(f"{output_folder}/N_{N}_sol.or", "w", encoding="utf-8") as f2:
            for t in itkns:
                if abs(t) > len(M):
                    continue
                if t > 0:
                    f2.write(f"A_{M[t]}\n")
                elif t < 0:
                    f2.write(f"B_{M[-t]}\n")

        if argparser.parse_args().counts:
            counts = count_polygons(f"{output_folder}/N_{N}_sol.or")
            content = []
            for c in counts:
                content.append(f"{counts[c]} {c}-gons")
            print(
                f"Solution has {sum(counts.values())} convex polygons: {', '.join(content)}"
            )

        if argparser.parse_args().realizer is None:
            print("No realizer path provided, skipping realizability test.")
        else:
            res = test_and_output(
                f"{output_folder}/N_{N}_sol.or",
                timeout=2,
                points_output_file=f"realizations/{N}.txt",
                realizer_path=argparser.parse_args().realizer,
            )
            if res:
                print("Realizable solution!")
            else:
                print("Solution is (seemingly) not realizable.")
