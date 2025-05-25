# Automated Symmetric Constructions in Discrete Geometry
Code corresponding to our CICM'2025 submission "Automated Symmetric Constructions in Discrete Geometry".
This is joint work of Bernardo Subercaseaux, Ethan Mackey, Long Qian, and Marijn Heule.

## Requirements

For SAT encodings in Python, we use the `eznf` and `PySAT` libraries. They can be installed by
```
pip install eznf
pip install python-sat
```

For enumerating solutions, we use `allsat-cadical` (https://github.com:jreeves3/allsat-cadical). When we refer to the `allsat` executable, we mean the file `cadical` in the `build` directory of the `allsat-cadical` repository after running `./configure && make` in that directory.

For the realizability problem we use the `Localizer` solver (https://github.com/bsubercaseaux/localizer). We will assume that `localizer` is the executable in your path.


## Symmetries in the Erdős-Szekeres problem

To generate the 66 pointsets with 16 points and no 6-point convex subset that are 4-fold symmetric, run
```
sh experiments/16-6-4sym.sh <path to allsat> <path to localizer>
```
This should take no more than a couple of minutes, and should leave the realizations in the `realizations` folder. The filenames indicate the number of convex 4-gons, convex 5-gons, and convex 6-gons in the pointset, and then the index, respectively.

For the 932 solutions with a 5-fold symmetry, run

```
sh experiments/16-6-5sym.sh <path to allsat> <path to localizer>
```

This should take under 20 minutes.

## Symmetries in the Everywhere-unbalanced-points problem

To generate an s-fold symmetric pointset with n points and an unbalance of k, run

```
python3 encoders/everywhere_unbalanced.py -n <n> -k <k> -s <s>
```

For example, running `python3 encoders/everywhere_unbalanced.py -n 21 -k 2 -s 3` will encode the existence of a set of 15 points with an unbalance of 2 and a 3-fold symmetry, as depicted in the paper. The solver (https://github.com/arminbiere/kissat) can solve this instance in about 10 seconds. In turn, the instance resulting from `python3 encoders/everywhere_unbalanced.py -n 21 -k 2 -s 1`, where no non-trivial symmetries are enforced (any pointset has a 1-fold symmetry), is much harder to solve. 

In terms of the minimality, we use a different encoder for obtaining UNSAT results. Namely,
 `python3 encoders/everywhere_unbalanced_unsat.py -n 11 -k 2 ` will generate an instance which can be proved unsatisfiable in under a second, while `-n 13` takes a few seconds. The instance for `-n 15` takes under a minute, while `-n 17` takes  The instance for `-n 19` took 110 CPU hours with `kissat`.



## Automated proof of Proposition 4.2

The file `scripts/axiom_proof.py` generates the SAT instance described in the proof of Proposition 4.2, which is unsatisfiable. To generate it simply run 

```
python3 scripts/axiom_proof.py
```
after which the file `axiom_proof_5.cnf` will be created. 