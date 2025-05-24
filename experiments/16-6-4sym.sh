ALLSAT=$1
REALIZER=$2

if [ -z "$ALLSAT" ] ; then
	echo "Error: Required parameters missing"
	echo "Usage: $0 <allsat_solver> <localizer>"
	exit 1
fi

# encode
python3 encoders/erdos_szekeres.py -n 16 -g 6 -s 4

# solve with allsat
$ALLSAT formulas/es_16_6_4sym.cnf --allsat --datavars=560 > solutions/es_16_6_4sym.sols
# (16 choose 3) = 560.

# decode
if [ -z "$REALIZER" ] ; then
    python3 scripts/decoder.py -f solutions/es_16_6_4sym.sols -n 16 -o orientations --allsat --counts
else
    python3 scripts/decoder.py -f solutions/es_16_6_4sym.sols -n 16 -o orientations --allsat --counts --realizer $REALIZER --fix realizations/helpers/16-6-4fold-fix.txt --sym realizations/helpers/16-6-4fold-sym.txt
fi


