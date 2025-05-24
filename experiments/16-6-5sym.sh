ALLSAT=$1
REALIZER=$2

if [ -z "$ALLSAT" ] ; then
	echo "Error: Required parameters missing"
	echo "Usage: $0 <allsat_solver> <localizer>"
	exit 1
fi

# encode
python3 encoders/erdos_szekeres.py -n 16 -g 6 -s 5

# solve with allsat
$ALLSAT formulas/es_16_6_5sym.cnf --allsat --datavars=560 > solutions/es_16_6_5sym.sols
# (16 choose 3) = 560.

# decode
if [ -z "$REALIZER" ] ; then
    python3 scripts/decoder.py -f solutions/es_16_6_5sym.sols -n 16 -o orientations --allsat --counts
else
    python3 scripts/decoder.py -f solutions/es_16_6_5sym.sols -n 16 -o orientations --allsat --counts --realizer $REALIZER --sym realizations/helpers/16-6-5fold-sym.txt
fi


