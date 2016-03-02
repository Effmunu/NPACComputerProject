make
# If it compiles: run, else error message
if [ $? -ne 0 ]; then
    echo "error while compiling"
    exit 1
else
    for i in `seq 2 2 48`; do
        root -b -l -q "plot_alphas.C(\"mc\", \"Z\", \"10000\", $i, 0)"
    done;
    exit 0
fi

