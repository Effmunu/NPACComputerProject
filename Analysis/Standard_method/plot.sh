for i in `seq 24 2 24`; do
    root -b -l -q "plot_alphas.C(\"data\", \"Z\", \"50000\", $i, 0)"
done;

