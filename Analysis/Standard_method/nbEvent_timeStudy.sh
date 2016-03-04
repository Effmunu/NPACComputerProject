for i in `seq 1000 1000 10000`; do
    ./bin/main mc Z 10000 24 0 $i
done

for i in `seq 1000 1000 50000`; do
    ./bin/main data Z 50000 24 0 $i
done
