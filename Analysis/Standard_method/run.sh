make;
for i in `seq 2 2 48`; do
    ./bin/main mc Z 10000 $i 0;
done;

