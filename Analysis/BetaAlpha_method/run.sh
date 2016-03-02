make
# If it compiles: run, else error message
if [ $? -ne 0 ]; then
    echo "error while compiling"
    exit 1
else
    for i in `seq 2 2 24`; do
        ./bin/main mc Z 10000 $i 0
    done;
    exit 0
fi

