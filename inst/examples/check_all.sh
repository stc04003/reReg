#! /bin/bash

rm .RData
rm .Rhistory
rm *.Rout
rm *.pdf

for f in ex_*
do 
    echo "library(reReg);source('$f')" > "$f"2.R
    R -d "valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all -s --log-file=$f"out --vanilla < "$f"2.R
    rm "$f"2.R
done
