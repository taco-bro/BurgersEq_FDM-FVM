#!/bin/bash

# 设置标志来检查是否已经输出了 "Finish"

FILE="burgersEq.f90"

Finish=false

rm -f *.dat

let nstart=1600
let relast=800

for ((j=1; j<=3; j+=1)); do

    let re=2*$relast

    sed -i "s/integer, parameter[[:space:]]*::[[:space:]]*inputRe[[:space:]]*=[[:space:]]*[0-9]\+/integer, parameter          :: inputRe = $re/" "$FILE"

    # Check
    if grep -q "integer, parameter          :: inputRe = $re" "$FILE"; then
        echo "成功替换Re为 $re"
    else
        echo "替换Re失败"
        exit 1
    fi
    echo "Current Re value: $re"

    let N=2*$nstart

    # Replace the value of N in fortran code
    sed -i "s/integer, parameter[[:space:]]*::[[:space:]]*N[[:space:]]*=[[:space:]]*[0-9]\+/integer, parameter          :: N = $N/" "$FILE"
    sed -i "s/integer, parameter[[:space:]]*::[[:space:]]*nn[[:space:]]*=[[:space:]]*[0-9]\+/integer, parameter          :: nn = $N/" "$FILE"

    # Check
    if grep -q "integer, parameter          :: N = $N" "$FILE"; then
        echo "成功替换N为 $N"
    else
        echo "替换N失败"
        exit 1
    fi

    if grep -q "integer, parameter          :: nn = $N" "$FILE"; then
        echo "成功替换nn为 $N"
    else
        echo "替换nn失败"
        exit 1
    fi

    # 编译Fortran程序
    make
    
    # 运行程序
    ./solveBurgersEq

    let nstart=$N
    let relast=$re

done
    echo "successful"

