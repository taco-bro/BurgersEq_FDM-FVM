#!/bin/bash

# 设置标志来检查是否已经输出了 "Finish"

FILE="burgersWENOEq.f90"

Finish=false

# rm -f *.dat
# rm -f *.bat

let nstart=100
let relast=800
# cfllast=$(echo "0.005" | bc)
let cfllast=95

for ((j=1; j<=3; j+=1)); do

    let re=2*$relast
    sed -i "s/integer, parameter[[:space:]]*::[[:space:]]*inputRe[[:space:]]*=[[:space:]]*[0-9]\+/integer, parameter          :: inputRe = $re/" "$FILE"

    # Check
    if grep -q "integer, parameter          :: inputRe = $re" "$FILE"; then
        echo "成功替换Re为 $re"
    else
        echo "替换N失败"
        exit 1
    fi
    echo "Current Re value: $re"

    for ((N=$nstart; N<=4000; N+=100)); do

        for ((i=1; i<=10000; i+=1)); do
            # cfl=$(echo "$cfllast + 0.005" | bc)
            let cfl=$cfllast+5
            # Replace the value of N in fortran code
            sed -i "s/integer, parameter[[:space:]]*::[[:space:]]*N[[:space:]]*=[[:space:]]*[0-9]\+/integer, parameter          :: N = $N/" "$FILE"
            sed -i "s/integer, parameter[[:space:]]*::[[:space:]]*nn[[:space:]]*=[[:space:]]*[0-9]\+/integer, parameter          :: nn = $N/" "$FILE"
            # sed -i "s/real*8, parameter[[:space:]]*::[[:space:]]*inputCFL[[:space:]]*=[[:space:]]*[0-9]*\.?[0-9]\+/real*8, parameter           :: inputCFL = $cfl/" "$FILE"
            sed -i "s/integer, parameter[[:space:]]*::[[:space:]]*inputCFL[[:space:]]*=[[:space:]]*[0-9]\+/integer, parameter          :: inputCFL = $cfl/" "$FILE"

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

            if grep -q "integer, parameter          :: inputCFL = $cfl" "$FILE"; then
                echo "成功替换CFL为 $cfl"
            else
                echo "替换CFL失败"
                exit 1
            fi
            echo "Current CFL value: $cfl"

            # 编译Fortran程序
            make
            
            # 运行程序
            ./solveBurgersWENOEq

            # 检查输出是否发散
            if grep -q "NaN" <<< "$(./solveBurgersWENOEq)"; then
                echo "$re, $N, $cfl, unstable" >> unstableCase.bat
                echo "$re, $N, $cfl" >> unstableOutput.dat
                Finish=true
                # cfllast=$(echo "$cfl" | bc)
                let cfllast=95
                break
            fi

            # cfllast=$(echo "$cfl" | bc)
            let cfllast=$cfl
        done

    done
    # let nstart=$N
    let relast=$re

done
    echo "successful"

