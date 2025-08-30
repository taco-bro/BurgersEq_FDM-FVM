#!/bin/bash

# 设置标志来检查是否已经输出了 "Finish"

FILE="burgersEq.f90"

Finish=false

# rm -f *.dat

let nstart=40000

for ((re=8000; re<=10000; re+=1000)); do

    sed -i "s/integer, parameter[[:space:]]*::[[:space:]]*inputRe[[:space:]]*=[[:space:]]*[0-9]\+/integer, parameter          :: inputRe = $re/" "$FILE"

    for ((N=$nstart; N<=90000; N+=4000)); do

        # 修改Fortran程序中的N值

        sed -i "s/integer, parameter[[:space:]]*::[[:space:]]*N[[:space:]]*=[[:space:]]*[0-9]\+/integer, parameter          :: N = $N/" "$FILE"
        sed -i "s/integer, parameter[[:space:]]*::[[:space:]]*nn[[:space:]]*=[[:space:]]*[0-9]\+/integer, parameter          :: nn = $N/" "$FILE"



        # 验证替换是否成功
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
        echo "Current N value: $N"

        # 编译Fortran程序
        make
        
        # 运行程序
        ./solveBurgersEq
        # 检查输出是否包含 "Finish" 
        if grep -q "This case is converged !!!" <<< "$(./solveBurgersEq)"; then
            Finish=true
            let nstart=$N
            break  # 如果包含 "Finish"，就退出循环
        fi
    done

done
    echo "successful"

