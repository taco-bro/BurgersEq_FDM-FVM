#!/bin/bash

# 设置标志来检查是否已经输出了 "Finish"

FILE="burgersEq.f90"

Finish=false

rm -f *.dat

let nstart=200
let relast=50

for ((j=1; j<=3; j+=1)); do

    let re=2*$relast
    let nstart=2*$nstart
    sed -i "s/integer, parameter[[:space:]]*::[[:space:]]*inputRe[[:space:]]*=[[:space:]]*[0-9]\+/integer, parameter          :: inputRe = $re/" "$FILE"

    for ((N=$nstart; N<=$nstart; N+=10)); do

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
            # let nstart=$N
            # break  # 如果包含 "Finish"，就退出循环
            echo '为收敛'
            exit 1
        fi
        let nstart=$N
        let relast=$re
    done

done
    echo "successful"

