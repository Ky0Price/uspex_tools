#!/bin/bash
for i in Li*;do echo -e $i "," $(grep 'enthalpy' $i/OUTCAR | tail -n 1 | awk '{print $5}');done > Enthalpy.csv
# 这里的 Li* 的意思是所有以 Li 开头的文件夹，如果你的文件夹不是以 Li 开头的，那么你需要修改这个参数。
# 高压下可以用这个脚本来批量获取能量，与常压下的脚本有所不同，因为高压下的能量是用 enthalpy 来表示的，enthalpy = energy + pressure * volume