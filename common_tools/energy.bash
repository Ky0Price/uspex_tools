#!/bin/bash
for i in Li*;do echo -e $i "," $(grep 'without' $i/OUTCAR | tail -n 1 | awk '{print $7}');done > Energy.csv
# 这里的 Li* 的意思是所有以 Li 开头的文件夹，如果你的文件夹不是以 Li 开头的，那么你需要修改这个参数。
# 这个脚本可以用于获取常压下的能量，与高压下的脚本有所不同，因为高压下的能量是用 enthalpy 来表示的，enthalpy = energy + pressure * volume