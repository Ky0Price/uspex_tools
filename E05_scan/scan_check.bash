#!/bin/bash
rm -f scan_check
for i in EA*
do
name=$i
state=`tail -n 1 $name/log`
#cd $i
echo name=$name    state=$state>>scan_check
#cd $OLDPWD
done
