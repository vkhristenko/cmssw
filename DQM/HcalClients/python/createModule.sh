#!/bin/sh

ref=$1
dest=$2

cp "Hcal$1Client.py" "Hcal$2Client.py"
cp "../interface/Hcal$1Client.h" "../interface/Hcal$2Client.h"
cp "../plugins/Hcal$1Client.cc" "../plugins/Hcal$2Client.cc"

sed -i "s/$1/$2/g" "Hcal$2Client.py" \
	"../interface/Hcal$2Client.h" \
	"../plugins/Hcal$2Client.cc"
sed -i "s/${1^^}/${2^^}/g" "Hcal$2Client.py" \
	"../interface/Hcal$2Client.h" \
	"../plugins/Hcal$2Client.cc"

