#!/bin/bash

PREFIX="/run/media/mschu/My_Book_3TB/SPEED"

#GSEs=$(cat index.txt | grep GPL570 | egrep -o "GSE[0-9]+" | sort | uniq)
#GSEs=$(cat index.txt | grep GPL571 | egrep -o "GSE[0-9]+" | sort | uniq)
#GSEs=$(cat index.txt | grep GPL96 | egrep -o "GSE[0-9]+" | sort | uniq)
GSEs=$(cat index.txt | grep GPL6244 | egrep -o "GSE[0-9]+" | sort | uniq)

for GSE in $GSEs; do

    GSEDIR=$(echo $GSE | sed "s/...$/nnn/")
    URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/$GSEDIR/$GSE/suppl"

    if [ ! -d "$PREFIX/$GSE" ]; then
        mkdir -p "$PREFIX/$GSE"
#        echo $GSE >> GPLcurrent.all
    fi

    if [ ! -f "$PREFIX/$GSE/${GSE}_RAW.tar" ]; then
#        wget -O "$PREFIX/$GSE/${GSE}_RAW.tar" $URL/${GSE}_RAW.tar
        wget -O "$PREFIX/$GSE/${GSE}_RAW.tar" "http://www.ncbi.nlm.nih.gov/geo/download/?acc=${GSE}&format=file"
#        tar -xvf "$PREFIX/$GSE/${GSE}_RAW.tar"
        sleep 10
    fi

done

