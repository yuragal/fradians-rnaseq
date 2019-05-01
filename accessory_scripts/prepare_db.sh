#!/bin/bash -vex
sp=$(df -h | grep sda1 | awk '{print $4}' | sed 's/G//')
tempd=/store/$USER
db=$1
f=${db##*/}
if [ -e ${tempd}/${f}.ok ]; then
    echo "Nothing to do, database exists."
    exit 0
else
    if [[ $# -eq 2 ]] && [[ $sp -lt $2 ]]; then
        if [[ -e ${tempd} ]]; then
            rm -rf ${tempd}
            mkdir $tempd
        fi
        sp=$(df -h | grep sda1 | awk '{print $4}' | sed 's/G//')
        if [[ $sp -lt $2 ]]; then
            echo `hostname`":there is no space on device to extract the database!!!"
            df -h
            echo "Exiting..."
            exit 1
        fi
    fi

    cd /dev/shm
    while [[ -e ${db}.lock ]]; do
        sleep 30
    done
    touch ${db}.lock
    rm -f ${f}* 2>/dev/null
    cp -f $db ${db}.md5 ./
    if [[ -e ${db}.lock ]]; then
        rm ${db}.lock
    fi
    md5sum -c ${f}.md5

    mkdir -p ${tempd} && cd ${tempd}
    if [[ $f =~ tar\.bz2$ ]]; then
        pbzip2 -dc /dev/shm/${f} | tar x
    else
        pbzip2 -dc /dev/shm/$f > ${f%*.bz2}
    fi
    rm -f /dev/shm/${f}* 2>/dev/null
    touch ${f}.ok
    cd - &>/dev/null
fi
