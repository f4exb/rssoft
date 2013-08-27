#!/bin/sh
for pi in $(ps -edf | grep ${1} | cut -c10-15); 
do 
    kill -9 $pi; 
done
