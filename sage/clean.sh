#!/bin/sh
find . ! -regex '^.*sage$' ! -regex '^.*sh$' -exec rm -f {} \; 
