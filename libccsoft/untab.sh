#!/bin/bash
sed 's/\t/    /g' ${1} > tmp_untab
cp tmp_untab ${1}
rm tmp_untab
