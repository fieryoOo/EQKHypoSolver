#!/bin/bash

for dir in `ls -d results_*/`; do
	echo $dir
	mkdir -p old_results/${dir}
	ls -d $dir* 2> /dev/null | xargs -I file mv file old_results/${dir}
	rm -r $dir
done
