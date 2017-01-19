#!/usr/bin/env bash

for nb in {1..68} 
	do 
		new=$(printf "batch%04d.sh" ${nb})
		cp template.sh ${new}
		echo "" >> ${new}
		echo "" >> ${new}
			for counter in {1..16}
				do
					let nmod=0
		        	let nmod=$nmod+$counter+$nb*16-1*16 
					a="matlab -nojvm -nodesktop -r \"getINITmodelFromPostProbs($nmod,0.99,0.05)\" > INITOutputLog$nmod 2> INITErrorLog$nmod &"
					echo $a >> ${new}
				done
				echo "wait" >> ${new}
		sbatch ${new}
	done