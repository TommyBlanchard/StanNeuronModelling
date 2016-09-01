#!/bin/bash
        for i in {1..20}
        do
                for j in {1..3}
		do
			sbatch --wrap="Rscript test_noise_increase.R $i $j" -p ncf --mem=10000 -t 600
		done
        done
