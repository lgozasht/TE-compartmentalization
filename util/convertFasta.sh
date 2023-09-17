bioawk -c fastx '{print $name"\t"length($seq)}' $1/$2 | awk '$2>50000' > $1/scaffolds.bed 
