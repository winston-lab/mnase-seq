#!/bin/bash

#SBATCH -p short
#SBATCH -t 8:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH -n 1
#SBATCH -e mogrify.err
#SBATCH -o mogrify.log
#SBATCH -J mogrify-netseq

find . -name "*.svg" ! -path "*.git*" ! -path "*.snakemake*" ! -name "rulegraph.svg" ! -name "dag.svg" | while read svg; do
    png=$(echo $svg | sed -e 's/.svg$/.png/g')
    if [ -e $png ]; then
        #if png exists but svg is newer, make new png
        if [ $svg -nt $png ]; then
            rm $png
            echo "updating png of $svg"
            dim=$(grep -oP "(?<=viewBox=['\"]0 0 ).*?(?=['\"])" $svg)
            width=$(echo $dim | cut -d ' ' -f1 | paste - <(echo \*326\/96) | bc -l)
            height=$(echo $dim | cut -d ' ' -f2 | paste - <(echo \*326\/96) | bc -l)
            svg2png -w $width -h $height -o $png $svg
            # sbatch -p short -t 60 --mem-per-cpu=4G -J "update $png" -o mogrify.log -e mogrify.err --open-mode=append convert.sh $svg $png
        fi
    else
        #if png doesn't exist, make new png
        # sbatch -p short -t 60 --mem-per-cpu=4G -J "create $png" -o mogrify.log -e mogrify.err --open-mode=append convert.sh $svg $png
        echo "creating png of $svg"
        dim=$(grep -oP "(?<=viewBox=['\"]0 0 ).*?(?=['\"])" $svg)
        width=$(echo $dim | cut -d ' ' -f1 | paste - <(echo \*326\/96) | bc -l)
        height=$(echo $dim | cut -d ' ' -f2 | paste - <(echo \*326\/96) | bc -l)
        svg2png -w $width -h $height -o $png $svg
    fi
done
echo "mogrification complete!"
