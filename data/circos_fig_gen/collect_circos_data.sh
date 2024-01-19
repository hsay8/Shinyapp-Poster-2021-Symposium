#Uses GTDBTK taxonomy prediction, mosdepth, prokka annotation data to produce circos plots.
#Run this AFTER you've done all the previous analysis
mkdir ./circos
#Get final size (after polishing) and mean depth for all contigs in the folder from mosdepth output.
for file in ./mosdepth/polished/*/*.mosdepth.summary.txt; do awk 'FNR==2 {printf $1 "\t"; printf $2 "\t" ; printf "not_real" "\t"; printf "not_real" "\n"}' ${file}; done > cytoband.txt

for file in ./mosdepth/polished/*/*.mosdepth.summary.txt; do awk 'FNR==2 {filename=${file}; printf $1 "\t"; printf $2 "\t" ;printf "not_real" "\t";printf "not_real" "\n"}' ${file}; done > ./circos/cytoband.txt

#Get GC content of final polished assemblies
for file in ./final_assemblies/*.fasta; do Rscript 0-scripts_tools/circlize_gc_information.R -i $file -o ./circos/${file##*/}.txt;done

#Convert prokka output to bed for plotting CDS+ and CDS-, tRNA and rRNA. Run each line separately.
mkdir ./prokka/beds
find ./prokka/* -name '*.gff' -exec cp "{}" ./prokka/beds \;
#remove sequence information
for files in ./prokka/beds/*; do sed -n '/##FASTA*/q;p' $files > $files.txt; done
#remove header lines
for files in ./prokka/beds/*.txt; do awk 'NR>2' $files > $files.txt; done

#Use bakta CDS+ and CDS- instead
mkdir ./bakta/beds
find ./bakta/* -name '*.gff3' -exec cp "{}" ./bakta/beds \;
for files in ./bakta/beds/*; do sed -n '/##FASTA*/q;p' $files > $files.txt; done
for files in ./bakta/beds/*.txt; do awk 'NR>9' $files > ${files}_2.txt; done