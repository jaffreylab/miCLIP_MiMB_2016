# This script uses bedtools to determine the 5-mer sequences of potential m6A sites.
# It also generates a summary of the DRACH sites detected. 
# The script requires two input arguments: 
# 1.) the bed file of truncations to work on 
# 2.) the fasta file of the genome used for alignment. This fasta file needs to be indexed with a "genome.fasta.fai" index file in the same directory. 
# Example usage: kmer.annotate.cits.sh my.truncations.bed genomes/genome.fasta

# determine file prefix for temp files
s=$1
pf=$(echo ${s%.*})
#as the truncations directly call the m6A, no offsect correction is needed

#from the m6A positions, get the position of the 5mer motif (+-2 of the m6A)
echo "getting positions of 5mers..."
bedtools slop -s -l 2 -r 2 -i $pf.bed -g /home/shareowner/genomes/human/all_chrom.fa.fai > $1.m6A.+-2.bed 

#get the sequence of the +-2 5mers !!! genomic coordinates !!!
echo "getting sequences of 5mers"
bedtools getfasta -s -tab -name -fi /home/shareowner/genomes/human/all_chrom.fa -bed $1.m6A.+-2.bed -fo $1.m6A.+-2.tab
#get the sequence of the +-30 region !!! genomic coordinates !!!
bedtools slop -s -l 1 -r -1 -i $pf.bed -g $2.fai > $pf.m6A.tmp
#from the m6A positions, get the position of the 5mer motif (+-2 of the m6A)
bedtools slop -s -l 2 -r 2 -i $pf.m6A.tmp -g $2.fai > $pf.m6A.+-2.tmp 
#get the sequence of the +-2 5mers !!! genomic coordinates !!!
bedtools getfasta -s -tab -name -fi $2 -bed $pf.m6A.+-2.tmp -fo $pf.m6A.+-2.tab
#extract k and m for calculations and filtering
awk -F "=" '{print $2 }' $pf.m6A.tmp | cut -d"]" -f1 > $pf.k.tab
# combine data 
paste $pf.m6A.tmp $pf.k.tab $pf.m6A.+-2.tab | awk '{print $1"\t"$2"\t"$3"\t"$4"_"$9"\t"$5/$7"\t"$6"\t"$7"\t"$5"\t"$9}' > $pf.m6A.txt
# get counts for unfiltered sites
nsites=$(awk 'IGNORECASE=1 $4 ~ /[GACT][GACT][GACT][GACT][GACT]$/' $pf.m6A.txt | wc -l)
Asites=$(awk 'IGNORECASE=1 $4 ~ /[GACT][GACT][A][GACT][GACT]$/' $pf.m6A.txt | wc -l)
DRACHsites=$(awk 'IGNORECASE=1 $4 ~ /[GAT][GA][A][C][ACT]$/' $pf.m6A.txt | wc -l)
ratioA=$(echo "scale=2; ($Asites*100/$nsites)" | bc )
ratioDRACH=$(echo "scale=2; ($DRACHsites*100/$Asites)" | bc)
# create summary file
printf "Unfiltered data:\n" > $pf.summary.txt 
printf "Total sites:\t$nsites\n" >> $pf.summary.txt
printf "A sites:\t$Asites\t$ratioA\t%% of total sites\n" >> $pf.summary.txt
printf "DRACH sites:\t$DRACHsites\t$ratioDRACH\t%% of A sites\n\n" >> $pf.summary.txt  

# remove additional columns to create standard bed6 file
awk 'IGNORECASE=1 $4 ~ /[GACT][GACT][A][GACT][GACT]$/ ' $pf.m6A.txt | cut -f1,2,3,4,5,6 > $pf.m6A.filtered.bed

#remove temp files
rm $pf.m6A.tmp
rm $pf.m6A.+-2.tmp
rm $pf.m6A.+-2.tab
rm $pf.m6A.txt
