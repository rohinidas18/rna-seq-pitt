#!/bin/bash

# this is the bash script for downloading the reference transcriptome followed by the primary assembly
# then activating salmon to create a decoy-aware index of the transcriptome
# and finally quantifying the reads per gene per sample (since we have only one) for the 2 species

#curl ftp.ensembl.org/pub/release-106/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz -o mmus.fa.gz

#curl ftp.ensembl.org/pub/release-106/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz -o mmus_pi.fa.gz

#grep "^>" <(gunzip -c mmus_pi.fa.gz) | cut -d " " -f 1 > decoys.txt
#sed -i.bak -e 's/>//g' decoys.txt

#cat mmus.fa.gz mmus_pi.fa.gz > mmus.fa.gz



#insert source of conda in the system if using miniconda to activate salmon
#source /home/rohini/miniconda3/etc/profile.d/conda.sh  
#conda activate salmon

#else run salmon using the full path to the binary executable file
export PATH="/usr/local/bin/salmon-1.8.0_linux_x86_64/bin:$PATH"
salmon index -t mmus.fa.gz -d decoys.txt -p 12 -i mmus_index --keepDuplicates

for fn in data/F2_KO_S23;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i mmus_index -l A \
         -1 ${fn}/${samp}_R1_001.fastq.gz \
         -2 ${fn}/${samp}_R2_001.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}_quant
done

for fn in data/F7_WT_S22;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i mmus_index -l A \
         -1 ${fn}/${samp}_R1_001.fastq.gz \
         -2 ${fn}/${samp}_R2_001.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}_quant
done
