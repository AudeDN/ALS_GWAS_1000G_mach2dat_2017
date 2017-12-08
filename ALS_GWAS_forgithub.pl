#!/usr/local/bin/perl


################################################################################
# This perl script is a pipeline, designed to submit swarm jobs to National Ins-
# titutes of Health HPC Biowulf cluster, written in December 2017 by Aude
# Nicolas (aude.nicolas.deydier@gmail.com) on behalf of Bryan Traynor,
# Laboratory of Neurogenetics, National Institute for Aging, 
# National Institutes of Health.
#
# This pipeline is divided into 3 steps:
# First step is the imputation (SHAPEIT, minimac3)
# Second Step is the association analysis (mach2dat)
# Third step is the Meta-analysis (METAL)

################################################################################

# Variables
 
$minimac3 = "/Minimac3/bin/Minimac3-omp"; #full path of minimac executable
$haps_location = "/1kG.phase3v5";
$shapeit_memory = "120";    # gigabytes of memory required for a single Mach job
$minimac_memory = "120"; # gigabytes of memory required for a single minimac-omp multi-threaded job


$|=1;
$PID = "$$";
$prefix = "genetic_map_chr";
$suffix = "_combined_b37.txt";
$path ="/pathtobinary";
$binaryfiles ="binaryfiles";
$pathtogeneticmaps="/pathtogeneticmaps";
$datfile ="/pathtodatfile";
$phenofile ="/pathtophenofile";
$metalscript ="/pathtoscriptwith commandsformeta-analysis.txt";


$swarmfile1 = "shapeit_phase.$PID.swarm";
$swarmfile2 = "shapeit_convert.$PID.swarm";
$swarmfile3 = "minimac3.$PID.swarm";
$swarmfile4a = "formatingforassociation.$PID.swarm";
$swarmfile4b = "tabix4b.$PID.swarm";
$swarmfile5 = "mach2dat3PCs.$PID.swarm";
$swarmfile6 = "formattingfilespart1formeta.$PID.swarm";
$swarmfile7 = "formattingfilespart2formeta.$PID.swarm";
$swarmfile8 = "formattingfilespart3formeta.$PID.swarm";
$swarmfile9 = "metal.$PID.swarm";

print "$PID";


################################################################################
# Step 1- Imputation which needs a reference panel, the genetic maps an the 
# binary files
################################################################################

## 1- Phasing step - Estimation of the haplotypes for all individuals - 
# using SHAPEIT version 2.r790


open (SWFILE1, ">$swarmfile1"); 

foreach $chr (1..22)
{

    print SWFILE1 "shapeit -T 10 --input-bed $path/$binaryfiles.chr$chr.bed $path/$binaryfiles.chr$chr.bim $path/$binaryfiles.chr$chr.chr$chr.fam --input-map /pathtogeneticmaps/$prefix$chr$suffix -O $path/$binaryfiles.chr$chr.phased --output-log $path/$binaryfiles.chr$chr.phased.log\n";

}

close (SWFILE1);

##############
$jobnum1 =`swarm -f $swarmfile1 --module shapeit -t 10 -g $shapeit_memory --partition=norm --time=10-00:00:00`;
$jobstr1 = substr($jobnum1, -9, 8);
print "=====Phasing haplotypes with JobID $jobstr1.===\n";


################################################################################

## 2- Conversion SHAPEIT's ouputs in vcf (Minimac3 needs an input in a vcf format)

open (SWFILE2, ">$swarmfile2");

foreach $chr (1..22)
{

	print SWFILE2 "shapeit -convert --input-haps $path/$binaryfiles.chr$chr.phased --output-vcf $path/$binaryfiles.chr$chr.phased.vcf\n";

}

close (SWFILE2);

##############
$jobnum2 =`swarm -f $swarmfile2 -t 1 -g $shapeit_memory --module shapeit --time=24:00:00 --dependency=afterany:$jobstr1`;
$jobstr2 = substr($jobnum2, -9, 8);
print "=====Converting to VCF with JobID $jobstr2.=======\n";

################################################################################

#3 - Imputation into phased haplotypes using minimac3 (version 1.0.11) and
# 1000 Genomes Project dataset (phase 3, version 5a, release 2013-05-02, 1000genomes.org)

open (SWFILE3, ">$swarmfile3");

foreach $chr (1..22)
{

    $haps = "$haps_location/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
    print SWFILE3 "$minimac3 --cpus 10 --refHaps $haps --rounds 5 --states 200 --haps $path/$binaryfiles.chr$chr.phased.vcf --doseOutput  --prefix $path/$binaryfiles.chr$chr.phased.imputed > $path/$binaryfiles.chr$chr.phased.imputed.minimac3.log 2>&1\n";

}

close (SWFILE3);

##############
$jobnum3 =`swarm -f $swarmfile3 -t 10 -g $minimac_memory --partition=norm --time=10-00:00:00 --dependency=afterany:$jobstr2`;
$jobstr3 = substr($jobnum3, -9, 8);
print "Imputing with JobID $jobstr3.\n"; 

################################################################################
# Second step - Association
################################################################################

## 1- Formatting the dosage files and info files

open (SWFILE4a, ">$swarmfile4a");

foreach $chr (1..22)
{

	print SWFILE4a "zcat $path/$binaryfiles.chr$chr.phased.imputed.dose.gz | awk -F\"\t\" '{print \$1\"->\"\$1\"\t\"\$0}' | cut --complement -f2 | bgzip -c > $path/$binaryfiles.chr$chr.phased.imputed.dose.mach2datcompatible.gz\n";

}

close (SWFILE4a);

##############
$jobnum4a =`swarm -f $swarmfile4a -g 120  --partition=norm  --time=01-00:00:00 --dependency=afterany:$jobstr3`;
$jobstr4a = substr($jobnum4a, -9, 8);
print "=====FOrmatting dosage files with JobID $jobstr4a.===\n";

################################################################################

open (SWFILE4b, ">$swarmfile4b");

foreach $chr (1..22)
{

	print SWFILE4b "tabix -p vcf $path/$binaryfiles.chr$chr.phased.imputed.dose.mach2datcompatible.gz\n";

}

close (SWFILE4b);

##############
$jobnum4b =`swarm -f $swarmfile4b -g 40  --partition=quick  --time=00-02:00:00 --module Tabix --dependency=afterany:$jobstr4a`;
$jobstr4b = substr($jobnum4b, -9, 8);
print "=====Splitting big dosage files with JobID $jobstr4b.===\n";

################################################################################
# Formating the info files > check the Rsq column is the good one

################################################################################

## 2- Association using mach2dat software (version 1.0.24)

open (SWFILE5, ">$swarmfile5"); 


foreach $chr (1..22)
{

	print SWFILE5 "mach2dat -d $datfile -p $phenofile -i $path/$binaryfiles.chr$chr.phased.imputed.dose.info_compmach2dat --dosefile $path/$binaryfiles.chr$chr.phased.imputed.dose.mach2datcompatible.gz --frequency --verbosesamplesize --rsqcutoff 0.3 > $path/$binaryfiles.chr$chr.mach2dat.out\n";


}

close (SWFILE5);

##############
$jobnum5 =`swarm -f $swarmfile5 -g 240 --partition=norm --time=10-00:00:00 --module mach2dat --dependency=afterany:$jobstr4e`;
$jobstr5 = substr($jobnum5, -9, 8);
print "=====Running mach2dat association with imputed variants JobID $jobstr5.===\n";


################################################################################
# 3a- FOrmatting the files for the meta-analysis


open (SWFILE6, ">$swarmfile6"); 

foreach $chr (1..22)
{

	print SWFILE6 "echo -e \"SNP\tMEANFREQ1\tN\tTRAIT\tCHR\tBP\tALLELE1\tALLELE2\tFREQCI1\tFREQCI2\tRSQR\tEFFECT1\tOR\tSTDERR\tWALDCHISQ\tPVALUE\tLRCHISQ\tLRPVAL\tNCASES\tNCONTROLS\" > $path/$binaryfiles.chr$chr.mach2dat.out.tomerge\n";
	print SWFILE6 "echo -e \"SNP\tMARKER\tMEANFREQ1\tN\tTRAIT\tCHR\tPOS\tALLELE1\tALLELE2\tFREQ1\tFREQCI1\tFREQCI2\tRSQR\tEFFECT1\tOR\tSTDERR\tWALDCHISQ\tP\tLRCHISQ\tLRPVAL\tNCASES\tNCONTROLS\" > $path/$binaryfiles.chr$chr.mach2dat.out.formeta\n";

}

close (SWFILE6);

##############
$jobnum6 =`swarm -f $swarmfile6 -g 60 --partition=quick  --time=00-00:05:00  --dependency=afterany:$jobstr5`;
$jobstr6 = substr($jobnum6, -9, 8);
print "=====Formattingfiles for GWAS part1 JobID $jobstr6.===\n";

################################################################################
# 3b- FOrmatting the files for the meta-analysis


open (SWFILE7, ">$swarmfile7"); 

foreach $chr (1..22)
{

	print SWFILE7 "cat $path/$binaryfiles.chr$chr.mach2dat.out | grep \"^ALS\" | sed -e 's|:| |' | sed -e 's|,| |g' | sed -e 's|\\.[0-9][0-9][0-9][0-9]|  &|'  |  sed -e 's|(| |' |  sed -e 's|)| |'  | sed -E \"s/[[:blank:]]+/\$(printf '\t')/g\" | awk -F\"\t\" '{print \$2\":\"\$3\"\t\"\$2\":\"\$3\":\"\$4\":\"\$5\"\t\"(\$7+\$8)\/2\"\t\"4\/((1\/\$17)+(1\/\$18))\"\t\"\$0}'  > $path/$binaryfiles.chr$chr.mach2dat.out.tomerge\n";

}

close(SWFILE7);

##############
$jobnum7 =`swarm -f $swarmfile7 -g 60 --partition=quick --time=00-00:30:00 --dependency=afterany:$jobstr6`;
$jobstr7 = substr($jobnum7, -9, 8);
print "=====FOrmatting files for GWAS part2 JobID $jobstr7.===\n";

################################################################################
# 3c- Formatting the files for the meta-analysis
# concatenation of all the chromosomes in one file and filtering the SNPs with a
# |Beta coeffincent| < 5


open (SWFILE8, ">$swarmfile8");

print SWFILE8 qq[cat $path/$binaryfiles.chr1.mach2dat.out.tomerge $path/$binaryfiles.chr2.mach2dat.out.tomerge $path/$binaryfiles.chr3.mach2dat.out.tomerge $path/$binaryfiles.chr4.mach2dat.out.tomerge $path/$binaryfiles.chr5.mach2dat.out.tomerge $path/$binaryfiles.chr6.mach2dat.out.tomerge $path/$binaryfiles.chr7.mach2dat.out.tomerge $path/$binaryfiles.chr8.mach2dat.out.tomerge $path/$binaryfiles.chr9.mach2dat.out.tomerge $path/$binaryfiles.chr10.mach2dat.out.tomerge $path/$binaryfiles.chr11.mach2dat.out.tomerge $path/$binaryfiles.chr12.mach2dat.out.tomerge $path/$binaryfiles.chr13.mach2dat.out.tomerge $path/$binaryfiles.chr14.mach2dat.out.tomerge $path/$binaryfiles.chr15.mach2dat.out.tomerge $path/$binaryfiles.chr16.mach2dat.out.tomerge $path/$binaryfiles.chr17.mach2dat.out.tomerge $path/$binaryfiles.chr18.mach2dat.out.tomerge $path/$binaryfiles.chr19.mach2dat.out.tomerge $path/$binaryfiles.chr20.mach2dat.out.tomerge $path/$binaryfiles.chr21.mach2dat.out.tomerge $path/$binaryfiles.chr22.mach2dat.out.tomerge | grep -v "^SNP" | awk -F\"\t\" '{if (\$14 < 5 && \$14 > -5 ) print \$0}' >> /data/ALS_50k/5ksamples/GWASTempfiles/Mach2datresults/Italians/4PCs/italians.allchr.imputed.mach2dat_Rsq0.0.3PCs.out.formeta\n];

}

close (SWFILE8);

##############
$jobnum8 =`swarm -f $swarmfile8 -g 60 --partition=quick  --time=00-00:05:00  --dependency=afterany:$jobstr7`;
$jobstr8 = substr($jobnum8, -9, 8);
print "=====Formattingfiles for GWAS part3 $jobstr8.===\n";


################################################################################
##3- Meta-analysis
################################################################################



open (SWFILE9, ">$swarmfile9"); 

print SWFILE9 "metal << SOURCE $metalscript \n";

close(SWFILE9);

##############
$jobnum9 =`swarm -f $swarmfile9 -g 10 --partition=norm --time=00-03:00:00 --module metal --dependency=afterany:$jobstr8`;
$jobstr9 = substr($jobnum9, -9, 8);
print "=====Running the Meta-analysis JobID $jobstr9.===\n";


