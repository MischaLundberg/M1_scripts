#!/bin/sh

if [ $# -ne 4 ]
then
    echo "usage: $0 <GWAS INPUT FILE> <BIN SIZE PER THREAD> <TYPE> <LINES_TO_USE>"
    echo "<GWAS INPUT FILE might be e.g.: phenotypes.both_sexes.tsv"
    echo "<TYPE> regards binary, categorial..."
    echo "<LINES_TO_USE> how many lines of the input file shall be downloaded and processed (0==all, 1==1, 2==2, ...). can also be a comma delimited list of phenotypes to download/process (e.g. C21,G54,F45)"
    echo "Walltime is calculated by the number of lines to check times the prest \
        value of 20 minutes per line plus in the end an additional 30 minutes buffer. \
        The preset 20 minutes can be easily changed with setting the calue of \"approx_time\" to a different value."
    exit 65
fi

gwas=$1
splits=$(($2))
typ=$3
lines_to_use=$4


if [[ $lines_to_use == *".txt" ]]; then
    ## or do ls ../passed/* | awk -F "\." '{ print $4 }' | uniq | awk -vORS=\| '{ print $1 }'
    lines=$(awk -vORS=\| '{ print $1 }' $lines_to_use | head -c -1 | awk '{ print "\"^"$0"\""}' | sed 's/,/\\|^/g')
    cat_command="grep "$lines" "
elif [[ $lines_to_use == *","* ]]; then
    lines=$(echo "$lines_to_use" | awk '{ print "\"^"$0"\""}' | sed 's/,/\\|^/g')
    cat_command="grep "$lines" "
elif [[ $lines_to_use -gt 0 ]]; then
    cat_command="head -"$lines_to_use" "
else
    cat_command="cat "
fi


log_dir=~/git-repos/ctg-vl-utils/pipeline_tools/.logs/
ldsc_dir=~/downloaded/ldsc
check_integrity_dir=~/git-repos/ctg-vl-utils/pipeline_tools/
#get_commands="get_commands.txt"
variant_info="variant.info"
get_commands="get_commands_UKBB_20170915.txt"
approx_time=90 ##needed time for one file
current_dir=$(pwd)
min_buffer=30 #minimum buffer time in minutes
buffer=0

gwas_head_option=$(head -1 $gwas | cut -f3,8 | grep -E -o 'PHESANT.notes|variable_type')

##Sort variant.info
head -n 1 $variant_info > $variant_info.sorted; tail -n +2 $variant_info | sort -n -k5 >> $variant_info.sorted; mv $variant_info.sorted $variant_info


if [[ $gwas_head_option == "variable_type" ]]; then
    if [ $splits -gt 1 ]; then
        echo "Binning into binsize of $splits"
        if [ ${typ^^} == CONTINUOUS ]; then
            #echo "grep $typ $gwas | cut -d'   ' -f1,5,3 | grep $typ - | cut -f1,3 "
            grep -e $typ -e ${typ^^} $gwas | cut -d'	' -f1,5,3 | grep -e $typ -e ${typ^^} - | cut -f1,3 | split -l ${splits} - $typ.splitted. --verbose | cut -d' ' -f3 | sed 's/^‘\(.*\)’$/\1/' > $typ.all.splitted.files.txt
        else
            grep -e $typ -e ${typ^^} $gwas | cut -d'	' -f1,7,8,3 | grep -e $typ -e ${typ^^} - | cut -f1,3,4 | split -l ${splits} - $typ.splitted. --verbose | cut -d' ' -f3 | sed 's/^‘\(.*\)’$/\1/' > $typ.all.splitted.files.txt
        fi
        wall_e=$(($splits*$approx_time))
        buffer=$[$( wc -l $typ.all.splitted.files.txt | awk '{print $1}' )*1]
    else
        echo "single run"
        if [ ${typ^^} == CONTINUOUS ]; then
            grep -e $typ -e ${typ^^} $gwas | cut -f1,5,3 | grep -e $typ -e ${typ^^} - | cut -f1,3 > $typ.txt
        else
            grep -e $typ -e ${typ^^} $gwas | cut -f1,7,8,3 | grep -e $typ -e ${typ^^} - | cut -f1,3,4 > $typ.txt
        fi
        echo $typ.txt > $typ.all.splitted.files.txt
        wall_e=$[$( wc -l $typ.txt | awk '{print $1}' )*$approx_time]
        buffer=$[$( wc -l $typ.all.splitted.files.txt | awk '{print $1}' )*1]
    fi
elif [[ $gwas_head_option == "PHESANT.notes" ]]; then
    if [ $splits -gt 1 ]; then
        echo "Binning into binsize of $splits"
        if [ ${typ^^} == CONTINUOUS ]; then
            #echo "grep -e $typ -e ${typ^^} $gwas | cut -f1,3,8 | grep -e $typ -e ${typ^^} - | cut -f1,2 | split -l ${splits} - $typ.splitted. --verbose | cut -d' ' -f3 | sed 's/^‘\(.*\)’$/\1/' > $typ.all.splitted.files.txt"
            grep -e $typ -e ${typ^^} $gwas | cut -f1,3,8 | grep -e $typ -e ${typ^^} - | cut -f1,2 | split -l ${splits} - $typ.splitted. --verbose | cut -d' ' -f3 | sed 's/^‘\(.*\)’$/\1/' > $typ.all.splitted.files.txt
        else
            grep -e $typ -e ${typ^^} $gwas | cut -f1,5,6,8 | grep -e $typ -e ${typ^^} - | cut -f1,2,3 | split -l ${splits} - $typ.splitted. --verbose | cut -d' ' -f3 | sed 's/^‘\(.*\)’$/\1/' > $typ.all.splitted.files.txt
        fi
        wall_e=$(($splits*$approx_time))
        buffer=$[$( wc -l $typ.all.splitted.files.txt | awk '{print $1}' )*1]
    else
        echo "single run"
        if [ ${typ^^} == CONTINUOUS ]; then
            grep -e $typ -e ${typ^^} $gwas | cut -f1,3,8 | grep -e $typ -e ${typ^^} - | cut -f1,2 > $typ.txt
        else
            grep -e $typ -e ${typ^^} $gwas | cut -f1,5,6,8 | grep -e $typ -e ${typ^^} - | cut -f1,2,3 > $typ.txt
        fi
        echo $typ.txt > $typ.all.splitted.files.txt
        wall_e=$[$( wc -l $typ.txt | awk '{print $1}' )*$approx_time]
        buffer=$[$( wc -l $typ.all.splitted.files.txt | awk '{print $1}' )*1]
    fi
else
    echo "Please check input file. Consider editing this download script"
fi

counter_files=$( cat $typ.all.splitted.files.txt )
if [[ $buffer -lt $min_buffer ]]; then
    buffer=$min_buffer
fi

wall_e=$(($(($wall_e+$buffer))*60))
((h=((${wall_e}-${wall_e}%3600)/3600)))
((m=(${wall_e}%3600)/60))
((s=0))

wall=$(echo "$h:$m:$s")
echo "Estimated runtime: $h:$m:$s"
DATE=`date '+%Y-%m-%d_%H-%M-%S'`

if [ -f $typ.download.log ]; then
    mv $typ.download.log $typ.$DATE.download.log 
fi

if [ -f passed ]; then
    mkdir passed    
fi
if [ -f failed ]; then
    mkdir failed    
fi
binary_categorical_dl (){
for file in $counter_files; do   
    echo "#PBS -A tb
#PBS -V
#PBS -j oe
#PBS -l select=1:ncpus=2:mem=4G
#PBS -l walltime=$wall
#PBS -N dl_`echo $file | cut -c1-3`_gwas
#PBS -o $log_dir`echo $file`_$DATE

cd $current_dir
$cat_command $file | while read pheno Ncase Ncon; do
    grep -w \$pheno $get_commands | while read com1 com2 com3 com4; do
        \$(\$com1 \$com2 \$com3 \$com4)
        zcat \$com4 > \$com4.1
        head -n 1 \$com4.1 > \$com4.1.sorted; tail -n +2 \$com4.1 | sort -n -k2 >> \$com4.1.sorted; mv \$com4.1.sorted \$com4.1
        var_lines=\$(wc -l $variant_info)
        phen_lines=\$(wc -l \$com4.1)
        if [[ \$phen_lines -ne \$var_lines ]]; then
            echo \"There is a inconsistency between your variant.info ($variant_info == \$var_lines lines) and your phenotype file (\$com4.1 == \$phen_lines lines). Please check this!\" >> $typ.download.log
            exit 0
        fi
        paste \$com4.1 $variant_info > \$com4
        header=\$(head -1 \$com4 | awk '{ if (\$10 == \"chr\") print \"10\";else if (\$12 == \"chr\") print \"12\"; else if (\$13 == \"chr\") print \"13\"; else print \"0\"}' -)
        if [[ \$header == 13 ]]; then
            grep -Fvw NaN \$com4 | grep -Fvw true | awk '{if(\$13==\"chr\") print \"CHR BP A2 A1 SNP FREQ BETA SE P\"; else print \$13,\$14,\$15,\$16,\$17,\$18,\$9,\$10,\$12}' > $typ.\${pheno}.txt
        elif [[ \$header == 12 ]]; then
            grep -Fvw NaN \$com4 | grep -Fvw true | awk '{if(\$12==\"chr\") print \"CHR BP A2 A1 SNP FREQ BETA SE P\"; else print \$12,\$13,\$14,\$15,\$16,\$17,\$8,\$9,\$11}' > $typ.\${pheno}.txt
        elif [[ \$header == 10 ]]; then
            grep -Fvw NaN \$com4 | grep -Fvw true | awk '{if(\$10==\"chr\") print \"CHR BP A2 A1 SNP FREQ BETA SE P\"; else print \$10,\$11,\$12,\$13,\$14,\$15,\$6,\$7,\$9}' > $typ.\${pheno}.txt
        else
            echo \"please check input and change script towards your input file\"
            echo \"current header (\$header) looks like: \$(head -1 \$com4)\"
            exit 0
        fi
        source activate ldsc
        python $ldsc_dir/munge_sumstats.py --sumstats $typ.\${pheno}.txt --out $typ.\${pheno} --merge-alleles $ldsc_dir/w_hm3.snplist --N-cas \$Ncase --N-con \$Ncon
        python $ldsc_dir/ldsc.py --h2 $typ.\${pheno}.sumstats.gz --ref-ld-chr $ldsc_dir/eur_w_ld_chr/ --w-ld-chr $ldsc_dir/eur_w_ld_chr/ --out $typ.\${pheno}
        awk '\$9<=0.05' $typ.\${pheno}.txt > $typ.\${pheno}.txt.significant
        heritability=\$(grep \"Heritability test\" $typ.\${pheno}.log)
        if ! grep -q \"Heritability test passed\" $typ.\${pheno}.log; then
            if [[ \${heritability} == \"\" ]] || [[ \${heritability} == \" \" ]]; then
                heritability=\"No heritability calculation possible\"
            fi
            echo \"\${heritability} for file $typ.\${pheno}.txt ... For more information check $typ.\${pheno}.log\" >> $typ.download.log
            #mv $typ.\${pheno}.*sumstats* failed/
            mv $typ.\${pheno}*txt.significant failed/
            mv $typ.\${pheno}.log failed/
        else
            echo \"\${heritability} for file $typ.\${pheno}.txt\" >> $typ.download.log
            sh ${check_integrity_dir}/check_integrity.sh $typ.\${pheno}.txt
            #if [ -f $typ.\${pheno}.log ]; then
            #    rm $typ.\${pheno}.log
            #fi
            mv $typ.\${pheno}.*checked* passed/
            mv $typ.\${pheno}.*sumstats* passed/
            mv $typ.\${pheno}*txt.significant passed/
            mv $typ.\${pheno}.log passed/
        fi
        #if [ -f $typ.\${pheno}.txt ]; then
        #    rm $typ.\${pheno}.txt
        #fi
        #rm \$com4.1
        #rm \$com4
        #rm $typ.\${pheno}.sumstats.gz
        #rm $typ.\${pheno}.*
        conda deactivate
    done
done
#rm $file
#rm $file.sh
2>&1" > $file.sh; qsub $file.sh;
done
}

others_dl () {
for file in $counter_files; do
    echo "#PBS -A tb
#PBS -V
#PBS -j oe
#PBS -l select=1:ncpus=2:mem=4G
#PBS -l walltime=$wall
#PBS -N dl_`echo $file | cut -c1-3`_gwas
#PBS -o $log_dir`echo $file`_$DATE
cd $current_dir
$cat_command $file | while read pheno N; do
    grep -w \$pheno $get_commands | while read com1 com2 com3 com4; do
        \$(\$com1 \$com2 \$com3 \$com4)
        zcat \$com4 > \$com4.1
        head -n 1 \$com4.1 > \$com4.1.sorted; tail -n +2 \$com4.1 | sort -n -k2 >> \$com4.1.sorted; mv \$com4.1.sorted \$com4.1
        var_lines=\$(wc -l $variant_info)
        phen_lines=\$(wc -l \$com4.1)
        if [[ \$phen_lines -ne \$var_lines ]]; then
            echo \"There is a inconsistency between your variant.info ($variant_info == \$var_lines lines) and your phenotype file (\$com4.1 == \$phen_lines lines). Please check this!\" >> $typ.download.log
            exit 0
        fi
        paste \$com4.1 $variant_info > \$com4
        header=\$(head -1 \$com4 | awk '{ if (\$10 == \"chr\") print \"10\"; else if (\$12 == \"chr\") print \"12\"; else if (\$13 == \"chr\") print \"13\"; else print \"0\"}' - )
        if [[ \$header == 13 ]]; then
            grep -Fvw NaN \$com4 | grep -Fvw true | awk '{if(\$13==\"chr\") print \"CHR BP A2 A1 SNP FREQ BETA SE P\"; else print \$13,\$14,\$15,\$16,\$17,\$18,\$9,\$10,\$12}' > $typ.\${pheno}.txt
        elif [[ \$header == 12 ]]; then
            grep -Fvw NaN \$com4 | grep -Fvw true | awk '{if(\$12==\"chr\") print \"CHR BP A2 A1 SNP FREQ BETA SE P\"; else print \$12,\$13,\$14,\$15,\$16,\$17,\$8,\$9,\$11}' > $typ.\${pheno}.txt
        elif [[ \$header == 10 ]]; then
            grep -Fvw NaN \$com4 | grep -Fvw true | awk '{if(\$10==\"chr\") print \"CHR BP A2 A1 SNP FREQ BETA SE P\"; else print \$10,\$11,\$12,\$13,\$14,\$15,\$6,\$7,\$9}' > $typ.\${pheno}.txt
        else
            echo \"please check input and change script towards your input file\"
            echo \"current header (\$header) looks like: \$(head -1 \$com4)\"
            exit 0
        fi
        source activate ldsc
        python $ldsc_dir/munge_sumstats.py --sumstats $typ.\${pheno}.txt --out $typ.\${pheno} --merge-alleles $ldsc_dir/w_hm3.snplist --N \${N}
        python $ldsc_dir/ldsc.py --h2 $typ.\${pheno}.sumstats.gz --ref-ld-chr $ldsc_dir/eur_w_ld_chr/ --w-ld-chr $ldsc_dir/eur_w_ld_chr/ --out $typ.\${pheno}
        awk '\$9<=0.05' $typ.\${pheno}.txt > $typ.\${pheno}.txt.significant
        heritability=\$(grep \"Heritability test\" $typ.\${pheno}.log)
        if ! grep -q \"Heritability test passed\" $typ.\${pheno}.log; then
            if [[ \${heritability} == \"\" ]] || [[ \${heritability} == \" \" ]]; then
                heritability=\"No heritability calculation possible\"
            fi
            echo \"\${heritability} for file $typ.\${pheno}.txt ... For more information check $typ.\${pheno}.log\" >> $typ.download.log
            #mv $typ.\${pheno}.*sumstats* failed/
            mv $typ.\${pheno}*txt.significant failed/
            mv $typ.\${pheno}.log failed/
        else
            echo \"\${heritability} for file $typ.\${pheno}.txt\" >> $typ.download.log
            sh ${check_integrity_dir}/check_integrity.sh $typ.\${pheno}.txt
            #if [ -f $typ.\${pheno}.log ]; then
            #    rm $typ.\${pheno}.log
            #fi
            mv $typ.\${pheno}.*checked* passed/
            mv $typ.\${pheno}.*sumstats* passed/
            mv $typ.\${pheno}*txt.significant passed/
            mv $typ.\${pheno}.log passed/
        fi
        #if [ -f $typ.\${pheno}.txt ]; then
        #    rm $typ.\${pheno}.txt
        #fi
        #rm \$com4.1
        #rm \$com4
        #rm $typ.\${pheno}.sumstats.gz
        #rm $typ.\${pheno}.*
        conda deactivate
    done
done
#rm $file
#rm $file.sh
2>&1" > $file.sh; qsub $file.sh;
done
}


if [[ $typ == binary ]] || [[ $typ == categorical ]]; then
    binary_categorical_dl
else
    others_dl
fi

#for file in $counter_files; do
#    rm $file*
#done
