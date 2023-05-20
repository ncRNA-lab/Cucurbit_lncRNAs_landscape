#!/bin/bash

####### FUNCTIONS

task_Repair(){
	echo -e "\n"$1"..."
	
	if [ -e $1"_1.fastq.gz" ] && [ -e $1"_2.fastq.gz" ]; then
		## Repair.
		repair.sh in1=$1_1.fastq.gz in2=$1_2.fastq.gz out1=$1_1_repaired.fastq.gz out2=$1_2_repaired.fastq.gz outs=$1_single.fastq.gz repair
		rm $1_single.fastq.gz
		rm $1_1.fastq.gz
		rm $1_2.fastq.gz
		mv $1_1_repaired.fastq.gz $1_1.fastq.gz
		mv $1_2_repaired.fastq.gz $1_2.fastq.gz
	elif [ -e $1".fastq.gz" ]; then
		## Repair.
		repair.sh in1=$1.fastq.gz out1=$1_repaired.fastq.gz outs=$1_single.fastq.gz repair
		rm $1_single.fastq.gz
		rm $1.fastq.gz
		mv $1_repaired.fastq.gz $1.fastq.gz
	else
		echo -e "\n$1 doesn't exist. Check it."
	fi
}

task_Fastqc1(){
	echo -e "\n$1..."
	
	if [ -e $1"_1.fastq.gz" ] && [ -e $1"_2.fastq.gz" ]; then
		rm ./fastqc/$1_1_fastqc.html
		rm ./fastqc/$1_2_fastqc.html
		rm ./fastqc/$1_1_fastqc.zip
		rm ./fastqc/$1_2_fastqc.zip
		## Execute fastqc.
		echo -e "\nExecute fastqc..."
		fastqc $1"_1.fastq.gz"
		fastqc $1"_2.fastq.gz"
		## Save fastqc results in fastqc folder.
		echo -e "\nMove fastqc files to fastqc folder..."
		mv $1"_1_fastqc.zip" ./fastqc
		mv $1"_1_fastqc.html" ./fastqc
		mv $1"_2_fastqc.zip" ./fastqc
		mv $1"_2_fastqc.html" ./fastqc
	elif [ -e $1".fastq.gz" ]; then
		rm ./fastqc/$1_fastqc.html
		rm ./fastqc/$1_fastqc.zip
		## Execute fastqc.
		echo -e "\nExecute fastqc..."
		fastqc $1".fastq.gz"
		## Save fastqc results in fastqc folder.
		echo -e "\nMove fastqc files to fastqc folder..."
		mv $1"_fastqc.zip" ./fastqc
		mv $1"_fastqc.html" ./fastqc
	else
		echo -e "\n$1 doesn't exist. Check it."
	fi
}

task_Trimming(){
	echo -e "\n"$1"..."
	
	if [ -e $1"_tr_1P.fastq.gz" ] && [ -e $1"_tr_2P.fastq.gz" ]; then
		echo -e "\n$1 exists but it's wrong. These files were created in step 2.1 before repairing the files. Removing..."
		rm $1"_tr_1P.fastq.gz"
		rm $1"_tr_2P.fastq.gz"
	elif [ -e $1"_tr.fastq.gz" ]; then
		echo -e "\n$1 exists but it's wrong. This file was created in step 2.1 before repairing the file. Removing..."
		rm $1"_tr.fastq.gz"
	else
		echo -e "\n$1 doesn't exist. Finally, these file/files weren't created in step 2.1 before repairing the files. Continue..."
	fi
	
	if [ -e ../01-Raw_data/$1"_1.fastq.gz" ] && [ -e ../01-Raw_data/$1"_2.fastq.gz" ]; then
		## Trimming the libraries.
		echo -e "\nTrimming..."
		fastp \
			--adapter_fasta $2/adapters.fa \
			--thread $3 \
			-i ../01-Raw_data/$1"_1.fastq.gz" \
			-I ../01-Raw_data/$1"_2.fastq.gz" \
			-o $1"_tr_1P.fastq.gz" \
			-O $1"_tr_2P.fastq.gz" \
			--correction \
			--cut_right \
			--cut_right_window_size 4 \
			--cut_right_mean_quality 20 \
			--cut_front \
			--cut_front_window_size 1 \
			--cut_front_mean_quality 3 \
			--cut_tail \
			--cut_tail_window_size 1 \
			--cut_tail_mean_quality 3 \
			--length_required 49 \
			--trim_poly_x \
			--poly_x_min_len 10
	elif [ -e ../01-Raw_data/$1".fastq.gz" ]; then
		## Trimming the libraries.
		echo -e "\nTrimming..."
		fastp \
			--adapter_fasta $2/adapters.fa \
			--thread $3 \
			-i ../01-Raw_data/$1".fastq.gz" \
			-o $1"_tr.fastq.gz" \
			--cut_right \
			--cut_right_window_size 4 \
			--cut_right_mean_quality 20 \
			--cut_front \
			--cut_front_window_size 1 \
			--cut_front_mean_quality 3 \
			--cut_tail \
			--cut_tail_window_size 1 \
			--cut_tail_mean_quality 3 \
			--length_required 49 \
			--trim_poly_x \
			--poly_x_min_len 10
	else
		echo -e "\n$1 doesn't exist. Check it."
	fi
}

"$@"

