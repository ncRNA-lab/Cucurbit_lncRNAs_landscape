#!/bin/bash

####### FUNCTIONS

task_Fastqc1(){
	echo -e "\n$1..."
	
	if [ -e $1"_1.fastq.gz" ] && [ -e $1"_2.fastq.gz" ]; then
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
	
	if [ -e $2/$1"_1.fastq.gz" ] && [ -e $2/$1"_2.fastq.gz" ]; then
		## Trimming the libraries.
		echo -e "\nTrimming..."
		fastp \
			--adapter_fasta $4/adapters.fa \
			--thread $5 \
			-i $2/$1"_1.fastq.gz" \
			-I $2/$1"_2.fastq.gz" \
			-o $3/$1"_tr_1P.fastq.gz" \
			-O $3/$1"_tr_2P.fastq.gz" \
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
	elif [ -e $2/$1".fastq.gz" ]; then
		## Trimming the libraries.
		echo -e "\nTrimming..."
		fastp \
			--adapter_fasta $4/adapters.fa \
			--thread $5 \
			-i $2/$1".fastq.gz" \
			-o $3/$1"_tr.fastq.gz" \
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

task_Fastqc2(){
	echo -e "\n"$1"..."
	
	if [ -e $1"_tr_1P.fastq.gz" ] && [ -e $1"_tr_2P.fastq.gz" ]; then
		## Execute fastqc.
		echo -e "\nExecute fastqc..."
		fastqc $1"_tr_1P.fastq.gz"
		fastqc $1"_tr_2P.fastq.gz"
		## Save fastqc results in fastqc folder.
		echo -e "\nMove fastqc files to fastqc folder..."
		mv $1"_tr_1P_fastqc.zip" ./fastqc
		mv $1"_tr_1P_fastqc.html" ./fastqc
		mv $1"_tr_2P_fastqc.zip" ./fastqc
		mv $1"_tr_2P_fastqc.html" ./fastqc
	elif [ -e $1"_tr.fastq.gz" ]; then
		## Execute fastqc.
		echo -e "\nExecute fastqc..."
		fastqc $1"_tr.fastq.gz"
		## Save fastqc results in fastqc folder.
		echo -e "\nMove fastqc files to fastqc folder..."
		mv $1"_tr_fastqc.zip" ./fastqc
		mv $1"_tr_fastqc.html" ./fastqc
	else
		echo -e "\n$1 doesn't exist. Check it."
	fi
}

task_Summary_table_trimming(){
	echo -e "\n"$1"..."
	
	if [ -e $1"_tr_1P.fastq.gz" ] && [ -e $1"_tr_2P.fastq.gz" ]; then
		## Calculate the library depth.
		echo -e "\nCalculate the library depth..."
		total_lines_1=$(zcat $1"_tr_1P.fastq.gz" | wc -l)
		num_sequences_1=$(($total_lines_1/4))
		total_lines_2=$(zcat $1"_tr_2P.fastq.gz" | wc -l)
		num_sequences_2=$(($total_lines_2/4))
		## Check if the libraries are empty files.
		echo -e "\nCheck if the libraries are empty files..."
		if [ $num_sequences_1 -ge 1000000 ] && [ $num_sequences_2 -ge 1000000 ]; then
			echo -e $1" passes the trimmomatic step."
			echo -e $1"\t"$2"\tpasses the trimmomatic step\tYES\t$num_sequences_1/$num_sequences_2" >> $3/$1.tsv
			echo -e $1 >> $4/$1.txt
		elif [ $num_sequences_1 -ge 1000000 ] || [ $num_sequences_2 -ge 1000000 ]; then
			echo -e $1" passes the filter only forward or reverse)."
			echo -e $1"\t"$2"\tpasses the filter only forward or reverse\tNO\t$num_sequences_1/$num_sequences_2" >> $3/$1.tsv
		else
			echo -e $1" doesn't pass the trimmomatic step."
			echo -e $1"\t"$2"\tdoesn't pass the trimmomatic step\tNO\t$num_sequences_1/$num_sequences_2" >> $3/$1.tsv
		fi
	elif [ -e $1"_tr.fastq.gz" ]; then
		## Calculate the library depth.
		echo -e "\nCalculate the library depth..."
		total_lines=$(zcat $1"_tr.fastq.gz" | wc -l)
		num_sequences=$(($total_lines/4))
		## Check if the library is an empty file.
		echo -e "\nCheck if the library is an empty file..."
		if [ $num_sequences -ge 1000000 ]; then
			echo -e $1" passes the trimmomatic step."
			echo -e $1"\t"$2"\tpasses the trimmomatic step\tYES\t$num_sequences" >> $3/$1.tsv
			echo -e $1 >> $4/$1.txt
		else
			echo -e $1" doesn't pass the trimmomatic step."
			echo -e $1"\t"$2"\tdoesn't pass the trimmomatic step\tNO\t$num_sequences" >> $3/$1.tsv
		fi
	else
		echo -e "\n$1 doesn't exist. Check it."
	fi
}

task_Pseudomapping(){
	echo -e "\n"$1"..."
	
	if [ -e $2/$1"_tr_1P.fastq.gz" ] && [ -e $2/$1"_tr_2P.fastq.gz" ]; then
		echo -e "\nPseudomapping (Libraries vs reference transcriptome)..."
		salmon quant -i $3/02-Index/$4 -l A \
		-1 $2/$1"_tr_1P.fastq.gz" -2 $2/$1"_tr_2P.fastq.gz" \
		-g $5/GFF3_genes/$4.gff3 \
		-p $6 -o $1 --gcBias
	elif [ -e $2/$1"_tr.fastq.gz" ]; then
		echo -e "\nPseudomapping (Library vs reference transcriptome)..."
		salmon quant -i $3/02-Index/$4 -l A \
	    	-r $2/$1"_tr.fastq.gz" \
	    	-g $5/GFF3_genes/$4.gff3 \
	    	-p $6 -o $1 --gcBias
	else
		echo -e "\n$1 doesn't exist. Check it."
	fi
}

task_Select_SS_libraries(){
	
	if [ -e $3/$1"_tr_1P.fastq.gz" ] && [ -e $3/$1"_tr_2P.fastq.gz" ]; then
		echo -e "Move the strand-specific library $1 to 04-Selected_data..."
		cp $3/$1"_tr_1P.fastq.gz" ./$2"_tr_1P.fastq.gz"
		cp $3/$1"_tr_2P.fastq.gz" ./$2"_tr_2P.fastq.gz"
	elif [ -e $3/$1"_tr.fastq.gz" ]; then
		echo -e "Move the strand-specific library $1 to 04-Selected_data..."
		cp $3/$1"_tr.fastq.gz" ./$2"_tr.fastq.gz"
	else
		echo -e "\n$1 doesn't exist. Check it."
	fi
}

"$@"

