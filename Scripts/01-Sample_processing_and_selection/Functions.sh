#!/bin/bash

####### FUNCTIONS

task_Prefetch(){
	echo -e "\n$1..."
	
	## Download the SRA file.
	echo -e "\nExecute prefetch..."
    	prefetch $1 --verify yes --check-all
	mv ./$1/$1.sra .
	rmdir $1
	## Extract the fastq file (SINGLE) or files (PAIRED).
	echo -e "\nExecute fasterq-dump..."
	fasterq-dump --split-files $1.sra
	rm $1.sra
}

task_Fastqc1(){
	echo -e "\n$1..."
	
	if [ -e $1"_1.fastq" ] && [ -e $1"_2.fastq" ]; then
		## Compress fastq files to fastq.gz.
		echo -e "\nCompress the libraries..."
		pigz -p $2 $1"_1.fastq"
		pigz -p $2 $1"_2.fastq"
		## Execute fastqc.
		echo -e "\nExecute fastqc..."
		fastqc $1"_1.fastq.gz"
		fastqc $1"_2.fastq.gz"
		## Save fastqc results in Fastqc folder.
		echo -e "\nMove fastqc files to Fastqc folder..."
		mv $1"_1_fastqc.zip" ./Fastqc
		mv $1"_1_fastqc.html" ./Fastqc
		mv $1"_2_fastqc.zip" ./Fastqc
		mv $1"_2_fastqc.html" ./Fastqc
	elif [ -e $1".fastq" ]; then
		## Compress fastq files to fastq.gz.
		echo -e "\nCompress the library..."
		pigz -p $2 $1".fastq"
		## Execute fastqc.
		echo -e "\nExecute fastqc..."
		fastqc $1".fastq.gz"
		## Save fastqc results in Fastqc folder.
		echo -e "\nMove fastqc files to Fastqc folder..."
		mv $1"_fastqc.zip" ./Fastqc
		mv $1"_fastqc.html" ./Fastqc
	else
		echo -e "\n$1 doesn't exist. Check it."
	fi
}

task_Trimming(){
	echo -e "\n"$1"..."
	
	if [ -e ../01-Raw_data/$1"_1.fastq.gz" ] && [ -e ../01-Raw_data/$1"_2.fastq.gz" ]; then
		## Trimming the libraries.
		echo -e "\nTrimming..."
		fastp \
			--adapter_fasta $2/TruSeq3-PE-2.fa \
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
			--adapter_fasta $2/TruSeq3-SE.fa \
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

task_Fastqc2(){
	echo -e "\n"$1"..."
	
	if [ -e $1"_tr_1P.fastq.gz" ] && [ -e $1"_tr_2P.fastq.gz" ]; then
		## Execute fastqc.
		echo -e "\nExecute fastqc..."
		fastqc $1"_tr_1P.fastq.gz"
		fastqc $1"_tr_2P.fastq.gz"
		## Save fastqc results in Fastqc folder.
		echo -e "\nMove fastqc files to Fastqc folder..."
		mv $1"_tr_1P_fastqc.zip" ./Fastqc
		mv $1"_tr_1P_fastqc.html" ./Fastqc
		mv $1"_tr_2P_fastqc.zip" ./Fastqc
		mv $1"_tr_2P_fastqc.html" ./Fastqc
	elif [ -e $1"_tr.fastq.gz" ]; then
		## Execute fastqc.
		echo -e "\nExecute fastqc..."
		fastqc $1"_tr.fastq.gz"
		## Save fastqc results in Fastqc folder.
		echo -e "\nMove fastqc files to Fastqc folder..."
		mv $1"_tr_fastqc.zip" ./Fastqc
		mv $1"_tr_fastqc.html" ./Fastqc
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
		if [ $num_sequences_1 -ge 750000 ] && [ $num_sequences_2 -ge 750000 ]; then
			echo -e $1" passes the trimming step."
			echo -e $1"\t$num_sequences_1/$num_sequences_2\tpasses the trimming step" >> $2/$1.tsv
			echo $1 >> $3/$1.txt
		elif [ $num_sequences_1 -ge 750000 ] || [ $num_sequences_2 -ge 750000 ]; then
			echo -e $1" passes the trimming step only forward or reverse)."
			echo -e $1"\t$num_sequences_1/$num_sequences_2\tpasses the trimming step only forward or reverse" >> $2/$1.tsv
		else
			echo -e $1" doesn't pass the trimming step."
			echo -e $1"\t$num_sequences_1/$num_sequences_2\tdoesn't pass the trimming step" >> $2/$1.tsv
		fi
	elif [ -e $1"_tr.fastq.gz" ]; then
		## Calculate the library depth.
		echo -e "\nCalculate the library depth..."
		total_lines=$(zcat $1"_tr.fastq.gz" | wc -l)
		num_sequences=$(($total_lines/4))
		## Check if the library is an empty file.
		echo -e "\nCheck if the library is an empty file..."
		if [ $num_sequences -ge 750000 ]; then
			echo -e $1" passes the trimming step."
			echo -e $1"\t$num_sequences\tpasses the trimming step" >> $2/$1.tsv
			echo $1 >> $3/$1.txt
		else
			echo -e $1" doesn't pass the trimming step."
			echo -e $1"\t$num_sequences\tdoesn't pass the trimming step" >> $2/$1.tsv
		fi
	else
		echo -e "\n$1 doesn't exist. Check it."
	fi
}

task_Pseudomapping(){
	echo -e "\n"$1"..."
	
	if [ -e $2/02-Trimmed_data/$1"_tr_1P.fastq.gz" ] && [ -e $2/02-Trimmed_data/$1"_tr_2P.fastq.gz" ]; then
		echo -e "\nPseudomapping (Libraries vs reference transcriptome)..."
		salmon quant \
		-i $2/03-Strand_detection/02-Index/$3 \
		-l A \
		-1 $2/02-Trimmed_data/$1"_tr_1P.fastq.gz" \
		-2 $2/02-Trimmed_data/$1"_tr_2P.fastq.gz" \
		-g $4/GFF3_genes/$3.gff3 \
		-p $5 \
		-o $1 \
		--gcBias
	elif [ -e $2/02-Trimmed_data/$1"_tr.fastq.gz" ]; then
		echo -e "\nPseudomapping (Library vs reference transcriptome)..."
		salmon quant \
		-i $2/03-Strand_detection/02-Index/$3 \
		-l A \
	    	-r $2/02-Trimmed_data/$1"_tr.fastq.gz" \
	    	-g $4/GFF3_genes/$3.gff3 \
	    	-p $5 \
	    	-o $1 \
	    	--gcBias
	else
		echo -e "\n$1 doesn't exist. Check it."
	fi
}

task_Select_SS_libraries(){
	
	if [ -e $2/02-Trimmed_data/$1"_tr_1P.fastq.gz" ] && [ -e $2/02-Trimmed_data/$1"_tr_2P.fastq.gz" ]; then
		echo -e "Move the strand-specific library $1 to 04-Selected_data..."
		cp $2/02-Trimmed_data/$1"_tr_1P.fastq.gz" ./
		cp $2/02-Trimmed_data/$1"_tr_2P.fastq.gz" ./
	elif [ -e $2/02-Trimmed_data/$1"_tr.fastq.gz" ]; then
		echo -e "Move the strand-specific library $1 to 04-Selected_data..."
		cp $2/02-Trimmed_data/$1"_tr.fastq.gz" ./
	else
		echo -e "\n$1 doesn't exist. Check it."
	fi
}

"$@"

