#!/bin/bash

####### FUNCTIONS

task_LncRNAs_prediction_STEP1(){
	
	####### VARIABLES
	P1=$2/STEP1
	G=$P1/Original_genes
	L=$P1/Potential_lncRNAs
	
	####### DIRECTORY
	mkdir -p $G
	mkdir -p $L
	
	####### PIPELINE
	
	### ORIGINAL GENES
	echo -e "\nORIGINAL GENES..."
	cd $G
	
	## Ids.
	echo -e "\nGet Ids..."
	gffread -w OG_temp.fasta -W -g $3/Genome/$4.fa $3/GTF_genes/$4.gtf
	sed '/^>/s/ .*//' OG_temp.fasta > OG_temp_mod.fasta 
	grep ">" OG_temp_mod.fasta | sed 's/>//g' > ORIGINAL_GENES_ids.txt
	rm OG_temp.fasta

	## Fasta.
	echo -e "\nGet Fasta..."
	seqtk subseq OG_temp_mod.fasta ORIGINAL_GENES_ids.txt > ORIGINAL_GENES.fasta
	rm OG_temp_mod.fasta

	## GTF.
	echo -e "\nGet GTF..."
	cp $3/GTF_genes/$4.gtf $G/ORIGINAL_GENES.gtf

	## TSV.
	echo -e "\nGet TSV..."
	GFF3orGTF2TABLE.py \
		--gtf $G/ORIGINAL_GENES.gtf \
		--fasta $G/ORIGINAL_GENES.fasta \
		--tab $G/ORIGINAL_GENES.tsv \
		--mode 'original' \
		--sep $'\t'


	### POTENTIAL LNCRNAS
	echo -e "\nPOTENTIAL LNCRNAS..."
	cd $L

	## Ids: Select the potential lncRNAs trascripts by class code and remove the unstranded transcripts.
	echo -e "\nGet Ids: Select the potential lncRNAs trascripts by class code and remove the unstranded transcripts..."
	GFF3orGTF2TABLE.py \
		--gtf $1/$4\_merged_compared_$5.annotated.gtf \
		--tab $L/POTENTIAL_LNCRNAS_temp.tsv \
		--mode 'assembly' \
		--sep $'\t'
	tail -n +2 POTENTIAL_LNCRNAS_temp.tsv | awk '$8 == "x" || $8 == "u" || $8 == "i" || $8 == "o" || $8 == "e" {print $0}' | awk '$5 != "." {print $6}' > POTENTIAL_LNCRNAS_ids.txt
	rm POTENTIAL_LNCRNAS_temp.tsv

	## GTF.
	echo -e "\nGet GTF..."
	Filter_GTF.py \
		--gtf-initial $1/$4\_merged_compared_$5.annotated.gtf \
		--gtf-final $L/POTENTIAL_LNCRNAS.gtf \
		--ids $L/POTENTIAL_LNCRNAS_ids.txt

	## Fasta.
	echo -e "\nGet Fasta..."
	gffread -w PL_temp.fasta -W -g $3/Genome/$4.fa POTENTIAL_LNCRNAS.gtf
	sed '/^>/s/ .*//' PL_temp.fasta > PL_temp_mod.fasta
	seqtk subseq PL_temp_mod.fasta POTENTIAL_LNCRNAS_ids.txt > POTENTIAL_LNCRNAS.fasta
	rm PL_temp.fasta PL_temp_mod.fasta

	## TSV.
	echo -e "\nGet TSV..."
	GFF3orGTF2TABLE.py \
		--gtf $L/POTENTIAL_LNCRNAS.gtf \
		--fasta $L/POTENTIAL_LNCRNAS.fasta \
		--tab $L/POTENTIAL_LNCRNAS.tsv \
		--mode 'assembly' \
		--sep $'\t'
}

task_LncRNAs_prediction_STEP2(){
	
	####### VARIABLES
	P1=$1/STEP1
	P2=$1/STEP2
	
	####### DIRECTORY
	mkdir -p $P2/CPC2
	mkdir -p $P2/FEELnc
	mkdir -p $P2/CPAT
	
	####### PIPELINE
	cd $P2
	
	########### Predict lncRNAs by CPC2.
	echo -e "\nPREDICTING LNCRNAS BY CPC2...\n"
	CPC2.py -i $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta -o ./CPC2/output_CPC2 --ORF
	awk '$9 == "noncoding" {print $1}' ./CPC2/output_CPC2.txt > ./CPC2/noncoding_ids_CPC2.txt

	########### Predict lncRNAs by FEELNC.
	echo -e "\n\n\nPREDICTING LNCRNAS BY FEELNC...\n"
	#--spethres=0.99,0.99 -> Esta opcion es mucho mas estricta.
	FEELnc_codpot.pl --outdir="./FEELnc" -o output_FEELnc -i $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta -a $P1/Original_genes/ORIGINAL_GENES.fasta --mode=shuffle --proc $3
	awk '$11 == 0 {print $1}' ./FEELnc/output_FEELnc_RF.txt > ./FEELnc/noncoding_ids_FEELnc.txt
	cat ./FEELnc/noncoding_ids_FEELnc.txt > ./FEELnc/noncoding_and_no_ORF_ids_FEELnc.txt
	if [ -f ./FEELnc/output_FEELnc.noORF.fa ]; then
	    grep ">" ./FEELnc/output_FEELnc.noORF.fa | sed 's/>//g' >> ./FEELnc/noncoding_and_no_ORF_ids_FEELnc.txt
	fi

	########### Predict the lncRNAs with CPAT.
	echo -e "\n\n\nPREDICTING LNCRNAS BY CPAT...\n"
	cpat.py -g $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta -o ./CPAT/output_cpat -x $2/cpat-3.0.2/crema_models/ath_hexamer.tsv -d $2/cpat-3.0.2/crema_models/ath_logit.RData
	awk '$11 < 0.5 {print $1}' ./CPAT/output_cpat.ORF_prob.best.tsv > ./CPAT/noncoding_ids_CPAT.txt
	cat ./CPAT/noncoding_ids_CPAT.txt > ./CPAT/noncoding_and_no_ORF_ids_CPAT.txt
	if [ -f ./CPAT/output_cpat.no_ORF.txt ]; then
	    cat ./CPAT/output_cpat.no_ORF.txt >> ./CPAT/noncoding_and_no_ORF_ids_CPAT.txt
	fi
	rm CPAT_run_info.log
}

task_LncRNAs_prediction_STEP3(){
	
	####### VARIABLES
	P1=$1/STEP1
	P3=$1/STEP3
	
	####### DIRECTORY
	mkdir -p $P3/SwissProt
	mkdir -p $P3/Pfam
	
	####### PIPELINE
	cd $P3
	
	########### Compare the potential lncRNAs with SwissProt database
	echo -e "\nFINDING CODING SEQUENCES IN SWISSPROT..."
	
	## Make the Swissprot database.
	#http://gensoft.pasteur.fr/docs/diamond/0.9.19/diamond_manual.pdf
	echo -e "\nMaking the SwissProt database...\n"
	diamond makedb --in $2/Swissprot/uniprot_sprot.fasta -d ./SwissProt/swissprot.dmnd
	## Execute blastx against Swissprot database.
	echo -e "\nExecuting blastx...\n"
	diamond blastx \
		-d ./SwissProt/swissprot.dmnd \
		-q $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta \
		-o ./SwissProt/diamond_output_temp.tsv \
		--un ./SwissProt/diamond_output_unaligned.fasta \
		--al ./SwissProt/diamond_output_aligned.fasta \
		--unfmt fasta \
		--alfmt fasta \
		--strand plus \
		--matrix BLOSUM62 \
		--gapopen 11 \
		--gapextend 1 \
		--more-sensitive \
		--top 5 \
		--evalue $3 \
		--threads $4 \
		-f 6 qseqid qlen sseqid slen qstart qend sstart send qframe pident gaps length evalue bitscore

	## Add info of alignment and the table header.
	echo -e "\nAdding info and the table header...\n"
	Add_info_to_Swissprot_blastx_tsv.py \
		--tsv-initial $P3/SwissProt/diamond_output_temp.tsv \
		--fasta $2/Swissprot/uniprot_sprot.fasta \
		--tsv-final $P3/SwissProt/diamond_output.tsv
	rm ./SwissProt/diamond_output_temp.tsv


	########### Compare the potential lncRNAs with Pfam database
	echo -e "\n\nFINDING CODING SEQUENCES IN PFAM-A..."

	## Trandecoder
	#https://github.com/TransDecoder/TransDecoder/wiki
	echo -e "\nFinding ORFs...\n"
	TransDecoder.LongOrfs -t $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta -m 20 -S --output_dir ./Pfam/Transdecoder
	rm *.cmds
	rm -r ./Pfam/*.__checkpoints_longorfs

	## hmmscan or hmmsearch
	#http://eddylab.org/software/Hmmer3/3.1b2/Userguide.pdf
	# hmmsearch is more faster than hmmscan. However hmmscan generates a table with the desirable characteristics.

	mkdir -p $P3/Pfam/Hmmer

	echo -e "\nCopying the PFAM-A database..."
	cp $2/PFAM/Pfam-A.hmm ./Pfam/Hmmer/

	echo -e "\nExecuting hmmsearch...\n"
	hmmsearch \
		--tblout ./Pfam/Hmmer/Results_PFAM_tblout_temp.tsv \
		--domtblout ./Pfam/Hmmer/Results_PFAM_domtblout_temp.tsv \
		--cpu $4 \
		-E $3 \
		--domE $3 \
		./Pfam/Hmmer/Pfam-A.hmm \
		./Pfam/Transdecoder/longest_orfs.pep 1> ./Pfam/Hmmer/stdout.txt

	## Convert to tabular file.
	echo -e "\nConverting result tables to tabular tables..."
	# FS = Full Sequence
	# D = Domain
	echo -e "target_name\ttaccession\ttlen\tquery_name\tqaccession\tqlen\te-value.FS\tscore.FS\tbias.FS\tn.D\ttotal.D\tc-evalue.D\ti-evalue.D\tscore.D\tbias.D\tfrom.hmm\tto.hmm\tfrom.ali\tto.ali\tfrom.env\tto.env\tacc\ttarget_description" > ./Pfam/Hmmer/Results_PFAM_domtblout.tsv
	tail -n +4 ./Pfam/Hmmer/Results_PFAM_domtblout_temp.tsv | head -n -10 | awk -F" " '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23" "$24" "$25}' >> ./Pfam/Hmmer/Results_PFAM_domtblout.tsv
	rm ./Pfam/Hmmer/Results_PFAM_domtblout_temp.tsv
	# FS = Full Sequence
	# BD = Best Domain
	# DNE = Domain Number Estimation
	echo -e "target_name\ttaccession\tquery_name\tqaccession\te-value.FS\tscore.FS\tbias.FS\te-value.BD\tscore.BD\tbias.BD\t\texp.DNE\treg.DNE\tclu.DNE\tov.DNE\tenv.DNE\tdom.DNE\trep.DNE\tinc.DNE\ttarget_description" > ./Pfam/Hmmer/Results_PFAM_tblout.tsv
	tail -n +4 ./Pfam/Hmmer/Results_PFAM_tblout_temp.tsv | head -n -10 | awk -F" " '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19" "$20" "$21}' >> ./Pfam/Hmmer/Results_PFAM_tblout.tsv
	rm ./Pfam/Hmmer/Results_PFAM_tblout_temp.tsv
}

task_LncRNAs_prediction_STEP4(){
	
	####### VARIABLES
	P1=$1/STEP1
	P4=$1/STEP4
	
	####### DIRECTORY
	mkdir -p $P4/ORF
	
	####### PIPELINE
	cd $P4

	########### Calculate the longest ORFs for each potential lncRNA
	#https://github.com/TransDecoder/TransDecoder/wiki
	
	## Limit 80 aa.
	echo -e "\nFINDING ORFS > 80 AA...\n"
	TransDecoder.LongOrfs -t $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta -m 80 -S --output_dir ./ORF/Transdecoder_80
	rm *.cmds
	rm -r ./ORF/*.__checkpoints_longorfs
	sed '/^>/s/ .*//' ./ORF/Transdecoder_80/longest_orfs.pep | grep ">" | sed 's/>//g' | sed 's/\.p[0-9]\+//g' | sort | uniq > ./ORF/Transdecoder_80/ids_more_than_80_aa.txt
	cat $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS_ids.txt ./ORF/Transdecoder_80/ids_more_than_80_aa.txt | sort | uniq -u > ./ORF/Transdecoder_80/noncoding_ids_ORF_length.txt

	## Limit 100 aa.
	echo -e "\nFINDING ORFS > 100 AA...\n"
	TransDecoder.LongOrfs -t $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta -m 100 -S --output_dir ./ORF/Transdecoder_100
	rm *.cmds
	rm -r ./ORF/*.__checkpoints_longorfs
	sed '/^>/s/ .*//' ./ORF/Transdecoder_100/longest_orfs.pep | grep ">" | sed 's/>//g' | sed 's/\.p[0-9]\+//g' | sort | uniq > ./ORF/Transdecoder_100/ids_more_than_100_aa.txt
	cat $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS_ids.txt ./ORF/Transdecoder_100/ids_more_than_100_aa.txt | sort | uniq -u > ./ORF/Transdecoder_100/noncoding_ids_ORF_length.txt

	## Limit 120 aa.
	echo -e "\nFINDING ORFS > 120 AA...\n"
	TransDecoder.LongOrfs -t $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta -m 120 -S --output_dir ./ORF/Transdecoder_120
	rm *.cmds
	rm -r ./ORF/*.__checkpoints_longorfs
	sed '/^>/s/ .*//' ./ORF/Transdecoder_120/longest_orfs.pep | grep ">" | sed 's/>//g' | sed 's/\.p[0-9]\+//g' | sort | uniq > ./ORF/Transdecoder_120/ids_more_than_120_aa.txt
	cat $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS_ids.txt ./ORF/Transdecoder_120/ids_more_than_120_aa.txt | sort | uniq -u > ./ORF/Transdecoder_120/noncoding_ids_ORF_length.txt
}

task_LncRNAs_prediction_STEP5(){
	
	####### VARIABLES
	P1=$1/STEP1
	P5=$1/STEP5
	
	####### DIRECTORY
	mkdir -p $P5/RNAcentral
	mkdir -p $P5/PmiREN
	mkdir -p $P5/miRBase
	
	####### PIPELINE
	cd $P5
	
	########### Annotate the sequences with the RNAcentral database (rRNA, tRNA, snRNA, snoRNA).
	#https://academic.oup.com/bioinformatics/article/29/19/2487/186765
	#https://rnacentral.org/help/sequence-search
	echo -e "\nCLASSIFYING POTENTIAL LNCRNAS WITH RNACENTRAL DATABASE (rRNA, tRNA, snRNA, snoRNA)...\n"
	
	## Make the RNAcentral database (rRNA, tRNA, snRNA, snoRNA).
	echo -e "Making the Blast database..."
	makeblastdb \
		-in $2/RNAcentral/ALL/cucurbitaceae_RNAcentral.fasta \
		-dbtype nucl \
		-parse_seqids \
		-out ./RNAcentral/RNAcentral \
		-title RNAcentral
	## Execute blastn against RNAcentral database (rRNA, tRNA, snRNA, snoRNA).
	echo -e "Executing the blastn..."
	blastn \
		-query $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta \
		-task blastn \
		-db ./RNAcentral/RNAcentral \
		-out ./RNAcentral/output_blastn_temp.tsv \
		-evalue $3\
		-num_threads $4 \
		-strand plus \
		-outfmt "6 delim= qseqid qlen sseqid slen qstart qend sstart send qframe sframe sstrand pident gaps length evalue bitscore"
	## Add info of alignment and the table header.
	echo -e "\nAdding info and the table header..."
	Add_info_to_RNAcentral_blastn_tsv.py \
		--tsv-initial $P5/RNAcentral/output_blastn_temp.tsv \
		--fasta $2/RNAcentral/ALL/cucurbitaceae_RNAcentral.fasta \
		--tsv-final $P5/RNAcentral/output_blastn.tsv
	rm ./RNAcentral/output_blastn_temp.tsv
	grep -v "bacterial" ./RNAcentral/output_blastn.tsv > ./RNAcentral/output_blastn_no_bacterial.tsv


	########### Annotate the sequences with the PmiREN database (miRNA-precursors and mature miRNAs).

	echo -e "\n\n\nSEARCHING FOR MIRNA-PRECURSORS USING PMIREN DATABASE...\n"
	
	## Make the miRNA precursors database (miRNA-precursors).
	echo -e "Making the Blast database..."
	makeblastdb \
		-in $2/miRNAs/PmiREN/ALL/cucurbitaceae_miRNA_precursors.fa \
		-dbtype nucl \
		-parse_seqids \
		-out ./PmiREN/Precursor_miRNAs \
		-title Precursor_miRNAs
	## Execute blastn against miRNA precursors database (pre-miRNAs).
	echo -e "Executing the blastn..."
	blastn \
		-query $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta \
		-task blastn \
		-db ./PmiREN/Precursor_miRNAs \
		-out ./PmiREN/output_blastn_prec_temp.tsv \
		-evalue $3 \
		-num_threads $4 \
		-strand plus \
		-outfmt "6 delim= qseqid qlen sseqid slen qstart qend sstart send qframe sframe sstrand pident gaps length evalue bitscore"
	## Add the table header. sseqid contains the info.
	echo -e "\nAdding the table header..."
	echo -e "qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tqframe\tsframe\tsstrand\tpident\tgaps\tlength\tevalue\tbitscore" > ./PmiREN/output_blastn_prec.tsv
	cat ./PmiREN/output_blastn_prec_temp.tsv >> ./PmiREN/output_blastn_prec.tsv
	rm ./PmiREN/output_blastn_prec_temp.tsv
	# Ids.
	echo -e "\nGet Ids..."
	tail -n +2 ./PmiREN/output_blastn_prec.tsv | awk '{print $1}' | sort | uniq > ./PmiREN/potential_prec_ID.txt
	# Fasta.
	echo -e "\nGet Fasta..."
	seqtk subseq $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta ./PmiREN/potential_prec_ID.txt > ./PmiREN/potential_prec.fasta

	echo -e "\n\n\nSEARCHING FOR MATURE MIRNAS USING PMIREN DATABASE...\n"

	## Make the mature miRNA database (mature miRNAs).
	echo -e "Making the Blast database..."
	makeblastdb \
		-in $2/miRNAs/PmiREN/ALL/cucurbitaceae_miRNA_mature.fa \
		-dbtype nucl \
		-parse_seqids \
		-out ./PmiREN/Mature_miRNAs \
		-title Mature_miRNAs
	## Execute blastn against mature miRNA database (mature miRNAs).
	echo -e "Executing the blastn..."
	blastn \
		-query ./PmiREN/potential_prec.fasta \
		-task blastn \
		-db ./PmiREN/Mature_miRNAs \
		-out ./PmiREN/output_blastn_mat_temp.tsv \
		-dust no \
		-evalue 0.05 \
		-word_size 7 \
		-num_threads $4 \
		-strand plus \
		-outfmt "6 delim= qseqid qlen sseqid slen qstart qend sstart send qframe sframe sstrand pident gaps mismatch length evalue bitscore"
	## Add the table header. sseqid contains the info.
	echo -e "\nAdding the table header..."
	echo -e "qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tqframe\tsframe\tsstrand\tpident\tgaps\tmismatch\tlength\tevalue\tbitscore" > ./PmiREN/output_blastn_mat.tsv
	cat ./PmiREN/output_blastn_mat_temp.tsv >> ./PmiREN/output_blastn_mat.tsv
	rm ./PmiREN/output_blastn_mat_temp.tsv

	echo -e "\n\n\nCHECKING THE POTENTIAL MIRNA PRECURSORS WITH MIRENA...\n"

	echo -e "Executing MIReNA with miRNA-precursor length = 100 nucleotides..."

	Create_fasta_for_MIReNA.py \
		--blastn-mat $P5/PmiREN/output_blastn_mat.tsv \
		--blastn-prec $P5/PmiREN/output_blastn_prec.tsv \
		--fasta-initial $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta \
		--fasta-MIReNA $P5/PmiREN/MIReNA-100.fasta \
		--length 100 1> ./PmiREN/Create_fasta_for_MIReNA-100.txt
		
	if [ -f "$P5/PmiREN/MIReNA-100.txt" ]; then
	    rm $P5/PmiREN/MIReNA-100.txt
	fi
	source activate ViennaRNA2.5.1
	MIReNA.sh -v -x -f $P5/PmiREN/MIReNA-100.fasta -o $P5/PmiREN/MIReNA-100.txt
	conda deactivate
	grep ">" $P5/PmiREN/MIReNA-100.txt | sed '/^>/ s/ .*//' | sed 's/>//g' | sort | uniq > ./PmiREN/prec_ID-100.txt
	seqtk subseq $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta ./PmiREN/prec_ID-100.txt > ./PmiREN/prec-100.fasta

	echo -e "\nExecuting MIReNA with miRNA-precursor length = 200 nucleotides..."

	Create_fasta_for_MIReNA.py \
		--blastn-mat $P5/PmiREN/output_blastn_mat.tsv \
		--blastn-prec $P5/PmiREN/output_blastn_prec.tsv \
		--fasta-initial $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta \
		--fasta-MIReNA $P5/PmiREN/MIReNA-200.fasta \
		--length 200 \
		--remove-IDs $P5/PmiREN/prec-100.fasta 1> ./PmiREN/Create_fasta_for_MIReNA-200.txt
		
	if [ -f "$P5/PmiREN/MIReNA-200.txt" ]; then
	    rm $P5/PmiREN/MIReNA-200.txt
	fi
	source activate ViennaRNA2.5.1
	MIReNA.sh -v -x -f $P5/PmiREN/MIReNA-200.fasta -o $P5/PmiREN/MIReNA-200.txt
	conda deactivate
	grep ">" $P5/PmiREN/MIReNA-200.txt | sed '/^>/ s/ .*//' | sed 's/>//g' | sort | uniq > ./PmiREN/prec_ID-200.txt
	seqtk subseq $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta ./PmiREN/prec_ID-200.txt > ./PmiREN/prec-200.fasta
	cat ./PmiREN/prec-100.fasta ./PmiREN/prec-200.fasta | sort | uniq > ./PmiREN/prec_ID-100-200.txt
		
	echo -e "\nExecuting MIReNA with miRNA-precursor length = 300 nucleotides..."

	Create_fasta_for_MIReNA.py \
		--blastn-mat $P5/PmiREN/output_blastn_mat.tsv \
		--blastn-prec $P5/PmiREN/output_blastn_prec.tsv \
		--fasta-initial $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta \
		--fasta-MIReNA $P5/PmiREN/MIReNA-300.fasta \
		--length 300 \
		--remove-IDs $P5/PmiREN/prec_ID-100-200.txt 1> ./PmiREN/Create_fasta_for_MIReNA-300.txt
		
	rm ./PmiREN/prec_ID-100-200.txt
		
	if [ -f "$P5/PmiREN/MIReNA-300.txt" ]; then
	    rm $P5/PmiREN/MIReNA-300.txt
	fi
	source activate ViennaRNA2.5.1
	MIReNA.sh -v -x -f $P5/PmiREN/MIReNA-300.fasta -o $P5/PmiREN/MIReNA-300.txt
	conda deactivate
	grep ">" $P5/PmiREN/MIReNA-300.txt | sed '/^>/ s/ .*//' | sed 's/>//g' | sort | uniq > ./PmiREN/prec_ID-300.txt
	seqtk subseq $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta ./PmiREN/prec_ID-300.txt > ./PmiREN/prec-300.fasta


	cat ./PmiREN/prec_ID-100.txt ./PmiREN/prec_ID-200.txt ./PmiREN/prec_ID-300.txt | sort | uniq > ./PmiREN/prec_ID.txt
	seqtk subseq $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta ./PmiREN/prec_ID.txt > ./PmiREN/prec.fasta


	########### Annotate the sequences with the miRBase database (miRNA-precursors and mature miRNAs).

	echo -e "\n\n\nSEARCHING FOR MIRNA-PRECURSORS USING MIRBASE DATABASE...\n"

	## Make the miRNA precursors database (miRNA-precursors).
	echo -e "Making the Blast database..."
	makeblastdb \
		-in $2/miRNAs/miRBase/ALL/cucurbitaceae_miRNA_precursors.fa \
		-dbtype nucl \
		-parse_seqids \
		-out ./miRBase/Precursor_miRNAs \
		-title Precursor_miRNAs
	## Execute blastn against miRNA precursors database (pre-miRNAs).
	echo -e "Executing the blastn..."
	blastn \
		-query $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta \
		-task blastn \
		-db ./miRBase/Precursor_miRNAs \
		-out ./miRBase/output_blastn_prec_temp.tsv \
		-evalue $3 \
		-num_threads $4 \
		-strand plus \
		-outfmt "6 delim= qseqid qlen sseqid slen qstart qend sstart send qframe sframe sstrand pident gaps length evalue bitscore"
	## Add the table header. sseqid contains the info.
	echo -e "\nAdding the table header..."
	echo -e "qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tqframe\tsframe\tsstrand\tpident\tgaps\tlength\tevalue\tbitscore" > ./miRBase/output_blastn_prec.tsv
	cat ./miRBase/output_blastn_prec_temp.tsv >> ./miRBase/output_blastn_prec.tsv
	rm ./miRBase/output_blastn_prec_temp.tsv
	# Ids.
	echo -e "\nGet Ids..."
	tail -n +2 ./miRBase/output_blastn_prec.tsv | awk '{print $1}' | sort | uniq > ./miRBase/potential_prec_ID.txt
	# Fasta.
	echo -e "\nGet Fasta..."
	seqtk subseq $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta ./miRBase/potential_prec_ID.txt > ./miRBase/potential_prec.fasta

	echo -e "\n\n\nSEARCHING FOR MATURE MIRNAS USING MIRBASE DATABASE...\n"

	## Make the mature miRNA database (mature miRNAs).
	echo -e "Making the Blast database..."
	makeblastdb \
		-in $2/miRNAs/miRBase/ALL/cucurbitaceae_miRNA_mature.fa \
		-dbtype nucl \
		-parse_seqids \
		-out ./miRBase/Mature_miRNAs \
		-title Mature_miRNAs
	## Execute blastn against mature miRNA database (mature miRNAs).
	echo -e "Executing the blastn..."
	blastn \
		-query ./miRBase/potential_prec.fasta \
		-task blastn \
		-db ./miRBase/Mature_miRNAs \
		-out ./miRBase/output_blastn_mat_temp.tsv \
		-dust no \
		-evalue 0.05 \
		-word_size 7 \
		-num_threads $4 \
		-strand plus \
		-outfmt "6 delim= qseqid qlen sseqid slen qstart qend sstart send qframe sframe sstrand pident gaps mismatch length evalue bitscore"
	## Add the table header. sseqid contains the info.
	echo -e "\nAdding the table header..."
	echo -e "qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tqframe\tsframe\tsstrand\tpident\tgaps\tmismatch\tlength\tevalue\tbitscore" > ./miRBase/output_blastn_mat.tsv
	cat ./miRBase/output_blastn_mat_temp.tsv >> ./miRBase/output_blastn_mat.tsv
	rm ./miRBase/output_blastn_mat_temp.tsv

	echo -e "\n\n\nCHECKING THE POTENTIAL MIRNA PRECURSORS WITH MIRENA...\n"

	echo -e "Executing MIReNA with miRNA-precursor length = 100 nucleotides..."

	Create_fasta_for_MIReNA.py \
		--blastn-mat $P5/miRBase/output_blastn_mat.tsv \
		--blastn-prec $P5/miRBase/output_blastn_prec.tsv \
		--fasta-initial $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta \
		--fasta-MIReNA $P5/miRBase/MIReNA-100.fasta \
		--length 100 1> ./miRBase/Create_fasta_for_MIReNA-100.txt
		
	if [ -f "$P5/miRBase/MIReNA-100.txt" ]; then
	    rm $P5/miRBase/MIReNA-100.txt
	fi
	source activate ViennaRNA2.5.1
	MIReNA.sh -v -x -f $P5/miRBase/MIReNA-100.fasta -o $P5/miRBase/MIReNA-100.txt
	conda deactivate
	grep ">" $P5/miRBase/MIReNA-100.txt | sed '/^>/ s/ .*//' | sed 's/>//g' | sort | uniq > ./miRBase/prec_ID-100.txt
	seqtk subseq $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta ./miRBase/prec_ID-100.txt > ./miRBase/prec-100.fasta

	echo -e "\nExecuting MIReNA with miRNA-precursor length = 200 nucleotides..."

	Create_fasta_for_MIReNA.py \
		--blastn-mat $P5/miRBase/output_blastn_mat.tsv \
		--blastn-prec $P5/miRBase/output_blastn_prec.tsv \
		--fasta-initial $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta \
		--fasta-MIReNA $P5/miRBase/MIReNA-200.fasta \
		--length 200 \
		--remove-IDs $P5/miRBase/prec-100.fasta 1> ./miRBase/Create_fasta_for_MIReNA-200.txt
		
	if [ -f "$P5/miRBase/MIReNA-200.txt" ]; then
	    rm $P5/miRBase/MIReNA-200.txt
	fi
	source activate ViennaRNA2.5.1
	MIReNA.sh -v -x -f $P5/miRBase/MIReNA-200.fasta -o $P5/miRBase/MIReNA-200.txt
	conda deactivate
	grep ">" $P5/miRBase/MIReNA-200.txt | sed '/^>/ s/ .*//' | sed 's/>//g' | sort | uniq > ./miRBase/prec_ID-200.txt
	seqtk subseq $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta ./miRBase/prec_ID-200.txt > ./miRBase/prec-200.fasta
	cat ./miRBase/prec-100.fasta ./miRBase/prec-200.fasta | sort | uniq > ./miRBase/prec_ID-100-200.txt

	echo -e "\nExecuting MIReNA with miRNA-precursor length = 300 nucleotides..."
		
	Create_fasta_for_MIReNA.py \
		--blastn-mat $P5/miRBase/output_blastn_mat.tsv \
		--blastn-prec $P5/miRBase/output_blastn_prec.tsv \
		--fasta-initial $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta \
		--fasta-MIReNA $P5/miRBase/MIReNA-300.fasta \
		--length 300 \
		--remove-IDs $P5/miRBase/prec_ID-100-200.txt 1> ./miRBase/Create_fasta_for_MIReNA-300.txt
		
	rm ./miRBase/prec_ID-100-200.txt
		
	if [ -f "$P5/miRBase/MIReNA-300.txt" ]; then
	    rm $P5/miRBase/MIReNA-300.txt
	fi
	source activate ViennaRNA2.5.1
	MIReNA.sh -v -x -f $P5/miRBase/MIReNA-300.fasta -o $P5/miRBase/MIReNA-300.txt
	conda deactivate
	grep ">" $P5/miRBase/MIReNA-300.txt | sed '/^>/ s/ .*//' | sed 's/>//g' | sort | uniq > ./miRBase/prec_ID-300.txt
	seqtk subseq $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta ./miRBase/prec_ID-300.txt > ./miRBase/prec-300.fasta


	cat ./miRBase/prec_ID-100.txt ./miRBase/prec_ID-200.txt ./miRBase/prec_ID-300.txt | sort | uniq > ./miRBase/prec_ID.txt
	seqtk subseq $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta ./miRBase/prec_ID.txt > ./miRBase/prec.fasta
}

task_LncRNAs_prediction_STEP6(){
	
	####### VARIABLES
	P1=$1/STEP1
	P6=$1/STEP6
	
	####### DIRECTORY
	mkdir -p $P6/CANTATAdb
	mkdir -p $P6/PLncDB
	mkdir -p $P6/GreeNC
	
	####### PIPELINE
	cd $P6
	
	########### Annotate the sequences with the known lncRNA of lncRNA databases.
	#http://nebc.nerc.ac.uk/bioinformatics/documentation/blast+/user_manual.pdf

	###### CANTATADB
	
	echo -e "\nANNOTATE THE SEQUENCES WITH KNOWN CUCURBITACEAE LNCRNAS STORED IN CANTATADB...\n"
	#http://cantata.amu.edu.pl/download.php

	## Make the LncRNAs database.
	echo -e "Making the Blast database..."
	makeblastdb \
		-in $2/LncRNAs_databases/CANTATAdb/ALL/cucurbitaceae_CANTATAdb.fa \
		-dbtype nucl \
		-parse_seqids \
		-out ./CANTATAdb/database \
		-title LncRNAs
	## Execute blastn against LncRNAs database.
	echo -e "Executing the blastn..."
	blastn \
		-query $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta \
		-task blastn \
		-db ./CANTATAdb/database \
		-out ./CANTATAdb/output_blastn_temp.tsv \
		-evalue $3 \
		-num_threads $4 \
		-strand plus \
		-outfmt "6 delim= qseqid qlen sseqid slen qstart qend sstart send qframe sframe sstrand pident gaps length evalue bitscore"
	## Add the table header. sseqid contains the info.
	echo -e "\nAdding the table header..."
	echo -e "qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tqframe\tsframe\tsstrand\tpident\tgaps\tlength\tevalue\tbitscore" > ./CANTATAdb/output_blastn.tsv
	cat ./CANTATAdb/output_blastn_temp.tsv >> ./CANTATAdb/output_blastn.tsv
	rm ./CANTATAdb/output_blastn_temp.tsv
	## Filter by coverage.
	echo -e "\nFiltering by coverage..."
	tail -n +2 ./CANTATAdb/output_blastn.tsv | awk '{$16 = ( ( ( $6 - $5 ) + 1 ) * 100 ) / $2}1' | column -t | awk '{$17 = ( ( ( $8 - $7 ) + 1 ) * 100 ) / $4}1' | column -t | awk '$16 >= 50 || $17 >= 50 {print $0}' > ./CANTATAdb/output_blastn_cov_temp.tsv
	echo -e "qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tqframe\tsframe\tsstrand\tpident\tgaps\tlength\tevalue\tbitscore\tqcov\tscov" > ./CANTATAdb/output_blastn_cov.tsv
	cat ./CANTATAdb/output_blastn_cov_temp.tsv >> ./CANTATAdb/output_blastn_cov.tsv
	rm ./CANTATAdb/output_blastn_cov_temp.tsv


	###### PLNCDB
	
	echo -e "\n\n\nANNOTATE THE SEQUENCES WITH KNOWN CUCURBITACEAE LNCRNAS STORED IN PLNCDB...\n"
	#https://www.tobaccodb.org/plncdb/searchList
	
	## Make the LncRNAs database.
	echo -e "Making the Blast database..."
	makeblastdb \
		-in $2/LncRNAs_databases/PLncDB/ALL/cucurbitaceae_PLncDB.fa \
		-dbtype nucl \
		-parse_seqids \
		-out ./PLncDB/database \
		-title LncRNAs
	## Execute blastn against LncRNAs database.
	echo -e "Executing the blastn..."
	blastn \
		-query $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta \
		-task blastn \
		-db ./PLncDB/database \
		-out ./PLncDB/output_blastn_temp.tsv \
		-evalue $3 \
		-num_threads $4 \
		-strand plus \
		-outfmt "6 delim= qseqid qlen sseqid slen qstart qend sstart send qframe sframe sstrand pident gaps length evalue bitscore"
	## Add the table header. sseqid contains the info.
	echo -e "\nAdding the table header..."
	echo -e "qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tqframe\tsframe\tsstrand\tpident\tgaps\tlength\tevalue\tbitscore" > ./PLncDB/output_blastn.tsv
	cat ./PLncDB/output_blastn_temp.tsv >> ./PLncDB/output_blastn.tsv
	rm ./PLncDB/output_blastn_temp.tsv
	## Filter by coverage.
	echo -e "\nFiltering by coverage..."
	tail -n +2 ./PLncDB/output_blastn.tsv | awk '{$16 = ( ( ( $6 - $5 ) + 1 ) * 100 ) / $2}1' | column -t | awk '{$17 = ( ( ( $8 - $7 ) + 1 ) * 100 ) / $4}1' | column -t | awk '$16 >= 50 || $17 >= 50 {print $0}' > ./PLncDB/output_blastn_cov_temp.tsv
	echo -e "qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tqframe\tsframe\tsstrand\tpident\tgaps\tlength\tevalue\tbitscore\tqcov\tscov" > ./PLncDB/output_blastn_cov.tsv
	cat ./PLncDB/output_blastn_cov_temp.tsv >> ./PLncDB/output_blastn_cov.tsv
	rm ./PLncDB/output_blastn_cov_temp.tsv

	###### GREENC

	echo -e "\n\n\nANNOTATE THE SEQUENCES WITH KNOWN CUCURBITACEAE LNCRNAS STORED IN GREENC...\n"
	#http://greenc.sequentiabiotech.com/wiki2/Species:Cucumis_melo_(Ensembl_Plants_51)
	
	## Make the LncRNAs database.
	echo -e "Making the Blast database..."
	makeblastdb \
		-in $2/LncRNAs_databases/GreeNC/ALL/cucurbitaceae_GreeNC.fa \
		-dbtype nucl \
		-parse_seqids \
		-out ./GreeNC/database \
		-title LncRNAs
	## Execute blastn against LncRNAs database.
	echo -e "Executing the blastn..."
	blastn \
		-query $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta \
		-task blastn \
		-db ./GreeNC/database \
		-out ./GreeNC/output_blastn_temp.tsv \
		-evalue $3 \
		-num_threads $4 \
		-strand plus \
		-outfmt "6 delim= qseqid qlen sseqid slen qstart qend sstart send qframe sframe sstrand pident gaps length evalue bitscore"
	## Add the table header. sseqid contains the info.
	echo -e "\nAdding the table header..."
	echo -e "qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tqframe\tsframe\tsstrand\tpident\tgaps\tlength\tevalue\tbitscore" > ./GreeNC/output_blastn.tsv
	cat ./GreeNC/output_blastn_temp.tsv >> ./GreeNC/output_blastn.tsv
	rm ./GreeNC/output_blastn_temp.tsv
	## Filter by coverage.
	echo -e "\nFiltering by coverage..."
	tail -n +2 ./GreeNC/output_blastn.tsv | awk '{$16 = ( ( ( $6 - $5 ) + 1 ) * 100 ) / $2}1' | column -t | awk '{$17 = ( ( ( $8 - $7 ) + 1 ) * 100 ) / $4}1' | column -t | awk '$16 >= 50 || $17 >= 50 {print $0}' > ./GreeNC/output_blastn_cov_temp.tsv
	echo -e "qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tqframe\tsframe\tsstrand\tpident\tgaps\tlength\tevalue\tbitscore\tqcov\tscov" > ./GreeNC/output_blastn_cov.tsv
	cat ./GreeNC/output_blastn_cov_temp.tsv >> ./GreeNC/output_blastn_cov.tsv
	rm ./GreeNC/output_blastn_cov_temp.tsv
}

task_LncRNAs_prediction_STEP-FINAL(){
	
	####### VARIABLES
	P1=$1/STEP1
	P7=$1/STEP-FINAL
	
	####### DIRECTORY
	mkdir -p $P7
	mkdir -p $P7/Database
	mkdir -p $P7/Database/Extra
	mkdir -p $P7/Figures
	mkdir -p $P7/Files
	mkdir -p $P7/Files/Genes
	mkdir -p $P7/Files/LncRNAs
	mkdir -p $P7/Files/LncRNAs/r
	mkdir -p $P7/Files/LncRNAs/nr
	mkdir -p $P7/Files/Joined
	mkdir -p $P7/Files/Joined/ALL
	mkdir -p $P7/Files/Joined/ALL/r
	mkdir -p $P7/Files/Joined/ALL/nr
	mkdir -p $P7/Files/Joined/FILTERED
	mkdir -p $P7/Files/Joined/FILTERED/r
	mkdir -p $P7/Files/Joined/FILTERED/nr
	mkdir -p $P7/Redundancy_analysis
	mkdir -p $P7/Redundancy_analysis/AGAT
	mkdir -p $P7/Redundancy_analysis/CGAT
	mkdir -p $P7/Redundancy_analysis/CDHIT

	####### PIPELINE
	cd $P7

	DB=$P7/Database
	FG=$P7/Figures
	FL=$P7/Files
	NRA=$P7/Redundancy_analysis

	#########################################################################

	####### CREATE THE LNCRNAS DATABASE.

	echo -e "\nCREATE THE LNCRNAS DATABASE...\n"
	Rscript $3/Create_customed_LncRNA_db.R $5 $1 $2 $4

	#########################################################################

	####### REDUNDANT: BUILD THE FOLDERS FOR DOWNSTREAM ANALYSIS.

	echo -e "\n\n\nREDUNDANT: BUILD THE FOLDERS FOR DOWNSTREAM ANALYSIS..."

	######################################

	#### GENES
	echo -e "\n\nGENES..."
	cd $FL/Genes

	## Ids, fasta, gtf and tsv
	echo -e "\nCopy from STEP1/Original_genes ids file..."
	echo -e "\nCopy from STEP1/Original_genes fasta file..."
	echo -e "\nCopy from STEP1/Original_genes gtf file..."
	echo -e "\nCopy from STEP1/Original_genes tsv file..."
	cp $P1/Original_genes/* ./

	######################################

	#### LNCRNAS
	echo -e "\n\nLNCRNAS..."
	cd $FL/LncRNAs/r

	## ids
	echo -e "\nCreate new ids file..."
	tail -n +2 $DB/Database_LncRNAs.tsv | awk '{print $1}' > POTENTIAL_LNCRNAS_ids_pred.txt
	## fasta
	echo -e "\nCreate new fasta file..."
	seqtk subseq $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta POTENTIAL_LNCRNAS_ids_pred.txt > POTENTIAL_LNCRNAS_pred.fasta
	## gtf
	echo -e "\nCreate new gtf file..."
	Filter_GTF.py \
		--gtf-initial $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.gtf \
		--gtf-final $FL/LncRNAs/r/POTENTIAL_LNCRNAS_pred.gtf \
		--ids $FL/LncRNAs/r/POTENTIAL_LNCRNAS_ids_pred.txt
	## tsv
	echo -e "\nCreate new tsv file..."
	GFF3orGTF2TABLE.py \
		--gtf $FL/LncRNAs/r/POTENTIAL_LNCRNAS_pred.gtf \
		--fasta $FL/LncRNAs/r/POTENTIAL_LNCRNAS_pred.fasta \
		--tab $FL/LncRNAs/r/POTENTIAL_LNCRNAS_pred.tsv \
		--mode 'assembly' \
		--sep $'\t'
		
	######################################

	#### JOINED: ALL (ALL lncRNAs and genes from original annotation file)
	echo -e "\n\nJOINED: ALL: Join genes and lncRNAs without any filter..."
	cd $FL/Joined/ALL/r

	## ids
	echo -e "\nCreate ids file..."
	cat $FL/Genes/ORIGINAL_GENES_ids.txt $FL/LncRNAs/r/POTENTIAL_LNCRNAS_ids_pred.txt > ALL_ids.txt
	## gtf
	echo -e "\nCreate gtf file..."
	cat $FL/Genes/ORIGINAL_GENES.gtf $FL/LncRNAs/r/POTENTIAL_LNCRNAS_pred.gtf > ALL.gtf
	## fasta
	echo -e "\nCreate fasta file..."
	gffread -w ALL_temp.fasta -W -g $2/Genome/$5.fa ALL.gtf
	sed '/^>/s/ .*//' ALL_temp.fasta > ALL.fasta 
	rm ALL_temp.fasta
	## tsv
	echo -e "\nCreate tsv file..."
	cp $FL/Genes/ORIGINAL_GENES.tsv ./
	cp $FL/LncRNAs/r/POTENTIAL_LNCRNAS_pred.tsv ./
	tail -n +2 POTENTIAL_LNCRNAS_pred.tsv > POTENTIAL_LNCRNAS_pred_mod.tsv
	cat ORIGINAL_GENES.tsv POTENTIAL_LNCRNAS_pred_mod.tsv > ALL.tsv
	rm POTENTIAL_LNCRNAS_pred_mod.tsv

	######################################

	#### JOINED: FILTERED (u, x (not overlap with coding transcripts) and genes from original annotation file)
	echo -e "\n\nJOINED: FILTERED: Join genes and lncRNAs but filter lncRNAs by class code (x and u) and overlap (x)..."
	cd $FL/Joined/FILTERED/r

	echo -e "\nFilter the LncRNAs..."
	## GTF/TSV LncRNAs
	cp $FL/LncRNAs/r/POTENTIAL_LNCRNAS_pred.gtf ./
	cp $FL/LncRNAs/r/POTENTIAL_LNCRNAS_pred.tsv ./
	## Select the u and x lncRNAs.
	tail -n +2 POTENTIAL_LNCRNAS_pred.tsv | awk '$8 == "u" || $8 == "x" {print $6}' > POTENTIAL_LNCRNAS_ids_pred_x_u.txt
	Filter_GTF.py \
		--gtf-initial $FL/Joined/FILTERED/r/POTENTIAL_LNCRNAS_pred.gtf \
		--gtf-final $FL/Joined/FILTERED/r/POTENTIAL_LNCRNAS_pred_x_u.gtf \
		--ids $FL/Joined/FILTERED/r/POTENTIAL_LNCRNAS_ids_pred_x_u.txt
	## Convert GTF to BED.
	cat POTENTIAL_LNCRNAS_pred_x_u.gtf | convert2bed -i gtf --attribute-key=transcript_id | awk '$8 == "transcript" {print $0}' > POTENTIAL_LNCRNAS_pred_x_u.bed
	rm POTENTIAL_LNCRNAS_pred.tsv

	## GTF/TSV Genes.
	cp $FL/Genes/ORIGINAL_GENES.gtf ./
	## Convert GTF to BED.
	cat ORIGINAL_GENES.gtf | convert2bed -i gtf --attribute-key=transcript_id | awk '$8 == "transcript" {print $0}' > ORIGINAL_GENES.bed

	## Intersect LncRNAs (u and x) BED file with Genes BED file.
	bedtools intersect -a POTENTIAL_LNCRNAS_pred_x_u.bed -b ORIGINAL_GENES.bed -s -wao -nonamecheck | awk -F"\t" '$21 == 0 {print $4}' > POTENTIAL_LNCRNAS_ids_pred_x_u_filtered.txt
	rm POTENTIAL_LNCRNAS_pred_x_u.bed ORIGINAL_GENES.bed POTENTIAL_LNCRNAS_ids_pred_x_u.txt

	## ids
	echo -e "\nCreate ids file..."
	cp $FL/Genes/ORIGINAL_GENES_ids.txt ./
	cat ORIGINAL_GENES_ids.txt POTENTIAL_LNCRNAS_ids_pred_x_u_filtered.txt > FILTERED_ids.txt 

	## gtf: Select lncRNAs which don't overlap with coding transcripts and join genes gtf and filtered lncRNAs gtf.
	echo -e "\nCreate gtf file..."
	Filter_GTF.py \
		--gtf-initial $FL/Joined/FILTERED/r/POTENTIAL_LNCRNAS_pred.gtf \
		--gtf-final $FL/Joined/FILTERED/r/POTENTIAL_LNCRNAS_pred_x_u_filtered.gtf \
		--ids $FL/Joined/FILTERED/r/POTENTIAL_LNCRNAS_ids_pred_x_u_filtered.txt
	cat ORIGINAL_GENES.gtf POTENTIAL_LNCRNAS_pred_x_u_filtered.gtf > FILTERED.gtf
	rm POTENTIAL_LNCRNAS_pred.gtf POTENTIAL_LNCRNAS_pred_x_u.gtf

	## fasta
	echo -e "\nCreate fasta file..."
	gffread -w FILTERED_temp.fasta -W -g $2/Genome/$5.fa FILTERED.gtf
	sed '/^>/s/ .*//' FILTERED_temp.fasta > FILTERED.fasta 
	rm FILTERED_temp.fasta

	## tsv
	echo -e "\nCreate tsv file..."
	cp $FL/Genes/ORIGINAL_GENES.tsv ./
	gffread -w PL_temp.fasta -W -g $2/Genome/$5.fa POTENTIAL_LNCRNAS_pred_x_u_filtered.gtf
	sed '/^>/s/ .*//' PL_temp.fasta > PL.fasta 
	rm PL_temp.fasta
	GFF3orGTF2TABLE.py \
		--gtf $FL/Joined/FILTERED/r/POTENTIAL_LNCRNAS_pred_x_u_filtered.gtf \
		--fasta $FL/Joined/FILTERED/r/PL.fasta \
		--tab $FL/Joined/FILTERED/r/POTENTIAL_LNCRNAS_pred_x_u_filtered.tsv \
		--mode 'assembly' \
		--sep $'\t'
	rm PL.fasta
	tail -n +2 POTENTIAL_LNCRNAS_pred_x_u_filtered.tsv > POTENTIAL_LNCRNAS_pred_x_u_filtered_mod.tsv
	cat ORIGINAL_GENES.tsv POTENTIAL_LNCRNAS_pred_x_u_filtered_mod.tsv > FILTERED.tsv
	rm POTENTIAL_LNCRNAS_pred_x_u_filtered_mod.tsv

	#########################################################################

	####### DETERMINE WHICH SEQUENCES ARE REDUNDANT.

	echo -e "\n\n\nDETERMINE WHICH SEQUENCES ARE REDUNDANT..."

	###########################################################################
	##
	## Determine which sequences are redundant.
	## 	- CD-HIT (https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide)
	## 	- AGAT (https://github.com/NBISweden/AGAT)
	##	- CGAT (https://github.com/cgat-developers/cgat-apps)
	##
	###########################################################################

	cd $NRA

	###### CD-HIT
	echo -e "\nExecuting the cd-hit (90% similarity)..."
	cd-hit-est -i $FL/LncRNAs/r/POTENTIAL_LNCRNAS_pred.fasta -o ./CDHIT/POTENTIAL_LNCRNAS_filt_90.fasta -c 0.90 -n 9 -d 999 -M 0 -T $6 1> ./CDHIT/stdout_90.txt
	grep ">" ./CDHIT/POTENTIAL_LNCRNAS_filt_90.fasta | sed 's/>//g' > ./CDHIT/POTENTIAL_LNCRNAS_filt_90_ids.txt

	echo -e "\nExecuting the cd-hit (95% similarity)..."
	cd-hit-est -i $FL/LncRNAs/r/POTENTIAL_LNCRNAS_pred.fasta -o ./CDHIT/POTENTIAL_LNCRNAS_filt_95.fasta -c 0.95 -n 10 -d 999 -M 0 -T $6 1> ./CDHIT/stdout_95.txt
	grep ">" ./CDHIT/POTENTIAL_LNCRNAS_filt_95.fasta | sed 's/>//g' > ./CDHIT/POTENTIAL_LNCRNAS_filt_95_ids.txt

	###### AGAT
	echo -e "\nExecuting the AGAT (keep longest isoform)..."
	source activate AGAT
	agat_sp_keep_longest_isoform.pl -gff $FL/LncRNAs/r/POTENTIAL_LNCRNAS_pred.gtf -o ./AGAT/POTENTIAL_LNCRNAS_filt.gff3 2> ./AGAT/progress_bar.txt 1> ./AGAT/stdout.txt
	gffread ./AGAT/POTENTIAL_LNCRNAS_filt.gff3 -T -F --keep-exon-attrs > ./AGAT/POTENTIAL_LNCRNAS_filt_sort.gtf
	gffread -w ./AGAT/POTENTIAL_LNCRNAS_filt_sort_temp.fasta -W -g $2/Genome/$5.fa ./AGAT/POTENTIAL_LNCRNAS_filt_sort.gtf
	sed '/^>/s/ .*//' ./AGAT/POTENTIAL_LNCRNAS_filt_sort_temp.fasta > ./AGAT/POTENTIAL_LNCRNAS_filt_sort_temp_mod.fasta 
	grep ">" ./AGAT/POTENTIAL_LNCRNAS_filt_sort_temp_mod.fasta | sed 's/>//g' > ./AGAT/POTENTIAL_LNCRNAS_filt_sort_ids.txt
	rm ./AGAT/POTENTIAL_LNCRNAS_filt_sort_temp.fasta
	rm ./AGAT/POTENTIAL_LNCRNAS_filt_sort_temp_mod.fasta
	conda deactivate

	###### CGAT
	echo -e "\nExecuting the CGAT (keep longest isoform)..."
	source activate cgat-a
	cgat gtf2gtf --method=sort --sort-order=contig+gene -I $FL/LncRNAs/r/POTENTIAL_LNCRNAS_pred.gtf | cgat gtf2gtf --method=filter --filter-method=longest-transcript -S ./CGAT/POTENTIAL_LNCRNAS_sort_filt.gtf 1> ./CGAT/stdout.txt
	gffread ./CGAT/POTENTIAL_LNCRNAS_sort_filt.gtf -T -F --keep-exon-attrs > ./CGAT/POTENTIAL_LNCRNAS_sort_filt_sort.gtf
	gffread -w ./CGAT/POTENTIAL_LNCRNAS_sort_filt_sort_temp.fasta -W -g $2/Genome/$5.fa ./CGAT/POTENTIAL_LNCRNAS_sort_filt_sort.gtf
	sed '/^>/s/ .*//' ./CGAT/POTENTIAL_LNCRNAS_sort_filt_sort_temp.fasta > ./CGAT/POTENTIAL_LNCRNAS_sort_filt_sort_temp_mod.fasta 
	grep ">" ./CGAT/POTENTIAL_LNCRNAS_sort_filt_sort_temp_mod.fasta | sed 's/>//g' > ./CGAT/POTENTIAL_LNCRNAS_sort_filt_sort_ids.txt
	rm ./CGAT/POTENTIAL_LNCRNAS_sort_filt_sort_temp.fasta
	rm ./CGAT/POTENTIAL_LNCRNAS_sort_filt_sort_temp_mod.fasta
	conda deactivate

	#########################################################################

	####### CREATE THE NON-REDUNDANT LNCRNAS DATABASE.

	echo -e "\nCREATE THE NON-REDUNDANT LNCRNAS DATABASE...\n"
	Rscript $3/Create_customed_LncRNA_db_NR.R $1

	#########################################################################

	####### NON-REDUNDANT: BUILD THE FOLDERS FOR DOWNSTREAM ANALYSIS.

	echo -e "\n\n\nNON-REDUNDANT: BUILD THE FOLDERS FOR DOWNSTREAM ANALYSIS..."

	######################################

	#### LNCRNAS
	echo -e "\n\nLNCRNAS..."
	cd $FL/LncRNAs/nr

	## ids
	echo -e "\nCreate new ids file..."
	tail -n +2 $DB/Database_LncRNAs_NR.tsv | awk '{print $1}' > POTENTIAL_LNCRNAS_ids_pred.txt
	## fasta
	echo -e "\nCreate new fasta file..."
	seqtk subseq $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.fasta POTENTIAL_LNCRNAS_ids_pred.txt > POTENTIAL_LNCRNAS_pred.fasta
	## gtf
	echo -e "\nCreate new gtf file..."
	Filter_GTF.py \
		--gtf-initial $P1/Potential_lncRNAs/POTENTIAL_LNCRNAS.gtf \
		--gtf-final $FL/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.gtf \
		--ids $FL/LncRNAs/nr/POTENTIAL_LNCRNAS_ids_pred.txt
	## tsv
	echo -e "\nCreate new tsv file..."
	GFF3orGTF2TABLE.py \
		--gtf $FL/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.gtf \
		--fasta $FL/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.fasta \
		--tab $FL/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.tsv \
		--mode 'assembly' \
		--sep $'\t'
		

	######################################

	#### JOINED: ALL (ALL lncRNAs and genes from original annotation file)
	echo -e "\n\nJOINED: ALL: Join genes and lncRNAs without any filter..."
	cd $FL/Joined/ALL/nr

	## ids
	echo -e "\nCreate ids file..."
	cat $FL/Genes/ORIGINAL_GENES_ids.txt $FL/LncRNAs/nr/POTENTIAL_LNCRNAS_ids_pred.txt > ALL_ids.txt
	## gtf
	echo -e "\nCreate gtf file..."
	cat $FL/Genes/ORIGINAL_GENES.gtf $FL/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.gtf > ALL.gtf
	## fasta
	echo -e "\nCreate fasta file..."
	gffread -w ALL_temp.fasta -W -g $2/Genome/$5.fa ALL.gtf
	sed '/^>/s/ .*//' ALL_temp.fasta > ALL.fasta 
	rm ALL_temp.fasta
	## tsv
	echo -e "\nCreate tsv file..."
	cp $FL/Genes/ORIGINAL_GENES.tsv ./
	cp $FL/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.tsv ./
	tail -n +2 POTENTIAL_LNCRNAS_pred.tsv > POTENTIAL_LNCRNAS_pred_mod.tsv
	cat ORIGINAL_GENES.tsv POTENTIAL_LNCRNAS_pred_mod.tsv > ALL.tsv
	rm POTENTIAL_LNCRNAS_pred_mod.tsv

	######################################

	#### JOINED: FILTERED (u, x (not overlap with coding transcripts) and genes from original annotation file)
	echo -e "\n\nJOINED: FILTERED: Join genes and lncRNAs but filter lncRNAs by class code (x and u) and overlap (x)..."
	cd $FL/Joined/FILTERED/nr

	echo -e "\nFilter the LncRNAs..."
	## GTF/TSV LncRNAs
	cp $FL/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.gtf ./
	cp $FL/LncRNAs/nr/POTENTIAL_LNCRNAS_pred.tsv ./
	## Select the u and x lncRNAs.
	tail -n +2 POTENTIAL_LNCRNAS_pred.tsv | awk '$8 == "u" || $8 == "x" {print $6}' > POTENTIAL_LNCRNAS_ids_pred_x_u.txt
	Filter_GTF.py \
		--gtf-initial $FL/Joined/FILTERED/nr/POTENTIAL_LNCRNAS_pred.gtf \
		--gtf-final $FL/Joined/FILTERED/nr/POTENTIAL_LNCRNAS_pred_x_u.gtf \
		--ids $FL/Joined/FILTERED/nr/POTENTIAL_LNCRNAS_ids_pred_x_u.txt
	## Convert GTF to BED.
	cat POTENTIAL_LNCRNAS_pred_x_u.gtf | convert2bed -i gtf --attribute-key=transcript_id | awk '$8 == "transcript" {print $0}' > POTENTIAL_LNCRNAS_pred_x_u.bed
	rm POTENTIAL_LNCRNAS_pred.tsv

	## GTF/TSV Genes.
	cp $FL/Genes/ORIGINAL_GENES.gtf ./
	## Convert GTF to BED.
	cat ORIGINAL_GENES.gtf | convert2bed -i gtf --attribute-key=transcript_id | awk '$8 == "transcript" {print $0}' > ORIGINAL_GENES.bed

	## Intersect LncRNAs (u and x) BED file with Genes BED file.
	bedtools intersect -a POTENTIAL_LNCRNAS_pred_x_u.bed -b ORIGINAL_GENES.bed -s -wao -nonamecheck | awk -F"\t" '$21 == 0 {print $4}' > POTENTIAL_LNCRNAS_ids_pred_x_u_filtered.txt
	rm POTENTIAL_LNCRNAS_pred_x_u.bed ORIGINAL_GENES.bed POTENTIAL_LNCRNAS_ids_pred_x_u.txt

	## ids
	echo -e "\nCreate ids file..."
	cp $FL/Genes/ORIGINAL_GENES_ids.txt ./
	cat ORIGINAL_GENES_ids.txt POTENTIAL_LNCRNAS_ids_pred_x_u_filtered.txt > FILTERED_ids.txt 

	## gtf: Select lncRNAs which don't overlap with coding transcripts and join genes gtf and filtered lncRNAs gtf.
	echo -e "\nCreate gtf file..."
	Filter_GTF.py \
		--gtf-initial $FL/Joined/FILTERED/nr/POTENTIAL_LNCRNAS_pred.gtf \
		--gtf-final $FL/Joined/FILTERED/nr/POTENTIAL_LNCRNAS_pred_x_u_filtered.gtf \
		--ids $FL/Joined/FILTERED/nr/POTENTIAL_LNCRNAS_ids_pred_x_u_filtered.txt
	cat ORIGINAL_GENES.gtf POTENTIAL_LNCRNAS_pred_x_u_filtered.gtf > FILTERED.gtf
	rm POTENTIAL_LNCRNAS_pred.gtf POTENTIAL_LNCRNAS_pred_x_u.gtf

	## fasta
	echo -e "\nCreate fasta file..."
	gffread -w FILTERED_temp.fasta -W -g $2/Genome/$5.fa FILTERED.gtf
	sed '/^>/s/ .*//' FILTERED_temp.fasta > FILTERED.fasta 
	rm FILTERED_temp.fasta

	## tsv
	echo -e "\nCreate tsv file..."
	cp $FL/Genes/ORIGINAL_GENES.tsv ./
	gffread -w PL_temp.fasta -W -g $2/Genome/$5.fa POTENTIAL_LNCRNAS_pred_x_u_filtered.gtf
	sed '/^>/s/ .*//' PL_temp.fasta > PL.fasta 
	rm PL_temp.fasta
	GFF3orGTF2TABLE.py \
		--gtf $FL/Joined/FILTERED/nr/POTENTIAL_LNCRNAS_pred_x_u_filtered.gtf \
		--fasta $FL/Joined/FILTERED/nr/PL.fasta \
		--tab $FL/Joined/FILTERED/nr/POTENTIAL_LNCRNAS_pred_x_u_filtered.tsv \
		--mode 'assembly' \
		--sep $'\t'
	rm PL.fasta
	tail -n +2 POTENTIAL_LNCRNAS_pred_x_u_filtered.tsv > POTENTIAL_LNCRNAS_pred_x_u_filtered_mod.tsv
	cat ORIGINAL_GENES.tsv POTENTIAL_LNCRNAS_pred_x_u_filtered_mod.tsv > FILTERED.tsv
	rm POTENTIAL_LNCRNAS_pred_x_u_filtered_mod.tsv
}

task_LncRNAs_prediction_STEP-FINAL_circos(){
	
	####### VARIABLES
	P7=$1/STEP-FINAL
	
	####### DIRECTORY
	mkdir -p $P7/Figures
	mkdir -p $P7/Figures/Circos
	mkdir -p $P7/Figures/Circos/r
	mkdir -p $P7/Figures/Circos/nr

	####### PIPELINE
	cd $P7/Figures/Circos

	## Get chromosome sizes.
	echo -e "\nGet chromosome sizes..."
	cp $2/Genome/$4.fa ./
	samtools faidx $4.fa
	cut -f1,2 $4.fa.fai > sizes_genome_temp.txt
	awk '{print $1"\t"0"\t"$2"\t"$1"\t"$1}' sizes_genome_temp.txt > sizes_genome.txt
	rm sizes_genome_temp.txt $4.fa $4.fa.fai

	## REDUNDANT: Visualize the lncRNAs distribution.
	echo -e "\nREDUNDANT: Visualize the lncRNAs distribution..."
	cp sizes_genome.txt ./r/
	LncRNAs_tab="$P7/Database/Database_LncRNAs.tsv"
	Genes_tab="$P7/Files/Genes/ORIGINAL_GENES.tsv"
	Rscript $3/Genome-wide_distribution.R $4 $LncRNAs_tab $Genes_tab $P7/Figures/Circos/r $2

	## NON-REDUNDANT: Visualize the lncRNAs distribution.
	echo -e "\nNON-REDUNDANT: Visualize the lncRNAs distribution..."
	cp sizes_genome.txt ./nr/
	LncRNAs_tab="$P7/Database/Database_LncRNAs_NR.tsv"
	Genes_tab="$P7/Files/Genes/ORIGINAL_GENES.tsv"
	Rscript $3/Genome-wide_distribution.R $4 $LncRNAs_tab $Genes_tab $P7/Figures/Circos/nr $2
}

"$@"


