# processBarcode

Detect and trim barcodes generated by BLESS method using Striped Smith Waterman (https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library).

#create an exception list for mappable reads
bowtie /data/store/yizhu/my_databases/mouse/mm10/mm10.bowtie2 -p 20 -v 2 -r B_J1_0_GTGGCC_R1.fq.gz.dump.fa > map_read_to_genome/B_J1_0_GTGGCC_R1-vs-genome.bt
cut -f1 map_read_to_genome/B_J1_0_GTGGCC_R1-vs-genome.bt | sort | uniq > map_read_to_genome/B_J1_0_GTGGCC_R1-vs-genome.bt.except
bowtie /data/store/yizhu/my_databases/mouse/mm10/mm10.bowtie2 -p 20 -v 2 -r B_J1_0_GTGGCC_R2.fq.gz.dump.fa > map_read_to_genome/B_J1_0_GTGGCC_R2-vs-genome.bt
cut -f1 map_read_to_genome/B_J1_0_GTGGCC_R2-vs-genome.bt | sort | uniq > map_read_to_genome/B_J1_0_GTGGCC_R2-vs-genome.bt.except

#process barcode sequences 
../process_barcode_pe -m 9 -M 1 -p test -s 1 -a map_read_to_genome/B_J1_0_GTGGCC_R1-vs-genome.bt.except -b map_read_to_genome/B_J1_0_GTGGCC_R2-vs-genome.bt.except R1_10000.fa R2_10000.fa > test.log
../process_barcode_se -m 9 -M 1 -p merge -s 1 -a map_read_to_genome/B_J1_0_GTGGCC_R1-vs-genome.bt.except merge.fa > merge.log
