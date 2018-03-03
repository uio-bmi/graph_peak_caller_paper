

./call_peaks.sh blank_config.sh ENCSR471GSA 1 DM_JRA chr3R,chr3L,chr2R,chr2L,chrX,chr4 https://www.encodeproject.org/files/ENCFF597CVA/@@download/ENCFF597CVA.bam . 97958719 /media/storage1/dm6/dm3_main_chromosomes.fasta /media/storage1/dm6/wg.xg /media/storage1/dm6/wg.gcsa /media/storage1/dm6/  > dm_jra_output.txt 2 >&1 &
./call_peaks.sh blank_config.sh ENCSR923VWW 1 DM_SQZ chr3R,chr3L,chr2R,chr2L,chrX,chr4 https://www.encodeproject.org/files/ENCFF255JBC/@@download/ENCFF255JBC.bam . 97958719 /media/storage1/dm6/dm3_main_chromosomes.fasta /media/storage1/dm6/wg.xg /media/storage1/dm6/wg.gcsa /media/storage1/dm6/ > dm_sqz_output.txt 2 >&1 &
./call_peaks.sh blank_config.sh ENCSR978WED 1 DM_JIM chr3R,chr3L,chr2R,chr2L,chrX,chr4 . . 97958719 /media/storage1/dm6/dm3_main_chromosomes.fasta /media/storage1/dm6/wg.xg /media/storage1/dm6/wg.gcsa /media/storage1/dm6/ > dm_jim_output.txt 2 >&1 &
./call_peaks.sh blank_config.sh ENCSR082RBU 1 DM_ANTP chr3R,chr3L,chr2R,chr2L,chrX,chr4 https://www.encodeproject.org/files/ENCFF713VJE/@@download/ENCFF713VJE.bam . 97958719 /media/storage1/dm6/dm3_main_chromosomes.fasta /media/storage1/dm6/wg.xg /media/storage1/dm6/wg.gcsa /media/storage1/dm6/ > dm_antp_output.txt 2 >&1 &

wait

./generate_paper_figures_and_tables.sh  ENCSR471GSA 1 DM_JRA chr3R,chr3L,chr2R,chr2L,chrX,chr4 . /media/storage1/dm6/ /media/storage1/dm6/dm3_main_chromosomes.fasta 97958719 &
./generate_paper_figures_and_tables.sh  ENCSR923VWW 1 DM_SQZ chr3R,chr3L,chr2R,chr2L,chrX,chr4 . /media/storage1/dm6/ /media/storage1/dm6/dm3_main_chromosomes.fasta 97958719 &
./generate_paper_figures_and_tables.sh  ENCSR082RBU 1 DM_ANTP chr3R,chr3L,chr2R,chr2L,chrX,chr4 . /media/storage1/dm6/ /media/storage1/dm6/dm3_main_chromosomes.fasta 97958719 &
./generate_paper_figures_and_tables.sh  ENCSR978WED 1 DM_JIM chr3R,chr3L,chr2R,chr2L,chrX,chr4 . /media/storage1/dm6/ /media/storage1/dm6/dm3_main_chromosomes.fasta 97958719 &

wait
