
vg snarls 1.vg > 1.snarls &
vg snarls 2.vg > 2.snarls &
vg snarls 3.vg > 3.snarls &
vg snarls 4.vg > 4.snarls &
vg snarls 5.vg > 5.snarls &
vg snarls 6.vg > 6.snarls &
vg snarls 7.vg > 7.snarls &
vg snarls 8.vg > 8.snarls &
vg snarls 9.vg > 9.snarls &
vg snarls 10.vg > 10.snarls &
vg snarls 11.vg > 11.snarls &
vg snarls 12.vg > 12.snarls &
vg snarls 13.vg > 13.snarls &
vg snarls 14.vg > 14.snarls &
vg snarls 15.vg > 15.snarls &
vg snarls 16.vg > 16.snarls &
vg snarls 17.vg > 17.snarls &
vg snarls 18.vg > 18.snarls &
vg snarls 19.vg > 19.snarls &
vg snarls 20.vg > 20.snarls &
vg snarls 21.vg > 21.snarls &
vg snarls 22.vg > 22.snarls &
vg snarls X.vg > X.snarls &
vg snarls Y.vg > Y.snarls &


vg view -Vj 1.vg > 1.json &
vg view -Vj 2.vg > 2.json &
vg view -Vj 3.vg > 3.json &
vg view -Vj 4.vg > 4.json &
vg view -Vj 5.vg > 5.json &
vg view -Vj 6.vg > 6.json &
vg view -Vj 7.vg > 7.json &
vg view -Vj 8.vg > 8.json &
vg view -Vj 9.vg > 9.json &
vg view -Vj 10.vg > 10.json &
vg view -Vj 11.vg > 11.json &
vg view -Vj 12.vg > 12.json &
vg view -Vj 13.vg > 13.json &
vg view -Vj 14.vg > 14.json &
vg view -Vj 15.vg > 15.json &
vg view -Vj 16.vg > 16.json &
vg view -Vj 17.vg > 17.json &
vg view -Vj 18.vg > 18.json &
vg view -Vj 19.vg > 19.json &
vg view -Vj 20.vg > 20.json &
vg view -Vj 21.vg > 21.json &
vg view -Vj 22.vg > 22.json &
vg view -Vj X.vg > X.json &
vg view -Vj X.vg > Y.json &




~/bin/vg1.5 view -Vj 14.vg > 14.json &
~/bin/vg1.5 view -Vj 17.vg > 17.json &
~/bin/vg1.5 view -Vj 12.vg > 12.json &
~/bin/vg1.5 view -Vj 10.vg > 10.json &
~/bin/vg1.5 view -Vj 8.vg > 8.json &
~/bin/vg1.5 view -Vj 3.vg > 3.json &


for chr in $(seq 16 22);
do
	graph_peak_caller create_ob_graph $chr.json $chr.nobg &
done

for chr in $(seq 16 21; echo X; echo Y);
do
	graph_peak_caller find_linear_path $chr.json $chr.nobg $chr ${chr}_linear_pathv2.interval &
done

for chr in $(seq 1 22; echo X; echo Y);
do
    git mv BATF_ENCSR000BGT/1/chr${chr}_sequences.fasta BATF_ENCSR000BGT/1/${chr}_sequences.fasta
done

for chr in $(seq 1 22; echo X; echo Y); do     git mv CTCF_ENCSR000DUB/1/chr${chr}_sequences.fasta CTCF_ENCSR000DUB/1/${chr}_sequences.fasta; done
for chr in $(seq 1 22; echo X; echo Y); do     git mv EBF_ENCSR000BGU/1/chr${chr}_sequences.fasta EBF_ENCSR000BGU/1/${chr}_sequences.fasta; done
for chr in $(seq 1 22; echo X; echo Y); do     git mv ELF1_ENCSR000BMZ/1/chr${chr}_sequences.fasta ELF1_ENCSR000BMZ/1/${chr}_sequences.fasta; done
for chr in $(seq 1 22; echo X; echo Y); do     git mv GATA3_ENCSR000EXZ/1/chr${chr}_sequences.fasta GATA3_ENCSR000EXZ/1/${chr}_sequences.fasta; done
for chr in $(seq 1 22; echo X; echo Y); do     git mv JUN_ENCSR996DUT/1/chr${chr}_sequences.fasta JUN_ENCSR996DUT/1/${chr}_sequences.fasta; done
for chr in $(seq 1 22; echo X; echo Y); do     git mv MAX_ENCSR521IID/1/chr${chr}_sequences.fasta MAX_ENCSR521IID/1/${chr}_sequences.fasta; done
for chr in $(seq 1 22; echo X; echo Y); do     git mv NR3C1_ENCSR070VXO/1/chr${chr}_sequences.fasta NR3C1_ENCSR070VXO/1/${chr}_sequences.fasta; done
for chr in $(seq 1 22; echo X; echo Y); do     git mv SP1_ENCSR000ECD/1/chr${chr}_sequences.fasta SP1_ENCSR000ECD/1/${chr}_sequences.fasta; done
for chr in $(seq 1 22; echo X; echo Y); do     git mv SRF_ENCSR000BIV/1/chr${chr}_sequences.fasta SRF_ENCSR000BIV/1/${chr}_sequences.fasta; done
for chr in $(seq 1 22; echo X; echo Y); do     git mv USF2_ENCSR000ECD/1/chr${chr}_sequences.fasta USF2_ENCSR000ECD/1/${chr}_sequences.fasta; done


# Create linear maps
for chr in $(seq 10 22; echo X; echo Y);
do
	python3 ~/dev/graph_peak_caller/graph_peak_caller.py create_linear_map $chr.obg $chr.snarls linear_map_$chr &
done


# Create heads/tails
for chr in $(seq 1 22; echo X; echo Y);
do
	vg stats -H $chr.vg > $chr_head_node.txt
	vg stats -T $chr.vg > $chr_tail_node.txt
done


# Get node range
for chr in $(seq 1 22; echo X; echo Y);
do
	vg stats -r $chr.vg  | awk '{print $2}' > node_range_$chr.txt
done


# Fix wrong format
for chr in $(seq 2 22; echo X; echo Y);
do
	cat \{$chr\}_node_range.txt  | awk '{print $2}' > node_range_$chr.txt
done


# DM
for chr in $(echo chr3R; echo chr3L; echo chr2R; echo chr2L; echo chrX; echo chr4);
do
    vg view -Vj $chr.vg > $chr.json &
done


for chr in $(echo chr3R; echo chr3L; echo chr2R; echo chr2L; echo chrX; echo chr4);
do
    vg snarls $chr.vg > $chr.snarls &
done


for chr in $(echo chr3R; echo chr3L; echo chr2R; echo chr2L; echo chrX; echo chr4);
do
	vg stats -r $chr.vg  | awk '{print $2}' > node_range_$chr.txt
done


for chr in $(echo chr3R; echo chr3L; echo chr2R; echo chr2L; echo chrX; echo chr4);
do
	graph_peak_caller create_ob_graph $chr.json $chr.nobg &
done


for chr in $(echo chr3R; echo chr3L; echo chr2R; echo chr2L; echo chrX; echo chr4);
do
	graph_peak_caller find_linear_path $chr.json $chr.nobg $chr ${chr}_linear_pathv2.interval &
done

# Create linear maps
for chr in $(echo chr3R; echo chr3L; echo chr2R; echo chr2L; echo chrX; echo chr4);
do
	graph_peak_caller create_linear_map $chr.nobg $chr.snarls linear_map_$chr &
done