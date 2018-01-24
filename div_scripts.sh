
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


for chr in $(seq 1 22; echo X; echo Y);
do
	#python3 ~/dev/graph_peak_caller/graph_peak_caller.py create_ob_graph $chr.json $chr.obg
	python3 ~/dev/graph_peak_caller/graph_peak_caller.py create_ob_numpy_graph $chr.json $chr
done

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
