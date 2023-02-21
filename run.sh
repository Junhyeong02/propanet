#!/bin/bash

# <Replace with the description and/or purpose of this shell script.>
expFile=$1 # "data/DEG.AtGenExpress.signed_zstats.heat_shoots"
seedFile=$2 # "data/DEG.AtGenExpress.signed_binary.heat_shoots"
TFliFile=$3 # "data/Ath_TF_list.gene"
gSet=$4 # (optional) additional gene list file if the user wants
resD="result"

if [ -z ${expFile} ]; then
    expFile="data/DEG.AtGenExpress.signed_zstats.heat_shoots"
fi
if [ -z ${seedFile} ]; then
    seedFile="data/DEG.AtGenExpress.signed_binary.heat_shoots"
fi
if [ -z ${TFliFile} ]; then
	TFliFile="data/Ath_TF_list.gene"
fi

prefix="PropaNet"
n=$(head -n 1 $seedFile| awk '{print NF}')

if [ -e result/intermediate_results ]; then
    echo
else
    mkdir result/intermediate_results
fi

# Weighted template network construction
# python network_weight.py -nwk data/Ara_GSE5628_2_normFalse_GRN_0.5_tp1.tsv -exp ${expFile} -o ${resD}/intermediate_results/${prefix}.Ara_GSE5628_2_normFalse_GRN_0.5_tp1.tsv -p 10 || exit 1
# python network_weight.py -nwk data/Ara_GSE5628_2_normFalse_GRN_0.5_tp2.tsv -exp ${expFile} -o ${resD}/intermediate_results/${prefix}.Ara_GSE5628_2_normFalse_GRN_0.5_tp2.tsv -p 10 || exit 1
# python network_weight.py -nwk data/Ara_GSE5628_2_normFalse_GRN_0.5_tp3.tsv -exp ${expFile} -o ${resD}/intermediate_results/${prefix}.Ara_GSE5628_2_normFalse_GRN_0.5_tp3.tsv -p 10 || exit 1
# python network_weight.py -nwk data/Ara_GSE5628_2_normFalse_GRN_0.5_tp4.tsv -exp ${expFile} -o ${resD}/intermediate_results/${prefix}.Ara_GSE5628_2_normFalse_GRN_0.5_tp4.tsv -p 10 || exit 1

wait
echo "Extract optimal TF list"

# Extract optimal TF list
python new_TF_adding_NP_noCtrl.py -TFliFile ${TFliFile} -nwkFile ${resD}/intermediate_results/${prefix}.Ara_GSE5628_2_normFalse_GRN_0.5_tp -expFile ${expFile} -binFile ${seedFile} -p 10 -cond ${prefix} -outD ${resD}/intermediate_results || exit 1
sleep 10

wait 
echo "Extract target genes of the optimal TFs"

# Extract target genes of the optimal TFs
# for ((i=1;i<=$[$n-1];i++));do
#     python Target_genes.py ${resD}/intermediate_results/${prefix}.nwk.t$i ${resD}/intermediate_results/${prefix}.DEGli.t$i $TFliFile ${resD}/intermediate_results/${prefix}.TF_rank.t$i.trim ${gSet} ${resD}/TG $i&
# done; wait

python Target_genes.py ${resD}/intermediate_results/${prefix}.nwk.t1 ${resD}/intermediate_results/${prefix}.DEGli.t1 $TFliFile ${resD}/intermediate_results/${prefix}.TF_rank.t1.trim ${gSet} ${resD}/TG 1 &
python Target_genes.py ${resD}/intermediate_results/${prefix}.nwk.t2 ${resD}/intermediate_results/${prefix}.DEGli.t2 $TFliFile ${resD}/intermediate_results/${prefix}.TF_rank.t2.trim ${gSet} ${resD}/TG 2 &
python Target_genes.py ${resD}/intermediate_results/${prefix}.nwk.t3 ${resD}/intermediate_results/${prefix}.DEGli.t3 $TFliFile ${resD}/intermediate_results/${prefix}.TF_rank.t3.trim ${gSet} ${resD}/TG 3 &
python Target_genes.py ${resD}/intermediate_results/${prefix}.nwk.t4 ${resD}/intermediate_results/${prefix}.DEGli.t4 $TFliFile ${resD}/intermediate_results/${prefix}.TF_rank.t4.trim ${gSet} ${resD}/TG 4

wait

echo "Networks are comprised of resulting TFs/TGs"

# Final Result : Networks are comprised of resulting TFs/TGs
python makeTGDesc.py ${resD}

wait
echo "Propanet done"

# cd ${resD}
# awk '{print}' *.trim|sort|uniq > ../propanet_${cond}_TFs
# awk '{print}' TG/*.TG|sort|uniq > ../propanet_${cond}_TGs
# cd ..
# cat propanet_${cond}_TFs propanet_${cond}_TGs|sort|uniq >propanet_${cond}_TFTGs

# # Exit with an explicit exit status.
# exit 0
