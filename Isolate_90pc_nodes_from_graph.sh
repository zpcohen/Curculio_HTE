awk -F"," '$17 == "TRUE" {print $1","$9","$21}' Ccaryae_node_dense_Mar6.tsv > CcaryaeMar6_node90.tsv

awk -F"," '{print $3}' CcaryaeMar6_node90.tsv | awk -F"e" '{print $1 * (10 ^ $2)}' > CcaryaeMar6_90pc_integerstop.list
awk -F"," '{print $2}' CcaryaeMar6_node90.tsv | awk -F"e" '{print $1 * (10 ^ $2)}' > CcaryaeMar6_90pc_integerstart.list
paste -d"," CcaryaeMar6_90pc_chrom.list CcaryaeMar6_90pc_integerstart.list CcaryaeMar6_90pc_integerstop.list > CcaryaeMar6_90pc_csv

# edit and run identify_hits_by_window.sh
bash identify_hits_by_window_Mar6.sh > CcaryaeMar6_pangenomenodes_90pc.tsv # copy this to server, run with species specific nodes using odgi -position
