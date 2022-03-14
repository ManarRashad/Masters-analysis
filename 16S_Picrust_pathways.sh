for i in *biom; do filter_otus_from_otu_table.py -i $i -o $i.out --negate_ids_to_exclude -e /home/mh_user/anaconda2/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta; done
# PICRUSt: normalize OTU table by predicted 16S copy numbers.
for i in *out; do normalize_by_copy_number.py -i $i -o $i.norm; done
#PICRUSt: Predict abundances of KEGG orthologs (KOs) in each sample.
for i in *norm; do predict_metagenomes.py -i $i -o $i.ko; done
#PICRUSt: Collapse KOs into KEGG Pathways.
for i in *ko; do categorize_by_function.py -i $i -l 3 -c KEGG_Pathways -o $i.L3; done
for i in *L3; do biom_to_stamp.py -m KEGG_Pathways $i > $i.spf; done
#convert out to txt
for i in *out; do biom convert -i $i -o $i.txt --to-tsv; done
#convert norm to txt
for i in *norm; do biom convert -i $i -o $i.txt --to-tsv; done
#convert ko to txt
for i in *ko; do biom convert -i $i -o $i.txt --to-tsv; done
#convert L3 to txt
for i in *L3; do biom convert -i $i -o $i.txt --to-tsv; done

