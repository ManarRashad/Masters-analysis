for i in *biom; do filter_otus_from_otu_table.py -i $i -o $i.out --negate_ids_to_exclude -e /home/mh_user/anaconda2/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta; done

# PICRUSt: normalize OTU table by predicted 16S copy numbers.

for i in *out; do normalize_by_copy_number.py -i $i -o $i.norm; done

# contribution of ko in norm biom

for i in *norm; do metagenome_contributions.py -i $i -o $i.contributions.tab; done

# convert each norm biom to taxonomy

for i in *norm; do biom convert -i $i -o $i.taxonomy.txt --to-tsv --header-key taxonomy; done
