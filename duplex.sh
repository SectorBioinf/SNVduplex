#!/bin/bash
cohort=kid
### set-based cohort testing
plink --noweb --file "$cohort"_10 --r2 --inter-chr --out LD_"$cohort"chrs --allow-no-sex
cat LD_"$cohort"chrs.ld | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' > tab_delim_LD_"$cohort"chrs.tsv
input=tab_delim_LD_"$cohort"chrs.tsv
output=set_file.set
awk -F '\t' '{printf "%s__%s__%s\n%s\n%s\nEND\n\n", $3, $6, $7, $3, $6}' "$input" > "$output"
plink --noweb --bfile "$cohort"_10 --set-test --set set_file.set --assoc --mperm 10000 --out "$cohort"chrs --set-max 99999 --set-p 1 --set-r2 1 --allow-no-sex
cat "$cohort"chrs.assoc.set.mperm | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > tab_del_"$cohort"chrs.assoc.set.mperm
row_res=tab_del_"$cohort"chrs.assoc.set.mperm
cat $row_res | awk '{print $6}' > 'res_only_6_col'$row_res'.txt'
sed -e '1,2d' 'res_only_6_col'$row_res'.txt' > 'res_o_6_col'$row_res'.txt'
while IFS= read line; do echo -e "$line" | sed "s/|/\t/g"; done < 'res_o_6_col'$row_res'.txt' > 'res_only_6_col'$row_res'_splited.txt'
while IFS= read line; do echo -e "$line" | sed "s/:/\t/g"; done < 'res_only_6_col'$row_res'_splited.txt' > all_dropped_cols
paste 'res_only_6_col'$row_res'_splited.txt' all_dropped_cols > 'all_necessary_info'$row_res
cat all_dropped_cols | awk '{print $1 "\t" $2 "\t" $2 }' > 'Segmented_file'$row_res
cat 'Segmented_file'$row_res| awk '{print $3}' > increasing
input_file=increasing
output_file=out_increasing
while read -r number; do  if [[ $number =~ ^[0-9]+$ ]]; then ((number++)); echo "$number" >> "$output_file"; else echo "Некорр. знач. $number"; fi; done < "$input_file"
cat Segmented_filetab_del_"$cohort"chrs.assoc.set.mperm | awk '{print $1 "\t" $2}' > 1_2_col_out_increasing
paste 1_2_col_out_increasing out_increasing > incr
head incr
### annot genes
sort-bed incr > new_1.bed
sort-bed ref.bed > new_2.bed
sort-bed --check-sort new_1.bed
sort-bed --check-sort new_2.bed
bedmap --echo --echo-map-id --delim '\t' new_1.bed new_2.bed > new_annot.bed
cat new_1.bed | wc -l
cat new_annot.bed | wc -l
sed -i.bkp '1i CHR\tPOS1_N1\tPOS2_N1\tGENE_N1' output_file_first_snp.xls
head -20 output_file_first_snp.xls incr
cat all_dropped_cols | awk '{print $5 "\t" $6}' > increasing_second_snp
cat increasing_second_snp | awk '{print $2}' > incr_only_pos_second_snp
input_f=incr_only_pos_second_snp
output_f=out_incr_only_pos_second_snp
while IFS= read -r line; do if [[ -z "$line" ]]; then echo -en "\n" >> "$output_f"; else sum=$(($line+1)); echo "$sum" >> "$output_f"; fi; done < "$input_f"
paste increasing_second_snp out_incr_only_pos_second_snp > second_part
Rscript narows.R
sort-bed na_droped_second_part > new_1.bed
sort-bed ref.bed > new_2.bed
sort-bed --check-sort new_1.bed
sort-bed --check-sort new_2.bed
bedmap --echo --echo-map-id --delim '\t' new_1.bed new_2.bed > new_annot.bed
cat new_1.bed | wc -l
cat new_annot.bed | wc -l
python merge_annot.py
head -n 40 second_part output_file_second_snp.xls
sed -i.bkp '1i CHR_N2\tPOS1_N2\tPOS2_N2\tGENE_N2' output_file_second_snp.xls
awk '{print $1 "\t" $5 "\t"}' tab_del_"$cohort"chrs.assoc.set.mperm > set_name_p_val
sed -i.bkp '1i SNP1\tSNP2' res_only_6_coltab_del_"$cohort"chrs.assoc.set.mperm_splited.txt
sed -i 's/SET\tNSNP\tNSIG\tISIG\tEMP1\tSNPS//g' tab_del_"$cohort"chrs.assoc.set.mperm
sed -i.bkp 's/\t0\t0\t0\t1\tNA//g' tab_del_"$cohort"chrs.assoc.set.mperm
sed -i.bkp '1i CHR_N2\tPOS1_N2\tPOS2_N2\tGENE_N2' tab_del_"$cohort"chrs.assoc.set.mperm
sed -i.bkp '2d; 3d' tab_del_"$cohort"chrs.assoc.set.mperm
sed -i.bkp '2d' set_name_p_val
paste set_name_p_val res_only_6_coltab_del_"$cohort"chrs.assoc.set.mperm_splited.txt > res_annot
### freq annot
cat "$cohort"chrs.assoc | awk '{print $2 "\t" $5 "\t" $6}' > all_freqAU.bed
cat res_annot | awk '{print $3}' > freqAU1.bed
cat res_annot | awk '{print $4}' > freqAU2.bed
python merge_annot_first_p.py
python merge_annot_second_p.py
sed -i.bkp '1i CHR_N1\tF_A_N1\tF_U_N1' FREQ_AU_1.xls
sed -i.bkp '1i CHR_N2\tF_A_N2\tF_U_N2' FREQ_AU_2.xls
mkdir results 
mv FREQ_AU_2.xls results/
mv FREQ_AU_1.xls results/
mv output_file_second_snp.xls results/
mv output_file_first_snp.xls results/
mv res_annot results/
for i in *; do cat $i | wc -l; done
