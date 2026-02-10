# Check if the file exists
if [ ! -f "cases_vs_controls.counts.tab.gz" ]; then
    echo "Error: cases_vs_controls.counts.tab.gz does not exist."
    exit 1
fi

## get dPSI of exons pass the min_reads cutoff
# col1 : exon, col2: inc-case, col3:exc-case, col4: inc-control, col5: exc-control, col6: sum-case, col7: sum-control, col8: dPSI
zcat cases_vs_controls.counts.tab.gz | awk -F" " '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$2+$3"\t"$4+$5}' | awk -F"\t" '$6>0 && $7>0 {print $0"\t"($2/$6)-($4/$7)}' | sort -u > cases_vs_controls.counts.dPSI
### get |dPSI|>=0.5 exons
cat cases_vs_controls.counts.dPSI | awk -F"\t" '$8<=-0.5 || $8>=0.5 {print $0}' | sort -u > cases_vs_controls.counts.dPSIabove0.5

#### get dPSI values for exons tested by chi-sqr
zcat cases_vs_controls.counts.pVal.FDR.LOR.tab.gz | awk -F"\t" '{print $0"\t"$2+$3"\t"$4+$5}' | awk -F"\t" '$9>0 && $10>0 {print $0"\t"($2/$9)-($4/$10)}' | sort -u > cases_vs_controls.counts.pVal.FDR.LOR.dPSI
### get |dPSI|>=0.5 tested by chi-sqr
cat cases_vs_controls.counts.pVal.FDR.LOR.dPSI | awk -F"\t" '$7<0.05 && ($11>=0.5 || $11<=-0.5)  {print $0}' | sort -u > cases_vs_controls.counts.pVal.FDR.LOR_sig_dPSIabove0.5

#exons not tested by chi-sqr OR tested but not significant
join -j1 -a1 -o1.1,1.8,2.1,2.11 cases_vs_controls.counts.dPSIabove0.5 cases_vs_controls.counts.pVal.FDR.LOR_sig_dPSIabove0.5 | column -t | awk -F" " '$1!=$3 {print $1}' | sort | uniq > cases_vs_controls.counts.dPSIabove0.5_candidates
