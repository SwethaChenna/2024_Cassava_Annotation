# convert bed file of Low_expressed_genes.bed into bed12 format and replace "Nascent" with "SMC_gene" & add _tx_01 to the end of transcript_id and save file with file name Low_expressed_transcripts.bed

sed -i -e 's/SMC_gene/Ctrl_SMC_gene/g' Low_expressed_transcripts_wt_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed
sed -i -e 's/SMC_gene/Cold_SMC_gene/g' Low_expressed_transcripts_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed
sed -i -e 's/HC_gene/Ctrl_HC_gene/g' -e 's/MC_gene/Ctrl_MC_gene/g' -e 's/LC_gene/Ctrl_LC_gene/g' Called_transcripts_wt_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed
sed -i -e 's/HC_gene/Cold_HC_gene/g' -e 's/MC_gene/Cold_MC_gene/g' -e 's/LC_gene/Cold_LC_gene/g' Called_transcripts_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed
cat Called_transcripts_wt_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed Low_expressed_transcripts_wt_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed >> Called_SMC_transcripts_wt_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed
cat Called_transcripts_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed Low_expressed_transcripts_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed >> Called_SMC_transcripts_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed
cat Called_SMC_transcripts_wt_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed Called_SMC_transcripts_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed >> Called_SMC_transcripts_ctrl_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed
bedtools sort -i Called_SMC_transcripts_ctrl_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50.bed > Called_SMC_transcripts_ctrl_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50_sorted.bed

# Add gene_id with ';' delimited 4th column(gene_id;transcript_id) - as required by TAMA and save with _forTAMA.bed suffix

grep -v -E "LC|SMC" Called_SMC_transcripts_ctrl_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50_sorted_forTAMA.bed > Called_SMC_transcripts_ctrl_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50_sorted_HC-MC_forTAMA.bed
grep -E "LC|SMC" Called_SMC_transcripts_ctrl_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50_sorted_forTAMA.bed > Called_SMC_transcripts_ctrl_cold_PASminTPM0.1_qTrimTSS0.5PAS0.8_PAS25_ReadSupport5_skipMinTx50_sorted_LC-SMC_forTAMA.bed

######### TAMA ################
conda activate tama

# TMA Merge
tama_merge.py -f filelist.txt -p tama_merge_mod -d merge_dup -a 100 -z 50 -s Ensembl,TR-HC,TR-LC #default parameters m 10
