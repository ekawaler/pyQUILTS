
cd /ifs/data/proteomics/tcga/databases/refseq_human_20170327
dos2unix *.gpff
dos2unix *.fasta
dos2unix *.txt

perl /ifs/data/proteomics/tcga/databases/../scripts/quilts/v3.0/check_and_clean_fasta.pl .
perl /ifs/data/proteomics/tcga/databases/../scripts/quilts/v2.0/refseq_gpff_to_bed.pl human.protein.gpff refSeqAli.txt
mv refSeqAli.txt.bed transcriptome.bed
mv human.protein.gpff.bed proteome-first.bed
mv human.protein.gpff-descriptions.txt proteome-first-descriptions.txt
mv human.protein.gpff-genes.txt proteome-first-genes.txt
/ifs/data/proteomics/tcga/databases/../scripts/quilts/v3.0/read_chr_bed proteome-first.bed /ifs/data/proteomics/tcga/databases/hg38
perl /ifs/data/proteomics/tcga/databases/../scripts/quilts/v2.0/translate.pl proteome-first.bed.dna
mv proteome-first.bed.dna.fasta proteome-first.fasta 
rm -f proteome-first.bed-mod.fasta
rm -f proteome-first.bed-mod.bed
rm -f proteome-first.bed-mod.log
rm -f proteome-first.bed-mod.stat
rm -f proteome-first.bed-frame-shift.fasta
perl /ifs/data/proteomics/tcga/databases/../scripts/quilts/v3.0/compare_seq_fasta.pl proteome-first.fasta human.protein.fasta.cleaned.fasta
mv proteome-first-corrected.bed proteome.bed
cp proteome-first-descriptions.txt proteome-descriptions.txt
cp proteome-first-genes.txt proteome-genes.txt
/ifs/data/proteomics/tcga/databases/../scripts/quilts/v3.0/read_chr_bed proteome.bed /ifs/data/proteomics/tcga/databases/hg38
perl /ifs/data/proteomics/tcga/databases/../scripts/quilts/v2.0/translate.pl proteome.bed.dna
mv proteome.bed.dna.fasta proteome.fasta 
rm -f proteome.bed-mod.fasta
rm -f proteome.bed-mod.bed
rm -f proteome.bed-mod.log
rm -f proteome.bed-mod.stat
rm -f proteome.bed-frame-shift.fasta
perl /ifs/data/proteomics/tcga/databases/../scripts/quilts/v3.0/compare_seq_fasta.pl proteome.fasta human.protein.fasta.cleaned.fasta
cp proteome.fasta proteome_original.fasta
cat human.protein.fasta.cleaned.fasta.extra.fasta >> proteome_original.fasta
python /ifs/data/proteomics/tcga/databases/../scripts/quilts/v2.0/../../pgx/v0.9/pgx_index.py /ifs/data/proteomics/tcga/databases/refseq_human_20170327

mkdir /ifs/data/proteomics/tcga/databases/refseq_human_original_20170327
cd /ifs/data/proteomics/tcga/databases/refseq_human_original_20170327
cp /ifs/data/proteomics/tcga/databases/refseq_human_20170327/proteome_original.fasta proteome.fasta
dos2unix proteome.fasta
python /ifs/data/proteomics/tcga/databases/../scripts/quilts/v2.0/../../pgx/v0.9/pgx_index.py /ifs/data/proteomics/tcga/databases/refseq_human_original_20170327
cd /ifs/data/proteomics/tcga/databases/refseq_human_20170327
			