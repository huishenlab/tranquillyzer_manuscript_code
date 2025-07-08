# Junction bed files

# Path to minimap2 directory
# (EDIT WITH YOUR PATH OR MODIFY IF paftools.js IN PATH)
MINIMAP=/path/to/minimap2-2.24

# Path to reference GTF
# (EDIT WITH YOUR PATH)
GTF=/path/to/hg38_gencode.gtf

${MINIMAP}/misc/paftools.js gff2bed ${GTF} > hg38_gencode.junction.bed

# refFlat files
# Assumes gtfToGenePred in PATH
gtfToGenePred \
    -genePredExt \
    -geneNameAsName2 \
    ${GTF} \
    hg38_gencode.TMP.refFlat

awk -F"\t" 'BEGIN{OFS="\t"} {
    print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
}' hg38_gencode.TMP.refFlat > hg38_gencode.refFlat
