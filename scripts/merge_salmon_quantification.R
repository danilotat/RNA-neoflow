
library("dplyr")
library("tximport")
library("tximeta")
library("tibble")

# This script will convert transcript level quantification into
# gene-level one. A single matrix will be retrieved, using TPM as unit
# Wildcards or list of samples?

# 

files <- file.path(snakemake@input[["quant"]])
patients <- snakemake@params[["patients"]]
coldata <- data.frame(files=files, names=patients)
# # build linked Txome
indexDir <- file.path(snakemake@input[["index"]])
fasta <- file.path(snakemake@input[["cdna_fasta"]])
gtf <- file.path(snakemake@input[["annotation"]])
makeLinkedTxome(indexDir=indexDir,source="Ensembl",
organism="Homo sapiens",release="105",genome="GRCh38",
fasta=fasta,gtf=gtf)
se <- tximeta(coldata, useHub=FALSE)
# transcript_level_TPM <- se@assays@data@listData[["abundance"]]
# write.csv(transcript_level_TPM, file=snakemake@output["transcript"])
# gse <- summarizeToGene(se)
# gene_level_TPM <- gse@assays@data@listData[["abundance"]]
# write.csv(gene_level_TPM, file=snakemake@output["gene"])
#samples = c("ERR3549177","ERR3549175","ERR3549179")
#files <- file.path(paste0("/hpc/scratch/danilo.tatoni/Neoepitopes/Neoepi_output/quantification/", samples, "/quant.sf"))
# build linked Txome
#indexDir <- file.path("/hpc/scratch/danilo.tatoni/Neoepitopes/resources/")
#fasta <- file.path("/hpc/scratch/danilo.tatoni/Neoepitopes/resources/Homo_sapiens.GRCh38.cdna.all.fa")
#gtf <- file.path("/hpc/scratch/danilo.tatoni/Neoepitopes/resources/Homo_sapiens.GRCh38.105.gtf")
#makeLinkedTxome(indexDir=indexDir,source="Ensembl",
#organism="Homo sapiens",release=105,genome="GRCh38",
#fasta=fasta,gtf=gtf,write=FALSE)
#se <- tximeta(files)


# report quantification for genes & transcripts

outdir <- file.path(snakemake@params[["outdir"]])

transcript_level_TPM <- se@assays@data@listData[["abundance"]] %>% as.data.frame() %>% tibble::rownames_to_column(var="transcripts")
write.table(transcript_level_TPM, sep="\t", file=paste0(outdir,"/transcript_expression.tsv"),row.names = F, quote=F)
gse <- summarizeToGene(se, countsFromAbundance="lengthScaledTPM")
gene_level_TPM <- gse@assays@data@listData[["abundance"]] %>% as.data.frame() %>% tibble::rownames_to_column(var="genes")
write.table(gene_level_TPM, sep="\t", file=paste0(outdir,'/gene_expression.tsv'),row.names = F, quote=F)

