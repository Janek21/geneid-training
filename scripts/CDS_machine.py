import pandas as pd
import os


def CDS_read(filepath):
    """Given the filepath to a gff file containing CDS, it transforms it to a pandas dataframe and formats"""
    cds_data=pd.read_csv(filepath, sep="\t", names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "geneId"]) #name columns

    cds_data["start"]=pd.to_numeric(cds_data["start"])
    cds_data["end"]=pd.to_numeric(cds_data["end"])

    return cds_data 


def gene_assign(df):
    """Given a dataframe from a gff file containing CDS, it assigns a gene to each group of CDS that belong to a geneId"""
    #get cds only dataframe
    cds_rows=df[df["feature"]=="CDS"]

    #minimum start and end value for each geneId group
    min_start=cds_rows.groupby("geneId")["start"].min().reset_index()
    max_end=cds_rows.groupby("geneId")["end"].max().reset_index()

    startEnd_rows=pd.merge(min_start, max_end, on="geneId")
    
    #get one complete row for each geneId
    gene_rows=cds_rows.drop_duplicates(subset="geneId").reset_index(drop=True)

    #Set feature to "gene"
    gene_rows["feature"]="gene"

    #set start and end of each gene(from the start of the first cds to the end of the last cds of the same geneId)
    gene_rows.drop(columns=["start", "end"], inplace=True)
    gene_rows_complete=pd.merge(gene_rows, startEnd_rows , on="geneId")[["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "geneId"]] #reorder
    
    #concatenate the new gene rows to the CDS dataframe    
    gene_newdata=pd.concat([gene_rows_complete, df]).sort_index(kind="merge").reset_index(drop=True)
    sorted_gene=gene_newdata.sort_values(["geneId", "feature"], ascending=[True, False]).reset_index(drop=True) #sort

    return sorted_gene


def CDS_2_exon(df):
    """Given a dataframe from a gff file containing CDS, it assigns an exon to each CDS"""
    cds_rows=df[df["feature"]=="CDS"].copy() #get duplicate cds_rows in dataframe
    cds_rows["feature"]="exon" #change CDS for exon (they now are exon rows)
    cds_newdata=pd.concat([cds_rows, df]).sort_index(kind="merge").reset_index(drop=True) #insert exon rows together with cds rows (all ordered)
    sorted_cds=cds_newdata.sort_values(["geneId", "feature"], ascending=[True, False]).reset_index(drop=True) #sort
    return sorted_cds


def hierarchy_setter(df):
    """Given a dataframe containing genes, exons and CDS(with a plain ID in the last column), it establishes the proper hierarchy"""

    n_df=df.copy()
    #set all rows (gene, exon, cds) with Parent genId
    n_df["geneId"]="Parent="+n_df["geneId"]
    
    #replace Parent with ID in genes (correct hierarchy)
    is_gene=n_df["feature"]=="gene"
    n_df.loc[is_gene, "geneId"]=n_df.loc[is_gene, "geneId"].str.replace("Parent", "ID")

    return n_df




def CDS_2_gffcomp(inFile, outFile):
    """Takes a CDS only gff file and an output path as an input, and converts the CDS file to gff-comparable"""
    #read file to pandas
    cds_data=CDS_read(inFile)

    #add a gene for each geneId
    gened_data=gene_assign(cds_data)

    #add exons for each cds
    exoned_data=CDS_2_exon(gened_data)

    #set cds>exon>gene hierarchy
    hdata=hierarchy_setter(exoned_data)

    #create folders if they dont exist and write file
    outFolder=outFile[:outFile.rfind("/")]
    os.makedirs(outFolder, exist_ok=True)
    hdata.to_csv(outFile, sep="\t", index=None, header=None)


def main():
    #cds_file="../data/species/Plasmodium_vivax/CDS_Pvivax_ann.gff"
    cds_file="/home/jj/Desktop/Data_science/CRG/TFM/projects/geneid-training/data/species/Drosophila_melanogaster/CDS_Dmelanogaster_ann.gff"
    outPath="../data/"
    CDS_2_gffcomp(cds_file, "../results/file.gff")

#main()
