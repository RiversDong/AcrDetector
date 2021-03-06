#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
AcrDetector: A tool to determine anti-CRISPR proteins in the whole genome-scale

Version 1.0
Author: Chuand Dong 
Email: chuand@uchicago.edu and chuand.dong@gmail.com

Please contact me if you find bugs or any other questions when using the script
'''

import pandas as pd
import numpy as np
from joblib import load
import sys,os
from Bio import SeqIO
import getopt
import re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

script_path = os.path.dirname(os.path.realpath(__file__))
modelHthdb=script_path + "/modelHthDb"

def optParse(argv):
    if len(argv)==0:
        print("-i and --infile are same. Specify the <inputfile>\n-o and --outfile are same.  Specify the <outputfile>")
        sys.exit()
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["help", "infile=", "outfile="])
    except getopt.GetoptError:
        print("main.py -i <inputfile> -o <outputfile> OR main.py -infile <inputfile> -outfile <outputfile>")
        sys.exit()
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("-i and --infile are same. Specify the <inputfile>\n-o and --outfile are same.  Specify the <outputfile>")
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile=arg
        elif opt in ("-o", "--outfile"):
            outfile=arg
    return infile, outfile

def codonFreq(target_seq):
    code_l = re.findall(r'(\w\w\w)',str(target_seq))
    sum_num = len(code_l)
    CodonsDict = dict([('TTT', 0), ('TTC', 0), ('TTA', 0), ('TTG', 0), ('CTT', 0), ('CTC', 0), ('CTA', 0),
                     ('CTG', 0), ('ATT', 0), ('ATC', 0),
                     ('ATA', 0), ('ATG', 0), ('GTT', 0), ('GTC', 0), ('GTA', 0), ('GTG', 0), ('TAT', 0),
                     ('TAC', 0), ('TAA', 0), ('TAG', 0),
                     ('CAT', 0), ('CAC', 0), ('CAA', 0), ('CAG', 0), ('AAT', 0), ('AAC', 0), ('AAA', 0),
                     ('AAG', 0), ('GAT', 0), ('GAC', 0),
                     ('GAA', 0), ('GAG', 0), ('TCT', 0), ('TCC', 0), ('TCA', 0), ('TCG', 0), ('CCT', 0),
                     ('CCC', 0), ('CCA', 0), ('CCG', 0),
                     ('ACT', 0), ('ACC', 0), ('ACA', 0), ('ACG', 0), ('GCT', 0), ('GCC', 0), ('GCA', 0),
                     ('GCG', 0), ('TGT', 0), ('TGC', 0),
                     ('TGA', 0), ('TGG', 0), ('CGT', 0), ('CGC', 0), ('CGA', 0), ('CGG', 0), ('AGT', 0),
                     ('AGC', 0), ('AGA', 0), ('AGG', 0),
                     ('GGT', 0), ('GGC', 0), ('GGA', 0), ('GGG', 0)])
    for item in CodonsDict.keys():
        if item in code_l:
            CodonsDict[item] = code_l.count(item)
    tensor = []
    for i in sorted(CodonsDict.keys()):
        tensor.append(float('%.4f' % (CodonsDict[i]/float(sum_num))))
    return np.array(tensor)

def codonDistance(gene2seq):
    genome_cds=''
    gene2distance={}
    proteins_in_genomes=gene2seq.keys()
    for i in proteins_in_genomes:
        sequence=gene2seq[i]
        genome_cds=genome_cds+sequence
    genome_codon = codonFreq(genome_cds)
    for i in proteins_in_genomes:
        sequence=gene2seq[i]
        i_codon=codonFreq(sequence)
        distance=np.linalg.norm(genome_codon-i_codon)
        gene2distance[i]=distance
    return gene2distance

def geneDownstreamHth(geneIDs, gene2strand, gene2chr, chr2gene, hthSearchResult):
    f=open(hthSearchResult).read().split("\n")
    f=f[3:-11]
    query2hth={}
    gene2hth={}
    gene2downstream={}
    del_=[]
    for i in f:
        i_info=i.split()
        query=i_info[2]
        accession=i_info[1]
        if query not in query2hth.keys():
            query2hth[query]=[accession]
        else:
            query2hth[query].append(accession)
    geneWithHth=query2hth.keys()
    for i in geneIDs:
        if gene2strand[i]=="+":
            chromosome=gene2chr[i]
            genes=chr2gene[chromosome]
            i_index=genes.index(i)
            neighbour_index=i_index+1
            downstream=genes[neighbour_index:neighbour_index+3]
            gene2downstream[i]=downstream
        else:
            chromosome=gene2chr[i]
            genes=chr2gene[chromosome]
            i_index=genes.index(i)
            neighbour_index=i_index-3
            if neighbour_index>=0:
                downstream=genes[neighbour_index:i_index]
            else:
                downstream=genes[0:i_index]
            gene2downstream[i]=downstream
        if downstream:
            tmp=[]
            for j in downstream:
                if j in geneWithHth:
                    tmp.append("1")
                else:
                    tmp.append(1)
            if "1" in tmp:
                gene2hth[i]=1
            else:
                gene2hth[i]=0
        else:
            gene2hth[i]=0
    return gene2hth

def deviation(gene2distance):
    gene2dev={}
    genes=gene2distance.keys()
    comparisions=len(genes)-1
    for i in genes:
        large_num=0
        for j in genes:
            if gene2distance[i]>gene2distance[j]:
                large_num=large_num+1
        dev=large_num/comparisions
        gene2dev[i]=dev
    return gene2dev

def extractFeatures(infile, outfile):
    OUT=open(outfile,"w")
    proteinOutFile=outfile+"_protein"
    outProSeq=open(proteinOutFile, "w")
    outHth=outfile+"_tbl"
    outHthDomain=outfile+"_domain"
    records=SeqIO.parse(infile, "fasta")
    gene2seq={}
    gene2len={}
    gene2chr={}
    gene2func={}
    gene2strand={}
    chr2gene={}
    geneIDs=[]
    gene2stuation={}
    for i in records:
        header=str(i.description)
        tmp_header=header
        seq=str(i.seq)
        if len(seq)%3==0:
            if re.search("[^ATGC]",seq):
                gene2stuation[str(i.id)]="SBS"
                pass
            else:
                if "complement" in header:
                    strand="-"
                else:
                    strand="+"
                header=header.split(" [")
                seqID=header[0]
                tmp_info=header[0].split("|")[1].split("_")
                chromosome=tmp_info[0]
                protein=seqID
                geneIDs.append(seqID)
                if "[protein=" in tmp_header:
                    proteinFunc=re.search(r"\[protein=.*\]*?",tmp_header)
                    proteinFunc=proteinFunc.group(0).split(" [")[0]
                else:
                    proteinFunc="unknown"
                gene2chr[protein]=chromosome
                gene2seq[protein]=seq
                if proteinFunc in ["[protein=hypothetical protein]","[protein=Uncharacterised protein]","[protein=putative uncharacterized protein]"]:
                    gene2func[protein]=1
                else:
                    gene2func[protein]=0
                gene2len[protein]=len(seq)
                gene2strand[protein]=strand
                if chromosome not in chr2gene.keys():
                    chr2gene[chromosome]=[protein]
                else:
                    chr2gene[chromosome].append(protein)
                coding_dna=Seq(seq, IUPAC.unambiguous_dna)
                proteinSeq=str(coding_dna.translate()).replace("*","")
                outProSeq.write(">{0}\n{1}\n".format(seqID, proteinSeq))
        else:
            gene2stuation[str(i.id)]="!3"
    outProSeq.close()
    gene2distance=codonDistance(gene2seq)
    hmmscan_cmd="hmmscan -o {0}  --tblout {1} --noali \
            -E 1e-5 {2}/hthpfam {3}".format(outHthDomain, outHth, modelHthdb, proteinOutFile)
    os.system(hmmscan_cmd)
    gene2hth=geneDownstreamHth(geneIDs, gene2strand, gene2chr, chr2gene, outHth)
    gene2dev=deviation(gene2distance)
    features=[]
    gene2feature={}
    for i in geneIDs:
        chromosome_i=gene2chr[i]
        chromosome_genes=chr2gene[chromosome_i]
        i_index=chromosome_genes.index(i)
        chrGenesNum=len(chromosome_genes)-1
        if i_index==0:
            if len(chromosome_genes)>=3:
                geneI5=chromosome_genes[1]
                geneI3=chromosome_genes[2]
        elif i_index==chrGenesNum:
            geneI5=chromosome_genes[i_index-1]
            geneI3=chromosome_genes[i_index-2]
        else:
            geneI5=chromosome_genes[i_index-1]
            geneI3=chromosome_genes[i_index+1]
        func=gene2func[i]
        dev_i5=gene2dev[geneI5]
        dev_i3=gene2dev[geneI3]
        len_i5=gene2len[geneI5]
        len_i3=gene2len[geneI3]
        strand_i=gene2strand[i]
        if strand_i==gene2strand[geneI5] or strand_i==gene2strand[geneI3]:
            codirection=1
        else:
            codirection=0
        i_feature=[codirection,len_i5,gene2len[i],len_i3, gene2func[i],\
                gene2distance[i],dev_i5, gene2dev[i],dev_i3, gene2hth[i]]
        features.append(i_feature)
        str_feature=map(str,i_feature)
        gene2feature[i]="\t".join(str_feature)
    features=np.array(features)
    rf = load('{0}/rf.joblib'.format(modelHthdb))
    y_prediction = rf.predict(features)
    y_prediction_probability=rf.predict_proba(features)
    y_prediction_probability=list(y_prediction_probability)
    probability_range=range(0,len(y_prediction_probability))
    acr_rank_tmp=[]
    for i_index in probability_range:
        acr_id=geneIDs[i_index]
        probability=y_prediction_probability[i_index][1]
        if probability>0:
            acr_rank_tmp.append({"acr_id":acr_id,"probability":probability})
    acr_rank=sorted(acr_rank_tmp,key=lambda x:x['probability'], reverse=True)
    title_info=['Accession', 'Codirection', 'LNL', 'TL', 'RNL',\
            'Function', 'codon distance', 'LDev', 'TDev', 'RDev', 'HTH', 'Score']
    OUT.write("\t".join(title_info)+"\n")
    if acr_rank_tmp:
        for i in acr_rank:
            acr_accession=i["acr_id"]
            result=[acr_accession,gene2feature[acr_accession],str(i["probability"])]
            result="\t".join(result)
            OUT.write(result+"\n")
    OUT.close()

if __name__ == "__main__":
    argv=sys.argv[1:]
    infile, outfile=optParse(argv)
    extractFeatures(infile, outfile)

