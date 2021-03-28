#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 15:01:37 2021
@author: axel
to visualise contact map and intermediars, loops and a chip seq 
"""
import scipy
from scipy.sparse import *
import hicstuff as hcs
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import cooler 
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
from scipy.signal import find_peaks
import random 

bin_matrice=2000

# contact data
chr1="chr7"

# npz format: 
bank="/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a91/cohesin_kochland_paper_data/all_cool_files_renamed/mitotic_dump/chr5-chr5_01_detrended.npz"
matraw = scipy.sparse.load_npz(bank)
matraw = matraw + np.transpose(matraw)
matraw = matraw.toarray()
matscn = matraw

# cool format
cool_file= "/home/axel/Bureau/micro-hackathon/all_cool_files_renamed_all_cool_files_koshland_res_2000_out3_SRR11893085_tmp_valid_idx_pcrfree.pairs.2000.cool"
c = cooler.Cooler(cool_file)
cooler.balance_cooler(c)   # Normalisation 

# loops file: 
df=pd.read_table('/home/axel/Bureau/micro-hackathon/SRR11893084and85_default.tsv',header=None, delimiter="\t", skiprows=1 )

# chip-seq data
chip_ip = pd.read_csv('/home/axel/Bureau/micro-hackathon/SRR2065097.fastq.sam.MQ30.FigS5_Scc1PK9_IP_G1_releasing_60min', 
                    sep=" ", header= None)
chip_input = pd.read_csv('/home/axel/Bureau/micro-hackathon/SRR2065092.fastq.sam.MQ30.FigS5_Scc1PK9_WCE_G1_releasing_60min', 
                    sep=" ", header= None)
name_prot="Cohesin"

# output files 
f_out = open("positive_seq.txt","w+")
f_out2 = open("negative_seq.txt","w+")

fh = open("positive_seq.fasta","w+")
fh2 = open("negative_seq.fasta","w+")

list_all_chrms= ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
                 "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16")

j=0
for chr1 in list_all_chrms :
# computataions 
    matscn = c.matrix().fetch(chr1, chr1)
    matscn[np.isnan(matscn)] = 0 
    coverage = matscn.sum(axis=0)
    coverage[coverage==0] = np.nan
    
    maxi= matscn.shape[0]*bin_matrice
    bin_histo = 100
    chip_ip_chr = chip_ip[chip_ip[0] ==chr1]
    chip_input_chr = chip_input[chip_input[0] ==chr1]
    
    pos_set = df.loc[(df[0] == chr1)]
    pos_set = np.array(pos_set)
    
    v_ip, b_ip = np.histogram(chip_ip_chr[1],
                              bins=range(0,maxi+bin_histo,bin_histo), density=True)
    v_input, b_input = np.histogram(chip_input_chr[1],
                            bins=range(0,maxi+bin_histo,bin_histo), density=True)
    
    b_ip = b_ip/bin_matrice
    b_ip = b_ip[:len(b_ip)-1]
    b_input = b_input/bin_matrice
    b_input = b_input[:len(b_input)-1]
    
    values_chip = v_ip / v_input
    #values_chip = v_ip  # for RNAseq or Exo ChIP
    #values_chip = np.log(values_chip)  # for RNAseq or other 
    
    # fin dpeak in chip-seq:
    #peak_chip = find_peaks(values_chip, 
    #                       height = 1.0,
    #                       width = 2.,
    #                       distance=50)  # for 200 bp resolution 
    
    peak_chip = find_peaks(values_chip, 
                           height = 1.5,
                           width = 1.,
                           distance=100)  # for 100 bp resolution 
    
    h=peak_chip[1]
    left_ips = h['left_ips']
    right_ips = h['right_ips']
    center_ips = (left_ips + right_ips)/2.0 
    height = h['peak_heights']
    
    # Plot: 
    j+=1
    plt.figure(j)
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(8, 11.8)
    gs = gridspec.GridSpec(3, 1,height_ratios=[9,1,1])
    
    ax1 = plt.subplot(gs[0])
    ax1.imshow( matscn**0.15,
               interpolation="none",
               cmap="afmhot_r",
               vmin=0.0, vmax=0.8,
               aspect="auto")
    plt.title(chr1+" ")
    ns =0
    list_loop_chr=[]
    for i in range(pos_set.shape[0] ) :
        site1 = int(pos_set[i,1])
        site2 = int(pos_set[i,4] )
        site1 = int(site1 / bin_matrice)
        site2 = int(site2 / bin_matrice)
        plt.scatter(site1, site2, s=20, facecolors='none', edgecolors='yellow')
        ns +=1
        pi =0
        list_loop_chr.append(site1)
        list_loop_chr.append(site2) 
    ax2 = plt.subplot(gs[1], sharex=ax1)
    ax2.plot(coverage,color='green')
    plt.title("Proportion of the signal in INTRA (in %)")
    
    ax4 = plt.subplot(gs[2], sharex=ax1)
    ax4.plot(b_ip, values_chip, color="royalblue")
    plt.xlabel("Position along the chromosome (bins"+ str(bin_histo) +"bp)")
    plt.title(name_prot+" chip-seq signal")
    
    plt.plot(left_ips*bin_histo/bin_matrice,left_ips*bin_histo/bin_matrice*0.0,'|')
    plt.plot(right_ips*bin_histo/bin_matrice,right_ips*bin_histo/bin_matrice*0.0,'|')
    plt.plot(center_ips*bin_histo/bin_matrice,height,'o',color="orange")
    
    center_ips2 = [int(c*bin_histo/bin_matrice) for c in center_ips]
    common = list(set(list_loop_chr).intersection(center_ips2))
    indices = [center_ips2.index(x) for x in common]
    center_ips_loop = center_ips[indices] 
    plt.plot(center_ips_loop*bin_histo/bin_matrice,np.array(center_ips_loop)*0.,'^',color="red")
    
    list_all = range(matscn.shape[0])
    negatif_group = list(set(list_all) - set(list_loop_chr))
    negatif_group = list(set(negatif_group) - set(center_ips2))
    negatif_group = list(set(negatif_group) - set([0,matscn.shape[0]]))
    negatif_group =random.choices(negatif_group, k=len(center_ips_loop))
    
    plt.plot(negatif_group,np.array(negatif_group)*0.,'s',color="yellow")
    
    
    # writting of genomic positions
    area = 100 # nb of base at left and right of the center 
    
    fasta_file="/home/axel/Bureau/YEAST/agnes_test/sacCer3_with_plasmid_2micron/"+chr1+".fa"
    record = SeqIO.read(open(fasta_file), "fasta")  
    
    for p in center_ips_loop :  #  to write positions of peaks
        start=int(p*bin_histo-area)
        end=int(p*bin_histo+area)
        f_out.write(chr1 +'\t'+ str(start) + '\t' + str(end) + '\n')
        
        fh.write(">loop "+str(chr1)+" "+ str(start) + " " + str(end) + "\n" )
        sequence = record.seq[start: end]
        fh.write(str(sequence)+"\n")
    
    for p in negatif_group :  #  to write positions of peaks
        start=int(p*bin_matrice-area)
        end=int(p*bin_matrice+area)
        f_out2.write(chr1 +'\t'+ str(start) + '\t' + str(end) + '\n')
        
        fh2.write(">neg_loop "+str(chr1)+" "+ str(start) + " " + str(end) + "\n" )
        sequence = record.seq[start: end]
        fh2.write(str(sequence)+"\n")
 
# close of files
pdf = matplotlib.backends.backend_pdf.PdfPages("all.pdf")
for fig in range(1, 16): 
    pdf.savefig( fig )
pdf.close()
plt.close('all')     
        
f_out.close()
f_out2.close()

fh.close()
fh2.close()


