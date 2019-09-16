# -*- coding: utf-8 -*- 

import ipyparallel
import math
import collections
import os
import gc
import numpy as np
from collections import Counter

#os.system("ipcluster start -n 8 &")   
#os.system("ipcluster stop  &")

binId = lambda : None
bin_size=-1
chrms = []
chrms_flat = []
chrms_length = {}
    

def init(bin_size_local):
    global bin_size, binId
    global chrms_length
    global chrms_stat, chrms, chrms_flat    

    bin_size = bin_size_local
    binId = lambda x: int(x/bin_size) 
    chrms_length_file_name = "./" + "hg19_chrms_length"
    chrms_length = {}
    chrms=[]
    chrms_flat=[] 
    
    with open(chrms_length_file_name, 'r') as chrms_length_file:
        for i, l in enumerate(chrms_length_file):
            ln = l.rstrip().split("\t")
            chrms_length[ln[0]] = int(ln[1])   
    
    for j in range (0, 5):    
        tmp=[]
        for i in range (1, 5):
            tmp.append('chr' + str(j*4 + i))
        chrms.append(tmp)
    
    tmp=[]   
    tmp.append('chr21')
    tmp.append('chr22')
    tmp.append('chrX')
    tmp.append('chrY')
    chrms.append(tmp)    
    chrms_flat=[chr for tmp in chrms for chr in tmp]
#
#
def load_stat(fl_in):
    global chrms_stat
    chrms_stat = {}
    
    with open(fl_in, 'r') as tab_fl:
        for i, l in enumerate(tab_fl):
            ln = l.rstrip().split("\t")
            if len(ln)>2:                
                chrms_stat[ln[0]] = [chrms_length[ln[0]], binId(chrms_length[ln[0]]), int(ln[2])]
    
    for chr in chrms_flat:
        if chr not in chrms_stat:
            chrms_stat[chr] = [chrms_length[chr], binId(chrms_length[chr]), 0]
#
# bin background
def bg_create_binned_profiles(track_name, annotation_table_file, background_file):
    print annotation_table_file
    print background_file
    print chrms
    chrs = {} 
    reads_count = 0
    annot = open(annotation_table_file, 'r')
    bg = open(background_file, 'w+')    
    bg.write("track type=bedGraph name=\"bg_" \
             + track_name + "\"description=\"bg track " \
             + track_name + " binned\"\n")         

    for chr in chrms_flat:
        chrs[chr] = collections.defaultdict(list)
        
    for line in annot:
        ln = line.rstrip().split("\t")        
        assert len (ln)==24, \
               "Warning: create_binnded_background_from_annot_tables \
               length of the input line is less than 24 elements %r" % ln        
        if ln[11]!="protein_coding":
            continue        
        if ln[7]==ln[13]:
            continue        
        chr = ln[13]      
        bin1 = binId(int(ln[14]))        
        reads_count+=1        
        if bin1 in chrs[chr]:
            chrs[chr][bin1]+=1
        else:
            chrs[chr][bin1] = 1     
 
    for chr in chrms_flat:
        for bn in sorted(chrs[chr].keys()):            
            bg.write(chr + "\t" + str(bn * bin_size) \
                     + "\t" + str(bn * bin_size+bin_size) 
                     + "\t" + str(chrs[chr][bn]) + "\n")    
    bg.close()
    annot.close()

#
# normalize background
def bg_normalize(fl_in, fl_norm_genome_mean_out, 
                        fl_norm_chrm_mean_out, fl_norm_no_out, 
                        fl_bg_stat):
    bg=open(fl_in, "r")
    bg_norm_genome_mean=open(fl_norm_genome_mean_out,'w+')
    bg_norm_chrm_mean=open(fl_norm_chrm_mean_out,'w+')
    bg_norm_no=open(fl_norm_no_out,'w+')
    bg_stat=open(fl_bg_stat, 'w+')
    
    chr=""
    chr_size=0
    total=0.0
    total_truncated=0.0
    total_genome=0.0
    lines=[]
    norm_bg=[]
    number_of_bins=0
    
    chrms_max={}
    chrms_filled={}
    chrms_sum={}
    genome_mean_total=0.0
    genome_sum=0.0
    bins_total=0.0
    bins_filled=0.0
    
    for l in bg:
        ln=l.rstrip().split("\t")  
        if "bedGraph" in l:
            continue                    
        if ln[0] in chrms_sum:
            chrms_sum[ln[0]]+=int(ln[3]) 
            chrms_filled[ln[0]]+=1
        else:
            chrms_sum[ln[0]]=int(ln[3])    
            chrms_filled[ln[0]]=1        
        
        if ln[0] in chrms_length and ln[0] not in chrms_max:
            chrms_max[ln[0]]=chrms_length[ln[0]]
        elif ln[0] not in chrms_length:
            print "Warning: simple_normalization_for_bg chr not in chrm_length dict " + ln[0]   
    
    print 'chr chrm_sum_signal chrms_filled_bins chrm_length chrm_mean' 
    for chr in sorted(chrms_sum.keys()):
        print chr + "\t" + str(chrms_sum[chr]) + "\t" + str (chrms_filled[chr]) \
              + "\t" + str (chrms_max[chr]) + "\t" \
              + str(((chrms_sum[chr] * 1.0)/(1.0 * binId(chrms_max[chr]))))
        bg_stat.write(chr + "\t" + str(chrms_sum[chr]) + "\t" \
                      + str(chrms_filled[chr]) + "\t" + str(chrms_max[chr]) \
                      + "\t" +  str((chrms_sum[chr] * 1.0)/(1.0 * binId(chrms_max[chr]))) + "\n")
        genome_sum+=chrms_sum[chr]
        bins_total+=binId(chrms_max[chr])
        bins_filled+=chrms_filled[chr]
    
    print 'total_bins_total ',
    print bins_total    
    bg_stat.write('bins_total\t' + str(bins_total)  + "\n")
    
    print 'total_bins_filled ',
    print bins_filled
    bg_stat.write('bins_filled\t' + str(bins_filled) + "\n" )
    
    print 'genome_mean bins total ',
    genome_mean_total = (genome_sum * 1.0)/(1.0 * bins_total)
    print genome_mean_total
    bg_stat.write('genome_mean bins total\t' + str(genome_mean_total) + "\n")
    
    print 'genome_mean bins filled ',
    genome_mean_filled = (genome_sum * 1.0)/(1.0 * bins_filled)
    print genome_mean_filled
    bg_stat.write('genome_mean bins filled\t' + str(genome_mean_filled) + "\n") 
    
    bg.close()
    bg = open(fl_in, "r")        
    prev_chrm = ""
    prev_end = 0
    for l in bg:         
        if "bedGraph" in l:
            continue
        ln = l.rstrip().split("\t") 
        bin1 = binId(int(ln[1]))
        bin2 = binId(int(ln[2]))        
        if chr=="" or prev_chrm!=ln[0]:
            chr = ln[0]
            for i in range (1, bin1 + 1):
                bg_norm_genome_mean.write(ln[0] + "\t" + str( bin_size*(i-1)) + "\t" + str(bin_size*(i)) + "\t" + "0\n")
                bg_norm_chrm_mean.write(ln[0] + "\t" + str( bin_size*(i-1)) + "\t" + str(bin_size*(i)) + "\t" + "0\n")
                bg_norm_no.write(ln[0] + "\t" + str( bin_size*(i - 1)) + "\t" + str(bin_size * (i)) + "\t" + "0\n")                
        if  ln[0]==chr :
            for i in range (1, 1 + bin1 - prev_end):
                bg_norm_genome_mean.write(ln[0] + "\t"+ str(prev_end * bin_size + bin_size*(i - 1)) \
                                          + "\t" + str(prev_end*bin_size+ bin_size*(i)) + "\t" + "0\n")
                bg_norm_chrm_mean.write(ln[0] + "\t" + str(prev_end * bin_size + bin_size * (i - 1)) \
                                        + "\t" + str(prev_end * bin_size + bin_size * (i)) + "\t" + "0\n")
                bg_norm_no.write(ln[0] + "\t" + str(prev_end * bin_size + bin_size * (i - 1)) \
                                 + "\t" + str(prev_end*bin_size+bin_size * (i)) + "\t" + "0\n")
            for i in range(1, (1 + bin2 - bin1)):
                bg_norm_genome_mean.write(ln[0] + "\t" + str(bin1 * bin_size + bin_size * (i - 1)) \
                                          + "\t" + str(bin1 * bin_size + bin_size * (i)) \
                                          + "\t" + str(float(ln[3])/genome_mean_total) + "\n")           
                bg_norm_chrm_mean.write(ln[0] + "\t" + str(bin1 * bin_size + bin_size * (i - 1)) \
                                        + "\t"+  str(bin1 * bin_size + bin_size * (i)) \
                                        + "\t" + str(float(ln[3])/((chrms_sum[chr] * 1.0)/(1.0 * binId(chrms_max[chr])))) + "\n") 
                bg_norm_no.write(ln[0] + "\t" + str(bin1 * bin_size + bin_size * (i - 1)) \
                                 + "\t" + str(bin1 * bin_size + bin_size * (i)) + "\t" + str(ln[3]) + "\n") 
        prev_end = bin2
        prev_chrm = ln[0]     
        
    bg.close()
    bg_norm_genome_mean.close()
    bg_norm_chrm_mean.close()
    bg_norm_no.close() 
    bg_stat.close()

#
#
def core_bg_smooth_with_window(lines_list, chr, smoothed):
    lines = np.asarray(lines_list, dtype=object) 
    proximity = [ float(x) for x in lines[0:11,3] ]     
    for j in range(len(lines)):
        if j<=5:
            mn = np.mean(proximity)
        elif j>len(lines) - 5 - 1:
            mn = np.mean(proximity)
        else:    
            del proximity[0]
            proximity.append(float(lines[j+5][3]))    
            mn = np.mean(proximity)    
        smoothed.write(lines[j][0] + "\t" + lines[j][1] + "\t" + lines[j][2] + "\t" +str(mn) + "\n")   
    del lines    
#
# smooth bg with 11-bin window
def bg_smooth_with_window(fl_bg, fl_bg_smoothed):
    smoothed = open(fl_bg_smoothed, 'w+')
    lines_list = []
    chr = ""
    with open(fl_bg, 'r') as bg:
        for i, l in enumerate(bg):
            ln = l.rstrip().split("\t")
            if chr=="":
                chr = ln[0]
                lines_list.append(ln)
            elif ln[0]==chr:
                lines_list.append(ln)
            else:
                core_bg_smooth_with_window(lines_list,chr, smoothed)                 
                chr = ln[0]
                del lines_list                    
                gc.collect()                    
                lines_list = []
                lines_list.append(ln)
        core_bg_smooth_with_window(lines_list,chr, smoothed)            
    smoothed.close()    
#
# filter annotation tabel files
def rna_filter(file_in, file_out, mask, name):
    fl_in = open (file_in, 'r')
    fl_out = open (file_out,'w+')
    if 'region' in mask:
        mk = mask.split("_")        
        chr = mk[1]
        start = int(mk[2])
        end = int(mk[3])
    else:
        mk = []
        start = -1
        end = -1
    print "filter_rna : " ,    
    print "mask = " ,
    print mk
    for ln in fl_in:
        lna = ln.rstrip().split("\t")        
        if lna[11]=='chr19' and int(ln[12])==27738479:            
            continue        
        if mask=='name' and lna[12]==name :
            fl_out.write(ln)
        if mask == 'pc' and 'protein_coding' in ln:
            fl_out.write(ln)
        if mask == 'pc,trans' and 'protein_coding' in ln and lna[7]!=lna[13]:
            fl_out.write(ln)  
        if mask=='all,trans' and lna[7]!=lna[13]:
            fl_out.write(ln)
        if 'region' in mask:
            if lna[0]==chr and int(lna[1])>=start and int(lna[2])<=end:
                fl_out.write(ln)            
    fl_in.close()
    fl_out.close()
#
# calculate rna coverage
def rna_calculate_coverage(fl_annot_tab, fl_out):
    total_coverage = 0
    rna = {}
    with open(fl_annot_tab, 'r') as tab_fl:
        for i, l in enumerate(tab_fl):
            ln = l.rstrip().split("\t")            
            bin1 = binId(int(ln[14]))
            bin2 = binId(int(ln[15]))            
            if bin2-bin1>1: 
                print "Warning: calculate_single_RNA_coverage very long RNA " + l            
            if ln[13] in rna:
                rna[ln[13]][bin1] = rna[ln[13]][bin1]+1 if bin1 in rna[ln[13]] else 1
                total_coverage+=1                
            else:
                rna[ln[13]] = {}
                rna[ln[13]][bin1] = 1
                total_coverage+=1
    print  fl_out ,  
    end = 0
    print ' total RNA coverage ',
    print total_coverage
    with open(fl_out, 'w+') as rna_bin:
        if ".tab" in fl_out:
            track_name = fl_out.split(".tab")[0]
        elif ".bedGraph" in fl_out:
            track_name = fl_out.split(".bedGraph")[0]
        else:
            track_name = ""
            print "Warning: normalize_single_rna problem with track name"         
        rna_bin.write("track type=bedGraph name=\"" + track_name + "\"description=\"" + track_name + "\"\n")   
        for chr in sorted(rna.keys()):
            for bin in sorted(rna[chr].keys()):
                if (bin+1)*bin_size > chrms_stat[chr][0] : 
                    end = (bin+1)*bin_size
                else:
                    end = (bin+1)*bin_size 
                rna_bin.write(chr+"\t"+str(bin*bin_size) +"\t"+ str(end) +"\t"+ str(rna[chr][bin]) + "\n")                
#
# normalize rna coverage
def rna_normalize(fl_rna_binned, fl_rna_norm):
    rna = {}
    total_coverage = []
    with open(fl_rna_binned, 'r') as rna_bin:
        for i, l in enumerate(rna_bin):
            ln = l.rstrip().split("\t")    
            if len(ln) <4:
                print "Warning: normalize_single_rna short line " + l
                continue
            if ln[0] not in rna: rna[ln[0]] = []
            rna[ln[0]].append(ln)
            total_coverage.append(float(ln[-1]))    
    print 'total coverage ',
    print sum(total_coverage)      
    number_of_bins = sum(chrms_stat[chr][1] for chr in chrms_stat)    
    mean_total = sum(total_coverage)/(number_of_bins*1.0)
    
    print 'mean by all bins in the genome',
    print mean_total    
    
    print 'mean by only filled bins ',
    mean_filled = np.mean(total_coverage)    
    print mean_filled
    
    with open(fl_rna_norm, 'w+') as rna_norm:
        if ".tab" in fl_rna_norm:
            track_name = fl_rna_norm.split(".tab")[0]
        elif ".bedGraph" in fl_rna_norm:
            track_name = fl_rna_norm.split(".bedGraph")[0]
        else:
            track_name = ""
            print "Warning: normalize_single_rna problem with track name"  
        rna_norm.write("track type=bedGraph name=\"" + track_name + "\"description=\"" + track_name + "\"\n")  
        for chr in sorted(rna.keys()):
            for ln in rna[chr]:                
                rna_norm.write(chr+"\t"+str(int(ln[1])) +"\t"+ str(int(ln[2])) +"\t"+ str((float(ln[3])*1.0)/mean_total) + "\n")
#
# calculate fold enrichment
def rna_calculate_fold_enrichment(fl_bg, fl_rna, fl_out):
    fold_enrichment =  open(fl_out, 'w+')
    if ".tab" in fl_out:
        track_name = fl_out.split(".tab")[0]
    elif ".bedGraph" in fl_out:
        track_name = fl_out.split(".bedGraph")[0]
    else:
        track_name = ""
        print "Warning: normalize_single_rna problem with track name"      
    fold_enrichment.write("track type=bedGraph name=\"" +track_name + "\"description=\"" + track_name + "\"\n")    
    bg = {}
    eps = 0.0001
    total_coverage = []
    with open(fl_bg, 'r') as f_bg:
        for i, l in enumerate(f_bg):
            ln = l.rstrip().split("\t")         
            if ln[0] not in bg: bg[ln[0]] = {}
            bg[ln[0]][int(ln[1])] = float(ln[-1])
    rna = {}
    count_zeros = 0
    with open(fl_rna, 'r') as f_rna:
        for i, l in enumerate(f_rna):
            ln = l.rstrip().split("\t")    
            if len(ln)<4:
                print "Warning: calculate_fold_enrichment_for_single_RNA short line " + l
                continue            
            bin1 = int(ln[1])
            if ln[0] not in rna: rna[ln[0]] = {}
            if bin1 not in rna[ln[0]] : rna[ln[0]][bin1] = []
            if bin1 not in bg[ln[0]]: 
                print "Warning: calculate_fold_enrichment_for_single_RNA bin not in bg"
                print ln[0] + "\t" + str(bin1) 
                continue
            fold = (float(ln[-1]))/(bg[ln[0]][bin1]) if bg[ln[0]][bin1]!=0  else (-1) * float(ln[-1])
            if bg[ln[0]][bin1]==0:
                count_zeros+=1
            fold_enrichment.write(ln[0] + "\t" + ln[1] + "\t" +ln[2] + "\t" +str(fold) + "\n")
    print 'zeros ',        
    print count_zeros
    fold_enrichment.close()   
#
# filter fold enrichment
def rna_filter_fold_enrichment(fl_RNA, fl_out):
    filtered = open(fl_out,"w+")   
    if ".tab" in fl_out:
        track_name = fl_out.split(".tab")[0]
    elif ".bedGraph" in fl_out:
        track_name = fl_out.split(".bedGraph")[0]
    else:
        track_name = ""
        print "Warning: normalize_single_rna problem with track name"  
    filtered.write("track type=bedGraph name=\"" +track_name + "\"description=\"" + track_name + "\"\n")                                             
    
    list_of_coordinates = {}
    with open(fl_RNA, 'r') as fl_RNA:
        for i, l in enumerate(fl_RNA): 
            ln = l.rstrip().split("\t")
            if len(ln)<4:                 
                print "Warning: filtering_enrichment_signal_for_RNA short line " + l
                continue            
            if ln[0] not in list_of_coordinates: list_of_coordinates[ln[0]] = []
            if abs(float(ln[3]))<2:
                continue
            list_of_coordinates[ln[0]].append(ln)
    for chr in sorted(list_of_coordinates.keys()):
        for i,ln in enumerate(list_of_coordinates[chr]):
            bin1 = binId(int(ln[1]))
            start = i-5 if i>4 else 0
            end = i+6 if i<len(list_of_coordinates[chr])-5 else len(list_of_coordinates[chr])
            count = 0
            for j in range(start,end):
                if abs(binId(int(list_of_coordinates[chr][j][1]))-bin1)<=5:
                    count+=1
            if count >=3:
                filtered.write(ln[0]+"\t"+ ln[1]+"\t"+ ln[2]+"\t"+ ln[3]+"\n")
    filtered.close()
#
#
def core_rna_smooth_with_window(lines_list, chr, smoothed):
    lines = np.asarray(lines_list, dtype=object)  
    checked_bins = []    
    for i, ln in enumerate(lines):
        bin = binId(int(ln[1]))
        proximity1 = [0] * 11
        proximity2 = [0] * 11
        proximity3 = [0] * 11
        for j in range(-5,6):
            if i+j>=0 and i+j<len(lines) and abs(binId(int(lines[i+j][1])) - binId(int(ln[1])))<=5:
                proximity1[5 + binId(int(lines[i + j][1])) - binId(int(ln[1]))] = (float(lines[i + j][3])) 
                proximity2[5 + binId(int(lines[i + j][1])) - binId(int(ln[1]))] = (float(lines[i + j][3])) 
                proximity3[5 + binId(int(lines[i + j][1])) - binId(int(ln[1]))] = (float(lines[i + j][3])) 
        mn2 = np.mean(proximity2)     
        tmp = {}                  
        for k in range(-1,-6, -1): 
            del proximity1[-1]
            proximity1 = [0.0] + proximity1
            mn1 = np.mean(proximity1)
            tmp[k] = mn1            
        for k in range(-5,0, +1):            
            if bin + k not in checked_bins and bin+k>=0:
                smoothed.write(ln[0] + "\t" + str((bin + k) * bin_size) \
                               + "\t" + str((bin + k + 1) * bin_size) + "\t" + str(tmp[k]) + "\n")
                checked_bins.append(bin + k)  
        if bin not in checked_bins:
            smoothed.write(ln[0] + "\t" + ln[1] + "" + "\t" + ln[2] + "\t" + str(mn2) + "\n")     
        checked_bins.append(bin)
        for k in range(1,6):
            del proximity3[0]
            found = 0
            for m in range(0, 11):
                if i + m<len(lines) and bin + 5 + k==binId(int(lines[i + m][1])):
                    proximity3.append(float(lines[i + m][3])) 
                    found = 1   
                    break
            if found==0:
                proximity3.append(0.0)                  
            if bin + k not in checked_bins and bin + k <= binId(chrms_length[ln[0]]):
                mn3 = np.mean(proximity3) 
                tmp_count = 0                
                smoothed.write(ln[0] + "\t" + str((bin+k)*bin_size) + "\t" + str((bin+k+1)*bin_size) + "\t" +str(mn3) + "\n")                
                checked_bins.append(bin+k)   
            if i+1 <len(lines) and bin+k+1 == binId(int(lines[i+1][1])) :
                break            
    del lines   
#
# smooth rna profile with 11-bin window
def rna_smooth_with_window(fl_rna, fl_rna_smoothed):
    smoothed = open(fl_rna_smoothed, 'w+')    
    if ".tab" in fl_rna_smoothed:        
        track_name = fl_rna_smoothed.split("/")[-1].split(".tab")[0]
    elif ".bedGraph" in fl_rna_smoothed:
        track_name = fl_rna_smoothed.split("/")[-1].split(".bedGraph")[0]
    else:
        track_name = ""
        print "Warning: normalize_single_rna problem with track name"  
    smoothed.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_name + "\"\n")     
    lines_list = []
    chr = ""
    checked_bins = []
    with open(fl_rna, 'r') as rna:
        for i, l in enumerate(rna):
            ln = l.rstrip().split("\t")
            if len(ln)<4:
                print "Warning: smooth_with_window_rna short line " + l
                #smoothed.write(l)
                continue
            if chr=="":
                chr = ln[0]
                lines_list.append(ln)                
            elif ln[0]==chr:
                lines_list.append(ln)
            else:
                core_rna_smooth_with_window(lines_list, chr, smoothed)                 
                chr = ln[0]
                del lines_list                    
                gc.collect()                    
                lines_list = []
                lines_list.append(ln)
                checked_bins = []
        if len(lines_list)!=0:
            core_rna_smooth_with_window(lines_list, chr, smoothed)        
    smoothed.close()        
#
# correct last bin in profiles to visualize in genome browser
def correct_last_bin_chmrs_end(fl_in_name, fl_out_name):
    with open(fl_in_name, 'r') as fl_in, \
         open(fl_out_name, 'w+') as fl_out:
        for i, l in enumerate(fl_in):
            ln=l.rstrip().split("\t")
            if "track type=" in l:
                fl_out.write(l)
            if len (ln)<4:
                print "Warning: correct_last_bin_chmrs_end short line " + l
                continue              
            if int(ln[2]) > chrms_length[ln[0]]:
                fl_out.write(ln[0]+ "\t" + ln[1] + "\t" + str(chrms_length[ln[0]]) + "\t" + ln[3] + "\n")
            else:
                fl_out.write(l)    
#
#
def Workflow():    
    init(-1)       
    rnas = { "AXIN1" : {"K562" : [1000000, 100000, 50000, 20000, 10000]}, 
             "MIR3648" : {"K562" : [1000000, 100000, 50000, 20000, 10000]}, 
             "MIR3687" : {"K562" : [1000000, 100000, 50000, 20000, 10000]},
             "GAPDH" : {"K562" : [1000000, 100000, 50000, 20000, 10000]}, 
             "MALAT1" : {"K562" : [1000000, 100000, 50000, 20000, 10000]}, 
             "XIST" : {"K562" : [1000000, 100000, 50000, 20000, 10000], 
                       "fibr" : [1000000, 100000, 50000, 20000, 10000]}
             }      
    bin_sizes = {}
    for rna in rnas:
        for cline in rnas[rna].keys():
            if cline not in bin_sizes:
                bin_sizes[cline] = []
            for bin_size in rnas[rna][cline]:    
                if bin_size not in bin_sizes[cline]:
                    bin_sizes[cline].append(bin_size)   
    #
    # background
    #
    for cline in bin_sizes.keys():
        for bin_size in bin_sizes[cline]:
            print "====================================================" ,
            print "background" ,
            print cline ,
            print bin_size ,
            print "===================================================="  
            
            bin_size_local = bin_size                
            init(bin_size_local)   
            bg_create_binned_profiles(cline + ".background", 
                                      cline + ".full.tab", 
                                      cline + ".background.binned" + str(bin_size_local) + ".tab")    
                                                 
            bg_normalize(cline + ".background.binned" + str(bin_size_local) + ".tab", 
                         cline + ".background.binned" + str(bin_size_local) + ".norm_genome_mean.tab", 
                         cline + ".background.binned" + str(bin_size_local) + ".norm_chrm_mean.tab", 
                         cline + ".background.binned" + str(bin_size_local) + ".norm_no.tab", 
                         cline + ".background.binned" + str(bin_size_local) + ".stat")
                        
            load_stat(cline + ".background.binned" + str(bin_size_local) + ".stat")         
            
            bg_smooth_with_window(cline + ".background.binned" + str(bin_size_local) + ".norm_genome_mean.tab", 
                               cline + ".background.binned" + str(bin_size_local) + ".norm_genome_mean.smoothed.tab")
        
            correct_last_bin_chmrs_end(cline + ".background.binned" + str(bin_size_local) + ".norm_genome_mean.smoothed.tab", 
                                       cline + ".background.binned" + str(bin_size_local) + ".norm_genome_mean.smoothed.corrected.tab") 
    #
    # rnas
    #    
    for rna in rnas:
        print rna
        print rnas[rna]
        for cline in rnas[rna].keys():
            for bin_size in rnas[rna][cline]:                
                bin_size_local = bin_size
                init(bin_size_local) 
                load_stat(cline + ".background.binned" + str(bin_size_local) + ".stat")   
                
                print "====================================================" ,
                print "rna" ,
                print rna ,
                print cline ,
                print bin_size ,
                print "===================================================="  
                
                rna_filter(cline + ".full.tab", 
                           cline + "." + rna + ".tab", 
                           "name", 
                           rna)                                
                rna_calculate_coverage(cline + "." + rna + ".tab", 
                                       cline + "." + rna + ".binned" + str(bin_size_local) + ".tab")   
                
                rna_normalize(cline + "." + rna + ".binned" + str(bin_size_local) + ".tab", 
                              cline + "." + rna + ".binned" + str(bin_size_local) + ".normalized.tab")
                
                rna_calculate_fold_enrichment(cline + ".background.binned" + str(bin_size_local) + ".norm_genome_mean.smoothed.tab", 
                                              cline + "." + rna + ".binned" + str(bin_size_local) + ".normalized.tab", 
                                              cline + "." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.tab")
                
                rna_filter_fold_enrichment(cline + "." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.tab",
                                           cline + "." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.filtered.tab")
                
                rna_smooth_with_window(cline + "." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.filtered.tab", 
                                       cline + "." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.filtered.smoothed.tab")            
                                                           
                correct_last_bin_chmrs_end(cline + "." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.filtered.smoothed.tab", 
                                           cline + "." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.filtered.smoothed.corrected.tab")                 
 
Workflow()

