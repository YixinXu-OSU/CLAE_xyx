#! /usr/bin/env python

#packages that come with python
import sys
import argparse
import os
import pkg_resources
import subprocess
#packages used for various purposes
import numpy as np; import pandas as pd; import matplotlib.pyplot as plt
#packages needed from biopython
from Bio import AlignIO; from Bio.Align import AlignInfo; from Bio.Align.Applications import MuscleCommandline;
from Bio.Align import MultipleSeqAlignment; from Bio.Seq import Seq; from Bio.SeqRecord import SeqRecord
#packages dealing with runtime of the algorithm
from tqdm import tqdm; from halo import Halo
#other
from io import StringIO; import collections

def long_con_start_end_trim(trial):
    #beginning
    if trial[0] == '-':
        while trial[0] == '-':
            trial = trial[1:]
    trial = trial[len(trial)::-1] #flips str
    #end
    if trial[0] == '-':
        while trial[0] == '-':
            trial = trial[1:]
    trial = trial[len(trial)::-1]
    return trial
    
def muscle_generation(file, thresholds=.6809677419354837):

    #align file using muscle msa
    spinner = Halo(text = "Aligning Sequences", spinner = {"interval": 70,"frames": ["_","_","_","-","`","`","'","´","-","_","_","_"]}, color='green')
    spinner.start()
    muscle_cline = MuscleCommandline(os.getcwd() + "/muscle",input= file,clw = True)
    child = subprocess.Popen(str(muscle_cline), 
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                             universal_newlines=True, shell=(sys.platform!="win32"))
    stdout, stderr = muscle_cline()
    aligned_file = AlignIO.read(child.stdout,'clustal')
    spinner.stop()
    print(aligned_file)

    con = []
    counter = 0
    #gets nt info at each index for the alinged sequnces
    all_info = {}
    sequences_column = []
    for i in range(len(aligned_file[0])):
        column = []
        for record in aligned_file:
            record = record.seq
            column.append(record[i])    
        sequences_column.append(column)
        
    #determining threshold for nt/dash
    for threshold in thresholds: #each threshold
        threshold = float(threshold)
        con_str = ''
        for i in range(len(sequences_column)): #go through each index for each threshold
            column = sequences_column[i]
            prop_dash = (collections.Counter(column)['-'])/len(sequences_column[0]) #tells how often a dash shows up
                
            #nt or dash
            #above threshold, a dash in consensus
            if prop_dash > threshold:
                con_str = con_str + '-'
                #just most frequent nt 
            else:
                freq = collections.Counter(column).most_common()
                freq = list(filter(lambda x: x[0] != '-', freq))
                con_str = con_str + str(freq[0][0])
        con.append(con_str)
        
        
    for c in con:
        noDash = c.replace('-','')
        write = open("Consensus_" + os.path.basename(file),'w')
        writeNoDash = open("Consensus_No_Dashes_" + os.path.basename(file),'w')
        
        #write to all dashes
        write.writelines(">Consensus_" + os.path.basename(file) + "\n")
        write.writelines(long_con_start_end_trim(c))
        
        #write to no dashes
        writeNoDash.writelines(">Consensus_No_Dashes_" + os.path.basename(file) + "\n")
        writeNoDash.writelines(noDash)
        
        write.close()
        writeNoDash.close()

        #add consensus to msa alignment
        new_align = []
        new_align.append(SeqRecord(Seq(con[0]), id="CON"))
        for a in aligned_file:
            new_align.append(a)
        new_align = MultipleSeqAlignment(new_align)
        AlignIO.write(new_align,"MSA_Consensus_" + os.path.basename(file),"clustal")
        
    print("Consensus Sequence Generated")

    
    