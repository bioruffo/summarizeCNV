#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 11:05:23 2021

@author: roberto

Summarize CNV data from CNVfinder
"""

import os
import openpyxl
from plotly import graph_objects as go
import copy
import glob

def to_tsv(path=None, outdir='out'):
    '''
    Extract data from CNVfinder's XLSX output and save as TSV, kinda like CNVfinder already does
    '''
    for filename in os.listdir(path=path):
        if filename.endswith('.xlsx'):
            print(filename)
            basename = filename.split('.')[0]
            workbook = openpyxl.load_workbook(filename=filename)
            # Assume that the sheet is named `<samplename>_cnvs-annotated` 
            sheet = workbook[basename+'_cnvs-annotated']
            with open(os.path.join(outdir, basename+"_cnvs-annotated.tsv"), "w") as f:
                for row in sheet.rows:
                    f.write('\t'.join(str(x.value) for x in row)+'\n')



def do_sum(path=None, recursive=True, bonferroni=0.05, minscore=0.5):
    '''
    Resume all data by chromosome and by gain/loss
    '''
    all_chromosomes = {'gain': dict(), 'loss': dict()}
    '''
    Final structure will be:
    all_chromosomes = {'gain': {
                                'chr1': {
                                         12522: [2, ['sample1', 'sample2']], ...
                                         },
                                'chr2': {} ...
                                },
                       'loss': <...>
                       }
        
    '''
    for filename in glob.glob(os.path.join(path, '*_cnvs-annotated.tsv'), recursive=recursive):
        if filename.endswith('_cnvs-annotated.tsv'):
            print(filename)
            basename = os.path.basename(filename).split('_cnvs-')[0]
            
            for line in open(filename, "r"):
                # Ignore the header
                if line.startswith('cnv id'):
                    continue
                else:
                    line = line.split('\t')
                    cnvname, chrom, start, end, *_ = line
                    # Need to pass through float because apparently sometimes it's represented as float
                    start = int(float(start))
                    end = int(float(end))
                    gainloss = line[15]
                    bonf = float(line[17])
                    score = float(line[11])
                    chromosomes = all_chromosomes[gainloss]
                    # Only pass calls with Bonferroni below threshold and score of at least threshold
                    if bonf < bonferroni and score >= minscore:
                        if chrom not in chromosomes:
                            chromosomes[chrom] = dict()
                            chromosomes[chrom][start] = [1, [basename]]
                            chromosomes[chrom][end+1] = [0, []]
                        else:
                            curr = sorted(list(chromosomes[chrom].keys()))
                            pre = None
                            first = None
                            last = None
                            for item in curr:
                                if item < start:
                                    pre = item
                                elif item >= start and item <= end:
                                    chromosomes[chrom][item][0] += 1
                                    chromosomes[chrom][item][1].append(basename)
                                    if first is None:
                                        first = item
                                    last = item
                            if start not in chromosomes[chrom]:
                                if pre is None:
                                    chromosomes[chrom][start] = [1, [basename]]
                                else:
                                    chromosomes[chrom][start] = [chromosomes[chrom][pre][0] + 1,
                                                                 chromosomes[chrom][pre][1] + [basename]]
                            if end+1 not in chromosomes[chrom]:
                                if pre is None or last is None:
                                    chromosomes[chrom][end+1] = [0, []]
                                else:
                                    chromosomes[chrom][end+1] = [chromosomes[chrom][last][0] - 1,
                                                                 chromosomes[chrom][last][1].copy()]
                                    chromosomes[chrom][end+1][1].remove(basename)
    # Prettify the data
    for gainloss in all_chromosomes:
        for chrom in all_chromosomes[gainloss]:
            # Set a base 0-point
            if 0 not in all_chromosomes[gainloss][chrom]:
                all_chromosomes[gainloss][chrom][0] = [0, []]
                last = 0
            # Add endpoints for each region
            for position, data in sorted(all_chromosomes[gainloss][chrom].items()):
                if position-1 not in all_chromosomes[gainloss][chrom]:
                    all_chromosomes[gainloss][chrom][position-1] = copy.deepcopy(all_chromosomes[gainloss][chrom][last])
                last = position
                
    return all_chromosomes
                            

def write_gainloss(data, outdir='out'): 
    '''
    Write on disk the data we generated
    '''                   
    for gainloss in data:
        for chromosome in data[gainloss]:
            with open(gainloss+'_'+chromosome+'.tsv', 'w') as f:
                for item in sorted(data[gainloss][chromosome].keys()):
                    f.write('{}\t{}\t{}\n'.format(item, data[gainloss][chromosome][item][0], ','.join(data[gainloss][chromosome][item][1])))
                    if data[gainloss][chromosome][item][0] > 6:
                        print('{}\t{}\t{}\t{}'.format(gainloss, chromosome, item, data[gainloss][chromosome][item][0]))            


def plot_figures(data, outdir='out', extratitletext=''):
    '''
    Plot data both as image and html
    '''
    global all_data, chrom, gainloss
    all_data = {}
    for gainloss in data:
        for chrom in data[gainloss]:
            if not chrom in all_data:
                all_data[chrom] = {'gain': dict(), 'loss': dict()}
            sorted_data = sorted(list(data[gainloss][chrom].items()), key = lambda x: x[0])
            all_data[chrom][gainloss]['x'] = [n[0] for n in sorted_data]
            all_data[chrom][gainloss]['y'] = [n[1][0]*[-1, 1][gainloss=='gain'] for n in sorted_data]
            all_data[chrom][gainloss]['samples'] = [','.join(sorted(['', n[1][1]][len(n[1][1])>0])) for n in sorted_data]
            
            
    for chrom in all_data:
        fig = go.Figure()
        fig.update_layout(title='{} ({})'.format(chrom, extratitletext))
        fig.update_yaxes(title='Cases with gain/loss')
        for gainloss in ['gain', 'loss']:
            if len(all_data[chrom][gainloss]) > 0:
                fig.add_trace(go.Scatter(
                        x = all_data[chrom][gainloss]['x'],
                        y = all_data[chrom][gainloss]['y'],
                        mode = "lines",
                        name = gainloss,
                        line=dict(color=['red', 'green'][gainloss=='loss']),
                        customdata = all_data[chrom][gainloss]['samples'],
                        hovertemplate='<br>%{customdata}</br>'
                 )
                )
        #fig.update_yaxes(range=[-24, 24])
        fig.write_image(os.path.join(outdir, chrom+'.jpeg'))
        fig.write_html(os.path.join(outdir, chrom+'.html'))
        
        
        


if __name__ == '__main__':
    outdir = 'out'
    bonferroni = 0.05
    minscore = 0.5
    
    #print('Extracting data from .xlsx')
    #to_tsv(path=None, outdir=outdir)
    print('Summarizing data from TSVs')
    all_chromosomes = do_sum(path=outdir, recursive=True, bonferroni=bonferroni, minscore=minscore)
    #print('writing count files')
    #write_gainloss(data=all_chromosomes, outdir=outdir)
    print('printing figures')
    plot_figures(data=all_chromosomes, outdir=outdir, 
                 extratitletext="Bonferroni < {}, score >= {}".format(bonferroni, minscore))
    
    
