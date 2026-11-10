# -*- coding: UTF-8 -*-
# Source: from the upstream variant-analysis pipeline scripts used for the
# China 100K Genomes Project cohort (internal deployment).

import json
import re
import os,sys
from lxml import etree
import csv
import pandas as pd
input_file = sys.argv[1]
out_file = sys.argv[2]
sample_name = sys.argv[3]
chrom_list = sys.argv[4].split(",")


def re_genome_results(patterns, line):
    m = re.match(pattern, line)
    if m:
        val = m.group(1)
        return val
   
if __name__ == '__main__':
    data = {}

    file_name = input_file
    text = open(file_name).readlines()
    html = etree.HTML('\n'.join(text))
    base = html.xpath('/html/body/div[2]/div/div[1]/div/div/div/div[2]/div')[0]
    groups = base.xpath('div')
    data[sample_name] = {}
    for group in groups:
        title = group.xpath('h3/text()')[0]
        trs = group.xpath('table/tr')
        headers = ''
        

        for tr in trs:
            if title != 'Chromosome stats':
                tds = tr.xpath('td/text()')
                #data[title][tds[0]] = tds[1]
                data[sample_name][title +" "+tds[0]] = tds[1]
            else:
                if headers == '':
                    headers = tr.xpath('td/b/text()')
                    continue
                tds = tr.xpath('td/text()')
                d = {}
                for i in range(len(tds)):
                    if tds[0] not in chrom_list:
                        break
                    
                    #d[title+"_"+headers[i]] = tds[i]
                    if headers[i] == "Name" or headers[i] == "Length" :
                        continue
                    data[sample_name][tds[0]+"_"+headers[i]] = tds[i]

    filter_file_dir  = os.path.dirname(file_name)
    # print filter_file_dir
    # filter_file = os.path.join(filter_file_dir, "qualimap_filter.txt")
    # with open(filter_file, "r") as filter_f:
    #     for line in filter_f:
    #         line = line.strip()
    #         if line.startswith("[qualimap]"):
    #             continue
    #         if line:
    #             #print line
    #             name, base, val, pass_ = line.strip().split("\t")
    #             data[sample_name][name] = val

    genome_results = os.path.join(filter_file_dir, "genome_results.txt")

    with open(genome_results, "r") as filter_f:
        total_reads_num = None
        sup_reads_num = None
        dup = None
        for line in filter_f:
            line = line.strip()
            
            pattern = r"number of reads = (.*)"
            if not total_reads_num:
                total_reads_num = re_genome_results(pattern, line)
                if total_reads_num:
                    total_reads_num = int("".join(total_reads_num.split(",")))
            pattern = r"number of supplementary alignments = (.*)\("
            if not sup_reads_num:
                sup_reads_num = re_genome_results(pattern, line)
            if total_reads_num and sup_reads_num :
                reads_num =  total_reads_num  -  int("".join(sup_reads_num.split(",")))
                data[sample_name]["Total_reads_num_in_bam"] = reads_num

            pattern = r"mean coverageData = (.*)X"
            Coverage_depth_after_deduplication = re_genome_results(pattern, line)
            if Coverage_depth_after_deduplication:
                data[sample_name]["Coverage_depth_after_deduplication"] = Coverage_depth_after_deduplication

            pattern = r"There is a (.*)\% of reference with a coverageData >= 10X"
            Coverage_depth_on_10x = re_genome_results(pattern, line)
            if Coverage_depth_on_10x:
                data[sample_name]["Coverage_depth_on_10x"] = Coverage_depth_on_10x

            pattern = r"There is a (.*)\% of reference with a coverageData >= 1X"
            Coverage_depth_on_1x = re_genome_results(pattern, line)
            if Coverage_depth_on_1x:
                data[sample_name]["Coverage_depth_on_1x"] = Coverage_depth_on_1x

            pattern = r"There is a (.*)\% of reference with a coverageData >= 20X"
            Coverage_depth_on_20x = re_genome_results(pattern, line)
            if Coverage_depth_on_20x:
                data[sample_name]["Coverage_depth_on_20x"] = Coverage_depth_on_20x

            pattern = r"number of duplicated reads \(flagged\) = (.*)"
            if not dup:
                dup = re_genome_results(pattern, line)

            if dup and total_reads_num:
                data[sample_name]["Duplication_rate"] = float("".join(dup.split(",")))/total_reads_num

            pattern = r"GC percentage = (.*)\%"
            Gc_content = re_genome_results(pattern, line)
            if Gc_content:
                data[sample_name]["Gc_content"] = Gc_content

            pattern = r"number of mapped reads = .*?\((.*)\%\)"
            Mapping_rate = re_genome_results(pattern, line)
            if Mapping_rate:
                data[sample_name]["Mapping_rate"] = Mapping_rate

        data[sample_name]["Average_read_length"] = data[sample_name]["Globals Read min/max/mean length"].split("/")[1]
        data[sample_name]["Target_region_size"] = "NA"

    #print data
    ofile = open(out_file, 'w')
    jstr = json.dumps(data)
    ofile.write(jstr)
    ofile.close()

    # header=data.keys()
    # headers={}
    # for key in header:
    #      headers[key] = key
    # data = pd.read_csv(out_file)
    # data.set_index("sample_name", inplace=True)
    # data.loc[["E100010645_L01_115"],["Insert size Mean"]]=0
    # print data.loc[["E100010645_L01_115"],["Insert size Mean"]]
    # data.to_csv(out_file)
    #with open(out_file,'a') as f:
    #    f_csv = csv.DictWriter(f,header)
        #f_csv.writeheader()
    #    f_csv.writerow(data)


    

