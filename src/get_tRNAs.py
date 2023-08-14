#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 10:45:48 2023

@author: micromol
"""

import csv, os
from itertools import islice
from Bio import Entrez, SeqIO
Entrez.email = "your.email@gmail.com"


##Get phage and host genomes from list and download from NCBI

# Downloads genome from NCBI in fasta format using Nucleotide accession
def acc_to_genome(nucacc, filename):
    '''
    Parameters
    ----------
    nucacc : nucleotide accession of the genome.
    filename : filename of the output.

    Returns
    -------
    Fasta file of the genome.

    '''
    handle = Entrez.efetch(db="nuccore", id=nucacc, rettype="fasta", retmode="text")
    record = handle.read()    
    out_handle = open(filename, "w")
    out_handle.write(record)
    out_handle.close()
    handle.close()


# get host genomes
# in this case, because we don't have host genome accessions we have to search for them. Because there can be several genomes for a single 
# host, we will prioritize the genomes according to the following dictionary, and will keep the one with the maximum score.
#dictionary used to assign priority to accession numbers
priority = {"NC_": '7', "AC_": '7', "AE": '6', "CP": '6', "CY": '6', "LT": "6", "AL": "6", "FN": "6", \
                    "NZ_": '5', "NT_": '5', "NW_": '5', "AAAA-AZZZ": '4',\
                    "U": '3', "AF": '3', "AY": '3', "DQ": '3'}

def host_reference_genomes(path_to_phages, out_file):
    '''
    Parameters
    ----------
    path_to_phages : path to the csv file with the phage and host information. In this case, phage acc is found in the first column
        and host name is found in the second column.
    out_file : csv file were the host name and accession will be stored

    Returns
    -------
    csv file were the host name and accession will be stored.

    '''
    path_to_phages = path_to_phages
    host_and_phage = {}
    with open(path_to_phages) as csv_file:
        for i in range(1):
            next(csv_file)
        for line in csv.reader(csv_file, delimiter=","):
            phage = line[0]
            host_without_strain = line[1].split()[:2]
            host = " ".join(host_without_strain)
            if host not in host_and_phage.keys():
                host_and_phage[str(host)] = []
            host_and_phage[str(host)].append(phage)
    
    sequenced_strains_accessions = {}
    host_accessions = []
    for host in host_and_phage:
        search_term = '((('+host+') AND complete) AND (genome OR chromosome) NOT (phage OR plasmid))'
        genome_record = Entrez.read(Entrez.esearch(db="nuccore", term=search_term, idtype="acc", field="title", retmax=10000000))
        if int(genome_record["Count"]) == 0:
            genome_record2 = Entrez.read(Entrez.esearch(db="nuccore", term=search_term, idtype="acc"))
            host_accessions.append(genome_record2["IdList"])
        else:
            host_accessions.append(genome_record["IdList"])
        acc_scr = 0
        for acc_1 in host_accessions:
            for acc_2 in acc_1:
                for key in priority:
                    if acc_2.startswith(key):
                        acc_scr = priority[key]
                        mylist = [acc_2, acc_scr]
                        if host not in sequenced_strains_accessions.keys():
                            sequenced_strains_accessions[str(host)] = []
                        sequenced_strains_accessions[str(host)].append(mylist)
    
    prioritized_sequenced_strains_accessions = {}        
    for host in sequenced_strains_accessions.keys():
        max_genome_score = 0
        selection_score = 7
        max_genome_acc = 0
        for record in sequenced_strains_accessions[host]:
            if int(record[1]) >= selection_score:
                if int(record[1]) == max_genome_score:
                    max_genome_score = record[1]
                    max_genome_acc = "%s, %s" % (max_genome_acc, record[0])
                else:
                    if int(record[1]) > int(max_genome_score):
                        max_genome_score = record[1]
                        max_genome_acc = record[0]
            elif int(record[1]) >= int(max_genome_score):
                max_genome_score = record[1]
                max_genome_acc = record[0]
        prioritized_sequenced_strains_accessions[str(host)] = max_genome_acc
    
    out_file = out_file
    out_handle = open(out_file, "a+")
    out_handle.write("Host,Accessions\n")
    for host in prioritized_sequenced_strains_accessions.keys():
        out_handle.write("%s," % (host))
        out_handle.write("%s," % (prioritized_sequenced_strains_accessions[host]))
        out_handle.write("\n")


def run_trnascan(filename, result, tRNAscan_results_path, genomes_path):
    '''
    Run tRNAscan-SE
    For each of the files, construct the command and run it. (https://seaphages.org/media/docs/Predicting_tRNA_and_tmRNA_genes_12-2-16.pdf)
    -B for bacterial tRNAs
    -I search using infernal
    -H show breakdown of primary and secondary structure components to covariance model bit scores.
    -D disable pseudogene checking.
    -o final results in table format
    -f tRNA secondary structures
    -X set cutoff score, in this case to 17.
    
    Parameters
    ----------
    filename : phage genome in fasta format
    result: name of filename without extension
    tRNAscan_results_path: path where tRNAscan results will be saved
    genomes_path = path where the phage genomes in fasta format can be found

    Returns
    -------
    Creates an output file for each phage genome in the specified folder.
    
    '''
     
    print(filename, result)
    outfile_path = tRNAscan_results_path + result
    infile_path = genomes_path + filename
    os.system("tRNAscan-SE -B -I -H -D -o %s.out -f %s.ss -X17 %s" % (outfile_path, outfile_path, infile_path))


def run_aragorn(filename, result, aragorn_results_path, genomes_path):
    '''
    Run ARAGORN
    For each of the files, construct the command and run it. (https://seaphages.org/media/docs/Predicting_tRNA_and_tmRNA_genes_12-2-16.pdf)
    -t look for tRNAs
    -m look for tmRNAs
    -gcbact use bacterial genetic code
    -c assume sequences are circular (because phage genomes circularize when infecting)
    -o path of outfile

    Parameters
    ----------
    filename : phage genome in fasta format
    result: name of filename without extension
    aragorn_results_path: path where Aragorn results will be saved
    genomes_path = path where the phage genomes in fasta format can be found

    Returns
    -------
    Creates an output file for each phage genome in the specified folder.

    '''
    
    print(filename, result)
    outfile_path = aragorn_results_path + result
    infile_path = genomes_path + filename
    os.system("aragorn -t -m -gcbact -c -o %s.aout %s" % (outfile_path, infile_path))



def genome_to_prot(genome, filename):
    '''
    Get the proteome of a provided genome.    

    Parameters
    ----------
    genome : gb file of the phage you want to obtain the proteome from
        
    filename : file name of the output

    '''
    
    record = SeqIO.read(genome, "gb")    
    out_handle = open(filename, "a+")
    for feature in record.features:
        if feature.type == "source":
            try:
                taxon = feature.qualifiers["db_xref"][0].split(":")[-1]
            except KeyError:
                taxon = "no info"
    for feature in record.features:
        if feature.type == "CDS":
            try:
                proteinID = feature.qualifiers["protein_id"][0]
            except KeyError:
                proteinID = "no_protein_ID"
            try:
                product = feature.qualifiers["product"][0]
            except KeyError:
                product = "no_product_info"
            try:
                locustag = feature.qualifiers["locus_tag"][0]
            except KeyError:
                locustag = "no_locus_tag"
            try: 
                geneID = feature.qualifiers["db_xref"][0]
            except KeyError:
                geneID = "no_GeneID"
            try: 
                translation = feature.qualifiers["translation"][0]
            except KeyError:
                continue                
            out_handle.write(">%s|%s|%s|%s|%s|locus_tag:%s|%s|txid:%s" % (proteinID.replace(" ", "_"), product.replace(" ", "_"), record.id.replace(" ", "_"), record.description.replace(" ", "_"), str(feature.location).replace(" ", "_"), locustag.replace(" ", "_"), geneID.replace(" ", "_"), taxon.replace(" ", "_")) + "\n"+"%s" % (translation) + "\n")
    out_handle.close()



def parse_trnas(result_files, accession, outfile, reference_proteome, tRNAscan_results_path, aragorn_results_path):
    '''
    Parses the output generated by tRNAscan_SE and Aragorn, to obtain a single csv file per phage that
    contains all the desired information of each of the tRNAs predicted by each software. 

    Parameters
    ----------
    result_files : tRNAscan or Aragorn output files. 
    
    accession: phage accession to analize.    
        
    outfile : output file
        
    reference_proteome : proteome of the phage. 
    
    tRNAscan_results_path : path to tRNAscan results
    
    aragorn_results_path : path to aragorn results

    '''
    infernal_score = 35.0
    
    for filename in result_files:
        if accession in filename:
            
            if filename.endswith(".out"):
                
                with open(tRNAscan_results_path+filename) as file:
                    tRNAscan_csv_dict = {}
                    i = 1
                    
                    #Skip the first three lines (description)
                    for j in range(3):
                        try:
                            next(file)
                        except StopIteration:
                            continue                            
                    
                    for line in csv.reader(file, delimiter="\t", skipinitialspace=True):
                        if line:
                            tRNAlen = float(line[3])-float(line[2])
                            
                            #If column 9 has a value of more than 35 (Infernal score), and the length of the tRNA is less than 90, create a info_dict and append all the information.
                            if float(line[8]) >= infernal_score and -90 < tRNAlen < 90:
                                csv_info_dict = {}
                                csv_info_dict[str("tRNA prediction")] = []
                                csv_info_dict[str("tRNA #")] = []
                                csv_info_dict[str("tRNA Begin")] = []
                                csv_info_dict[str("tRNA End")] = []
                                csv_info_dict[str("tRNA Strand")] = []
                                csv_info_dict[str("tRNA Type")] = []
                                csv_info_dict[str("Anticodon")] = []
                                csv_info_dict[str("tRNA prediction")].append("tRNAscan")
                                csv_info_dict[str("tRNA #")].append(line[1])
                                csv_info_dict[str("tRNA Begin")].append(line[2])
                                csv_info_dict[str("tRNA End")].append(line[3])
                                csv_info_dict[str("tRNA Type")].append(line[4])
                                csv_info_dict[str("Anticodon")].append(line[5])
                                if -90 < tRNAlen < 0:
                                    csv_info_dict[str("tRNA Strand")].append("-")
                                if 0 <= tRNAlen < 90:
                                    csv_info_dict[str("tRNA Strand")].append("+")
                                                                
                            else:
                                continue
                            
                        # Append all the tRNA information compiled in csv_info_dict to csv_dict
                        tRNAscan_csv_dict[str(i)] = []
                        tRNAscan_csv_dict[str(i)].append(csv_info_dict)
                        i = i + 1
                    
                
            if filename.endswith(".aout"):
                with open(aragorn_results_path+filename) as file:
                    for k in range(20):
                        next(file)
                    aragorn_csv_dict = {}
                    previous_tRNA_num = 0
                    m = 1
                    for line in file:
                        if line.startswith("tRNA-", 4, 9):
                            info = []
                            info.append(''.join(islice(file, 6)))
                            tRNAlen = int(info[0].split(" bases")[0])
                            if tRNAlen < 90:
                                csv_info_dict = {}
                                csv_info_dict[str("tRNA prediction")] = []
                                csv_info_dict[str("tRNA #")] = []
                                csv_info_dict[str("tRNA Begin")] = []
                                csv_info_dict[str("tRNA End")] = []
                                csv_info_dict[str("tRNA Strand")] = []
                                csv_info_dict[str("tRNA Type")] = []
                                csv_info_dict[str("Anticodon")] = []
                                csv_info_dict[str("tRNA prediction")].append("Aragorn")
                                try:
                                    csv_info_dict[str("tRNA #")].append(int(info[0].split(".\n")[0].split("]\n")[1]) - 1)
                                    previous_tRNA_num = int(info[0].split(".\n")[0].split("]\n")[1]) - 1
                                except ValueError:
                                    csv_info_dict[str("tRNA #")].append(int(previous_tRNA_num) + 1)
                                csv_info_dict[str("tRNA Begin")].append(info[0].split("[")[1].split(",")[0])
                                csv_info_dict[str("tRNA End")].append(info[0].split("[")[1].split(",")[1].split("]")[0])
                                csv_info_dict[str("tRNA Type")].append(line.split("-")[1].split("(")[0])
                                if line.split("-")[1].split("(")[0] == "?":
                                    csv_info_dict[str("Anticodon")].append(line.split("(")[2].split(")")[0])
                                else:
                                    csv_info_dict[str("Anticodon")].append(line.split("(")[1].split(")")[0])
                                if info[0].split("[")[0].split("Sequence")[1] == " c":
                                    csv_info_dict[str("tRNA Strand")].append("-")
                                if info[0].split("[")[0].split("Sequence")[1] == " ":
                                    csv_info_dict[str("tRNA Strand")].append("+")
                                                               
                            else:
                                continue
                            
                        else:
                            continue
                    
                        aragorn_csv_dict[str(m)] = []
                        aragorn_csv_dict[str(m)].append(csv_info_dict)
                        m = m + 1
                

    if bool(tRNAscan_csv_dict) or bool(aragorn_csv_dict):
        with open(outfile, mode="w") as csv_file:
            fieldnames = ["tRNA prediction", "tRNA #", "tRNA Begin", "tRNA End", "tRNA Strand", "tRNA Type", "Anticodon"]
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            #For each of the keys in csv dict (1, 2, 3, 4, etc), write their value in each of the rows of the csv file. 
            for key in tRNAscan_csv_dict.keys():
                for item in tRNAscan_csv_dict[key]:
                    writer.writerow(item)
            for key in aragorn_csv_dict.keys():
                for item in aragorn_csv_dict[key]:
                    writer.writerow(item) 
                    

def filter_tRNAs(parsed_path, filename, filtered_path):
    """
    Filters the parsed tRNAscan and aragorn results. 
    Also determines if a tRNA is predicted by both tRNAscan and aragorn 
    (same anticodon, location about 10bp different max).
    If there is a match, it will keep the aragorn location. In the case where 
    the anticodon is CAT, aragorn always determines the type as Met, but tRNAscan 
    determines Ile2 or fMET, so it keeps the tRNAscan determined type. Also if the 
    tRNA is undetermined by tRNAscan, it will keep the aragorn tRNA type.

    Parameters
    ----------
    parsed_path : path where all the parsed tRNAscan and aragorn results can be found
    filename : filename of each file to filter.
    filtered_path : path where the filtered results should be saved. 

    Returns
    -------
    None. Outputs csv files with filtered results. 

    """
    tRNA_dict = {}
    with open(parsed_path+filename) as csv_file:
        for i in range(1):
            next(csv_file)
        tRNA_num = 1
        for line in csv.reader(csv_file, delimiter="\t"):
            tRNA_info_dict = {}
            tRNA_prediction = line[0].split("'")[1].split("'")[0]
            tRNA_in_dict = "Not yet"
            
            # Determine all the info on the tRNA and append it to a dictionary
            tRNA_strand = line[4].split("'")[1].split("'")[0]
            tRNA_type = line[5].split("'")[1].split("'")[0]
            tRNA_anticodon = line[6].split("'")[1].split("'")[0].lower()
            tRNA_info_dict[str("tRNA strand")] = []
            tRNA_info_dict[str("tRNA type")] = []
            tRNA_info_dict[str("tRNA anticodon")] = []
            tRNA_info_dict[str("tRNA strand")].append(tRNA_strand)
            tRNA_info_dict[str("tRNA type")].append(tRNA_type)
            tRNA_info_dict[str("tRNA anticodon")].append(tRNA_anticodon)
            
            if tRNA_strand == "+":
                tRNA_begin = line[2].split("'")[1].split("'")[0]
                tRNA_end = line[3].split("'")[1].split("'")[0]
                tRNA_info_dict[str("tRNA begin")] = []
                tRNA_info_dict[str("tRNA end")] = []
                tRNA_info_dict[str("tRNA begin")].append(tRNA_begin)
                tRNA_info_dict[str("tRNA end")].append(tRNA_end)
            if tRNA_strand == "-":
                if tRNA_prediction == "tRNAscan":
                    tRNA_begin = line[3].split("'")[1].split("'")[0]
                    tRNA_end = line[2].split("'")[1].split("'")[0]
                    tRNA_info_dict[str("tRNA begin")] = []
                    tRNA_info_dict[str("tRNA end")] = []
                    tRNA_info_dict[str("tRNA begin")].append(tRNA_begin)
                    tRNA_info_dict[str("tRNA end")].append(tRNA_end)
                if tRNA_prediction == "Aragorn":
                    tRNA_begin = line[2].split("'")[1].split("'")[0]
                    tRNA_end = line[3].split("'")[1].split("'")[0]
                    tRNA_info_dict[str("tRNA begin")] = []
                    tRNA_info_dict[str("tRNA end")] = []
                    tRNA_info_dict[str("tRNA begin")].append(tRNA_begin)
                    tRNA_info_dict[str("tRNA end")].append(tRNA_end)
            
            if tRNA_prediction == "tRNAscan":
                if bool(tRNA_dict):
                    for tRNA in tRNA_dict.keys():
                        if tRNA_dict[tRNA][0][str("tRNA anticodon")][0] == tRNA_anticodon and -10 < int(tRNA_dict[tRNA][0][str("tRNA begin")][0]) - int(tRNA_begin) < 10 and -10 < int(tRNA_dict[tRNA][0][str("tRNA end")][0]) - int(tRNA_end) < 10:
                            tRNA_in_dict = "Yes"
                            if tRNA_anticodon == "cat":
                                if tRNA_type == "Met":
                                    # If the anticodon is CAT and the tRNA type determined by aragorn is Met, keep the tRNAscan anticodon
                                    tRNA_dict[tRNA][0][str("tRNA type")] = tRNA_type
                        else:
                            continue
                    
                    if tRNA_in_dict == "Not yet":
                        tRNA_dict[str(tRNA_num)] = []
                        tRNA_dict[str(tRNA_num)].append(tRNA_info_dict)
                        tRNA_num = tRNA_num + 1
                else:
                    tRNA_dict[str(tRNA_num)] = []
                    tRNA_dict[str(tRNA_num)].append(tRNA_info_dict)
                    tRNA_num = tRNA_num + 1
            
            if tRNA_prediction == "Aragorn":
                if bool(tRNA_dict):
                    for tRNA in tRNA_dict.keys():
                        if tRNA_dict[tRNA][0][str("tRNA anticodon")][0] == tRNA_anticodon and -10 < int(tRNA_dict[tRNA][0][str("tRNA begin")][0]) - int(tRNA_begin) < 10 and -10 < int(tRNA_dict[tRNA][0][str("tRNA end")][0]) - int(tRNA_end) < 10:
                            # If the tRNA is already in the dictionary we will update the location of the tRNA, because aragorn tends to be more accurate than tRNAscan in this. 
                            tRNA_in_dict = "Yes"
                            tRNA_dict[tRNA][0][str("tRNA begin")] = tRNA_begin
                            tRNA_dict[tRNA][0][str("tRNA end")] = tRNA_end
                            if tRNA_dict[tRNA][0][str("tRNA type")] == "Undet":
                                # If the tRNA type already in the dictionary (tRNAscan) is Undet, change it for the one specified by aragorn
                                tRNA_dict[tRNA][0][str("tRNA type")] = tRNA_type
                        else:
                            continue
                    
                    if tRNA_in_dict == "Not yet":
                        tRNA_dict[str(tRNA_num)] = []
                        tRNA_dict[str(tRNA_num)].append(tRNA_info_dict)
                        tRNA_num = tRNA_num + 1
                else:
                    tRNA_dict[str(tRNA_num)] = []
                    tRNA_dict[str(tRNA_num)].append(tRNA_info_dict)
                    tRNA_num = tRNA_num + 1
                    
    with open(filtered_path + filename, "w+") as out_file:
        fieldnames = ["tRNA begin", "tRNA end", "tRNA strand", "tRNA type", "tRNA anticodon"]
        writer = csv.DictWriter(out_file, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for key in tRNA_dict.keys():
            for item in tRNA_dict[key]:
                writer.writerow(item)
           

def tRNA_count(filtered_path, filename, tRNA_count_path):
    '''
    Parameters
    ----------
    filtered_path : path to the filtered results from tRNA prediction
    filename : file of the filtered results from tRNA prediction
    tRNA_count_path : path to the directory where the tRNA count files have to be saved

    Returns
    -------
    None.

    '''
    with open(filtered_path+filename) as csv_file:
        for i in range(1):
            next(csv_file)
        counter = {}
        for line in csv.reader(csv_file, delimiter="\t"):
            if line[4] not in counter:
                count = 0
                anticodon = line[4].split("']")[0].split("['")[1]
                counter[anticodon] = str(count)
            count +=1
            counter[anticodon] = str(count)
            
    
    with open(tRNA_count_path + filename, "w+") as outfile:
        outfile.write("tRNA anticodon, count\n")
        w = csv.writer(outfile)
        w.writerows(counter.items())
   
