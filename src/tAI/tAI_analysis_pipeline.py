
"""

Main pipeline to get the results on tAI for the paper.

"""

from Bio import SeqIO
from tai import TAI
import warnings
import time
import os

import matplotlib.pyplot as plt
# import seaborn as sns
import pandas as pd
import copy
from scipy.stats import wilcoxon, mannwhitneyu


tRNA_data_dir = "../../datasets/tRNAs/"
output_dir = "../../results/tAI/" + time.strftime("%Y%m%d%H%M%S") + "/"
os.mkdir(output_dir)




def load_tGCN_dict(filepath):
    '''
    Reads the table of anticodon - tRNA counts values, and returns a dictionary
    where the keys are the tRNA anticodons and the values are corresponding
    tRNA gene counts. It also sets to 0 the counts for any missing triplet.
    '''
    triplets = [x+y+z for x in 'tcag' for y in 'tcag' for z in 'tcag']
    
    tGCN_dict = {}
    with open(filepath, "r") as f:
        for line in f:
            anticodon, count = line.rstrip().split(",")
            # Ignore non-standard anticodons / column names line
            if anticodon not in triplets:
                continue
            tGCN_dict[anticodon.lower().replace("u", "t")] = float(count)
    
    # Set to a count of 0 all the missing anticodons
    for triplet in triplets:
        if triplet not in tGCN_dict.keys():
            tGCN_dict[triplet] = 0
    return tGCN_dict

def is_complete_CDS(sequence, gtg_start=False, ttg_start=False):
    '''
    Checks that the input sequence is oriented correctly (5' to 3') and that it
    is a complete coding sequence (CDS). To be complete (not truncated), it has
    to start with s START codon and end with a STOP codon. Returns a boolean.
    '''
    start_codons = ['ATG']
    stop_codons = ["TAA", "TGA", "TAG"]
    if gtg_start:
        start_codons.append("GTG")
    if ttg_start:
        start_codons.append("TTG")
    
    if sequence[:3].upper() not in start_codons:
        return False
    elif sequence[-3:].upper() not in stop_codons:
        return False
    else:
        return True

def load_gene_seq(filepath):
    '''
    Reads a fasta file and returns the sequence as a lowercase string. If the
    number of records in the file is not 1, a warning will be generated.
    '''
    records = SeqIO.parse(open(filepath),'fasta')
    seq_list = []
    for record in records:
        seq_list.append(record.seq)
    if len(seq_list) > 1:
        raise ValueError("More than one sequence record in fasta file")
    elif len(seq_list) == 0:
        raise ValueError("No sequence record found in fasta file")
    else:
        # check_CDS will generate a warning if the CDS is not complete
        if not is_complete_CDS(seq_list[0], gtg_start=True, ttg_start=True):
            warnings.warn("The sequence in " + filepath + "is an incomplete CDS") 
        return str(seq_list[0]).lower()

def frange(minimum, maximum, step=1.0):
    '''
    Makes an iterator that can be used to iterate over float values from
    `minimum` to `maximum` (included) with a float step (default step is 1.0).
    '''
    out_of_bounds = False
    count = 0
    while not out_of_bounds:
        temp = float(minimum + count * step)
        if step > 0 and temp >= maximum:
            out_of_bounds = True
        elif step < 0 and temp <= maximum:
            out_of_bounds = True
        yield temp
        count += 1

def mixture_dict(dict_1, dict_2, alpha=1):
    '''
    Adds the values of the two dictionaries in a proportion of 1:alpha and
    returns the mixture dictionary. I.e., in the returned dictionary, for every
    key k, the corresponding value will be
    v_1 + (alpha * v_2)
    where v_1 is the value from dict_1 and v_2 is the value from dict_2. An
    error is raised if the two input dictionaries don't have the same keys.
    '''
    # Check keys
    if set(dict_1.keys()) != set(dict_2.keys()):
        raise ValueError("The two input dictionaries don't have the same keys")
    
    for k in dict_2.keys():
        dict_1[k] += dict_2[k] * alpha
    return dict_1

def selective_mixture_dict(dict_1, dict_2, alpha=1):
    '''
    Returns a combination of dict_1 and dict_2, like in the mixture_dict
    function, except that values from dict_2 can only contribute with the
    starting value in dict_1 is not 0. '''
    # Check keys
    if set(dict_1.keys()) != set(dict_2.keys()):
        raise ValueError("The two input dictionaries don't have the same keys")
    
    for k in dict_1.keys():
        if dict_1[k] == 0:
            continue
        dict_1[k] += dict_2[k] * alpha
    return dict_1

def make_phage_host_dict():
    '''
    Reads the table of phage accession - host accession IDs, and returns a
    dictionary where the keys are phage accessions IDs and the values are the
    accession IDs of the corresponding hosts.
    '''
    filepath = tRNA_data_dir + "phage_to_host_table.csv"
    phage_host_dict = {}
    with open(filepath, "r") as f:
        for line in f:
            phage, host = line.rstrip().split(",")
            # Skip line with column names, if present
            if (phage, host) == ('Phage_accession', 'Host_accession'):
                continue
            phage_host_dict[phage] = host.split(".")[0]
    return phage_host_dict

def make_phage_acc_phage_name_dict():
    '''
    Reads the table with phages info, and returns a dictionary where the keys
    are phage accessions IDs and the values are the corresponding phage names.
    '''
    df = pd.read_csv(tRNA_data_dir + "BE_phages.csv")
    acc_name_dict = {}
    for idx in df.index:
        name = df.loc[idx, 'Phage Name']
        acc = df.loc[idx, 'Accession #']
        if pd.isna(acc):
            continue
        acc_name_dict[acc] = name
    
    # XXX MAybe update this part with GenBank accessions of MiniFlayer and MulchRoom
    acc_name_dict['MiniFlayer'] = 'MiniFlayer'
    acc_name_dict['MulchRoom'] = 'MulchRoom'
    
    return acc_name_dict

def phage_tGCN_file(acc):
    ''' Returns the path of the tGCN file for the input phage accession. '''
    return tRNA_data_dir + "phage_tRNAs/tRNA_count/" + acc + '.csv'

def host_tGCN_file(acc):
    ''' Returns the path of the tGCN file for the input host accession. '''
    return tRNA_data_dir + "host_tRNAs/tRNA_count/" + acc + '.csv'

def phage_MCP_file(acc):
    ''' Returns the path of the Major Capsid Protein file for the input phage accession. '''
    return tRNA_data_dir + "BE_phage_MCPs/" + acc + "_MCP.fas"

def satellite_MCP_file(acc):
    ''' Returns the path of the Major Capsid Protein file for the input satellite accession. '''
    return tRNA_data_dir + "satellite_MCPs/" + acc + "_MCP.fas"

def host_Ribosomal_file(acc):
    ''' Returns the path of the ribosomal proteins file for the input host accession. '''
    return tRNA_data_dir + "BE_host_Ribosomal_CDS/" + acc + "_Ribosomal_CDS.fas"

def host_All_CDS_file(acc):
    ''' Returns the path of the all-CDS file for the input host. '''
    if acc == 'CP032543':
        acc = 'NC_010572'
    return tRNA_data_dir + "BE_host_All_CDS/" + acc + "_All_CDS.fas"

def analyze_host_phage_systems(phage_host_dict):
    '''
    phage_host_dict : dict
        Dictionary where each key is a phage genome accession, and the
        corresponding value is the genome accession of the host it infects.
    
    Returns a dictionary where for each phage in the keys of `phage_host_dict`
    the tAI of the major capsid protein is computed using the tRNA pool of the
    host ("tai_h" key of the output dictionary), or using the tRNA pool of the
    phage ("tai_p" key of the output dictionary). The "host_acc" and the
    "phage_acc" store the host accessions and the phage accessions.
    '''
    
    tai_h_list = []
    tai_p_list = []
    phage_acc_list = []
    host_acc_list = []
    
    for phage_acc in phage_host_dict.keys():
        
        # Host
        host_acc = phage_host_dict[phage_acc]
        if host_acc == 'NZ_AMLP00000000':
            continue
        
        # Load host tGCN
        host_tGCN = load_tGCN_dict(host_tGCN_file(host_acc))
        
        # Load phage tGCN
        phage_tGCN = load_tGCN_dict(phage_tGCN_file(phage_acc))
        
        # Load MCP sequence
        MCP_seq = load_gene_seq(phage_MCP_file(phage_acc))
        
        # Compute tAI for MCP with host tGCN
        tai_calc = TAI(host_tGCN)
        tai_h = tai_calc.get_tai(MCP_seq)
        
        # Compute tAI for MCP with phage tGCN
        tai_calc.update(phage_tGCN)
        tai_p = tai_calc.get_tai(MCP_seq)
        
        # Store results
        tai_h_list.append(tai_h)
        tai_p_list.append(tai_p)
        phage_acc_list.append(phage_acc)
        host_acc_list.append(host_acc)
    
    return {"tai_h": tai_h_list,
            "tai_p": tai_p_list,
            "phage_acc": phage_acc_list,
            "host_acc": host_acc_list}

def host_CDS_control(phage_host_dict, control):
    '''
    phage_host_dict : dict
        Dictionary where each key is a phage genome accession, and the
        corresponding value is the genome accession of the host it infects.
    
    control : str
        String that can be "ribosomal" or "all".
    
    Returns a dictionary where for each host in the values of `phage_host_dict`
    the average tAI of some host genes is computed using the tRNA pool of the
    host ("tai_h" key of the output dictionary), or using the tRNA pool of the
    phage ("tai_p" key of the output dictionary). The host genes to be tested
    can be used as a control for the `analyze_host_phage_systems` analysis.
    If `control` is "ribosomal", the host genes (CDS) that will be used are the
    ribosomal proteins of the host. If `control` is "all", all the host CDSs
    will be used. The "host_acc" and the "phage_acc" store the host accessions
    and the phage accessions.
    '''
    
    hosts = list(set(phage_host_dict.values()))
    hosts.sort()
    
    avg_tai_h_list = []
    avg_tai_p_list = []
    host_acc_list = []
    phage_acc_list = []
    
    for host_acc in hosts:
        if host_acc == 'NZ_AMLP00000000':
            continue
        
        # All phages that infect that host
        phages = [k for k in phage_host_dict.keys() if phage_host_dict[k] == host_acc]
        
        # Control set of CDS of the host
        if control.lower() == "ribosomal":
            records = SeqIO.parse(open(host_Ribosomal_file(host_acc)), 'fasta')
        elif control.lower() == "all":
            records = SeqIO.parse(open(host_All_CDS_file(host_acc)), 'fasta')
        else:
            raise ValueError("'control' should be 'ribosomal' or 'all'.")
        
        seq_list = []
        for record in records:
            if is_complete_CDS(record.seq, gtg_start=True, ttg_start=True):
                seq_list.append(str(record.seq))
        print("{}: {} sequences".format(host_acc, len(seq_list)))
        
        for phage_acc in phages:
            
            tai_h_list = []
            tai_p_list = []
            
            for sequence in seq_list:
                
                # Load host tGCN
                host_tGCN = load_tGCN_dict(host_tGCN_file(host_acc))
                
                # Load phage tGCN
                phage_tGCN = load_tGCN_dict(phage_tGCN_file(phage_acc))
                
                # Compute tAI with host tGCN
                tai_calc = TAI(host_tGCN)
                tai_h = tai_calc.get_tai(sequence)
                
                # Compute tAI with phage tGCN
                tai_calc.update(phage_tGCN)
                tai_p = tai_calc.get_tai(sequence)
                
                # Store results
                tai_h_list.append(tai_h)
                tai_p_list.append(tai_p)
            
            avg_tai_h_list.append(sum(tai_h_list)/len(tai_h_list))
            avg_tai_p_list.append(sum(tai_p_list)/len(tai_p_list))
            host_acc_list.append(host_acc)
            phage_acc_list.append(phage_acc)
    
    return {"tai_h": avg_tai_h_list,
            "tai_p": avg_tai_p_list,
            "phage_acc": phage_acc_list,
            "host_acc": host_acc_list}

def analyze_helper_satellite_systems(satellite_helper_dict, phage_host_dict):
    '''
    satellite_helper_dict : dict
        Dictionary where each key is a satellite identifier, and the
        corresponding value is the genome accession of its helper phage.
    
    phage_host_dict : dict
        Dictionary where each key is a phage genome accession, and the
        corresponding value is the genome accession of the host it infects.
    
    Returns a dictionary where for each phage in the keys of `phage_host_dict`
    the tAI of the major capsid protein is computed using the tRNA pool of the
    host ("tai_h" key of the output dictionary), or using the tRNA pool of the
    phage ("tai_p" key of the output dictionary). The "host_acc" and the
    "phage_acc" store the host accessions and the phage accessions.
    '''
    tai_host_list = []
    tai_helper_list = []
    host_acc_list = []
    helper_acc_list = []
    satellite_acc_list = []
    
    for satellite in satellite_helper_dict.keys():
        
        # Helper
        helper_acc = satellite_helper_dict[satellite]
        
        # Host
        host_acc = phage_host_dict[helper_acc]
        if host_acc == 'NZ_AMLP00000000':
            continue
        
        # Load host tGCN
        host_tGCN = load_tGCN_dict(host_tGCN_file(host_acc))
        
        # Load helper phage tGCN
        helper_phage_tGCN = load_tGCN_dict(phage_tGCN_file(helper_acc))
        
        # Load satellite MCP sequence
        MCP_seq = load_gene_seq(satellite_MCP_file(satellite))
        
        # Compute tAI for MCP with host tGCN
        tai_calc = TAI(host_tGCN)
        tai_host = tai_calc.get_tai(MCP_seq)
        
        # Compute tAI for MCP with helper phage tGCN
        tai_calc.update(helper_phage_tGCN)
        tai_helper = tai_calc.get_tai(MCP_seq)
        
        # Store results
        tai_host_list.append(tai_host)
        tai_helper_list.append(tai_helper)
        host_acc_list.append(host_acc)
        helper_acc_list.append(helper_acc)
        satellite_acc_list.append(satellite)
    
    return {"tai_host": tai_host_list,
            "tai_helper": tai_helper_list,
            "host_acc": host_acc_list,
            "helper_acc": helper_acc_list,
            "satellite_acc": satellite_acc_list}

def barplot_hosts_phages(results_dict, title):
    '''
    Make a bar plot to visualize the results from `analyze_host_phage_systems`.
    '''
    xlabels = [phage_acc_name_dict[acc] for acc in results_dict['phage_acc']]
    
    # xlabels = []
    # for i in range(len(results_dict['phage_acc'])):
    #     xlabels.append(results_dict['host_acc'][i] + " - " + results_dict['phage_acc'][i])
    
    index = list(range(len(xlabels)))
    bar_width = 0.35
    
    fig, ax = plt.subplots()
    ax.bar(index, results_dict['tai_h'], bar_width, label="Host tRNA")
    ax.bar([i+bar_width for i in index], results_dict['tai_p'], bar_width, label="Phage tRNA")
    
    ax.set_xlabel('Host-phage systems')
    ax.set_ylabel('tAI')
    ax.set_title(title)
    ax.set_xticks([i + bar_width / 2 for i in index])
    ax.set_xticklabels(xlabels, rotation = 90, fontsize=5)
    ax.legend()
    plt.ylim((0,1))
    plt.show()
    filepath = output_dir + 'barplot_' + title.replace(" ", "_").lower() + '.png'
    fig.savefig(filepath, dpi=600, bbox_inches='tight')
    plt.close()

def barplot_hosts_phages_satellites(mcp_results, satellite_mcp_results):
    '''
    Make a bar plot to visualize the results from `analyze_host_phage_systems`
    and `analyze_helper_satellite_systems`.
    '''
    
    # Prepare data
    results_dict = copy.deepcopy(mcp_results)
    results_dict['tai_h'] += satellite_mcp_results['tai_host']
    results_dict['tai_p'] += satellite_mcp_results['tai_helper']
    results_dict['phage_acc'] += satellite_mcp_results['satellite_acc']
    results_dict['host_acc'] += satellite_mcp_results['host_acc']
    
    # Make plot
    
    xlabels = [phage_acc_name_dict[acc] for acc in results_dict['phage_acc']]
    
    # xlabels = []
    # for i in range(len(results_dict['phage_acc'])):
    #     xlabels.append(results_dict['host_acc'][i] + "\n" + results_dict['phage_acc'][i])

    index = list(range(len(xlabels)))
    bar_width = 0.3

    fig, ax = plt.subplots()
    ax.bar(index, results_dict['tai_h'], bar_width, label="Host tRNA")
    ax.bar([i+bar_width for i in index[:-2]], results_dict['tai_p'][:-2], bar_width, label="Phage tRNA")
    ax.bar([i+bar_width for i in index[-2:]], results_dict['tai_p'][-2:], bar_width, label="Helper phage tRNA")

    #ax.set_xlabel('Phage      Satellite\nMCP          MCP')
    ax.set_ylabel('tAI')
    ax.set_title("Major capsid protein (MCP)")
    ax.set_xticks([i + bar_width / 2 for i in index])
    ax.set_xticklabels(xlabels, rotation = 90, fontsize=5)
    for i in range(len(xlabels)):
        if xlabels[i] in ['MindFlayer', 'MulchMansion', 'MiniFlayer', 'MulchRoom']:
            ax.get_xticklabels()[i].set_color('red')
    
    ax.legend()
    plt.ylim((0,1))
    plt.show()
    
    # Save plot
    filepath = output_dir + "barplot_MCP_phages_and_satellites.png"
    fig.savefig(filepath, dpi=600, bbox_inches='tight')
    plt.close()

# -----------------------------------------------------------------------------

# Input data
phage_host_dict = make_phage_host_dict()
satellite_helper_dict = {"MiniFlayer": "MW291014", "MulchRoom": "MT897905"}
phage_acc_name_dict = make_phage_acc_phage_name_dict()



# Analize pahges MCP tAI
mcp_results = analyze_host_phage_systems(phage_host_dict)
barplot_hosts_phages(mcp_results, "Major Capsid Protein")

# Control experiment with all host CDS
host_all_CDS_results = host_CDS_control(phage_host_dict, 'all')
barplot_hosts_phages(host_all_CDS_results, "All host CDS")

# Control experiment with host ribosomal proteins
host_ribosomal_results = host_CDS_control(phage_host_dict, 'ribosomal')
barplot_hosts_phages(host_ribosomal_results, "Host ribosomal proteins")

# Analize satellite-helper systems
satellite_mcp_results = analyze_helper_satellite_systems(satellite_helper_dict, phage_host_dict)
barplot_hosts_phages_satellites(mcp_results, satellite_mcp_results)
# violinplot_hosts_phages_satellites(mcp_results, satellite_mcp_results)



# Significance

# Wilcoxon signed-rank test
res = wilcoxon(mcp_results['tai_p'], mcp_results['tai_h'], alternative='greater')
print(res)

# Mann-Whitney-U
res = mannwhitneyu(mcp_results['tai_p'], mcp_results['tai_h'], alternative='greater')
print(res)





