
'''

Definition of functions and variables that model biochemical relationships to
be used by tai.py.

'''


import utils as ut
import warnings


# Load settings
settings = ut.read_json_file("tAI_settings.json")


bases = "tcag"
codons = [x+y+z for x in bases for y in bases for z in bases]
complement = {'t':'a', 'u':'a', 'a':'t', 'c':'g', 'g':'c'}


def rev_compl(sequence, use_uracil=False):
    ''' Returns the reverse complement of the input sequence. '''
    rev_compl_string = ""
    for x in sequence[::-1].lower():
        rev_compl_string += complement[x]
    if use_uracil:
        return rev_compl_string.replace("t", "u")
    else:
        return rev_compl_string

def anticodons(codon):
    ''' Returns the list of anticodons that pair perfectly with the first two
    bases of the input codon. '''
    return [b + rev_compl(codon[:2]) for b in bases]

def get_sequence_codons(seq):
    ''' Returns the list of codons composing the input sequence `seq`. The
    codons are returned in the order in which they appear. However, methionine
    codons (included start codon) and stop codons are ignored (skipped).
    '''
    if int(len(seq)) % 3 != 0:
        warnings.warn("The input sequence is " + str(int(len(seq))) +
                      " bp long, which is not a multiple of 3.")
        # raise ValueError("The input sequence is " + str(int(len(seq))) +
        #                  " bp long, which is not a multiple of 3.")
    seq = seq.lower().replace("u", "t")
    codons_list = [seq[i:i+3] for i in range(0, (len(seq)//3)*3, 3)]
    # Ignore Methionine and STOP codons
    codons_to_remove = {'atg', 'taa', 'tga', 'tag'}
    return [codon for codon in codons_list if codon not in codons_to_remove]

def load_s_dict(filepath):
    ''' Makes a dictionary from the input CSV file for s values.
    In the dictionary, each key is a codon-anticodon pairing at the
    third position of the codon, and each corresponding value is the s value
    that will be used to model affinity. '''
    s_dict = {}
    with open(filepath, "r") as f:
        for line in f:
            try:
                pairing, s = line.rstrip().split(",")
            except:
                raise ValueError("The file " + filepath + " should be a " +
                                 "2-columns CSV, where each line follows " +
                                 "the format xy,0.12")
            pairing = pairing.lower().replace("u", "t")
            s_dict[pairing] = float(s)
    for codon_base in bases:
        for anticodon_base in bases + 'il':
            if codon_base + anticodon_base not in s_dict.keys():
                s_dict[codon_base + anticodon_base] = 1
    return s_dict



s_dict = load_s_dict("../../datasets/" + settings['s_table_filename'])
use_naive_s_values = settings['use_naive_s_values']







