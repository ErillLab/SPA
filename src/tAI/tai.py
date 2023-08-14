
'''

Definition of the TAI class.
Objects belonging to this class are tAI calculators.
The calculator can be used to calculate tAI for an input CDS, by assuming the
tRNA counts (gene copy numbers) specified when initializing the object.
It is also possible to update the object with a new set of tRNA counts.

'''


import biochem as bc
import utils as ut


'''
Lysidine (L, N*, or k2C) (see Figure 1C) is a lysine (Lys)-containing cytidine derivative that occurs at the wobble position (position 34) of the isoleucine tRNAs (tRNAIle) that bear the CAT anticodon in their tDNA sequence and are specific for the AUA codon

Notably, the CAT anticodon also occurs in tRNAMet responsible for AUG decoding. The lysidine modification converts the codon specificity of a tRNA with the CAT anticodon from AUG to AUA and the amino acid specificity from methionine (Met) to isoleucine (Ile), thus preventing the misreading of AUG as Ile and AUA as Met.

https://doi.org/10.1016/j.molcel.2005.06.007
'''


class TAI(object):

    def __init__(self, tRNA_counts_dict):
        
        self.tGCN_dict = None
        self.w_dict = None
        
        # Set attributes
        self.update(tRNA_counts_dict)
    
    def set_tGCN_dict(self, tRNA_counts_dict):
        ''' Sets the `tRNA_counts_dict` attribute. '''
        
        # Check data type
        if not isinstance(tRNA_counts_dict, dict):
            raise ValueError("tRNA_counts_dict type should be dict, not " +
                             type(tRNA_counts_dict).__name__)
        
        # Make a dictionary of tRNA gene copy numbers
        tGCN_dict = {}
        for k in tRNA_counts_dict.keys():
            tGCN_dict[k.lower().replace("u","t")] = tRNA_counts_dict[k]
        # Set to a count of 0 all the missing anticodons
        for triplet in bc.codons:
            if triplet not in tGCN_dict.keys():
                tGCN_dict[triplet] = 0
        
        # Set attribute
        self.tGCN_dict = tGCN_dict
    
    def get_s(self, codon, anticodon):
        ''' Returns the s value (affinity) of the codon-anticodon pairing. '''
        
        if anticodon not in bc.anticodons(codon):
            raise ValueError(anticodon, "is not  an anticodon recognizing", codon)
        
        # The 3rd base of the codon (`third_base`) pairs with the first base of
        # the anticodon (`recognizer`)
        third_base = codon[2]
        recognizer = anticodon[0]
        
        # Modification of Adenine to Inosine is assumed to be the standard in
        # Bacteria for anticodons that have an Adenine in 1st position (which
        # pairs with the 3rd position of the codon).
        if recognizer == 'a':
            recognizer = 'i'
            
        # The ATA codon for Isoleucine can be recognized by a modified CAT
        # anticodon, where Cytosine is modified to Lysidine (L). The L in the
        # LAT anticodon will pair with the third base (Adenine) of the ATA
        # codon.
        if (codon, anticodon) == ('ata', 'cat'):
            recognizer = 'l'
        
        return bc.s_dict[third_base + recognizer]
    
    def get_naive_s(self, third_codon_base, first_anticodon_base):
        ''' Returns a naive s value (affinity) given the recognition of the 3rd
        codon base by the first anticodon base. Perfect matching on the first
        two codon positions is assumed. Returns 0 or 0.5. '''
        
        # Perfect match
        if bc.complement[third_codon_base] == first_anticodon_base:
            return 0
        # Mismatch at 3rd position (wobble)
        else:
            return 0.5
    
    def get_W(self, codon):
        ''' Returns the W value for the input codon.
        
        Don't confuse W with w.
        w = W / max(W)  when W is not 0, otherwise w is the geometric mean of
        all the W values that are not 0.
        '''
        
        # Old code reflecting R code behavior
        # # "Check Methionine"
        # # As in lines 125,126 in https://github.com/mariodosreis/tai/blob/master/R/tAI.R
        # if codon == 'atg':
        #     return (1-self.get_s(codon, 'cat')) * self.tGCN_dict['cat']
        
        
        # No W for Methionine or STOP codons.
        if codon in ['atg', 'tga', 'taa', 'tag']:
            return None
        
        # In bacteria, the recognition of the 'ata' (AUA) codon for Isoleucine
        # occurs through a modified CAT anticodon (see comments in the get_s
        # function). The abundance of that anticodon is unknown because there
        # is no automatic way to differentiate between 'START' Met-tRNA genes
        # and normal Met-tRNAs in any genome. Moreover, the frequency of
        # CAT -> LAT modifications is another unknown variable. Following the
        # original R code [https://github.com/mariodosreis/tai/blob/master/R/tAI.R
        # lines 128, 129], we assign a token tGCN count of 1 for this pairing.
        elif codon == 'ata':
            return (1-self.get_s(codon, 'cat')) * 1
        
        else:
            W = 0
            for anticodon in bc.anticodons(codon):
                # print("\t+ (1 - s({}))\t* tGCN({}) =".format(anticodon[::-1], anticodon[::-1]))
                # print("\t  (1 - {})   \t* {}".format(self.get_s(codon, anticodon),
                #                                  self.tGCN_dict[anticodon]))
                # print("")
                W += (1-self.get_s(codon, anticodon)) * self.tGCN_dict[anticodon]
            # print("W =", W)
            # print("")
            return W
    
    def set_w_dict(self):
        ''' Sets the `w_dict` attribute. First computes W_dict, then computes
        w_dict by deviding all the W values by the largest W value. '''
        
        # Make codon -> W dictionary
        W_dict = {}
        for codon in bc.codons:
            # Ignore Methionine and STOP codons
            if codon not in ['atg', 'taa', 'tag', 'tga']:
                W_dict[codon] = self.get_W(codon)
        
        # Make codon -> w dictionary
        max_W = max(W_dict.values())
        w_dict = {}
        for k in W_dict.keys():
            w_dict[k] = W_dict[k] / max_W
        
        # Substitute zeros with geometric mean of non-zeros
        if len([w for w in w_dict.values() if w == 0]) > 0:
            gm_non_zeros = ut.geo_mean([w for w in w_dict.values() if w != 0])
            for k in w_dict.keys():
                if w_dict[k] == 0:
                    w_dict[k] = gm_non_zeros
        
        # Set the w_dict attribute
        self.w_dict = w_dict
    
    def get_tai(self, seq):
        ''' Returns the tAI index of the input coding sequence. '''
        seq_codons = [codon.lower() for codon in bc.get_sequence_codons(seq)]
        ws = [self.w_dict[codon] for codon in seq_codons]
        return ut.geo_mean(ws)
    
    def update(self, tRNA_counts_dict):
        ''' Update the TAI object according to a new input set of tGCN values,
        (the `tRNA_counts_dict` dictionary). '''
        self.set_tGCN_dict(tRNA_counts_dict)
        self.set_w_dict()
    
    def add_tRNA_genes(self, tRNA_counts_dict, expr_level=1):
        ''' Add extra tRNA gene copy numbers (specified by the `tRNA_counts_dict`
        dictionary) to the already existing ones.
        The `expr_level` parameter can be used to tune the expression level of the
        added genes. This is useful to model the introduction of viral tRNA genes,
        that may have a different gene expression level than the host's. The tGCN
        values (gene copy numbers) to be added are first multiplied by expr_level.
        The default expression level is 1, meaning that the added genes are
        supposed to be expressed as much as the already present genes. '''
        
        # Check data type
        if not isinstance(tRNA_counts_dict, dict):
            raise ValueError("tRNA_counts_dict type should be dict, not " +
                             type(tRNA_counts_dict).__name__)
        
        # Add to present tRNA genes new genes in a proportion of 1 to expr_level
        for k in tRNA_counts_dict.keys():
            self.tGCN_dict[k.lower().replace("u","t")] += tRNA_counts_dict[k] * expr_level
        
        
    






