#!/usr/bin/env python
'''
Setting global variables.
'''

#####################
# IMPORT OPERATIONS #
#####################

import MyExceptions as ME

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2018 Michael Gruenstaeudl'
__info__ = 'annonex2embl'
__version__ = '2018.03.26.2000'

#############
# DEBUGGING #
#############

import pdb
# pdb.set_trace()


#############
# VARIABLES #
#############

warnings = []

allowed_marker_abbrev = ['ITS','18S','28S','trnK','matK','IGS','genomic_CDS']
minimum_rrna_marker_abbrev = ['18S','28S']
genomic_CDS_marker_abbrev = ['gene','CDS']

#############
# FUNCTIONS #
#############

def getOutlist(checklist_type):
    if checklist_type == 'ITS':
        return ['entrynumber','organism_name','isolate','env_sam','country','spec_vouch','RNA_18S','ITS1_feat','RNA_58S','ITS2_feat','RNA_26S','sequence']
    elif checklist_type == 'rRNA':
        return ['entrynumber','organism_name','sediment','isolate','isol_source','country','lat_lon','collection_date','sequence']
    elif checklist_type == 'trnK_matK':
        return ['entrynumber','organism_name','fiveprime_cds','threeprime_cds','fiveprime_partial','threeprime_partial','trnK_intron_present','isolate','spec_vouch','country','ecotype','sequence']
    elif checklist_type == 'IGS':
        return ['entrynumber','organism_name','env_sam','gene1','g1present','gene2','g2present','isolate','spec_vouch','country','sequence']
    elif checklist_type == 'genomic_CDS':
        return ['entrynumber','organism_name','env_sam','gene_symbol','product_name','transl_table','fiveprime_cds','threeprime_cds','fiveprime_partial','threeprime_partial','read_frame','isolate','spec_vouch','country','ecotype','sequence']
    elif checklist_type == 'ETS':
        return ['entrynumber','organism_name','ets_type','isolate','clone','strain','variety','cultivar','breed','ecotype','mating_type','sex','isolation_source','host','tissue','country','area','locality','lat_lon','coldate','col_by','cult_coll','spech_vouch','bio_mat','fwd_name1','fwd_seq1','rev_seq1','rev_name1','fwd_name2','fwd_seq2','sequence']
    else:
        raise ME.CheckListTypeNotKnownError(checklist_type + ' is not a known checklist type')

def getOutdict(checklist_type):
    outlist = getOutlist(checklist_type)
    outdict = {}
    for elem in outlist:
        outdict.update({elem:''})

    return outdict
