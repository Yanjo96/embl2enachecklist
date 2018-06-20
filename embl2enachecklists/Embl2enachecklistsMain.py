#!/usr/bin/env python2.7

''' Main operations in embl2enachecklist '''

#####################
# IMPORT OPERATIONS #
#####################

import MyExceptions as ME
import ChecklistOps as ClOps
import sys
import os

from Bio import SeqIO
from termcolor import colored

# Add specific directory to sys.path in order to import its modules
# NOTE: THIS RELATIVE IMPORTING IS AMATEURISH.
# NOTE: COULD THE FOLLOWING IMPORT BE REPLACED WITH 'import annonex2embl'?

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'embl2enachecklists'))

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2018 Michael Gruenstaeudl'
__info__ = 'embl2enachecklists'
__version__ = '2018.05.23.2000'

#############
# DEBUGGING #
#############

import pdb
# pdb.set_trace()

###########
# CLASSES #
###########

#############
# FUNCTIONS #
#############


def embl2enachecklists(path_to_embl,
                       path_to_outfile,
                       checklist_type=None):

########################################################################

# 1. OPEN OUTFILE
    outp_handle = open(path_to_outfile, 'a')

# 1.1 WRITE HEADER
    ClOps.Writer().header(checklist_type, outp_handle)

########################################################################

# 2. PARSE DATA FROM EMBL-FILE
    try:
        seq_records = SeqIO.parse(open(path_to_embl, "r"), "embl")
    except:
        sys.exit('%s annonex2embl ERROR: Cannot parse `%s`' % ('\n',
            colored(path_to_embl, 'red')))

########################################################################

# 3. CONVERSION TO CHECKLISTS
    for counter, seq_record in enumerate(seq_records): # Looping through records
        try:

    # 3.1. Extraction of charset symbols
            sym_keywrds = ['gene', 'note']
            charset_syms = ClOps.Parser().parse_charset_sym(seq_record, sym_keywrds)

    # 3.2. Conversion to checklist format
            if checklist_type == 'ITS':
                if checklist_type in charset_syms:
                    ClOps.Writer().ITS(seq_record, counter, outp_handle)
                else:
                    sys.exit('%s annonex2embl ERROR: Checklist_type not present in '
                             'charset symbols: `%s`' % ('\n', colored(charset_syms, 'red')))

            elif checklist_type == 'rRNA':
                keyw = ['18S', '28S']
                if any(elem in keyw for elem in charset_syms):
                    ClOps.Writer().rRNA(seq_record, counter, charset_syms,
                                        outp_handle)
                else:
                    sys.exit('%s annonex2embl ERROR: None of the keywords (`%s`) present in '
                             'charset symbols: `%s`' % ('\n', colored(keyw, 'red'), colored(charset_syms, 'red')))

            elif checklist_type == 'trnK_matK':
                keyw = ['trnK', 'matK']
                if any(elem in keyw for elem in charset_syms):
                    ClOps.Writer().trnK_matK(seq_record, counter,
                                             outp_handle)
                else:
                    sys.exit('%s annonex2embl ERROR: None of the keywords (`%s`) present in '
                             'charset symbols: `%s`' % ('\n', colored(keyw, 'red'), colored(charset_syms, 'red')))

            elif checklist_type == 'IGS':
                if charset_syms:
                    ClOps.Writer().IGS(seq_record, counter,
                                       charset_syms, outp_handle)
                else:
                    sys.exit('%s annonex2embl ERROR: To few charset symbols: `%s`' % ('\n', colored(charset_syms, 'red')))

            elif checklist_type == 'genomic_CDS':
                if len(charset_syms)>=2:
                    ClOps.Writer().genomic_CDS(seq_record, counter,
                                           charset_syms, outp_handle)
                else:
                    sys.exit('%s annonex2embl ERROR: To few charset symbols: `%s`' % ('\n', colored(charset_syms, 'red')))

            else:
                sys.exit('%s annonex2embl ERROR: Checklist type `%s` '
                         'not recognized.' % ('\n',
                         colored(checklist_type, 'red')))

        except:
            raise ME.MyException('%s annonex2embl ERROR: Processing of record `%s` failed.' % ('\n', colored(seq_record.name, 'red')))

########################################################################

# 4. CLOSE OUTFILE
    outp_handle.close()

# 5. Return True if no errors occurred
    return True
