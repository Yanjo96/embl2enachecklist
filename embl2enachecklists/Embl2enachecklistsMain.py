#!/usr/bin/env python2.7

''' Main operations in embl2enachecklist '''

#####################
# IMPORT OPERATIONS #
#####################

import MyExceptions as ME
import ChecklistOps as ClOps
import PrerequisiteOps as PreOps
import globalVariables as GlobVars
import sys
import os

from Bio import SeqIO
from termcolor import colored

# Add specific directory to sys.path in order to import its modules
# NOTE: THIS RELATIVE IMPORTING IS AMATEURISH.
# NOTE: COULD THE FOLLOWING IMPORT BE REPLACED WITH 'import embl2enachecklist'?

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'embl2enachecklists'))

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2018 Michael Gruenstaeudl'
__info__ = 'embl2enachecklists'
__version__ = '2018.06.27.2030'

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



def embl2enachecklists(path_to_embl, path_to_outfile, checklist_type=None):

########################################################################

# 1. OPEN OUTFILE
    outp_handle = []

########################################################################

# 2. PARSE DATA FROM EMBL-FILE
    try:
        seq_records = SeqIO.parse(open(path_to_embl, "r"), "embl")
    except:
        raise ME.FileNotExist(path_to_embl + " does not exists")

########################################################################

# 3. CONVERSION TO CHECKLISTS
    for counter, seq_record in enumerate(seq_records): # Looping through records
        try:
            outdict = {}
    # 3.1. Extraction of marker abbreviations
            target_qualifiers = ['gene','note','standard_name']
            marker_abbrev = ClOps.Parser().parse_marker_abbrevs(seq_record, target_qualifiers)

    # 3.2. Check if marker abbreviation has implemented checklist type if not it skip this seq_record
            if PreOps.checkMinimalPrerequisites(checklist_type, marker_abbrev):
                pass
            else:
                raise ME.MinmalPrerequisitesNotMet('The minimal prerequisites are not met')

    # 3.3. Check if marker abbreviation has implemented checklist type if not it skip this seq_record
            if PreOps.checkFeaturePrerequisites(checklist_type, seq_record.features):
                pass
            else:
                raise ME.FeaturePrerequisutesNotMet('Youre features: "' + ' '.join(seq_record.features) + '" are not part of the allowed Marker abbreviations')

    # 3.4. Check if marker abbreviation has implemented checklist type if not it skip this seq_record
            if not checkCorrectCheckListType(checklist_type, marker_abbrev):
                print checkCorrectCheckListType(checklist_type, marker_abbrev)
                GlobVars.warnings.append('Warning: Checklist type not found as marker abbreviation: `%s`' % (marker_abbrev))
                continue

    # 3.5. Conversion to checklist format
            if checklist_type == 'ITS':
                outdict = ClOps.Writer().ITS(seq_record, counter, GlobVars.getOutdict(checklist_type))

            elif checklist_type == 'rRNA':
                outdict = ClOps.Writer().rRNA(seq_record, counter, marker_abbrev,  GlobVars.getOutdict(checklist_type))

            elif checklist_type == 'trnK_matK':
                outdict = ClOps.Writer().trnK_matK(seq_record, counter, GlobVars.getOutdict(checklist_type))

            elif checklist_type == 'IGS':
                outdict = ClOps.Writer().IGS(seq_record, counter, marker_abbrev, GlobVars.getOutdict(checklist_type))

            elif checklist_type == 'genomic_CDS':
                outdict = ClOps.Writer().genomic_CDS(seq_record, counter, GlobVars.getOutdict(checklist_type))

            outp_handle.append(outdict)

        except ME.ParserError as parserError:
            raise parserError

        except Exception as warning:
            GlobVars.warnings.append('Warning: Processing of record `%s` failed.' % (seq_record.id))
            continue

########################################################################

# 4. CLOSE OUTFILE
    outp_file = open(path_to_outfile,"w")
    ClOps.Writer().writer(checklist_type, outp_handle, outp_file)
    outp_file.close()

# 5. Return True if no errors occurred
    return True
