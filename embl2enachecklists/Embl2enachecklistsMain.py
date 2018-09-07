#!/usr/bin/env python2.7

''' Main operations in embl2enachecklist '''

#####################
# IMPORT OPERATIONS #
#####################

import MyExceptions as ME
import ChecklistOps as ClOps
import PrerequisiteOps as PreOps
import globalVariables as GlobVars
import Tkinter as tk
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
def embl2enachecklists(path_to_embl, path_to_outfile, inputFormat, checklist_type=None, env_sample='no'):
    '''Main functions

    This function does the pipeline job
    It execute all functions of the programm in the right order

    Args:
        path_to_embl (String)
        path_to_outfile (String)
        checklist_type (String)

    Returns:
        bool : True if nithing goes wrong
    Raises:
        Give all the errors whicxh occur to the GUI or CMD function so it can
        drop the error for the user
    '''

    globalvariables = GlobVars.GlobalVariables()
########################################################################

# 1. OPEN OUTFILE
    outp_handle = []

########################################################################

# 2. PARSE DATA FROM EMBL-FILE
    try:
        seq_records = SeqIO.parse(open(path_to_embl, "r"), inputFormat)
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

    # 3.2.1. Check if marker abbreviation has implemented checklist type if not it skip this seq_record
            PreOps.Checker().checkMinimalPrerequisites(checklist_type, marker_abbrev)

    # 3.2.2. Check if marker abbreviation has implemented checklist type if not it skip this seq_record
            PreOps.Checker().checkFeaturePrerequisites(checklist_type, seq_record.features)

    # 3.2.3. Check if marker abbreviation has implemented checklist type if not it skip this seq_record
            if not PreOps.Checker().checkCorrectCheckListType(checklist_type, marker_abbrev):
                GlobVars.warnings.append('Warning: Checklist type not found as marker abbreviation: `%s`' % (marker_abbrev))
                continue

    # 3.3. Conversion to checklist format
            outdict = ClOps.Parser().mandatoryQualifiers(seq_record,
                                                         marker_abbrev,
                                                         counter,
                                                         checklist_type,
                                                         env_sample)


    # 3.4. When outdict wasnt False it will add the optional qualifiers
            if outdict:
                optional_qualifiers = ClOps.Parser().optionalQualifiers(seq_record, globalvariables.getQualifiers(checklist_type,'o'))
                outdict.update(optional_qualifiers)
                outp_handle.append(outdict)

        except Exception as error:
            raise error

########################################################################

# 4. CLOSE OUTFILE
    outp_file = open(path_to_outfile,"w")
    ClOps.Writer().writer(checklist_type, outp_handle, outp_file)
    outp_file.close()

# 5. Return True if no errors occurred
    return True
