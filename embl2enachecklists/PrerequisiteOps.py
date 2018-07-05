#!/usr/bin/env python
'''
Custom operations to generate ENA checklists
'''

#####################
# IMPORT OPERATIONS #
#####################

import Bio
import globalVariables as GlobVars
import MyExceptions as ME

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2018 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2018.06.27.2030'

#############
# DEBUGGING #
#############

import pdb
# pdb.set_trace()

###########
# CLASSES #
###########


class Parser:
    ''' This class contains functions to parse out information from
    sequence records.
    Args:
        [specific to function]
    Returns:
        [specific to function]
    Raises:
        -
    '''

    def __init__(self):
        pass


    def checkMinimalPrerequisites(self, checklist_type, marker_abbrev):
        if checklist_type in GlobVars.allowed_checklists:
            return True
        else:
            raise ME.MinmalPrerequisitesNotMet('The minimal prerequisites are not met')
        if any(elem in GlobVars.allowed_marker_abbrev for elem in marker_abbrev):
            return True
        else:
            return False

    def checkFeaturePrerequisites(self, checklist_type, features):
            if checklist_type == 'genomic_CDS':
                if any(elem in features for elem in GlobVars.genomic_CDS_marker_abbrev):
                    return True
                else:
                    return False

    def checkCorrectCheckListType(self, checklist_type, marker_abbrev):
        if checklist_type == 'ETS' and 'ETS' in marker_abbrev:
            return True
        if checklist_type == 'ITS' and 'ITS' in marker_abbrev:
            return True
        if checklist_type == 'rRNA' and any(elem in GlobVars.allowed_rrna_marker_abbrev for elem in marker_abbrev):
            return True
        if checklist_type == 'trnK_matK' and any(elem in GlobVars.allowed_marker_abbrev for elem in marker_abbrev):
            return True
        if checklist_type == 'IGS' and marker_abbrev:
            return True
        if checklist_type == 'genomic_CDS' and len(marker_abbrev)>=2:
            return True
        return False
