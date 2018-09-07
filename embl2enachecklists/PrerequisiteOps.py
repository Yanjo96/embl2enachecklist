#!/usr/bin/env python
'''
Custom operations to generate ENA checklists
'''

 ###                                      #######
  #  #    # #####   ####  #####  #####    #     # #####  ###### #####    ##   ##### #  ####  #    #  ####
  #  ##  ## #    # #    # #    #   #      #     # #    # #      #    #  #  #    #   # #    # ##   # #
  #  # ## # #    # #    # #    #   #      #     # #    # #####  #    # #    #   #   # #    # # #  #  ####
  #  #    # #####  #    # #####    #      #     # #####  #      #####  ######   #   # #    # #  # #      #
  #  #    # #      #    # #   #    #      #     # #      #      #   #  #    #   #   # #    # #   ## #    #
 ### #    # #       ####  #    #   #      ####### #      ###### #    # #    #   #   #  ####  #    #  ####

import Bio
import globalVariables as GlobVars
import MyExceptions as ME

    #    #     # ####### #     # ####### ######     ### #     # ####### #######
   # #   #     #    #    #     # #     # #     #     #  ##    # #       #     #
  #   #  #     #    #    #     # #     # #     #     #  # #   # #       #     #
 #     # #     #    #    ####### #     # ######      #  #  #  # #####   #     #
 ####### #     #    #    #     # #     # #   #       #  #   # # #       #     #
 #     # #     #    #    #     # #     # #    #      #  #    ## #       #     #
 #     #  #####     #    #     # ####### #     #    ### #     # #       #######

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2018 Michael Gruenstaeudl'
__info__ = 'embl2enachecklist'
__version__ = '2018.07.17'


 ######  ####### ######  #     #  #####   #####  ### #     #  #####
 #     # #       #     # #     # #     # #     #  #  ##    # #     #
 #     # #       #     # #     # #       #        #  # #   # #
 #     # #####   ######  #     # #  #### #  ####  #  #  #  # #  ####
 #     # #       #     # #     # #     # #     #  #  #   # # #     #
 #     # #       #     # #     # #     # #     #  #  #    ## #     #
 ######  ####### ######   #####   #####   #####  ### #     #  #####

import pdb
# pdb.set_trace()

  #####  #          #     #####   #####  #######  #####
 #     # #         # #   #     # #     # #       #     #
 #       #        #   #  #       #       #       #
 #       #       #     #  #####   #####  #####    #####
 #       #       #######       #       # #             #
 #     # #       #     # #     # #     # #       #     #
  #####  ####### #     #  #####   #####  #######  #####

class Checker:
    """
    This class contains functions to check the marker abbreviations from a sequence record.

    Attributes:
        -
    """

    def __init__(self):
        pass


     #######
     #       #    # #    #  ####  ##### #  ####  #    #  ####
     #       #    # ##   # #    #   #   # #    # ##   # #
     #####   #    # # #  # #        #   # #    # # #  #  ####
     #       #    # #  # # #        #   # #    # #  # #      #
     #       #    # #   ## #    #   #   # #    # #   ## #    #
     #        ####  #    #  ####    #   #  ####  #    #  ####

    def checkMinimalPrerequisites(self, checklist_type, marker_abbrev):
        """Check minimal prerequisites

        Function proofs if the choosen checklist type is implemented and
        proof if the checklist type is correct for the specific

        Args:
            checklist_type [string]: user input which contains the type of the checklist
            marker_abbrev  [list]  : list of the abbreviations of the sequence record

        Returns:
            -

        Raises:
            ME.MinimalPrerequisitesNotMet: If the checklist type is not part of implemented checklist types
        """
        if not checklist_type in GlobVars.allowed_checklists:
            raise ME.MinimalPrerequisitesNotMet(checklist_type + ' is not an implemented checklist type')
        #NACHFRAGEN wegen IGS marker abbreviation sind nicht fest
        if not any(elem in " ".join(marker_abbrev) for elem in GlobVars.allowed_marker_abbrev) and checklist_type != 'IGS':
            raise ME.MinimalPrerequisitesNotMet('The minimal prerequisites are not met')

    def checkFeaturePrerequisites(self, checklist_type, features):
        """Check feature prerequisites

        Are the marker abbreviation part of the allowed marker abbreviations
        Args:
            checklist_type [string]: user input which contains the type of the checklist
            marker_abbrev  [list]  : list of the abbreviations of the sequence record

        Returns:
            -

        Raises:
            ME.FeaturePrerequisutesNotMet: If the marker abbreviations are not part of the allowed marker abbreviations
        """
        if checklist_type == 'genomic_CDS':
            types = [elem.type for elem in features]
            if not any(elem in types for elem in GlobVars.genomic_CDS_marker_abbrev):
                raise ME.FeaturePrerequisutesNotMet('Youre features: "' + ', '.join(types) + '" are not part of the allowed Marker abbreviations')


    def checkCorrectCheckListType(self, checklist_type, marker_abbrev):
        """Check if the checklist type matched to the marker abbreviations

        Args:
            checklist_type [string]: user input which contains the type of the checklist
            marker_abbrev  [list]  : list of the abbreviations of the sequence record

        Returns:
            Boolean: True if it matched False if not
            Maybe switch to raises so the code in the main class become more readable

        Raises:
            -
        """
        if checklist_type == 'ETS' and 'ETS' in marker_abbrev:
            return True
        if checklist_type == 'ITS' and any(elem in " ".join(marker_abbrev) for elem in GlobVars.allowed_its_marker_abbrev):
            return True
        if checklist_type == 'rRNA' and any(elem in " ".join(marker_abbrev) for elem in GlobVars.allowed_rrna_marker_abbrev):
            return True
        if checklist_type == 'trnK_matK' and any(elem in " ".join(marker_abbrev) for elem in GlobVars.allowed_marker_abbrev):
            return True
        if checklist_type == 'IGS' and marker_abbrev:
            return True
        if checklist_type == 'genomic_CDS' and len(marker_abbrev)>=2:
            return True
        if checklist_type == 'gene_intron':
            return True
        return False
