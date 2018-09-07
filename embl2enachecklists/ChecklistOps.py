#!/usr/bin/env python
'''
Custom operations to generate ENA checklists
'''

# Find all exercises: TOEDIT



#### ##     ## ########   #######  ########  ########     #######  ########  ######## ########     ###    ######## ####  #######  ##    ##  ######
 ##  ###   ### ##     ## ##     ## ##     ##    ##       ##     ## ##     ## ##       ##     ##   ## ##      ##     ##  ##     ## ###   ## ##    ##
 ##  #### #### ##     ## ##     ## ##     ##    ##       ##     ## ##     ## ##       ##     ##  ##   ##     ##     ##  ##     ## ####  ## ##
 ##  ## ### ## ########  ##     ## ########     ##       ##     ## ########  ######   ########  ##     ##    ##     ##  ##     ## ## ## ##  ######
 ##  ##     ## ##        ##     ## ##   ##      ##       ##     ## ##        ##       ##   ##   #########    ##     ##  ##     ## ##  ####       ##
 ##  ##     ## ##        ##     ## ##    ##     ##       ##     ## ##        ##       ##    ##  ##     ##    ##     ##  ##     ## ##   ### ##    ##
#### ##     ## ##         #######  ##     ##    ##        #######  ##        ######## ##     ## ##     ##    ##    ####  #######  ##    ##  ######

import Bio
import globalVariables as GlobVars
import MyExceptions as ME



   ###    ##     ## ######## ##     ##  #######  ########     #### ##    ## ########  #######
  ## ##   ##     ##    ##    ##     ## ##     ## ##     ##     ##  ###   ## ##       ##     ##
 ##   ##  ##     ##    ##    ##     ## ##     ## ##     ##     ##  ####  ## ##       ##     ##
##     ## ##     ##    ##    ######### ##     ## ########      ##  ## ## ## ######   ##     ##
######### ##     ##    ##    ##     ## ##     ## ##   ##       ##  ##  #### ##       ##     ##
##     ## ##     ##    ##    ##     ## ##     ## ##    ##      ##  ##   ### ##       ##     ##
##     ##  #######     ##    ##     ##  #######  ##     ##    #### ##    ## ##        #######

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2018 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2018.06.27.2030'



########  ######## ########  ##     ##  ######    ######   #### ##    ##  ######
##     ## ##       ##     ## ##     ## ##    ##  ##    ##   ##  ###   ## ##    ##
##     ## ##       ##     ## ##     ## ##        ##         ##  ####  ## ##
##     ## ######   ########  ##     ## ##   #### ##   ####  ##  ## ## ## ##   ####
##     ## ##       ##     ## ##     ## ##    ##  ##    ##   ##  ##  #### ##    ##
##     ## ##       ##     ## ##     ## ##    ##  ##    ##   ##  ##   ### ##    ##
########  ######## ########   #######   ######    ######   #### ##    ##  ######

import pdb
# pdb.set_trace()



 ######  ##          ###     ######   ######  ########  ######
##    ## ##         ## ##   ##    ## ##    ## ##       ##    ##
##       ##        ##   ##  ##       ##       ##       ##
##       ##       ##     ##  ######   ######  ######    ######
##       ##       #########       ##       ## ##             ##
##    ## ##       ##     ## ##    ## ##    ## ##       ##    ##
 ######  ######## ##     ##  ######   ######  ########  ######

class Writer:
    ''' This class contains functions to write tab-separated
    spreadsheets for submission via the WEBIN checklist submission system.
    Args:
        [specific to function]
    '''

    def __init__(self):
        pass

    def deleteEmptyKeys(self, checklist_type, outp_handle):
        ''' This function delete the not necessary keys
            So the output is dynamical
        Args:
            keys (list)
            outp_handle (list)
        Returns:
            keys (list)
        Raises:
            -
        '''

        keys = GlobVars.GlobalVariables().getQualifiers(checklist_type, 'om')
        toDelete = []
        for i in keys:
            prev = ''
            isEmpty = False
            for j in range(len(outp_handle)):
                try:
                    if prev == outp_handle[j][i]:
                        prev = outp_handle[j][i]
                        isEmpty = True
                    else:
                        isEmpty = False
                        break
                except:
                    isEmpty = True

            if isEmpty:
                toDelete.append(i)

        for delete in toDelete:
            j = 0
            for i in range(len(keys)):
                if keys[j] == delete:
                    del keys[j]
                else:
                    j = j + 1

        return keys



    def writer(self, checklist_type, outp_handle, outp_file):
        ''' This function writes a TSV spreadsheet for submission via
            the WEBIN checklist submission system.
        Args:
            checklist_type (string)
            outp_handle (list)
            outp_file (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''

        translator = GlobVars.GlobalVariables().getTranslator(checklist_type)
        keys = self.deleteEmptyKeys(checklist_type, outp_handle)
        for out_list in outp_handle:
            for key in keys:
                try:
                    out_list[translator[key]] = out_list.pop(key)
                except:
                    pass

        for i in range(len(keys)):
            keys[i] = translator[keys[i]]

        out_string = '\t'.join(keys) + '\n'
        outp_file.write(out_string)

        for out_list in outp_handle:
            out_array = []
            for key in keys:
                    try:
                        out_array.append(out_list[key])
                    except:
                        out_array.append('')

            out_string = '\t'.join(out_array) + '\n'
            outp_file.write(out_string)


class Parser:
    ''' This class contains functions to parse out information from
    sequence records.
    Args:
        -
    '''

    def __init__(self):
        pass

     #####    ##   #####   ####  ######         #    #   ##   #####  #    # ###### #####            ##   #####  #####  #####  ###### #    #  ####
     #    #  #  #  #    # #      #              ##  ##  #  #  #    # #   #  #      #    #          #  #  #    # #    # #    # #      #    # #
     #    # #    # #    #  ####  #####          # ## # #    # #    # ####   #####  #    #         #    # #####  #####  #    # #####  #    #  ####
     #####  ###### #####       # #              #    # ###### #####  #  #   #      #####          ###### #    # #    # #####  #      #    #      #
     #      #    # #   #  #    # #              #    # #    # #   #  #   #  #      #   #          #    # #    # #    # #   #  #       #  #  #    #
     #      #    # #    #  ####  ######         #    # #    # #    # #    # ###### #    #         #    # #####  #####  #    # ######   ##    ####
                                        #######                                           #######

    def parse_marker_abbrevs(self, seq_record, target_qualifiers):
        ''' This function extracts the marker abbreviations from a sequence record.
        Args:
            seq_record (SeqIO)
            target_qualifiers (list)
        Returns:
            marker_abbrev (list)
        Raises:
            ParserError: When there are no marker abbreviations it would raise an error
        '''
        gene_qualifiers = [f.qualifiers for f in seq_record.features if not f.type=='source'] # Produces a list of dictionaries
        marker_abbrev = []
        # Extract marker abbreviations
        for keyw in target_qualifiers:
            for dct in gene_qualifiers:
                try:
                    marker_abbrev.extend(dct[keyw])
                except KeyError:
                    marker_abbrev = marker_abbrev
        # Parse raw list of marker abbreviations
        if marker_abbrev:
            # 3.1.1. Extract all unique values in list, keep order of original list
            seen = set()
            marker_abbrev = [elem for elem in marker_abbrev if elem not in seen and not seen.add(elem)]
            # 3.1.2. Remove any multi-word elements
            #marker_abbrev = [elem for elem in marker_abbrev if len(elem.split(" ")) == 1]
        else:
            raise ME.ParserError('ERROR: No marker abbreviation parsed successfully for record `%s`' % (seq_record.id))
        return marker_abbrev

                                                                   #####
     #    #   ##   #    # #####    ##   #####  ####  #####  #   # #     # #    #   ##   #      # ###### # ###### #####   ####
     ##  ##  #  #  ##   # #    #  #  #    #   #    # #    #  # #  #     # #    #  #  #  #      # #      # #      #    # #
     # ## # #    # # #  # #    # #    #   #   #    # #    #   #   #     # #    # #    # #      # #####  # #####  #    #  ####
     #    # ###### #  # # #    # ######   #   #    # #####    #   #   # # #    # ###### #      # #      # #      #####       #
     #    # #    # #   ## #    # #    #   #   #    # #   #    #   #    #  #    # #    # #      # #      # #      #   #  #    #
     #    # #    # #    # #####  #    #   #    ####  #    #   #    #### #  ####  #    # ###### # #      # ###### #    #  ####

    def mandatoryQualifiers(self, seq_record, marker_abbrev, counter, checklist_type, env_sample):
        '''fill the mandatory qualifiers for specific checklist_type

        At first it fill the mandatory qualifiers where part of each checklist
        After that it go further for the specific manatory qualifiers and fill the outdict
        with information. It fill the warnings if it cant fill an information to the
        dictionary cause of it is mandatory at all qualifiers have to be filled

        Args:
            -
        Returns:
            outdict (dict) : if information can be filled
            False          : when something went wrong so the main function can skip
                             this seq_record
        Raises:
            it will not raise any errors cause it only give warnings for the specific
            seq_record
        '''

        outdict = GlobVars.GlobalVariables().getOutdict(GlobVars.GlobalVariables().getQualifiers(checklist_type, 'm'))
                                                         #####  ######   #####
          ####  ###### #    #  ####  #    # #  ####     #     # #     # #     #
         #    # #      ##   # #    # ##  ## # #    #    #       #     # #
         #      #####  # #  # #    # # ## # # #         #       #     #  #####
         #  ### #      #  # # #    # #    # # #         #       #     #       #
         #    # #      #   ## #    # #    # # #    #    #     # #     # #     #
          ####  ###### #    #  ####  #    # #  ####      #####  ######   #####

        if checklist_type == 'genomic_CDS':
            # ENTRYNUMBER
            outdict["entrynumber"] = str(counter+1)  # enumerate counter starts counting at 0

            # ORGANISM_NAME
            try:
                outdict["organism"] = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('organism', seq_record.id))
                return False


            # SEQUENCE
            try:
                outdict["sequence"] = str(seq_record.seq)
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('sequence', seq_record.id))
                return False

          ####  ###### #    # ######    # #    # ##### #####   ####  #    #
         #    # #      ##   # #         # ##   #   #   #    # #    # ##   #
         #      #####  # #  # #####     # # #  #   #   #    # #    # # #  #
         #  ### #      #  # # #         # #  # #   #   #####  #    # #  # #
         #    # #      #   ## #         # #   ##   #   #   #  #    # #   ##
          ####  ###### #    # ######    # #    #   #   #    #  ####  #    #

        if checklist_type == 'gene_intron':
            # ENTRYNUMBER
            outdict["entrynumber"] = str(counter+1)  # enumerate counter starts counting at 0
            # ORGANISM_NAME
            try:
                outdict["organism"] = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('organism', seq_record.id))
                return False

            # ENV_SAMPLE
            try:
                outdict["env_sample"] = env_sample
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('env_sample', seq_record.id))
                return False

            # GENE
            try:
                outdict["gene"] = [f.qualifiers['gene'] for f in seq_record.features if f.type=='gene'][0][0]
            except:
                try:
                    outdict["gene"] = [f.qualifiers['gene'] for f in seq_record.features if f.type=='CDS'][0][0]
                except:
                    GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('gene', seq_record.id))
                    return False

            try:
                intron = [f for f in seq_record.features if f.type=='intron'][0]
            except:
                return False

            # 5' PARTIAL
            if type(intron.location.start) == Bio.SeqFeature.ExactPosition:
                outdict["fiveprime_partial"] = 'no'
            if type(intron.location.start) == Bio.SeqFeature.BeforePosition:
                outdict["fiveprime_partial"] = 'yes'

            # 3' PARTIAL
            if type(intron.location.end) == Bio.SeqFeature.ExactPosition:
                outdict["threeprime_partial"] = 'no'
            if type(intron.location.end) == Bio.SeqFeature.AfterPosition:
                outdict["threeprime_partial"] = 'yes'

            # 5' INTRON
            outdict["fiveprime_cds"] = str(intron.location.start.position+1)

            # 3' INTRON
            outdict["threeprime_cds"] = str(intron.location.end.position)

            # NUMBER
            try:
                outdict["number"] = [f.qualifiers['number'] for f in seq_record.features if f.type=='intron'][0][0]
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('intron', seq_record.id))
                return False

            # SEQUENCE
            try:
                outdict["sequence"] = str(seq_record.seq)
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('sequence', seq_record.id))
                return False


                             #    #                        #    #
         ##### #####  #    # #   #     #    #   ##   ##### #   #
           #   #    # ##   # #  #      ##  ##  #  #    #   #  #
           #   #    # # #  # ###       # ## # #    #   #   ###
           #   #####  #  # # #  #      #    # ######   #   #  #
           #   #   #  #   ## #   #     #    # #    #   #   #   #
           #   #    # #    # #    #    #    # #    #   #   #    #

        elif checklist_type == 'trnK_matK':
            # ENTRYNUMBER
            outdict["entrynumber"] = str(counter+1)  # enumerate counter starts counting at 0

            # ORGANISM_NAME
            try:
                outdict["organism"] = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('organism', seq_record.id))
                return False

            try:
                gene_features = [f for f in seq_record.features if not f.type=='source']
                # trnK_intron
                try:
                    trnK_intron_present = 'no'
                    for t in gene_features:
                        if t.type == "intron":
                            if 'trnK' in t.qualifiers['gene']:
                                trnK_intron_present = 'yes'
                    outdict['trnK_intron_present'] = trnK_intron_present
                except:
                    outdict['trnK_intron_present'] = 'no'
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('trnK_intron_present', seq_record.id))
                return False

            # matK
            matK_gene = None
            try:
                for gene in [f for f in gene_features if f.type=="gene"]:
                    if "gene" in gene.qualifiers:
                        if 'matK' in gene.qualifiers['gene']:
                            matK_gene = gene
                    if "note" in gene.qualifiers:
                        if 'matK' in gene.qualifiers['note']:
                            matK_gene = gene
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('matK', seq_record.id))
                return False

            # 5'_CDS and 5'_PARTIAL
            # 5'_CDS: Start of the matK coding region relative to the submitted sequence. For a full length CDS this is the position of the first base of the start codon.
            # NOTE: One nucleotide position has to be added to the start position to make it correct.

            outdict["fiveprime_cds"] = str(matK_gene.location.start.position+1)
            # 5'_PARTIAL: cds partial at 5'? (yes/no) For an incomplete CDS with the start codon upstream of the submitted sequence.
            if type(matK_gene.location.start) == Bio.SeqFeature.ExactPosition:
                outdict["fiveprime_partial"] = 'no'
            if type(matK_gene.location.start) == Bio.SeqFeature.BeforePosition:
                outdict["fiveprime_partial"] = 'yes'

            # 3'_CDS and 3'_PARTIAL
            # 3'_CDS: End of the matK coding region relative to the submitted sequence. For a full length CDS this is the position of the last base of the stop codon.

            outdict["threeprime_cds"] = str(matK_gene.location.end.position)
            # 3'_PARTIAL: cds partial at 3'? (yes/no) For an incomplete CDS with the stop codon downstream of the submitted sequence.
            if type(matK_gene.location.end) == Bio.SeqFeature.ExactPosition:
                outdict["threeprime_partial"] = 'no'
            if type(matK_gene.location.end) == Bio.SeqFeature.AfterPosition:
                outdict["threeprime_partial"] = 'yes'

            # SEQUENCE
            try:
                outdict["sequence"] = str(seq_record.seq)
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('sequence', seq_record.id))
                return False

                ######  #     #    #
         #####  #     # ##    #   # #
         #    # #     # # #   #  #   #
         #    # ######  #  #  # #     #
         #####  #   #   #   # # #######
         #   #  #    #  #    ## #     #
         #    # #     # #     # #     #

        elif checklist_type == 'rRNA':

            # ENTRYNUMBER
            outdict["entrynumber"] = str(counter+1)  # enumerate counter starts counting at 0

            # ORGANISM_NAME
            try:
                outdict["organism"] = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('organism', seq_record.id))
                return False

            # ENV_SAMPLE
            outdict["env_sample"] = env_sample

            # SEDIMENT
            # sediment can also occur in the qualifiers: sediment, gene and note
            try:
                outdict["sediment"] = [f.qualifiers['product'] for f in seq_record.features if f.type == 'rRNA'][0][0]
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('sediment', seq_record.id))
                return False

            # SEQUENCE
            try:
                outdict["sequence"] = str(seq_record.seq)
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('sequence', seq_record.id))
                return False

         ### #######  #####
          #     #    #     #
          #     #    #
          #     #     #####
          #     #          #
          #     #    #     #
         ###    #     #####

        elif checklist_type == 'ITS':

            # ENTRYNUMBER
            outdict["entrynumber"] = str(counter+1)  # enumerate counter starts counting at 0

            # ORGANISM_NAME
            try:
                outdict["organism"] = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('organism', seq_record.id))
                return False

            # ENV_SAMPLE
            try:
                outdict["env_sample"] = env_sample
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('env_sam', seq_record.id))
                return False

            # ISOLATION_SOURCE
            if env_sample == 'yes':
                try:
                    outdict["isolation_source"] = [f.qualifiers['isolation_source'] for f in seq_record.features if f.type=='source'][0][0]
                except:
                    GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('isolation_source', seq_record.id))
                    return False

            # 18S, ITS1, ITS2, 28S
            # NACHFRAGEN
            try:
                gene_seqrec_features = " ".join([f.qualifiers['gene'] for f in seq_record.features if f.type=='rRNA'][0])
                # 18S ITS1
                if '18S' in gene_seqrec_features:
                    outdict["RNA_18S"] = 'partial'
                    outdict["ITS1_feat"] = 'complete'
                else:
                    outdict["RNA_18S"] = 'no'
                    outdict["ITS1_feat"] = 'partial'
                # 28S ITS2
                if '28S' in gene_seqrec_features or '26S' in gene_seqrec_features:
                    outdict["RNA_28S"] = 'partial'
                    outdict["ITS2_feat"] = 'complete'
                else:
                    outdict["RNA_28S"] = 'no'
                    outdict["ITS2_feat"] = 'partial'
            except:
                outdict["RNA_18S"] = 'no'
                outdict["ITS1_feat"] = 'partial'
                outdict["RNA_28S"] = 'no'
                outdict["ITS2_feat"] = 'partial'

            # 58S
            try:
                note_seqrec_features = " ".join([f.qualifiers['note'] for f in seq_record.features if f.type=='misc_feature'][0])
                if 'ITS1' in note_seqrec_features and 'ITS2' in note_seqrec_features:
                    outdict["RNA_58S"] = 'complete'
                else:
                    outdict["RNA_58S"] = 'partial'
            except:
                outdict["RNA_58S"] = 'partial'

            # SEQUENCE
            try:
                outdict["sequence"] = str(seq_record.seq)
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('sequence', seq_record.id))
                return False

         ###  #####   #####
          #  #     # #     #
          #  #       #
          #  #  ####  #####
          #  #     #       #
          #  #     # #     #
         ###  #####   #####

        elif checklist_type == 'IGS':
            # ENTRYNUMBER
            outdict["entrynumber"] = str(counter+1)  # enumerate counter starts counting at 0

            # ORGANISM_NAME
            try:
                outdict["organism"] = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('organism', seq_record.id))
                return False

            # ENV_SAMPLE
            try:
                outdict["env_sample"] = env_sample
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('env_sam', seq_record.id))
                return False

            # GENE1 and G1PRESENT
            try:
                try:
                    outdict["gene1"] = [f.qualifiers['gene'] for f in seq_record.features if f.type=='gene'][0][0]
                    outdict["g1present"] = 'yes'
                except:
                    outdict['gene1'] = 'placeholder'
                    outdict["g1present"] = 'no'
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('gene1, g1present', seq_record.id))
                return False

            # GENE2 and G2PRESENT
            try:
                try:
                    outdict["gene2"] = [f.qualifiers['gene'] for f in seq_record.features if f.type=='gene'][1][0]
                    outdict["g2present"] = 'yes'
                except:
                    outdict['gene2'] = 'placeholder'
                    outdict["g2present"] = 'no'
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('gene2, g2present', seq_record.id))
                return False

            # SEQUENCE
            try:
                outdict["sequence"] = str(seq_record.seq)
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('sequence', seq_record.id))
                return False

         ####### #######  #####
         #          #    #     #
         #          #    #
         #####      #     #####
         #          #          #
         #          #    #     #
         #######    #     #####

        elif checklist_type == 'ETS':
            # ENTRYNUMBER
            outdict["entrynumber"] = str(counter+1)  # enumerate counter starts counting at 0

            # ORGANISM_NAME
            try:
                outdict["organism"] = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('organism', seq_record.id))
                return False

            #TOEDIT
            try:
                outdict["ets_type"] = "3'"
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('ets_type', seq_record.id))
                return False

            # SEQUENCE
            try:
                outdict["sequence"] = str(seq_record.seq)
            except:
                GlobVars.warnings.append('Warning: The mandatory feature "%s" in "%s" is missing.' % ('sequence', seq_record.id))
                return False

         ######  ####### #     # #######
         #     # #     # ##    # #
         #     # #     # # #   # #
         #     # #     # #  #  # #####
         #     # #     # #   # # #
         #     # #     # #    ## #
         ######  ####### #     # #######

        return outdict


                                                          #######
      ####  #####  ##### #  ####  #    #   ##   #         #       ######   ##   ##### #    # #####  ######  ####
     #    # #    #   #   # #    # ##   #  #  #  #         #       #       #  #    #   #    # #    # #      #
     #    # #    #   #   # #    # # #  # #    # #         #####   #####  #    #   #   #    # #    # #####   ####
     #    # #####    #   # #    # #  # # ###### #         #       #      ######   #   #    # #####  #           #
     #    # #        #   # #    # #   ## #    # #         #       #      #    #   #   #    # #   #  #      #    #
      ####  #        #   #  ####  #    # #    # ######    #       ###### #    #   #    ####  #    # ######  ####

    def optionalQualifiers(self, seq_record, optional_qualifiers):
        '''fill the optional qualifiers

        Will add all information in the qualifiers to the outdict
        if they are part of the optional qualifiersfor the specific checklist_type

        Args:
            seq_record (SeqIO)
            optional_qualifiers (list)

        Returns:
            outdict (dict) : if information can be filled
            False          : when something went wrong so the main function can skip
                             this seq_record
        Raises:
            it will not raise any errors cause it only give warnings for the specific
            seq_record
        '''
        optional_dictionary = {}
        for i in seq_record.features:
            for j in [quali for quali in i.qualifiers if quali in optional_qualifiers]:
                optional_dictionary.update({j:i.qualifiers[j][0]})

            # if statement fue pcr_primers erstmal ignorieren

        if 'breed' in optional_qualifiers:
            try:
                notes = [f.qualifiers['note'] for f in seq_record.features if f.type=='source'][0][0].split(';')
                for note in notes:
                    if 'breed' in note:
                        optional_dictionary.update({'breed':note.split('breed:')[-1].strip()})
                        break
            except:
                pass

        # PCR primer
        try:
            pcrPrimers = [ i.split(",]") for i in [f.qualifiers['PCR_primers'] for f in seq_record.features if f.type=='source'][0][0].replace(" ","").split(",[") ]
            tmp = []
            for t in pcrPrimers:
                for i in t:
                    tmp.append(i.split(":"))
            pcrPrimers = tmp
            print pcrPrimers
            fwd_name = []
            fwd_seq = []
            rev_name = []
            rev_seq = []
            for i in pcrPrimers:
                if i[0].replace("[","").replace("]","") == 'fwd_name':
                    fwd_name.append(i[1])
                if i[0].replace("[","").replace("]","") == 'fwd_seq':
                    fwd_seq.append(i[1])
                if i[0].replace("[","").replace("]","") == 'rev_name':
                    rev_name.append(i[1])
                if i[0].replace("[","").replace("]","") == 'rev_seq':
                    rev_seq.append(i[1])


            for counter, elem in enumerate(fwd_name):
                optional_dictionary.update({'fwd_name'+str(counter+1):elem})

            for counter, elem in enumerate(fwd_seq):
                optional_dictionary.update({'fwd_seq'+str(counter+1):elem})

            for counter, elem in enumerate(rev_name):
                optional_dictionary.update({'rev_name'+str(counter+1):elem})

            for counter, elem in enumerate(rev_seq):
                optional_dictionary.update({'rev_seq'+str(counter+1):elem})

        except:
            pass

        return optional_dictionary
