#!/usr/bin/env python
'''
Custom operations to generate ENA checklists
'''

# Find all exercises: TOEDIT

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

    def parse_marker_abbrevs(self, seq_record, target_qualifiers):
        ''' This function extracts the marker abbreviations from a sequence record.
        Args:
            seq_record (obj)
            target_qualifiers (list)
        Returns:
            marker_abbrev (list)
        Raises:
            -
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
            marker_abbrev = [elem for elem in marker_abbrev if len(elem.split(" ")) == 1]
        else:
            raise ME.ParserError('ERROR: No marker abbreviation parsed successfully for record `%s`' % (seq_record.id))
        return marker_abbrev


class Writer:
    ''' This class contains functions to write tab-separated
    spreadsheets for submission via the WEBIN checklist submission system.
    Args:
        [specific to function]
    Returns:
        [specific to function]
    Raises:
        -
    '''

    def __init__(self):
        pass

    def deleteEmptyKeys(self, checklist_type, outp_handle):
        ''' This function delete the not necessary keys
            So the ouput is dynamical
        Args:
            keys (list)
            outp_handle (list)
        Returns:
            keys (list)
        Raises:
            -
        '''

        keys = GlobVars.getOutlist(checklist_type)
        toDelete = []
        for i in keys:
            prev = ''
            isEmpty = False
            for j in range(len(outp_handle)):
                if prev == outp_handle[j][i]:
                    prev = outp_handle[j][i]
                    isEmpty = True
                else:
                    isEmpty = False
                    break

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

        keys = self.deleteEmptyKeys(checklist_type, outp_handle)
        out_string = '\t'.join(keys) + '\n'
        outp_file.write(out_string)

        for out_list in outp_handle:
            out_array = []
            for key in keys:
                    out_array.append(out_list[key])
            out_string = '\t'.join(out_array) + '\n'
            outp_file.write(out_string)

    def genomic_CDS(self, seq_record, counter, outdict):
        ''' This function writes a TSV spreadsheet for submission via
            the WEBIN checklist submission system.
        Args:
            seq_record (obj)
            counter (int)
            outp_handle (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''

        # ENTRYNUMBER
        outdict["entrynumber"] = str(counter+1) # enumerate counter starts counting at 0
        # ORGANISM_NAME
        outdict["organism_name"] = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]
        # ENV_SAMPLE
        outdict["env_sam"] = 'no'
        # GENE               # Symbol of the gene corresponding to a sequence region; example: RdRp, sigA, inv
        outdict["gene_symbol"] = "foo bar"
        # PRODUCT            # Name of the product associated with the feature; example: RNA dependent RNA polymerase, sigma factor A
        outdict["product_name"] = "foo bar"
        # TRANSLATION TABLE  # Translation table for this organism. Chose from a drop-down list; example: 1, 2, 3, 5, 11
        outdict["transl_table"] = "12345"

        # the gene
        the_gene = [f for f in seq_record.features
                    if f.type == 'gene']
        try:
            the_gene = the_gene[0]
        except:
            try:
                the_gene = [f for f in seq_record.features if f.type=='CDS']
            except:
                GlobVars.warnings.append('Warning: Problem with `%s`. %s gene not found.' % (seq_name, 'The gene'))

        # 5' CDS LOCATION and 5'_PARTIAL
            # 5' CDS LOCATION   # Start of the coding region relative to the submitted sequence. For a full length CDS this is the position of the first base of the start codon.
        outdict["fiveprime_cds"] = str(the_gene.location.start.position)
        # PARTIAL AT 5'? (yes/no)  # For an incomplete CDS with the start codon upstream of the submitted sequence.
        if type(the_gene.location.start) == Bio.SeqFeature.ExactPosition:
            outdict["fiveprime_partial"] = 'no'
        if type(the_gene.location.start) == Bio.SeqFeature.BeforePosition:
            outdict["fiveprime_partial"] = 'yes'
        # 3' CDS LOCATION and 3'_PARTIAL
            # 3' CDS LOCATION # End of the coding region relative to the submitted sequence. For a full length CDS this is the position of the last base of the stop codon.
        outdict["threeprime_cds"] = str(the_gene.location.end.position)
        # PARTIAL AT 3'? (yes/no) # For an incomplete CDS with the stop codon downstream of the submitted sequence.
        if type(the_gene.location.end) == Bio.SeqFeature.ExactPosition:
            outdict["threeprime_partial"] = 'no'
        if type(the_gene.location.end) == Bio.SeqFeature.AfterPosition:
            outdict["threeprime_partial"] = 'yes'
        # READING FRAME  # Mandatory if your CDS is 5' partial as it defines the reading frame. Location of the first base of the first fully-encoded amino acid., Example: 1,2 or 3
        outdict["read_frame"] = "12345"

        source_qualifiers = [f.qualifiers for f in seq_record.features if f.type=='source'][0]
        # ISOLATE
        try:
            outdict["isolate"] = source_qualifiers['isolate'][0]
        except:
            pass

        # SPEC_VOUCH
        try:
            outdict["spec_vouch"] = source_qualifiers['specimen_voucher'][0]
        except:
            pass

        # LOCALITY
        try:
            outdict["country"] = source_qualifiers['country'][0]
        except:
            pass

        # ECOTYPE
        try:
            outdict["ecotype"] = source_qualifiers['ecotype'][0]
        except:
            pass

        # SEQUENCE
        outdict["sequence"] = str(seq_record.seq)

        return outdict

    def trnK_matK(self, seq_record, counter, outdict):
        ''' This function writes a TSV spreadsheet for submission via
            the WEBIN checklist submission system.
        Args:
            seq_record (obj)
            counter (int)
            outp_handle (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''
        # ENTRYNUMBER
        outdict["entrynumber"] = str(counter+1)  # enumerate counter starts counting at 0
        # ORGANISM_NAME is mandatory
        if [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]:
            outdict["organism_name"] = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]
        else:
            GlobVars.warnings.append('Warning: ' + seq_record.id + ' is missing a mandatory feature: name of the organism')
            # TOEDIT: Right here the program should stop with this seq_record, delete this on and continue

        gene_features = [f for f in seq_record.features if not f.type=='source']
        # trnK_intron
        try:
            trnK_intron = [f for f in gene_features
                           if 'trnK' in f.qualifiers['gene'] and f.type=='intron']
            trnK_intron_present = 'yes'
        except:
            trnK_intron_present = 'no'

        # matK
        try:
            matK_gene = [f for f in gene_features
                         if 'matK' in f.qualifiers['gene'] and
                         (f.type=='gene' or f.type=='CDS')]
            matK_gene = matK_gene[0]
        except:
            GlobVars.warnings.append('Warning: Qualifiers for gene `%s` not found.' % ('\n', 'matK'))

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

        source_qualifiers = [f.qualifiers for f in seq_record.features if f.type=='source'][0]
        # ISOLATE
        try:
            outdict["isolate"] = source_qualifiers['isolate'][0]
        except:
            pass

        # SPEC_VOUCH
        try:
            outdict["spec_vouch"] = source_qualifiers['specimen_voucher'][0]
        except:
            pass

        # LOCALITY
        try:
            outdict["country"] = source_qualifiers['country'][0]
        except:
            pass

        # ECOTYPE
        try:
            outdict["ecotype"] = source_qualifiers['ecotype'][0]
        except:
            pass

        # SEQUENCE
        outdict["sequence"] = str(seq_record.seq)

        return outdict

    def rRNA(self, seq_record, counter, marker_abbrev, outdict):
        ''' This function writes a TSV spreadsheet for submission via
            the WEBIN checklist submission system.
        Args:
            seq_record (obj)
            counter (int)
            marker_abbrev (list)
            outp_handle (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''

        # ENTRYNUMBER
        outdict["entrynumber"] = str(counter+1)  # enumerate counter starts counting at 0
        # ORGANISM_NAME
        outdict["organism_name"] = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]
        # SEDIMENT
        outdict["sediment"] = '_'.join(marker_abbrev)

        source_qualifiers = [f.qualifiers for f in seq_record.features if f.type=='source'][0]
        # ISOLATE
        try:
            outdict["isolate"] = source_qualifiers['isolate'][0]
        except:
            pass

        # ISOLATION_SOURCE
        try:
            outdict["isol_source"] = source_qualifiers['isolation_source']
        except:
            pass

        # COUNTRY
        try:
            outdict["country"] = source_qualifiers['country'][0]
        except:
            pass

        # ECOTYPE
        try:
            outdict["lat_lon"] = source_qualifiers['lat_lon'][0]
        except:
            pass

        # COLLECTION_DATE
        try:
            outdict["collection_date"] = source_qualifiers['collection_date']
        except:
            pass


        # SEQUENCE
        outdict["sequence"] = str(seq_record.seq)

        return outdict

    def ITS(self, seq_record, counter, outdict):
        ''' This function writes a TSV spreadsheet for submission via
            the WEBIN checklist submission system.
        Args:
            seq_record (obj)
            counter (int)
            outp_handle (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''

        # ENTRYNUMBER
        outdict["entrynumber"] = str(counter+1)  # enumerate counter starts counting at 0

        # ORGANISM_NAME
        outdict["organism_name"] = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]

        source_qualifiers = [f.qualifiers for f in seq_record.features if f.type=='source'][0]
        # ISOLATE
        try:
            outdict["isolate"] = source_qualifiers['isolate'][0]
        except:
            pass

        # ENV_SAMPLE
        outdict["env_sam"] = 'no'

        # COUNTRY
        try:
            outdict["country"] = source_qualifiers['country'][0]
        except:
            pass

        # SPEC_VOUCH
        try:
            outdict["spec_vouch"] = source_qualifiers['specimen_voucher'][0]
        except:
            pass

        all_seqrec_features = [f.qualifiers['gene'] for f in seq_record.features]
        # 18S
        if '18S' in all_seqrec_features:
            outdict["RNA_18S"] = 'partial'
        else:
            outdict["RNA_18S"] = 'no'
        # 26S
        if '26S' in all_seqrec_features:
            outdict["RNA_26S"] = 'partial'
        else:
            outdict["RNA_26S"] = 'no'
        # ITS1
        if '18S' in all_seqrec_features:
            outdict["ITS1_feat"] = 'complete'
        else:
            outdict["ITS1_feat"] = 'partial'
        # ITS2
        if '26S' in all_seqrec_features:
            outdict["ITS2_feat"] = 'complete'
        else:
            outdict["ITS2_feat"] = 'partial'
        # 58S
        if 'ITS1' in all_seqrec_features and 'ITS2' in all_seqrec_features:
            outdict["RNA_58S"] = 'complete'
        else:
            outdict["RNA_58S"] = 'partial'

        # SEQUENCE
        outdict["sequence"] = str(seq_record.seq)

        return outdict

    def IGS(self, seq_record, counter, marker_abbrev, outdict):
        ''' This function writes a TSV spreadsheet for submission via
            the WEBIN checklist submission system.
        Args:
            seq_record (obj)
            counter (int)
            marker_abbrev (list)
            outp_handle (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''

        # ENTRYNUMBER
        outdict["entrynumber"] = str(counter+1)  # enumerate counter starts counting at 0
        # ORGANISM_NAME
        outdict["organism_name"] = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]

        # ENV_SAMPLE
        outdict["env_sam"] = 'no'
        # GENE1 and G1PRESENT
        try:
            outdict["gene1"] = marker_abbrev[0]
            outdict["g1present"] = 'yes'
        except:
            #gene1 = 'placeholder'
            outdict["g1present"] = 'no'
        # GENE2 and G2PRESENT
        try:
            outdict["gene2"] = marker_abbrev[1]
            outdict["g2present"] = 'yes'
        except:
            #gene2 = 'placeholder'
            outdict["g2present"] = 'no'

        source_qualifiers = [f.qualifiers for f in seq_record.features if f.type=='source'][0]
        # ISOLATE
        try:
            outdict["isolate"] = source_qualifiers['isolate'][0]
        except:
            pass
        # SPEC_VOUCH
        try:
            outdict["spec_vouch"] = source_qualifiers['specimen_voucher'][0]
        except:
            pass
        # COUNTRY
        try:
            outdict["country"] = source_qualifiers['country'][0]
        except:
            pass

        # SEQUENCE
        outdict["sequence"] = str(seq_record.seq)

        return outdict
