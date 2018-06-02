#!/usr/bin/env python
'''
Custom operations to generate ENA checklists
'''

#####################
# IMPORT OPERATIONS #
#####################

import Bio
import MyExceptions as ME

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2018 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2018.05.23.2000'

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

    def parse_charset_sym(self, seq_record, sym_keywrds):
        ''' This function extracts the charset symbols from a sequence record.
        Args:
            seq_record (obj)
            sym_keywrds (list)
        Returns:
            charset_syms (list)
        Raises:
            -
        '''
        gene_quals = [f.qualifiers for f in seq_record.features if not f.type=='source'] # Produces a list of dictionaries
        charset_syms = []
        for keyw in sym_keywrds:
            for dct in gene_quals:
                try:
                    charset_syms.extend(dct[keyw])
                except KeyError:
                    charset_syms = charset_syms
        if charset_syms:
            # 3.1.1. Extract all unique values in list, keep order of original list
            seen = set()
            charset_syms = [e for e in charset_syms if e not in seen and not seen.add(e)]
            # 3.1.2. Remove any multi-word elements
            charset_syms = [e for e in charset_syms if len(e.split(" "))==1]
        if not charset_syms:
            sys.exit('%s annonex2embl ERROR: Parsing of charset symbol '
                     'unsuccessful')
        return charset_syms



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

    def header(self, checklist_type, outp_handle):
        out_list = []

        if checklist_type == 'ITS':
            out_list = ['entrynumber',
                        'organism_name',
                        'isolate',
                        'env_sam',
                        'country',
                        'spec_vouch',
                        'RNA_18S',
                        'ITS1_feat',
                        'RNA_58S',
                        'ITS2_feat',
                        'RNA_26S',
                        'sequence'
                        ]

        elif checklist_type == 'rRNA':
            out_list = ['entrynumber',
                        'organism_name',
                        'sediment',
                        'isolate',
                        'isol_source',
                        'country',
                        'lat_lon',
                        'collection_date',
                        'sequence'
                        ]

        elif checklist_type == 'trnK_matK':
            out_list = ['entrynumber',
                        'organism_name',
                        'fiveprime_cds',
                        'threeprime_cds',
                        'fiveprime_partial',
                        'threeprime_partial',
                        'trnK_intron_present',
                        'isolate',
                        'spec_vouch',
                        'country',
                        'ecotype',
                        'sequence'
                        ]

        elif checklist_type == 'IGS':
            out_list = ['entrynumber',
                        'organism_name',
                        'env_sam',
                        'gene1',
                        'g1present',
                        'gene2',
                        'g2present',
                        'isolate',
                        'spec_vouch',
                        'country',
                        'sequence'
                        ]

        elif checklist_type == 'genomic_CDS':
            out_list = ['entrynumber',
                        'organism_name',
                        'env_sam',
                        'gene_symbol',
                        'product_name',
                        'transl_table',
                        'fiveprime_cds',
                        'threeprime_cds',
                        'fiveprime_partial',
                        'threeprime_partial',
                        'read_frame',
                        'isolate',
                        'spec_vouch',
                        'country',
                        'ecotype',
                        'sequence'
                        ]

        else:
            pass

        out_string = '\t'.join(out_list) + '\n'
        outp_handle.write(out_string)

    def genomic_CDS(self, seq_record, counter, charset_syms, outp_handle):
        ''' This function writes a TSV spreadsheet for submission via
            the WEBIN checklist submission system.
        Args:
            seq_record (obj)
            counter (int)
            charset_syms (list)
            outp_handle (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''

        # ENTRYNUMBER
        entrynumber = str(counter+1)  # enumerate counter starts counting at 0
        # ORGANISM_NAME
        organism_name = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]
        # ENV_SAMPLE
        env_sam = 'no'
        # GENE               # Symbol of the gene corresponding to a sequence region; example: RdRp, sigA, inv
        gene_symbol = "foo bar"
        # PRODUCT            # Name of the product associated with the feature; example: RNA dependent RNA polymerase, sigma factor A
        product_name = "foo bar"
        # TRANSLATION TABLE  # Translation table for this organism. Chose from a drop-down list; example: 1, 2, 3, 5, 11
        transl_table = "12345"

        # the gene
        the_gene = [f for f in seq_record.features
                    if f.type == 'gene']
        try:
            the_gene = the_gene[0]
        except:
            try:
                the_gene = [f for f in seq_record.features
                            if f.type=='CDS']
            except:
                raise ME.MyException('%s annonex2embl ERROR: Problem \
                    with `%s`. %s gene not found.' % ('\n', seq_name, 'The gene'))

        # 5' CDS LOCATION and 5'_PARTIAL
            # 5' CDS LOCATION   # Start of the coding region relative to the submitted sequence. For a full length CDS this is the position of the first base of the start codon.
        fiveprime_cds = str(the_gene.location.start.position)
        # PARTIAL AT 5'? (yes/no)  # For an incomplete CDS with the start codon upstream of the submitted sequence.
        if type(the_gene.location.start) == Bio.SeqFeature.ExactPosition:
            fiveprime_partial = 'no'
        if type(the_gene.location.start) == Bio.SeqFeature.BeforePosition:
            fiveprime_partial = 'yes'
        # 3' CDS LOCATION and 3'_PARTIAL
            # 3' CDS LOCATION # End of the coding region relative to the submitted sequence. For a full length CDS this is the position of the last base of the stop codon.
        threeprime_cds = str(the_gene.location.end.position)
        # PARTIAL AT 3'? (yes/no) # For an incomplete CDS with the stop codon downstream of the submitted sequence.
        if type(the_gene.location.end) == Bio.SeqFeature.ExactPosition:
            threeprime_partial = 'no'
        if type(the_gene.location.end) == Bio.SeqFeature.AfterPosition:
            threeprime_partial = 'yes'
        # READING FRAME  # Mandatory if your CDS is 5' partial as it defines the reading frame. Location of the first base of the first fully-encoded amino acid., Example: 1,2 or 3
        read_frame = "12345"

        source_qualifiers = [f.qualifiers for f in seq_record.features if f.type=='source'][0]
        # ISOLATE
        try:
            isolate = source_qualifiers['isolate'][0]
        except:
            isolate = ''
        # SPEC_VOUCH
        try:
            spec_vouch = source_qualifiers['specimen_voucher'][0]
        except:
            spec_vouch = ''
        # LOCALITY
        try:
            country = source_qualifiers['country'][0]
        except:
            country = ''
        # ECOTYPE
        try:
            ecotype = source_qualifiers['ecotype'][0]
        except:
            ecotype = ''

        # SEQUENCE
        sequence = str(seq_record.seq)

        out_list = [entrynumber,
                    organism_name,
                    env_sam,
                    gene_symbol,
                    product_name,
                    transl_table,
                    fiveprime_cds,
                    threeprime_cds,
                    fiveprime_partial,
                    threeprime_partial,
                    read_frame,
                    isolate,
                    spec_vouch,
                    country,
                    ecotype,
                    sequence
                    ]
        out_string = '\t'.join(out_list) + '\n'
        outp_handle.write(out_string)

    def trnK_matK(self, seq_record, counter, outp_handle):
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
        entrynumber = str(counter+1)  # enumerate counter starts counting at 0
        # ORGANISM_NAME
        organism_name = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]

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
            raise ME.MyException('%s annonex2embl ERROR: Qualifiers for gene `%s` not found.' % ('\n', 'matK'))

        # 5'_CDS and 5'_PARTIAL
            # 5'_CDS: Start of the matK coding region relative to the submitted sequence. For a full length CDS this is the position of the first base of the start codon.
            # NOTE: One nucleotide position has to be added to the start position to make it correct.
        fiveprime_cds = str(matK_gene.location.start.position+1)
        # 5'_PARTIAL: cds partial at 5'? (yes/no) For an incomplete CDS with the start codon upstream of the submitted sequence.
        if type(matK_gene.location.start) == Bio.SeqFeature.ExactPosition:
            fiveprime_partial = 'no'
        if type(matK_gene.location.start) == Bio.SeqFeature.BeforePosition:
            fiveprime_partial = 'yes'
        # 3'_CDS and 3'_PARTIAL
            # 3'_CDS: End of the matK coding region relative to the submitted sequence. For a full length CDS this is the position of the last base of the stop codon.
        threeprime_cds = str(matK_gene.location.end.position)
        # 3'_PARTIAL: cds partial at 3'? (yes/no) For an incomplete CDS with the stop codon downstream of the submitted sequence.
        if type(matK_gene.location.end) == Bio.SeqFeature.ExactPosition:
            threeprime_partial = 'no'
        if type(matK_gene.location.end) == Bio.SeqFeature.AfterPosition:
            threeprime_partial = 'yes'

        source_qualifiers = [f.qualifiers for f in seq_record.features if f.type=='source'][0]
        # ISOLATE
        try:
            isolate = source_qualifiers['isolate'][0]
        except:
            isolate = ''
        # SPEC_VOUCH
        try:
            spec_vouch = source_qualifiers['specimen_voucher'][0]
        except:
            spec_vouch = ''
        # LOCALITY
        try:
            country = source_qualifiers['country'][0]
        except:
            country = ''
        # ECOTYPE
        try:
            ecotype = source_qualifiers['ecotype'][0]
        except:
            ecotype = ''

        # SEQUENCE
        sequence = str(seq_record.seq)

        out_list = [entrynumber,
                    organism_name,
                    fiveprime_cds,
                    threeprime_cds,
                    fiveprime_partial,
                    threeprime_partial,
                    trnK_intron_present,
                    isolate,
                    spec_vouch,
                    country,
                    ecotype,
                    sequence
                    ]
        out_string = '\t'.join(out_list) + '\n'
        outp_handle.write(out_string)

    def rRNA(self, seq_record, counter, charset_syms, outp_handle):
        ''' This function writes a TSV spreadsheet for submission via
            the WEBIN checklist submission system.
        Args:
            seq_record (obj)
            counter (int)
            charset_syms (list)
            outp_handle (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''

        # ENTRYNUMBER
        entrynumber = str(counter+1)  # enumerate counter starts counting at 0
        # ORGANISM_NAME
        organism_name = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]
        # SEDIMENT
        sediment = '_'.join(charset_syms)

        source_qualifiers = [f.qualifiers for f in seq_record.features if f.type=='source'][0]
        # ISOLATE
        try:
            isolate = source_qualifiers['isolate'][0]
        except:
            isolate = ''
        # ISOLATION_SOURCE
        try:
            isol_source = source_qualifiers['isolation_source']
        except:
            isol_source = ''
        # COUNTRY
        try:
            country = source_qualifiers['country'][0]
        except:
            country = ''
        # ECOTYPE
        try:
            lat_lon = source_qualifiers['lat_lon'][0]
        except:
            lat_lon = ''
        # COLLECTION_DATE
        try:
            collection_date = source_qualifiers['collection_date']
        except:
            collection_date = ''

        # SEQUENCE
        sequence = str(seq_record.seq)

        out_list = [entrynumber,
                    organism_name,
                    sediment,
                    isolate,
                    isol_source,
                    country,
                    lat_lon,
                    collection_date,
                    sequence
                    ]
        out_string = '\t'.join(out_list) + '\n'
        outp_handle.write(out_string)

    def ITS(self, seq_record, counter, outp_handle):
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
        entrynumber = str(counter+1)  # enumerate counter starts counting at 0

        # ORGANISM_NAME
        organism_name = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]

        source_qualifiers = [f.qualifiers for f in seq_record.features if f.type=='source'][0]
        # ISOLATE
        try:
            isolate = source_qualifiers['isolate'][0]
        except:
            isolate = ''

        # ENV_SAMPLE
        env_sam = 'no'

        # COUNTRY
        try:
            country = source_qualifiers['country'][0]
        except:
            country = ''

        # SPEC_VOUCH
        try:
            spec_vouch = source_qualifiers['specimen_voucher'][0]
        except:
            spec_vouch = ''

        all_seqrec_features = [f.qualifiers['gene'] for f in seq_record.features]
        # 18S
        if '18S' in all_seqrec_features:
            RNA_18S = 'partial'
        else:
            RNA_18S = 'no'
        # 26S
        if '26S' in all_seqrec_features:
            RNA_26S = 'partial'
        else:
            RNA_26S = 'no'
        # ITS1
        if '18S' in all_seqrec_features:
            ITS1_feat = 'complete'
        else:
            ITS1_feat = 'partial'
        # ITS2
        if '26S' in all_seqrec_features:
            ITS2_feat = 'complete'
        else:
            ITS2_feat = 'partial'
        # 58S
        if 'ITS1' in all_seqrec_features and 'ITS2' in all_seqrec_features:
            RNA_58S = 'complete'
        else:
            RNA_58S = 'partial'

        # SEQUENCE
        sequence = str(seq_record.seq)

        out_list = [entrynumber,
                    organism_name,
                    isolate,
                    env_sam,
                    country,
                    spec_vouch,
                    RNA_18S,
                    ITS1_feat,
                    RNA_58S,
                    ITS2_feat,
                    RNA_26S,
                    sequence
                    ]
        out_string = '\t'.join(out_list) + '\n'
        outp_handle.write(out_string)

    def IGS(self, seq_record, counter, charset_syms, outp_handle):
        ''' This function writes a TSV spreadsheet for submission via
            the WEBIN checklist submission system.
        Args:
            seq_record (obj)
            counter (int)
            charset_syms (list)
            outp_handle (obj)
        Returns:
            currently nothing; writes string to file
        Raises:
            -
        '''

        # ENTRYNUMBER
        entrynumber = str(counter+1)  # enumerate counter starts counting at 0
        # ORGANISM_NAME
        organism_name = [f.qualifiers['organism'] for f in seq_record.features if f.type=='source'][0][0]

        # ENV_SAMPLE
        env_sam = 'no'
        # GENE1 and G1PRESENT
        try:
            gene1 = charset_syms[0]
            g1present = 'yes'
        except:
            gene1 = 'placeholder'
            g1present = 'no'
        # GENE2 and G2PRESENT
        try:
            gene2 = charset_syms[1]
            g2present = 'yes'
        except:
            gene2 = 'placeholder'
            g2present = 'no'

        source_qualifiers = [f.qualifiers for f in seq_record.features if f.type=='source'][0]
        # ISOLATE
        try:
            isolate = source_qualifiers['isolate'][0]
        except:
            isolate = ''
        # SPEC_VOUCH
        try:
            spec_vouch = source_qualifiers['specimen_voucher'][0]
        except:
            spec_vouch = ''
        # COUNTRY
        try:
            country = source_qualifiers['country'][0]
        except:
            country = ''

        # SEQUENCE
        sequence = str(seq_record.seq)

        out_list = [entrynumber,
                    organism_name,
                    env_sam,
                    gene1,
                    g1present,
                    gene2,
                    g2present,
                    isolate,
                    spec_vouch,
                    country,
                    sequence
                    ]

        out_string = '\t'.join(out_list) + '\n'
        outp_handle.write(out_string)
