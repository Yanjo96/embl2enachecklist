#!/usr/bin/env python2.7
'''
annonex2embl wrapper
'''

#####################
# IMPORT OPERATIONS #
#####################

# Add specific directory to sys.path in order to import its modules
# NOTE: THIS RELATIVE IMPORTING IS AMATEURISH.
# NOTE: COULD THE FOLLOWING IMPORT BE REPLACED WITH 'import embl2enachecklists'?
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'embl2enachecklists'))

import Embl2enachecklistsMain as EMBL2ENAclMain
import globalVariables as GlobVars
import MyExceptions as ME

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2018 Michael Gruenstaeudl'
__info__ = 'embl2enachecklists'
__version__ = '2018.03.26.2000'

#############
# DEBUGGING #
#############

import pdb
#pdb.set_trace()

####################
# GLOBAL VARIABLES #
####################

########
# TODO #
########

############
# ARGPARSE #
############

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))

    # Required
    parser.add_argument('-e',
                        '--embl',
                        help='absolute path to infile; infile in EMBL format; Example: /path_to_input/test.embl',
                        default='/home/username/Desktop/test.embl',
                        required=True)

    parser.add_argument('-o',
                        '--outfile',
                        help='absolute path to outfile; outfile in ENA checklist format (tsv-format); Example: /path_to_output/test.tsv',
                        default='/home/username/Desktop/test.tsv',
                        required=True)

    parser.add_argument('-c',
                        '--cltype',
                        help='Any of the currently implemented checklist types (i.e., `ITS`, `rRNA`, `trnK_matK`, `IGS`, `genomic_CDS`)',
                        default=None,
                        required=True)

    parser.add_argument('--version',
                        help='Print version information and exit',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()

########
# MAIN #
########

    try:
        if args.embl.split('.')[-1] == 'embl':
            EMBL2ENAclMain.embl2enachecklists(args.embl, "output/" + args.outfile + ".enachecklist" ,args.cltype )
            if len(GlobVars.warnings) != 0:
                for warning in GlobVars.warnings:
                    print warning
            print 'Done! :D You can find your output here: ' + "output/" + args.outfile + ".enachecklist"
        else:
            raise ME.WrongInputFile('Error: ' + args.embl + ' is not in embl format')
    except Exception as error:
        try:
            print error.value
        except:
            print error
