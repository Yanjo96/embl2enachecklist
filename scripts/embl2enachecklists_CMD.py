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
__version__ = '2018.09.01'

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
    parser.add_argument('-i',
                        '--infile',
                        help='absolute path to infile; infile in EMBL or GenBank format; Example: /path_to_input/test.embl',
                        default='/home/username/Desktop/test.embl',
                        required=True)

    parser.add_argument('-o',
                        '--outfile',
                        help='absolute path to outfile; outfile in ENA checklist format (tsv-format); Example: /path_to_output/test.tsv',
                        default='/home/username/Desktop/test.tsv',
                        required=True)

    parser.add_argument('-c',
                        '--cltype',
                        help='Any of the currently implemented checklist types:' + ", ".join(GlobVars.allowed_checklists),
                        default=None,
                        required=True)

    parser.add_argument('-e',
    		            '--environmental',
    		            help="Is your organism from an environmental/uncultured sample? yes/no",
    		            required=True)

    parser.add_argument('--version',
                        help='Print version information and exit',
                        action='version',
                        version=__info__ + ", version: " + __version__)

    args = parser.parse_args()

########
# MAIN #
########

    try:
        if args.infile.split('.')[-1] == 'embl':
            EMBL2ENAclMain.embl2enachecklists(args.infile, args.outfile, 'embl', args.cltype, args.environmental)
            if len(GlobVars.warnings) != 0:
                for warning in GlobVars.warnings:
                    print warning
            print 'Done! :D You can find your output here: ' + args.outfile
        elif args.infile.split('.')[-1] == 'gb':
            EMBL2ENAclMain.embl2enachecklists(args.infile, args.outfile, 'gb', args.cltype, args.environmental)
            if len(GlobVars.warnings) != 0:
                for warning in GlobVars.warnings:
                    print warning
            print 'Done! :D You can find your output here: ' + args.outfile
        else:
            raise ME.WrongInputFile('Error: ' + args.infile + ' is not in embl format')
    except Exception as error:
        try:
            print error.value
        except:
            print error
