# -*- coding: utf-8 -*-
"""Shared global constants and functions"""
from math import floor, log10

PROGRAM_NAME = 'biofletch'
#
# Code tables for amino-acid and nucleic-acid encodings.
#
# Amino acids are encoded in frequency order for better
# bit-string compression.  Frequencies do vary somewhat
# among species, but the design judgement was made that
# variations aren't large enough to warrant the complication
# of making a custom table just for mildly-better compression.
# We use the overall frequency found in UNIPROT taken from
# https://www.ebi.ac.uk/uniprot/TrEMBLstats in May 2019.
#
# Nucleic-acid encoding order was arbitrarily chosen to be
# in alphabetical order.  While this encoding is not identical
# with others in common use, it does have the salutary properties
# of being memorable and that biological complements are also bit
# complements, enabling rapid bit-parallel calculation of a
# frequent downstream operation.
#
CODE_TABLE = {'AA':{'bits':5,
                    'order': ['L', 'A', 'G', 'V',
                              'S', 'E', 'R', 'I',
                              'T', 'D', 'K', 'P',
                              'F', 'N', 'Q', 'Y',
                              'M', 'H', 'W', 'C',
                              'X']
                    },
              'NA':{'bits':2,
                    'order':['A', 'C', 'G', 'T']
                    }
              }
_SUFFIXES = ['', 'K', 'M', 'G', 'T', 'P']
#
# Shared functions
#
def engr_notation(qty, powers_of_2=True, digits=None):

    if powers_of_2:
        divisor = 1024.
    else:
        divisor = 1000.
    exponent = 0
    while qty >= divisor and exponent < len(_SUFFIXES) - 1:
        qty /= divisor
        exponent += 1
    if digits:
        fmt_string = '%.'+ '%d'%(digits) + 'f'
    else:
        fmt_string = '%f'
    if qty and digits:
        round_factor = 10**(-int(floor(log10(abs(qty)))) + digits-1)
        qty = round(qty*round_factor)/round_factor
    f = (fmt_string %qty).rstrip('0').rstrip('.')
    return '%s %s' % (f, _SUFFIXES[exponent])