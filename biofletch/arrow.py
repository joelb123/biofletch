# -*- coding: utf-8 -*-
"Schema for Arrow object"
import json
import sys
from collections import OrderedDict
#
# third-party imports
#
import arrow    # this package is about time, nothing to do wity pyarrow
import bitarray
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import bcolz
#
# package imports
#
from .version import version as VERSION
from .common import *
#
# global variables
#
COMPRESSION_TYPES = set(['none'])
#
# Functions
#
def get_bit_dict(code_type='AA'):
    '''Get a dictionary of codes to bit-array elements'''
   # print(CODE_TABLE[code_type])
    n_bits = CODE_TABLE[code_type]['bits']
    code_order = CODE_TABLE[code_type]['order']
    bit_dict = OrderedDict([(code, [bool((i >> b)&1) for b in range(n_bits)])
                             for i, code in enumerate(code_order)])
#    if code_type == 'AA':
#        bit_dict['X'] = [None]*n_bits
#    else:
#        bit_dict['N'] = [None]*n_bits
    return bit_dict

def get_int_dict(code_type='AA'):
    code_order = CODE_TABLE[code_type]['order']
    return OrderedDict([(code, i) for i, code in enumerate(code_order)])

def transpose(l):
    return list(map(list, zip(*l)))


class SequenceArrayDict(OrderedDict):
    """Construct an ordered dictionary of sequence-derived arrays

    This dictionary will later be made into a table after
    possible addition of other arrays such as masks,
    quality scores, or alignments.  All arrays must therefore be of
    the same length as the sequence from which it is derived.

    The arrays are generally derived from code tables that are
    specific to either nucleic acids or amino acids. Anything not
    found in these code tables (e.g., 'X', or 'N') will be stored
    as null values.  This choice was made to enable a 33% savings
    in storage space for nucleic-acid strings.  Null values are
    part of the Apache Arrow architecture, so we might as well
    use them."""


    def __init__(self, seqs, logger, seq_type='AA', encoding=None):

        if encoding is None: # default encoding is shortest for type
            if seq_type == 'AA':
                encoding = 'string'
            elif seq_type == 'NA':
                encoding = 'bitstring'
        if logger:
            logger.debug('encoding %d sequences with %s encoder',
                         len(seqs), encoding)
        if encoding == 'string':
            seq_list = [('%s'%seq_type, pa.array([str(seq)
                                                  for seq in seqs]))]
        elif encoding == 'charlist':
            seq_list = [('%s'%seq_type, pa.array([list(str(seq))
                                                  for seq in seqs]))]
        elif encoding == 'bytearray':
            seq_list = [('%s'%seq_type, pa.array([bytearray(str(seq).encode('UTF-8'))
                                                  for seq in seqs],
                                                 type=pa.list_(pa.uint8())))]
        elif encoding == 'int8list':
            int_dict = get_int_dict(code_type=seq_type)
            seq_list = [('%s'%seq_type, pa.array([[int_dict[c] for c in seq]
                                                  for seq in seqs],
                                                 type=pa.list_(pa.uint8())))]
        elif encoding == 'packed':
            int_dict = get_int_dict(code_type=seq_type)
            n_bits = CODE_TABLE[seq_type]['bits']
            bit_arr = []
            for seq in seqs:
                bytestring = bitarray.bitarray(endian='little')
                packed = bitarray.bitarray(endian='little')
                bytestring.frombytes(bytes([int_dict[c] for c in seq.upper()]))
                for i in range(len(seq)):
                    packed[i*n_bits:] = bytestring[i*8:(i*8)+n_bits]
                bit_arr.append(packed.tobytes())
            seq_list = [('%s'%seq_type, pa.array(bit_arr,
                                                 type=pa.list_(pa.uint8())))]
        elif encoding == 'bitstring':
            self.bit_dict = get_bit_dict(code_type=seq_type)
            self.bits = CODE_TABLE[seq_type]['bits']
            bitlists = transpose([transpose([self.bit_dict[c] for c in seq]) for seq in seqs])
            seq_list = [('%s%d'%(seq_type, b), pa.array(arr, type=pa.list_(pa.bool_())))
                          for b, arr in enumerate(bitlists)]
        else:
            if logger:
                logger.error('unknown encoding type "%s"', encoding)
            else:
                print('unknown encoding type "%s"'%encoding)
            sys.exit(1)
        super().__init__(seq_list)


class BioFletchTableMaker(object):
    def __init__(self,
                 logger,
                 set_id,
                 seqs,
                 protein=True,
                 ids=None,
                 alignment_columns=0,
                 dups=None,
                 substr_ranges=None,
                 compression=None,
                 encoding=None,
                 storage=None,
                 extra_metadata=None):
        self.logger = logger
        self.set_id = set_id
        self.rows = len(seqs)
        self.storage = storage
        if protein:
            self.seq_type = 'AA'
        else:
            self.seq_type = 'NA'
        self.columns = SequenceArrayDict(seqs,
                                         self.logger,
                                         seq_type=self.seq_type,
                                         encoding=encoding)
        #
        # If alignment, use a fixed number of columns
        #
        self.alignment_columns = alignment_columns
        if not self.alignment_columns:
            schema_name = PROGRAM_NAME + '.alignment'
        else:
            schema_name = PROGRAM_NAME + '.genome'
        #
        # If no ID's provided, number them starting at zero
        #
        if ids is None:
            ids = ['%d' % i for i in range(self.rows)]
        #self.columns.update([('id', pa.array(ids))])
        #
        # Calculate length-related parameters
        #
        lengths = np.array([len(s) for s in seqs])
        len_max = int(max(lengths))
        len_min = int(min(lengths))
        len_total = int(lengths.sum())
        len_mean = round(lengths.mean(),1)
        len_std = lengths.std()
        if np.isnan(len_std):
            len_std = 0.
        else:
            len_std = round(len_std,2)
        #self.columns.update([('len', pa.array(lengths))])
        del lengths
        #
        # Deal with duplications
        #
        if dups:
            self.columns.update(('duplicate_of', pa.array(dups, type=pa.int32())))
            dedup_type = 'duplicate'
            if substr_ranges:
                self.columns.update(('substr_start', pa.array([r[0] for r in substr_ranges],
                                                              type=pa.int32)))
                self.columns.update(('substr_stop', pa.array([r[1] for r in substr_ranges],
                                                              type=pa.int32)))
                dedup_type = 'substring'
        else:
            dedup_type = 'none'
        #
        # Compression type
        #
        if compression is None:
            self.compression = 'none'
        else:
            self.compression = compression
        #else:
        #    if compression not in COMPRESSION_TYPES:
        #        if logger:
        #            logger.error('Unknown compression type "%s"', compression)
        #        sys.exit(1)

        #
        # Metadata.  Note this has to all be byte strings, so anything non-string
        # should go into the JSON-encoded 'extra' dictionary.
        #
        self.extra_metadata = OrderedDict([('seq_type', self.seq_type),
                                          ('len_min', len_min),
                                          ('len_max', len_max),
                                          ('len_mean', len_mean),
                                          ('len_std', len_std),
                                          ('len_total', len_total)
                                          ])
        if extra_metadata:
            self.extra_metadata.update(extra_metadata)
        self.metadata_dict = OrderedDict([
            ('set_id', set_id),
            ('schema_name', schema_name),
            ('version', VERSION),
            ('sequence_compression', self.compression),
            ('deduplication', dedup_type),
            ('extra', json.dumps(self.extra_metadata))
        ])
        #
        # ID string fields
        #
        #'organism_id', pa.string())
        #'genome_id', pa.string())
        #'fragment_id', pa.string())
        #'alt_id_list', pa.list(pa.string()))
        #
        # Origin fields
        #
        #start_loc_fld = pa.field('start_loc', pa.uint32())
        #positive_dir_fld = pa.field('positive_dir', pa.bool_())
        #plastid_fragment = pa.field('plastid_fragment', pa.bool_())
        #plastid_gene = pa.field('plastid_gene', pa.bool())
        #
        # Length-related fields
        #
        #n_masked_fld = pa.field('n_masked', pa.uint32())
        #n_ambig_fld = pa.field('n_ambig', pa.uint32())


    def add_array(self, name, values, arr_type=None):
        if arr_type is None:
            arr_type = pa.uint32()
        self.columns.update([(name, pa.array(values, type=arr_type))])

    def add_dictionary_array(self, name, values):
        if arr_type is None:
            arr_type = pa.string()
        self.columns.update([(name, pa.array(values))]).dictionary_encode()

    def create_table(self):
        if self.logger:
            self.logger.debug('creating table with %d columns',
                              len(self.columns))
        if self.storage == 'arrow':
            self.table = pa.Table.from_arrays([arr for arr in self.columns.values()],
                                          [key for key in self.columns.keys()],
                                          #metadata=self.metadata_dict
                                          )
        elif self.storage == 'bcolz':
            self.table = bcolz.
        del self.columns
        return self.table

    def write_table(self,filename):
        if self.logger:
            self.logger.info('Writing to file "%s"', filename)
        pq.write_table(self.table, filename, #use_dictionary=False,
                       compression=self.compression,
                       version='2.0'
                    )


class BioFletchReader(object):
    def __init__(self, filename, logger=None):
        self.filename = filename
        self.logger = logger

    def read(self):
        start_time = arrow.utcnow()
        if self.logger:
            self.logger.debug('Reading %s', self.filename)
        table = pq.read_table(self.filename)
        read_time = (arrow.utcnow()-start_time).total_seconds()
        n_bytes = pa.total_allocated_bytes()
        if self.logger:
            self.logger.debug('Reading %s sequences (%sB) took %.1f s (%sB/s).',
                              engr_notation(len(table), powers_of_2=False),
                              engr_notation(n_bytes, digits=2),
                              read_time,
                              engr_notation(n_bytes/read_time, digits=2))
        return table