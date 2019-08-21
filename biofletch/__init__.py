# -*- coding: utf-8 -*-
"Convert sequence file to Apache Arrow object"
#
# stdlib imports
#
import logging
import zlib
from collections import defaultdict, OrderedDict
from operator import itemgetter
from pathlib import Path
# third-party imports
from Bio import SeqIO
from Bio.Data import IUPACData
# package imports
from .arrow import BioFletchTableMaker, BioFletchReader
from .common import *
#
# global variables
#
LOGINT_FMT = '%s\t%s\t%s\t%d'
LOGFLOAT_FMT = '%s\t%s\t%s\t%f'
LOGSTR_FMT = '%s\t%s\t%s\t%s'
DUP_NAME = 'Dups'
SUBSTRING_NAME = 'Substr'
AMBIG_NAME = 'Ambig%'
SHORT_NAME = 'Short'
FILESMALL_NAME = 'Small'
LENGTH_NAME = 'Length'
DEFAULT_SEQ_FILE_TYPE = 'fasta'
#
# Classes
#
class ShortCounter(object):
    """count sequences too short
    """

    def __init__(self, length=0, log=False, name=None):
        """
        :param length: minimum length
        :param log: logger object
        :param name: file name
        """
        self.count = 0
        self.log = log
        self.minlen = length
        self.name = name

    def test(self, s, id):
        """test and count if too short
        :param s: sequence
        :param id: ID of sequence for log
        :return: True if s is too short
        """
        if len(s) < self.minlen:
            self.count += 1
            if self.log:
                self.log.debug(LOGINT_FMT,
                               self.name,
                               id,
                               SHORT_NAME,
                               len(s))
            return True
        else:
            return False


class AmbigCounter(object):
    """count sequences with too high a fraction ambiguous
    """
    def __init__(self, frac=0.0, log=False, name=None):
        """
        :param frac: maximum fraction ambiguous
        :param log: logger object
        :param name: file name for log
        """
        self.count = 0
        self.log = log
        self.frac = frac
        self.name = name

    def test(self, s, id):
        """test for fraction ambiguous and count
        :param s: sequence
        :param id: ID of sequence for log
        :return: True if s is above fraction ambiguous
        """
        ambig = sum([i == 'X' for i in s])
        if ambig > self.frac:
            self.count += 1
            if self.log:
                self.log.debug(LOGFLOAT_FMT,
                               self.name,
                               id,
                               AMBIG_NAME,
                               ambig * 100. / len(s))
            return True
        else:
            return False


class DuplicateCounter(object):
    """count sequences too short
    """

    class DuplicateIDDict(defaultdict):
        """A dictionary of lists with a get() that returns the surviving ID
        """

        def __init__(self, concat=False):
            self.concat = concat
            super().__init__(list)

        def get(self, k):
            if self.concat:
                return '|'.join([k] + self[k])
            else:
                return k

    def __init__(self, log=False, name=None, concat_names=False):
        self.exact_count = 0
        self.substring_count = 0
        self.log = log
        self.name = name
        self.hash_dict = {}
        self.dupDict = self.DuplicateIDDict(concat=concat_names)

    def exact(self, s, id):
        "test and count if exact duplicate"
        seq_hash = zlib.adler32(bytearray(str(s), 'utf-8'))
        if seq_hash not in self.hash_dict:
            self.hash_dict[seq_hash] = id
            return False
        else:
            self.exact_count += 1
            first_ID = self.hash_dict[seq_hash]
            self.dupDict[first_ID].append(id)
            if self.log:
                self.log.debug(LOGSTR_FMT,
                               self.name,
                               id,
                               DUP_NAME,
                               first_ID)
            return True

    def index_by_length(self, items):
        length_idx = [(i, len(item))
                      for i, item in enumerate(items)]
        length_idx.sort(key=itemgetter(1))
        return [idx for idx, length in length_idx]

    def substring(self, seqs, remove=True):
        "find and optionally remove exact substring matches in set of seqs"
        removal_list = []
        ascending = self.index_by_length([r.seq for r in seqs])
        for item_num, idx in enumerate(ascending):
            test_seq = seqs[idx].seq
            test_id = seqs[idx].id
            # traverse from biggest to smallest to find the biggest match
            for record in [seqs[i] for i in
                           reversed(ascending[item_num+1:])]:
                if str(test_seq) in str(record.seq):
                    self.substring_count += 1
                    self.dupDict[record.id].append(test_id)
                    removal_list.append(idx)
                    if self.log:
                        self.log.debug(LOGSTR_FMT,
                                       self.name,
                                       test_id,
                                       SUBSTRING_NAME,
                                       record.id)
                    break
        # Optionally remove exact substring matches
        if remove and len(removal_list) > 0:
            removal_list.sort()
            for item_num, idx in enumerate(removal_list):
                seqs.pop(idx-item_num)
        return seqs


class Sanitizer(object):
    """clean up and count potential problems with sequence
       potential problems are:
          dashes:    (optional, removed if remove_dashes=True)
          alphabet:  if not in IUPAC set, changed to 'X'
    """

    def __init__(self, protein=True, remove_dashes=False):
        self.remove_dashes = remove_dashes
        self.seqs_sanitized = 0
        self.chars_in = 0
        self.chars_removed = 0
        self.chars_fixed = 0
        self.endchars_removed = 0
        if protein:
            self.alphabet = IUPACData.protein_letters + IUPACData.protein_letters.lower()\
                    + 'X' +'x' + '-'
        else:
            self.alphabet = IUPACData.unambiguous_dna_letters + 'N' + '-'

    def char_remover(self, s, character):
        """remove positions with a given character
        :param s: mutable sequence
        :return: sequence with characters removed
        """
        removals = [i for i, j in enumerate(s) if j == character]
        self.chars_removed += len(removals)
        [s.pop(pos - k) for k, pos in enumerate(removals)]
        return s

    def fix_alphabet(self, s):
        """replace everything out of alphabet with 'X'
        :param s: mutable sequence
        :return: fixed sequence
        """
        fix_positions = [pos for pos, char in enumerate(s)
                         if char not in self.alphabet]
        self.chars_fixed = len(fix_positions)
        [s.__setitem__(pos, 'X') for pos in fix_positions]
        return s

    def remove_char_on_ends(self, s, character):
        """remove leading/trailing ambiguous residues
        :param s: mutable sequence
        :return: sequence with characterss removed from ends
        """
        in_len = len(s)
        while s[-1] == character:
            s.pop()
        while s[0] == character:
            s.pop(0)
        self.endchars_removed += in_len - len(s)
        return s

    def sanitize(self, s):
        """sanitize alphabet use while checking lengths
        :param s: mutable sequence
        :return: sanitized sequence
        """
        self.seqs_sanitized += 1
        self.chars_in += len(s)
        if len(s) and self.remove_dashes:
            s = self.char_remover(s, '-')
        if len(s):
            s = self.fix_alphabet(s)
        if len(s):
            s = self.remove_char_on_ends(s, 'X')
            s = self.remove_char_on_ends(s, 'x')
        return s


def read_seq(in_file,
             logger=None,
             seq_file_type=DEFAULT_SEQ_FILE_TYPE,
             protein=True,
             out_path=None,
             min_len=0,
             storage='arrow',
             max_ambiguous=100.0,
             remove_dashes=False,
             remove_dups=False,
             remove_substrings=False,
             compression=None,
             encoding=None,
             lengths=False):
    """convert sequence file to arrow object (and optionally parquet file)
    :param in_file: pathlib path to input file
    :param logger: if True, logger object
    :param seq_file_type: a valid SeqIO type, [default FASTA]
    :param protein: if True, protein sequence, else DNA
    :param write:  if True, write parquet file of arrow object
    :param min_len: minimum length of sequences after trims
    :param max_ambiguous: maximum fraction of sequence that can be 'X'
    :param remove_dashes: if True, remove '-' characters
    :param remove_duplicates: if True, remove exact matches
    :param remove_substrings: if True, remove substring matches
    :param lengths: if True, write sequence lengths to log file
    :return:
    """
    if logger:
        logger.debug('processing %s', in_file)
    out_sequences = []
    sanitizer = Sanitizer(protein=protein,
                          remove_dashes=remove_dashes)
    #oo_short = ShortCounter(length=min_len,
    #                        log=logger,
    #                        name=in_file.name)
    #oo_ambig = AmbigCounter(frac=max_ambiguous,
    #                        log=logger,
    #                        name=in_file.name)
    #uplicated = DuplicateCounter(log=logger,
    #                             name=in_file.name)
    nseqs = 0
    with in_file.open('rU') as handle:
        for record in SeqIO.parse(handle, seq_file_type):
            #logger.debug('in[%d]=%s', nseqs, record.seq)
            #seq = record.seq.tomutable()
            #seq = sanitizer.sanitize(seq)
            #logger.debug('out[%d]=%s', nseqs, seq)
            #if not len(seq):
            #    # zero-length string after sanitizing
            #    continue
            #if too_short.test(seq, record.id):
            #   # delete sequences too short
            #    continue
            #if too_ambig.test(seq, record.id):
            #    # delete sequences too ambig
            #    continue
            #if duplicated.exact(seq, record.id) and remove_dups:
            #    # delete exact-duplicate sequences
            #    continue
            #record.seq = seq.toseq()
            out_sequences.append(record)
            nseqs += 1
    # Search for exact substring matches in the set
    #if remove_dups and remove_substrings:
    #    out_sequences = duplicated.substring(out_sequences)
    # Indicate (in ID) records deduplicated
    #if (remove_dups or remove_substrings) \
    #       and len(duplicated.dupDict) > 0:
    #    for record in out_sequences:
    #        if record.id in duplicated.dupDict:
    #            record.id = duplicated.dupDict.get(record.id)
    #            record.description = ''
    #if lengths and logger:
    #   for record in out_sequences:
    #       logger.debug(LOGINT_FMT,
    #                    in_file.name,
    #                    record.id,
    #                    LENGTH_NAME,
    #                    len(record.seq))


#    file_dict = {'name': file.name,
#                     'in_int': sanitizer.chars_in,
#                     'removed_int': sanitizer.chars_removed,
#                     'fixed_int': sanitizer.chars_fixed
#                 'seqs': { # sequence stats
#                     'total': len(out_sequences),
#                     'in': sanitizer.seqs_sanitized,
#                     'ambig': too_ambig.count,
#                     'short': too_short.count,
#                     'dedup': { #deduplicated stats
#                    'dup': duplicated.exact_count,
#                    'substr': duplicated.substring_count,
#                 }
#                 }
    logger.debug('Done reading %s sequences', len(out_sequences))
    seqs = [record.seq for record in out_sequences]
    extra_metadata = OrderedDict([('something',1)])
    tablemaker = BioFletchTableMaker(logger,
                 'test',
                 seqs,
                 protein=False,
                 ids=None,
                 alignment_columns=0,
                 dups=None,
                 storage=storage,
                 substr_ranges=None,
                 compression=compression,
                 encoding=encoding,
                 extra_metadata=extra_metadata)
    table = tablemaker.create_table()
    if out_path:
        tablemaker.write_table(out_path)
    return table


def test():
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(PROGRAM_NAME)
    #testpath = Path('/scratch/joelb/biofletch/nr-smallids-first1M.faa')
    filestem = '/scratch/joelb/biofletch/human-nonambig'
    testpath = Path(filestem+'.fna')
    #compression = 'snappy'
    #compression = 'brotli'
    compression = 'LZ4'
    #compression = 'ZSTD'
    #compression = 'gzip'
    #compression = 'none'
    #encoding = 'bitstring'
    encoding = 'packed'
    #encoding = 'string'
    #encoding = 'bytearray'
    #encoding = 'charlist'
    #encoding = 'int8list'
    outpath = Path(filestem + '-' + encoding
                   + '-' + compression + '.biofl')
    read_seq(testpath,
             logger,
             protein=False,
             storage='bcolz',
             out_path=outpath,
             compression=compression,
             encoding=encoding)
    print('written')
    reader = BioFletchReader(outpath, logger=logger)
    table = reader.read()
    #print(table[0][-1])
    print('done')