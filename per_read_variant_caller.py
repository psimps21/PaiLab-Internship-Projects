import pysam
import json
import csv
import gzip
import os
import datetime
from Bio.SeqIO.FastaIO import SimpleFastaParser
from argparse import ArgumentParser

#
class VCTotalStats:
    """
    A class used to track data from the VariantCaller class

    Attributes:
    _json_list (list): .
    _read_count (int): number of reads.
    _total_read_len (int): running sum of read lengths.
    _total_aligned_bases (int): running sum of aligned bases.
    _total_inserts (int): running sum of inserted bases.
    _total_dels (int): running sum of deleted bases.
    _total_introns (int): running sum of introns.
    _total_sclips (int): running sum of .
    _total_snps (dict): running sum of SNPs.
    """

    def __init__(self):
        """
        The constructor of VCTotal Stats class.

        Parameters:
        _json_list (list): .
        _read_count (int): number of reads.
        _total_read_len (int): running sum of read lengths.
        _total_aligned_bases (int): running sum of aligned bases.
        _total_inserts (int): running sum of inserted bases.
        _total_dels (int): running sum of deleted bases.
        _total_introns (int): running sum of introns.
        _total_sclips (int): running sum of .
        _total_snps (dict): running sum of SNPs.
        """

        self._json_list = []
        self._read_count = 0
        self._total_read_len = 0
        self._total_aligned_bases = 0
        self._total_inserts = 0
        self._total_dels = 0
        self._total_introns = 0
        self._total_sclips = 0
        self._total_snps = {
            'A/C': 0,
            'A/G': 0,
            'A/T': 0,
            'A/N': 0,
            'C/A': 0,
            'C/G': 0,
            'C/T': 0,
            'C/N': 0,
            'G/A': 0,
            'G/C': 0,
            'G/T': 0,
            'G/N': 0,
            'T/A': 0,
            'T/C': 0,
            'T/G': 0,
            'T/N': 0

        }

    @property
    def json_list(self):
        return self._json_list

    @property
    def total_snps(self):
        return self._total_snps

    @property
    def total_inserts(self):
        return self._total_inserts

    @property
    def total_aligned_beses(self):
        return self._total_aligned_bases

    def gen_snp_attribute_dicts(self, vc):
        vc_attr = vars(vc)
        del vc_attr['_read']
        del vc_attr['_ref_seq']
        del vc_attr['_read_seq']
        del vc_attr['_insert_dict']
        del vc_attr['_del_dict']
        del vc_attr['_sclip_dict']
        del vc_attr['_intron_dict']
        del vc_attr['_chr']
        self._json_list.append(vc_attr)

    def gen_full_vc_attribute_dicts(self, vc):
        vc_attr = vars(vc)
        del vc_attr['_read']
        del vc_attr['_ref_seq']
        del vc_attr['_read_seq']
        self._json_list.append(vc_attr)

    def increment_total_snps(self, vc):
        for k, v in vc.snp_dict.items():
            self._total_snps[k] += v

    def gen_percent_snp_type(self):
        new_dict = {}
        num_snps = sum(self._total_snps.values())
        for k, v in self._total_snps.items():
            perc = int(round(v/num_snps, 2) * 100)
            new_dict[k] = (v, perc)

        return new_dict

    def gen_perc_snp_in_read(self):
        return int(round(sum(self._total_snps.values()) / self._total_read_len, 20) * 100)

    def add_snp_vc(self, vc):
        self.increment_total_snps(vc)
        self._total_aligned_bases += vc.read.query_alignment_length
        self.gen_snp_attribute_dicts(vc)
        self._total_read_len += vc.read_len
        self._read_count += 1

    def add_full_vc(self, vc):
        self.increment_total_snps(vc)
        self._total_aligned_bases += vc.read.query_alignment_length
        self._total_inserts += len(vc.insert_dict)
        self._total_dels += len(vc.del_dict)
        self._total_sclips = len(vc.sclip_dict)
        self._total_introns += len(vc.intron_dict)
        self.gen_full_vc_attribute_dicts(vc)
        self._total_read_len += vc.read_len
        self._read_count += 1


class ReferenceGen:
    """
    A class to hold the reference file

    Attributes:
        _ref_file: path to reference genome file
        _chr_dict: dictionary holding the sequences for each chromosome

    """
    def __init__(self, filepath):
        """
        The constructor of the ReferenceGen class

        Parameters:
            _ref_file: path to reference genome file
            _chr_dict: dictionary holding the sequences for each chromosome
        """

        self._ref_file = filepath
        self._chr_dict = self.make_chr_strings()

    @property
    def chr_dict(self):
        return self._chr_dict

    def make_chr_strings(self):
        """
        Generates dictionary with the chromosome header as key and the sequence as a value

        :return: A dictionary of chromosome sequences
        """
        count = 0
        curr_line = None
        chrdict = {}
        full_ref_seq = ''

        if self._ref_file.endswith('.gz'):
            open_ref = gzip.open(self._ref_file, 'rt')
        else:
            open_ref = open(self._ref_file, 'r')

        with open_ref as file:
            for line in file:
                if line[0] == ">":
                    if curr_line:
                        chrdict[chrm] = full_ref_seq
                        count += 1
                        full_ref_seq = ""
                    chrm = line.replace('\n', '').split(' ')[0].replace('>', '')
                    curr_line = line
                    continue
                full_ref_seq += line.strip("\n")

            if len(full_ref_seq) != 0:
                chrdict[chrm] = full_ref_seq

        return chrdict

    def get_chr_seq(self, ref_name):
        """
        Retrieves the sequence of chromosome

        :param ref_name: The header of the chomosome sequence
        :return: The chomosome sequence
        """
        return self._chr_dict[ref_name]



class VariantCaller:
    """
    A class that perfomrs variant calling

    Attributes:
        _scvs (csv.writer): a csv.riter object storing data
        _lcsv csv.writer(): a csv.riter object storing data
        _read (pysam.AlignedSegment): the AlignedSegment object for a read
        _read_len (int): length of a read
        _read_name (str): name of a read
        _ref_seq (str): the reference sequence
        _read_seq (str): the read sequence
        _chr (str): the chromosome header
        _read_start (int): the read start position in the reference sequence
        _snp_dict (dict): a dictionary containing  all present types of SNPs
        _insert_dict (dict): a dictionary containing all insertions
        _sgap_dict (dict): a dictionary containing  all gaps
        _sclip_dict (dict): a dictionary containing all short clippings
        _snp_total (int): total number of SNPs
        _percent_snps (int): the percentage of bases that are SNPs
        _type_snp_dict (dict): a dictionary containing the count of all types of SNPs

    """

    def __init__(self, read, reference, short_csv, long_csv):
        """
        The constructor for the Variant Caller class.

        :param read: an Aligned Segment object
        :param reference: the ReferenceGen object
        :param short_csv: csv.writer object for storing data to be written to file
        :param long_csv: csv.writer object for storing data to be written to file
        """

        self._scsv = short_csv
        self._lcsv = long_csv
        self._read = read
        self._read_len = read.query_length
        self._read_name = read.query_name
        self._ref_seq = reference.get_chr_seq(read.reference_name)
        self._read_seq = read.query_sequence
        self._chr = read.reference_name
        self._read_start = read.reference_start
        self._snp_dict = {}
        self._insert_dict = {}
        self._del_dict = {}
        self._sgap_dict = {}
        self._sclip_dict = {}
        self._snp_total = 0
        self._percent_snps = 0
        self._type_snp_dict = {
            'A/C': 0,
            'A/G': 0,
            'A/T': 0,
            'C/A': 0,
            'C/G': 0,
            'C/T': 0,
            'G/A': 0,
            'G/C': 0,
            'G/T': 0,
            'T/A': 0,
            'T/C': 0,
            'T/G': 0
        }

    @property
    def read(self):
        return self._read

    @property
    def read_len(self):
        return self._read_len

    @property
    def snp_dict(self):
        return self._snp_dict

    @property
    def insert_dict(self):
        return self._insert_dict

    @property
    def sgap_dict(self):
        return self._sgap_dict

    @property
    def sclip_dict(self):
        return self._sclip_dict

    @property
    def del_dict(self):
        return self._del_dict

    @property
    def snp_total(self):
        return self._snp_total

    def gen_dict_str(self, a_dict):
        """

        :param a_dict:
        :return:
        """
        new_str = ''
        for k, v in a_dict.items():
            new_str += k + ':' + v + ';'

        return new_str[:-1]

    def csv_generator(self):
         """ Writes the date of a read to csv """

         self._scsv.writerow([self._read_name, self._read_len, self._chr, self._read_start, self._read.reference_end,
                              self._read.is_reverse, len(self._snp_dict), len(self._del_dict), len(self._insert_dict),
                              len(self._sgap_dict)])

         self._lcsv.writerow([self._read_name, self._read_len, self._chr, self._read_start, self._read.reference_end,
                              self._read.is_reverse, self.gen_dict_str(self._snp_dict),
                              self.gen_dict_str(self._del_dict), self.gen_dict_str(self._insert_dict),
                              self.gen_dict_str(self._sgap_dict)])

    def get_snps(self, pairs):
        """
        Update class attributes with snp data

        :param pairs: a list of tuples (read_pos, ref_pos)
        :return: void
        """
        for read_pos, ref_pos in pairs:
            if ref_pos >= len(self._ref_seq):
                break
            if read_pos and ref_pos:
                if self._ref_seq[ref_pos] == 'N':
                    continue
                if self._ref_seq[ref_pos] != self._read_seq[read_pos]:
                    self._snp_dict[str(ref_pos) + ':' + str(read_pos)] = str(self._ref_seq[ref_pos]) + ':' + \
                                                                              str(self._read_seq[read_pos])
                    self._type_snp_dict[self._ref_seq[ref_pos] + '/' + self._read_seq[read_pos]] += 1

        self._snp_total = sum(self._type_snp_dict.values())
        self._percent_snps = int(round(self._snp_total/self.read_len, 2) * 100)

    def analyze_cig_str(self, cig_tups):
        """

        :param cig_tups: A list of tuples (CIGAR operation, length)
        :return: void
        """

        cigar_ops = {
            # Match
            0: self.cig_op_0,
            # Insertion
            1: self.cig_op_1,
            # Deletion
            2: self.cig_op_2,
            # Skip from reference
            3: self.cig_op_3,
            # Soft clipping
            4: self.cig_op_4,
            # Hard clipping
            5: self.cig_op_5,
            # Padding
            6: self.cig_op_6,
            # Sequence math
            7: self.cig_op_7,
            # SNP
            8: self.cig_op_8
        }
        pairs = self._read.get_aligned_pairs()
        pair_pos = 0
        for cig_tup in cig_tups:
            if cig_tup[0] == 5 or cig_tup[0] == 6:
                continue
            op_pairs = pairs[pair_pos:pair_pos+cig_tup[1]]
            cigar_ops[cig_tup[0]](op_pairs)
            pair_pos += cig_tup[1]

    def cig_op_0(self, pairs):
        """
        Calls function to find all SNPs in read

        :param pairs: a tuple of (read position, reference position) For inserts, deletions, skipping either query or
                        reference position may be None
        :return: void
        """
        self.get_snps(pairs)

    def cig_op_1(self, pairs):
        """
        Find all inserts for the read

        :param pairs: a tuple of (read position, reference position) For inserts, deletions, skipping either query or
                        reference position may be None
        :return: void
        """

        self._insert_dict[str(pairs[0][0])] = ''
        for pair in pairs:
            if isinstance(pair[1], int) and pair[1] >= len(self._ref_seq):
                break
            self._insert_dict[str(pairs[0][0])] += self._read_seq[pair[0]]

    def cig_op_2(self, pairs):
        """
        Find all deletions for the read

        :param pairs: a tuple of (read position, reference position) For inserts, deletions, skipping either query or
                        reference position may be None
        :return: void
        """
        self._del_dict[str(pairs[0][1])] = ''
        for pair in pairs:
            if pair[1] >= len(self._ref_seq):
                break
            self._del_dict[str(pairs[0][1])] += self._ref_seq[pair[1]]

    def cig_op_3(self, pairs):
        """
        Find all gaps in the read

        :param pairs: a tuple of (read position, reference position) For inserts, deletions, skipping either query or
                        reference position may be None
        :return: void
        """

        self._sgap_dict[str(pairs[0][1])] = ''
        for pair in pairs:
            if pair[1] >= len(self._ref_seq):
                break
            self._sgap_dict[str(pairs[0][1])] += self._ref_seq[pair[1]]

    def cig_op_4(self, pairs):
        self._sclip_dict[str(pairs[0][0])] = ''
        for pair in pairs:
            if isinstance(pair[1], int) and pair[1] >= len(self._ref_seq):
                break
            self._sclip_dict[str(pairs[0][0])] += self._read_seq[pair[0]]

    def cig_op_5(self, pairs):
        pass

    def cig_op_6(self, pairs):
        pass

    def cig_op_7(self, pairs):
        for pair in pairs:
            if pair[1] >= len(self._ref_seq):
                break

    def cig_op_8(self, pairs):
        self.get_snps(pairs)


def snp_vc_json_generator(vctotal, data_fname):
    """
    Generate a json file containing the dictionaries of SNP data
    
    :param vctotal: a VCTotalStats object
    :param data_fname: name of data file
    :return: void
    """
    fname = data_fname + '_PRVCSNP_Data.js'
    cfname = data_fname + 'Concise_PRVCSNP_Data.js'
    with open(fname, 'w') as jfile:
        json.dump(vctotal.json_list, jfile)
        # json.dump(time, jfile)

    atts = vars(vctotal)
    del atts['_json_list']
    atts['_percent_snps_in_file'] = vctotal.gen_percent_snp_type()
    atts['_percent_snps_in_read'] = vctotal.gen_perc_snp_in_read()
    atts['_percent_snps_in_aligned'] = int(round(sum(vctotal.total_snps.values())/vctotal.total_aligned_beses, 2) * 100)
    with open(cfname, 'w') as cjfile:
        json.dump(atts, cjfile)


def full_vc_json_generator(vctotal, data_fname):
    """
    Generate a file containing the dictionaries from the VCTotalStats attributes
    
    :param vctotal: a VCTotoalStats object
    :param data_fname: name of data file
    :return: void
    """
    fname = data_fname + '_PRVC_Data.js'
    cfname = data_fname + '_Concise_PRVC_Data.js'
    with open(fname, 'w') as long_doc:
        json.dump(vctotal.json_list, long_doc)

    atts = vars(vctotal)
    del atts['_json_list']
    atts['_percent_snps_in_file'] = vctotal.gen_percent_snp_type()
    atts['_percent_snps_in_read'] = vctotal.gen_perc_snp_in_read()
    atts['_percent_snps_in_aligned'] = int(round(sum(vctotal.total_snps.values())/vctotal.total_aligned_beses, 2) * 100)
    with open(cfname, 'w') as short_doc:
        json.dump(atts, short_doc)


def process_sam_snps(reffname, *argv):
    """
    Generate SNP data for all reads
    
    :param reffname: a reference file path
    :param argv: command line argumants
    :return: void
    """
    reference = ReferenceGen(reffname)
    for arg in argv:
        data_total = VCTotalStats()
        open_sam = pysam.AlignmentFile(arg, 'rb')

        for read in open_sam:
            if len(read.get_reference_positions()) == 0:
                continue
            if read.query_sequence is None:
                continue

            read_vc = VariantCaller(read, reference)
            read_vc.get_snps(read.get_aligned_pairs(matches_only=True))
            data_total.add_snp_vc(read_vc)

        data_fname = arg.split('.')[0]
        snp_vc_json_generator(data_total, data_fname)


def get_args():
    """
    Defines command line interface instructions

    :return: namespace of command line arguments
    """
    parser = ArgumentParser(
        description='Parse and store data from cigar strings.')
    parser.add_argument('-o', '--outdir',
                        help='Path to the directory in which the output files are to be created. Defaults to present '
                             'directory',
                        default='.',
                        nargs='?')
    parser.add_argument('-b', '--bam',
                        help='A path to each BAM file to be analyzed.',
                        required=True,
                        nargs='+',
                        )
    parser.add_argument('-r', '--reference',
                        help='Path to a reference file in fasta format.',
                        required=True,
                        nargs=1
                        )

    return parser.parse_args()


def process_sam_cigstring(args):
    """
    Executes the calls to analyze bam file and write data to file

    :param args: commandline arguments
    :return: void
    """

    reference = ReferenceGen(args.reference[0])

    for arg in args.bam:
        open_sam = pysam.AlignmentFile(arg, 'rb')


        data_fname = os.path.basename(arg).split('.')[0]
        fbasename = data_fname + '_PRVC_Data.csv.gz'
        cfbasename = data_fname + '_Concise_PRVC_Data.csv.gz'
        ffbasename = data_fname + '_Failed_Reads.csv.gz'
        fname = os.path.join(args.outdir, fbasename)
        cfname = os.path.join(args.outdir, cfbasename)
        ffname = os.path.join(args.outdir, ffbasename)

        file = gzip.open(fname, 'wt')
        cfile = gzip.open(cfname, 'wt')
        ffile = gzip.open(ffname, 'wt')
        long_csv = csv.writer(file)
        long_csv.writerow(['readname', 'readlength', 'readchromosome', 'readstart', 'readend',' reverse', 'snps', 'deletions',
                           'insertions', 'splicinggap'])
        short_csv = csv.writer(cfile)
        short_csv.writerow(['readname', 'readlength', 'readchromosome', 'readstart', 'readend', 'reverse', 'snps', 'deletions',
                           'insertions', 'splicinggap'])
        failed_csv = csv.writer(ffile)
        failed_csv.writerow(['readname', 'readchromosome'])

        count = 0
        for read in open_sam:
            if count % 1000 == 0:
                file.flush()
                cfile.flush()
                ffile.flush()

            if len(read.get_reference_positions()) == 0 or read.query_sequence is None:
                continue
            elif read.reference_name not in reference.chr_dict.keys():
                failed_csv.writerow([read.query_name, read.reference_name])
                continue

            read_vc = VariantCaller(read, reference, short_csv, long_csv)
            read_vc.analyze_cig_str(read.cigartuples)
            read_vc.csv_generator()
            count += 1

        file.close()
        cfile.close()
        ffile.close()


if __name__ == '__main__':
    print('Started:', datetime.datetime.now())
    clargs = get_args()
    process_sam_cigstring(clargs)
    print('Finished:', datetime.datetime.now())

