#!/usr/bin/env python3
# the above sources Python from $PATH

##!/usr/local/bin/python3
##!/usr/bin/python3
# the above uses specific Python version; allows script name in top

"""
usage: intronIC [-h] [-g GENOME] [-a ANNOTATION] [-b BED_FILE] -n SPECIES_NAME
                [-q SEQUENCE_FILE] [-f {cds,exon}] [-s] [--no_nc] [-i] [-v]
                [--pwms {PWM files)} [{PWM file(s} ...]]
                [--reference_u12s {reference U12 intron sequences}]
                [--reference_u2s {reference U2 intron sequences}] [--no_plot]
                [--format_info] [-d] [-u] [--no_abbreviate] [-t 0-100]
                [--no_sequence_output] [--five_score_coords start stop]
                [--three_score_coords start stop]
                [--branch_point_coords start stop]
                [-r {five,bp,three} [{five,bp,three} ...]]
                [--abbreviate_filenames] [--recursive]
                [--n_subsample N_SUBSAMPLE] [--cv_processes CV_PROCESSES]
                [-p PROCESSES] [--matrix_score_info] [-C HYPERPARAMETER_C]
                [--min_intron_len MIN_INTRON_LEN] [--pseudocount PSEUDOCOUNT]
                [--exons_as_flanks]
"""

__author__ = 'Graham E. Larue'
__maintainer__ = "Graham E. Larue"
__email__ = 'egrahamlarue@gmail.com'
__license__ = 'GPL v3.0'

# imports
import argparse
import copy
import logging
import math
import os
import re
import random
import sys
import time
import gzip
import numpy as np
import warnings
import types
import pkg_resources

# hacky way to ignore annoying sklearn warnings
# (https://stackoverflow.com/a/33616192/3076552)
def warn(*args, **kwargs):
    pass
warnings.warn = warn

from multiprocessing.pool import Pool
from multiprocessing import get_all_start_methods, set_start_method
from scipy import stats as pystats
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict, deque, namedtuple
from itertools import islice, repeat, chain
from operator import attrgetter
from functools import partial
# from hashlib import sha1
from sklearn.metrics import precision_recall_curve, auc, classification_report
# from sklearn.cluster import SpectralClustering
from sklearn import svm, preprocessing
from sklearn.base import clone
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.metrics import f1_score
from sklearn import linear_model
from biogl import fasta_parse, get_runtime, rev_comp, flex_open, GxfParse

try:
    import matplotlib
    matplotlib.use('Agg') # allow to run without X display server
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from mpl_toolkits.axes_grid1 import make_axes_locatable
except ModuleNotFoundError:
    pass

# nuclear option to ignore sklearn warnings during fitting process
# does not appear to persist as env variable after script execution
# if not sys.warnoptions:
#     warnings.simplefilter("ignore")
#     os.environ["PYTHONWARNINGS"] = "ignore" # Also affect subprocesses

# improve sklearn parallel performance/stability
# this option is limited by certain OS constraints,
# and may cause problems outside Unix...?
# os.environ['JOBLIB_START_METHOD'] = 'forkserver'

# Classes ####################################################################

class GenomeFeature(object):
    """
    Features that all genomic entities should
    have.

    >start< and >stop< are always relative to positive strand, i.e.
    >start< is always less than >stop<.

    """
    count = 1

    # __slots__ prevents objects from adding new attributes, but it
    # significantly reduces the memory footprint of the objects
    # in use. Idea from http://tech.oyster.com/save-ram-with-python-slots/

    __slots__ = [
        'line_number', 'region', 'start', 'stop', 'parent_type',
        'strand', 'name', 'parent', 'seq', 'flank', 'feat_type', 'phase',
        'upstream_flank', 'downstream_flank', 'family_size', 'unique_num'
    ]

    def __init__(
            self, line_number=None, region=None,
            start=None, stop=None, parent_type=None,
            strand=None, name=None, parent=None,
            seq=None, flank=0, feat_type=None,
            upstream_flank=None, downstream_flank=None,
            family_size=0, phase=None
    ):
        self.region = region
        self.start = start
        self.stop = stop
        self.strand = strand
        self.name = name
        self.unique_num = self.__class__.count
        self.feat_type = feat_type
        self.parent = parent
        self.parent_type = parent_type
        self.family_size = family_size
        self.line_number = line_number
        self.phase = phase
        self.seq = seq
        self.flank = flank
        self.upstream_flank = upstream_flank
        self.downstream_flank = downstream_flank

    @property
    def length(self):
        """
        Returns the length of the object, preferentially
        inferring it from the start and stop coordinates,
        then from the length of the full sequence.

        """
        if not (self.start and self.stop):
            if not self.seq:
                return None
            else:
                return len(self.seq)
        return abs(self.start - self.stop) + 1

    def get_coding_length(self, child_type="cds"):
        """
        Returns an integer value of the aggregate
        length of all children of type child_type.

        If child_type is None, returns aggregate
        length of all children

        """
        total = 0
        while True:
            # while children themselves have children, recurse
            try:
                children = [c for c in self.children if c.children]
                for child in children:
                    total += child.get_coding_length(child_type)
            except AttributeError:
                try:
                    children = self.get_children(child_type)
                    total += sum(c.length for c in children)
                except AttributeError:  # if was called on exon or cds objects
                    total += self.length
                break
            break
        return total

    def set_family_size(self):
        """
        Assigns family_size attribute to all children.

        """
        try:
            all_children = self.children
        except AttributeError:
            return
        child_types = set([c.feat_type for c in all_children])
        for ct in child_types:
            siblings = [c for c in all_children if c.feat_type == ct]
            sibling_count = len(siblings)
            for sib in siblings:
                sib.family_size = sibling_count
                sib.set_family_size()

    def compute_name(self, parent_obj=None, set_attribute=True):
        """
        Returns a generated unique ID for the object
        in question. Passing a parent object should
        be preferred, as the ID is otherwise not
        particularly informative (but is unique!)

        """
        parent = self.parent
        parent_type = self.parent_type
        # Make a coord_id relative to parent's contents if possible
        try:
            # 1-based indexing
            unique_num = parent_obj.children.index(self) + 1
        # Otherwise, make relative to all made objects of same type
        except AttributeError:
            # every feature type should have this
            unique_num = "u{}".format(self.unique_num)
        coord_id = ("{}:{}_{}:{}".format(parent_type,
                                    parent,
                                    self.feat_type,
                                    unique_num))
        if set_attribute:
            setattr(self, "name", coord_id)
        return coord_id

    def update(self, other, attrs=None):
        """
        Updates attributes based on another object and
        optional attribute filter list

        """
        attrs_to_use = vars(other)  # use all if not attrs
        if attrs:
            attrs_to_use = {a: v for a, v in attrs_to_use.items() if a in attrs}
        for key, value in attrs_to_use.items():
            if hasattr(self, key):  # don't make new attributes
                setattr(self, key, value)


    def get_seq(self, region_seq=None, start=None, stop=None,
                flank=None, strand_correct=True):
        """
        Retrieves object's sequence from its parent
        sequence, with optional flanking sequence.

        If >strand_correct<, will return strand-
        corrected (e.g. reverse-complemented) sequence.

        """
        if region_seq is None:
            return self.seq

        # Assign class defaults if not values
        if start is None:
            start = self.start
        if stop is None:
            stop = self.stop

        # Correct for 1-based indexing in start and stop
        start -= 1
        # avoid negative indices
        start = max(start, 0)
        # Pull sequence and reverse if necessary
        seq = region_seq[start:stop]
        # TODO reverse the region seq once outside of the 
        # function to speed things up?
        if strand_correct and self.strand == '-':
            seq = rev_comp(seq, use_lower=False)
        if flank:
            try:
                up_flank_length, down_flank_length = flank
            except TypeError:  # same length for both
                up_flank_length = flank
                down_flank_length = flank
            up_flank = self.upstream_seq(region_seq, up_flank_length)
            down_flank = self.downstream_seq(region_seq, down_flank_length)
            seq = list(map(str.upper, [up_flank, seq, down_flank]))

        return seq

    def upstream_seq(self, region_seq, n, strand_correct=True):
        """
        Get sequence of n length from upstream of
        feature start, relative to coding direction.

        """
        if self.strand == "-":
            start = self.stop
            stop = start + n
        else:
            stop = self.start - 1
            start = max(stop - n, 0) # no negative indices
        seq = region_seq[start:stop]
        if strand_correct and self.strand == '-':
            seq = rev_comp(seq, use_lower=False)

        return seq

    def downstream_seq(self, region_seq, n, strand_correct=True):
        """
        Get sequence of n length from downstream of
        feature start, relative to coding direction.

        """
        if self.strand == "-":
            stop = self.start - 1
            start = max(stop - n, 0) # no negative indices
        else:
            start = self.stop
            stop = start + n
        seq = region_seq[start:stop]
        if strand_correct and self.strand == '-':
            seq = rev_comp(seq, use_lower=False)

        return seq


class Parent(GenomeFeature):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.children = []
        self._intronator = intronator

    def get_children(self, child_type=None):
        """
        Returns a list of all children of type >child_type<.

        If >child_type< not specified, returns all children.

        """
        if not child_type:
            selected = self.children
        else:
            selected = [c for c in self.children if c.feat_type == child_type]
        return selected


    def get_introns(self, child_types, flat_annots):
        """
        Returns all introns based on >child_type<,
        including any duplicates across children.

        """
        introns = []
        non_redundant = []
        filtered_children = []
        intron_count = 1
        try:
            children = [child for child in self.children if
                                 child.feat_type in child_types]
        except AttributeError:
            return introns
        if not children:
            try:
                for child in self.children:
                    introns += child.get_introns(child_types, flat_annots)
                return introns
            except AttributeError:
                return introns
        # track coding lengths calculated under different
        # feature types to preference CDS-based lengths at the end
        coding_lengths = {}
        for ct in child_types:  # CDS > exons depends on this list's order?
            tmp_introns = []
            filtered_children = [
                c for c in children if c.feat_type == ct]
            if not filtered_children:
                continue
            coding_lengths[ct] = self.get_coding_length(ct)
            for indx, intron in enumerate(
                self._intronator(filtered_children), start=1):
                intron.defined_by = ct
                intron.index = indx
                tmp_introns.append(intron)
            if not non_redundant:
                non_redundant = tmp_introns
            else:
                existing_coords = sorted(
                    (i.start, i.stop) for i in non_redundant)
                for i in tmp_introns:
                    coords = (i.start, i.stop)
                    if not coord_overlap(
                        coords, existing_coords, presorted=True):
                        non_redundant.append(i)
        if not non_redundant:
            return non_redundant
        
        # prioritize protein-coding transcript lengths over others
        try:
            coding_length = coding_lengths['cds']
        except KeyError:
            coding_length = coding_lengths['exon']
        family_size = len(non_redundant)
        # get the lengths of all intron-forming features to allow
        # calculation of fractional position of introns

        # TODO: figure out how to correctly get sorted list on
        # this call, given that the features being sorted are
        # introns - partially resolved by using CDS info if present
        non_redundant = sorted_in_coding_direction(non_redundant)
        exon_lengths = [flat_annots[i.upstream_exon].length for i in non_redundant]
        last_exon = non_redundant[-1].downstream_exon
        exon_lengths.append(flat_annots[last_exon].length)
        exon_cumsum = np.array(exon_lengths)[:-1].cumsum()
        aggregate_length = sum(exon_lengths)
        frac_positions = ((exon_cumsum / aggregate_length) * 100).round(3)
        
        for index, i in enumerate(
            non_redundant, start=1):
            i.index = index
            i.family_size = family_size
            i.parent_length = coding_length
            i.fractional_position = frac_positions[index - 1]
            if i.defined_by == 'exon':
                i.dynamic_tag.add('[e]')
            introns.append(i)

        return introns


class Gene(Parent):
    def __init__(self, parent=None, **kwargs):
        super().__init__(**kwargs)
        self.__class__.count += 1
        self.feat_type = "gene"
        self.parent = parent
        self.parent_type = None


class Transcript(Parent):
    def __init__(self, parent=None, **kwargs):
        super().__init__(**kwargs)
        self.__class__.count += 1
        self.feat_type = "transcript"
        self.parent_type = "gene"
        if not parent:
            self.parent = self.name
        else:
            self.parent = parent


class Exon(GenomeFeature):
    def __init__(self, feat_type="exon", parent=None,
                 grandparent=None, **kwargs):
        super().__init__(**kwargs)
        self.__class__.count += 1
        self.parent_type = "transcript"
        self.feat_type = feat_type
        self.parent = parent
        self.grandparent = grandparent


class Intron(GenomeFeature):
    __slots__ = [
        '__dict__', 'bp_raw_score', 'bp_region_seq', 'bp_seq', 'bp_start',
        'bp_stop', 'bp_z_score', 'corrected', 'dnts', 'downstream_exon',
        'downstream_flank', 'duplicate', 'dynamic_tag', 'family_size',
        'feat_type', 'five_display_seq', 'five_raw_score',
        'five_score_coords', 'five_seq', 'five_start', 'five_stop',
        'five_z_score', 'flank', 'fractional_position', 'grandparent',
        'index', 'line_number', 'longest_isoform', 'name', 'noncanonical',
        'omitted', 'overlap', 'parent', 'parent_type', 'parent_length', 
        'phase', 'region', 'seq', 'start', 'stop', 'strand', 'type_id',
        'three_display_seq', 'relative_score', 'defined_by',
        'three_score_coords', 'three_seq', 'three_z_score',
        'three_raw_score', 'u12_matrix', 'u2_matrix', 'unique_num',
        'upstream_exon', 'upstream_flank', 'svm_score', 'matrices'
    ]

    def __init__(self, parent=None, grandparent=None, **kwargs):
        # Set certain attrs from parent class
        super().__init__(**kwargs)
        self.__class__.count += 1
        # Set other intron-only attrs
        # Inherit from transcript
        self.feat_type = "intron"
        self.parent_type = "transcript"
        self.parent = parent
        self.parent_length = 0
        self.grandparent = grandparent
        self.five_raw_score = None
        self.five_z_score = None
        self.bp_raw_score = None
        self.bp_z_score = None
        self.index = None  # in coding orientation
        self.defined_by = None
        self.five_seq = None
        self.three_seq = None
        self.bp_seq = None
        self.bp_region_seq = None
        self.five_start = None
        self.five_stop = None
        self.bp_start = None
        self.bp_stop = None
        self.omitted = False
        self.corrected = False
        self.duplicate = False
        self.overlap = None
        self.longest_isoform = None
        self.dynamic_tag = set()
        self.noncanonical = False
        self.upstream_exon = None
        self.downstream_exon = None
        self.upstream_flank = None
        self.downstream_flank = None
        self.fractional_position = '.'
        self.five_score_coords = None
        self.three_score_coords = None
        self.five_display_seq = None
        self.three_display_seq = None
        self.u2_matrix = None
        self.u12_matrix = None
        self.matrices = set()
        self.three_raw_score = None
        self.three_z_score = None
        self.dnts = None
        self.svm_score = None
        self.relative_score = None
        self.type_id = None

    @classmethod
    def from_exon_pair(cls, up_ex, down_ex):
        """
        Takes a pair of Exon objects and builds an intron
        based on their information.

        """
        # Infer intron coordinates
        start = min(up_ex.stop, down_ex.stop) + 1
        stop = max(up_ex.start, down_ex.start) - 1
        # Get applicable attributes from one of the defining coding objects
        strand = up_ex.strand
        parent = up_ex.parent
        grandparent = up_ex.grandparent
        region = up_ex.region
        # derive phase from upstream exon (CDS) phase annotation
        # (if available)
        if up_ex.phase != '.':  # default phase value if not present
            phase = (up_ex.length - up_ex.phase) % 3
        else:
            phase = '.'
        fam = up_ex.family_size - 1  # intron number is exon number - 1

        # average the line numbers of the children that define each intron
        # to enable tie-breaking downstream when deciding which duplicates
        # to exclude (for those whose parents have equal length)
        line_number = sum([x.line_number for x in (up_ex, down_ex)]) / 2

        return cls(start=start, stop=stop, strand=strand, family_size=fam,
                   parent=parent, grandparent=grandparent, region=region,
                   line_number=line_number, phase=phase)


    def get_rel_coords(self, relative_to, relative_range):
        """
        Calculates and retrieves a pair of genomic sequence coordinates
        based on input relative to the five-prime or three-prime
        end of the sequence. Returns a set of adjusted coords.

        e.g. get_rel_coords("five", "five", (0, 12)),
             get_rel_coords("bp", "three", (-45, -5))
        """

        def __upstream(ref, x, strand):
            x = abs(x)
            if strand == "-":
                return ref + x
            else:
                return ref - x

        def __downstream(ref, x, strand):
            x = abs(x)
            if strand == "-":
                return ref - x
            else:
                return ref + x

        start = self.start
        stop = self.stop
        strand = self.strand
        # Begin orientation gymnastics
        if strand == "-":
            start, stop = stop, start
        if relative_to == "five":
            ref_point = start
        else:
            ref_point = stop
        new_coords = []
        for n in relative_range:  # each number can be + or -
            if n < 0:
                new_coords.append(__upstream(ref_point, n, strand))
            else:
                new_coords.append(__downstream(ref_point, n, strand))
        new_coords = tuple(sorted(new_coords))
        # setattr(self, seq_name, new_coords)
        return new_coords

    def get_name(self, spcs, simple=False, special='?'):
        """
        Build a unique name (not including score) from metadata.

        """
        if self.name is not None:
            return self.name
        if self.omitted:
            omit_tag = ';[o:{}]'.format(self.omitted)
        else:
            omit_tag = ''
        if self.dynamic_tag:
            dyn_tag = ';{}'.format(';'.join(sorted(self.dynamic_tag)))
        else:
            dyn_tag = ''
        if simple is True:
            return '{}-i_{}{}{}'.format(spcs, self.unique_num, omit_tag, dyn_tag)
        elements = [
            self.grandparent, self.parent,
            self.index, self.family_size, omit_tag]
        tags = [e if e is not None else special for e in elements]
        tags.append(dyn_tag)  # compatibility with earlier Python 3s
        name = "{}-{}@{}-intron_{}({}){}{}".format(spcs, *tags)
        # setattr(self, "name", name)
        return name

    def get_label(self, spcs, simple=False, special='?'):
        """
        Builds a unique intron label from metadata
        """
        # if self.omitted:
        #     setattr(self, 'relative_score', 0)
        if self.relative_score is not None:
            # clipping prevents rounding from pushing introns over the
            # u12 boundary
            # *clip* float to 3 places (without rounding)
            truncated = math.floor(self.relative_score * 1000) / 1000
            score = '{}%'.format(truncated)
            # *round* float to 4 places
            # rel_score = '{:.4f}%'.format(self.relative_score)
        else:
            score = None
        if not self.name:
            self.name = self.get_name(spcs, simple)
        label = ('{};{}'
                 .format(*[e if e is not None else special for e in
                           [self.name, score]]))
        return label


    def omit_check(
        self, min_length, bp_matrix_length, allow_noncanon=False,
        allow_overlap=False, longest_only=True):
        #TODO make .omitted be a list containing all tags
        # instead of just a single string
        """
        Checks an intron object for omission criteria, and sets
        the >omitted< attribute accordingly.

        """
        omit_tags = {
            'short': 's',
            'ambiguous sequence': 'a',
            'noncanonical': 'n',
            'coordinate overlap': 'v',
            'not in longest isoform': 'i'
        }
        scoring_regions = ['five_seq', 'three_seq']
        omission_reason = None
        if self.length < min_length:
            omission_reason = 'short'
        elif any(valid_chars(getattr(self, region)) is False
                 for region in scoring_regions):
            omission_reason = 'ambiguous sequence'
        # check if there is sufficiently long sequence in the
        # bp region to score at least one bp motif
        elif longest_match(self.bp_region_seq) < bp_matrix_length:
                omission_reason = 'short'
        elif not allow_noncanon and self.noncanonical:
            omission_reason = 'noncanonical'
        elif longest_only and self.longest_isoform is False:
            omission_reason = 'not in longest isoform'
        elif allow_overlap is False and self.overlap:
            omission_reason = 'coordinate overlap'
        if omission_reason:
            setattr(self, 'omitted', omit_tags[omission_reason])
    
    def motif_string(self, exonic=3):
        """
        Returns a schematic of all scored motifs, and the BPS context as strings

        """
        five_boundary = '{}|{}'.format(
            self.upstream_flank[-exonic:], self.five_display_seq)
        three_boundary = '{}|{}'.format(
            self.three_display_seq, self.downstream_flank[:exonic])
        schematic_bits = [five_boundary, self.bp_seq, three_boundary]
        schematic_bits = [b for b in schematic_bits if b is not None]
        schematic_string = '...'.join(schematic_bits)
        
        return schematic_string

    def bps_context(self):
        try:
            context = annotate(
                self.bp_region_seq,
                *self.bp_relative_coords) + self.three_display_seq
        except:
            context = None
        
        return context
    
    def model_score(
        self, 
        models, 
        scoring_regions, 
        THRESHOLD,
        model_weights=None, 
        add_type=True
    ):
        """
        Assigns probability scores to self using supplied
        models; if multiple models, will assign average of all
        scores.

        """
        score_vector = np.array([[getattr(self, r) for r in scoring_regions]])
        probabilities, labels, distances = svm_predict(score_vector, models)
        probabilities = list(chain.from_iterable(probabilities))
        labels = list(chain.from_iterable(labels))
        distances = list(chain.from_iterable(distances))

        info = average_svm_score_info(
            probabilities, labels, distances, model_weights)
        
        self.svm_score = info['u12_avg'] * 100
        self.score_distance = info['avg_distance']
        
        # self.error = avg_u12_probabilities[idx]['u12_sem'] * 100

        # relative as percentage of threshold
        # self.relative_score = self.svm_score - THRESHOLD) / THRESHOLD * 100
        # relative as probability
        self.relative_score = self.svm_score - THRESHOLD
        label_ratio = u12_label_ratio(info['labels'])
        if add_type is True:
            if label_ratio > 0.5:
                type_id = 'u12'
            else:
                type_id = 'u2'
            self.type_id = type_id
        self.label_ratio = label_ratio

# /Classes ###################################################################

# Functions ##################################################################

def make_parser():
    parser = argparse.ArgumentParser(
        description='intronIC (intron Interrogator and Classifier) is a '
        'script which collects all of the annotated introns found in a '
        'genome/annotation file pair, and produces a variety of output '
        'files (*.iic) which describe the annotated introns and (optionally) '
        'their similarity to known U12 sequences.\n'
        'Without the \'-m\' flag, there MUST exist a matrix file in the '
        '\'data\' subdirectory in the same parent directory '
        'as intronIC.py, with filename \'scoring_matrices.fasta.iic\'. '
        'In the same data directory, there must also be a pair of sequence '
        'files (see --format_info) with reference intron sequences named '
        '\'[u2, u12]_reference_set.introns.iic\'',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    req_parse_grp = parser.add_argument_group(
        title='required arguments (-g, [-a, -b] | -q)')
    req_parse_grp.add_argument(
        '-g',
        '--genome',
        help='Genome file in FASTA format (gzip compatible)')
    req_parse_grp.add_argument(
        '-a',
        '--annotation',
        help='Annotation file in gff/gff3/gtf format (gzip compatible)')
    req_parse_grp.add_argument(
        '-b',
        '--bed',
        metavar='BED_FILE',
        help='Supply intron coordinates in BED format'
    )
    req_parse_grp.add_argument(
        '-n',
        '--species_name',
        type=str,
        help='Binomial species name, used in output file and intron label '
        'formatting. It is recommended to include at least the first letter '
        'of the species, and the full genus name since intronIC (by default) '
        'abbreviates the provided name in its output '
        '(e.g. Homo_sapiens --> HomSap)',
        required=True)
    req_parse_grp.add_argument(
        '-q',
        '--sequence_file',
        help='Provide intron sequences directly, rather than using a '
        'genome/annotation combination. Must follow the introns.iic '
        'format (see README for description)')
    parser.add_argument(
        '-f',
        '--feature',
        help='Specify feature to use to define introns. By default, '
        'intronIC will identify all introns uniquely defined by '
        'both CDS and exon features. Under the default mode, introns '
        'defined by exon features only will be demarcated by an '
        '\'[e]\' tag',
        choices=['cds', 'exon'],
        default=None)
    parser.add_argument(
        '-s',
        '--sequences_only',
        action='store_true',
        help='Bypass the scoring system and simply report the intron '
        'sequences present in the annotations')
    parser.add_argument(
        '--no_nc',
        action='store_true',
        help='Omit introns with non-canonical terminal dinucleoties from scoring')
    parser.add_argument(
        '-i',
        '--allow_multiple_isoforms',
        action='store_true',
        help='Include non-duplicate introns from isoforms other than '
        'the longest in the scored intron set (including those with alt. splice '
        'boundaries unless also -v')
    parser.add_argument(
        '-v',
        '--no_intron_overlap',
        action='store_true',
        help='Requires -i (or weird annotations). Exclude any introns with '
        'boundaries that overlap other introns from higher-priority transcripts '
        '(longer coding length, etc.). This will exclude, for example, introns '
        'with alternative 5′/3′ boundaries',
        default=False
    )
    parser.add_argument(
        '--pwms',
        metavar='{PWM file(s)}',
        help='One or more PWMs to use in place of the defaults. '
        'Must follow the formatting described by the --format_info '
        'option',
        nargs='+'
    )
    parser.add_argument(
        '--reference_u12s',
        '--r12',
        metavar='{reference U12 intron sequences}',
        help='introns.iic file with custom reference introns to be used '
        'for setting U12 scoring expectation, including flanking regions')
    parser.add_argument(
        '--reference_u2s',
        '--r2',
        metavar='{reference U2 intron sequences}',
        help='introns.iic file with custom reference introns to be used '
        'for setting U12 scoring expectation, including flanking regions')
    parser.add_argument(
        '--no_plot',
        action='store_true',
        help='Do not output illustrations of intron scores/distributions'
        '(plotting requires matplotlib)')
    parser.add_argument(
        '--format_info',
        action='store_true',
        help='Print information about the system '
        'files required by this script')
    parser.add_argument(
        '-d',
        '--include_duplicates',
        action='store_true',
        help='Include introns with duplicate '
        'coordinates in the intron seqs file')
    parser.add_argument(
        '-u',
        '--uninformative_naming',
        action='store_true',
        help='Use a simple naming scheme for introns instead of the '
        'verbose, metadata-laden default format'
    )
    parser.add_argument(
        '--no_abbreviate',
        '--na',
        action='store_true',
        help='Use the provided species name in full within the output files'
    )
    parser.add_argument(
        '-t',
        '--threshold',
        metavar='0-100',
        type=float,
        default=90,
        help='Threshold value of the SVM-calculated probability of being a U12 to '
        'determine output statistics')
    parser.add_argument(
        '--no_sequence_output',
        '--ns',
        action='store_true',
        help='Do not create a file with the full intron sequences '
        'of all annotated introns')
    parser.add_argument(
        '--five_score_coords',
        '--5c',
        default=(-3, 9),
        metavar=('start', 'stop'),
        nargs=2,
        type=int,
        help=(
            'Coordinates describing the 5\' sequence to be scored, relative to '
            'the 5\' splice site (e.g. position 0 is the first base of the '
            'intron); half-closed interval [start, stop)'
        )
    )
    parser.add_argument(
        '--three_score_coords',
        '--3c',
        default=(-10, 4),
        metavar=('start', 'stop'),
        nargs=2,
        type=int,
        help=(
            'Coordinates describing the 3\' sequence to be scored, relative to '
            'the 3\' splice site (e.g. position -1 is the last base of the '
            'intron); half-closed interval (start, stop]'
        )
    )
    parser.add_argument(
        '--branch_point_coords',
        '--bpc',
        default=(-55, -5),
        metavar=('start', 'stop'),
        nargs=2,
        type=int,
        help=(
            'Coordinates describing the region to search for branch point '
            'sequences, relative to the 3\' splice site (e.g. position -1 is the '
            'last base of the intron); half-closed interval [start, stop).'
        )
    )
    parser.add_argument(
        '-r',
        '--scoring_regions',
        help='Intron sequence regions to include in intron score calculations.',
        default=('five', 'bp'),
        choices=('five', 'bp', 'three'),
        nargs='+'
    )
    parser.add_argument(
        '--abbreviate_filenames',
        '--afn',
        action='store_true',
        help='Use abbreviated species name when creating '
        'output filenames.'
    )
    parser.add_argument(
        '--recursive',
        action='store_true',
        help=(
            'Generate new scoring matrices and training data using '
            'confident U12s from the first scoring pass. This option may '
            'produce better results in species distantly related to the '
            'species upon which the training data/matrices are based, though '
            'beware accidental training on false positives. Recommended only '
            'in cases where clear separation between types is seen with default '
            'data.'
        )
    )
    parser.add_argument(
        '--n_subsample',
        default=0,
        type=int,
        help=(
            'Number of sub-samples to use to generate SVM classifiers; 0 uses the '
            'entire training set and should provide the best results; otherwise, '
            'higher values will better approximate the entire set at the expense '
            'of speed.'
        )
    )
    parser.add_argument(
        '--cv_processes',
        # default=1,
        type=int,
        help=(
            'Number of parallel processes to use during cross-validation')
    )
    parser.add_argument(
        '-p',
        '--processes',
        default=1,
        type=int,
        help=(
            'Number of parallel processes to use for scoring (and '
            'cross-validation, unless --cv_processes is also set)'
        )
    )
    parser.add_argument(
        '--matrix_score_info',
        action='store_true',
        help='Produce additional per-matrix raw score information for each intron'
    )
    parser.add_argument(
        '-C',
        '--hyperparameter_C',
        default=None,
        help=(
            'Provide the value for hyperparameter C directly (bypasses optimized '
            'parameter search)')
    )
    parser.add_argument(
        '--min_intron_len',
        default=30,
        type=int,
        help='Minimum intron length to consider for scoring'
    )
    parser.add_argument(
        '--pseudocount',
        default=0.0001,
        type=float,
        help='Pseudocount value to add to each matrix value to avoid 0-div errors'
    )
    parser.add_argument(
        '--exons_as_flanks',
        action='store_true',
        help='Use entire up/downstream exonic sequence as flank sequence in output'
    )

    return parser


def check_thresh_arg(t):
    """
    Used in argument parsing to reject malformed threshold values

    """
    t = float(t)
    if not 0 <= t <= 100:
        raise argparse.ArgumentTypeError("'{}' is not within the range 0-100".
                                         format(t))
    return t


def abbreviate(species_name, n=3, separator=''):
    """
    Make a truncated binomial name from a longer
    version. If <n> is an int, both genus and 
    species will be truncated by the same amount.
    Otherwise, <n> is assumed to be an iterable of 
    integers for (genus, species).

    abbreviate(Homo_sapiens, n=3) --> HomSap
    abbreviate(Homo_sapiens, n=(3, 4)) --> HomSapi

    If n=None, will keep full species suffix and
    only abbreviate first name, e.g.

    Gallus_gallus --> GGallus

    """
    species_name = species_name.lower()
    bits = re.split(r"\W|_", species_name)
    if len(bits) == 1:  # no special character in string
        bit = bits[0]
        bits = [bit[0], bit[1:]]
    if type(n) == int:
        genus_n = n
        sp_n = n
    elif n == None:
        genus_n = 1
        sp_n = None
    else:
        genus_n, sp_n = n
    genus = bits[0][:genus_n]
    genus = genus[0].upper() + genus[1:]
    sp = bits[1]
    species = sp[0].upper() + sp[1:sp_n]
    abbreviated = separator.join([genus, species])

    return abbreviated


def load_external_matrix(matrix_file):
    """
    Makes matrices using the following file format:
    >maxtrix_name   start={n}
    A       C       G       T
    0.29	0.24	0.26	0.21
    0.22	0.33	0.19	0.26
    0.27	0.18	0.04	0.51
    0	    0	    1	    0
    0	    0	    0	    1
    0.99	0.004	0.001	0.005
    0.004	0.001	0.005	0.99
    0.004	0.99	0.001	0.005
    0.004	0.99	0.001	0.005
    0.009	0.02	0.001	0.97
    0.01	0.06	0.02	0.91
    0.05	0.22	0.08	0.65
    0.19	0.29	0.2	    0.32

    Returns a dictionary of the format
    {matrix_classification: {base: [frequencies]}}

    """
    def __name_parser(matrix_name):
        """
        Will attempt to interpret a matrix name within the matrix
        file header for categorization into a tuple, describing
        region and intron type (e.g. ("u12", "bp", "gtag")

        Returns a tuple of (subtype, region, boundaries)

        """
        subtypes = {
            "u12": ["u12", "12", "minor"],
            "u2": ["u2", "major"]
        }
        regions = {
            "five": ["five", "5"],
            "bp": ["bp", "branch-point"],
            "three": ["three", "3"]
        }
        boundaries = {
            "atac": ["at-ac", "atac"],
            "gtag": ["gt-ag", "gtag"],
            "gcag": ["gc-ag", "gcag"]
        }
        name_bits = []
        for cat in [subtypes, boundaries, regions]:
            bit = next(k for k, v in cat.items() if
                       any(subv in matrix_name.lower() for subv in v))
            name_bits.append(bit)

        # add an optional version tag to the name if present
        matrix_version = re.findall('[^A-Za-z]v\.?([^_\s]+)', matrix_name)
        if matrix_version:
            name_bits.append(matrix_version[0])

        return tuple(name_bits)

    matrices = {}
    for name, rows in fasta_parse(
        matrix_file,
        separator=None,
        trim_header=False
    ):
        try:
            start_bit = next(e for e in name.split() if 'start=' in e)
            start_index = int(start_bit.split('=')[1])
        except StopIteration:
            start_index = 0
        formatted_name = __name_parser(name.split()[0])
        matrices[formatted_name] = defaultdict(dict)
        # first row is bases in order
        bases = [b for b in rows.pop(0).split() if b in 'AGCT']
        base_index = {}
        for i, r in enumerate(rows, start=start_index):
            freqs = [float(f) for f in r.split()]
            for base, freq in zip(bases, freqs):
                matrices[formatted_name][base][i] = freq

    return matrices


def add_pseudos(matrix, pseudo=0.0001):
    """
    Apply a pseudo-count of ^pseudo to every frequency in a
    frequency matrix (to every number in each keys' values)

    """
    with_pseudos = {}
    for key, value in matrix.items():
        if isinstance(value, dict):
            with_pseudos[key] = add_pseudos(value, pseudo)
        else:
            # with_pseudos[key] = [(float(f) + pseudo) for f in value]
            with_pseudos[key] = float(value) + pseudo
            
    return with_pseudos


def average_matrices(a, b):
    outerdict = {}
    for key, a_bases in a.items():
        if key not in b:
            outerdict[key] = a_bases
            continue
        matrix = defaultdict(dict)
        for base, a_freqs in a_bases.items():
            b_freqs = b[key][base]
            for position, a_freq in sorted(a_freqs.items()):
                try:
                    b_freq = b_freqs[position]
                    avg_freq = (a_freq + b_freq) / 2
                    matrix[base][position] = avg_freq
                except KeyError:
                    continue
        outerdict[key] = matrix

    return outerdict


def format_matrix(matrix, label="frequencies", precision=None):
    """
    Formats the contents of a matrix in FASTA
    format, with the first line being the
    order of characters, {index_order}, and the
    following lines containing the frequency of
    each character; each position in the sequence
    has its own line in the matrix.

    {label} is used as the header for the frequencies
    entry.

    example output:

    >{index_order}
    A C G T
    >{label}
    0.25 0.5 0.25 0.0
    0.0 0.75 0.25 0.0
    ...

    """
    string_list = []
    characters = sorted(matrix.keys())
    freq_index = defaultdict(list)
    for character, frequencies in sorted(matrix.items()):
        for i, e in sorted(frequencies.items()):
            if precision:
                e = round(e, precision)
            freq_index[i].append(str(e))
    start_index = min(freq_index.keys())
    char_order = '\t'.join(characters)
    # string_list.append(">index_order\n{}".format(character_order))
    string_list.append(">{}\tstart={}\n{}".format(label, start_index, char_order))
    for i, freqs in sorted(freq_index.items()):
        string_list.append('\t'.join(freqs))
    return '\n'.join(string_list)


def introns_from_flatfile(
    flatfile,
    args,
    hashgen=False,
    type_id=None
):
    """
    Build a list of Intron objects from a reference file in
    introns.iic format

    """
    FIVE_SCORE_COORDS = args['FIVE_SCORE_COORDS']
    THREE_SCORE_COORDS = args['THREE_SCORE_COORDS']
    MIN_INTRON_LENGTH = args['MIN_INTRON_LENGTH']
    allow_noncanon = args['ALLOW_NONCANONICAL']
    allow_overlap = args['ALLOW_OVERLAP']
    BP_COORDS = args['BP_REGION_COORDS']
    BP_MATRIX_LENGTH = args['BP_MATRIX_LENGTH']

    ref_introns = []
    auto_name = 'auto_intron_'
    auto_int = 0
    with flex_open(flatfile) as flat:
        for line in flat:
            if line.startswith('#'):
                continue
            bits = line.strip().split('\t')
            name = bits[0]
            # ignore any columns that might be old format-style with
            # scores/lengths in them
            others = [b for b in bits[1:] if not b[0].isdigit()]
            try:
                five, int_seq, three = others[:3]
            except:
                sys.exit('Malformed line in file:\n{}'.format(line))
            if not name or name in '-.*_':
                name = auto_name + str(auto_int)
                auto_int += 1
            # make intron coordinates relative to the combined seq
            intronic_start = len(five) + 1
            intronic_stop = len(int_seq) + len(five)
            flank_size = max(len(five), len(three))
            seq = five + int_seq + three
            new_intron = Intron(
                name=name,
                start=intronic_start,
                stop=intronic_stop
            )
            # Set sub-sequence attributes for each intron object
            new_intron = assign_seqs(
                new_intron,
                seq.upper(),
                flank_size,
                FIVE_SCORE_COORDS,
                THREE_SCORE_COORDS,
                # five_score_length,
                BP_COORDS
            )
            new_intron.omit_check(
                MIN_INTRON_LENGTH,
                BP_MATRIX_LENGTH,
                allow_noncanon=allow_noncanon,
                allow_overlap=allow_overlap
            )
            # if hashgen:
            #     new_intron.sha1 = sha1(
            #         new_intron.seq.encode('utf-8')).digest()
            new_intron.seq = None

            if type_id:
                new_intron.type_id = type_id

            yield new_intron


def get_reference_introns(ref_file, args, type_id=None):
    """
    Build a list of Intron objects from a reference file in
    introns.iic format

    """
    refs = introns_from_flatfile(
        ref_file,
        args,
        type_id=type_id)
    
    refs = list(refs)

    ref_introns = [i for i in refs if not i.omitted]

    omitted_refs = len(refs) - len(ref_introns)

    if omitted_refs:
        write_log(
            '{} reference introns omitted from {}',
            omitted_refs, ref_file)

    return ref_introns


def make_feat_instance(line_info, feat_type=None):
    """
    Takes a GxfParse instance and returns a feature-specific object

    """
    if feat_type is None:
        feat_type = line_info.feat_type.lower()
    containers = {
        "gene": Gene,
        "transcript": Transcript,
        "exon": Exon,
        "cds": Exon
    }
    # default to Transcript class if it's not an obvious feature
    if feat_type not in containers:
        feat_type = "transcript"
    feats = []
    # Get standard feature info
    for p in line_info.parent:
        info = {
            "name": line_info.name,
            "feat_type": feat_type,
            "parent": p,
            "region": line_info.region,
            "strand": line_info.strand,
            "start": line_info.start,
            "stop": line_info.stop,
            "line_number": line_info.line_number,
            "phase": line_info.phase
        }
        if line_info.phase is None:
            info['phase'] = '.'
        # Put each type of data in the right container
        C = containers[feat_type]  # Gene, Transcript, Exon
        # Initialize new instance
        new_feat = C(**info)  # dict unpacking FTW
        feats.append(new_feat)

    return feats


def has_feature(annotation, check_types, column=2):
    """
    Checks annotation file for presence of check_type
    in column position of each line

    Returns True at first match found, False if
    check_type not found on any line

    """
    with flex_open(annotation) as f:
        for l in f:
            if l.startswith("#"):
                continue
            bits = l.strip().split("\t")
            if len(bits) < 9:
                continue
            feat_type = bits[column].lower()
            if feat_type in check_types:
                return True
    return False


def consolidate(*levels):
    """
    Arranges list of object dictionaries levels into a
    hierarchical data structure, based upon the attribute
    "parent"

    Must pass arguments in ascending order (e.g. children,
    parents, grandparents)

    """
    # consolidated = list(levels)  # modify copy
    consolidated = levels
    for i, lvl in enumerate(levels):
        try:
            parents = consolidated[i + 1]
        except IndexError:
            return consolidated[i]
        for name, obj in lvl.items():
            p = obj.parent
            if p in parents:
                parents[p].children.append(obj)
                obj.grandparent = parents[p].parent

    return consolidated


# Build a hierarchy of objects to allow for object functions to operate
# correctly (won't work if everything is flattened)
def annotation_hierarchy(annotation_file, child_feats):
    """
    Build an object heirarchy of gene/transcript/child_feats
    entries from an annotation file and chosen feature(s).
    Examines entire file in case features are not already
    grouped in hierarchical structure.

    child_feats may be multiple features, which will all be
    assigned as children of Transcript entries.


    Returns a list of aggregate objects, as well as a
    dictionary of all objects indexed by name.

    """
    genes = {}
    transcripts = {}
    children = defaultdict(list)
    destination = {"gene": genes, "transcript": transcripts}
    destination.update({c: children for c in child_feats})
    # to_collect = parent_types + list(child_feats)
    parent_containers = {"gene": Gene, "transcript": Transcript}
    parent_containers.update({c: Transcript for c in child_feats})
    # Track parents that don't meet criteria to avoid adding their children
    # to list of objects
    # # Keep record of unique coordinates to allow setting of dupe flag
    # # if dupe children are detected
    unique_coords = defaultdict(set)
    with flex_open(annotation_file) as f:
        # Track line numbers for tagging downstream objects with
        # their source lines
        parent_map = {}
        for ln, l in enumerate(f):
            try:
                line_info = GxfParse(l, ln)
            except TypeError:  # not a proper annotation line
                continue
            # will return as many objects as there are parents
            all_feats = make_feat_instance(line_info)
            for new_feat in all_feats:
                parent = new_feat.parent
                # check to ignore duplicate child entries
                if new_feat.feat_type not in child_feats:  # can be None
                    # if new_feat.feat_type.lower() in ('cds', 'exon'):  # why?
                    #     continue
                    parent_map[new_feat.name] = line_info
                    # parent_map[new_feat.name] = new_feat  ###???
                    if parent is not None and parent not in parent_map:
                        parent_map[parent] = Gene(name=parent, parent=[None])
                        # parent_map[parent] = line_info
                    continue
                elif parent:
                    check_coords = (
                        new_feat.feat_type,
                        new_feat.start,
                        new_feat.stop)
                    if check_coords not in unique_coords[parent]:
                        children[new_feat.compute_name()] = new_feat
                        unique_coords[parent].add(check_coords)

    parent_dest = transcripts
    grandparent_dest = genes
    # now make parent objs
    for name, child in children.items():
        parent = child.parent
        # grandparent = child.parent
        if parent not in parent_dest:
            try:
                parent_info = parent_map[parent]
                parent_objs = make_feat_instance(parent_info, 'transcript')
                for p in parent_objs:
                    parent_name = p.name
                    parent_dest[parent_name] = p
                    grandparent = p.parent
                    gp_info = parent_map[grandparent]
                    gp_objs = make_feat_instance(gp_info, 'gene')
                    for gp in gp_objs:
                        gp_name = gp.name
                        grandparent_dest[gp_name] = gp

            except KeyError:  # there was no parent line in gff
                # without parent line, use transcript name for
                # both gene and transcript
                parent_obj = Transcript(name=parent)
                gp_obj = Gene(name=parent)
                grandparent = parent
                # make parent and grandparent objs in respective
                # containers
                parent_dest[parent] = parent_obj
                grandparent_dest[grandparent] = gp_obj

    # Now collapse everything down to topmost objects
    collected = [e for e in [children, transcripts, genes] if e]
    consolidated = consolidate(*collected)
    if not consolidated:
        write_log(
            '[!] ERROR: could not establish parent-child relationships '
            'among feature set. Check annotation file format. Exiting now.'
        )
        sys.exit()
    top_level_objs = list(consolidated.values())
    # generate sibling numbers
    for obj in top_level_objs:
        obj.set_family_size()
        obj.coding_length = obj.get_coding_length()

    return top_level_objs


#TODO finish this function to allow creation of transcript/gene info file
# def high_level_metadata(obj_hierarchy, file_name):
#     """
#     Gets metadata about high-level genome features (genes, transcripts).
#     Returns a generator of formatted info strings.

#     """
#     for obj in obj_hierarchy:


def flatten(obj_list, all_objs=None, feat_list=None):
    """
    Takes a list of top-level objects and creates a dictionary
    of all feature objects in the list, including children.

    If feat_list, only returns objects of with >feat_type<s
    present in >feat_list<.

    """
    if all_objs is None:
        all_objs = {}
    for obj in obj_list:
        try:
            all_objs.update(
                flatten(obj.children, all_objs, feat_list)) # oooh recursion
            if feat_list:
                if obj.feat_type not in feat_list:
                    continue
            all_objs[obj.name] = obj
        except AttributeError:
            if feat_list:
                if obj.feat_type not in feat_list:
                    continue
            all_objs[obj.name] = obj

    return all_objs


def collect_introns(objs, feat_types, flat_annots):
    """
    Creates a dictionary of all intron objects from a
    list of higher-level objects.

    """
    introns = defaultdict(list)
    total_count = 0
    intron_defining_feats = ('cds', 'exon')  # this order matters in get_introns()
    for o in objs:
        # use all features to get introns initially, then filter down
        # to desired set afterward to allow family size calculation
        # to be consistent regardless of feature type being used
        new_introns = o.get_introns(intron_defining_feats, flat_annots)
        new_introns = [ni for ni in new_introns if ni.defined_by in feat_types]
        for i in new_introns:
            total_count += 1
            introns[i.region].append(i)

    return introns, total_count


def coord_overlap(coords, coord_list, presorted=False):
    """
    Check a list of integer coordinate tuples
    >coord_list< for overlap with the tuple
    >coords<. Returns the first overlapping
    tuple from >coord_list<, or False if no
    overlapping tuples found.

    """
    if not coord_list:
        return False
    if presorted is False:
        coord_list = sorted(coord_list)
    starts, stops = zip(*coord_list)
    c_start, c_stop = coords
    start_idx = bisect_left(stops, c_start)
    stop_idx = bisect_right(starts, c_stop)
    if start_idx != stop_idx:
        overlap_idx = min(start_idx, stop_idx)
        return coord_list[overlap_idx]
    else:
        return False


def add_tags(
    intron, intron_index,
    longest_isoforms, allow_overlap=False, longest_only=True
):
    """
    Add >duplicate<, >overlap< and >longest_isoform< tags
    to an intron.

    """
    # use intron omitted status to build parallel dictionary 
    # of only omitted introns
    region_id = (intron.region, intron.strand, not intron.omitted)

    i_coords = tuple([getattr(intron, f) for f in ['start', 'stop']])
    region_idx = intron_index[region_id]
    seen_coords = list(region_idx.keys())
    # is it a duplicate?
    if i_coords not in region_idx:
        intron.duplicate = False
        region_idx[i_coords]['max'] = intron.parent_length
        # region_idx[i_coords]['index'] = index
        region_idx[i_coords]['fam'] = intron.family_size
        region_idx[i_coords]['unique_num'] = intron.unique_num
    else:
        intron.duplicate = region_idx[i_coords]['unique_num']
        intron.overlap = intron.duplicate
    # is it in longest isoform?
    parent = intron.parent
    gene = intron.grandparent
    if gene not in longest_isoforms:
        longest_isoforms[gene] = parent
        intron.longest_isoform = True
    # elif intron.duplicate or longest_isoforms[gene] != parent:
    elif longest_isoforms[gene] != parent:
        intron.longest_isoform = False
    else:
        intron.longest_isoform = True

    # only check for overlap if intron is not a duplicate
    if intron.duplicate is False and not intron.omitted:
        if (intron.longest_isoform is not True and
            allow_overlap is False and longest_only is False):
            # iff the intron isn't in the longest isoform, check if its
            # coords overlap any other intron's coords and tag it
            # accordingly
            overlap = coord_overlap(i_coords, seen_coords)
            if overlap:
                overlap = region_idx[overlap]['unique_num']
                intron.overlap = overlap
            else:
                intron.overlap = False

    return intron, intron_index, longest_isoforms


def intronator(exons):
    """
    Builds introns from pairs of Exon objects
    in >exons<.

    sort_attr is used for ordering of introns.

    """

    def _window(seq):
        """
        Taken from https://docs.python.org/release/2.3.5/lib/itertools-example.html

        Returns a sliding window (of width n) over data from the iterable
        s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
        """
        n = 2
        it = iter(seq)
        result = tuple(islice(it, n))
        if len(result) == n:
            yield result
        for elem in it:
            result = result[1:] + (elem,)
            yield result

    exons = sorted_in_coding_direction(exons)

    for index, pair in enumerate(_window(exons)):
        # edge case where exons might overlap in same transcript/gene;
        # causes problems downstream but parsing later requires more
        # work and it's rare enough that warning user should suffice
        pair_coords = [(p.start, p.stop) for p in pair]
        if overlap_check(*pair_coords):
            overlap_log(pair)
            # don't create intron from overlapping features
            continue
        new_intron = Intron.from_exon_pair(*pair)
        # Record up- and downstream exon names in case needed
        # for annotation file modification later
        us_ex, ds_ex = [ex.name for ex in pair]
        new_intron.upstream_exon = us_ex
        new_intron.downstream_exon = ds_ex

        # moved to get_introns()
        # new_intron.fractional_position = frac_positions[index]

        yield new_intron


def overlap_log(objs):
    parent = objs[0].parent
    feature = objs[0].feat_type
    coords = [(o.start, o.stop) for o in objs]
    write_log(
        '[!] WARNING: overlapping {} features found in {}: {} - skipping',
        feature, parent, coords)


def overlap_check(a, b):
   """
   Check to see if the two ordered tuples, >a< and >b<,
   overlap with one another.

   By way of Jacob Stanley <3
   """

   val = (a[0] - b[1]) * (a[1] - b[0])

   if val < 0:
       return True
   else:
       return False


# def overlap_check(a, b):
#     """
#     >a< and >b< are sorted coordinate tuples;
#     returns True if they overlap, False otherwise

#     """
#     lowest = min([a, b], key=lambda x: min(x))
#     highest = max([a, b], key=lambda x: min(x))
#     if min(highest) <= max(lowest):
#         return True
#     else:
#         return False


# TODO use desired feature type (if specified) to decide the
# sort order of the object list
def sorted_in_coding_direction(obj_list):
    """
    Sorts a list of GenomeFeature objects by their
    stop/start attribute, depending on the value of
    their strand attributes.

    """
    strands = set([o.strand for o in obj_list])
    rev = False
    order_attr = None
    if len(strands) > 1:  # can't pick single orientation
        # prioritize CDS features for strand info if mixed strand
        # features present. If no CDS, default to file order
        cds_strands = set(
            [o.strand for o in obj_list if 
            any(a == 'cds' for a in (o.feat_type, o.defined_by))])
        if cds_strands and len(cds_strands) == 1:
            strand = cds_strands.pop()
            order_choice = (
                'strand of CDS features ({}) in parent feature'.format(strand))
        else:
            strand = False
            order_choice = 'order of features in source file (line number)'
        write_log(
            "WARNING: mixed strands found in provided feature list "
            "(parent feature: '{}'); "
            "defaulting to {}", 
            obj_list[0].parent, order_choice, wrap_chars=None)
    else:
        strand = strands.pop()
    if strand == "-":
        # Want to sort by "first" coord, which for reverse-
        # strand features is the stop
        order_attr = "stop"
        rev = True
    elif strand == False:
        order_attr = 'line_number'
    else:  # this defaults any non-(+,-) features to sort by start coord
        order_attr = "start"

    return sorted(obj_list, key=attrgetter(order_attr), reverse=rev)


def sliding_window(seq, n):
    """
    Taken from https://docs.python.org/release/2.3.5/lib/itertools-example.html

    Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...

    """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def seq_score(seq, matrix, start_index=0):
    """
    Score >seq< using values from >matrix<.

    Returns a float.

    """
    score = None
    for i, e in enumerate(seq, start=start_index):
        if score is None:
            score = matrix[e][i]
        else:
            score *= matrix[e][i]

    return score


# def bp_score(seq, matrix):
#     """
#     Score every sub-sequence of >seq< with length equal
#     to the value of the keys in >matrix<.

#     Returns the highest score achieved by any sub-sequence,
#     the relative coords of that sub-sequence within >seq<,
#     and the sub-sequence itself.

#     """
#     # If the matrix has different lengths for the value of any key,
#     # use the shortest
#     window_size = matrix_length(matrix)
#     start = 0
#     stop = window_size
#     best_score = None
#     best_coords = None
#     best_seq = None
#     for sub_seq in sliding_window(seq, window_size):
#         # convert from tuple to string
#         sub_seq = ''.join(sub_seq)
#         new_score = seq_score(sub_seq, matrix)
#         new_coords = (start, stop)
#         if best_score is None:
#             best_score = new_score
#             best_coords = new_coords
#             best_seq = sub_seq
#         else:
#             if new_score > best_score:
#                 best_score = new_score
#                 best_coords = new_coords
#                 best_seq = sub_seq
#         # Adjust window coordinates for next round
#         start += 1
#         stop += 1
#     return best_score, best_coords, best_seq


# bp score version with position-specific score weighting
# for secret species
def bp_score(seq, matrix, use_bpx=False, BPX=None, matrix_tag='TTTGA'):
    """
    Score every sub-sequence of >seq< with length equal
    to the value of the keys in >matrix<.

    Returns the highest score achieved by any sub-sequence,
    the relative coords of that sub-sequence within >seq<,
    and the sub-sequence itself.

    """
    # If the matrix has different lengths for the value of any key,
    # use the shortest
    window_size = matrix_length(matrix)
    start_index = min(matrix['A'].keys())
    start = 0
    stop = window_size
    best_score = None
    best_coords = None
    best_seq = None
    best_bpx_mod = None
    for sub_seq in sliding_window(seq, window_size):
        # convert from tuple to string
        sub_seq = ''.join(sub_seq)
        if not valid_chars(sub_seq):
            start += 1
            stop += 1
            continue
        # flag to send back if branch point score was modified
        # using BPX
        bpx_mod = None
        new_score = seq_score(sub_seq, matrix, start_index=start_index)
        # calculate the distance of the end of the motif from
        # the 3' end of the full intron using the location of
        # the end of the bp region and the window's current
        # end coordinate
        if BPX and use_bpx is True and 'TTGA' in matrix_tag:
            # TODO: fix this global reference
            dist_from_3 = (len(seq) - stop) + abs(BP_REGION_COORDS[1])
            # perform multiplier if present in dictionary
            multiplier = BPX.get(dist_from_3)
            if multiplier is not None:
                # bpx_mod = ( (multiplier - 1) / 2 ) + 1
                bpx_mod = multiplier
                # bp score adjustment accounting for neg. initial scores
                delta = abs(new_score * (bpx_mod - 1))
                new_score += delta
        new_coords = (start, stop)
        if best_score is None or new_score > best_score:
            best_bpx_mod = bpx_mod
            best_score = new_score
            best_coords = new_coords
            best_seq = sub_seq
        # Adjust window coordinates for next round
        start += 1
        stop += 1

    return best_score, best_coords, best_bpx_mod, best_seq


def matrix_from_seqs(seqs, start_index=0):
    """
    Constructs a position-weight matrix from one or
    more sequences.

    Returns a dictionary of {character: [f1, f2...fn]},
    where fn is the character frequency at the nth
    position in the sequence.

    """
    matrix = defaultdict(list)
    characters = set(['G', 'T', 'A', 'C'])
    lengths = set()
    n_seqs = 0

    for s in seqs:
        lengths.add(len(s))
        n_seqs += 1
        for i, e in enumerate(s, start=start_index):
            matrix[i].append(e)
            characters.add(e)
    
    if n_seqs == 0:
        return {}, 0

    character_order = sorted(characters)
    seq_length = min(lengths)

    frequencies = defaultdict(dict)

    for i in range(start_index, seq_length - abs(start_index)):
        chars = matrix[i]
        freqs = []
        for c in character_order:
            n = chars.count(c)
            f = n / n_seqs
            freqs.append(f)
            frequencies[c][i] = f

    return frequencies, n_seqs


def max_min_matrix_values(matrix):
    max_by_pos = defaultdict(float)
    min_by_pos = defaultdict(float)
    for char, pos_dict in matrix.items():
        for pos, v in pos_dict.items():
            if v >= max_by_pos.get(pos, 0):
                max_by_pos[pos] = v
            if v <= min_by_pos.get(pos, 1):
                min_by_pos[pos] = v

    return max_by_pos, min_by_pos


def mark_seq_score(seq, matrix, start_index=0):
    max_values, min_values = max_min_matrix_values(matrix)
    relative_scores = []
    for i, c in enumerate(seq, start=start_index):
        score = matrix[c][i]
        max_score = max_values.get(i)
        min_score = min_values.get(i)
        if score == max_score:
            rel_score = '*'
        elif score == min_score:
            rel_score = '/'
        else:
            rel_score = math.floor((score / max_values.get(i)) * 10)
        relative_scores.append(str(rel_score))

    return relative_scores


def valid_chars(seq, bad_chars=re.compile("[^ACTG]")):
    """
    Returns True if >seq< contains only characters in the set
    [ACTG], False otherwise.

    """
    if bad_chars.search(seq):
        return False
    else:
        return True


def canonical_bounds(intron):
    """
    Checks an intron sequence's dinucleotide boundaries
    for canonical splice sites. Returns False if boundaries
    are non-canonical.

    Splice site variants derived from table 1 in
    DOI: 10.1093/nar/gku744

    """
    canonical = {
        "AT": ["AC"], #"AC", "AG", "AA", "AT"
        "GT": ["AG"], #"GG", "AA", "TG", "AT"
        "GC": ["AG"],
        # "GG": ["AG"],
        # "GA": ["AG"],
        # "TT": ["AG"]
    }
    # Get terminal dinucleotides
    five, three = intron.dnts

    if five not in canonical:
        return False
    elif three not in canonical[five]:
        return False
    else:
        return True


def u12_correction(intron):
    """
    Checks >intron< for the presence of misannotated
    U12 boundaries. This requires a shift by the same number
    of nucleotides on both ends of the intron.

    If a correction can be made to produce a strong 5' boundary,
    the integer value of the shift will be placed in the intron's
    >corrected< attribute. Otherwise, >intron< is returned unchanged.

    """
    def _shift_phase(phase, shift):
        phases = deque([0, 1, 2])
        try:
            index = phases.index(int(phase))
        except ValueError:  # e.g. '.' for exons
            return phase
        phases.rotate(-shift)

        return phases[index]

    up_n = 5
    down_n = 12
    # strict_motif = re.compile(r'[AG]TATC[CT]{2}')
    strict_motif = re.compile(r'[AG]TATC[CT]([ACTG]T|T[ACTG])')
    lax_motif = re.compile(r'[AG]TATC[CT]')
    # relax constraints if we're correcting a non-canonical intron
    # if canonical_bounds(intron):
    #     motif = strict_motif
    # else:
    #     motif = lax_motif
    motif = strict_motif
    region = intron.upstream_flank[-up_n:] + intron.seq[:down_n]
    match = motif.search(region)
    if not match:
        return False
    match_index = match.start()
    shift = match_index - up_n
    if shift == 0:
        return False
    intron.corrected = shift
    intron.phase = _shift_phase(intron.phase, shift)
    intron.dynamic_tag.add('[c:{}]'.format(shift))
    if intron.strand == '-':  # reverse adjustment for neg. strand
        shift *= -1
    intron.start += shift
    intron.stop += shift

    return True


def correct_annotation(
    introns, 
    flattened, 
    annotation_file, 
    RUN_DIR, 
    FN_ANNOT
):
    """
    Adjusts intron-defining entries in >annotation_file<
    based upon the intron attribute >corrected<.

    If an intron's coordinates can be corrected with a -2 shift,
    for example, this would require shifting the stop coord
    of the 5' exon and the start coord of the 3' exon by -2 each.

    A ';shifted:[start, stop]:[1, 2, -1, -2]' tag is added to each
    line where such adjustments are made.

    >flattened< is a dictionary index of all objects in the annotation
    file indexed by name.

    Returns the modified filename.

    """

    # two functions to correct the phase of a sequence
    # whose start coordinate is changed by >shift< number
    # of nt
    def _shift_phase(phase, shift):
        phases = deque([0, 1, 2])
        try:
            index = phases.index(int(phase))
        except ValueError:  # e.g. '.' for exons
            return phase
        phases.rotate(-shift)

        return phases[index]

    def _correct_exon(exon, shift, side):
        """
        Corrects one of an exon object's coordinates,
        by amount >shift< and value of >side< (either
        5 or 3).

        Intended to function in the service of correcting
        exons to recapitulate canonical intron boundaries.

        """
        # Coordinate choice depends on which side of intron it's on
        if side == 5:
            target = 'stop'
            phase_shift = False
        else:
            target = 'start'
            # only flag for phase correction if start coord is modified
            phase_shift = True
        # Change goes the other way for negative strand coordinates
        if exon.strand == '-':
            shift = shift * -1
            # coords are always relative to + strand
            target = next(e for e in ('start', 'stop') if e != target)
        old_coord = getattr(exon, target)
        new_coord = old_coord + shift
        setattr(exon, target, new_coord)
        # Return extra values for use in _change_coords
        return exon, target, shift, phase_shift

    def _change_coords(line, coords, coord_tag, shift_tag, phase_shift):
        bits = line.strip().split('\t')
        insert = ';'
        if bits[8].endswith(';'):
            insert = ''
        mod_tag = '{}shift:{}:{}'.format(insert, coord_tag, shift_tag)
        bits[8] = bits[8] + mod_tag
        bits[3], bits[4] = map(str, coords)
        # correct phase if applicable
        if phase_shift is True:
            if bits[6] == '+':
                shift_tag *= -1
            # bits[7] = str(int(bits[7]) + (shift_tag * -1))
            bits[7] = str(_shift_phase(bits[7], shift_tag))
        return '\t'.join(bits)

    # Iterate over only those introns with non-zero values for {corrected}
    corrected_count = 0
    corrected_dupes = 0
    corrected_exon_coords = {}
    corrected_introns = [i for i in introns if i.corrected]
    corrected_count = len(corrected_introns)
    if corrected_count == 0:
        return
    # for intron in filter(attrgetter('corrected'), introns):
    for intron in corrected_introns:
        corrected_count += 1
        if intron.duplicate:
            corrected_dupes += 1
        shift = intron.corrected
        five_exon = flattened[intron.upstream_exon]
        three_exon = flattened[intron.downstream_exon]
        # Build index of corrected coords by line number
        for ex, intron_side in zip([five_exon, three_exon], [5, 3]):
            (cor_ex, coord_tag,
            shift_tag, phase_shift) = _correct_exon(ex, shift, intron_side)
            coords = (cor_ex.start, cor_ex.stop)
            corrected_exon_coords[ex.line_number] = (
                coords, coord_tag, shift_tag, phase_shift)

    # if corrected_count == 0:
    #     return
    # We now have an index of lines to modify in the annotation file
    modded_filename = FN_ANNOT
    modded_filepath = os.path.join(
        RUN_DIR, modded_filename)
    with flex_open(annotation_file) as infile, \
    open(modded_filepath, 'w') as outfile:
        for ln, l in enumerate(infile):
            if ln not in corrected_exon_coords:
                outfile.write(l)
                continue
            (new_coords, coord_tag,
            shift_tag, phase_shift) = corrected_exon_coords[ln]
            new_line = _change_coords(
                l, new_coords, coord_tag, shift_tag, phase_shift)
            outfile.write(new_line + '\n')
    write_log(
        '{} ({} unique, {} redundant) putatively misannotated U12 introns '
        'corrected in {}',
        corrected_count,
        corrected_count - corrected_dupes,
        corrected_dupes,
        modded_filename
    )


def write_log(string, *variables, wrap_chars='[]', level='critical'):  # level was info
    """
    Prints to screen and writes to log file a
    formatted string with >variables< surrounded
    by >wrap_chars<.

    e.g. write_log('The answer is {}', 42) --> 'The answer is [42]'

    If >wrap_chars< is None, variables are printed without
    additional formatting.

    """
    if wrap_chars is None:
        formatted_vars = variables
    else:
        open_char, close_char = wrap_chars
        formatted_vars = [
            '{}{}{}'.format(open_char, v, close_char) for v in variables
        ]
    formatted_string = string.format(*formatted_vars)
    logger = getattr(logging, level)
    logger(formatted_string)


def assign_seqs(
    intron,
    region_seq,
    int_flank_size,
    five_score_coords,
    three_score_coords,
    bp_coords
):
    """
    Assign sub-sequences to intron objects based on object metadata
    and provided variables in argument.

    Returns a modified intron object.

    """
    def _short_bp_adjust(intron, bpc, fsl):
        """
        Check for and adjust bp coords to correct
        for smaller initial bp search region in
        short introns (where original bp region coords
        might extend upstream of the 5' end of the intron)

        """
        ic = (intron.start, intron.stop)
        rev = False
        if intron.strand == '-':
            rev = True
            ic = ic[::-1]
            bpc = bpc[::-1]
            fsl *= -1
        bp_start, bp_stop = bpc
        int_start = ic[0]
        shift = abs((int_start + fsl) - bp_start)
        if rev:
            bp_start -= shift
            bpc = (bp_stop, bp_start)
        else:
            bp_start += shift
            bpc = (bp_start, bp_stop)
        return bpc

    upf, intron.seq, downf = intron.get_seq(region_seq, flank=int_flank_size)
    intron.upstream_flank = upf
    intron.downstream_flank = downf

    # scoring sequences
    us_length = len(intron.upstream_flank)
    ds_length = len(intron.downstream_flank)

    # adjust coords to account for flanking sequence
    scoring_seq = upf + intron.seq + downf
    five_rel_coords = [c + us_length for c in five_score_coords]
    three_rel_coords = [c - ds_length for c in three_score_coords]

    # remove zeroes to avoid indexing errors if range ends in 0
    five_rel_coords = [c if c != 0 else None for c in five_rel_coords]
    three_rel_coords = [c if c != 0 else None for c in three_rel_coords]

    # pull sequence corresponding to relative coordinate ranges
    intron.five_seq = scoring_seq[slice(*five_rel_coords)]
    intron.three_seq = scoring_seq[slice(*three_rel_coords)]

    # fixed five and three sequences for display purposes
    intron.three_display_seq = intron.seq[bp_coords[1]:]
    intron.five_display_seq = intron.seq[:10]
    intron.dnts = (intron.seq[:2], intron.seq[-2:])
    if not canonical_bounds(intron):
        intron.noncanonical = True
    else:
        intron.noncanonical = False

    # account for the 1-based indexing adjustment in get_seqs()
    # which should not apply to these kinds of relative coords
    bp_coords = (bp_coords[0] + 1, bp_coords[1])
    bp_region_coords = intron.get_rel_coords('three', bp_coords)

    # correct bp region if intron is short
    five_score_length = len([e for e in range(*five_score_coords) if e >=0])
    if intron.length < abs(bp_coords[0]) + five_score_length:
        bp_region_coords = _short_bp_adjust(
            intron, bp_region_coords, five_score_length)

    intron.bp_region_seq = intron.get_seq(region_seq, *bp_region_coords)

    return intron


def longest_match(seq, pattern=r'[ATCG]+'):
    """
    Takes a string, and returns the length of
    the longest stretch of characters matching
    {pattern}.

    """
    matches = re.findall(pattern, seq)
    if matches:
        return len(max(matches, key=len))
    else:
        return 0


def heirarchical_sort_attrs(intron):
    sort_features = [
        intron.defined_by,  # CDS before exon
        intron.parent_length * -1,
        intron.parent,
        intron.family_size * -1,
        intron.index,
        intron.line_number]

    return tuple(sort_features)


def get_sub_seqs(
    introns_by_region,
    flat_objs,
    genome,
    int_flank_size,
    five_score_coords,
    three_score_coords,
    bp_coords,
    exons_as_flanks=False
):
    """
    Generator that populates objects in >introns< with short
    sub-sequences using >genome<.

    >introns_by_region< is a dictionary of intron objects keyed
    with the >region< attribute.

    """
    # get the total number of FASTA headers we need to consider
    # to allow for early exit if there are a bunch of headers
    # without introns
    total_regions = len(introns_by_region)
    intron_index = defaultdict(lambda: defaultdict(dict))
    for region_name, region_seq in fasta_parse(flex_open(genome)):
        if region_name not in introns_by_region:
            continue
        longest_isoforms = {}
        region_seq = region_seq.upper()
        for intron in sorted(
            introns_by_region[region_name],
            key=heirarchical_sort_attrs):
            if exons_as_flanks is True:
                # determine length of upstream and downstream
                # defining features to use as flanking sequence
                up_ex = flat_objs[intron.upstream_exon]
                down_ex = flat_objs[intron.downstream_exon]
                up_length = abs(up_ex.stop - up_ex.start) + 1
                down_length = abs(down_ex.stop - down_ex.start) + 1
                int_flank_size = (up_length, down_length)
            intron = assign_seqs(
                intron,
                region_seq,
                int_flank_size,
                five_score_coords,
                three_score_coords,
                bp_coords)
            if intron.noncanonical:
                if u12_correction(intron):  # coords have changed
                    intron = assign_seqs(
                        intron,
                        region_seq,
                        int_flank_size,
                        # five_flank,
                        # five_score_length,
                        five_score_coords,
                        three_score_coords,
                        bp_coords)

            yield intron

        total_regions -= 1
        if total_regions == 0:  # we've used all the headers we need
            break


def write_format(obj, *attribs, fasta=True, separator='\t', null='.'):
    """
    Formats a set of object attributes into a string
    for writing to a file.

    If >fasta< is True, will use the first attribute
    as a header in FASTA format.

    Can retrieve method attributes, but not ones which
    require arguments.

    No trailing newline in returned string.

    """
    attribs = list(attribs)
    values = []
    for atr in attribs:
        try:
            value = getattr(obj, atr)
        except AttributeError:  # is simply a variable being passed in
            values.append(atr)
            continue
            # try them as functions; if that doesn't work, they're
            # already values
        try:
            value = value()
        except TypeError:
            pass
        values.append(value)
    if fasta is True:
        header = '>{}'.format(values.pop(0))
    if null is not None:
        values = [v if v is not None else null for v in values]
    if separator is not None:
        content = separator.join([str(v) for v in values])
        if fasta is True:
            content = '\n'.join([header, content])
    else:
        content = values

    return content


def counter_format_top(cntr, num=None):
    """
    Formats a Counter object >cntr< for printing its
    values in a numbered list, along with percentage
    statistics for each value.

    >num< is the number of entries that will be printed.
    If >num< is None, all values will be printed.

    Returns an iterator of formatted strings.

    """
    top = cntr.most_common(num)
    total = sum(cntr.values())
    for element, count in top:
        fraction = (count / total) * 100
        count_info = "* {} ({}/{}, {:.2f}%)".format(
            element, count, total, fraction)
        yield count_info


def build_u2_bp_matrix(introns, u12_matrices, spcs, simple_name, dnt_list=None):
    """
    Builds a matrix for branch point sequences
    based on highest-scoring sequences when
    scored using a u12 bp matrix.

    Returns a matrix dictionary.

    """
    def _iter_bps(ints, matrices, dnt_list=None):
        bp_seqs = []
        for intron in ints:
            if dnt_list is not None:
                if intron.dnts not in dnt_list:
                    continue
            bp_region_seq = intron.bp_region_seq
            
            ###!!!
            best_score = 0
            best_seq = None
            for name, matrix in matrices.items():
                m_score, *_, seq = bp_score(
                    bp_region_seq, matrix, use_bpx=True, matrix_tag=name[-1])
                if m_score > best_score:
                    best_seq = seq
                    best_score = m_score
            seq = best_seq
            ###!!!

            if not seq:
                print('NO BP SEQ: ', intron.get_name(spcs, simple_name))
                sys.exit(0)
            yield seq

    bp_seqs = _iter_bps(introns, u12_matrices, dnt_list)

    return matrix_from_seqs(bp_seqs)


def matrix_length(matrix):
    """
    Returns the shortest value length found in >matrix<

    """
    return min(set([len(vals) for vals in matrix.values()]))


def get_score_bounds(matrix):
    """
    Retrieves the highest and lowest possible scores for >matrix<,
    e.g. the score achieved if a sequence matched the highest-
    or lowest-frequency character at every position in the matrix.

    """
    position_indexed = defaultdict(list)
    # Collect frequencies at each position in lists
    for character, freqs in matrix.items():
        for i, f in freqs.items():
            position_indexed[i].append(f)
    # Calculate max and min values for each position
    min_freqs, max_freqs = [], []
    for position, freqs in position_indexed.items():
        min_freqs.append(min(freqs))
        max_freqs.append(max(freqs))

    # Total scores are products of each list
    min_score = np.prod(min_freqs)
    max_score = np.prod(max_freqs)

    return min_score, max_score


def multi_matrix_score(
    intron,
    matrices,
    FIVE_SCORE_COORDS,
    THREE_SCORE_COORDS,
    PSEUDOCOUNT,
    regions=('five', 'bp', 'three'),
    matrix_tags=None,
    use_bpx=False):
    """
    Finds the highest-scoring matrix key and value for the
    specified intron, using the sum of multiple regions for
    scoring.

    Returns a tuple of the score and associated matrix key.

    """
    region_map = {
        'five': 'five_seq',
        'bp': 'bp_region_seq',
        'three': 'three_seq'
    }
    score_funcs = {
        'five': partial(seq_score, start_index=FIVE_SCORE_COORDS[0]),
        'bp': partial(bp_score, use_bpx=use_bpx),
        'three': partial(seq_score, start_index=THREE_SCORE_COORDS[0])}
    score_info = defaultdict(lambda: defaultdict(dict))
    if matrix_tags is not None:
        matrices = {
            k: v for k, v in matrices.items()
            if all(t in k for t in matrix_tags)}
    else:
        matrices = {
            k: v for k, v in matrices.items()
            if any(r in k for r in regions)}
    for matrix_key, matrix in matrices.items():
        subtype, dnts, region, *_ = matrix_key
        matrix_category = (subtype, dnts)
        score_function = score_funcs[region]
        seq = getattr(intron, region_map[region])
        if region == 'bp':
            score = score_function(seq, matrix, matrix_tag=matrix_key[-1])
        else:
            score = score_function(seq, matrix)
        info = []
        if score is None:  # edge case where sequence is at end of contig
            s = PSEUDOCOUNT * len(matrix)
        elif region == 'bp':
            # bp score function returns additional information
            s = score[0]
            info = score[1:]
        else:
            s = score
        target = score_info[matrix_category][region]
        # check for existing score and only overwrite if current
        # score is better (multiple scores possible with versioned
        # matrices, e.g. VA9, VA10 bp matrices)
        if target.get('score', 0) > s:
            continue
        else:
            target['score'] = s
            target['info'] = info
            target['name'] = matrix_key

    return score_info


def best_matrix(score_info, scoring_regions, priority_tag=None):
    # summary_score_regions = ['five', 'bp', 'three']
    category_scores = []
    for matrix_category, regions in score_info.items():
        region_scores = []
        for r, r_score in regions.items():
            score = r_score['score']
            if r in scoring_regions:
                region_scores.append(score)
        summary_score = summarize(region_scores)
        category_scores.append((summary_score, matrix_category))

    # filter best scores by tag if present, otherwise simply return
    # highest-scoring matrix
    category_scores.sort(reverse=True)
    best_score_info = category_scores[0]
    if priority_tag:
        try:
            best_score_info = next(
                s for s in category_scores if priority_tag in s[1])
        except StopIteration:
            pass

    # Return top score and matrix category
    return best_score_info  # e.g. (score, (u12, gtag))


def log_ratio(a, b):
    """
    Computes the log base 2 score of a over b.

    """
    return math.log2(a / b)


# def get_raw_scores(introns, matrices, scoring_regions, multiplier_d=None):
#     """
#     Assigns scores to each intron in >introns<,
#     using available matrices in >matrices<.

#     Returns a generator of introns with raw scores added

#     """
#     scored_introns = []
#     for intron in introns:
#         # Determine the best type of matrix to used based on 5' score
#         # e.g. ('u12', 'gtag')
#         u2_matrix_info = multi_matrix_score(
#             intron, matrices, matrix_tags=['u2'])
#         u12_matrix_info = multi_matrix_score(
#             intron, matrices, matrix_tags=['u12'], use_bpx=True)
#         dnts = ''.join(intron.dnts).lower()
#         u2_score, best_u2_key = best_matrix(
#             u2_matrix_info, scoring_regions, dnts)
#         u12_score, best_u12_key = best_matrix(
#             u12_matrix_info, scoring_regions, dnts)
#         u12_bp_score = u12_matrix_info[best_u12_key]['bp']['score']
#         u12_bp_info = u12_matrix_info[best_u12_key]['bp']['info']
#         bp_rel_coords, bpm, u12_bp_seq = u12_bp_info

#         if bpm is not None:
#             intron.dynamic_tag.add('bpm={}'.format(bpm))

#         # Get log scores for each region
#         score_regions = ['five', 'bp', 'three']
#         for region, attr in zip(
#             score_regions,
#             ['five_raw_score', 'bp_raw_score', 'three_raw_score']
#         ):
#             try:
#                 u12_region_matrix = u12_matrix_info[best_u12_key][region]
#                 u2_region_matrix = u2_matrix_info[best_u2_key][region]
#                 ratio_score = log_ratio(
#                     u12_region_matrix['score'],
#                     u2_region_matrix['score'])
#                 setattr(intron, attr, ratio_score)
#                 intron.matrices.add(u12_region_matrix['name'])
#                 intron.matrices.add(u2_region_matrix['name'])
#             except:
#                 continue

#         intron.bp_seq = u12_bp_seq
#         intron.bp_relative_coords = bp_rel_coords

#         # record the name of the matrix used for each score type
#         intron.u2_matrix = best_u2_key  # e.g. ('u12', 'gtag')
#         intron.u12_matrix = best_u12_key

#         scored_introns.append(intron)
#         # yield intron
#     return scored_introns

def assign_raw_score(
    intron, 
    matrices, 
    scoring_regions, 
    five_score_coords, 
    three_score_coords, 
    pseudocount, 
    multiplier_d=None
):
    """
    Assigns raw scores to an Intron object based on the supplied matrices

    """
    # Determine the best type of matrix to used based on 5' score
    # e.g. ('u12', 'gtag')
    u2_matrix_info = multi_matrix_score(
        intron, 
        matrices, 
        five_score_coords, 
        three_score_coords, 
        pseudocount, 
        matrix_tags=['u2']
    )
    u12_matrix_info = multi_matrix_score(
        intron, 
        matrices, 
        five_score_coords, 
        three_score_coords, 
        pseudocount, 
        matrix_tags=['u12'],
        use_bpx=True
    )
    dnts = ''.join(intron.dnts).lower()
    u2_score, best_u2_key = best_matrix(
        u2_matrix_info, scoring_regions, dnts)
    u12_score, best_u12_key = best_matrix(
        u12_matrix_info, scoring_regions, dnts)
    u12_bp_score = u12_matrix_info[best_u12_key]['bp']['score']
    u12_bp_info = u12_matrix_info[best_u12_key]['bp']['info']
    bp_rel_coords, bpm, u12_bp_seq = u12_bp_info

    if bpm is not None:
        intron.dynamic_tag.add('bpm={}'.format(bpm))

    # Get log scores for each region
    score_regions = ['five', 'bp', 'three']
    for region, attr in zip(
        score_regions,
        ['five_raw_score', 'bp_raw_score', 'three_raw_score']
    ):
        try:
            u12_region_matrix = u12_matrix_info[best_u12_key][region]
            u2_region_matrix = u2_matrix_info[best_u2_key][region]
            ratio_score = log_ratio(
                u12_region_matrix['score'],
                u2_region_matrix['score'])
            setattr(intron, attr, ratio_score)
            intron.matrices.add(u12_region_matrix['name'])
            intron.matrices.add(u2_region_matrix['name'])
        except:
            continue

    intron.bp_seq = u12_bp_seq
    intron.bp_relative_coords = bp_rel_coords

    # record the name of the matrix used for each score type
    intron.u2_matrix = best_u2_key  # e.g. ('u12', 'gtag')
    intron.u12_matrix = best_u12_key

    return intron


def get_raw_scores(
    introns, 
    matrices, 
    scoring_regions, 
    five_score_coords, 
    three_score_coords, 
    pseudocount, 
    multiplier_d=None, 
    processes=1
):
    with Pool(processes=processes) as pool:
        raw_introns = pool.starmap(
            assign_raw_score, zip(
                introns, 
                repeat(matrices), 
                repeat(scoring_regions),
                repeat(five_score_coords),
                repeat(three_score_coords),
                repeat(pseudocount),
                repeat(multiplier_d)))
    return raw_introns


def unique_attributes(introns, a_list):
    seen_attrs = set()
    for i in introns:
        i_attrs = get_attributes(i, a_list)
        for ia in i_attrs:
            if ia not in seen_attrs:
                seen_attrs.add(ia)
                yield i


def make_z_score(raw_score, mean, stdev):
    try:
        z_score = (raw_score - mean) / stdev
        return z_score
    except RuntimeWarning:  # if stdev = 0
        return raw_score - mean


def mutate(intron, five_score_coords, three_score_coords):
    intron.dnts = ('GT', 'AG')
    five_before = intron.five_seq[:abs(five_score_coords[0])]
    five_after = intron.five_seq[abs(five_score_coords[0]) + 2:]
    intron.five_seq = five_before + 'GT' + five_after
    three_before = intron.three_seq[:abs(three_score_coords[0]) - 2]
    three_after = intron.three_seq[abs(three_score_coords[0]):]
    intron.three_seq = three_before + 'AG' + three_after

    return intron


def write_matrix_file(matrices, fn, precision=None):
    """
    Writes one or more matrices to a file, using
    keys as headers.

    """
    with open(fn, 'w') as outfile:
        for name, matrix in sorted(matrices.items()):
            label = '-'.join(name)
            # outfile.write('[#] {}\n'.format('-'.join(name)))
            outfile.write(format_matrix(matrix, label, precision) + '\n')


def min_max_from_bounds(u12_bounds, u2_bounds):
    u12_min, u12_max = u12_bounds
    u2_min, u2_max = u2_bounds

    min_score = log_ratio(u12_min, u2_max)
    max_score = log_ratio(u12_max, u2_min)

    return min_score, max_score


def rescale(old_score, oldmin, oldmax, newmin=1e-5, newmax=1):
    scaled_score = (
        ((old_score - oldmin) * (newmax - newmin)) /
        ((oldmax - oldmin)) + newmin)

    return scaled_score


def annotate(string, start, stop):
    annotated = string[:start] + '[' + string[start:stop] + ']' + string[stop:]

    return annotated


def motif_strings(intron):
    """
    Returns a schematic of all scored motifs, and the BPS context as strings

    """
    schematic = '{}|{}...{}...{}|{}'.format(
        intron.upstream_flank[-3:],
        intron.five_display_seq, intron.bp_seq,
        intron.three_display_seq, intron.downstream_flank[:3])
    bp_context = annotate(
        intron.bp_region_seq,
        *intron.bp_relative_coords) + intron.three_display_seq
    
    return schematic, bp_context


def dupe_list_format(intron, spcs, simple_name, dupe_index):
    """
    Check an intron against an index of duplicates, and
    return a string of all duplicate entries in list format,
    or None.

    """
    intron_uid = (intron.region, intron.start, intron.stop)
    if intron_uid in dupe_index:
        list_bits = [intron.region, intron.strand, intron.start, intron.stop]
        dupe_strings = []
        for dupe in dupe_index[intron_uid]:
            # dupe.svm_score = intron.svm_score
            dupe_name = '{}-{}'.format(spcs, dupe.get_label(simple_name))
            dupe_bits = list_bits + [dupe_name]
            dupe_string = '\t'.join(map(str, dupe_bits))
            dupe_strings.append(dupe_string)
        dupe_output = '\n'.join(dupe_strings)
    else:
        dupe_output = None

    return dupe_output


def flip_check(intron, flip_dict):
    """
    Checks an intron against a dictionary of introns with scores
    that did not survive boundary switching.

    If present, the score resulting from the boundary switch
    will be returned, otherwise None.

    """
    # name = intron.get_name()
    uid = intron.unique_num
    if uid in flip_dict:
    # if name in flip_dict:
        return flip_dict[uid]
        # return flip_dict[name]
    else:
        return None


def bp_multiplier(intron, multiplier_dict):
    """
    Checks the position of the branch point seq and adds a
    multiplier from multiplier_dict if the bp seq distance
    from the 3' end is in multiplier_dict.

    Returns None, or the appropriate multiplier.

    """
    bpr_length = len(intron.bp_region_seq)
    bp_stop = intron.bp_relative_coords[1]
    # Calculate distance from the 3' end of the intron
    d_from_3 = (bpr_length - bp_stop) + abs(BP_REGION_COORDS[1])

    return multiplier_dict.get(d_from_3)


def add_scores(introns, intron_seq_file, spcs):
    """
    Add scores to previously-made intron sequences file

    """
    tmp_file = '{}.tmp'.format(intron_seq_file)
    score_dict = {
        intron.name.split(';')[0]: intron.svm_score
        for intron in introns}
    with open(intron_seq_file) as inf, open(tmp_file, 'w') as outf:
        for l in inf:
            bits = l.strip().split('\t')
            # bits[0] is intron name
            name = bits[0].split(';')[0]
            if name in score_dict:
                score = str(round(score_dict[name], 4))
                bits[1] = score
            outf.write('\t'.join(bits) + '\n')
    os.remove(intron_seq_file)
    os.rename(tmp_file, intron_seq_file)


def summarize(scores):
    summary = pystats.gmean(scores)  # all scores are in [0-1]
    # summary = np.cbrt(np.prod(scores))
    # summary = np.mean(scores)
    # summary = np.sqrt(sum([s ** 2 for s in scores]))

    return summary


def summary_score(region_thresholds, weights=None):
    if weights is not None:
        region_thresholds = {r: t * weights[r] for r, t in region_thresholds.items()}

    return summarize(list(region_thresholds.values()))


def progress_bar(current, total, done_char='#', fill_char=' ', length=20):
    fraction_done = current / total
    n_done = math.floor(fraction_done * length)
    n_remaining = length - n_done
    prog_bar = (done_char * n_done) + (fill_char * n_remaining)
    prog_bar = '[{}]'.format(prog_bar)

    return prog_bar


def make_matrices(
    introns, 
    u12_threshold,
    five_start_coord,
    three_start_coord,
    regions=['five', 'bp'], 
    min_seqs=5):
    seq_attrs = {
        'five': 'five_seq',
        'bp': 'bp_seq',
        'three': 'three_seq'
    }
    region_starts = {
        'five': five_start_coord,
        'three': three_start_coord,
        'bp': 0
    }
    u12_atac = re.compile('^AT[ACT]{3}')
    u12_gtag = re.compile('^GT[ACT]{3}')
    u2_gtag = re.compile('^GT')
    u2_gcag = re.compile('^GC')
    u12_motifs = {
        ('u12', 'atac'): {
            'pattern': u12_atac,
            'dnts': ('AT', 'AC')},
        ('u12', 'gtag'): {
            'pattern': u12_gtag,
            'dnts': ('GT', 'AG')}
    }
    u2_motifs = {
        ('u2', 'gtag'): {
            'pattern': u2_gtag,
            'dnts': ('GT', 'AG')},
        ('u2','gcag'): {
            'pattern': u2_gcag,
            'dnts': ('GC', 'AG')}
    }
    # define generators to feed sequences
    matrices = {}
    for r in regions:
        attr_name = seq_attrs[r]
        attr_start = region_starts[r]
        for subtype, info in u12_motifs.items():
            subtype_key = subtype + (r,)
            subtype_motif = info['pattern']
            dnts = info['dnts']
            u12_seqs = (
                getattr(i, attr_name) for i in introns if 
                i.svm_score > u12_threshold and 
                i.dnts == dnts and
                subtype_motif.match(i.five_display_seq)
            )
            matrix, n_seqs = matrix_from_seqs(
                u12_seqs, start_index=attr_start)
            if n_seqs > min_seqs:
                matrices[subtype_key] = matrix
        for subtype, info in u2_motifs.items():
            subtype_key = subtype + (r,)
            subtype_motif = info['pattern']
            dnts = info['dnts']
            u2_seqs = (
                getattr(i, attr_name) for i in introns if 
                i.type_id == 'u2' and 
                i.dnts == dnts and
                subtype_motif.match(i.five_display_seq)
            )
            matrix, n_seqs = matrix_from_seqs(
                u2_seqs, start_index=attr_start)
            if n_seqs > min_seqs:
                matrices[subtype_key] = matrix
    
    return matrices


def add_u2_matrix(
    introns, 
    args,
    min_u2_count=100
):
    matrices = args['MATRICES']
    five_score_coords = args['FIVE_SCORE_COORDS']
    three_score_coords = args['THREE_SCORE_COORDS']
    pseudocount = args['PSEUDOCOUNT']
    spcs = args['SPCS']
    simple_name = args['SIMPLE_NAME']
    U2_BPS_MATRIX_FILE = args['U2_BPS_MATRIX_FILE']

    u2_bp_key = ('u2', 'gtag', 'bp')
    u12_five_key = ('u12', 'gtag', 'five')

    if u2_bp_key not in matrices:
        # get the Nth percentile of 5' scores to cull putative U12s
        five_scores = []
        for i in introns:
            dnt_string = ''.join((e.lower() for e in i.dnts))
            tag = ('u12', dnt_string)
            score_info = multi_matrix_score(
                i, 
                matrices,
                five_score_coords,
                three_score_coords,
                pseudocount,
                regions=('five'),
                matrix_tags=['u12'])
            # try:
            best_tag = best_matrix(score_info, ['five'])[1]
            best_score = score_info[best_tag]['five']['score']
            five_scores.append(best_score)

        u2_five_percentile = 95
        u2_threshold = np.percentile(five_scores, u2_five_percentile)
        u2_bp_introns = (
            i for i, five in zip(introns, five_scores) 
            if five < u2_threshold
        )
        u12_bp_matrices = {
            k: v for k, v in matrices.items() if
            all(x in k for x in ['u12', 'bp', 'gtag'])}
        u2_bp_matrix, n_introns_used = build_u2_bp_matrix(
            u2_bp_introns,
            u12_bp_matrices,
            spcs,
            simple_name,
            dnt_list=[('GT', 'AG'), ('GC', 'AG')]
        )
        # use conserved u2 bp matrix as fallback if insufficient
        # u2s in provided data
        if n_introns_used < min_u2_count:
            u2_bp_matrix = load_external_matrix(U2_BPS_MATRIX_FILE)
            u2_bp_matrix = u2_bp_matrix[u2_bp_key]
            write_log(
                ('Insufficient U2 introns available to build BPS matrix; '
                 'using {} instead'), U2_BPS_MATRIX_FILE
            )
        else:
            write_log(
            '{} introns used to build U2 branch point matrix '
            '(5\'SS in bottom {}th percentile)',
            n_introns_used, u2_five_percentile)
        matrices[u2_bp_key] = add_pseudos(u2_bp_matrix, pseudo=pseudocount)
        matrices[('u2', 'gcag', 'bp')] = matrices[u2_bp_key]
    if ('u2', 'gcag', 'bp') not in matrices:
        matrices[('u2', 'gcag', 'bp')] = matrices[u2_bp_key]
        
    return matrices


def get_flipped(
    introns, 
    model, 
    scaler,
    scoring_region_labels,
    args
):
    SPCS = args['SPCS']
    N_PROC = args['N_PROC']
    MATRICES = args['MATRICES']
    SCORING_REGIONS = args['SCORING_REGIONS']
    FIVE_SCORE_COORDS = args['FIVE_SCORE_COORDS']
    THREE_SCORE_COORDS = args['THREE_SCORE_COORDS']
    PSEUDOCOUNT = args['PSEUDOCOUNT']
    THRESHOLD = args['THRESHOLD']
    SIMPLE_NAME = args['SIMPLE_NAME']
    flipped = {}
    swap_introns = copy.deepcopy(
        [
            i for i in introns if 
            (
                i.dnts == ('AT', 'AC') and not 
                i.five_display_seq.startswith('ATATC')
            )
            or i.noncanonical == True
        ]
    )
    # Only perform these calculation if we actually found applicable introns
    if swap_introns:
        raw_swaps = get_raw_scores(
            swap_introns, 
            MATRICES, 
            SCORING_REGIONS, 
            FIVE_SCORE_COORDS, 
            THREE_SCORE_COORDS, 
            PSEUDOCOUNT, 
            processes=N_PROC
        )
        scored_swaps = scale_scores(raw_swaps, scaler)
        mutant_swaps = [
            mutate(intron, FIVE_SCORE_COORDS, THREE_SCORE_COORDS)
            for intron in copy.deepcopy(swap_introns)]
        raw_mutants = get_raw_scores(
            mutant_swaps, 
            MATRICES, 
            SCORING_REGIONS, 
            FIVE_SCORE_COORDS, 
            THREE_SCORE_COORDS, 
            PSEUDOCOUNT, 
            processes=N_PROC
        )
        scored_mutants = scale_scores(raw_mutants, scaler)

        # scored_swaps = assign_svm_scores(
        #     scored_swaps, model, scoring_region_labels)

        # scored_mutants = assign_svm_scores(
        #     scored_mutants, model, scoring_region_labels)
        scored_swaps = parallel_svm_score(
            scored_swaps, model, scoring_region_labels, THRESHOLD, processes=N_PROC)
        scored_mutants = parallel_svm_score(
            scored_mutants, model, scoring_region_labels, THRESHOLD, processes=N_PROC)

        # for old, new in zip(
        #     sorted(scored_swaps, key=lambda x: x.unique_num), 
        #     sorted(scored_mutants, key=lambda x: x.unique_num)):
        for old, new in zip(scored_swaps, scored_mutants):
            if old.svm_score > THRESHOLD and new.svm_score <= THRESHOLD:
                name = old.get_name(SPCS, SIMPLE_NAME)
                flipped[old.unique_num] = new
                # flipped[name] = new

    return flipped
    

def demote(introns, flipped, spcs, simple_name, THRESHOLD):
    # For putative U12 introns whose scores don't survive dnt switching
    # demoted_swaps = []
    for i in introns:
        if (i.svm_score < THRESHOLD or (
            i.dnts != ('AT', 'AC') and i.noncanonical is False)):
            yield i
            continue
        if '[d]' in i.dynamic_tag:
            i.dynamic_tag.remove('[d]')
        flipped_i = flip_check(i, flipped)
        if flipped_i is not None:
            old_score = i.svm_score
            old_five = i.five_z_score
            old_bp = i.bp_z_score
            old_three = i.three_z_score
            flip_relative = flipped_i.relative_score
            flip_score = flipped_i.svm_score
            flip_five = flipped_i.five_z_score
            flip_bp = flipped_i.bp_z_score
            flip_three = flipped_i.three_z_score
            flip_dynamic_tag = flipped_i.dynamic_tag
            flip_type_id = flipped_i.type_id
            flip_matrices = flipped_i.matrices
            flip_info = [
                i.get_name(spcs, simple_name), 
                old_score, 
                old_five, 
                old_bp, 
                old_three,
                flip_score, 
                flip_five, 
                flip_bp, 
                flip_three,
                flipped_i.five_seq,
                flipped_i.bp_seq,
                flipped_i.three_seq
            ]
            i.demote_info = flip_info
            # intron = flipped_intron   # produces incorrect display seqs

            # Set original intron's score to new adjusted score
            demote_attrs = {
                'svm_score': flip_score,
                'relative_score': flip_relative,
                'five_z_score': flip_five,
                'bp_z_score': flip_bp,
                'three_z_score': flip_three,
                'type_id': flip_type_id,
                'matrices': flip_matrices
            }
            for da, v in demote_attrs.items():
                setattr(i, da, v)

            i.dynamic_tag.add('[d]')

        yield i


def set_attributes(objs, attr_list, attr_names):
    new_objs = []
    o_append = new_objs.append  # optimize function call
    for o, attrs in zip(objs, attr_list):
        for a, name in zip(attrs, attr_names):
            setattr(o, name, a)
        o_append(o)
    
    return new_objs


def scale_scores(
    introns,
    scaler,
    get_names=['five_raw_score', 'bp_raw_score', 'three_raw_score'], 
    set_names=['five_z_score', 'bp_z_score', 'three_z_score']):
    s_vect = get_score_vector(
        introns, get_names)
    s_vect = scaler.transform(s_vect)
    introns = set_attributes(introns, s_vect, set_names)

    return introns


def apply_scores(
    ref_set, 
    exp_set,
    matrices,
    scoring_region_labels,
    args,
    log=True
):
    SCORING_REGIONS = args['SCORING_REGIONS']
    N_PROC = args['N_PROC']
    OPTIMIZE_N = args['OPTIMIZE_N']
    SVM_ITER = args['SVM_ITER']
    CV_JOBS = args['CV_JOBS']
    HYPER_C = args['HYPER_C']
    SPECIES = args['SPECIES']
    FIG_DPI = args['FIG_DPI']
    THRESHOLD = args['THRESHOLD']
    SPCS = args['SPCS']
    SIMPLE_NAME = args['SIMPLE_NAME']
    FIVE_SCORE_COORDS = args['FIVE_SCORE_COORDS']
    THREE_SCORE_COORDS = args['THREE_SCORE_COORDS']
    PSEUDOCOUNT = args['PSEUDOCOUNT']

    # Get the raw log ratio scores of each scoring region in each intron
    # pool = Pool(processes=4)
    # raw_introns = pool.starmap(
    #     get_raw_scores, zip(exp_set, repeat(matrices), repeat(SCORING_REGIONS)))
    raw_introns = get_raw_scores(
        exp_set, 
        matrices, 
        SCORING_REGIONS,
        FIVE_SCORE_COORDS,
        THREE_SCORE_COORDS,
        PSEUDOCOUNT, 
        processes=N_PROC)
    # raw_introns = get_raw_scores(exp_set, matrices, SCORING_REGIONS)
    # raw_introns = list(raw_introns)

    # Same for the reference sequences
    # raw_refs = pool.starmap(
    #     get_raw_scores, zip(ref_set, repeat(matrices), repeat(SCORING_REGIONS)))
    raw_refs = get_raw_scores(
        ref_set, 
        matrices, 
        SCORING_REGIONS,
        FIVE_SCORE_COORDS,
        THREE_SCORE_COORDS,
        PSEUDOCOUNT, 
        processes=N_PROC)
    # raw_refs = get_raw_scores(ref_set, matrices, SCORING_REGIONS)
    # raw_refs = list(raw_refs)

    raw_score_names = ['five_raw_score', 'bp_raw_score', 'three_raw_score']

    # scale based on the training set only, to avoid
    # test set data leak
    scale_vector = get_score_vector(
        raw_refs, score_names=raw_score_names)
    
    # make a scaler to adjust raw scores --> z-scores
    score_scaler = preprocessing.StandardScaler().fit(scale_vector)

    scored_refs = scale_scores(raw_refs, score_scaler)
    ref_u12s = [i for i in scored_refs if i.type_id == 'u12']
    ref_u2s = [i for i in scored_refs if i.type_id == 'u2']
    if log is True:
        write_log(
            'Raw scores calculated for {} U2 and {} U12 reference introns',
            len(ref_u2s), len(ref_u12s)
        )
    scored_introns = scale_scores(raw_introns, score_scaler)
    if log is True:
        write_log(
            'Raw scores calculated for {} experimental introns',
            len(scored_introns)
        )
    
    # make score vectors for training data
    ref_u12_vector = get_score_vector(
        ref_u12s, score_names=scoring_region_labels)
    ref_u2_vector = get_score_vector(
        ref_u2s, score_names=scoring_region_labels)

    # NOTE: in this application with a "soft-margin" (C > 0) SVM,
    # redundant data arguably should be included, if present.
    # In this case, that distinct introns with identical 
    # sequence in their scoring regions, so here we should not filter
    # out redundant score vectors so long as we've checked that our 
    # input full sequences are unique.
    # sources: 
    # https://martin-thoma.com/svm-with-sklearn/, 
    # https://stats.stackexchange.com/questions/45045/regarding-redundant-training-data-in-building-svm-based-classifier,
    # https://stats.stackexchange.com/questions/23143/remove-duplicates-from-training-set-for-classification)
    # https://stats.stackexchange.com/questions/310005/should-we-remove-duplicates-when-training-a-svm
    # ref_u12_vector, unique_u12_index = np.unique(
        # ref_u12_vector, axis=0, return_index=True)
    # ref_u2_vector, unique_u2_index = np.unique(
        # ref_u2_vector, axis=0, return_index=True)

    # unique_ref_u12s = [
        # e for i, e in enumerate(ref_u12s) if i in unique_u12_index]
    # unique_ref_u2s = [
        # e for i, e in enumerate(ref_u2s) if i in unique_u2_index]

    # scored_refs = unique_ref_u2s + unique_ref_u12s

    if log is True:
        write_log(
            'Training set score vectors constructed: {} U2, {} U12',
            len(ref_u2_vector), len(ref_u12_vector))
        write_log('Training SVM using reference data')

    svm_start = time.time()

    ###!!!!! reinstate
    model = svm.SVC(probability=True, kernel='linear', cache_size=1000)
    # model = svm.SVC(probability=True, kernel='rbf', gamma='scale', cache_size=1000)
    model, model_performance = optimize_svm(
        model,
        ref_u12_vector,
        ref_u2_vector,
        n_optimize=OPTIMIZE_N,
        iterations=SVM_ITER,
        cv_jobs=CV_JOBS,
        hyper_C=HYPER_C
    )

    svm_train_time = get_runtime(svm_start)

    # PERFORMANCE METRICS
    if log is True:
        auc = round(np.mean(model_performance['pr_auc']), 6)
        f1 = round(np.mean(model_performance['f1']), 6)
        all_reports = model_performance['classification_report']
        report = '\n'.join(model_performance['classification_report'])

        write_log(
            'Average classifier performance on training data:'
            '\n\tF1\t{}\n\tP-R AUC\t{}',
            f1, auc
        )
        if not SVM_ITER:
            write_log(
                'Classifier performance details:\n{}', report, wrap_chars=None)
        else:
            stats_outfile = '{}.classifier_performance.iic'.format(SPECIES)
            with open(stats_outfile, 'w') as f:
                f.write(report)
            write_log(
                'Per-subsample classifier performance info written to {}', 
                stats_outfile)

        pr_curves = model_performance['pr_curve']
        for curve in pr_curves:
            precision, recall, _ = curve
            plt.plot(recall, precision)
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('{} precision-recall AUC: {:.3f}'.format(SPECIES, auc))
        plt.tight_layout()
        plt.savefig('{}.AUC.iic.png'.format(SPECIES), dpi=FIG_DPI)

    #/ PERFORMANCE METRICS

    if len(scoring_region_labels) > 1:
        ref_scatter(
            ref_u2_vector, 
            ref_u12_vector, 
            model[0], 
            SCORING_REGIONS, 
            SPECIES
        )

        density_hexplot(
            np.concatenate((ref_u2_vector, ref_u12_vector)), 
            title='{}.ref_hex'.format(SPECIES),
            xlab='5\' z-score', ylab='BPS z-score', fig_dpi=FIG_DPI)

    # run trained SVM on experimental data
    # scored_introns = assign_svm_scores(
    #     scored_introns, model, scoring_regions, weights=model_performance['f1'])

    scored_introns = parallel_svm_score(
        scored_introns, model, scoring_region_labels, THRESHOLD,
        weights=model_performance['f1'], processes=N_PROC)

    flipped = get_flipped(
        scored_introns, model, score_scaler, scoring_region_labels, args)
    if log is True:
        write_log(
            '{} putative U12 scores were not robust to boundary switching',
            len(flipped.keys()))

    finalized_introns = []
    u12_count = 0
    atac_count = 0
    demoted_swaps = []

    for i in demote(scored_introns, flipped, SPCS, SIMPLE_NAME, THRESHOLD):
        if '[d]' in i.dynamic_tag:
            demoted_swaps.append(i.demote_info)
        if i.svm_score > THRESHOLD:
            u12_count += 1
            if i.dnts == ('AT', 'AC'):
                atac_count += 1
        finalized_introns.append(i)
    
    return finalized_introns, model, u12_count, atac_count, demoted_swaps


def recursive_scoring(
    finalized_introns,
    refs,
    model,
    args,
    raw_score_names,
    z_score_names
):

    FIVE_SCORE_COORDS = args['FIVE_SCORE_COORDS']
    THREE_SCORE_COORDS = args['THREE_SCORE_COORDS']
    matrices = args['MATRICES']
    threshold = args['THRESHOLD']
    SCORING_REGIONS = args['SCORING_REGIONS']
    PSEUDOCOUNT = args['PSEUDOCOUNT']
    U12_COUNT = args['U12_COUNT']
    N_PROC = args['N_PROC']
    SCORING_REGION_LABELS = args['SCORING_REGION_LABELS']
    # use introns from first round to create new matrices
    # for second round
    write_log('Updating scoring matrices using empirical data')
    new_matrices = make_matrices(
        finalized_introns, 
        threshold, 
        FIVE_SCORE_COORDS[0], 
        THREE_SCORE_COORDS[0],
        min_seqs=5)

    new_matrices = add_pseudos(new_matrices, pseudo=PSEUDOCOUNT)
    if U12_COUNT < 100:
        mod_type = 'Averaging'
        matrices = average_matrices(matrices, new_matrices)
    # replace old matrices with new ones
    else:
        mod_type = 'Replacing'
        matrices.update(new_matrices)
    write_log('{} matrices using empirically-derived data'.format(mod_type))

    # re-run previously-trained SVM models on introns scored with new 
    # matrices to filter out any introns whose scored changed dramatically 
    # before picking a species-specific reference set
    raw_introns = get_raw_scores(
        finalized_introns, 
        matrices, 
        SCORING_REGIONS,
        FIVE_SCORE_COORDS,
        THREE_SCORE_COORDS,
        PSEUDOCOUNT,
        processes=N_PROC
    )
    raw_refs = get_raw_scores(
        refs, 
        matrices, 
        SCORING_REGIONS,
        FIVE_SCORE_COORDS,
        THREE_SCORE_COORDS,
        PSEUDOCOUNT,
        processes=N_PROC
    )
    scale_vector = get_score_vector(
        raw_refs, score_names=raw_score_names)
    recursive_scaler = preprocessing.RobustScaler().fit(scale_vector)
    scored_introns = scale_scores(raw_introns, recursive_scaler)
    # scored_introns = assign_svm_scores(
    #     scored_introns, model, scoring_region_labels)
    scored_introns = parallel_svm_score(
        scored_introns, 
        model, 
        SCORING_REGION_LABELS, 
        threshold, 
        processes=N_PROC
    )

    # filter reference introns to non-redundant set
    # before training SVM
    ref_u12_threshold = 95
    ref_u2_threshold = 5
    unambiguous = [
        copy.deepcopy(i) for i in scored_introns if 
        i.svm_score > ref_u12_threshold or i.svm_score < ref_u2_threshold]
    recursive_refs = get_attributes(
        unambiguous, ['five_seq', 'bp_region_seq', 'three_seq'])

    scale_vector = get_score_vector(
        recursive_refs, score_names=raw_score_names)

    recursive_scaler = preprocessing.RobustScaler().fit(scale_vector)

    scored_introns = scale_scores(raw_introns, recursive_scaler)
    scored_refs = scale_scores(recursive_refs, recursive_scaler)

    recursive_ref_u12s = [
        i for i in scored_refs if i.svm_score > ref_u12_threshold]
    recursive_ref_u2s = [
        i for i in scored_refs if i.type_id == 'u2']

    write_log(
        'Empirically-derived training data: {} U2, {} U12', 
        len(recursive_ref_u2s), len(recursive_ref_u12s))
    # ensure that the training set is sufficiently large
    # adding the reference introns with scores derived from the
    # initial matrices ensures that changes to the matrices do
    MIN_REF_U2 = 5000
    MIN_REF_U12 = 50

    ref_u2_delta = MIN_REF_U2 - len(recursive_ref_u2s)
    ref_u12_delta = MIN_REF_U12 - len(recursive_ref_u12s)
    random.seed(42)
    if ref_u2_delta > 0 or ref_u12_delta > 0:
        raw_default_refs = get_raw_scores(
            refs, 
            matrices, 
            SCORING_REGIONS,
            FIVE_SCORE_COORDS,
            THREE_SCORE_COORDS,
            PSEUDOCOUNT,
            processes=N_PROC
        )
        scored_default_refs = scale_scores(raw_default_refs, recursive_scaler)
        # working_default_refs = assign_svm_scores(
        #     scored_default_refs, model, scoring_region_labels)
        working_default_refs = parallel_svm_score(
            scored_default_refs, 
            model, 
            SCORING_REGION_LABELS, 
            threshold, 
            processes=N_PROC
        )
        scored_ref_u2s = [
            i for i in working_default_refs if i.svm_score < ref_u2_threshold]
        scored_ref_u12s = [
            i for i in working_default_refs if i.svm_score > ref_u12_threshold]
    if ref_u2_delta > 0:
        write_log('[!] Adding {} reference U2s to meet minimum', ref_u2_delta)
        recursive_ref_u2s += random.sample(scored_ref_u2s, ref_u2_delta)
    if ref_u12_delta > 0:
        write_log('[!] Adding {} reference U12s to meet minimum', ref_u12_delta)
        recursive_ref_u12s += random.sample(scored_ref_u12s, ref_u12_delta)

    recursive_refs = recursive_ref_u2s + recursive_ref_u12s

    return recursive_refs, scored_introns, matrices


def add_custom_matrices(custom_matrix_files, default_matrices):
    '''
    Integrates user-supplied matrices with the default matrices,
    overwriting those in the default set as needed.
    
    '''
    tmp_custom_matrices = {}
    all_custom_keys = set()
    for m in custom_matrix_files:
        custom_matrix = load_external_matrix(m)
        custom_keys = ['-'.join(key) for key in custom_matrix.keys()]
        for k in custom_keys:
            all_custom_keys.add(k)
        tmp_custom_matrices.update(custom_matrix)
    # in special case where there are diff. versions of a matrix,
    # e.g. u12-gtag-bp-vA10
    versioned_keys = [k for k in tmp_custom_matrices if len(k) > 3]
    for k in versioned_keys:
        category = tuple(k[:3])
        matches = [mk for mk in default_matrices.keys() if all(e in mk for e in category)]
        for m in matches:
            del default_matrices[m]
            # del MATRICES[m]
    default_matrices.update(tmp_custom_matrices)
    write_log('Custom matrices used: {}', ','.join(sorted(all_custom_keys)))

    return default_matrices


def introns_from_seqs(
    source_file, 
    five_coords, 
    three_coords, 
    bp_region_coords, 
    allow_noncanonical):
    # source_file = args.sequence_file
    all_introns = introns_from_flatfile(
        source_file,
        five_coords,
        three_coords,
        bp_region_coords,
        allow_noncanonical,
        hashgen=False,
        allow_overlap=False)
    all_introns = list(all_introns)
    final_introns = [i for i in all_introns if not i.omitted]
    total_count = len(final_introns)
    omit_count = len(all_introns) - total_count
    if omit_count > 0:
        write_log('{} introns omitted', omit_count)

    if total_count == 0:
        write_log(
            '[!] ERROR: No intron sequences found. Exiting now.'
        )
        sys.exit(1)

    return final_introns, total_count


def introns_from_annotation(annotation, feature):
    # source_file = ANNOTATION
    feats_to_try = (feature, ('cds',), ('exon',))
    feature_found = False
    # Check to make sure annotation file has desired feature type
    for feat in feats_to_try:
        if not has_feature(annotation, feat):
            write_log(
                '[!] No {} features in annotation.', ','.join(feat))
        else:
            feature_found = True
            feature = feat
            break
    if feature_found is False:
        write_log(
            '[!] ERROR: No intron-defining features found in annotation. '
            'Exiting now.')       
        sys.exit(1)

    write_log('Using {} features to define introns', ','.join(feature))

    # set a temporary feature list here to allow all possible
    # introns to be gathered initially, then pruned in
    # collect_introns() if needed
    working_feature = ('cds', 'exon')   ###!!!

    # original functionality
    # working_feature = FEATURE

    # Pull all annotation entries into a hierarchy of objects
    top_level_annots = annotation_hierarchy(
        annotation, working_feature)

    # Make a dictionary index for all intron-defining objects in the
    # annotated heirarchy. This index will allow easy access to parent
    # objects when mucking about with the introns

    # Add transcripts here to allow longest-isoform tagging
    feat_list = working_feature + ('transcript',)
    flat_annots = flatten(top_level_annots, feat_list=feat_list)

    # Make intron object dictionary from the collected top-level objects,
    # including whatever duplicates might exist in the annotation
    all_introns, total_count = collect_introns(
        top_level_annots, feature, flat_annots)

    if total_count == 0:
        write_log(
            '[!] ERROR: No intron sequences found. '
            'Annotation <-> genome ID mismatch, perhaps? Exiting now.')
        sys.exit(1)
    
    return all_introns, total_count, flat_annots


def introns_from_bedfile(genome, bed):
    """
    chr start   stop    name    score   strand
    
    """
    bed_introns = []
    seen_coords = set()
    intron_index = defaultdict(list)
    intron_count = 0
    with flex_open(bed) as b:
        for line_number, l in enumerate(b, start=1):
            if l.startswith('#'):
                continue
            try:
                loc, start, stop, name, _, strand = l.strip().split('\t')[:6]
            except:
                continue
            try:
                # ensure start is lowest coord
                start, stop = sorted(map(int, [start, stop]))
            except:
                continue
            if strand not in ('+', '-'):
                continue
            if not name or name in ('-.*'):  # Null indicators
                name = 'i_{}'.format(line_number)
            start += 1  # BED files are 0-indexed
            intron = Intron(
                name=name, 
                start=start, 
                stop=stop, 
                strand=strand, 
                region=loc
            )
            intron_index[loc].append(intron)
            intron_count += 1
    
    return intron_index, intron_count


def output_format(
    intron, 
    out_type, 
    spcs, 
    simple_name=False, 
    null='.', 
    separator='\t'
):
    attrib_index = {
        'BED': [
            'region',
            'start',
            'stop',
            # 'get_label',
            'svm_score',
            'strand'
        ],
        'META': [
            # 'get_name',
            'motif_string',
            'bps_context',
            'length',
            'parent',
            'grandparent',
            'index',
            'family_size',
            'fractional_position',
            'phase',
            'type_id',
            'defined_by'
        ],
        'SCORE': [
            # 'get_name',
            'relative_score',
            'svm_score',
            'score_distance',
            'five_seq',
            'five_raw_score',
            'five_z_score',
            'bp_seq',
            'bp_raw_score',
            'bp_z_score',
            'three_seq',
            'three_raw_score',
            'three_z_score'
        ],
        'SEQ': [
            # 'get_name',
            'svm_score',
            'upstream_flank',
            'seq',
            'downstream_flank'
        ]
    }
    to_get = attrib_index[out_type]
    # attribs = write_format(intron, to_get, separator=None, null=null)
    attribs = [getattr(intron, a) for a in to_get]
    # run any methods that were returned along with other static types
    attribs = [
        a() if type(a)==types.MethodType else a for a in attribs]
    attribs = [null if a is None else a for a in attribs]
    name = intron.get_name(spcs, simple_name)
    label = intron.get_label(spcs, simple_name)
    if out_type in ('META', 'SEQ', 'SCORE'):
        attribs.insert(0, name)
        if out_type == 'META':
            try:
                rounded_score = round(intron.relative_score, 4)
            except:
                rounded_score = null
            dnts = '-'.join(intron.dnts)
            attribs.insert(1, rounded_score)
            attribs.insert(2, dnts)
    elif out_type == 'BED':
        attribs.insert(3, label)
    if separator is not None:
        attribs = separator.join([str(a) for a in attribs])
    
    return attribs

def get_args(argv, arg_parser):
    format_info_message = (
    """
    [ data formatting information ]

    # Matrix files ################################################################

    /// Description ///

    The matrix file describes the frequency of each nucleotide at each position
    in scoring region for different intron types, in FASTA format with specific
    header and first line of each record requirements (see below).

    /// Naming ///

    The matrix file must be named 'scoring_matrices.fasta.iic', and placed in a
    directory called 'data' located in the main intronIC directory.

    /// Format ///

    Each matrix in the file consists of a FASTA header line, indicating
    which type of intron the matrix represents ('u2', 'u12'), which region
    it describes ('five', 'bp', 'three') and which intron subtype it applies to
    ('gtag', 'atac', etc). Optionally, this header may also include a special
    keyword phrase 'start={integer}', where {integer} is the position of
    the first matrix entry relative to the corresponding splice site. If the
    5' matrix begins with three exonic positions upstream of the beginning of 
    the intron, for example, then the header for that matrix would include
    'start=-3'. If this keyword is omitted, the matrix is assumed to start at 
    position 0, the first base of corresponding scoring region. Because the branch 
    point motif occurs at variable locations, branch point matrices do not use this 
    convention (e.g. all BPS matrices start at 0).

    The line after the header must be the whitespace-separated order of the bases 
    as they appear in the subsequent lines of the matrix. The rest of the lines
    under the header are tab-separated values representing the frequency of each 
    base (columns) at each position in the target sequence (rows).

    Example formatting ('scoring_matrices.fasta.iic'):

    >u12_gtag_five  start=-3
    A   C   G   T
    0.293478260869565	0.239130434782609	0.239130434782609	0.228260869565217
    0.271739130434783	0.326086956521739	0.152173913043478	0.25
    0.217391304347826	0.184782608695652	0.0543478260869565	0.543478260869565
    0	0	1	0
    0	0	0	1
    [...]
    >u12_gtag_bp
    A   C   G   T
    0.12	0.15	0.21	0.52
    0.12	0.17	0.18	0.53
    0.13	0.2 0.12	0.55
    [...]

    # Reference introns ###########################################################

    /// Description ///

    The reference intron set is a collection of introns with evidence supporting 
    classification into one of the two types. Each intron will be scored against 
    the same matrices used for the experimental dataset and fed into the classifier
    as training data.

    /// Naming ///

    The files needs to be named '[u2, u12]_reference_set.introns.iic[.gz]', and
    be located in the data directory, unless specified via the --r[2, 12]
    command line argument.

    /// Format ///

    Each line in the file contains the
    following columns:

    1. Identifier
    2. 5' exonic sequence (>= the length used for scoring, if any)
    3. Intronic sequence
    4. 3' exonic sequence (>= the length used for scoring, if any)

    """)
    if len(argv) == 1:
        arg_parser.print_usage()
        sys.exit(0)

    if '--format_info' in argv:  # can't use argparse unless other args
        sys.exit(format_info_message)

    args = arg_parser.parse_args()

    if args.format_info:
        sys.exit(format_info_message)

    if not (args.genome and (args.annotation or args.bed)) or args.sequence_file:
        arg_parser.error(
            '{}: error: must be run with either a genome and annotation/BED, '
            'or intron sequences file'.format(os.path.basename(sys.argv[0]))
        )
    if args.sequence_file and (args.genome or args.annotation or args.bed):
        arg_parser.error(
            'Must specify either direct sequence input (via -q) or a '
            'genome/annotation combination'
        )
    if args.bed and args.exons_as_flanks:
        arg_parser.error(
            'Cannot use --exons_as_flanks in combination with BED input'
        )
    if args.bed and args.annotation:
        arg_parser.error(
            'Use either -a or -b, but not both'
        )
    
    return args

def get_custom_args(args, argv):
    custom_args = {}
    matrix_filename = "scoring_matrices.fasta.iic"
    reference_u12s_filename = "u12_reference.introns.iic"
    reference_u2s_filename = "u2_reference.introns.iic"
    backup_u2_bps_filename = 'u2.conserved_empirical_bp_matrix.iic'

    # Get script directory and external file paths
    DATA_DIR = pkg_resources.resource_filename(__name__, 'data/')
    DATA_DIR = os.path.realpath(DATA_DIR)

    # HOME = os.path.dirname(os.path.realpath(sys.argv[0]))
    # custom_args['HOME'] = HOME
    # DATA_DIR = os.path.join(HOME, "data")

    custom_args['DATA_DIR'] = DATA_DIR
    RUN_DIR = os.getcwd()
    custom_args['RUN_DIR'] = os.getcwd()
    MATRIX_FILE = os.path.join(DATA_DIR, matrix_filename)
    custom_args['MATRIX_FILE'] = MATRIX_FILE
    U2_BPS_MATRIX_FILE = os.path.join(DATA_DIR, backup_u2_bps_filename)
    custom_args['U2_BPS_MATRIX_FILE'] = U2_BPS_MATRIX_FILE
    # This is inelegant, but to do otherwise would require reorg of args
    # and constants (defaulting the arg would require knowing the data dir
    # and matrix filename, for example)
    if args.reference_u12s:
        REFERENCE_U12_FILE = args.reference_u12s
    else:
        ref_u12_fns = [
            os.path.join(DATA_DIR, r) for r in
            [reference_u12s_filename, reference_u12s_filename + '.gz']]
        REFERENCE_U12_FILE = next(r for r in ref_u12_fns if os.path.isfile(r))
    if args.reference_u2s:
        REFERENCE_U2_FILE = args.reference_u2s
    else:
        ref_u2_fns = [
            os.path.join(DATA_DIR, r) for r in
            [reference_u2s_filename, reference_u2s_filename + '.gz']]
        REFERENCE_U2_FILE = next(r for r in ref_u2_fns if os.path.isfile(r))

    custom_args['REFERENCE_U12_FILE'] = REFERENCE_U12_FILE
    custom_args['REFERENCE_U2_FILE'] = REFERENCE_U2_FILE

    custom_args['SVM_ITER'] = args.n_subsample
    custom_args['N_PROC'] = args.processes
    custom_args['CV_JOBS'] = args.cv_processes
    # set cross-validation parallel processing to same as main
    # parallel arg unless explicitly specified
    if not custom_args['CV_JOBS']:
        custom_args['CV_JOBS'] = custom_args['N_PROC']
    custom_args['RECURSIVE'] = args.recursive
    # scoring region coordinates are relative to the 5' and 3' ends of the intron
    custom_args['FIVE_SCORE_COORDS'] = tuple(args.five_score_coords)
    custom_args['THREE_SCORE_COORDS'] = tuple(args.three_score_coords)
    custom_args['BP_REGION_COORDS'] = tuple(args.branch_point_coords)
    custom_args['OPTIMIZE_N'] = 5  # rounds of SVM optimization
    custom_args['INTRON_FLANK_SIZE'] = 200  # for output file
    custom_args['MIN_INTRON_LENGTH'] = args.min_intron_len
    custom_args['NUM_NC_USER'] = 5 # number of nc splice sites to print to screen
    custom_args['NUM_NC_LOG'] = 25  # number of nc splice sites to print to log
    custom_args['PSEUDOCOUNT'] = args.pseudocount # value to add to matrix values to avoid div 0 errors
    custom_args['EXONS_AS_FLANKS'] = args.exons_as_flanks
    custom_args['FIG_DPI'] = 300

    custom_args['GENOME'] = args.genome
    custom_args['ANNOTATION'] = args.annotation
    custom_args['SEQUENCE_INPUT'] = args.sequence_file
    custom_args['SPECIES_FULL'] = args.species_name  # used in log file header
    custom_args['SPECIES'] = args.species_name  # used to name files
    if not args.no_abbreviate:
        SPCS = abbreviate(custom_args['SPECIES'])  # used within files
        SPECIES_NAME_INFO = '{} ({})'.format(custom_args['SPECIES_FULL'], SPCS)
    else:
        SPCS = args.species_name
        SPECIES_NAME_INFO = args.species_name
    custom_args['SPCS'] = SPCS
    custom_args['SPECIES_NAME_INFO'] = SPECIES_NAME_INFO
    if not args.feature:  # to define introns
        FEATURE = ('cds', 'exon')
    else:
        FEATURE = (args.feature.lower(),)
    custom_args['FEATURE'] = FEATURE
    if args.no_nc:
        ALLOW_NONCANONICAL = False
    else:
        ALLOW_NONCANONICAL = True
    custom_args['ALLOW_NONCANONICAL'] = ALLOW_NONCANONICAL
    custom_args['ALLOW_OVERLAP'] = not args.no_intron_overlap
    custom_args['LONGEST_ONLY'] = not args.allow_multiple_isoforms
    custom_args['THRESHOLD'] = args.threshold
    custom_args['INCLUDE_DUPES'] = args.include_duplicates
    custom_args['ONLY_SEQS'] = args.sequences_only
    custom_args['NO_SEQS'] = args.no_sequence_output
    custom_args['CUSTOM_MATRICES'] = args.pwms
    custom_args['SIMPLE_NAME'] = args.uninformative_naming
    custom_args['SCORING_REGIONS'] = args.scoring_regions
    custom_args['MATRIX_SCORE_INFO'] = args.matrix_score_info
    if args.hyperparameter_C:
        HYPER_C = float(args.hyperparameter_C)
    else:
        HYPER_C = None
    custom_args['HYPER_C'] = HYPER_C

    custom_args['WANT_PLOT'] = not args.no_plot
    custom_args['PLOT'] = not args.no_plot
    if custom_args['PLOT'] and check_plot() is False:
        custom_args['PLOT'] = False

    # Change file-naming variable if specified
    if args.abbreviate_filenames:
        SPECIES_FN = abbreviate(custom_args['SPECIES'])
    else:
        SPECIES_FN = custom_args['SPECIES']

    # Output filenames
    custom_args['FN_SEQS'] = '{}.introns.iic'.format(SPECIES_FN)
    custom_args['FN_RANKINGS'] = '{}.rankings.iic'.format(SPECIES_FN)
    custom_args['FN_BED'] = '{}.bed.iic'.format(SPECIES_FN)
    custom_args['FN_DUPE_MAP'] = '{}.dupe_map.iic'.format(SPECIES_FN)
    custom_args['FN_OVERLAP_MAP'] = '{}.overlap.iic'.format(SPECIES_FN)
    custom_args['FN_MATRICES'] = '{}.pwms.iic'.format(SPECIES_FN)
    custom_args['FN_META'] = '{}.meta.iic'.format(SPECIES_FN)
    custom_args['FN_SCORE'] = '{}.score_info.iic'.format(SPECIES_FN)
    custom_args['FN_SWAP'] = '{}.demoted.iic'.format(SPECIES_FN)
    custom_args['FN_ANNOT'] = '{}.annotation.iic'.format(SPECIES_FN)
    custom_args['FN_LOG'] = '{}.log.iic'.format(SPECIES_FN)

    # convert arg index to a namedtuple for prettier referencing
    # custom_args = namedtuple(
        # 'CustomArgs', custom_args.keys())(*custom_args.values())
    
    return custom_args

def add_pwm_args(args):
    MATRICES = load_external_matrix(args['MATRIX_FILE'])

    if args['CUSTOM_MATRICES']:
        MATRICES = add_custom_matrices(args['CUSTOM_MATRICES'], MATRICES)

    # Pseudocounts added to avoid division by 0 during scoring
    MATRICES = add_pseudos(MATRICES, pseudo=args['PSEUDOCOUNT'])

    # Determine length of 5' sequence region from supplied matrices
    # and use if shorter than >FIVE_LENGTH<
    five_matrix_length = matrix_length(MATRICES[('u12', 'gtag', 'five')])
    three_matrix_length = matrix_length(MATRICES[('u12', 'gtag', 'three')])

    FIVE_SCORE_COORDS = args['FIVE_SCORE_COORDS']
    THREE_SCORE_COORDS = args['THREE_SCORE_COORDS']

    ARG_FIVE_LENGTH = abs(FIVE_SCORE_COORDS[0] - FIVE_SCORE_COORDS[1])
    ARG_THREE_LENGTH = abs(THREE_SCORE_COORDS[0] - THREE_SCORE_COORDS[1])

    FIVE_LENGTH = min(ARG_FIVE_LENGTH, five_matrix_length)
    THREE_LENGTH = min(ARG_THREE_LENGTH, three_matrix_length)

    # Notify user if desired length is incompatible with matrix
    for l, arg, tag in zip(
        [FIVE_LENGTH, THREE_LENGTH],
        [ARG_FIVE_LENGTH, ARG_THREE_LENGTH],
        ['5\'', '3\'']
        ):
        if l < arg:
            write_log(
                '[!] Length of {} region limited to {} by scoring matrices',
                tag, l
            )
    args['FIVE_LENGTH'] = FIVE_LENGTH
    args['THREE_LENGTH'] = THREE_LENGTH

    # check min intron length against total size
    # of scoring matrices and adjust up if necessary
    bp_key = next(k for k in MATRICES.keys() if 'u12' in k and 'bp' in k)
    args['BP_MATRIX_LENGTH'] = matrix_length(MATRICES[bp_key])
    # BP_MATRIX_LENGTH = matrix_length(MATRICES[('u12', 'gtag', 'bp', 'vA9')])

    # TODO adjust the bp margin to account for change in interval openness
    # This seems to have been done...the empirical data min/max BPS
    # sequence locations match the coords listed in BP_REGION_COORDS
    bp_margin = abs(args['BP_REGION_COORDS'][1])

    # calculate length required by matrices at either end of intron
    intronic_five = len([e for e in range(*FIVE_SCORE_COORDS) if e >= 0])
    intronic_three = len([e for e in range(*THREE_SCORE_COORDS) if e < 0])
    intronic_three = max(
        (args['BP_MATRIX_LENGTH'] + bp_margin), intronic_three)

    # scoring_region_size = FIVE_LENGTH + BP_MATRIX_LENGTH + bp_margin
    scoring_region_size = intronic_five + intronic_three

    MIN_INTRON_LENGTH = args['MIN_INTRON_LENGTH']

    if MIN_INTRON_LENGTH < scoring_region_size and not args['ONLY_SEQS']:
        MIN_INTRON_LENGTH = scoring_region_size
        write_log(
            '[!] Minimum intron length set to {} due to size of scoring regions',
            MIN_INTRON_LENGTH
        )
    args['MIN_INTRON_LENGTH'] = MIN_INTRON_LENGTH
    args['MATRICES'] = MATRICES

    return args

def introns_from_args(args, ext_args):
    # return a trio of values regardless of method
    introns = []
    total_count = None
    flat_annots = {}
    if ext_args.sequence_file:
        source_file = ext_args.sequence_file
        introns, total_count = introns_from_seqs(
            source_file,
            args['FIVE_SCORE_COORDS'],
            args['THREE_SCORE_COORDS'],
            args['BP_REGION_COORDS'],
            args['ALLOW_NONCANONICAL']
        )
    elif ext_args.bed:
        source_file = ext_args.bed
        introns, total_count = introns_from_bedfile(args['GENOME'], source_file)
    else:
        source_file = args['ANNOTATION']
        introns, total_count, flat_annots = introns_from_annotation(
            source_file, args['FEATURE'])
    
    write_log('{} introns found in {}', total_count, source_file)
    
    return introns, total_count, flat_annots


def log_omitted(stats_omitted_counts, min_intron_length):
    write_log(
        '{} introns omitted from scoring based on the following criteria:',
        sum(stats_omitted_counts.values()))
    write_log(
        '* short (<{} nt): {}', min_intron_length, stats_omitted_counts['s'],
        wrap_chars=None)
    write_log(
        '* ambiguous nucleotides in scoring regions: {}',
        stats_omitted_counts['a'],
        wrap_chars=None)
    write_log(
        '* non-canonical boundaries: {}', stats_omitted_counts['n'],
        wrap_chars=None)
    write_log(
        '* overlapping coordinates: {}', stats_omitted_counts['v'],
        wrap_chars=None)
    write_log(
        '* not in longest isoform: {}', stats_omitted_counts['i'],
        wrap_chars=None)


def filter_introns_write_files(
    all_introns, 
    flat_annots, 
    total_count, 
    args,
    ext_args
):
    # Populate introns with sub-sequences and write full seqs to file
    # Also, tag introns with >omitted< if meet omission criteria
    # Omitted introns are written to same file as all other introns
    stats_omitted_counts = Counter()
    stats_corrected_count = 0
    stats_nc_types = Counter()  # omitted for non-canon bounds
    final_introns = []
    omitted_introns = []
    corrected_duplicates = []
    duplicate_count = 0
    # Only make seqs file if not argument preventing it
    opened_files = []
    if not args['NO_SEQS']:
        seq_file = open(args['FN_SEQS'], 'w')
        opened_files.append(seq_file)
    # Write the omitted introns to bed file first
    meta_file = open(args['FN_META'], 'w')
    opened_files.append(meta_file)
    bed_file = open(args['FN_BED'], 'w')
    opened_files.append(bed_file)
    # build an index of duplicate introns keyed by (region, start, stop) of
    # the intron which superceded them to allow score propagation
    dupe_intron_index = defaultdict(set)
    overlap_index = defaultdict(set)
    # map unique_num attributes to final intron names for use in
    # dupe_map output
    intron_name_index = {}

    # these two dictionaries enable communication between get_sub_seqs()
    # and add_tags()
    intron_index = defaultdict(lambda: defaultdict(dict))
    longest_isoforms = {}

    GENOME = args['GENOME']
    INTRON_FLANK_SIZE = args['INTRON_FLANK_SIZE']
    FIVE_SCORE_COORDS = args['FIVE_SCORE_COORDS']
    THREE_SCORE_COORDS = args['THREE_SCORE_COORDS']
    BP_REGION_COORDS = args['BP_REGION_COORDS']
    EXONS_AS_FLANKS = args['EXONS_AS_FLANKS']
    MIN_INTRON_LENGTH = args['MIN_INTRON_LENGTH']
    ALLOW_NONCANONICAL = args['ALLOW_NONCANONICAL']
    ALLOW_OVERLAP = args['ALLOW_OVERLAP']
    LONGEST_ONLY = args['LONGEST_ONLY']
    ONLY_SEQS = args['ONLY_SEQS']
    INCLUDE_DUPES = args['INCLUDE_DUPES']
    NO_SEQS = args['NO_SEQS']
    FN_ANNOT = args['FN_ANNOT']
    RUN_DIR = args['RUN_DIR']
    BP_MATRIX_LENGTH = args['BP_MATRIX_LENGTH']
    SIMPLE_NAME = args['SIMPLE_NAME']
    SPCS = args['SPCS']
    FN_DUPE_MAP = args['FN_DUPE_MAP']
    NUM_NC_USER = args['NUM_NC_USER']
    NUM_NC_LOG = args['NUM_NC_LOG']
    FN_SEQS = args['FN_SEQS']
    FN_OVERLAP_MAP = args['FN_OVERLAP_MAP']
    START_TIME = args['START_TIME']

    # Iterate over generator with transient full sequences
    # Keep your wits about you here, given the number of flags at play
    for intron in get_sub_seqs(
        all_introns,
        flat_annots,
        GENOME,
        INTRON_FLANK_SIZE,
        FIVE_SCORE_COORDS,
        THREE_SCORE_COORDS,
        BP_REGION_COORDS,
        EXONS_AS_FLANKS
    ):
        # Set omission status before generating headers
        intron.omit_check(
            MIN_INTRON_LENGTH, 
            BP_MATRIX_LENGTH,
            ALLOW_NONCANONICAL,
            ALLOW_OVERLAP, 
            LONGEST_ONLY
        )

        #TODO the problem is that for a duplicate short intron, say, 
        # since it is omitted at this point (which we want) it 
        # won't go into add_tags() and therefore won't be tagged
        # as a duplicate. But, sending omitted introns to add_tags()
        # won't work as-is because we don't want to make decisions
        # on other introns based on omitted intron values.
        # One way would be to run a series of checks within add_tags()
        # to make sure the intron isn't omitted, but that's hacky.
        # Unsure what to do...

        # UPDATE: modified add_tags() to partition things by 
        # omission status — seems to be working.

        # add tags with filtering internal to add_tags
        intron, intron_index, longest_isoforms = add_tags(
            intron, 
            intron_index, 
            longest_isoforms,
            ALLOW_OVERLAP, 
            LONGEST_ONLY)

        # this step seems redundant but it's not!
        intron.omit_check(
            MIN_INTRON_LENGTH, 
            ALLOW_NONCANONICAL,
            ALLOW_OVERLAP, 
            LONGEST_ONLY)

        if intron.noncanonical:
            # tag non-canonical introns in label even if being scored
            intron.dynamic_tag.add('[n]')
        if not intron.longest_isoform:
            intron.dynamic_tag.add('[i]')
        if not ONLY_SEQS:
            if not intron.omitted and intron.duplicate is False:
                final_introns.append(intron)
                # compute hash of intron sequence to identify introns
                # across different annotations
            # if intron.duplicate is False:
                    # intron.sha1 = sha1(intron.seq.encode('utf-8')).digest()
        intron_name_index[intron.unique_num] = intron.get_name(SPCS, SIMPLE_NAME)
        if intron.duplicate is not False:
            duplicate_count += 1
            dupe_intron_index[intron.duplicate].add(intron.get_name(SPCS, SIMPLE_NAME))
            if intron.corrected:
                corrected_duplicates.append(intron)
            if not INCLUDE_DUPES:
                continue
        # keep a index of introns which are being omitted due to
        # overlap
        if intron.omitted == 'v':
            overlap_index[intron.overlap].add(intron.get_name(SPCS, SIMPLE_NAME))
        if intron.corrected:
            stats_corrected_count += 1
        # if intron is non-canonical (omitted or not), add to nc stats
        if intron.noncanonical and intron.omitted not in ('s', 'a'):
        # if intron.noncanonical and not intron.omitted:
            dnts = '-'.join(intron.dnts)
            stats_nc_types[dnts] += 1
        if ONLY_SEQS: # no scoring info in intron label
            intron.omitted = False
            intron.name = intron.get_name(SPCS, SIMPLE_NAME)
            for f, tag in zip([bed_file, meta_file], ['BED', 'META']):
                out_string = output_format(intron, tag, SPCS, SIMPLE_NAME)
                f.write(out_string + '\n')

        # write omitted introns to list file (unscored)
        elif intron.omitted and intron.duplicate is False:
            stats_omitted_counts[intron.omitted] += 1
            if intron.omitted == 'n':
                dnts = '-'.join(intron.dnts)
                intron.dynamic_tag.add(dnts)
            for f, tag in zip([bed_file, meta_file], ['BED', 'META']):
                out_string = output_format(intron, tag, SPCS, SIMPLE_NAME)
                f.write(out_string + '\n')

            # odds of two introns producing the same 128-bit hash in a 
            # set of 1 trillion introns is ~1.44e-15
            # (via https://www.ilikebigbits.com/2018_10_20_estimating_hash_collisions.html#toc1.7)
            # intron.sha1.hex()

        if not NO_SEQS:
            # Write all introns to file
            seq_string = output_format(intron, 'SEQ', SPCS, SIMPLE_NAME)
            seq_file.write(seq_string + '\n')

        intron.seq = None  # reduce memory footprint

    if dupe_intron_index:
        with open(FN_DUPE_MAP, 'w') as dupe_map:
            for chosen, dupes in dupe_intron_index.items():
                name = intron_name_index[chosen]
                for d in dupes:
                    dupe_map.write('{}\t{}\n'.format(name, d))

    if not ALLOW_OVERLAP and overlap_index:
        with open(FN_OVERLAP_MAP, 'w') as ol_map:
            for chosen, overlapping in overlap_index.items():
                name = intron_name_index[chosen]
                for o in overlapping:
                    ol_map.write('{}\t{}\n'.format(name, o))
    
    # close any files that were opened
    for f in opened_files:
        f.close()

    if not INCLUDE_DUPES:
        write_log(
            '{} introns with redundant coordinates excluded', 
            duplicate_count)

    # If they only want the intron sequences, exit after writing seq file
    if ONLY_SEQS:
        if not INCLUDE_DUPES:
            total_count = total_count - duplicate_count
        # info_file.close()
        # meta_file.close()
        # seq_file.close()
        # bed_file.close()
        write_log(
            '{} intron sequences written to {}',
            (total_count),
            FN_SEQS)
        run_time = get_runtime(START_TIME)
        write_log('Run finished in {}', run_time)
        sys.exit(0)

    if stats_omitted_counts:
        log_omitted(stats_omitted_counts, MIN_INTRON_LENGTH)

    if stats_nc_types:
        write_log(
            'Most common non-canonical splice sites:',
        )
        for stat in counter_format_top(stats_nc_types, NUM_NC_USER):
            print('[#] {}'.format(stat))
        for stat in counter_format_top(stats_nc_types, NUM_NC_LOG):
            write_log(stat, level='warning')  # was debug

    if not final_introns:
        write_log(
            '[!] ERROR: No intron sequences found. '
            'Check for annotation/genome headers mismatch. Exiting now.')
        sys.exit(1)

    # Correct coordinates of defining features in annotation if any
    # misannotated AT-AC introns are found
    correct_annotation(
        final_introns + corrected_duplicates, 
        flat_annots, 
        ext_args.annotation, 
        RUN_DIR, 
        FN_ANNOT
    )
    
    return final_introns

def plot_figures(finalized_introns, args, logging):
    # build lists of component scores for making plots
    score_vector = np.array(
        [(i.five_z_score, i.bp_z_score) for i in finalized_introns])

    write_log('Generating figures')
    # disable logging temporarily to avoid matplotlib writing
    # to root logging stream
    logging.disable(logging.INFO)
    # plot intron component scores with and without marked U12s
    density_hexplot(
        score_vector,
        '{}.plot.hex'.format(args['SPECIES_FULL']),
        xlab='5\' z-score',
        ylab='BPS z-score',
        fig_dpi=args['FIG_DPI']
    )
    scatter_plot(
        finalized_introns,
        score_vector,
        '{}.plot.scatter'.format(args['SPECIES_FULL']),
        xlab='5\' z-score',
        ylab='BPS z-score',
        threshold=args['THRESHOLD'],
        fig_dpi=args['FIG_DPI']
    )
    summary_scores = [i.svm_score for i in finalized_introns]
    hist_title = '{}.plot.score_histogram'.format(args['SPECIES'])
    histogram(
        summary_scores,
        args['THRESHOLD'],
        bins=100, 
        title=hist_title, 
        fig_dpi=args['FIG_DPI']
    )
    # re-enable logging
    logging.disable(logging.NOTSET)

def write_pwm_info(finalized_introns, args):

    SPECIES = args['SPECIES']
    MATRICES = args['MATRICES']
    SIMPLE_NAME = args['SIMPLE_NAME']
    SPCS = args['SPCS']

    matrix_score_info_fn = '{}.pwm_score_info.iic'.format(SPECIES)
    matrix_score_file = open(matrix_score_info_fn, 'w')
    scored_seqs = ['five_seq', 'bp_seq', 'three_seq']
    region_tags = ['five', 'bp', 'three']
    starts = [args['FIVE_SCORE_COORDS'][0], 0, args['THREE_SCORE_COORDS'][0]]
    # map to demote_info attribute indices for different scored seqs
    demote_seq_map = {
        'five': -3,
        'bp': -2,
        'three': -1
    }
    for index, intron in enumerate(
        sorted(finalized_introns, key=attrgetter('svm_score')), start=1):
        u12_key = intron.u12_matrix  # e.g. ('u12', 'gtag')
        u2_key = intron.u2_matrix
        if getattr(intron, 'demote_info', None) is not None:
            u12_key = ('u12', 'gtag')
            u2_key = ('u2', 'gtag')
        score_strings = []
        scoring_matrices = {
            u12_key: [],
            u2_key: []
        }
        name = intron.get_name(SPCS, SIMPLE_NAME)
        header = '\t'.join([name, str(intron.svm_score)])
        matrix_score_file.write('>{}\n'.format(header))
        for seq_name, tag, start in zip(scored_seqs, region_tags, starts):
            try:
                seq = intron.demote_info[demote_seq_map[tag]]
            except AttributeError:
                seq = getattr(intron, seq_name)
            u12_region_key = u12_key + (tag,)
            u2_region_key = u2_key + (tag,)

            if u12_region_key not in MATRICES:
                u12_region_key = next(
                    m for m in intron.matrices if 
                    all(e in m for e in u12_region_key)
                )
            if u2_region_key not in MATRICES:
                u2_region_key = next(
                    m for m in intron.matrices if 
                    all(e in m for e in u2_region_key)
                )
            scoring_matrices[u12_key].append('-'.join(u12_region_key))
            scoring_matrices[u2_key].append('-'.join(u2_region_key))
            u12_m = MATRICES[u12_region_key]
            u2_m = MATRICES[u2_region_key]
            u12_m_scores = ''.join(mark_seq_score(seq, u12_m, start))
            u2_m_scores = ''.join(mark_seq_score(seq, u2_m, start))
            output_string = (u12_m_scores, seq, u2_m_scores)
            score_strings.append(output_string)
        
        score_strings.append(
            (
                ','.join(scoring_matrices[u12_key]), 
                '', 
                ','.join(scoring_matrices[u2_key])
            )
        )
        
        for line_bits in zip(*score_strings):
            matrix_score_file.write('\t'.join(line_bits) + '\n')
    matrix_score_file.close()

def main():
    # Capture start time of run to facilitate calculation of total runtime at end
    START_TIME = time.time()
    # set better forking method for Unix-based OSes
    fork_types = get_all_start_methods()
    if 'forkserver' in fork_types:
        set_start_method('forkserver')
    elif 'spawn' in fork_types:
        set_start_method('spawn')
    warnings.filterwarnings('ignore')
    parser = make_parser()
    EXT_ARGS = get_args(sys.argv, parser)
    ARGS = get_custom_args(EXT_ARGS, sys.argv)
    ARGS['START_TIME'] = START_TIME
    ARGS = add_pwm_args(ARGS)

    # Logging setup ###########################################################

    # Default config used as a base
    # Debug messasges go only to log file, not screen
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=ARGS['FN_LOG'],
                        format='[#] %(asctime)s | %(message)s',
                        datefmt='%Y.%m.%d-%H.%M.%S',
                        level=logging.WARNING,  # was DEBUG
                        filemode='w')
    # Add logging module for printing to stdout
    screenlogger = logging.StreamHandler(stream=sys.stdout)
    screenlogger.setLevel(logging.CRITICAL)  # was INFO
    # Set less-verbose print syntax for screen vs. log file
    screenformatter = logging.Formatter('[#] %(message)s')
    screenlogger.setFormatter(screenformatter)
    # Add screen logger to root logger
    logging.getLogger('').addHandler(screenlogger)

    # /Logging setup ##########################################################

    write_log('Starting run on {}', ARGS['SPECIES_NAME_INFO'])

    # Determine whether we can plot figures, assuming it is asked for
    if ARGS['WANT_PLOT'] and not ARGS['PLOT']:
        write_log('Matplotlib not detected; plotting disabled')

    # Start processing files
    full_path_args = [
        os.path.abspath(a) if os.path.exists(a) else a for a in sys.argv
    ]
    write_log('Run command: {}', ' '.join(full_path_args))


    # /Global constants ##########################################################

    # Initial data collection ####################################################

    # Collect reference introns, with flanking exonic region
    # on either end
    ARGS['REF_U2S'] = get_reference_introns(
        ref_file=ARGS['REFERENCE_U2_FILE'],
        args=ARGS,
        type_id='u2'
    )
    ARGS['REF_U12S'] = get_reference_introns(
        ref_file=ARGS['REFERENCE_U12_FILE'],
        args=ARGS,
        type_id='u12'
    )

    # remove any reference data with redundant scores
    # filter_attributes = ['five_seq', 'bp_region_seq', 'three_seq']
    # REF_U12S = list(unique_attributes(REF_U12S, filter_attributes))
    # REF_U2S = list(unique_attributes(REF_U2S, filter_attributes))

    ARGS['REFS'] = ARGS['REF_U12S'] + ARGS['REF_U2S']

    all_introns, total_count, flat_annots = introns_from_args(ARGS, EXT_ARGS)

    if not ARGS['SEQUENCE_INPUT']:
        final_introns = filter_introns_write_files(
            all_introns, 
            flat_annots, 
            total_count, 
            ARGS, 
            EXT_ARGS
        )
    else:
        final_introns = all_introns

    FINAL_INTRON_COUNT = len(final_introns)

    write_log('{} introns included in scoring analysis', FINAL_INTRON_COUNT)

    # /Initial data collection ###################################################

    # Scoring ####################################################################
    #TODO make cmdline arg to require bp be built from all introns (in case is
    # u2 matrix in matrices file but don't want to use it)

    BPX = None

    # [ v1 ]
    # secret species introns conserved as U12s in other species, using the
    # pattern '.{4}TTTGA.{3}.{6,8}$'
    # BPX = {
    #     6: 1.371429,
    #     7: 1.428571,
    #     8: 1.2
    # }
    # [ v2 ]
    # secret species introns conserved as U12s in other species plus ATACs
    # from conserved regions using the pattern '.{4}TTTGA.{3}.{6,8}$'
    # BPX = {
    #     6: 1.354167, # 17/48
    #     7: 1.4375, # 21/48
    #     8: 1.208333  # 10 /48
    # }
    # [ v3 ]
    # secret species introns conserved as U12s in others and BUSCO matches,
    # plus ATACs in conserved regions
    # BPX = {
    #     6: 1.352941, # 54/153
    #     7: 1.431373, # 66/153
    #     8: 1.215686  # 33/153
    # }

    ###!!! FIRST ROUND OF SCORING
    MATRICES = add_u2_matrix(
        final_introns,
        ARGS
    )

    # Get the maximum and minimum score possible for each matrix, to be used
    # to scale scores
    # MATRIX_BOUNDS = defaultdict(dict)
    # for k, m in MATRICES.items():
    #     category = k[:2]
    #     sub_category = k[-1]
    #     MATRIX_BOUNDS[category][sub_category] = get_score_bounds(m)
    # ARGS['MATRIX_BOUNDS'] = MATRIX_BOUNDS

    ARGS['MATRICES'] = MATRICES

    write_log(
        'Scoring introns using the following regions: {}',
        ', '.join(ARGS['SCORING_REGIONS'])
    )

    score_label_table = {
        'five': 'five_z_score',
        'bp': 'bp_z_score',
        'three': 'three_z_score'
    }

    # ensure that scoring regions are sorted in the correct order for plotting
    scoring_region_labels = []
    label_order = ['five', 'bp', 'three']
    for label in label_order:
        if label in ARGS['SCORING_REGIONS']:
            scoring_region_labels.append(score_label_table[label])
    ARGS['SCORING_REGION_LABELS'] = scoring_region_labels
    raw_score_names = ['five_raw_score', 'bp_raw_score', 'three_raw_score']
    z_score_names = ['five_z_score', 'bp_z_score', 'three_z_score']

    finalized_introns, model, u12_count, atac_count, demoted_swaps = apply_scores(
        ARGS['REFS'], final_introns, ARGS['MATRICES'], scoring_region_labels, ARGS)

    ###!!! /FIRST ROUND OF SCORING

    if ARGS['RECURSIVE'] and u12_count > 5:
        ARGS['U12_COUNT'] = u12_count
        r_refs, r_introns, r_matrices = recursive_scoring(
            finalized_introns, 
            ARGS['REFS'], 
            model, 
            ARGS,
            raw_score_names, 
            z_score_names)
        
        (
            finalized_introns, 
            model, 
            u12_count, 
            atac_count, 
            demoted_swaps
        ) = apply_scores(
            r_refs, 
            r_introns, 
            r_matrices,
            scoring_region_labels,
            ARGS
        )

    write_log(
        '{} putative AT-AC U12 introns found.', atac_count)
    write_log(
        '{} putative U12 introns found with scores > {}%', 
        u12_count, 
        ARGS['THRESHOLD']
    )

    # /Scoring ###################################################################

    # format and write output files

    if demoted_swaps:
        with open(ARGS['FN_SWAP'], 'w') as swaps:
            for info in demoted_swaps:
                swaps.write(
                    '{}\t{:.5f} ({:.5f}, {:.5f}, {:.5f}) --> '
                    '{:.5f} ({:.5f}, {:.5f}, {:.5f})\n'.format(*info))

    # Write matrices to file for reference,
    # removing added pseudocounts first
    print_matrices = {
        k: add_pseudos(v, pseudo=-ARGS['PSEUDOCOUNT'])
        for k, v in ARGS['MATRICES'].items()}
    write_matrix_file(print_matrices, ARGS['FN_MATRICES'], precision=15)
    
    SPCS = ARGS['SPCS']
    SIMPLE_NAME = ARGS['SIMPLE_NAME']
    
    if not ARGS['SEQUENCE_INPUT']:
        bed_file = open(ARGS['FN_BED'], 'a')
        for i in finalized_introns:
            bed_string = output_format(i, 'BED', SPCS, SIMPLE_NAME)
            # bed_string = write_format(
            #     i, 'region', 'strand', 'start', 'stop', 'get_label', fasta=False)
            bed_file.write(bed_string + '\n')
        bed_file.close()

    # ranking_file = open(FN_RANKINGS, 'w')
    score_file = open(ARGS['FN_SCORE'], 'w')
    meta_file = open(ARGS['FN_META'], 'a')

    if ARGS['SEQUENCE_INPUT']:
        ref_format = True
    else:
        ref_format=False

    #TODO add headers to each output file with column names
    for index, intron in enumerate(sorted(
            finalized_introns, 
            key=attrgetter('svm_score', 'parent', 'index')), start=1
        ):

        score_string = output_format(
            intron, 'SCORE', SPCS, SIMPLE_NAME)    
        score_file.write(score_string + '\n')
        meta_string = output_format(
            intron, 'META', SPCS, SIMPLE_NAME)
        meta_file.write(meta_string + '\n')

    score_file.close()
    meta_file.close()
    # ranking_file.close()

    if not ARGS['NO_SEQS'] and not ARGS['SEQUENCE_INPUT']:
        write_log('Adding scores to intron sequences file')
        add_scores(finalized_introns, ARGS['FN_SEQS'], SPCS)

    if ARGS['PLOT']:
        plot_figures(finalized_introns, ARGS, logging)

    # per-matrix raw scoring information for debugging
    if ARGS['MATRIX_SCORE_INFO']:
        write_pwm_info(finalized_introns, ARGS)

    run_time = get_runtime(START_TIME)

    write_log('Run finished in {}', run_time)

    logging.shutdown()

    sys.exit(0)

# SVM functions ##############################################################

def rank_ones(cv_results, key):
    ranks = cv_results['rank_test_score']
    params = cv_results['params']
    rank_one_params = [e[0][key] for e in zip(params, ranks) if e[1] == 1]

    return rank_one_params


def index_of_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()

    return idx


## v1 (working, custom rank-1 parameter search)
def rank1_param_avg(models):
    """
    Returns geometric mean of rank-1 parameter values across all models
    """
#     params = defaultdict(Counter)
    params = defaultdict(list)
    for m in models:
        for p, v in m.best_params_.items():
            best_params = rank_ones(m.cv_results_, p)
            params[p].extend(best_params)
    refined_params = {}
    for p, v in params.items():
        # use geometric mean to find a value in-between the bounds of the range
        best_avg_param = pystats.gmean(v)
        refined_params[p] = best_avg_param

    return refined_params

# ## v2, uses default sklearn infrastructure and suffers from returing
# the first rank-1 result
# def rank1_param_avg(models):
#     """
#     Returns geometric mean of rank-1 parameter values across all models
#     """
# #     params = defaultdict(Counter)
#     params = defaultdict(list)
#     for m in models:
#         for p, v in m.best_params_.items():
#             # best_params = rank_ones(m.cv_results_, p)
#             params[p].append(v)
#     refined_params = {}
#     for p, v in params.items():
#         # use geometric mean to find a value in-between the bounds of the range
#         best_avg_param = pystats.gmean(v)
#         refined_params[p] = best_avg_param

#     return refined_params

# iterations
def train_svm(
    base_model,
    u2_vector,
    u12_vector,
    iterations,
    cv_parameters=None,
    train_fraction=0.80,  # fraction of the training data used to train SVM
    seed=42):
    u12s = np.array(u12_vector)
    u2s = np.array(u2_vector)
    np.random.seed(seed)
    if iterations:
        subset_size = len(u12_vector)
        possible_iterations = math.floor(len(u2s) / len(u12s))
        if possible_iterations < iterations:
            iterations = possible_iterations
            write_log(
                '[!] Training iterations limited to {} '
                'due to size of reference U2 set', iterations)
        u2_sets = np.random.choice(
            len(u2s), size=(iterations, subset_size), replace=False)
        u2_sets = u2s[u2_sets]
        zero_labels = np.zeros(subset_size)
    else:
        base_model.class_weight = 'balanced'
        # base_model.class_weight = 'balanced_subsample'  # random forest
        u2_sets = [u2s]
        zero_labels = np.zeros(len(u2s))
    trained_models = []
    performance_scores = {
        'f1': [],
        'pr_auc': [],
        'pr_curve': [],
        'classification_report': []
    }
    input_labels = np.concatenate([zero_labels, np.ones(len(u12s))])
    for iter_n, u2_subset in enumerate(u2_sets, start=1):
        if iterations:
            prog_bar = progress_bar(iter_n, iterations)
            percent_complete = round(iter_n / iterations * 100, 2)
            sys.stdout.write(
                '\rTraining progress: {} - {}% ({}/{}) complete'.format(
                    prog_bar, percent_complete, iter_n, iterations))
        if cv_parameters is not None:
            model = GridSearchCV(
                base_model,
                error_score=np.nan,
                **cv_parameters
            )
        else:
            model = clone(base_model)  # normal
            # model = base_model   # linearSVC
        input_feats = np.concatenate([u2_subset, u12s])
        (train_scores, test_scores,
         train_labels, test_labels) = train_test_split(
            input_feats,
            input_labels,
            train_size=train_fraction,
            stratify=input_labels,
            random_state=seed)
        
        model.fit(train_scores, train_labels)
        predict_labels = model.predict(test_scores)
        subset_f1 = f1_score(test_labels, predict_labels)
        # subset_auc = roc_auc_score(test_labels, predict_labels)
        # roc_curve_data = roc_curve(test_labels, predict_labels)
        prob_labels = model.predict_proba(test_scores)[:,1]
        pr_curve_data = precision_recall_curve(
            test_labels, prob_labels)
        precision, recall = pr_curve_data[:2]
        subset_auc = auc(recall, precision)        
        class_report = classification_report(
            test_labels, predict_labels, target_names=['U2', 'U12'])

        performance_scores['f1'].append(subset_f1)
        performance_scores['pr_auc'].append(subset_auc)
        performance_scores['pr_curve'].append(pr_curve_data)
        performance_scores['classification_report'].append(class_report)
        trained_models.append(model)
    
    if iterations:
        print()  # leave space after last progress bar update

    return trained_models, performance_scores

# iterations
def optimize_svm(
    base_model,
    u12_vector,
    u2_vector,
    scorer='balanced_accuracy',
    # scorer='f1',  # performs same as balanced_accuracy but has warnings
    n_optimize=5,
    range_subdivisions=100,  # used by np.geomspace to build parameter ranges
    iterations=0,
    seed=42,
    cv_jobs=1,
    hyper_C=None):
    if hyper_C is not None:
        avg_hyper = hyper_C
    else:
        parameter_range_start = -6
        parameter_range_stop = 6
        log_intervals = np.logspace(
            parameter_range_start,
            parameter_range_stop,
            num=abs(parameter_range_stop - parameter_range_start) + 1)

        initial_parameters = {'C': log_intervals}

        # make things run faster for testing
        # u2_vector = random.sample(list(u2_vector), 5000)   ###!!!

        cv_params = {
            'cv': 5,
            'iid': False,
            'scoring': scorer,
            'n_jobs': cv_jobs,
            'param_grid': initial_parameters
        }
        refined_params = initial_parameters
        for search_round in range(1, n_optimize + 1):
            round_seed = seed + search_round
            print('Starting optimization round {}/{}'.format(
                    search_round, n_optimize))
            search_model, performance = train_svm(
                base_model,
                u2_vector,
                u12_vector,
                iterations=iterations,
                cv_parameters=cv_params,
                seed=round_seed)
            
            # print out cross-validation performance stats for C
            # import pandas as pd
            # with open('{}.r{}_stats.txt'.format(SPECIES, search_round), 'w') as f:
            #     for m in search_model:
            #         df = pd.DataFrame(m.cv_results_)
            #         df.to_string(f)

            # sklearn's default behavior for returning the "best" params
            # for a model (via the best_params_ attribute) will simply return 
            # the first parameter value out of a set of equally-well-performing
            # values. For example, if we are optimizing <C>, and values 1, 10,
            # 100, and 1000 are all rank-1 according to the performance metrics,
            # best_params_ will return 1 as the "best" value for <C>. However,
            # a lower value for <C> is less conservative than a larger value. 
            # Therefore, in this specific application (with a single 
            # hyperparameter) taking the average of all equally-well-
            # performing values seems appropriate. 
            best_first_params = rank1_param_avg(search_model)

            for p, v in best_first_params.items():
                best_index = index_of_nearest(refined_params[p], v)
                low_bound = refined_params[p][max(best_index - 1, 0)]
                high_bound = refined_params[p][min(
                    best_index + 1, len(refined_params[p]) - 1)]
                p_range = np.geomspace(low_bound, high_bound, range_subdivisions)
                refined_params[p] = p_range
            cv_params['param_grid'] = refined_params
        
        if iterations:
            print()
        min_hyper = min(cv_params['param_grid']['C'])
        max_hyper = max(cv_params['param_grid']['C'])
        write_log(
            'Range for \'C\' after {} rounds of optimization: {}-{}',
            n_optimize, min_hyper, max_hyper
        )
        avg_hyper = pystats.gmean(cv_params['param_grid']['C'])
        cv_params['param_grid']['C'] = avg_hyper

    write_log(
        'Set classifier value for \'C\': {}', 
        avg_hyper
    )
    write_log('Training classifier with optimized hyperparameters')
    base_model.C = avg_hyper
    trained_model, performance = train_svm(
        base_model,
        u2_vector,
        u12_vector,
        iterations=iterations)

    return trained_model, performance

# single-process
# def svm_predict(score_vectors, model_list):
#     probability_ledger = defaultdict(list)
#     for m in model_list:
#         probabilities = m.predict_proba(score_vectors)
#         labels = m.predict(score_vectors)
#         distances = m.decision_function(score_vectors)
#         for i, lp in enumerate(zip(probabilities, labels, distances)):
#             probability_ledger[i].append(lp)

#     return probability_ledger

# multi-process WIP
def svm_predict(score_vectors, models):
    probabilities, labels, distances = [], [], []
    for m in models:
        p = m.predict_proba(score_vectors)[:,1]
        l = m.predict(score_vectors)
        d = m.decision_function(score_vectors)
        probabilities.append(p)
        labels.append(l)
        distances.append(d)

    return probabilities, labels, distances

# multi-process
# def svm_predict(score_vector, model_list):
#     class_info = []
#     for m in model_list:
#         probability = m.predict_proba(score_vector)[0]
#         label = m.predict(score_vector)[0]
#         distance = m.decision_function(score_vector)[0]
#         class_info.append((probability, label, distance))

#     return class_info


# TODO reformat to accept vector input, incorporate into main script
def assign_clusters(introns, attributes=['five', 'bp']):
    def _get_attrs(i, attr_list):
        return tuple([i.a for a in attr_list])
    attr_map = {
        'five': 'five_z_score',
        'bp': 'bp_z_score',
        'three': 'three_z_score'
    }
    attrs = [v for k, v in attr_map if k in attributes]
    score_vectors = [_get_attrs(i, attrs) for i in introns]
    labels = SpectralClustering(
        n_clusters=2, 
        assign_labels='discretize', 
        n_init=42, 
        gamma=5).fit_predict(score_vectors)
    
    return labels

# # single-process
# def average_svm_score(probability_dict, weights=None):
#     average_probs = defaultdict(dict)
#     if weights is not None:
#         weights = np.array(weights)
#     for k, v in probability_dict.items():
#         u2_probs = np.array([e[0][0] for e in v])
#         u12_probs = np.array([e[0][1] for e in v])
#         distances = np.array([e[2] for e in v])
#         if weights is not None:
#             # weight each average by the performance of its source classifier
#             u2_probs = u2_probs * weights
#             u12_probs = u12_probs * weights
#             distances = distances * weights
#         label_count = Counter([e[1] for e in v])
#         # avg_u2_prob = np.mean(u2_probs)
#         avg_u12_prob = np.mean(u12_probs)
#         # u2_prob_stdev = np.std(u2_probs)
#         # u12_prob_stdev = np.std(u12_probs)
#         # u12_prob_sem = pystats.sem(u12_probs)
#         avg_distance = np.mean(distances)
#         average_probs[k] = {
#             'u12_avg': avg_u12_prob,
#             # 'u12_sem': u12_prob_sem,
#             # 'u2_avg': avg_u2_prob,
#             # 'u12_std': u12_prob_stdev,
#             # 'u2_std': u2_prob_stdev,
#             'avg_distance': avg_distance,
#             'labels': label_count
#         }

#     return average_probs

# multi-process
# def average_svm_score(class_info, weights=None):
#     if weights is not None:
#         weights = np.array(weights)
#     u2_probs = np.array([e[0][0] for e in class_info])
#     u12_probs = np.array([e[0][1] for e in class_info])
#     distances = np.array([e[2] for e in class_info])
#     if weights is not None:
#         # weight each average by the performance of its source classifier
#         u2_probs = u2_probs * weights
#         u12_probs = u12_probs * weights
#         distances = distances * weights
#     label_count = Counter([e[1] for e in class_info])
#     # avg_u2_prob = np.mean(u2_probs)
#     avg_u12_prob = np.mean(u12_probs)
#     # u2_prob_stdev = np.std(u2_probs)
#     # u12_prob_stdev = np.std(u12_probs)
#     # u12_prob_sem = pystats.sem(u12_probs)
#     avg_distance = np.mean(distances)
#     average_probs = {
#         'u12_avg': avg_u12_prob,
#         # 'u12_sem': u12_prob_sem,
#         # 'u2_avg': avg_u2_prob,
#         # 'u12_std': u12_prob_stdev,
#         # 'u2_std': u2_prob_stdev,
#         'avg_distance': avg_distance,
#         'labels': label_count
#         }

#     return average_probs

# multi-process WIP
def average_svm_score_info(probabilities, labels, distances, weights=None):
    if weights:
        weights = np.array(weights)
    if len(probabilities) == 1:
        avg_u12_prob = probabilities[0]
        avg_distance = distances[0]
        # if weights:
        #     w = weights[0]
        #     avg_u12_prob *= w
        #     avg_distance *= w
        label_count = Counter(labels)
    else:
        u12_probs = np.array(probabilities)
        distances = np.array(distances)
        if weights is not None:
            # weight each average by the performance of its source classifier
            # u2_probs = u2_probs * weights
            u12_probs = u12_probs * weights
            distances = distances * weights
        label_count = Counter(labels)
        # avg_u2_prob = np.mean(u2_probs)
        avg_u12_prob = np.mean(u12_probs)
        # u2_prob_stdev = np.std(u2_probs)
        # u12_prob_stdev = np.std(u12_probs)
        # u12_prob_sem = pystats.sem(u12_probs)
        avg_distance = np.mean(distances)
    average_probs = {
        'u12_avg': avg_u12_prob,
        # 'u12_sem': u12_prob_sem,
        # 'u2_avg': avg_u2_prob,
        # 'u12_std': u12_prob_stdev,
        # 'u2_std': u2_prob_stdev,
        'avg_distance': avg_distance,
        'labels': label_count
        }

    return average_probs


def get_attributes(objs, attr_names):
    if type(objs) is not list:
        objs = [objs]
    return [tuple([getattr(o, a) for a in attr_names]) for o in objs]


def get_score_vector(introns, score_names):
    vect = [[getattr(i, n) for n in score_names] for i in introns]
    
    return np.asarray(vect)


def u12_label_ratio(label_dict):
    u12 = label_dict[1]
    total = sum(label_dict.values())
    
    return u12 / total

# single-process
# def assign_svm_scores(
#     introns, 
#     models, 
#     scoring_region_labels, 
#     weights=None,
#     add_type=True):
#     intron_score_vector = get_score_vector(
#             introns, score_names=scoring_region_labels)
#     u12_probability_index = svm_predict(intron_score_vector, models)
#     avg_u12_probabilities = average_svm_score(
#         u12_probability_index, weights)

#     id_map = {
#         0: 'u2',
#         1: 'u12'
#     }

#     for idx, intron in enumerate(introns):
#         intron.svm_score = avg_u12_probabilities[idx]['u12_avg'] * 100

#         intron.score_distance = avg_u12_probabilities[idx]['avg_distance']
        
#         # intron.error = avg_u12_probabilities[idx]['u12_sem'] * 100
#         # relative as percentage of threshold
#         # intron.relative_score = (intron.svm_score - THRESHOLD) / THRESHOLD * 100
#         # relative as probability
#         intron.relative_score = intron.svm_score - THRESHOLD
#         label_ratio = u12_label_ratio(avg_u12_probabilities[idx]['labels'])
#         if add_type is True:
#             if label_ratio > 0.5:
#                 type_id = 'u12'
#             else:
#                 type_id = 'u2'
#             intron.type_id = type_id
#         intron.label_ratio = label_ratio

#     return introns

# def apply_svm_score(introns, info, add_type=True):
#     intron.svm_score = info['u12_avg'] * 100
#     intron.score_distance = info['avg_distance']
    
#     # intron.error = avg_u12_probabilities[idx]['u12_sem'] * 100

#     # relative as percentage of threshold
#     # intron.relative_score = (intron.svm_score - THRESHOLD) / THRESHOLD * 100
#     # relative as probability
#     intron.relative_score = intron.svm_score - THRESHOLD
#     label_ratio = u12_label_ratio(info['labels'])
#     if add_type is True:
#         if label_ratio > 0.5:
#             type_id = 'u12'
#         else:
#             type_id = 'u2'
#         intron.type_id = type_id
#     intron.label_ratio = label_ratio

#     return intron

# multi-process WIP
# def assign_svm_score(
#     introns, 
#     models, 
#     scoring_region_labels,
#     weights=None,
#     add_type=True,
#     processes=1
# ):
#     intron_score_vector = get_score_vector(
#             introns, score_names=scoring_region_labels)
#     probabilities, labels, distances = svm_predict(intron_score_vector, models)

#     p_group = zip(*probabilities)
#     l_group = zip(*labels)
#     d_group = zip(*distances)

#     info_groups = zip(*[p_group, l_group, d_group])
#     score_info = [
#         average_svm_score_info(p, l, d, weights) for (p, l, d) in info_groups
#     ]

#     pool = Pool(processes=processes)
#     scored_introns = pool.starmap(
#         apply_svm_score, zip(
#             introns,
#             score_info,
#             repeat(add_type)
#         )
#     )

#     return scored_introns


def assign_svm_score(
    intron,
    models, 
    scoring_region_labels,
    weights,
    threshold,
    add_type
):
    intron.model_score(
        models, scoring_region_labels, threshold, weights, add_type)
    
    return intron


def parallel_svm_score(
    introns, 
    models, 
    scoring_region_labels,
    threshold,
    weights=None, 
    add_type=True, 
    processes=1
):
    with Pool(processes=processes) as pool:
    # score_func = partial(
    #     assign_svm_score, 
    #     models=models,
    #     scoring_region_labels=scoring_region_labels,
    #     weights=weights,
    #     add_type=add_type)
    # classified_introns = pool.imap(
    #     score_func, introns)
        classified_introns = pool.starmap(
            assign_svm_score, zip(
                introns,
                repeat(models),
                repeat(scoring_region_labels),
                repeat(weights),
                repeat(threshold),
                repeat(add_type)
            )
        )

    return classified_introns


# # multi-process
# def assign_svm_score(
#     intron, 
#     models, 
#     scoring_region_labels, 
#     weights=None,
#     add_type=True):
#     intron_score_vector = get_score_vector(
#             [intron], score_names=scoring_region_labels)
#     u12_probability_index = svm_predict(intron_score_vector, models)
#     avg_u12_probabilities = average_svm_score(
#         u12_probability_index, weights)

#     intron.svm_score = avg_u12_probabilities['u12_avg'] * 100

#     intron.score_distance = avg_u12_probabilities['avg_distance']
    
#     # intron.error = avg_u12_probabilities[idx]['u12_sem'] * 100

#     # relative as percentage of threshold
#     # intron.relative_score = (intron.svm_score - THRESHOLD) / THRESHOLD * 100
#     # relative as probability
#     intron.relative_score = intron.svm_score - THRESHOLD
#     label_ratio = u12_label_ratio(avg_u12_probabilities['labels'])
#     if add_type is True:
#         if label_ratio > 0.5:
#             type_id = 'u12'
#         else:
#             type_id = 'u2'
#         intron.type_id = type_id
#     intron.label_ratio = label_ratio

#     return intron

# def parallel_svm_score(
#     introns, models, scoring_region_labels, weights=None, add_type=True, processes=1):
#     pool = Pool(processes=processes)
#     classified_introns = pool.starmap(
#         assign_svm_score, zip(
#             introns, 
#             repeat(models), 
#             repeat(scoring_region_labels),
#             repeat(weights),
#             repeat(add_type)
#         )
#     )

#     return classified_introns

# /SVM functions ######################################################

# Plotting functions #########################################################

def check_plot():
    # check for the plotting library required to produce
    # optional figures
    try:
        import matplotlib
        can_plot = True
    except ModuleNotFoundError:
        can_plot = False

    return can_plot

def density_hexplot(
    scores, title, xlab=None, ylab=None, outfmt='png', fsize=14, fig_dpi=300):
    # plt.rcParams['figure.figsize'] = [12, 10.2]
    plt.figure(figsize=(8, 8))
    ax = plt.gca()
    hx = ax.hexbin(
            *scores.T, mincnt=1, cmap='inferno', 
            bins='log', linewidths=0)
    title_with_n = '{} (n={})'.format(title, len(scores))
    if xlab:
        plt.xlabel(xlab, fontsize=fsize)
    if ylab:
        plt.ylabel(ylab, fontsize=fsize)
    plt.title(title_with_n, fontsize=fsize)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(hx, cax=cax)
    cb.set_label('Bin density (log10(n))')
    title = '_'.join(title.split())
    plt.savefig(
        "{}.iic.{}".format(title, outfmt), 
        format=outfmt, dpi=fig_dpi, bbox_inches='tight')
    ###!!!
    # outfmt = 'pdf'
    # plt.savefig(
    #     "{}.iic.{}".format(title, outfmt), 
    #     format=outfmt, dpi=fig_dpi, bbox_inches='tight')
    ###!!!
    plt.close()


def scatter_plot(
    introns, 
    scores, 
    title, 
    xlab, 
    ylab,
    threshold,
    fsize=14, 
    outfmt='png', 
    fig_dpi=300):
    plt.figure(figsize=(8, 8))
    cluster_colors = []
    u2_count, u12_low, u12_med, u12_high = [0] * 4
    score_stdev = np.std([i.svm_score for i in introns])
    high_val = threshold
    med_val = threshold - score_stdev
    # low_val = THRESHOLD - (score_stdev * 2)
    for i in introns:
        itype = i.type_id
        p = i.svm_score
        if itype == 'u2':
            u2_count += 1
            color = 'xkcd:medium grey'
        elif p > high_val:
            u12_high += 1
            color = 'xkcd:green'
        elif med_val < p <= high_val:
            u12_med += 1
            color = 'xkcd:orange'
        elif p <= med_val:
            u12_low += 1
            color = 'xkcd:red'
        cluster_colors.append(color)

    legend_colors = [
        'xkcd:medium grey', 'xkcd:red', 'xkcd:orange', 'xkcd:green']
    # legend_labels = ['U2', 'U12<=68', '68<U12<=95', 'U12>95']
    round_vals = map(round, [high_val, med_val])
    legend_labels = [
        'U2', 
        'U12<={}'.format(int(med_val)), 
        '{}<U12<={}'.format(int(med_val), int(high_val)), 
        'U12>{}'.format(int(high_val))
    ]
    legend_counts = [u2_count, u12_low, u12_med, u12_high]
    legend_patches = []
    for label, count, color in zip(legend_labels, legend_counts, legend_colors):
        label = '{} ({})'.format(label, count)
        patch = mpatches.Patch(color=color, label=label)
        legend_patches.append(patch)
    plt.scatter(
        *scores[:,:2].T, s=20, c=cluster_colors, alpha=0.5, rasterized=True)
    # plt.axes().set_aspect('equal')   ###!!!!
    plt.legend(handles=legend_patches)
    plt.xlabel(xlab, fontsize=fsize)
    plt.ylabel(ylab, fontsize=fsize)
    plt.title(title, fontsize=fsize)
    plt.tight_layout()
    plt.savefig(
        '{}.iic.{}'.format(title, outfmt), 
        format=outfmt, dpi=fig_dpi, bbox_inches='tight')

    # ###!!!
    # outfmt = 'pdf'
    # plt.savefig(
    #     '{}.iic.{}'.format(title, outfmt), 
    #     format=outfmt, dpi=fig_dpi, bbox_inches='tight')
    plt.close()


def ref_scatter(
    u2_vector, 
    u12_vector, 
    clf, 
    scoring_regions, 
    species, 
    fsize=14, 
    fig_dpi=300
):
    plt.figure(figsize=(8, 8))
    plt.scatter(
        *u2_vector[:,:2].T,
        c='xkcd:medium grey',
        alpha=0.5,
        s=42,
        label='U2 (n={})'.format(len(u2_vector)), 
        rasterized=True)
    plt.scatter(
        *u12_vector[:,:2].T, 
        c='xkcd:green',
        alpha=0.5,
        s=42,
        label='U12 (n={})'.format(len(u12_vector)),
        rasterized=True)
    plt.xlabel('5\' z-score', fontsize=fsize)
    plt.ylabel('BPS z-score', fontsize=fsize)

    # only plot decision function if 2 scoring regions used
    if len(scoring_regions) == 2:
        # plot the decision function
        ax = plt.gca()
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        # create grid to evaluate model
        xx = np.linspace(xlim[0], xlim[1], 30)
        yy = np.linspace(ylim[0], ylim[1], 30)
        YY, XX = np.meshgrid(yy, xx)
        xy = np.vstack([XX.ravel(), YY.ravel()]).T
        Z = clf.decision_function(xy).reshape(XX.shape)

        # plot decision boundary and margins
        ax.contour(XX, YY, Z, colors='k', levels=[-1, 0, 1], alpha=0.5,
                linestyles=['--', '-', '--'])
        # plot support vectors
        ax.scatter(clf.support_vectors_[:, 0], clf.support_vectors_[:, 1], s=100,
                linewidth=1, facecolors='none', edgecolors='k', rasterized=True)
    plt.axes().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(
        '{}.plot.training_scatter.iic.png'.format(species), 
        format='png', 
        dpi=fig_dpi)

    ###!!!
    # plt.savefig(
        # '{}.plot.training_scatter.iic.pdf'.format(SPECIES), format='pdf', dpi=fig_dpi)
    plt.close()


def histogram(data_list, threshold, title=None, grid=True, bins='auto', log=True, fig_dpi=300):
    plt.figure(figsize=(10, 6))
    if log is True:
        plt.yscale('log')
    plt.hist(data_list, bins=bins)
    if grid:
        plt.grid(True, which="both", ls="--", alpha=0.7)
    if title is not None:
        plt.title(title, fontsize=14)
    plt.xlabel('U12 score', fontsize=14)
    plt.ylabel('Number of introns', fontsize=14)
    plt.axvline(
        threshold, color='orange', linestyle='--',
        label='U12 threshold: {}'.format(threshold))
    plt.legend()
    plt.tight_layout()
    # plt.savefig('{}.iic.png'.format(title.replace(' ', '_')), dpi=600)
    plt.savefig('{}.iic.png'.format(title), dpi=fig_dpi)
    plt.close()

# /Plotting functions ########################################################

# /Functions #################################################################

# Main #######################################################################

if __name__ == '__main__':
    main()

# TODO report types of features that were skipped
# TODO modify functions to return dictionaries for extensibility?
