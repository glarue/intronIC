import sys
import re
from collections import defaultdict
from biogl import flex_open
import argparse

def format_matrix(matrix, label="frequencies", precision=None):
    """
    Formats the contents of a matrix in FASTA
    format, with the first line being the
    order of characters, and the
    following lines containing the frequency of
    each character; each position in the sequence
    has its own line in the matrix.

    {label} is used as the header for the frequencies
    entry.

    example output:

    >{label}
    A C G T
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

def canonical_bounds(dnts, subtype):
    canonical = {
        'u2': [
            ('GT', 'AG'),
            ('GC', 'AG'),
        ],
        'u12': [
            ('GT', 'AG'),
            ('GC', 'AG'),
            ('AT', 'AC')
        ]
    }
    if dnts not in canonical[subtype]:
        return False
    else:
        return True


def valid_seq(seq, invalid=re.compile('[^ACTGN]')):
    if invalid.search(seq):
        return False
    else:
        return True


def seq_field(s, seq_pattern=re.compile(r'[^\W0-9_]')):
    if seq_pattern.fullmatch(s):
        return True
    else:
        return False

def score_field(s, pattern=re.compile(r'[.\-*\s0-9]')):
    if pattern.match(s):
        return True
    else:
        return False

def aggregate_seqs(
    seqfile,
    subtype,
    allow_noncanonical=False,
    five_bounds=(-20, 20), 
    three_bounds=(-20, 20), 
    blank=re.compile(r'[.\-*\s]')
):
    seq_index = {
        'five': defaultdict(list),
        'three': defaultdict(list)
    }
    fives = []
    threes = []
    five_start, five_stop = five_bounds
    three_start, three_stop = three_bounds
    five_length = abs(five_bounds[0] - five_bounds[1])
    three_length = abs(three_bounds[0] - three_bounds[1])
    with flex_open(seqfile) as f:
        for i, l in enumerate(f):
            if l.startswith('#'):
                continue
            bits = l.strip().split('\t')
            # uid = bits[0]
            seq_bits = [b for b in bits[1:5] if not score_field(b)] # skip score if present
            five, seq, three = seq_bits
            five = five.upper()
            three = three.upper()
            seq = seq.upper()
            dnts = (seq[:2], seq[-2:])
            if not canonical_bounds(dnts, subtype) and not allow_noncanonical:
                continue
            # if blank.match(uid):
            #     uid = 'training_intron_{}'.format(i)
            five_seq = five[five_start:] + seq[:five_stop]
            three_seq = seq[three_start:] + three[:three_stop]
            if valid_seq(five_seq) and len(five_seq) == five_length:
                seq_index['five'][dnts].append(five_seq)
            if valid_seq(three_seq) and len(three_seq) == three_length:
                seq_index['three'][dnts].append(three_seq)

    return seq_index

parser = argparse.ArgumentParser(
    description='Generate PWMs from intron sequence file'
)
parser.add_argument(
    'sequence_file',
    help=(
        'Tab-delimited sequence file with columns 1: name (or blank), 2: 5\' seq, '
        '3: intron seq, 4: 3\' seq '
        '(see https://github.com/glarue/intronIC/wiki/Training-data-and-PWMs#reference_u2-u12_setintronsiicgz) '
        'for details'
    )
)
parser.add_argument(
    'intron_subtype',
    choices=('u2', 'u12'),
    type=str
)
parser.add_argument(
    '--three_coords',
    '--3c',
    default=(-20, 20),
    metavar=('start', 'stop'),
    nargs=2,
    type=int,
    help=(
        'Coordinates describing the 3\' sequence to be used, relative to '
        'the 3\' splice site (e.g. position -1 is the last base of the '
        'intron); 0-indexed half-closed interval (start, stop]'
    )
)
parser.add_argument(
    '--five_coords',
    '--5c',
    default=(-20, 20),
    metavar=('start', 'stop'),
    nargs=2,
    type=int,
    help=(
        'Coordinates describing the 5\' sequence to be used, relative to '
        'the 3\' splice site (e.g. position -1 is the last base of the '
        'intron); 0-indexed half-closed interval (start, stop]'
    )
)
parser.add_argument(
    '--non_canonical_bounds',
    '--nc',
    help='allow non-canonical bounds during PWM creation',
    default=False,
    action='store_true'
)

args = parser.parse_args()

seqfile = args.sequence_file
subtype_tag = args.intron_subtype
ALLOW_NC = args.non_canonical_bounds
# seqfile = sys.argv[1]
# subtype_tag = sys.argv[2]  # u2 or u12

five_bounds = args.five_coords
three_bounds = args.three_coords
bounds = {
    'five': five_bounds,
    'three': three_bounds
}

seqs = aggregate_seqs(
    seqfile, 
    subtype_tag, 
    allow_noncanonical=ALLOW_NC,
    five_bounds=bounds['five'],
    three_bounds=bounds['three'])

for region, seq_dict in seqs.items():
    region_bounds = bounds[region]
    for dnts, seq_list in seq_dict.items():
        dnt_tag = ''.join(dnts)
        label = '_'.join([subtype_tag, dnt_tag, region])
        pwm, n_seqs = matrix_from_seqs(seq_list, start_index=region_bounds[0])
        if n_seqs < 1:
            print('ERROR: no seqs found for {}'.format(label), file=sys.stderr)
            continue
        label = '{} (n={})'.format(label.lower(), n_seqs)
        formatted = format_matrix(pwm, label=label)
        print(formatted)
