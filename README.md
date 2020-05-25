![intronIC_logo](https://user-images.githubusercontent.com/6827531/82829967-62872480-9e69-11ea-94e9-fa7306c7df1b.png)

# (__intron__ <ins>I</ins>nterrogator and <ins>C</ins>lassifier)

## __Quick start__

First, install dependencies:

* `python -m pip install numpy scipy matplotlib scikit-learn biogl`

Change to the `test_data` subdirectory, which contains Ensembl annotations for chr 19 of the human genome.

To identify putative U12-type introns:

* `../intronIC -g Homo_sapiens.Chr19.Ensembl_91.fa.gz -a Homo_sapiens.Chr19.Ensembl_91.gff3.gz -n homo_sapiens`

To extract all annotated intron sequences:

* `../intronIC -g Homo_sapiens.Chr19.Ensembl_91.fa.gz -a Homo_sapiens.Chr19.Ensembl_91.gff3.gz -n homo_sapiens -s`

For many more details about options, output file, etc, please read on. And if you find this tool useful, please cite it! http://dx.doi.org/10.1093/nar/gkaa464

## **Overview**

`intronIC` has two primary uses, both of which require a genome and corresponding annotation file (or BED file of intron coordinates):

1. Score all annotated introns against expectations for U12 introns, and classify introns as either U2- or U12-type using a support-vector machine (SVM)-based approach.
2. Retrieve all annotated intron sequences and associated metadata (using the `-s` flag).

Importantly, by default `intronIC` _only processes introns with unique coordinates from the longest annotated isoform for each gene_ (though this behavior is adjustable). Therefore, the same intron from multiple isoforms will only be included once, and named based upon the longest isoform. See `Details` for additional caveats.

`intronIC` should be able to process most annotations (even kinda crappy ones) provided they roughly adhere to [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format) formatting standards, and produces a set of output files described in detail below. At a high level, it works by aggregating all of the CDS/exon sequences under their parent transcripts and/or genes based on the parent-child relationships given by the last column in the annotation file. 

Then, assuming the scoring mode is being used, it will score the introns against configurable position-weight matrices and reference introns, feed that information into an SVM training and optimization routine, and output a set of files containing intron sequences and score information.

A typical default scoring run, which will use both `CDS` and `exon` features to define introns, and will include introns with non-canonical splice boundaries, might go like this:

`intronIC -g {genome} -a {annotation} -n {binomial_name}`

To concretize with real file names, this could be:

`intronIC -g Homo_sapiens.Chr19.Ensembl_91.fa.gz -a Homo_sapiens.Chr19.Ensembl_91.gff3.gz -n homo_sapiens`

#### A note on the `-n` (name) argument:

By default, `intronIC` expects names in binomial (genus, species) form separated by a non-alphanumeric character, e.g. 'homo_sapiens', 'homo.sapiens', etc. `intronIC` then formats that name internally into a tag that it uses to label all output intron IDs, ignoring anything past the second non-alphanumeric character. Output *files*, on the other hand, are named using the full name supplied via `-n`. If you'd prefer to have it leave whatever argument you supply to `-n` unmodified, use the `--na` flag.

If you are running multiple versions of the same species and would like to keep the same species abbreviations in the output intron data, simply add a tag to the end of the name, e.g. "homo_sapiens.v2".

### _Dependencies_

* [Python >=3.3](https://www.python.org/downloads/)
* [numpy & scipy](https://www.scipy.org/scipylib/download.html)
* [scikit-learn >=0.20.1](http://scikit-learn.org/stable/index.html)
* [biogl](https://github.com/glarue/biogl)
* [matplotlib](https://matplotlib.org/) (optional, required for plotting)

To install prerequisites using `pip`:

`python -m pip install numpy scipy matplotlib scikit-learn biogl`

`intronIC` was built and tested on Linux, but should run on Windows or Mac OSes without much trouble (I say that now...).

#### Note on resource usage

For genomes with a large number of annotated introns, memory usage can be on the order of gigabytes. This should rarely be a problem even for most modern personal computers, however. For reference, the Ensembl 95 release of the human genome requires ~5 GB of memory.

For many non-model genomes, `intronIC` should run fairly quickly (e.g. tens of minutes). For human and other very well annotated genomes, runtime may be longer (the human Ensembl release 95 takes ~20-35 minutes in testing); run time scales relatively linearly with the total number of annotated introns, and can be improved by using parallel processes via `-p`.

## **Details**

### _Output files_

`intronIC` will automatically generate a set of output files. A brief description of the contents of each file follows (numbered lists represent the columns within each file).

#### `annotation.iic`

If putatively misannotated introns are found, `intronIC` will correct their coordinates by adjusting the features that define them in this file. These entries will contain a 'shift' tag that indicates the change made to either their start or stop coordinate (or both).

#### `demoted.iic`

Scoring information for putative U12s whose scores fell below the threshold after their boundaries were switched to GT-AG (as a check to avoid scoring non-canonical introns as U12-type by way of superficial similarities to U12-type motifs):

1. label
2. initial score (with five and bp scores in parentheses), followed by reduced score after boundary switching

#### `dupe_map.iic`

A mapping table of unique, scored intron labels and their corresponding duplicate intron labels:

1. scored intron label
2. duplicate intron label

#### `introns.iic`

All of the annotated intron sequences, including any introns not meeting scoring criteria (e.g. too short, including non-ATCG characters in scoring regions, etc.; including duplicate introns if `-d`):

1. label (without score)
2. score (or '`-`' if run with `-s`)
3. upstream (5′ exon) sequence (200 nt by default)
4. intron sequence
5. downstream (3′ exon) sequence (200 nt by default)

#### `bed.iic`

A BED format file of intron coorindates with U12 probabilty scores and labels:

1. genomic region (e.g. 'chr1')
2. start coordinate (0-indexed)
3. stop coordinate (1-indexed)
4. label (including rounded score in label)
5. U12 probability score
6. strand

#### `log.iic`

A log of all of the information generated during operation, including total number of introns processed, omitted, etc., and total number of U12 introns identified.

#### `pwms.iic`

A FASTA file of the PWMs used, including those built from the experimental dataset.

#### `meta.iic`

A hodgepodge of other intron data (ordered by increasing U12 score):

1. intron label
2. relative score based upon score threshold (U12-type: >0)
3. terminal dinucleotides (e.g. 'GT-AG', 'AT-AC')
4. intron motif string (5', BPS, 3')
5. BPS location in relation to the 3' end of the intron
6. intron length
7. parent transcript
8. parent gene
9. ordinal intron position in transcript
10. total introns in transcript
11. fractional position of intron in transcript as a percentage of the coding length, e.g. 50.0 for an intron that interrupts the coding sequence between codons 15 and 16 out of 30.
12. intron phase
13. binary classification made by the classifier ('u12' or 'u2'), which may include many lower-probability U12-type introns

#### `score_info.iic`

Various scoring information (in order of increasing score).

1. intron label
2. relative score (maximum precision); U2 <= 0 < U12
3. SVM-assigned U12 probability score (0-100) (averaged across `N` SVM classifiers if `--subsample_n N`)
4. 5′ sequence used for scoring
5. 5′ raw score
6. 5′ z-score
7. branch point sequence used for scoring
8. branch point raw score
9. branch point z-score
10. 3' sequence used for scoring
11. 3' raw score
12. 3' z-score
13. distance from hyperplane, e.g. the raw classifier output prior to scikit-learn's implementation of [Platt scaling](https://en.wikipedia.org/wiki/Platt_scaling) to convert distances to probabilities

### _Configuration_

Configuration files can be found in the `intronIC_data` directory. They consist of the following set of files:

* `reference_[u2, u12]_set.introns.iic[.gz]` — Intron sequences (**including any flanking exon sequence required for scoring**) to establish a score threshold. These should contain representative introns of both types for whatever species is being analyzed to produce the best results (although the default set should work well in many cases).
  * The default reference introns are from a set found to be conserved across one or more groupings of species. In the U2 set, all introns are conserved between human, mouse and zebrafish, and also between human, macaque and marmoset. The U12 set contains human U12 introns conserved in the previous groups, as well as between human, mouse, chicken and hagfish. In addition, conserved human U12 introns found in the [U12DB](http://genome.crg.es/datasets/u12/) ([Alioto 2007](https://www.ncbi.nlm.nih.gov/pubmed/17082203)) and [SpliceRack](http://katahdin.mssm.edu/splice/splice_matrix.cgi?database=spliceNew) ([Sheth et al. 2006](https://www.ncbi.nlm.nih.gov/pubmed/16914448)), and those with evidence of retention due to U12 spliceosome knockdown ([Madan et al. 2015](https://www.ncbi.nlm.nih.gov/pubmed/25586593), [Niemelä et al. 2014](https://www.ncbi.nlm.nih.gov/pubmed/24848017)) are also included. The identifiers for each intron follow either the naming convention used by intronIC, or begin with the U12DB identifier, and each line ends with a tag indicating the source(s) of their supporting evidence. Multiple independent sources of evidence were required for introns to be included in the U12 training set. See the files themselves for additional details.

* `scoring_matrices.fasta.iic` - A set of position-weight matrices (PWMs) representing different intron motifs for the scoring regions (five prime and branch point), e.g.

   ```
    >u12_atac_five  start=-20       (n=61)
    A       C       G       T
    0.262295082     0.2459016393    0.3278688525    0.1639344262
    0.2295081967    0.262295082     0.1967213115    0.3114754098
    0.2459016393    0.262295082     0.262295082     0.2295081967
    ```

  * The headers of the FASTA records must contain information about the type ("U2", "U12"), sub-type ("GTAG", "ATAC", etc.) and scoring region ("five", "bp", "three"). Order of these tags is not important.
  * The `start=` tag signifies the start location  of the PWM (relative to the first, `0`, or last, `-1`, base of the intron). For example, if the 5' scoring matrix begins with 5 nt of exonic sequence, the header should be tagged with `start=-5`. The first line of each entry must consist only of the order of the bases in the matrix (see default file for a clearer idea).
  * The default scoring PWMs are a combination of introns  supported by multiple publications, as well as introns with clear conservation as U12-type based on internal analyses. Details of the filtering critera is documented within the default matrices file.
* `u2.conserved_empirical_bp_matrix.iic` - A matrix derived from branch point sequences in conserved U2 introns (based on data from [Pineda and Bradley 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5959240)). Used in cases where the number of experimental introns is too low to construct a robust U2 branch point matrix.

#### Note: Any/all of the above files can be replaced as the defaults either permanently by modifying the existing files, or in one-off fashion using the '-m' and '-r[2, 12]' command line arguments

### _Data filtering notes_

1. When scoring introns, `intronIC` only processes introns with unique coordinates, and by default only includes introns from the longest isoform for a given gene. If run with `-i` (to include multiple isoforms), introns with duplicate coordinates are still excluded; in such cases, introns from the **longest isoform** (computed as the sum of the component coding sequences) will be preferentially included over introns with identical coordinates from shorter isoforms. Duplicate introns may optionally be included in the sequences output file using `-d`.
2. There are a number of criteria by which introns may be omitted from the processed data, depending on run options. These introns will be included in the `bed.iic` and `introns.iic` files (and summarized in `log.iic`), however, tagged with `[o:x]` where `x` is one of the following:
    * `s` | short: Introns that are shorter than 30 nt (by default) cannot be scored, due to length requirements for the scored sub-sequences.
    * `n` | non-canonical: Introns without terminal dinucleotides in the set [`GT-AG`, `GC-AG`, `AT-AC`] are excluded when run with `--no_nc`
    * `a` | ambiguous characters: Introns with ambiguous characters (e.g. 'N') in scoring regions cannot be properly scored and are therefore excluded.
    * `i` | short isoform: If run without `-i`, introns not present in the longest isoform are excluded.
3. Non-canonical introns with very strong U12-like 5′ motifs near their annotated start will have their start and stop coordinates corrected (by equal amounts) to reflect the more U12-like splicing boundaries. These introns are tagged with `[c:x]`, where `x` is the relative coordinate shift applied (the total number of corrected introns is also summarized in `log.iic`).

## **Complete usage info**

```
usage: intronIC [-h] [-g GENOME] [-a ANNOTATION] -n SPECIES_NAME
                [-q SEQUENCE_FILE] [-f {cds,exon}] [-s] [--no_nc] [-i] [-v]
                [-m {matrix file} [{matrix file} ...]]
                [--r12 {reference U12 intron sequences}]
                [--r2 {reference U2 intron sequences}] [--no_plot]
                [--format_info] [-d] [-u] [--na] [-t 0-100] [--ns]
                [--5c start stop] [--3c start stop] [--bpc start stop]
                [-r {five,bp,three} [{five,bp,three} ...]] [--afn]
                [--recursive] [--n_subsample N_SUBSAMPLE]
                [--cv_processes CV_PROCESSES] [-p PROCESSES]
                [--matrix_score_info] [-C HYPERPARAMETER_C]
                [--min_intron_len MIN_INTRON_LEN] [--pseudocount PSEUDOCOUNT]
                [--exons_as_flanks] [-b BED_FILE]

intronIC (intron Interrogator and Classifier) is a script which collects all
of the annotated introns found in a genome/annotation file pair, and produces
a variety of output files (*.iic) which describe the annotated introns and
(optionally) their similarity to known U12 sequences. Without the '-m' flag,
there MUST exist a matrix file in the 'intronIC_data' subdirectory in the same
parent directory as intronIC.py, with filename 'scoring_matrices.fasta.iic'.
In the same data directory, there must also be a pair of sequence files (see
--format_info) with reference intron sequences named '[u2,
u12]_reference_set.introns.iic'

optional arguments:
  -h, --help            show this help message and exit
  -f {cds,exon}, --feature {cds,exon}
                        Specify feature to use to define introns. By default,
                        intronIC will identify all introns uniquely defined by
                        both CDS and exon features. Under the default mode,
                        introns defined by exon features only will be
                        demarcated by an '[e]' tag (default: None)
  -s, --sequences_only  Bypass the scoring system and simply report the intron
                        sequences present in the annotations (default: False)
  --no_nc               Omit introns with non-canonical terminal dinucleoties
                        from scoring (default: False)
  -i, --allow_multiple_isoforms
                        Include non-duplicate introns from isoforms other than
                        the longest in the scored intron set (default: False)
  -v, --allow_intron_overlap
                        Allow introns with boundaries that overlap other
                        introns from higher-priority transcripts (longer
                        coding length, etc.) to be included. This will
                        include, for instance, introns with alternative 5′/3′
                        boundaries (default: False)
  -m {matrix file} [{matrix file} ...], --matrices {matrix file} [{matrix file} ...]
                        One or more matrices to use in place of the defaults.
                        Must follow the formatting described by the
                        --format_info option (default: None)
  --r12 {reference U12 intron sequences}, --reference_u12s {reference U12 intron sequences}
                        introns.iic file with custom reference introns to be
                        used for setting U12 scoring expectation, including
                        flanking regions (default: None)
  --r2 {reference U2 intron sequences}, --reference_u2s {reference U2 intron sequences}
                        introns.iic file with custom reference introns to be
                        used for setting U12 scoring expectation, including
                        flanking regions (default: None)
  --no_plot             Do not output illustrations of intron
                        scores/distributions(plotting requires matplotlib)
                        (default: False)
  --format_info         Print information about the system files required by
                        this script (default: False)
  -d, --include_duplicates
                        Include introns with duplicate coordinates in the
                        intron seqs file (default: False)
  -u, --uninformative_naming
                        Use a simple naming scheme for introns instead of the
                        verbose, metadata-laden default format (default:
                        False)
  --na, --no_abbreviate
                        Use the provided species name in full within the
                        output files (default: False)
  -t 0-100, --threshold 0-100
                        Threshold value of the SVM-calculated probability of
                        being a U12 to determine output statistics (default:
                        90)
  --ns, --no_sequence_output
                        Do not create a file with the full intron sequences of
                        all annotated introns (default: False)
  --5c start stop, --five_score_coords start stop
                        Coordinates describing the 5' sequence to be scored,
                        relative to the 5' splice site (e.g. position 0 is the
                        first base of the intron); half-closed interval
                        [start, stop) (default: (-3, 9))
  --3c start stop, --three_score_coords start stop
                        Coordinates describing the 3' sequence to be scored,
                        relative to the 3' splice site (e.g. position -1 is
                        the last base of the intron); half-closed interval
                        (start, stop] (default: (-10, 4))
  --bpc start stop, --branch_point_coords start stop
                        Coordinates describing the region to search for branch
                        point sequences, relative to the 3' splice site (e.g.
                        position -1 is the last base of the intron); half-
                        closed interval [start, stop). (default: (-55, -5))
  -r {five,bp,three} [{five,bp,three} ...], --scoring_regions {five,bp,three} [{five,bp,three} ...]
                        Intron sequence regions to include in intron score
                        calculations. (default: ('five', 'bp'))
  --afn, --abbreviate_filenames
                        Use abbreviated species name when creating output
                        filenames. (default: False)
  --recursive           Generate new scoring matrices and training data using
                        confident U12s from the first scoring pass. This
                        option may produce better results in species distantly
                        related to the species upon which the training
                        data/matrices are based, though beware accidental
                        training on false positives. Recommended only in cases
                        where clear separation between types is seen with
                        default data. (default: False)
  --n_subsample N_SUBSAMPLE
                        Number of sub-samples to use to generate SVM
                        classifiers; 0 uses the entire training set and should
                        provide the best results; otherwise, higher values
                        will better approximate the entire set at the expense
                        of speed. (default: 0)
  --cv_processes CV_PROCESSES
                        Number of parallel processes to use during cross-
                        validation (default: None)
  -p PROCESSES, --processes PROCESSES
                        Number of parallel processes to use for scoring (and
                        cross-validation, unless --cv_processes is also set)
                        (default: 1)
  --matrix_score_info   Produce additional per-matrix raw score information
                        for each intron (default: False)
  -C HYPERPARAMETER_C, --hyperparameter_C HYPERPARAMETER_C
                        Provide the value for hyperparameter C directly
                        (bypasses optimized parameter search) (default: None)
  --min_intron_len MIN_INTRON_LEN
                        Minimum intron length to consider for scoring
                        (default: 30)
  --pseudocount PSEUDOCOUNT
                        Pseudocount value to add to each matrix value to avoid
                        0-div errors (default: 0.0001)
  --exons_as_flanks     Use entire up/downstream exonic sequence as flank
                        sequence in output (default: False)
  -b BED_FILE, --bed BED_FILE
                        Supply intron coordinates in BED format (default:
                        None)

required arguments (-g, -a | -q):
  -g GENOME, --genome GENOME
                        Genome file in FASTA format (gzip compatible)
                        (default: None)
  -a ANNOTATION, --annotation ANNOTATION
                        Annotation file in gff/gff3/gtf format (gzip
                        compatible) (default: None)
  -n SPECIES_NAME, --species_name SPECIES_NAME
                        Binomial species name, used in output file and intron
                        label formatting. It is recommended to include at
                        least the first letter of the species, and the full
                        genus name since intronIC (by default) abbreviates the
                        provided name in its output (e.g. Homo_sapiens -->
                        HomSap) (default: None)
  -q SEQUENCE_FILE, --sequence_file SEQUENCE_FILE
                        Provide intron sequences directly, rather than using a
                        genome/annotation combination. Must follow the
                        introns.iic format (see README for description)
                        (default: None)
```

## __Example usage__

A sample genome and annotation file can be found in `test_data`. To run `intronIC` on the sample data as organized in the repository simply do the following from within the `test_data` directory:

`$ python3 ../intronIC -g Homo_sapiens.Chr19.Ensembl_91.fa.gz -a Homo_sapiens.Chr19.Ensembl_91.gff3.gz -n homo_sapiens`

Information about the run will be printed to the screen; this same information (plus some additional details) can be found in the `log.iic` file:

```
python3 ../intronIC -g Homo_sapiens.Chr19.Ensembl_91.fa.gz -a Homo_sapiens.Chr19.Ensembl_91.gff3.gz -n homo_sapiens
[#] Starting run on [homo_sapiens (HomSap)]
[#] Run command: [/home/glarue/Documents/Coding/Python/Research/intronIC/intronIC -g /home/glarue/Documents/Coding/Python/Research/intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz -a /home/glarue/Documents/Coding/Python/Research/intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz -n homo_sapiens]
[#] Using [cds,exon] features to define introns
[#] [58933] introns found in [Homo_sapiens.Chr19.Ensembl_91.gff3.gz]
[#] [38681] introns with duplicate coordinates excluded
[#] [8178] introns omitted from scoring based on the following criteria:
[#] * short (<30 nt): 66
[#] * ambiguous nucleotides in scoring regions: 0
[#] * non-canonical boundaries: 0
[#] * overlapping coordinates: 0
[#] * not in longest isoform: 8112
[#] Most common non-canonical splice sites:
[#] * AT-AG (17/328, 5.18%)
[#] * GT-TG (12/328, 3.66%)
[#] * GG-AG (12/328, 3.66%)
[#] * GA-AG (11/328, 3.35%)
[#] * AG-AG (10/328, 3.05%)
[#] [12] ([3] unique, [9] redundant) putatively misannotated U12 introns corrected in [homo_sapiens.annotation.iic]
[#] [12074] introns included in scoring analysis
[#] [11272] introns used to build U2 branch point matrix (5'SS in bottom [95]th percentile)
[#] Scoring introns using the following regions: [five, bp]
[#] Raw scores calculated for [20689] U2 and [387] U12 reference introns
[#] Raw scores calculated for [12074] experimental introns
[#] Non-redundant training sets: [20556] U2, [387] U12
[#] Training SVM using reference data
Starting optimization round 1/5
Starting optimization round 2/5
Starting optimization round 3/5
Starting optimization round 4/5
Starting optimization round 5/5
[#] Range for 'C' after [5] rounds of optimization: [976.5411685881514]-[976.5419176464368]
[#] Set classifier value for 'C': [976.5415431172212]
[#] Training classifier with optimized hyperparameters
[#] Average classifier performance on training data:
	F1	[1.0]
	P-R AUC	[1.0]
[#] Classifier performance details:
              precision    recall  f1-score   support

          U2       1.00      1.00      1.00      4112
         U12       1.00      1.00      1.00        77

    accuracy                           1.00      4189
   macro avg       1.00      1.00      1.00      4189
weighted avg       1.00      1.00      1.00      4189

[#] [1] putative U12 scores were not robust to boundary switching
[#] [10] putative AT-AC U12 introns found.
[#] [31] putative U12 introns found with scores > [90]%
[#] Adding scores to intron sequences file
[#] Generating figures
[#] Run finished in [7.161 minutes]
```

If only the intron sequences are desired, scoring can be bypassed using the `-s` flag which will significantly reduce the processing time and produce only a subset of the output files:

```
➜ python3 ../intronIC -g /home/glarue/Documents/Coding/Python/Research/intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz -a /home/glarue/Documents/Coding/Python/Research/intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz -n homo_sapiens -s
[#] Starting run on [homo_sapiens (HomSap)]
[#] Run command: [/home/glarue/Documents/Coding/Python/Research/intronIC/intronIC -g /home/glarue/Documents/Coding/Python/Research/intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz -a /home/glarue/Documents/Coding/Python/Research/intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz -n homo_sapiens -s]
[#] Using [cds,exon] features to define introns
[#] [58933] introns found in [/home/glarue/Documents/Coding/Python/Research/intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz]
[#] [38681] introns with duplicate coordinates excluded
[#] [20252] intron sequences written to [homo_sapiens.introns.iic]
[#] Run finished in [27.41 seconds]
```

Many additional options exist for a variety of use cases. Run `intronIC --help` for additional details.

## __[background]__