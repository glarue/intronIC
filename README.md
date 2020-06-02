![intronIC_logo](https://user-images.githubusercontent.com/6827531/82829967-62872480-9e69-11ea-94e9-fa7306c7df1b.png)

# (intron <ins>I</ins>nterrogator and <ins>C</ins>lassifier)

`intronIC` is a program that can be used to classify intron sequences as minor (U12-type) or major (U2-type), using a genome and annotation or the sequences themselves.

First, clone the repo to your local machine:

```console
$ git clone https://github.com/glarue/intronIC.git
$ cd intronIC
```

You may also wish to add `intronIC` to your system path (how you do this is platform-dependent).

See the [wiki](https://github.com/glarue/intronIC/wiki) for more detail information about configuration/run options.

## Cite

If you find this tool useful, please cite https://doi.org/10.1093/nar/gkaa464

## Dependencies

* [Python >=3.3](https://www.python.org/downloads/)
* [numpy & scipy](https://www.scipy.org/scipylib/download.html)
* [scikit-learn >=0.20.1](http://scikit-learn.org/stable/index.html)
* [biogl](https://github.com/glarue/biogl)
* [matplotlib](https://matplotlib.org/) (optional, required for plotting)

To install dependencies using `pip`, do

`python3 -m pip install numpy scipy matplotlib scikit-learn biogl`

`intronIC` was built and tested on Linux, but should run on Windows or Mac OSes without much trouble (I say that now...).

## Useful arguments

The required arguments for any classification run include a name (`-n`), along with:

1. Genome (`-g`) and annotation/BED (`-a`, `-b`) files or,
2. Intron sequences file (`-q`) (see [Training-data-and-PWMs](https://github.com/glarue/intronIC/wiki/Training-data-and-PWMs) for formatting information, which matches the reference sequence format)

By default, `intronIC` **includes non-canonical introns**, and considers **only the longest isoform of each gene**. Helpful arguments may include:

* `-p`  parallel processes, which can significantly reduce runtime

* `-f cds`  use only `CDS` features to identify introns (by default, uses both `CDS` and `exon` features)

* `--no_nc` exclude introns with non-canonical (non-`GT-AG`/`GC-AG`/`AT-AC`) boundaries

* `-i`  include introns from multiple isoforms of the same gene (default: longest isoform only)

## Running on test data

To test the installation, change to the `test_data` subdirectory, which contains Ensembl annotations for chromosome 19 of the human genome.

### Classify annotated introns

* `../intronIC -g Homo_sapiens.Chr19.Ensembl_91.fa.gz -a Homo_sapiens.Chr19.Ensembl_91.gff3.gz -n homo_sapiens`

### Extract all annotated intron sequences

If you just want to retrieve all annotated intron sequences, add the `-s` flag:

* `../intronIC -g Homo_sapiens.Chr19.Ensembl_91.fa.gz -a Homo_sapiens.Chr19.Ensembl_91.gff3.gz -n homo_sapiens -s`

See the rest of the [wiki](https://github.com/glarue/intronIC/wiki) for more details about [output files](https://github.com/glarue/intronIC/wiki/Output-files), etc.

## A note on the `-n` (name) argument

By default, `intronIC` expects names in binomial (genus, species) form separated by a non-alphanumeric character, e.g. 'homo_sapiens', 'homo.sapiens', etc. `intronIC` then formats that name internally into a tag that it uses to label all output intron IDs, ignoring anything past the second non-alphanumeric character.

Output *files*, on the other hand, are named using the full name supplied via `-n`. If you'd prefer to have it leave whatever argument you supply to `-n` unmodified, use the `--na` flag.

If you are running multiple versions of the same species and would like to keep the same species abbreviations in the output intron data, simply add a tag to the end of the name, e.g. "homo_sapiens.v2"; the tags within files will be consistent ("HomSap"), but the file names across runs will be distinct.

## Resource usage

For genomes with a large number of annotated introns, memory usage can be on the order of gigabytes. This should rarely be a problem even for most modern personal computers, however. For reference, the Ensembl 95 release of the human genome requires ~5 GB of memory.

For many non-model genomes, intronIC should run fairly quickly (e.g. tens of minutes). For human and other very well annotated genomes, runtime may be longer (the human Ensembl 95 release takes ~20-35 minutes in testing); run time scales relatively linearly with the total number of annotated introns, and can be improved by using parallel processes via `-p`.

See the rest of the wiki for more detailed instructions.

## About

`intronIC` was written to provide a customizable, open-source method for identifying minor (U12-type) spliceosomal introns from annotated intron sequences. Minor introns usually represent ~0.5% (at most) of a given genome's introns, and contain distinct splicing motifs which make them amenable to bioinformatic identification.

Earlier minor intron resources (U12DB, SpliceRack, ERISdb, etc.), while important contributions to the field, are static by design. As such, these databases fail to reflect the dramatic increase in available genome sequences and annotation quality of the last decade.

In addition, other published identification methods employ a certain amount of heuristic fuzziness in defining the classification criteria of their U12-type scoring systems (i.e how "U12-like" does an intron need to look before being called a U12-type intron). `intronIC` relegates this decision to the well-established support-vector machine (SVM) classification method, which produces an easy-to-interpret "probability of being U12-type" score for each intron.

Furthermore, `intronIC` provides researchers the opportunity to tailor the underlying training data/position-weight matrices, should they have species-specific data to take advantage of.

Finally, `intronIC` performs a fair amount of bookkeping during the intron collection process, resulting in (potentially) useful metadata about each intron including parent gene/transcript, ordinal index and phase, information which (as far as I'm aware) is otherwise non-trivial to acquire.
