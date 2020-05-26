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
2. Intron sequences file (`-q`) (see [[Configuration]] for formatting information, which matches the reference sequence format)

By default, `intronIC` **includes non-canonical introns**, and considers **only the longest isoform of each gene**. Helpful arguments may include:

* `-p` | parallel processes, which can significantly reduce runtime

* `--no_nc` | exclude introns with non-canonical (non-GT-AG, GC-AG, AT-AC) boundaries

* `-i` | include introns from multiple isoforms of the same gene (default: longest isoform only)

* `-v` | include introns with overlapping boundaries (e.g. alt. 5'/3' boundaries) across multiple isoforms (if `-i`)

## Running on test dataset

To test the installation, change to the `test_data` subdirectory, which contains Ensembl annotations for chr 19 of the human genome.

### Classify annotated introns

* `../intronIC -g Homo_sapiens.Chr19.Ensembl_91.fa.gz -a Homo_sapiens.Chr19.Ensembl_91.gff3.gz -n homo_sapiens`

### Extract all annotated intron sequences

If you just want to retrieve all annotated intron sequences, add the `-s` flag:

* `../intronIC -g Homo_sapiens.Chr19.Ensembl_91.fa.gz -a Homo_sapiens.Chr19.Ensembl_91.gff3.gz -n homo_sapiens -s`

See the rest of the Wiki for more extensive details about [[output files|Output-files]], [[usage info|Usage-info]], etc.

## Resource usage

For genomes with a large number of annotated introns, memory usage can be on the order of gigabytes. This should rarely be a problem even for most modern personal computers, however. For reference, the Ensembl 95 release of the human genome requires ~5 GB of memory.

For many non-model genomes, intronIC should run fairly quickly (e.g. tens of minutes). For human and other very well annotated genomes, runtime may be longer (the human Ensembl 95 release takes ~20-35 minutes in testing); run time scales relatively linearly with the total number of annotated introns, and can be improved by using parallel processes via `-p`.

See the rest of the wiki for more detailed instructions.

## About

`intronIC` was written to provide a customizable, open-source method for identifying minor (U12-type) spliceosomal introns from annotated intron sequences. Minor introns usually represent at most ~0.5% of a given genome's intron complement, but contain distinct splicing motifs which make them amenable to bioinformatic identification.

Earlier minor intron resources (U12DB, SpliceRack, ERISdb, etc.), while hugely important to the field, are by design static. As such, these databases fail to reflect the dramatic increase in available genome sequences and annotation quality of the last decade.

In addition, other published identification methods employ a certain amount of heuristic fuzziness to define their U12-type scoring system. `intronIC` relegates this decision to the well-established support-vector machine (SVM) classification approach, which produces an easily-interpretable "probability of being U12-type" score for each intron.

Furthermore, `intronIC` provides researchers the opportunity to tailor the underlying training data/position-weight matrices, should they have species-specific data that they can take advantage of.

Finally, `intronIC` performs a fair amount of bookkeping during the intron collection process, resulting in (potentially) useful metadata about each intron including parent gene/transcript, ordinal index and phase, information for which there is no other striaghtforward pipeline (as far as I'm aware) through which to easily acquire.

## Cite

If you find this tool useful, please cite http://dx.doi.org/10.1093/nar/gkaa464
