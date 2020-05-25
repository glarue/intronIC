![intronIC_logo](https://user-images.githubusercontent.com/6827531/82829967-62872480-9e69-11ea-94e9-fa7306c7df1b.png)

# (intron <ins>I</ins>nterrogator and <ins>C</ins>lassifier)

See the [[Quick start|Quick-start]] section and the rest of the wiki for installation and other instructions.

## About

`intronIC` was written to provide a customizable, open-source method for identifying minor (U12-type) spliceosomal introns from annotated intron sequences. Minor introns usually represent at most ~0.5% of a given genome's intron complement, but contain distinct splicing motifs which make them amenable to bioinformatic identification.

Earlier minor intron resources (U12DB, SpliceRack, ERISdb, etc.), while hugely important to the field, are by design static. As such, these databases fail to reflect the dramatic increase in available genome sequences and annotation quality of the last decade.

In addition, other published identification methods employ a certain amount of heuristic fuzziness to define their U12-type scoring system. `intronIC` relegates this decision to the well-established support-vector machine (SVM) classification approach, which produces an easily-interpretable "probability of being U12-type" score for each intron.

Furthermore, `intronIC` provides researchers the opportunity to tailor the underlying training data/position-weight matrices, should they have species-specific data that they can take advantage of.

Finally, `intronIC` performs a fair amount of bookkeping during the intron collection process, resulting in (potentially) useful metadata about each intron including parent gene/transcript, ordinal index and phase, information for which there is no other striaghtforward pipeline (as far as I'm aware) to acquire.

## Cite

If you find this tool useful, please cite http://dx.doi.org/10.1093/nar/gkaa464
