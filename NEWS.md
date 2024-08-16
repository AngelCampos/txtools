# txtools 0.1.5

* Added the `loadSecondaryAligns` argument to `tx_load_bam()` so that secondary 
alignments can be easily discarded by setting it to FALSE. The same effect can be
accomplished by inputing a scanFlag built with the `Rsamtools::scanBamFlag()` 
helper, but this way it is a more streamlined feature.

# txtools 0.1.4

* Added tx_combineTxReads() to combine reads processed by tx_reads().
* Setup pkgdown [website](https://angelcampos.github.io/txtools).
* Bug fix: Loading Bed-6 gene annotations broke due to a dependency update.

# txtools 0.1.2

* Bug fix: stringr::str_sub() update broke functions.

# txtools 0.1.1

* Bug fix: Error processing single-end reads

# txtools 0.1.0

* Major: 
    * Added ignore.strand functionality to load reads regardless if the RNA-seq 
    library is strand-aware.
    * New plotting function tx_plot_numeric(): Plots numeric variables along
    the transcriptomic space.
    * tx_metagene_*() functions have a normalization argument which yields the
    are under the curve to be proportional to 1*(transcriptomic window length).
* Plus minor fixes.

# txtools 0.0.8

* Fixed loading paired end BAM files when strand sign is defined by read2 
("last").
* Bugs fixed
* Minor improvements:
    * Loading genome name as first word separated by spaces
    * Edited main vignette tutorial "txtools"
    * Added bam2txDT.R script in inst/ dir
* Patch 0.0.8.1:
    * Added "splicingSite" option to tx_get_metageAtCDS() and 
    tx_plot_metageneAtCDS().

# txtools 0.0.7

* General Update
* Added vignette "txtools_user_guide"
* Changed tx_add_diffNucToRef() for tx_add_misincCount()
* Changed tx_add_diffNucToRefRatio() for tx_add_misincRate()
* Fix minor bugs
* **0.0.7.4 patch**: Fixed a bug that incorrectly included alignments with N 
operations in their CIGAR strings (gaps), that didn't match the gene structure.
* Removed plotly conversion in tx_plot_*() functions. As can be easily done by
the user with plotly::ggplotly()
* **0.0.7.8 patch**: Fixed tx_orderDT() so that order is not affected by levels 
in the 'gene' factor; accordingly fix tx_unifyTxDTL() so DTs have the same 
txcoor order. Changing name of function tx_add_misincorpRateNucSpec()
in favor of tx_add_misincRateNucSpec(), so that it conforms to other function 
names, former can be used as alias but has been deprecated.

# txtools 0.0.6

* Fixed a bug detected when loading BAM files with mappings with empty 
sequence field
* 0.0.6.2 patch: Fixed bug related with spurious paired-end alignments, R1 and
R2 order is wrong.

# txtools 0.0.5

* Update Bioconductor remotes.
* Adding tx_core_mc() which is deprecated as alias of tx_reads().
* Fix bug in tx_complete_DT() which added refSeq before necessary.

# txtools 0.0.4

* Uni and Multi-core functions merged into unique functions 
capable of using multi-cores with the parallel package.

# txtools 0.0.3

* Bug fixed when using BED6 annotation

# txtools 0.0.2

* Updating plotting function tx_plot_staEndCov() to remove 
coverage counts.
* Other minor changes

# txtools 0.0.1 

* Initial draft of package
