# Updates during development

## v0.0.6

* Fixed a bug detected when loading BAM files with mappings with empty 
sequence field
* 0.0.6.2 patch: Fixed bug related with spurious paired-end alignments, R1 and
R2 order is wrong.

## v0.0.5

* Update Bioconductor remotes.
* Adding tx_core_mc() which is deprecated as alias of tx_reads().
* Fix bug in tx_complete_DT() which added refSeq before necessary.

## v0.0.4

* Uni and Multi-core functions merged into one function
capable of running in one or multi-cores with the parallel package.

## v0.0.3

* Bug fixed when using BED6 annotation

## v0.0.2

* Updating plotting function tx_plot_staEndCov() to remove 
coverage counts.
* Other minor changes

## v0.0.1 

* Initial draft of package
