## v0.0.7
* General Update
* Added vignette "txtools_user_guide"
* Changed tx_add_diffNucToRef() for tx_add_misincCount()
* Changed tx_add_diffNucToRefRatio() for tx_add_misincRate()
* Fix minor bugs
 
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

* Uni and Multi-core functions merged into unique functions 
capable of using multi-cores with the parallel package.

## v0.0.3

* Bug fixed when using BED6 annotation

## v0.0.2

* Updating plotting function tx_plot_staEndCov() to remove 
coverage counts.
* Other minor changes

## v0.0.1 

* Initial draft of package