# GRPutils

This repository contains utilities to manipulate Artemis compatible coverage **gr**a**p**hs (denote *.grp format), and a set of commands developed based on this file format. The grp file can be visualized using Artemis directly by clicking on `Graphs => Add User Plot`. The majority of the grp manipulating commands can be found within `grptools`.

* `wig2grp.py`: convert Wiggle files generated by `bam2wig.py` from the RseQC package to grp format.  
* `merge_grp.py`: given a list of grp files from the same genome, this script will average/add-up the coverage at each genome coordinate, then return a new grp file.
* `grptools`: a suite of subcommands for grp file manipulations.
* `tptools`: a suite of subcommands for TSS prediction and classification.
* `tss_analysis_pipline.sh`: an example shell script utilizing `grptools` and `tptools` for grp file processing, TSS prediction and classification. 


