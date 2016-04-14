Detection of Duplicates on Patterned Flowcells
==============================================

THE PLN: Directly detect local duplicates in patterned flow cells without alignment or demultiplexing.

Overview
--------

We sample a number targets on each tile of an Illumina Hiseq 4000 (or Hiseq X) flowcell to scan for duplicates.  A target is a specific well and the surrounding wells out to a pre-defined distance (say, five layers out in the honeycomb).  For each target we read a substring of the sequence for every well directly from the raw BCL files.  Then we count the duplicates found by looking for wells that match the centre sequence.

Usage
-----

```prepare_cluster_indexes.py``` will come up with a list of cluster locations (targets) to be sampled, and work out the co-ordinates of all the surrounding wells.  It parses the standard .locs file found in the Data directory for every Illumina run.  Note that the layout of wells is specific to the generation of flowcell rather than being specific to the machine, so watch out if you are planning to use the same locations file for scanning multiple flowcells - check that the .locs files are indeed the same.

```count_well_duplicates.py``` will read the data from your BCL files and output duplication stats.  It needs to be supplied with a run to be analysed and also a targets file produced with the ```prepare_cluster_indexes.py``` script.

BCL Direct Reader
-----------------

The file ```bcl_direct_reader.py``` contains pure Python code for retrieving sequence from raw BCL files.  For our purposes there is little to be gained from porting this to C as most of the time is spent Gunzipping the data.

Run ```pydoc ./bcl_direct_reader.py``` for more info.

Warning
-------

This code is not yet well tested, specifically it has not been tested on Hiseq X data, though it is designed to work on those flowcells.
