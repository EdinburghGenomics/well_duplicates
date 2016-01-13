THE PLN
=======

Detect well duplicates in patterned flow cells without alignment

Setup stage
-----------

input:

* s.locs file
* sample size N
        
output: file with cluster_index => x.coord, y.coord

    open s.locs
    generate N random numbers between 0..max_num of clusters
    for (1..N)
        find cluster_index entry in s.locs
        convert to x.coord, y.coord
        find adjacent 6,12,18 clusters
        get their cluster_index
        convert to x.coord, y.coord
        add to output

Run stage
---------

on each run in the QC pipeline stage

    get bcl reader
    for each lane
        for each tile
            for each cluster_index in file
                get central read (nuc_seq, pass)
                get 6,12,18 surrounding reads (nuc_seq, pass)
                check_edit_distance (LD<=2 -> duplicate read)
                calculate %duplicates


