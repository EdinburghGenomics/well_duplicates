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

Hiseq 4000 s.locs coordinates:

| pixel_distance | x | y | rownr | row_distance | max_pix_dist | pixel_distance | x | y | rownr | row_distance | max_pix_dist | pixel_distance | x | y | rownr | row_distance | max_pix_dist |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 9536 | 6589 | 500000 | 0 |  | 0 | 18101 | 12181 | 1000000 | 0 |  | 0 | 21268 | 998 | 1000 | 0 |
| 20.5913 | 9546 | 6607 | 501571 | -1571 |  | 19.7231 | 18111 | 12198 | 1001571 | -1571 |  | 19.7231 | 21278 | 1015 | 2571 | -1571 |
| 20.5913 | 9526 | 6607 | 501570 | -1570 |  | 19.7231 | 18091 | 12198 | 1001570 | -1570 |  | 20.2485 | 21257 | 1015 | 2570 | -1570 |
| 21 | 9557 | 6589 | 500001 | -1 |  | 21 | 18122 | 12181 | 1000001 | -1 |  | 20 | 21288 | 998 | 1001 | -1 |
| 20 | 9516 | 6589 | 499999 | 1 |  | 20 | 18081 | 12181 | 999999 | 1 |  | 21 | 21247 | 998 | 999 | 1 |
| 19.7231 | 9546 | 6572 | 498429 | 1571 |  | 20.5913 | 18111 | 12163 | 998429 | 1571 |  |  |  |  |  |  |
| 19.7231 | 9526 | 6572 | 498428 | 1572 | 22 | 20.5913 | 18091 | 12163 | 998428 | 1572 |  |  |  |  |  |  |
| 40.8167 | 9557 | 6624 | 503143 | -3143 |  | 40.8167 | 18122 | 12216 | 1003143 | -3143 |  | 40.3113 | 21288 | 1033 | 4143 | -3143 |
| 35 | 9536 | 6624 | 503142 | -3142 |  | 35 | 18101 | 12216 | 1003142 | -3142 |  | 35 | 21268 | 1033 | 4142 | -3142 |
| 40.3113 | 9516 | 6624 | 503141 | -3141 |  | 40.3113 | 18081 | 12216 | 1003141 | -3141 |  | 40.8167 | 21247 | 1033 | 4141 | -3141 |
| 35.8469 | 9567 | 6607 | 501572 | -1572 |  | 35.3553 | 18132 | 12198 | 1001572 | -1572 |  | 34.4819 | 21298 | 1015 | 2572 | -1572 |
| 34.9857 | 9506 | 6607 | 501569 | -1569 |  | 34.4819 | 18071 | 12198 | 1001569 | -1569 |  | 35.3553 | 21237 | 1015 | 2569 | -1569 |
| 41 | 9577 | 6589 | 500002 | -2 |  | 41 | 18142 | 12181 | 1000002 | -2 |  | 40 | 21308 | 998 | 1002 | -2 |
| 40 | 9496 | 6589 | 499998 | 2 |  | 40 | 18061 | 12181 | 999998 | 2 |  | 41 | 21227 | 998 | 998 | 2 |
| 35.3553 | 9567 | 6572 | 498430 | 1570 |  | 35.8469 | 18132 | 12163 | 998430 | 1570 |  |  |  |  |  |  |
| 34.4819 | 9506 | 6572 | 498427 | 1573 |  | 34.9857 | 18071 | 12163 | 998427 | 1573 |  |  |  |  |  |  |
| 40.8167 | 9557 | 6554 | 496859 | 3141 |  | 40.8167 | 18122 | 12146 | 996859 | 3141 |  |  |  |  |  |  |
| 35 | 9536 | 6554 | 496858 | 3142 |  | 35 | 18101 | 12146 | 996858 | 3142 |  |  |  |  |  |  |
| 40.3113 | 9516 | 6554 | 496857 | 3143 | 42 | 40.3113 | 18081 | 12146 | 996857 | 3143 |  |  |  |  |  |  |
| 61.4003 | 9567 | 6642 | 504714 | -4714 |  | 61.4003 | 18132 | 12234 | 1004714 | -4714 |  | 60.0333 | 21298 | 1050 | 5714 | -4714 |
| 53.9351 | 9546 | 6642 | 504713 | -4713 |  | 53.9351 | 18111 | 12234 | 1004713 | -4713 |  | 52.9528 | 21278 | 1050 | 5713 | -4713 |
| 53.9351 | 9526 | 6642 | 504712 | -4712 |  | 53.9351 | 18091 | 12234 | 1004712 | -4712 |  | 53.1507 | 21257 | 1050 | 5712 | -4712 |
| 60.9016 | 9506 | 6642 | 504711 | -4711 |  | 60.9016 | 18071 | 12234 | 1004711 | -4711 |  | 60.5392 | 21237 | 1050 | 5711 | -4711 |
| 53.9073 | 9577 | 6624 | 503144 | -3144 |  | 53.9073 | 18142 | 12216 | 1003144 | -3144 |  | 53.1507 | 21308 | 1033 | 4144 | -3144 |
| 53.1507 | 9496 | 6624 | 503140 | -3140 |  | 53.1507 | 18061 | 12216 | 1003140 | -3140 |  | 53.9073 | 21227 | 1033 | 4140 | -3140 |
| 54.0833 | 9587 | 6607 | 501573 | -1573 |  | 53.7587 | 18152 | 12198 | 1001573 | -1573 |  | 52.811 | 21318 | 1015 | 2573 | -1573 |
| 54.0833 | 9485 | 6607 | 501568 | -1568 |  | 52.811 | 18051 | 12198 | 1001568 | -1568 |  | 53.7587 | 21217 | 1015 | 2568 | -1568 |
| 61 | 9597 | 6589 | 500003 | -3 |  | 61 | 18162 | 12181 | 1000003 | -3 |  | 60 | 21328 | 998 | 1003 | -3 |
| 61 | 9475 | 6589 | 499997 | 3 |  | 61 | 18040 | 12181 | 999997 | 3 |  | 61 | 21207 | 998 | 997 | 3 |
| 53.7587 | 9587 | 6572 | 498431 | 1569 |  | 54.0833 | 18152 | 12163 | 998431 | 1569 |  |  |  |  |  |  |
| 53.7587 | 9485 | 6572 | 498426 | 1574 |  | 53.1413 | 18051 | 12163 | 998426 | 1574 |  |  |  |  |  |  |
| 53.9073 | 9577 | 6554 | 496860 | 3140 |  | 53.9073 | 18142 | 12146 | 996860 | 3140 |  |  |  |  |  |  |
| 53.1507 | 9496 | 6554 | 496856 | 3144 |  | 53.1507 | 18061 | 12146 | 996856 | 3144 |  |  |  |  |  |  |
| 61.4003 | 9567 | 6536 | 495288 | 4712 |  | 61.4003 | 18132 | 12128 | 995288 | 4712 |  |  |  |  |  |  |
| 53.9351 | 9546 | 6536 | 495287 | 4713 |  | 53.9351 | 18111 | 12128 | 995287 | 4713 |  |  |  |  |  |  |
| 53.9351 | 9526 | 6536 | 495286 | 4714 |  | 53.9351 | 18091 | 12128 | 995286 | 4714 |  |  |  |  |  |  |
| 60.9016 | 9506 | 6536 | 495285 | 4715 | 62 | 60.9016 | 18071 | 12128 | 995285 | 4715 |  |  |  |  |  |  |

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


