# Test extendGOAnnosWithParentsRcpp
print("Testing extendGOAnnosWithParentsRcpp(...)")
goa <- read.table( stringsAsFactors=FALSE, text=
"GO:0008734 IMP P10902
GO:0034628 EXP P10902
GO:0034628 IDA P66699" )
go.prnts <- read.table( stringsAsFactors=FALSE, header=TRUE, text=
'"id" "name" "term_type" "acc" "is_obsolete" "is_root" "is_relation" "relation_distance" "child_acc"
7226 "L-aspartate oxidase activity" "molecular_function" "GO:0008734" 0 0 0 7 "GO:0008734"
11107 "oxidoreductase activity, acting on the CH-NH2 group of donors" "molecular_function" "GO:0016638" 0 0 0 4 "GO:0008734"
1227 "L-amino-acid oxidase activity" "molecular_function" "GO:0001716" 0 0 0 6 "GO:0008734"
18612 "de novo NAD biosynthetic process from aspartate" "biological_process" "GO:0034628" 0 0 0 9 "GO:0034628"
23453 "single-organism cellular process" "biological_process" "GO:0044763" 0 0 0 3 "GO:0034628"
24954 "nicotinamide nucleotide metabolic process" "biological_process" "GO:0046496" 0 0 0 7 "GO:0034628"' )

res.extendGOAnnosWithParentsRcpp <- extendGOAnnosWithParentsRcpp( goa, go.prnts )
exp.extendGOAnnosWithParentsRcpp <- read.table( stringsAsFactors=FALSE, header=TRUE, text=
"acc ec prot.acc term_type
GO:0008734 IMP P10902 molecular_function
GO:0016638 IMP P10902 molecular_function
GO:0001716 IMP P10902 molecular_function
GO:0034628 EXP P10902 biological_process
GO:0044763 EXP P10902 biological_process
GO:0046496 EXP P10902 biological_process
GO:0034628 IDA P66699 biological_process
GO:0044763 IDA P66699 biological_process
GO:0046496 IDA P66699 biological_process" )
checkEquals( res.extendGOAnnosWithParentsRcpp, exp.extendGOAnnosWithParentsRcpp ) 
