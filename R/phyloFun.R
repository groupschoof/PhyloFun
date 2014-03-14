library(gRain)

roundBranchLengths <- function( phyl.tree, digits=2 ) {
  phyl.tree$edge.length <- round( phyl.tree$edge.length, digits=digits )
  phyl.tree
}

# Returns unique node, that has no edge leading to it. Thus, it must be a node
# in the first column whereas mustn't appear in the second column of the
# phylTree's edge matrix.
get.root.node <- function(phylTree) {
  as.character(
    unique(phylTree$edge[!(phylTree$edge[, 1] %in% phylTree$edge[, 2]), ][, 1])[[1]]
    )
}

edgeIndex <- function( desc.node, phyl.tree ) {
  # Returns the index of the tree edge of which the 'desc.node' is the
  # descendant.
  #
  # Args:
  #  desc.node : The index of the descendant node. 
  #  phyl.tree : An object of class phylo representing the phylogenetic tree.
  #
  # Returns: The edge's index, of which 'desc.node' is the descendant.
  #   
  which( phyl.tree$edge[ , 2] == desc.node )
}

parentNode <- function( desc.node, phyl.tree ) {
  # Finds the index of the node being the parent of 'desc.node'.
  #
  # Args:
  #  desc.node : The index of the node for which to find its parents. 
  #  phyl.tree : An object of class phylo representing the phylogenetic tree.
  #
  # Returns: The index of 'desc.node's parent node.
  #   
  phyl.tree$edge[ which( phyl.tree$edge[ , 2] == desc.node ), 1 ][[ 1 ]]
}

cumulativeBranchLengthsToRoot <- function( curr.node, phyl.tree,
  root.node=get.root.node(phyl.tree), curr.dist=0.0 ) {
  # Computes the cumulative branch lengths from curr.node to the root.node of
  # the phylogenetic tree 'phyl.tree'.
  #
  # Args:
  #  curr.node : Phylogenetic node to measure its distance to the root node
  #              for.
  #  phyl.tree : Object of class phylo representing the phylogenetic tree.
  #  root.node : The phylogenetic tree's root node.
  #  curr.dist : As this is recursive function, this is the start value for the
  #              distance measuring.
  #
  # Returns: The sum of branch lengths on the path from curr.node to the tree's
  # root.
  #   
  cum.dist <- curr.dist +
    phyl.tree$edge.length[[ edgeIndex( curr.node, phyl.tree ) ]]
  prnt.node <- parentNode( curr.node, phyl.tree )
  if ( prnt.node != root.node ) {
    cumulativeBranchLengthsToRoot( prnt.node, phyl.tree, root.node,
      cum.dist )
  } else {
    cum.dist
  }
}

maxDistanceToRoot <- function( phyl.tree ) {
  # Computes the maximum phylogenetic distance of any tip to the tree's
  # root node.
  #
  # Args:
  #  phyl.tree : Object of class phylo representing the phylogenetic tree. 
  #
  # Returns: The numeric maximum phylogenetic distance, as cumulative
  # branch lengths on the path from any tree's tip to tree's root.
  #   
  max(
    as.numeric(
      lapply(
        1:length( phyl.tree$tip.label ),
        cumulativeBranchLengthsToRoot, phyl.tree=phyl.tree
      )
    )
  )
}

getDescendantNodes <- function( phylo.tree, node.index ) {
  # Looks up all direct descendants of 'node.index'.
  #
  # Args:
  #  phylo.tree : An object of class phylo as returned by read.tree.
  #  node.index : The integer index of the node to find direct descendants for.
  #               See phylo.tree$edge for more details.
  #
  # Returns: The integer vector of those node indices being direct descendants
  # of argument 'node.index', or NULL if no descendants can be found.
  #   
  desc.nds <- phylo.tree$edge[ phylo.tree$edge[ , 1 ] == node.index, 2 ]
  if ( length( desc.nds ) > 0 ) desc.nds
}

get.node.label <- function(phylTree, node.index) {
  if ( node.index <= length( phylTree$tip.label ) )
    phylTree$tip.label[[ node.index ]]
  else
    as.character( node.index )
}

edge.to.formula <- function(phyloTree, edge.index) {
  prnt <- get.node.label(phyloTree,
    phyloTree$edge[[edge.index, 1]])
  child <- get.node.label(phyloTree,
    phyloTree$edge[[edge.index, 2]])
  eval(bquote(~ .(child) | .(prnt) ))
}

goTypeAnnotationMatrices <- function( annotation.df,
  valid.go.terms=names( GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE ) ) {
  # For each type of Gene Ontology terms, 'molecular_function',
  # 'biological_process', and 'cellular_component', a separate annotation
  # matrix is generated. These matrices hold the set of GO term annotations for
  # each protein as in the column names of argument 'annotation.df'. If a
  # non NULL and non empty valid.go.terms list is provided only GO term
  # annotations appearing in this list will appear.
  #
  # Args:
  #  annotation.df  : A Gene Ontology (GO) term annotation data frame as
  #                   returned by 'extendGOAnnosWithParents( …,
  #                   append.term.type=TRUE )'.
  #  valid.go.terms : If non NULL and non empty only terms appearing in
  #                   this list will be assigned.
  #
  # Returns: A named list with three annotation matrices 'biological_process',
  # 'cellular_component', and 'molecular_function'. Each of these annotation
  # data frames is the subset of argument 'annotation.df' where the fourth
  # column matches the respective GO type. 
  #   
  go.types <- c( 'biological_process', 'cellular_component', 'molecular_function' ) 
  setNames(
    lapply( go.types, function( go.tp ) {
      rows.inds <- annotation.df[ , 4 ] == go.tp &
        if ( ! is.null( valid.go.terms ) ) {
          annotation.df[ , 1 ] %in% valid.go.terms
        } else {
          TRUE
        }
      annotation.df[ which( rows.inds ), ]
    }),
    go.types
  )
}

eliminateUnreachableStates <- function( conditional.probs.tbl ) {
  # gRain throws errors if conditional probability tables have unreachable
  # stated. This function identifies them by looking up indices of those
  # columns summing up to zero and eliminating columns and rows of these
  # indices.
  #
  # Args:
  #  conditional.probs.tbl : A probability matrix as generated by
  #                          conditionalProbabilityTable.
  #
  # Returns: The conditional.probs.tbl without rows and columns of unreachable
  # states.
  #   
  zero.cols <- apply( conditional.probs.tbl, 2, function( x ) sum(x) == 0 )
  if ( any( zero.cols ) ) {
    cpt <- conditional.probs.tbl[ , - which( zero.cols ), drop=F ]
  } else {
    conditional.probs.tbl # unchanged
  }
}

annotationToString <- function( compound.annotation ) {
  # Throughout PhyloFun this is used to generate a unique string from a
  # character vector of protein annotations. - We could of course use the
  # toString method, but this function produces a more readable output when
  # more compound annotations appear as random variable values.
  #
  # Args:
  #  compound.annotation : a character vector, i.e. c( "GO_AGO_B", "GO_C" ).
  #
  # Returns: A character vector of length one, i.e. "GO_A & GO_B & GO_C"
  #   
  paste( compound.annotation, collapse=" & " )
}

annotationToCharacterVector <- function( compound.annotation ) {
  # Deconstructs a compound annotation into a character vector.
  #
  # Args:
  #  compound.annotation : Result of annotationToString
  #
  # Returns: A split character vector. 'GO_A & GO_B' -> c( 'GO_A', 'GO_B' )
  #   
  strsplit( compound.annotation, '\\s*&\\s*' )[[ 1 ]]
}

mostAppropriateAnnotation <- function( named.annotation.score.vector,
  unknown.annot='unknown' ) {
  # Returns the most appropriate annotation, fulfilling the following three
  # criteria:
  # 1. It is one of the highest scoring ( There might be more than one! )
  # 2. Preferably it is not the unknown.annot
  # 3. Preferably it is the longest annotation, that is it has the highest
  #    number of terms annotated.
  #
  # Args:
  #  named.annotation.score.vector : A vector of posterior probabilities whose
  #     names are the actual compound annotations. This vector is the single entry
  #     of a typical PhyloFun result list for a given Gene Ontology term type (
  #     i.e. 'molecular_function' ).  
  #  unknown.annot : The string used to represent the empty annotation
  #     'UNKNOWN'.
  #
  # Returns: The selected annotation as character vector.
  #   
  score.srt <- sort( named.annotation.score.vector )
  max.ind <- length( score.srt )
  max.scr <- score.srt[[ max.ind ]]
  # Check for equally scoring annots, prefer any other than 'unknown' in that
  # case:
  best.annos <- score.srt[ which( score.srt == max.scr ) ]
  selected.ind <- if ( length( best.annos ) > 1 ) {
    no.unknown <- setdiff( names( best.annos ), unknown.annot )
    lns <- sort(
      setNames( as.integer( lapply( no.unknown, nchar ) ), no.unknown )
    )
    names( lns[ length(lns) ] )
  } else
    max.ind
  # Return annotated terms as a character vector:
  annotationToCharacterVector(
    names( score.srt[ selected.ind ] )
  )
}

highScoringAnnotations <- function( phyloFun.annos, query.accession,
  unknown.annot='unknown', exclude.unknown=TRUE ) {
  # Finds those PhyloFun Gene Ontology Annotations whose probability is higher
  # or equal than equal distribution.
  #
  # Args:
  #  phyloFun.annos  : The result of PhyloFun as list of gRain predictions,
  #                    where the list names are the three GO types.
  #  query.accession : The Query Protein's Accession surrounded with double
  #                    quotes.
  #  unknown.annot   : The string used to describe the UNKNOWN annotation.
  #  exclude.unknown : If set to TRUE the unknown annotation will not be
  #                    included in the results.
  #
  # Returns: A named list of character vectors of "high scoring" Gene Ontology
  # Term accessions, in which the names are the three Gene Ontology type.
  #   
  if ( ! is.null( phyloFun.annos ) ) {
    setNames(
      lapply( names( phyloFun.annos ), function( gt ) {
        gt.annos <- phyloFun.annos[[ gt ]]
        rslt <- if ( ! is.null( gt.annos ) ) {
          annos <- gt.annos[[ 'pred' ]][[ query.accession ]]
          co <- 1/length(annos)
          ba <- colnames(annos[ 1, which( annos[1,] >= co ), drop=F ])
          if ( ! is.null( ba ) ) 
            unlist( lapply( ba, annotationToCharacterVector ) )
          else
            ba
        } else {
          NA
        }
        unique(
          ( if ( exclude.unknown ) setdiff( rslt, unknown.annot ) else rslt )
        )
      } ),
      names( phyloFun.annos )
    )
  } else {
    NA
  }
}

predictionsToCharacterVector <- function( go.prediction.list,
  query.protein.accession, unknown.annot='unknown' ) {
  # For each predicted Gene Ontology annotation the most appropriate is
  # selected and returned as a character vector of the anntotated GO terms. 
  #
  # Args:
  #  go.prediction.list : PhyloFun result list with entries for each GO type
  #                       'biological_process', 'cellular_component', and
  #                       'molecular_function'
  #  query.protein.accession : The quoted query protein accession, as returned
  #                       by calling function surroundEachWithQuotes(
  #                       original.protein.accession ). If this protein
  #                       accession does not start and end with '"' it will be
  #                       passed to function surroundEachWithQuotes to quote
  #                       it.
  #  unknown.annot      : The string used as the empty annotation.
  #
  # Returns: A character vector of annotated GO terms that matched the most
  # appropriate requirements. The string used to express PhyloFun could not
  # annotate this query protein ( unknown.annot ) is returned, if and only if
  # not a single term could be annotated.
  #   
  quoted.acc <- if ( grepl( '^".+"$', query.protein.accession ) )
    query.protein.accession
  else
    surroundEachWithQuotes( query.protein.accession )
  annots <- setdiff(
    unique(
      unlist(
        lapply( names( go.prediction.list ), function( go.type ) {
          agps <- go.prediction.list[[ go.type ]]
          if ( ! is.null( agps ) ) {
            preds <- agps$pred[[ quoted.acc ]][ 1, ]
            if ( ! is.null( preds ) && length( preds ) > 0 ) {
              mostAppropriateAnnotation( preds )
            }
          }
        })
      )
    ),
    unknown.annot
  )
  if ( is.null( annots ) || length( annots ) == 0 ) unknown.annot else annots 
}

goTermPredictionTable <- function( go.prediction.list, query.protein.accession,
  unknown.annot='unknown', go.con=connectToGeneOntology(),
  close.db.connection=TRUE ) {
  # Finds for each category of the Gene Ontology 'biological_process',
  # 'cellular_component', and 'molecular_function' the most appropriate
  # annotation in the PhyloFun result 'go.prediction.list'. Most appropriate is
  # defined in function mostAppropriateAnnotation(…). The resulting GO term
  # accessions are used to query the Gene Ontology database and all available
  # information is compiled into a table. See function
  # goTermsForAccessionWithLevel(…) for more details.
  #
  # Args:
  #  go.prediction.list : PhyloFun result list with entries for each GO type
  #                       'biological_process', 'cellular_component', and
  #                       'molecular_function'
  #  query.protein.accession : The query protein's accession.
  #  unknown.annot      : The string used as the empty annotation.
  #  go.con : The database connection to the Gene Ontology as returned by
  #                       function connectToGeneOntology.
  #  close.db.connection : If set to TRUE, the database connection go.con will
  #                       be closed automatically.
  #
  # Returns: A matrix with columns "id", "name", "term_type", "acc",
  # "is_obsolete", "is_root", "is_relation", "relation_distance" in which each
  # annotated Gene Ontology (GO) term fills up a row. If and only if no GO
  # terms are annotated the UNKNOWN matrix is returned.
  #   
  go.accs <- predictionsToCharacterVector( go.prediction.list,
    query.protein.accession, unknown.annot=unknown.annot )
  go.tbl <- if ( length( go.accs ) == 1 && go.accs == unknown.annot ) 
    matrix(
      c( NA, 'unknown', NA, NA, NA, NA, NA, NA ), byrow=T, ncol=8,
      dimnames=list( c(), c( "id", "name", "term_type", "acc", "is_obsolete",
        "is_root", "is_relation", "relation_distance" ) )
    )
  else 
    goTermsForAccessionWithLevel( go.accs, con=go.con )
  # Close DB connection?
  if ( close.db.connection ) dbDisconnect( go.con )
  # return
  go.tbl
}

goAnnotationSpaceList <- function( go.type.annotation.matrices,
  unknown.annot='unknown' ) {
  # Generates the spaces of unique annotations for each GO term type
  # 'biological_process', 'cellular_component', and 'molecular_function'.
  #
  # Args:
  #  go.type.annotation.matrices : list of GO annotation matrices as returned
  #                                by 'goTypeAnnotationMatrices'.
  #  unknown.annot               : Set to NULL, if the UNKNOWN annotation
  #                                should be excluded from each annotation
  #                                space. Set to any other convenient
  #                                character, otherwise.
  #
  # Returns: A list of lists. Names are the GO term types and contained lists
  # are the unique GO term annotations as found in the argument
  # 'go.type.annotation.matrices'.
  #   
  setNames(
    lapply( names( go.type.annotation.matrices ), function( go.type ) {
      as <- annotationSpace( go.type.annotation.matrices[[ go.type ]] )
      if ( ! is.null( unknown.annot ) )
        as <- c( as, unknown.annot )
      as
    }),
    names( go.type.annotation.matrices )
 )
}

annotationSpace <- function( annotation.df ) {
  # Generates the set of unique anntations found in the column 1 of argument
  # 'annotation.df'.
  #
  # Args:
  #  annotation.df : The matrix of protein annotations as returned by
  #                  retrieveGOAnnotations(…)
  #
  # Returns: A character vector of unique annotations.
  #   
  if ( nrow( anno.df ) > 0 ) {
    unique( annotation.df[ , 1 ] )
  } else {
    ''
  }
}

bayesNode <- function( phylo.tree, annotation.space, node.index, cond.prob.tbls,
  mutation.probability.tables.list=GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE,
  unknown.annot='unknown' ) {
  # Compiles the conditional probability tables for argument phylogenetic tree
  # node 'node.index'. 
  #
  # Args:
  #  phylo.tree                       : An object of class phylo representing
  #                                     the phylogenetic tree to translate into
  #                                     a Bayesian network.
  #  annotation.space                 : The set of annotations as found in the
  #                                     parent node of siblings 'node.indices'.
  #  node.index                       : The index of the node to compile the
  #                                     CPT for.
  #  cond.prob.tbls                   : The named list of annotation mutation
  #                                     conditional probability tables as
  #                                     returned by function
  #                                     conditionalProbabilityTables. 
  #  mutation.probability.tables.list : The list of mutation probabilities for
  #                                     measured branch lengths.
  #  unknown.annot                    : The label of the unknown annotation.
  #                                     Default is 'unknown'.
  #
  # Returns: The compiled conditional probability table (CPT) for argument node
  # 'node.index' as a list.
  #   
  edge.ind <- as.integer( which( phylo.tree$edge[ , 2 ] == node.index ) )
  is.root.node <- node.index == get.root.node( phylo.tree )
  annotations <- as.character( lapply( annotation.space, annotationToString ) )
  # conditional probability table:
  cond.prob.mtrx <- if ( is.root.node ) {
    matrix( 1.0, ncol=1, nrow=length( annotation.space ), 
      dimnames=list( annotations, c() )
    )
  } else {
    cond.prob.tbls[[ as.character( phylo.tree$edge.length[[ edge.ind ]] ) ]]
  }
  # Generate the conditional probability table for current Bayesian node:
  edge.formula <- if ( is.root.node ) {
    eval( bquote( ~ .( get.root.node( phylo.tree ) ) ) )
  } else {
    edge.to.formula( phylo.tree, edge.ind )
  }
  list(
    cptable(
      edge.formula,
      values=cond.prob.mtrx,
      levels=annotations
    )
  )
}

bayesNodes <- function( phylo.tree, annotation.space,
  mutation.probability.tables.list=GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE,
  unknown.annot='unknown' ) {
  # Computes the conditional probability tables (CPT) for each node in the
  # argument 'phylo.tree' using function 'bayesNode'.
  #
  # Args:
  #  phylo.tree                       : An object of class phylo as returned by
  #                                     function 'read.tree' of library ape.
  #  annotation.space                 : The set of annotations as found in the
  #                                     parent node of siblings 'node.indices'.
  #  mutation.probability.tables.list : The list of mutation probabilities for
  #                                     measured branch lengths.
  #  unknown.annot                    : The label of the unknown annotation.
  #                                     Default is 'unknown'.
  #
  # Returns: A list of CPTs, one for each phylogenetic node in argument
  # 'phylo.tree'.
  #   
  uniq.branch.lengths <- unique( phylo.tree$edge.length )
  annos.as.strs <- as.character( lapply( annotation.space,
    annotationToString ) )
  cpts <- conditionalProbabilityTables( uniq.branch.lengths,
    annotation.space, annos.as.strs,
    mutation.probability.tables.list, mutTblLengthColIndx=4
  )
  phylo.nodes <- c( get.root.node( phylo.tree ), phylo.tree$edge[ , 2 ] )
  unlist(
    lapply( phylo.nodes, function( phylo.node ) {
      bayesNode( phylo.tree, annotation.space, phylo.node, cpts,
        mutation.probability.tables.list, unknown.annot
      )
    }),
    recursive=F
  )
}

# Returns the labels of phylo.tree's tips, that have no annotation of specified
# annotation.type, or those that have non NA annotations, if negate=TRUE.
getTipsWithNaAnnotation <- function(phylo.tree, annotation.matrix, annotation.type='GO', negate=F) {
  bool.indxs <- is.na(annotation.matrix[annotation.type, ])
  surroundEachWithQuotes(
    colnames(annotation.matrix)[which( if(negate) ! bool.indxs else bool.indxs)]
  )
}

# Each entry of char.vector is surrounded by "\"".
surroundEachWithQuotes <- function(char.vector) {
  sapply(char.vector,
    function(x){
      paste("\"", x, "\"", sep='')
    },
    USE.NAMES=F)
}

annotationMatrixForBayesNetwork <- function( annotation.matrix,
  all.accessions=colnames( annotation.matrix ), annotation.type='GO',
  unknown.annot='unknown' ) {
  # Prepares an annotation.matrix to be used i.e. as diagnostic evidence inside
  # an Bayesian independence network.
  #
  # Args:
  #  annotation.matrix : A matrix of function annotations for some leaves in a
  #                      phylogenetic tree. Use for example all experimentally
  #                      verified molecular function annotations. The matrix's
  #                      columns should be the protein accessions and the row
  #                      the type of annotation, i.e. 'GO'.
  #  all.accessions    : The accessions of the query protein's homologs. If
  #                      more than in colnames( annotation.matrix ) those
  #                      without an experimentally verified function annotation
  #                      will be annotated with 'unknown.annot'.
  #  annotation.type   : The row of the 'annotation.matrix' to select.
  #  unknown.annot     : The unkown function annotation, set to any convenient
  #                      string.
  #
  # Returns: An annotation matrix, in which each protein accession is
  # surrounded by escaped quotes (surroundEachWithQuotes) and each cell is the
  # compound and sorted set of annotation terms (annotationToString). If
  # requested unkown annotations are added for protein accessions missing
  # annotations in argument 'annotation.matrix'.
  #   
  if ( ! is.null( annotation.matrix ) ) {
    am <- do.call( 'cbind', setNames(
      lapply( annotation.matrix[ annotation.type, ],
        function( a ) annotationToString( sort( a ) ) ),
      surroundEachWithQuotes( colnames( annotation.matrix ) )
      )
    )
    rownames( am ) <- annotation.type
    # Add 'unkown' annotations, if needed:
    accs.without.annos <- setdiff( all.accessions,
      colnames( annotation.matrix )
    )
    if ( length( accs.without.annos ) > 0 ) {
      am <- cbind( am,
        matrix( unknown.annot, ncol=length( accs.without.annos ), nrow=1,
          dimnames=list( annotation.type, surroundEachWithQuotes(
           accs.without.annos ) )
        )
      )
    }
    # return
    am
  }
}

diagnosticEvidence <- function( uniprot.accessions,
  valid.go.annos=names(GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE)
  ) {
  # Uses function 'retrieveExperimentallyVerifiedGOAnnotations' to download
  # experimentally verified Gene Ontology annotations for the argument
  # 'uniprot.accessions'. 
  #
  # Args:
  #  uniprot.accessions : A character vector of Uniprot accessions.
  #  valid.go.annos     : Those GO terms, that are valid in the model, others
  #                       will have to be excluded from the diagnostic
  #                       evidence, otherwise gRain will raise a 'Finding for
  #                       row 1 has probability 0 in then model.  Exiting...'
  #                       error message.
  #
  # Returns: A data.frame with column names being the argument Uniprot
  # accessions and a single row 'GO'. Each cell holds the list of annotated GO
  # terms. Accessions without experimentally verified GO annotations will be
  # ommited.  
  #   
  diag.evdnc <- as.data.frame(
    retrieveExperimentallyVerifiedGOAnnotations( uniprot.accessions )
  )
  colnames( diag.evdnc ) <- surroundEachWithQuotes( colnames( diag.evdnc ) )
  # Exclude annotations not present in our model:
  for ( i in 1:ncol(diag.evdnc) ) {
    diag.evdnc[[ 1, i ]] <- intersect( diag.evdnc[[ 1, i ]], valid.go.annos )
  }
  diag.evdnc
}

queryPhylBayesNetwork <- function( phylo.tree, response, annotation.matrix,
  annotation.space, diagnostic.evidence, annotation.type='GO',
  bayes.netw=grain( compileCPT( bayesNodes( phylo.tree,
        annotation.matrix, annotation.space, annotation.type ) )
  ), type='distribution' ) {
  # DEPRECATED!
  #
  # Creates and queries the Bayesian independence network. Argument
  # 'diagnostic.evidence' holds the predictor variables and their states those
  # tips of the phylogenetic tree, that are not referenced in the evidence will
  # be response variables. Meaning the probability distribution of their states
  # will be computed and returned.
  #
  # Args:
  #  phylo.tree          : The phylogenetic tree as returned by read.tree from
  #                        package 'ape'.
  #  annotation.matrix   : The matrix of tree tips annotations as returned by
  #                        either function
  #                        'retrieveExperimentallyVerifiedGOAnnotations' or
  #                        'retrieveAnnotationsBiomart', the former gives more
  #                        trustworthy results, the latter might yield more
  #                        sensitivity. For GO prediction we strongly suggest
  #                        to pass the annotation.matrix through function
  #                        'goTypeAnnotationMatrices' and invoke this function
  #                        for each type of GO term separatly. I.e. with
  #                        <your.go.type.anno.matrices>$molecular_function.
  #  annotation.space    : The set of unique compound annotations,
  #                        i.e. entry 'molecular_function' of the
  #                        list returned by function
  #                        'goAnnotationSpaceList'. Choose the same GO type as
  #                        in argument 'annotation.matrix', i.e.
  #                        'molecular_function'.
  #  diagnostic.evidence : The evidence matrix as returned by function
  #                        'retrieveExperimentallyVerifiedGOAnnotations'. If
  #                        you choose to get very trustworthy results (see
  #                        comment on argument 'annotation.matrix'), this
  #                        argument might be the identical object as
  #                        'annotation.matrix'.
  #  response            : The list of nodes in 'phylo.tree' to predict states
  #                        for. It is crucial to name the response nodes with
  #                        escaped quotes. So "MyProtein" should be
  #                        "\"MyProtein\"", use function surroundEachWithQuotes
  #                        to achieve this.
  #  annotation.type     : The type of annotation to process, default is 'GO'.
  #  bayes.netw          : The Bayesian network to generate and query.
  #                        Typically generated by the other arguments.
  #  type                : The type of prediction to return, either 'state' or
  #                        'distribution'. See ?predict.grain for more details.
  #
  # Returns: The result of calling either querygrain or predict.grain. A named
  # list with a data.frame for each node in argument 'response'. 
  #
  warning( "queryPhylBayesNetwork is DEPRECATED" )
  if ( ncol(diagnostic.evidence) == 0 && nrow(diagnostic.evidence) == 0 ) {
    # Just propagate the probability distributions from root to leaves, without
    # any diagnostic.evidence:
    querygrain( bayes.netw, nodes=response, result="data.frame" )
  } else {
    # Feed diagnostic.evidence into the network and then compute the
    # probability distributions based on both the prior and posterior
    # knowledge:
    predict(
      bayes.netw, response=response, newdata=diagnostic.evidence, type=type
    )
  }
}
