library(gRain)

# Returns unique node, that has no edge leading to it. Thus, it must be a node
# in the first column whereas mustn't appear in the second column of the
# phylTree's edge matrix.
get.root.node <- function(phylTree) {
  as.character(
    unique(phylTree$edge[!(phylTree$edge[, 1] %in% phylTree$edge[, 2]), ][, 1])[[1]]
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
  res <- try(phylTree$tip.label[[node.index]], silent=T)
  if(identical("try-error", class(res))) as.character(node.index) else res
}

edge.to.formula <- function(phyloTree, edge.index) {
  prnt <- get.node.label(phyloTree,
    phyloTree$edge[[edge.index, 1]])
  child <- get.node.label(phyloTree,
    phyloTree$edge[[edge.index, 2]])
  eval(bquote(~ .(child) | .(prnt) ))
}

findMatchingColumn <- function( mutation.probability.table, value, column.index ) {
  # Looks up the _last_ row of argument 'mutation.probability.table' who's
  # entry in column 'column.index' is smaller or equal to argument 'value'. If
  # none can be found returns the last row of mutation.probability.table. 
  #
  # Args:
  #  mutation.probability.table : The matrix holding mutation probabilities as
  #                               returned i.e. by function
  #                               'mutationProbabilityDistribution'.
  #  value                      : The current distance value, i.e. the branch
  #                               length to find the correct mutation
  #                               probability for.
  #  column.index               : The index of the column holding the distance
  #                               measure to compare argument value with.
  #
  # Returns: The argument mutation.probability.table's matching row as matrix.
  #   
  hits <- mutation.probability.table[
    which( mutation.probability.table[ , column.index ] >= value ), , drop=F
  ]
  if ( nrow( hits ) > 0 ) {
    hits[ 1, , drop=F ]
  } else {
    mutation.probability.table[ nrow(mutation.probability.table), , drop=F ]
  }
}

goTypeAnnotationMatrices <- function( annotation.matrix,
  valid.go.terms=names( GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE ),
  go.con=connectToGeneOntology(), exclude.empty.cols=T ) {
  # For each type of Gene Ontology terms, 'molecular_function',
  # 'biological_process', and 'cellular_component', a separate annotation
  # matrix is generated. These matrices hold the set of GO term annotations for
  # each protein as in the column names of argument 'annotation.matrix'. If a
  # non NULL and non empty valid.go.terms list is provided only GO term
  # annotations appearing in this list will appear.
  #
  # Args:
  #  annotation.matrix  : A funtion annotation matrix as returned by
  #                       'retrieveAnnotationsBiomart'.
  #  valid.go.terms     : If non NULL and non empty only terms appearing in
  #                       this list will be assigned.
  #  go.con             : database connection to the gene-ontology-MySQL-db
  #  exclude.empty.cols : If set to FALSE proteins with no annotations in the
  #                       current GO term type will have character(0) entries.
  #                       If this switch is set to TRUE, such columns will
  #                       simply be excluded from the resulting matrices. 
  #
  # Returns: A named list with three annotation matrices 'biological_process',
  # 'cellular_component', and 'molecular_function'. Each annotation matrix has
  # a single row and the same columns as argument 'annotation.matrix'. Each
  # cell contains the filtered GO annotations of the corresponding GO term
  # type.
  #   
  all.distinct.annos <- uniq.annotations( annotation.matrix, 'GO', exclude.NAs=T )
  annos <- goTermsForAccessionWithLevel( all.distinct.annos )
  go.types <- c( 'biological_process', 'cellular_component', 'molecular_function' ) 
  setNames(
    lapply( go.types, function( go.tp ) {
      gos.of.type <- annos[ which( annos$term_type == go.tp ), 'acc' ]
      if ( ! is.null( valid.go.terms ) && length( valid.go.terms ) > 0 )
        gos.of.type <- intersect( gos.of.type, valid.go.terms )
      do.call( 'cbind',
        lapply( colnames( annotation.matrix ), function( prot ) {
          prot.mtrx <- matrix( list(), ncol=1, nrow=1, dimnames=list( 'GO', prot ) )
          prot.mtrx[[ 'GO', prot ]] <- sort( intersect( annotation.matrix[[ 'GO', prot ]], gos.of.type ) )
          if ( ! is.null(prot.mtrx[[ 'GO', prot ]]) &&
              length(prot.mtrx[[ 'GO', prot ]]) > 0 ||
              ! exclude.empty.cols ) 
            prot.mtrx
          else
            NULL
        })
      )
    }),
    go.types
  )
}

mutationProbability <- function( annotation, branch.length,
  annots.mut.prob.table.list=GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE,
  distance.column.index, select.funk=max ) {
  # Returns the maximum mutation probability found for any annotation in the
  # argument list 'annotation'.
  #
  # Args:
  #  annotation                 : A character vector of function annotations,
  #                               i.e. c( 'GO_A', 'GO_B', 'GO_C' )
  #  branch.length              : The length of the current phylogenetic
  #                               branch.
  #  annots.mut.prob.table.list : The list to lookup each single annotations's
  #                               mutation probability table.
  #  distance.column.index      : The index of the column in which above
  #                               probability tables hold the maximum distance
  #                               measures mapped to their corresponding
  #                               mutation probabilities.
  #  select.funk                : Set to min, mean or max to define which
  #                               mutation probability should be returned.
  #                               Default is max.
  #
  # Returns: The numeric maximum of all matching mutation probabilities.
  #   
  select.funk(
    as.numeric(
      lapply( annotation, function( singl.anno ) {
        findMatchingColumn( annots.mut.prob.table.list[[ singl.anno ]],
          branch.length, distance.column.index )[[ 1, 1 ]]
      })
    )
  )
}

conditionalProbsTbl <- function( edge.length, annos,
  annots.mut.prob.table.list, mut.tbl.length.col.indx=5, p.mut.col.indx=1,
  unknown.annot='unknown' ) {
  # Looks up the mutation probability tables for each annotation and constructs
  # the transition probabilities based on the current branch's length.
  # Introduces a new annotation 'UNKNOWN' which can mutate to every other
  # annotation in argument 'uniq.annotations' and does have a zero probability
  # of self retainment.
  #
  # Args:
  #  edge.length : sequence distance as the current branch's length.
  #  annos : The unique anntotations as they appear in the current phylogenetic
  #          tree. Should NOT contain NULL or NA values. 
  #  annots.mut.prob.table.list : The list of empirical mutation probabilities
  #                               for each annotation.
  #  mut.tbl.length.col.indx : The index of the mutation probability table's
  #                            column in which to find the corresponding
  #                            sequence distance.
  #  p.mut.col.indx : The mutation probability column's index.
  #  unknown.annot : The UNKNOWN annotation, set to any convenient name.
  #
  # Returns: 
  #   A length(annos)*length(annos) probability matrix with
  #   the m[i, j] := P(j | i), i = child function, j = parent function
  #   
  do.call( 'cbind', 
    lapply( annos, function( anno ) {
      a <- annotationToString( anno )
      rownms <- lapply( annos, annotationToString )
      if ( identical( a, unknown.annot ) ) {
        trans.probs <- matrix( 1 / ( length( annos ) - 1 ),
          ncol=1, nrow=length( annos ),
          dimnames=list( rownms, a )
        )
        # unknown retainment probability:
        trans.probs[[ a, a ]] <- 0
        trans.probs
      } else {
        p.mut <- mutationProbability( anno, edge.length,
          annots.mut.prob.table.list, mut.tbl.length.col.indx
        )
        trans.probs <- matrix( ( p.mut / ( length( annos ) - 1 ) ),
          ncol=1, nrow=length( annos ),
          dimnames=list( rownms, a )
        )
        # annotation retainment probability:
        trans.probs[[ a, a ]] <- 1 - p.mut
        trans.probs
      }
    })
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
  #                          conditionalProbsTbl.
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

annotationSpace <- function( annotation.matrix, annotation.type='GO' ) {
  # Generates the set of unique anntations found in the columns of argument
  # 'annotation.matrix' in row 'annotation.type'.
  #
  # IMPORTANT NOTE: This function expects each proteins annotations vector to
  # be alphabetically sorted! - See function goTypeAnnotationMatrices for
  # further details.
  #
  # Args:
  #  annotation.matrix : The matrix of protein annotations as returned by
  #                      retrieveAnnotationsBiomart.
  #  annotation.type   : The type of annotations to select, default 'GO'.
  #
  # Returns: A character vector of unique annotations.
  #   
  unique( annotation.matrix[ annotation.type, ] )
}

bayesNode <- function( phylo.tree, annotation.space, node.index,
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
    conditionalProbsTbl( phylo.tree$edge.length[[ edge.ind ]],
      annotation.space, mutation.probability.tables.list,
      unknown.annot=unknown.annot )
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
  unknown.annot='unknown', lapply.funk=mclapply ) {
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
  #  lapply.funk                      : Set to 'mclapply' if parallel
  #                                     computation of CPTs is wanted, set to
  #                                     lapply, if serial computation is
  #                                     preferred. 'mclapply' is *strongly*
  #                                     recommended when used on trees with more
  #                                     than 250 nodes. Default is 'mclapply'.
  #
  # Returns: A list of CPTs, one for each phylogenetic node in argument
  # 'phylo.tree'.
  #   
  phylo.nodes <- c( get.root.node( phylo.tree ), phylo.tree$edge[ , 2 ] )
  unlist(
    lapply.funk( phylo.nodes, function( phylo.node ) {
      bayesNode( phylo.tree, annotation.space, phylo.node,
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
