join.funk <- function( ... ) { paste( ..., sep=", " ) }

connectToGeneOntology <- function( driver=MySQL(), user="go_select",
  password="amigo", dbname="go_latest", host="mysql.ebi.ac.uk", port=4085 ) {
  # Connects to the latest gene ontology database using DBI and MySQL.
  #
  # Args:
  #  driver   : The database driver to use, default MySQL
  #  user     : The database user, default go_select
  #  password : The database password, default amigo
  #  dbname   : The database name to use, default go_latest
  #  host     : The database host to connect to, default mysql.ebi.ac.uk
  #  port     : The port on the host to connect to, default 4085
  #
  # Returns: A database connection as created by DBI using the arguments.
  #   
  dbConnect( driver, user=user, password=password, dbname=dbname,
  host=host, port=port )
}

goTermForAccession <- function( accession, con=connectToGeneOntology() ) {
  dbGetQuery( con,
    paste( "SELECT * FROM term WHERE acc = '",
      accession, "'", sep="" )
  )
}

goTermsForAccession <- function( accessions, con=connectToGeneOntology() ) {
  dbGetQuery( con,
    paste( "SELECT * FROM term WHERE acc in (",
      do.call( 'join.funk', as.list( paste( "'", accessions, "'", sep="" ) ) ),
      ")", sep=""
    )
  )
}

goTermsForAccessionWithDefinition <- function( accessions,
  con=connectToGeneOntology() ) {
  # Selects all Gene Ontology (GO) terms whose accession equals argument
  # 'accessions'. Also looks up the terms definitions, longer texts explaining
  # the GO terms in more detail.
  #
  # Args:
  #  accessions : GO term accessions as to be found in table 'term.acc'
  #  con        : A valid and active database connection to an instance of the
  #               GO, default is obtained from calling function
  #               connectToGeneOntology(…)
  #
  # Returns: A data frame with the query results.
  #   
  dbGetQuery( con,
    paste( "SELECT t.*, d.term_definition FROM term t ",
      "LEFT JOIN term_definition d ON t.id = d.term_id WHERE t.acc in (",
      paste( paste( '"', accessions, '"', sep='' ), collapse=',' ),
      ")", sep=""
    )
  )
}

goTermsForAccessionWithDefinitionAndLevel <- function( accessions,
  con=connectToGeneOntology() ) {
  # Selects all Gene Ontology (GO) terms whose accession equals argument
  # 'accessions'. Also looks up the terms distance to the GO directed acyclic
  # graph's (DAG) root, as well as the terms' definitions, longer texts
  # explaining the GO terms in more detail.
  #
  # Args:
  #  accessions : GO term accessions as to be found in table 'term.acc'
  #  con        : A valid and active database connection to an instance of the
  #               GO, default is obtained from calling function
  #               connectToGeneOntology(…)
  #
  # Returns: A data frame with the query results.
  #   
  go.accs <- paste( paste( '"', accessions, '"', sep='' ), collapse=',' )
  dbGetQuery( con,
    paste( "SELECT t.*, d.term_definition, p.relation_distance, s.term_synonym ",
      "FROM term t LEFT JOIN graph_path p ON t.id = p.term2_id AND ",
      "p.term1_id = (SELECT r.id FROM term r WHERE r.is_root = 1) ",
      "LEFT JOIN term_definition d ON t.id = d.term_id ",
      "LEFT JOIN term_synonym s ON t.id = s.term_id ",
      "WHERE t.acc in (", go.accs, ") OR ",
      "s.term_synonym in (", go.accs, ") ",
      "GROUP BY t.id ORDER BY p.relation_distance", sep=""
    )
  )
}

goTermForAccessionOrSynonym <- function( acc.or.synonym,
  con=connectToGeneOntology() ) {
  go.term <- goTermForAccession( acc.or.synonym )
  if ( nrow( go.term ) == 0 && ncol( go.term ) == 0 ) {
    go.term <- dbGetQuery( con, paste(
        "SElECT t.* FROM term t LEFT JOIN term_synonym s ON ",
        "t.id = s.term_id WHERE s.term_synonym = '",
        acc.or.synonym,
        "'", sep=""
      )
    )
  }
  go.term
}

goTermsForAccessionOrSynonymWithLevel <- function( accessions,
  con=connectToGeneOntology() ) {
  # Looks up Gene Ontology (GO) terms in the term and term_synonym table, thus
  # finding also outdated GO terms.
  #
  # Args:
  #  accessions : A character vector of GO term accessions
  #  con        : A valid and active database connection
  #
  # Returns: A data frame
  #   
  accs <- unique( accessions )
  dbGetQuery( con, paste(
      "SElECT t.*, g.relation_distance FROM term t LEFT JOIN term_synonym s ON ",
      "t.id = s.term_id LEFT JOIN graph_path g ON t.id = g.term2_id WHERE ",
      "g.term1_id = ( SELECT r.id FROM term r WHERE r.is_root = 1 ) AND ",
      "t.acc in (",
      paste( paste( "'", accs, "'", sep='' ), collapse=',' ),
      ") OR ",
      "s.term_synonym in (",
      paste( paste( "'", accs, "'", sep='' ), collapse=',' ),
      ") GROUP BY t.id", sep=""
    )
  )
}

parentGoTerms <- function( go.term.id, con=connectToGeneOntology() ) {
  dbGetQuery( con,
    paste( "SELECT t.* FROM graph_path res ",
      "LEFT JOIN term t ON t.id = res.term1_id ",
      "WHERE res.relationship_type_id = 1 ",
      "AND res.term1_id != (SELECT r.id FROM term r WHERE r.is_root = 1) ",
      "AND res.term2_id = ",
      go.term.id,
      " AND t.id != ",
      go.term.id,
      " ORDER BY res.distance DESC", 
      sep=""
    )
  )
}

ancestralGoTerms <- function( go.term.id, con=connectToGeneOntology() ) {
  # Selects all Gene Ontology (GO) terms that have any ancestral
  # 'gene_ontology' relationship to the argument descendant term 'go.term.id'.
  # All available 'gene_ontology' relationships are selected from the term
  # table itself.
  #
  # Args:
  #  go.term.id : The descendant GO term to find ancestral GO terms for with
  #               any 'gene_ontology' type relationship. 
  #  con        : A valid and active database connection to an instance of the
  #               Gene Ontology mysql database.
  #
  # Returns: A data frame, excluding both the 'all' GO DAG root node as well as
  # the self match 'go.term.id'.
  #   
  dbGetQuery( con,
    paste( "SELECT t.*, r.acc AS relationship FROM term t ",
      "LEFT JOIN graph_path p ",
      "ON t.id = p.term1_id LEFT JOIN term r ON p.relationship_type_id = r.id ",
      "WHERE p.term2_id = ", go.term.id,
      " AND NOT t.is_root AND t.id != ", go.term.id, " GROUP BY t.id",
      sep="" )
  )
}

parentGoTermsForAccession <- function( go.term.accs, include.selves=FALSE,
  relationship.type.id=1, con=connectToGeneOntology() ) {
  # Finds the Gene Ontology (GO) terms that are parent to the argument
  # 'go.term.accs' and are NOT obsolete.
  #
  # Args:
  #  go.term.accs         : The GO term accessions of the GO terms to find the
  #                         parents for.
  #  include.selves       : If set to TRUE the result includes the GO terms
  #                         given by argument 'go.term.accs'. The default is
  #                         FALSE.
  #  relationship.type.id : The type of parent child relationship to restrict
  #                         the query to. The default '1' refers to 'is_a'. To
  #                         select all available relationships provide this
  #                         argument as NULL.
  #  con                  : A valid and active database connection to an
  #                         instance of the Gene Ontology relation database.
  #
  # Returns: A data frame.
  #
  gta <- paste( paste( "'", go.term.accs, "'", sep='' ), collapse=',' )
  incl.selves.sql <- if ( include.selves ) {
    ''
  } else {
    paste( "AND t.acc not in (", gta, ") ", sep='' )
  }
  rel.sql <- if ( ! is.null( relationship.type.id ) ) {
    paste( 'res.relationship_type_id = ', relationship.type.id,
      ' AND ', sep='' )
  } else {
    ''
  }
  sql <- paste( "SELECT t.*, to_root.relation_distance, child.acc as child_acc ",
      "FROM graph_path res LEFT JOIN term t ON t.id = res.term1_id ",
      "LEFT JOIN graph_path to_root ON t.id = to_root.term2_id ",
      "LEFT JOIN term child ON child.id = res.term2_id ",
      "WHERE ", rel.sql,
      "res.term1_id != (SELECT r.id FROM term r WHERE r.is_root = 1) ",
      "AND child.acc in (", gta, ") ",
      incl.selves.sql,
      "AND to_root.term1_id = (SELECT r.id FROM term r WHERE r.is_root = 1) ",
      "AND t.is_obsolete = 0 ",
      "GROUP BY t.id ORDER BY to_root.relation_distance ASC",
      sep=""
    )
  dbGetQuery( con, sql )
}

ancestorGoTermsForAccession <- function( go.term.accs,
  con=connectToGeneOntology() ) {
  # Selects those Gene Ontology (GO) terms that are ancestors ( parents ) of
  # the GO terms whose accessions are given in argument 'go.term.accs'. In
  # contrast to the other 'parent' functions, this function includes the "all"
  # root node and does not return the distance to the root.
  #
  # Args:
  #
  #  go.term.accs : The GO term accessions of the GO terms to find the parents
  #                 for.
  #  con          : A valid and active database connection to an instance of
  #                 the Gene Ontology relation database.
  #
  # Returns: A data frame with results from the GO term table.
  #   
  gta <- paste( paste( "'", go.term.accs, "'", sep='' ), collapse=',' )
  sql <- paste( "SELECT t.* FROM term_ancestor a ",
    "LEFT JOIN term t ON a.ancestor_id = t.id WHERE ",
    "a.acc IN (", gta, ") AND a.distance > 0 GROUP BY t.id" )
  dbGetQuery( con, sql )
}

parentGoTermsOfLevel <- function( go.term.id, go.level=3,
  con=connectToGeneOntology() ) {
  # Finds the Gene Ontology (GO) terms that are parent to the argument
  # 'go.term.id' and have argument 'go.level' distance to the root of the GO
  # directed acyclic graph. 
  #
  # Args:
  #  go.term.id : The database identifier of the GO term to find the parent
  #               for.
  #  go.level   : The distance of the looked for parent term to the GO DAG's
  #               root node.
  #  con        : A valid and active database connection to an instance of the
  #               Gene Ontology relation database.
  #
  # Returns: A data frame.
  #
  dbGetQuery( con,
    paste( "SELECT t.*, to_root.relation_distance FROM graph_path res ",
      "LEFT JOIN term t ON t.id = res.term1_id ",
      "LEFT JOIN graph_path to_root ON t.id = to_root.term2_id ",
      "WHERE res.relationship_type_id = 1 ",
      "AND res.term1_id != (SELECT r.id FROM term r WHERE r.is_root = 1) ",
      "AND res.term2_id = ",
      go.term.id, " ",
      "AND t.id != ",
      go.term.id, " ",
      "AND to_root.term1_id = (SELECT r.id FROM term r WHERE r.is_root = 1) ",
      "AND to_root.relation_distance = ", go.level, " ",
      "ORDER BY to_root.relation_distance ASC",
      sep=""
    )
  )
}


parentGoTermsOfLevelForAccessions <- function( go.term.accs, go.level=3,
  con=connectToGeneOntology() ) {
  # Finds the Gene Ontology (GO) terms that are parent to the argument
  # 'go.term.accs' and have argument 'go.level' distance to the root of the GO
  # directed acyclic graph. 
  #
  # Args:
  #  go.term.accs : The GO term accessions of the GO terms to find the parents
  #                 for.
  #  go.level     : The distance of the looked for parent term to the GO DAG's
  #                 root node.
  #  con          : A valid and active database connection to an instance of
  #                 the Gene Ontology relation database.
  #
  # Returns: A data frame.
  #
  gta <- paste( paste( "'", go.term.accs, "'", sep='' ), collapse=',' )
  sql <- paste( "SELECT t.*, to_root.relation_distance FROM graph_path res ",
      "LEFT JOIN term t ON t.id = res.term1_id ",
      "LEFT JOIN graph_path to_root ON t.id = to_root.term2_id ",
      "LEFT JOIN term child ON child.id = res.term2_id ",
      "WHERE res.relationship_type_id = 1 ",
      "AND res.term1_id != (SELECT r.id FROM term r WHERE r.is_root = 1) ",
      "AND child.acc in (", gta, ") ",
      "AND t.acc not in (", gta, ") ",
      "AND to_root.term1_id = (SELECT r.id FROM term r WHERE r.is_root = 1) ",
      "AND to_root.relation_distance = ", go.level, " ",
      "GROUP BY t.id ORDER BY to_root.relation_distance ASC",
      sep=""
    )
  dbGetQuery( con, sql )
}

parentChildRelation <- function( go.1, go.2, con=connectToGeneOntology() ) {
  # Returns TRUE if either of the argument Gene Ontology (GO) term accessions
  # is an ancestor of the other one or if the accessions are identical; FALSE
  # is returned otherwise.
  #
  # Args:
  #  go.1 : The first GO term accesion
  #  go.2 : The second GO term accession
  #  con  : A valid and active database connection to an instance of
  #         the Gene Ontology relation database.
  #
  # Returns: A boolean value explained above.
  #   
  go.1.prnts <- ancestorGoTermsForAccession( go.1, con=con )
  go.2.prnts <- ancestorGoTermsForAccession( go.2, con=con )
  (
    is.element( go.1, go.2.prnts$acc ) || is.element( go.2, go.1.prnts$acc ) ||
    identical( go.1, go.2 )
  )
}


goProfile <- function( accessions, go.level=3, con=connectToGeneOntology(),
  close.db.con=TRUE ) {
  # Measures the annotation frequencies of parent Gene Ontology (GO) terms that
  # are parents of the argument 'accessions' terms and have argument 'go.level'
  # distance to the root of the GO directed acyclic graph (GO DAG). 
  #
  # Args:
  #  accessions   : The annotated GO term accessions
  #  go.level     : The distance of parent GO terms to the GO DAG's root node
  #  con          : A valid and active database connection
  #  close.db.con : If TRUE the data connection 'con' is automatically closed
  #                 at the of this function's execution.
  #
  # Returns: A data frame of parent GO terms with their annotation frequencies.
  #   
  go.terms <- goTermsForAccessionOrSynonymWithLevel( accessions, con=con )
  go.prnts <- setNames(
    lapply( as.integer( go.terms$id ), function( go.id ) {
      parentGoTermsOfLevel( go.id, go.level=go.level, con=con )
    }), 
    as.character( go.terms$id )
  )
  setFrequeny <- function( go.profile, prnt.go.trm ) {
    if (
        is.null( go.profile ) ||
        length( intersect( go.profile$id, prnt.go.trm$id ) ) == 0
    ) {
      go.profile <- rbind( go.profile,
        cbind( prnt.go.trm, list( 'frequency'=1 ) )
      )
    } else {
      go.profile[ which( go.profile$id == prnt.go.trm$id ), ]$frequency <-
        go.profile[ which( go.profile$id == prnt.go.trm$id ), ]$frequency + 1
    }
    go.profile
  }
  go.profile <- NULL
  lapply( as.character( go.terms$acc ), function( go.acc ) {
    go.id <- go.terms[ which( go.terms$acc == go.acc ), ]$id
    prnt.go.trms <- go.prnts[[ as.character( go.id ) ]]
    if ( nrow( prnt.go.trms ) > 0 ) {
      for ( i in 1:nrow( prnt.go.trms ) ) {
        prnt.go.trm <- prnt.go.trms[ i, , drop=FALSE ]
        go.profile <<- setFrequeny( go.profile, prnt.go.trm )
      }
    } else if (
      go.terms[ which( go.terms$id == go.id ), ]$relation_distance ==
        go.level
    ) { 
      go.profile <<- setFrequeny( go.profile,
        go.terms[ which( go.terms$id == go.id ), ] )
    }
  })
  if ( close.db.con ) {
    dbDisconnect( con )
  }
  go.profile
}

spawnedGoTerms <- function( go.term.id, relationship.type.id=1,
  include.selve=TRUE, con=connectToGeneOntology() ) {
  # Queries the Gene Ontology ( GO ) database to retrieve all GO terms placed
  # in the descending sub-graph spawned by GO term with ID 'go.term.id'.
  #
  # Args:
  #  go.term.id           : The parent GO term's database ID
  #  relationship.type.id : The type of parent child relationship to restrict
  #                         the query to. The default '1' refers to 'is_a'. To
  #                         select all available relationships provide this
  #                         argument as NULL.
  #  include.selve        : If TRUE, the default, the term specified by
  #                         argument 'go.term.id' will be included in the
  #                         result.
  #  con                  : A valid and alive database connection an instance
  #                         of the Gene Ontology
  #
  # Returns: A table with all GO terms being direct and indirect children of
  # 'go.term.id'
  #   
  rel.sql <- if ( ! is.null( relationship.type.id ) ) {
    paste( 'p.relationship_type_id = ', relationship.type.id,
      ' AND ', sep='' )
  } else {
    ''
  }
  incl.slv.sql <- if ( include.selve ) {
    ''
  } else {
    paste( " t.id !=", go.term.id, 'AND ' )
  }
  dbGetQuery( con,
    paste( "SELECT t.* FROM graph_path p ",
      "LEFT JOIN term t ON t.id = p.term2_id WHERE ",
      rel.sql, incl.slv.sql, 
      "p.term1_id = ", go.term.id,
      " GROUP BY t.id ORDER BY p.distance",
      sep=''
    )
  )
}

goTermsOfLevelAndType <- function( level, term.type,
  con=connectToGeneOntology() ) {
  dbGetQuery( con,
    paste(
      "SELECT t.* FROM graph_path g ",
      "LEFT JOIN term t ON g.term2_id = t.id ",
      "WHERE t.term_type = '", term.type ,"' AND ",
      "g.relation_distance = ", level ," AND ",
      "g.relationship_type_id = 1 AND ",
      "g.term1_id = ( SELECT r.id FROM term r WHERE r.is_root = 1 )",
      sep=""
    )
  )
}

goTermsForAccessionWithLevel <- function( accessions,
  con=connectToGeneOntology() ) {
  # For argument Gene Ontology (GO) term accession all data in the SQL term
  # table and additionaly their respective distance to the GO directed acyclic
  # graph's root node is looked up.
  #
  # Args:
  #  accessions         : The GO terms' accessions ( column 'acc' in table
  #                       'term' )
  #  con                : A valid and active database connection
  #
  # Returns: A data frame.
  #   
  dbGetQuery( con, paste(
      "SELECT t.*, g.relation_distance FROM term t LEFT JOIN graph_path g ON ",
      "t.id = g.term2_id WHERE g.term1_id = ( SELECT r.id FROM term r WHERE r.is_root = 1 ) ",
      "AND t.acc in (",
      do.call( 'join.funk',
        as.list( paste( "'", unique( accessions ), "'", sep="" ) ) ),
      ") GROUP BY t.id"
    ) 
  )
}

goTermsForProteinAccession <- function( prot.acc, con=connectToGeneOntology() ) {
  dbGetQuery( con,
    paste( "select t.* from term t left join association a on a.term_id = t.id ",
           "left join gene_product gp on gp.id = a.gene_product_id ",
           "left join dbxref d on d.id = gp.dbxref_id ",
           "where d.xref_key ='", prot.acc, "'", sep=""
    )
  )
}

goTermsForProteinAccessionsAndEvidenceCodes <- function( prot.accs,
  ec=EVIDENCE.CODES, con=connectToGeneOntology() ) {
  # Queries the Gene Ontology (GO) database accessible via DB connection 'con'
  # for all GO term annotations available for protein accessions 'prot.accs'.
  # If argument evidence codes 'ec' are provided only those annotations
  # matching any of the provided evidence codes will be returned.
  #
  # Args:
  #  prot.accs : The protein accessions to find GO term annotations for.
  #  ec        : A character vector of evidence codes, default is the constant
  #              'EVIDENCE.CODES'. Set this argument to NULL or 'ALL', if no
  #              restrictions should be made.
  #  con       : A valid and active MySQL database connection to an instance of
  #              the Gene Ontology database.
  #
  # Returns: A data frame with the found GO term annotations. Columns are
  # 'xref_key', 'term_type', 'acc', 'is_obsolete', 'is_root', 'is_relation',
  # 'code'.
  #   
  ec.sql <- if ( is.null( ec ) || is.na( ec ) || ec == 'ALL' || ec == 'all' ) {
    ''
  } else {
    paste( " AND ev.code in ('", paste( ec, collapse="','" ), "')", sep="" )
  }
  dbGetQuery( con,
    paste( "select d.xref_key, t.*, ev.code from term t ",
           "left join association a on a.term_id = t.id ",
           "left join gene_product gp on gp.id = a.gene_product_id ",
           "left join dbxref d on d.id = gp.dbxref_id ",
           "left join evidence ev on ev.association_id = a.id ",
           "where d.xref_key in ('",
           paste( prot.accs, collapse="','" ), "')", ec.sql , sep=""
    )
  )
}

extendGOAnnosWithParents <- function( go.anno.df, con=connectToGeneOntology(),
  close.db.con=TRUE ) {
  # For every GO term a protein is annotated with it will be annotated with all
  # ancestors of this GO term no matter what kind of relationship the ancestor
  # has to its descendant. The evidence code will be "reused". If requested,
  # the GO term type ('biological_process', 'cellular_component',
  # 'molecular_function') is also added to the annotation matrix.
  #
  # Args:
  #  go.anno.df    : A data frame of GO term annotations for proteins.  Columns
  #                  are expected to be 1. GO term accession, 2.  Evidence
  #                  Code, and 3. Protein accession. It is required to be
  #                  cleaned up by
  #                  uniqueGOAnnotationsWithMostSignificantEvidenceCodes(…).
  #  con           : A valid and active MySQL connection to an instance of the
  #                  GO database - default is connectToGeneOntology(…)
  #  close.db.con  : If set to TRUE, the database connection 'con' will be
  #                  closed automatically before this function returns
  #
  # Returns: A data frame, the extension of the argument 'go.anno.df', if
  # requested with an additional column holding the GO term types.
  #   
  unq.gos <- unique( go.anno.df[ , 1 ] )
  go.prnts <- parentGoTermsForAccession( unq.gos, include.selves=TRUE,
    relationship.type.id=NULL, con=con )
  if ( close.db.con ) dbDisconnect( con )

  d.f <- as.data.frame( matrix( vector(), nrow=0, ncol=4,
    dimnames=list( c(), c('acc', 'ec', 'prot.acc', 'term_type') ) ),
    stringsAsFactors=FALSE
  )
  for ( i in 1:nrow(go.prnts) ) {
    # For each prot select a df with cols
    # 'acc', 'ec', 'prot.acc', 'term_type'.
    # evidence.code should be inherited by original annotation!
    go.row <- go.prnts[ i, ]
    prot.annos <- go.anno.df[ which( go.anno.df[ , 1 ] ==
      go.row[[1,'child_acc']] ), ]
    d.f <- rbind( d.f, data.frame( list(
      acc = rep( go.row[[1, "acc"]], nrow(prot.annos) ), 
      ec = prot.annos[, 2],
      prot.acc = prot.annos[, 3],
      term_type = rep( go.row[[1, "term_type"]], nrow(prot.annos) )
    ), stringsAsFactors = FALSE ) ) 
  }
  d.f
}

selectMostSignificantEvidenceCode <- function( evidence.codes,
  experiment=c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP"),
  experiment.general="EXP",
  author=c( "TAS", "NAS"), curator=c( "IC", "ND"), automatic= "IEA" ) {
  # Identifies the most trustworthy evidence.code in a set of many.
  #
  # Args:
  #  evidence.codes     : The evidence.codes assigned to a single protein GO
  #                       term pair
  #  experiment         : The evidence.codes categorized as proved by
  #                       laboratory experiment
  #  experiment.general : The evidence.code to return, if more than one
  #                       experimental evidence codes are in argument
  #                       'evidence.codes'
  #  author             : The evidence.codes indicating a publication or
  #                       scientist statement yielding the GO term
  #                       annotation(s)
  #  curator            : The evidence.codes indicating a curator decision
  #                       yielding the GO term annotation(s)
  #  automatic          : The default evidence.code to return if none of the
  #                       above categories are matched
  #
  # Returns: The single evidence.code interpreted most significant and
  # trustworthy.
  #   
  if ( length( evidence.codes ) == 1 ) {
    evidence.codes
  } else if ( any( experiment %in% evidence.codes ) ) {
    experiment.general
  } else if ( any( author %in% evidence.codes ) ) {
    intersect( author, evidence.codes )[[1]]
  } else if ( any( curator %in% evidence.codes ) ) {
    intersect( curator, evidence.codes )[[1]]
  } else {
    automatic[[1]] 
  }
}
