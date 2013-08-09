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

goTermForAccessionOrSynonym <- function( acc.or.synonym,
  con=connectToGeneOntology() ) {
  go.term <- goTermForAccession( acc.or.synonym )
  if ( nrow( go.term ) == 0 && ncol( go.term ) == 0 ) {
    go.term <- dbGetQuery( con, paste(
        "SElECT t.* FROM term_synonym s LEFT JOIN term t ON ",
        "s.term_id = t.id WHERE s.term_synonym = '",
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
      "SElECT t.*, g.relation_distance FROM term_synonym s LEFT JOIN term t ON ",
      "s.term_id = t.id LEFT JOIN graph_path g ON t.id = g.term2_id WHERE ",
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

goProfile <- function( accessions, go.level=3, con=connectToGeneOntology() ) {
  # Measures the annotation frequencies of parent Gene Ontology (GO) terms that
  # are parents of the argument 'accessions' terms and have argument 'go.level'
  # distance to the root of the GO directed acyclic graph (GO DAG). 
  #
  # Args:
  #  accessions : The annotated GO term accessions
  #  go.level   : The distance of parent GO terms to the GO DAG's root node
  #  con        : A valid and active database connection
  #
  # Returns: A data frame of parent GO terms with their annotation frequencies.
  #   
  go.terms <- goTermsForAccessionOrSynonymWithLevel( accessions, con=go.con )
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
  lapply( accessions, function( go.acc ) {
    go.id <- go.terms[ which( go.terms$acc == go.acc ), ]$id
    prnt.go.trms <- go.prnts[[ as.character( go.id ) ]]
    if ( nrow( prnt.go.trms ) > 0 ) {
      for ( i in 1:nrow( prnt.go.trms ) ) {
        prnt.go.trm <- prnt.go.trms[ i, , drop=FALSE ]
        go.profile <<- setFrequeny( go.profile, prnt.go.trm )
      }
    } else if (
      go.terms[ which( go.terms$id == go.id ), ]$relation_distance == go.level
    ) { 
      go.profile <<- setFrequeny( go.profile,
        go.terms[ which( go.terms$id == go.id ), ] )
    }
  })
  go.profile
}

spawnedGoTerms <- function( go.term.id, con=connectToGeneOntology() ) {
  # Queries the Gene Ontology ( GO ) database to retrieve all GO terms placed
  # in the descending sub-graph spawned by GO term with ID 'go.term.id'.
  #
  # Args:
  #  go.term.id : The parent GO term's database ID 
  #  con        : A valid and alive database connection an instance of the Gene
  #               Ontology
  #
  # Returns: A table with all GO terms being direct and indirect children of
  # 'go.term.id'
  #   
  dbGetQuery( con,
    paste( "SELECT t.* FROM graph_path p ",
      "LEFT JOIN term t ON t.id = p.term2_id ",
      "where p.relationship_type_id = 1 ",
      "AND p.term1_id = ", go.term.id,
      " AND t.id != ", go.term.id,
      " ORDER BY p.distance",
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

goTermsForProteinAccessionAndEvidenceCodes <- function( prot.acc,
  ec=EVIDENCE.CODES, con=connectToGeneOntology() ) {
  dbGetQuery( con,
    paste( "select t.*, ev.code from term t left join association a on a.term_id = t.id ",
           "left join gene_product gp on gp.id = a.gene_product_id ",
           "left join dbxref d on d.id = gp.dbxref_id ",
           "left join evidence ev on ev.association_id = a.id ",
           "where d.xref_key ='", prot.acc, "' ",
           "and ev.code in ('",
           paste( ec, collapse="','" ),
           "')", sep=""
    )
  )
}
