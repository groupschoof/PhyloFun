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
  # Finds the parent GO term for argument term 'accession', where the parent
  # has distance 'go.level' to the GO directed acyclic graph (GO-DAG) root.
  #
  # Args:
  #  go.term.id : The database identifier of the argument GO term  
  #  go.level  : The distance of the requested parent GO term to the GO-DAG's
  #              root node
  #  con       : A valid and active database connection to the Gene Ontology
  #              realtional database.
  #
  # Returns: A data.frame holding the requested row, if found. Column names are
  # those of the returned SQL set. 
  #   
  dbGetQuery( con,
    paste( "SELECT t.*, res.relation_distance FROM graph_path res ",
      "LEFT JOIN term t ON t.id = res.term1_id ",
      "WHERE res.relationship_type_id = 1 ",
      "AND res.relation_distance = ",
      go.level,
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

goTermsForAccessionWithLevel <- function( accessions, con=connectToGeneOntology() ) {
  dbGetQuery( con, paste(
      "SELECT t.*, g.relation_distance FROM term t LEFT JOIN graph_path g ON ",
      "t.id = g.term2_id WHERE g.term1_id = ( SELECT r.id FROM term r WHERE r.is_root = 1 ) ",
      "AND t.acc in (",
      do.call( 'join.funk', as.list( paste( "'", accessions, "'", sep="" ) ) ),
      ") GROUP BY t.id"
    ) 
  )
}

isConnectionAlive <- function( go.con ) {
  if ( class( try( dbGetInfo( go.con ), silent=T ) ) == 'try-error' )
    FALSE
  else
    TRUE
}

reConnectIfExpired <- function( go.con ) {
  # If database connection go.con has expired returns a new connection to Gene
  # Ontology, retruns argument go.con otherwise.
  #
  # Args:
  #  go.con : A database connection as obtainable by connectToGeneOntology()
  #
  # Returns: An active database connection to the Gene Ontology database.
  #   
  if ( isConnectionAlive( go.con ) ) {
    go.con
  } else {
    connectToGeneOntology()
  }
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
