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
  join.funk <- function( ... ) { paste( ..., sep=", " ) }
  dbGetQuery( con, paste(
      "SELECT t.*, g.relation_distance FROM term t LEFT JOIN graph_path g ON ",
      "t.id = g.term2_id WHERE g.term1_id = ( SELECT r.id FROM term r WHERE r.is_root = 1 ) ",
      "AND t.acc in (",
      do.call( 'join.funk', as.list( paste( "'", accessions, "'", sep="" ) ) ),
      ")"
    ) 
  )
}
