REPORT.CSS <- 
'<style type="text/css">
  .table {
    min-width: 820px;
    width: 820px;
    border: 1px solid #e6e6e6;
    border-collapse: collapse;
    -moz-border-radius: 4px;
    -webkit-border-radius: 4px;
    border-radius: 4px;
    color: #171717;
    font-size: 9pt;
  }

  .table td, .table th, .tablescroll_head th {
    padding: 0.5em;
    font-size: 9pt;
    border-bottom: 1px solid #e6e6e6;
  }

  .table th, .tablescroll_head th {
    background-color: #e6e6e6;
  }

  .table td {
    border-right: 1px solid #e6e6e6;
  }
</style>'

JS.PHYLO.SVG.RENDER.TEMPLATE <-
  '<render>
    <parameters>
      <circular>
        <bufferRadius>0.5</bufferRadius>
      </circular>
      <rectangular>
        <alignRight>1</alignRight>
        <bufferX>300</bufferX>
      </rectangular>
    </parameters>
    <charts>
      <biological_process type="binary" thickness="10" />
      <cellular_component type="binary" thickness="10" />
      <molecular_function type="binary" thickness="10" />
    </charts>
    <styles>
      <query fill="#A93" stroke="#DDD" />
      <% anno.colors <- rainbow( sum( as.integer( lapply( names( go.anno.spaces ), function( go.type ) length( go.anno.spaces[[ go.type ]] ) ) ) ) ) %>
      <% i <- 1 %>
      <% for ( go.type in names( go.anno.spaces ) ) { %>
        <% for ( go.anno in go.anno.spaces[[ go.type ]] ) { %> 
          <<%= paste( go.type, i, sep="" ) %> fill="<%= anno.colors[[ i ]] %>" stroke="#DDD" />
          <% i <- i + 1 %>
        <% } %>
      <% } %>
    </styles>	
  </render>'

LEGEND.TEMPLATE <- 
'<div id="legend">
  <table id="go_annotation_legend">
    <tr><th>GO accession</th><th>GO name</th><th>distance to root</th></tr>
    <% i <- 1 %>
    <% for ( go.type in names( go.anno.spaces ) ) { %>
      <tr><th colspan="3"><%= go.type %></th></tr>
      <% for ( go.anno in go.anno.spaces[[ go.type ]] ) { %> 
        <% for ( go.acc in go.anno ) { %>
          <% go.term <- go.terms.tbl[ which( go.terms.tbl[ , "acc" ] == go.acc ), ] %>
          <tr style="background-color:rgb(<%= paste( anno.colors[[ i ]], collapse="," ) %>);">
            <% if ( nrow( go.term ) > 0 ) { %>
              <td><%= go.term[[ "acc" ]] %></td>
              <td><%= go.term[[ "name" ]] %></td>
              <td><%= go.term[[ "relation_distance" ]] %></td>
            <% } else { %>
              <td colspan="3">unknown</td>
            <% } %>
          </tr>
        <% } %>
        <% i <- i + 1 %>
      <% } %>
    <% } %>
  </table>
</div>'

PHYLO.XML.TEMPLATE <- 
'<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns="http://www.phyloxml.org">
  <phylogeny rooted="false">
    <%= render.tag %>
    <%= convertToJsPhyloXML( phylo.tree, annotation.matrix, query.accession ) %>
  </phylogeny>
</phyloxml>'

PHYLO.XML.CLADE.TEMPLATE <- 
'<clade>
  <% if ( confidence != "" ) { %>
    <confidence type="bootstrap"><%= confidence %></confidence>
  <% } %>
  <% if ( branch.length != "" ) { %>
    <branch_length><%= branch.length %></branch_length>
  <% } %>
  <% if ( is.leaf.node ) { %>
    <% prot.acc <- phylo.tree$tip.label[[ node.index ]] %>
    <name <%= if ( is.query.node ) \'bgStyle="query"\' %>><%= prot.acc %></name>
    <% if ( ! is.query.node ) { %>
      <annotation>
        <% if ( ! is.query.node ) { %>
          <uri>http://www.uniprot.org/uniprot/<%= prot.acc %></uri>
        <% } %>
        <desc><%= popUpAnnotations( go.terms.tbl, annotation.matrix, prot.acc )%></desc>
      </annotation>
    <% } %>
 </clade>
  <% } %>'

PHYLO.XML.NAMESPACE <- c( xmlns="http://www.phyloxml.org" )

annoColors <- function( go.anno.spaces ) {
  setNames(
    as.character(
      lapply(
        rainbow(
          sum( as.integer(
            lapply( names( go.anno.spaces ),
              function( go.type ) length( go.anno.spaces[[ go.type ]] )
            )
          ) )
        ),
        col2rgb
      )
    ),
    as.character( unlist( go.anno.spaces, recursive=F ) )
  )
}

goTermsTable <- function( go.anno.spaces, go.con=connectToGeneOntology() ) {
  goTermsForAccessionWithLevel(
    setdiff(
      unique( as.character( unlist( go.anno.spaces ) ) ),
      'unknown'
    ),
    go.con=go.con
  )
}

popUpAnnotations <- function( go.terms.tbl, annotation.matrix, prot.acc ) {
  prot.annos <- annotation.matrix[[ 'GO', prot.acc ]]
  paste( as.character(
      lapply( prot.annos, function( go.acc ) {
        gt <- go.terms.tbl[ which( go.terms.tbl[ , 'acc' ] == go.acc ), , drop=F ]
        if ( nrow( gt ) > 0 ) {
          paste( go.acc, '(', gt[[ 1, 'name' ]], ')' )
        }
      })
    ),
    collapse='\r'
  )
}

renderHTMLReport <- function( phylo.tree, annotation.matrix, annotation.space,
  query.accession ) {
  
}

renderHeader <- function( annotation.space ) {
  
}

convertToJsPhyloXML <- function( phylo.tree,
  annotation.matrix, annotation.space, query.accession, 
  go.terms.tbl, node.index=as.integer( get.root.node( phylo.tree  ) ) ) {
  desc.nds <- getDescendantNodes( phylo.tree, node.index )
  is.leaf.node <- is.null( desc.nds )
  is.root.node <- node.index == get.root.node( phylo.tree ) 
  is.query.node <- is.leaf.node &&
    phylo.tree$tip.label[[ node.index ]] == query.accession
  edge.ind <- as.integer( which( phylo.tree$edge[ , 2 ] == node.index ) )
  confidence <- if ( ! is.leaf.node ) {
    phylo.tree$node.label[[
      node.index - length( phylo.tree$tip.label )
    ]]
  } else {
    ''
  }
  branch.length <- if ( ! is.root.node ) {
    phylo.tree$edge.length[[ edge.ind ]]
  } else {
    ''
  }
  phylo.clade <- capture.output( brew( text=PHYLO.XML.CLADE.TEMPLATE ) )
  if ( ! is.leaf.node ) {
    phylo.xml <- c(
      phylo.clade,
      unlist(
        lapply( desc.nds, function( desc.node ) {
          convertToJsPhyloXML( phylo.tree, annotation.matrix, annotation.space,
            query.accession, go.terms.tbl, desc.node
          )
        }),
        recursive=F
      ),
      '</clade>'
    )
    paste( phylo.xml, collapse='' )
  } else {
    phylo.clade
  }
}


