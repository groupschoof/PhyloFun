#include "rcpp_phylo_fun.h"
using namespace Rcpp ;

SEXP findMatchingRow( SEXP sTable, SEXP val, SEXP colInd ){
  BEGIN_RCPP

    NumericMatrix table( sTable );
    NumericVector value = NumericVector( val );
    NumericVector columnIndex = NumericVector( colInd );
    NumericVector clmn = table( _ , columnIndex( 0 ) );
    int indx = clmn.size() - 1;
    for( int i=0; i<clmn.size(); ++i ) {
      if( clmn( i ) >= value( 0 ) ) {
        indx = i;
        break;
      }
    }
    NumericVector rw = table( indx, _ );
    return( wrap( rw ) );

  END_RCPP
}

SEXP mutationProbability( SEXP compositeAnnotation, SEXP branchLength, SEXP
    annotsMutProbTables, SEXP distanceColumnIndx ) {
  BEGIN_RCPP

    CharacterVector compAnnos( compositeAnnotation );
    List anMutProbTbls = List( annotsMutProbTables );
    double mutProb = 0.0;
    for ( int i = 0; i < compAnnos.size(); ++i ) {
      std::string singlAnno( compAnnos( i ) );
      NumericMatrix pMutTbl = anMutProbTbls( singlAnno );
      NumericVector pMutRow = findMatchingRow( pMutTbl, branchLength, distanceColumnIndx );
      if ( pMutRow( 0 ) > mutProb ) {
        mutProb = pMutRow( 0 );
      }
    }
    return( wrap( NumericVector( 1, mutProb ) ) );

  END_RCPP
}

SEXP conditionalProbabilityTable( SEXP branchLength, SEXP annos,
    SEXP stringifiedAnnotations, SEXP annotsMutationProbTables, SEXP
    mutTblLengthColIndx ) {
  BEGIN_RCPP

    List annosLst( annos );
    CharacterVector annosAsStrs( stringifiedAnnotations );
    NumericMatrix cpt = NumericMatrix( annosLst.size(), annosLst.size() );
    cpt.attr( "dimnames" ) = List::create( annosAsStrs, annosAsStrs );

    for ( int i = 0; i < annosLst.size(); ++i ) {
      CharacterVector compositeAnnotation = annosLst( i );
      double compAnnoMutProb = 1.0;
      std::string ua = "unknown";
      std::string caFirst = as<std::string>( compositeAnnotation( 0 ) );
      if ( ua != caFirst ) {
        compAnnoMutProb = ( (NumericVector) mutationProbability(
              compositeAnnotation, branchLength, annotsMutationProbTables,
              mutTblLengthColIndx ) )( 0 );
      }
      double mutToOtherAnnosProb = compAnnoMutProb / ( annosLst.size() - 1 );
      NumericVector colmn( annosLst.size(), mutToOtherAnnosProb );
      colmn( i ) = 1.0 - compAnnoMutProb;
      cpt( _, i ) = colmn;
    }

    return( wrap( cpt ) );

  END_RCPP
}

SEXP conditionalProbabilityTables( SEXP uniqueEdgeLengths, SEXP annos, SEXP
    stringifiedAnnotations, SEXP annotsMutationProbTableList, SEXP
    mutTblLengthColIndx ) {
  BEGIN_RCPP

    NumericVector edgeLengths = NumericVector( uniqueEdgeLengths );
    CharacterVector edgeLengthsAsStrs = as<CharacterVector>( edgeLengths );
    List cpts = List( 0 );

    for ( int i = 0; i < edgeLengths.size(); ++i ) {
      NumericVector currBranchLength( 1, edgeLengths( i ) );
      NumericMatrix cpt = conditionalProbabilityTable( currBranchLength, annos,
          stringifiedAnnotations, annotsMutationProbTableList,
          mutTblLengthColIndx );
      cpts.push_back( cpt, std::string( edgeLengthsAsStrs( i ) ) );
    }

    return( wrap( cpts ) );

  END_RCPP
}

SEXP proteinsWithGOtermAnnotation( SEXP goTermAccession, SEXP
    proteinGOannotationMatrix, SEXP proteinAccessionColIndx, SEXP goTermColInd
    ) {
  BEGIN_RCPP

    // Initializations:
    std::string goAcc( CharacterVector( goTermAccession )( 0 ) );
    
    int colProt( NumericVector( proteinAccessionColIndx )( 0 ) );
    int colGO( NumericVector( goTermColInd )( 0 ) );
    CharacterMatrix protGOannos( proteinGOannotationMatrix );
    CharacterVector annotatedGOs( protGOannos( _, colGO ) );
    CharacterVector annotatedProts( protGOannos( _, colProt ) );
    CharacterVector matchingProts;

    // Find accessions of proteins annotated with the GO term
    // 'goTermAccession':
    for (int i = 0; i < annotatedGOs.size(); ++i) {
      if ( std::string( annotatedGOs( i ) ) == goAcc ) {
        matchingProts.push_back( annotatedProts( i ) );
      }
    }

    return( wrap( matchingProts ) );

  END_RCPP
}

// Functor used for the creation of a set of _unordered_ vectors of
// std::strings and length 2. Used in function uniqueProteinPairs(…)
struct Comparator {
  bool operator()(const std::vector<std::string> & a, const std::vector<std::string> & b) {
    const bool swapA = a[0] < a[1];
    const std::string & al = swapA ? a[0] : a[1];
    const std::string & ar = swapA ? a[1] : a[0];
    const bool swapB = b[0] < b[1];
    const std::string & bl = swapB ? b[0] : b[1];
    const std::string & br = swapB ? b[1] : b[0];
    return al < bl || (al == bl && ar < br);
  }
};

// Extracts the subset of unique protein pairs from argument 'proteinPairsTbl'.
// Here each row is considered a symmetric pair of strings. Uses above functor
// to construct an un-ordered set.
SEXP uniqueProteinPairs( SEXP proteinPairsTbl, SEXP pairFirstMemberColIndex,
    SEXP pairSecondMemberColIndex ) {
  BEGIN_RCPP

    StringMatrix tbl( proteinPairsTbl );
    int colA( NumericVector( pairFirstMemberColIndex )( 0 ) );
    int colB( NumericVector( pairSecondMemberColIndex )( 0 ) );
    std::set<std::vector<std::string>, Comparator> prs;

    for ( int i=0; i<tbl.nrow(); ++i ) {
      CharacterVector r( tbl( i, _ ) );
      std::vector< std::string > p;
      p.push_back( std::string( r( colA ) ) );
      p.push_back( std::string( r( colB ) ) );
      prs.insert( p );
    }

    return( wrap( prs ) );

  END_RCPP
}
