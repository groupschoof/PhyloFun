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

SEXP extractProteinPairs( SEXP proteinAccessionMatrix, SEXP accession, SEXP
    pairFirstMemberColIndex, SEXP pairSecondMemberColIndex ) {
  StringMatrix tbl( proteinAccessionMatrix );
  StringVector rAcc( accession );
  std::string acc( rAcc( 0 ) );
  NumericVector rColA( pairFirstMemberColIndex );
  int colA( rColA( 0 ) );
  NumericVector rColB( pairSecondMemberColIndex );
  int colB( rColB( 0 ) );

  List prs( 0 );
  for ( int i=0; i<tbl.nrow(); ++i ) {
    if ( acc == std::string( tbl( i, 0 ) ) &&
      acc != std::string( tbl( i, 1 ) ) )
    {
      prs.push_back( tbl( i, _ ) );
    }
  }
  return( wrap( prs ) );
}

SEXP countProteinPairs( SEXP proteinAccessionMatrix, SEXP accession, SEXP
    pairFirstMemberColIndex, SEXP pairSecondMemberColIndex ) {
  StringMatrix tbl( proteinAccessionMatrix );
  StringVector rAcc( accession );
  std::string acc( rAcc( 0 ) );
  NumericVector rColA( pairFirstMemberColIndex );
  int colA( rColA( 0 ) );
  NumericVector rColB( pairSecondMemberColIndex );
  int colB( rColB( 0 ) );

  int cnt = 0;
  for ( int i=0; i<tbl.nrow(); ++i ) {
    if ( acc == std::string( tbl( i, 0 ) ) &&
      acc != std::string( tbl( i, 1 ) ) )
    {
      ++cnt;
    }
  }
  return( wrap( cnt ) );
}

SEXP countOrExtractProteinPairs( SEXP proteinAccessionMatrix, SEXP accessions,
    SEXP pairFirstMemberColIndex, SEXP pairSecondMemberColIndex, SEXP
    funcAbbrev ) {
  CharacterVector rFuncAbbrev( funcAbbrev );
  bool countPairs = std::string( rFuncAbbrev( 0 ) ) == "count";
  CharacterVector protAccs( accessions );

  List rslts( 0 );
  for ( int i = 0; i < protAccs.size(); ++i ) {
    std::string protAcc( protAccs( i ) );
    if ( countPairs ) {
      rslts.push_back( countProteinPairs( proteinAccessionMatrix, protAccs( i ),
            pairFirstMemberColIndex, pairSecondMemberColIndex ), protAcc );
    } else {
      rslts.push_back( extractProteinPairs( proteinAccessionMatrix, protAccs( i ),
            pairFirstMemberColIndex, pairSecondMemberColIndex ), protAcc );
    }
  }

  return( wrap( rslts ) );
}