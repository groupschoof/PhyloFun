#ifndef _PhyloFun_RCPP_PHYLO_FUN_H
#define _PhyloFun_RCPP_PHYLO_FUN_H

#include <Rcpp.h>

using namespace Rcpp ;

/*
 * Function returns the row of matrix 'table' in which column 'colInd's value
 * is >= argument 'value'. Exported to be callable from R.
 */
RcppExport SEXP findMatchingRow( SEXP sTable, SEXP val, SEXP colInd );

/*
 * Finds the maximum mutation probability for any atomic annotation in argument
 * 'compositeAnnotation'. Exported to be callable from R.
 */
RcppExport SEXP mutationProbability( SEXP compositeAnnotation, SEXP
    branchLength, SEXP annotsMutProbTables, SEXP distanceColumnIndx );

/*
 * Generates the conditional probability table which holds the mutation
 * probability for each composite annotation in argument 'annos' to each
 * composite annotation. The table is returned as a NumericMatrix with row and
 * column names being the respective stringified annotations as in argument
 * 'stringifiedAnnotations'. Exported to be callable from R.
 */
RcppExport SEXP conditionalProbabilityTable( SEXP branchLength, SEXP annos,
    SEXP stringifiedAnnotations, SEXP annotsMutationProbTables, SEXP
    mutTblLengthColIndx );

/*
 * Generates one mutation probability table for each unique phylogenetic branch
 * length in argument 'uniqueEdgeLengths'. Invokes function
 * conditionalProbabilityTable with each such branch length. Returns a named
 * list of numeric matrices, in which the names are the unique branch lengths
 * as strings. Exported to be callable from R.
 */
RcppExport SEXP conditionalProbabilityTables( SEXP uniqueEdgeLengths, SEXP
    annos, SEXP stringifiedAnnotations, SEXP annotsMutationProbTableList, SEXP
    mutTblLengthColIndx );

/*
 * Iterates through the argument proteinAccessionMatrix looking up rows in
 * which the accession stored in column 'pairFirstMemberColIndex' equals
 * argument 'accession'. If the row is a self match, that is the identity pair,
 * it is excluded from the returned list of pairs ( StringVectors of length 2
 * ).
 */
RcppExport SEXP extractProteinPairs( SEXP proteinAccessionMatrix, SEXP
    accession, SEXP pairFirstMemberColIndex, SEXP pairSecondMemberColIndex );

/*
 * Iterates through the argument proteinAccessionMatrix looking up rows in
 * which the accession stored in column 'pairFirstMemberColIndex' equals
 * argument 'accession'. If the row is a self match, that is the identity pair,
 * it is not counted. Instead of returning a list of protein pairs this
 * function just counts the number of matching pairs.
 */
RcppExport SEXP countProteinPairs( SEXP proteinAccessionMatrix, SEXP
    accession, SEXP pairFirstMemberColIndex, SEXP pairSecondMemberColIndex );

/*
 * Iterates over the vector 'accessions' and invokes countProteinPairs(…) or
 * extractProteinPairs(…), respectively. Counts are measured if and only if
 * argument 'funcAbbrev' equals "count". Returns a named list of results from
 * the respective function calls.
 */
RcppExport SEXP countOrExtractProteinPairs( SEXP proteinAccessionMatrix, SEXP
    accessions, SEXP pairFirstMemberColIndex, SEXP pairSecondMemberColIndex,
    SEXP funcAbbrev );

#endif
