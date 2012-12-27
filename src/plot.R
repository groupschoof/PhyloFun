as.df.for.box.plot <- function( p.mut.mtrxs.list, distance.column.index,
  p.column.index, lapply.funk=lapply ) {
  # Prepares mutation probability matrices in argument list 'p.mut.mtrxs.list'
  # for a boxplot with ggplot2.
  #
  # Example plot:
  #--------------
  # ggplot( as.df.for.box.plot( p.mut.mtrxs.list, 3, 1, mclapply ), aes(
  # factor( p.mutation.Sequence.Distance.grid ), min.Sequence.Distance ) ) +
  # geom_boxplot()
  #
  # Args:
  #  p.mut.mtrxs.list      : The list of mutation probability matrices as
  #                          result of multiple function calls of function
  #                          'gridPMutation'
  #  distance.column.index : The integer index of the column to be box plotted,
  #                          i.e. 'min.Sequence.Distance'.
  #  p.column.index        : The integer index of the column holding the
  #                          gridded mutation probabilities, i.e.
  #                          'p.mutation.Sequence.Distance.grid'.
  #  lapply.funk           : Set to 'mclapply' if parallel execution is wanted.
  #
  # Returns: A data.frame that can be fed into a ggplot call. See example above.
  #   
  as.data.frame( 
    do.call( 'rbind', 
      lapply.funk( p.mut.mtrxs.list, function( mtrx ) {
        mtrx[ , c( p.column.index, distance.column.index ), drop=F ]
      })
    )
  )
}
