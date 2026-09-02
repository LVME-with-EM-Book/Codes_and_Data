################################################################################
F_Matrix2Tabular <- function(A, labelRowNames='', envir=TRUE){
   # begin
   if(envir){
      cat('\\begin{tabular}{')
      if(!is.null(rownames(A))){cat('l|')}
      cat(rep('r', ncol(A)), sep='')
      cat('} \n')
   }
   # headline
   if(!is.null(colnames(A))){
      if(!is.null(rownames(A))){cat(labelRowNames, '& ')}
      if(ncol(A)>1){for(j in 1:(ncol(A)-1)){cat(colnames(A)[j], '& ')}}
      cat(colnames(A)[ncol(A)], '\\\\ \n')
      cat('\\hline \n')
   }
   # lines
   for(i in 1:nrow(A)){
      if(!is.null(rownames(A))){cat(rownames(A)[i], '& ')}
      if(ncol(A) > 1){for(j in 1:(ncol(A)-1)){cat(as.character(A[i, j]), '& ')}}
      if(i<nrow(A)){cat(as.character(A[i, ncol(A)]), '\\\\ \n')}else{cat(as.character(A[i, ncol(A)]), '\n')}
   }
   # end
   if(envir){cat('\\end{tabular} \n')}
}
