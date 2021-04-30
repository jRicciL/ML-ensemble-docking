library(scmamp)
library(scales)

full_names <- c('LR', 'GBT', 'DClf', 'csAVG', 'csGEO', 'csMIN')
cbbPalette <- c( '#785EF0', '#3F93D2', '#44AA99', '#FE6100', '#DC267F', '#FFB000')
names(cbbPalette) <- full_names

getNemenyiCD <- function (alpha = 0.05, 
                          num.alg, 
                          num.problems) {
  #' Auxiliar function to compute the critical difference for Nemenyi test
  #' @Parms:
  #'   alpha:        Alpha for the test
  #'   num.alg:      Number of algorithms tested
  #'   num.problems: Number of problems where the algorithms have been tested
  #'
  #' @Returns:
  #'   Corresponding critical difference
  
  df <- num.alg * (num.problems - 1)
  qa <- qtukey(p=1 - alpha, nmeans=num.alg, df=df)/sqrt(2)
  cd <- qa * sqrt((num.alg * (num.alg + 1)) / (6 * num.problems))
  return(cd)
}

plot_Critical_Difference_color <- function (results.matrix, alpha=0.05, cex=0.75, colPalette=cbbPalette, 
            side_marging=0, line.spacing=0.9, labels.cex=0.6, ...) {
  
  opar <- par(mai = c(0,0,0,0))
  on.exit(par(opar))
  
  k <- dim(results.matrix)[2]
  N <- dim(results.matrix)[1]
  cd <- getNemenyiCD(alpha=alpha, num.alg=k, num.problems=N)
    
  ## Get color palette
#   colPalette <- hue_pal()(k)
#   names(cbbPalette) <- rev(colnames(results.matrix))
  
  mean.rank <- sort(colMeans(rankMatrix(results.matrix, ...)))
  
  # Separate the algorithms in left and right parts
  lp <- floor(k/2)
  left.algs <- mean.rank[1:lp]
  left.colors <- colPalette[names(left.algs)]
  right.algs <- mean.rank[(lp+1):k]
  right.colors <- colPalette[names(right.algs)]
  max.rows <- ceiling(k/2)
  
  # Basic dimensions and definitions
  char.size    <- 0.001  # Character size
  line.spacing <- line.spacing   # Line spacing for the algorithm name
  m            <- floor(min(mean.rank))
  M            <- ceiling(max(mean.rank))
  max.char     <- max(sapply(colnames(results.matrix), FUN = nchar))  # Longest length of a label
  text.width   <- (max.char + 8) * char.size
  w            <- (M-m) + 2 * text.width + side_marging
  h.up         <- 2.5 * line.spacing  # The upper part is fixed. Extra space is for the CD
  h.down       <- (max.rows + 2.25) * line.spacing # The lower part depends on the no. of algorithms. 
  # The 2 extra spaces are for the lines that join algorithms
  tick.h       <- 0.25 * line.spacing 
  
  label.displacement <- 1    # Displacement of the label with respect to the axis
  line.displacement  <- 0.025  # Displacement for the lines that join algorithms
  
  # Background of the plot
  plot(0, 0, type="n", xlim=c(m - w / (M - m), M + w / (M - m)), 
       ylim=c(-h.down, h.up), xaxt="n", yaxt="n", xlab= "", ylab="", bty="n")
  
  # Draw the axis
  lines (c(m,M), c(0,0), lwd=2)
  dk <- sapply(m:M, 
               FUN=function(x) {
                 lines(c(x,x), c(0, tick.h))
                 text(x, 3*tick.h, labels=x, cex=cex)
               })
  
  # Draw the critical difference
  lines(c(m, m + cd), c(1.75 * line.spacing, 1.75 * line.spacing), lwd=1)
  text(m + cd / 2, 2.5 * line.spacing, paste("CD:", round(cd,2)), cex=cex)
  lines(c(m, m), c(1.75 * line.spacing - tick.h, 
                   1.75 * line.spacing + tick.h))
  lines(c(m + cd, m + cd), c(1.75 * line.spacing - tick.h, 
                             1.75 * line.spacing + tick.h))
  
  # Left part, labels
  dk <- sapply (1:length(left.algs), 
                FUN=function(x) {
                  line.h <- -line.spacing * (x + 2)
                  text(x=m - label.displacement, y=line.h, 
                       labels=names(left.algs)[x], cex=labels.cex, adj=1)
                  lines(c(m - label.displacement*0.75, left.algs[x]), 
                        c(line.h, line.h), col=left.colors[x], lwd=2)
                  lines(c(left.algs[x], left.algs[x]), c(line.h, 0),
                       col=left.colors[x], lwd=2)
                })
  
  # Right part, labels
  dk <- sapply (1:length(right.algs), 
                FUN=function(x) {
                  line.h <- -line.spacing * (x + 2)
                  text(x=M + label.displacement, y=line.h, 
                       labels=names(right.algs)[x], cex=labels.cex, adj=0)
                  lines(c(M + label.displacement*0.75, right.algs[x]), 
                        c(line.h, line.h), col=right.colors[x], lwd=2)
                  lines(c(right.algs[x], right.algs[x]), c(line.h, 0),
                       col=right.colors[x], lwd=2)
                })
  
  # Draw the lines to join algorithms
  getInterval <- function (x) {
    from <- mean.rank[x]
    diff <- mean.rank - from
    ls <- which(diff > 0 & diff < cd)
    if (length(ls) > 0) {
      c(from, mean.rank[max(ls)])
    }
  }
  
  intervals <- mapply (1:k, FUN=getInterval)
  aux <- do.call(rbind, intervals)
  if(NROW(aux) > 0) {
    # With this strategy, there can be intervals included into bigger ones
    # We remove them in a sequential way
    to.join <- aux[1,]
    if(nrow(aux) > 1) {  
      for (r in 2:nrow(aux)) {
        if (aux[r - 1, 2] < aux[r, 2]) {
          to.join <- rbind(to.join, aux[r, ])
        }
      }
    }

    row <- c(1)
    # Determine each line in which row will be displayed
    if (!is.matrix(to.join)) {  # To avoid treating vector separately
      to.join <- t(as.matrix(to.join))
    }
    nlines <- dim(to.join)[1]

    for(r in 1:nlines) {
      id <- which(to.join[r, 1] > to.join[, 2])
      if(length(id) == 0) {
        row <- c(row, tail(row, 1) + 1)
      } else {
        row <- c(row, min(row[id]))
      }
    }
    
    step <- max(row) / 2

    # Draw the line
    dk <- sapply (1:nlines, 
                  FUN = function(x) {
                    y <- -line.spacing * (0.5 + row[x] / step)
                    lines(c(to.join[x, 1] - line.displacement, 
                            to.join[x, 2] + line.displacement), 
                          c(y, y), lwd=4)
                  })
  }
}