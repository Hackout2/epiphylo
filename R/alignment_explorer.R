#' Extracts and plots information on DNA alignments
#'
#' This function extracts information from the DNA sequences in an
#' \linkS4class{obkData} object and prints or plots them.
#'
#' @param x an \linkS4class{obkData} object.
#' @return nothing, the results are printed on the current graphical device.
#' @author Emmanuel Paradis
#'
#' @details The functions scan the data for the different genes and
#' plots the output successively (the user is asked to type enter).
#'
#' @examples
#' \dontrun{
#' require(OutbreakTools)
#' data(ToyOutbreak)
#' alignment_explorer(ToyOutbreak)
#' }

alignment_explorer <- function(x)
{
    require(pegas)

    op <- options(warn = -1)
    on.exit(options(op))

    DNA <- x@dna@dna
    LOCI <- names(DNA)

    cat("okbData object has", length(LOCI), "genetic data set(s)\n")
    for (i in seq_along(LOCI)) {
        cat("Processing data from", LOCI[i], "\n")
        dna <- DNA[[i]]
        s <- seg.sites(dna)
        ss <- site.spectrum(dna)
        more.than.two.states <- 0L
        ## lines taken from pegas::site.spectrum()
        for (j in s) {
            bfi <- base.freq(dna[, j], freq = TRUE)
            bfi <- bfi[bfi > 0]
            if (length(bfi) > 2)
                more.than.two.states <- more.than.two.states + 1L
        }
        cat("Type enter to plot\n")
        readLines(n = 1)
        layout(matrix(c(1, 2, 1, 3), 2, 2))
        image(dna, xlab = LOCI[i])
        rug(s, -0.05, 3, 1)
        mtext("Tick marks show segregating sites", at = 0, adj = 0, line = 2, font = 3)
        plot(ss)

        BF <-  base.freq(dna, TRUE, TRUE)
        nb <- sum(BF)
        ROW <- c("ambiguous bases", "alignment gaps",
                  "A, C, G or T", "", "TOTAL")
        nambi <- sum(BF[5:15])
        ngaps <- BF[16]
        nacgt <- sum(BF[1:4])
        COL1 <- c(nambi, ngaps, nacgt, NA, nb)
        COL2 <- 100 * COL1/nb
        COL2 <- paste(COL2[], "%")
        COL1[4] <- COL2[4] <- ""
        COL1 <- paste(COL1, collapse = "\n")
        COL2 <- paste(COL2, collapse = "\n")

        plot(NA, type = "n", xlim = c(0, 100), ylim = c(0, 100), bty = "n",
             xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        rect(0, 0, 100, 100, col = "blue", border = NA)
        textcol <- "yellow"
        cx <- 1.2
        ft <- 2
        yy <- 60
        xx1 <- 10
        xx2 <- xx1 + strwidth(COL1[1], cex = cx, font = ft) + 20
        xx3 <- xx2 + strwidth(COL2[4], cex = cx, font = ft) + 15
        text(xx1, yy, paste(ROW, collapse = "\n"), adj = c(0, 0.5),
             col = textcol, font = ft, cex = cx)
        text(xx2, yy, COL1, adj = c(1, 0.5), col = textcol, font = ft, cex = cx)
        text(xx3, yy, COL2, adj = c(1, 0.5), col = textcol, font = ft, cex = cx)
        if (more.than.two.states)
            msg <- paste0(more.than.two.states, " sites with more than two states were ignored in the spectrum")
        text(xx1, 30, msg, adj = c(0, 0.5), col = textcol, font = ft, cex = cx)
    }
}
