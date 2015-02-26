#' Plot phylogenies with frequencies of unique sequences
#'
#' This function plots a phylogeny built from the unique sequences
#' together with a representation of the abundance (frequencies) of
#' these in the data.
#'
#' @param x an \linkS4class{obkData} object.
#' @param ncatmax an integer giving the number of categories of
#' haplotype frequency
#' @param colfun a function to define the colours
#' @return nothing, the results are printed on the current graphical device.
#' @author Emmanuel Paradis
#'
#' @details First, the unique sequences are extracted using the
#' function \code{haplotype} (in \pkg{pegas}). Then, a
#' neighbor-joigning (NJ) tree is built. The tree is plotted on the
#' left-hand side of the graph, and the frequencies of each unique
#' sequence is represented in two ways: coloured labels at the tips of
#' the tree (with the colour scale drawn at the top), and a horizontal
#' barplot on the right-hand side of the graph.
#'
#' The functions scan the data for the different genes and plots the
#' output successively (the user is asked to type enter)
#'
#' @examples
#' \dontrun{
#' require(OutbreakTools)
#' data(ToyOutbreak)
#' tree_okbData(ToyOutbreak)
#' }

tree_okbData <- function(x, ncatmax = 10, colfun = topo.colors)
{
    require(pegas)
    DNA <- x@dna@dna
    LOCI <- names(DNA)
    cat("okbData object has", length(LOCI), "genetic data set(s)\n")
    for (i in seq_along(LOCI)) {
        cat("Processing data from", LOCI[i], "\n")

        dna <- DNA[[i]]
        h <- haplotype(dna)
        nh <- nrow(h)
        attr(h, "index") <- lapply(attr(h, "index"), function(x) rownames(dna)[x])
        names(attr(h, "index")) <- rownames(h) <- paste0("uniqseqID", seq_len(nh))
        phy <- nj(dist.dna(h, "N", pairwise.deletion = TRUE))

        hfreq <- sapply(attr(h, "index"), length)
        rangehfreq <- range(hfreq)
        if (rangehfreq[2] < 10) ncatmax <- rangehfreq[2]
        sq <- seq(rangehfreq[1] - 1, rangehfreq[2], length.out = ncatmax)
        cat <- cut(hfreq, sq)
        co <- colfun(ncatmax)
        n <- Ntip(phy)
        cat("Type enter to continue\n")
        readLines(n = 1)
        lttcoords <- ltt.plot.coords(phy, backward = FALSE)
        xlm <- lttcoords[nrow(lttcoords), 1] + max(hfreq)
        plot(phy, show.tip.label = FALSE, x.lim = xlm)
        axisPhylo()
        phydataplot(hfreq, phy, "b")
        tiplabels(text = "   ", bg = co[cat], adj = 0)
        mtext(LOCI[i], at = 0, font = 2)

        psr <- par("usr")
        xx <- psr[2]/2
        yy <- psr[4] * (0.5 + 0.5/par("plt")[4])
        legend(xx, yy, legend = round(sq[-1], 1), pch = 22, pt.bg = co,
               pt.cex = 2, bty = "n", xjust = 0.5, yjust = 0.5,
               horiz = TRUE, xpd = TRUE)
    }
}
