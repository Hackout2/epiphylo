#' Plot Phylogenies and Time
#'
#' These two functions plot a tree built from the sequences available
#' in the \linkS4class{obkData} object together with the time of
#' sampling of each sequence (also taken from the same object).
#'
#' @param x an \linkS4class{obkData} object.
#' @return nothing, the results are plotted on the current graphical device.
#' @author Emmanuel Paradis
#'
#' @details The first function \code{tree_time_1} plots the NJ tree
#' growing upwards, without the branch length information, and the
#' time axis is represented on top with from left to right. The second
#' function \code{tree_time_2} plots the NJ tree in the usual way,
#' from left to right, with the branch length information, and the
#' time axis is represented on bottom also from left to right.
#'
#' Coloured lines show the link between the tips of the tree (i.e.,
#' unique sequences) and the time of sampling of the individuals. The
#' colour scale goes from blue (oldest dates) to red (youngest dates).
#'
#' The functions scan the data for the different genes and plots the
#' output successively (the user is asked to type enter).
#'
#' @note The neighbor-joigning (NJ) method reconstructs unrooted trees
#' whereas the present functions represent them as rooted.
#'
#' @examples
#' \dontrun{
#' require(OutbreakTools)
#' data(ToyOutbreak)
#' tree_time_1(ToyOutbreak)
#' tree_time_2(ToyOutbreak)
#' }

tree_time_1 <- function(x)
{
    require(pegas)
    DNA <- x@dna@dna
    LOCI <- names(DNA)
    META <- x@dna@meta

    cat("okbData object has", length(LOCI), "genetic data set(s)\n")
    for (i in seq_along(LOCI)) {
        cat("Processing data from", LOCI[i], "\n")

        dna <- DNA[[i]]

        h <- haplotype(dna)
        nh <- nrow(h)
        attr(h, "index") <- lapply(attr(h, "index"), function(x) rownames(dna)[x])
        names(attr(h, "index")) <- rownames(h) <- paste0("uniqseqID", seq_len(nh))
        phy <- nj(dist.dna(h, "N", pairwise.deletion = TRUE))
        n <- Ntip(phy)

        sel <- META$locus == LOCI[i]
        submeta <- META[sel, ]
        dates <- as.integer(submeta$date)
        range.dates <- range(dates)

        footrans <- function(x) (n - 1) * (x - range.dates[1])/diff(range.dates) + 1
        transfdates <- footrans(dates)
        x1 <- transfdates
        names(x1) <- rownames(submeta)
        x2 <- x1
        cat("Type enter to plot\n")
        readLines(n = 1)

        plot(phy, "p", FALSE, show.tip.label = FALSE, direction = "u",
             font = 1, label.offset = 1, node.depth = 2, y.lim = 40)
        psr <- par("usr")
        nms <- rownames(submeta)
        index <- attr(h, "index")
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        for (j in 1:n) {
            s <- index[[phy$tip.label[j]]]
            x2[s] <- lastPP$xx[j]
        }
        y2 <- psr[4] - 0.1
        y1 <- max(lastPP$yy) + 0.1
        segments(x2, y1, x1, y2, lwd = 1, col = rgb(x1/n, 0, 1 - x1/n, alpha = .5))
        at <- pretty(dates)
        axis(3, at = footrans(at), labels = as.Date(at, origin = "1970-01-01"))
    }
}

#old code:
#
#nms <- rownames(submeta)
#index <- attr(h, "index")
#lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
#for (i in 1:n) {
#    s <- index[[phy$tip.label[i]]]
#    x2[s] <- lastPP$xx[i]
#}
#
#y1 <- psr[3]
#y2 <- 5
#segments(x1, y1, x2, y2, lwd = .5)
#at <- pretty(dates)
#axis(1, at = footrans(at), labels = as.Date(at, origin = "1970-01-01"))

# o <- order(dates)
# phy2 <- rotateConstr(phy, rownames(submeta)[o])


tree_time_2 <- function(x)
{
    require(pegas)
    DNA <- x@dna@dna
    LOCI <- names(DNA)
    META <- x@dna@meta

    cat("okbData object has", length(LOCI), "genetic data set(s)\n")
    for (i in seq_along(LOCI)) {
        cat("Processing data from", LOCI[i], "...\n")

        dna <- DNA[[i]]

        h <- haplotype(dna)
        nh <- nrow(h)
        attr(h, "index") <- lapply(attr(h, "index"), function(x) rownames(dna)[x])
        names(attr(h, "index")) <- rownames(h) <- paste0("uniqseqID", seq_len(nh))
        phy <- nj(dist.dna(h, "N", pairwise.deletion = TRUE))
        n <- Ntip(phy)

        sel <- META$locus == LOCI[i]
        submeta <- META[sel, ]
        dates <- as.integer(submeta$date)
        range.dates <- range(dates)

        cat("Type enter to plot\n")
        readLines(n = 1)

        plot(phy, "p", show.tip.label = FALSE, font = 1, label.offset = 1, node.depth = 2, y.lim = c(-100, n), edge.width = 2)
        psr <- par("usr")

        footrans2 <- function(x) psr[2] * (x - range.dates[1])/diff(range.dates)

        x1 <- footrans2(dates)

        y1 <- psr[3]
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        index <- attr(h, "index")
        x2 <- y2 <- numeric(length(dates))
        names(x2) <- names(y2) <- rownames(submeta)
        for (j in 1:n) {
            s <- index[[phy$tip.label[j]]]
            x2[s] <- lastPP$xx[j]
            y2[s] <- lastPP$yy[j]
        }
        segments(x1, y1, x2, y2, lwd = 1, col = rgb(x1/max(x1), 0, 1 - x1/max(x1), alpha = .5))
        at <- pretty(dates)
        axis(1, at = footrans2(at), labels = as.Date(at, origin = "1970-01-01"))
        add.scale.bar(psr[2] - 3, psr[4] - 5, 1)
        mtext(LOCI[i])
    }
}
