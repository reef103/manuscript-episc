# Functions required for the manuscript code

#' Simple one-tail rank based enrichment analysis sREA
#'
#' This function performs simple 1-tail rank based enrichment analysis
#'
#' @param signatures Numeric matrix of signatures
#' @param groups List containing the groups as vectors of sample names
#' @param scaled Logical, whether the enrichment score should be scaled to the possible maximum
#' @return Matrix of Normalized Enrichment Zcores
#' @export
sREA <- function(signatures, groups, scaled=FALSE) {
    if (is.null(nrow(signatures))) signatures <- matrix(signatures, length(signatures), 1, dimnames=list(names(signatures), "sample1"))
    sig <- qnorm(apply(signatures, 2, rank)/(nrow(signatures)+1))
    gr <- sapply(groups, function(x, samp) {
        samp %in% x
    }, samp=rownames(sig))
    gr <- t(gr)
    nn <- rowSums(gr)
    gr <- gr/nn
    es <- gr %*% sig
    if (scaled) {
        gr1 <- t(apply(gr, 1, sort))
        es1 <- as.vector(gr1 %*% matrix(sort(sig[, 1]), nrow(sig), 1))
        return(es/es1)
    }
    return(es*sqrt(nn))
}

#' Z-score Stouffer's integration with Brown's correction
#'
#' This function integrate the z-scores per row using the Stouffer method
#'
#' @param x Matrix of Z-scores
#' @param brown Character string indicating the correlation method used for the Brown's correction, either none, pearson, spearman or kendall
#' @return Vector of integrated scores
#' @export

rowStouffer <- function(x, brown=c("none", "pearson", "spearman", "kendall")) {
    brown <- match.arg(brown)
    if (brown=="none") bc <- 0
    else {
        tmp <- abs(cor(x, method=brown, use="pairwise.complete.obs"))
        bc <- mean(tmp[upper.tri(tmp)], na.rm=TRUE)
        bc[is.na(bc)] <- 0
    }
    xsum <- rowSums(x, na.rm=TRUE)
    wsum <- ncol(x)
    bc * xsum / wsum + (1-bc) * xsum / sqrt(wsum)
}

# Plotting functions
plothm. <- function(x, color=c("royalblue","firebrick2"), gama=1, grid=T, scmax=0, box=TRUE, ...) {
    coli <- colorScale(x=filterRowMatrix(x, nrow(x):1), color=color, gama=gama, scmax=scmax)
    image(1:ncol(x), 1:nrow(x), t(matrix(1:(ncol(x)*nrow(x)), nrow(x), ncol(x))), col=coli, ylab="", xlab="", axes=F, ...)
    if (box) box()
    if (grid) grid(ncol(x), nrow(x), col="lightgrey", lty=1)
}

#' colorScale
#' 
#' This function generates a color scale
#' 
#' @param x Vector or matrix of numeric values
#' @param color Vector of character strings indicating the colors for the scale. Up to three colors can be defined. While is used for the missing color
#' @param gama Number indicating the gama transformation
#' @param alpha Number between 0 and 1 indicating the transparency of the color (1 for absolute color)
#' @param scmax Number indicating the maximum value for the scale
#' @param nacol Character string indicating the color for missing values
#' @return Vector of colors
#' @export
colorScale <- function(x, color=c("royalblue","firebrick2"), gama=1, alpha=1, scmax=0, nacol="grey80") {
    if (length(color)==1) color <- c(color, "white", color)
    if (length(color)==2) color <- c(color[1], "white", color[2])
    if (scmax==0) scmax <- max(abs(x), na.rm=T)
    pos <- which(abs(x) > scmax)
    if (length(pos)>0) x[pos] <- scmax*sign(x[pos])
    x <- abs(x/scmax)^gama*sign(x)
    color <- t(col2rgb(color))
    col <- sapply(x, function(x, color) {
        colSums(color*c(abs(x)*(x<0), 1-abs(x), x*(x>0)))
    }, color=color/255)
    pos <- which(colSums(is.na(col))>0)
    col[is.na(col)] <- 0
    col <- apply(col, 2, function(x, alpha) rgb(x[1], x[2], x[3], alpha=alpha), alpha=alpha)
    col[pos] <- nacol
    return(col)
}

#' Plot heatmap
#'
#' This function produce a heatmap plot from a numerical matrix
#'
#' @param x Numerical matrix
#' @param color Two character strings vector describing the colors for the heatmap
#' @param gama Number, indicating the exponential transformation for the color scale
#' @param cex Number indicating the magnification factor for the labels
#' @param grid Logical, whether a grid should be ploted
#' @param scale Number between 0 and .9 indicating the proportion of vertical space used to draw the color scale
#' @param scmax Optional number indicating the maximum value to be allowed for the heatmap
#' @param box Logical, whether to draw a box around the plot
#' @param ... Additional parameters to pass to the plot function
#' @return Nothing, a heatmap is produced in the default output device
#' @export

plothm <- function(x, color=c("royalblue","firebrick2"), gama=1, cex=1, grid=T, scale=F, scmax=0, box=TRUE, ...) {
    if (scale>0) {
        if (scale==1) ff <- 6/(nrow(x)+5)
        else ff <- scale
        pari <- par("mai")
        layout(matrix(1:2, 2, 1), h=c(1-ff, ff))
        if (round(sum(pari-c(1.02, .82, .82, .42)), 2)==0) pari <- c(.2, .2, 1.2, 1.2)
        par(mai=pari)
        plothm.(x, color=color, gama=gama, scmax=scmax, box=box, ...)
        axis(4, nrow(x):1, rownames(x), tick=F, line=0, las=2, adj=0, cex.axis=cex)
        axis(3, 1:ncol(x), colnames(x), tick=F, line=0, las=2, adj=0, cex.axis=cex)
        ra <- seq(-1, 1, length=100)
        coli <- colorScale(x=ra, color=color, gama=gama, scmax=scmax)
        par(mai=pari*c(0, 1, 0, 3)+c(.5, 0, .1, 0))
        image(1:length(ra), 1, matrix(1:length(ra), length(ra), 1), col = coli, ylab = "", xlab = "", axes = F)
        if (scmax==0) scmax <- max(abs(x), na.rm=T)
        axis(1, seq(1, length(ra), length=5), round(seq(-scmax, scmax, length=5), 1), cex.axis=cex)
    }
    else plothm.(x=x, color=color, gama=gama, grid=grid, scmax=scmax, box=box, ...)
}

#' Capitalize first letters
#' 
#' This function capitalizes first letter of words
#' 
#' @param x Character strig to capitalize
#' @param all Logical, whether all words should be capitalized
#' @return String
#' @export
capitalize <- function(x, all=FALSE) {
    if (all) x <- strsplit(x, " ")[[1]]
    paste(toupper(substring(x, 1,1)), substring(x, 2), sep="", collapse=" ")
}

#' Nice Exponential representations of scientific notation
#' 
#' This function generates a plotmath or latex representation of scientific notation
#' 
#' @param x Numeric vector
#' @param drop.1 Logical, whether 1 in 1 x type of representatons should be dropped
#' @param sub10 Either logical, "10", a non-negative integer or a length 2 integer vector, indicating if some expression should be formatted traditionally, when integer, all expression before the integer are simplified. when a 2 elements vector, all between the indicated range are simplified
#' @param digits Number of significant digits
#' @param lab.type Character string indicating how the result should look like, either plotmath or latex
#' @param lab.sep Character separator between mantissa and exponent
#' @return Vector of formated numbers
#' @export
niceExponent <- function(x, drop.1 = TRUE, sub10 = "10", digits = 2, digits.fuzz, lab.type = c("plotmath", "latex"), lab.sep = c("cdot", "times"))
{
    lab.type <- match.arg(lab.type)
    lab.sep <- match.arg(lab.sep)
    eT <- floor(log10(abs(x)) + 10^-digits)
    mT <- signif(x / 10^eT, digits)
    ss <- vector("list", length(x))
    if(sub.10 <- !identical(sub10, FALSE)) {
        if(identical(sub10, TRUE))
            sub10 <- c(0,0)
        else if(identical(sub10, "10"))
            sub10 <- 0:1
        sub10 <- as.integer(sub10)
        noE <-
            if(length(sub10) == 1) {
                if(sub10 < 0)
                    stop("'sub10' must not be negative if a single number")
                eT <= sub10
            } else if(length(sub10) == 2) {
                stopifnot(sub10[1] <= sub10[2])
                sub10[1] <= eT & eT <= sub10[2]
            } else stop("invalid 'sub10'")
        mT[noE] <- mT[noE] * 10^eT[noE]
    }
    if (lab.type == "plotmath") {
        for(i in seq(along = x))
            ss[[i]] <-
                if(x[i] == 0) quote(0)
        else if(sub.10 &&  noE[i]    ) substitute( A, list(A = mT[i]))
        else if(drop.1 && mT[i] ==  1) substitute( 10^E, list(E = eT[i]))
        else if(drop.1 && mT[i] == -1) substitute(-10^E, list(E = eT[i]))
        else substitute(A %*% 10^E, list(A = mT[i], E = eT[i]))
        do.call("expression", ss)
    }
    else {
        mTf <- format(mT)
        eTf <- format(eT)
        for(i in seq(along = x))
            ss[[i]] <-
            if(x[i] == 0) ""
        else if(sub.10 &&  noE[i]    ) mTf[i]
        else if(drop.1 && mT[i] ==  1) sprintf("$10^{%s}$", eTf[i])
        else if(drop.1 && mT[i] == -1) sprintf("$-10^{%s}$",eTf[i])
        else sprintf("$%s \\%s 10^{%s}$", mTf[i], lab.sep,  eTf[i])
        unlist(ss, use.names=FALSE)  ## perhaps unlist(ss) ?
    }
}

#' Concatenate text
#' 
#' This function concatenates text using , and and
#' 
#' @param x Vector of character strings
#' @return Character string
#' @export
textConcatenate <- function(x) {
    if (length(x)<2) return(x)
    if (length(x)==2) {
        tmp <- paste(x, collapse=" and ")
        return(tmp)
    }
    tmp <- paste(paste(x[1:(length(x)-1)], collapse=", "), x[length(x)], sep=" and ")
    return(tmp)
}


#' Variance stabilization transformation for RNAseq data
#'
#' This function stabilizes the variance, transform the data and add shot noise to the data
#'
#' @param x CountDataSet or matrix containing the raw counts, with genes in rows and samples in columns
#' @param method Character string indicating the method for estimating the dispersion (see DESeq::estimateDispersions)
#' @param fitType Character string indicating the type of fit for the dispersion (see DESeq::estimateDispersions)
#' @param seed Integer indicating the fixed seed for random numbers, 0 for not setting the seed
#' @return Expression matrix
#' @export
DEtransform <- function(x, method=c("blind", "pooled", "pooled-CR", "per-condition"), fitType=c("parametric", "local"), noise=TRUE, seed=1) {
    if (seed>0) set.seed(seed)
    method <- match.arg(method)
    fitType <- match.arg(fitType)
    cnames <- NULL
    if (!("CountDataSet" %in% class(x))) {
        if (!("matrix" %in% class(x))) stop("x must be a CountDataSet or integer matrix object", call.=F)
        if (length(which(duplicated(colnames(x))))>0) {
            cnames <- colnames(x)
            colnames(x) <- 1:ncol(x)
        }
        x <- newCountDataSet(x, factor(colnames(x)))
    }
    x <- estimateSizeFactors(x)
    x <- estimateDispersions(x, method=method, fitType=fitType)
    x <- getVarianceStabilizedData(x)
    tmp <- x
    if (noise) {
        tmp <- unlist(apply(x, 2, function(x) {
            x <- sort(unique(x))
            x <- cbind(x[1:(length(x)-1)], x[2:length(x)])
            x <- cbind(x[, 1], sqrt(frvarna(x)))
            return(list(x))
        }), recursive=FALSE)
        tmp <- cbind(unlist(lapply(tmp, function(x) x[, 1]), use.names=F), unlist(lapply(tmp, function(x) x[, 2]), use.names=F))
        tmp1 <- smooth.spline(tmp[, 1], tmp[, 2], spar=.5)
        tmp[tmp[, 1]>tmp1$x[which.min(tmp1$y)], 2] <- 0
        tmp1 <- smooth.spline(tmp[, 1], tmp[, 2], spar=.5)
        tmp <- x+rnorm(length(x))*predict(tmp1, x)$y
    }
    if (any(!is.null(cnames))) colnames(tmp) <- cnames
    return(tmp)
}

#' Mean for groups
#'
#' This function computes the mean for groups by rows for \code{colnames}-defined groups
#'
#' @param dset Numeric matrix
#' @return Numeric matrix
#' @export
meanGroup <- function(dset) {
    groups <- colnames(dset)
    glev <- unique(groups)
    repit <- tapply(groups, groups, length)
    res <- dset[,match(names(repit)[repit==1],colnames(dset))]
    for (i in names(repit)[repit>1]) res <- cbind(res,as.matrix(rowMeans(as.matrix(dset[,groups==i]))))
    res <- res[,match(glev,c(names(repit)[repit==1],names(repit)[repit>1]))]
    if (is.null(dim(res))) dim(res) <- c(length(res),1)
    colnames(res) <- glev
    rownames(res) <- rownames(dset)
    res
}

#' PCA plot for pilot perturbation experiment
#' 
#' @param pca List of PCA results with slots x and pvar
#' @annot annot Matrix of 7 columns indicating the metadata for the pilot experiment
plotPCApilotPerturbations <- function(pca, annot) {
    drugs <- sort(unique(annot[, 4]))
    drugs <- c("none", drugs[!(drugs %in% c("none", "CM"))])
    col <- c("black", rainbow(length(drugs)-1, s=.8, v=.9, start=0, end=.7))
    pt <- rep(1, nrow(annot))
    pt[annot[, 5]=="24"] <- 2
    pt[annot[, 5]=="36" & annot[, 3]=="RA"] <- 5
    pt1 <- pt+15
    pt1[pt1==20] <- 18
    pt.cex <- rep(1, length(pt1))
    pt.cex[pt1==18] <- 1.2
    par(mai=c(.8, .8, .2, .2))
    plot(pca$x[1, ], pca$x[2, ], xlab="", ylab="", pch=pt1, col=col[match(annot[, 4], drugs)], cex=pt.cex, xlim=c(min(pca$x[1, ]), max(pca$x[1, ]*1.7)), axes=FALSE)
    points(pca$x[1, ], pca$x[2, ], pch=pt)
    axis(2)
    axis(1, axisTicks(range(pca$x[1, ]*1.1), FALSE))
    axis(1, mean(range(pca$x[1, ])), paste0("PC1 - Var: ", round(pca$pvar[1]*100, 1), "%"), tick=FALSE, line=1.5)
    axis(2, mean(range(pca$x[2, ])), paste0("PC2 - Var: ", round(pca$pvar[2]*100, 1), "%"), tick=FALSE, line=1.5)
    legend("topright", c("mock", drugs[-1], "cpmd 24h", "cpmd 36h + mrf 24h", "mrf 36h + cpmd 24h"), bty="n", col=c(col, rep("black", 3)), pch=c(rep(15, length(drugs)), 1, 2, 5), cex=.8)
}

#' PCA plot for pilot perturbation experiment
#' 
#' @param pca List of PCA results with slots x and pvar
#' @annot annot Matrix of 7 columns indicating the metadata for the pilot experiment
plotPCApilotPerturbationsGG <- function(pca, annot, pc = c(1, 2)) {
    tmp <- tibble::tibble(x=pca$x[pc[1], ], y=pca$x[pc[2], ], Perturbagen=annot[, "perturbagen"], Condition="cpmd 24h")
    tmp$Perturbagen[tmp$Perturbagen=="none"] <- "Mock"
    pert <- sort(unique(tmp$Perturbagen))
    pert <- c("Mock", pert[pert != "Mock"])
    tmp$Perturbagen <- factor(tmp$Perturbagen, levels=pert)
    tmp$Condition[annot[, "morphogen"]=="RA" & annot[, "m_time"]=="36"] <- "cpmd 24h + mrf 36h"
    tmp$Condition[annot[, "morphogen"]=="RA" & annot[, "m_time"]=="24"] <- "cpmd 36h + mrf 24h"
    ggplot(tmp, aes(x=x, y=y)) + theme_classic(base_size = 12) +
        xlab(paste0("PC", pc[1], " - Var: ",  round(pca$pvar[pc[1]]*100, 1), "%")) +
        ylab(paste0("PC", pc[2], " - Var: ", round(pca$pvar[pc[2]]*100, 1), "%")) +
        geom_point(size=4, aes(colour=Perturbagen, shape=Condition)) +
        scale_color_manual(values=c("#666666", RColorBrewer::brewer.pal(10, "Paired")))
}


#' Matrix summary for the tume-course differentiation experiment
#' 
#' @param Matrix of 2 columns with treatment and time-point
timeCourseTreatmentMatrix <- function(diff_samp) {
    diff_arms_ref <- c("Ctrl" = "Ctrl", "RA" = "Retinoic acid", "SB" = "SB431542", "B4" = "BMP4", "Meso" = "Mesoderm", "Endo" = "Endoderm")
    tmp <- as.numeric(sub("H", "", diff_samp[, 2]))
    tmp <- paste(diff_samp[, 1], tmp, sep="-x-")[order(match(diff_samp[, 1], names(diff_arms_ref)), tmp)]
    replics <- table(tmp)
    replics <- replics[match(unique(tmp), names(replics))]
    tmp <- cbind(do.call(rbind, strsplit(names(replics), "-x-")), replics)
    tmp[, 1] <- diff_arms_ref[match(tmp[, 1], names(diff_arms_ref))]
    tmp[, 2] <- paste0(tmp[, 2], "h")
    colnames(tmp) <- c("Treatment", "Time", "Replicates")
    tmp
}

#' Plor PCA for differentiation time course experiments
#' 
#' @param pca List of principal components matrix and proportion of variance
plotPCAdiffTimeCourse <- function(pca) {
    diff_arms_ref <- c("Ctrl" = "Ctrl", "RA" = "Retinoic acid", "SB" = "SB431542", "B4" = "BMP4", "Meso" = "Mesoderm", "Endo" = "Endoderm")
    diff_samp <- vapply(strsplit(colnames(pca$x), "-"), function(x) x[1], character(1))
    par(mai=c(.8, .8, .2, .2))
    plot(pca$x[1, ], pca$x[2, ], xlab=paste0("PC-1 Var = ", round(pca$pvar[1]*100, 1), "%"), ylab=paste0("PC-2 Var = ", round(pca$pvar[2]*100, 1), "%"), axes=FALSE, type="n", xlim=c(0, max(pca$x[1, ])*1.15), ylim=range(pca$x[2, ])*1.05)
    axis(1)
    axis(2)
    d3 <- meanGroup(pca$x)
    tt <- c("RA", "B4", "Endo", "Meso", "SB")
    col <- c("#0097CE", "#FF7F40", "#007167", "#B68150", "#A87BC9")
    for (i in 1:length(tt)) {
        tmp <- d3[, grep(paste0(tt[i], "-"), colnames(d3))]
        tmp <- tmp[, order(as.numeric(sub("H", "", sapply(strsplit(colnames(tmp), "-"), function(x) x[2]))))]
        tmp <- cbind(d3[, "Ctrl-0H"], tmp)
        lines(tmp[1, ], tmp[2, ], col=col[i], lwd=2)
    }
    col1 <- c("black", col)[match(diff_samp, c("Ctrl", "RA", "B4", "Endo", "Meso", "SB"))]
    points(pca$x[1, ], pca$x[2, ], pch=20, cex=.8, col=col1)
    points(pca$x[1, ], pca$x[2, ], pch=1, cex=.9, col="black")
    #col <- hsv(c(.05, .2, .3, .5, .6, .8), .7, .8)
    gap <- diff(range(d3[2, ]))*.05
    gap2 <- diff(range(d3[1, ]))*.05
    for (i in 1:length(tt)) {
        tmp <- d3[, grep(paste0(tt[i], "-"), colnames(d3))]
        tmp <- tmp[, order(as.numeric(sub("H", "", sapply(strsplit(colnames(tmp), "-"), function(x) x[2]))))]
        tmp <- cbind(d3[, "Ctrl-0H"], tmp)
        lines(tmp[1, ], tmp[2, ], col=col[i], lwd=2)
#        text(tmp[1, ][-c(1, ncol(tmp))], tmp[2, ][-c(1, ncol(tmp))]+gap, sub("H", "", sapply(strsplit(colnames(tmp)[-c(1, length(tmp))], "-"), function(x) x[2])), col=col[i])
        text(tmp[1, ], tmp[2, ]+gap, sub("H", "", sapply(strsplit(colnames(tmp), "-"), function(x) x[2])), col=col[i])
        text(tmp[1, ncol(tmp)]+gap2, tmp[2, ncol(tmp)]+gap, diff_arms_ref[match(tt[i], names(diff_arms_ref))], col=col[i], adj=0)
    }
}


#' rpkm
#' 
#' This function computes rpkm for human transcripts
#' 
#' @param dset Numeric matrix of rawcounts
#' @param log Logical, whether the data should be log10 transformed
#' @param offset Number to offset zero values when applying log transformation
#' @param genome Filename with transcript-size data
#' @return Numeric matrix of rpkm data
rpkm <- function(dset, log = TRUE, offset = 1, genome) {
    load(genome)
    dset <- dset[rownames(dset) %in% names(geneLength), , drop=FALSE]
    tmp <- 1e9 *
        t(t(dset) /colSums(dset, na.rm = TRUE)) / geneLength[match(rownames(dset), names(geneLength))]
    if (!log)
        return(tmp)
    log10(tmp + offset)
}

#' Find the modes for a multimodal distribution
#' 
#' This function returns the modes of a multimodal distribution
#' 
#' @param x Numeric vector
#' @param adj Number indicating the adjust parameter for bandwidth of the density estimation
#' @param thr Threshold for lambda, distrivutions with lambda below this threshold will be discarded
#' @return list of three elements: mean, sd and lambda
getPeaks3 <- function(x, adj=1.2, thr=1e-2) {
    den <- density(x, adj=adj, na.rm=TRUE, n=512*50)
    sp <- smooth.spline(den)
    sp2 <- predict(sp, den$x, deriv=2)$y
    pos <- which(sp2[-length(sp2)]*sp2[-1]<0)
    x2 <- (den$x[pos]+den$x[pos+1])/2
    sp3 <- predict(sp, x2, deriv=3)$y
    posi <- sapply(which(sp3<0), function(i, x2, posi) {
        if (length(x2)<(i+1)) return(NULL)
        pos <- which(posi>x2[i] & posi<x2[i+1])
        if (length(pos)==0) return((x2[i+1]+x2[i])/2)
        return(posi[pos[1]])
    }, x2=x2, posi=getPeaks(x, adj=adj))
    posi <- unlist(posi[sapply(posi, length)>0], use.names=FALSE)
    posi1 <- sapply(which(sp3>0), function(i, x2) {
        if (length(x2)<(i+1)) return(NULL)
        (x2[i]+x2[i+1])/2
    }, x2=x2)
    posi1 <- unlist(posi1[sapply(posi1, length)>0], use.names=FALSE)
    posi <- sapply(posi, function(posi, x) which(x>posi)[1], x=den$x)
    posi1 <- sapply(posi1, function(posi, x) which(x>posi)[1], x=den$x)
    m <- sort(den$x[posi])
    if (length(posi1)==0) b <- c(min(den$x)-diff(range(den$x)),  max(den$x)+diff(range(den$x)))
    else b <- c(min(den$x)-diff(range(den$x)), sort(den$x[posi1]), max(den$x)+diff(range(den$x)))
    sigma <- sapply(m, function(x, den, b) {
        ymax <- approx(den, xout=x)$y
        x1 <- den$x[den$y<(ymax/2)]
        dd <- x1-x
        dd1 <- b-x
        x2 <- c(x1[dd<0][which.max(dd[dd<0])], x1[dd>0][which.min(dd[dd>0])])
        x3 <- c(b[dd1<0][which.max(dd1[dd1<0])], b[dd1>0][which.min(dd1[dd1>0])])
        dd2 <- x2-x
        dd3 <- x3-x
        opt <- c(dd2[1]-dd3[1], dd3[2]-dd2[2])
        pos <- which(opt>0)
        if (length(pos)>0) tmp <- dd2[pos][which.min(abs(dd2[pos]))]
        else tmp <- dd3[which.min((dd2-dd3)^2)]
        res <- -3.934e-6+.8493*abs(tmp)
        if (length(res)==0) res <- NA
        res
    }, den=den, b=b)
    lambda <- sapply(1:length(m), function(i, m, sigma, den) {
        if (is.na(sigma[i])) return(NA)
        x2 <- seq(m[i]-sigma[i], m[i]+sigma[i], length=100)
        integrateTZ(x2, approx(den, xout=x2)$y)
    }, m=m, sigma=sigma, den=den)
    lambda[is.na(lambda)] <- 0
    lambda <- lambda/sum(lambda)
    pos <- which(lambda>thr)
    list(m=m[pos], sigma=sigma[pos], lambda=lambda[pos])
}

#' Find the modes for a multimodal distribution
#' 
#' This function returns the modes of a multimodal distribution
#' 
#' @param x Numeric vector
#' @param adj Number indicating the adjust parameter for bandwidth of the density estimation
#' @return Numeric vector of modes
getPeaks <- function(x, adj=1.5) {
    den <- density(x, adj=adj, na.rm=TRUE, n=512*50)
    posi <- which(c(FALSE, diff(diff(den$y)>0)<0, FALSE))
    posi <- posi[order(den$y[posi], decreasing=TRUE)]
    sort(den$x[posi])
}

#' Fit a mixture of gaussian curves
#' 
#' This function fits a mixture of gaussians to the distribution of any population
#' 
#' @param x Numeric vectors
#' @param thr Minimum lambda (proportion of the distribution) to include in the analysis
#' @param min Optional minimum number of distributions
#' @param max Optional maximum number of distributions
#' @param adj Optional vector of 2 components defining the range of values for the density adj parameter
#' @return Object of class mgfit. List of fitted parameters
mixGaussianFit <- function(x, thr=1e-2, min=1, max=1e6, adj=c(.8, 2)) {
    x <- x[is.finite(x)]
    fit <- lapply(seq(min(adj), max(adj), length=10), function(adj, x, thr, den, min, max) {
        param <- getPeaks3(x, adj=adj, thr=thr)
        param <- lapply(param, function(x) c(x, 0))
        if ((length(param$m)-1)<min | (length(param$m)-1)>max) return(NULL)
        while(any(param$lambda<thr)) {
            pos <- which(param$lambda<thr)
            param <- lapply(param, function(x, pos) x[-pos], pos=pos)
            suppressMessages(fit <- normalmixEM(x, lambda=param$lambda, mu=param$m, sd=param$sigma, epsilon=1e-50))
            ye <- sapply(1:length(fit$lambda), function(i, x, fit) {
                dnorm(x, fit$mu[i], fit$sigma[i])*fit$lambda[i]
            }, x=den$x, fit=fit)
            if (is.null(dim(ye))) ye <- matrix(ye, length(ye), 1)
            mse <- mean((rowSums(ye)-den$y)^2)
            param <- list(mu=fit$mu, sigma=fit$sigma, lambda=fit$lambda)
        }
        list(mu=fit$mu, sigma=fit$sigma, lambda=fit$lambda, loglik=fit$loglik, mse=mse)
    }, x=x, thr=thr, den=density(x, n=512*20), min=min, max=max)
    fit <- fit[sapply(fit, length)>0]
    if (length(fit)==0) stop("No sucessful fit with provided parameters", call.=FALSE)
    mse <- sapply(fit, function(x) x$mse)
    mse[!is.finite(mse)] <- 1e6
    res <- fit[[which.min(mse)]]
    class(res) <- "mgfit"
    return(res)
}

#' Predict relative likelihood for mixGaussianFit
#' 
#' This function computes the relative likelihood based on parameters fitted with mixGaussianFit function
#' 
#' @param fit Object generated by mixGaussianFit function
#' @param x Numerical vector of observations
#' @param k Optional vector of integers indicating the distributions to use
#' @return Matrix of relative likelihood, observations in rowa and distributions in columns
predict.mgfit <- function(fit, x, k=NULL) {
    # If k is not specified then fit for all gaussians
    if (is.null(k)) k <- 1:length(fit$mu)
    # Ensure 1 <= k <= gaussians
    k <- k[k %in% (1:length(fit$mu))]
    # Check there are still usable k
    if (length(k)==0) stop("Selected gaussians not present in the fit object", call.=FALSE)
    res <- sapply(k, function(n, fit, x) {
        res1 <- res2 <- NULL
        pos <- x<fit$mu[n]
        if (length(which(pos))>0) {
            tmp1 <- pnorm(x[pos], fit$mu[n], fit$sigma[n], lower.tail=TRUE)
            tmp2 <- 0
            if (n>1) tmp2 <- sapply(x[pos], function(x, fit, n) sum(pnorm(x, fit$mu[1:(n-1)], fit$sigma[1:(n-1)], lower.tail=FALSE)), fit=fit, n=n)
            res1 <- tmp1/(tmp1+tmp2)
        }
        if (length(which(!pos))>0) {
            tmp1 <- pnorm(x[!pos], fit$mu[n], fit$sigma[n], lower.tail=FALSE)
            tmp2 <- 0
            if (n<length(fit$mu)) tmp2 <- sapply(x[!pos], function(x, fit, n) sum(pnorm(x, fit$mu[(n+1):length(fit$mu)], fit$sigma[(n+1):length(fit$mu)], lower.tail=TRUE)), fit=fit, n=n)
            res2 <- tmp1/(tmp1+tmp2)
        }
        tmp <- 1:length(x)
        c(res1, res2)[order(c(tmp[pos], tmp[!pos]))]
    }, fit=fit, x=x)
    if(!is.null(nrow(res))) return(res)
    rownames(res) <- names(x)
    colnames(res) <- 1:ncol(res)
    return(res)
}

#' Plot mixGaissuinFit objects
#' 
#' Plot a mixture of distributions fitted with mixGaussianFit
#' 
#' @param fit Object generated by the function mixGaussianFit
#' @param x Vector of values to plot the distribution
#' @param col Color for the background distribution
#' @param lwd Line weight for the fitted distributions
#' @param fitCol Optional vector of colors for the fitted distributions
#' @param main Character string for the main title
#' @param scaled Whether scaled density lines should be added
#' @param nonscaled Whetehr non-scaled density lines should be added
#' @param ... Additional parameters to pass to the plot function
#' @return Nothing, a plot is generated
plot.mgfit <- function(fit, x, col="grey75", lwd=2, fitCol=NULL, main="", ylim=NULL, scaled=TRUE, nonscaled=FALSE, ...) {
    if (is.null(fitCol)) fitCol <- rainbow(length(fit$mu), .8, .8, end=.8)
    den <- density(x, from=min(x, na.rm=TRUE), to=max(x, na.rm=TRUE), na.rm=TRUE)
    denmax <- 0
    for (i in 1:length(fit$mu)) denmax <- max(denmax, dnorm(den$x, fit$mu[i], fit$sigma[i])*fit$lambda[i])
    denmax <- max(denmax, den$y)
    if (is.null(ylim))
        ylim <- c(0, denmax)
    plot(den, type="n", axes=FALSE, main=main, ylim=ylim, ...)
    axis(1)
    axis(2)
    polygon(c(min(den$x), den$x, max(den$x)), c(0, den$y, 0), col=col, border=NA)
    if (scaled) {
        for (i in 1:length(fit$mu)) {
            lines(den$x, dnorm(den$x, fit$mu[i], fit$sigma[i])*fit$lambda[i], lwd=lwd, col=fitCol[i])
        }
    }
    if (nonscaled) {
        datax <- seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length=100)
        dataymax <- vapply(fit, function(x, datax) {
            max(dnorm(datax, fit$mu[i], fit$sigma[i]))
        }, numeric(1), datax=datax)
        f <- denmax/dataymax/3
        for (i in 1:length(fit$mu)) {
            lines(datax, dnorm(datax, fit$mu[i], fit$sigma[i])*f, lwd=lwd, col=fitCol[i], lty=3)
        }
    }
}

#' Integration with trapezoid method
#' 
#' This function integrate over a numerical range using the trapezoid method
#' 
#' @param x Numeric vector of x values
#' @param y Numeric vector of y values
#' @return Number
integrateTZ <- function(x, y) {
    pos <- order(x)
    x <- x[pos]
    y <- y[pos]
    idx = 2:length(x)
    return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
}

#' Compute VIPER results for the differentiation experiment
#' 
#' @param x Matrix of normalized gene expression
#' @return Viper matrix
viperEpiSCdiff <- function(x, regul) {
    samp_ref <- c("RA", "SB", "B4", "Meso", "Endo")
    treat <- unique(colnames(x)[colnames(x) != "Ctrl-0H"])
    tmp <- do.call(rbind, strsplit(treat, "-"))
    tmp[, 2] <- sub("H", "", tmp[, 2])
    treat <- treat[order(match(tmp[, 1], samp_ref), as.numeric(tmp[, 2]))]
    vp <- mclapply(treat, function(treat1, x, regul, ctrl, treat) {
        test <- x[, colnames(x)==treat1, drop=FALSE]
            ges <- rowTtest(test, ctrl)
        ges <- qnorm(ges$p.value/2, lower.tail=FALSE) * sign(ges$statistic)
        if (ncol(test)<5) {
            pos <- which(treat==treat1)
            if (pos>1 && length(grep(strsplit(treat1, "-")[[1]][1], treat[pos-1]))==1) {
                test1 <- x[, colnames(x)==treat[pos-1], drop=FALSE]
                test <- cbind(test, test1[, sample(ncol(test1), min(ncol(test1), 5-ncol(test))), drop=FALSE])
            }
            if (ncol(test)<5) {
                if (pos==length(treat)) {
                    pos <- pos-2
                } else {
                    pos <- pos+1
                }
                test1 <- x[, colnames(x)==treat[pos], drop=FALSE]
                test <- cbind(test, test1[, sample(ncol(test1), min(ncol(test1), 5-ncol(test))), drop=FALSE])
            }
        }
        dnull <- ttestNull(test, ctrl, per=1000, verbose=FALSE)
        msviper(ges, regul, dnull, minsize=20, adaptive.size=20, verbose=FALSE)$es$nes
    }, x=x, regul=regul, ctrl=x[, colnames(x)=="Ctrl-0H", drop=FALSE], treat=treat, mc.cores=3)
    names(vp) <- treat
    vp <- do.call(cbind, vp)
    return(vp)
}

#' Compute integrated differentiation signatures for eahc time point
#' 
#' @param x Matrix of normalized gene expression
#' @return GES matrix
differentiationSignature <- function(x) {
    samp_ref <- c("RA", "SB", "B4", "Meso", "Endo")
    treat <- unique(colnames(x)[colnames(x) != "Ctrl-0H"])
    tmp <- do.call(rbind, strsplit(treat, "-"))
    tmp[, 2] <- sub("H", "", tmp[, 2])
    treat <- treat[order(match(tmp[, 1], samp_ref), as.numeric(tmp[, 2]))]
    diff_ges <- mclapply(treat, function(treat1, x, ctrl) {
        test <- x[, colnames(x)==treat1, drop=FALSE]
            ges <- rowTtest(test, ctrl)
        ges <- qnorm(ges$p.value[, 1]/2, lower.tail=FALSE) * sign(ges$statistic[, 1])
    }, x=x, ctrl=x[, colnames(x)=="Ctrl-0H", drop=FALSE], mc.cores=3)
    names(diff_ges) <- gsub("H", "", treat)
    diff_ges <- do.call(cbind, diff_ges)
    diff_ges[, colSums(!is.na(diff_ges))>0, drop=FALSE]
}


#' Size Factor
#'
#' This function estimates the size factor for a matrix of raw-counts
#'
#' @param x Matrix of raw-counts
#' @param ref Optional vector of raw-counts to be used as reference
#' @return Vector of size factors
sfactor <- function(x, ref=NULL) {
    if (is.null(nrow(x))) x <- matrix(x, length(x), 1, dimnames=list(names(x), "Sample1")) # convert to matrix if vector
    x1 <- x # copy values to working variable
    x1[is.na(x1)] <- 0 # replace NAs with 0
    x1 <- x1[rowSums(x1==0)==0, , drop=FALSE] # Keep only genes detected (at least 1 count) in all samples
    if (is.null(ref)) ref <- r <- exp(rowMeans(log(x1))) # If not specified, define the reference as the geometric mean by rows
    genes <- intersect(rownames(x1), names(ref)) # common genes between input data and reference
    x1 <- x1[match(genes, rownames(x1)), , drop=FALSE]
    ref <- ref[match(genes, names(ref))] # sort data and reference to be compatible at the gene level
    apply(x1, 2, function(x, ref) 2^median(log2(x)-log2(ref)), ref=ref) # Compute the size factors
}

#' Plot size factor vs sequencing depth and infer outliers from residuals
#'
#' @param x Matrix of raw counts
#' @param outlier number indicating the threshold percentage for identifying outliers
#' @return outlier samples
plotDepthvsSizeFactor <- function(x, outlier=1e6) {
    resid <- sf1 <- depth1 <- samp <- NULL
    sf <- sfactor(x)
    depth <- colSums(x)
    flag <- TRUE
    while(flag && length(sf)>10) {
        flag <- FALSE
        fit <- lm(y~x, data=list(y=sf, x=depth))
        res <- abs(residuals(fit))*100/predict(fit, xout=depth)
        if (max(res)>outlier) {
            flag <- TRUE
            pos <- which.max(res)
            resid <- c(resid, res[pos])
            sf1 <- c(sf1, sf[pos])
            depth1 <- c(depth1, depth[pos])
            samp <- c(samp, names(sf)[pos])
            sf <- sf[-pos]
            depth <- depth[-pos]
        }
    }
    par(mai=c(.8, .8, .2, .2))
    plot(colSums(x), sfactor(x), type="n", xlab="Sequencing Depth", ylab="Size Factor", axes=FALSE)
    axis(1)
    axis(2)
    points(depth, sf, pch=20, cex=.8, col="grey")
    points(depth1, sf1, pch=20, col="red")
    abline(fit)    
    text(depth1, sf1+max(sf)*.02, samp, col="red", adj=0)
    return(samp)
}

#' Enrichment plot for a subset of genes
#' 
#' @param signature Named vector for the signature
#' @param regulon Regulon object
#' @param output String indicating the name of the output file
#' @param genes Vector of genes
#' @param xlab String indicating the xlab
#' @param margin Vector of 4 numbers indicating the margins (bottom, left, top, right)
#' 
#' @return Nothing, plot generated in a pdf output
episcEnrichmentPlot <- function(signature, regulon, output="", genes, xlab="", margin=c(.1, .1, .1, .7*1.5)) {
    genes <- gene2entrez(genes)
    if (length(which(is.na(genes)))>0)
        warning(paste0(paste(names(genes)[is.na(genes)], collapse=", "), " were not found"))
    genes <- genes[!is.na(genes)]
    pos <- which(!(genes %in% names(regulon)))
    if (length(pos)>0) {
        warning(paste0(paste(names(genes)[pos], collapse=", "), " are not represented in the regulons"))
        genes <- genes[-pos]
    }
    regulon <- regulon[match(genes, names(regulon))]
    if (output != "")
        pdf(output, w=1.4*2+sum(margin[c(2, 4)]), h=.22*1.5*length(regulon)+sum(margin[c(1, 3)]), pointsize=15, useD=FALSE)
    par(mai=margin)
    enrichmentPlot(signature, lapply(regulon, function(x) x[[1]]), col=c("darkviolet", "chocolate1"))
    if (margin[4]>.5)
        axis(4, length(regulon):1, names(genes), tick=FALSE, las=2, line=-.5)
    if (margin[1]>.2 & xlab!="")
        axis(1, length(signature)/2, xlab, tick=FALSE, las=1, line=-.5)
    if (output != "")
        dev.off()
}

#' entrez2gene
#' 
#' This function converts human entrezID to gene symbol
#' 
#' @param genes Vector of character strings of geneIDs
#' @return Vector of gene symbols
#' @export
entrez2gene <- function(genes) {
    entrez2gene <- getEntrez2gene()
    tmp <- entrez2gene[match(genes, names(entrez2gene))]
    names(tmp) <- genes
    tmp
}

#' gene2entrez
#' 
#' This function converts human gene symbol to geneIDs
#' 
#' @param genes Vector of human gene symbols
#' @return Vector of gene IDs
#' @export
gene2entrez <- function(genes) {
    entrez2gene <- getEntrez2gene()
    tmp <- names(entrez2gene)[match(genes, entrez2gene)]
    names(tmp) <- genes
    tmp
}

#' get entrez2gene metadata
#' 
#' @return Named vector of gene symbol qith entrezID in names
getEntrez2gene <- function() {
    if (exists("_entrez_2_gene", envir=globalenv())) {
        entrez2gene <- get("_entrez_2_gene", envir=globalenv())
    } else {
        entrez2gene <- readRDS(file.path(data_dir, "entrez2gene.rds"))
        assign("_entrez_2_gene", entrez2gene, envir=globalenv())
    }
    return(entrez2gene)
}

#' get candidates MR info
#' 
#' @return Matrix of entrezID and type (MR or POS)
getCandidatesInfo <- function() {
    if (exists("_candidates_info", envir=globalenv())) {
        cand <- get("_candidates_info", envir=globalenv())
    } else {
        cand <- do.call(rbind, strsplit(readLines(file.path(data_dir, "candidates.tsv")), "\t"))
        assign("_candidates_info", cand, envir=globalenv())
    }
    return(cand)
}

#' enrichmentPlot
#' 
#' Generates bar-code like plots to show enrichment of a set of features on a signature
#' 
#' @param signature Vector of character strings indicating the IDs of the elements in the signature
#' @param regulon List of vectors containing the IDs of the elements of the signature to highlight
#' @param col Character string indicating the color for the barcode plot
#' @param bins Number indicating the number of elements per bin to smooth. Default to 1/20 of the signature length
#' @param offset Number between 0 and 0.5 indicating the separation between sets
#' @param sep Character string indicating the type of separation between sets, either box, line or none
#' @param ... Additional graphic paramiters to pass to the plot function
#' @description This function generates a barcode-like plot based on the signature. This function does not sort the signature, it has to be sorted before passing it to this function.
#' @export
enrichmentPlot <- function(signature, regulon, col="cornflowerblue", bins=length(signature)/20, offset=.1, sep=c("box", "line", "none"), xlab="", ylab="", ylim=NULL, ...) {
    sep <- match.arg(sep)
    if (is.numeric(signature) & !is.null(names(signature))) signature <- names(signature)[order(signature)]
    if (is.null(ylim)) ylim <- c(.5, length(regulon)+.5)
    plot(0, 0, type="n", xlab=xlab, ylab=ylab, xlim=c(0, length(signature)+1), axes=FALSE, ylim=ylim, yaxs="i", xaxs="i", ...)
    regulon <- rev(regulon)
    if (any(sapply(regulon, function(x) {
        if (is.numeric(x)) return(prod(range(x))<0)
        return(FALSE)
    }))) {
        if (length(col)==1) col <- rep(col, 2)
        col <- sapply(col, col2hsv)
        for (i in 1:length(regulon)) {
            reg <- list(names(regulon[[i]])[regulon[[i]]<0], names(regulon[[i]])[regulon[[i]]>=0])
            for (ii in 1:2) {
                densi <- rep(0, length(signature))
                x <- which(signature %in% reg[[ii]])
                if (length(x)>0) {
                    densi[x] <- 1
                    denStep <- round(length(densi)/bins)
                    x1 <- x[x<denStep]
                    x2 <- x[x>=denStep & x <= (length(signature)-denStep)]
                    x3 <- x[x>(length(signature)-denStep)]
                    densiRes <- sapply(x2, function(i, densi, denStep) sum(densi[(i-denStep):(i+denStep)]), densi=densi, denStep=denStep)
                    densiRes <- densiRes/max(densiRes)
                    temp <- hsv(col[1, ii], densiRes, 1-(1-col[3, ii])*densiRes)
                    if (ii==1) for (iii in order(densiRes)) lines(c(x[iii], x[iii]), c(i-.5+offset, i), col=temp[iii])
                    if (ii==2) for (iii in order(densiRes)) lines(c(x[iii], x[iii]), c(i, i+.5-offset), col=temp[iii])
                }
            }
            switch(sep, 
                   box={lines(c(0, length(signature)+1, length(signature)+1, 0, 0), c(i-.5+offset, i-.5+offset, i+.5-offset, i+.5-offset, i-.5+offset))},
                   line={
                       lines(c(0, length(signature)+1), c(i+.5, i+.5))
                       if (i==1) lines(c(0, length(signature)+1), c(.5, .5))
                   })
        }
    }
    else {
        col <- col2hsv(col)
        for (i in 1:length(regulon)) {
            densi <- rep(0, length(signature))
            x <- which(signature %in% regulon[[i]])
            if (length(x)>0) {
                densi[x] <- 1
                denStep <- round(length(densi)/bins)
                x1 <- x[x<denStep]
                x2 <- x[x>=denStep & x <= (length(signature)-denStep)]
                x3 <- x[x>(length(signature)-denStep)]
                densiRes <- sapply(x2, function(i, densi, denStep) sum(densi[(i-denStep):(i+denStep)]), densi=densi, denStep=denStep)
                densiRes <- densiRes/max(densiRes)
                temp <- hsv(col[1], densiRes, 1-(1-col[3])*densiRes)
                for (ii in order(densiRes)) lines(c(x[ii], x[ii]), c(i-.5+offset, i+.5-offset), col=temp[ii])
            }
            switch(sep, 
                   box={lines(c(0, length(signature)+1, length(signature)+1, 0, 0), c(i-.5+offset, i-.5+offset, i+.5-offset, i+.5-offset, i-.5+offset))},
                   line={
                       lines(c(0, length(signature)+1), c(i+.5, i+.5))
                       if (i==1) lines(c(0, length(signature)+1), c(.5, .5))
                   })
        }
    }
}

col2hsv <- function(color) {
    tmp <- col2rgb(color)
    rgb2hsv(tmp[1], tmp[2], tmp[3])
}

#' Plot group viper matrix
#' 
#' @param vp Viper matrix
#' @param genes Vector of genes
plotHeatmapGroups <- function(vp, genes, gama=1.5) {
    genes1 <- gene2entrez(genes)
    annot <- getPositionsFromString(colnames(vp), "-", 1:2)
    pos <- c("RA", "SB", "B4", "Meso", "Endo")
    pos <- order(match(annot[, 1], pos), as.numeric(annot[, 2]))
    annot <- annot[pos, ]
    vp <- vp[, pos]
    vp1 <- lapply(unique(annot[, 1]), function(treat, annot, vp, genes1) {
        vp[match(genes1, rownames(vp)), , drop=FALSE][, annot[, 1]==treat, drop=FALSE]
    }, annot=annot, vp=vp, genes1=genes1)
    names(vp1) <- unique(annot[, 1])
    diff_arms_ref <- c("RA" = "Retinoic acid", "SB" = "SB431542", "B4" = "BMP4", "Meso" = "Mesoderm", "Endo" = "Endoderm")
    layout(matrix(seq_len(length(vp1)*2), 2, length(vp1), byrow=TRUE), widths = c(rep(1.6, length(vp1)-1), 2.2), heights=c(.9+.22*length(genes), 1))
    for (i in 1:length(vp1)) {
        if (i < length(vp1)) {
            par(mai=c(.05, .05, .5, .05))
            plothm(vp1[[i]], gama=gama, scmax=round(max(abs(vp))))
            axis(3, seq_len(ncol(vp1[[i]])), paste0(getPositionsFromString(colnames(vp1[[i]]), "-", 2)[, 1], "h"), tick=FALSE, las=2, line=-.5)
            axis(3, ncol(vp1[[i]])/2+.5, diff_arms_ref[names(vp1)[i]], tick=FALSE, line=1.5, cex.axis=1.2)
        } else {
            par(mai=c(.05, .05, .5, .6))
            plothm(vp1[[i]], gama=gama, scmax=round(max(abs(vp))))
            axis(3, seq_len(ncol(vp1[[i]])), paste0(getPositionsFromString(colnames(vp1[[i]]), "-", 2)[, 1], "h"), tick=FALSE, las=2, line=-.5)
            axis(3, ncol(vp1[[i]])/2+.5, diff_arms_ref[names(vp1)[i]], tick=FALSE, line=1.5, cex.axis=1.2)
            axis(4, nrow(vp1[[i]]):1, entrez2gene(rownames(vp1[[i]])), tick=FALSE, las=2, line=-.5)
        }
    }
    par(mai=c(.5, .2, .05, .2))
    sc <- round(max(vp))
    plothm(matrix(seq(-sc, sc, length=100), 1, 100), grid=FALSE, gama=gama)
    axis(1, seq(1, 100, length=sc+1), seq(-sc, sc, length=sc+1))
    axis(1, 50, "NES", tick=FALSE, line=1.5)
}

#' Plot group of enrichment analysis
#' 
#' @param diff_ges Differential gene expression signature
#' @param regul Interactome
#' @param genes Vector of gene names to include
plotEnrichmentGroups <- function(diff_ges, regul, genes=c("Dnmt3b", "Bcl11a", "Otx2", "Nfyb", "Bcl11b")) {
    pos <- paste0(c("RA", "SB", "B4", "Meso", "Endo"), "-72")
    diff_arms_ref <- c("RA" = "Retinoic acid", "SB" = "SB431542", "B4" = "BMP4", "Meso" = "Mesoderm", "Endo" = "Endoderm")
    layout(matrix(seq_len(length(pos)), 1, length(pos)), widths = c(rep(1.6, length(pos)-1), 2.2))
    for (i in 1:length(pos)) {
        if (i==length(pos)) {
            episcEnrichmentPlot(diff_ges[, pos[i]], regul, genes=genes, margin=c(.4, .1, .1, .7), xlab=diff_arms_ref[i])
        }
        else {
            episcEnrichmentPlot(diff_ges[, pos[i]], regul, genes=genes, margin=c(.4, .1, .1, .1), xlab=diff_arms_ref[i])
        }
    }
}

#' Plot principal component trajectories
#' 
#' @param x Matrix of protein activity with colnames as treatment-time
#' @param n Vector of PCs to plot
#' @param treat String indicating the treatment to plot
#' @param ... Additional parameters to pass to plot function
#' @return Plot of PC trajectories
plotPCAtrajectory <- function(x, n=1:3, treat="", ...) {
    if (treat != "") {
        x <- x[, grep(treat, colnames(x)), drop=FALSE]
    }
    tmp <- prcomp(t(x), center = TRUE, scale. = TRUE)
    pca <- t(tmp$x)
    pvar <- tmp$sdev^2/sum(tmp$sdev^2)
    n <- n[n<nrow(pca)]
    time_points <- as.numeric(sapply(strsplit(colnames(x), "-"), function(x) x[2]))
    pca <- pca[, order(time_points), drop=FALSE]
    time_points <- sort(time_points)
    pca <- pca[n, , drop=FALSE]
    par(mai=c(.8, .8, .5, .2))
    plot(0, 0, type="n", axes=FALSE, xlab="Time points (h)", ylab="Principal components", xlim=c(0, max(time_points)), ylim=range(pca), ...)
    axis(1)
    axis(2)
    col <- rainbow(nrow(pca), s=.8, v=.8, end=.8)
    for (i in 1:nrow(pca)) {
        lines(time_points, pca[i, ], col=col[i], lwd=seq(2, .5, length=nrow(pca))[i])
    }
    legend("bottomright", paste0("PC-", n, ": ", round(pvar[n]*100, 1), "%"), lwd=seq(2, .5, length=nrow(pca)), col=col, bty="n")
}

#' Integrate zscores for the differentiation experiment
#' 
#' @param x VIPER matrix for the differentiation experiment
#' @return Vector of z-scores
integrateEpiSCviperZ <- function(x) {
#    x1 <- stouffers4treatment(x)
#    rowWStouffer(x1, c(1, 1, 1, .5, .5))
    w <- rep(1, ncol(x))
    w[getPositionsFromString(colnames(x), "-", 1) %in% c("Endo", "Meso")] <- .5
    rowWStouffer(x, w)
}

#' Integrate zscores for the differentiation experiment
#' 
#' @param x VIPER matrix for the differentiation experiment
#' @return Vector of z-scores
integrateEpiSCviperP <- function(x) {
#    x1 <- stouffers4treatment(x)
#    x1 <- pnorm(abs(x1), lower.tail=FALSE)*2
#    rowWFisher(x1, c(1, 1, 1, .5, .5))
    x <- pnorm(abs(x), lower.tail=FALSE)*2
    w <- rep(1, ncol(x))
    w[getPositionsFromString(colnames(x), "-", 1) %in% c("Endo", "Meso")] <- .5
    rowWFisher(x, w)
}

#' Stouffers integration for each treatment of the differentiation experiment
#' 
#' @param x VIPER matrix for the differentiation experiment
#' @return Matrix of stoufer integrated samples across time points
stouffers4treatment <- function(x) {
    tmp <- tapply(1:ncol(x), getPositionsFromString(colnames(x), sep="-", pos=1)[, 1], function(i, x) {
        rowSums(x[, i, drop=FALSE], na.rm=TRUE)/sqrt(rowSums(!is.na(x[, i, drop=FALSE])))
    }, x=x)
    do.call(cbind, tmp)
}

#' Get positions from string
#'
#' @param x Vector of strings
#' @param sep String indicaitng the separation characters
#' @param pos Vector of integers indicating the positions to take
#' @return Matrix
getPositionsFromString <- function(x, sep="_", pos=1) {
    tmp <- lapply(strsplit(x, sep), function(x, pos) {
        x[pos]
    }, pos=pos)
    do.call(rbind, tmp)
}

#' Weighted version of Z-score Stouffer's integration
#'
#' This function integrate the z-scores per row using a weighted version of the Stouffer method
#' 
#' @param x Matrix of Z-scores
#' @param w Vector of weights with length = ncol(x) or number indicating the exponent for the weights
#' @param brown Character string indicating the correlation method used for the Brown's correction, either none, pearson, spearman or kendall
#' @param tail Character string indicating which tail to use for the weights, either both, less or greater
#' @return Vector of integrated scores
#' @export

rowWStouffer <- function(x, w=2, brown=c("none", "pearson", "spearman", "kendall"), tail=c("both", "less", "greater")) {
    brown <- match.arg(brown)
    tail=match.arg(tail)
    if (length(w)==1) {
        w <- abs(x)^w
    }
    else {
        if (length(w) != ncol(x)) stop("Length of weight vector should match columns in input matrix x")
        w <- matrix(abs(w), nrow(x), ncol(x), byrow=TRUE)
    }
    if (tail=="less") w[x>0] <- 0
    if (tail=="greater") w[x<0] <- 0
    w[is.na(x)] <- 0
    w <- w/apply(w, 1, max, na.rm=TRUE)
    if (brown=="none") bc <- 0
    else {
        tmp <- abs(cor(x, method=brown, use="pairwise.complete.obs"))
        bc <- mean(tmp[upper.tri(tmp)], na.rm=TRUE)
    }
    xsum <- rowSums(x*w, na.rm=TRUE)
    wsum <- rowSums(w)
    bc * xsum / wsum + (1-bc) * xsum / sqrt(wsum)
}

#' Weighted version of Fisher integration of p-values
#' 
#' @param x Matrix of p-values
#' @param w Vector of weights with length = ncol(x)
#' @return Vector of integrated p-values
rowWFisher <- function(x, w=rep(1, ncol(x))) {
    if (length(w) != ncol(x)) stop("Length of weight vector should match columns in input matrix x")
    w <- matrix(abs(w), nrow(x), ncol(x), byrow=TRUE)
    x <- -log(x)
    w[is.na(x)] <- 0
    w <- w/apply(w, 1, max, na.rm=TRUE)
    x <- x * w
    tmp <- 2 * rowSums(x, na.rm=TRUE)
    pchisq(tmp, 2 * rowSums(w), lower.tail = F)
}    
    
#' Plot a heatmap fo the differentiation signatures
#' 
#' @param x Viper matrix
#' @param gamma Number indicating the exponential transformation for the color scale
#' @return Plot
plotHeatmapDiffViper <- function(x, gamma=1.5) {
    defin <- c("RA"="Retinoic acid", "SB"="SB431542", "B4"="BMP4", "Meso"="Mesoderm", "Endo"="Endoderm")
    layout(matrix(seq_len(2), 2, 1), heights=c(10, 1.5))
    par(mai=c(.05, .2, 1, .2))
    hc <- hclust(dist(x))
    sc <- round(max(abs(x)))
    plothm(x[hc$order, ], grid=FALSE, scmax=sc, gama=gamma)
    annot <- getPositionsFromString(colnames(x), "-", 1:2)
    axis(3, 1:ncol(x), annot[, 2], las=2, tick=FALSE, line=-.5)
    pos1 <- table(annot[, 1])
    pos1 <- cumsum(pos1[match(unique(annot[, 1]), names(pos1))])+.5
    axis(3, c(.5, pos1), rep("", length(pos1)+1), line=2, tck=.02)
    axis(3, pos1-c(pos1[1], diff(pos1))/2, defin[match(names(pos1), names(defin))], tick=FALSE, line=1.5)
    abline(v=pos1[-length(pos1)])
    par(mai=c(.8, 3, .05, 3))
    plothm(matrix(seq(-sc, sc, length=100), 1, 100), grid=FALSE, gama=gamma)
    axis(1, seq(1, 100, length=sc+1), seq(-sc, sc, length=sc+1))
    title(xlab="NES")
}

#' Plot PCA
pcplot <- function (dset, clases = colnames(dset), pc = 1:2, center = T, scale. = F, labels = TRUE, ...) {
    tmp <- prcomp(t(dset), center = center, scale. = scale.)
    d2 <- t(tmp$x)
    pvar <- tmp$sdev^2/sum(tmp$sdev^2)
    col <- c("#000000FF", rainbow(length(unique(clases)) - 1, s = 1, v = 0.8, 0, 0.66))
    col <- col[match(clases, unique(clases))]
    tt <- "p"
    if (labels) 
        tt <- "n"
    plot(d2[pc[1], ], d2[pc[2], ], type = tt, xlim = range(d2[pc[1], ]) + c(-0.3, 0.2) * max(abs(range(d2[pc[1], ]))), ylim = range(d2[pc[2], ]) + c(-0.1, 0.1) * max(abs(range(d2[pc[2], ]))), xlab = paste("PC", "-", pc[1], ", var: ", signif(pvar[pc[1]] * 100, 2), "%", sep = ""), ylab = paste("PC", "-", pc[2], ", var: ", signif(pvar[pc[2]] * 100, 2), "%", sep = ""), col = col, ...)
    if (labels) 
        text(d2[pc[1], ], d2[pc[2], ], colnames(d2), adj = 0.5, col = col, font = 2)
}

#' Plot a heatmap for selected candidates
#' 
#' @param x Viper matrix
#' @param genes Vector of genes for incell
#' @param ps Vector of genes for plateseq
#' @param gamma Number indicating exponential transformaiton for the color scale
#' @return Plot
plotHeatmapCandidateViper <- function(x, genes, ps, gamma=1.5) {
    defin <- c("RA"="Retinoic acid", "SB"="SB431542", "B4"="BMP4", "Meso"="Mesoderm", "Endo"="Endoderm")
    layout(matrix(seq_len(2), 2, 1), heights=c(10, .6))
    par(mai=c(.05, .8, 1, .4))
    x <- x[rownames(x) %in% genes, , drop=FALSE]
    hc <- hclust(dist(x))
    sc <- round(max(abs(x)))
    plothm(x[hc$order, ], scmax=sc, gama=gamma)
    axis(2, nrow(x):1, entrez2gene(rownames(x))[hc$order], tick=FALSE, las=2, line=-.5, cex.axis=.7)
    pos <- which(rev(rownames(x)[hc$order]) %in% ps)
    axis(4, pos, rep("<=", length(pos)), tick=FALSE, las=2, line=-.5, cex.axis=.7)
    annot <- getPositionsFromString(colnames(x), "-", 1:2)
    axis(3, 1:ncol(x), annot[, 2], las=2, tick=FALSE, line=-.5)
    pos1 <- table(annot[, 1])
    pos1 <- cumsum(pos1[match(unique(annot[, 1]), names(pos1))])+.5
    axis(3, c(.5, pos1), rep("", length(pos1)+1), line=2, tck=.02)
    axis(3, pos1-c(pos1[1], diff(pos1))/2, defin[match(names(pos1), names(defin))], tick=FALSE, line=1.5)
    abline(v=pos1[-length(pos1)])
    par(mai=c(.8, 3, .05, 3))
    plothm(matrix(seq(-sc, sc, length=100), 1, 100), grid=FALSE, gama=gamma)
    axis(1, seq(1, 100, length=sc+1), seq(-sc, sc, length=sc+1))
    title(xlab="NES")
}

#' Incell remove outliers
#' 
#' @param incell List containing incell data
#' @param iqrx Threshold as times of IQR
#' @return List od incell data
incellRemoveOutliers <- function(incell, iqrx=6) {
    ctrlnames <- c("Empty-Sh001", "TRC1-Sh002", "TRC2-Sh202", "TRC2-Shc002", "TRC2-Shc202", "Empty-Shc001", "TRC1-Shc002", "Scramble-", "Mock-")
    ctrl <- lapply(incell, function(x, ctrl) {
        unlist(x[which(names(x) %in% ctrl)], use.names=FALSE)
    }, ctrl=ctrlnames)
    thr <- vapply(ctrl, function(x, iqrx) median(x)+diff(quantile(x, c(.5, .75)))*iqrx, numeric(1), iqrx=iqrx)
    for (i in 1:length(incell)) incell[[i]] <- lapply(incell[[i]], function(x, top) x[x<top], top=thr[i])
    incell
}

#' Approximate empirical commulative distribution function
#'
#' This function generates an empirical null model that computes a normalized statistics and p-value
#' 
#' @param dnull Numerical vector representing the null model
#' @param symmetric Logical, whether the distribution should betreated as symmetric around zero and only one tail should be approximated
#' @param k Integer indicating the number of points to evaluate the empirical cummulative probability function
#' @return function with two parameters, \code{x} and \code{alternative}
aecdf <- function(dnull, symmetric=FALSE, k=100) {
    dnull <- dnull[is.finite(dnull)]
    if (symmetric) {
        tmp <- sort(abs(dnull), decreasing=T)
        i <- 4
        n <- 4
        while(n<14) {
            i <- i+1
            n <- length(unique(tmp[1:i]))
            if (n==5) iq1 <- i
        }
        tl1 <- i
        iqr <- quantile(abs(dnull), c(.5, 1-iq1/length(dnull)))
        epd <- ecdf(abs(dnull))
        a <- list(x=abs(dnull), y=epd(abs(dnull)))
        fit <- lm(y~0+x, data=list(x=a$x[length(a$x)-(tl1:iq1)+1]-iqr[2], y=log(1-epd(iqr[2]))-log(1-a$y[length(a$x)-(tl1:iq1)+1])))
        val <- seq(0, iqr[2], length=k)
        pd <- approxfun(val, epd(val), method="linear", yleft=0, rule=2)
        dnull <- function(x, alternative=c("two.sided", "greater", "less")) {
            alternative <- match.arg(alternative)
            x1 <- abs(x)
            p <- exp(log(1-pd(iqr[2]))-predict(fit, list(x=x1-iqr[2])))
            p[!is.finite(p)] <- 1
            p <- p * (x1>iqr[2]) + (1-pd(x1)) * (x1<=iqr[2])
            nes <- qnorm(p/2, lower.tail=F)*sign(x)
            switch(alternative,
                   two.sided={p <- p},
                   greater={p <- p/2; p[x<0] <- 1-p[x<0]},
                   less={p <- p/2; p[x>0] <- 1-p[x>0]}
            )
            names(nes) <- names(p) <- names(x)
            list(nes=nes, p.value=p)
        }
        return(dnull)
    }
    tmp <- sort(dnull, decreasing=FALSE)
    i <- 4
    n <- 4
    while(n<14) {
        i <- i+1
        n <- length(unique(tmp[1:i]))
        if (n==5) iq1 <- i
    }
    tl1 <- i
    tmp <- sort(dnull, decreasing=TRUE)
    i <- 4
    n <- 4
    while(n<14) {
        i <- i+1
        n <- length(unique(tmp[1:i]))
        if (n==5) iq2 <- i
    }
    tl2 <- i
    iqr <- quantile(dnull, c(iq1/length(dnull), .5, 1-iq2/length(dnull)))
    epd <- ecdf(dnull)
    #    a <- list(x=knots(epd), y=epd(knots(epd)))
    a <- list(x=dnull, y=epd(dnull))
    fit1 <- lm(y~0+x, data=list(x=a$x[iq1:tl1]-iqr[1], y=log(epd(iqr[1]))-log(a$y[iq1:tl1])))
    fit2 <- lm(y~0+x, data=list(x=a$x[length(a$x)-(tl2:iq2)+1]-iqr[3], y=log(1-epd(iqr[3]))-log(1-a$y[length(a$x)-(tl2:iq2)+1])))
    val <- seq(iqr[1], iqr[3], length=k)
    pd <- approxfun(val, epd(val), method="linear", rule=2)
    dnull <- function(x, alternative=c("two.sided", "greater", "less")) {
        alternative <- match.arg(alternative)
        p1 <- exp(log(pd(iqr[1]))-predict(fit1, list(x=x-iqr[1])))
        p2 <- exp(log(1-pd(iqr[3]))-predict(fit2, list(x=x-iqr[3])))
        p1[!is.finite(p1)] <- 1
        p2[!is.finite(p2)] <- 1
        p <- p1*(x<iqr[1]) + p2*(x>iqr[3]) + pd(x)*(x>=iqr[1] & x<iqr[2]) + (1-pd(x))*(x>=iqr[2] & x<=iqr[3])
        nes <- qnorm(p, lower.tail=F)*sign(x-iqr[2])
        switch(alternative,
               two.sided={p <- p*2},
               greater={p[x<iqr[2]] <- 1-p[x<iqr[2]]},
               less={p[x>=iqr[2]] <- 1-p[x>=iqr[2]]}
        )
        names(nes) <- names(p) <- names(x)
        list(nes=nes, p.value=p)
    }
    return(dnull)
}

#' Incell intra-plate normalization
#' 
#' @param incell List of incell data
#' @return Normalized results as a list
incellNorm <- function(incell) {
    ctrlnames <- c("Empty-Sh001", "TRC1-Sh002", "TRC2-Sh202", "TRC2-Shc002", "TRC2-Shc202", "Empty-Shc001", "TRC1-Shc002", "Scramble-", "Mock-")
    incell_norm <- mclapply(incell, function(x, ctrlnames) {
        tmp <- aecdf(unlist(x[names(x) %in% ctrlnames], use.names=FALSE))
        lapply(x, function(x, tmp) {
            x <- tmp(x)$nes
            tmp <- max(abs(x[is.finite(x)]))
            x[!is.finite(x)] <- tmp*sign(x[!is.finite(x)])
            x
        }, tmp=tmp)
    }, ctrlnames=ctrlnames, mc.cores=3)
    incell_norm
}

#' Integration based on CDF (AOC)
#' 
#' This function integrates a distribution of scores based on the area over the CDF curve
#' 
#' @param x Numeric vector, matrix or list of vectors or matrixes
#' @param xlim Numeric vector of 2 elements indicating the range where to perform the integration
#' @details This function computes the area over the curve for the vector o columns of the matrix provided as input
cdfInteg <- function(x, xlim=NULL) {
    if (is.null(xlim)) xlim <- range(unlist(x, use.names=FALSE), na.rm=TRUE)
    if (is.list(x)) return(sapply(x, cdfInteg, xlim=xlim))
    if (is.matrix(x)) return(apply(x, 2, cdfInteg, xlim=xlim))
    1 - integrateFunction(ecdf(x), xlim[1], xlim[2], steps=1000)/diff(xlim)
}

#' Numerical integration of functions
#' 
#' Integrates numerically a function over a range using the trapezoid method
#' 
#' @param f Function of 1 variable (first argument)
#' @param xmin Number indicating the min x value
#' @param xmax Number indicating the max x value
#' @param steps Integer indicating the number of steps to evaluate
#' @param ... Additional arguments for \code{f}
#' @return Number
integrateFunction <- function(f, xmin, xmax, steps=100, ...) {
    x <- seq(xmin, xmax, length=steps)
    y <- f(x, ...)
    integrateTZ(x, y)
}

#' Violin plot for Incell results
#' 
#' @param x list of incell results
#' @param ... Additional parameters for title
plotIncellViolin <- function(x, ...) {
    ctrlnames <- c("Empty-Sh001", "TRC1-Sh002", "TRC2-Sh202", "TRC2-Shc002", "TRC2-Shc202", "Empty-Shc001", "TRC1-Shc002", "Scramble-", "Mock-")
    ctrl <- which(names(x) %in% ctrlnames)
    sc <- cdfInteg(x)
    sc <- sc-mean(sc[ctrl])
    pos <- order(sc, decreasing=TRUE)
    pos <- c(ctrl, pos[!(pos %in% ctrl)])
    col <- colorScale(sc[pos])
    if (max(unlist(x, use.names=FALSE))>20) {
        violin(x[pos], col=col, horiz=TRUE)
    } else {
        violin(x[pos], col=col, horiz=TRUE, ylim=c(-6, 6))
    }
    title(...)
}

#' Integrate and compute z-score for incell data
#' 
#' @param x list of incell data
#' @return Vector of z-scores
incellZscore <- function(x) {
    ctrlnames <- c("Empty-Sh001", "TRC1-Sh002", "TRC2-Sh202", "TRC2-Shc002", "TRC2-Shc202", "Empty-Shc001", "TRC1-Shc002", "Scramble-", "Mock-")
    # Merge plates
    x <- unlist(x, recursive=FALSE)
    names(x) <- getPositionsFromString(names(x), "\\.", 2)
    # Integration
    x <- cdfInteg(x, xlim=c(-6, 6))
    # Aggregating replicates
    names(x)[names(x) %in% ctrlnames] <- "ctrl"
    xag <- split(x, names(x))
    # Computing z-score
    zs <- vapply(xag[names(xag) != "ctrl"], function(x, y) {
        if (length(x)<2) {
            tmp <- t.test(x-y)
            return(qnorm(tmp$p.value/2, lower.tail=FALSE)*sign(tmp$statistic))
        }
        tmp <- t.test(x, y)
        return(qnorm(tmp$p.value/2, lower.tail=FALSE)*sign(tmp$statistic))
    }, numeric(1), y=xag[[which(names(xag)=="ctrl")]])
    zs
}

#' Table of incell results
#' 
#' @param incell List of incell results
#' @return data.frame
incellTable <- function(incell) {
    ctrlnames <- c("Empty-Sh001", "TRC1-Sh002", "TRC2-Sh202", "TRC2-Shc002", "TRC2-Shc202", "Empty-Shc001", "TRC1-Shc002", "Scramble-", "Mock-")
    cand <- getCandidatesInfo()
    incell <- unlist(incell, recursive=FALSE)
    names(incell) <- getPositionsFromString(names(incell), "\\.", 2)
    incell_size <- vapply(incell, length, numeric(1))
    incell_size <- tapply(incell_size, names(incell_size), sum)
    tmp <- getPositionsFromString(names(incell_size), "-", 1:2)
    entrezid <- tmp[, 1]
    genes <- entrez2gene(tmp[, 1])
    genes[is.na(genes)] <- tmp[is.na(genes), 1]
    tmp[, 1] <- genes
    pos <- !(names(incell_size) %in% ctrlnames)
    ctrl <- rep("Ctrl", length(incell_size))
    ctrl[pos] <- cand[match(entrezid[pos], cand[, 1]), 2]
    pos <- order(pos, tmp[, 1])
    res <- data.frame(Type=ctrl[pos], Target=tmp[pos, 1], Hairpin=tmp[pos, 2], Cells=incell_size[pos])
    res[is.na(res)] <- ""
    res[order(res[, 1]), ]
}

#' violin
#' 
#' Plot density distributions using a violin plot
#' 
#' @param data List or matrix containing the univariate distributions to plot
#' @param col Vector of fill colors. If NA no violin is ploted
#' @param border Vector of line colors for the violin borders
#' @param trim Logical, whether the density estimation should be trimmed to match the data
#' @param scaled Number between 0 and 1 indicating the scaling level.
#' @param adjust Number indicating the adjust parameter for the density estimation
#' @param axes Logical, whether plot axis should be added
#' @param ylim Optional vector of two elements indicating the limits for the y axis
#' @param add Logical, whether the plot should be added on top of hte current default plot
#' @param lty Number indicating the style for the violin border
#' @param lwd Number indicating the width for the violin border
#' @param hline optional vector indicating the horizontal positions to add lines to the plot
#' @param center optional character string indicating the type of center statistic to plot, either mean, median or mode
#' @param horiz Logical, whether the plot should be horizontal
#' @param pch Integer indicating the style of the data points, 0 to not plot them
#' @param pt.cex Size of the points
#' @param pt.col Vector of colors for the points
#' @param pt.disp Number between 0 and 1 indicating the dispersion factor for the points
#' @param sep Number indicating the minimum allowed distance beteen points
#' @return Nothing, a violin plot is created in the default output device
violin <- function(data, col="grey", border="black", trim=TRUE, scaled=0, adjust=1, axes=TRUE, ylim=NULL, add=FALSE, lty=1, lwd=1, hline=NULL, center=c("none", "mean", "median", "mode"), horiz=FALSE, pch=0, pt.cex=1, pt.col="black", pt.disp=1, sep=.01) {
    if (pt.disp>1) pt.disp <- 1
    center <- match.arg(center)
    if (is(data, "matrix")) {
        den <- apply(data, 2, function(x, adjust, trim, ylim) {
            if (trim) return(density(x, adjust=adjust, from=max(min(x, na.rm=TRUE), ylim[1]), to=min(max(x, na.rm=TRUE), ylim[2]), na.rm=TRUE))
            return(density(x, adjust=adjust, na.rm=TRUE))
        }, adjust=adjust, trim=trim, ylim=ylim)
    }
    else if (is.list(data)) {
        den <- lapply(data, function(x, adjust, trim, ylim) {
            if (length(x)==1) return(NA)
            if (trim) return(density(x, adjust=adjust, from=max(min(x, na.rm=TRUE), ylim[1]), to=min(max(x, na.rm=TRUE), ylim[2]), na.rm=TRUE))
            return(density(x, adjust=adjust, na.rm=TRUE))
        }, adjust=adjust, trim=trim, ylim=ylim)
    }
    else stop("data must be a list or matrix", call.=FALSE)
    maxy <- sapply(den, function(x) {
        if (length(x)==1) return(NA)
        max(x$y, na.rm=TRUE)
    })
    den <- lapply(den, function(x, scale) {
        if (length(x)<3) return(NA)
        x$y <- x$y/max(scale, max(x$y, na.rm=TRUE))
        return(x)
    }, scale=max(maxy, na.rm=TRUE)*(1-scaled))
    if (is.null(ylim)) {
        ylim <- sapply(den, function(x) {
            if (length(x)<3) return(c(0, 0))
            range(x$x, na.rm=TRUE)
        })
        ylim <- c(min(ylim[1, ]), max(ylim[2, ]))
    }
    den <- lapply(den, function(x, top) {
        if (length(x)<3) return(NA)
        x$y <- x$y/top*.45
        return(x)
    }, top=max(sapply(den, function(x) {
        if (length(x)<3) return(NA)
        max(x$y, na.rm=TRUE)
    }), na.rm=TRUE))
    if (is(data, "matrix")) {
        switch(center,
               none={center <- NULL},
               mean={center <- colMeans(data, na.rm=TRUE)},
               median={center <- apply(data, 2, median, na.rm=TRUE)}, 
               mode={center <- apply(data, 2, distMode)})
    }
    else {
        switch(center,
               none={center <- NULL},
               mean={center <- sapply(data, mean, na.rm=TRUE)},
               median={center <- sapply(data, median, na.rm=TRUE)}, 
               mode={center <- sapply(data, distMode)})
    }    
    if (!add) plot.new()
    if (horiz) {
        plot.window(ylim=c(.5, length(den)+.5), xlim=ylim, yaxs="i")
        if (axes) {
            axis(1, line=.5)
            axis(2, 1:length(den), names(den), las=2)
        }
        if (length(col)<length(den)) col <- rep(col[1], length(den))
        if (length(pt.col)<length(den)) pt.col <- rep(pt.col[1], length(den))
        if (length(pt.cex)<length(den)) pt.cex <- rep(pt.cex[1], length(den))
        if (length(pch)<length(den)) pch <- rep(pch[1], length(den))
        if (length(border)<length(den)) border <- rep(border[1], length(den))
        if (!is.null(hline)) abline(h=hline)
        for (i in 1:length(den)) {
            if (!is.na(col[i])) {
                if (length(den[[i]])>2) polygon(c(den[[i]]$x, rev(den[[i]]$x)), c(i+den[[i]]$y, i-rev(den[[i]]$y)), col=col[i], border=border[i], lty=lty, lwd=lwd)
            }
            if (pch[i]>0) {
                y <- data[[i]]
                xl <- rep(0, length(y))
                if (length(unique(y))>1) xl <- approx(den[[i]], xout=y)$y
                x1 <- runif(1, -xl[1], xl[1])*pt.disp
                while(length(x1) < length(xl)) {
                    pos <- length(x1)+1
                    x2 <- runif(100, -xl[pos], xl[pos])*pt.disp
                    pos1 <- which.min((y[1:length(x1)]-y[pos])^2)
                    dd <- abs(x2-x1[pos1])/length(den)
                    pos2 <- which(dd>(sep*pt.disp))
                    if (length(pos2)==0) pos2 <- which.max(dd)
                    x1 <- c(x1, x2[pos2[1]])
                }
                points(y, i+x1, pch=pch[i], col=pt.col[i], cex=pt.cex[i])
            }
        }
        if (length(center)>0) {
            for (i in 1:length(center)) {
                ll <- .5
                if (length(unique(data[[i]]))>2) 
                    ll <- approx(den[[i]]$x, den[[i]]$y, xout=center[i])$y
                lines(rep(center[i], 2), c(i-ll*.9, i+ll*.9))
            }
        }        
    }
    else {
        plot.window(xlim=c(.5, length(den)+.5), ylim=ylim, xaxs="i")
        if (axes) {
            axis(2, line=.5)
            axis(1, 1:length(den), names(den), las=2)
        }
        if (length(col)<length(den)) col <- rep(col[1], length(den))
        if (length(pt.col)<length(den)) pt.col <- rep(pt.col[1], length(den))
        if (length(pt.cex)<length(den)) pt.cex <- rep(pt.cex[1], length(den))
        if (length(pch)<length(den)) pch <- rep(pch[1], length(den))
        if (length(border)<length(den)) border <- rep(border[1], length(den))
        if (!is.null(hline)) abline(h=hline)
        for (i in 1:length(den)) {
            if (!is.na(col[i])) {
                if (length(den[[i]])>2) polygon(c(i+den[[i]]$y, i-rev(den[[i]]$y)), c(den[[i]]$x, rev(den[[i]]$x)), col=col[i], border=border[i], lty=lty, lwd=lwd)
            }
            if (pch[i]>0) {
                y <- data[[i]]
                xl <- rep(0, length(y))
                if (length(unique(y))>1) xl <- approx(den[[i]], xout=y)$y
                x1 <- runif(1, -xl[1], xl[1])*pt.disp
                while(length(x1) < length(xl)) {
                    pos <- length(x1)+1
                    x2 <- runif(100, -xl[pos], xl[pos])*pt.disp
                    pos1 <- which.min((y[1:length(x1)]-y[pos])^2)
                    dd <- abs(x2-x1[pos1])/length(den)
                    pos2 <- which(dd>(sep*pt.disp))
                    if (length(pos2)==0) pos2 <- which.max(dd)
                    x1 <- c(x1, x2[pos2[1]])
                }
                points(i+x1, y, pch=pch[i], col=pt.col[i], cex=pt.cex[i])
            }
        }
        if (length(center)>0) {
            for (i in 1:length(center)) {
                ll <- .5
                if (length(unique(data[[i]]))>2)
                    ll <- approx(den[[i]]$x, den[[i]]$y, xout=center[i])$y
                lines(c(i-ll*.9, i+ll*.9), rep(center[i], 2))
            }
        }
    }
}

#' Plot distributions and proportion of fit
#' 
#' @param x Vector of values to plot the distribution
#' @param min Number indicating the minimum number of distros
#' @param max Number indicating the maximum number of distros
#' @param prop Logical, whether the proportions should be displayed
#' @param ... Additional parameters to plot
#' @return mgfit object
plotGaussianFit <- function(x, min=1, max=10, prop=TRUE, ...) {
    fit <- mixGaussianFit(x, min=min, max=max)
    tmp <- approx(density(x), xout=fit$mu)$y
    plot(fit, x, ylim=c(0, max(tmp)*1.2), ...)
    if (prop) {
        text(fit$mu, tmp+(max(tmp)*.1), paste0(round(fit$lambda*100), "%"), col=rainbow(length(fit$mu), .8, .8, end=.8))
    }
    legend("topright", paste0("log(likelihood): ", round(fit$loglik)), bty="n")
    invisible(fit)
}

#' Table fro plateseq samples
#' 
#' @param x matrix of plateseq experiment metadata
#' @return Matrix
plateseqSamplesTable <- function(x) {
    cand <- getCandidatesInfo()
    x[x[, 1] %in% c("Mock", "NT"), 2] <- "Control"
    tmp <- table(paste(x[, 1], x[, 2], sep="-"))
    tmp1 <- cbind(getPositionsFromString(names(tmp), "-", 1:2), tmp)
    tmp <- entrez2gene(tmp1[, 1])
    tmp[is.na(tmp)] <- names(tmp)[is.na(tmp)]
    res <- cbind(EntrezID=tmp1[, 1], Symbol=tmp, Type=cand[match(tmp1[, 1], cand[, 1]), 2], Clone=tmp1[, 2], Replicates=tmp1[, 3])
    res[is.na(res)] <- ""
    res[order(as.numeric(res[, 4] != "Control"), res[, 3], res[, 2]), ]
}

#' RNAseq data quantile normalization
#' 
#' This function performs quantile normalizartion and missing data imputation for RNAseq data
#' 
#' @param rawcounts Matrix of raw counts
#' @param ref Matrix of raw counts to use as reference
#' @param thr Threshold for gene detection
#' @param cores Number of CPU cores to use
#' @return Normalized counts matrix
rnaseqQNorm <- function(rawcounts, ref=NULL, thr=1, cores=1, seed=1) {
    if (seed>0) {
        set.seed(seed)
        if (cores>1) {
            RNGkind("L'Ecuyer-CMRG")
            mc.reset.stream()
        }
    }
    if (is.null(nrow(rawcounts))) rawcounts <- matrix(rawcounts, length(rawcounts), 1, dimnames=list(names(rawcounts), "Sample1"))
    # Global distribution
    if (!is.null(ref)) {
        if (is.null(nrow(ref))) ref <- matrix(ref, length(ref), 1, dimnames=list(names(refg), "Sample1"))
        grc <- rowSums(ref)
    }
    else grc <- rowSums(rawcounts)
    grc[grc<=thr] <- NA
    grr <- rank(grc, na.last="keep", ties.method="random")
    # Rank transformation and imputation    
    d1 <- mclapply(1:ncol(rawcounts), function(i, rawcounts, thr, grr, grc) {
        tmp <- rank(rawcounts[, i], na.last="keep", ties.method="random")
        tmp[rawcounts[, i] <= thr] <- NA
        pos <- which(is.na(tmp) & !is.na(grr))
        tmp[pos] <- rpois(length(pos), rank(grr[pos]))
        r1 <- rank(tmp, na.last="keep", ties.method="random")
        grc[match(round(r1/max(r1, na.rm=TRUE)*max(grr, na.rm=TRUE)), grr)]
    }, rawcounts=rawcounts, thr=thr, grr=grr, grc, mc.cores=min(cores, ncol(rawcounts)))
    d1 <- do.call(cbind, d1)
    colnames(d1) <- colnames(rawcounts)
    rownames(d1) <- rownames(rawcounts)
    d1[is.na(d1)] <- 0
    d1
}

#' Gene expression signatures for the plateseq experiment
#' 
#' @param expmat Plateseq expression matrix
#' @return Matrix with gene expression signatures
plateseqGES <- function(expmat) {
    samp <- getPositionsFromString(colnames(expmat), "-", 1:3)
    pos <- which(samp[, 1] %in% c("Mock", "NT"))
    samp[pos, 1] <- "ctrl"
    samp[pos, 2] <- "xx"
    treat <- paste(samp[, 1], samp[, 2], sep="-")
    batch <- samp[, 3]
    # Computing GES
    ges <- mclapply(unique(treat[treat != "ctrl-xx"]), function(x, expmat, samp, treat, batch) {
        pos1 <- which(treat==x)
        batch1  <- batch[pos1]
        if (length(pos1)==1) {
            ctrl <- which(batch==batch1 & treat=="ctrl-xx")
            if (length(ctrl)==0) return(NULL)
            if (length(ctrl)==1) {
                tmp <- expmat[, pos1] - expmat[, ctrl]
                warning(paste0("No replicates for ", x), call.=FALSE)
                return(tmp/max(abs(tmp)))
            }
            dtmp <- expmat[, pos1] - expmat[, ctrl]
            design <- rep(1, ncol(dtmp))
            tmp <- eBayes(lmFit(dtmp, design))
            return(qnorm(tmp$p.value[, 1]/2, lower.tail=FALSE)*sign(tmp$t[, 1]))
        }
        pos2 <- lapply(unique(batch1), function(x, batch, treat) which(batch==x & treat=="ctrl-xx"), batch=batch, treat=treat)
        names(pos2) <- unique(batch1)
        pos2 <- pos2[lapply(pos2, length)>0]
        pos3 <- batch1 %in% names(pos2)
        pos1 <- pos1[pos3]
        batch1 <- batch1[pos3]
        treat1 <- factor(c(rep("treat", length(pos1)), rep("ctrl", length(unlist(pos2, use.names=FALSE)))), levels=c("ctrl", "treat"))
        if (length(pos2)==1) {
            dtmp <- expmat[, c(pos1, unlist(pos2, use.names=FALSE))]
            design <- model.matrix(~treat1)
            tmp <- eBayes(lmFit(dtmp, design))
            return(qnorm(tmp$p.value[, 2]/2, lower.tail=FALSE)*sign(tmp$t[, 2]))
        }
        batch1 <- factor(c(batch1, rep(names(pos2), sapply(pos2, length))))
        dtmp <- expmat[, c(pos1, unlist(pos2, use.names=FALSE))]
        contrasts(treat1) <- contr.sum(2)
        contrasts(batch1) <- contr.sum(length(levels(batch1)))
        design <- model.matrix(~treat1+batch1)
        tmp <- eBayes(lmFit(dtmp, design))
        return(qnorm(tmp$p.value[, 2]/2, lower.tail=FALSE)*sign(-tmp$t[, 2]))
    }, expmat=expmat, samp=samp[, 1], treat=treat, batch=batch, mc.cores=3)
    names(ges) <- unique(treat[treat != "ctrl-xx"])
    # Collapsing into matrix
    ges <- ges[sapply(ges, length)>0]
    ges <- lapply(ges, function(x, genes) {
        x <- x[match(genes, names(x))]
        x[is.na(x)] <- 0
        names(x) <- genes
        x
    }, genes=unique(unlist(lapply(ges, names), use.names=FALSE)))
    do.call(cbind, ges)
}

#' Integration of hairpin effect per protein to identify inverted regulons
#' 
#' @param vp Matrix of protein activity per silencing clone
#' @param plot Logical whether a plot should be generated
#' @return Plot and inversion likelihood ratio
invertedRegulonLikelihood_old <- function(vp, plot=TRUE) {
    # Hairpins in cis
    genes <- getPositionsFromString(colnames(vp), "-", 1)[, 1]
    vp_cis <- lapply(unique(genes), function(x, vp, genes) {
        vp[, genes==x, drop=FALSE][rownames(vp)==x, ]
    }, vp=vp, genes=genes)
    names(vp_cis) <- unique(genes)
    # Compute z-score for differential activity across hairpins != 0
    prot_z <- vapply(vp_cis, function(x) {
        if (length(x)==1) return(NA)
        tmp <- t.test(x, rep(0, length(x)), paired=TRUE)
        -qnorm(tmp$p.value/2, lower.tail=FALSE)*sign(tmp$statistic)
    }, numeric(1))
    prot_z <- prot_z[!is.na(prot_z)]
    # Fit 2 gaussians to the z-score of differential activity != 0
    set.seed(1)
    fit <- list(mu = c(-1, 1), sigma = c(1, 1), lambda = c(.5, .5), all.loglik = rep(0, 1001))
    count <- 0
    while (length(fit$all.loglik) > 1000 & count < 3) {
        fit <- normalmixEM(prot_z, mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda, maxit = 1000, verb = FALSE)
        count <- count + 1
    }
    class(fit) <- "mgfit"
    a <- pnorm(prot_z, mean=fit$mu[2], sd=fit$sigma[2], lower.tail=TRUE)
    b <- pnorm(prot_z, mean=fit$mu[1], sd=fit$sigma[1], lower.tail=FALSE)
    rl <- a/b # likelihood ratio
    if (plot) {
        par(mai=c(.8, .8, .2, .2))
        plot(fit, prot_z, xlab="Reduction in protein activity (z-score)")
        legend("topleft", paste("log-likelihood: ", round(fit$loglik), sep=""), cex=.9, bty="n")
        thr <- approx(rl, prot_z, xout=.25)$y
        points(thr, 0, pch=6, lwd=2)
    }
    return(rl)
}

#' Integration of hairpin effect per protein to identify inverted regulons
#' 
#' @param vp Matrix of protein activity per silencing clone
#' @param plot Logical whether a plot should be generated
#' @return Plot and inversion likelihood ratio
invertedRegulonLikelihood <- function(vp, plot=TRUE) {
    # Hairpins in cis
    genes <- getPositionsFromString(colnames(vp), "-", 1)[, 1]
    vp_cis <- lapply(unique(genes), function(x, vp, genes) {
        vp[, genes==x, drop=FALSE][rownames(vp)==x, ]
    }, vp=vp, genes=genes)
    names(vp_cis) <- unique(genes)
    # Compute z-score for differential activity across hairpins != 0
    prot_z <- vapply(vp_cis, function(x) {
        if (length(x)==1) return(NA)
        tmp <- t.test(x, rep(0, length(x)), paired=TRUE)
        -qnorm(tmp$p.value/2, lower.tail=FALSE)*sign(tmp$statistic)
    }, numeric(1))
    prot_z <- prot_z[!is.na(prot_z)]
    # Fit 2 gaussians to the z-score of differential activity != 0
    set.seed(1)
    fit <- mixGaussianFit(prot_z)
    a <- pnorm(prot_z, mean=fit$mu[2], sd=fit$sigma[2], lower.tail=TRUE)
    b <- pnorm(prot_z, mean=fit$mu[1], sd=fit$sigma[1], lower.tail=FALSE)
    rl <- a/b # likelihood ratio
    if (plot) {
        par(mai=c(.8, .8, .2, .2))
        plot(fit, prot_z, xlab="Reduction in protein activity (z-score)")
        legend("topleft", paste("log-likelihood: ", round(fit$loglik), sep=""), cex=.9, bty="n")
        thr <- approx(rl, prot_z, xout=.25)$y
        points(thr, 0, pch=6, lwd=2)
    }
    return(rl)
}


#' Correct direction for viper matrix
#' 
#' @param x Viper matrix
#' @param rl Named vector of relative likelihood
#' @param thr Threshold for inverting the regulon
#' @return Corrected viper matrix
viperMatrixCorrectSign <- function(x, rl, thr=.25) {
    rl <- rl[rl<thr]
    if (is.matrix(x)) {
        rl <- rownames(x) %in% names(rl)
        names(rl) <- rownames(x)
    } else {
        rl <- names(x) %in% names(rl)
        names(rl) <- names(x)
    }
    f <- -2*as.numeric(rl)+1
    x * f
}

#' Sigmoid transformation
#' 
#' @param x Numeric vector
#' @param slope Number indicating the slope
#' @param inflection Number indicating the inflection point
#' @return Numeric vector
sigT <- function (x, slope = 20, inflection = 0.5) {
    1 - 1/(1 + exp(slope * (x - inflection)))
}

#' Get Protein activity for the silenced proteins
#' 
#' @param x viper matrix
#' @rerturn Named vector
proteinActivityCis <- function(x) {
    vp_cis <- tapply(seq_len(ncol(x)), getPositionsFromString(colnames(x), "-", 1)[, 1], function(i, vp) {
        vp[, i, drop=FALSE][rownames(vp)==getPositionsFromString(colnames(vp)[i], "-", 1)[1], , drop=FALSE]
    }, vp=x, simplify=FALSE)
    res <- unlist(vp_cis, use.names=FALSE)
    names(res) <- unlist(lapply(vp_cis, colnames), use.names=FALSE)
    res[match(colnames(x), names(res))]
}

#' Integrate shRNA clones
#' 
#' @param x vector or matrix of clones (columns)
#' @param ws Weighting score
#' @param sig Vector of slope and inflection for a sigmoid transformation of the weights
#' @return Integrated vector or matrix by proteins
integrateshRNAclones <- function(x, ws=NULL, sig=NULL) {
    mat <- !is.matrix(x)
    if (mat) {
        x <- matrix(x, 1, length(x), dimnames=list("readout", names(x)))
    }
    if (is.null(ws)) {
        ws <- rep(1, ncol(x))
        names(ws) <- colnames(x)
    }
    if (is.null(sig)) {
        ws <- abs(ws)
        ws <- ws/max(ws)
    } else {
        ws <- sigT(ws, sig[1], sig[2])
    }
    ws <- ws[match(colnames(x), names(ws))]
    names(ws) <- colnames(x)
    ws <- imputeShRNAweight(ws)
    ws[is.na(ws)] <- 1
    genes <- getPositionsFromString(colnames(x), "-", 1)[, 1]
    integ <- tapply(seq_len(ncol(x)), genes, function(i, x, ws) {
        ws <- ws[i]
        if (sum(ws)==0)
            ws <- rep(1, length(ws))
        rowWStouffer(x[, i, drop=FALSE], ws)
    }, x=x, ws=ws, simplify=FALSE)
    integ <- do.call(cbind, integ)
    if (mat) integ <- integ[1, ]
    integ
}

#' Input weight for shRNA integration
#' 
#' @param ws Vector of weighting scores
#' @return Vector of weighting score
imputeShRNAweight <- function(ws) {
    tmp <- tapply(ws, getPositionsFromString(names(ws), "-", 1)[, 1], function(x) {
        x[is.na(x)] <- mean(x, na.rm=TRUE)
        x
    })
    res <- unlist(tmp, use.names=FALSE)
    names(res) <- unlist(lapply(tmp, names), use.names=FALSE)
    res[is.na(res)] <- 1
    res[match(names(ws), names(res))]
}

#' Get full lineage trajectory from acronym
#' 
#' @param x Vector of acronyms
#' @return Vector of lineages
acro2lineage <- function(x) {
    ref <- c("RA"="Retinoic acid", "B4"="BMP4", "SB"="SB431542", "Endo"="Endoderm", "Meso"="Mesoderm")
    ref[match(x, names(ref))]
}

#' Waterfall plot with IDs
#'
#' @param x Vector of values to plot
#' @param pval Number indicating the p-value threshold
#' @param genes selected genes to highlight
#' @return Plot
waterfallPlot <- function(x, pval=.01, genes) {
    if (length(which(genes %in% names(x)))==0) {
        names(x) <- entrez2gene(names(x))
    }
    x <- sort(x)
    suppressWarnings(thr <- approx(p.adjust(pnorm(abs(x), lower.tail=FALSE)*2, "fdr"), abs(x), xout=pval)$y)
    dat <- data.frame(genes=names(x), x=seq_len(length(x)), y=x, color="grey", label="", lab.col="black")
    genes <- genes[genes %in% dat$genes]
    pos <- which(dat$genes %in% genes)
    dat$label[pos] <- dat$genes[pos]
    dat$color[pos] <- hsv(c(.6, .05)[as.numeric(dat$y[pos]>0)+1], .8, .8)
    dat$color[pos][-thr < dat$y[pos] & dat$y[pos] < thr] <- "black"
    dat$lab.col[dat$y<(-thr)] <- "darkblue"
    dat$lab.col[dat$y>thr] <- "darkred"
    dat <- dat[order(dat$label != ""), ]
    ggplot(dat, aes(x, y, label=label)) + geom_hline(yintercept=0, color="grey") + geom_hline(yintercept=c(-thr, thr), color="grey", linetype="dashed")+ geom_text_repel(max.overlaps=Inf, color=dat$lab.col) + geom_point(color=dat$color) + theme_classic(base_size = 16) + labs(x="", y="Pou5f1 protein levels (z-score)") + theme(axis.line.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
}

#' Plot Pou5f1 score vs Differentiation score
#' 
#' @param incell Vector on incell score
#' @param diff Vector of differentiation score
#' @param genes Optional vector of gene names to highlight
#' @param log Logical, whether axis should be log-transformed
#' @param thr P-value threshold
#' @param method p-value correction method
#' @return plot
plotPou5f1vsDiffScore <- function(incell, diff, genes="", log=TRUE, thr=.05, method="fdr") {
    genescom <- intersect(names(incell), names(diff))
    incell <- incell[match(genescom, names(incell))]
    diff <- diff[match(genescom, names(diff))]
    corr <- cor(incell, diff, method="spearman")
    p <- cortest(corr, length(incell))$p.value
    dat <- data.frame(genes=entrez2gene(names(incell)), x=diff, y=incell, col=hsv(0, 0, 0, .5), label="", lab.col="grey")
    pos1 <- which(p.adjust(pnorm(abs(dat$x), lower.tail=FALSE)*2, method)<thr & dat$x>0)
    dat$col[pos1] <- hsv(.3, .8, .6, .5)
    dat$lab.col[pos1] <- hsv(.3, .8, .6)
    pos2 <- which(dat$genes %in% genes)
    dat$label[pos2] <- dat$genes[pos2]
    pval <- p.adjust(pnorm(abs(dat$x), lower.tail=FALSE)*2, method)
    pval_thr <- suppressWarnings(approx(pval, abs(dat$x), xout=thr)$y)
    if (log) {
        #tmp <- log2(abs(dat$x))
        #tmp[tmp<0] <- 0
        #dat$x <- tmp * sign(dat$x)
        dat$x <- log2(abs(dat$x)+1) * sign(dat$x)
        #tmp <- log2(abs(dat$y))
        #tmp[tmp<0] <- 0
        #dat$y <- tmp * sign(dat$y)
        dat$y <- log2(abs(dat$y)+1) * sign(dat$y)
        pval_thr <- log2(pval_thr+1)
    }
    ylim <- range(dat$y)
    ylim <- ylim + c(-diff(ylim)*.05, diff(ylim)*.14)
    xleft <- (min(dat$x)+pval_thr)/2 - pval_thr
    xright <- (max(dat$x)-pval_thr)/2 + pval_thr
    if (any(genes != "")) {
        p1 <- ggplot(dat, aes(x, y, label=label)) + geom_hline(yintercept=0, color="grey", linetype="dashed") +
            geom_vline(xintercept=0, color="grey", linetype="dashed") + geom_vline(xintercept=c(-1, 1) * pval_thr, color=hsv(.3, .8, .6),
            linetype="dotted") + geom_text_repel(max.overlaps=Inf, color=dat$lab.col) + geom_point(color=dat$col) + theme_classic(base_size = 14)
        if (log) {
            p1 <- p1 + labs(x=expression(Differentiation~score~Log[2]('z-score + 1')), y=expression(Pou5f1~protein~score~Log[2]('z-score + 1')))
        } else {
            p1 <- p1 + labs(x="Differentiation score (z-score)", y="Pou5f1 protein score (z-score)")
        }
        p1 <- p1 + annotate(geom="text", x=xleft, y=ylim[2], label=paste0(nrow(dat)-length(pos1), " MRs not validated"), size=4.5, fontface="bold") +
            annotate(geom="text", x=xleft, y=ylim[2]-diff(ylim)*.1, label=paste0("Spearman's Rho: ", round(corr, 2), "\np-value: ",
            signif(p, 2)), size=3.5) + annotate(geom="text", x=xright, y=ylim[2], label=paste0(length(pos1),
            " MRs validated"), color=hsv(.3, .8, .6), size=4.5, fontface="bold") + ylim(ylim[1], ylim[2])
    } else {
        dat <- dat[, 1:4]
        p1 <- ggplot(dat, aes(x, y, label=genes)) + geom_hline(yintercept=0, color="grey", linetype="dashed") +
            geom_vline(xintercept=0, color="grey", linetype="dashed") + geom_vline(xintercept=c(-1, 1) * pval_thr, color=hsv(.3, .8, .6), linetype="dotted") + geom_point(color=dat$col) + theme_classic(base_size = 14)
        if (log) {
            p1 <- p1 + labs(x="Differentiation score Log2(z-score + 1)", y="Pou5f1 protein score Log2(z-score + 1)")
        } else {
            p1 <- p1 + labs(x="Differentiation score (z-score)", y="Pou5f1 protein score (z-score)")
        }
        p1 <- p1 + annotate(geom="text", x=xleft, y=ylim[2], label=paste0(nrow(dat)-length(pos1), " MRs not validated"), size=4.5, fontface="bold") + annotate(geom="text", x=xleft, y=ylim[2]-diff(ylim)*.1, label=paste0("Spearman's Rho: ", round(corr, 2), "\np-value: ", signif(p, 2)), size=3.5) + annotate(geom="text", x=xright, y=ylim[2], label=paste0(length(pos1), " MRs validated"), color=hsv(.3, .8, .6), size=4.5, fontface="bold") + ylim(ylim[1], ylim[2])
    }
    p1
}

#' Plot Pou5f1 score vs Differentiation score 90 rotated
#' 
#' @param incell Vector on incell score
#' @param diff Vector of differentiation score
#' @param genes Optional vector of gene names to highlight
#' @param log Logical, whether axis should be log-transformed
#' @param thr P-value threshold
#' @param method p-value correction method
#' @return plot
plotPou5f1vsDiffScoreRotated <- function(incell, diff, genes="", log=TRUE, thr=.05, method="fdr") {
    genescom <- intersect(names(incell), names(diff))
    incell <- incell[match(genescom, names(incell))]
    diff <- diff[match(genescom, names(diff))]
    corr <- cor(incell, diff, method="spearman")
    p <- cortest(corr, length(incell))$p.value
    dat <- data.frame(genes=entrez2gene(names(incell)), y=diff, x=incell, col=hsv(0, 0, 0, .5), label="", lab.col="grey")
    pos1 <- which(p.adjust(pnorm(abs(dat$y), lower.tail=FALSE)*2, method)<thr & dat$y>0)
    dat$col[pos1] <- hsv(.3, .8, .6, .5)
    dat$lab.col[pos1] <- hsv(.3, .8, .6)
    pos2 <- which(dat$genes %in% genes)
    dat$label[pos2] <- dat$genes[pos2]
    pval <- p.adjust(pnorm(abs(dat$y), lower.tail=FALSE)*2, method)
    pval_thr <- suppressWarnings(approx(pval, abs(dat$y), xout=thr)$y)
    if (log) {
        dat$y <- log2(abs(dat$y)+1) * sign(dat$y)
        dat$x <- log2(abs(dat$x)+1) * sign(dat$x)
        pval_thr <- log2(pval_thr+1)
    }
    ylim <- range(dat$y)
    ylim <- ylim + c(-diff(ylim)*.12, diff(ylim)*.12)
    xleft <- min(dat$x)/2
    xleft <- 0
    xright <- max(dat$x)
    if (any(genes != "")) {
        p1 <- ggplot(dat, aes(x, y, label=label)) + geom_hline(yintercept=0, color="grey", linetype="dashed") +
            geom_hline(yintercept=c(-1, 1) * pval_thr, color=hsv(.3, .8, .6), linetype="dotted") +
            geom_text_repel(max.overlaps=Inf, color=dat$lab.col) + geom_point(color=dat$col) + theme_classic(base_size = 14)
#        p1 <- p1 + geom_vline(xintercept=0, color="grey", linetype="dashed")
        if (log) {
            p1 <- p1 + labs(y=expression(Differentiation~score~Log[2]('z-score + 1')), x=expression(Pou5f1~protein~score~Log[2]('z-score + 1')))
        } else {
            p1 <- p1 + labs(y="Differentiation score (z-score)", x="Pou5f1 protein score (z-score)")
        }
        p1 <- p1 + annotate(geom="text", x=xleft, y=ylim[1], label=paste0(nrow(dat)-length(pos1), " MRs not validated"), size=4.5, fontface="bold") +
            annotate(geom="text", x=xright, y=ylim[2], label=paste0("Spearman's Rho: ", round(corr, 2), "\np-value: ", signif(p, 2)), size=3.5, hjust=1) +
            annotate(geom="text", x=xleft, y=ylim[2], label=paste0(length(pos1), " MRs validated"), color=hsv(.3, .8, .6), size=4.5, fontface="bold") + ylim(ylim[1], ylim[2])
    } else {
        dat <- dat[, 1:4]
        p1 <- ggplot(dat, aes(x, y, label=genes)) + geom_hline(yintercept=0, color="grey", linetype="dashed") +
            geom_hline(yintercept=c(-1, 1) * pval_thr, color=hsv(.3, .8, .6), linetype="dotted") +
            geom_point(color=dat$col) + theme_classic(base_size = 14)
#        p1 <- p1 +  geom_vline(xintercept=0, color="grey", linetype="dashed")
        if (log) {
            p1 <- p1 + labs(y="Differentiation score Log2(z-score + 1)", x="Pou5f1 protein score Log2(z-score + 1)")
        } else {
            p1 <- p1 + labs(y="Differentiation score (z-score)", x="Pou5f1 protein score (z-score)")
        }
        p1 <- p1 + annotate(geom="text", x=xleft, y=ylim[1], label=paste0(nrow(dat)-length(pos1), " MRs not validated"), size=4.5, fontface="bold") +
            annotate(geom="text", x=xright, y=ylim[2], label=paste0("Spearman's Rho: ", round(corr, 2), "\np-value: ", signif(p, 2)), size=3.5, hjust=1) +
            annotate(geom="text", x=xleft, y=ylim[2], label=paste0(length(pos1), " MRs validated"), color=hsv(.3, .8, .6), size=4.5, fontface="bold") +
            ylim(ylim[1], ylim[2])
    }
    p1
}


#' Test for significance of correlation coefficient
#' 
#' This function test whether the correlation coefficient differs from zero
#' 
#' @param r Correlation between x and y
#' @param n sample size
#' @param alternative character string indicating the alternative hypothesis
#' @return z-score and p-value
cortest <- function(r, n, alternative="two.sided") {
    z <- log((1+r)/(1-r))/2*sqrt(n-3)
    switch(match.arg(alternative, c("two.sided", "less", "greater")),
           two.sided={p <- pnorm(abs(z), lower.tail=F)*2},
           less={p <- pnorm(z, lower.tail=T)},
           greater={p <- pnorm(z, lower.tail=F)})
    return(list(z=z, p.value=p))
}

#' Enrichment plot for the differentiation score
#' 
#' @param plateseq_vp VIPER matrix with regulators in rows and shRNA clones in columns
#' @param reg_diff Regulon with the top 100 MRs for each differentiation path at 72h
#' @param gene String indicating the gene symbol to plot
plotDifferentiationEnrichment <- function(plateseq_vp, reg_diff, gene="Pou5f1") {
    diff_treat <- c("RA-72"="Retinoic acid", "SB-72"="SB431542", "B4-72"="BMP4", "Meso-72"="Mesoderm", "Endo-72"="Endoderm")
    diff_treat <- diff_treat[match(names(reg_diff), names(diff_treat))]
    if (!is.na(gene2entrez(gene))) {
        vpmat <- plateseq_vp[, getPositionsFromString(colnames(plateseq_vp), "-")[, 1]==gene2entrez(gene), drop=FALSE]
        if (length(vpmat)==0) {
            vpmat <- plateseq_vp[, getPositionsFromString(colnames(plateseq_vp), "-")[, 1]==gene, drop=FALSE]
        }
    } else {
        vpmat <- plateseq_vp[, colnames(plateseq_vp)==gene, drop=FALSE]
    }
    layout(matrix(seq_len(ncol(vpmat)), 1, ncol(vpmat)), widths=rep(4, ncol(vpmat))+c(.8, rep(0, ncol(vpmat)-1)))
    nes <- aREA(vpmat, reg_diff)$nes
    for (i in seq_len(ncol(nes))) {
        if (i==1) {
            par(mai=c(.5, 1.1, .2, .5))
        } else {
            par(mai=c(.5, .2, .2, .5))
        }
        enrichmentPlot(vpmat[, i], lapply(reg_diff, function(x) x$tfmode),
                       col=c("royalblue","firebrick2"), bins=500,
                       xlab="")
        axis(1, nrow(vpmat)/2, getPositionsFromString(colnames(vpmat), "-", 2)[i], tick=FALSE, line=-.5)
        if (i==1)
            axis(2, length(diff_treat):1, diff_treat, tick=FALSE, las=2, line=-.5)
        axis(4, nrow(nes):1, round(nes[, i], 2), tick=FALSE, las=2, line=-.5)
    }
}


#' Differentiation score
#' 
#' @param nes_plateseq_diff Matrix of enrichment of the differentiation signatures on the silencing data
#' @param plateseq_vp Plateseq viper matrix
#' @param diff_vp Optional matrix of differentiation signatures
#' @return Vector of differentiation scores
differentiationScore <- function(nes_plateseq_diff, plateseq_vp, diff_vp=NULL) {
    genes <- getPositionsFromString(colnames(nes_plateseq_diff), "-", 1)[, 1]
    shrna_z <- proteinActivityCis(plateseq_vp)
    wt <- matrix(-1, length(genes), nrow(nes_plateseq_diff), dimnames=list(genes, rownames(nes_plateseq_diff)))
    if (!is.null(diff_vp)) {
        wt <- sign(diff_vp[match(genes, rownames(diff_vp)), , drop=FALSE][, match(rownames(nes_plateseq_diff), colnames(diff_vp)), drop=FALSE])
    }
    wt <- t(t(wt)/(as.numeric(rownames(nes_plateseq_diff) %in% c("Endo-72", "Meso-72"))+1))
    wt <- wt * shrna_z
    diff_score <- tapply(1:length(genes), genes, function(i, x, ws) {
        tmp <- x[, i, drop=FALSE]
        ws1 <- ws[, i, drop=FALSE]
        sum(tmp*ws1)/sqrt(sum(ws1^2))
    }, x=nes_plateseq_diff, ws=t(wt))
    diff_score
}

#' Tables of differentiation score
#' 
#' @param nes_plateseq_diff Matrix of enrichment of the differentiation signatures on the silencing data
#' @param plateseq_vp Plateseq viper matrix
#' @return Vector of differentiation scores
differentiationScoreTables <- function(nes_plateseq_diff, plateseq_vp) {
    genes <- entrez2gene(getPositionsFromString(colnames(nes_plateseq_diff), "-", 1)[, 1])
    shrna_z <- proteinActivityCis(plateseq_vp)
    # Resorting data
    pos <- order(genes, shrna_z)
    genes <- genes[pos]
    shrna_z <- shrna_z[pos]
    nes_plateseq_diff <- nes_plateseq_diff[, pos, drop=FALSE]
    # Export hairpin table
    exportshRNAtable(t(nes_plateseq_diff), file="Differentiation-shRNA")
    # Integrate differentiation paths
    wt <- matrix(1, length(genes), nrow(nes_plateseq_diff), dimnames=list(genes, rownames(nes_plateseq_diff)))
    wt <- t(t(wt)/(as.numeric(rownames(nes_plateseq_diff) %in% c("Endo-72", "Meso-72"))+1))
    tmp <- rowSums(t(nes_plateseq_diff) * wt) / sqrt(rowSums(wt^2))
    tmp <- matrix(tmp, length(tmp), 1, dimnames=list(names(tmp), "Differentiation-score"))
    exportshRNAtable(tmp, file="Differentiation-integrated-shRNA")
    # Integrate hairpins
    diff_score <- tapply(1:length(genes), genes, function(i, x, ws) {
        tmp <- t(x[, i, drop=FALSE])
        ws1 <- ws[i]
        colSums(tmp*ws1)/sqrt(sum(ws1^2))
    }, x=nes_plateseq_diff, ws=-shrna_z)
    diff_score <- do.call(rbind, diff_score)
    tmp <- cbind(Gene=rownames(diff_score), round(diff_score, 3))
    file <- gzfile(file.path(working_dir, "Differentiation-integrated-gene.tsv.gz"), open="w")
    write.table(tmp, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    close(file)
    # Integrate all
    diff_score <- tapply(1:length(genes), genes, function(i, x, ws) {
        tmp <- x[, i, drop=FALSE]
        ws1 <- ws[, i, drop=FALSE]
        sum(tmp*ws1)/sqrt(sum(ws1^2))
    }, x=nes_plateseq_diff, ws=t(-wt*shrna_z))
    tmp <- cbind(Gene=rownames(diff_score), "Differentiation-score"=round(diff_score, 3))
    file <- gzfile(file.path(working_dir, "Differentiation-integrated.tsv.gz"), open="w")
    write.table(tmp, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    close(file)
}

#' Table of differentiation score
#' 
#' @param nes_plateseq_diff Matrix of enrichment of the differentiation signatures on the silencing data
#' @param plateseq_vp Plateseq viper matrix
#' @param diff_mrs Viper differentiation matrix
#' @return Vector of differentiation scores
differentiationScoreTable <- function(nes_plateseq_diff, plateseq_vp, diff_mrs) {
    diff_score <- differentiationScore(nes_plateseq_diff, plateseq_vp, diff_mrs)
    genes <- getPositionsFromString(colnames(nes_plateseq_diff), "-", 1)[, 1]
    shrna_z <- proteinActivityCis(plateseq_vp)
    # Resorting data
    pos <- order(entrez2gene(genes), shrna_z)
    genes <- genes[pos]
    shrna_z <- shrna_z[pos]
    nes_plateseq_diff <- nes_plateseq_diff[, pos, drop=FALSE]
    platesq_vp <- plateseq_vp[pos, , drop=FALSE]
    # Integrate all
    ref <- c("RA-72"="Retinoic acid", "B4-72"="BMP4", "SB-72"="SB431542", "Endo-72"="Endoderm", "Meso-72"="Mesoderm")
    nes_tmp <- t(nes_plateseq_diff)
    wt <- sign(diff_mrs[match(genes, rownames(diff_mrs)), , drop=FALSE][, match(rownames(nes_plateseq_diff), colnames(diff_mrs)), drop=FALSE])
    pos <- colnames(wt) %in% c("Meso-72", "Endo-72")
    wt[, pos] <- wt[, pos]/2
    ds <- diff_score[match(genes, names(diff_score))]
    tmp <- data.frame(Gene=entrez2gene(genes), shRNA_clone=getPositionsFromString(rownames(nes_tmp), "-", 2)[, 1],
                      round(nes_tmp, 2), shRNA_clone_score=round(shrna_z, 2), wt,
                      "Differentiation-score"=round(ds, 2), DS=round(log2(abs(ds)+1)*sign(ds), 2))
    # Formating for xlsx
    # Creating the workbook
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName="Differentiation score")
    openxlsx::freezePane(wb, sheet=1, firstActiveRow=3, firstActiveCol=3)
    openxlsx::writeData(wb, sheet=1, tmp, startCol=1, startRow=3, colNames=FALSE, rowNames=FALSE)
    # Titles
    titles <- rbind(c("", "", rep("Differentiation enrichment score", 5), "", rep("MR direction and weight", 5), "", ""),
                    c("Gene", "shRNA clone", ref[match(colnames(nes_tmp), names(ref))], "Silencing efficiency",
                      ref[match(colnames(wt), names(ref))], "Differentiation score", "log2(DS+1)"))
    openxlsx::writeData(wb, sheet=1, as.data.frame(titles), startCol=1, startRow=1, colNames=FALSE, rowNames=FALSE)
    openxlsx::mergeCells(wb, sheet=1, cols=3:7, rows=1)
    openxlsx::mergeCells(wb, sheet=1, cols=9:13, rows=1)
    # Merges
    pos <- split(seq_len(nrow(tmp))+2, factor(tmp$Gene, levels=unique(tmp$Gene)))
    borders <- openxlsx::createStyle(border="bottom", borderStyle="thin")
    for (i in 1:length(pos)) {
        openxlsx::mergeCells(wb, sheet=1, cols=1, rows=pos[[i]])
        openxlsx::mergeCells(wb, sheet=1, cols=14, rows=pos[[i]])
        openxlsx::mergeCells(wb, sheet=1, cols=15, rows=pos[[i]])
        openxlsx::addStyle(wb, sheet=1, borders, rows=pos[[i]][length(pos[[i]])], cols=1:15, stack=TRUE)
    }
    openxlsx::saveWorkbook(wb, file=file.path(working_dir, "Differentiation-score-table.xlsx"), overwrite=TRUE)
}


#' Export shRNA table
#' 
#' @param x Matrix with samples in columns and hairpins in rows
#' @param file Filename
#' @return Nothing
exportshRNAtable <- function(x, file="shRNA") {
    tmp <- getPositionsFromString(rownames(x), "-", 1:2)
    tmp[, 1] <- entrez2gene(tmp[, 1])
    tmp <- cbind(Gene=tmp[, 1], shRNA=tmp[, 2], round(x, 3))
    file <- gzfile(file.path(working_dir, paste0(file, ".tsv.gz")), open="w")
    write.table(tmp, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    close(file)
}


#' Export the Differentiation normalized expression matrix
#' 
#' @param x Normalized expression matrix
#' @param file Filename
#' @param gz Logical whether to gzip the output
#' @return None
exportDiffExpmat <- function(x, file, gz=TRUE) {
    ref <- c("RA"="Retinoic acid", "B4"="BMP4", "SB"="SB431542", "Endo"="Endoderm", "Meso"="Mesoderm")
    annot <- getPositionsFromString(colnames(x), "-", 1:2)
    pos <- order(match(annot[, 1], names(ref)), as.numeric(sub("H", "", annot[, 2])))
    x <- x[, pos]
    x <- cbind(EntrezID=rownames(x), Symbol=entrez2gene(rownames(x)), round(x, 5))
    x[is.na(x)] <- ""
    if (gz) file <- gzfile(paste0(file, ".gz"), open="w")
    write.table(x, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    if (gz) close(file)
}

#' Export gene signature
#' @param x Vector of gene signature
#' @param file Filename
#' @param n Number of genes to include
#' @param gz Logical, whether to gzip the output
#' @return None
exportGESsignature <- function(x, file, n = Inf, gz = TRUE) {
    tmp <- data.frame(EntrezID = names(x), Symbol = entrez2gene(names(x)),
        Score = round(x, 4))
    tmp <- tmp[!is.na(tmp$Symbol), , drop = FALSE]
    if (n == Inf) {
        n <- nrow(tmp)
    }
    tmp <- tmp[order(abs(tmp$Score), decreasing = TRUE)[seq_len(n)], , drop = FALSE]
    tmp <- tmp[order(tmp$Score), , drop = FALSE]
    tmp[is.na(tmp)] <- ""
    if (gz) {
        file <- gzfile(paste0(file, ".gz"), open="w")
    }
    write.table(tmp, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    if (gz) {
        close(file)
    }
}

#' Table of genes in a signature
#' 
#' @param x Vector of signature genes
#' @param n Number indicating the nuber of genes to include
#' @return data.frame
tableGESsignature <- function(x, n = Inf) {
    tmp <- data.frame(EntrezID = names(x), Symbol = entrez2gene(names(x)),
        Score = round(x, 4))
    tmp <- tmp[!is.na(tmp$Symbol), , drop = FALSE]
    if (n == Inf) {
        n <- nrow(tmp)
    }
    tmp <- tmp[order(abs(tmp$Score), decreasing = TRUE)[seq_len(n)], , drop = FALSE]
    tmp <- tmp[order(tmp$Score), , drop = FALSE]
    tmp[is.na(tmp)] <- ""
    return(tmp)
}

#' Export the Differentiation viper matrix
#' 
#' @param x Viper matrix
#' @param file Filename
#' @param gz Logical, whether to gzip the output file
#' @return None
exportDiffVPmat <- function(x, file, gz=TRUE) {
    x <- cbind(EntrezID=rownames(x), Symbol=entrez2gene(rownames(x)), round(x, 5))
    x[is.na(x)] <- ""
    if (gz) file <- gzfile(paste0(file, ".gz"), open="w")
    write.table(x, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    if (gz) close(file)
}

#' Heatmap for MRs specific of lineage
#' 
#' @param mr_mat Viper matrix
#' @param ztest
#' @param genes Integer indicating the number of genes
#' @param labels Vector of labels
#' @param gamma Number
#' @return Plot
complexHeatmap <- function(mr_mat, ztest, genes=5, labels=c(RA="Retinoic Acid", "B4"="BMP4", "SB"="SB431542", "Meso"="Mesoderm", "Endo"="Endoderm", "MEndo"="Meso + Endo"), gamma=1) {
    scmax <- max(abs(mr_mat))
    treat <- getPositionsFromString(colnames(mr_mat), "-", 1:2)
    pos <- order(match(treat[, 1], colnames(ztest)), as.numeric(treat[, 2]))
    treat <- treat[pos, ]
    mr_mat <- mr_mat[, pos]
    w <- table(treat[, 1])
    w <- w/w+1
    w <- w[match(colnames(ztest)[colnames(ztest) %in% names(w)], names(w))]
    w <- .1*w+.04
    w[1] <- w[1] + .28
    w[length(w)] <- w[length(w)] + .4
    h <- rep(.1*genes*2, ncol(ztest))+.04
    h[1] <- h[1] + .28
    h[length(h)] <- h[length(h)] + .38
    layout(matrix(seq_len(length(w)*length(h)), length(h), length(w), byrow=TRUE), widths=w, heights=h)
    for (i in seq_len(length(h))) {
        for (ii in seq_len(length(w))) {
            if (i==1) {
                if (ii==1) {
                    par(mai=c(.02, .3, .3, .02))
                } else if (ii==length(w)) {
                    par(mai=c(.02, .02, .3, .6))
                } else {
                    par(mai=c(.02, .02, .3, .02))
                }
            } else if (i==length(h)) {
                if (ii==1) {
                    par(mai=c(.4, .3, .02, .02))
                } else if (ii==length(w)) {
                    par(mai=c(.4, .02, .02, .6))
                } else {
                    par(mai=c(.4, .02, .02, .02))
                }
            } else {
                if (ii==1) {
                    par(mai=c(.02, .3, .02, .02))
                } else if (ii==length(w)) {
                    par(mai=c(.02, .02, .02, .6))
                } else {
                    par(mai=c(.02, .02, .02, .02))
                }
            }
            tmp <- mr_mat[, treat[, 1]==(colnames(ztest)[ii]), drop=FALSE]
            tmp <- tmp[order(ztest[, i], decreasing=TRUE)[c(1:genes, (nrow(tmp)-genes+1):nrow(tmp))], , drop=FALSE]
            tmp <- matrix(rowMeans(tmp, na.rm=TRUE), nrow(tmp), 1, dimnames=list(rownames(tmp), colnames(ztest)[ii]))
            plothm(tmp, scmax=scmax, gama=gamma)
            abline(h=genes+.5, lwd=3)
            abline(h=genes+.5, lwd=1, col="white")
            if (i==ii)
                box(lwd=2)
            if (i==length(h) & ii==(length(w)-1)) {
                abline(v=.5, h=c(0, nrow(tmp))+.5, lwd=3)
            }
            if (i==length(h) & ii==length(w)) {
                abline(v=ncol(tmp)+.5, h=c(0, nrow(tmp))+.5, lwd=3)
            }
            if (i==1) {
                axis(3, ncol(tmp)/2+.5, labels[colnames(ztest)[ii]], tick=FALSE, line=-.5)
            }
            if (i==length(h)) {
                #axis(1, seq_len(ncol(tmp)), paste0(sapply(strsplit(colnames(tmp), "-"), function(x) x[2]), "h"), tick=FALSE, las=2, line=-.5)
            }
            if (ii==1) {
                axis(2, nrow(tmp)/2+.5, labels[colnames(ztest)[i]], tick=FALSE, line=-.5)
            }
            if (ii==length(w)) {
                axis(4, nrow(tmp):1, entrez2gene(rownames(tmp)), tick=FALSE, las=2, line=-.5)
            }
        }
    }
}

#' Master regulators specific of lineage
#' 
#' @param res VIPER matrix
#' @param genes Number of genes
#' @param gamma Gamma transformation for the colors
#' 
mrSpecificLineage <- function(res, genes=3, gamma=1.5) {
    groups <- factor(vapply(strsplit(colnames(res), "-"), function(x) x[1], character(1)), levels=c("RA", "SB", "B4", "Meso", "Endo"))
    groups <- split(seq_len(ncol(res)), groups)
    groups[["MEndo"]] <- c(groups[["Meso"]], groups[["Endo"]])
    # T-test for each group vs. rest
    t_test <- vapply(groups, function(i, res) {
        tmp <- rowTtest(res[, i, drop=FALSE], res[, -i, drop=FALSE])
        qnorm(tmp$p.value[, 1]/2, lower.tail=FALSE) * sign(tmp$statistic[, 1])
    }, numeric(nrow(res)), res=res)
    complexHeatmap(res, t_test, gamma=gamma, genes=genes)
}

#' Table of MRs - S6
#' 
#' @param rpkm Vector of rpkm with entrezID as names
#' @param exp_lh Vector of expression likelihood with EntrezID as names
#' @param regulon Vector of EntrezIDs for regulators
#' @param inversion_lr Inversion likelihood ratio
#' @param msvier_vp Viper matrix computed with msviper
#' @param viper_vp Viper matrix computed with viper
#' @param incell_sc Incell validation score
#' @param diff_sc Differentiation validation score
MRtable <- function(rpkm=NULL, exp_lh=NULL, regulon=NULL, inversion_lr=NULL, msviper_vp=NULL, viper_vp=NULL, incell_sc=NULL, diff_sc=NULL) {
    genes <- do.call(rbind, strsplit(readLines(file.path(data_dir, "candidates.tsv")), "\t"))
    tmp <- cbind(EntrezID=genes[, 1], Symbol=entrez2gene(genes[, 1]))
    if (!is.null(rpkm)) {
        tmp <- cbind(tmp, RPKM=round(rpkm[match(genes[, 1], names(rpkm))], 3))
    }
    if (!is.null(exp_lh)) {
        tmp <- cbind(tmp, Expression_lh=round(exp_lh[match(genes[, 1], names(exp_lh))], 2))
    }
    tmp <- cbind(tmp, MR=genes[, 2]=="MR")
    if (!is.null(regulon)) {
        tmp <- cbind(tmp, Regulon=genes[, 1] %in% regulon)
    }
    incell <- readLines(file.path(data_dir, "incell.txt"))
    plateseq <- readLines(file.path(data_dir, "plateseq.txt"))
    tmp <- cbind(tmp, Pou5f1_Incell=genes[, 1] %in% incell, PLATEseq=genes[, 1] %in% plateseq)
    if (!(is.null(inversion_lr) | is.null(msviper_vp) | is.null(viper_vp))) {
        ref <- c("RA"="Retinoic Acid", "B4"="BMP4", "SB"="SB431542", "Meso"="Mesoderm", "Endo"="Endoderm")
        inversion_lr <- inversion_lr[match(genes[, 1], names(inversion_lr))]
        inversion_lr[is.na(inversion_lr)] <- 1
        msviper_vp <- msviper_vp[match(genes[, 1], rownames(msviper_vp)), , drop=FALSE]
        viper_vp <- viper_vp[match(genes[, 1], rownames(viper_vp)), , drop=FALSE]
        msviper_vp <- msviper_vp[, match(paste0(names(ref), "-72"), colnames(msviper_vp)), drop=FALSE]
        viper_vp <- viper_vp[, match(paste0(names(ref), "-72"), colnames(viper_vp)), drop=FALSE]
        sf1 <- sign(rowMeans(msviper_vp, na.rm=TRUE))
        sf2 <- sign(rowMeans(viper_vp, na.rm=TRUE))
        sf1[sf1==0] <- 1
        sf2[sf2==0] <- 1
        msviper_vp <- msviper_vp * sf1 * sf2 * (-as.numeric(inversion_lr<.25)*2+1)
        colnames(msviper_vp) <- ref[match(getPositionsFromString(colnames(msviper_vp), "-", 1)[, 1], names(ref))]
        tmp <- cbind(tmp, Inverted=inversion_lr < .25 & (sf1*sf2) > 0)
        ws <- rep(1, ncol(msviper_vp))
        ws[colnames(msviper_vp) %in% c("Endoderm", "Mesoderm")] <- .5
        act_integ <- rowWMeans(msviper_vp, ws)
        act_integ_sign <- sign(act_integ)
        act_integ_sign[act_integ_sign==0] <- 1
        tmp <- cbind(tmp, round(msviper_vp, 2), Mean_activity=round(act_integ, 2))
        if (!(is.null(incell_sc) | is.null(diff_sc))) {
            incell_sc <- incell_sc[match(genes[, 1], names(incell_sc))]
            diff_sc <- diff_sc[match(genes[, 1], names(diff_sc))]
            tmp <- cbind(tmp, Pou5f1=round(incell_sc * act_integ_sign, 2), Diff_effect=round(-diff_sc * act_integ_sign, 2))
            tmp <- cbind(tmp, Pou5f1_score=round(incell_sc, 2), Diff_score=round(diff_sc, 2))
            tmp <- cbind(tmp, FDR.incell=signif(p.adjust(pnorm(incell_sc, lower.tail=FALSE)*2, "fdr"), 3))
            tmp <- cbind(tmp, FDR.diff=signif(p.adjust(pnorm(diff_sc, lower.tail=FALSE)*2, "fdr"), 3))
            tmp <- tmp[order(as.numeric(tmp[, "Expression_lh"])<.5, -diff_sc), , drop=FALSE] 
        }
    }
    tmp[is.na(tmp)] <- ""
    return(tmp)
}

#' Weighted means by columns
#' 
#' This function computes weighted means for the columns of a matrix
#' 
#' @param x Numeric matrix
#' @param w Numeric vector of weights of length equal to matrix \code{x} rows
#' @return Vector of weighted means
colWMeans <- function(x, w) {
    w <- w/sum(w, na.rm=TRUE)
    colSums(x*w, na.rm=TRUE)
}

#' Weighted variance by columns
#' 
#' This function computes weighted variance for the columns of a matrix
#' 
#' @param x Numeric matrix
#' @param w Numeric vector of weights of length equal to matrix \code{x} rows
#' @return Vector of weighted variance
colWVars <- function(x, w) {
    w <- w/sum(w, na.rm=TRUE)
    m <- colSums(x*w, na.rm=TRUE)
    colSums(w * t((t(x)-m)^2), na.rm=TRUE)
}

#' Weighted means by rows
#' 
#' This function computes weighted means for the rows of a matrix
#' 
#' @param x Numeric matrix
#' @param w Numeric vector of weights of length equal to matrix \code{x} columns
#' @return Vector of weighted means
rowWMeans <- function(x, w) colWMeans(t(x), w)

#' Weighted variance by rows
#' 
#' This function computes weighted variance for the rows of a matrix
#' 
#' @param x Numeric matrix
#' @param w Numeric vector of weights of length equal to matrix \code{x} columns
#' @return Vector of weighted variance
rowWVars <- function(x, w) colWVars(t(x), w)

#' Silencing efficiency by integrating genes expression and viper activity
#' 
#' @param ges Cis gene expression signature
#' @param vp Cis viper signature
silencingEfficiency <- function(ges, vp) {
    tapply(1:length(ges), getPositionsFromString(names(ges), "-", 1)[, 1], function(i, x) {
        x <- c(-x[[1]][i], abs(x[[2]][i]))
        sum(x, na.rm=T)/sqrt(length(which(is.finite(x))))
    }, x=lapply(list(ges, vp), function(x) x/var(x)))
}

#' Interactome based on ARACNe and silencing experimental data
#' 
#' @param regul ARACNe interactome
#' @param inv_lr Regulon inversion likelihood ratio
#' @param ps_gen Gene expression signatures in response to gene silencing
#' @param ps_vp_cis Protein activity of the targeted protein corrected by inversion of regulon
#' @param sief Silencing efficiency
silencingInteractome <- function(regul, inv_lr, ps_ges, ps_vp_cis, sief) {
    # Evaluable genes
    genes <- getPositionsFromString(colnames(ps_ges), "-", 1)[, 1]
    # Correct interactome for inverted regulons
    regul <- regul[names(regul) %in% unique(genes)]
    # Regulons to invert
    pos <- which(names(regul) %in% (names(inv_lr)[inv_lr < .25]))
    regul[pos] <- lapply(regul[pos], function(x) {
        x[["tfmode"]] <- -x[["tfmode"]]
        x
    })
    class(regul) <- "regulon"
    # Integration of the plateseq gene expression signature
    ps_ges_integ <- tapply(seq_len(ncol(ps_ges)), genes, function(i, ws, ps_ges) {
        tmp <- filterColMatrix(ps_ges, i)
        rowSums(t(t(tmp)*ws[i]))/sum(ws[i])
    }, ws=sigT(ps_vp_cis, -2, -2), ps_ges=ps_ges, simplify=FALSE)
    ps_ges_integ <- do.call(cbind, ps_ges_integ)
    # Getting the MI
    mi <- readRDS(file.path(data_dir, "episc-mi.rds"))
    # Generating the regulon
    regulon <- lapply(names(regul), function(x, reg1, mi, ges, genes, ws) {
        mi <- rank(mi[[x]])/(length(mi[[x]])+1)
        ges <- rank(abs(ges[, x]))/(nrow(ges)+1)*sign(ges[, x])
        tfmode <- reg1[[x]][["tfmode"]]
        tfmode <- tfmode[match(genes, names(tfmode))]
        mi <- mi[match(genes, names(mi))]
        ges <- ges[match(genes, names(ges))]
        names(tfmode) <- names(mi) <- names(ges) <- genes
        int <- cbind(mi*sign(tfmode), ges*ws[x])
        int <- rowSums(int, na.rm=TRUE)/sqrt(rowSums(matrix(c(1, ws[x]), nrow(int), 2, byrow=TRUE)^2 * !is.na(int))) 
        int.a <- cbind(mi, abs(ges)*ws[x])
        int.a <- rowSums(int.a)/sqrt(rowSums(matrix(c(1, ws[x]), nrow(int.a), 2, byrow=TRUE)^2 * !is.na(int.a)))
        tfmode1 <- tfmode
        tfmode1[is.na(tfmode1)] <- 0
        likelihood <- abs(int) * abs(tfmode1) + int.a * (1-abs(tfmode))
        tfmode <- rowWMeans(cbind(tfmode, ges), c(1, ws[x]))
        likelihood <- likelihood/max(likelihood, na.rm=TRUE)
        likelihood[is.na(likelihood)] <- 0
        pos <- order(likelihood, decreasing=TRUE)[1:200]
        list(tfmode=tfmode[pos], likelihood=likelihood[pos])
    }, reg1=regul, mi=mi, ges=-ps_ges_integ, genes=unique(c(names(mi[[1]]), rownames(ps_ges_integ))), ws=sigT(sief, 3, 1))
    names(regulon) <- names(regul)
    class(regulon) <- "regulon"
    return(regulon)
}

#' Scatter-plot for the comparison between viper matrices generated with different interactomes
#' 
#' @param int1 Viper matrix for first interactome
#' @param int2 Viper matrix for second interactome
scatterPlotVPcomparison <- function(int1, int2, xlab="", ylab="") {
    genes <- intersect(rownames(int1), rownames(int2))
    treat <- intersect(colnames(int1), colnames(int2))
    int1 <- int1[match(genes, rownames(int1)), , drop=FALSE][, match(treat, colnames(int1)), drop=FALSE]
    int2 <- int2[match(genes, rownames(int2)), , drop=FALSE][, match(treat, colnames(int2)), drop=FALSE]
    dat <- data.frame(x=as.vector(int1), y=as.vector(int2))
    ggplot(dat, aes(x, y)) + theme_classic(base_size = 16) + labs(x=xlab, y=ylab) + geom_hline(yintercept=0, color="grey", linetype="dotted") + geom_vline(xintercept=0, color="grey", linetype="dotted") + geom_point(color=hsv(0, 0, 0, .3), size=.2) + geom_density_2d(color="grey") + annotate("text", x=min(dat$x), y=max(dat$y), label=paste0("Pearson's R: ", round(cor(dat$x, dat$y), 3)), hjust=0)
}

#' Barplot for the correlations between regulons
#' 
#' @param int1 Viper matrix for first interactome
#' @param int2 Viper matrix for second interactome
barplotVPcomparison <- function(int1, int2) {
    genes <- intersect(rownames(int1), rownames(int2))
    treat <- intersect(colnames(int1), colnames(int2))
    int1 <- int1[match(genes, rownames(int1)), , drop=FALSE][, match(treat, colnames(int1)), drop=FALSE]
    int2 <- int2[match(genes, rownames(int2)), , drop=FALSE][, match(treat, colnames(int2)), drop=FALSE]
    corr <- diag(cor(t(int1), t(int2)))
    corr <- sort(corr)
    names(corr) <- entrez2gene(names(corr))
    dat <- data.frame(gene=names(corr), cor=corr)
    ggplot(dat, aes(x=reorder(gene, cor), y=cor)) + theme_classic(base_size=16) + geom_bar(stat="identity", width=.7) + labs(x="", y="Pearson's correlation") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

#' Export MR-crosstalk table
#' 
#' @param x VIPER matrix
exportTablePsVpInteg <- function(x, adjust="fdr") {
    rownames(x) <- colnames(x) <- entrez2gene(colnames(x))
    tmp <- cbind(MR=rownames(x), round(x, 3))
    file <- gzfile(file.path(working_dir, "MR_crosstalk.tsv.gz"), open="w")
    write.table(tmp, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    close(file)
    pval <- apply(x, 2, function(x, method) {
        x <- p.adjust(pnorm(abs(x), lower.tail=FALSE)*2, method)
        -log10(x)
    }, method=adjust)
    tmp <- cbind(MR=rownames(pval), round(pval, 2))
    file <- gzfile(file.path(working_dir, "MR_crosstalk_log10FDR.tsv.gz"), open="w")
    write.table(tmp, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    close(file)
}

#' Heatmap for the MR protein activity
#' 
#' @param ps_vp_integ Viper matrix
heatmapMReffect <- function(ps_vp_integ) {
    max_val <- round(max(abs(ps_vp_integ)))
    layout(matrix(1:2, 1, 2), widths=c(1, .03))
    par(mai=c(.1, .8, .8, .05))
    plothm(ps_vp_integ, grid=FALSE, scmax=max_val, gama=1)
    axis(2, nrow(ps_vp_integ):1, entrez2gene(rownames(ps_vp_integ)), las=2, tick=FALSE, line=-.5)
    axis(3, 1:ncol(ps_vp_integ), entrez2gene(colnames(ps_vp_integ)), las=2, tick=FALSE, line=-.5)
    par(mai=c(25, .05, .8, .5))
    plothm(matrix(seq(max_val, -max_val, length=100), 100, 1), grid=FALSE, scmax=max_val, gama=1)
    axis(4, seq(1, 100, length=6), seq(-max_val, max_val, length=6), las=2)
    axis(3, 1, "NES", tick=FALSE, las=1, line=-.5)
}

#' Heatmap for the MR protein activity
#' 
#' @param x Viper matrix
heatmapMReffect2 <- function(x) {
    colnames(x) <- entrez2gene(colnames(x))
    rownames(x) <- entrez2gene(rownames(x))
    max_val <- round(max(abs(x)))
    col_fun <- circlize::colorRamp2(c(-max_val, 0, max_val), c("royalblue", "white", "firebrick2"))
    dd <- as.dist(minimumEuclideanDistance(x))
    hc <- hclust(dd)
    hcser <- seriate(dd, "OLO_complete")
    Heatmap(x, name="NES", col=col_fun, rect_gp = gpar(col = "white", lwd = .5), cluster_rows=hcser[[1]], cluster_columns=hcser[[1]], row_dend_reorder=FALSE, column_dend_reorder=FALSE, heatmap_legend_param = list(title = "NES", at = seq(-30, 30, length=5)))
}

#' Analysis of networks of different size
#' 
#' @param vp matrix with perturbations in columns and protein activity in rows
#' @param len Integer indicating the number of points to scan
#' @param gamma Number indicating the exponential transformation for the edges sequence
#' @param edges Optional maximum number of edges
#' @param adjust Multiple hypothesis correction method
#' @return Diagnostic matrix
networkSizeScan <- function(vp, len=1000, gamma=2, edges=NULL, adjust="fdr") {
    tbl <- data.frame(Regulator=rep(colnames(vp), each=nrow(vp)), Target=rep(rownames(vp), ncol(vp)), z.score=as.vector(vp))
    pval <- apply(vp, 2, function(x) {
        p.adjust(pnorm(abs(x), lower.tail=FALSE)*2, method=adjust)
    })
    tbl$p.value <- as.vector(pval)
    # remove self
    tbl <- tbl[tbl$Regulator != tbl$Target, ]
    # Limiting the number of edges
    if (is.null(edges))
        edges <- nrow(tbl)
    # Sort by z.score
    tbl <- tbl[order(abs(tbl$z.score), decreasing=TRUE), ]
    tmp <- lapply(unique(round(seq(0, 1, length=len)^gamma*(edges-1)))+1, function(i, x) {
        ntw <- igraph::graph.data.frame(x[seq_len(i), ], directed=TRUE)
        ntw1 <- as.undirected(ntw)
        clus <- igraph::cluster_louvain(ntw1)
        c(Edges=length(E(ntw)), Nodes=length(V(ntw)), Communities=length(clus), Modularity=modularity(ntw1, membership(clus)))
    }, x=tbl[, 1:2])
    tmp <- do.call(rbind, tmp)
    cbind(tbl[tmp[, 1], , drop=FALSE], tmp)
}

#' Network optimization figure
#' 
#' @param vp Viper matrix
#' @param edges Vector of maximum number of edges for the networks
#' @param adjust Multiple hypothesis correction method
#' @return Plot of Nodes vs. Edges
networkOptimizationFigure <- function(vp, edges, show=c("nodes", "NOE", "community", "modularity"), adjust="fdr") {
    show <- match.arg(show)
    if (length(edges)>1) {
        par(mfrow=c(1, length(edges)))
        tmp <- lapply(edges, function(x, vp, show) {
            networkOptimizationFigure(vp, x, show, adjust=adjust)
        }, vp=vp, show=show)
        return(NULL)
    }
    par(mai=c(.8, .8, .8, .2))
    res <- networkSizeScan(vp, edges=edges, adjust=adjust)
    switch(show,
        nodes={
            plot(res[, 5], res[, 6], pch=20, cex=.5, col=hsv(0, 0, 0, .5), xlab="Edges", ylab="Nodes", axes=FALSE, ylim=c(0, nrow(vp)))
            abline(h=nrow(vp))
            axis(3, 0, nrow(vp), tick=FALSE, line=-1.5, hadj=1)
        },
        NOE={
            plot(res[, 5], log(res[, 6]/res[, 5]), pch=20, cex=.5, col=hsv(0, 0, 0, .5), xlab="Edges", ylab="log(Nodes / Edges) ratio", axes=FALSE)
            x <- round(edges*.75):edges
            fit <- lm(y~x, list(x=x, y=approx(res[, 5], log(res[, 6]/res[, 5]), xout=x)$y))
            abline(fit, col="blue")
        },
        community={
            plot(res[, 5], res[, 7], pch=20, cex=.5, col=hsv(0, 0, 0, .5), xlab="Edges", ylab="Communities", axes=FALSE)
        },
        modularity={
            plot(res[, 5], res[, 8], pch=20, cex=.5, col=hsv(0, 0, 0, .5), xlab="Edges", ylab="Modularity", axes=FALSE)
        })
    axis(1)
    axis(2)
    pvthr <- 2:10
    if (edges>5000)
        pvthr <- pvthr[pvthr != 9]
    suppressWarnings(pos <- approx(-log10(res[, 4]), res[, 5], xout=pvthr)$y)
    abline(v=pos[!is.na(pos)], col="grey", lty=3)
    for (i in which(!is.na(pos))) {
        axis(3, pos[i], niceExponent(10^(-(pvthr[i]))), las=2, line=-.5)
    }
    # pos <- pos[pos != 0]
    # pval <- approx(res[, 5], res[, 4], xout=pos)$y
    # abline(v=pos, col="grey", lty=3)
    # for (i in seq_len(length(pos))) {
    #     axis(3, pos[i], niceExponent(pval[i]), tick=FALSE, line=-.5)
    # }
    axis(3, max(res[, 5])/2, "Regulation FDR", tick=FALSE, line=1.5)
}

#' Nice Exponential representations of scientific notation
#' 
#' This function generates a plotmath or latex representation of scientific notation
#' 
#' @param x Numeric vector
#' @param drop.1 Logical, whether 1 in 1 x type of representatons should be dropped
#' @param sub10 Either logical, "10", a non-negative integer or a length 2 integer vector, indicating if some expression should be formatted traditionally, when integer, all expression before the integer are simplified. when a 2 elements vector, all between the indicated range are simplified
#' @param digits Number of significant digits
#' @param lab.type Character string indicating how the result should look like, either plotmath or latex
#' @param lab.sep Character separator between mantissa and exponent
#' @return Vector of formated numbers
niceExponent <- function(x, drop.1 = TRUE, sub10 = "10", digits = 2, digits.fuzz, lab.type = c("plotmath", "latex"), lab.sep = c("cdot", "times"))
{
    lab.type <- match.arg(lab.type)
    lab.sep <- match.arg(lab.sep)    
    eT <- floor(log10(abs(x)) + 10^-digits)
    mT <- signif(x / 10^eT, digits)
    ss <- vector("list", length(x))
    if(sub.10 <- !identical(sub10, FALSE)) {
        if(identical(sub10, TRUE))
            sub10 <- c(0,0)
        else if(identical(sub10, "10"))
            sub10 <- 0:1
        sub10 <- as.integer(sub10)
        noE <-
            if(length(sub10) == 1) {
                if(sub10 < 0)
                    stop("'sub10' must not be negative if a single number")
                eT <= sub10
            } else if(length(sub10) == 2) {
                stopifnot(sub10[1] <= sub10[2])
                sub10[1] <= eT & eT <= sub10[2]
            } else stop("invalid 'sub10'")
        mT[noE] <- mT[noE] * 10^eT[noE]
    }
    if (lab.type == "plotmath") {
        for(i in seq(along = x))
            ss[[i]] <-
                if(x[i] == 0) quote(0)
        else if(sub.10 &&  noE[i]    ) substitute( A, list(A = mT[i]))
        else if(drop.1 && mT[i] ==  1) substitute( 10^E, list(E = eT[i]))
        else if(drop.1 && mT[i] == -1) substitute(-10^E, list(E = eT[i]))
        else substitute(A %*% 10^E, list(A = mT[i], E = eT[i]))
        do.call("expression", ss)
    }
    else { 
        mTf <- format(mT)
        eTf <- format(eT)
        for(i in seq(along = x))
            ss[[i]] <-
            if(x[i] == 0) ""
        else if(sub.10 &&  noE[i]    ) mTf[i]
        else if(drop.1 && mT[i] ==  1) sprintf("$10^{%s}$", eTf[i])
        else if(drop.1 && mT[i] == -1) sprintf("$-10^{%s}$",eTf[i])
        else sprintf("$%s \\%s 10^{%s}$", mTf[i], lab.sep,  eTf[i])
        unlist(ss, use.names=FALSE)  ## perhaps unlist(ss) ?
    }
}


#' Generate network from perturbation matrix
#' 
#' @param vp Viper matrix with perturbations in columns and protein activity in rows
#' @param pval P-value threshold
#' @param method Multiple testing correction method
#' @param w Logical, whether the network should be weighted
#' @param directed Logical, whether the network is directed
#' @return igraph netowrk object
perturbationNetwork <- function(vp, pval=1e-5, method="bonferroni", w=FALSE, directed=TRUE) {
    ntw <- apply(vp[, match(rownames(vp), colnames(vp))], 2, function(x, pval, method) {
        x <- p.adjust(pnorm(abs(x), lower.tail=FALSE)*2, method)
        res <- -log(x)
        names(res) <- names(x)
        res[x<pval]
    }, pval=pval, method=method)
    tmp <- cbind(rep(names(ntw), sapply(ntw, length)), unlist(lapply(ntw, names), use.names=FALSE), round(unlist(ntw, use.names=FALSE), 3))
    tmp <- tmp[tmp[, 1] != tmp[, 2], ]
    colnames(tmp) <- c("from", "to", "weight")
    if (!directed) {
        tmp <- tmp[order(as.numeric(tmp[, 3]), decreasing=TRUE), ]
        tmp <- tmp[!duplicated(apply(tmp[, 1:2], 1, function(x) paste(sort(x), collapse="-x-"))), ]
    }
    if (w) tmp <- data.frame(from=tmp[, 1], to=tmp[, 2], weight=as.numeric(tmp[, 3]))
    else tmp <- data.frame(from=tmp[, 1], to=tmp[, 2])
    res <- NULL
    try(res <- igraph::graph.data.frame(as.data.frame(tmp), directed=directed))
    return(res)
}

#' Degree table
#' 
#' @param ntw Network as a directed igraph-class object
#' @erturn Table of 4 columns, with degree, out-degree, in-degree and regularized out-degree
degreeTable <- function(ntw) {
    deg <- sapply(c("all", "out", "in"), function(x, ntw) igraph::degree(ntw, mode=x), ntw=ntw)
    colnames(deg) <- c("Degree", "Out-degree", "In-degree")
    #weight for the prioritization based on degree difference
    deg_diff_weight <- deg[, 2]-deg[, 3]
    deg_diff_weight <- deg_diff_weight-min(deg_diff_weight)
    deg_diff_weight <- deg_diff_weight/max(deg_diff_weight)
    #regularized out-degree
    deg <- cbind(deg, "Reg_out-degree"=deg[, 2] * deg_diff_weight)
    #weight for the prioritization based on degree difference
    deg_diff_weight <- deg[, 3]-deg[, 2]
    deg_diff_weight <- deg_diff_weight-min(deg_diff_weight)
    deg_diff_weight <- deg_diff_weight/max(deg_diff_weight)
    #regularized out-degree
    deg <- cbind(deg, "Reg_in-degree"=deg[, 3] * deg_diff_weight)
    # Closeness
    suppressWarnings(closeness <- sapply(c("all", "out", "in"), function(x, ntw) igraph::closeness(ntw, mode=x), ntw=ntw))
    colnames(closeness) <- c("Closeness", "Out-closeness", "In-Closeness")
    deg <- cbind(deg, closeness)
    #Betweness
    betweenness <- igraph::betweenness(ntw)
    deg <- cbind(deg, Betweenness=betweenness)
    #Eigencentrality
    eigen <- eigen_centrality(ntw, directed=TRUE)$vector
    deg <- cbind(deg, Eigen_centrality=eigen)
    #Hub
    hub <- hub_score(ntw)$vector
    deg <- cbind(deg, Hub_score=hub)
    #authority
    authority <- authority_score(ntw)$vector
    deg <- cbind(deg, Authority=authority)
    deg
}

#' Barplots for degree analysis
#' 
#' @param ntw Network
#' @return list with 4 ggplots
barplotDegree <- function(ntw) {
    dg_tbl <- degreeTable(ntw)
    dg_tbl <- dg_tbl[order(dg_tbl[, "Betweenness"], dg_tbl[, "Reg_out-degree"]), , drop=FALSE]
    rownames(dg_tbl) <- entrez2gene(rownames(dg_tbl))
    p <- lapply(colnames(dg_tbl), function(x, dg_tbl) {
        tmp <- data.frame(gene=factor(rownames(dg_tbl), levels=rownames(dg_tbl)), y=dg_tbl[, x])
        ggplot(tmp, aes(x=gene, y=y)) + theme_classic(base_size=16) + geom_bar(stat="identity", width=.7) + labs(x="", y=x) + coord_flip()
    }, dg_tbl=dg_tbl)
    names(p) <- colnames(dg_tbl)
    return(p)
}

#' Heatmap for the degree analysis
#' 
#' @param ntw Network
heatmapDegree <- function(ntw) {
    dg_tbl <- degreeTable(ntw)
    dg_tbl <- dg_tbl[, c("Degree", "Out-degree", "In-degree", "Reg_out-degree", "Reg_in-degree", "Betweenness")]
    dg_tbl <- t(apply(dg_tbl, 2, function(x) x/max(x)))
    dg_tbl <- dg_tbl[, order(dg_tbl[4, ], decreasing=TRUE), drop=FALSE]
    colnames(dg_tbl) <- entrez2gene(colnames(dg_tbl))
    layout(matrix(1:2, 2, 1), heights=c(1, .6))
    par(mai=c(.05, 1.5, .8, .1))
    plothm(dg_tbl, color=c("darkgreen", "darkgreen"))
    axis(2, nrow(dg_tbl):1, rownames(dg_tbl), tick=FALSE, las=2, line=-.5)
    axis(3, 1:ncol(dg_tbl), colnames(dg_tbl), tick=FALSE, las=2, line=-.5)
    par(mai=c(.8, 1.5, .05, 20))
    plothm(matrix(seq(0, 1, length=100), 1, 100), grid=FALSE, col=c("darkgreen", "darkgreen"))
    axis(1, seq(1, 100, length=6), seq(0, 100, length=6))
    axis(1, 50, "Degree (%)", tick=FALSE, line=1.2)
}

#' Scatter-plot for the correlation between silencing efficiency and degree
#' 
#' @param viper matrix
#' @param network
#' @param silencing efficiency
scatterDegreeSief <- function(ps_vp_integ, ntw, sief) {
    dg_tbl <- degreeTable(ntw)
    dg_tbl <- dg_tbl[, c("Out-degree", "Reg_out-degree", "Betweenness")]
    sief <- sief[names(sief) %in% rownames(dg_tbl)]
    diag(ps_vp_integ) <- NA
    tbl <- cbind(x=sief, "Average protein activity"=colMeans(abs(ps_vp_integ), na.rm=TRUE)[match(names(sief), colnames(ps_vp_integ))], dg_tbl[match(names(sief), rownames(dg_tbl)), , drop=FALSE])
    tmp <- lapply(2:ncol(tbl), function(i, tbl) {
        scatterPlotCorrelation(tbl[, 1], tbl[, i], method="spearman", xlab="Silencing Efficiency", ylab=colnames(tbl)[i])
    }, tbl=tbl)
    tmp
}

#' Scatterplot for correlation
#' 
#' @param x Numeric vector
#' @param y Numeric vector
#' @param method Character string indicating the correlation method (pearson or spearman)
#' @param xlab Character string
#' @param ylab Character string
#' @param match_names Logical, whether names should be matched
#' @return ggplot
scatterPlotCorrelation <- function(x, y, method=c("pearson", "spearman"), xlab="", ylab="", match_names=FALSE) {
    method <- match.arg(method)
    if (match_names) {
        genes <- intersect(names(x), names(y))
        x <- x[match(genes, names(x))]
        y <- y[match(genes, names(y))]
    }
    tmp <- data.frame(x=x, y=y)
    corr <- round(cor(x, y, method=method), 2)
    tit <- paste0("Pearson's R: ", corr)
    if (method=="spearman")
        tit <- paste0("Spearman's Rho: ", corr)
    ggplot(tmp, aes(x=x, y=y)) + theme_classic(base_size=16) + geom_point(alpha=5) + labs(x=xlab, y=ylab, subtitle=tit)
}

#' Community stability analysis performed by FET
#' 
#' @param ntw Network in igraph format
#' @param adjust Method for multiple hypothesis correction
#' @return Vector of community stability per node
communityStabilityFET <- function(ntw, adjust="fdr") {
    # Reference community
    #comm_ref <- igraph::cluster_louvain(as.undirected(ntw))
    #mem_ref <- igraph::membership(comm_ref)
    mem_ref <- episcCommunityAnalysis(ntw)
    # Get edges
    res <- sapply(names(mem_ref), communityFET, member=mem_ref, ntw=ntw)
    return(-log10(p.adjust(res, adjust)))
#    res/max(res, na.rm=TRUE)
}

#' Community stability by FET per gene
#'
#' @param gene Gene name
#' @param member Vector of cluster membership
#' @param ntw Network in igraph format
#' @return Score
communityFET <- function(gene, member, ntw) {
    gr <- igraph::get.edgelist(ntw)
    pos <- which(gr[, 2]==gene)
    if (length(pos)>0) {
        tmp <- gr[pos, , drop=FALSE][, 2:1, drop=FALSE]
        gr <- rbind(gr[-pos, , drop=FALSE], tmp)
    }
    mem <- member[match(gr[, 2], names(member))]
    res <- fisher.test(gr[, 1]==gene, mem==member[gene], alternative="greater")
    res$p.value
}

#' Community stability analysis performed by FET for all communities
#' 
#' @param ntw Network in igraph format
#' @param adjust Multiple hypothesis correction
#' @return Vector of community stability per node
allCommunityStabilityFET <- function(ntw, adjust="fdr") {
    # Reference community
    #comm_ref <- igraph::cluster_louvain(as.undirected(ntw))
    #mem_ref <- igraph::membership(comm_ref)
    mem_ref <- episcCommunityAnalysis(ntw)
    # Get edges
    res <- sapply(names(mem_ref), allCommunityFET, member=mem_ref, ntw=ntw)
    res <- t(res)
    colnames(res) <- seq_len(ncol(res))
    matrix(-log10(p.adjust(res, adjust)), nrow(res), ncol(res), dimnames=list(rownames(res), colnames(res)))
}

#' Community stability by FET per gene for all communities
#'
#' @param gene Gene name
#' @param member Vector of cluster membership
#' @param ntw Network in igraph format
#' @return Score
allCommunityFET <- function(gene, member, ntw) {
    gr <- igraph::get.edgelist(ntw)
    pos <- which(gr[, 2]==gene)
    if (length(pos)>0) {
        tmp <- gr[pos, , drop=FALSE][, 2:1, drop=FALSE]
        gr <- rbind(gr[-pos, , drop=FALSE], tmp)
    }
    mem <- member[match(gr[, 2], names(member))]
    sapply(sort(unique(member)), function(x, gr, gene, mem) {
        res <- fisher.test(gr[, 1]==gene, mem==x, alternative="greater")
        res$p.value
    }, gr=gr, gene=gene, mem=mem)
}

#' Heatmap plot for the FET analysis of communities
#' 
#' @param ntw Network
#' @param mem Optional community membership
heatmapCommStbl <- function (ntw, mem=NULL) {
    if (is.null(mem)) {
        mem <- episcCommunityAnalysis(ntw)
    }
    all_fet <- allCommunityStabilityFET(ntw)
    com_stbl <- communityStabilityFET(ntw)
    mem <- mem[match(names(com_stbl), names(mem))]
    pos <- order(mem, -com_stbl)
    tmp <- all_fet[pos, , drop=FALSE]
    com_stbl <- com_stbl[pos]
    mem <- mem[pos]
    layout(matrix(1:2, 1, 2), widths=c(1.67+.85, .33+.45))
    par(mai=c(.1, .8, 1.2, .05))
    plothm(tmp, col=rep("darkgreen", 2))
    axis(2, nrow(tmp):1, entrez2gene(rownames(tmp)), tick=FALSE, las=2, line=-.5)
    axis(3, 1:ncol(tmp), colnames(tmp), tick=FALSE, line=-.5)
    text(rep(1:ncol(tmp), each=nrow(tmp)), rep(nrow(tmp):1, ncol(tmp)), round(tmp, 2), cex=.7)
    par(mai=c(.1, .05, 1.2, .4))
    plothm(matrix(com_stbl, length(com_stbl), 1), col=rep("brown", 2))
    axis(4, nrow(tmp):1, mem, tick=FALSE, las=2, line=-.5)
    axis(3, 1, "Community\nscore", tick=FALSE, las=2, line=-.5)
    text(rep(1, nrow(tmp)), nrow(tmp):1, round(com_stbl, 2), cex=.7)
}

#' FET analysis of communities
#' @param ntw Network
#' @param mem Optional community membership
commStblData <- function(ntw, mem) {
    if (is.null(mem)) {
        mem <- episcCommunityAnalysis(ntw)
    }
    all_fet <- allCommunityStabilityFET(ntw)
    com_stbl <- communityStabilityFET(ntw)
    mem <- mem[match(names(com_stbl), names(mem))]
    pos <- order(mem, -com_stbl)
    tmp <- all_fet[pos, , drop=FALSE]
    com_stbl <- com_stbl[pos]
    mem <- mem[pos]
    res <- data.frame(MR=entrez2gene(rownames(tmp)), round(tmp, 2), Community.score=round(com_stbl, 2), Community=mem)
    colnames(res) <- c("MR", colnames(tmp), "Community.score", "Community")
    return(res)
}

#' Violin plot for communities
#' 
#' @param x Variable to plot
#' @param group cluster membership
#' @param xlab String
#' @param ylab String
#' @param match logical, whether to match names
#' @param ... Additional arguments for violin()
#' @return ggplot
violinCommunities <- function(x, group, xlab="", ylab="", match=TRUE, ...) {
    if (match) {
        tmp <- intersect(names(x), names(group))
        x <- x[match(tmp, names(x))]
        group <- group[match(tmp, names(group))]
    }
    tmp <- split(x, group)
    violin(tmp, axes=FALSE, center="median", ...)
    if (is.numeric(group)) {
        axis(1, las=1)
    } else {
        axis(1, 1:length(tmp), names(tmp), tick=FALSE, las=2, line=-.5)
    }
    axis(2)
    axis(3, seq_len(length(tmp)), sapply(tmp, length), tick=FALSE, las=1, line=0)
    title(xlab=xlab, ylab=ylab)
}

violinCommunitiesTest <- function(x, group, comm_colors, ...) {
    data <- tibble::tibble(x=x, Community=group)
    comp <- combn(unique(group), 2)
    comp <- lapply(seq_len(ncol(comp)), function(i, comp) {
        as.character(comp[, i])
    }, comp=comp)
    p1 <- ggviolin(data, x="Community", y="x", fill="Community", palette=as.vector(comm_colors),
                   trim=TRUE, add="dotplot", add.params=list(fill="black", dotsize=.5), ...) +
        theme(legend.position="none")
#    p1 <- p1 + stat_compare_means(comparisons=comp)
    tests <- compare_means(x ~ Community, comparisons = comp, p.adjust.method = "fdr",
                           method='wilcox.test', data = data)
    p1 <- p1 + stat_pvalue_manual(tests, label = "p.adj", y.position=max(x)*1.1, step.increase=.15)
    posy <- ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range
    p1 + stat_compare_means(label.y=posy[2]*1.05)
}

#' Enrichment analysis for community members
#' 
#' @param x Vector of values
#' @param groups community groups
#' @param title String

communityEnrichmentAnalysis <- function(x, groups, title) {
    gr <- split(names(groups), groups)
    names(gr) <- paste0("Comm-", names(gr))
    pval <- sREA(x, gr) %>% abs() %>% pnorm(lower.tail=FALSE) %>% p.adjust(method="fdr")
    par(mai=c(.6, 1, .2, .8))
    enrichmentPlot(x, gr, col="black", bins=100)
    axis(2, length(gr):1, names(gr), tick=FALSE, las=2, line=-.5)
    axis(4, length(gr):1, signif(pval, 3), tick=FALSE, las=2, line=-.5)
    axis(1, length(x)/2, title, tick=FALSE, las=1, line=-.5)
}

#' Table of network parameters
#' 
#' @param degree_tbl Table of degree results
#' @param com_membership Community membership vector
#' @param com_scr Community score vector
#' @param sief Silencing efficiency vector
#' @return Matrix
tableNtwParams <- function(degree_tbl, com_membership, com_scr, sief, spk_cat) {
    pos <- !(colnames(degree_tbl) %in% c("Degree", "Out-degree", "In-degree"))
    degree_tbl[, pos] <- signif(degree_tbl[, pos], 3)
    sief <- sief+(min(sief))
    sief <- round(100*sief/max(sief), 1)
    com_scr <- round(com_scr, 2)
    tmp <- cbind(EntrezID=rownames(degree_tbl), Symbol=entrez2gene(rownames(degree_tbl)), "Silencing efficiency"=sief[match(rownames(degree_tbl), names(sief))], degree_tbl, "Spk category"=spk_cat[match(rownames(degree_tbl), names(spk_cat))], Community=com_membership[match(rownames(degree_tbl), names(com_membership))], "Comm score"=com_scr[match(rownames(degree_tbl), names(com_scr))])
    tmp[order(as.numeric(tmp[, "Community"]), -as.numeric(tmp[, "Comm score"])), ]
}

#' Regularized out-degree distribution
#' 
#' @param x Degree table
#' @return category vector
fitRegOutDegree <- function(x) {
    x <- x[, "Reg_out-degree"]
    set.seed(1)
    fit <- mixGaussianFit(x, min=4, max=5)
    fit <- lapply(fit, function(x, pos) {
        if (length(x)!=length(pos)) return(x)
        x[pos]
    }, pos=order(fit$mu))
    class(fit) <- "mgfit"
    fit$mu <- fit$mu[1:3]
    fit$sigma <- fit$sigma[1:3]
    fit$lambda <- fit$lambda[1:3]
    par(mai=c(.8, .8, .2, .1))
    plot(fit, x, xlab="Regularized out-degree", nonscaled=TRUE)
    #legend("topright", c("Listeners", "Comunicators-1", "Communicators-2", "Speakers"), bty="n", lwd=2, col=rainbow(length(fit$mu), .8, .8, end=.8))
    lr <- log(pnorm(x, fit$mu[1], fit$sigma[1], lower.tail=FALSE)/pnorm(x, fit$mu[2], fit$sigma[2]))
    suppressWarnings(thr <- c(approx(lr, x, xout=0)$y, qnorm(.05, fit$mu[3], fit$sigma[3], lower.tail=FALSE)))
    abline(v=thr, col="grey", lty=3)
    text(c(0, thr)+diff(c(0, thr, max(x)))/2, mean(axTicks(2)), c("Listeners", "Communicators", "Speakers"), font=2)
    tmp <- predict(fit, x)
    tmp <- apply(tmp, 1, which.max)
    tmp[tmp==3] <- 2
    tmp[pnorm(x, fit$mu[3], fit$sigma[3], lower.tail=FALSE)<.05] <- 3
    names(tmp) <- names(x)
    tmp
}

#' Regularized out-degree distribution
#' 
#' @param x Degree table
#' @return category vector
fitRegOutDegreePlot <- function(x) {
    x <- x[, "Reg_out-degree"]
    set.seed(1)
    fit <- mixGaussianFit(x, min=4, max=5)
    fit <- lapply(fit, function(x, pos) {
        if (length(x)!=length(pos)) return(x)
        x[pos]
    }, pos=order(fit$mu))
    class(fit) <- "mgfit"
    fit$mu <- fit$mu[1:3]
    fit$sigma <- fit$sigma[1:3]
    fit$lambda <- fit$lambda[1:3]
    lr <- log(pnorm(x, fit$mu[1], fit$sigma[1], lower.tail=FALSE)/pnorm(x, fit$mu[2], fit$sigma[2]))
    suppressWarnings(thr <- c(approx(lr, x, xout=0)$y, qnorm(.05, fit$mu[3], fit$sigma[3], lower.tail=FALSE)))
    
    data <- data.frame(x=x, names=entrez2gene(names(x)), y=0)
    pos <- which(data$x<thr[2])
    data$names[pos] <- ""
    pos <- which(data$x>=thr[2])
    p1 <- ggplot(data=data, aes(x=x)) + theme_classic(base_size = 14) +
        geom_density(fill="grey80", color="grey80") +
        geom_vline(xintercept=thr, linetype="dashed", color="grey50") +
        annotate(geom="point", x=data$x[pos], y=rep(0, length(pos))) +
        geom_text_repel(data=data, aes(x=x, y=y, label=names)) +
        xlab("Regularized out-degree") +
        ylab("Density")
    posy <- max(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range)*.8
    p1 + annotate(geom="text", x=c(0, thr)+diff(c(0, thr, max(x)))/2, y=rep(posy, 3), label=c("Listeners", "Communicators", "Speakers"))
}    


#' Regularized in-degree distribution
#' 
#' @param x Degree table
#' @return category vector
fitRegInDegree <- function(x) {
    x <- x[, "Reg_in-degree"]
    set.seed(1)
    fit <- mixGaussianFit(x[x<40])
    fit <- lapply(fit, function(x, pos) {
        if (length(x)!=length(pos)) return(x)
        x[pos]
    }, pos=order(fit$mu))
    class(fit) <- "mgfit"
    par(mai=c(.8, .8, .2, .1))
    plot(fit, x, xlab="Regularized in-degree")
    #legend("topright", c("Listeners", "Comunicators-1", "Communicators-2", "Speakers"), bty="n", lwd=2, col=rainbow(length(fit$mu), .8, .8, end=.8))
    lr <- log(pnorm(x, fit$mu[3], fit$sigma[3], lower.tail=FALSE)/pnorm(x, fit$mu[4], fit$sigma[4]))
    suppressWarnings(thr <- approx(lr, x, xout=0)$y)
    abline(v=thr, col="grey", lty=3)
    text(c(0, thr)+diff(c(0, thr, max(x)))/2, mean(axTicks(2)), c("", "Integrators"), font=2)
    tmp <- predict(fit, x)
    tmp <- apply(tmp, 1, which.max)
    tmp[tmp<4] <- 1
    tmp[tmp==4] <- 2
    names(tmp) <- names(x)
    tmp
}

#' Betweenness distribution
#' 
#' @param x Degree table
#' @return category vector
fitBetweenness <- function(x) {
    x <- x[, "Betweenness"]
    x1 <- sort(x)
    # set.seed(1)
    # thr=450
    # fit <- mixGaussianFit(x[x<thr])
    # fit <- lapply(fit, function(x, pos) {
    #     if (length(x)!=length(pos)) return(x)
    #     x[pos]
    # }, pos=order(fit$mu))
    # class(fit) <- "mgfit"
    # par(mai=c(.8, .8, .2, .1))
    # plot(fit, x, xlab="Network centrality (Betweenness)")
    # abline(v=thr, col="grey", lty=3)
    # text(c(0, thr)+diff(c(0, thr, max(x)))/2, mean(axTicks(2)), c("", "Integrators"), font=2)
    # tmp <- x>=thr
    # names(tmp) <- names(x)
    # tmp
    fit <- fitdistrplus::fitdist(round(x1), "geom")
    pval <- pgeom(round(x1), fit$estimate, lower.tail=FALSE)
    col <- rep("grey70", length(x1))
    pos <- which(p.adjust(pval, "fdr")<.05)
    col[pos] <- "darkred"
    label <- rep("", length(x1))
    label[pos] <- entrez2gene(names(x1)[pos]) 
    p1 <- fitdistrplus::qqcomp(fit, fitpch=20, plotstyle="ggplot") +
        ggplot2::theme_classic(base_size=14) +
        ggplot2::geom_point(size=3, col=col) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle("")
    tmp <- ggplot_build(p1)$data[[1]]
    tmp$label <- label
    tmp <- tmp[, c("x", "y", "label")]
    p1 <- p1 + geom_text_repel(data=subset(tmp, label != ""), aes(x=x, y=y, label=label), inherit.aes=FALSE)
    list(plot=p1, pval=pval[match(names(x), names(pval))], estimate=fit$estimate)
}

#' Betweenness distribution
#' 
#' @param x Degree table
#' @param thr Threshold for betweeness
#' @return category vector
fitBetweennessPlot <- function(x, thr) {
    x <- x[, "Betweenness"]
    set.seed(1)

    data <- data.frame(x=x, names=entrez2gene(names(x)), y=0)
    pos <- which(data$x<thr)
    data$names[pos] <- ""
    pos <- which(data$x>=thr)
    p1 <- ggplot(data=data, aes(x=x)) + theme_classic(base_size = 14) +
        geom_density(fill="grey80", color="grey80") +
        geom_vline(xintercept=thr, linetype="dashed", color="grey50") +
        annotate(geom="point", x=data$x[pos], y=rep(0, length(pos))) +
        geom_text_repel(data=subset(data, names != ""), aes(x=x, y=y, label=names)) +
        xlab("Network centrality (Betweenness)") +
        ylab("Density")
    posy <- max(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range)*.8
    p1 + annotate(geom="text", x=thr+diff(c(thr, max(x)))/2, y=posy, label="Integrators")
}


#' Export interactome
#' 
#' @param reg Interactome
#' @param file Filename
#' @param gz Logical, whether gzip
#' @param genesymbol Logical, whether to transform entrezIDs to gene symbol
exportInteractome <- function(reg, file, gz=TRUE, genesymbol=FALSE) {
    tmp <- cbind(Regulator=rep(names(reg), sapply(reg, function(x) length(x[["tfmode"]]))), Target=unlist(lapply(reg, function(x) names(x[["tfmode"]])), use.names=FALSE), Mode=round(unlist(lapply(reg, function(x) x[["tfmode"]]), use.names=FALSE), 5), Likelihood=round(unlist(lapply(reg, function(x) x[["likelihood"]]), use.names=FALSE), 5))
    if (genesymbol) {
        tmp[, 1] <- entrez2gene(tmp[, 1])
        tmp[, 2] <- entrez2gene(tmp[, 2])
        tmp <- tmp[!is.na(tmp[, 2]), , drop=FALSE]
    }
    if (gz) file <- gzfile(paste0(file, ".gz"), open="w")
    write.table(tmp, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    if (gz) close(file)
}

#' Export igraph network
#' 
#' @param ntw Interactome
#' @param file Filename
#' @param gz Logical, whether gzip
exportNetwork <- function(ntw, ntw1, file, gz=FALSE) {
    tmp1 <- igraph::as_data_frame(ntw)
    tmp2 <- igraph::as_data_frame(ntw1)
    index1 <- paste(tmp1[, 1], tmp1[, 2], sep="x")
    index2 <- paste(tmp2[, 1], tmp2[, 2], sep="x")
    tmp <- cbind("From GeneID"=tmp1[, 1], "From Symbol"=entrez2gene(tmp1[, 1]), "To GeneID"=tmp1[, 2], "To Symbol"=entrez2gene(tmp1[, 2]), "Community network"=index1 %in% index2)
    if (gz) file <- gzfile(paste0(file, ".gz"), open="w")
    write.table(tmp, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    if (gz) close(file)
}

#' Layout based on the communities
#' 
#' @param ntw Network
#' @param stbl Community stability matrix
#' @return layout
radialLayout <- function(ntw, stbl) {
    stbl <- stbl[match(names(V(ntw)), rownames(stbl)), , drop=FALSE]
    r <- seq(0, 2 * pi, length=ncol(stbl)+1)[seq_len(ncol(stbl))]
    cord <- cbind(x=cos(r), y=sin(r))
    #mm <- max(stbl)*.2
    #stbl <- (stbl+mm)/max(stbl+mm)
    #tmp <- sigT(stbl, 10, .2) %*% cord
    stbl1 <- stbl + abs(rnorm(length(stbl), 0, .8))
    stbl1 <- stbl1/max(stbl1)
    tmp <- sigT(stbl1, 15, .2)
    tmp %*% cord
}

#' Layout based on umap
#' 
#' @param ntw Network
#' @param k Integer indicating the nearest neighbors for umap
#' @param x Optional viper matrix
#' @return layout
umapLayout <- function(ntw, k=15, x=NULL) {
    if (is.null(x)) {
        x <- as.matrix(as.undirected(ntw)[])
        d <- 1-jaccard(x)
    }
    else {
        pos <- match(V(ntw)$name, rownames(x))
        x <- x[pos, , drop=FALSE][, pos, drop=FALSE]
        d <- minimumEuclideanDistance(x)
    }
    umap_res <- umap(d, n_neighbors=k, input="dist")$layout
}

#' Minimum Euclidean distance
#' 
#' @param x Matrix
#' @return Distance matrix
minimumEuclideanDistance <- function(x) {
    d1 <- as.matrix(dist(x))
    d2 <- as.matrix(dist(t(x)))
    tmp <- apply(cbind(as.vector(d1), as.vector(d2)), 1, min, na.rm=TRUE)
    matrix(tmp, nrow(d1), ncol(d1), dimnames=list(rownames(d1), colnames(d1)))
}

#' jaccard distance
#' 
#' This function computes the jaccard distance between two vectors or the columns of a matrix
#' 
#' @param x Vector or matrix
#' @param y Optional vector
#' @return Jaccard distance or vector of jaccard distances
jaccard <- function(x, y=NULL) {
    if (is.matrix(x)) {
        tmp <- lapply(1:(ncol(x)-1), function(i, x) {
            apply(filterColMatrix(x, (i+1):ncol(x)), 2, function(x1, x2) {
                jaccard(x1, x2)
            }, x2=x[, i])
        }, x=x)
        tmp <- unlist(tmp, use.names=FALSE)
        ji <- matrix(1, ncol(x), ncol(x))
        ji[lower.tri(ji)] <- tmp
        ji <- t(ji)
        ji[lower.tri(ji)] <- tmp
        rownames(ji) <- colnames(ji) <- colnames(x)
        ji[is.na(ji)] <- 0
        return(ji)
    }
    sum(x&y)/sum(x|y)
}

#' Prune communities object
#' 
#' @param x communities-class object
#' @param genes vector of genes to keep
pruneCommunity <- function(x, genes) {
    if (length(which(x$names %in% genes)) < length(genes)) {
        x$names <- entrez2gene(x$names)
    }
    pos <- match(genes, x$names)
    x$membership <- x$membership[pos]
    x$memberships <- x$memberships[, pos, drop=FALSE]
    x$names <- x$names[pos]
    x$vcount <- length(x$names)
    x
}

#' Export network data for cytoscape
#' 
#' @param ntw1 igraph network object
#' @param layout Graph layout as a 2-columns matrix
#' @param filename Prefix for the export file names
exportCytoscape <- function(ntw1, layout, filename="Network") {
    # Exporting network
    tmp <- as.matrix(as_data_frame(ntw1))
    tmp <- cbind(from=tmp[, 1], to=tmp[, 2], weight=match(tmp[, 3], sort(unique(tmp[, 3]))))
    write.table(tmp, file=file.path(working_dir, paste0(filename, "-network.tsv")),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    # Exporting node attributes
    tmp <- igraph::vertex_attr(ntw1)
    tmp <- cbind(name=tmp$name, size=tmp$size, color=tmp$color, community=tmp$community)
    write.table(tmp, file=file.path(working_dir, paste0(filename, "-nodeAttr.tsv")),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    # Exporting for coordinates layout
    tmp <- cbind(Number=seq_len(nrow(layout)),
                 Name=entrez2gene(rownames(layout)),
                 Cluster=igraph::vertex_attr(ntw1, "community"),
                 layout/max(abs(layout)))
    write.table(tmp, file=file.path(working_dir, paste0(filename, "-coordLayout-nodes.csv")),
                quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
    tmp <- as.matrix(as_data_frame(ntw1))[, 1:3, drop=FALSE]
    tmp[, 1] <- match(tmp[, 1], igraph::vertex_attr(ntw1, "name"))
    tmp[, 2] <- match(tmp[, 2], igraph::vertex_attr(ntw1, "name"))
    write.table(tmp, file=file.path(working_dir, paste0(filename, "-coordLayout-edges.csv")),
                quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
}

#' Integrate GO enrichment results based on the communities
#' 
#' @param go_res Matrix of GO p-values
#' @param ntw Network
#' @param adjust Multiple hypothesis adjustment method
#' @param method String indicating the integration method
#' @return Integrated GO matrix
integrateGOcommunity <- function(go_res, ntw, adjust="none", method=c("mean", "stouffer")) {
    method <- match.arg(method)
    # Community membership
    com <- episcCommunityAnalysis(ntw)
    # Community stability
    stbFET <- communityStabilityFET(ntw)
    # Correct p-values
    go_res <- apply(go_res, 2, p.adjust, method=adjust)
    # p-value to Z-score
    z <- qnorm(go_res, lower.tail=FALSE)
    # Replace -Inf 
    z[z==(-Inf)] <- min(z[is.finite(z)])
    switch(method,
        mean={
            # WStouffer per community
            z_com <- tapply(names(com), com, WRowStoufferSubset, z=z, stb=stbFET)
        },
        stouffer={
            # WMean per community
            z_com <- tapply(names(com), com, WRowMeanSubset, z=z, stb=stbFET)
        }
    )
    z_com <- do.call(cbind, z_com)
    pnorm(z_com, lower.tail=FALSE)
}

#' Weighted Stouffer's integration for a subset of a matrix
#' 
#' @param genes Columns of z (colnames) to include in the analysis
#' @param z Matrix of z scores
#' @param stb Community stability score as z-score to use as weight
#' @return Vector of integrated z-scores
WRowStoufferSubset <- function(genes, z, stb) {
    z <- z[, match(genes, colnames(z)), drop=FALSE]
    stb <- stb[match(genes, names(stb))]
    rowWStouffer(z, sigT(stb, 2, 2), brown="pearson")
}

#' Weighted Mean integration for a subset of a matrix
#' 
#' @param genes Columns of z (colnames) to include in the analysis
#' @param z Matrix of z scores
#' @param stb Community stability score as z-score to use as weight
#' @return Vector of integrated z-scores
WRowMeanSubset <- function(genes, z, stb) {
    z <- z[, match(genes, colnames(z)), drop=FALSE]
    stb <- stb[match(genes, names(stb))]
    rowWMeans(z, sigT(stb, 2, 2))
}

#' Analysis of GOBP categories using topGO
#' 
#' @param plateseq_ges_integ Matrix of gene expression signatures
#' @param myGENE2GO Gene2go mapping
#' @param gobp_mm List of GO sets
#' 
#' @return Matrix of NES
episcTopGOanalysis <- function(plateseq_ges_integ, myGENE2GO = NULL, gobp_mm = NULL) {
    if (length(myGENE2GO) > 0) {
        gobp_mm <- myGENE2GO@BP
    }
    if (length(gobp_mm) == 0) {
        stop("No GO-BP terms have been provided")
    }
    if (is.null(dim(plateseq_ges_integ)) | length(dim(plateseq_ges_integ)) == 1) {
        plateseq_ges_integ <- matrix(plateseq_ges_integ, length(plateseq_ges_integ), 1, dimnames = list(names(plateseq_ges_integ), "GES"))
    }
    go_res <- apply(plateseq_ges_integ, 2, function(x, gobp_mm) {
        sel_genes <- function(x) p.adjust(x, "fdr")<1e-3
        all_genes <- pnorm(abs(x), lower.tail=FALSE)*2
        tgo <- new("topGOdata", ontology="BP", allGenes=all_genes, geneSel=sel_genes, annot=annFUN.gene2GO, gene2GO=gobp_mm, nodeSize=10)
        res <- topGO::runTest(tgo, algorithm ="weight01", statistic = "ks")
        score(res)
    }, gobp_mm=gobp_mm)
    go_res
}

#' Analysis of GOBP categories using aREA and shadow analysis
#' 
#' @param plateseq_ges_integ Matrix of gene expression signatures
#' @param myGENE2GO Gene2go mapping
#' @param gobp_mm List of GO sets
#' @param abs Logical, whether the enrichment should be on the absolute value of the signature
#' 
#' @return Matrix of NES
episcGOBPshadow <- function(plateseq_ges_integ, myGENE2GO = NULL, gobp_mm = NULL, abs = TRUE) {
    if (length(myGENE2GO) > 0) {
        gobp_mm <- myGENE2GO@BP
    }
    if (length(gobp_mm) == 0) {
        stop("No GO-BP terms have been provided")
    }
    if (is.null(dim(plateseq_ges_integ)) | length(dim(plateseq_ges_integ)) == 1) {
        plateseq_ges_integ <- matrix(plateseq_ges_integ, length(plateseq_ges_integ), 1, dimnames = list(names(plateseq_ges_integ), "GES"))
    }
    gobp_mm <- gobp_mm[names(gobp_mm) %in% rownames(plateseq_ges_integ)]
    gobp_mm <- cbind(unlist(gobp_mm, use.names=FALSE), rep(names(gobp_mm), sapply(gobp_mm, length)))
    # Generate the gene-sets
    gobp_reg <- tapply(gobp_mm[, 2], gobp_mm[, 1], function(x) {
        tfmode <- rep(1, length(x))
        names(tfmode) <- x
        list(tfmode=tfmode, likelihood=rep(1, length(tfmode)))
    })
    # Remove small gene-sets
    gobp_reg <- gobp_reg[sapply(gobp_reg, function(x) {
        length(x[["tfmode"]])
    })>=15]
    expmat <- plateseq_ges_integ
    if (abs) {
        expmat <- abs(plateseq_ges_integ)
    }
    #save(expmat, gobp_reg, version=2, file="rawdata4gobp.rda")
    # Enrichment analysis
    res <- apply(expmat, 2,  function(x, gobp_reg) {
        res <- viper::msviper(x, gobp_reg, minsize=15, ges.filter=FALSE, verbose=FALSE)
        res <- viper::shadow(res, verbose=FALSE)$es$nes
        res <- res[match(names(gobp_reg), names(res))]
        names(res) <- names(gobp_reg)
        res
    }, gobp_reg=gobp_reg)
    if (is.null(ncol(res))) {
        res <- res[!is.na(res)]
    } else {
        res <- res[rowSums(is.na(res))==0, , drop=FALSE]
    }
    return(res)
}


#' Integrate GO enrichment results based on speaker, communicator, listener
#' 
#' @param go_res Matrix of GO p-values
#' @param spkr speaker category
#' @param adjust p-value correction method
#' @param method String indicating the integration method
#' @return Integrated GO matrix
integrateGOspeaker <- function(go_res, spkr, adjust="none", method=c("mean", "stouffer")) {
    method <- match.arg(method)
    # p-value correction
    go_res <- apply(go_res, 2, p.adjust, method=adjust)
    # p-value to Z-score
    z <- qnorm(go_res, lower.tail=FALSE)
    # Replace -Inf 
    z[z==(-Inf)] <- min(z[is.finite(z)])
    # filter and sort z
    z <- z[, match(names(spkr), colnames(z))]
    switch(method,
        mean={
            # WMean per community
            z_spkr <- tapply(names(spkr), spkr, WRowMeanSubset, z=z, stb=spkr/spkr)
        },
        stouffer={
            # WStouffer per community
            z_spkr <- tapply(names(spkr), spkr, WRowStoufferSubset, z=z, stb=spkr/spkr)
        }
    )
    z_spkr <- do.call(cbind, z_spkr)
    pnorm(z_spkr, lower.tail=FALSE)
}

#' Get GO terms from GOID
#'
#' @param goid Vector of GO IDs
#' @return Vector of terms
goID2Term <- function(goid) {
    if (length(goid)==0) return(NULL)
    term <- Term(goid)
    vapply(term, capitalize, character(1))
}

#' Export GO table
#' 
#' @param go_res Matrix of GO results (p-value) to export
#' @param go_metadata Data.frame with Category, Size and Level
#' @param file Filename without extension
#' @param pval P-value threshold
#' @param adjust Multiple hypothesis correction method
#' @return Nothing a file is saved to the data dir
exportGOtable <- function(go_res, go_metadata, file, pval=.01, adjust="fdr") {
    go_res_sel <- filterGOmatrix(go_res, pval=pval, adjust=adjust)
    go_res_sel <- sortGOtable(go_res_sel, pval=pval)
    tmp <- go_metadata[match(rownames(go_res_sel), go_metadata$Category), , drop=FALSE][, -1, drop=FALSE]
    rownames(go_res_sel) <- goID2Term(rownames(go_res_sel))
    tmp <- cbind(GOBP=rownames(go_res_sel), tmp, signif(go_res_sel, 3))
    write.table(tmp, file=file.path(working_dir, paste0(file, ".tsv")), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}

#' Filter GO matrix to keep significant results only
#' 
#' @param go_res GO results matrix (p-values)
#' @param pval Number indicating the p-value threshold
#' @param adjust Method for multiple hypothesis correction

filterGOmatrix <- function(go_res, pval=1e-2, adjust="fdr") {
    res <- apply(go_res, 2, p.adjust, method=adjust)
    res[rowSums(res < pval)>0, , drop=FALSE]
}

#' Sort a matrix of GO results
#' 
#' @param x Matrix of GO results (p-values)
#' @param pval P-value threshold
sortGOtable <- function(x, pval) {
    tmp <- apply(x, 2, function(x, pval) {
        x <- x*(x<pval)
        x[x==0] <- 1
        x
    }, pval=pval)
    pos <- do.call(order, data.frame(mult=-rowSums(x<pval), tmp))
    x[pos, , drop=FALSE]
}

#' Get MsigDB annotation
getMsigDB <- function() {
    if (exists("_episc_msigdb", envir = globalenv())) {
        msigdb <- get("_episc_msigdb", envir=globalenv())
    } else {
        msigdb <- readRDS(file.path(data_dir, "msig_mouse_entrez_v7.rds"))
        assign("_episc_msigdb", msigdb, envir=globalenv())
    }
    msigdb
}

#' Enrichment of MsigDB gene-sets
#' 
#' @param expmat Expression or viper matrix
#' @param minsize minimum size for the gene-sets
msigdbAnalysis <- function(expmat, minsize=10) {
    msigdb <- getMsigDB()
    # Filter msigdb
    msigdb <- lapply(msigdb, function(x, genes) {
        x[x %in% genes]
    }, genes=rownames(expmat))
    # Remove small gene-sets
    msigdb <- msigdb[sapply(msigdb, length) >= minsize]
    # Running enrichment analysis
    sREA(expmat, msigdb)
}

#' Enrichment of MsigDB gene-sets with shadow analysis
#' 
#' @param expmat Expression or viper matrix
#' @param minsize minimum size for the gene-sets
msigdbAnalysisShadow <- function(expmat, minsize=10) {
    res <- readRDS(file.path(data_dir, "msigdb_results.rds"))
    return(res)
    # Following code is for illustrative purposes
    msigdb <- getMsigDB()
    # Filter msigdb
    msigdb <- lapply(msigdb, function(x, genes) {
        tmp <- x[x %in% genes]
        tfmode <- rep(1, length(tmp))
        names(tfmode) <- tmp
        list(tfmode=tfmode, likelihood=rep(1, length(tfmode)))
    }, genes=rownames(expmat))
    # Remove small gene-sets
    msigdb <- msigdb[sapply(msigdb, function(x) {
        length(x[["tfmode"]])
    }) >= minsize]
    # Running enrichment analysis with shadow correction
    res <- apply(expmat, 2,  function(x, msigdb) {
        res <- viper::msviper(x, msigdb, minsize=minsize, ges.filter=FALSE, verbose=FALSE)
        res <- viper::shadow(res, verbose=FALSE)$es$nes
        res <- res[match(names(msigdb), names(res))]
        names(res) <- names(msigdb)
        res
    }, msigdb=msigdb)
    res[rowSums(is.na(res))==0, , drop=FALSE]
}

#' Integrate MsigDB enrichment results based on the communities
#' 
#' @param nes Matrix of MsigDB NES
#' @param ntw Network
#' @param adjust Multiple hypothesis adjust method
#' @param method String indicating the integration method
#' @return Integrated NES matrix
integrateMsigDBcommunity <- function(nes, ntw, adjust="fdr", method=c("mean", "stouffer")) {
    method <- match.arg(method)
    # Multiple-hypothesis correction
    if (adjust != "none") {
        pval <- pnorm(abs(nes), lower.tail=FALSE)*2
        pval <- p.adjust(pval, method=adjust)
        nes <- qnorm(pval/2, lower.tail=FALSE) * sign(nes)
    }
    # Community membership
    com <- episcCommunityAnalysis(ntw)
    # Community stability
    stbFET <- communityStabilityFET(ntw)
    switch(method,
        mean={
            # RowMeans
            z_com <- tapply(names(com), com, WRowMeanSubset, z=nes, stb=stbFET)
        },
        stouffer={
            # Stouffer
            z_com <- tapply(names(com), com, WRowStoufferSubset, z=nes, stb=stbFET)
        }
    )
    do.call(cbind, z_com)
}

#' Integrate MsigDB enrichment results based on speaker, communicator, listener
#' 
#' @param nes Matrix of MsigDB NES
#' @param spkr speaker category
#' @param adjust P-value correction method
#' @param method String indicating the integration method
#' @return Integrated NES matrix
integrateMsigDBspeaker <- function(nes, spkr, adjust="fdr", method=c("mean", "stouffer")) {
    method <- match.arg(method)
    # Multiple-hypothesis correction
    if (adjust != "none") {
        pval <- pnorm(abs(nes), lower.tail=FALSE)*2
        pval <- p.adjust(pval, method=adjust)
        nes <- qnorm(pval/2, lower.tail=FALSE) * sign(nes)
    }
    # filter and sort z
    nes <- nes[, match(names(spkr), colnames(nes))]
    switch(method,
        mean={
            # WMean per community
            z_spkr <- tapply(names(spkr), spkr, WRowMeanSubset, z=nes, stb=spkr/spkr)
        },
        stouffer={
            # WStouffer per community
            z_spkr <- tapply(names(spkr), spkr, WRowStoufferSubset, z=nes, stb=spkr/spkr)
        }
    )
    do.call(cbind, z_spkr)
}


#' Export MsigDB table
#' 
#' @param nes Matrix of GO results (p-value) to export
#' @param metadata Data.frame with Category and Size
#' @param file Filename without extension
#' @param thr P-value threshold
#' @param adjust Multiple hypothesis correction method
#' @param direction Whether the enrichment is 2-tailed or 1-tailed
#' @return Nothing a file is saved to the data dir
exportMsigDBtable <- function(nes, metadata, file, thr=.01, adjust="fdr", direction=TRUE) {
    if (direction) {
        pval <- pnorm(abs(nes), lower.tail=FALSE)*2
    } else {
        pval <- pnorm(nes, lower.tail=FALSE)
    }
    pval <- filterGOmatrix(pval, pval=thr, adjust=adjust)
    res_sel <- sortMsigDBtable(nes[match(rownames(pval), rownames(nes)), , drop=FALSE], pval=pval, thr=thr)
    if (direction) {
        tmp <- cbind(getPositionsFromString(rownames(res_sel), "_", 1:2), metadata$Size[match(rownames(res_sel), metadata$Category)], round(res_sel[, seq_len(ncol(nes)), drop=FALSE], 2),
                 signif(res_sel[, -seq_len(ncol(nes)), drop=FALSE], 3))
        colnames(tmp) <- c("Database", "Pathway", "Size", paste0("NES (", colnames(nes), ")"), paste0("FDR (", colnames(nes), ")"))
    } else {
        tmp <- cbind(getPositionsFromString(rownames(res_sel), "_", 1:2), metadata$Size[match(rownames(res_sel), metadata$Category)], signif(res_sel[, -seq_len(ncol(nes)), drop=FALSE], 3))
        colnames(tmp) <- c("Database", "Pathway", "Size", paste0("FDR (", colnames(nes), ")"))
    }
    write.table(tmp, file=file.path(working_dir, paste0(file, ".tsv")), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}

#' Sort a matrix of NES and p-values for MsigDB results
#' 
#' @param nes Matrix of NES
#' @param pval Matrix o p-values
#' @param thr P-value threshold
sortMsigDBtable <- function(nes, pval, thr) {
    tmp <- apply(pval, 2, function(x, thr) {
        x <- x*(x<thr)
        x[x==0] <- 1
        x
    }, thr=thr)
    pos <- do.call(order, data.frame(mult=-rowSums(pval<thr), tmp, rowMeans(nes, na.rm=TRUE)))
    cbind(nes[pos, , drop=FALSE], pval[pos, , drop=FALSE])
}

#' Plot UMAP representation of the expression matrix
plotExpmatUmap <- function(expmat, annot) {
    set.seed(1)
    pos <- umap::umap(scale(t(expmat)))$layout
    col1 <- rainbow(length(levels(annot$Morphogen)), s=.7, v=.8, start=0, end=.8, alpha=.8)
    col <- col1[match(annot$Morphogen, levels(annot$Morphogen))]
    par(mai=c(.5, .5, .2, .2))
    plot(pos[, 1], pos[, 2], col=col, pch=20, axes=FALSE, xlab="Dimension 1", ylab="Dimension 2")
    tmp <- levels(annot$Morphogen)
    tmp[tmp == ""] <- "Mock"
    legend("topleft", tmp, fill=col1, bty="n")
}

#' Plot UMAP representation of the expression matrix
plotExpmatUmapGG <- function(expmat, annot) {
    set.seed(pi)
    pos <- umap::umap(scale(t(expmat)))$layout
    morph <- as.vector(annot$Morphogen)
    morph[morph==""] <- "Mock"
    tmp <- unique(morph)
    tmp <- c("Mock", "CMAF", tmp[!(tmp %in% c("Mock", "CMAF"))])
    tmp <- tibble::tibble(x=pos[, 1], y=pos[, 2], Morphogen=factor(morph, levels=tmp),
                          Strain = annot$Strain)
    ggplot(tmp, aes(x=x, y=y, colour=Morphogen, shape = Strain)) + theme_classic(base_size = 12) +
        xlab("UMAP dimension 1") + ylab("UMAP dimension 2") + geom_point(alpha=.5, size=2)
}

#' Heatmap for the differentiation trajectories
#' 
#' @param vpmat Viper matrix
#' @param output Optional string indicating the output file for the data
differentiationHeatmap <- function(vpmat, output=NULL) {
    annot <- getPositionsFromString(colnames(vpmat), "-", 1:2)
    annot[, 2] <- paste0(annot[, 2], "h")
    tmp <- vpmat
    rownames(tmp) <- entrez2gene(rownames(tmp))
    colnames(tmp) <- annot[, 2]
    col_fun = circlize::colorRamp2(c(-10, 0, 10), c("royalblue", "grey94", "firebrick2"))
    # Clustering
    hc <- hclust(dist(tmp)) %>%
        as.dendrogram() %>%
        stats::reorder(wts=rowMeans(tmp, na.rm=TRUE))
    fig <- ComplexHeatmap::Heatmap(tmp, name="NES", col=col_fun,
                                   rect_gp = gpar(col = "white", lwd = 1),
                                   cluster_columns=FALSE, cluster_rows=hc,
                                   column_split=factor(annot[, 1], levels=unique(annot[, 1])),
                                   column_names_side="top")
    if (!is.null(output)) {
        tmp1 <- tmp[as.hclust(hc)$order, , drop=FALSE]
        tmp1 <- cbind(rownames(tmp1), round(tmp1, 3))
        colnames(tmp1) <- c("Symbol", colnames(vpmat))
        write.table(tmp1, file=output, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    }
    fig
}

#' Heatmap for the differentiation trajectories for each community
#' 
#' @param vpmat Viper matrix
#' @param ntw Network
#' @param output Optional file name to output data
differentiationHeatmapCommunity <- function(vpmat, ntw, output=NULL) {
    # Community membership
    com <- episcCommunityAnalysis(ntw)
    # Community stability
    stbFET <- communityStabilityFET(ntw)
    # WMean per community
    z_com <- tapply(names(com), com, WRowMeanSubset, z=t(vpmat), stb=stbFET)
    tmp <- do.call(rbind, z_com)
    annot <- getPositionsFromString(colnames(vpmat), "-", 1:2)
    annot[, 2] <- paste0(annot[, 2], "h")
    colnames(tmp) <- annot[, 2]
    rownames(tmp) <- paste0("Community ", rownames(tmp))
    col_fun = circlize::colorRamp2(c(-10, 0, 10), c("royalblue", "grey94", "firebrick2"))
    if (!is.null(output)) {
        tmp1 <- cbind(rownames(tmp), round(tmp, 3))
        colnames(tmp1) <- c("Community", colnames(vpmat))
        write.table(tmp1, file=output, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    }
    ComplexHeatmap::Heatmap(tmp, name="NES", col=col_fun,
                            rect_gp = gpar(col = "white", lwd = 1),
                            cluster_columns=FALSE, cluster_rows=FALSE,
                            column_split=factor(annot[, 1], levels=unique(annot[, 1])),
                            row_names_side="left", column_names_side="top")
}

#' Clean interactome to contain only known genes
#' 
#' @param reg Interactome
#' @return Clean interactome
cleanInteractomeKnownGenes <- function(reg) {
    reg <- reg[!is.na(entrez2gene(names(reg)))]
    targets <- unique(unlist(lapply(reg, function(x) names(x[["tfmode"]])), use.names=FALSE))
    targets <- targets[!is.na(entrez2gene(targets))]
    reg <- lapply(reg, function(x, targets) {
        pos <- which(names(x[["tfmode"]]) %in% targets)
        list(tfmode=x[["tfmode"]][pos], likelihood=x[["likelihood"]][pos])
    }, targets=targets)
    reg <- reg[sapply(reg, function(x) length(x[["tfmode"]]))>0]
    reg
}

#' Map Illumina probes to entrezID
#' 
#' @param x Vector of probe IDs
#' @param data_dir String indicating the data directory
#' @return Vector of entrezIDs
probe2entrez <- function(x, data_dir) {
    ref <- read.table(file.path(data_dir, "mWG6v2.txt.gz"), header = TRUE, sep = "\t")
    ref$ENTREZ_GENE_ID[match(x, ref$ProbeID)]
}

#' Community analysis
#' 
#' @param ntw igraph network object
#' @return Vector of community assignments
episcCommunityAnalysis <- function(ntw) {
    community <- igraph::cluster_louvain(igraph::as.undirected(ntw))
    com_membership <- igraph::membership(community)
    # Getting the saved community info
#    tmp <- read.table(file.path(data_dir, "communities.tsv"), sep="\t", header=TRUE, as.is=TRUE)
#    com_membership <- tmp$Community
#    names(com_membership) <- tmp$EntrezID
    com_membership
}

#' MsigDB metadata
#' 
#' @param categ Vector of MsigDB category names
#' @param genes Vector of genes in the signatures
#' 
#' @return data.frame with size per category
msigdbMetadata <- function(categ, genes) {
    msigdb <- getMsigDB()
    # Filter categories
    msigdb <- msigdb[match(categ, names(msigdb))]
    # Filter msigdb
    msigdb <- lapply(msigdb, function(x, genes) {
        x[x %in% genes]
    }, genes=genes)
    size <- vapply(msigdb, length, numeric(1))
    data.frame(Category=names(msigdb), Size=size)
}

GOtable <- function(gores, ges, thr = 0.01, adjust = "fdr") {
    pval <- p.adjust(pnorm(abs(gores), lower.tail = FALSE) * 2, adjust)
    go_md <- getGOmetadata(rownames(gores), names(ges))
    tbl <- data.frame(GOID = go_md$Category, GOterm = goID2Term(go_md$Category),
        Size = go_md$Size, Level = go_md$Level, NES = round(gores[match(go_md$Category,
        rownames(gores)), ], 2), FDR = signif(pval[match(go_md$Category, rownames(gores))], 3))
    tbl <- tbl[order(tbl$NES), ]
    tbl[tbl$FDR < thr, ]
}

#' GOBP metadata
#' 
#' @param categ Vector of GOIDs
#' @param genes Vector of genes in the signatures
#'
#' @return data.frame with Category, Size and Level
getGOmetadata <- function(categ, genes) {
#    EntrezGene<-ViSEAGO::EntrezGene2GO()
#    myGENE2GO<-ViSEAGO::annotate("10090", EntrezGene)
#    gobp_mm <- myGENE2GO@BP
#    save(gobp_mm, version=2, file=file.path(data_dir, "gobp_mm.rda"))
    load(file.path(data_dir, "gobp_mm.rda"))
    gobp_mm <- gobp_mm[names(gobp_mm) %in% genes]
    gobp_mm <- cbind(unlist(gobp_mm, use.names=FALSE), rep(names(gobp_mm), sapply(gobp_mm, length)))
    gobp_mm <- gobp_mm[gobp_mm[, 1] %in% categ, , drop=FALSE]
    size <- table(gobp_mm[, 1])
    level <- sapply(names(size), goLevel, parents=as.list(GO.db::GOBPPARENTS, use.names=FALSE))
    data.frame(Category=names(size), Size=as.vector(size), Level=level)
}

goLevel <- function(x, parents=NULL, sc=0) {
    if (is.null(parents))
        parents <- as.list(GO.db::GOBPPARENTS)
    pos <- which(names(parents) %in% x)
    while(length(pos)>0) {
        x <- unlist(parents[pos], use.names=FALSE)
        sc <- sc + 1
        pos <- which(names(parents) %in% x)
    }
    sc
}

#' CRISPR-mediated knock-out gene expression signatures
#' 
#' @param data_dir String indicating the path to the data
koGES <- function(data_dir) {
    rc <- readRDS(file.path(data_dir, "KO_counts.rds"))
    expmat <- DEtransform(rc)
    ges <- cbind(expmat[, 1:4]-rowMeans(expmat[, 5:6], na.rm=TRUE))
    ges <- cbind(ges, Empty=rowMeans(expmat[, 5:6], na.rm=TRUE)-expmat[, "Cas9"])
    ges <- cbind(ges, Cas9=expmat[, "Cas9"]-expmat[, "Mock"])
    return(ges)
}

#' Integrate list of matrices
#' 
#' @param x List of matrices
#' @return Matrix
#' @export
concatenateMatrixList <- function(x, method=c("union", "intersection")) {
    # Check argument
    checkmate::assertList(x, types="matrix", min.len=1)
    method <- match.arg(method)
    switch(method,
           union={genes <- unique(unlist(lapply(x, rownames), use.names=FALSE))},
           intersection={
               genes <- table(unlist(lapply(x, rownames), use.names=FALSE))
               genes <- names(genes)[genes==max(genes)]
           })
    do.call(cbind, lapply(x, matrixOrderRowsByName, names=genes))
}

#' Matrix order rows by name
#' 
#' @param x Matrix
#' @param names Vector of strings
#' @return Matrix with ordered rows
#' @export
matrixOrderRowsByName <- function(x, names) {
    # Assert arguments
    checkmate::assertMatrix(x, mode="numeric", row.names="named", min.rows=1)
    checkmate::assertCharacter(names, min.len=1, any.missing=FALSE, min.chars=1)
    # Order
    x <- x[match(names, rownames(x)), , drop=FALSE]
    rownames(x) <- names
    return(x)
}

heatmapCorrelation <- function(x) {
    x <- cor(x, use="pairwise.complete.obs")
    col <- circlize::colorRamp2(c(-1, 0, 1), c("slateblue3", "grey93", "orange3"))
    Heatmap(x, name="R", col=col, row_order=seq_len(nrow(x)),
            column_order=seq_len(ncol(x)), row_names_side="left", column_names_side="top",
            rect_gp = gpar(col = "white", lwd = 2))
}

heatmapViperSimilarity <- function(x, nn=100) {
    x <- viperSimilarity(x, nn=nn) %>% scale
    class(x) <- "matrix"
    col <- circlize::colorRamp2(c(-1, 0, 1), c("royalblue", "grey93", "firebrick2"))
    Heatmap(x, name="RE", col=col, row_order=seq_len(nrow(x)),
            column_order=seq_len(ncol(x)), row_names_side="left", column_names_side="top",
            rect_gp = gpar(col = "white", lwd = 2))
}

btwPvalTable <- function(x, pval) {
    x <- x[, "Betweenness"]
    tbl <- cbind(Gene=entrez2gene(names(x)), Betweenness=round(x),
                 "p-value"=signif(pval, 3), FDR=signif(p.adjust(pval, "fdr"), 3))
    tbl[order(pval)[1:10], ]
}

btwThr <- function(x, pval, adjust="none", pthr=.05) {
    x <- x[, "Betweenness"]
    tbl <- cbind(round(x), p.adjust(pval, adjust))
    suppressWarnings(btw <- approx(tbl[, 2], tbl[, 1], xout=pthr)$y)
    btw
}

#' Read community membership from file
#' 
#' @param path String indicating the directory
readCommunityFromFile <- function(path) {
    tmp <- strsplit(readLines(file.path(path, "communities.tsv"))[-1], "\t")
    tmp <- do.call(rbind, tmp)
    comm <- as.numeric(tmp[, 2])
    names(comm) <- tmp[, 1]
    return(comm)
}
