
#' @title Computing quantile normalised MAD for gene expression
#'To produce a quantile normalised MAD values for genes expression
#' @param geneMat Expression matrix where columns are cells and rows are genes
#' @param celltype Vector of cell type labels, must correspond to columns in geneMat
#' @param zeroPerc Percentage of lowly expressed cells that need to be removed
#' @param OutlierRemoval Removal of outlier gene expression set at 4x of MAD values
#' Default set to TRUE
#' @param QuantNo Total number of quantiles the gene expression distribution is to be spli into
#' @param method Different analysis building upon the quantile methodology:
#' "std", returns the quantile NMAD values for each gene, for each cell type and a list of data frames.
#' "SG" calculates genes with and ranks genes with based on how stable the gene expressions are, genes are ranked from most to least stable, where there higher the ranking the more stable the genes are.
#' "DS" calculates genes that are differentially stably expressed between cell types, using two-sample Wilcoxon test.
#' "NoPval", calculates differentially stable genes by comparing the ranking between genes expression from different cell types.
#' "SV" calculate and returns a ranking of genes that are stably variable either between cell types or between datasets
#' @return with std returns the quantile NMAD values
#'
#' @import data.table
#' @import Hmisc
#' @import Seurat
#' @import pbapply
#' @import dplyr
#' @import progress
#'
#'

scQuant = function(geneMat = geneMat, celltype = celltype,
                   zeroPerc = zeroPerc,
                   OutlierRemoval = TRUE,
                   QuantNo =  QuantNo,
                   method = c("std", "SG", "SV", "DS",
                              "DS.Pval")){


    message("Removing cell types with fewer than 50 cells")

    CT = as.data.frame(table(celltype))

    CT = CT[CT$Freq >= 50 ,]

    celltype = as.character(celltype)
    celltype = celltype[celltype %in% CT$celltype]

    message("Subsetting data by Cell Type")
    geneMat_list = pblapply(names(table(celltype)), function(x) {
        (geneMat[, celltype == x] )})

    names(geneMat_list) = names(table(celltype))
    #
    #print message to show removal of genes with zeros
    message("removing genes with excessive zeros")
    zeRM = pblapply(geneMat_list, function(x){
        RMze = x[rowMeans(x == 0) < zeroPerc ,]})

    if (OutlierRemoval == TRUE) {
        message("removing outliers")
        OutRM = pblapply(zeRM, function(x){
            apply(x, 1, function(w){
                wt = RMOutlier(w, STD = 4)
                wt = as.list(wt)
                return(wt)})
        })
    }

    if(OutlierRemoval == FALSE) {
        OutRM = lapply(zeRM, function (x){
            apply(x, 1, function(w){
                wt = w
                wt = as.list(wt)
                return(wt)})
        })
    }

    OutRM.unlist = lapply(OutRM, lapply, unlist)

    message("Running scQuant")
    SetOut = pblapply(OutRM.unlist, lapply, function(x){
        dt = as.data.frame(table(x))
        if(nrow(dt) == 1){
            t24 = as.data.frame(t(c(0,0,0,0)))
            colnames(t24) = c(seq(1:ncol(t24)))
            return(t24)
        }else{
            t22 = split(x, cut(x, unique(quantile(x, (0:4/4),
                                                  names = FALSE),
                                         include = TRUE)))
            t21 = as.data.frame(t(unlist(lapply(t22, NMAD))))
            colnames(t21) = c(seq(1:ncol(t21)))
            return(t21)
        }
    })

    Quant_out.tb = lapply(SetOut, function(x){
        list_data = as.data.frame(data.table::rbindlist(x, fill = T))
        list_data[is.na(list_data)] = 0
        list_data$genes = names(x)
        return(list_data)
    })

    if (method == 'std') {
        return(Quant_out.tb)
    }

    if (method == "SG"){
        SG = scSG (dat.list = Quant_out.tb, QuantNo = QuantNo)
        return(SG)
    }
    if (method == "SV"){
        SV = scSV(dat.list = Quant_out.tb, QuantNo = QuantNo)
        return(SV)
    }
    if (method == "DS.Pval"){
        DS = scDS.Pval(dat.list = Quant_out.tb)
        return(DS)
    }
    if (method == "DS"){
        DS = scDS(dat.list = Quant_out.tb, QuantNo = QuantNo)
    }
}

#' scSG an method for scQuant to calculate and rank stable genes. This function can be run within scQuant or be used as a standalone, to be run using the results from scQuant (method = "std")
#' @param dat.list Output from scQuant (method = "std")
#' @param QuantNo Number of quantiles that were used to generate dat.list
#' @return ranking of stably expressed genes across different cell types


#For calculating SG genes
scSG = function (dat.list = dat.list, QuantNo = QuantNo){
    message("Performing SG calculation")
    scSGenens = pblapply(dat.list, function(v){
            gen.name = v$genes
            v$genes = NULL
            v$DG.mean = base::rowMeans(v[,1:ncol(v)])
            v$DG.sum = base::rowSums(v[,1:ncol(v)])
            v$genes = gen.name
            v$DG.Rank = 1 - rank(v$DG.mean)/(length(v$DG.mean) + 1)
        return(v)
    })
}

#' scSV an method for scQuant to calculate and rank genes that are similarly variable. This function can be run within scQuant or be used as a standalone, to be run using the results from scQuant (method = "SV")
#' @param dat.list Output from scQuant (method = "std")
#' @param QuantNo Number of quantiles that were used to generate dat.list
#' @import dplyr
#' @import purrr
#' @import data.table
#' @import pbapply
#' @importFrom philentropy distance
scSV = function (dat.list = dat.list, QuantNo = QuantNo){

    dtD = dat.list %>% purrr::reduce(inner_join, by ='genes')
    rownames(dtD) = dtD$genes
    dtD$genes = NULL
    dtD1 = dtD[rowMeans(dtD == 0) <= 0.5 ,]
    colnames(dtD1) = NULL
    scSV_split = apply(dtD1, 1 , function(x){
        scSV_split = as.data.frame(split(x, ceiling(seq_along(x) / QuantNo)))
    })
    message("Perfoming SV")
    SV.l = pblapply(scSV_split, function(d){
        xw= t(d)
        d1 = philentropy::distance(xw,
                                   method = "euclidean", mute.message = T)
        d1[lower.tri(d1)] = NA
        SV.d = mean(d1, na.rm = TRUE)
        return(SV.d)
    })
    SV = as.data.frame(t(rbindlist(list(SV.l))))
    colnames(SV) = "SV.values"
    SV$Rank = 1 - rank(SV$SV.values)/(length(SV$SV.values) + 1)
    return(SV)
}

#' scDS.Pval an method for scQuant to calculate genes that are differentially stable (DS) between cell types. This function can be run within scQuant or be used as a standalone, to be run using the results from scQuant (method = "DS.Pval")
#' @param dat.list Output from scQuant (method = "std")
#' @import progress
#' @import data.table
#' @importFrom stats complete.cases
#For calculating DS genes with p.value
scDS.Pval = function(dat.list = dat.list){
    SetUni = lapply(dat.list, function(unidat){
        genenames = unidat$genes
        unidat$genes = NULL
        vars = combn(colnames(unidat), 2)
        vars = cbind(vars, vars[2:1,])
        done = sapply(seq_len(ncol(vars)),
                      function(x) unidat[, vars[1, x]] / unidat[, vars[2, x]])

        done = abs(1 - done)
        done[!is.finite(done)] = 0
        done = as.data.frame(done)

        colnames(done) = paste(vars[1, ], vars[2, ], sep="/")
        done$genes = genenames
        return(done)
    })


    vars2 = combn(names(SetUni), 2)


    pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                           total = ncol(vars2),
                           complete = "+",
                           incomplete = "_",
                           current = "+",
                           clear = FALSE,
                           width = 100)

    DSout = data.frame()
    message("Performing DS")
    for(i in 1:ncol(vars2)){
        pb$tick()
        t1 = SetUni[[vars2[1, i]]]
        t2 = SetUni[[vars2[2, i]]]
        t3 = merge(t1, t2, by = "genes", all.x = TRUE)
        rownames(t3) = t3$genes
        t3$genes = NULL
        t3 = t3[complete.cases(t3), ]
        t3_split = apply(t3, 1 , function(x){
            t3_split = as.data.frame(split(x, ceiling(seq_along(x) / (ncol(t1)-1))))
        })

        y3 = lapply(t3_split, function(w){
            w1 = wilcox.test(w$X1, w$X2,
                             paired = TRUE,
                             exact = FALSE,
                             alternative = "two.sided")$p.value
        })

        y3 = t(as.data.frame(y3))
        colnames(y3) = "p.value"
        y3 = as.data.frame(y3)
        y3[is.nan(y3$p.value),] <- 1
        y3$celltype.1 = vars2[1,i]
        y3$celltype.2 = vars2[2, i]
        y3$genes = rownames(y3)
        rownames(y3) = NULL
        y3$p.adj = p.adjust(y3$p.value, method = "BH")
        DSout = rbindlist(list(DSout, y3))
    }
    DSout = DSout[complete.cases(DSout), ]

    return(DSout)
}

#' scDS.Pval an method for scQuant to calculate genes that are differentially stable (DS) between cell types. This function can be run within scQuant or be used as a standalone, to be run using the results from scQuant (method = "DS")
#' @param dat.list Output from scQuant (method = "std")
#' @param QuantNo Number of quantiles that were used to generate dat.list
#' @import progress
#' @import data.table

#For calcualting DS genes
scDS = function(dat.list, QuantNo){
    xt = scSG(dat.list, QuantNo)

    ntsel = lapply(xt, function(x){
        sel = x %>% dplyr::select(DG.Rank, genes)
    })


    vars2 = combn(names(ntsel), 2)
    vars3 = rbind(vars2[2,], vars2[1,])

    vars4 = cbind(vars2, vars3)

    pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                          total = ncol(vars4),
                          complete = "+",
                          incomplete = "_",
                          current = "+",
                          clear = FALSE,
                          width = 100)

    DSout = data.frame()
    message("Performing DS")
    for(i in 1:ncol(vars4)){
        #pb$tick()
        t1 = ntsel[[vars4[1, i]]]
        t2 = ntsel[[vars4[2, i]]]
        t3 = merge(t1, t2, by = "genes", all.x = TRUE)

        t3$diff.abs = t3$DG.Rank.x - t3$DG.Rank.y
        t3$celltype.1 = vars4[1, i]
        t3$celltype.2 = vars4[2, i]
        t4 = t3 %>% top_frac(0.2, diff.abs)
        DSout = rbindlist(list(DSout, t4))
    }
    return(DSout)
}


#' NMAD is for calculating normalised MAD, where the MAD value is divided by the mean. This allows for comparison of NMAD results between different datasets.
#' @param x Vector of values
#'
#For calcualting normalised MAD values where x is a vector
NMAD = function(x) {
    mad(x[x!= 0])/mean(x[x!= 0])}

#' RMOutlier is for the removal of outliers within gene expression datasets. Users can alter the cutoffs for what is deemed to be an outlier using the STD variable. The STD variable denotes how many standard deviations from teh median.
#' @param x Vector of values
#' @param STD number of deviations from median. Default set at 4
#For the removal of outliers, where x is a vector and STD is the median deviations. Default at 4.
RMOutlier = function(x, STD = 4){
    x2 = x
    x1 = as.vector(x2[!is.na(x2)])
    w2 = as.vector(x2[is.na(x2)])
    m <- median(x1)
    if (m == 0) {
        x2[is.na(x2)] = 0
        return(x2)
    }else if (m == as.numeric(names(which.max(table(x2))))){
        x2[is.na(x2)] = 0
        return(x2)
    }else{
        deviations <- abs(x1 - m)
        lowerMAD <- 1.4826 * median(deviations[x1 <= m])
        upperMAD <- 1.4826 * median(deviations[x1 >= m])

        xt = (x1[x1 >= m - STD * lowerMAD & x1 <= m + STD* upperMAD])
        xt1 = c(xt, w2)
        xt1[is.na(xt1)] = 0
        return(xt)}
}

#'
#' scQuant.single is scQuant for a single dataset without celltypes, All functions are similar to scQuant. The matrix variable is to return the results as a matrix or as a method list. Default is set to NULL
#' @param geneMat Expression matrix where columns are cells and rows are genes
#' @param zeroPerc Percentage of lowly expressed cells that need to be removed
#' Default set to TRUE
#' @param QuantNo Total number of quantiles the gene expression distribution is to be spli into
#' @param method Different analysis building upon the quantile methodology:
#' "std", returns the quantile NMAD values for each gene, for each cell type and a list of data frames.
#' "SG" calculates genes with and ranks genes with based on how stable the gene expressions are, genes are ranked from most to least stable, where there higher the ranking the more stable the genes are.
#' "DS" calculates genes that are differentially stably expressed between cell types, using two-sample Wilcoxon test.
#' "NoPval", calculates differentially stable genes by comparing the ranking between genes expression from different cell types.
#' "SV" calculate and returns a ranking of genes that are stably variable either between cell types or between datasets
#' @param Matrix whether the results needs to be returned as a matrix or as a list. Default set to list
scQuant.single = function(geneMat = geneMat,
                          zeroPerc = zeroPerc, QuantNo = QuantNo,
                          method = c("SG", "SV", "DS", "std"),
                          Matrix = NULL){

    geneMat_list = list(geneMat, geneMat)

    names(geneMat_list) = c('dat1', 'dat2')

    message("Reorganising data")
    geneMat_list = pblapply(geneMat_list, as_matrix)

    message("Removing genes with excess zeros")

    zeRM = pblapply(geneMat_list, function(x){
        RMze = x[rowMeans(x == 0) < zeroPerc ,]})

    message("Removing outliers")
    OutRM = pblapply(zeRM, function(x){
        apply(x, 1, function(w){
            wt = RMOutlier(w)
            wt = as.list(wt)
            return(wt)})
    })



    OutRM.unlist = lapply(OutRM, lapply, unlist)

    message("Running scQuant")
    SetOut = pblapply(OutRM.unlist, lapply, function(x){
        dt = as.data.frame(table(x))
        if(nrow(dt) == 1){
            t24 = as.data.frame(t(c(0,0,0,0)))
            colnames(t24) = c(seq(1:ncol(t24)))
            return(t24)
        }else{
            t22 = split(x, cut(x, unique(quantile(x, (0:QuantNo/QuantNo),
                                                  names = FALSE),
                                         include = TRUE)))
            t21 = as.data.frame(t(unlist(lapply(t22, NMAD))))
            colnames(t21) = c(seq(1:ncol(t21)))
            return(t21)
        }
    })

    Quant_out.tb = lapply(SetOut, function(x){
        list_data = as.data.frame(data.table::rbindlist(x, fill = T))
        list_data[is.na(list_data)] = 0
        list_data$genes = names(x)
        return(list_data)
    })

    if (method == 'std') {
        return(Quant_out.tb)
    }

    if (method == "SG"){
        SG = scSG (dat.list = Quant_out.tb, QuantNo = QuantNo)
        return(SG)

        if (method == "DS"){
            DS = scDS(dat.list = Quant_out.tb)
            return(DS)
        }

    }
    if (method == "SV"){
        Quant_out.tb = Quant_out.tb[[1]]
        gen.name = Quant_out.tb$genes
        Quant_out.tb$genes = NULL

        datL <- lapply(seq_len(nrow(Quant_out.tb)),
                       function(i) lapply(Quant_out.tb, "[", i))
        names(datL) = gen.name

        f1 <- function(x, y) philentropy::distance(rbind(
            x, y), method = "euclidean", mute.message = TRUE)

        message("Performing pairwise comparisons")
        out <- outer(datL, datL, FUN = Vectorize(f1))

        out[lower.tri(out)] = NA

        if(isTRUE(Matrix)){
            return(out)
        }else{
            df1 = reshape2::melt(out)
            df1 = df1[!is.na(df1$value),]

            df1$Gene_pair = paste0(df1$Var1, "_", df1$Var2)

            df1$rank = rank(df1$value)
            return(df1)
        }
    }
}


#' For converting a sparse matrix into a dense matrix. Use this when after obtaining log normalised count data from either singlecelldata objects or Seurat objects
#' @param mat a sparse matrix that is to be converted to a dense matrix
#
as_matrix = function(mat){

    tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])

    row_pos <- mat@i+1
    col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
    val <- mat@x

    for (i in seq_along(val)){
        tmp[row_pos[i],col_pos[i]] <- val[i]
    }

    row.names(tmp) <- mat@Dimnames[[1]]
    colnames(tmp) <- mat@Dimnames[[2]]
    return(tmp)
}








