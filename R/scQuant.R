#Final functions for scQuant
library(data.table)
library(Hmisc)
library(Seurat)
library(pbapply)
library(dplyr)
library(progress)

scQuant = function(geneMat = geneMat, celltype = celltype, 
                   zeroPerc = zeroPerc, 
                   OutlierRemoval = TRUE,
                   QuantNo =  QuantNo,
                   method = c("SG", "SV", "DS", "std",
                              "NoPval")){
    
    
    message("Removing cell types with fewer than 50 cells")
    
    CT = as.data.frame(table(celltype))
    
    CT = CT[CT$Freq >= 50 ,]
    
    celltype = as.character(celltype)
    celltype = celltype[celltype %in% CT$celltype]
    
    message("Subsetting data by Cell Type")
    geneMat_list = pblapply(names(table(celltype)), function(x) {
        (geneMat[, celltype == x] )
    })
    
   
    names(geneMat_list) = names(table(celltype))
    

    message("removing genes with excessive zeros")
    zeRM = pblapply(geneMat_list, function(x){
        RMze = x[rowMeans(x == 0) < zeroPerc ,]})
    
    
    if (OutlierRemoval == TRUE) {
        message("removing outliers")
        OutRM = pblapply(zeRM, function(x){
            apply(x, 1, function(w){
                wt = RMOutlier(w)
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
    if (method == "DS"){
        DS = scDS(dat.list = Quant_out.tb)
        return(DS)
    }
    if (method == "NoPval"){
        DS = OldscSG(dat.list = Quant_out.tb, QuantNo = QuantNo)
    }
}  
   

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

#For calculating SV genes
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
scDS = function(dat.list = dat.list){
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

OldscSG = function(dat.list, QuantNo){
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
#NMAD
NMAD = function(x) {
    mad(x[x!= 0])/mean(x[x!= 0])}
RMOutlier = function(x){
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
        
        xt = (x1[x1 >= m - 4 * lowerMAD & x1 <= m + 4* upperMAD])
        xt1 = c(xt, w2)
        xt1[is.na(xt1)] = 0
        return(xt)}
}
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

as_matrix <- function(mat){
    
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


      
    



