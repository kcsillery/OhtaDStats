## by Katalin Csillery, kati.csillery@gmail.com
## modified after OhtaDStats, https://github.com/pfpetrowski/OhtaDStats

dstat.ph <- function(index, data, phased, tot_maf = 0, pop_maf = 0){
    tot_maf = tot_maf * 2
    pop_maf = pop_maf * 2
    tot_max_thresh = 2 - tot_maf
    pop_max_thresh = 2 - pop_maf

    ## KATI
    if(phased){
        geno.phased <- data[,c(index[1],index[2])]
        geno <- data[,c(index[1],index[2])]
        ## for the geno, simplify: confound the AB and BA
        geno[geno==2] <- 1
        geno[geno==3] <- 2
        ##print(table(geno.phased, geno))
    }
    else geno <- data[,c(index[1],index[2])]
    
    if(mean(geno[,1],na.rm=T)>=tot_maf & mean(geno[,2],na.rm=T)>=tot_maf & mean(geno[,1],na.rm=T)<=tot_max_thresh & mean(geno[,2],na.rm=T)<=tot_max_thresh){
        freqs1 <- unlist(by(geno[,1],rownames(geno),mean,na.rm=T))  # unlist/by is essentially the same thing as tapply
        freqs2 <- unlist(by(geno[,2],rownames(geno),mean,na.rm=T))
        rm <- c(names(freqs1)[which(freqs1<pop_maf | freqs1>pop_max_thresh)], names(freqs2)[which(freqs2<pop_maf | freqs2>pop_max_thresh)])
        if(length(rm)>0) geno <- geno[-which(rownames(geno)%in%rm),]
        if(nrow(geno)>0){
            
            ## Compute n.pops
            n.pops <- length(unique(rownames(geno)))

            ## PER POP: Compute expected haplo freq 
            exp.haplo.pop <-  by(geno,rownames(geno),P)
            exp.haplo.pop <- as.numeric(unlist(exp.haplo.pop))

            ## PER POP: Compute observed (estimated) haplo freq
            if(phased) obs.haplo.pop <- by(geno.phased, rownames(geno.phased), T.phase.known)
            else obs.haplo.pop <- by(geno,rownames(geno),T.phase.unknown)
            obs.haplo.pop <- as.numeric(unlist(obs.haplo.pop))
                        
            ## MEAN: Compute expected haplo freq 
            exp.haplo.mean <- P(geno)

            ## MEAN: Compute observed (estimated) haplo freq
            if(phased) obs.haplo.mean <- T.phase.known(geno.phased)
            else obs.haplo.mean <- T.phase.unknown(geno)
            obs.haplo.mean <- as.numeric(unlist(obs.haplo.mean))

            ## Compute D statistics
            D2is <- sum((obs.haplo.pop - exp.haplo.pop)^2)/n.pops
            D2it <- sum((obs.haplo.pop - exp.haplo.mean)^2)/n.pops
            D2st <- sum((exp.haplo.pop - exp.haplo.mean)^2)/n.pops
            Dp2st <- sum((obs.haplo.mean - exp.haplo.mean)^2)/n.pops
            Dp2is <- sum((obs.haplo.pop - obs.haplo.mean)^2)/n.pops
        }
        else{
            D2it <- NA; D2is <- NA; D2st <- NA; Dp2st <- NA; Dp2is <- NA; n.pops<- NA
        }
    }
    else{
        D2it <- NA; D2is <- NA; D2st <- NA; Dp2st <- NA; Dp2is <- NA; n.pops<- NA
    }

    return(round(c(n.pops, D2it, D2is, D2st, Dp2st, Dp2is),6))
    
}

## KATI: original T
T.phase.unknown <- function(geno){ ### T is a function to compute Tij.s, for any specific subpopulation.
    geno <- geno[which(is.na(geno[,1])==F & is.na(geno[,2])==F),] #Remove any rows with NA values
    length <- nrow(geno)
    s.00 <- length(which(geno[,1]==0 & geno[,2]==0))  # Getting haplotype counts for the designated loci
    s.01 <- length(which(geno[,1]==0 & geno[,2]==1))
    s.02 <- length(which(geno[,1]==0 & geno[,2]==2))
    s.10 <- length(which(geno[,1]==1 & geno[,2]==0))
    s.11 <- length(which(geno[,1]==1 & geno[,2]==1))
    s.12 <- length(which(geno[,1]==1 & geno[,2]==2))
    s.20 <- length(which(geno[,1]==2 & geno[,2]==0))
    s.21 <- length(which(geno[,1]==2 & geno[,2]==1))
    s.22 <- length(which(geno[,1]==2 & geno[,2]==2))
    T00 <- (2 * s.00 + s.01 + s.10 + 0.5 * s.11)/(2*length)
    T02 <- (2 * s.02 + s.01 + s.12 + 0.5 * s.11)/(2*length)
    T20 <- (2 * s.20 + s.10 + s.21 + 0.5 * s.11)/(2*length)
    T22 <- (2 * s.22 + s.21 + s.12 + 0.5 * s.11)/(2*length)
    return(c(T00,T02,T20,T22))
}

## KATI: modified T, for phased data, ## coded 0 (homAA), 1 (hetAB), 2 (hetBA), 3 (homBB)
T.phase.known <- function(x){ ### T is a function to compute Tij.s, for any specific subpopulation.
    x <- x[which(is.na(x[,1])==F & is.na(x[,2])==F),] #Remove any rows with NA values
    length <- nrow(x)
    s.00 <- length(which(x[,1]==0 & x[,2]==0))  # Getting haplotype counts for the designated loci
    s.01 <- length(which(x[,1]==0 & x[,2]==1))
    s.03 <- length(which(x[,1]==0 & x[,2]==3))
    s.10 <- length(which(x[,1]==1 & x[,2]==0))
    s.11 <- length(which(x[,1]==1 & x[,2]==1))
    s.13 <- length(which(x[,1]==1 & x[,2]==3))
    s.30 <- length(which(x[,1]==3 & x[,2]==0))
    s.31 <- length(which(x[,1]==3 & x[,2]==1))
    s.33 <- length(which(x[,1]==3 & x[,2]==3))
    ## new
    s.20 <- length(which(x[,1]==2 & x[,2]==0))
    s.21 <- length(which(x[,1]==2 & x[,2]==1))
    s.23 <- length(which(x[,1]==2 & x[,2]==3))
    s.02 <- length(which(x[,1]==0 & x[,2]==2))
    s.12 <- length(which(x[,1]==1 & x[,2]==2))
    s.32 <- length(which(x[,1]==3 & x[,2]==2))
    s.22 <- length(which(x[,1]==2 & x[,2]==2))
    
    T00 <- (2 * s.00 + s.01 + s.10 + s.11 + s.20 + s.02 + s.22)/(2*length)
    T02 <- (s.01 + 2 * s.03 + s.13 + s.21 + s.23 + s.02 + s.12)/(2*length)
    T20 <- (s.10 + 2 * s.30 + s.31 + s.20 + s.21 + s.12 + s.32)/(2*length)
    T22 <- (s.11 + s.13 + s.31 + 2 * s.33  + s.23 + s.32 + s.22)/(2*length)
    return(c(T00,T02,T20,T22))
}

P <- function(x){ ## for a given sub population calculate expected haplo frequencies
    x <- x[which(is.na(x[,1])==F & is.na(x[,2])==F),]
    length <- nrow(x)
    ## loc 1
    p1 <- 1-sum(x[,1],na.rm=T)/(2*length) ## 1- because 0 is AA
    q1 <- 1-p1
    ## loc 2
    p2 <- 1-sum(x[,2],na.rm=T)/(2*length) ## 
    q2 <- 1-p2
    return(c(p1*p2, p1*q2, q1*p2, q1*q2))
} 
