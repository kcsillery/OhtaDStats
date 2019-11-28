## Katalin Csillery, kati.csillery@gmail.com
## modified after OhtaDStats, https://github.com/pfpetrowski/OhtaDStats
## see original code below
## the modified code is wrapped in a new set of functions: dstat.ph, T.phase.known, T.phase.unknown

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


## original version by Paul Petrowski
#' Tomoka Ohta's D Statistics
#' 
#' Implements Ohta's D statistics for a pair of loci. Statistics are returned in a vector in the following order:
#' Number of populations, D2it, D2is, D2st, D'2st, D'2is.
#' 
#' @param index A two-element vector of column names or numbers for which
#' Ohta's D Statistics will be computed.
#' @param data_set Matrix containing genotype data with individuals as rows and
#' loci as columns. Genotypes should be coded as 0 (homozygous), 1 (heterozygous),
#' or 2 (homozygous). Rownames must be subpopulation names and column names
#' should be marker names.
#' @param tot_maf Minimum minor allele frequency across the total population
#' for a marker to be included in the analysis.
#' @param pop_maf Minimum minor allele frequency across a subpopulation for
#' that subpopulation to be included in analysis.
#' 
#' @return nPops Number of subpopulations used for computation, after filtering.
#' @return D2it A measure of the correlation of alleles at two loci on the same
#' gametes in a subpopulation relative to their expectation according to allele
#' frequencies in the total population.
#' @return D2is Expected variance of LD for subpopulations.
#' @return D2st Expected correlation of alleles in a subpopulation relative to their
#' expected correlation in the total population.
#' @return Dp2st Variance of LD for the total population computed over alleles only.
#' @return Dp2is Correlation of alleles at two loci on the same gamete in subpopulations
#' relative to their expected correlation in the total population.
#'
#' @details When the loci being evaluated fail to pass the filtering thresholds determined by
#' tot_maf and pop_maf, NAs are returned.
#'
#' @references
#' Beissinger et al. (2016) Heredity. (https://www.nature.com/articles/hdy201581) &
#' Ohta. (1982) Proc. Natl. Acad. Science. (http://www.pnas.org/content/79/6/1940)
#' 
#' @examples
#' data(beissinger_data)
#' dstat(index = c(5,6), data_set = beissinger_data)
#' 
#' @export

dstat <- function(index, data_set, tot_maf = 0.1, pop_maf = 0.05){
    tot_maf = tot_maf * 2
    pop_maf = pop_maf * 2
    tot_max_thresh = 2 - tot_maf
    pop_max_thresh = 2 - pop_maf 
    geno <- data_set[,c(index[1],index[2])]
    if(mean(geno[,1],na.rm=T)>=tot_maf & mean(geno[,2],na.rm=T)>=tot_maf & mean(geno[,1],na.rm=T)<=tot_max_thresh & mean(geno[,2],na.rm=T)<=tot_max_thresh){
        freqs1 <- unlist(by(geno[,1],rownames(geno),mean,na.rm=T))  # unlist/by is essentially the same thing as tapply
        freqs2 <- unlist(by(geno[,2],rownames(geno),mean,na.rm=T))
        rm <- c(names(freqs1)[which(freqs1<pop_maf | freqs1>pop_max_thresh)], names(freqs2)[which(freqs2<pop_maf | freqs2>pop_max_thresh)])
        if(length(rm)>0) geno <- geno[-which(rownames(geno)%in%rm),]
        if(nrow(geno)>0){
            # Count number of pops
            nPops <- length(table(rownames(geno)))
        
            T <- function(geno){ ### T is a function to compute Tij.s, for any specific subpopulation.
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
                T00 <- (2 * s.00 + s.01 + s.10 + 0.5 * s.11)/length
                T02 <- (2 * s.02 + s.01 + s.12 + 0.5 * s.11)/length
                T20 <- (2 * s.20 + s.10 + s.21 + 0.5 * s.11)/length
                T22 <- (2 * s.22 + s.21 + s.12 + 0.5 * s.11)/length
                return(c(T00,T02,T20,T22))
            }
          
        
            P <- function(geno,out){
                geno <- geno[which(is.na(geno[,1])==F & is.na(geno[,2])==F),]
                length <- nrow(geno)
                p1.0 <- 1-sum(geno[,1],na.rm=T)/(2*length)
                p1.2 <- 1-p1.0
                p2.0 <- 1-sum(geno[,2],na.rm=T)/(2*length)
                p2.2 <- 1-p2.0
                if(out=="freq")  return(c(p1.0, p1.2, p2.0, p2.2))
                if(out=="manip1") return(c(2*p1.0*p2.0, 2*p1.0*p2.2, 2*p1.2*p2.0, 2*p1.2*p2.2))
                if(out=="manip2") return(c(p1.0*p2.0, p1.0*p2.2, p1.2*p2.0, p1.2*p2.2))
            } #P returns allele frequencies or manipulated frequencies
        
            ### Compute Tijs
            Tijs <- by(geno,rownames(geno),T)
            Tijs.n <- as.numeric(unlist(Tijs))
            
            ### Compute P manipulated1
            P1m <-  by(geno,rownames(geno),P,"manip1")
            P1m.n <- as.numeric(unlist(P1m))
            
            ### Compute P manipulated2
            P2m <-  by(geno,rownames(geno),P,"manip2")
            P2m.n <- as.numeric(unlist(P2m))
            
            ### Compute P freq
            Pf <- by(geno,rownames(geno),P,"freq")
            Pf.n <- as.numeric(unlist(Pf))
            
            ### Compute n.pops
            n.pops <- length(levels(as.factor(rownames(geno))))
            
            ### Compute D2is
            D2is <- sum((Tijs.n-P1m.n)^2)/n.pops
            
            ### Compute mean allele frequencies
            p1.0 <- mean(Pf.n[seq(1,length(Pf.n),4)]) # Compute mean frequencies
            p1.2 <- mean(Pf.n[seq(2,length(Pf.n),4)]) #
            p2.0 <- mean(Pf.n[seq(3,length(Pf.n),4)]) #
            p2.2 <- mean(Pf.n[seq(4,length(Pf.n),4)]) # done
            
            ### Create vector of mean allele frequencies (that is comparable to Tijs.n)
            Pf.v <- c()
            Pf.v[seq(1,length(Pf.n),4)] <- 2*p1.0*p2.0
            Pf.v[seq(2,length(Pf.n),4)] <- 2*p1.0*p2.2
            Pf.v[seq(3,length(Pf.n),4)] <- 2*p1.2*p2.0
            Pf.v[seq(4,length(Pf.n),4)] <- 2*p1.2*p2.2
            
            ### Compute mean T
            T00 <- mean(Tijs.n[seq(1,length(Tijs.n),4)]) # Compute mean T
            T02 <- mean(Tijs.n[seq(2,length(Tijs.n),4)]) #
            T20 <- mean(Tijs.n[seq(3,length(Tijs.n),4)]) #
            T22 <- mean(Tijs.n[seq(4,length(Tijs.n),4)]) # done
            
            ### Create vector of mean T (that is comparable to Tijs.n)
            Tijs.n2 <- c()
            Tijs.n2[seq(1,length(Tijs.n),4)] <- T00
            Tijs.n2[seq(2,length(Tijs.n),4)] <- T02
            Tijs.n2[seq(3,length(Tijs.n),4)] <- T20
            Tijs.n2[seq(4,length(Tijs.n),4)] <- T22
            
            ### Compute D2it
            D2it <- sum((Tijs.n-Pf.v)^2)/n.pops
            
            ### Compute D2st
            D2st <- sum((P2m.n-Pf.v/2)^2)/n.pops
            
            ### Compute Dp2st
            p0p0 <- 2*p1.0*p2.0
            p0p2 <- 2*p1.0*p2.2
            p2p0 <- 2*p1.2*p2.0
            p2p2 <- 2*p1.2*p2.2
            Dp2st <- (T00-p0p0)^2+(T02-p0p2)^2+(T20-p2p0)^2+(T22-p2p2)^2
            
            ### Compute Dp2is
            Dp2is <- sum((Tijs.n-Tijs.n2)^2)/n.pops
        }
        else{
            D2it <- NA; D2is <- NA; D2st <- NA; Dp2st <- NA; Dp2is <- NA; nPops <- NA
        }
    }
    else{
        D2it <- NA; D2is <- NA; D2st <- NA; Dp2st <- NA; Dp2is <- NA; nPops <- NA
    }
return(round(c(nPops, D2it, D2is, D2st, Dp2st, Dp2is),6))
}



