#!/usr/bin/Rscript

# Copyright (C) 2016-2018 University of Southern California
#                         Andrew D Smith
# Author: Andrew D. Smith, Xiaojing Ji
# 
#  This is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
#   This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
# along with this software; if not, write to the Free Software
# Foundation, Incpg., 51 Franklin St, Fifth Floor, Boston, MA
# 02110-1301 USA

reaction <- function(dnas, Nt, ncpg) {
  for (k in 1:dim(dnas)[1]) {
    dna <- dnas[k, ]
    for (i in 1:Nt*ncpg) {
      idx <- sample(1:length(params[, "rate"]), 1, prob=params[, "rate"])
      idx.cpg <- sample(1:length(dna), 2)
      target <- idx.cpg[1]
      mediator <- idx.cpg[2]
      if(is.na(params[idx, "med"]) &&  dna[target] == params[idx, "from"]) {
        dna[target] = params[idx, "to"]
      } else if (!is.na(params[idx, "med"])) {
        if (dna[target] == params[idx, "from"] && dna[mediator] == params[idx, "med"]) {
          dna[target] = params[idx, "to"]
        }
      }
    }
    dnas[k, ] <- dna
  }
  return(dnas)
}


replicate <- function(dnas) {
    dnas.new <- dnas
    hemi.idx <- dnas.new[dnas.new[,]==0.5]
    dnas.new[hemi.idx] <- rbinom(n=length(hemi.idx), size=1, prob=0.5)*0.5
    dnas.new[dnas.new[,]==1] <- 0.5
  return(dnas.new)
}

#############################################
## Below are simulation parameters for testing

# Number of CpGs
ncpg <- 80

# Number of average reactions per CpG per generation
Nt <- 100

# Number of generations
Ng <- 200


# Initial proportion of hemimethylated CpGs
p0.hemi <- 0.15

# Initial proportion of unmethylated CpGs
p0.u <- 0.07 


# Reaction parameters
params <- matrix(0, nrow=7, ncol=4)
rownames(params) <- c("uh", "hu", "hm", "mh", "uhm", "hmh", "hmm")
colnames(params) <- c("rate", "from", "to", "med")
params[, "from"] <- c(0, 0.5, 0.5, 1, 0, 0.5, 0.5)
params[, "to"] <- c(0.5, 0, 1, 0.5, 0.5, 1, 1)
params[, "med"] <- c(NA, NA, NA, NA, 1, 0.5, 1)
params[, "rate"] <- c(0.01, 0.005, 0.01, 0.01, 0.2, 0.3, 0.2)
params[, "rate"] <- params[, "rate"] / sum(params[, "rate"])


cpgs <- rep(1, ncpg)
cpgs[1:round(ncpg*p0.hemi)] <- 0.5
cpgs[round(ncpg*p0.hemi+1):round(ncpg*p0.hemi+ncpg*p0.u+1)] <- 0
cpgs <- sample(cpgs)
dnas <- matrix(cpgs, nrow=1, ncol=ncpg)


meth <- rep(0, Ng+1)
meth[1] <- mean(dnas)
#############################################

for (i in 1:Ng) {
  dnas <- replicate(dnas)
  dnas <- reaction(dnas, Nt, ncpg)
  meth[i+1] <- mean(dnas)
}

plot(0:Ng, meth, xlab="Generation", ylab="Methylation", type="l", ylim=c(0, 1))

cpgs <- rep(0, ncpg)
cpgs[1:round(ncpg*p0.hemi)] <- 0.5
cpgs[round(ncpg*p0.hemi+1):round(ncpg*p0.hemi+ncpg*p0.u+1)] <- 1
cpgs <- sample(cpgs)
dnas <- matrix(cpgs, nrow=1, ncol=ncpg)


meth <- rep(0, Ng+1)
meth[1] <- mean(dnas)

for (i in 1:Ng) {
  dnas <- replicate(dnas)
  dnas <- reaction(dnas, Nt, ncpg)
  meth[i+1] <- mean(dnas)
}

lines(0:Ng, meth, xlab="Generation", ylab="Methylation", type="l", ylim=c(0, 1),
      col="blue")