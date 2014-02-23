#
source("utilities.R")



n.tax <- 50
tree.age <- 50
old.cal <- 20

sim.het.tree <- function(n.tax = 50, tree.age = 50, old.cal = 10){
    tip.lengths <- old.cal /2
    n <- 0
    while(all(tip.lengths < old.cal)){
        tr <- rcoal(n.tax)
        init.ages <- allnode.times(tr)
        root.age <- max(init.ages)
        scale.factor <- root.age / tree.age
        tr$edge.length <- tr$edge.length / scale.factor
        ages.tips <- init.ages[1:n.tax]
        root.node <- (1:nrow(tr$edge))[which.max(node.depth(tr))]
        ascendants.tips <- list()
        for(i in 1:length(tr$tip.label)){
            tip.temp <- i
            ascendants <- vector()
            br.temp <- tip.temp
            while(tip.temp != root.node){
                br.temp <- which(tr$edge[, 2] == tip.temp)
                tip.temp <- tr$edge[, 1][br.temp]
                ascendants <- c(ascendants, br.temp)
            }
            ascendants.tips[[i]] <- ascendants
        }
        tip.edges <- sapply(ascendants.tips, function(x) return(x[1]))
        tip.lengths <- sapply(tip.edges, function(x) tr$edge.length[x])
#if any tip.length has the desired depth (old.cal)
        n <- n+1
    }
    longest.tip <- (1:length(tr$tip.label))[which.max(tip.lengths)]
    trim.tips <- runif(n=length(tip.edges), min = 0, max = old.cal / tree.age)
    tip.trimmed <- tip.lengths - (tip.lengths*trim.tips)
    tip.trimmed[longest.tip] <- tip.lengths[longest.tip] - old.cal
    tr$edge.length[tip.edges] <- tip.trimmed
    print(paste("Finished with", n, "iterations"))
    return(tr)
}


for(i in 1:20){
tr <- sim.het.tree(n.tax = 20, tree.age = 50, old.cal = 30)
plot(tr, show.tip.label = F)
tiplabels(round(allnode.times(tr, tipsonly = T), 1), cex = 0.7)
system("sleep 2")
}
