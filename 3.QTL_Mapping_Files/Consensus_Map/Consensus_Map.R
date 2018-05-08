library(tidyverse)
library(LPmerge)
##I editted the code to the LPmerge function.
##The only edits I did was had LPmerge function report back the statistics for each interval so I could chose which interval to use for the chromosome
LPmerge <- function (Maps, max.interval = 1:3, weights = NULL) 
{
  n.maps <- length(Maps)
  if (n.maps < 2) {
    print("Error.  Must have at least two maps.")
    stop()
  }
  if (is.null(weights)) {
    weights <- rep(1, n.maps)
  }
  stopifnot(length(weights) == n.maps)
  map.names <- attributes(Maps)$names
  if (is.null(map.names)) {
    map.names <- 1:n.maps
  }
  num.mark <- rep(NA, n.maps)
  num.unique.mark <- rep(NA, n.maps)
  for (i in 1:n.maps) {
    map <- Maps[[i]]
    map <- map[order(map[, 2]), ]
    map[, 2] <- round(map[, 2], 2)
    map[, 1] <- as.character(map[, 1])
    num.unique.mark[i] <- length(unique(map[, 1]))
    num.mark[i] <- nrow(map)
    Maps[[i]] <- map
  }
  errs <- which(num.unique.mark != num.mark)
  if (length(errs) > 0) {
    print("Error. Redundant markers present in following maps:")
    print(paste(map.names[errs], collapse = " "))
    stop()
  }
  map <- Maps[[1]]
  markers <- map[, 1]
  bins <- unique(map[, 2])
  mark.bins <- match(map[, 2], bins)
  for (i in 2:n.maps) {
    map <- Maps[[i]]
    mark.i <- map[, 1]
    bins.i <- unique(map[, 2])
    mark.bins.i <- match(map[, 2], bins.i)
    mo <- length(markers)
    m <- length(union(mark.i, markers))
    new.mark.bins <- rep(0, m)
    shared <- which(is.element(markers, mark.i))
    not.shared <- setdiff(1:mo, shared)
    shared.i <- match(markers[shared], mark.i)
    not.shared.i <- setdiff(1:length(mark.bins.i), shared.i)
    if (length(shared) > 0) {
      new.mark.bins[shared] <- paste(mark.bins[shared], 
                                     mark.bins.i[shared.i], sep = ".")
    }
    if (length(not.shared) > 0) {
      new.mark.bins[not.shared] <- paste(mark.bins[not.shared], 
                                         0, sep = ".")
    }
    if (length(not.shared.i) > 0) {
      if (m > mo) {
        new.mark.bins[(mo + 1):m] <- paste(0, mark.bins.i[not.shared.i], 
                                           sep = ".")
      }
      markers <- c(markers, mark.i[not.shared.i])
    }
    mark.bins <- new.mark.bins
  }
  n.mark <- length(markers)
  bins <- unique(mark.bins)
  n.bin <- length(bins)
  bin.list <- list()
  for (j in 1:n.bin) {
    bin.list[[j]] <- markers[which(mark.bins == bins[j])]
  }
  constraints <- numeric(0)
  for (i in 1:n.maps) {
    map <- Maps[[i]]
    map <- data.frame(bin = match(mark.bins[match(map[, 
                                                      1], markers)], bins), pos = map[, 2])
    uniq <- unique(map$bin)
    map <- map[match(uniq, map$bin), ]
    m <- nrow(map)
    j <- which(map[1, 2] == map[, 2])
    while (max(j) < m) {
      k <- which(map[max(j) + 1, 2] == map[, 2])
      d <- map[k[1], 2] - map[j[1], 2]
      z <- expand.grid(j, k)
      constraints <- rbind(constraints, t(apply(z, 1, 
                                                function(x) {
                                                  return(c(i, map[x, 1]))
                                                })))
      j <- k
    }
  }
  n.constraint <- nrow(constraints)
  print(paste("# markers:", n.mark))
  print(paste("# bins:", n.bin))
  print(paste("# constraints:", n.constraint))
  A <- Matrix(0, nrow = 0, ncol = n.bin, sparse = TRUE)
  for (i in 1:n.constraint) {
    v <- rep(0, n.bin)
    v[constraints[i, 2]] <- -1
    v[constraints[i, 3]] <- 1
    A <- rBind(A, v)
  }
  maxFS <- function(A) {
    print("Finding maximum feasible subsystem.")
    n.constraint <- nrow(A)
    n.mark <- ncol(A)
    B <- cBind(A, Diagonal(n.constraint))
    f <- c(rep(0, n.mark), rep(1, ncol(B) - n.mark))
    ans <- Rglpk_solve_LP(f, B, rep(">=", nrow(B)), rep(1, 
                                                        nrow(B)))
    if (ans$status != 0) {
      stop("Error in LP solver.")
    }
    else {
      if (ans$optimum < 1e-06) {
        return(integer(0))
      }
      elastic.var <- ans$solution[n.mark + 1:nrow(B)]
      sorted <- sort(elastic.var, decreasing = TRUE, index.return = TRUE)
      HoldSet <- sorted$ix[which(sorted$x > 1e-04)]
      eliminate <- integer(0)
      if (length(HoldSet) == 1) {
        return(HoldSet)
      }
      else {
        repeat {
          candidates <- HoldSet
          min.SINF <- Inf
          for (i in 1:length(candidates)) {
            B <- cBind(A[-c(eliminate, candidates[i]), 
                         ], Diagonal(n.constraint - length(eliminate) - 
                                       1))
            constraint.id <- (1:n.constraint)[-c(eliminate, 
                                                 candidates[i])]
            f <- c(rep(0, n.mark), rep(1, ncol(B) - 
                                         n.mark))
            ans <- Rglpk_solve_LP(f, B, rep(">=", nrow(B)), 
                                  rep(1, nrow(B)))
            if (ans$status != 0) {
              stop("Error in LP solver.")
            }
            else {
              if (ans$optimum < 1e-06) {
                eliminate <- c(eliminate, candidates[i])
                return(eliminate)
              }
              else {
                if (ans$optimum < min.SINF) {
                  winner <- candidates[i]
                  min.SINF <- ans$optimum
                  elastic.var <- ans$solution[n.mark + 
                                                1:nrow(B)]
                  sorted <- sort(elastic.var, decreasing = TRUE, 
                                 index.return = TRUE)
                  HoldSet <- constraint.id[sorted$ix[which(sorted$x > 
                                                             1e-04)]]
                  if (length(HoldSet) == 1) {
                    next.winner <- HoldSet
                  }
                  else {
                    next.winner <- NULL
                  }
                }
              }
            }
          }
          eliminate <- c(eliminate, winner)
          if (!is.null(next.winner)) {
            return(c(eliminate, next.winner))
          }
        }
      }
    }
  }
  eliminate <- maxFS(A)
  n.bad <- length(eliminate)
  if (n.bad > 0) {
    print("Eliminated following constraints to resolve marker order conflicts: ")
    for (i in 1:n.bad) {
      mark1 <- paste(bin.list[[constraints[eliminate[i], 
                                           2]]], collapse = " ")
      mark2 <- paste(bin.list[[constraints[eliminate[i], 
                                           3]]], collapse = " ")
      print(paste("Map ", map.names[constraints[eliminate[i], 
                                                1]], ": ", mark1, " < ", mark2, sep = ""))
    }
    A <- A[-eliminate, ]
  }
  else {
    print("Linkage maps had no ordering conflicts.")
  }
  n.composite.maps <- length(max.interval)
  result <- list()
  for (p in seq(2,2 * n.composite.maps,2)) {
    error.terms <- numeric(0)
    n.terms <- rep(0, n.maps)
    for (i in 1:n.maps) {
      map <- Maps[[i]]
      map <- data.frame(bin = match(mark.bins[match(map[, 
                                                        1], markers)], bins), pos = map[, 2])
      uniq <- unique(map$bin)
      map <- map[match(uniq, map$bin), ]
      m <- nrow(map)
      for (q in 1:max.interval[p / 2]) {
        for (j in 1:(m - q)) {
          n.terms[i] <- n.terms[i] + 1
          d <- map[j + q, 2] - map[j, 2]
          error.terms <- rbind(error.terms, c(map[j, 
                                                  1], map[j + q, 1], d, i))
        }
        for (j in (m - q + 1):m) {
          n.terms[i] <- n.terms[i] + 1
          d <- map[j, 2] - map[(j + q)%%m, 2]
          error.terms <- rbind(error.terms, c(map[(j + 
                                                     q)%%m, 1], map[j, 1], d, i))
        }
      }
    }
    n.error.terms <- nrow(error.terms)
    print("Generating consensus map.")
    N <- n.bin + n.error.terms
    G <- cBind(A, Matrix(0, nrow = nrow(A), ncol = n.error.terms, 
                         sparse = TRUE))
    b <- rep(0, nrow(G))
    for (i in 1:n.error.terms) {
      v <- rep(0, N)
      v[error.terms[i, 1]] <- -1
      v[error.terms[i, 2]] <- 1
      v[n.bin + i] <- 1
      b <- c(b, error.terms[i, 3])
      G <- rBind(G, v)
      v <- rep(0, N)
      v[error.terms[i, 1]] <- 1
      v[error.terms[i, 2]] <- -1
      v[n.bin + i] <- 1
      b <- c(b, -error.terms[i, 3])
      G <- rBind(G, v)
    }
    f <- c(rep(0, n.bin), weights[error.terms[, 4]]/n.terms[error.terms[, 
                                                                        4]])
    ans <- Rglpk_solve_LP(f, G, rep(">=", nrow(G)), b)
    if (ans$status != 0) {
      stop("Error in LP solver.")
    }
    else {
      map <- data.frame(marker = markers, position = ans$solution[match(mark.bins, 
                                                                        bins)], stringsAsFactors = F)
      composite.map <- map[order(map$position), ]
      row.names(composite.map) <- NULL
      print(paste("Max.Interval = ", max.interval[p / 2], 
                  sep = ""))
      print(paste("Consensus map length:", max(composite.map$position)))
      link.maps <- numeric(0)
      for (k in 1:n.maps) {
        map <- Maps[[k]]
        ix <- match(composite.map$marker, map[, 1])
        link.maps <- cbind(link.maps, ifelse(is.na(ix), 
                                             NA, map[ix, 2]))
      }
      RMSE <- apply(link.maps, 2, function(x) {
        sqrt(mean((composite.map$position - x)^2, na.rm = TRUE))
      })
      print(data.frame(map = c(map.names, "mean", "sd"), 
                       RMSE = round(c(RMSE, mean(RMSE), sd(RMSE)), 
                                    2)))
      colnames(link.maps) <- map.names
      result[[p - 1]] <- data.frame(marker = composite.map$marker, 
                                consensus = composite.map$position, link.maps, 
                                stringsAsFactors = F)
      result[[p]] <- data.frame(k = rep(p / 2,4),
                                map = c(map.names, "mean", "sd"), 
                                RMSE = round(c(RMSE, mean(RMSE), sd(RMSE)), 
                                             2))
    }
  }
  return(result)
}


##This reads in the genetic maps that are ordered long format.
geneticMap <- read_tsv("Genetic_Maps.txt") %>%
  group_by(Population,Chromosome) %>%
  mutate(Max = max(Position)) %>%
  ungroup()

##This builds the consensus map for each chromosome. It uses the LPmerge function that I editted above to construct it
##We wanted to minimize RMSE and pulled the first ranked one as the map of choice
consensusMap <- NULL
for(i in 1:10) {
  maxK <- 6
  mapList <- list("B73" = geneticMap %>%
                           filter(Population == "B73",
                                  Chromosome == i) %>%
                           select(Marker,Position)%>%
                           as.data.frame(),
                         "Mo17" = geneticMap %>%
                           filter(Population == "Mo17",
                                  Chromosome == i) %>%
                           select(Marker,Position)%>%
                           as.data.frame())
  makeConMap <- LPmerge(mapList,max.interval = 1:maxK)
  mapTally <- NULL
  for (i in 1:maxK){
    mapTally <- bind_rows(mapTally,
                          as.data.frame(makeConMap[i * 2]))
  }
  map_to_select <- mapTally %>% 
    filter(map == "sd") %>%
    arrange(RMSE) %>% 
    mutate(Rnk = rank(RMSE,ties.method = "first")) %>% 
    filter(Rnk == 1) %>% 
    pull(k)
  consensusMap <- bind_rows(consensusMap,
                            as.data.frame(makeConMap[(map_to_select * 2) - 1]))
}

write_tsv(consensusMap,
          "ConsensusMap.txt")

#####Consensus Summary######
##This is used to summarize the genetic maps for each one to report back in the manuscript
##cM per Chr
consensusMap %>% 
  separate(marker,into = c("Chr","BP"), sep  = "_",
           remove = FALSE) %>%
  mutate(Chr = as.numeric(Chr),
         BP = as.numeric(BP)) %>% 
  group_by(Chr) %>% 
  summarise(Length = max(consensus))

##Total genetic distance
consensusMap %>% 
  separate(marker,into = c("Chr","BP"), sep  = "_",
           remove = FALSE) %>%
  mutate(Chr = as.numeric(Chr),
         BP = as.numeric(BP)) %>% 
  group_by(Chr) %>% 
  summarise(Length = max(consensus)) %>% 
  summarise(Total = sum(Length)) %>% 
  pull(Total)

##Genetic distance of B73 Population
##Summarizes the two maps:
geneticMap %>%
  select(Population,Max) %>%
  unique() %>% 
  group_by(Population) %>% 
  summarise(Total = sum(Max)) %>% 
  pull(Total)

##Number of bins used for each population
geneticMap %>% 
  group_by(Population) %>% 
  summarise(n = n())
