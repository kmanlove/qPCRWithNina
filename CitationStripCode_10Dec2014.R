# 1) store a couple functions that will be used to clean things up later
trim.leading <- function (x)  sub("^\\s+", "", x)

# 2) read in data.frame (which is csv of citations as exported from Web of Science, but you could modify that)
data.frame <- read.csv(file.choose())
data.frame$CitedRefs <- data.frame$CR # renames Web of Science export column "CR" in data.frame "CitedRefs"
data.frame$Authors <- data.frame$AU # renames Web of Science export column "AU" to "Authors"
data.frame$PubYear <- data.frame$PY # renames Web of Science export column "PY" to "PubYear"

# 3) build some storage objects to stick the citations in
citation.list <- citation.frame <- citation.frame.small <- vector("list", dim(data.frame)[1])
first.author <- pub.year <- rep(NA, dim(data.frame)[1])

# 4) strip out any papers that don't have references (for me, those were the 1st, 484th, 591st and 843rd entried in my data.frame)
papers.with.cites <- c(1:dim(data.frame)[1])[-c(1, 484, 591, 843)]

# 5) loop over all papers in the papers.with.cites object
for(i in papers.with.cites){
  citation.list[[i]] <- strsplit(x = as.character(data.frame$CitedRefs)[i], split = ";")[[1]]
  citation.frame[[i]] <- matrix(NA, nrow = length(citation.list[[i]]), ncol = 15)
  citation.frame.small[[i]] <-  matrix(NA, nrow = length(citation.list[[i]]), ncol = 12)
  first.author[i] <- strsplit(as.character(data.frame$Authors[i]), split = ";")[[1]][1]
  pub.year[i] <- data.frame$PubYear[i]
  for(j in 1:length(citation.list[[i]])){
    citation.frame.small[[i]][j, ] <- trim.leading(c(strsplit(citation.list[[i]][j], split = ",")[[1]], rep(NA, 12 - length(strsplit(citation.list[[i]][j], split = ",")[[1]]))))
    if(is.na(as.numeric(citation.frame.small[[i]][j, 1])) == F){
      citation.frame.small[[i]][j, ] <- c(NA, citation.frame.small[[i]][j, -12])
    }
    citation.frame[[i]][j, ] <- c(as.character(data.frame$DOI)[i], as.character(first.author[i]), pub.year[i], citation.frame.small[[i]][j, ])
    citation.frame[[i]][j, 1] <- paste("DOI ", citation.frame[[i]][j, 1], sep = "")    
    citation.frame[[i]][j, 2] <- tolower(trim.leading(citation.frame[[i]][j, 2]))
    citation.frame[[i]][j, 4] <- tolower(trim.leading(citation.frame[[i]][j, 4]))
  }
}

full.citation.frame <- do.call("rbind", citation.frame) # builds data frame that contains ALL cited refs
table(full.citation.frame$SO) # tables source journals for all papers cited in original search
