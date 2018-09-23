options(warn=-1)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop("Wrong parameters. Run as: Rscript EntropyFilter.R <max row % ?> <max col % ?> <min entropy> <input file> <output file> <species list>", call.=FALSE)
}

max_rq <- as.numeric(args[1])
max_cq <- as.numeric(args[2])
min_ent <- as.numeric(args[3])
infile <- args[4]
outfile <- args[5]
species_file <- args[6]

logab <- function(a, b) {
  return(log(b)/log(a))
}

entropia <- function(vals, max_cq) {
  tbl <- NULL
  tbl <- as.matrix(table(vals))
  q <- 0
  if ('?' %in% rownames(tbl)) {
    q <- tbl['?',1]
    tbl['?',1] <- 0
  }
  tbl <- as.matrix(tbl[tbl > 0,])
  len_tbl <- length(tbl)
  n <- sum(tbl)
  ent <- 0
  if (len_tbl > 0) {
    for (i in 1:len_tbl) {
      nn <- as.integer(tbl[i,1])
      p <- nn/n
      h <- (-1*p*logab(len_tbl,p))
      ent <- ent + h
    }
  }
  if ((n == 0) | (q/n > max_cq)) {ent <- 0}
  return(ent)
}

numq <- function(val) {
  n <- length(val)
  nq <- table(val=="?")[2]
  return(nq/n)
}

data <- read.table(infile,sep="\t",row.names=1, header=F)

# remove species with too many ?
print("Removing species with too many undefined characters...")
nqs <- apply(data,1,numq)
data_filtered <- data[which(nqs<=max_rq),]

# get species in list
print("Retrieving submatrix with specified species...")
species <- as.matrix(read.table(species_file,sep="\t",header=F))
data_filtered2 <- data_filtered[rownames(data_filtered) %in% species,]

# do entropy filter
print("Performing entropy filter...")
cut <- min_ent
H <- apply(data_filtered2,2,entropia,max_cq=max_cq)
H[is.na(H)] <- 0

new_data <- data_filtered[,which(H > cut)]

write.table(new_data,outfile,sep="\t",quote=F, col.names=F)

print("Filtered data set is ready!...")