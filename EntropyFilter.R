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
n_species <- dim(data)[1]
n_chars <- dim(data)[2]
q_data <- table(data == "?")
q <- q_data[2]/(q_data[1]+q_data[2])

# calculate pre-filter entropy
H_pre <- apply(data,2,entropia,max_cq=1)
H_pre[is.na(H_pre)] <- 0
H_pre_mean <- sum(H_pre)/length(H_pre)

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
n_species_filtered <- dim(new_data)[1]
n_chars_filtered <- dim(new_data)[2]
q_data_filt <- table(new_data == "?")
q_filt <- q_data_filt[2]/(q_data_filt[1]+q_data_filt[2])

# calculate post-filter entropy
H_post <- apply(new_data,2,entropia,max_cq=1)
H_post[is.na(H_post)] <- 0
H_post_mean <- sum(H_post)/length(H_post)

red <- (n_species * n_chars)/(n_species_filtered * n_chars_filtered)

write.table(new_data,outfile,sep="\t",quote=F, col.names=F)

print("Filtered data set is ready!...")

# write stats file
header <- paste("####### ","Analysis done on ",date()," #######",sep="")
write.table(header,'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('Input file: ',infile,sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('Species file: ',species_file,sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('Output file: ',outfile,sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('Max row undefined %: ',toString(max_rq),sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('Max column undefined %: ',toString(max_cq),sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('Min entropy: ',toString(min_ent),sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('no. species before filter: ',toString(n_species),sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('no. characters before filter: ',toString(n_chars),sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('% undefined before filter: ',toString(q),sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('no. species after filter: ',toString(n_species_filtered),sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('no. characters after filter: ',toString(n_chars_filtered),sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('% undefined after filter: ',toString(q_filt),sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('% reduction in data: ',toString(red),sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('% Mean pre-filter entropy: ',toString(H_pre_mean),sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table(paste('% Mean post-filter entropy: ',toString(H_post_mean),sep=''),'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)
write.table("\n",'stats.txt',append=TRUE,quote=F,col.names=F,row.names=F)