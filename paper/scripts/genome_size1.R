##USAGE
##Rscript genome_size1.R  AMtb_1.5000.fq AMtb_1.5000.paf AMtb_1.gsize AMtb_1.5000.fq 5000 false

library(ShortRead)
library(jsonlite)

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  args = fromJSON('["AMtb_1.5000.fq",
                  "AMtb_1.5000.paf",
                  "AMtb_1.gsize",
                      "AMtb_1.5000.fq",
                  5000,
                  "false"
                  ] ')
}


fq_file = args[[1]]   ##overlap reads
paf_in = args[[2]]  ##PAF
gsize_file= args[[3]]
fq_file_long =args[[4]] 
num_long_reads =as.numeric(args[[5]])
separate = args[6]=="true"  ## whether to separate the overlap reads from long reads

overlap_thresh=100

##genomesize for comparison
gsize = scan(gsize_file)


###read PAF
paf1 = (read.delim(paf_in, head=F))
names(paf1) = fromJSON('["qname","qlen","qstart","qend","strand","tname","tlen","tstart","tend","nmatch","alen","mapq","tp","cm","s1","dv","rl"]')


##get rid of overlaps and duplicates
self_overl = apply(paf1,1,function(v)v[1]==v[2])
paf1 = paf1[!self_overl,,drop=F]
ch1 = apply(paf1[,c(1,6)],1,function(v) paste(sort(v), collapse="."))
dupls = duplicated(ch1)
paf1 = paf1[!dupls,,drop=F]


##readfastq

##function to get read lengths from fastq
#returns a vector of lengths which names as read ids
.getLens<-function(fq,num=length(fq@sread)){
  lens = unlist(lapply(fq@sread[1:num], function(x) x@length))
  ids = id(fq)[1:num]
  names(lens) = lapply(ids, function(id1) strsplit(as.character(id1)," ")[[1]][1])
  lens
}

overlap_reads = .getLens( fq = readFastq(fq_file))
if(fq_file_long==fq_file){
  long_reads = overlap_reads
}else{
  long_reads = .getLens( fq = readFastq(fq_file_long), num_long_reads)
}
#SORT LONG READS
long_reads = sort(long_reads, decreasing = T)
long_reads = long_reads[1:num_long_reads]
if(separate){
  overlap_reads = overlap_reads[!(names(overlap_reads) %in% names(long_reads))]
}

##FIGURE OUT WHICH READS WE NEED TO SUBTRACT THE LENGTH
subtract = !is.na(match( names(long_reads),names(overlap_reads)))
names(subtract) = names(long_reads)


### SUBSET THE PAF TO THOSE WHICH INCLUDE OVERLAP READS
paf2 = paf1[paf1$qname %in% names(overlap_reads) | paf1$tname %in% names(overlap_reads),c(1,6),drop=F]

tot_read_length = sum(overlap_reads)
num_reads = length(overlap_reads)


##THIS BIT ITERATES THROUGH LONG READS TO CALCULATE ESTIAMATE
gs=unlist(lapply(names(long_reads), function(lr){
  n_overlaps = length(which(paf2$qname==lr | paf2$tname==lr))
  read_length=long_reads[[lr]]
  if(subtract[[lr]]){
    avg_read_length = (tot_read_length-read_length )/(num_reads-1)
    x = (num_reads-1)/n_overlaps
    
  }else{
    avg_read_length = (tot_read_length )/(num_reads)
    x = (num_reads)/n_overlaps
  }
  (x*(read_length+avg_read_length-2*overlap_thresh))
}))

##print results for different numbers of long reads, in steps of 200
for(k in c(seq(100,length(long_reads),100), length(long_reads))){
  result = round(median(gs[1:k]))
  print(paste0(k," ",length(overlap_reads)," ",result/1e6, "mb ",round(100*result/gsize,2),"% "))
}


