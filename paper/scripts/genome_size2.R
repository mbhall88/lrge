##USAGE
##cd /data/scratch/projects/punim2009/genome_size/

##Rscript genome_size1.R  fq_O=AMtb_1.5000.fq paf=AMtb_1.5000.paf gsize=AMtb_1.gsize fq_L=AMtb_1.5000.fq olen=3000,5000 max_l_len=5000 sep=false

library(ShortRead)
library(jsonlite)

args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  args = fromJSON('["fq_O=AMtb_1.5000.fq",
                  "paf=AMtb_1.5000.paf",
                  "gsize=AMtb_1.gsize",
                      "fq_L=AMtb_1.5000.fq",
                  "olen=3000,5000",
                  "max_l_len=5000", 
                  "sep=false"
                  ] ')
}
##just setting option values
df=data.frame(lapply(args, function(x) v=strsplit(x,"=")[[1]]))
opts = as.list(df[2,])
names(opts) = df[1,]
opts$olen=as.numeric(strsplit(opts$olen,",")[[1]])
opts$max_l_len =as.numeric(opts$max_l_len)
opts$sep = opts$sep=="true"|| opts$sep=="TRUE"

overlap_thresh=100

##genomesize for comparison
gsize = if(file.exists(opts$gsize)) scan(opts$gsize) else NA


###read PAF
paf1 = (read.delim(opts$paf, head=F))
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

overlap_reads = .getLens( fq = readFastq(opts$fq_O))

if(opts$fq_L==opts$fq_O){
  long_reads = overlap_reads[1:opts$max_l_len]
}else{
  long_reads = .getLens( fq = readFastq(opts$fq_L), opts$max_l_len)
}
if(opts$sep){
  #MAKE OVERLAP READS DISTINCT FROM LONG READS
  overlap_reads = overlap_reads[!(names(overlap_reads) %in% names(long_reads))]
}

###this loop is over number of short reads

print("num_long_reads num_short_reads shared_reads est_size_mb rel_est_size_percent")

for(num_short_read in opts$olen){
      overlap_reads1 = overlap_reads[1:num_short_read]

      ##FIGURE OUT WHICH READS WE NEED TO SUBTRACT THE LENGTH
      subtract = !is.na(match( names(long_reads),names(overlap_reads1)))
      names(subtract) = names(long_reads)
      shared_reads = length(which(subtract))
      
      ### SUBSET THE PAF TO THOSE WHICH INCLUDE OVERLAP READS
      paf2 = paf1[paf1$qname %in% names(overlap_reads1) | paf1$tname %in% names(overlap_reads1),c(1,6),drop=F]
      
      tot_read_length = sum(overlap_reads1)
      num_reads = length(overlap_reads1)
      
      
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
        print(paste0(k," ",length(overlap_reads1)," ",shared_reads," ",result/1e6," " ,round(100*result/gsize,2)))
      }
}

