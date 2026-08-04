rm(list =ls())
basedir="/data/scratch/projects/punim2009/genome_size"
library(jsonlite); 
#libraries for plots
library(ggridges); library(ggplot2); library(cowplot)
args = commandArgs(trailingOnly=TRUE)
params_file = "params.json"
if(length(args)>0){
  params_file = args[[1]]
}
params=fromJSON("params.json"); options(params)
library(pafr)

 ##get genome size based on all reads 
  getGenomeSize<-function(cum_read_length, num_reads, num_overlaps,overlap_thresh ){
    read_length = cum_read_length/num_reads
    potential_overlaps = (num_reads *(num_reads-1)/2)
    #    print(potential_overlaps)
    x = (potential_overlaps / num_overlaps)
    print(x)
    genome_length = 2*x*(read_length-overlap_thresh)
    return(genome_length)
  }
  ##per read
  
  
#  getGenomeSize(4.4e3,10000,118e3,100)
  
  ## get genome size for individual read
  getGenomeSizePerRead<-function(read_length, cum_read_length, num_reads, num_overlaps, overlap_thresh){
    potential_overlaps = num_reads-1
    x = potential_overlaps/num_overlaps
    avg_read_length = (cum_read_length - read_length)/(num_reads-1)
    genome_length = x*(read_length+avg_read_length-2*overlap_thresh)
    return(genome_length)
  }
  # this is probability two reads overlap
#  p(overlap(r_2,r_1)) = (rl_1 +rl_2)/gs  
#  expected_num_overlaps (r1) = (n-1)*(rl_1  + mean(rl_2))/gs
 # getGenomeSizePerRead(4.4e3,4.4e3*10000,10000, 20, 100)


  ## this gets a data frame with genome size per read  
  .analyse<-function(basedir, genome, thresh_end=NA, len_thresh = 0){
      mtb=grep(genome,dir(basedir),v=T)
      paf = grep("paf",mtb,v=T)
      paf1 = read_paf(paste(basedir,paf,sep="/"))
      gsize = scan(paste(basedir, grep("gsize", mtb,v=T),sep="/"))
      paf2 = subset(paf1, !is.na(tp))
      gc_mat = read.delim(paste(basedir,grep("gc",mtb,v=T),sep="/"), head=F)
      dimnames(gc_mat)[[1]] = gc_mat[,1]
      gc_mat = gc_mat[,-1]
      
      names(gc_mat) = c("readlength","gc")
      reads_incl = dimnames(gc_mat)[[1]][gc_mat$readlength>=len_thresh]
      paf2 =paf2[paf2$qname %in% reads_incl & paf2$tname %in% reads_incl,,drop=F]
      gc_mat=subset(gc_mat, readlength>=len_thresh)
      ids = dimnames(gc_mat)[[1]]
      if(!is.na(thresh_end)){
        
            q_contain =  paf2$qstart <thresh_end & (paf2$qlen - paf2$qend) < thresh_end
            t_contain=paf2$tstart <thresh_end & (paf2$tlen - paf2$tend) < thresh_end
            q_mid =  paf2$qstart >thresh_end & (paf2$qlen - paf2$qend) > thresh_end
            t_mid=paf2$tstart >thresh_end & (paf2$tlen - paf2$tend) > thresh_end
          
            to_remove = q_contain & t_mid | t_contain & q_mid
            print(paste("dropped",length(which(to_remove))))
            paf2 = paf2[!to_remove ,,drop=F]
      }
    
  #    gc = gc_mat$gc
#      gc_z = (gc -mean(gc))/sd(gc)
#      gc2=gc_z^2
      num_reads = nrow(gc_mat)
      overlap_thresh=100
      overlaps=unlist(lapply(ids, function(id){
        inds1 = which(paf2$qname==id) 
        inds2 = which(paf2$tname==id)
        length(inds1)+length(inds2)
       
      }))
      df =  data.frame(cbind(gc_mat,  overlaps))
      cum_read_length = sum(df$readlength)
      num_reads = nrow(df)
      mi = match(c("readlength", "overlaps"), names(df))
      genomsize=apply(df, 1, function(v){
          getGenomeSizePerRead(v[mi[1]], cum_read_length, num_reads, v[mi[2]], overlap_thresh)
      })
      df2 = cbind(df, genomsize)
      total_num_overlaps = nrow(paf2)
      gs1=getGenomeSize(cum_read_length, num_reads, total_num_overlaps,overlap_thresh )
      attr(df2,"overall")=gs1
      attr(df2,"correct")=gsize
    print(paste("overall",genome,gs1/1e6, "correct",gsize/1e6))  
    df2      
}

  ##this gets plots based on data frame from .analyse
  .getPlots<-function(genome,df2, br=c(300,1000,2000,3000), gc_br = c(20,40,60,80,100)){
    correct=attr(df2,"correct")
    if(length(correct)==0)correct=0
    br = c(br, max(df2$readlength))
    rl =factor(unlist( lapply(df2$readlength, function(x) br[which(br>=x)[1]])))
    gcl =factor(unlist( lapply(df2$gc, function(x) gc_br[which(gc_br>=x)[1]])))
    df2 = cbind(df2[,!(names(df2)%in% c("rl","gl"))],rl, gcl)
    df3 = subset(df2,!is.infinite(genomsize))
    
   # ggplot(df3, aes(x=gc,y=readlength,z=genomsize))+geom_contour()
    #+geom_point( size=0.1,alpha = 0.1)
    plots=list(
      "rl_size" =ggplot(df2,aes(x=readlength, y=genomsize,  color=gc))+scale_color_gradient(low="blue", high="green")+geom_hline(yintercept=correct),
      "rl_overlaps"=ggplot(df2,aes(x=readlength, y=overlaps,  color=gc))+scale_color_gradient(low="blue", high="green"),
      "gc_size" =ggplot(df2,aes(color=readlength, y=genomsize, x=gc))+scale_color_gradient(low="blue", high="green")+geom_hline(yintercept=correct),
      "gc_size2" =ggplot(df2,aes(color=rl, y=genomsize, x=gc))+geom_hline(yintercept=correct),
      "histlen"=ggplot(df3, aes(x=genomsize, y=rl)) + geom_density_ridges() +geom_vline(xintercept=correct),
      "histgc"=ggplot(df3, aes(x=genomsize, y=gcl)) + geom_density_ridges() +geom_vline(xintercept=correct)
      
    )
    #geom_density()
    plots1=lapply(plots, function(p) p<-p+geom_point(size=1)+theme_minimal()+ggtitle(genome))
    combined=plot_grid(plots1[[1]], plots1[[2]], plots1[[3]], plots1[[4]], plots1[[5]], plots1[[6]])
  #  df3 = subset(df2, !is.infinite(genomsize))
  #  a=lm(df3$genomsize~df3$gc+df3$readlength+df3$gc2) #+df1$gc2)
  #  coeff= summary(a)
    combined
  }
  
## gets several estimates of genome size from dataframe
.calcGenomeSize<-function(df2){
  h=hist(df2$genomsize, br=100)
  mode=h$mids[ which.max(h$counts)]
  correct=attr(df2,"correct")
  if(length(correct)==0)correct=NA
  result = list(mode=mode,median=median(df2$genomsize),  overall=attr(df2,"overall"),
               correct=correct)
  o = order(df2$genomsize)
  df3 = subset(df2,overlaps>0)
  #lm1 = lm(overlaps~readlength,data=df3)
  #coeff = summary(lm1)$coeff
  cums=cumsum(df2$readlength[o] )
  midp=which(cums>=max(cums)/2)[1]
  result[["midp"]]=df2$genomsize[o][midp]
  #df3$readlength[o]
  
  lens = c(2000,3000,4000)
  names(lens)=paste("median",lens,sep="_")
  res2=lapply(lens, function(l){
    median(df2$genomsize[df2$readlength>l])
  })
  result2 = c(result ,res2)  
  unlist( result2)
}


### this is where main program starts

genomes = sub(".5000.paf","",grep("paf",dir(basedir),v=T))
names(genomes)=genomes  
allres = lapply(genomes, function(genome) .analyse(basedir ,genome,thresh_end=NA, len_thresh=0))
plots_all =lapply(names(allres), function(n) .getPlots(n, allres[[n]], br=c(0,200,1000,2000,3000)))
names(plots_all) = names(allres)

df4=data.frame(t(data.frame(lapply(allres, function(df2){
  .calcGenomeSize(df2)
}))))

nme=dimnames(df4)[[1]]
df5=cbind(df4,nme)

sze=3
plots2 = list(
ggplot(df5, aes(x=correct, y=median))+geom_point()+geom_text(aes(x=correct, y=median, label=nme),size=sze),
ggplot(df5, aes(x=correct, y=overall))+geom_point()+geom_text(aes(x=correct, y=overall, label=nme),size=sze),
ggplot(df5, aes(x=correct, y=midp))+geom_point()+geom_text(aes(x=correct, y=midp, label=nme),size=sze),
ggplot(df5, aes(x=correct, y=median_2000))+geom_point()+geom_text(aes(x=correct, y=median_2000, label=nme),size=sze),
ggplot(df5, aes(x=correct, y=median_3000))+geom_point()+geom_text(aes(x=correct, y=median_3000, label=nme),size=sze),
ggplot(df5, aes(x=correct, y=median_4000))+geom_point()+geom_text(aes(x=correct, y=median_4000, label=nme),size=sze)
)
plots2 = lapply(plots2, function(ggp)ggp+geom_abline(intercept = 0, slope = 1))
combined=plot_grid(plots2[[1]], plots2[[2]], plots2[[3]], plots2[[4]], plots2[[5]], plots2[[6]])

pdf("out4.pdf", width=20, height=10)
print(combined)
for(k in 1:length(plots_all)){
  print(plots_all[[k]])
}
dev.off()


df6 = subset(df5, !is.infinite(median))

l = lm(df5$correct~df5$median_2000)
summary(l)
