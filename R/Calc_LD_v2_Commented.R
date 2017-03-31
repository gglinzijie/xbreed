#'  Calculate linkage disequilibrium
#'
#' Different measures of linkage disequilibrium (LD) such as \eqn{D},\eqn{r} and \eqn{r^2} are calculated for phased genotypes. LD measurements can be calculated both for adjacent and pairwise loci. Decay of LD between marker pairs can be assessed as well.


################################################
#################PARAMETERS#####################
#################################################

#' @param mat (\code{matrix}) Matrix with phased genotypes of individuals coded as {11,12,21,22} for genotypes {AA, Aa, aA and AA} respectively. Dimension is \eqn{n*p} where n is number of individuals and p is twice the number of loci.
#' @param MAF  \emph{Optional} Minor allele frequency threshold for LD calculation. Loci with minor allele frequency less than \code{MAF} will be excluded. Default: {0.05}.

#' @param method \emph{Optional} (\code{character}) Method to be used for calculation of LD. Possible options are:
	#' \itemize{
	#'\item{"adjacent"}  {LD will be calculated for adjacent loci only}.
	#'\item{"pairwise"}  {LD will be calculated for all pairwise loci}.
	#' }
#'	Default: "adjacent"

#' @param LD_summary \emph{Optional} (\code{Logical}) Display LD calculation summary if is not \code{FALSE}. Default: \code{True}. 

#' @param saveAt \emph{Optional} (\code{character}). Name to be used to save output files.

#' @param linkage_map \emph{Optional} (\code{data.frame}) or (\code{vector}) Linkage map of the loci for the measurement of LD decay. \bold{Note:} both argument \code{linkage_map} and \code{interval} should be present in the function and also argument \code{method} should be set to "pairwise" for the calculation of LD decay.

#' @param interval \emph{Optional} Interval to be used for grouping loci by their pairwise distance. See \code{details}. 

################################################
#################RETURN/KEY/EXPORT##############
################################################

#' @return \code{list} with data of LD calculations.\cr
	#' \describe{
	#'\item{$Mean_r2}{Mean \eqn{r^2} for provided genotypes based on the method specified}. 
	#'\item{$ld_data}{Data frame with 12 columns including pair id, frequencies of alleles and haplotypes for each pair as well as measurements of LD}.  
	#'\item{$ld_decay}{Data for LD decay including the average \eqn{r^2} for loci pairs in each interval}.  
	#' }
	
#' 
#' @export calc_LD

################################################
#################EXAMPLS########################
################################################

#' @examples 
#' # Calculate mean r2 and LD decay.
#'
#'genome<-data.frame(matrix(NA, nrow=1, ncol=6))
#'names(genome)<-c("chr","len","nmrk","mpos","nqtl","qpos")
#'genome$chr<-c(1)
#'genome$len<-c(100)	
#'genome$nmrk<-c(100)
#'genome$mpos<-c("rnd")	
#'genome$nqtl<-c(50)
#'genome$qpos<-c("even")	
#'genome
#'
 #'hp<-make_hp(hpsize=100,
   #' ng=10,h2=0.3,phen_var=1 ,genome=genome,
   #' mutr=2.5e-4)
   #'
#'# Mean r2
#'
 #'mat<-hp$hp_mrk[,-1]
#'rLD<-calc_LD(mat=mat,MAF=0.1,method='adjacent',LD_summary=TRUE)
#' 
#' # LD decay
#'
#'  linkage_map<-hp$linkage_map_mrk[,3]
#'  rLD<-calc_LD(mat=mat,MAF=0.1,method='pairwise'
#' ,LD_summary=TRUE,linkage_map=linkage_map,interval=5)
#'
#'rLD$ld_decay

################################################
#################DETAILS########################
################################################

#' @details	
#' The extent of LD is an important factor both in association studies and genomic selection. Commonly used measure to calculate LD between loci A and B is Pearson coefficient (\eqn{r}) of correlation as: \cr 

#' \deqn{r = D/ \sqrt{(p_1p_2q_1q_2)}}  \cr 
#' where \eqn{D} is \cr 
#' \deqn{D_{ij}=p(A_iB_j)-p(A_i)p(B_j)} \cr 
#' However, squared coefficient of correlation \eqn{r^2} is often used to remove the arbitrary sign introduced: \cr 
 #' \deqn{r^2_{ij} = D^2_{ij}/ (p(A_i)(1-p(A_i))p(B_j)(1-p(B_j)))}
#' \cr
#' To determine the decay of LD with increasing distance between loci (SNPs), the average \eqn{r^2} can be expressed as a function of distance between SNPs. SNP pairs are grouped by their pairwise distance into intervals defined by the user in argument \code{interval}. The average \eqn{r^2} for SNP pairs in each interval are estimated as the mean of all \eqn{r^2} within that interval.



calc_LD<-function(mat,MAF,method,LD_summary,saveAt,linkage_map,interval) {

  # # loading .dll
# dyn.load('ld.dll')
# is.loaded("ld")

# dyn.load('cf.dll')
# is.loaded("cf")


# Note:in function mitavand pairwise ra baraye hameye chromosome yekja hesab konad be sharti ke ld_decay darkhast nashavad. agar ld decay bayad mohasebe shavad dar in sorat faghat yek choromosome info bayesti vared shavad
# ---------------------START-----------------
# INPUT by user control section
# ---------------------START-----------------

# mat
	if(missing(mat)) {
	stop('\n','--- Argument "mat" is missing')
	}
	
	mat<-as.matrix(mat)	
		
	if(dim(mat)[2]%%2!=0) {
	stop('\n','--- Error in number of colomns in argument "mat"')
	}
	
	if(length(which(mat!=1 & mat!=2))>0){
	stop('\n','--- Error in argument "mat". Wrong values in mat')
	}


# MAF
	if(missing(MAF)) {
	cat('MAF is missing, it has been set to default value of 0.05.',fill=TRUE)
	MAF<-0.05
	}
	
	if(!missing(MAF)) {
		if(MAF>0.49999 | MAF<0){
		stop('\n','--- Error in argument "MAF". Possible range: 0<=x<=0.49999')
		}
	}
	
# method
	# method
	if(missing(method)) {
	cat('Method to calculate LD is missing, LD will be calculated for all adjacent loci.',fill=TRUE)
	method<-'adjacent'
	}

	test <- c('adjacent','pairwise')
	if(method%in%test==FALSE){
	stop('\n','--- Possible options in argument "method" are "adjacent" and "pairwise"')
	}
	
	
# LD_summary
	if(missing(LD_summary)) {
	LD_summary<-TRUE
	}
	
	if(!missing(LD_summary)) {
			if(is.logical(LD_summary)==FALSE){
			stop('\n','--- Argument "LD_summary" should be type "logical"')
			}
	}

# saveAt
	if(!missing(saveAt)){
		if(is.character(saveAt)==FALSE){
		stop('Define a name for the saveAt argument as type "character"')
		}
	outFile<-paste(saveAt,'.txt',sep='')
	}
	

	# linkage_map,interval
	
	if(missing(linkage_map) & missing(interval)) {
	Ld_decay<-FALSE	
	}
	
		if(!missing(linkage_map) & missing(interval)){
		stop('\n','--- For calculation of LD decay, arguments "linkage_map" and "interval" should be provided in the function.')
		}
		
		if(missing(linkage_map) & !missing(interval)){
		stop('\n','--- For calculation of LD decay, arguments "linkage_map" and "interval" should be provided in the function.')
		}
		
	
	if(!missing(linkage_map)) {

		if(length(linkage_map)!=(dim(mat)[2]/2)){
		stop('\n','--- Error in linkage map. Linkage map does not match number of loci')
		}
		if(length(linkage_map)!=length(unique(linkage_map))){
		stop('\n','--- Error in linkage map. Some loci have the same position in linkage map')
		}
		if(length(which(linkage_map<0)>0)){
		stop('\n','--- Error in linkage map. Negative value in linkage map')
		}
		
	}
	
	# interval
		if(!missing(linkage_map) & !missing(interval) & method=='pairwise'){
		Ld_decay<-TRUE
		chrom_length<-max(linkage_map)
		intervals<-seq(0,chrom_length,interval)
			if(length(intervals)==1){
			stop('\n','--- Wrong input for argument "interval"')
			}
		}
		
	# test whether 	method is correct for Ld_decay
		if(!missing(linkage_map) & !missing(interval) & method!='pairwise'){
		stop('\n','--- For the calculation of LD decay argument "method" should be set to "pairwise"')
		}
		
	
# ---------------------FINISH-----------------
# INPUT by user control section
# ---------------------FINISH-----------------
# Necessary internal functions---START
		Calc_freq<-function(Loci_Mat){
		freq1<-as.numeric()
		freq2<-as.numeric()
		
		s1<-seq(1,length(Loci_Mat[1,]),2)
		s2<-seq(2,length(Loci_Mat[1,]),2)
		a1<-Loci_Mat[,s1]+Loci_Mat[,s2]
		a1[a1==3]=1
		a1[a1==4]=0
		Loci_Mat2<-a1
		dunkan<-function(vec){
			sum(vec)/(length(vec)*2)
			}
		freq1<-apply(Loci_Mat2,2,dunkan)
		freq2<-1-freq1
		return(freq1)
	    }
		
# Necessary internal functions---FINISH

if(Ld_decay==TRUE){
chrom_length<-max(linkage_map)
intervals<-seq(0,chrom_length,interval)
intervals
}


mat<-as.matrix(mat)		
# freq1<-Calc_freq(mat)

# Passing to .Fortran for calc Freq--start
loci_mat<-mat
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1<-as.numeric(outfreq[[5]])
 # Passing to .Fortran for calc Freq--finish

freq2<-1-freq1
# -----FOR MAF
freqmatrix<-matrix(0,ncol=2,nrow=length(freq1))
freqmatrix[,1]<-freq1
freqmatrix[,2]<-freq2
MAFv<-apply(freqmatrix,1,min)
MAFv1<-MAFv[MAFv>=MAF]
index<-which(MAFv%in%MAFv1)
# -------------------------------

if(Ld_decay==TRUE){
linkage_map1<-linkage_map[index]
}

# extract loci passed MAF
new_index<-index*2
index<-c(new_index-1,new_index)
index<-sort(index)
mat_after_maf<-mat[,index]

# args for .Fortran
r_size<-nrow(mat_after_maf)
c_size<-ncol(mat_after_maf)
outi2<-unlist(mat_after_maf)
outi2<-as.integer(mat_after_maf)

if (method=='pairwise'){
npair1<-combn(1:(ncol(mat_after_maf)/2),2)
npair<-length(npair1[1,])
arg5<-npair*12
argmet<-1 #two options 1:pairwise 2:adjacent
}
if (method=='adjacent'){
npair<-(ncol(mat_after_maf)/2)-1
arg5<-npair*12
argmet<-2 #two options 1:pairwise 2:adjacent
}

LD_Out<-.Fortran("ld",r_size= as.integer(r_size), c_size = as.integer(c_size),loci_mat=outi2,npair=as.integer(npair),method=as.integer(argmet),ld_data=double(arg5))

ld_data<-matrix(as.numeric(LD_Out[[6]]),nrow = npair,ncol = 12)
dim(ld_data)
Mean_r2<-mean(ld_data[,12])


if(Ld_decay==TRUE){
	npair1<-t(npair1)
	mydif<-linkage_map1[npair1[,2]]-linkage_map1[npair1[,1]]
	length(mydif)
	if(length(which(mydif<0))){
	stop('--- Internal error in linkage map')
	}
	ld_data2<-matrix(ncol=2,nrow=length(mydif))
	ld_data2[,1]<-ld_data[,12]
	ld_data2[,2]<-mydif

	r2_A<-matrix(nrow=length(intervals)-1,ncol=2)
	r_breedA<-ld_data2
	r_breedA<-r_breedA[order(r_breedA[,2]),] 

	 for (i in 1:(length(intervals)-1)){
	 insideA<-subset(r_breedA,r_breedA[,2]>intervals[i] & r_breedA[,2]<=intervals[i+1])
	 insideA<-insideA[,1]
	 r2_A[i,1]<-mean(insideA)
	 r2_A[i,2]<-(intervals[i]+intervals[i+1])/2
	 }
	 		 
	if(length(which(is.nan(r2_A[,1])==TRUE))>0){
	warning('\n','Mean LD is not available for one or more intervals due to lack of loci on those intervals.')
	}
	 	 
}

	if(LD_summary==TRUE){
	Mean_D<-mean(ld_data[,10])
	Mean_r<-mean(ld_data[,11])
	Mean_r2<-mean(ld_data[,12])
	cat('\n','****  Linkage disequilibrium output **** ',fill=TRUE,'\n')
	cat('Method:',method,fill = TRUE)
	cat('No. marker:',(ncol(mat_after_maf)/2),fill = TRUE)
	cat('No. marker pairs:',npair,fill = TRUE)
	cat('\n','"D" summary:',fill=TRUE)
	print(summary(ld_data[,10]))
	cat('\n','"r" summary:',fill=TRUE)
	print(summary(ld_data[,11]))
	cat('\n','"r2" summary:',fill=TRUE)
	print(summary(ld_data[,12]))
	}



	#write to output
	if(!missing(saveAt)){
	# ld_data[,2:length(ld_data[1,])]<-formatC(ld_data[,2:length(ld_data[1,])],digits=4, format="f")
	colnames(ld_data) <- c("Pair No.","--p1--","--q1--","--p2--","--q2--","--Hap11--","--Hap12--","--Hap21--","--Hap22--","--D--","--r--","--r2--")
	 write.table(formatC(ld_data,digits=4, format="f") ,file=outFile,row.names = FALSE,quote=FALSE)
	  cat('\n','LD data has been written to the output file:',outFile,fill = TRUE)
	}
 


# RETURN SECTION	
if(Ld_decay==TRUE){
Final<-list(Mean_r2=Mean_r2,ld_data=ld_data,ld_decay=r2_A)
} else 
Final<-list(Mean_r2=Mean_r2,ld_data=ld_data)

return(Final)

}

