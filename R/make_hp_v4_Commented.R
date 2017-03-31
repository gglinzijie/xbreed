#' Create historical population
#'
#' Simulates historical generations to establish mutation-drift equilibrium and create linkage disequilibrium.

################################################
#################PARAMETERS#####################
#################################################


#' @param hpsize Size of historical population. Range: \eqn{0<\code{hpsize} \leq 10000}.
#' @param ng  Number of generations. Range: \eqn{0 \leq \code{ng} \leq 10000}.

#////////////////////////////////////////////////////////
#' @param genome \code{data.frame} specifying genome parameters with dimension \eqn{n*p}  where n is number of chromosomes and p={6} are columns defined as following:\cr
#'   "chr" Chromosome id starting from {1} to the number of specified chromosomes (Max={100}). \cr
#'   "len" Chromosome length in cM in range \eqn{0< \code{len} \leq 10000}.\cr
#'   "nmrk" Number of markers in each chromosome in range \eqn{1\leq \code{nmrk} \leq 10000}.\cr
#'   "mpos" Marker position along chromosome with options:
	#' \itemize{
	#'\code{'rnd'} - samples marker positions from uniform distribution.\cr
	#'\code{'even'} - markers are evenly spaced.
	#' }
#'   "nqtl" Number of qtl in each chromosome in range \eqn{1\leq \code{nqtl} \leq 2000}.\cr
#'   "qpos" QTL position along chromosome with options: \code{'rnd'} and \code{'even'} similar to "mpos".\cr

#////////////////////////////////////////////////////////

#' @param h2 Narrow sense heritability. Range: \eqn{0< \code{h2} \leq 1}.

#' @param d2 Ratio of dominance variance to phenotypic variance. Range: \eqn{0\leq \code{d2} \leq 1}.

#' @param phen_var Phenotypic variance. Range: \eqn{0< \code{x} \leq 10000}.

#' @param mutr Mutation rate. Range: \eqn{0\leq \code{mutr} \leq 0.01}.


#' @param laf \emph{Optional} Loci (marker and qtl) allele frequencies in the first historical generation with options:
	#' \itemize{
	#'\item{"rnd"}  {where allele frequencies will be sampled from uniform distribution}.
	#'\item{0 < \code{laf} <1}  {where all loci allele frequencies will be equal to \code{laf}}.
	#' }
#'	Default: "rnd"

#' @param sel_seq_qtl \emph{Optional} Select segregating qtl in the last historical generation. QTL with minor allele frequency larger than or equal to \code{sel_seq_qtl}  will be selected.  Range: \eqn{0\leq \code{x} \leq 0.4999}. Default: {0}

#' @param sel_seq_mrk \emph{Optional} Select segregating markers in the last historical generation. Markers with minor allele frequency larger than or equal to \code{sel_seq_mrk}  will be selected. Range: \eqn{0\leq \code{x} \leq 0.4999}. Default: {0}

#' @param saveAt \emph{Optional} (\code{character}). Name to be used to save output files.

################################################
#################RETURN/KEY/EXPORT##############
################################################

#' @return \code{list} with data of last generation of historical population.\cr
	#' \describe{
	#'\item{$hploci}{Genotype (both marker (SNP) and QTL) of individuals. Coded as {11,12,21,22} where first allele is always paternal allele}.
	#'\item{$hp_mrk}{Marker genotye of individuals}.
	#'\item{$hp_qtl}{QTL genotye of individuals}.
	#'\item{$freqQTL}{QTL allele frequency}.
	#'\item{$freqMrk}{Marker allele frequency}.
	#'\item{$linkage_map_qtl}{Linkage map for qtl}.
	#'\item{$linkage_map_mrk}{Linkage map for marker}.
	#'\item{$linkage_map_qtl_mrk}{Integrated linkage map for both marker and qtl}.
	#'\item{$allele_effcts}{QTL allelic effects}.
	#'\item{$hp_data}{Individuals data except their genotypes}.
	#'\item{$trait}{Trait specifications}.
	#'\item{$mut_data}{Mutation data.}
	#' }

#'
#' @export make_hp


################################################
#################EXAMPLS########################
################################################

#' @examples
#' # # # EXAMPLE 1 Simulation of a historical population
#' #for an additive trait (h2=0.3) for 10 generations.
#' # Two  chromosome with different parameters
#'
 #'genome<-data.frame(matrix(NA, nrow=2, ncol=6))
 #'names(genome)<-c("chr","len","nmrk","mpos","nqtl","qpos")
 #'genome$chr<-c(1,2)
 #'genome$len<-c(100,200)
 #'genome$nmrk<-c(100,100)
 #'genome$mpos<-c("rnd","even")
 #'genome$nqtl<-c(50,50)
 #'genome$qpos<-c("even","rnd")
 #'genome
 #'
 #'hp<-make_hp(hpsize=100,
 #'  ng=10,h2=0.3,phen_var=1 ,genome=genome,
 #'  mutr=2.5e-4,saveAt="hp1")
#'
#' head(hp$hp_data)
#' head(hp$freqQTL)
#' head(hp$linkage_map_qtl_mrk)

#'
#'
#' # # # EXAMPLE 2 Simulation of a historical population for a trait with both additive and
#' # dominance effects (h2=0.3, d2=0.1).
#' # All loci will have the same allele frequencies in the first generation.
#' # Segregating markers and qtl with MAF>0.1 will be selected in the last historical population.
#'
#' genome<-data.frame(matrix(NA, nrow=3, ncol=6))
#' names(genome)<-c("chr","len","nmrk","mpos","nqtl","qpos")
#' genome$chr<-c(1,2,3)
#' genome$len<-c(12,8,11)
#' genome$nmrk<-c(140,80,73)
#' genome$mpos<-c("rnd","even","even")
#' genome$nqtl<-c(40,65,24)
#' genome$qpos<-rep("rnd",3)
#' genome
#'
 #'hp2<-make_hp(hpsize=100,
 #'   ng=10,h2=0.3,d2=0.1,phen_var=1 ,genome=genome,
 #'   mutr=2.5e-4,sel_seq_qtl=0.1,sel_seq_mrk=0.1,
 #'   laf=0.1,saveAt="hp2")
#'
#'head(hp2$hp_data)
#'head(hp2$freqQTL)
#'head(hp2$linkage_map_qtl_mrk)

################################################
#################DETAILS########################
################################################

#' @details
#' \bold{Historical population} \cr
#'In order to create initial linkage disequilibrium (LD)  and to establish mutation-drift equilibrium, a historical population is simulated by considering only two evolutionary forces: mutation and drift. Mutation constantly introduces new variation and genetic drift shifts the variation to fixation. Offspring are produced by random union of gametes, each from the male and female gametic pools. Population size is constant over discrete generations with equal number of males and females. \cr
#'\cr
#' \bold{Genome} \cr
#' A wide range of parameters can be specified for simulating the genome, such as: number of chromosomes, markers and QTL, location of markers and QTL, mutation rate and initial allelic frequencies. This flexibility permits for a wide variety of genetic architectures to be considered. No allelic effects are simulated for markers, so they are treated as neutral. For QTL, additive allelic effects are sampled from gamma distribution with shape and scale parameters of {0.4} and {1.66}, respectively. This provided an L-shaped distribution of QTL effects. To simulate dominance effects, first dominance degrees \eqn{h_i}  are sampled from normal distribution (\eqn{N(0.5,1)}) then absolute dominance effects are considered as \eqn{d_i=h_i.|a_i|} where \eqn{|a_i|} is the absolute value of the additive effect. Thus, additive and dominance effects are dependent. Next, additive and dominance effects are scaled such that user defined \eqn{h^2} and \eqn{d^2} are met. Trait phenotypes are simulated by adding a standard normal residual effect to the genotypic value of each individual. \cr
#'\cr
#' One important aspect of genome simulation is to model the recombination appropriately to produce realistic level of LD, given the recent and past population structures. \code{make_hp} models crossover process, using a Poisson model. This is done by sampling the number of crossovers from a Poisson distribution and then the crossovers are located randomly across the chromosome. Because the input map is in centiMorgan it is straightforward to take into account the pattern of recombination hotspots and cold spots along the genome by adjusting the distances between markers. To establish mutation-drift equilibrium in the historical generations recurrent mutation model is used. The recurrent mutation model assumes that a mutation alters an allelic state to another and does not create a new allele. In the recurrent model, transition probabilities from one allelic state to another are assumed equal. Different mutation rates for simulated loci can be specified. The number of mutations is sampled from a Poisson distribution and it is assumed that mutation rates are equal for all loci. \cr
#'\cr
#' In conclusion the main features for \code{make_hp} are as following:
#' \itemize{
#'  \item{}{Multiple chromosomes with similar or different genome length in cM, each with different or similar density of markers and QTL maps, can be generated.}
#'  \item{}{Trait of interest can be controlled by additive or both additive and dominance effects.}
#' }

# keep wrapping for later on
make_hp<-function(hpsize,ng,genome,h2,d2,phen_var,mutr,laf,sel_seq_qtl,sel_seq_mrk,saveAt) {

  # # loading .dll this will be removed in the package
# dyn.load('hp.dll')
# is.loaded("hp")
# library(stats)
# dyn.load('cf.dll')
# is.loaded("cf")

# library(data.table)

if (is.loaded("hp")==FALSE) {
cat('HP subroutine from Fortran has not been loaded',fill=TRUE)
}

# ---------------------START-----------------
# INPUT by user control section
# ---------------------START-----------------

# hpsize
	if(missing(hpsize)) {
	stop('---hpsize is missing')
	}

	# check later on
	if(hpsize<1 | hpsize>10000 | floor(hpsize)!=hpsize) {
	stop('---hpsize should be a positive integer in range 1-10000')
	}

	# check later on
	if(hpsize%%2!=0) {
	hpsize<-hpsize+1
	}

# ng
	if(missing(ng)) {
	stop('---ng is missing')
	}

	if(ng<1 | ng>20000 | floor(ng)!=ng) {
	stop('---Number of generations for historical population should be an integer in range 1-20000')
	}

# genome
	if(missing(genome)) {
	stop('---argument "genome" is missing',fill=TRUE)
	}

	control_names<-c("chr","len","nmrk","mpos","nqtl","qpos")
	missingCols <- setdiff(control_names, colnames(genome))
	if(length(missingCols)>0){
			for (er_test in 1:length(missingCols)){
			cat(missingCols[er_test],'in argument "genome" is missing',fill=TRUE)
			}
			stop('---Error in argument "genome" ')
	}

	# genome$chr
		if(all(genome$chr == floor(genome$chr))==FALSE){
		stop('---genome$chr in argument "genome" should be integer')
		}

		t1<-genome$chr==c(1:length(genome$chr))
		if(all(t1)!=TRUE){
		stop('---Chromosome ID in argument "genome" should start from 1 to the number of chromosomes')
		}

		if(length(genome$chr)>100){
		stop('---Maximum number of chromosomes allowed in argument "genome" is 100')
		}

	# genome$len
		tst<-genome$len
		if(any(tst<=0) | any(tst>10000)){
		stop('---Chromosome length (genome$len) in argument "genome" shoule be in range 0<R<10000')
		}

	# genome$nmrk
		tst<-genome$nmrk
		if(any(tst<1) | any(tst>10000)){
		stop('---Number of marker (genome$nmrk) in argument "genome" shoule be in range 1-10000')
		}

	# genome$mpos
		if(length(setdiff(genome$mpos, c('rnd','even')))>0){
		stop('---Marker position (genome$mpos) in argument "genome" shoule be either "rnd" or "even".')
		}

	# genome$nqtl
		tst<-genome$nqtl
		if(any(tst<1) | any(tst>2000)){
		stop('---Number of qtl (genome$nqtl) in argument "genome" shoule be in range 1-2000')
		}

	# genome$qpos
		if(length(setdiff(genome$qpos, c('rnd','even')))>0){
		stop('---QTL position (genome$qpos) in argument "genome" shoule be either "rnd" or "even".')
		}

	# make correct order
	genome<-data.frame(genome$chr,genome$len,genome$nmrk,genome$mpos,genome$nqtl,genome$qpos)
	colnames(genome)<-c("chr","len","nmrk","mpos","nqtl","qpos")
	genome
	genome[,4]<-as.character(genome[,4])
	genome[,6]<-as.character(genome[,6])

	# class(genome[,1])
	# class(genome[,2])
	# class(genome[,3])
	# class(genome[,4])
	# class(genome[,5])
	# class(genome[,6])

	# mutr
	if(missing(mutr)) {
	cat('---mutr is missing, it has been set to default value of 2.5e-5',fill=TRUE)
	mutr<- 2.5e-5
	}
	if(mutr<0 | 0.01<mutr) {
	stop('---Mutation rate should be in range 0-0.01')
	}

	# sel_seq_qtl
	if(missing(sel_seq_qtl)) {
	cat('---sel_seq_qtl is missing, it has been set to default value of 0',fill=TRUE)
	sel_seq_qtl<- 0
	}
	if(sel_seq_qtl<0 | 0.49999<sel_seq_qtl) {
	stop('---sel_seq_qtl should be in range 0-0.4999')
	}

	# sel_seq_mrk
	if(missing(sel_seq_mrk)) {
	cat('---sel_seq_mrk is missing, it has been set to default value of 0',fill=TRUE)
	sel_seq_mrk<- 0
	}
	if(sel_seq_mrk<0 | 0.4999<sel_seq_mrk) {
	stop('---sel_seq_mrk should be in range 0-0.4999')
	}

	# laf
	if(missing(laf)) {
	cat('---laf is missing, loci allele frequency in the first historical generation will be sampled from uniform distribution',fill=TRUE)
	laf<- 'rnd'
	}
	if(class(laf)=="character"){
		if( laf!='rnd') {
		stop('---argument "laf" can be defined as "rnd" or in a range 0<laf<1')
		}
	}

	if(class(laf)=="numeric"){
		if( laf<0 | 1<laf) {
		stop('---argument "laf" can be defined as "rnd" or in a range 0<=laf<1')
		}
	}

	# Trait
		# h2
		if(missing(h2)) {
		stop('---h2 is missing')
		}
		if(h2<=0 | 1<h2) {
		stop('---h2 should be in range 0<h2<=1')
		}

		# d2
		if(missing(d2)) {
		d2<-0
		}

		if(missing(d2) | d2==0) {
		addTra<-TRUE
		} else {
		addTra<-FALSE
		}

		if(!missing(d2)) {
			if(d2<0 | 1<d2) {
			stop('---d2 should be in range 0<=d2<=1')
			}
		}

		# phen_var
		if(missing(phen_var)) {
		cat('---phen_var is missing, it has been set to default value of 1',fill=TRUE)
		}
		if(phen_var<=0 | 100000<phen_var) {
		stop('---phen_var should be in range 0<phen_var<=100000')
		}


	# saveAt
	if(!missing(saveAt)){
		if(is.character(saveAt)==FALSE){
		stop('---Define a name for the saveAt argument as type "character"')
		}
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

	bin_snp<-function(mat){
	s1<-seq(1,ncol(mat),2)
	s2<-seq(2,ncol(mat),2)
	a1<-mat[,s1]+mat[,s2]
	a1[a1==3]=1
	a1[a1==4]=0
	snp_code<-a1
	return(snp_code)
	 }

	 # function for TGV
	calc_TGV<-function(mat,add_eff,dom_eff){
	mat<-bin_snp(mat)
	xprogeny<-mat
	xprogeny<-as.matrix(xprogeny)
	q1<-xprogeny
	for (i in 1:length(xprogeny[,1])){
	ti<-xprogeny[i,]
	two<-which(ti==2)
	one<-which(ti==1)
	zero<-which(ti==0)
	q1[i,two]<-add_eff_1[two]
	q1[i,one]<-dom_eff[one]
	q1[i,zero]<-add_eff_2[zero]
	}
	xprogeny<-q1
	tgv<- rowSums(xprogeny)
	return(tgv)
	}

# Necessary internal functions---Finish



  # genome
nqtl_input<-sum(genome[,5])
nmarker_input<-sum(genome[,3])
nchr<-length(genome[,1])
nloci=nqtl_input+nmarker_input
# arg4<-sum(genome[,c(3,6)])
arg4<-(nqtl_input+nmarker_input)*2
arg5<-hpsize*arg4
poiss<-rpois((hpsize*ng),lambda=(2*(nqtl_input+nmarker_input)*mutr))
# index of loci to be used for crossover
	index_loc<-matrix(nrow=nchr,ncol=2)
	a<-genome[,3]+genome[,5]
	a1<-cumsum(a)
	a2<-c(0,a1)
	a3<-a2+1
	index_loc[,1]<-a3[-c(length(a3))]
	index_loc[,2]<-a1
	index_loc
	index_loc<-index_loc*2
	index_loc[,1]<-index_loc[,1]-1
	# cat(index_loc,fill=TRUE)
	rcs<-nrow(index_loc)
	ccs<-ncol(index_loc)
	index_loc<-as.integer(unlist(index_loc))
	arg8<-index_loc


if(any(poiss>0)){
	a<-which(poiss>0)
	number_mut<-poiss[a]
	we<-seq(1,(length(poiss)+hpsize),hpsize)
	naslina<-findInterval(a,  we,all.inside = FALSE)
	ng1<-which(1<=a & a<=hpsize)
	# peyda kardan heyvan
	if(length(ng1)>0){
	heyvan<-a[ng1]
	a_d<-a[-ng1]
	naslina_d<-naslina[-ng1]
	a_d<-a_d-(we[naslina_d-1]+hpsize)
	indi<-a_d+1
	indi<-c(heyvan,indi)
	outi<-cbind(naslina,indi,number_mut)
	outi<-as.data.frame(outi)
	}
		# peyda kardan heyvan
	if(length(ng1)==0){
	a_d<-a
	naslina_d<-naslina
	a_d<-a_d-(we[naslina_d-1]+hpsize)
	indi<-a_d+1
	outi<-cbind(naslina,indi,number_mut)
	outi<-as.data.frame(outi)
	}
}   else {
outi<-matrix(0,3,3)
outi<-as.data.frame(outi)
# cat('hech')
}
outi

if (laf=='rnd'){
arglaf<-1 #two options 1:unif(rnd) 2:equal
arglaf_F<-0.01
} else {
arglaf<-2 #two options 1:unif(rnd) 2:equal
arglaf_F<-laf
}

# -------------
# outi[31,2]=25
# outi[31,3]=5
# outi
# --------------
# ------
arg_r_mut<-nrow(outi)
arg_c_mut<-ncol(outi)
outi2<-unlist(outi)
outi2<-as.integer(outi2)
arg7<-outi2
  # Passing to .Fortran
  cat('Historical pop is initialized...',fill=TRUE)

# hp_str<-.Fortran("hp",hpsize= as.integer(hpsize), nmarker = as.integer(nmarker_input),nqtl=as.integer(nqtl_input),loc_c=as.integer(arg4),ng=as.integer(ng),in_r_mu=as.integer(arg_r_mut),in_c_mu=as.integer(arg_c_mut),in_pois=arg7,hp_loci=integer(arg5),method_laf=as.integer(arglaf),arglaf_F=as.double(arglaf_F))


hp_str<-.Fortran("hp",hpsize= as.integer(hpsize), nmarker = as.integer(nmarker_input),nqtl=as.integer(nqtl_input),loc_c=as.integer(arg4),ng=as.integer(ng),in_r_mu=as.integer(arg_r_mut),in_c_mu=as.integer(arg_c_mut),in_pois=arg7,hp_loci=integer(arg5),method_laf=as.integer(arglaf),arglaf_F=as.double(arglaf_F),rcs=as.integer(rcs),cross=arg8)


# After Fortran calcs
lmi<-matrix(as.numeric(hp_str[[9]]),nrow = hpsize,ncol = arg4)
# cat('HP is established',fill=TRUE)

#-----------------------------------
# After getting last generation of hp
#-----------------------------------
# Assign QTL and marker
index1<-sample((nloci/2),nqtl_input,replace=FALSE)
index1<-sort(index1)
index2<-index1*2
index3<-index2-1
indexQTL<-c(index2,index3)
indexQTL<-sort(indexQTL)
qtlLoci<-lmi[,indexQTL]
dim(qtlLoci)
mrkLoci<-lmi[,-indexQTL]
dim(mrkLoci)

# freq1Mrk<-Calc_freq(mrkLoci)
	 # Passing to .Fortran for calc Freq--start
loci_mat<-mrkLoci
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1Mrk<-as.numeric(outfreq[[5]])
	 # Passing to .Fortran for calc Freq--finish

freq2Mrk<-1-freq1Mrk

# SAMPLING POSITION FOR MARKER AND QTL
# # original
# linkage_map_qtl<-list()
# for (i in 1:length(genome[,1])){
	# if(genome[i,6]=='even'){
	 # dis<-genome[i,2]/genome[i,5]
	# linkage_map_qtl[[i]]<-seq(dis,genome[i,2],dis)
	# }
    # else
	# linkage_map_qtl[[i]]<-sort(runif(genome[i,5],1,genome[i,2]))
# }

linkage_map_qtl<-list()
for (i in 1:length(genome[,1])){
	if(genome[i,6]=='even'){
	 dis<-genome[i,2]/genome[i,5]
	linkage_map_qtl[[i]]<-seq(dis,genome[i,2],dis)
		# check 1
		if(length(linkage_map_qtl[[i]])<genome[i,5]){
		cat('#pos in linkage_map_qtl in chromosome',i,'is smaller than genome[,5]')
		stop('Internal error in linkage_map_qtl')
		}
		# check 2
		if(length(unique(linkage_map_qtl[[i]]))<genome[i,5]){
		cat('some qtl in chromosome',i,'have the same position')
		stop('Solutions:change qpos or length of chromosome in argument"genome"')
		}

	}
    else
	linkage_map_qtl[[i]]<-sort(runif(genome[i,5],0,genome[i,2]))
		# check 1
		if(length(linkage_map_qtl[[i]])<genome[i,5]){
		cat('#pos in linkage_map_qtl in chromosome',i,'is smaller than genome[,5]')
		stop('Internal error in linkage_map_qtl')
		}
		# check 2
		if(length(unique(linkage_map_qtl[[i]]))<genome[i,5]){
		cat('some qtl in chromosome',i,'have the same position')
		stop('Solutions:change qpos or length of chromosome in argument"genome"')
		}

}

# # original
# linkage_map_mrk<-list()
# for (i in 1:length(genome[,1])){
	# if(genome[i,4]=='even'){
	# dis<-genome[i,2]/genome[i,3]
	# linkage_map_mrk[[i]]<-seq(dis,genome[i,2],dis)
	# }
    # else
	# linkage_map_mrk[[i]]<-sort(runif(genome[i,3],1,genome[i,2]))
# }

linkage_map_mrk<-list()
for (i in 1:length(genome[,1])){
	if(genome[i,4]=='even'){
	dis<-genome[i,2]/genome[i,3]
	linkage_map_mrk[[i]]<-seq(dis,genome[i,2],dis)
		# check 1
		if(length(linkage_map_mrk[[i]])<genome[i,3]){
		cat('#pos in linkage_map_mrk in chromosome',i,'is smaller than genome[,3]')
		stop('Internal error in linkage_map_mrk')
		}
		# check 2
		if(length(unique(linkage_map_mrk[[i]]))<genome[i,3]){
		cat('some marker in chromosome',i,'have the same position')
		stop('Solutions:change mpos or length of chromosome in argument"genome"')
		}
	}

    else
	linkage_map_mrk[[i]]<-sort(runif(genome[i,3],0,genome[i,2]))
	   # check 1
		if(length(linkage_map_mrk[[i]])<genome[i,3]){
		cat('#pos in linkage_map_mrk in chromosome',i,'is smaller than genome[,3]')
		stop('Internal error in linkage_map_mrk')
		}
		# check 2
		if(length(unique(linkage_map_mrk[[i]]))<genome[i,3]){
		cat('some marker in chromosome',i,'have the same position')
		stop('Solutions:change mpos or length of chromosome in argument"genome"')
		}
}

# # Test of linkage map
# for (i in 1:length(genome[,1])){
	# if(genome[i,3]==genome[i,5] & genome[i,4]==genome[i,6]){
	# checker<-linkage_map_mrk[[i]]==linkage_map_qtl[[i]]
	# internal_lm_mrk<-linkage_map_mrk[[i]]
	# inde<-which(checker==TRUE)

		# if (length(inde)>0){
		# internal_lm_mrk[inde]<-internal_lm_mrk[inde]+0.0001
		# linkage_map_mrk[[i]]<-internal_lm_mrk
	# }
# }
# }

# cat('Sampling position for Marker and QTL ...',fill=TRUE)

# FINAL CHECK FOR POSITIONS START
	 all_poses<-c()
	 sizez<-c()
	 for (i in 1:nchr){
	 ap1<-matrix(as.numeric(unlist(linkage_map_qtl[[i]])))
	 ap2<-matrix(as.numeric(unlist(linkage_map_mrk[[i]])))
	 sizez[i]<-length(ap1)+length(ap2)
	 all_poses<-rbind(all_poses,ap1,ap2)
	 }
	 # length(all_poses)
	 # length(unique(all_poses))
	 all_poses2<-unique(all_poses)

	 for_mess<-1
	while(length(all_poses)!=length(all_poses2)){
	indexp<-match(all_poses2,all_poses)
	all_poses[-indexp]<-all_poses[-indexp]-0.0001
	# dak<-1:length(all_poses)
	# changed_pos<-dak[-indexp]
		sizez<-cumsum(sizez)
		sizez<-c(0,sizez)
		sizez
		a1<-1
		b1<-2
		all_poses_t<-list()
		for (i in 1:nchr){
		all_poses_t[[i]]<-all_poses[(sizez[a1]+1):sizez[b1]]
		a1<-1+a1
		b1<-1+b1
		if(a1==nchr+1) break
		}

		for (i in 1:nchr){
		uv<-c(genome[i,5],genome[i,3])
		uv<-cumsum(uv)
		uv<-c(0,uv)
		a1<-1
		b1<-2
		linkage_map_qtl[[i]]<-all_poses_t[[i]][(uv[a1]+1):uv[b1]]
		a1<-1+a1
		b1<-1+b1
		linkage_map_mrk[[i]]<-all_poses_t[[i]][(uv[a1]+1):uv[b1]]
		if(a1==nchr+1) break
		}

	 all_poses<-c()
	 sizez<-c()
	 for (i in 1:nchr){
	 ap1<-matrix(as.numeric(unlist(linkage_map_qtl[[i]])))
	 ap2<-matrix(as.numeric(unlist(linkage_map_mrk[[i]])))
	 sizez[i]<-length(ap1)+length(ap2)
	 all_poses<-rbind(all_poses,ap1,ap2)
	 }
	 all_poses2<-unique(all_poses)
	 if(for_mess==1){
	 # message('\n','"loci positions have been modified slightly to avoid identical positions along the genome"')
	 }
	  for_mess<-1+ for_mess
} #end while

# FINAL CHECK FOR POSITIONS FINISH



qtlLoci_chri<-list()
mrkLoci_chri<-list()

# create names for QTL and marker
QTL_Label<-paste('Q',1:nqtl_input,sep='')
MRK_Label<-paste('M',1:nmarker_input,sep='')



if (sel_seq_qtl==0){
		dim(qtlLoci)

	 # index of qtl in each chromosome
	  selcted_qtlindex_chri<-list()
		for (i in 1:length(genome[,1])){
	selcted_qtlindex_chri[[i]]<-1:genome[i,5]
		}
	   	# extracting QTL number in each chr
		nqtl_chr<-genome[,5]
		# extracting QTL loci for each chromosome
		index_qtl_chr1<-cumsum(nqtl_chr*2)
		index_qtl_chr2<-c(0,index_qtl_chr1)
		a1<-1
		b1<-2
		for (i in 1:nchr){
		qtlLoci_chri[[i]]<-qtlLoci[,(index_qtl_chr2[a1]+1):index_qtl_chr2[b1]]
		a1<-1+a1
		b1<-1+b1
		if(a1==nchr+1) break
		}
		dim(qtlLoci_chri[[i]])
		# create names for QTL and marker
		# nqtl after maf
		Finall_nqtl<-sum(nqtl_chr)
	    QTL_Label<-paste('Q',1:Finall_nqtl,sep='')
	   nqtl<-Finall_nqtl
	}

if (sel_seq_qtl>0){
cat('Extracting segregating qtl loci ...',fill=TRUE)
		# freq1<-Calc_freq(qtlLoci)

# Passing to .Fortran for calc Freq--start
loci_mat<-qtlLoci
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish

        freq2<-1-freq1
		summary(freq1)
        freqmatrix<-matrix(0,ncol=2,nrow=length(freq1))
        freqmatrix[,1]<-freq1
        freqmatrix[,2]<-freq2
		MAFv<-apply(freqmatrix,1,min)
		MAFv1<-MAFv[MAFv>=sel_seq_qtl]
		index<-which(MAFv%in%MAFv1)
		selcted_qtl_index<-index
		cat('----No. segregating QTL:',length(index),'out of',length(freq1),fill=TRUE)
		if(length(index)==0){
		stop('ERROR: There is no segregating QTL')
		}

		# extract loci passed seq
		new_index<-index*2
		index<-c(new_index-1,new_index)
		index<-sort(index)
		mat_after_maf<-qtlLoci[,index]
		# summary(Calc_freq(mat_after_maf))
		qtlLoci<-mat_after_maf
		# dim(qtlLoci)

	 # index of qtl in each chromosome
	  selcted_qtlindex_chri<-list()
	  	for_chr<-c(0,genome[,5])
		for_chr<-cumsum(for_chr)
		a1<-1
		b1<-2
		for (i in 1:length(for_chr)){
		# cat(a1,fill=TRUE)
		# cat(b1,fill=TRUE)
		internal_index<-which(selcted_qtl_index>for_chr[a1] & selcted_qtl_index<=for_chr[b1])
	selcted_qtlindex_chri[[i]]<-selcted_qtl_index[internal_index]
		a1<-1+a1
		b1<-1+b1
		}
	selcted_qtlindex_chri
	for (i in 1:length(selcted_qtlindex_chri)){
	selcted_qtlindex_chri[[i]]<-selcted_qtlindex_chri[[i]]-for_chr[i]
	}
	selcted_qtlindex_chri<-selcted_qtlindex_chri[-length(selcted_qtlindex_chri)]


    	# extracting QTL number in each chr
		for_chr<-c(0,genome[,5])
		for_chr<-cumsum(for_chr)
		a1<-1
		b1<-2
		nqtl_chr<-c()
		for (i in 1:length(for_chr)){
		# cat(a1,fill=TRUE)
		# cat(b1,fill=TRUE)
	nqtl_chr[i]<-sum(selcted_qtl_index>for_chr[a1] & selcted_qtl_index<=for_chr[b1])
		a1<-1+a1
		b1<-1+b1
		}
		nqtl_chr<-nqtl_chr[-length(nqtl_chr)]

		# extracting QTL loci for each chromosome
		index_qtl_chr1<-cumsum(nqtl_chr*2)
		index_qtl_chr2<-c(0,index_qtl_chr1)
		a1<-1
		b1<-2
		for (i in 1:nchr){
			# cat(a1,fill=TRUE)
		# cat(b1,fill=TRUE)
		qtlLoci_chri[[i]]<-qtlLoci[,(index_qtl_chr2[a1]+1):index_qtl_chr2[b1]]
		a1<-1+a1
		b1<-1+b1
		if(a1==nchr+1) break
		}
		dim(qtlLoci_chri[[i]])
		# create names for QTL and marker
		# nqtl after maf
		Finall_nqtl<-sum(nqtl_chr)
	    QTL_Label<-paste('Q',1:Finall_nqtl,sep='')
	   nqtl<-Finall_nqtl
	}



if (sel_seq_mrk==0){
		dim(mrkLoci)

	 # index of qtl in each chromosome
	  selcted_mrkindex_chri<-list()
		for (i in 1:length(genome[,1])){
	selcted_mrkindex_chri[[i]]<-1:genome[i,3]
		}
		# extracting marker number in each chr
		nmrk_chr<-genome[,3]

		# extracting MRK loci for each chromosome
		index_mrk_chr1<-cumsum(nmrk_chr*2)
		index_mrk_chr2<-c(0,index_mrk_chr1)
		a1<-1
		b1<-2
		for (i in 1:nchr){
		mrkLoci_chri[[i]]<-mrkLoci[,(index_mrk_chr2[a1]+1):index_mrk_chr2[b1]]
		a1<-1+a1
		b1<-1+b1
		if(a1==nchr+1) break
		}
		dim(mrkLoci_chri[[i]])
		# create names for marker
		# nqtl after maf
		Finall_nmarker<-sum(nmrk_chr)
		MRK_Label<-paste('M',1:Finall_nmarker,sep='')
	nmarker<-Finall_nmarker
	}

if (sel_seq_mrk>0){
 cat('Extracting segregating markers ...',fill=TRUE)
		freq1<-Calc_freq(mrkLoci)

# Passing to .Fortran for calc Freq--start
loci_mat<-mrkLoci
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish


        freq2<-1-freq1
         # summary(freq1)
        freqmatrix<-matrix(0,ncol=2,nrow=length(freq1))
        freqmatrix[,1]<-freq1
        freqmatrix[,2]<-freq2
		MAFv<-apply(freqmatrix,1,min)
		MAFv1<-MAFv[MAFv>=sel_seq_mrk]
		index<-which(MAFv%in%MAFv1)
		selcted_mrk_index<-index
		cat('----No. segregating markers:',length(index),'out of',length(freq1),fill=TRUE)
		if(length(index)==0){
		stop('ERROR: There is no segregating marker')
		}

		# extract loci passed seq
		new_index<-index*2
		index<-c(new_index-1,new_index)
		index<-sort(index)
		mat_after_maf<-mrkLoci[,index]
		# summary(Calc_freq(mat_after_maf))
		mrkLoci<-mat_after_maf
		dim(mrkLoci)

		# index of mrk in each chromosome
	  selcted_mrkindex_chri<-list()
	  	for_chr<-c(0,genome[,3])
		for_chr<-cumsum(for_chr)
		a1<-1
		b1<-2
		for (i in 1:length(for_chr)){
		# cat(a1,fill=TRUE)
		# cat(b1,fill=TRUE)
		internal_index<-which(selcted_mrk_index>for_chr[a1] & selcted_mrk_index<=for_chr[b1])
	 selcted_mrkindex_chri[[i]]<-selcted_mrk_index[internal_index]
		a1<-1+a1
		b1<-1+b1
		}

 	selcted_mrkindex_chri
	for (i in 1:length(selcted_mrkindex_chri)){
	selcted_mrkindex_chri[[i]]<-selcted_mrkindex_chri[[i]]-for_chr[i]
	}
	selcted_mrkindex_chri<-selcted_mrkindex_chri[-length(selcted_mrkindex_chri)]

		# extracting MRK in each chr
		for_chr<-c(0,genome[,3])
		for_chr<-cumsum(for_chr)
		a1<-1
		b1<-2
		nmrk_chr<-c()
		for (i in 1:length(for_chr)){
		# cat(a1,fill=TRUE)
		# cat(b1,fill=TRUE)
	nmrk_chr[i]<-sum(selcted_mrk_index>for_chr[a1] & selcted_mrk_index<=for_chr[b1])
		a1<-1+a1
		b1<-1+b1
		}
		nmrk_chr<-nmrk_chr[-length(nmrk_chr)]

		index_mrk_chr1<-cumsum(nmrk_chr*2)
		index_mrk_chr2<-c(0,index_mrk_chr1)
		a1<-1
		b1<-2
		for (i in 1:nchr){
		mrkLoci_chri[[i]]<-mrkLoci[,(index_mrk_chr2[a1]+1):index_mrk_chr2[b1]]
		a1<-1+a1
		b1<-1+b1
		if(a1==nchr+1) break
		}

		# create names for QTL and marker
		# nqtl after maf
		Finall_nmarker<-sum(nmrk_chr)
		MRK_Label<-paste('M',1:Finall_nmarker,sep='')
	nmarker<-Finall_nmarker
}



	qtl_pos_chri<-list()
	mrk_pos_chri<-list()
	for (i in 1:nchr) {
	  # QTL linkage map
		qtl_pos_chri[[i]]<-as.numeric(unlist(linkage_map_qtl[[i]][selcted_qtlindex_chri[[i]]]))
	  # MRK linkage map
		mrk_pos_chri[[i]]<-as.numeric(unlist(linkage_map_mrk[[i]][selcted_mrkindex_chri[[i]]]))
		}


	# length(qtl_pos_chri[[1]])
	# length(qtl_pos_chri[[2]])
	# length(mrk_pos_chri[[1]])
	# length(mrk_pos_chri[[2]])


# sequence loci in each chromosome
seque_loci<-list()
for (i in 1:nchr){
	posqtl1<-qtl_pos_chri[[i]]
	posqtl2<-posqtl1
	posqtldo<-c(posqtl1,posqtl2)
	posqtl<-sort(posqtldo)
	length(posqtl)
	qtlseq_bapos<-rbind(posqtl,qtlLoci_chri[[i]])
	dim(qtlseq_bapos)
	qtlseq_bapos<-as.data.frame(qtlseq_bapos)

	posmarker1<-mrk_pos_chri[[i]]
	posmarker2<-posmarker1
	posmarkerdo<-c(posmarker1,posmarker2)
	posmarker<-sort(posmarkerdo)
	length(posmarker)
	mrkseq_bapos<-rbind(posmarker,mrkLoci_chri[[i]])
	dim(mrkseq_bapos)
	mrkseq_bapos<-as.data.frame(mrkseq_bapos)

	pos<-c(posmarker,posqtl)
	pos<-sort(pos)
	length(pos)

	sequence_loci<-matrix(ncol=length(posqtl)+length(posmarker),nrow=(dim(qtlLoci_chri[[i]])[1])+1)
	dim(sequence_loci)
	sequence_loci<-as.data.frame(sequence_loci)

	sequence_loci[1,]<-pos
	sequence_loci2<-sequence_loci
	sequence_loci2[sequence_loci2[1,]%in%posqtl]<-qtlseq_bapos
	sequence_loci2[sequence_loci2[1,]%in%posmarker]<-mrkseq_bapos
	sequence_loci<-sequence_loci2[-1,]
	sequence_loci<-as.matrix(sequence_loci)
	seque_loci[[i]]<-sequence_loci

	# sequence_loci2[1:10,1:5]
	# mrkseq_bapos[1:10,1:5]

	# dim(sequence_loci)
	# sequence_loci[1:10,1:20]
}


# QTL and Mrk and sequence loci in all chromosomes
qtlLoci<-c()
mrkLoci<-c()
seqLoci<-c()
for (i in 1:nchr){
	qseq<-qtlLoci_chri[[i]]
	qtlLoci<-cbind(qtlLoci,qseq)
	mseq<-mrkLoci_chri[[i]]
	mrkLoci<-cbind(mrkLoci,mseq)
	all_seq<-seque_loci[[i]]
	seqLoci<-cbind(seqLoci,all_seq)
}
dim(qtlLoci)
dim(mrkLoci)
dim(seqLoci)
class(seqLoci)

# To be used for output  DO NOT DELETE
# freq1Mrk<-Calc_freq(mrkLoci)

# Passing to .Fortran for calc Freq--start
loci_mat<-mrkLoci
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1Mrk<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish
freq2Mrk<-1-freq1Mrk
summary(freq1Mrk)

# # ORIGINAL
# # Create Trait
# cat('Simulating trait ...',fill=TRUE)
# var_add<-h2*phen_var
# var_dom	<-d2*phen_var
# add_eff_1<-rgamma(nqtl, shape=0.4,1.66)
# e<-sample(length(add_eff_1),length(add_eff_1)/2,replace=FALSE)
# add_eff_1[e]<-add_eff_1[e]*-1
# add_eff_2<-add_eff_1*(-1)
# dom_degree<-rnorm(nqtl,mean=0.5)
# dom_eff<-dom_degree*abs(add_eff_1)

# Create Trait
cat('Simulating trait ...',fill=TRUE)
var_add<-h2*phen_var
if(addTra==FALSE){
var_dom	<-d2*phen_var
} else {
var_dom<-0
d2<-0
}
add_eff_1<-rgamma(nqtl, shape=0.4,1.66)
e<-sample(length(add_eff_1),length(add_eff_1)/2,replace=FALSE)
add_eff_1[e]<-add_eff_1[e]*-1
add_eff_2<-add_eff_1*(-1)
dom_degree<-rnorm(nqtl,mean=0.5)
if(addTra==FALSE){
dom_eff<-dom_degree*abs(add_eff_1)
} else {
dom_eff<-rep(0,length(dom_degree))
}



#TBV based on allele frequency
# freq1<-Calc_freq(qtlLoci)

# Passing to .Fortran for calc Freq--start
loci_mat<-qtlLoci
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish

freq2<-1-freq1
freq1qtl<-freq1
freq2qtl<-freq2
Snp_BreedB<-bin_snp(qtlLoci)
alpha=abs(add_eff_1)+((freq2-freq1)*dom_eff)
#Breeding Value
xprogeny<-Snp_BreedB
xprogeny<-as.matrix(xprogeny)
q1<-xprogeny
for (i in 1:length(xprogeny[,1])){
ti<-xprogeny[i,]
two<-which(ti==2)
one<-which(ti==1)
zero<-which(ti==0)
q1[i,two]<-2*freq2[two]*alpha[two]
q1[i,one]<-(freq2[one]-freq1[one])*alpha[one]
q1[i,zero]<--2*freq1[zero]*alpha[zero]
}
xprogeny<-q1
tbv<- rowSums(xprogeny)
var(tbv)


#Dominance Deviation
xprogeny<-Snp_BreedB
xprogeny<-as.matrix(xprogeny)
q1<-xprogeny
for (i in 1:length(xprogeny[,1])){
ti<-xprogeny[i,]
two<-which(ti==2)
one<-which(ti==1)
zero<-which(ti==0)
q1[i,two]<-(freq2[two]**2)*-2*dom_eff[two]  #A1A1  -2q^2d   A1A2 2pqd   A2A2 -2p^2d
q1[i,one]<-2*freq1[one]*freq2[one]*dom_eff[one]
q1[i,zero]<--2*(freq1[zero]**2)*dom_eff[zero]
}
xprogeny<-q1
dom_dev<- rowSums(xprogeny)
var(dom_dev)

# scaling
ratio_add<-sqrt(var_add)/sqrt(var(tbv))
if(addTra==FALSE){
ratio_dom<-sqrt(var_dom)/sqrt(var(dom_dev))
} else {
ratio_dom<-0
}
dom_eff_new<-dom_eff*ratio_dom
dom_eff<-dom_eff_new
alpha_new<-alpha*ratio_add
alpha<-alpha_new
add_eff_1<-alpha-((freq2-freq1)*dom_eff)
add_eff_2<-add_eff_1*(-1)

#Breeding Value
alpha=abs(add_eff_1)+((freq2-freq1)*dom_eff)
xprogeny<-Snp_BreedB
xprogeny<-as.matrix(xprogeny)
q1<-xprogeny
for (i in 1:length(xprogeny[,1])){
ti<-xprogeny[i,]
two<-which(ti==2)
one<-which(ti==1)
zero<-which(ti==0)
q1[i,two]<-2*freq2[two]*alpha[two]
q1[i,one]<-(freq2[one]-freq1[one])*alpha[one]
q1[i,zero]<--2*freq1[zero]*alpha[zero]
}
xprogeny<-q1
tbv<- rowSums(xprogeny)
var(tbv)

#Dominance Deviation
xprogeny<-Snp_BreedB
xprogeny<-as.matrix(xprogeny)
q1<-xprogeny
for (i in 1:length(xprogeny[,1])){
ti<-xprogeny[i,]
two<-which(ti==2)
one<-which(ti==1)
zero<-which(ti==0)
q1[i,two]<-(freq2[two]**2)*-2*dom_eff[two]  #A1A1  -2q^2d   A1A2 2pqd   A2A2 -2p^2d
q1[i,one]<-2*freq1[one]*freq2[one]*dom_eff[one]
q1[i,zero]<--2*(freq1[zero]**2)*dom_eff[zero]
}
xprogeny<-q1
dom_dev<- rowSums(xprogeny)
var(dom_dev)

# ---------------------START-----------------
# DATA RETURN BACK BY FUNCTION
# ---------------------START-----------------

#Fill in data of last hp
cat('Output data preparation ...',fill=TRUE)
ID<-1:hpsize
# Sex<-sample(c('F','M'),hpsize,replace=T)
Sex1<-rep('F',hpsize/2)
Sex2<-rep('M',hpsize/2)
Sex<-sample(c(Sex1,Sex2),hpsize,replace=FALSE)
	vv<-phen_var-(var_add+var_dom)
Res<-rnorm(length(ID),mean=0,sd=sqrt(vv))
	TGV<-calc_TGV(qtlLoci,add_eff_1,dom_eff)
Phen<-TGV+Res
TBV<-tbv
data_hp<-data.frame(ID,Sex,Phen,Res,TBV)
names(data_hp)<-c('ID','Sex','Phen','Res','TBV')



# genotypes sequence
all_loci<-data.frame(ID,seqLoci)
all_loci<-as.matrix(all_loci)

# QTL Loci
QTL_loci<-data.frame(ID,qtlLoci)

# MRK Loci
MRK_loci<-data.frame(ID,mrkLoci)

# freq QTL
ID_qtl<-QTL_Label
ID_qtl_chr<-rep(genome[,1],nqtl_chr)
freqQTL<-data.frame(ID_qtl,ID_qtl_chr,freq1qtl,freq2qtl)
names(freqQTL)<-c('ID','Chr','Freq.Allele1','Freq.Allele2')

# freq MRK
ID_mrk<-MRK_Label
ID_mrk_chr<-rep(genome[,1],nmrk_chr)
freqMRK<-data.frame(ID_mrk,ID_mrk_chr,freq1Mrk,freq2Mrk)
names(freqMRK)<-c('ID','Chr','Freq.Allele1','Freq.Allele2')

 # QTL linkage map
ID_qtl<-QTL_Label
ID_qtl_chr<-rep(genome[,1],nqtl_chr)
internal_map<-c()
	 for (i in 1:length(linkage_map_qtl)){
	a<-matrix(as.numeric(unlist(qtl_pos_chri[[i]])))
	internal_map<-rbind(internal_map,a)
	}
	internal_map_QTL<-data.frame(ID_qtl,ID_qtl_chr,internal_map)
	names(internal_map_QTL)<-c('ID','Chr','Position')

	# MRK linkage map
ID_mrk<-MRK_Label
ID_mrk_chr<-rep(genome[,1],nmrk_chr)
	internal_map<-c()
	 for (i in 1:length(linkage_map_mrk)){
	a<-matrix(as.numeric(unlist(mrk_pos_chri[[i]])))
	internal_map<-rbind(internal_map,a)
	}
	internal_map_MRK<-data.frame(ID_mrk,ID_mrk_chr,internal_map)
	names(internal_map_MRK)<-c('ID','Chr','Position')

 # QTL/MRK linkage map
 qtl_mrk_lm<-list()
 	for (i in 1:length(genome[,1])){
labelqtl<-subset(internal_map_QTL,internal_map_QTL[,2]==genome[i,1])
labelmrk<-subset(internal_map_MRK,internal_map_MRK[,2]==genome[i,1])
QtlMrk_lm<-rbind(labelqtl,labelmrk )
QtlMrk_lm<-QtlMrk_lm[order(QtlMrk_lm$Position),]
qtl_mrk_lm[[i]]<-QtlMrk_lm
	}
internal_map<-c()
	 for (i in 1:length(qtl_mrk_lm)){
	a<-qtl_mrk_lm[[i]]
	internal_map<-rbind(internal_map,a)
	}
	internal_map_qtl_mrk<-internal_map
	names(internal_map_qtl_mrk)<-c('ID','Chr','Position')


#Allele effects
ID_qtl<-QTL_Label
ID_qtl_chr<-rep(genome[,1],nqtl_chr)

allele_effcts<-data.frame(ID_qtl,ID_qtl_chr,add_eff_1,add_eff_2,dom_eff)
names(allele_effcts)<-c('ID','Chr','Add.Allele1','ADD.Allele2','Dominance')

# Trait
trait<-data.frame(h2=h2,dom_her=d2,phen_var=phen_var)
# Genome
genome_l<-genome
for(g in 1:length(genome_l[,1])){
genome_l[g,3]<-length(which(internal_map_MRK[,2]==genome_l[g,1]))
genome_l[g,5]<-length(which(internal_map_QTL[,2]==genome_l[g,1]))
}


# ---------------------FINISH-----------------
# DATA RETURN BACK BY FUNCTION
# ---------------------FINISH-----------------


#WRITE TO OUTPUT
	if(!missing(saveAt)){
	cat('Writing output files ...',fill=TRUE)
	# ALL Loci
	outFile_allLoci<-paste(saveAt,'_genome.txt',sep='')
	writeLines(c('ID  Genotypes (paternal allele, maternal allele) ...'),outFile_allLoci)
	write.table(all_loci,file=outFile_allLoci,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)


	# QTL Loci
	outFile_QTL_loci<-paste(saveAt,'_qtl.txt',sep='')
	writeLines(c('ID  Genotypes (paternal allele, maternal allele) ...'),outFile_QTL_loci)
	write.table(QTL_loci,file=outFile_QTL_loci,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)

	# MRK Loci
	outFile_MRK_loci<-paste(saveAt,'_mrk.txt',sep='')
	writeLines(c('ID  Genotypes (paternal allele, maternal allele) ...'),outFile_MRK_loci)
	write.table(MRK_loci,file=outFile_MRK_loci,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)

	# freq QTL
	dom<-format(freqQTL,  justify = "right")
	outFileFreqQTL<-paste(saveAt,'_freq_qtl.txt',sep='')
    write.table(dom,file=outFileFreqQTL,row.names=FALSE,col.names=TRUE,quote = FALSE)

	# freq MRK
	dom<-format(freqMRK,  justify = "right")
	outFileFreqMRK<-paste(saveAt,'_freq_mrk.txt',sep='')
	write.table(dom,file=outFileFreqMRK,row.names=FALSE,col.names=TRUE,quote = FALSE)

    # QTL linkage map
	dom<-format(internal_map_QTL,  justify = "right")
	outFilelm_QTL<-paste(saveAt,'_lm_qtl.txt',sep='')
    write.table(dom,file=outFilelm_QTL,row.names=FALSE,col.names=TRUE,quote = FALSE)

	# MRK linkage map
	dom<-format(internal_map_MRK,  justify = "right")
	outFilelm_MRK<-paste(saveAt,'_lm_mrk.txt',sep='')
	write.table(dom,file=outFilelm_MRK,row.names=FALSE,col.names=TRUE,quote = FALSE)

	# QTL/MRK linkage map
	dom<-format(internal_map_qtl_mrk,  justify = "right")
	outFilelm_QTL_MRK<-paste(saveAt,'_lm_qtl_mrk.txt',sep='')
	write.table(dom,file=outFilelm_QTL_MRK,row.names=FALSE,col.names=TRUE,quote = FALSE)


	# QTL allele effects
	dom<-format(allele_effcts,  justify = "right")
	outFile_Alele_eff<-paste(saveAt,'_effect_qtl.txt',sep='')
	write.table(dom,file=outFile_Alele_eff,row.names=FALSE,col.names=TRUE,quote = FALSE)


	# DATA of last hp
	dom<-format(data_hp,  justify = "right")
	outFile_data_hp<-paste(saveAt,'_last_hp_data.txt',sep='')
	write.table(dom,file=outFile_data_hp,row.names=FALSE,col.names=TRUE,quote = FALSE)


}

# RETURN SECTION
MakeHP<-list(hpsize=hp_str[[1]],ng=hp_str[[5]],nqtl=nqtl,nmarker=nmarker,mutr=mutr,mut_data=outi,hploci=all_loci,hp_mrk=MRK_loci,hp_qtl=QTL_loci,freqQTL=freqQTL,freqMrk=freqMRK,linkage_map_qtl=internal_map_QTL,linkage_map_mrk=internal_map_MRK,linkage_map_qtl_mrk=internal_map_qtl_mrk,allele_effcts=allele_effcts,hp_data=data_hp,trait=trait,genome=genome_l)

cat('Establishment of historical population completed',fill=TRUE)

return(MakeHP)

}


.onAttach = function(library, pkg){
packageStartupMessage('("|-----------------------------------------------------|")', appendLF=TRUE)
packageStartupMessage('("|                      xbreed                         |")', appendLF=TRUE)
packageStartupMessage('("|    Genomic simulation of purebreds and crossbreds   |")', appendLF=TRUE)
packageStartupMessage('("|               March 2017 Version 1.0.0              |")', appendLF=TRUE)
packageStartupMessage('("|                                                     |")', appendLF=TRUE)
packageStartupMessage('("|             H.Esfandyari,A.C.Sorensen               |")', appendLF=TRUE)
packageStartupMessage('("| Center for Quantitative Qenetics and Genomics (QGG) |")', appendLF=TRUE)
packageStartupMessage('("|             Aarhus University,Denmark               |")', appendLF=TRUE)
packageStartupMessage('("|                                                     |")', appendLF=TRUE)
packageStartupMessage('("|-----------------------------------------------------|")', appendLF=TRUE)
packageStartupMessage('("|Questions and bugs: esfandyari.hadi@gmail.com        |")', appendLF=TRUE)
packageStartupMessage('("|Development of xbreed was supported by GenSAP.       |")', appendLF=TRUE)
packageStartupMessage('("|-----------------------------------------------------|")', appendLF=TRUE)
  invisible()
}

