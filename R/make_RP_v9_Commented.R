#' Make recent population
#'
#' Creates recent population similar to function \code{\link{sample_hp}} with some modifications.

################################################
#################PARAMETERS#####################
#################################################

#' @param sh_out (\code{list}) Output of function  \code{\link{sample_hp}} or \code{\link{make_rp}}.

# -------------male_founders ---------------
#' @param Male_founders  (\code{data.frame}) Data frame with {1} row and {3} colomns as following:\cr
#'   Colomn 1) "number" is the number of male individuals to be selected. \cr

#'   Colomn 2) "generation" is the generation number from which males will be selected.\cr

#'   Colomn 3) "select" indicates the type of selection with options:
	#' \itemize{
	#'\item{"rnd"}  {Select individuals randomly}.
	#'\item{"phen"}  {Select individuals based on their phenotypes}.
	#'\item{"tbv"}  {Select individuals based on their true breeding value (tbv)}.
	#'\item{"gebv"}  {Select individuals based on their genomic estimated breeding value (gebv). \bold{Note:} If "gebv" was not selection criteria for the population from which founders are being selected, then "gebv" as selection criteria will be ignored}.
	#' }

#'   Colomn 4) "value" Indicates to select high: "h" or low: "l" values. Note: This colomn is ignored if indviduals are selected randomly.\cr
# -------------male_founders finish---------------
# -------------Female_founders---------------

#' @param Female_founders  (\code{data.frame}) Data frame with {1} row and {3} columns containing information to select female founders. Details are similar to argument \code{Male_founders}.

# -------------Female_founders finish---------------
#' @param ng  Number of generations. Range: \eqn{1 \leq \code{ng} \leq 500}.

#' @param litter_size  Litter size or the number of progeny per dam. Range: \eqn{1 \leq \code{x} \leq 200}.


# -------------selection start---------------
#' @param Selection (\code{data.frame}) Data frame with {2} rows and {3} columns. First row is for the selection design of males and second row is for the selection design of females. The columns are as following:\cr
#'   Column 1) "size" is the number of individuals to be selected as sires/dams. \cr
#'   Column 2) "type" indicates the type of selection with options:
	#' \itemize{
	#'\item{"rnd"}  {Select individuals randomly}.
	#'\item{"phen"}  {Select individuals based on their phenotypes}.
	#'\item{"tbv"}  {Select individuals based on their true breeding value (tbv)}.
	#'\item{"gebv"}  {Select individuals based on their genomic estimated breeding value (gebv)}.
	#'}

#'   Column 3) "value" Indicates to select "h" or "l" values. Note: This colomn is ignored if indviduals are selected randomly.\cr
# -------------selection finished---------------

# -------------training start---------------
#' @param Training \emph{Optional} (\code{data.frame}) Data frame with {1} row and {8} columns. The columns are as following:\cr
#'  Column 1) "size" is the number of individuals to be selected for training. \cr
#'  Column 2) "sel" \emph{Optional} (\code{character}) Indicates the type of the selection of individuals for training. The possible options are:
	#' \itemize{
	#'\item{"rnd"}  {Select individuals for training randomly}.
	#'\item{"min_rel_mrk"}  {Select individuals for training, where genomic relationship among individuals based on marker information is minimum}.
	#'\item{"max_rel_mrk"}  {Select individuals for training, where genomic relationship among individuals based on marker information is maximum}.
	#'\item{"min_rel_qtl"}  {Select individuals for training, where genomic relationship among individuals based on qtl information is minimum}.
	#'\item{"max_rel_qtl"}  {Select individuals for training, where genomic relationship among individuals based on qtl information is maximum}.
	#'}
#'	Default: "rnd" \cr

#'   Column 3) "method" \emph{Optional} (\code{character}) Method used for the estimation of marker effects. The possible options are:
	#' \itemize{
	#'\item{"BRR"}  {Gaussian prior}.
	#'\item{"BayesA"}  {scaled-t prior}.
	#'\item{"BL"}  {Double-Exponential prior}.
	#'\item{"BayesB"}  {two component mixture prior with a point of mass at zero and a sclaed-t slab}.
	#'\item{"BayesC"}  {two component mixture prior with a point of mass at zero and a Gaussian slab}.
	#'}
	#'	Default: "BRR" \cr

#'  Column 4) "nIter" \emph{Optional} The number of iterations. Default: \eqn{1500} \cr
#'   Column 5) "burnIn" \emph{Optional} The number of burn-in. Default: \eqn{500} \cr
#'   Column 6) "thin" \emph{Optional} The number of thinning. Default: \eqn{5} \cr
#'   Column 7) "save" \emph{Optional} This may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs. Default:"Out_BGLR" \cr
#'   Column 8) "show" \emph{Optional} (\code{Logical}) if TRUE the iteration history is printed. Default: \code{TRUE}. \cr
#' \bold{Note:} This argument is compulsory if \code{"type"} in argument \code{Selection} is "gebv".  More details about the argument can be found in package \pkg{BGLR}.
# -------------training finished---------------

#' @param saveAt \emph{Optional} (\code{character}). Name to be used to save output files.

# -------------sh_output starts---------------
#' @param rp_output \emph{Optional} (\code{data.frame}). Data frame to specify generations indexs and type of data to be written to output files. User can define which type of data and which generation to be written to output files. The possible options are:\cr
#'\cr
#'   "data"      Individuals data except their genotypes. \cr
#'   "qtl"       QTL genotye of individuals coded as {11,12,21,22}. \cr
#'   "marker"    Marker genotye of individuals. \cr
#'   "seq"       Genotype (both marker (SNP) and QTL) of individuals.\cr
#'   "freq_qtl"  QTL allele frequency. \cr
#'   "freq_mrk"  Marker allele frequency. \cr
#' \bold{Note:} Both arguments \code{rp_output} and \code{saveAt} should present in the function in order to write the output files.
# -------------sh_output finished---------------

#' @param Display \emph{Optional} (\code{Logical}) Display summary of the simulated generations if is not \code{FALSE}. Default: \code{TRUE}.

################################################
#################RETURN/KEY/EXPORT##############
################################################

#' @return \code{list} with all data of simulated generations.\cr
	#' \describe{
	#'\item{$output}{(\code{list}) Two-level list  (\code{$output[[x]][[y]]}) containing information about simulated generations. First index (x) indicates generation number. It should be noted that as data for base generation (0) is also stored by the function, to retrive data for a specific generation, index should be equal to generation number plus one. As an example to observe data for generation 2 index should be 3 i.e, \code{$output[[3]]$data}. Second index (y) that ranges from {1} to {6} contain the information as following:
		#' \itemize{
		#'\item{\code{$output[[x]]$data}}  {Individuals data except their genotypes. Here x is the generation index}
		#'\item{\code{$output[[x]]$qtl}}  {QTL genotye of individuals.}.
		#'\item{\code{$output[[x]]$mrk}}   {Marker genotye of individuals}.
		#'\item{\code{$output[[x]]$sequ}}  {Genotype (both marker (SNP) and QTL) of individuals}.
		#'\item{\code{$output[[x]]$freqQTL}}   {QTL allele frequency}.
		#'\item{\code{$output[[x]]$freqMRK}}   {Marker allele frequency}.
		#'}
	#'}
	#'\item{$summary_data}{Data frame with summary of simulated generations}.
	#'\item{$linkage_map_qtl}{Linkage map for qtl}.
	#'\item{$linkage_map_mrk}{Linkage map for marker}.
	#'\item{$linkage_map_qtl_mrk}{Integrated linkage map for both marker and qtl}.
	#'\item{$allele_effcts}{QTL allele effects}.
	#'\item{$trait}{Trait specifications}.
	#'\item{$genome}{Genome specifications}.
	#'}

#' @useDynLib xbreed
#' @export make_rp
#' @import utils
#' @import BGLR
#' @import stats
#' @seealso \code{\link{sample_hp}}

################################################
#################EXAMPLS########################
################################################

#' @examples
#' # # # Simulation of a population where founders are from a population created by function sample_hp.
#'
#'# CREATE HISTORICAL POPULATION
#'
#' genome<-data.frame(matrix(NA, nrow=2, ncol=6))
#' names(genome)<-c("chr","len","nmrk","mpos","nqtl","qpos")
#' genome$chr<-c(1,2)
#' genome$len<-c(12,8)
#' genome$nmrk<-c(140,80)
#' genome$mpos<-c("rnd","even")
#' genome$nqtl<-c(25,25)
#' genome$qpos<-rep("rnd",2)
#' genome
#'
#'hp<-make_hp(hpsize=100
#',ng=10,h2=0.3,phen_var=1
#',genome=genome,mutr=5*10**-4,sel_seq_qtl=0.05,sel_seq_mrk=0.05,laf=0.5)
#'
#'# # MAKE FIRST RECENT POPULATION USING FUNCTION sample_hp
#'
#'Male_founders<-data.frame(number=50,select='rnd')
#'Female_founders<-data.frame(number=50,select='rnd')
#'
#'# Selection scheme in each generation of recent population
#'Selection<-data.frame(matrix(NA, nrow=2, ncol=2))
#'names(Selection)<-c('Number','type')
#'Selection$Number[1:2]<-c(50,50)
#'Selection$type[1:2]<-c('rnd','rnd')
#'Selection
#'
#' RP_1<-sample_hp(hp_out=hp,Male_founders=
#' Male_founders,Female_founders=Female_founders,
#' ng=4,Selection=Selection,Training=Training,
#' litter_size=3,Display=TRUE)
#'
#'# # MAKE SECOND RP (RP2) USING FUNCTION make_hp
#'   # Select founders
#'   # Select 30 males based on 'tbv' from generation  2 of RP1.
#'   # Select 40 females based on 'phen' from generation 4 of RP1.
#'
#'Males<-data.frame(number=30,generation=2,select='tbv',value='h')

#'Females<-data.frame(number=40,generation=4,select='phen',value='l')
#'
#'# Selection scheme for RP2
#'
#'	# Selection of 20 sires and 50 dam
#'	# Selection criteria is "tbv" for sires and "phen" for dams
#'
 #'Selection<-data.frame(matrix(NA, nrow=2, ncol=3))
#' names(Selection)<-c('Number','type','Value')
#' Selection$Number[1:2]<-c(20,50)
#' Selection$type[1:2]<-c('tbv','phen')
#' Selection$Value[1:2]<-c('h','h')
#' Selection
 #'
#'	# Save "data" and "qtl" for first and last generation of RP1
#'
#'rp2_output<-data.frame(matrix(NA, nrow=2, ncol=2))
#'names(rp2_output)<-c("data","qtl")
#'rp2_output[,1]<-c(1,4) # Save data for generations 1 and 4
#'rp2_output[,2]<-c(1,4) # Save qtl genotype for generations 1 and 4
#'rp2_output
#'
#'RP_2<-make_rp(sh_out=RP_1,Male_founders=Males,
#'Female_founders=Females,Selection=Selection,
#'ng=4,litter_size=4,saveAt='RP2',
#'rp_output=rp2_output)
#'
#' # Some output display
#'
#' RP_2$summary_data
#' RP_2$output[[1]]$data # Data for base Generation
#' RP_2$output[[2]]$freqQTL # qtl frequencies for 1st Generation
#' RP_2$output[[4]]$freqMRK # Marker frequencies for 3rd Generation
#' RP_2$linkage_map_qtl
#' RP_2$allele_effcts

################################################
#################DETAILS########################
################################################

#' @details
#' Function \code{make_rp} is used to create recent population(s) similar to function \code{sample_hp}. The difference between these two functions is that, for function \code{sample_hp}, male and female founders are always from the last generation of historical population. However, in function \code{make_rp}, male and female founders can be selected from any generation of population created by function \code{sample_hp} or in the second usage, founders can be selected from any generation of population created by the function \code{make_rp} itself. So, basically function \code{make_rp} can be used multiple times to sample individuals from the desired generation of population created by function \code{sample_hp} or function \code{make_rp}. See details for function \code{\link{sample_hp}} and the package vignette for more clarification. \cr


#sh_out2 for female founder
# keep wrapping for later on
make_rp_half<-function(sh_out,sh_out2,Male_founders,Female_founders,ng,litter_size,Selection,Training,saveAt,rp_output,Display) {

# dyn.load('sh.dll')
# is.loaded("sh")

# dyn.load('cf.dll')
# is.loaded("cf")



# ---------------------START-----------------
# INPUT by user control section
# ---------------------START-----------------
cat('Controlling input data ...',fill=TRUE)

# sh_out
	if(missing(sh_out)) {
	stop('--- argument "sh_out" is missing')
	}

	if(class(sh_out)!='list') {
	stop('--- argument "sh_out" should be the output of function sample_hp')
	}

# Male_founders
	if(missing(Male_founders)) {
	stop('--- argument "Male_founders" is missing')
	}


	outforLD<-sh_out


	tst<-outforLD$output[[(Male_founders[,2]+1)]]$data
	t_mad<-subset(tst,tst[,5]=='M')
	t_mad<-length(t_mad[,5])




	test_Males <- c('rnd','tbv','phen','gebv')
	if(Male_founders[,3]%in%test_Males==FALSE){
	stop('--- poosible options in argument "Male_founders" are "phen","tbv","rnd","gebv"')
	}

	test_Males <- c('tbv','phen','gebv')
	if(Male_founders[,3]%in%test_Males==TRUE){
	if(length(Male_founders[1,])!=4){
	stop('\n','--- If males in argument "Male_founders" are selected based on "phen","tbv" or "gebv", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}

	test_Males <- c('tbv','phen','gebv')
	if(Male_founders[,3]%in%test_Males==TRUE){
	if(Male_founders[,4]!='h' & Male_founders[,4]!='l' ){
	stop('\n','--- If males in argument "Male_founders" are selected based on "phen","tbv" or "gebv", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}

	test_Selection <- c('gebv')
	tst_gebv<-colnames(outforLD$summary_data)
	if(Male_founders[2]%in%test_Selection==TRUE){
	if(is.na(match('GEBV',tst_gebv))){
	stop('\n','--- "gebv" as selection criteria for male founders in argument "Male_founders" is not allowed due to different user input in function sample_hp/make_rp')
	}
	}

# Female_founders
	if(missing(Female_founders)) {
	stop('--- argument "Female_founders" is missing')
	}


	outforLD2<-sh_out2


	tst<-outforLD$output[[(Female_founders[,2]+1)]]$data
	t_mad<-subset(tst,tst[,5]=='F')
	t_mad<-length(t_mad[,5])



	test_Males <- c('rnd','tbv','phen','gebv')
	if(Female_founders[,3]%in%test_Males==FALSE){
	stop('--- poosible options in argument "Female_founders" are "phen","tbv","rnd","gebv"')
	}

	test_Males <- c('tbv','phen','gebv')
	if(Female_founders[,3]%in%test_Males==TRUE){
	if(length(Female_founders[1,])!=4){
	stop('\n','--- If females in argument "Female_founders" are selected based on "phen","tbv" or "gebv", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}


	test_Males <- c('tbv','phen','gebv')
	if(Female_founders[,3]%in%test_Males==TRUE){
	if(Female_founders[,4]!='h' & Female_founders[,4]!='l' ){
	stop('\n','--- If females in argument "Female_founders" are selected based on "phen","tbv" or "gebv", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}

	test_Selection <- c('gebv')
	tst_gebv<-colnames(outforLD$summary_data)
	if(Female_founders[2]%in%test_Selection==TRUE){
	if(is.na(match('GEBV',tst_gebv))){
	stop('\n','--- "gebv" as selection criteria for female founders in argument "Female_founders" is not allowed due to different user input in function sample_hp/make_rp')
	}
	}


	# ng
	if(missing(ng)) {
	stop('\n','--- ng is missing')
	}

	if(ng<1 | ng>500 | floor(ng)!=ng) {
	stop('\n','--- Number of generations should be an integer in range 1-500')
	}

	# litter_size
	if(missing(litter_size)) {
	stop('\n','--- litter_size is missing')
	}

	if(litter_size<1 | litter_size>200 | floor(litter_size)!=litter_size) {
	stop('\n','--- argument "litter_size" should be an integer in range 1<=x<=200')
	}

	# Selection
	if(missing(Selection)) {
	stop('\n','--- argument "Selection" is missing')
	}

		# males
		if(Selection[1,1]<1 | floor(Selection[1,1])!=Selection[1,1] ) {
		stop('\n','--- Number of selected sires in argument "Selection" should be an integer grater than 1')
		}

		# females
		if(Selection[2,1]<1 | floor(Selection[2,1])!=Selection[2,1] ) {
		stop('\n','--- Number of selected dams in argument "Selection" should be an integer grater than 1')
		}

		test_Selection <- c('rnd','tbv','phen','gebv')
		if(any(Selection[,2]%in%test_Selection==FALSE)){
		stop('\n','---  selection criteria for selected sires/dams in argument "Selection" can be  "rnd","phen","tbv" or "gebv"')
		}

		test_Selection <- c('tbv','phen','gebv')
		if(any(Selection[,2]%in%test_Selection==TRUE)){
		if(length(Selection[1,])!=3){
		stop('\n','--- if selection criteria for selected sires/dams in argument "Selection" is  "phen", "tbv" or "gebv", it should be defined to select "h" (high) or "l" (low) values.')
		}
		}

		test_Selection <- c('tbv','phen','gebv')
		if(any(Selection[,2]%in%test_Selection==TRUE)){
		if(any(Selection[,3]!='h') & any(Selection[,3]!='l' )){
		stop('\n','--- if selection criteria for selected sires/dams in argument "Selection" is  "phen" or "tbv" or "gebv", it should be defined to select "h" (high) or "l" (low) values.')
		}
		}



		# test_Selection <- c('gebv')
		# tst_gebv<-colnames(outforLD$summary_data)
		# if(any(Selection[,2]%in%test_Selection==TRUE)){
		# if(is.na(match('GEBV',tst_gebv))){
		# stop('\n','---"gebv" as selection criteria for selected sires/dams in argument "Selection" is not allowed due to different user input in function Sample_hp/Make_rp')
		# }
		# }

	# Training
	test_Selection <- c('gebv')
	if(any(Selection[,2]%in%test_Selection==TRUE)){
	if(missing(Training)){
	stop('\n','--- If selection criteria for selected sires/dams in argument "Selection" is  "gebv", argument "Training" should be defined.')
	}
	}

	if(any(Selection[,2]%in%test_Selection==TRUE) & !missing(Training)){
	# control of training
	control_names<-c('size','method','sel','nIter','burnIn','thin','save','show')
	test_w<-intersect(names(Training),control_names)

	if(any(test_w=='size')==FALSE){
		stop('\n','--- size in argument "Training" is  missing.')
	}

	if(any(test_w=='method')==FALSE){
		cat('\n','method in argument "Training" is set to default method:"BRR"')
	}

	# size
	if(any(test_w=='size')){
	tra_size<-Training$size
		if(floor(tra_size)!=tra_size | tra_size<1){
		stop('\n','--- size in argument "Training" should be positive integer.')
		}


		if(tra_size>(litter_size*Female_founders[1])){
		avai<-litter_size*as.numeric(Female_founders[1])
		cat('\n','--- Training size:',tra_size,' in argument "Training" is grater than number of available individuals:',avai,fill=TRUE)
		stop("\n",'Solution 1: Increase litter size.',"\n",'Solution 2: Decrease size of training.',"\n",'Solution 3: Increase size in Female_founders.')
		}
	}

	# method
	test_method <- c('BRR', 'BayesA', 'BL', 'BayesB', 'BayesC')
	if(any(test_w=='method')){
	tra_method<-Training$method
		if(tra_method%in%test_method==FALSE){
		stop('\n','--- possible options for method in argument "Training" are:"BRR", "BayesA", "BL", "BayesB", "BayesC"')
		}
	} else {
	tra_method<-'BRR'
	}

	# sel
	test_sel <- c('rnd', 'min_rel_mrk', 'max_rel_mrk', 'min_rel_qtl', 'max_rel_qtl')
	if(any(test_w=='sel')){
	tra_sel<-Training$sel
		if(tra_sel%in%test_sel==FALSE){
		stop('\n','--- possible options for sel in argument "Training" are:"rnd", "min_rel_mrk", "max_rel_mrk", "min_rel_qtl", "max_rel_qtl"')
		}
	} else {
	tra_sel<-'rnd'
	}

	if(any(test_w=='nIter')){
	tra_nIter<-Training$nIter
		if(floor(tra_nIter)!=tra_nIter | tra_nIter<1){
		stop('\n','--- tra_nIter in argument "Training" should be positive integer.')
		}
	} else {
	tra_nIter<-1500
	}

	if(any(test_w=='burnIn')){
	tra_burnIn<-Training$burnIn
		if(floor(tra_burnIn)!=tra_burnIn | tra_burnIn<1){
		stop('\n','--- tra_burnIn in argument "Training" should be positive integer.')
		}
	} else {
	tra_burnIn<-500
	}

	if(any(test_w=='thin')){
	tra_thin<-Training$thin
		if(floor(tra_thin)!=tra_thin | tra_thin<1){
		stop('\n','--- tra_thin in argument "Training" should be positive integer.')
		}
	} else {
	tra_thin<-5
	}

	if(any(test_w=='save')){
	tra_save<-Training$save
		if(is.character(tra_save)==FALSE){
		stop('\n','--- Define a name for tra_save in argument "Training" as type "character"')
		}
	} else {
	tra_save<-'Out_BGLR'
	}

		if(any(test_w=='show')){
		tra_show<-Training$show
			if(is.logical(tra_show)==FALSE){
			stop('\n','--- tra_show in argument "Training" should be type "logical"')
			}
		} else {
		tra_show<-TRUE
		}
	}


	if(all(Selection[,2]%in%test_Selection==FALSE) & !missing(Training)){
	warning('\n','argument "Training" is ignored as selection criteria is not "gebv"')
	}


	# saveAt
	if(!missing(saveAt) & missing(rp_output)){
		stop('\n','--- argument "rp_output" is missing')
	}
	if(missing(saveAt) & !missing(rp_output)){
		stop('\n','--- argument "saveAt" is missing')
	}

	if(!missing(saveAt) & !missing(rp_output)){

			if(is.character(saveAt)==FALSE){
			stop('\n','--- Define a name for the saveAt argument as type "character"')
			}

		control_names<-c("data","qtl","marker","seq","freq_qtl","freq_mrk")
		test_w<-intersect(names(rp_output),control_names)
		if(length(test_w)==0){
			stop('\n','--- possible output files in argument "rp_output" are:"data","qtl","marker","seq","freq_qtl","freq_mrk"')
		}

		if(any(floor(rp_output)==rp_output)==FALSE){
			stop('\n','--- generation indexes in argument "rp_output" for writting output files should be integer')
		}

		if(any(floor(rp_output)>ng)==TRUE){
			stop('\n','--- there is error in generation indexes in argument "rp_output" for writting output files')
		}
	}


	# Display
	if(missing(Display)) {
	Display<-TRUE
	}
	if(!missing(Display)) {
			if(is.logical(Display)==FALSE){
			stop('\n','--- argument "Display" should be type "logical"')
			}
	}

# cat('Controlling input data done.',fill=TRUE)
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

	funi<-function(vec){
L<-sample(vec,1)
return(L)
}

list_fun<-function() {
    list(list(), list(),list(), list(), list(),list())
}
# Necessary internal functions---Finish


# Total_Total<-list(list(),list(),list(),list())
Total_Total<-replicate((ng+1), list_fun(), simplify=FALSE)

cat('Intializing base population ...',fill=TRUE)
# calculation
Males=Male_founders
Females=Female_founders


# effectsi
outforLD<-sh_out
add_eff_1<-outforLD$allele_effcts[,3]
add_eff_2<-outforLD$allele_effcts[,4]
dom_eff<-outforLD$allele_effcts[,5]
	if(outforLD$trait[2]==0){
	addTra<-TRUE
	} else {
	addTra<-FALSE
	}


# genome
# postions
	posmarker_org<-outforLD$linkage_map_mrk[,3]
	posmarker1<-posmarker_org
	posmarker2<-posmarker_org
	posmarkerdo<-c(posmarker1,posmarker2)
	posmarker<-sort(posmarkerdo)

	posqtl_org<-outforLD$linkage_map_qtl[,3]
	posqtl1<-posqtl_org
	posqtl2<-posqtl_org
	posqtldo<-c(posqtl1,posqtl2)
	posqtl<-sort(posqtldo)
	length(posqtl)

nchr<-length(outforLD$genome[,1])
chrom_length<-outforLD$genome[,2]

pos<-sort(rep(outforLD$linkage_map_qtl_mrk[,3],2))
pos_each_chr<-list()
for (i in 1:nchr){
pos_each_chr[[i]]<-subset(outforLD$linkage_map_qtl_mrk,outforLD$linkage_map_qtl_mrk[,2]==i)
}

sumqtl_mrk_chr<-table(outforLD$linkage_map_qtl_mrk[,2])
sumqtl_mrk_chr<-as.numeric(sumqtl_mrk_chr)
#sumqtl_mrk_chr
nqtl_allele<-sum(outforLD$genome[,5])*2
# nqtl_allele<-length(outforLD$freqQTL[,1])*2
nmarker_allele<-sum(outforLD$genome[,3])*2
# nmarker_allele<-length(outforLD$freqMrk[,1])*2
#nmarker_allele

# trait
trait_spec<-as.numeric(unlist(outforLD$trait))
#trait_spec
h2<-trait_spec[1]
d2<-trait_spec[2]
phen_var<-trait_spec[3]
var_add<-h2*phen_var
var_dom	<-d2*phen_var
vv<-phen_var-(var_add+var_dom)



# sire data
my_internal_data1<-outforLD$output[[Males[,2]+1]]
my_internal_data<-my_internal_data1[[1]]
my_internal_data
provided_sire_data<-subset(my_internal_data,my_internal_data[,5]=='M')
provided_sire_data
length(provided_sire_data[,1])
index_sires<-provided_sire_data[,1]
index_sires_1<-match(index_sires,my_internal_data1[[1]][,1])
index_sires_1
provided_sire_seq<-my_internal_data1[[4]][index_sires_1,]
dim(provided_sire_seq)
provided_sire_qtl<-my_internal_data1[[2]][index_sires_1,]
dim(provided_sire_qtl)
provided_sire_mrk<-my_internal_data1[[3]][index_sires_1,]
dim(provided_sire_mrk)

# sires_toC
if(Males[,3]=='rnd'){
ID_sires_toC<-sample(index_sires,Males[,1],replace=FALSE)
index<-match(ID_sires_toC,index_sires)
} else if(Males[,3]=='phen'){

	if(Males[,4]=='h'){
	sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,6]), ]
	ID_sires_toC<-sorted_sire_data[1:Males[,1],]
	index<-match(ID_sires_toC[,1],index_sires)

	}
	if(Males[,4]=='l'){
	sorted_sire_data<-provided_sire_data[order(provided_sire_data[,6]), ]
	ID_sires_toC<-sorted_sire_data[1:Males[,1],]
	index<-match(ID_sires_toC[,1],index_sires)
	}

} else if (Males[,3]=='tbv'){
	if(Males[,4]=='h'){
	sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,8]), ]
	ID_sires_toC<-sorted_sire_data[1:Males[,1],]
	index<-match(ID_sires_toC[,1],index_sires)
	}
		if(Males[,4]=='l'){
	sorted_sire_data<-provided_sire_data[order(provided_sire_data[,8]), ]
	ID_sires_toC<-sorted_sire_data[1:Males[,1],]
	index<-match(ID_sires_toC[,1],index_sires)
	}

}

sires_toC_data<-provided_sire_data[index,]
sires_toC_seq<-provided_sire_seq[index,]
# sires_toC_data
sires_toC_qtl<-provided_sire_qtl[index,]
sires_toC_mrk<-provided_sire_mrk[index,]



# dam data
my_internal_data1<-outforLD2$output[[Females[,2]+1]]
my_internal_data<-my_internal_data1[[1]]
my_internal_data
provided_dam_data<-subset(my_internal_data,my_internal_data[,5]=='F')
provided_dam_data
length(provided_dam_data[,1])
index_dams<-provided_dam_data[,1]
index_dams_1<-match(index_dams,my_internal_data1[[1]][,1])
index_dams_1
provided_dam_seq<-my_internal_data1[[4]][index_dams_1,]
dim(provided_dam_seq)
provided_dam_qtl<-my_internal_data1[[2]][index_dams_1,]
dim(provided_dam_qtl)
provided_dam_mrk<-my_internal_data1[[3]][index_dams_1,]
dim(provided_dam_mrk)

# dams_toC
if(Females[,3]=='rnd'){
ID_dams_toC<-sample(index_dams,Females[,1],replace=FALSE)
index<-match(ID_dams_toC,index_dams)
} else if(Females[,3]=='phen'){

	if(Females[,4]=='h'){
	sorted_sire_data<-provided_dam_data[order(-provided_dam_data[,6]), ]
	ID_dams_toC<-sorted_sire_data[1:Females[,1],]
	index<-match(ID_dams_toC[,1],index_dams)

	}
	if(Females[,4]=='l'){
	sorted_sire_data<-provided_dam_data[order(provided_dam_data[,6]), ]
	ID_dams_toC<-sorted_sire_data[1:Females[,1],]
	index<-match(ID_dams_toC[,1],index_dams)
	}

} else if (Females[,3]=='tbv'){
	if(Females[,4]=='h'){
	sorted_sire_data<-provided_dam_data[order(-provided_dam_data[,8]), ]
	ID_dams_toC<-sorted_sire_data[1:Females[,1],]
	index<-match(ID_dams_toC[,1],index_dams)
	}
		if(Females[,4]=='l'){
	sorted_sire_data<-provided_dam_data[order(provided_dam_data[,8]), ]
	ID_dams_toC<-sorted_sire_data[1:Females[,1],]
	index<-match(ID_dams_toC[,1],index_dams)
	}

}

dams_toC_data<-provided_dam_data[index,]
dams_toC_seq<-provided_dam_seq[index,]
# dams_toC_data
dams_toC_qtl<-provided_dam_qtl[index,]
dams_toC_mrk<-provided_dam_mrk[index,]
# dim(dams_toC_mrk)

# tamirat----------
start_data<-rbind(sires_toC_data,dams_toC_data)
total_B<-list()
 total_B$data<-start_data
 total_B$data[,4]<-0
  total_B$data

total_B$qtl<-rbind(sires_toC_qtl,dams_toC_qtl)
total_B$qtl[,2]<-0
 dim(total_B$qtl)
 class(total_B$qtl)
total_B$qtl<-as.matrix(total_B$qtl)
# class(total_B$qtl)

total_B$mrk<-rbind(sires_toC_mrk,dams_toC_mrk)
total_B$mrk[,2]<-0
# dim(total_B$mrk)
total_B$mrk<-as.matrix(total_B$mrk)

# seq
total_B$sequ<-rbind(sires_toC_seq,dams_toC_seq)
dim(total_B$sequ)
total_B$sequ<-as.matrix(total_B$sequ)
total_B$sequ[,2]<-0
# tamirat----------//////////


gene_counter<-0
 for (gene_counter in 0:ng){
cat('Generation',gene_counter,'started ')
 start_ge <- Sys.time()
  # ----------------------------
   Total_Total[[gene_counter+1]][[1]]<-total_B$data
   Total_Total[[gene_counter+1]][[2]]<-total_B$qtl
   Total_Total[[gene_counter+1]][[3]]<-total_B$mrk
   Total_Total[[gene_counter+1]][[4]]<-total_B$sequ
   # Freq QTL
	# ID_qtl<-outforLD$freqQTL[,1]
	ID_qtl<-outforLD$output[[1]][[5]][,1]
	ID_qtl_gen<-rep(gene_counter,length(ID_qtl))
	# ID_qtl_chr<-outforLD$freqQTL[,2]
	ID_qtl_chr<-outforLD$output[[1]][[5]][,3]

	#QTL allele freq
    # freq1qtl_B<-Calc_freq(total_B$qtl[,-c(1,2)])
# Passing to .Fortran for calc Freq--start
loci_mat<-total_B$qtl[,-c(1,2)]
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1qtl_B<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish

	freq2qtl_B<-1-freq1qtl_B
	freqQTL<-data.frame(ID_qtl,ID_qtl_gen,ID_qtl_chr,freq1qtl_B,freq2qtl_B)
	names(freqQTL)<-c('ID','Generation','Chr','Freq.Allele1','Freq.Allele2')
   Total_Total[[gene_counter+1]][[5]]<-freqQTL

  # Freq MRK
   	# ID_mrk<-outforLD$freqMrk[,1]
	ID_mrk<-outforLD$output[[1]][[6]][,1]
	ID_mrk_gen<-rep(gene_counter,length(ID_mrk))
	# ID_mrk_chr<-outforLD$freqMrk[,2]
	ID_mrk_chr<-outforLD$output[[1]][[6]][,3]
	#QTL allele freq
    # freq1mrk_B<-Calc_freq(total_B$mrk[,-c(1,2)])

# Passing to .Fortran for calc Freq--start
loci_mat<-total_B$mrk[,-c(1,2)]
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1mrk_B<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish

cat('.')

	freq2mrk_B<-1-freq1mrk_B
	freqMRK<-data.frame(ID_mrk,ID_mrk_gen,ID_mrk_chr,freq1mrk_B,freq2mrk_B)
	names(freqMRK)<-c('ID','Generation','Chr','Freq.Allele1','Freq.Allele2')
    Total_Total[[gene_counter+1]][[6]]<-freqMRK
 # ---------------------------------



  # cat('Getting subset for generations...',fill=TRUE)
# start.time <- Sys.time()
bb1<-subset(total_B$data,total_B$data[,4]==gene_counter)
bb2<-subset(total_B$qtl,total_B$qtl[,2]==gene_counter)
bb3<-subset(total_B$mrk,total_B$mrk[,2]==gene_counter)
bb4<-subset(total_B$sequ,total_B$sequ[,2]==gene_counter)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# cat('time.taken for subseting',time.taken,fill=TRUE)

	 if (gene_counter ==0){
	# selection of males and then top males
	males<-subset(bb1,bb1[,5]=='M')
	males_selected<-males
	# selection of females and then top females
	females<-subset(bb1,bb1[,5]=='F')
    females_selected<-females
	x<-length(females_selected[,1])*litter_size
	}

	 if (gene_counter >0){
	# selection of males and then top males
	males<-subset(bb1,bb1[,5]=='M')

		if(Selection[1,2]=='rnd'){
		ID_sires<-sample(males[,1],Selection[1,1],replace=FALSE)
		males<-males[match(ID_sires,males[,1]),]
		males_selected<-males[1:Selection[1,1],]
		}

		if(Selection[1,2]=='phen'){
			if(Selection[1,3]=='h'){
		males<-males[order(-males[,6]),]  #selection based on phen
		males_selected<-males[1:Selection[1,1],]
			}
			if(Selection[1,3]=='l'){
		males<-males[order(males[,6]),]  #selection based on phen
		males_selected<-males[1:Selection[1,1],]
			}
		}

	   if(Selection[1,2]=='tbv'){ #selection based on tbv
	   		if(Selection[1,3]=='h'){
		males<-males[order(-males[,8]),]  #selection based on phen
		males_selected<-males[1:Selection[1,1],]
			}
			if(Selection[1,3]=='l'){
		males<-males[order(males[,8]),]  #selection based on phen
		males_selected<-males[1:Selection[1,1],]
			}
    	}

		if(Selection[1,2]=='gebv'){ #selection based on tbv
	   		if(Selection[1,3]=='h'){
		males<-males[order(-males[,9]),]  #selection based on phen
		males_selected<-males[1:Selection[1,1],]
			}
			if(Selection[1,3]=='l'){
		males<-males[order(males[,9]),]  #selection based on phen
		males_selected<-males[1:Selection[1,1],]
			}
    	}


	# selection of females and then top females
	    females<-subset(bb1,bb1[,5]=='F')

		if(Selection[2,2]=='rnd'){
		ID_dams<-sample(females[,1],Selection[2,1],replace=FALSE)
		females<-females[match(ID_dams,females[,1]),]
		females_selected<-females[1:Selection[2,1],]
		}

		if(Selection[2,2]=='phen'){
			if(Selection[2,3]=='h'){
		females<-females[order(-females[,6]),]  #selection based on phen
		females_selected<-females[1:Selection[2,1],]
			}
			if(Selection[2,3]=='l'){
		females<-females[order(females[,6]),]  #selection based on phen
		females_selected<-females[1:Selection[2,1],]
			}
		}

	   if(Selection[2,2]=='tbv'){ #selection based on tbv
	   		if(Selection[2,3]=='h'){
		females<-females[order(-females[,8]),]  #selection based on phen
		females_selected<-females[1:Selection[2,1],]
			}
			if(Selection[2,3]=='l'){
		females<-females[order(females[,8]),]  #selection based on phen
		females_selected<-females[1:Selection[2,1],]
			}
    	}

		if(Selection[2,2]=='gebv'){ #selection based on tbv
	   		if(Selection[2,3]=='h'){
		females<-females[order(-females[,9]),]  #selection based on phen
		females_selected<-females[1:Selection[2,1],]
			}
			if(Selection[2,3]=='l'){
		females<-females[order(females[,9]),]  #selection based on phen
		females_selected<-females[1:Selection[2,1],]
			}
    	}
		x<-Selection[2,1]*litter_size
}

cat('.')
#mating
# x<-ndam_selected*litter_size
sires<-males_selected[,1]
if (anyNA(sires)==TRUE){
stop("\n",'Number of selected sires is less than available number of male progeny.'
,"\n",'Solution 1: Increase litter size.'
,"\n",'Solution 2: Decrease number of selected sires.')
}

dams<-females_selected[,1]
if (anyNA(dams)==TRUE){
stop("\n",'Number of selected dams is less than available number of female progeny.'
,"\n",'Solution 1: Increase litter size.'
,"\n",'Solution 2: Decrease number of selected dams.')
}



# Create id,sire,dam,generation,sex,env
 # cat('Create id,sire,dam...',fill=TRUE)
	my<-c()
	for (klm in 1:(gene_counter+1)){
	my[klm]<-length(Total_Total[[klm]][[1]][,1])
	}
counter_id<- sum(my)
id1 <- (counter_id+1)
id2 <- (counter_id+x)
id  <- id1:id2
# counter_id<-length(total_B$data[,1])
# id1 <- (counter_id+1)
# id2 <- (counter_id+x)
# id  <- id1:id2
# sire <- rep(sires,each=No_off_per_sire)
No_mat<-length(dams)*length(dams)

sire<-rep(rep(sires, each=(length(dams))),litter_size)
dam<-rep(rep(dams,length(sires)),litter_size)
generation <-gene_counter +1
sex <-sample(c('F','M'),x,replace=T)
env <-rnorm(x,mean=0,sd=sqrt(vv))
# sirem <- bb3[sire,]
sirem <- bb3[match(sire,bb3[,1]),]
# sirem[1:10,1:10]
# bb3[1:5,1:10]
# dim(sirem)
damm <- bb3[match(dam,bb3[,1]),]
# damm[1:10,1:10]
# bb3[1:5,1:10]
# dim(damm)
qtlseq_sire<-bb2[match(sire,bb2[,1]),]
# qtlseq_sire[1:10,1:10]
# bb2[1:5,1:10]
# dim(qtlseq_sire)
qtlseq_dam<-bb2[match(dam,bb2[,1]),]
# qtlseq_dam[1:10,1:10]
# bb2[1:5,1:10]
# dim(qtlseq_dam)
total_seq_sires<-bb4[match(sire,bb4[,1]),]
total_seq_sires<-total_seq_sires[,-c(1,2)]
total_seq_sires<-rbind(pos,total_seq_sires)
locii_sires<-total_seq_sires
# dim(locii_sires)
total_seq_dams<-bb4[match(dam,bb4[,1]),]
total_seq_dams<-total_seq_dams[,-c(1,2)]
total_seq_dams<-rbind(pos,total_seq_dams)
locii_dams<-total_seq_dams
 dim(locii_dams)


# ////////////////////////////////////
# COMBINED 1to2 and CROSSOVER FORTRAN///////START////
# ////////////////////////////////////
cat('.')
# New 1-----2-------------START-///////////
 start.time <- Sys.time()


# ONE TO TWO
		# Males
		locii_sires_npos<-locii_sires[-1,]
		locii_sires_npos<-as.matrix(locii_sires_npos)

		# For Fortan
		in_r_mu<-length(locii_sires_npos[,1])
		in_c_mu<-length(locii_sires_npos[1,])
		outi2<-as.integer(unlist(locii_sires_npos))
		arg7_m<-outi2
		arg4_m<-length(locii_sires_npos[1,])/2
		arg5_m<-(length(locii_sires_npos[,1])*2)*arg4_m

		# Females
		locii_dams_npos<-locii_dams[-1,]
		locii_dams_npos<-as.matrix(locii_dams_npos)

		# For Fortan
		# in_r_mu_f<-length(locii_dams_npos[,1])
		# in_c_mu_f<-length(locii_dams_npos[1,])
		outi2<-as.integer(unlist(locii_dams_npos))
		arg7_f<-outi2
		arg4_f<-length(locii_dams_npos[1,])/2
		arg5_f<-(length(locii_dams_npos[,1])*2)*arg4_f

# CROSSOVER
		# MALES
		For_recom_males<-list()
    	 for (chr_counter in 1:nchr){
			nrec<-rpois(length(locii_sires_npos[,1]),chrom_length[chr_counter]/100)
			nrec[nrec==0]=1
			nrec
			matri_rec<-matrix(-1,nrow=length(locii_sires_npos[,1]),ncol=(max(nrec)))
			matri_rec

			for(qw in 1:length(locii_sires_npos[,1])){
			pos_rec<-sample(pos_each_chr[[chr_counter]][,3],nrec[qw])
			pos_rec<-sort(pos_rec)
			# 4444444444
			qq<-match(pos_rec,pos_each_chr[[chr_counter]][,3])
			pos_rec<-qq
			# 444444444444444
			matri_rec[qw,1:length(pos_rec)]<-pos_rec
			}
			matri_rec
			For_recom_males[[chr_counter]]<-matri_rec
		 }

		 for_dim<-sapply(For_recom_males, dim)
		 new_For_recom_males<-matrix(-1,nrow=sum(for_dim[1,]),ncol=max(for_dim[2,]))
		 dim(new_For_recom_males)
		 for_c<-for_dim[2,]
		 for_r<-c(0,cumsum(for_dim[1,]))
		 a1=1
		 b1=2
		for (i in 1:length(For_recom_males)){
		new_For_recom_males[(for_r[a1]+1):for_r[b1],1:for_c[i]]<-as.matrix(For_recom_males[[i]])
		a1=a1+1
		b1=b1+1
		 if(a1==(length(For_recom_males)+1)) break
		}
		dim(new_For_recom_males)
		funi2<-function(vec){
		length(which(vec>=0))
		}
		 L1<-rep(1:for_dim[1,1],nchr)
		 L2<-rep(1:nchr,for_dim[1,])
		 L3<-apply(new_For_recom_males,1,funi2)
		 For_recom_males<-data.frame(L1,L2,L3,new_For_recom_males)
		 dim(For_recom_males)
		 For_recom_males<-as.matrix(For_recom_males)
		dim(For_recom_males)
		max(For_recom_males)
		# For Fortan CROSSOVER
		r_rec<-dim(For_recom_males)[1]
		c_rec<-dim(For_recom_males)[2]
		r_rec
		c_rec
		outi_rec_m<-as.integer(unlist(For_recom_males))
		arg_rec_m<-outi_rec_m
		class(arg_rec_m)
		index_chr1<-cumsum(sumqtl_mrk_chr)
		index_chr2<-c(0,index_chr1)
		arg8_m=as.integer(unlist(index_chr2))

		# FEMALES
		For_recom_females<-list()
		 for (chr_counter in 1:nchr){
			nrec<-rpois(length(locii_dams_npos[,1]),chrom_length[chr_counter]/100)
			nrec[nrec==0]=1
			nrec
			matri_rec<-matrix(-1,nrow=length(locii_dams_npos[,1]),ncol=(max(nrec)))

			for(qw in 1:length(locii_dams_npos[,1])){
			pos_rec<-sample(pos_each_chr[[chr_counter]][,3],nrec[qw])
			pos_rec<-sort(pos_rec)
				# 4444444444
			qq<-match(pos_rec,pos_each_chr[[chr_counter]][,3])
			pos_rec<-qq
			# 444444444444444

			matri_rec[qw,1:length(pos_rec)]<-pos_rec
			}
			matri_rec
			For_recom_females[[chr_counter]]<-matri_rec

		 }

		 for_dim<-sapply(For_recom_females, dim)
		 new_For_recom_females<-matrix(-1,nrow=sum(for_dim[1,]),ncol=max(for_dim[2,]))
		 dim(new_For_recom_females)
		 for_c<-for_dim[2,]
		 for_c
		 for_r<-c(0,cumsum(for_dim[1,]))
		 for_r
		 a1=1
		 b1=2
		for (i in 1:length(For_recom_females)){
		new_For_recom_females[(for_r[a1]+1):for_r[b1],1:for_c[i]]<-as.matrix(For_recom_females[[i]])
		a1=a1+1
		b1=b1+1
		 if(a1==(length(For_recom_females)+1)) break
		}

		dim(new_For_recom_females)
		funi2<-function(vec){
		length(which(vec>=0))
		}

		 L1<-rep(1:for_dim[1,1],nchr)
		 L2<-rep(1:nchr,for_dim[1,])
		 L3<-apply(new_For_recom_females,1,funi2)
		 For_recom_females<-data.frame(L1,L2,L3,new_For_recom_females)
		 dim(For_recom_females)
		 For_recom_females<-as.matrix(For_recom_females)
		dim(For_recom_females)
		# For Fortan CROSSOVER
		r_f<-dim(For_recom_females)[1]
		c_f<-dim(For_recom_females)[2]
		r_f
		c_f
		outi_rec_f<-as.integer(unlist(For_recom_females))
		arg_rec_f<-outi_rec_f
		class(arg_rec_f)


  # Passing to .Fortran
hp_str<-.Fortran("sh",in_r_mu= as.integer(in_r_mu),in_c_mu=as.integer(in_c_mu),in_pois=arg7_m,hp_loci=integer(arg5_m),in_pois_f=arg7_f,hp_loci_f=integer(arg5_f),r_rec=as.integer(r_rec),c_rec=as.integer(c_rec),nchr=as.integer(nchr),szch=arg8_m,rec_m=arg_rec_m,r_f=as.integer(r_f),c_f=as.integer(c_f),rec_f=arg_rec_f
)

 # pos_total<-pos[seq(1,length(pos),2)]
# After Fortran calcs
lmi_m<-matrix(as.numeric(hp_str[[4]]),nrow = (length(locii_sires_npos[,1])*2),ncol = arg4_m)
tempisirem<-lmi_m
# tempisirem[,1:5]

lmi_f<-matrix(as.numeric(hp_str[[6]]),nrow = (length(locii_dams_npos[,1])*2),ncol = arg4_f)
tempidamm<-lmi_f
# dim(tempidamm)
# tempidamm[,1:5]
 # end.time <- Sys.time()
# time.taken <- end.time - start.time
# cat('COMBINED 1to2 and CROSSOVER ',time.taken,fill=TRUE)

cat('.')
# ////////////////////////////////////
# COMBINED 1to2 and CROSSOVER FORTRAN///////FINISH////
# ////////////////////////////////////


 # dams/////////////////////////

#sampling markers for offspring from parents
# cat('sampling markers for offspring from parents',fill=TRUE)
 # start.time <- Sys.time()

# # changed and worked
# offsire<-matrix(ncol=length(tempisirem[1,]),nrow=x)
# offdam<-matrix(ncol=length(tempidamm[1,]),nrow=x)
# offsire<-as.data.frame(offsire)
# offdam<-as.data.frame(offdam)

# -------------------------------------------
dada<-matrix(ncol=length(tempisirem[1,]),nrow=x)
nana<-matrix(ncol=length(tempidamm[1,]),nrow=x)
# dada
vant1<-sample(c(1,2),x,replace=TRUE)
vant2<-c(0:(x-1))
vant3<-vant2*2
vant_1<-vant3+vant1
z1<-vant_1
dada<-tempisirem[z1,]
# nana
vant1<-sample(c(1,2),x,replace=TRUE)
vant2<-c(0:(x-1))
vant3<-vant2*2
vant_1<-vant3+vant1
z1<-vant_1
nana<-tempidamm[z1,]

offsire<-dada
offdam<-nana

# cat(dim(dada),'dim(dada)',fill=TRUE)
# ----------------------------------


# seq1<-seq(1,x*2,2)
# seq2<-seq(2,x*2,2)
# dada<-tempisirem
# nana<-tempidamm
# offsire<-dada
# offdam<-nana

# vant1<-sample(c(1,2),x,replace=TRUE)
# vant2<-c(0:(x-1))
# vant3<-vant2*2
# vant_1<-vant3+vant1
# vant_temp<-vant1
# vant_temp[vant_temp==1]=5
# vant_temp[vant_temp==2]=1
# vant_temp[vant_temp==5]=2
# vant_2<-vant3+vant_temp
# z1<-vant_1
# z3<-vant_2
# offsire[seq1,]<-dada[z1,]
# offsire[seq2,]<-dada[z3,]

# vant1<-sample(c(1,2),x,replace=TRUE)
# vant2<-c(0:(x-1))
# vant3<-vant2*2
# vant_1<-vant3+vant1
# vant_temp<-vant1
# vant_temp[vant_temp==1]=5
# vant_temp[vant_temp==2]=1
# vant_temp[vant_temp==5]=2
# vant_2<-vant3+vant_temp
# z1<-vant_1
# z3<-vant_2
# offdam[seq1,]<-nana[z1,]
# offdam[seq2,]<-nana[z3,]

# dada<-matrix(ncol=length(tempisirem[1,]),nrow=x)
# nana<-matrix(ncol=length(tempidamm[1,]),nrow=x)
# # # changed and worked
 # # dada<-as.data.frame(dada)
 # # nana<-as.data.frame(nana)

# w1<-sample(c(1,2),x,replace=TRUE)
# w2<-sample(c(1,2),x,replace=TRUE)
# # j<-1
# seq1<-seq(1,x,2)
# seq2<-seq(2,x,2)
# dada[seq1,]<-tempisirem[seq1,]
# dada[seq2,]<-tempisirem[seq2,]
# nana[seq1,]<-tempidamm[seq1,]
# nana[seq2,]<-tempidamm[seq2,]
# offsire<-dada
# offdam<-nana
# in_data<-data.frame(seq1,seq2)
# # funi<-function(vec){
# # L<-sample(vec,1)
# # return(L)
# # }
# z1<-apply(in_data,1,funi)
# zi<-seq(1,x,1)
# z2<-match(z1,zi)
# z3<-zi[-z2]
# offsire[seq1,]<-dada[z1,]
# offsire[seq2,]<-dada[z3,]

# z1<-apply(in_data,1,funi)
# zi<-seq(1,x,1)
# z2<-match(z1,zi)
# z3<-zi[-z2]
# offdam[seq1,]<-nana[z1,]
# offdam[seq2,]<-nana[z3,]

offspring<-matrix(ncol=length(tempisirem[1,]),nrow=x*2)
offspring<-as.data.frame(offspring)
seq1<-seq(1,x*2,2)
seq2<-seq(2,x*2,2)
# offspring[seq1,]<-offsire[1:x,]
# offspring[seq2,]<-offdam[1:x,]
offspring[seq1,]<-offsire
offspring[seq2,]<-offdam
offspringasli<-matrix(ncol=length(tempisirem[1,])*2,nrow=x)
offspringasli<-as.data.frame(offspringasli)
seq1<-seq(1,length(tempisirem[1,])*2,2)
seq2<-seq(2,length(tempisirem[1,])*2,2)


	sseq1<-seq(1,x*2,2)
	sseq2<-seq(2,x*2,2)
	offspringasli[1:x,seq1]<-offspring[sseq1,]
	offspringasli[1:x,seq2]<-offspring[sseq2,]
	offspringasli[1:5,1:5]
	offspringasli<-as.matrix(offspringasli)
# whole sequence of animals
sequence_sel<-offspringasli
#Extract qtl_genome for each offspring from total sequence based on positions
qtlsel1 <- matrix(ncol=nqtl_allele,nrow = x)
offspringasli_1<-rbind(pos,offspringasli)
index<-posqtl
index_qtls <-which(offspringasli_1[1,]%in%index)
qtlsel1<-offspringasli_1[,index_qtls]
qtlsel1<-qtlsel1[-1,]
# qtlsel1[1:5,1:5]
qtlsel1<-as.matrix(qtlsel1)

# delete row below if doesnt work
  offspringasli_1<-as.matrix(offspringasli_1)

msel1<- matrix(ncol=nmarker_allele,nrow = x)
index<-posmarker
index_marker <-which(offspringasli_1[1,]%in%index)
msel1<-offspringasli_1[,index_marker]
msel1<-msel1[-1,]
msel<-msel1
# msel[1:5,1:5]
msel<-as.matrix(msel)

# end.time <- Sys.time()
# time.taken <- end.time - start.time
# cat('sampling markers for offspring from parents',time.taken,fill=TRUE)
cat('.')
#QTL allele freq
# freq1qtl_B<-Calc_freq(qtlsel1)

# Passing to .Fortran for calc Freq--start
loci_mat<-qtlsel1
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1qtl_B<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish


freq2qtl_B<-1-freq1qtl_B
# freq_qtl<-data.frame(freq1qtl_B,freq2qtl_B)
# freq_qtl

#Marker allele freq
# freq1mrk_B<-Calc_freq(msel)
# Passing to .Fortran for calc Freq--start
loci_mat<-msel
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1mrk_B<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish

freq2mrk_B<-1-freq1mrk_B
# freq_mrk<-data.frame(freq1mrk_B,freq2mrk_B)
# freq_mrk
# # dim(freq_mrk)


id_B<-id
sire_B <-sire
dam_B<-dam
generation_B<-generation
sex_B<-sex
env_B<-env
msel_B<-msel
qtlsel1_B<-qtlsel1

#breed B starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id<-id_B
sire<-sire_B
dam<-dam_B
generation<-generation_B
sex<-sex_B
env<-env_B
msel<-msel_B
qtlsel1<-qtlsel1_B

#TBV
freq1<-freq1qtl_B
freq2<-freq2qtl_B
snp_validation<-qtlsel1
 s1<-seq(1,nqtl_allele,2)
 s2<-seq(2,nqtl_allele,2)
 a1<-snp_validation[,s1]+snp_validation[,s2]
 ##a1[a1==2]=2
 a1[a1==3]=1
 a1[a1==4]=0
 Snp_BreedB<-a1
xprogeny<-Snp_BreedB
xprogeny<-as.matrix(xprogeny)
q1<-xprogeny
for (i in 1:length(xprogeny[,1])){
ti<-xprogeny[i,]
two<-which(ti==2)
one<-which(ti==1)
zero<-which(ti==0)
q1[i,two]<-((freq1[two])*add_eff_1[two])+
         ((freq2[two])*dom_eff[two])
q1[i,one]<-((1/2*freq1[one])*add_eff_1[one])+
         ((1/2*freq2[one])*dom_eff[one])+
		 ((1/2*freq1[one])*dom_eff[one])+
		 ((1/2*freq2[one])*add_eff_2[one])
q1[i,zero]<-((freq2[zero])*add_eff_2[zero])+
         ((freq1[zero])*dom_eff[zero])
}
xprogeny<-q1
tbvp<- rowSums(xprogeny)
cat('.')

#True genetic value 4
tgv4<-calc_TGV(qtlsel1,add_eff_1,dom_eff)
# mean(tgv4)
phen <- tgv4 + env
# class(phen)
phen<-as.numeric(phen)

# TRAINING START-------------------
if(Selection[1,2]=='gebv' | Selection[2,2]=='gebv'){

	if(gene_counter ==0){
	 # library(BGLR)
	}

# # originak
	# if(tra_sel=='rnd'){
	# index_training<-sample(dim(msel)[1],tra_size,replace=FALSE)
	# # index_training<-sample(dim(msel)[1],200,replace=FALSE)
	# snp_reference<-msel[index_training,]
	# phen_reference<-phen[index_training]
	# }

	if(tra_sel=='rnd'){
	index_training<-sample(dim(msel)[1],tra_size,replace=FALSE)
	# index_training<-sample(dim(msel)[1],200,replace=FALSE)
	snp_reference<-msel[index_training,]
	phen_reference<-phen[index_training]
	} else if (tra_sel=='max_rel_mrk'){
		# MRK based
		# Genomic Relationship Matrix
		M_matrix<-msel
		# freq1_t<-Calc_freq(M_matrix)
# Passing to .Fortran for calc Freq--start
loci_mat<-M_matrix
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1_t<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish

		M_matrix1<-bin_snp(M_matrix)
		M_matrix1<-as.matrix(M_matrix1)
		M<-M_matrix1
		Z<-M
		for (i in 1:length(M[1,])){
		Z[,i]<-M[,i]-(2*(freq1_t[i]))
		}
		sumfreq<-2*sum(freq1_t*(1-freq1_t))
		Zprime<-t(Z)
		G_mat1<-(Z%*%Zprime)
		G_mat<-G_mat1/sumfreq
		ini<-which(G_mat>=min(G_mat), arr.ind = T)
		tn<-G_mat[ini]
		tn2<-cbind(ini,tn)
		tn2<-as.data.frame(tn2)
	# max to min
		tn3<-tn2[order(-tn2[,3]), ]
		y2<-c()
		for (i in 1:length(tn3[,1])){
		y1<-as.numeric(tn3[i,1:2])
		y2<-c(y1,y2)
		if(length(unique(y2))>=tra_size)  break
		}
		index_training<-unique(y2)
		index_training<-sample(index_training,tra_size,replace=FALSE)
		snp_reference<-msel[index_training,]
		phen_reference<-phen[index_training]
	} else if (tra_sel=='min_rel_mrk'){
		# MRK based
		# Genomic Relationship Matrix
		M_matrix<-msel
		dim(M_matrix)
		# freq1_t<-Calc_freq(M_matrix)
# Passing to .Fortran for calc Freq--start
loci_mat<-M_matrix
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1_t<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish
cat('.')


		M_matrix1<-bin_snp(M_matrix)
		M_matrix1<-as.matrix(M_matrix1)
		M<-M_matrix1
		Z<-M
		for (i in 1:length(M[1,])){
		Z[,i]<-M[,i]-(2*(freq1_t[i]))
		}
		sumfreq<-2*sum(freq1_t*(1-freq1_t))
		Zprime<-t(Z)
		G_mat1<-(Z%*%Zprime)
		G_mat<-G_mat1/sumfreq
		ini<-which(G_mat>=min(G_mat), arr.ind = T)
		tn<-G_mat[ini]
		tn2<-cbind(ini,tn)
		tn2<-as.data.frame(tn2)
		# min to max
		tn3<-tn2[order(tn2[,3]), ]
		y2<-c()
		for (i in 1:length(tn3[,1])){
		y1<-as.numeric(tn3[i,1:2])
		y2<-c(y1,y2)
		if(length(unique(y2))>tra_size)  break
		}
		index_training<-unique(y2)
		index_training<-sample(index_training,tra_size,replace=FALSE)
		snp_reference<-msel[index_training,]
		phen_reference<-phen[index_training]
	} else if (tra_sel=='max_rel_qtl'){
		# QTL based
		# Genomic Relationship Matrix
		M_matrix<-qtlsel1
		# freq1_t<-Calc_freq(M_matrix)

# Passing to .Fortran for calc Freq--start
loci_mat<-M_matrix
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1_t<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish

		freq2_t<-1-freq1_t
		M_matrix1<-bin_snp(M_matrix)
		M_matrix1<-as.matrix(M_matrix1)
		M<-M_matrix1
		Z<-M
		for (i in 1:length(M[1,])){
		Z[,i]<-M[,i]-(2*(freq1_t[i]))
		}
		sumfreq<-2*sum(freq1_t*(1-freq1_t))
		Zprime<-t(Z)
		G_mat1<-(Z%*%Zprime)
		G_mat<-G_mat1/sumfreq
		ini<-which(G_mat>=min(G_mat), arr.ind = T)
		tn<-G_mat[ini]
		tn2<-cbind(ini,tn)
		tn2<-as.data.frame(tn2)
	   # max to min
		tn3<-tn2[order(-tn2[,3]), ]
		y2<-c()
		for (i in 1:length(tn3[,1])){
		y1<-as.numeric(tn3[i,1:2])
		y2<-c(y1,y2)
		if(length(unique(y2))>=tra_size)  break
		}
		index_training<-unique(y2)
		index_training<-sample(index_training,tra_size,replace=FALSE)
		snp_reference<-msel[index_training,]
		phen_reference<-phen[index_training]
	} else if (tra_sel=='min_rel_qtl'){
		# QTL based
		# Genomic Relationship Matrix
		M_matrix<-qtlsel1
		dim(M_matrix)
		# freq1_t<-Calc_freq(M_matrix)

# Passing to .Fortran for calc Freq--start
loci_mat<-M_matrix
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1_t<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish

		M_matrix1<-bin_snp(M_matrix)
		M_matrix1<-as.matrix(M_matrix1)
		M<-M_matrix1
		Z<-M
		for (i in 1:length(M[1,])){
		Z[,i]<-M[,i]-(2*(freq1_t[i]))
		}
		sumfreq<-2*sum(freq1_t*(1-freq1_t))
		Zprime<-t(Z)
		G_mat1<-(Z%*%Zprime)
		G_mat<-G_mat1/sumfreq
		ini<-which(G_mat>=min(G_mat), arr.ind = T)
		tn<-G_mat[ini]
		tn2<-cbind(ini,tn)
		tn2<-as.data.frame(tn2)
	   #  min to max
		tn3<-tn2[order(tn2[,3]), ]
		y2<-c()
		for (i in 1:length(tn3[,1])){
		y1<-as.numeric(tn3[i,1:2])
		y2<-c(y1,y2)
		# cat('(length(unique(y2))',length(unique(y2)))
		# if(length(unique(y2))==tra_size)  break
		if(length(unique(y2))>=tra_size)  break
		}
		index_training<-unique(y2)
		index_training<-sample(index_training,tra_size,replace=FALSE)
		snp_reference<-msel[index_training,]
		phen_reference<-phen[index_training]
	} #'min_rel_qtl'

snp_reference<-bin_snp(snp_reference)
snp_reference<-as.matrix(snp_reference)
y<-as.matrix(phen_reference)


if(addTra==TRUE){
	X1<-snp_reference
	X1[X1==1]=0
	X1[X1==2]=1

	X3<-snp_reference
	X3[X3==0]=5
	X3[X3==1]=0
	X3[X3==2]=0
	X3[X3==5]=1

	ETA<-list(
	list(X=X1, model=tra_method),
	list(X=X3, model=tra_method)
	)

	answers<-BGLR(y=y,ETA=ETA, nIter=tra_nIter, burnIn=tra_burnIn,thin=tra_thin,rmExistingFiles=TRUE,verbose = tra_show,saveAt=tra_save)

	G11<-answers$ETA[[1]]$b
	G22<-answers$ETA[[2]]$b

	ghat1<-G11
	dhat<-rep(0,length(ghat1))
	ghat2<-G22
	ghat<-cbind(ghat1,ghat2)

} else if (addTra==FALSE){

	X1<-snp_reference
	X1[X1==1]=0
	X1[X1==2]=1
	X2<-snp_reference
	X2[X2==2]=0
	X3<-snp_reference
	X3[X3==0]=5
	X3[X3==1]=0
	X3[X3==2]=0
	X3[X3==5]=1

	ETA<-list(
	list(X=X1, model=tra_method),
	list(X=X2, model=tra_method),
	list(X=X3, model=tra_method)
	)

	answers<-BGLR(y=y,ETA=ETA, nIter=tra_nIter, burnIn=tra_burnIn,thin=tra_thin,rmExistingFiles=TRUE,verbose = tra_show,saveAt=tra_save)

	G11<-answers$ETA[[1]]$b
	G12<-answers$ETA[[2]]$b
	G22<-answers$ETA[[3]]$b

	ghat1<-G11
	dhat<-G12
	ghat2<-G22
	ghat<-cbind(ghat1,ghat2)

}




# write.table(ghat,file=fileghatBsc1[ii],row.names=F,col.names=F)
# write.table(dhat,file=filedhatBsc1[ii],row.names=F,col.names=F)
# write.table(ghat,file=fileghatAsc1[ii],row.names=F,col.names=F)
# write.table(dhat,file=filedhatAsc1[ii],row.names=F,col.names=F)
# TRAINING FINISH-------------------
}

cat('.')
#Gebv
if(Selection[1,2]=='gebv' | Selection[2,2]=='gebv'){
freq1<-freq1mrk_B
freq2<-freq2mrk_B
snp_validation<-msel
s1<-seq(1,nmarker_allele,2)
s2<-seq(2,nmarker_allele,2)
a1<-snp_validation[,s1]+snp_validation[,s2]
a1[a1==3]=1
a1[a1==4]=0
# if (gene_counter==1){
# ghat<-read.table(file=fileghatBsc1[ii])
# dhat<-read.table(file=filedhatBsc1[ii])
# }
# if (gene_counter>1){
# ghat<-read.table(file=fileghat_re_B[gene_counter])
# dhat<-read.table(file=filedhat_re_B[gene_counter])
# }
# ghat<-as.matrix(ghat)
dhat<-as.matrix(dhat)
ghat1<-ghat[,1]
ghat2<-ghat[,2]

snp_validation<-a1
xprogeny<-snp_validation
xprogeny<-as.matrix(xprogeny)
q1<-xprogeny
for (i in 1:length(xprogeny[,1])){
ti<-xprogeny[i,]
two<-which(ti==2)
one<-which(ti==1)
zero<-which(ti==0)
q1[i,two]<-((freq1[two])*ghat1[two])+
         ((freq2[two])*dhat[two])
q1[i,one]<-((1/2*freq1[one])*ghat1[one])+
         ((1/2*freq2[one])*dhat[one])+
		 ((1/2*freq1[one])*dhat[one])+
		 ((1/2*freq2[one])*ghat2[one])
q1[i,zero]<-((freq2[zero])*ghat2[zero])+
         ((freq1[zero])*dhat[zero])
}
xprogeny<-q1
gebvp<- rowSums(xprogeny)
# cor(gebvp,tbvp)
} else {
gebvp<-rep(0,length(tbvp))
}


qtlsel<-qtlsel1
newgen <- data.frame(id, sire, dam, generation,sex, phen,env,tbvp,gebvp)
# newgen[1:10,1:9]

 start.time <- Sys.time()

   # Total_Total[[gene_counter+1]][[1]]<-total_B$data
   # Total_Total[[gene_counter+1]][[2]]<-total_B$qtl
   # Total_Total[[gene_counter+1]][[3]]<-total_B$mrk
   # Total_Total[[gene_counter+1]][[4]]<-total_B$sequ
   # # Freq QTL
   	# ID_qtl<-my_internal_data1[[5]][,1]
	# ID_qtl_gen<-rep(gene_counter,length(ID_qtl))
	# ID_qtl_chr<-my_internal_data1[[5]][,3]
	# freqQTL<-data.frame(ID_qtl,ID_qtl_gen,ID_qtl_chr,freq1qtl_B,freq2qtl_B)
	# names(freqQTL)<-c('ID','Generation','Chr','Freq.Allele1','Freq.Allele2')
   # Total_Total[[gene_counter+1]][[5]]<-freqQTL

   # # Freq MRK
   # # outforLD$freqMrk
   	# ID_mrk<-my_internal_data1[[6]][,1]
	# ID_mrk_gen<-rep(gene_counter,length(ID_mrk))
	# ID_mrk_chr<-my_internal_data1[[6]][,3]
	# freqMRK<-data.frame(ID_mrk,ID_mrk_gen,ID_mrk_chr,freq1mrk_B,freq2mrk_B)
	# names(freqMRK)<-c('ID','Generation','Chr','Freq.Allele1','Freq.Allele2')
   # Total_Total[[gene_counter+1]][[6]]<-freqMRK


colnames(newgen) <- colnames(total_B$data)
# DATA
total_B$data <-newgen

# QTL
qtlsel_finall<-cbind(id,generation,qtlsel)
total_B$qtl<-qtlsel_finall

# Marker
msel_finall<-cbind(id,generation,msel)
total_B$mrk<-msel_finall

# SEQ
sequence_sel_finall<-cbind(id,generation,sequence_sel)
sequence_sel_finall<-as.matrix(sequence_sel_finall)
total_B$sequ<-sequence_sel_finall

cat('.','\n')
 # cat('Generation',gene_counter,'is finished',fill=TRUE)
 # end.time <- Sys.time()
# time.taken <- end.time - start.time
# cat('*******Time for Rbindin',time.taken,fill=TRUE)

end_ge <- Sys.time()
time_gen <- end_ge - start_ge
cat('Generation',gene_counter,'is finished.','Time taken:',time_gen,fill=TRUE)


} # end of do loop for generations



# ---------------------START-----------------
# DATA RETURN BACK BY FUNCTION
# ---------------------START-----------------


#Fill in data of last hp
cat('Output data preparation ...',fill=TRUE)

# AFTER GENERATIONS FINISHED FOR SUMMARY
    for_summary<-data.frame()
	for (i in 1:length(Total_Total)){
	ali<-Total_Total[[i]][[1]]
	for_summary<-rbind(for_summary,ali)
	}

mean_pheno<-as.numeric()
mean_tbv<-as.numeric()
mean_gebv<-as.numeric()
var_tbv<-as.numeric()
var_pheno<-as.numeric()
h2p<-as.numeric()
accuracy<-as.numeric()
M_accuracy<-as.numeric()
F_accuracy<-as.numeric()

for (i in 1:ng){
a<-subset(for_summary,for_summary[,4]==i)
mean_pheno[i]<-mean(a[,6])
var_pheno[i]<-var(a[,6])
mean_tbv[i]<-mean(a[,8])
var_tbv[i]<-var(a[,8])
mean_gebv[i]<-mean(a[,9])
# ACCURACY
	# Males
	b2<-subset(a,a[,5]=='M')
	if(Selection[1,2]=='rnd'){
		M_accuracy[i]<-cor(b2[,8],b2[,6])
		} else if (Selection[1,2]=='phen') {
		M_accuracy[i]<-cor(b2[,8],b2[,6])
		} else if (Selection[1,2]=='tbv') {
		M_accuracy[i]<-cor(b2[,8],b2[,8])
		} else if (Selection[1,2]=='gebv') {
		M_accuracy[i]<-cor(b2[,8],b2[,9])
	}
	# Females
	b2<-subset(a,a[,5]=='F')
	if(Selection[2,2]=='rnd'){
		F_accuracy[i]<-cor(b2[,8],b2[,6])
		} else if (Selection[2,2]=='phen') {
		F_accuracy[i]<-cor(b2[,8],b2[,6])
		} else if (Selection[2,2]=='tbv') {
		F_accuracy[i]<-cor(b2[,8],b2[,8])
		} else if (Selection[2,2]=='gebv') {
		F_accuracy[i]<-cor(b2[,8],b2[,9])
	}
h2p[i]<-var(a[,8])/var(a[,6])
}


Generation<-1:ng
Phenotype<-mean_pheno
TrueBV<-mean_tbv
GEBV<-mean_gebv
M_accuracy<-M_accuracy
F_accuracy<-F_accuracy
heritability<-h2p
heritability<-heritability*4


if(Selection[1,2]=='gebv' | Selection[2,2]=='gebv'){
summary_data<-data.frame(Generation,Phenotype,TrueBV,GEBV,M_accuracy,F_accuracy,heritability)
} else {
summary_data<-data.frame(Generation,Phenotype,TrueBV,M_accuracy,F_accuracy,heritability)
}

if(Display==TRUE){
print(summary_data)
}
internal_map_QTL<-outforLD$linkage_map_qtl
internal_map_MRK<-outforLD$linkage_map_mrk
internal_map_qtl_mrk<-outforLD$linkage_map_qtl_mrk
allele_effcts<-outforLD$allele_effcts
trait<-outforLD$trait
genome<-outforLD$genome

# ---------------------FINISH-----------------
# DATA RETURN BACK BY FUNCTION
# ---------------------FINISH-----------------

#WRITE TO OUTPUT
	if(!missing(saveAt) & !missing(rp_output)){
	cat('Writing output files ...',fill=TRUE)

control_names<-c("data","qtl","marker","seq","freq_qtl","freq_mrk")
test_w<-intersect(names(rp_output),control_names)

		# data to output
	if(any(test_w=='data')){
		in_w<-rp_output$data
		in_w2<-in_w+1
		for(P1 in 1:length(in_w)){
		outFile_data<-paste(saveAt,'_data_',sep='')
		outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
		dom<-format(Total_Total[[in_w2[P1]]][[1]],  justify = "right")
		write.table(dom,file=outFile_data,row.names=FALSE,col.names=TRUE,quote = FALSE)
		}
	}

	    # qtl to output
	if(any(test_w=='qtl')){
    in_w<-rp_output$qtl
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(saveAt,'_qtl_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-Total_Total[[in_w2[P1]]][[2]]

	writeLines(c('ID  Generaion  Genotypes (paternal allele, maternal allele) ...'),outFile_data)
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)
	}
    }

		# mrk to output
	if(any(test_w=='marker')){
    in_w<-rp_output$marker
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(saveAt,'_mrk_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-Total_Total[[in_w2[P1]]][[3]]

	writeLines(c('ID  Generaion  Genotypes (paternal allele, maternal allele) ...'),outFile_data)
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)

	}
	}

    # seq to output
	if(any(test_w=='seq')){
    in_w<-rp_output$seq
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(saveAt,'_seq_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-Total_Total[[in_w2[P1]]][[4]]
	writeLines(c('ID  Generaion  Genotypes (paternal allele, maternal allele) ...'),outFile_data)
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)

	}
	}

	    # freq qtl to output
	if(any(test_w=='freq_qtl')){
    in_w<-rp_output$freq_qtl
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(saveAt,'_freq_qtl_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-format(Total_Total[[in_w2[P1]]][[5]],  justify = "right")
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=TRUE,quote = FALSE)
	}
	}

		    # freq qtl to output
	if(any(test_w=='freq_mrk')){
    in_w<-rp_output$freq_mrk
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(saveAt,'_freq_mrk_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-format(Total_Total[[in_w2[P1]]][[6]],  justify = "right")
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=TRUE,quote = FALSE)
	}
	}


} #end loop for writing

userin<-list()
userin[[1]]<-Male_founders
userin[[2]]<-Female_founders
userin[[3]]<-ng
userin[[4]]<-litter_size
userin[[5]]<-Selection

# detach("package:BGLR")
# loadNamespace('BGLR')

for (g_index in 1:length(Total_Total)){
names(Total_Total[[g_index]])<-c('data','qtl','mrk','sequ','freqQTL','freqMRK')
}

cat('Making RP is done!',fill=TRUE)

# RETURN SECTION
Final<-list(output=Total_Total,summary_data=summary_data,linkage_map_qtl=internal_map_QTL,linkage_map_mrk=internal_map_MRK,linkage_map_qtl_mrk=internal_map_qtl_mrk,allele_effcts=allele_effcts,trait=trait,genome=genome,user_input=userin)
return(Final)

} # end of function
