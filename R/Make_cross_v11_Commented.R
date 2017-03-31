#' Create crossbred population
#'
#'  This function can be used for crossing between populations. Different crossbreeding schemes such as two-way, three-way and four-way crossbreeding schemes can be simulated.

################################################
#################PARAMETERS#####################
#################################################
#' @param pop1 (\code{list}) Output of function   \code{\link{make_rp}}, \code{\link{sample_hp}} or \code{\link{xbreed}}.   

#' @param pop2 (\code{list}) Output of function   \code{\link{make_rp}}, \code{\link{sample_hp}} or \code{\link{xbreed}}.  


# -------------founder_pop1 --------------- 
#' @param founder_pop1  (\code{data.frame}) Data frame with {2} rows and {3} columns. First row is for the selection of male founders and second row is for the selection of female founders from \code{pop1}. The columns are as following:\cr
  
#'   Colomn 1) "size" is the number of individuals to be selected. \cr

#'   Column 2) "generation" is the generation number from which males/females will be selected.\cr

#'   Column 3) "select" indicates the type of selection with options:
	#' \itemize{
	#'\item{"rnd"}  {Select individuals randomly}.
	#'\item{"phen"}  {Select individuals based on their phenotypes}.
	#'\item{"tbv"}  {Select individuals based on their true breeding value (tbv)}.
	#'\item{"gebv"}  {Select individuals based on their genomic estimated breeding value (gebv)}.
	#' }
	
#'   Column 4) "value" Indicates to select high: "h" or low: "l" values. Note: This colomn is ignored if individuals are selected randomly.\cr

# -------------founder_pop2 --------------- 
#' @param founder_pop2  (\code{data.frame}) Data frame with {2} rows and {3} columns. First row is for the selection of male founders and second row is for the selection of female founders from \code{pop2}. Details are similar to argument \code{founder_pop1}.
 
# -------------founder_cross --------------- 
#' @param founder_cross  (\code{data.frame}) Data frame with {2} rows and {4} columns that contains information about creating base generation of crossbreds.  First row is for the selection of male founders and second row is for the selection of female founders in order to create crossbred population. The columns are as following:\cr
#'   Column 1) "pop" Indication of the population name from which males/females for crossbreeding will be selected. The two possible options are:
	#' \itemize{
	#'\item{"pop1"}  {Select individuals are from pop1}.
	#'\item{"pop2"}  {Select individuals are from pop2}.
	#' }

#'   Column 2) "size" is the number of individuals to be selected. \cr

#'   Column 3) "select" indicates the type of selection with options:
	#' \itemize{
	#'\item{"rnd"}  {Select individuals randomly}.
	#'\item{"phen"}  {Select individuals based on their phenotypes}.
	#'\item{"tbv"}  {Select individuals based on their true breeding value (tbv)}.
	#'\item{"gebv"}  {Select individuals based on their genomic estimated breeding value (gebv)}.
	#' }
	
#'   Column 4) "value" Indicates to select high: "h" or low: "l" values. This column is ignored if individuals are selected randomly.\cr

#' \bold{Note:} After selecting founders for pop1 and pop2, by arguments \code{founder_pop1} and \code{founder_pop2}, user can select individuals from these founders as parents of crossbreds for the base generation of crossbreds. \cr

# -------------Selection_pop1 --------------- 
#' @param Selection_pop1 (\code{data.frame}) Selection design for pop1 (Breed {1}). Data frame with {2} rows and {3} columns. First row is for the selection design of males and second row is for the selection design of females. The columns are as following:\cr  
#'   Column 1) "Number" is the number of individuals to be selected as sires/dams. \cr
#'   Column 2) "type" indicates the type of selection with options:
	#' \itemize{
	#'\item{"rnd"}  {Select individuals randomly}.
	#'\item{"phen"}  {Select individuals based on their phenotypes}.
	#'\item{"tbv"}  {Select individuals based on their true breeding value (tbv)}.
	#'\item{"tbvc"}  {Select individuals based on true breeding value for crossbred performance(tbvc)}.
	#'\item{"gebv"}  {Select individuals based on their genomic estimated breeding value (gebv)}.
	#'\item{"gebvc"}  {Select individuals based on genomic estimated breeding value for crossbred performance (gebvc)}.
	#'}
	
#'   Column 3) "value" Indicates to select high: "h" or low: "l" values. Note: This colomn is ignored if individuals are selected randomly.\cr


# -------------Selection_pop2 --------------- 
#' @param Selection_pop2 (\code{data.frame}) Selection design for pop2 (Breed {2}). Details are similar to argument \code{Selection_pop1}.

 # -------------Cross_design --------------- 
#' @param Cross_design (\code{data.frame}) Data frame containing information on how to select individuals from pop1 and pop2 as parents of crossbreds over generations. This argument is a data frame with {2} rows and {4} columns. First row is for the selection of males as sires of crossbreds and second row is for the selection of females as dams of crossbreds. The columns are as following:\cr  

#'   Column 1) "pop" Indication of the population name from which males/females for crossing will be selected. The two possible options are:
	#' \itemize{
	#'\item{"pop1"}  {Selected individuals are from pop1}. 
	#'\item{"pop2"}  {Selected individuals are from pop2}.
	#' }
#'	As an example if row {1} and column {1} of argument \code{Cross_design} is equal to "pop2", this means that sires of crossbred are from pop2.\cr
#'   Column 2) "size" is the number of individuals to be selected as sires/dams. \cr
#'   Column 3) "select" indicates the type of selection with options:
	#' \itemize{
	#'\item{"rnd"}  {Select individuals randomly}.
	#'\item{"phen"}  {Select individuals based on their phenotypes}.
	#'\item{"tbv"}  {Select individuals based on their true breeding value (tbv)}.
	#'\item{"tbvc"}  {Select individuals based on true breeding value for crossbred performance(tbvc)}.
	#'\item{"gebv"}  {Select individuals based on their genomic estimated breeding value (gebv)}.
	#'\item{"gebvc"}  {Select individuals based on genomic estimated breeding value for crossbred performance (gebvc)}.
	#'}
	
#'   Column 4) "value" Indicates to select high: "h" or low: "l" values. Note: This column is ignored if individuals are selected randomly.\cr
 # -------------Cross_design finished--------------- 
 
#' @param ng  Number of generations. Range: \eqn{1 \leq \code{ng} \leq 200}.

#' @param litter_size  Litter size or the number of progeny per dam. Range: \eqn{1 \leq \code{x} \leq 200}.

 # -------------train_type --------------- 
#' @param train_type \emph{Optional} (\code{character}) Type of training for the estimation of marker effects. The two possible options are:
	#' \itemize{
	#'\item{"purebred"} {Training will be on purebreds. This means that training for each population is done separately. So, there will be two reference population; one for pop1 (Breed 1) and one for pop2 (Breed 2)}.
	#'\item{"crossbred"} {Training will be on crossbreds. So, estimated marker effects will be the same for both pop1 and pop2}.
	#'}
#'	\bold{Note:} If selection criteria for any population is defined as gebv/gebvc, argument \code{train_type} should be defined.\cr

 # -------------training pop 1 start---------------
#' @param train_pop1 \emph{Optional} (\code{data.frame}) Data frame with {1} row and {8} columns. The columns are as following:\cr  
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
	#'\item{"BayesB"}  {two component mixture prior with a point of mass at zero and a scaled-t slab}.
	#'\item{"BayesC"}  {two component mixture prior with a point of mass at zero and a Gaussian slab}.
	#'}
	#'	Default: "BRR" \cr

#'  Column 4) "nIter" \emph{Optional} The number of iterations. Default: \eqn{1500} \cr
#'   Column 5) "burnIn" \emph{Optional} The number of burn-in. Default: \eqn{500} \cr
#'   Column 6) "thin" \emph{Optional} The number of thinning. Default: \eqn{5} \cr
#'   Column 7) "save" \emph{Optional} This may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs.  Default:"Out_BGLR" \cr 
#'   Column 8) "show" \emph{Optional} (\code{Logical}) if TRUE the iteration history is printed. Default: \code{TRUE}. \cr
#' \bold{Note:} This argument is compulsory if argument \code{train_type} is "purbred". \cr
#' More details about the argument can be found in package \pkg{BGLR}.


 # -------------train_pop2---------------

#' @param train_pop2 \emph{Optional} (\code{data.frame}) Data frame with {1} row and {8} columns similar to argument \code{train_pop1}.  \cr
#' \bold{Note:} This argument is compulsory if argument \code{train_type} is "purbred". \cr

# -------------training finished---------------		

#' @param train_cross \emph{Optional} (\code{data.frame}) Data frame with {1} row and {8} columns similar to argument \code{train_pop1}.  \cr
#' \bold{Note:} This argument is compulsory if argument \code{train_type} is "crossbred".\cr
# -------------training finished---------------		


#' @param saveAt \emph{Optional} (\code{character}). Name to be used to save output files. 

 # -------------output_pop1 starts---------------
#' @param output_pop1 \emph{Optional} (\code{data.frame}). Data frame to specify generation indexs and type of data to be written to output files for "pop1". User can define which type of data and which generations to be written to output files. The possible options are:\cr  
#'\cr
#'   "data"      Individuals data except their genotypes. \cr
#'   "qtl"       QTL genotye of individuals coded as {11,12,21,22}. \cr
#'   "marker"    Marker genotye of individuals. \cr
#'   "seq"       Genotype (both marker (SNP) and QTL) of individuals.\cr
#'   "freq_qtl"  QTL allele frequency. \cr
#'   "freq_mrk"  Marker allele frequency. \cr
#' \bold{Note:} Both arguments \code{output_pop1} and \code{saveAt} should present in the function in order to write the output files for "pop1".\cr


 # -------------output_pop2 starts---------------
#' @param output_pop2 \emph{Optional} (\code{data.frame}). Data frame to specify generations indexs and type of data to be written to output files for pop2. Details are similar to argument \code{output_pop1}.\cr
#' \bold{Note:} Both arguments \code{output_pop2} and \code{saveAt} should present in the function in order to write the output files for "pop2".\cr


 # -------------output_cross starts---------------
#' @param output_cross \emph{Optional} (\code{data.frame}). Data frame to specify generations indexs and type of data to be written to output files for crossbred population. Details are similar to argument \code{output_pop1}.\cr
#' \bold{Note:} Both arguments \code{output_cross} and \code{saveAt} should present in the function in order to write the output files for crossbred population.

#' @param Display \emph{Optional} (\code{Logical}) Display summary of the simulated generations if is not \code{FALSE}. Default: \code{TRUE}. 


################################################
#################RETURN/KEY/EXPORT##############
################################################

# RETURN SECTION	

#' @return \code{list} with all data of simulated populations.\cr
	#' \describe{
	#'\item{$pop1}{(\code{list}) Two-level list  (\code{$pop1[[]][[]]}) containing information about simulated generations. First index (x) indicates generation number. It should be noted that as data for base generation (0) is also stored by the function, to retrive data for a specific generation, index should be equal to generation number plus one. As an example to observe data for generation 2 index should be 3 i.e, \code{$pop1[[3]]$data}. Second index (y) that ranges from {1} to {6} contain the information as following: 
	#' \itemize{
		#'\item{\code{$pop1[[x]]$data}}  {Individuals data except their genotypes. Here x is the generation index}.
		#'\item{\code{$pop1[[x]]$qtl}}  {QTL genotye of individuals}.
		#'\item{\code{$pop1[[x]]$mrk}}   {Marker genotye of individuals}.
		#'\item{\code{$pop1[[x]]$sequ}}  {Genotype (both marker (SNP) and QTL) of individuals}.
		#'\item{\code{$pop1[[x]]$freqQTL}}   {QTL allele frequency}.
		#'\item{\code{$pop1[[x]]$freqMRK}}   {Marker allele frequency}.
		#'}
	#'}

	#'\item{$pop2}{Similar to $pop1}
	
	
	#'\item{$cross}{(\code{list}) Two-level list (\code{$cross$output[[x]][[y]]}) containing information about simulated crossbred generations. Details are similar to \code{$pop1} such as:	
		#' \itemize{
		#'\item{\code{$cross$output[[x]]$data}}  {Crossbred individuals data except their genotypes. Here x is the generation index}.
		#'\item{\code{$cross$output[[x]]$qtl}}  {QTL genotye of crossbred individuals}.
		#'\item{\code{$cross$output[[x]]$mrk}}  {Marker genotye of individuals}.
		#'\item{\code{$cross$output[[x]]$sequ}}  {Genotype (both marker (SNP) and QTL) of individuals}.
		#'\item{\code{$cross$output[[x]]$freqQTL}}  {QTL allele frequency}.
		#'\item{\code{$cross$output[[x]]$freqMRK}}  {Marker allele frequency}.
		#'}
	#'}
	#'\item{$cross$summary_data_cross}{Data frame with summary of simulated crossbreds}. 	
	#'\item{$summary_data_pop1}{Data frame with summary of simulated generations for pop1 (Breed 1)}. 
	#'\item{$summary_data_pop2}{Data frame with summary of simulated generations for pop2 (Breed 2)}. 
	#'\item{$linkage_map_qtl}{Linkage map for qtl}.  
	#'\item{$linkage_map_mrk}{Linkage map for marker}.
	#'\item{$linkage_map_qtl_mrk}{Integrated linkage map for both marker and qtl}.
	#'\item{$allele_effcts}{QTL allele effects}.
	#'\item{$trait}{Trait specifications}.
	#'\item{$genome}{Genome specifications}.
	#'}
 
 
#' @export xbreed
#' @seealso \code{\link{make_rp}}, \code{\link{sample_hp}}

#' @references {Esfandyari H., A.C Sorensen and P. Bijma. 2015. A crossbred reference population can improve the response to genomic selection for crossbred performance. \emph{Genetics Selection Evolution} \bold{47}: 76}

#' @references {Esfandyari H., A.C Sorensen and P. Bijma. 2015. Maximizing crossbred performance through purebred genomic selection. \emph{Genetics Selection Evolution} \bold{47}: 16}



################################################
#################EXAMPLS########################
################################################
#' @examples 
#'\dontrun{

#'# # Simulation of a two-way crossbreeding program.
#' # The crossbreeding scheme in this example involves three steps:    
#'#Step 1: Historical population is created.      
#'#Step 2: Two recent populations named as Breed A and B
#'# are created by sampling individuals from historical population.       
#'#Step 3: Breed A and B are crossed.  
#'
#' #--------------------------------------
#'# # STEP 1: CREATE HISTORICAL POPULATION
#'
#'	# Genome consisted of 3 chromosomes
#'	genome<-data.frame(matrix(NA, nrow=3, ncol=6))
#'	names(genome)<-c("chr","len","nmrk","mpos","nqtl","qpos")
#'	genome$chr<-c(1:3)
#'	genome$len<-rep(100,3)	
#'	genome$nmrk<-rep(100,3)
#'	genome$mpos<-rep('rnd',3)	
#'	genome$nqtl<-rep(25,3)
#'	genome$qpos<-rep('rnd',3)	
#'	genome
#'
#'	historical<-make_hp(hpsize=300
#'	,ng=10,h2=0.25,d2=0.10,phen_var=1
#'	,genome=genome,mutr=5*10**-4,laf=0.5)
#'
#'# # STEP 2: MAKE BREED A AND B 
#'	# BREED A
#'	Breed_A_Male_fndrs<-data.frame(number=50,select='rnd') 
#'	Breed_A_Female_fndrs<-data.frame(number=50,select='rnd') 
#'
#'	# Selection and matings in Breed A
#'		# Selection of 50 sires and 50 dam			
#'		# Selection criteria is "rnd" for both sires and dams
#'	 Selection<-data.frame(matrix(NA, nrow=2, ncol=2))
#'	 names(Selection)<-c("Number","type")
#'	 Selection$Number[1:2]<-c(50,50)	
#'	 Selection$type[1:2]<-c("rnd","rnd")	
#'	 Selection
#'
#'	Breed_A<-sample_hp(hp_out=historical,Male_founders=
#'	Breed_A_Male_fndrs,Female_founders=Breed_A_Female_fndrs,
#'	ng=5,Selection=Selection,
#'	litter_size=3,Display=TRUE)
#'	 
#'	# BREED B
#'	Breed_B_Male_fndrs<-data.frame(number=50,select="rnd") 
#'	Breed_B_Female_fndrs<-data.frame(number=50,select="rnd") 
#'
#'	# Selection and matings in Breed B
#'		# Selection of 50 sires and 50 dam			
#'		# Selection criteria is "phen" for both sires and dams
#'
#'	Selection<-data.frame(matrix(NA, nrow=2, ncol=3))
#'	names(Selection)<-c("Number","type","Value")
#'	Selection$Number[1:2]<-c(50,50)	
#'	Selection$type[1:2]<-c("phen","phen")	
#'	Selection$Value[1:2]<-c("h","h") 
#'	Selection
#'
#'	Breed_B<-sample_hp(hp_out=historical,Male_founders=
#'	Breed_B_Male_fndrs,Female_founders=Breed_B_Female_fndrs,
#'	ng=5,Selection=Selection,
#'	litter_size=3,Display=TRUE)
#'
#'# # STEP 3: CROSSING BETWEEN BREED A AND B
#'
 #'   # Selection of founders in crossbreeding for Breed A 
 #' 		# Selection of 50 sires and 50 dams 
#'          # from last generation of pop in step 2. 			
#'		# Selection criteria is "rnd" for both sires and dams
#'	founder_pop1<-data.frame(matrix(NA, nrow=2, ncol=3))
#'	names(founder_pop1)<-c("size","generation","select")
#'	founder_pop1[1,]<-c(50,5,"rnd") 
#'	founder_pop1[2,]<-c(50,5,"rnd")
#'	founder_pop1
#'	  
#'	# Selection of founders in crossbreeding for Breed B 
#'		# Selection of 40 sires and 40 dams			
#'		# Selection criteria is "phen" for sires 
#'		# Selection criteria is "rnd" for dams
#'
#'	founder_pop2<-data.frame(matrix(NA, nrow=2, ncol=4))
#'	names(founder_pop2)<-c("size","generation","select","value")
#'	founder_pop2[1,]<-c(40,5,"phen","h") 
#'	founder_pop2[2,]<-c(40,5,"rnd","h") # "h" will be ignored as SC is "rnd"
#'	founder_pop2
#'
#'	# Selection of animals from founder_pop1 and founder_pop2 to be crossed
#'	founder_cross<-data.frame(matrix(NA, nrow=2, ncol=4))
#'	names(founder_cross)<-c("pop","size","select","value")
#'	founder_cross[1,]<-c("pop1",35,"tbv","h") # Select males from Breed A
#'	founder_cross[2,]<-c("pop2",40,"phen","h")  # Select females from Breed B
#'	founder_cross
#'
#'	# Selection scheme in Breed A to produce purebred replacement animals 
#'	Selection_pop1<-data.frame(matrix(NA, nrow=2, ncol=3))
#'	names(Selection_pop1)<-c("Number","type","Value")
#'	Selection_pop1$Number[1:2]<-c(70,70)				
#'	Selection_pop1$type[1:2]<-c("tbv","tbv")	
#'	Selection_pop1$Value[1:2]<-c("h","h")
#'	Selection_pop1
	#'
#'# Selection scheme in Breed B to produce purebred replacement animals 
#'Selection_pop2<-data.frame(matrix(NA, nrow=2, ncol=3))
#'names(Selection_pop2)<-c("Number","type","Value")
#'Selection_pop2$Number[1:2]<-c(40,40)				
#'	Selection_pop2$type[1:2]<-c("phen","phen")	
#'	Selection_pop2$Value[1:2]<-c("h","h")
#'	Selection_pop2
	#'
#'	# Selection scheme for crossing between A and B
#'	Cross_design<-data.frame(matrix(NA, nrow=2, ncol=4))
#'	names(Cross_design)<-c("pop","size","select","value")
#'	Cross_design[1,]<-c("pop1",50,"phen","h") 
#'	Cross_design[2,]<-c("pop2",100,"phen","h")
#'	Cross_design
#'
		 
#'	# Save data for crossbred AB
#'	output_cross<-data.frame(matrix(NA, nrow=1, ncol=1))
#"	names(output_cross)<-c("data")
#'	output_cross[,1]<-c(1)
#'	output_cross
	
#'
#'cross_AB<-xbreed(pop1=Breed_A,pop2=Breed_B,founder_pop1=
#'founder_pop1,founder_pop2=founder_pop2,
#'founder_cross=founder_cross,
#'Selection_pop1=Selection_pop1,Selection_pop2=Selection_pop2,
#'Cross_design=Cross_design,ng=2,litter_size=4,
#'saveAt="cross_pop",output_cross=output_cross,Display=TRUE)
#'}

################################################
#################DETAILS########################
################################################

#' @details	
#' Function \code{xbreed} is used for crossing between populations. These populations can be the ones created by functions \code{sample_hp} and \code{make_rp}. Also, if user would like to have multi-way crossbreeding schemes such as three-way crossbreeding, then function \code{xbreed} can be used to get the output of itself as input data in order to create the multi-way crossbreed populations. Simulations of two-way and multi-way crossbreeding schemes are presented in the package vignette.


# keep wrapping for later on
xbreed<-function(pop1,pop2,founder_pop1,founder_pop2,founder_cross,Selection_pop1,Selection_pop2,Cross_design,ng,litter_size,train_type,train_pop1,train_pop2,train_cross,saveAt,output_pop1,output_pop2,output_cross,Display) {

# dyn.load('sh.dll')
# is.loaded("sh")
# library(BGLR)
# dyn.load('cf.dll')
# is.loaded("cf")

	 
	 # NOTE
	 # ma founder_cross ra moshakhas kardim be in dalil ke bad az entekhab founder_pop1/founder_pop2 user betavanad az beyne inha entekhab konad ke kodam baham talaghi dahad. in option mahsob mishavad vali mitavand gij konande bashad.
	 # rahe digar in hast ke heyvanate entekhab shode ba founder_pop1/founder_pop1 be onvane founder baraye cross 0 estefade shavand va argument founder_cross
	 # ra delete konim. ya inke bemanad va optional bashad.
	 # bayesti rajebesh fekr konam
# ---------------------START-----------------
# INPUT by user control section
# ---------------------START-----------------
cat('Controlling input data ...',fill=TRUE)

# pop1
	if(missing(pop1)) {
	stop('--- argument "pop1" is missing')
	}

	if(class(pop1)!='list') {
	stop('--- argument "pop1" should be the output of function Sample_hp/Make_rp')
	}
	
# pop2
	if(missing(pop2)) {
	stop('--- argument "pop2" is missing')
	}

	if(class(pop2)!='list') {
	stop('--- argument "pop2" should be the output of function Sample_hp/Make_rp')
	}
	
# founder_pop1
	if(missing(founder_pop1)) {
	stop('--- argument "founder_pop1" is missing')
	}
	
	founder_pop1[,1]<-as.numeric(founder_pop1[,1])
	founder_pop1[,2]<-as.numeric(founder_pop1[,2])
	if(	length(which(founder_pop1[,1]<1))>0 | any(floor(founder_pop1[,1])!=founder_pop1[,1]))  {
	stop('\n','--- Number of selected males/females in argument "founder_pop1" should be an integer grater than 1')
	}

    outforLD<-pop1	
	if(founder_pop1[1,2]<0 | floor(founder_pop1[1,2])!=founder_pop1[1,2] | founder_pop1[1,2]>(length(outforLD$output)-1)) {
	stop('\n','--- Generation index for male founders in argument "founder_pop1" should be an integer in range 0<=x<=number of generations from pop1')
	}
	
	if(founder_pop1[2,2]<0 | floor(founder_pop1[2,2])!=founder_pop1[2,2] | founder_pop1[2,2]>(length(outforLD$output)-1)) {
	stop('\n','--- Generation index for female founders in argument "founder_pop1" should be an integer in range 0<=x<=number of generations from pop1')
	}
	
	tst<-outforLD$output[[(founder_pop1[1,2]+1)]]$data
	t_mad<-subset(tst,tst[,5]=='M')
	t_mad<-length(t_mad[,5])
	
	if(founder_pop1[1,1]>t_mad) {
	stop('\n','--- Number of males in argument "founder_pop1" is more than number of males in generation: ',founder_pop1[1,2],' of previous population.',"\n",'Decrease number of selected males')
	}
	
	tst<-outforLD$output[[(founder_pop1[2,2]+1)]]$data
	t_mad<-subset(tst,tst[,5]=='F')
	t_mad<-length(t_mad[,5])
	
	if(founder_pop1[2,1]>t_mad) {
	stop('\n','--- Number of females in argument "founder_pop1" is more than number of females in generation: ',founder_pop1[2,2],' of previous population.',"\n",'Decrease number of selected females')
	}
	
	
	
	test_Males <- c('rnd','tbv','phen','gebv')
	if(any(founder_pop1[,3]%in%test_Males==FALSE)){
	stop('--- possible options for selection of founders in argument "founder_pop1" are "phen","tbv","rnd","gebv"')
	}
	
	
	test_Males <- c('tbv','phen','gebv')
	hl<-c('h','l')
	if(any(founder_pop1[,3]%in%test_Males==TRUE)){
	if(length(founder_pop1[1,])!=4){
	stop('\n','--- If males/females in argument "founder_pop1" are selected based on "phen","tbv" or "gebv", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}
	
	test_Males <- c('tbv','phen','gebv')
	hl<-c('h','l')
	if(any(founder_pop1[,3]%in%test_Males==TRUE)){
	if(any(founder_pop1[,4]%in%hl==FALSE)){
	stop('\n','--- If males/females in argument "founder_pop1" are selected based on "phen","tbv" or "gebv", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}
				
	test_Selection <- c('gebv')
	tst_gebv<-colnames(outforLD$summary_data)
	if(any(founder_pop1[,3]%in%test_Selection==TRUE)){
	if(is.na(match('GEBV',tst_gebv))){
	stop('\n','--- "gebv" as selection criteria for selected sires/dams in argument "founder_pop1" is not allowed as data for this selection criteria is not provided by argument "pop1"')
	}
	}

	
# founder_pop2
	if(missing(founder_pop2)) {
	stop('\n','--- argument "founder_pop2" is missing')
	}
	
	founder_pop2[,1]<-as.numeric(founder_pop2[,1])
	founder_pop2[,2]<-as.numeric(founder_pop2[,2])
	if(	length(which(founder_pop2[,1]<1))>0 | any(floor(founder_pop2[,1])!=founder_pop2[,1]))  {
	stop('\n','--- Number of selected males/females in argument "founder_pop2" should be an integer grater than 1')
	}

    outforLD<-pop2	
	if(founder_pop2[1,2]<0 | floor(founder_pop2[1,2])!=founder_pop2[1,2] | founder_pop2[1,2]>(length(outforLD$output)-1)) {
	stop('\n','--- Generation index for male founders in argument "founder_pop2" should be an integer in range 0<=x<=number of generations from pop2')
	}
	
	if(founder_pop2[2,2]<0 | floor(founder_pop2[2,2])!=founder_pop2[2,2] | founder_pop2[2,2]>(length(outforLD$output)-1)) {
	stop('\n','--- Generation index for female founders in argument "founder_pop2" should be an integer in range 0<=x<=number of generations from pop2')
	}
	
	
	tst<-outforLD$output[[(founder_pop2[1,2]+1)]]$data
	t_mad<-subset(tst,tst[,5]=='M')
	t_mad<-length(t_mad[,5])
	
	if(founder_pop2[1,1]>t_mad) {
	stop('\n','--- Number of males in argument "founder_pop2" is more than number of males in generation: ',founder_pop2[1,2],' of previous population.',"\n",'Decrease number of selected males')
	}
	
	tst<-outforLD$output[[(founder_pop2[2,2]+1)]]$data
	t_mad<-subset(tst,tst[,5]=='F')
	t_mad<-length(t_mad[,5])
	
	if(founder_pop2[2,1]>t_mad) {
	stop('\n','--- Number of females in argument "founder_pop2" is more than number of females in generation: ',founder_pop2[2,2],' of previous population.',"\n",'Decrease number of selected females')
	}
	
		
	test_Males <- c('rnd','tbv','phen','gebv')
	if(any(founder_pop2[,3]%in%test_Males==FALSE)){
	stop('\n','--- Poosible options for selection of founders in argument "founder_pop2" are "phen","tbv","rnd","gebv"')
	}
	
	test_Males <- c('tbv','phen','gebv')
	hl<-c('h','l')
	if(any(founder_pop2[,3]%in%test_Males==TRUE)){
	if(length(founder_pop2[1,])!=4){
	stop('\n','--- If males/females in argument "founder_pop2" are selected based on "phen","tbv" or "gebv", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}
	
	test_Males <- c('tbv','phen','gebv')
	hl<-c('h','l')
	if(any(founder_pop2[,3]%in%test_Males==TRUE)){
	if(any(founder_pop2[,4]%in%hl==FALSE)){
	stop('\n','--- If males/females in argument "founder_pop2" are selected based on "phen","tbv" or "gebv", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}
	
	test_Selection <- c('gebv')
	tst_gebv<-colnames(outforLD$summary_data)
	if(any(founder_pop2[,3]%in%test_Selection==TRUE)){
	if(is.na(match('GEBV',tst_gebv))){
	stop('\n','--- "gebv" as selection criteria for selected sires/dams in argument "founder_pop2" is not allowed  as data for this selection criteria is not provided by argument "pop2"')
	}
	}
	
	# founder_cross
	if(missing(founder_cross)) {
	stop('\n','--- argument "founder_cross" is missing')
	}
	
	test <- c('pop1','pop2')
	if(any(founder_cross[,1]%in%test==FALSE)){
	stop('\n','--- Poosible options for cross founders in argument "founder_cross" are "pop1","pop2"')
	}
	
	founder_cross[,2]<-as.numeric(founder_cross[,2])
	if(length(which(founder_cross[,2]<1))>0 | any(floor(founder_cross[,2])!=founder_cross[,2]))  {
	stop('\n','--- Number of selected males/females in argument "founder_cross" should be an integer grater than 1')
	}
	
	if(founder_cross[1,1]=='pop1'){ #Males
		if(founder_cross[1,2]>founder_pop1[1,1]){
		stop('\n','--- Number of selected males in argument "founder_cross" should not be grater than number of sires in "founder_pop1"')
		} else if (founder_cross[1,1]=='pop2'){
		if(founder_cross[1,2]>founder_pop2[1,1]){
		stop('\n','--- Number of selected males in argument "founder_cross" should not be grater than number of sires in "founder_pop2"')
		}
		}
	}
	
	if(founder_cross[2,1]=='pop1'){ #Females
		if(founder_cross[2,2]>founder_pop1[2,1]){
		stop('\n','--- Number of selected females in argument "founder_cross" should not be grater than number of dams in "founder_pop1"')
		} else if (founder_cross[2,1]=='pop2'){
		if(founder_cross[2,2]>founder_pop2[2,1]){
		stop('\n','--- Number of selected females in argument "founder_cross" should not be grater than number of dams in "founder_pop2"')
		}
		}
	}

	test_Males <- c("rnd","tbv","phen","gebv")
	if(any(founder_cross[,3]%in%test_Males==FALSE)){
	stop('\n','---Possible options for cross founders in argument "founder_cross" are "rnd","tbv","phen","gebv"')
	}
	
		if(founder_cross[1,1]=='pop1' & founder_cross[1,3]=='gebv'){ 
		tst<-pop1$user_input[[5]][1,2]	
		if(tst!='gebv'){
			stop('\n','--- Sires in argument "founder_cross" can not be slected based on "gebv" as data for this selection criteria is not provided by argument "pop1"')
		}
		}
		
		if(founder_cross[1,1]=='pop2' & founder_cross[1,3]=='gebv'){ 
		tst<-pop2$user_input[[5]][1,2]	
		if(tst!='gebv'){
			stop('\n','--- Sires in argument "founder_cross" can not be slected based on "gebv" as data for this selection criteria is not provided by argument "pop2"')
		}
		}
		
		if(founder_cross[2,1]=='pop1' & founder_cross[2,3]=='gebv'){ 
		tst<-pop1$user_input[[5]][2,2]	
		if(tst!='gebv'){
			stop('\n','--- Dams in argument "founder_cross" can not be slected based on "gebv" as data for this selection criteria is not provided by argument "pop1"')
		}
		}
		
		if(founder_cross[2,1]=='pop2' & founder_cross[2,3]=='gebv'){ 
		tst<-pop2$user_input[[5]][2,2]	
		if(tst!='gebv'){
			stop('\n','--- Dams in argument "founder_cross" can not be slected based on "gebv" as data for this selection criteria is not provided by argument "pop2"')
		}
		}


	test_Males <- c('tbv','phen','gebv')
	hl<-c('h','l')
	if(any(founder_cross[,3]%in%test_Males==TRUE)){
	if(length(founder_cross[1,])!=4){
	stop('\n','--- If males/females in argument "founder_cross" are selected based on "phen","tbv" or "gebv", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}
			
	test_Males <- c('tbv','phen','gebv')
	hl<-c('h','l')
	if(any(founder_cross[,3]%in%test_Males==TRUE)){
	if(any(founder_cross[,4]%in%hl==FALSE)){
	stop('\n','--- If males/females in argument "founder_cross" are selected based on "phen","tbv" or "gebv", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}
	
	
	# Selection_pop1
	if(missing(Selection_pop1)) {
	stop('\n','--- argument "Selection_pop1" is missing')
	}
	
		# males
		if(Selection_pop1[1,1]<1 | floor(Selection_pop1[1,1])!=Selection_pop1[1,1] ) {
		stop('\n','--- Number of selected sires in argument "Selection_pop1" should be an integer grater than 1')
		}
		
		# females
		if(Selection_pop1[2,1]<1 | floor(Selection_pop1[2,1])!=Selection_pop1[2,1] ) {
		stop('\n','--- Number of selected dams in argument "Selection_pop1" should be an integer grater than 1')
		}
		
		test_Selection <- c('rnd','phen','tbv','gebv','tbvc','gebvc')
		if(any(Selection_pop1[,2]%in%test_Selection==FALSE)){
		stop('\n','--- Selection criteria for selected sires/dams in argument "Selection_pop1" can be  "rnd","phen","tbv","gebv","tbvc" or "gebvc"')
		}
		
		test_Selection <- c('phen','tbv','gebv','tbvc','gebvc')
		hl<-c('h','l')
		if(any(Selection_pop1[,2]%in%test_Selection==TRUE)){
		if(length(Selection_pop1[1,])!=3){
		stop('\n','--- If selection criteria for selected sires/dams in argument "Selection_pop1" is  "phen",  "tbv","gebv" "tbvc" or "gebvc", it should be defined to select "h" (high) or "l" (low) values.')
		}
		}
				
		test_Selection <- c('phen','tbv','gebv','tbvc','gebvc')
		hl<-c('h','l')
		if(any(Selection_pop1[,2]%in%test_Selection==TRUE)){
		if(any(Selection_pop1[,3]%in%hl==FALSE)){
		stop('\n','--- If selection criteria for selected sires/dams in argument "Selection_pop1" is  "phen",  "tbv","gebv" "tbvc" or "gebvc", it should be defined to select "h" (high) or "l" (low) values.')
		}
		}
		
		
		if(Selection_pop1[1,1]>0.5*(Selection_pop1[2,1]*litter_size)){
		stop("\n",'---  Number of selected males in argument "Selection_pop1" is greater than number of available male progeny',"\n",'Solution 1: Increase litter size.'
		,"\n",'Solution 2: Decrease number of selected sires. in argument "Selection_pop1"',"\n",'Solution 3: Increase number of selected dams in argument "Selection_pop1".')
		}
		
		
		
# Selection_pop2
	if(missing(Selection_pop2)) {
	stop('\n','--- argument "Selection_pop2" is missing')
	}
	
	# males
	if(Selection_pop2[1,1]<1 | floor(Selection_pop2[1,1])!=Selection_pop2[1,1] ) {
	stop('\n','--- Number of selected sires in argument "Selection_pop2" should be an integer grater than 1')
	}
		
	# females
	if(Selection_pop2[2,1]<1 | floor(Selection_pop2[2,1])!=Selection_pop2[2,1] ) {
	stop('\n','--- Number of selected dams in argument "Selection_pop2" should be an integer grater than 1')
	}
		
	test_Selection <- c('rnd','phen','tbv','gebv','tbvc','gebvc')
	if(any(Selection_pop2[,2]%in%test_Selection==FALSE)){
	stop('\n','--- Selection criteria for selected sires/dams in argument "Selection_pop2" can be  "rnd","phen","tbv","gebv","tbvc" or "gebvc"')
	}
	
	test_Selection <- c('phen','tbv','gebv','tbvc','gebvc')
	hl<-c('h','l')
	if(any(Selection_pop2[,2]%in%test_Selection==TRUE)){
	if(length(Selection_pop2[1,])!=3){
	stop('\n','--- If selection criteria for selected sires/dams in argument "Selection_pop2" is  "phen",  "tbv","gebv" "tbvc" or "gebvc", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}
		
	test_Selection <- c('phen','tbv','gebv','tbvc','gebvc')
	hl<-c('h','l')
	if(any(Selection_pop2[,2]%in%test_Selection==TRUE)){
	if(any(Selection_pop2[,3]%in%hl==FALSE)){
	stop('\n','--- if selection criteria for selected sires/dams in argument "Selection_pop2" is  "phen",  "tbv","gebv" "tbvc" or "gebvc", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}	
	
	
	if(Selection_pop2[1,1]>0.5*(Selection_pop2[2,1]*litter_size)){
	stop("\n",'--- Number of selected males in argument "Selection_pop2" is greater than number of available male progeny',"\n",'Solution 1: Increase litter size.'
	,"\n",'Solution 2: Decrease number of selected sires. in argument "Selection_pop2"',"\n",'Solution 3: Increase number of selected dams in argument "Selection_pop2".')
	}
		
		
		
	
		# ng
	if(missing(ng)) {
	stop('\n','--- ng is missing')
	}

	if(ng<1 | ng>200 | floor(ng)!=ng) {
	stop('\n','--- Number of generations should be an integer in range 1-200')
	}
	
	# litter_size
	if(missing(litter_size)) {
	stop('\n','--- litter_size is missing')
	}

	if(litter_size<1 | litter_size>200 | floor(litter_size)!=litter_size) {
	stop('\n','--- argument "litter_size" should be an integer in range 1<=x<=200')
	}
		
	# Cross_design
	if(missing(Cross_design)) {
	stop('\n','--- argument "Cross_design" is missing')
	}
	
	test <- c('pop1','pop2')
	if(any(Cross_design[,1]%in%test==FALSE)){
	stop('\n','--- Possible options for cross founders in argument "Cross_design" are "pop1","pop2"')
	}
	
	# Cross_design[1,2]<-as.numeric(Cross_design[1,2])
	Cross_design[,2]<-as.numeric(Cross_design[,2])
	# cat(class(Cross_design[1,2]),fill=TRUE)
	# cat(class(Cross_design[2,2]),fill=TRUE)
	
	if(length(which(Cross_design[,2]<1))>0 | any(floor(Cross_design[,2])!=Cross_design[,2]))  {
	stop('\n','--- Number of selected males/females in argument "Cross_design" should be an integer grater than 1')
	}
	
		if(Cross_design[1,1]=='pop1'){ #Males
				if(Cross_design[1,2]>floor(0.5*Selection_pop1[2,1]*litter_size)){
				stop('\n','--- Number of selected males: ',Cross_design[1,2],' in argument "Cross_design" is grater than number of available sires: ',floor(0.5*Selection_pop1[2,1]*litter_size),'  provided by "Selection_pop1"')
				}
		}
		
		if (Cross_design[1,1]=='pop2'){
			if(Cross_design[1,2]>floor(0.5*Selection_pop2[2,1]*litter_size)){
			stop('\n','--- Number of selected males: ',Cross_design[1,2],' in argument "Cross_design" is grater than number of available sires: ',floor(0.5*Selection_pop2[2,1]*litter_size),'  provided by "Selection_pop2"')
			}
		}
	
	
	if(Cross_design[2,1]=='pop1'){ #Females
		if(Cross_design[2,2]>floor(0.5*Selection_pop1[2,1]*litter_size)){
		stop('\n','--- Number of selected females: ',Cross_design[2,2],' in argument "Cross_design" is grater than number of available dams: ',floor(0.5*Selection_pop1[2,1]*litter_size),'  provided by "Selection_pop1"')
			}
	}
			
   if (Cross_design[2,1]=='pop2'){
		if(Cross_design[2,2]>floor(0.5*Selection_pop2[2,1]*litter_size)){
		stop('\n','--- Number of selected females: ',Cross_design[2,2],' in argument "Cross_design" is grater than number of available dams: ',floor(0.5*Selection_pop2[2,1]*litter_size),'  provided by "Selection_pop2"')
		}
		}

	
	test_Selection <- c('rnd','phen','tbv','gebv','tbvc','gebvc')
	if(any(Cross_design[,3]%in%test_Selection==FALSE)){
	stop('\n','--- Selection criteria for selected sires/dams in argument "Cross_design" can be  "rnd","phen","tbv","gebv","tbvc" or "gebvc"')
	}
	
	
	test_Selection <- c('phen','tbv','gebv','tbvc','gebvc')
	hl<-c('h','l')
	if(any(Cross_design[,3]%in%test_Selection==TRUE)){
	if(length(Cross_design[1,])!=4){
	stop('\n','--- if selection criteria for selected sires/dams in argument "Cross_design" is  "phen",  "tbv","gebv" "tbvc" or "gebvc", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}	
	
	test_Selection <- c('phen','tbv','gebv','tbvc','gebvc')
	hl<-c('h','l')
	if(any(Cross_design[,3]%in%test_Selection==TRUE)){
	if(any(Cross_design[,4]%in%hl==FALSE)){
	stop('\n','--- if selection criteria for selected sires/dams in argument "Cross_design" is  "phen",  "tbv","gebv" "tbvc" or "gebvc", it should be defined to select "h" (high) or "l" (low) values.')
	}
	}	
	
	# train_type
	test_Selection <- c('gebv','gebvc')
	if(any(Cross_design[,3]%in%test_Selection==TRUE) | any(	Selection_pop1[,2]%in%test_Selection==TRUE)| any(Selection_pop2[,2]%in%test_Selection==TRUE) ){
		if(missing(train_type)){
		stop('\n','--- if selection criteria for any population is defined as gebv/gebvc, argument "train_type" should be defined.')
	}
		if(!missing(train_type)){
		tst<-c('crossbred','purebred') 
		if(train_type%in%tst==FALSE)
		stop('\n','--- argument "train_type" can be either "crossbred" or "purebred"')
	}
	}
	
 # train_pop1
  	test_Selection <- c('gebv','gebvc')
	if(any(Cross_design[,3]%in%test_Selection==TRUE) | any(	Selection_pop1[,2]%in%test_Selection==TRUE)){
	
	 if(train_type=='purebred'){
			if(missing(train_pop1)){
			stop('\n','--- train_type is "purebred" but argument "train_pop1" is missing')}
			if(!missing(train_pop1)){
		control_names<-c('size','method','sel','nIter','burnIn','thin','save','show')
		test_w<-intersect(names(train_pop1),control_names)
		
		if(any(test_w=='size')==FALSE){
		stop('\n','--- size in argument "train_pop1" is  missing.')}
		if(any(test_w=='method')==FALSE){
			cat('\n','method in argument "train_pop1" is set to default method:"BRR"')}
			
	# size
	if(any(test_w=='size')){
		tra_size<-train_pop1$size
		if(floor(tra_size)!=tra_size | tra_size<1){
		stop('\n','--- size in argument "train_pop1" should be positive integer.')
		}
		if(tra_size>(litter_size*Selection_pop1[2,1])){
		avai<-litter_size*Selection_pop1[2,1]
		cat('\n','--- Training size:',tra_size,' in argument "train_pop1" is grater than number of available individuals:',avai,fill=TRUE)
		stop("\n",'Solution 1: Increase litter size.',"\n",'Solution 2: Decrease size of training.',"\n",'Solution 3: Increase number of dmas selected in argument "Selection_pop1" if possible.')		
		}
		
		if(tra_size>(litter_size*founder_pop1[2,1])){
		avai<-litter_size*as.numeric(founder_pop1[2,1])
		cat('\n','--- Training size:',tra_size,' in argument "train_pop1" is grater than number of available individuals:',avai,fill=TRUE)
		stop("\n",'Solution 1: Increase litter size.',"\n",'Solution 2: Decrease size of training.',"\n",'Solution 3: Increase number of dmas selected in argument "founder_pop1" if possible.')		
		}
	}
	
	# method
	test_method <- c('BRR', 'BayesA', 'BL', 'BayesB', 'BayesC')
	if(any(test_w=='method')){
	tra_method<-train_pop1$method
		if(tra_method%in%test_method==FALSE){
		stop('\n','--- possible options for method in argument "train_pop1" are:"BRR", "BayesA", "BL", "BayesB", "BayesC"')
		}
	} else {
	tra_method<-'BRR'
	}
	
	# sel
	test_sel <- c('rnd', 'min_rel_mrk', 'max_rel_mrk', 'min_rel_qtl', 'max_rel_qtl')
	if(any(test_w=='sel')){
	tra_sel<-train_pop1$sel
		if(tra_sel%in%test_sel==FALSE){
		stop('\n','--- possible options for sel in argument "train_pop1" are:"rnd", "min_rel_mrk", "max_rel_mrk", "min_rel_qtl", "max_rel_qtl"')
		}
	} else {
	tra_sel<-'rnd'
	}
	
	# nIter
	if(any(test_w=='nIter')){
	tra_nIter<-train_pop1$nIter
		if(floor(tra_nIter)!=tra_nIter | tra_nIter<1){
		stop('\n','--- tra_nIter in argument "train_pop1" should be positive integer.')
		}
	} else {
	tra_nIter<-1500
	}
	
	# burnIn
	if(any(test_w=='burnIn')){
	tra_burnIn<-train_pop1$burnIn
		if(floor(tra_burnIn)!=tra_burnIn | tra_burnIn<1){
		stop('\n','--- tra_burnIn in argument "train_pop1" should be positive integer.')
		}
	} else {
	tra_burnIn<-500
	}
	
	# thin
	if(any(test_w=='thin')){
	tra_thin<-train_pop1$thin
		if(floor(tra_thin)!=tra_thin | tra_thin<1){
		stop('\n','--- tra_thin in argument "train_pop1" should be positive integer.')
		}
	} else {
	tra_thin<-5
	}
		
	# save
	if(any(test_w=='save')){
	tra_save<-train_pop1$save
		if(is.character(tra_save)==FALSE){
		stop('\n','--- Define a name for tra_save in argument "train_pop1" as type "character"')
		}
	} else {
	tra_save<-'Out_BGLR'
	}
	
	# show
	if(any(test_w=='show')){
	tra_show<-train_pop1$show
		if(is.logical(tra_show)==FALSE){
		stop('\n','--- tra_show in argument "train_pop1" should be type "logical"')
		}
	} else {
	tra_show<-TRUE
	}	
		
		} #end loop !missing(train_pop1)
	  } #end loop train_type=='purebred'
	} #end loop any(Cross_design[,3]%in%test_Selection==TRUE) | any(	Selection_pop1[,2]%in%test_Selection==TRUE))
 
 # train_pop2
  	test_Selection <- c('gebv','gebvc')
	if(any(Cross_design[,3]%in%test_Selection==TRUE) | any(	Selection_pop2[,2]%in%test_Selection==TRUE)){
	
	 if(train_type=='purebred'){
			if(missing(train_pop2)){
			stop('\n','--- train_type is "purebred" but argument "train_pop2" is missing')}
			if(!missing(train_pop2)){
		control_names<-c('size','method','sel','nIter','burnIn','thin','save','show')
		test_w<-intersect(names(train_pop2),control_names)
		
		if(any(test_w=='size')==FALSE){
		stop('\n','--- size in argument "train_pop2" is  missing.')}
		if(any(test_w=='method')==FALSE){
			cat('\n','method in argument "train_pop2" is set to default method:"BRR"')}
			
	# size
	if(any(test_w=='size')){
		tra_size<-train_pop2$size
		if(floor(tra_size)!=tra_size | tra_size<1){
		stop('\n','--- size in argument "train_pop2" should be positive integer.')
		}
		if(tra_size>(litter_size*Selection_pop2[2,1])){
		avai<-litter_size*Selection_pop2[2,1]
		cat('\n','--- Training size:',tra_size,' in argument "train_pop2" is grater than number of available individuals:',avai,fill=TRUE)
		stop("\n",'Solution 1: Increase litter size.',"\n",'Solution 2: Decrease size of training.',"\n",'Solution 3: Increase number of dmas selected in argument "Selection_pop2" if possible.')		
		}
		
		if(tra_size>(litter_size*founder_pop2[2,1])){
		avai<-litter_size*as.numeric(founder_pop2[2,1])
		cat('\n','--- Training size:',tra_size,' in argument "train_pop2" is grater than number of available individuals:',avai,fill=TRUE)
		stop("\n",'Solution 1: Increase litter size.',"\n",'Solution 2: Decrease size of training.',"\n",'Solution 3: Increase number of dmas selected in argument "founder_pop2" if possible.')		
		}
	}
	
	# method
	test_method <- c('BRR', 'BayesA', 'BL', 'BayesB', 'BayesC')
	if(any(test_w=='method')){
	tra_method<-train_pop2$method
		if(tra_method%in%test_method==FALSE){
		stop('\n','--- possible options for method in argument "train_pop2" are:"BRR", "BayesA", "BL", "BayesB", "BayesC"')
		}
	} else {
	tra_method<-'BRR'
	}
	
	# sel
	test_sel <- c('rnd', 'min_rel_mrk', 'max_rel_mrk', 'min_rel_qtl', 'max_rel_qtl')
	if(any(test_w=='sel')){
	tra_sel<-train_pop2$sel
		if(tra_sel%in%test_sel==FALSE){
		stop('\n','--- possible options for sel in argument "train_pop2" are:"rnd", "min_rel_mrk", "max_rel_mrk", "min_rel_qtl", "max_rel_qtl"')
		}
	} else {
	tra_sel<-'rnd'
	}
	
	# nIter
	if(any(test_w=='nIter')){
	tra_nIter<-train_pop2$nIter
		if(floor(tra_nIter)!=tra_nIter | tra_nIter<1){
		stop('\n','--- tra_nIter in argument "train_pop2" should be positive integer.')
		}
	} else {
	tra_nIter<-1500
	}
	
	# burnIn
	if(any(test_w=='burnIn')){
	tra_burnIn<-train_pop2$burnIn
		if(floor(tra_burnIn)!=tra_burnIn | tra_burnIn<1){
		stop('\n','--- tra_burnIn in argument "train_pop2" should be positive integer.')
		}
	} else {
	tra_burnIn<-500
	}
	
	# thin
	if(any(test_w=='thin')){
	tra_thin<-train_pop2$thin
		if(floor(tra_thin)!=tra_thin | tra_thin<1){
		stop('\n','--- tra_thin in argument "train_pop2" should be positive integer.')
		}
	} else {
	tra_thin<-5
	}
		
	# save
	if(any(test_w=='save')){
	tra_save<-train_pop2$save
		if(is.character(tra_save)==FALSE){
		stop('\n','--- Define a name for tra_save in argument "train_pop2" as type "character"')
		}
	} else {
	tra_save<-'Out_BGLR'
	}
	
	# show
	if(any(test_w=='show')){
	tra_show<-train_pop2$show
		if(is.logical(tra_show)==FALSE){
		stop('\n','--- tra_show in argument "train_pop2" should be type "logical"')
		}
	} else {
	tra_show<-TRUE
	}	
		
		} #end loop !missing(train_pop2)
	  } #end loop train_type=='purebred'
	} #end loop any(Cross_design[,3]%in%test_Selection==TRUE) | any(Selection_pop2[,2]%in%test_Selection==TRUE))
 

# train_cross
  	test_Selection <- c('gebv','gebvc')
	if(any(Cross_design[,3]%in%test_Selection==TRUE) | any(	Selection_pop1[,2]%in%test_Selection==TRUE) | any(	Selection_pop2[,2]%in%test_Selection==TRUE)){
	
	 if(train_type=='crossbred'){
			if(missing(train_cross)){
			stop('\n','--- train_type is "crossbred" but argument "train_cross" is missing')}
			if(!missing(train_cross)){
		control_names<-c('size','method','sel','nIter','burnIn','thin','save','show')
		test_w<-intersect(names(train_cross),control_names)
		
		if(any(test_w=='size')==FALSE){
		stop('\n','--- size in argument "train_cross" is  missing.')}
		if(any(test_w=='method')==FALSE){
			cat('\n','method in argument "train_cross" is set to default method:"BRR"')}
			
	# size
	if(any(test_w=='size')){
		tra_size<-train_cross$size
		if(floor(tra_size)!=tra_size | tra_size<1){
		stop('\n','--- size in argument "train_cross" should be positive integer.')
		}
			if(tra_size>(litter_size*Cross_design[2,2])){
			avai<-litter_size*Cross_design[2,2]
			cat('\n','--- Training size:',tra_size,' in argument "train_cross" is grater than number of available individuals:',avai,fill=TRUE)
			stop("\n",'Solution 1: Increase litter size.',"\n",'Solution 2: Decrease size of training.',"\n",'Solution 3: Increase number of dmas selected in argument "Cross_design" if possible.')		
			}
		
			if(tra_size>(litter_size*founder_cross[2,2])){
			avai<-litter_size*as.numeric(founder_cross[2,2])
			cat('\n','--- Training size:',tra_size,' in argument "train_cross" is grater than number of available individuals:',avai,fill=TRUE)
			stop("\n",'Solution 1: Increase litter size.',"\n",'Solution 2: Decrease size of training.',"\n",'Solution 3: Increase number of dmas selected in argument "founder_cross" if possible.')		
			}
	}
	
	# method
	test_method <- c('BRR', 'BayesA', 'BL', 'BayesB', 'BayesC')
	if(any(test_w=='method')){
	tra_method<-train_cross$method
		if(tra_method%in%test_method==FALSE){
		stop('\n','--- possible options for method in argument "train_cross" are:"BRR", "BayesA", "BL", "BayesB", "BayesC"')
		}
	} else {
	tra_method<-'BRR'
	}
	
	# sel
	test_sel <- c('rnd', 'min_rel_mrk', 'max_rel_mrk', 'min_rel_qtl', 'max_rel_qtl')
	if(any(test_w=='sel')){
	tra_sel<-train_cross$sel
		if(tra_sel%in%test_sel==FALSE){
		stop('\n','--- possible options for sel in argument "train_cross" are:"rnd", "min_rel_mrk", "max_rel_mrk", "min_rel_qtl", "max_rel_qtl"')
		}
	} else {
	tra_sel<-'rnd'
	}
	
	# nIter
	if(any(test_w=='nIter')){
	tra_nIter<-train_cross$nIter
		if(floor(tra_nIter)!=tra_nIter | tra_nIter<1){
		stop('\n','--- tra_nIter in argument "train_cross" should be positive integer.')
		}
	} else {
	tra_nIter<-1500
	}
	
	# burnIn
	if(any(test_w=='burnIn')){
	tra_burnIn<-train_cross$burnIn
		if(floor(tra_burnIn)!=tra_burnIn | tra_burnIn<1){
		stop('\n','--- tra_burnIn in argument "train_cross" should be positive integer.')
		}
	} else {
	tra_burnIn<-500
	}
	
	# thin
	if(any(test_w=='thin')){
	tra_thin<-train_cross$thin
		if(floor(tra_thin)!=tra_thin | tra_thin<1){
		stop('\n','--- tra_thin in argument "train_cross" should be positive integer.')
		}
	} else {
	tra_thin<-5
	}
		
	# save
	if(any(test_w=='save')){
	tra_save<-train_cross$save
		if(is.character(tra_save)==FALSE){
		stop('\n','--- Define a name for tra_save in argument "train_cross" as type "character"')
		}
	} else {
	tra_save<-'Out_BGLR'
	}
	
	# show
	if(any(test_w=='show')){
	tra_show<-train_cross$show
		if(is.logical(tra_show)==FALSE){
		stop('\n','--- tra_show in argument "train_cross" should be type "logical"')
		}
	} else {
	tra_show<-TRUE
	}	
		
		} #end loop !missing(train_cross)
	  } #end loop train_type=='crossbred'
	} #any(Cross_design[,3]%in%test_Selection==TRUE) | any(	Selection_pop1[,2]%in%test_Selection==TRUE) | any(	Selection_pop2[,2]%in%test_Selection==TRUE)
	
	
	# saveAt and output_pop1
	# if(!missing(saveAt) & missing(output_pop1)){
		# stop('\n','--- argument "output_pop1" is missing')
	# }
	
	if(missing(saveAt) & !missing(output_pop1)){
		stop('\n','--- argument "saveAt" is missing')
	}
	
	if(!missing(saveAt) & !missing(output_pop1)){
			if(is.character(saveAt)==FALSE){
			stop('\n','--- Define a name for the saveAt argument as type "character"')
			}
		
	control_names<-c("data","qtl","marker","seq","freq_qtl","freq_mrk")
	test_w<-intersect(names(output_pop1),control_names)
	if(length(test_w)==0){
		stop('\n','--- possible output files in argument "output_pop1" are:"data","qtl","marker","seq","freq_qtl","freq_mrk"')
		}

	if(any(floor(output_pop1)==output_pop1)==FALSE){
		stop('\n','--- generation indexes in argument "output_pop1" for writting output files should be integer')
		}
		
	if(any(floor(output_pop1)>ng)==TRUE){
		stop('\n','--- Input error(s) in generation indexes in argument "output_pop1" for writting output files')
		}
	} # end loop !missing(saveAt) & !missing(output_pop1)

	# saveAt and output_pop2
	# if(!missing(saveAt) & missing(output_pop2)){
		# stop('\n','--- argument "output_pop2" is missing')
	# }
	
	if(missing(saveAt) & !missing(output_pop2)){
		stop('\n','--- argument "saveAt" is missing')
	}
	
	if(!missing(saveAt) & !missing(output_pop2)){
			if(is.character(saveAt)==FALSE){
			stop('\n','--- Define a name for the saveAt argument as type "character"')
			}
		
	control_names<-c("data","qtl","marker","seq","freq_qtl","freq_mrk")
	test_w<-intersect(names(output_pop2),control_names)
	if(length(test_w)==0){
		stop('\n','--- possible output files in argument "output_pop2" are:"data","qtl","marker","seq","freq_qtl","freq_mrk"')
		}

	if(any(floor(output_pop2)==output_pop2)==FALSE){
		stop('\n','--- generation indexes in argument "output_pop2" for writting output files should be integer')
		}
		
	if(any(floor(output_pop2)>ng)==TRUE){
		stop('\n','--- Input error(s) in generation indexes in argument "output_pop2" for writting output files')
		}
	} # end loop !missing(saveAt) & !missing(output_pop2)
	

	# saveAt and output_cross
	# if(!missing(saveAt) & missing(output_cross)){
		# stop('\n','--- argument "output_cross" is missing')
	# }
	
	if(missing(saveAt) & !missing(output_cross)){
		stop('\n','--- argument "saveAt" is missing')
	}
	
	if(!missing(saveAt) & !missing(output_cross)){
			if(is.character(saveAt)==FALSE){
			stop('\n','--- Define a name for the saveAt argument as type "character"')
			}
		
	control_names<-c("data","qtl","marker","seq","freq_qtl","freq_mrk")
	test_w<-intersect(names(output_cross),control_names)
	if(length(test_w)==0){
		stop('\n','--- possible output files in argument "output_cross" are:"data","qtl","marker","seq","freq_qtl","freq_mrk"')
		}

	if(any(floor(output_cross)==output_cross)==FALSE){
		stop('\n','--- generation indexes in argument "output_cross" for writting output files should be integer')
		}
		
	if(any(floor(output_cross)>ng)==TRUE){
		stop('\n','--- Input error(s) in generation indexes in argument "output_cross" for writting output files')
		}
	} # end loop !missing(saveAt) & !missing(output_cross)
	
	# Display
	if(missing(Display)) {
	Display<-TRUE
	}
	if(!missing(Display)) {
			if(is.logical(Display)==FALSE){
			stop('\n','--- argument "Display" should be type "logical"')
			}
	}
	
	
		
	# if(all(Selection[,2]%in%test_Selection==FALSE) & !missing(Training)){
	# warning('\n','argument "Training" is ignored as selection criteria is not "gebv"')
	# }
	
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

	funi2<-function(vec){
		length(which(vec>=0))
		}
		

# Necessary internal functions---Finish

# Total_Total<-list(list(),list(),list(),list())
Total_Total_A<-replicate((ng+1), list_fun(), simplify=FALSE)
Total_Total_B<-replicate((ng+1), list_fun(), simplify=FALSE)
Total_Total_C<-replicate((ng+1), list_fun(), simplify=FALSE)

cat('Intializing base population ...',fill=TRUE)
# calculation
flag_sel_A<-intersect(c('gebv','gebvc'),c(Selection_pop1[,2],Cross_design[,3]))
if(length(flag_sel_A)>0){
flag_sel_A<-TRUE
} else
flag_sel_A<-FALSE
# cat(flag_sel_A,"flag_sel_A",fill=TRUE)

flag_sel_B<-intersect(c('gebv','gebvc'),c(Selection_pop2[,2],Cross_design[,3]))
if(length(flag_sel_B)>0){
flag_sel_B<-TRUE
} else
flag_sel_B<-FALSE
flag_sel_B
# cat(flag_sel_B,"flag_sel_B",fill=TRUE)
# effectsi
outforLD<-pop1
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
nmarker_allele<-sum(outforLD$genome[,3])*2
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


# POPULATION 1 DATA PREPARATION START-------------------
cat('Population 1 data preparation ...',fill=TRUE)
	# SIRE DATA POP1
	founder_pop1[,1]<-as.numeric(founder_pop1[,1])
	founder_pop1[,2]<-as.numeric(founder_pop1[,2])
	
	my_internal_data1_A<-pop1$output[[founder_pop1[1,2]+1]]
	my_internal_data<-my_internal_data1_A[[1]]
	my_internal_data
	provided_sire_data<-subset(my_internal_data,my_internal_data[,5]=='M')
	provided_sire_data
	length(provided_sire_data[,1])
	index_sires<-provided_sire_data[,1]
	index_sires_1<-match(index_sires,my_internal_data1_A[[1]][,1])
	index_sires_1
	provided_sire_seq<-my_internal_data1_A[[4]][index_sires_1,]
	# dim(provided_sire_seq)
	provided_sire_qtl<-my_internal_data1_A[[2]][index_sires_1,]
	# dim(provided_sire_qtl)
	provided_sire_mrk<-my_internal_data1_A[[3]][index_sires_1,]
	# dim(provided_sire_mrk)

	#rnd,phen,tbv,gebv
	
	# sires_toC
	if(founder_pop1[1,3]=='rnd'){
	ID_sires_toC<-sample(index_sires,founder_pop1[1,1],replace=FALSE)
	index<-match(ID_sires_toC,index_sires)
	} else if(founder_pop1[1,3]=='phen'){

		if(founder_pop1[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,6]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_pop1[1,1],]
		index<-match(ID_sires_toC[,1],index_sires)

		}
		if(founder_pop1[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,6]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_pop1[1,1],]
		index<-match(ID_sires_toC[,1],index_sires)
		}

	} else if (founder_pop1[1,3]=='tbv'){
		if(founder_pop1[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,8]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_pop1[1,1],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
			if(founder_pop1[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,8]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_pop1[1,1],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
		
	} else if (founder_pop1[1,3]=='gebv'){
		if(founder_pop1[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,9]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_pop1[1,1],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
			if(founder_pop1[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,9]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_pop1[1,1],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
		
	}

	sires_toC_data<-provided_sire_data[index,]
	sires_toC_seq<-provided_sire_seq[index,]
	sires_toC_qtl<-provided_sire_qtl[index,]
	sires_toC_mrk<-provided_sire_mrk[index,]

    # DAM DATA POP1
	my_internal_data1_A<-pop1$output[[founder_pop1[2,2]+1]]
	my_internal_data<-my_internal_data1_A[[1]]
	my_internal_data
	provided_dam_data<-subset(my_internal_data,my_internal_data[,5]=='F')
	provided_dam_data
	length(provided_dam_data[,1])
	index_dams<-provided_dam_data[,1]
	index_dams_1<-match(index_dams,my_internal_data1_A[[1]][,1])
	index_dams_1
	provided_dam_seq<-my_internal_data1_A[[4]][index_dams_1,]
	# dim(provided_dam_seq)
	provided_dam_qtl<-my_internal_data1_A[[2]][index_dams_1,]
	# dim(provided_dam_qtl)
	provided_dam_mrk<-my_internal_data1_A[[3]][index_dams_1,]
	# dim(provided_dam_mrk)

	# dams_toC
	if(founder_pop1[2,3]=='rnd'){
	ID_dams_toC<-sample(index_dams,founder_pop1[2,1],replace=FALSE)
	index<-match(ID_dams_toC,index_dams)
	} else if(founder_pop1[2,3]=='phen'){

		if(founder_pop1[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,6]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_pop1[2,1],]
		index<-match(ID_dams_toC[,1],index_dams)

		}
		if(founder_pop1[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,6]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_pop1[2,1],]
		index<-match(ID_dams_toC[,1],index_dams)
		}

	} else if (founder_pop1[2,3]=='tbv'){
		if(founder_pop1[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,8]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_pop1[2,1],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
			if(founder_pop1[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,8]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_pop1[2,1],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
		
	} else if (founder_pop1[2,3]=='gebv'){
		if(founder_pop1[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,9]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_pop1[2,1],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
			if(founder_pop1[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,9]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_pop1[2,1],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
		
	}
	

	dams_toC_data<-provided_dam_data[index,]
	dams_toC_seq<-provided_dam_seq[index,]
	# dim(dams_toC_data)
	dams_toC_qtl<-provided_dam_qtl[index,]
	dams_toC_mrk<-provided_dam_mrk[index,]
	# dim(dams_toC_mrk)
	

	start_data<-rbind(sires_toC_data,dams_toC_data)
	addeded_colmns<-matrix(0,nrow=length(start_data[,1]),ncol=2)
	total_A<-list()
	total_A$data<-start_data
	total_A$data[,4]<-0
	total_A$data
	total_A$data<-cbind(total_A$data,addeded_colmns) 
	# dim(total_A$data)
	colnames(total_A$data) <- c('id','sire','dam','generation','sex','phen','env','tbv','gebv','tbv_c','gebv_c')
    # total_A$data

	total_A$qtl<-rbind(sires_toC_qtl,dams_toC_qtl)
	total_A$qtl[,2]<-0
	# dim(total_A$qtl)
	# class(total_A$qtl)
	total_A$qtl<-as.matrix(total_A$qtl)
	# class(total_A$qtl)

	total_A$mrk<-rbind(sires_toC_mrk,dams_toC_mrk)
	total_A$mrk[,2]<-0
	# dim(total_A$mrk)
	total_A$mrk<-as.matrix(total_A$mrk)

	# seq
	total_A$sequ<-rbind(sires_toC_seq,dams_toC_seq)
	# dim(total_A$sequ)
	total_A$sequ<-as.matrix(total_A$sequ)
	total_A$sequ[,2]<-0
	
# POPULATION 1 DATA PREPARATION FINISH--------------
# POPULATION 2 DATA PREPARATION START-------------------
cat('Population 2 data preparation ...',fill=TRUE)
	# SIRE DATA POP2
	founder_pop2[,1]<-as.numeric(founder_pop2[,1])
	founder_pop2[,2]<-as.numeric(founder_pop2[,2])
	founder_pop2
	my_internal_data1_B<-pop2$output[[founder_pop2[1,2]+1]]
	my_internal_data<-my_internal_data1_B[[1]]
	my_internal_data
	provided_sire_data<-subset(my_internal_data,my_internal_data[,5]=='M')
	provided_sire_data
	length(provided_sire_data[,1])
	index_sires<-provided_sire_data[,1]
	index_sires_1<-match(index_sires,my_internal_data1_B[[1]][,1])
	index_sires_1
	provided_sire_seq<-my_internal_data1_B[[4]][index_sires_1,]
	# dim(provided_sire_seq)
	provided_sire_qtl<-my_internal_data1_B[[2]][index_sires_1,]
	# dim(provided_sire_qtl)
	provided_sire_mrk<-my_internal_data1_B[[3]][index_sires_1,]
	# dim(provided_sire_mrk)

	# sires_toC
	if(founder_pop2[1,3]=='rnd'){
	ID_sires_toC<-sample(index_sires,founder_pop2[1,1],replace=FALSE)
	index<-match(ID_sires_toC,index_sires)
	} else if(founder_pop2[1,3]=='phen'){

		if(founder_pop2[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,6]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_pop2[1,1],]
		index<-match(ID_sires_toC[,1],index_sires)

		}
		if(founder_pop2[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,6]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_pop2[1,1],]
		index<-match(ID_sires_toC[,1],index_sires)
		}

	} else if (founder_pop2[1,3]=='tbv'){
		if(founder_pop2[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,8]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_pop2[1,1],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
			if(founder_pop2[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,8]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_pop2[1,1],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
		
	} else if (founder_pop2[1,3]=='gebv'){
		if(founder_pop2[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,9]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_pop2[1,1],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
			if(founder_pop2[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,9]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_pop2[1,1],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
		
	}

	sires_toC_data<-provided_sire_data[index,]
	sires_toC_seq<-provided_sire_seq[index,]
	sires_toC_qtl<-provided_sire_qtl[index,]
	sires_toC_mrk<-provided_sire_mrk[index,]

    # DAM DATA POP2
	my_internal_data1_B<-pop2$output[[founder_pop2[2,2]+1]]
	my_internal_data<-my_internal_data1_B[[1]]
	my_internal_data
	provided_dam_data<-subset(my_internal_data,my_internal_data[,5]=='F')
	provided_dam_data
	length(provided_dam_data[,1])
	index_dams<-provided_dam_data[,1]
	index_dams_1<-match(index_dams,my_internal_data1_B[[1]][,1])
	index_dams_1
	provided_dam_seq<-my_internal_data1_B[[4]][index_dams_1,]
	# dim(provided_dam_seq)
	provided_dam_qtl<-my_internal_data1_B[[2]][index_dams_1,]
	# dim(provided_dam_qtl)
	provided_dam_mrk<-my_internal_data1_B[[3]][index_dams_1,]
	# dim(provided_dam_mrk)

	# dams_toC
	if(founder_pop2[2,3]=='rnd'){
	ID_dams_toC<-sample(index_dams,founder_pop2[2,1],replace=FALSE)
	index<-match(ID_dams_toC,index_dams)
	} else if(founder_pop2[2,3]=='phen'){

		if(founder_pop2[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,6]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_pop2[2,1],]
		index<-match(ID_dams_toC[,1],index_dams)

		}
		if(founder_pop2[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,6]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_pop2[2,1],]
		index<-match(ID_dams_toC[,1],index_dams)
		}

	} else if (founder_pop2[2,3]=='tbv'){
		if(founder_pop2[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,8]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_pop2[2,1],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
			if(founder_pop2[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,8]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_pop2[2,1],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
		
	} else if (founder_pop2[2,3]=='gebv'){
		if(founder_pop2[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,9]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_pop2[2,1],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
			if(founder_pop2[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,9]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_pop2[2,1],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
		
	}
	

	dams_toC_data<-provided_dam_data[index,]
	dams_toC_seq<-provided_dam_seq[index,]
	# dim(dams_toC_data)
	dams_toC_qtl<-provided_dam_qtl[index,]
	dams_toC_mrk<-provided_dam_mrk[index,]
	# dim(dams_toC_mrk)
	

	start_data<-rbind(sires_toC_data,dams_toC_data)
	addeded_colmns<-matrix(0,nrow=length(start_data[,1]),ncol=2)
	total_B<-list()
	total_B$data<-start_data
	total_B$data[,4]<-0
	total_B$data
	total_B$data<-cbind(total_B$data,addeded_colmns) 
	# dim(total_B$data)
	colnames(total_B$data) <- colnames(total_A$data)
    total_B$data
	
	total_B$qtl<-rbind(sires_toC_qtl,dams_toC_qtl)
	total_B$qtl[,2]<-0
	# dim(total_B$qtl)
	# class(total_B$qtl)
	total_B$qtl<-as.matrix(total_B$qtl)
	# class(total_B$qtl)

	total_B$mrk<-rbind(sires_toC_mrk,dams_toC_mrk)
	total_B$mrk[,2]<-0
	# dim(total_B$mrk)
	total_B$mrk<-as.matrix(total_B$mrk)

	# seq
	total_B$sequ<-rbind(sires_toC_seq,dams_toC_seq)
	# dim(total_B$sequ)
	total_B$sequ<-as.matrix(total_B$sequ)
	total_B$sequ[,2]<-0
	
# POPULATION 2 DATA PREPARATION FINISH--------------

# POPULATION CROSS DATA PREPARATION START--------------
cat('Population cross data preparation ...',fill=TRUE)
# SIRE DATA CROSS
	founder_cross[,2]<-as.numeric(founder_cross[,2])
if(founder_cross[1,1]=='pop1'){

	my_internal_data1<-pop1$output[[founder_pop1[1,2]+1]]
	my_internal_data<-my_internal_data1[[1]]
	my_internal_data

} else if (founder_cross[1,1]=='pop2'){
	my_internal_data1<-pop2$output[[founder_pop2[1,2]+1]]
	my_internal_data<-my_internal_data1[[1]]
	my_internal_data
	
}
# SIRE CROSS
# cat('# SIRE CROSS',fill=TRUE)
	provided_sire_data<-subset(my_internal_data,my_internal_data[,5]=='M')
	provided_sire_data
	length(provided_sire_data[,1])
	index_sires<-provided_sire_data[,1]
	index_sires_1<-match(index_sires,my_internal_data1[[1]][,1])
	index_sires_1
	provided_sire_seq<-my_internal_data1[[4]][index_sires_1,]
	# dim(provided_sire_seq)
	provided_sire_qtl<-my_internal_data1[[2]][index_sires_1,]
	# dim(provided_sire_qtl)
	provided_sire_mrk<-my_internal_data1[[3]][index_sires_1,]
	# dim(provided_sire_mrk)

	# sires_toC
	if(founder_cross[1,3]=='rnd'){
	ID_sires_toC<-sample(index_sires,founder_cross[1,2],replace=FALSE)
	index<-match(ID_sires_toC,index_sires)
	} else if(founder_cross[1,3]=='phen'){

		if(founder_cross[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,6]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_cross[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)

		}
		if(founder_cross[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,6]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_cross[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}

	} else if (founder_cross[1,3]=='tbv'){
		if(founder_cross[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,8]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_cross[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
			if(founder_cross[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,8]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_cross[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
		
	} else if (founder_cross[1,3]=='gebv'){
		if(founder_cross[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,9]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_cross[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
			if(founder_cross[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,9]), ]
		ID_sires_toC<-sorted_sire_data[1:founder_cross[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
		
	}
	

	sires_toC_data<-provided_sire_data[index,]
	sires_toC_seq<-provided_sire_seq[index,]
	sires_toC_qtl<-provided_sire_qtl[index,]
	sires_toC_mrk<-provided_sire_mrk[index,]


# DAM DATA CROSS
# cat('# # DAM DATA CROSS',fill=TRUE)

if(founder_cross[2,1]=='pop1'){

	my_internal_data1<-pop1$output[[founder_pop1[2,2]+1]]
	my_internal_data<-my_internal_data1[[1]]
	my_internal_data

} else if (founder_cross[2,1]=='pop2'){
	my_internal_data1<-pop2$output[[founder_pop2[2,2]+1]]
	my_internal_data<-my_internal_data1[[1]]
	my_internal_data
	
}
# DAM CROSS

	provided_dam_data<-subset(my_internal_data,my_internal_data[,5]=='F')
	provided_dam_data
	length(provided_dam_data[,1])
	index_dams<-provided_dam_data[,1]
	index_dams_1<-match(index_dams,my_internal_data1[[1]][,1])
	index_dams_1
	provided_dam_seq<-my_internal_data1[[4]][index_dams_1,]
	# dim(provided_dam_seq)
	provided_dam_qtl<-my_internal_data1[[2]][index_dams_1,]
	# dim(provided_dam_qtl)
	provided_dam_mrk<-my_internal_data1[[3]][index_dams_1,]
	# dim(provided_dam_mrk)

	# dams_toC
	if(founder_cross[2,3]=='rnd'){
	ID_dams_toC<-sample(index_dams,founder_cross[2,2],replace=FALSE)
	index<-match(ID_dams_toC,index_dams)
	} else if(founder_cross[2,3]=='phen'){

		if(founder_cross[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,6]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_cross[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)

		}
		if(founder_cross[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,6]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_cross[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}

	} else if (founder_cross[2,3]=='tbv'){
		if(founder_cross[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,8]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_cross[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
			if(founder_cross[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,8]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_cross[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
		
	} else if (founder_cross[2,3]=='gebv'){
		if(founder_cross[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,9]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_cross[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
			if(founder_cross[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,9]), ]
		ID_dams_toC<-sorted_dam_data[1:founder_cross[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
		
	}

	dams_toC_data<-provided_dam_data[index,]
	dams_toC_seq<-provided_dam_seq[index,]
	# dim(dams_toC_data)
	dams_toC_qtl<-provided_dam_qtl[index,]
	dams_toC_mrk<-provided_dam_mrk[index,]
	# dim(dams_toC_mrk)
	
# cat('# # ghabl az rbind',fill=TRUE)
names(sires_toC_data) <- names(dams_toC_data) 
	start_data<-rbind(sires_toC_data,dams_toC_data)
	total_C<-list()
	total_C$data<-start_data
	total_C$data[,4]<-0
	total_C$data
	# dim(total_C$data)
	# cat('# # bad az rbind data',fill=TRUE)

	total_C$qtl<-rbind(sires_toC_qtl,dams_toC_qtl)
	total_C$qtl[,2]<-0
	# dim(total_C$qtl)
	# class(total_C$qtl)
	total_C$qtl<-as.matrix(total_C$qtl)
	# class(total_C$qtl)
	# cat('# # bad az rbind qtl',fill=TRUE)

	total_C$mrk<-rbind(sires_toC_mrk,dams_toC_mrk)
	total_C$mrk[,2]<-0
	# dim(total_C$mrk)
	total_C$mrk<-as.matrix(total_C$mrk)
	# cat('# # bad az rbind mrk',fill=TRUE)

	# seq
	total_C$sequ<-rbind(sires_toC_seq,dams_toC_seq)
	# dim(total_C$sequ)
	total_C$sequ<-as.matrix(total_C$sequ)
	total_C$sequ[,2]<-0
	# cat('# # bad az rbind sequ',fill=TRUE)



# POPULATION CROSS DATA PREPARATION FINISH--------------

# source('masalan.r')
gene_counter<-0
 for (gene_counter in 0:ng){ 

 start_ge <- Sys.time()
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CROSSBRED START
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# cat('CROSSBRED START ',fill=TRUE)
cat('Generation',gene_counter,'started ')

bb1<-subset(total_C$data,total_C$data[,4]==gene_counter) 
bb2<-subset(total_C$qtl,total_C$qtl[,2]==gene_counter) 
bb3<-subset(total_C$mrk,total_C$mrk[,2]==gene_counter) 
bb4<-subset(total_C$sequ,total_C$sequ[,2]==gene_counter) 


	 # if (gene_counter ==0){ 
	# selection of males and then top males
	males<-subset(bb1,bb1[,5]=='M')
	males_selected<-males
	# selection of females and then top females
	females<-subset(bb1,bb1[,5]=='F')
    females_selected<-females
	x<-length(females_selected[,1])*litter_size
	# }


#mating
# x<-ndam_selected*litter_size
sires<-males_selected[,1]
if (anyNA(sires)==TRUE){
stop("\n",'Number of required sires for crossing is less than available number of male progeny.'
,"\n",'Solution 1: Increase litter size.'
,"\n",'Solution 2: Decrease number of selected sires.')
}
		 
dams<-females_selected[,1]
if (anyNA(dams)==TRUE){
stop("\n",'Number of required dams for crossing is less than available number of female progeny.'
,"\n",'Solution 1: Increase litter size.'
,"\n",'Solution 2: Decrease number of selected dams.')
}

No_off_per_sire<-x/length(sires)
No_off_per_sire
# No_mat_per_sire<-No_off_per_sire/litter_size
# No_mat_per_sire
No_off_per_dam<-litter_size
No_off_per_dam
No_mat_per_dam<-1
No_mat_per_dam	

# Create id,sire,dam,generation,sex,env
# counter_id<-length(total_C$data[,1])
# id1 <- (counter_id+1)
# id2 <- (counter_id+x)
id1 <- 1
id2 <- x
id  <- id1:id2
if(No_off_per_sire==floor(No_off_per_sire)){
sire <- rep(sires,each=No_off_per_sire)
} else
sire<-sample(sires,x,replace=TRUE)
dam <- rep(dams,each=No_off_per_dam)
generation <-gene_counter +1 
sex <-sample(c('F','M'),x,replace=T)
env <-rnorm(x,mean=0,sd=sqrt(vv))

sirem <- bb3[match(sire,bb3[,1]),]
damm <- bb3[match(dam,bb3[,1]),]
qtlseq_sire<-bb2[match(sire,bb2[,1]),]
qtlseq_dam<-bb2[match(dam,bb2[,1]),]
total_seq_sires<-bb4[match(sire,bb4[,1]),]
total_seq_sires<-total_seq_sires[,-c(1,2)]
total_seq_sires<-rbind(pos,total_seq_sires)
locii_sires<-total_seq_sires
total_seq_dams<-bb4[match(dam,bb4[,1]),]
total_seq_dams<-total_seq_dams[,-c(1,2)]
total_seq_dams<-rbind(pos,total_seq_dams)
locii_dams<-total_seq_dams
 # dim(locii_dams)
cat('.')

# ////////////////////////////////////
# COMBINED 1to2 and CROSSOVER FORTRAN///////START////
# ////////////////////////////////////
 
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
		 # dim(new_For_recom_males)
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
		# dim(new_For_recom_males)
		funi2<-function(vec){
		length(which(vec>=0))
		}
		 L1<-rep(1:for_dim[1,1],nchr)
		 L2<-rep(1:nchr,for_dim[1,])
		 L3<-apply(new_For_recom_males,1,funi2)
		 For_recom_males<-data.frame(L1,L2,L3,new_For_recom_males)	
		 # dim(For_recom_males)
		 For_recom_males<-as.matrix(For_recom_males)
		# dim(For_recom_males)
		# For_recom_males[,1:5]
		# max(For_recom_males)
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
		 # dim(new_For_recom_females)
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

		# dim(new_For_recom_females)
		# funi2<-function(vec){
		# length(which(vec>=0))
		# }
		 
		 L1<-rep(1:for_dim[1,1],nchr)
		 L2<-rep(1:nchr,for_dim[1,])
		 L3<-apply(new_For_recom_females,1,funi2)
		 For_recom_females<-data.frame(L1,L2,L3,new_For_recom_females)	
		 # dim(For_recom_females)
		 For_recom_females<-as.matrix(For_recom_females)
		# dim(For_recom_females)
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
 cat('.')
 # pos_total<-pos[seq(1,length(pos),2)]
# After Fortran calcs
lmi_m<-matrix(as.numeric(hp_str[[4]]),nrow = (length(locii_sires_npos[,1])*2),ncol = arg4_m)
tempisirem<-lmi_m
# cat(dim(tempisirem),'dim(tempisirem)',fill=TRUE)

lmi_f<-matrix(as.numeric(hp_str[[6]]),nrow = (length(locii_dams_npos[,1])*2),ncol = arg4_f)
tempidamm<-lmi_f


# ////////////////////////////////////
# COMBINED 1to2 and CROSSOVER FORTRAN///////FINISH////
# ////////////////////////////////////


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
# ----------------------------------

# dada<-matrix(ncol=length(tempisirem[1,]),nrow=x)
# nana<-matrix(ncol=length(tempidamm[1,]),nrow=x)
# w1<-sample(c(1,2),x,replace=TRUE)
# w2<-sample(c(1,2),x,replace=TRUE)
# seq1<-seq(1,x,2)
# seq2<-seq(2,x,2)
# dada[seq1,]<-tempisirem[seq1,]
# dada[seq2,]<-tempisirem[seq2,]
# nana[seq1,]<-tempidamm[seq1,]  
# nana[seq2,]<-tempidamm[seq2,]
# offsire<-dada
# offdam<-nana
# in_data<-data.frame(seq1,seq2)
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
msel<-as.matrix(msel)
cat('.')
# cat('sampling markers for offspring from parents',time.taken,fill=TRUE)

#QTL allele freq
# freq1qtl_C<-Calc_freq(qtlsel1)

# Passing to .Fortran for calc Freq--start
loci_mat<-qtlsel1
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1qtl_C<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish

freq2qtl_C<-1-freq1qtl_C

#Marker allele freq
# freq1mrk_C<-Calc_freq(msel)

# Passing to .Fortran for calc Freq--start
loci_mat<-msel
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1mrk_C<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish

freq2mrk_C<-1-freq1mrk_C



id_C<-id  
sire_C <-sire
dam_C<-dam 
generation_C<-generation 
sex_C<-sex 
env_C<-env 
msel_C<-msel
qtlsel1_C<-qtlsel1
sequence_sel_C<-sequence_sel

# TBV of CROSS
snp_validation<-qtlsel1
 s1<-seq(1,nqtl_allele,2)
 s2<-seq(2,nqtl_allele,2)
 a1<-snp_validation[,s1]+snp_validation[,s2]
 a1[a1==3]=1
 a1[a1==4]=0
 Snp_BreedBB<-a1
freq1<-as.numeric()
freq2<-as.numeric()
for (i in 1:length(Snp_BreedBB[1,])){
freq1[i]<-sum(Snp_BreedBB[,i])/(length(Snp_BreedBB[,i])*2)
}
freq2<-1-freq1
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
tbvp_C<- rowSums(xprogeny)


# TGV of CROSS
tgv4<-calc_TGV(qtlsel1,add_eff_1,dom_eff)
mean(tgv4)
phen_C <- tgv4 + env
phen_C<-as.numeric(phen_C)

cat('.')
# GEBVP CROSS
# TRAINING CROSS START-------------------

if(flag_sel_A==TRUE | flag_sel_B==TRUE){
if(train_type=='crossbred'){

# if(gene_counter ==0){
# # # control of training
	
# }
	

	# if(tra_sel=='rnd'){
	# index_training<-sample(dim(msel_C)[1],tra_size,replace=FALSE)
	# # index_training<-sample(dim(msel_C)[1],200,replace=FALSE)
	# snp_reference1<-msel_C[index_training,]
	# phen_reference<-phen_C[index_training]
	# }
	
			
	if(tra_sel=='rnd'){
	index_training<-sample(dim(msel_C)[1],tra_size,replace=FALSE)
	# index_training<-sample(dim(msel_C)[1],200,replace=FALSE)
	snp_reference<-msel_C[index_training,]
	phen_reference<-phen_C[index_training]
	} else if (tra_sel=='max_rel_mrk'){
		# MRK based
		# Genomic Relationship Matrix
		M_matrix<-msel_C
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
		snp_reference<-msel_C[index_training,]
		phen_reference<-phen_C[index_training]
	} else if (tra_sel=='min_rel_mrk'){
		# MRK based
		# Genomic Relationship Matrix
		M_matrix<-msel_C
		# dim(M_matrix)
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
		snp_reference<-msel_C[index_training,]
		phen_reference<-phen_C[index_training]
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
		snp_reference<-msel_C[index_training,]
		phen_reference<-phen_C[index_training]
	} else if (tra_sel=='min_rel_qtl'){
		# QTL based
		# Genomic Relationship Matrix
		M_matrix<-qtlsel1
		# dim(M_matrix)
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
		snp_reference<-msel_C[index_training,]
		phen_reference<-phen_C[index_training]
	} #'min_rel_qtl'

snp_reference<-bin_snp(snp_reference)
snp_reference<-as.matrix(snp_reference)
y<-as.matrix(phen_reference)


# snp_reference<-bin_snp(snp_reference1)
# snp_reference<-as.matrix(snp_reference)
# y<-as.matrix(phen_reference)

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
		} # end loop else if (addTra==FALSE)


# write.table(ghat,file=fileghatBsc1[ii],row.names=F,col.names=F) 
# write.table(dhat,file=filedhatBsc1[ii],row.names=F,col.names=F) 
# write.table(ghat,file=fileghatAsc1[ii],row.names=F,col.names=F) 
# write.table(dhat,file=filedhatAsc1[ii],row.names=F,col.names=F) 
# TRAINING FINISH-------------------

	} #end loop train_type=='crossbred'

} # end loop flag



# GEBVP CROSS
if(flag_sel_A==TRUE | flag_sel_B==TRUE){
	if(train_type=='crossbred'){
		snp_validation<-msel_C
	s1<-seq(1,nmarker_allele,2)
	s2<-seq(2,nmarker_allele,2)
	a1<-snp_validation[,s1]+snp_validation[,s2]
	a1[a1==3]=1
	a1[a1==4]=0
	dhat<-as.matrix(dhat)
	ghat1<-ghat[,1]
	ghat2<-ghat[,2]
	Snp_BreedB<-a1
	freq1<-as.numeric()
	freq2<-as.numeric()
	for (i in 1:length(Snp_BreedB[1,])){
	freq1[i]<-sum(Snp_BreedB[,i])/(length(Snp_BreedB[,i])*2)
	}
	freq2<-1-freq1
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
	gebvp_C<- rowSums(xprogeny)
	  }
	if(train_type!='crossbred'){
	gebvp_C<-rep(0,length(tbvp_C))
	  }  
    } 
	
	# if(missing(train_type)){
	# gebvp_C<-rep(0,length(tbvp_C))
	# }
	
	if(flag_sel_A==FALSE & flag_sel_B==FALSE){
	gebvp_C<-rep(0,length(tbvp_C))
	}
	
	
# # GEBVP CROSS
# if(train_type=='crossbred'){
# snp_validation<-msel_C
# s1<-seq(1,nmarker_allele,2)
# s2<-seq(2,nmarker_allele,2)
# a1<-snp_validation[,s1]+snp_validation[,s2]
# a1[a1==3]=1
# a1[a1==4]=0

# dhat<-as.matrix(dhat)
# ghat1<-ghat[,1]
# ghat2<-ghat[,2]

# Snp_BreedB<-a1
# freq1<-as.numeric()
# freq2<-as.numeric()
# for (i in 1:length(Snp_BreedB[1,])){
# freq1[i]<-sum(Snp_BreedB[,i])/(length(Snp_BreedB[,i])*2)
# }
# freq2<-1-freq1
# snp_validation<-a1
# xprogeny<-snp_validation
# xprogeny<-as.matrix(xprogeny)
# q1<-xprogeny

# for (i in 1:length(xprogeny[,1])){
# ti<-xprogeny[i,]
# two<-which(ti==2)
# one<-which(ti==1)
# zero<-which(ti==0)
# q1[i,two]<-((freq1[two])*ghat1[two])+
         # ((freq2[two])*dhat[two])
# q1[i,one]<-((1/2*freq1[one])*ghat1[one])+
         # ((1/2*freq2[one])*dhat[one])+
		 # ((1/2*freq1[one])*dhat[one])+
		 # ((1/2*freq2[one])*ghat2[one])
# q1[i,zero]<-((freq2[zero])*ghat2[zero])+
         # ((freq1[zero])*dhat[zero])	
# }
# xprogeny<-q1
# gebvp_C<- rowSums(xprogeny)
# } else {
# gebvp_C<-rep(0,length(tbvp_C))
# }

	newgen_cross <- data.frame(id_C, sire_C, dam_C, gene_counter,sex_C, phen_C,env_C,tbvp_C,gebvp_C)
    gen0<-newgen_cross[,4]
	names(newgen_cross)<-c('ID', 'sire', 'dam', 'generation','sex', 'phen','env','tbvp','gebvp')
	
# QTL
qtlsel_finall<-cbind(id_C,gen0,qtlsel1_C)
# Marker
msel_finall<-cbind(id_C,gen0,msel_C)
# SEQ
sequence_sel_finall<-cbind(id_C,gen0,sequence_sel_C)
sequence_sel_finall<-as.matrix(sequence_sel_finall)
cat('.')
	# FILL IN DATA CROSS
   Total_Total_C[[gene_counter+1]][[1]]<-newgen_cross
   Total_Total_C[[gene_counter+1]][[2]]<-qtlsel_finall
   Total_Total_C[[gene_counter+1]][[3]]<-msel_finall
   Total_Total_C[[gene_counter+1]][[4]]<-sequence_sel_finall
   # Freq QTL
   	ID_qtl<-my_internal_data1_A[[5]][,1]
	ID_qtl_gen<-rep(gene_counter,length(ID_qtl))
	ID_qtl_chr<-my_internal_data1_A[[5]][,3]
	freqQTL<-data.frame(ID_qtl,ID_qtl_gen,ID_qtl_chr,freq1qtl_C,freq2qtl_C)
	names(freqQTL)<-c('ID','Generation','Chr','Freq.Allele1','Freq.Allele2')
   Total_Total_C[[gene_counter+1]][[5]]<-freqQTL

   # Freq MRK
   # outforLD$freqMrk
   	ID_mrk<-my_internal_data1_A[[6]][,1]
	ID_mrk_gen<-rep(gene_counter,length(ID_mrk))
	ID_mrk_chr<-my_internal_data1_A[[6]][,3]
	freqMRK<-data.frame(ID_mrk,ID_mrk_gen,ID_mrk_chr,freq1mrk_C,freq2mrk_C)
	names(freqMRK)<-c('ID','Generation','Chr','Freq.Allele1','Freq.Allele2')
    Total_Total_C[[gene_counter+1]][[6]]<-freqMRK
cat('.')
# cat("\n",'CROSSBRED FINISH',fill=TRUE)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CROSSBRED FINISH
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# BREED A STARTS
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# cat('BREED A STARTS ...',fill=TRUE)
cat('.')
  # ----------------------------
   Total_Total_A[[gene_counter+1]][[1]]<-total_A$data
   Total_Total_A[[gene_counter+1]][[2]]<-total_A$qtl
   Total_Total_A[[gene_counter+1]][[3]]<-total_A$mrk
   Total_Total_A[[gene_counter+1]][[4]]<-total_A$sequ
   # Freq QTL
   	ID_qtl<-my_internal_data1_A[[5]][,1]
	ID_qtl_gen<-rep(gene_counter,length(ID_qtl))
	ID_qtl_chr<-my_internal_data1_A[[5]][,3]
	
	#QTL allele freq
    # freq1qtl_A<-Calc_freq(total_A$qtl[,-c(1,2)])
	
	# Passing to .Fortran for calc Freq--start
loci_mat<-total_A$qtl[,-c(1,2)]
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1qtl_A<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish

	freq2qtl_A<-1-freq1qtl_A
	freqQTL<-data.frame(ID_qtl,ID_qtl_gen,ID_qtl_chr,freq1qtl_A,freq2qtl_A)
	names(freqQTL)<-c('ID','Generation','Chr','Freq.Allele1','Freq.Allele2')
    Total_Total_A[[gene_counter+1]][[5]]<-freqQTL
      
  
  # Freq MRK
   	ID_mrk<-my_internal_data1_A[[6]][,1]
	ID_mrk_gen<-rep(gene_counter,length(ID_mrk))
	ID_mrk_chr<-my_internal_data1_A[[6]][,3]
    # freq1mrk_A<-Calc_freq(total_A$mrk[,-c(1,2)])
	
		# Passing to .Fortran for calc Freq--start
loci_mat<-total_A$mrk[,-c(1,2)]
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1mrk_A<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish

	freq2mrk_A<-1-freq1mrk_A
	freqMRK<-data.frame(ID_mrk,ID_mrk_gen,ID_mrk_chr,freq1mrk_A,freq2mrk_A)
	names(freqMRK)<-c('ID','Generation','Chr','Freq.Allele1','Freq.Allele2')
    Total_Total_A[[gene_counter+1]][[6]]<-freqMRK
	  # ----------------------------

bb1<-subset(total_A$data,total_A$data[,4]==gene_counter) 
bb2<-subset(total_A$qtl,total_A$qtl[,2]==gene_counter) 
bb3<-subset(total_A$mrk,total_A$mrk[,2]==gene_counter) 
bb4<-subset(total_A$sequ,total_A$sequ[,2]==gene_counter) 


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
	
		if(Selection_pop1[1,2]=='rnd'){
		ID_sires<-sample(males[,1],Selection_pop1[1,1],replace=FALSE)
		males<-males[match(ID_sires,males[,1]),]  
		males_selected<-males[1:Selection_pop1[1,1],]		
		}
		
		if(Selection_pop1[1,2]=='phen'){
			if(Selection_pop1[1,3]=='h'){
		males<-males[order(-males[,6]),]  #selection based on phen
		males_selected<-males[1:Selection_pop1[1,1],]
			}
			if(Selection_pop1[1,3]=='l'){
		males<-males[order(males[,6]),]  #selection based on phen
		males_selected<-males[1:Selection_pop1[1,1],]
			}
		}
	
	   if(Selection_pop1[1,2]=='tbv'){ #selection based on tbv
	   		if(Selection_pop1[1,3]=='h'){
		males<-males[order(-males[,8]),]  #selection based on phen
		males_selected<-males[1:Selection_pop1[1,1],]
			}
			if(Selection_pop1[1,3]=='l'){
		males<-males[order(males[,8]),]  #selection based on phen
		males_selected<-males[1:Selection_pop1[1,1],]
			}
    	} 
		
		if(Selection_pop1[1,2]=='gebv'){ #selection based on gebv
	   		if(Selection_pop1[1,3]=='h'){
		males<-males[order(-males[,9]),]  #selection based on phen
		males_selected<-males[1:Selection_pop1[1,1],]
			}
			if(Selection_pop1[1,3]=='l'){
		males<-males[order(males[,9]),]  #selection based on phen
		males_selected<-males[1:Selection_pop1[1,1],]
			}
    	} 
 
 	   if(Selection_pop1[1,2]=='tbvc'){ #selection based on tbvc
	   		if(Selection_pop1[1,3]=='h'){
		males<-males[order(-males[,10]),]  #selection based on phen
		males_selected<-males[1:Selection_pop1[1,1],]
			}
			if(Selection_pop1[1,3]=='l'){
		males<-males[order(males[,10]),]  #selection based on phen
		males_selected<-males[1:Selection_pop1[1,1],]
			}
    	} 
		
		if(Selection_pop1[1,2]=='gebvc'){ #selection based on gebvc
	   		if(Selection_pop1[1,3]=='h'){
		males<-males[order(-males[,11]),]  #selection based on phen
		males_selected<-males[1:Selection_pop1[1,1],]
			}
			if(Selection_pop1[1,3]=='l'){
		males<-males[order(males[,11]),]  #selection based on phen
		males_selected<-males[1:Selection_pop1[1,1],]
			}
    	} 
		

	# selection of females and then top females
	    females<-subset(bb1,bb1[,5]=='F')
	
		if(Selection_pop1[2,2]=='rnd'){
		ID_dams<-sample(females[,1],Selection_pop1[2,1],replace=FALSE)
		females<-females[match(ID_dams,females[,1]),]
		females_selected<-females[1:Selection_pop1[2,1],]
		}
		
		if(Selection_pop1[2,2]=='phen'){
			if(Selection_pop1[2,3]=='h'){
		females<-females[order(-females[,6]),]  #selection based on phen
		females_selected<-females[1:Selection_pop1[2,1],]
			}
			if(Selection_pop1[2,3]=='l'){
		females<-females[order(females[,6]),]  #selection based on phen
		females_selected<-females[1:Selection_pop1[2,1],]
			}
		}
	
	   if(Selection_pop1[2,2]=='tbv'){ #selection based on tbv
	   		if(Selection_pop1[2,3]=='h'){
		females<-females[order(-females[,8]),]  #selection based on phen
		females_selected<-females[1:Selection_pop1[2,1],]
			}
			if(Selection_pop1[2,3]=='l'){
		females<-females[order(females[,8]),]  #selection based on phen
		females_selected<-females[1:Selection_pop1[2,1],]
			}
    	} 
		
		if(Selection_pop1[2,2]=='gebv'){ #selection based on tbv
	   		if(Selection_pop1[2,3]=='h'){
		females<-females[order(-females[,9]),]  #selection based on phen
		females_selected<-females[1:Selection_pop1[2,1],]
			}
			if(Selection_pop1[2,3]=='l'){
		females<-females[order(females[,9]),]  #selection based on phen
		females_selected<-females[1:Selection_pop1[2,1],]
			}
    	} 
		
		if(Selection_pop1[2,2]=='tbvc'){ #selection based on tbvc
	   		if(Selection_pop1[2,3]=='h'){
		females<-females[order(-females[,10]),]  #selection based on phen
		females_selected<-females[1:Selection_pop1[2,1],]
			}
			if(Selection_pop1[2,3]=='l'){
		females<-females[order(females[,10]),]  #selection based on phen
		females_selected<-females[1:Selection_pop1[2,1],]
			}
    	}

		if(Selection_pop1[2,2]=='gebvc'){ #selection based on gebvc
	   		if(Selection_pop1[2,3]=='h'){
		females<-females[order(-females[,11]),]  #selection based on phen
		females_selected<-females[1:Selection_pop1[2,1],]
			}
			if(Selection_pop1[2,3]=='l'){
		females<-females[order(females[,11]),]  #selection based on phen
		females_selected<-females[1:Selection_pop1[2,1],]
			}
    	} 
		x<-Selection_pop1[2,1]*litter_size
}
	

#mating
# x<-ndam_selected*litter_size
sires<-males_selected[,1]
if (anyNA(sires)==TRUE){
stop("\n",'Number of required sires in pop1 (Breed A) is less than available number of male progeny.'
,"\n",'Solution 1: Increase litter size.'
,"\n",'Solution 2: Decrease number of selected sires in argument "Selection_pop1"')
}
		 
dams<-females_selected[,1]
if (anyNA(dams)==TRUE){
stop("\n",'Number of required dams in pop1 (Breed A) is less than available number of female progeny.'
,"\n",'Solution 1: Increase litter size.'
,"\n",'Solution 2: Decrease number of selected dams in argument "Selection_pop1"')
}

No_off_per_sire<-x/length(sires)
No_off_per_sire
# No_mat_per_sire<-No_off_per_sire/litter_size
# No_mat_per_sire
No_off_per_dam<-litter_size
No_off_per_dam
No_mat_per_dam<-1
No_mat_per_dam	

# Create id,sire,dam,generation,sex,env
# counter_id<-length(total_A$data[,1])
# id1 <- (counter_id+1)
# id2 <- (counter_id+x)
# id  <- id1:id2
	my<-c()
	for (klm in 1:(gene_counter+1)){
	my[klm]<-length(Total_Total_A[[klm]][[1]][,1])
	}
counter_id<- sum(my)
id1 <- (counter_id+1)
id2 <- (counter_id+x)
id  <- id1:id2
if(No_off_per_sire==floor(No_off_per_sire)){
sire <- rep(sires,each=No_off_per_sire)
} else
sire<-sample(sires,x,replace=TRUE)
dam <- rep(dams,each=No_off_per_dam)
generation <-gene_counter +1 
sex <-sample(c('F','M'),x,replace=T)
env <-rnorm(x,mean=0,sd=sqrt(vv))

sirem <- bb3[match(sire,bb3[,1]),]
damm <- bb3[match(dam,bb3[,1]),]
qtlseq_sire<-bb2[match(sire,bb2[,1]),]
qtlseq_dam<-bb2[match(dam,bb2[,1]),]
total_seq_sires<-bb4[match(sire,bb4[,1]),]
total_seq_sires<-total_seq_sires[,-c(1,2)]
total_seq_sires<-rbind(pos,total_seq_sires)
locii_sires<-total_seq_sires
total_seq_dams<-bb4[match(dam,bb4[,1]),]
total_seq_dams<-total_seq_dams[,-c(1,2)]
total_seq_dams<-rbind(pos,total_seq_dams)
locii_dams<-total_seq_dams
 # dim(locii_dams)


# ////////////////////////////////////
# COMBINED 1to2 and CROSSOVER FORTRAN///////START////
# ////////////////////////////////////
 
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
		 # dim(new_For_recom_males)
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
		# dim(new_For_recom_males)
		# funi2<-function(vec){
		# length(which(vec>=0))
		# }
		 L1<-rep(1:for_dim[1,1],nchr)
		 L2<-rep(1:nchr,for_dim[1,])
		 L3<-apply(new_For_recom_males,1,funi2)
		 For_recom_males<-data.frame(L1,L2,L3,new_For_recom_males)	
		 # dim(For_recom_males)
		 For_recom_males<-as.matrix(For_recom_males)
		# dim(For_recom_males)

		# max(For_recom_males)
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
		 # dim(new_For_recom_females)
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

		# dim(new_For_recom_females)
		# funi2<-function(vec){
		# length(which(vec>=0))
		# }
		 
		 L1<-rep(1:for_dim[1,1],nchr)
		 L2<-rep(1:nchr,for_dim[1,])
		 L3<-apply(new_For_recom_females,1,funi2)
		 For_recom_females<-data.frame(L1,L2,L3,new_For_recom_females)	
		 # dim(For_recom_females)
		 For_recom_females<-as.matrix(For_recom_females)
		# dim(For_recom_females)
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
dim(tempidamm)

# ////////////////////////////////////
# COMBINED 1to2 and CROSSOVER FORTRAN///////FINISH////
# ////////////////////////////////////


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
# ----------------------------------

# dada<-matrix(ncol=length(tempisirem[1,]),nrow=x)
# nana<-matrix(ncol=length(tempidamm[1,]),nrow=x)
# w1<-sample(c(1,2),x,replace=TRUE)
# w2<-sample(c(1,2),x,replace=TRUE)
# seq1<-seq(1,x,2)
# seq2<-seq(2,x,2)
# dada[seq1,]<-tempisirem[seq1,]
# dada[seq2,]<-tempisirem[seq2,]
# nana[seq1,]<-tempidamm[seq1,]  
# nana[seq2,]<-tempidamm[seq2,]
# offsire<-dada
# offdam<-nana
# in_data<-data.frame(seq1,seq2)
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
msel<-as.matrix(msel)

# cat('sampling markers for offspring from parents',time.taken,fill=TRUE)

#QTL allele freq
# freq1qtl_A<-Calc_freq(qtlsel1)

# Passing to .Fortran for calc Freq--start
loci_mat<-qtlsel1
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1qtl_A<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish




freq2qtl_A<-1-freq1qtl_A
# freq_qtl<-data.frame(freq1qtl_A,freq2qtl_A)
# freq_qtl

#Marker allele freq
# freq1mrk_A<-Calc_freq(msel)

# Passing to .Fortran for calc Freq--start
loci_mat<-msel
r_size<-dim(loci_mat)[1]
c_size<-dim(loci_mat)[2]
nloci_freq<-dim(loci_mat)[2]/2
loci_mat<-as.integer(unlist(loci_mat))
outfreq<-.Fortran("cf",r_size= as.integer(r_size),c_size = as.integer(c_size),loci_mat=loci_mat,nloci=as.integer(nloci_freq),freq1_main=double(nloci_freq))
freq1mrk_A<-as.numeric(outfreq[[5]])
# Passing to .Fortran for calc Freq--finish

freq2mrk_A<-1-freq1mrk_A
# freq_mrk<-data.frame(freq1mrk_A,freq2mrk_A)
# freq_mrk
# # dim(freq_mrk)


id_A<-id  
sire_A <-sire
dam_A<-dam 
generation_A<-generation 
sex_A<-sex 
env_A<-env 
msel_A<-msel
qtlsel1_A<-qtlsel1
sequence_sel_A<-sequence_sel

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# BREED B STARTS
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# cat('BREED B STARTS ...',fill=TRUE)
cat('.')
  # ----------------------------
   Total_Total_B[[gene_counter+1]][[1]]<-total_B$data
   Total_Total_B[[gene_counter+1]][[2]]<-total_B$qtl
   Total_Total_B[[gene_counter+1]][[3]]<-total_B$mrk
   Total_Total_B[[gene_counter+1]][[4]]<-total_B$sequ
   # Freq QTL
   	ID_qtl<-my_internal_data1_B[[5]][,1]
	ID_qtl_gen<-rep(gene_counter,length(ID_qtl))
	ID_qtl_chr<-my_internal_data1_B[[5]][,3]
	
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
   Total_Total_B[[gene_counter+1]][[5]]<-freqQTL
      
  
  # Freq MRK
   	ID_mrk<-my_internal_data1_B[[6]][,1]
	ID_mrk_gen<-rep(gene_counter,length(ID_mrk))
	ID_mrk_chr<-my_internal_data1_B[[6]][,3]
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
	
	
	freq2mrk_B<-1-freq1mrk_B
	freqMRK<-data.frame(ID_mrk,ID_mrk_gen,ID_mrk_chr,freq1mrk_B,freq2mrk_B)
	names(freqMRK)<-c('ID','Generation','Chr','Freq.Allele1','Freq.Allele2')
    Total_Total_B[[gene_counter+1]][[6]]<-freqMRK
	  # ----------------------------


bb1<-subset(total_B$data,total_B$data[,4]==gene_counter) 
bb2<-subset(total_B$qtl,total_B$qtl[,2]==gene_counter) 
bb3<-subset(total_B$mrk,total_B$mrk[,2]==gene_counter) 
bb4<-subset(total_B$sequ,total_B$sequ[,2]==gene_counter) 


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
	
		if(Selection_pop2[1,2]=='rnd'){
		ID_sires<-sample(males[,1],Selection_pop2[1,1],replace=FALSE)
		males<-males[match(ID_sires,males[,1]),]  
		males_selected<-males[1:Selection_pop2[1,1],]		
		}
		
		if(Selection_pop2[1,2]=='phen'){
			if(Selection_pop2[1,3]=='h'){
		males<-males[order(-males[,6]),]  #selection based on phen
		males_selected<-males[1:Selection_pop2[1,1],]
			}
			if(Selection_pop2[1,3]=='l'){
		males<-males[order(males[,6]),]  #selection based on phen
		males_selected<-males[1:Selection_pop2[1,1],]
			}
		}
	
	   if(Selection_pop2[1,2]=='tbv'){ #selection based on tbv
	   		if(Selection_pop2[1,3]=='h'){
		males<-males[order(-males[,8]),]  #selection based on phen
		males_selected<-males[1:Selection_pop2[1,1],]
			}
			if(Selection_pop2[1,3]=='l'){
		males<-males[order(males[,8]),]  #selection based on phen
		males_selected<-males[1:Selection_pop2[1,1],]
			}
    	} 
		
		if(Selection_pop2[1,2]=='gebv'){ #selection based on tbv
	   		if(Selection_pop2[1,3]=='h'){
		males<-males[order(-males[,9]),]  #selection based on phen
		males_selected<-males[1:Selection_pop2[1,1],]
			}
			if(Selection_pop2[1,3]=='l'){
		males<-males[order(males[,9]),]  #selection based on phen
		males_selected<-males[1:Selection_pop2[1,1],]
			}
    	}
		
	   if(Selection_pop2[1,2]=='tbvc'){ #selection based on tbvc
	   		if(Selection_pop2[1,3]=='h'){
		males<-males[order(-males[,10]),]  #selection based on phen
		males_selected<-males[1:Selection_pop2[1,1],]
			}
			if(Selection_pop2[1,3]=='l'){
		males<-males[order(males[,10]),]  #selection based on phen
		males_selected<-males[1:Selection_pop2[1,1],]
			}
    	}

		if(Selection_pop2[1,2]=='gebvc'){ #selection based on gebvc
	   		if(Selection_pop2[1,3]=='h'){
		males<-males[order(-males[,11]),]  #selection based on phen
		males_selected<-males[1:Selection_pop2[1,1],]
			}
			if(Selection_pop2[1,3]=='l'){
		males<-males[order(males[,11]),]  #selection based on phen
		males_selected<-males[1:Selection_pop2[1,1],]
			}
    	}
		
		

	# selection of females and then top females
	    females<-subset(bb1,bb1[,5]=='F')
	
		if(Selection_pop2[2,2]=='rnd'){
		ID_dams<-sample(females[,1],Selection_pop2[2,1],replace=FALSE)
		females<-females[match(ID_dams,females[,1]),]
		females_selected<-females[1:Selection_pop2[2,1],]
		}
		
		if(Selection_pop2[2,2]=='phen'){
			if(Selection_pop2[2,3]=='h'){
		females<-females[order(-females[,6]),]  #selection based on phen
		females_selected<-females[1:Selection_pop2[2,1],]
			}
			if(Selection_pop2[2,3]=='l'){
		females<-females[order(females[,6]),]  #selection based on phen
		females_selected<-females[1:Selection_pop2[2,1],]
			}
		}
	
	   if(Selection_pop2[2,2]=='tbv'){ #selection based on tbv
	   		if(Selection_pop2[2,3]=='h'){
		females<-females[order(-females[,8]),]  #selection based on phen
		females_selected<-females[1:Selection_pop2[2,1],]
			}
			if(Selection_pop2[2,3]=='l'){
		females<-females[order(females[,8]),]  #selection based on phen
		females_selected<-females[1:Selection_pop2[2,1],]
			}
    	} 
		
		if(Selection_pop2[2,2]=='gebv'){ #selection based on tbv
	   		if(Selection_pop2[2,3]=='h'){
		females<-females[order(-females[,9]),]  #selection based on phen
		females_selected<-females[1:Selection_pop2[2,1],]
			}
			if(Selection_pop2[2,3]=='l'){
		females<-females[order(females[,9]),]  #selection based on phen
		females_selected<-females[1:Selection_pop2[2,1],]
			}
    	} 
		
		if(Selection_pop2[2,2]=='tbvc'){ #selection based on tbvc
	   		if(Selection_pop2[2,3]=='h'){
		females<-females[order(-females[,10]),]  #selection based on phen
		females_selected<-females[1:Selection_pop2[2,1],]
			}
			if(Selection_pop2[2,3]=='l'){
		females<-females[order(females[,10]),]  #selection based on phen
		females_selected<-females[1:Selection_pop2[2,1],]
			}
    	} 
		
		if(Selection_pop2[2,2]=='gebvc'){ #selection based on gebvc
	   		if(Selection_pop2[2,3]=='h'){
		females<-females[order(-females[,11]),]  #selection based on phen
		females_selected<-females[1:Selection_pop2[2,1],]
			}
			if(Selection_pop2[2,3]=='l'){
		females<-females[order(females[,11]),]  #selection based on phen
		females_selected<-females[1:Selection_pop2[2,1],]
			}
    	} 
		
		
		
		x<-Selection_pop2[2,1]*litter_size
}

	

#mating
# x<-ndam_selected*litter_size
sires<-males_selected[,1]
if (anyNA(sires)==TRUE){
stop("\n",'Number of required sires in pop2 (Breed B) is less than available number of male progeny.'
,"\n",'Solution 1: Increase litter size.'
,"\n",'Solution 2: Decrease number of selected sires in argument "Selection_pop2"')
}
		 
dams<-females_selected[,1]
if (anyNA(dams)==TRUE){
stop("\n",'Number of required dams in pop2 (Breed B) is less than available number of female progeny.'
,"\n",'Solution 1: Increase litter size.'
,"\n",'Solution 2: Decrease number of selected dams in argument "Selection_pop2"')
}

No_off_per_sire<-x/length(sires)
No_off_per_sire
# No_mat_per_sire<-No_off_per_sire/litter_size
# No_mat_per_sire
No_off_per_dam<-litter_size
No_off_per_dam
No_mat_per_dam<-1
No_mat_per_dam	

# Create id,sire,dam,generation,sex,env
# counter_id<-length(total_B$data[,1])
# id1 <- (counter_id+1)
# id2 <- (counter_id+x)
# id  <- id1:id2

	my<-c()
	for (klm in 1:(gene_counter+1)){
	my[klm]<-length(Total_Total_B[[klm]][[1]][,1])
	}
counter_id<- sum(my)
id1 <- (counter_id+1)
id2 <- (counter_id+x)
id  <- id1:id2
if(No_off_per_sire==floor(No_off_per_sire)){
sire <- rep(sires,each=No_off_per_sire)
} else
sire<-sample(sires,x,replace=TRUE)
dam <- rep(dams,each=No_off_per_dam)
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
 # dim(locii_dams)
# ////////////////////////////////////
# COMBINED 1to2 and CROSSOVER FORTRAN///////START////
# ////////////////////////////////////

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
		 # dim(new_For_recom_males)
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
		# dim(new_For_recom_males)
		# funi2<-function(vec){
		# length(which(vec>=0))
		# }
		 L1<-rep(1:for_dim[1,1],nchr)
		 L2<-rep(1:nchr,for_dim[1,])
		 L3<-apply(new_For_recom_males,1,funi2)
		 For_recom_males<-data.frame(L1,L2,L3,new_For_recom_males)	
		 # dim(For_recom_males)
		 For_recom_males<-as.matrix(For_recom_males)
		# dim(For_recom_males)
		# For_recom_males[,1:5]
		# max(For_recom_males)
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
		 # dim(new_For_recom_females)
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

		# dim(new_For_recom_females)
		# funi2<-function(vec){
		# length(which(vec>=0))
		# }
		 
		 L1<-rep(1:for_dim[1,1],nchr)
		 L2<-rep(1:nchr,for_dim[1,])
		 L3<-apply(new_For_recom_females,1,funi2)
		 For_recom_females<-data.frame(L1,L2,L3,new_For_recom_females)	
		 # dim(For_recom_females)
		 For_recom_females<-as.matrix(For_recom_females)
		# dim(For_recom_females)
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
# ////////////////////////////////////
# COMBINED 1to2 and CROSSOVER FORTRAN///////FINISH////
# ////////////////////////////////////
 
#sampling markers for offspring from parents

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
# ----------------------------------

# dada<-matrix(ncol=length(tempisirem[1,]),nrow=x)
# nana<-matrix(ncol=length(tempidamm[1,]),nrow=x)
# w1<-sample(c(1,2),x,replace=TRUE)
# w2<-sample(c(1,2),x,replace=TRUE)
# seq1<-seq(1,x,2)
# seq2<-seq(2,x,2)
# dada[seq1,]<-tempisirem[seq1,]
# dada[seq2,]<-tempisirem[seq2,]
# nana[seq1,]<-tempidamm[seq1,]  
# nana[seq2,]<-tempidamm[seq2,]
# offsire<-dada
# offdam<-nana
# in_data<-data.frame(seq1,seq2)
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
msel<-as.matrix(msel)

# TBVC of breed B 
freq1<-freq1qtl_A
freq2<-freq2qtl_A
Snp_BreedB<-bin_snp(qtlsel1)
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
tbvc_B<- rowSums(xprogeny)

# TGV of breed B
tgv4<-calc_TGV(qtlsel1,add_eff_1,dom_eff)
mean(tgv4)
phen_B <- tgv4 + env
phen_B<-as.numeric(phen_B)

# TBVP of breed B
snp_validation<-qtlsel1
 s1<-seq(1,nqtl_allele,2)
 s2<-seq(2,nqtl_allele,2)
 a1<-snp_validation[,s1]+snp_validation[,s2]
 a1[a1==3]=1
 a1[a1==4]=0
 Snp_BreedBB<-a1
freq1<-as.numeric()
freq2<-as.numeric()
for (i in 1:length(Snp_BreedBB[1,])){
freq1[i]<-sum(Snp_BreedBB[,i])/(length(Snp_BreedBB[,i])*2)
}
freq2<-1-freq1
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
tbvp_B<- rowSums(xprogeny)
cor(tbvp_B,tbvc_B)

# GEBVP and GEBVC BREED B
# TRAINING BREED B START-------------------

if(flag_sel_B==TRUE){
if(train_type=='purebred'){

	# if(gene_counter ==0){

# }
	

	# if(tra_sel=='rnd'){
	# index_training<-sample(dim(msel)[1],tra_size,replace=FALSE)
	# # index_training<-sample(dim(msel)[1],200,replace=FALSE)
	# snp_reference1<-msel[index_training,]
	# phen_reference<-phen_B[index_training]
	# }
	
				
	if(tra_sel=='rnd'){
	index_training<-sample(dim(msel)[1],tra_size,replace=FALSE)
	snp_reference<-msel[index_training,]
	phen_reference<-phen_B[index_training]
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
		phen_reference<-phen_B[index_training]
	} else if (tra_sel=='min_rel_mrk'){
		# MRK based
		# Genomic Relationship Matrix
		M_matrix<-msel
		# dim(M_matrix)
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
		phen_reference<-phen_B[index_training]
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
		phen_reference<-phen_B[index_training]
	} else if (tra_sel=='min_rel_qtl'){
		# QTL based
		# Genomic Relationship Matrix
		M_matrix<-qtlsel1
		# dim(M_matrix)
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
		phen_reference<-phen_B[index_training]
	} #'min_rel_qtl'

snp_reference<-bin_snp(snp_reference)
snp_reference<-as.matrix(snp_reference)
y<-as.matrix(phen_reference)


# snp_reference<-bin_snp(snp_reference1)
# snp_reference<-as.matrix(snp_reference)
# y<-as.matrix(phen_reference)
 
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
		} # end loop else if (addTra==FALSE)



# write.table(ghat,file=fileghatBsc1[ii],row.names=F,col.names=F) 
# write.table(dhat,file=filedhatBsc1[ii],row.names=F,col.names=F) 
# write.table(ghat,file=fileghatAsc1[ii],row.names=F,col.names=F) 
# write.table(dhat,file=filedhatAsc1[ii],row.names=F,col.names=F) 
# TRAINING FINISH-------------------
} #end loop train_type=='purebred'

} #end loop flag B

# GEBVC BREED B
if(flag_sel_B==TRUE){
freq1<-freq1mrk_A
freq2<-freq2mrk_A
snp_validation<-msel
s1<-seq(1,nmarker_allele,2)
s2<-seq(2,nmarker_allele,2)
a1<-snp_validation[,s1]+snp_validation[,s2]
a1[a1==3]=1
a1[a1==4]=0

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
gebvc_B<- rowSums(xprogeny)
} else {
gebvc_B<-rep(0,length(tbvp_B))
}

# GEBVP BREED B
if(flag_sel_B==TRUE){
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

Snp_BreedB<-a1
freq1<-as.numeric()
freq2<-as.numeric()
for (i in 1:length(Snp_BreedB[1,])){
freq1[i]<-sum(Snp_BreedB[,i])/(length(Snp_BreedB[,i])*2)
}
freq2<-1-freq1
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
gebvp_B<- rowSums(xprogeny)
} else {
gebvp_B<-rep(0,length(tbvp_B))
}
# cor(gebvp_B,gebvc_B)

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

id_B<-id  
sire_B <-sire
dam_B<-dam 
generation_B<-generation 
sex_B<-sex 
env_B<-env 
msel_B<-msel
qtlsel1_B<-qtlsel1
sequence_sel_B<-sequence_sel

# # FILL IN DATA BREED B
   # Total_Total_B[[gene_counter+1]][[1]]<-total_B$data
   # Total_Total_B[[gene_counter+1]][[2]]<-total_B$qtl
   # Total_Total_B[[gene_counter+1]][[3]]<-total_B$mrk
   # Total_Total_B[[gene_counter+1]][[4]]<-total_B$sequ
   # # Freq QTL
   	# ID_qtl<-my_internal_data1_B[[5]][,1]
	# ID_qtl_gen<-rep(gene_counter,length(ID_qtl))
	# ID_qtl_chr<-my_internal_data1_B[[5]][,3]
	# freqQTL<-data.frame(ID_qtl,ID_qtl_gen,ID_qtl_chr,freq1qtl_B,freq2qtl_B)
	# names(freqQTL)<-c('ID','Generation','Chr','Freq.Allele1','Freq.Allele2')
   # Total_Total_B[[gene_counter+1]][[5]]<-freqQTL

   # # Freq MRK
   # # outforLD$freqMrk
   	# ID_mrk<-my_internal_data1_B[[6]][,1]
	# ID_mrk_gen<-rep(gene_counter,length(ID_mrk))
	# ID_mrk_chr<-my_internal_data1_B[[6]][,3]
	# freqMRK<-data.frame(ID_mrk,ID_mrk_gen,ID_mrk_chr,freq1mrk_B,freq2mrk_B)
	# names(freqMRK)<-c('ID','Generation','Chr','Freq.Allele1','Freq.Allele2')
    # Total_Total_B[[gene_counter+1]][[6]]<-freqMRK
	
	
   # DATA OF Next GENERATION BREED B
newgen_B <- data.frame(id_B, sire_B, dam_B, generation_B,sex_B, phen_B,env_B,tbvp_B,gebvp_B,tbvc_B,gebvc_B)
names(newgen_B)<-c('ID', 'sire', 'dam', 'generation','sex', 'phen','env','tbv','gebv','tbvc','gebvc')

# DATA
total_B$data <-newgen_B

# QTL
qtlsel_finall<-cbind(id_B,generation,qtlsel1_B)
total_B$qtl<-qtlsel_finall

# Marker
msel_finall<-cbind(id_B,generation,msel_B)
total_B$mrk<-msel_finall

# SEQ
sequence_sel_finall<-cbind(id_B,generation,sequence_sel_B)
sequence_sel_finall<-as.matrix(sequence_sel_finall)
total_B$sequ<-sequence_sel_finall

# cat('BREED B FINISH ...',fill=TRUE)
cat('.')
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# BREED B FINISH
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# BREED A RESTART
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id<-id_A  
sire<-sire_A 
dam<-dam_A 
generation<-generation_A  
sex<-sex_A  
env<-env_A  
msel<-msel_A 
qtlsel1<-qtlsel1_A 


# TBVC of breed A
freq1<-freq1qtl_B
freq2<-freq2qtl_B
Snp_BreedB<-bin_snp(qtlsel1)
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
tbvc_A<- rowSums(xprogeny)

# TGV of breed A
tgv4<-calc_TGV(qtlsel1,add_eff_1,dom_eff)
mean(tgv4)
phen_A <- tgv4 + env
phen_A<-as.numeric(phen_A)

# TBVP of breed A
snp_validation<-qtlsel1
 s1<-seq(1,nqtl_allele,2)
 s2<-seq(2,nqtl_allele,2)
 a1<-snp_validation[,s1]+snp_validation[,s2]
 a1[a1==3]=1
 a1[a1==4]=0
 Snp_BreedBB<-a1
freq1<-as.numeric()
freq2<-as.numeric()
for (i in 1:length(Snp_BreedBB[1,])){
freq1[i]<-sum(Snp_BreedBB[,i])/(length(Snp_BreedBB[,i])*2)
}
freq2<-1-freq1
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
tbvp_A<- rowSums(xprogeny)
cor(tbvp_A,tbvc_A)

# GEBVP and GEBVC BREED A
# TRAINING BREED A START-------------------

if(flag_sel_A==TRUE){
if(train_type=='purebred'){

	# if(gene_counter ==0){
		
# }
	

	# if(tra_sel=='rnd'){
	# index_training<-sample(dim(msel)[1],tra_size,replace=FALSE)
	# # index_training<-sample(dim(msel)[1],200,replace=FALSE)
	# snp_reference1<-msel[index_training,]
	# phen_reference<-phen_A[index_training]
	# }
			
	if(tra_sel=='rnd'){
	index_training<-sample(dim(msel)[1],tra_size,replace=FALSE)
	snp_reference<-msel[index_training,]
	phen_reference<-phen_A[index_training]
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
		phen_reference<-phen_A[index_training]
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
		phen_reference<-phen_A[index_training]
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
		phen_reference<-phen_A[index_training]
	} else if (tra_sel=='min_rel_qtl'){
		# QTL based
		# Genomic Relationship Matrix
		M_matrix<-qtlsel1
		# dim(M_matrix)
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
		phen_reference<-phen_A[index_training]
	} #'min_rel_qtl'

snp_reference<-bin_snp(snp_reference)
snp_reference<-as.matrix(snp_reference)
y<-as.matrix(phen_reference)

# snp_reference<-bin_snp(snp_reference1)
# snp_reference<-as.matrix(snp_reference)
# y<-as.matrix(phen_reference)
 
 
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
		} # end loop else if (addTra==FALSE)


 

# write.table(ghat,file=fileghatBsc1[ii],row.names=F,col.names=F) 
# write.table(dhat,file=filedhatBsc1[ii],row.names=F,col.names=F) 
# write.table(ghat,file=fileghatAsc1[ii],row.names=F,col.names=F) 
# write.table(dhat,file=filedhatAsc1[ii],row.names=F,col.names=F) 
# TRAINING FINISH-------------------
}

} # end loop flag

# GEBVC BREED A
if(flag_sel_A==TRUE){
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
gebvc_A<- rowSums(xprogeny)
} else {
gebvc_A<-rep(0,length(tbvp_A))
}

# GEBVP BREED A
if(flag_sel_A==TRUE){
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

Snp_BreedB<-a1
freq1<-as.numeric()
freq2<-as.numeric()
for (i in 1:length(Snp_BreedB[1,])){
freq1[i]<-sum(Snp_BreedB[,i])/(length(Snp_BreedB[,i])*2)
}
freq2<-1-freq1
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
gebvp_A<- rowSums(xprogeny)
} else {
gebvp_A<-rep(0,length(tbvp_A))
}
# cor(gebvp_A,gebvc_A)

# # FILL IN DATA BREED A
   # Total_Total_A[[gene_counter+1]][[1]]<-total_A$data
   # Total_Total_A[[gene_counter+1]][[2]]<-total_A$qtl
   # Total_Total_A[[gene_counter+1]][[3]]<-total_A$mrk
   # Total_Total_A[[gene_counter+1]][[4]]<-total_A$sequ
   # # Freq QTL
   	# ID_qtl<-my_internal_data1_A[[5]][,1]
	# ID_qtl_gen<-rep(gene_counter,length(ID_qtl))
	# ID_qtl_chr<-my_internal_data1_A[[5]][,3]
	# freqQTL<-data.frame(ID_qtl,ID_qtl_gen,ID_qtl_chr,freq1qtl_A,freq2qtl_A)
	# names(freqQTL)<-c('ID','Generation','Chr','Freq.Allele1','Freq.Allele2')
   # Total_Total_A[[gene_counter+1]][[5]]<-freqQTL

   # # Freq MRK
   # # outforLD$freqMrk
   	# ID_mrk<-my_internal_data1_A[[6]][,1]
	# ID_mrk_gen<-rep(gene_counter,length(ID_mrk))
	# ID_mrk_chr<-my_internal_data1_A[[6]][,3]
	# freqMRK<-data.frame(ID_mrk,ID_mrk_gen,ID_mrk_chr,freq1mrk_A,freq2mrk_A)
	# names(freqMRK)<-c('ID','Generation','Chr','Freq.Allele1','Freq.Allele2')
    # Total_Total_A[[gene_counter+1]][[6]]<-freqMRK
	
	   # DATA OF Next GENERATION BREED A
newgen_A <- data.frame(id_A, sire_A, dam_A, generation_A,sex_A, phen_A,env_A,tbvp_A,gebvp_A,tbvc_A,gebvc_A)
names(newgen_A)<-c('ID', 'sire', 'dam', 'generation','sex', 'phen','env','tbv','gebv','tbvc','gebvc')

# DATA
total_A$data <-newgen_A

# QTL
qtlsel_finall<-cbind(id_A,generation,qtlsel1_A)
total_A$qtl<-qtlsel_finall

# Marker
msel_finall<-cbind(id_A,generation,msel_A)
total_A$mrk<-msel_finall

# SEQ
sequence_sel_finall<-cbind(id_A,generation,sequence_sel_A)
sequence_sel_finall<-as.matrix(sequence_sel_finall)
total_A$sequ<-sequence_sel_finall
# cat('BREED A FINISH ...',fill=TRUE)
cat('.')
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# BREED A FINISH
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SELECT FROM BREED A AND B AS PARENTS OF CROSBREDS--START
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# cat('Select from breed A and B as parents of crosbreds ...',fill=TRUE)
cat('.','\n')
# SIRE DATA CROSS
Cross_design[,2]<-as.numeric(Cross_design[,2])

if(Cross_design[1,1]=='pop1'){

	my_internal_data1<-total_A$data
	my_internal_data<-my_internal_data1
	my_internal_data
	
	# SIRE CROSS
	provided_sire_data<-subset(my_internal_data,my_internal_data[,5]=='M')
	provided_sire_data
	length(provided_sire_data[,1])
	index_sires<-provided_sire_data[,1]
	index_sires_1<-match(index_sires,my_internal_data1[,1])
	index_sires_1
	
	provided_sire_seq<-total_A$sequ[index_sires_1,]
	# dim(provided_sire_seq)
	provided_sire_qtl<-total_A$qtl[index_sires_1,]
	# dim(provided_sire_qtl)
	provided_sire_mrk<-total_A$mrk[index_sires_1,]
	# dim(provided_sire_mrk)


	

} else if (Cross_design[1,1]=='pop2'){
	my_internal_data1<-total_B$data
	my_internal_data<-my_internal_data1
	my_internal_data
	
	# SIRE CROSS
	provided_sire_data<-subset(my_internal_data,my_internal_data[,5]=='M')
	provided_sire_data
	length(provided_sire_data[,1])
	index_sires<-provided_sire_data[,1]
	index_sires_1<-match(index_sires,my_internal_data1[,1])
	index_sires_1
	
	provided_sire_seq<-total_B$sequ[index_sires_1,]
	# dim(provided_sire_seq)
	provided_sire_qtl<-total_B$qtl[index_sires_1,]
	# dim(provided_sire_qtl)
	provided_sire_mrk<-total_B$mrk[index_sires_1,]
	# dim(provided_sire_mrk)
	
}
#tbv,gebv,phengebvc,tbvc
	# sires_toC
	if(Cross_design[1,3]=='rnd'){
	ID_sires_toC<-sample(index_sires,Cross_design[1,2],replace=FALSE)
	index<-match(ID_sires_toC,index_sires)
	} else if(Cross_design[1,3]=='phen'){

		if(Cross_design[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,6]), ]
		ID_sires_toC<-sorted_sire_data[1:Cross_design[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)

		}
		if(Cross_design[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,6]), ]
		ID_sires_toC<-sorted_sire_data[1:Cross_design[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}

	} else if (Cross_design[1,3]=='tbv'){
		if(Cross_design[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,8]), ]
		ID_sires_toC<-sorted_sire_data[1:Cross_design[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
			if(Cross_design[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,8]), ]
		ID_sires_toC<-sorted_sire_data[1:Cross_design[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
		
	} else if (Cross_design[1,3]=='gebv'){
		if(Cross_design[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,9]), ]
		ID_sires_toC<-sorted_sire_data[1:Cross_design[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
			if(Cross_design[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,9]), ]
		ID_sires_toC<-sorted_sire_data[1:Cross_design[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
		
	} else if (Cross_design[1,3]=='tbvc'){
		if(Cross_design[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,10]), ]
		ID_sires_toC<-sorted_sire_data[1:Cross_design[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
			if(Cross_design[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,10]), ]
		ID_sires_toC<-sorted_sire_data[1:Cross_design[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
		
	} else if (Cross_design[1,3]=='gebvc'){
		if(Cross_design[1,4]=='h'){
		sorted_sire_data<-provided_sire_data[order(-provided_sire_data[,11]), ]
		ID_sires_toC<-sorted_sire_data[1:Cross_design[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
			if(Cross_design[1,4]=='l'){
		sorted_sire_data<-provided_sire_data[order(provided_sire_data[,11]), ]
		ID_sires_toC<-sorted_sire_data[1:Cross_design[1,2],]
		index<-match(ID_sires_toC[,1],index_sires)
		}
		
	}

	sires_toC_data<-provided_sire_data[index,]
	sires_toC_seq<-provided_sire_seq[index,]
	sires_toC_qtl<-provided_sire_qtl[index,]
	sires_toC_mrk<-provided_sire_mrk[index,]


# DAM DATA CROSS

if(Cross_design[2,1]=='pop1'){

    my_internal_data1<-total_A$data
	my_internal_data<-my_internal_data1
	my_internal_data
	# DAM CROSS
	provided_dam_data<-subset(my_internal_data,my_internal_data[,5]=='F')
	provided_dam_data
	length(provided_dam_data[,1])
	index_dams<-provided_dam_data[,1]
	index_dams_1<-match(index_dams,my_internal_data1[,1])
	index_dams_1
	provided_dam_seq<-total_A$sequ[index_dams_1,]
	# dim(provided_dam_seq)
	provided_dam_qtl<-total_A$qtl[index_dams_1,]
	# dim(provided_dam_qtl)
	provided_dam_mrk<-total_A$mrk[index_dams_1,]
	# dim(provided_dam_mrk)
	
	

} else if (Cross_design[2,1]=='pop2'){
	my_internal_data1<-total_B$data
	my_internal_data<-my_internal_data1
	my_internal_data
	
	# DAM CROSS
	provided_dam_data<-subset(my_internal_data,my_internal_data[,5]=='F')
	provided_dam_data
	length(provided_dam_data[,1])
	index_dams<-provided_dam_data[,1]
	index_dams_1<-match(index_dams,my_internal_data1[,1])
	index_dams_1
	provided_dam_seq<-total_B$sequ[index_dams_1,]
	# dim(provided_dam_seq)
	provided_dam_qtl<-total_B$qtl[index_dams_1,]
	# dim(provided_dam_qtl)
	provided_dam_mrk<-total_B$mrk[index_dams_1,]
	# dim(provided_dam_mrk)
	
	
}

	# dams_toC
	if(Cross_design[2,3]=='rnd'){
	ID_dams_toC<-sample(index_dams,Cross_design[2,2],replace=FALSE)
	index<-match(ID_dams_toC,index_dams)
	} else if(Cross_design[2,3]=='phen'){

		if(Cross_design[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,6]), ]
		ID_dams_toC<-sorted_dam_data[1:Cross_design[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)

		}
		if(Cross_design[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,6]), ]
		ID_dams_toC<-sorted_dam_data[1:Cross_design[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}

	} else if (Cross_design[2,3]=='tbv'){
		if(Cross_design[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,8]), ]
		ID_dams_toC<-sorted_dam_data[1:Cross_design[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
			if(Cross_design[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,8]), ]
		ID_dams_toC<-sorted_dam_data[1:Cross_design[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
		
	} else if (Cross_design[2,3]=='gebv'){
		if(Cross_design[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,9]), ]
		ID_dams_toC<-sorted_dam_data[1:Cross_design[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
			if(Cross_design[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,9]), ]
		ID_dams_toC<-sorted_dam_data[1:Cross_design[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
		
	} else if (Cross_design[2,3]=='tbvc'){
		if(Cross_design[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,10]), ]
		ID_dams_toC<-sorted_dam_data[1:Cross_design[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
			if(Cross_design[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,10]), ]
		ID_dams_toC<-sorted_dam_data[1:Cross_design[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
		
	} else if (Cross_design[2,3]=='gebvc'){
		if(Cross_design[2,4]=='h'){
		sorted_dam_data<-provided_dam_data[order(-provided_dam_data[,11]), ]
		ID_dams_toC<-sorted_dam_data[1:Cross_design[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
			if(Cross_design[2,4]=='l'){
		sorted_dam_data<-provided_dam_data[order(provided_dam_data[,11]), ]
		ID_dams_toC<-sorted_dam_data[1:Cross_design[2,2],]
		index<-match(ID_dams_toC[,1],index_dams)
		}
		
	}

	dams_toC_data<-provided_dam_data[index,]
	dams_toC_seq<-provided_dam_seq[index,]
	# dim(dams_toC_data)
	dams_toC_qtl<-provided_dam_qtl[index,]
	dams_toC_mrk<-provided_dam_mrk[index,]
	# dim(dams_toC_mrk)
	

	start_data<-rbind(sires_toC_data,dams_toC_data)
	total_C<-list()
	total_C$data<-start_data
	total_C$data
	# dim(total_C$data)

	total_C$qtl<-rbind(sires_toC_qtl,dams_toC_qtl)
	# dim(total_C$qtl)
	# class(total_C$qtl)
	total_C$qtl<-as.matrix(total_C$qtl)
	# class(total_C$qtl)

	total_C$mrk<-rbind(sires_toC_mrk,dams_toC_mrk)
	# dim(total_C$mrk)
	total_C$mrk<-as.matrix(total_C$mrk)

	# seq
	total_C$sequ<-rbind(sires_toC_seq,dams_toC_seq)
	# dim(total_C$sequ)
	total_C$sequ<-as.matrix(total_C$sequ)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SELECT FROM BREED A AND B AS PARENTS OF CROSBREDS--FINISH
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end_ge <- Sys.time()
time_gen <- end_ge - start_ge

cat('Generation',gene_counter,'is finished.','Time taken:',time_gen,fill=TRUE)

 
} # end of do loop for generations



# ---------------------START-----------------
# DATA RETURN BACK BY FUNCTION
# ---------------------START-----------------

# #Fill in data of last hp
 cat('Output data preparation ...',fill=TRUE) 

# SUMMARY BREED A START
    for_summary<-data.frame()
	for (i in 1:length(Total_Total_A)){
	ali<-Total_Total_A[[i]][[1]]
	names(ali)<-c(1:11)
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

# without generation zero
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
	if(Selection_pop1[1,2]=='rnd'){
		M_accuracy[i]<-cor(b2[,8],b2[,6])
		} else if (Selection_pop1[1,2]=='phen') {
		M_accuracy[i]<-cor(b2[,8],b2[,6])
		} else if (Selection_pop1[1,2]=='tbv') {
		M_accuracy[i]<-cor(b2[,8],b2[,8])
		} else if (Selection_pop1[1,2]=='gebv') {
		M_accuracy[i]<-cor(b2[,8],b2[,9])
	    } else if (Selection_pop1[1,2]=='tbvc') {
		M_accuracy[i]<-cor(b2[,8],b2[,10])
	    } else if (Selection_pop1[1,2]=='gebvc') {
		M_accuracy[i]<-cor(b2[,8],b2[,11])
	    }
	# Females
	b2<-subset(a,a[,5]=='F')
	if(Selection_pop1[2,2]=='rnd'){
		F_accuracy[i]<-cor(b2[,8],b2[,6])
		} else if (Selection_pop1[2,2]=='phen') {
		F_accuracy[i]<-cor(b2[,8],b2[,6])
		} else if (Selection_pop1[2,2]=='tbv') {
		F_accuracy[i]<-cor(b2[,8],b2[,8])
		} else if (Selection_pop1[2,2]=='gebv') {
		F_accuracy[i]<-cor(b2[,8],b2[,9])
	    } else if (Selection_pop1[1,2]=='tbvc') {
		F_accuracy[i]<-cor(b2[,8],b2[,10])
	    } else if (Selection_pop1[1,2]=='gebvc') {
		F_accuracy[i]<-cor(b2[,8],b2[,11])
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

# if(Selection_pop1[1,2]=='gebv' | Selection_pop1[2,2]=='gebv' | Selection_pop1[1,2]=='gebvc' | Selection_pop1[2,2]=='gebvc'){
# summary_data_A<-data.frame(Generation,Phenotype,TrueBV,GEBV,M_accuracy,F_accuracy,heritability)
# } else {
# summary_data_A<-data.frame(Generation,Phenotype,TrueBV,M_accuracy,F_accuracy,heritability)
# }

if(any(Selection_pop1[,2]=='gebv') | any(Selection_pop1[,2]=='gebvc')){
summary_data_A<-data.frame(Generation,Phenotype,TrueBV,GEBV,M_accuracy,F_accuracy,heritability)
} else {
summary_data_A<-data.frame(Generation,Phenotype,TrueBV,M_accuracy,F_accuracy,heritability)
}

if(Display==TRUE){
print(summary_data_A)
}
#  SUMMARY BREED A FINISH

# SUMMARY BREED B START
    for_summary<-data.frame()
	for (i in 1:length(Total_Total_B)){
	ali<-Total_Total_B[[i]][[1]]
	names(ali)<-c(1:11)
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

# without generation zero
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
	if(Selection_pop2[1,2]=='rnd'){
		M_accuracy[i]<-cor(b2[,8],b2[,6])
		} else if (Selection_pop2[1,2]=='phen') {
		M_accuracy[i]<-cor(b2[,8],b2[,6])
		} else if (Selection_pop2[1,2]=='tbv') {
		M_accuracy[i]<-cor(b2[,8],b2[,8])
		} else if (Selection_pop2[1,2]=='gebv') {
		M_accuracy[i]<-cor(b2[,8],b2[,9])
	    } else if (Selection_pop2[1,2]=='tbvc') {
		M_accuracy[i]<-cor(b2[,8],b2[,10])
	    } else if (Selection_pop2[1,2]=='gebvc') {
		M_accuracy[i]<-cor(b2[,8],b2[,11])
	    }
	# Females
	b2<-subset(a,a[,5]=='F')
	if(Selection_pop2[2,2]=='rnd'){
		F_accuracy[i]<-cor(b2[,8],b2[,6])
		} else if (Selection_pop2[2,2]=='phen') {
		F_accuracy[i]<-cor(b2[,8],b2[,6])
		} else if (Selection_pop2[2,2]=='tbv') {
		F_accuracy[i]<-cor(b2[,8],b2[,8])
		} else if (Selection_pop2[2,2]=='gebv') {
		F_accuracy[i]<-cor(b2[,8],b2[,9])
	    } else if (Selection_pop2[1,2]=='tbvc') {
		F_accuracy[i]<-cor(b2[,8],b2[,10])
	    } else if (Selection_pop2[1,2]=='gebvc') {
		F_accuracy[i]<-cor(b2[,8],b2[,11])
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

if(any(Selection_pop2[,2]=='gebv') | any(Selection_pop2[,2]=='gebvc')){
summary_data_B<-data.frame(Generation,Phenotype,TrueBV,GEBV,M_accuracy,F_accuracy,heritability)
} else {
summary_data_B<-data.frame(Generation,Phenotype,TrueBV,M_accuracy,F_accuracy,heritability)
}

if(Display==TRUE){
print(summary_data_B)
}
#  SUMMARY BREED B FINISH


# # SUMMARY CROSS START
    # for_summary<-data.frame()
	# for (i in 1:length(Total_Total_C)){
	# ali<-Total_Total_C[[i]][[1]]
	# for_summary<-rbind(for_summary,ali)
	# }
	
# mean_pheno<-as.numeric()
# mean_tbv<-as.numeric()
# mean_gebv<-as.numeric()
# var_tbv<-as.numeric()
# var_pheno<-as.numeric()
# h2p<-as.numeric()
# corr_tbv_gebv<-as.numeric()

# # without generation zero
# for (i in 1:ng){
# a<-subset(for_summary,for_summary[,4]==i)
# mean_pheno[i]<-mean(a[,6])
# var_pheno[i]<-var(a[,6])
# mean_tbv[i]<-mean(a[,8])
# var_tbv[i]<-var(a[,8])
# mean_gebv[i]<-mean(a[,9])
# if(Cross_design[1,3]=='gebv' | Cross_design[2,3]=='gebv'){
# corr_tbv_gebv[i]<-cor(a[,8],a[,9])
# } else {
# corr_tbv_gebv<-rep('Not available',ng)
# }
# h2p[i]<-var(a[,8])/var(a[,6])
# }

# Generation<-1:ng
# Phenotype<-mean_pheno
# TrueBV<-mean_tbv
# GEBV<-mean_gebv
# Corelation<-corr_tbv_gebv
# heritability<-h2p

# heritability<-heritability*4
# summary_data_cross<-data.frame(Generation,Phenotype,TrueBV,GEBV,Corelation,heritability)


# if(Display==TRUE){
# print(summary_data_cross)
# }

# # #  SUMMARY CROSS FINISH

# SUMMARY CROSS START
    for_summary<-data.frame()
	for (i in 1:length(Total_Total_C)){
	ali<-Total_Total_C[[i]][[1]]
	for_summary<-rbind(for_summary,ali)
	}
	
mean_pheno<-as.numeric()
mean_tbv<-as.numeric()
mean_gebv<-as.numeric()
var_tbv<-as.numeric()
var_pheno<-as.numeric()
h2p<-as.numeric()


# without generation zero
for (i in 1:ng){
a<-subset(for_summary,for_summary[,4]==i)
mean_pheno[i]<-mean(a[,6])
var_pheno[i]<-var(a[,6])
mean_tbv[i]<-mean(a[,8])
var_tbv[i]<-var(a[,8])
mean_gebv[i]<-mean(a[,9])
h2p[i]<-var(a[,8])/var(a[,6])
}

Generation<-1:ng
Phenotype<-mean_pheno
TrueBV<-mean_tbv
GEBV<-mean_gebv
heritability<-h2p
heritability<-heritability*4

summary_data_cross<-data.frame(Generation,Phenotype,TrueBV,heritability)
	
if(flag_sel_A==TRUE | flag_sel_B==TRUE){
	if(train_type=='crossbred'){
	summary_data_cross<-data.frame(Generation,Phenotype,TrueBV,GEBV,heritability)
	}
}


if(Display==TRUE){
print(summary_data_cross)
}

# #  SUMMARY CROSS FINISH


internal_map_QTL<-outforLD$linkage_map_qtl
internal_map_MRK<-outforLD$linkage_map_mrk
internal_map_qtl_mrk<-outforLD$linkage_map_qtl_mrk
allele_effcts<-outforLD$allele_effcts
trait<-outforLD$trait
genome<-outforLD$genome

# ---------------------FINISH-----------------
# DATA RETURN BACK BY FUNCTION
# ---------------------FINISH-----------------

#WRITE TO OUTPUT POP 1
	if(!missing(saveAt) & !missing(output_pop1)){
	cat('Writing output files of pop 1 ...',fill=TRUE)
	
control_names<-c("data","qtl","marker","seq","freq_qtl","freq_mrk")
test_w<-intersect(names(output_pop1),control_names)

		namepop_1<-paste(saveAt,'_pop1_',sep='')
		# data to output
	if(any(test_w=='data')){
		in_w<-output_pop1$data
		in_w2<-in_w+1
		for(P1 in 1:length(in_w)){
		outFile_data<-paste(namepop_1,'_data_',sep='')
		outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
		dom<-format(Total_Total_A[[in_w2[P1]]][[1]],  justify = "right")
		write.table(dom,file=outFile_data,row.names=FALSE,col.names=TRUE,quote = FALSE) 
		}
	}

	    # qtl to output
	if(any(test_w=='qtl')){
    in_w<-output_pop1$qtl
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_1,'_qtl_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-Total_Total_A[[in_w2[P1]]][[2]]
	writeLines(c('ID  Generaion  Genotypes (paternal allele, maternal allele) ...'),outFile_data)
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)
	
	}
    }	
		
		# mrk to output
	if(any(test_w=='marker')){
    in_w<-output_pop1$marker
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_1,'_mrk_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-Total_Total_A[[in_w2[P1]]][[3]]	
	writeLines(c('ID  Generaion  Genotypes (paternal allele, maternal allele) ...'),outFile_data)
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)
	}
	}
	
    # seq to output
	if(any(test_w=='seq')){
    in_w<-output_pop1$seq
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_1,'_seq_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-Total_Total_A[[in_w2[P1]]][[4]]
	writeLines(c('ID  Generaion  Genotypes (paternal allele, maternal allele) ...'),outFile_data)
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)

	}
	}
	
	    # freq qtl to output
	if(any(test_w=='freq_qtl')){
    in_w<-output_pop1$freq_qtl
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_1,'_freq_qtl_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-format(Total_Total_A[[in_w2[P1]]][[5]],  justify = "right")
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=TRUE,quote = FALSE) 
	}
	}
	
		    # freq qtl to output
	if(any(test_w=='freq_mrk')){
    in_w<-output_pop1$freq_mrk
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_1,'_freq_mrk_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-format(Total_Total_A[[in_w2[P1]]][[6]],  justify = "right")
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=TRUE,quote = FALSE) 
	}
	}
	
	
} #end loop for writing

#WRITE TO OUTPUT POP 2
	if(!missing(saveAt) & !missing(output_pop2)){
	cat('Writing output files of pop 2 ...',fill=TRUE)
	
control_names<-c("data","qtl","marker","seq","freq_qtl","freq_mrk")
test_w<-intersect(names(output_pop2),control_names)
		namepop_2<-paste(saveAt,'_pop2_',sep='')
		# data to output
	if(any(test_w=='data')){
		in_w<-output_pop2$data
		in_w2<-in_w+1
		for(P1 in 1:length(in_w)){
		outFile_data<-paste(namepop_2,'_data_',sep='')
		outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
		dom<-format(Total_Total_B[[in_w2[P1]]][[1]],  justify = "right")
		write.table(dom,file=outFile_data,row.names=FALSE,col.names=TRUE,quote = FALSE) 
		}
	}

	    # qtl to output
	if(any(test_w=='qtl')){
    in_w<-output_pop2$qtl
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_2,'_qtl_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-Total_Total_B[[in_w2[P1]]][[2]]
	writeLines(c('ID  Generaion  Genotypes (paternal allele, maternal allele) ...'),outFile_data)
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)
	}
    }	
		
		# mrk to output
	if(any(test_w=='marker')){
    in_w<-output_pop2$marker
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_2,'_mrk_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-Total_Total_B[[in_w2[P1]]][[3]]
	writeLines(c('ID  Generaion  Genotypes (paternal allele, maternal allele) ...'),outFile_data)
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)
	}
	}
	
    # seq to output
	if(any(test_w=='seq')){
    in_w<-output_pop2$seq
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_2,'_seq_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-Total_Total_B[[in_w2[P1]]][[4]]
	writeLines(c('ID  Generaion  Genotypes (paternal allele, maternal allele) ...'),outFile_data)
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)
	}
	}
	
	    # freq qtl to output
	if(any(test_w=='freq_qtl')){
    in_w<-output_pop2$freq_qtl
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_2,'_freq_qtl_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-format(Total_Total_B[[in_w2[P1]]][[5]],  justify = "right")
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=TRUE,quote = FALSE) 
	}
	}
	
		# freq qtl to output
	if(any(test_w=='freq_mrk')){
    in_w<-output_pop2$freq_mrk
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_2,'_freq_mrk_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-format(Total_Total_B[[in_w2[P1]]][[6]],  justify = "right")
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=TRUE,quote = FALSE) 
	}
	}
	
	
} #end loop for writing

#WRITE TO OUTPUT CROSS
	if(!missing(saveAt) & !missing(output_cross)){
	cat('Writing output files of crossbreds ...',fill=TRUE)
	
control_names<-c("data","qtl","marker","seq","freq_qtl","freq_mrk")
test_w<-intersect(names(output_cross),control_names)
		namepop_cross<-paste(saveAt,'_cross_',sep='')
		# data to output
	if(any(test_w=='data')){
		in_w<-output_cross$data
		in_w2<-in_w+1
		for(P1 in 1:length(in_w)){
		outFile_data<-paste(namepop_cross,'_data_',sep='')
		outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
		dom<-format(Total_Total_C[[in_w2[P1]]][[1]],  justify = "right")
		write.table(dom,file=outFile_data,row.names=FALSE,col.names=TRUE,quote = FALSE) 
		}
	}

	    # qtl to output
	if(any(test_w=='qtl')){
    in_w<-output_cross$qtl
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_cross,'_qtl_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-Total_Total_C[[in_w2[P1]]][[2]]
	writeLines(c('ID  Generaion  Genotypes (paternal allele, maternal allele) ...'),outFile_data)
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)
	}
    }	
		
		# mrk to output
	if(any(test_w=='marker')){
    in_w<-output_cross$marker
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_cross,'_mrk_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-Total_Total_C[[in_w2[P1]]][[3]]
	writeLines(c('ID  Generaion  Genotypes (paternal allele, maternal allele) ...'),outFile_data)
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)
	}
	}
	
    # seq to output
	if(any(test_w=='seq')){
    in_w<-output_cross$seq
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_cross,'_seq_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-Total_Total_C[[in_w2[P1]]][[4]]
	writeLines(c('ID  Generaion  Genotypes (paternal allele, maternal allele) ...'),outFile_data)
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=FALSE,quote = FALSE,append=TRUE)
	}
	}
	
	    # freq qtl to output
	if(any(test_w=='freq_qtl')){
    in_w<-output_cross$freq_qtl
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_cross,'_freq_qtl_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-format(Total_Total_C[[in_w2[P1]]][[5]],  justify = "right")
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=TRUE,quote = FALSE) 
	}
	}
	
		# freq qtl to output
	if(any(test_w=='freq_mrk')){
    in_w<-output_cross$freq_mrk
	in_w2<-in_w+1
	for(P1 in 1:length(in_w)){
	outFile_data<-paste(namepop_cross,'_freq_mrk_',sep='')
	outFile_data<-paste(outFile_data,paste(in_w[P1],'.txt',sep=''),sep='')
	dom<-format(Total_Total_C[[in_w2[P1]]][[6]],  justify = "right")
	write.table(dom,file=outFile_data,row.names=FALSE,col.names=TRUE,quote = FALSE) 
	}
	}
	
	
} #end loop for writing

userin<-list()
userin[[1]]<-founder_cross[1,]
userin[[2]]<-founder_cross[2,]
userin[[3]]<-ng
userin[[4]]<-litter_size
tob<-Cross_design[,-1]
tob[,2]<-'rnd'
if(flag_sel_A==TRUE | flag_sel_B==TRUE){
	if(train_type=='crossbred'){
	tob[,2]<-'gebv'
	}
}
userin[[5]]<-tob


for (g_index in 1:length(Total_Total_A)){
names(Total_Total_A[[g_index]])<-c('data','qtl','mrk','sequ','freqQTL','freqMRK')
}

for (g_index in 1:length(Total_Total_B)){
names(Total_Total_B[[g_index]])<-c('data','qtl','mrk','sequ','freqQTL','freqMRK')
}	  

for (g_index in 1:length(Total_Total_C)){
names(Total_Total_C[[g_index]])<-c('data','qtl','mrk','sequ','freqQTL','freqMRK')
}	


cat('xbreed is done!',fill=TRUE)	

# RETURN SECTION	
cross<-list(output=Total_Total_C,summary_data=summary_data_cross,linkage_map_qtl=internal_map_QTL,linkage_map_mrk=internal_map_MRK,linkage_map_qtl_mrk=internal_map_qtl_mrk,allele_effcts=allele_effcts,trait=trait,genome=genome,user_input=userin)

Final<-list(pop1=Total_Total_A,pop2=Total_Total_B,cross=cross,summary_data_pop1=summary_data_A,summary_data_pop2=summary_data_B,summary_data_cross=summary_data_cross,linkage_map_qtl=internal_map_QTL,linkage_map_mrk=internal_map_MRK,linkage_map_qtl_mrk=internal_map_qtl_mrk,allele_effcts=allele_effcts,trait=trait,genome=genome)
return(Final)

} # end of function


