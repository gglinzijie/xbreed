
SUBROUTINE hp(hpsize,nmarker,nqtl,loc_c,ng,in_r_mu,in_c_mu,in_pois,hp_loci,method_laf,arglaf_F,rcs,cross)

!Input official
INTEGER,INTENT(IN)::hpsize,nmarker,nqtl,loc_c,ng,method_laf
INTEGER,INTENT(IN)::in_r_mu,in_c_mu,rcs
DOUBLE PRECISION ,INTENT(IN) :: arglaf_F
!Output official

!Local variables
INTEGER::nloci,i,allele,in_nloci,j,k,counter,popsze,nalle,nsire,ad1,ad2,gen_counter,sh1,sh2,sh3,sh4,nchr
REAL :: start, finish
!INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT) :: hp_loci

INTEGER,DIMENSION(1:hpsize,1:loc_c),INTENT(INOUT) :: hp_loci
INTEGER,INTENT(IN), dimension(1:in_r_mu,1:in_c_mu):: in_pois
INTEGER,INTENT(IN), dimension(1:rcs,1:2):: cross
INTEGER,ALLOCATABLE,DIMENSION(:,:):: Males_M,Females_M,Start_pop
INTEGER,ALLOCATABLE,DIMENSION(:,:):: Tempo_Males_M,Tempo_Females_M
INTEGER,ALLOCATABLE,DIMENSION(:,:):: Males_M_Hap,Females_M_Hap

INTEGER,ALLOCATABLE,DIMENSION(:,:):: Haplo_prod_M
INTEGER,ALLOCATABLE,DIMENSION(:,:):: Haplo_prod_F
INTEGER,ALLOCATABLE,DIMENSION(:):: seq1,seq2,TemC2,TemC1
INTEGER,ALLOCATABLE,DIMENSION(:):: seq1_mtg,seq2_mtg
INTEGER,ALLOCATABLE,DIMENSION(:):: S_1,S_2,S_3,S_4
INTEGER,DIMENSION(:,:), ALLOCATABLE:: Offspring_2_str
INTEGER,DIMENSION(:,:), ALLOCATABLE:: Offspring_1_str,Offscross_1_str,in_cross,out_cross
INTEGER,ALLOCATABLE,DIMENSION(:):: M_shuf_array,F_shuf_array
INTEGER,DIMENSION(:,:), ALLOCATABLE:: in_pois_internal
INTEGER,DIMENSION(:,:), ALLOCATABLE:: Tempo_pois

!integer, dimension(1:nn1,1:3), intent(inout):: totalsub

REAL,ALLOCATABLE,DIMENSION(:) :: prob
REAL::r,RandomReal
CHARACTER (50)::name_1,name_2

!Declare new type anim
TYPE First
INTEGER::Idi
INTEGER::Generation
INTEGER::sex			
INTEGER,DIMENSION(:), ALLOCATABLE:: anim_loci  
END type First

TYPE (First), ALLOCATABLE,DIMENSION(:) :: Anim
!local variables for test, they will be omited

!Prelimary calcs
nloci=nmarker+nqtl
in_nloci=nloci*2
name_1='HP_Loci.txt'
name_2='HP_Loci_Prob.txt'
nsire=hpsize/2


!TEST start	
CALL init_random_seed()

!Allocate section
ALLOCATE(prob(1:nloci))
ALLOCATE(Start_pop(1:hpsize,1:loc_c))



!Fill Anim



!Create Loci
SELECT CASE (method_laf) 
CASE (1) !#two options 1:unif 2:equal
DO i=1,nloci
call RANDOM_NUMBER(r)
prob(i)=r 
!	print*,'prob(i)',prob(i)
END DO

counter=0	
DO i=1,in_nloci,2
counter=counter+1
DO j=1,hpsize
CALL RANDOM_NUMBER(r)
IF (r .LE. prob(counter)) THEN
Start_pop(j,i)=1
ELSE 
Start_pop(j,i)=2
END IF	   
IF (1-r .GT. prob(counter)) THEN
Start_pop(j,i+1)=1
ELSE 
Start_pop(j,i+1)=2		   
END IF
END DO
END DO
CASE (2) !#two options 1:unif 2:equal
DO i=1,nloci
prob(i)=arglaf_F 
!	print*,'prob(i)',prob(i)
END DO

counter=0	
DO i=1,in_nloci,2
counter=counter+1
DO j=1,hpsize
CALL RANDOM_NUMBER(r)
IF (r .LE. prob(counter)) THEN
Start_pop(j,i)=1		   
ELSE 
Start_pop(j,i)=2
END IF	   
IF (1-r .LE. prob(counter)) THEN
Start_pop(j,i+1)=1		   
ELSE 
Start_pop(j,i+1)=2		   
END IF
END DO
END DO

END  SELECT	

!/////////////////////////////////


ALLOCATE(S_1(1:nloci))
ALLOCATE(S_2(1:nloci))
ALLOCATE(S_3(1:hpsize/2))
ALLOCATE(S_4(1:hpsize/2))

counter=1
do i=1,loc_c,2
S_1(counter)=i
counter=1+counter
end do

counter=1
do i=2,loc_c,2
S_2(counter)=i
counter=1+counter
end do

counter=1
do i=1,hpsize,2
S_3(counter)=i
counter=1+counter
end do

counter=1
do i=2,hpsize,2
S_4(counter)=i
counter=1+counter
end do


!TEST FOR MUTATION--- START

!          DO i=1,3
!            DO j=1,2   
!           WRITE(*,'(I4,X)',advance='NO') cross(i,j)
!            END DO
!           write(*, *) '' 
!          END DO

!TEST FOR MUTATION--- FINISH

!Start generations
Do gen_counter=1, ng
CALL init_random_seed()


ALLOCATE(Males_M(1:hpsize/2,1:loc_c))
ALLOCATE(Females_M(1:hpsize/2,1:loc_c))
ALLOCATE(Tempo_Males_M(1:hpsize/2,1:loc_c))
ALLOCATE(Tempo_Females_M(1:hpsize/2,1:loc_c))

!	        ALLOCATE(Males_M_Hap(1:hpsize,1:loc_c/2))
!		    ALLOCATE(Females_M_Hap(1:hpsize,1:loc_c/2))
!			   ALLOCATE(seq1(1:nloci))
!			   ALLOCATE(seq2(1:nloci))
ALLOCATE(Tempo_pois(1:hpsize,1:loc_c))


!---------------------------------------------
!MUTATION  START-------------------------   
!---------------------------------------------

!     call cpu_time(start)
! put code to test here

!POINT	  
!when there is no mutation matrix(0,3*3) then this part should escape

Tempo_pois=Start_pop			  




!TEST FOR MUTATION--- START
!           IF (gen_counter==5) THEN
!           DO i=1,1
!           DO j=1,size(Tempo_pois,2)       
!           WRITE(61,'(I1,X)',advance='NO') Tempo_pois(25,j)
!          END DO
!           write(61, *) '' 
!           END DO
!          END IF
!TEST FOR MUTATION--- FINISH

counter=0
do i=1,size(in_pois,1)
if(in_pois(i,1) .eq. gen_counter) THEN
counter=counter+1
!                  print*,'shomaresh',counter
END IF
END DO

IF     (counter>0) THEN          
ALLOCATE(in_pois_internal(1:counter,1:3))
in_pois_internal(1:counter,1)=pack(in_pois(:,1),in_pois(:,1)==gen_counter)
in_pois_internal(1:counter,2)=pack(in_pois(:,2),in_pois(:,1)==gen_counter)
in_pois_internal(1:counter,3)=pack(in_pois(:,3),in_pois(:,1)==gen_counter)

Do i=1,size(in_pois_internal,1)
DO j=1,in_pois_internal(i,3)
CALL RANDOM_NUMBER(R)
ad1=1+int(in_nloci*R) !Random allele pos
!          print*,'Random allele pos',ad1
IF(Tempo_pois(in_pois_internal(i,2),ad1) .EQ. 1) THEN 
Tempo_pois(in_pois_internal(i,2),ad1)=2
ELSE 
IF (Tempo_pois(in_pois_internal(i,2),ad1) .EQ. 2)THEN 
Tempo_pois(in_pois_internal(i,2),ad1)=1
END IF
END IF
END DO
END DO


DEALLOCATE(in_pois_internal) 

END IF


!      call cpu_time(finish)
!   print '("Time mutation= ",f6.3," seconds.")',finish-start

!---------------------------------------------
!MUTATION  FINISH---------------    
!---------------------------------------------
!          call cpu_time(start)

Males_M=Tempo_pois(S_3,:)
Females_M=Tempo_pois(S_4,:)
!	   call cpu_time(finish)
!         print '("Time Filling= ",f6.3," seconds.")',finish-start	

!---------------------------------------------
!SUHFFLING--------START-----------------------   
!---------------------------------------------
!Shuffling males
!          call cpu_time(start)

ALLOCATE(M_shuf_array(1:size(Males_M,1)))
Do i=1,size(Males_M,1)
M_shuf_array(i)=i
END DO

ad1=size(Males_M,1)
CALL shuf(M_shuf_array,1,ad1)
Tempo_Males_M=Males_M
Males_M=Tempo_Males_M(M_shuf_array,:)
DEALLOCATE(M_shuf_array)


!Shuffling Females
ALLOCATE(F_shuf_array(1:size(Females_M,1)))
Do i=1,size(Females_M,1)
F_shuf_array(i)=i
END DO

ad2=size(Females_M,1)
CALL shuf(F_shuf_array,1,ad2)
Tempo_Females_M=Females_M
Females_M=Tempo_Females_M(F_shuf_array,:)
DEALLOCATE(F_shuf_array)
!
!	   call cpu_time(finish)
!          print '("Time shuffling= ",f6.3," seconds.")',finish-start			   

!---------------------------------------------
!SUHFFLING--------FINISH--------------------   
!---------------------------------------------
!********************************************
!NEW strategy start**************************		   !********************************************
!    call cpu_time(start)  	   
ALLOCATE(Offspring_1_str(1:hpsize,1:size(Males_M,2)))
ALLOCATE(Offscross_1_str(1:hpsize,1:size(Males_M,2)))



counter=0
DO  j=1,2 ! because two offspring per parent
DO i=1,size(Males_M,1)
CALL RANDOM_NUMBER(R)
counter=1+counter
IF (R .GE. 0.5) THEN
Offspring_1_str(counter,S_1)=Males_M(i,S_1)
Offspring_1_str(counter,S_2)=Females_M(i,S_1)
ELSE
Offspring_1_str(counter,S_1)=Males_M(i,S_2)
Offspring_1_str(counter,S_2)=Females_M(i,S_2)
END IF
END DO 
END DO


!	      call cpu_time(finish)
!         print '("Time NEW STR= ",f6.3," seconds.")',finish-start	

!********************************************
!NEW strategy FINISH**************************		   !********************************************	

!********************************************
!Cross over START****************************		   !********************************************

!       call cpu_time(start)  
nchr=size(cross,1)
Do i=1,nchr
ALLOCATE(in_cross(1:hpsize,cross(i,1):cross(i,2)))
ALLOCATE(out_cross(1:hpsize,cross(i,1):cross(i,2)))

in_cross=Offspring_1_str(:,cross(i,1):cross(i,2))
nalle=size(in_cross,2)
call crsovr(in_cross,out_cross,hpsize,nalle)
Offscross_1_str(:,cross(i,1):cross(i,2))=out_cross  

DEALLOCATE(in_cross)
DEALLOCATE(out_cross) 
End Do
!	      call cpu_time(finish)
!      print '("Time cross= ",f6.3," seconds.")',finish-start	

!********************************************
!Cross over FINISH***************************		   !********************************************

Start_pop=Offscross_1_str


!FILLING hp_LOCI TO GO BACK TO R
IF (gen_counter==ng) THEN
hp_loci=Offscross_1_str
END IF



!DEALLOCATION SECTION START-----------
DEALLOCATE(Males_M)
DEALLOCATE(Females_M)
DEALLOCATE(Tempo_Males_M)
DEALLOCATE(Tempo_Females_M)
DEALLOCATE(Offspring_1_str)
DEALLOCATE(Offscross_1_str)
DEALLOCATE(Tempo_pois)

!DEALLOCATION SECTION FINISH------------------			   
IF (MOD(gen_counter,51)== 0) THEN
! print*, 'Generation No.:',gen_counter,' is finished.'
call intpr ('Generation No.:',-1,gen_counter,1)
END IF



END DO   !Do main



END SUBROUTINE hp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!EXTERNAL SUBROUTINES!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!My Subroutines
! 1-SUBROUTINE FOR RANOMIZATION
SUBROUTINE init_random_seed()
INTEGER :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)
END SUBROUTINE init_random_seed

!  2-SUBROUTINE FOR SHUFFLING

SUBROUTINE shuf(SIREIDSH,n1,n2)
INTEGER :: ISH, RANDPOS, TEMP,n1,n2                 
INTEGER, INTENT(INOUT) :: SIREIDSH(n1:n2)
REAL :: R
DO ISH = SIZE(SIREIDSH), 2, -1
CALL RANDOM_NUMBER(R)
RANDPOS = INT(R * ISH) + 1
TEMP = SIREIDSH(RANDPOS)
SIREIDSH(RANDPOS) = SIREIDSH(ISH)
SIREIDSH(ISH) = TEMP
END DO
END SUBROUTINE shuf


!  2-SUBROUTINE FOR Gamete_union
SUBROUTINE gam_uni(pat_hap,mat_hap,off_hap,size1,size2)
INTEGER :: ISH,ls
INTEGER, INTENT(IN) ::size1,size2
INTEGER, INTENT(IN) :: pat_hap(1:2,1:size1),mat_hap(1:2,1:size1)
INTEGER, INTENT(INOUT) :: off_hap(1:4,1:size2)

REAL :: R
ls=2 !2 to keep pop size constant
CALL RANDOM_NUMBER(R)
IF (R .GE. 0.5) THEN
off_hap(1,:)=pat_hap(1,:) 
ELSE
off_hap(1,:)=pat_hap(2,:) 
END IF

CALL RANDOM_NUMBER(R)
IF (R .GE. 0.5) THEN
off_hap(2,:)=mat_hap(1,:) 
ELSE
off_hap(2,:)=mat_hap(2,:) 
END IF

CALL RANDOM_NUMBER(R)
IF (R .GE. 0.5) THEN
off_hap(3,:)=pat_hap(1,:) 
ELSE
off_hap(3,:)=pat_hap(2,:) 
END IF

CALL RANDOM_NUMBER(R)
IF (R .GE. 0.5) THEN
off_hap(4,:)=mat_hap(1,:) 
ELSE
off_hap(4,:)=mat_hap(2,:) 
END IF
END SUBROUTINE gam_uni


!  2-SUBROUTINE FOR Gamete_union
SUBROUTINE crsovr(recieve,output,size1,nallele)
INTEGER, INTENT(IN) ::size1,nallele !number of allele per chro
INTEGER, INTENT(IN) :: recieve(1:size1,1:nallele)
INTEGER, INTENT(INOUT) :: output(1:size1,1:nallele)
!Local 
REAL :: R
INTEGER :: sh1,nmrkpc,counter,i,j
INTEGER,ALLOCATABLE,DIMENSION(:):: seq1,seq2,TemC2,TemC1

nmrkpc=nallele/2

Do i=1,size1
CALL RANDOM_NUMBER(R)
sh1=1 + FLOOR((nmrkpc-1)*R) !Position of crossover
ALLOCATE(TemC1(1:(nallele-(sh1*2))))
ALLOCATE(TemC2(1:(nallele-(sh1*2))))
TemC1=recieve(i,sh1*2+1:nallele)
ALLOCATE(seq1(1:size(TemC1)/2))
ALLOCATE(seq2(1:size(TemC1)/2))

counter=1
Do j=1,size(TemC1),2
seq1(counter)=j
counter=1+counter
END DO
seq2=seq1+1

TemC2(seq1)=TemC1(seq2)
TemC2(seq2)=TemC1(seq1)
output(i,1:sh1*2)=recieve(i,1:sh1*2)	
output(i,sh1*2+1:nallele)=TemC2

!DEALLOCATE CROSS  
DEALLOCATE(seq1)
DEALLOCATE(seq2)
DEALLOCATE(TemC1)
DEALLOCATE(TemC2)
END DO
END SUBROUTINE crsovr


subroutine intpr(label, nchar, da, ndata)
integer:: nchar, ndata
character*(*) label
integer:: da
integer nc
nc = nchar
if(nc .lt. 0) nc = len(label)
end subroutine intpr











