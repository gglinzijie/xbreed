SUBROUTINE sh(in_r_mu,in_c_mu,in_pois,hp_loci,in_pois_f,hp_loci_f,r_rec,c_rec,nchr,szch,rec_m,rf,cf,rec_f)


!Input official
INTEGER,INTENT(IN)::in_r_mu,in_c_mu
INTEGER,INTENT(IN), dimension(1:in_r_mu,1:in_c_mu):: in_pois
INTEGER,INTENT(IN), dimension(1:in_r_mu,1:in_c_mu):: in_pois_f

! For crossover
!MALES
INTEGER,INTENT(IN)::r_rec,c_rec,nchr
INTEGER,INTENT(IN), DIMENSION(1:r_rec,1:c_rec):: rec_m
INTEGER,INTENT(IN), DIMENSION(1:(nchr+1))::szch

!FEMALES
INTEGER,INTENT(IN)::rf,cf

INTEGER,INTENT(IN), DIMENSION(1:rf,1:cf):: rec_f


!Output official
INTEGER,DIMENSION(1:(in_r_mu*2),1:(in_c_mu/2)),INTENT(INOUT) :: hp_loci
INTEGER,DIMENSION(1:(in_r_mu*2),1:(in_c_mu/2)),INTENT(INOUT) :: hp_loci_f

!Local variables
INTEGER::i,counter,j
REAL :: start, finish
INTEGER,ALLOCATABLE,DIMENSION(:,:):: Males_M
INTEGER,ALLOCATABLE,DIMENSION(:,:):: Males_M_Hap
INTEGER,ALLOCATABLE,DIMENSION(:,:):: Females_M
INTEGER,ALLOCATABLE,DIMENSION(:,:):: Females_M_Hap
INTEGER,ALLOCATABLE,DIMENSION(:):: seq1,seq2
INTEGER,ALLOCATABLE,DIMENSION(:):: S_1,S_2

!CROSSOVER
INTEGER::a1,b1,T1,T2,k1,shom,nrec,ad1,ad2,tcount,ali
INTEGER,ALLOCATABLE,DIMENSION(:,:):: Chrom,twoin,twoout
INTEGER,ALLOCATABLE,DIMENSION(:,:):: temp_config_chr	 
INTEGER,ALLOCATABLE,DIMENSION(:):: my_loc,posrec



ALLOCATE(Males_M(1:size(in_pois,1),1:size(in_pois,2)))
ALLOCATE(Females_M(1:size(in_pois_f,1),1:size(in_pois_f,2)))


ALLOCATE(Males_M_Hap(1:(in_r_mu*2),1:(in_c_mu/2)))
ALLOCATE(Females_M_Hap(1:(in_r_mu*2),1:(in_c_mu/2)))

!TEST----------------------
!        DO i=1,10
!        DO j=1,4
!        WRITE(*,'(I4,X)',advance='NO')  rec_m(i,j)
!        END DO
!        WRITE(*, *) '' 
!        END DO
!TEST--------------------FINISH

!---------------------------------------------
!ONE TO TWO STRING---MALES----START----------------   
!---------------------------------------------
ALLOCATE(seq1(1:size(Males_M,2)/2))
ALLOCATE(seq2(1:size(Males_M,2)/2))
ALLOCATE(S_1(1:size(Males_M,1)))
ALLOCATE(S_2(1:size(Males_M,1)))


Males_M=in_pois

counter=1
do i=1,size(Males_M,2),2
seq1(counter)=i
counter=1+counter
end do
counter=1
do i=2,size(Males_M,2),2
seq2(counter)=i
counter=1+counter
end do


counter=1
do i=1,size(Males_M,1)*2,2
S_1(counter)=i
counter=1+counter
end do

counter=1
do i=2,size(Males_M,1)*2,2
S_2(counter)=i
counter=1+counter
end do

Do i=1,size(Males_M,1)
Males_M_Hap(S_1(i),:)=Males_M(i,seq1)
Males_M_Hap(S_2(i),:)=Males_M(i,seq2)
END DO

hp_loci=Males_M_Hap  



DEALLOCATE(S_1)
DEALLOCATE(S_2)
DEALLOCATE(seq1)
DEALLOCATE(seq2)

!---------------------------------------------
!ONE TO TWO STRING-------FINISH---------------   
!---------------------------------------------



!---------------------------------------------
!CROSSOVER MALES----START----------------   
!---------------------------------------------
a1=1
b1=2
DO i=1,nchr
ALLOCATE(Chrom(1:size(hp_loci,1),(szch(a1)+1):szch(b1)))

Chrom=hp_loci(:,(szch(a1)+1):szch(b1))
T1=size(pack(rec_m(:,2),rec_m(:,2)==i))
ALLOCATE(my_loc(1:T1))

tcount=1
Do k1=1,size(rec_m,1)
IF(rec_m(k1,2)==i) THEN
my_loc(tcount)=k1
tcount=tcount+1
END IF
END DO

ALLOCATE(temp_config_chr(1:T1,1:size(rec_m,2)))
temp_config_chr=rec_m(my_loc,:)


T2=size(rec_m,2)-3
shom=1
DO j=1,size(Chrom,1),2
ALLOCATE(twoin(1:2,1:size(Chrom,2)))
ALLOCATE(posrec(1:T2))
ALLOCATE(twoout(1:2,1:size(Chrom,2)))
twoin(1:2,:)=Chrom(j:(j+1),:)
nrec=temp_config_chr(shom,3)
posrec(1:T2)=temp_config_chr(shom,4:size(rec_m,2))
ad1=size(twoin,2)
ad2=T2
CALL cross_over(twoin,nrec,posrec,twoout,ad1,ad2)
Chrom(j:(j+1),:)=twoout(1:2,:)
shom=1+shom
DEALLOCATE(twoin)
DEALLOCATE(posrec)
DEALLOCATE(twoout)

END DO

hp_loci(:,(szch(a1)+1):szch(b1))=Chrom
a1=a1+1
b1=b1+1
DEALLOCATE(Chrom)
DEALLOCATE(my_loc)
DEALLOCATE(temp_config_chr)

END DO
!---------------------------------------------
!CROSSOVER MALES-------FINISH---------------   
!---------------------------------------------


!---------------------------------------------
!ONE TO TWO STRING----Females---START----------------   
!---------------------------------------------
ALLOCATE(seq1(1:size(Females_M,2)/2))
ALLOCATE(seq2(1:size(Females_M,2)/2))
ALLOCATE(S_1(1:size(Females_M,1)))
ALLOCATE(S_2(1:size(Females_M,1)))


Females_M=in_pois_f

counter=1
do i=1,size(Females_M,2),2
seq1(counter)=i
counter=1+counter
end do
counter=1
do i=2,size(Females_M,2),2
seq2(counter)=i
counter=1+counter
end do


counter=1
do i=1,size(Females_M,1)*2,2
S_1(counter)=i
counter=1+counter
end do

counter=1
do i=2,size(Females_M,1)*2,2
S_2(counter)=i
counter=1+counter
end do

Do i=1,size(Females_M,1)
Females_M_Hap(S_1(i),:)=Females_M(i,seq1)
Females_M_Hap(S_2(i),:)=Females_M(i,seq2)
END DO

hp_loci_f=Females_M_Hap    


DEALLOCATE(S_1)
DEALLOCATE(S_2)
DEALLOCATE(seq1)
DEALLOCATE(seq2)

!---------------------------------------------
!ONE TO TWO STRING-------FINISH---------------   
!---------------------------------------------

!---------------------------------------------
!CROSSOVER FEMALES----START----------------   
!---------------------------------------------
a1=1
b1=2
DO i=1,nchr
ALLOCATE(Chrom(1:size(hp_loci_f,1),(szch(a1)+1):szch(b1)))

Chrom=hp_loci_f(:,(szch(a1)+1):szch(b1))
T1=size(pack(rec_f(:,2),rec_f(:,2)==i))
ALLOCATE(my_loc(1:T1))

tcount=1
Do k1=1,size(rec_f,1)
IF(rec_f(k1,2)==i) THEN
my_loc(tcount)=k1
tcount=tcount+1
END IF
END DO

ALLOCATE(temp_config_chr(1:T1,1:size(rec_f,2)))
temp_config_chr=rec_f(my_loc,:)


T2=size(rec_f,2)-3
shom=1
DO j=1,size(Chrom,1),2
ALLOCATE(twoin(1:2,1:size(Chrom,2)))
ALLOCATE(posrec(1:T2))
ALLOCATE(twoout(1:2,1:size(Chrom,2)))
twoin(1:2,:)=Chrom(j:(j+1),:)
nrec=temp_config_chr(shom,3)
posrec(1:T2)=temp_config_chr(shom,4:size(rec_f,2))
ad1=size(twoin,2)
ad2=T2
CALL cross_over(twoin,nrec,posrec,twoout,ad1,ad2)
Chrom(j:(j+1),:)=twoout(1:2,:)
shom=1+shom
DEALLOCATE(twoin)
DEALLOCATE(posrec)
DEALLOCATE(twoout)

END DO

hp_loci_f(:,(szch(a1)+1):szch(b1))=Chrom
a1=a1+1
b1=b1+1
DEALLOCATE(Chrom)
DEALLOCATE(my_loc)
DEALLOCATE(temp_config_chr)

END DO
!---------------------------------------------
!CROSSOVER FEMALES-------FINISH---------------   
!---------------------------------------------





!DEALLOCATION SECTION START-----------
DEALLOCATE(Males_M)
DEALLOCATE(Males_M_Hap)
DEALLOCATE(Females_M)
DEALLOCATE(Females_M_Hap)



END SUBROUTINE sh


!  1-SUBROUTINE FOR CROSSOVER

SUBROUTINE cross_over(twoin_rec,nrec_rec,posrec_rec,twoout_rec,size1,size2)

INTEGER, INTENT(IN) ::size1,size2,nrec_rec
INTEGER, INTENT(IN) :: twoin_rec(1:2,1:size1)
INTEGER, INTENT(IN) :: posrec_rec(1:size2)  
!OUT
INTEGER, INTENT(INOUT) :: twoout_rec(1:2,1:size1)
!LOCAL
INTEGER,ALLOCATABLE,DIMENSION(:):: s1,s2,s1new,s2new
INTEGER ::L1
!MOSAHEBAT
ALLOCATE(s1(1:size(twoin_rec,2)))
ALLOCATE(s2(1:size(twoin_rec,2)))
ALLOCATE(s1new(1:size(twoin_rec,2)))
ALLOCATE(s2new(1:size(twoin_rec,2)))

s1=twoin_rec(1,:)
s2=twoin_rec(2,:)
s1new=s1
s2new=s2

Do L1=1,nrec_rec
s1new(posrec_rec(L1):size(s1))=s2(posrec_rec(L1):size(s1))
s2new(posrec_rec(L1):size(s1))=s1(posrec_rec(L1):size(s1))
s1=s1new
s2=s2new
END DO

twoout_rec(1,:)=s1
twoout_rec(2,:)=s2

DEALLOCATE(s1)
DEALLOCATE(s2)
DEALLOCATE(s1new)
DEALLOCATE(s2new)

END SUBROUTINE cross_over