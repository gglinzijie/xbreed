SUBROUTINE ld(r_size,c_size,loci_mat,npair,method,ld_data)

!Input official
INTEGER,INTENT(IN)::r_size,c_size,npair,method
INTEGER,INTENT(IN), DIMENSION(1:r_size,1:c_size):: loci_mat
!Output official
DOUBLE PRECISION ,DIMENSION(1:npair,1:12),INTENT(INOUT) :: ld_data
!Local variables
INTEGER::i,j,nloci,counter,a,b,c,t,sth,paircounter
REAL :: outfreq,mean_LD
REAL :: p1,q1,p2,q2,hap_11,hap_12,hap_21,hap_22,r,D,r2
INTEGER,ALLOCATABLE,DIMENSION(:)::seq1,seq2
REAL,ALLOCATABLE,DIMENSION(:)::freq1_main,freq2_main
INTEGER,ALLOCATABLE,DIMENSION(:,:):: pair_1,pair_2
REAL,DIMENSION(1:4):: jj



SELECT CASE (method) 
CASE (1) !two options 1:pairwise 2:adjacent   
nloci=size(loci_mat,2)/2
!	print*,'nloci',nloci

ALLOCATE(freq1_main(1:nloci))
ALLOCATE(freq2_main(1:nloci))
ALLOCATE(seq1(1:nloci))
ALLOCATE(seq2(1:nloci))
ALLOCATE(pair_1(1:size(loci_mat,1),1:2))
ALLOCATE(pair_2(1:size(loci_mat,1),1:2))

!calc frequency

!			counter=1
!   			do i=1,nloci*2,2
!			seq1(counter)=i
!			counter=1+counter
!  			end do

!			counter=1
!		   	do i=2,nloci*2,2
!			seq2(counter)=i
!			counter=1+counter
!		   	end do

!   	   do i=1,nloci
!	     call calc_freq(loci_mat(:,seq1(i):seq2(i)),size(loci_mat,1),outfreq)
!        freq1_main(i)=outfreq
!	     END DO

!	        do i=1,10
!	  print*,'freq1(i)',freq1_main(i)
!            end do


a=1
b=2
sth=1 
paircounter=1
do i=1,nloci-1
pair_1=loci_mat(:,a:b)
c=a+2
t=b+2
call calc_freq(pair_1,size(loci_mat,1),outfreq)
p1=outfreq
q1=1-p1
do j=1,nloci-sth
pair_2=loci_mat(:,c:t)
call calc_freq(pair_2,size(loci_mat,1),outfreq)
p2=outfreq
q2=1-p2
call  calc_hap_freq(pair_1,pair_2,size(loci_mat,1),jj)
hap_11=jj(1)
hap_12=jj(2)
hap_21=jj(3)
hap_22=jj(4)
D=REAL(hap_11*hap_22)-REAL(hap_12*hap_21)
r=REAL(D)/REAL(sqrt(p1*q1*p2*q2))
r2=REAL(D**2)/REAL(p1*q1*p2*q2)
ld_data(paircounter,1)=	paircounter		 

ld_data(paircounter,2)=	p1		 
ld_data(paircounter,3)=	q1	
ld_data(paircounter,4)=	p2		 
ld_data(paircounter,5)=	q2	

ld_data(paircounter,6)=	hap_11	
ld_data(paircounter,7)=	hap_12		 
ld_data(paircounter,8)=	hap_21		 
ld_data(paircounter,9)=	hap_22	

ld_data(paircounter,10)=D		 
ld_data(paircounter,11)=r		 
ld_data(paircounter,12)=r2		 

c=c+2
t=t+2
paircounter=1+paircounter
END DO
a=a+2
b=b+2
sth=sth+1
END DO


!	 mean_LD=REAL(sum(ld_data(:,12)))/REAL(size(ld_data,1))

!     print*,'size(ld_data,1)',size(ld_data,1)
!      print*,'paircounter',paircounter
!	  print*,'mean_LD',mean_LD

CASE (2) !two options 1:pairwise 2:adjacent 
nloci=size(loci_mat,2)/2
!	   	print*,'nloci',nloci

ALLOCATE(freq1_main(1:nloci))
ALLOCATE(freq2_main(1:nloci))
ALLOCATE(seq1(1:nloci))
ALLOCATE(seq2(1:nloci))
ALLOCATE(pair_1(1:size(loci_mat,1),1:2))
ALLOCATE(pair_2(1:size(loci_mat,1),1:2))

!calc frequency

!			counter=1
!   			do i=1,nloci*2,2
!			seq1(counter)=i
!			counter=1+counter
!  			end do

!			counter=1
!		   	do i=2,nloci*2,2
!			seq2(counter)=i
!			counter=1+counter
!		   	end do

!   	   do i=1,nloci
!	     call calc_freq(loci_mat(:,seq1(i):seq2(i)),size(loci_mat,1),outfreq)
!        freq1_main(i)=outfreq
!	     END DO

!	        do i=1,10
!	  print*,'freq1(i)',freq1_main(i)
!            end do


a=1
b=2
c=3
t=4			
paircounter=1
do i=1,nloci-1
pair_1=loci_mat(:,a:b)
pair_2=loci_mat(:,c:t)
call calc_freq(pair_1,size(loci_mat,1),outfreq)
p1=outfreq
q1=1-p1
call calc_freq(pair_2,size(loci_mat,1),outfreq)
p2=outfreq
q2=1-p2

call  calc_hap_freq(pair_1,pair_2,size(loci_mat,1),jj)
hap_11=jj(1)
hap_12=jj(2)
hap_21=jj(3)
hap_22=jj(4)

D=(hap_11*hap_22)-(hap_12*hap_21)
r=D/sqrt(p1*q1*p2*q2)
r2=(D**2)/(p1*q1*p2*q2)

ld_data(paircounter,1)=	paircounter		 

ld_data(paircounter,2)=	p1		 
ld_data(paircounter,3)=	q1	
ld_data(paircounter,4)=	p2		 
ld_data(paircounter,5)=	q2	

ld_data(paircounter,6)=	hap_11	
ld_data(paircounter,7)=	hap_12		 
ld_data(paircounter,8)=	hap_21		 
ld_data(paircounter,9)=	hap_22	

ld_data(paircounter,10)=D		 
ld_data(paircounter,11)=r		 
ld_data(paircounter,12)=r2	

a=a+2
b=b+2
c=c+2
t=t+2			   
paircounter=1+paircounter
END DO


!	 mean_LD=REAL(sum(ld_data(:,12)))/REAL(size(ld_data,1))

!     print*,'size(ld_data,1)',size(ld_data,1)
!      print*,'paircounter',paircounter
!	  print*,'mean_LD',mean_LD




END  SELECT
END SUBROUTINE ld


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!EXTERNAL SUBROUTINES!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!My Subroutines

!  2-SUBROUTINE FOR FReq_calc
SUBROUTINE calc_freq(locus,size1,freq1)
INTEGER :: count1,count2
INTEGER, INTENT(IN) ::size1
INTEGER, INTENT(IN) :: locus(1:size1,1:2)
REAL, INTENT(INOUT) ::freq1

count1= count(locus(:,1)==1)
count2= count(locus(:,2)==1)
freq1=REAL(count1+count2)/REAL(2*size1)


END SUBROUTINE calc_freq


!  3-SUBROUTINE FOR Haplo_Freq
SUBROUTINE calc_hap_freq(loc1,loc2,size1,Freq_hap)
INTEGER :: i
INTEGER, INTENT(IN) ::size1
INTEGER, INTENT(IN) :: loc1(1:size1,1:2),loc2(1:size1,1:2)
REAL, INTENT(INOUT),DIMENSION(1:4) ::Freq_hap
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: Mat_1,Mat_2,haplolist
INTEGER,ALLOCATABLE,DIMENSION(:) :: haplolist2

ALLOCATE(Mat_1(1:size(loc1,1),1:2))
ALLOCATE(Mat_2(1:size(loc1,1),1:2))
ALLOCATE(haplolist(1:size(loc1,1)*2,1:2))
ALLOCATE(haplolist2(1:size(loc1,1)*2))

Mat_1(1:size1,1)=loc1(1:size1,1)
Mat_1(1:size1,2)=loc2(1:size1,1)

Mat_2(1:size1,1)=loc1(1:size1,2)
Mat_2(1:size1,2)=loc2(1:size1,2)

haplolist(1:size1,:)=Mat_1(:,:)
haplolist(size1+1:size(haplolist,1),:)=Mat_2(:,:)

DO i=1,size(haplolist,1)
IF (haplolist(i,1)==1) THEN
haplolist(i,1)=5
END IF
END DO

haplolist2= haplolist(:,1)+haplolist(:,2)


Freq_hap(1)=REAL(count(haplolist2==6))/REAL(size1*2) !11
Freq_hap(2)=REAL(count(haplolist2==7))/REAL(size1*2) !12
Freq_hap(3)=REAL(count(haplolist2==3))/REAL(size1*2) !21
Freq_hap(4)=REAL(count(haplolist2==4))/REAL(size1*2) !22

DEALLOCATE(Mat_1)
DEALLOCATE(Mat_2)
DEALLOCATE(haplolist)
DEALLOCATE(haplolist2)

END SUBROUTINE calc_hap_freq		