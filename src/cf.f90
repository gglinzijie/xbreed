SUBROUTINE cf(r_size,c_size,loci_mat,nloci,freq1_main)

!Input official
INTEGER,INTENT(IN)::r_size,c_size,nloci
INTEGER,INTENT(IN), DIMENSION(1:r_size,1:c_size):: loci_mat
!Output official
DOUBLE PRECISION ,DIMENSION(1:nloci),INTENT(INOUT) :: freq1_main
!Local variables
INTEGER::i,counter,hsize
REAL :: outfreq
INTEGER,ALLOCATABLE,DIMENSION(:)::seq1,seq2


ALLOCATE(seq1(1:nloci))
ALLOCATE(seq2(1:nloci))
!calc frequency

counter=1
do i=1,nloci*2,2
seq1(counter)=i
counter=1+counter
end do

counter=1
do i=2,nloci*2,2
seq2(counter)=i
counter=1+counter
end do

hsize=size(loci_mat,1)

do i=1,nloci
call freq_calc(loci_mat(:,seq1(i):seq2(i)),hsize,outfreq)
freq1_main(i)=outfreq
END DO

DEALLOCATE(seq1)
DEALLOCATE(seq2)

END SUBROUTINE cf
!My Subroutines


!  2-SUBROUTINE FOR freq_calc
SUBROUTINE freq_calc(locus,size1,freq1)
INTEGER :: count1,count2
INTEGER, INTENT(IN) ::size1
INTEGER, INTENT(IN) :: locus(1:size1,1:2)
REAL, INTENT(INOUT) ::freq1

count1= count(locus(:,1)==1)
count2= count(locus(:,2)==1)
freq1=REAL(count1+count2)/REAL(2*size1)


END SUBROUTINE freq_calc

