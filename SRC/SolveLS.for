CExample BLAS, LAPACK thread-safe version
implicit real*8 (a-h,o-z)
parameter(maxn=200,m=80,mblk=4,k=maxn)
parameter(zero=0.0d0,one=1.0d0)
real*8 a(k,maxn),aa(k,maxn),x(k,m),b(k,m)
integer ip(maxn)
C
===========================================================
C
Define the matrix
C
===========================================================
n=maxn
call inita(a,k,n)
do i=1,n
do j=1,n
aa(j,i)=a(j,i)
end do
end do
C
===========================================================
C
LU decomposition
C
===========================================================
call dgetrf(n,n,a,k,ip,info)
!$OMP PARALLEL PRIVATE(mb,info)
!$OMP DO
do i=1,m,mblk
mb=min(mblk,m-i+1)
C
===========================================================
C
Defeine the vectors
C
===========================================================
do jm=1,mb
do jn=1,n
x(jn,jm+i-1)=jn+i+jm-1
end do
end do
call dgemm('N','N',n,mb,n,one,aa,k,x(1,i),k,zero
&
,b(1,i),k)
C
===========================================================
C
Solution
C
===========================================================
call dgetrs('N',n,mb,a,k,ip,b(1,i),k,info)
11if(info.ne.0) then
write(6,*) 'error in dgetrs
info = ',info
stop
end if
end do
!$OMP END PARALLEL
C
===========================================================
C
Check result
C
===========================================================
call check(a,b,k,n,m)
end
