
subroutine pm(mat,sp,si,rep)
integer si,sp,rep
integer mat(sp,si)
double precision rn
jj=0
kswap=0
do while (kswap.lt.rep)
  jj=jj+1
  if(jj.gt.100000)goto 31
  call rng(rn)
  k1=int(sp*rn)+1
  call rng(rn)
  k2=int(si*rn)+1
  35 call rng(rn)
  k3=int(sp*rn)+1
  if(k3.eq.k1)goto 35
  36 call rng(rn)
  k4=int(si*rn)+1
  if(k4.eq.k2)goto 36
  if(mat(k1,k2).gt.0.and.mat(k3,k4).gt.0)then
    kswap=kswap+1
    a=mat(k1,k2)
    mat(k1,k2)=mat(k3,k4)
    mat(k3,k4)=a
  endif
enddo
31 return
end

subroutine rng(rn)
  double precision rn
  external rndstart, rndend
  double precision frunif
  external frunif

call rndstart()
rn = frunif()
call rndend()

return
end

