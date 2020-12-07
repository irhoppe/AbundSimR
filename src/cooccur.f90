
SUBROUTINE pm(mat, sp, si, rep)
  INTEGER:: si, sp, rep, jj, kswap
  INTEGER:: k1, k2, k3, k4
  DOUBLE PRECISION:: mat(sp,si), rn, a

  jj=0
  kswap=0

  DO WHILE(kswap < rep)
    jj=jj+1
    IF(jj > 100000) GOTO 3
    CALL rng(rn)
    k1=INT(sp*rn)+1
    CALL rng(rn)
    k2=INT(si*rn)+1
  1 CALL rng(rn)
    k3=INT(sp*rn)+1
    IF(k3 == k1) GOTO 1
  2 CALL rng(rn)
    k4=INT(si*rn)+1
    IF(k4 == k2) GOTO 2
    IF(mat(k1,k2) > 0 .AND. mat(k3,k4) > 0) THEN
      kswap=kswap+1
      a=mat(k1,k2)
      mat(k1,k2)=mat(k3,k4)
      mat(k3,k4)=a
    END IF
  END DO
3 RETURN
END

SUBROUTINE pc(mat, sp, si, rep)
  INTEGER:: si, sp, spp, rep, jj, kswap
  INTEGER:: j, k1, k3
  DOUBLE PRECISION:: mat(sp,si), rn, a

  spp=10*sp

  DO 1 j=1, si
    kswap=0
    DO WHILE(kswap < spp)
      jj=0
    2 jj=jj+1
      IF(jj > 10000) GOTO 1
      CALL rng(rn)
      k1=INT(sp*rn)+1
      IF(mat(k1,j) == 0) GOTO 2
    3 jj=jj+1
      IF(jj > 10000) GOTO 1
      CALL rng(rn)
      k3=INT(sp*rn)+1
      IF(mat(k3,j) == 0) GOTO 3
      kswap=kswap+1
      a=mat(k1,j)
      mat(k1,j)=mat(k3,j)
      mat(k3,j)=a
    END DO
1 CONTINUE
  RETURN
END

SUBROUTINE pr(mat, sp, si, rep)
  INTEGER:: si, sp, sii, rep, jj, kswap
  INTEGER:: i, k1, k3
  DOUBLE PRECISION:: mat(sp,si), rn, a

  sii=10*si

  DO 1 i=1, sp
    kswap=0
    DO WHILE(kswap < sii)
      jj=0
    2 jj=jj+1
      IF(jj > 10000) GOTO 1
      CALL rng(rn)
      k1=INT(si*rn)+1
      IF(mat(i,k1) == 0) GOTO 2
    3 jj=jj+1
      IF(jj > 10000) GOTO 1
      CALL rng(rn)
      k3=INT(si*rn)+1
      IF(mat(i,k3) == 0) GOTO 3
      kswap=kswap+1
      a=mat(i,k1)
      mat(i,k1)=mat(i,k3)
      mat(i,k3)=a
    END DO
1 CONTINUE
  RETURN
END

!SUBROUTINE o1(mat, sp, si, indiv)
!  INTEGER:: si, sp, jj, kswap
!  INTEGER:: k1, k2, k3, k4
!  DOUBLE PRECISION:: mat(sp,si), indiv
!  DOUBLE PRECISION:: rn
!
!  jj=0
!  kswap=0
!
!  DO WHILE(kswap < indiv)
!    IF(jj > 1000000) GOTO 4
!  1 CALL rng(rn)
!    k1=INT(sp*rn)+1
!    CALL rng(rn)
!    k2=INT(si*rn)+1
!  2 CALL rng(rn)
!    k3=INT(sp*rn)+1
!    IF(k3 == k1) GOTO 2
!  3 CALL rng(rn)
!    k4=INT(si*rn)+1
!    IF(k4 == k2) GOTO 3
!    jj=jj+1
!    IF(mat(k1,k2) <= 1 .OR. mat(k3,k2) <= 1 .OR. mat(k1,k4) <= 1 .OR. mat(k3,k4) <= 1) GOTO 1
!    kswap=kswap+1
!    IF(a < 0.5) THEN                      ! where does 'a' come from???
!      mat(k1,k2)=mat(k1,k2)+1
!      mat(k1,k4)=mat(k1,k4)-1
!      mat(k3,k2)=mat(k3,k2)-1
!      mat(k3,k4)=mat(k3,k4)+1
!    ELSE
!      mat(k1,k2)=mat(k1,k2)-1
!      mat(k1,k4)=mat(k1,k4)+1
!      mat(k3,k2)=mat(k3,k2)+1
!      mat(k3,k4)=mat(k3,k4)-1
!    END IF
!  END DO
!4 RETURN
!END

SUBROUTINE oa(mat, sp, si)
  INTEGER:: sp, si, i, j, k1
  DOUBLE PRECISION:: mat(sp,si), arow(sp), acol(si)
  DOUBLE PRECISION:: row(sp), col(si), item

  row=0
  col=0
  arow=0
  acol=0
  item=0

  DO 1 j=1, si
    DO 2 i=1, sp
      IF(mat(i,j) > 0) THEN
        acol(j)=acol(j)+mat(i,j)       ! compute row sums
        item=item+mat(i,j)             ! compute total abundance
      END IF
  2 CONTINUE
1 CONTINUE
  k1=0
  DO 3 i=1, sp
    DO 4 j=1, si
      IF(mat(i,j) > 0) THEN
        arow(i)=arow(i)+mat(i,j)       ! compute column sums
        mat(i,j)=1.0                   ! convert abundance matrix to presence-absence matrix
        k1=k1+1                        ! compute total abundance
      END IF
  4 CONTINUE
3 CONTINUE
  DO WHILE(k1 < item)
    CALL prob(arow,sp,ipos)
    CALL prob(acol,si,jpos)
    IF(mat(ipos,jpos) > 0) THEN
      mat(ipos,jpos)=mat(ipos,jpos)+1
      k1=k1+1
    END IF
  END DO
  RETURN
END

SUBROUTINE of(out, sp, si)
  INTEGER:: sp, si, i, j, k1, k2, ipos, jpos
  DOUBLE PRECISION:: mat(sp,si), arow(sp), acol(si)
  DOUBLE PRECISION:: row(sp), col(si), item

  row=0
  col=0
  arow=0
  acol=0
  item=0
  k1=0
  k2=0

  DO 1 j=1, si
    DO 2 i=1, sp
      IF(mat(i,j) > 0) THEN
        acol(j)=acol(j)+mat(i,j)
        item=item+mat(i,j)
      END IF
  2 CONTINUE
1 CONTINUE
  DO 3 i=1, sp
    DO 4 j=1,si
      IF(mat(i,j) > 0) THEN
        arow(i)=arow(i)+mat(i,j)
        mat(i,j)=1.0
        k1=k1+1
      END IF
  4 CONTINUE
3 CONTINUE
  DO 5 i=1, sp
    DO 6 j=1, si
      IF(mat(i,j) > 0) row(i)=row(i)+1
  6 CONTINUE
5 CONTINUE
  DO 7 j=1, si
    DO 8 i=1, sp
      IF(mat(i,j) > 0) col(j)=col(j)+1
    8 CONTINUE
7 CONTINUE
  DO WHILE(k1 < item)
  9 CALL prob(arow,sp,ipos)
    CALL prob(acol,si,jpos)
    IF(mat(ipos,jpos) > 0) THEN
      k2=k2+1
      IF(k2 > 10*item .OR. k1 == item-1) GOTO 10
      IF(row(ipos) >= arow(ipos)) GOTO 9
      IF(col(jpos) >= acol(jpos)) GOTO 9
      row(ipos)=row(ipos)+1
      col(jpos)=col(jpos)+1
      mat(ipos,jpos)=mat(ipos,jpos)+1
      k1=k1+1
    END IF
  END DO
  RETURN
10 CALL labelpr('Trials exeeded maximum',-1)
  RETURN
END

SUBROUTINE itt(mat, sp, si)
  INTEGER:: sp, si, i, j, k1, ipos, jpos
  DOUBLE PRECISION:: mat(sp,si), arow(sp), acol(si)
  DOUBLE PRECISION:: row(sp), col(si), item

  row=0
  col=0
  arow=0
  acol=0
  item=0

  DO 1 j=1, si
    DO 2 i=1, sp
    IF(mat(i,j) > 0) THEN
      acol(j)=acol(j)+mat(i,j)
      item=item+mat(i,j)
    END IF
  2 CONTINUE
1 CONTINUE
  DO 3 i=1, sp
    DO 4 j=1, si
      IF(mat(i,j) > 0) THEN
        arow(i)=arow(i)+mat(i,j)
      END IF
  4 CONTINUE
3 CONTINUE
  mat=0
  k1=0
  DO WHILE(k1 < item)
  5 CALL prob(arow,sp,ipos)
    IF(row(ipos) >= arow(ipos)) GOTO 5
  6 CALL prob(acol,si,jpos)
    IF(col(jpos) >= acol(jpos)) GOTO 6
    row(ipos)=row(ipos)+1
    col(jpos)=col(jpos)+1
    mat(ipos,jpos)=mat(ipos,jpos)+1
    k1=k1+1
  END DO
  RETURN
END

SUBROUTINE ir(mat, sp, si)
  INTEGER:: sp, si, i, j, k1, ipos, jpos
  DOUBLE PRECISION:: mat(sp,si), arow(sp), acol(si)
  DOUBLE PRECISION:: row(sp), item

  arow=0
  acol=0
  row=0
  item=0

  DO 1 j=1, si
    DO 2 i=1, sp
      IF(mat(i,j) > 0) THEN
        acol(j)=acol(j)+mat(i,j)
      END IF
  2 CONTINUE
1 CONTINUE
  DO 3 i=1, sp
    DO 4 j=1, si
      IF(mat(i,j) > 0) THEN
        arow(i)=arow(i)+mat(i,j)
      END IF
  4 CONTINUE
    DO 5 j=1, si
      IF(mat(i,j) > 0) THEN
        item=item+1
        GOTO 3
      END IF
  5 CONTINUE
3 CONTINUE
  mat=0
  k1=0
  DO WHILE(1 > -1)
    CALL prob(arow,sp,ipos)
    CALL prob(acol,si,jpos)
    mat(ipos,jpos)=mat(ipos,jpos)+1
    IF(row(ipos) == 0) THEN
      row(ipos)=1
      k1=k1+1
      IF(k1 >= item) RETURN
    END IF
  END DO
  RETURN
END

SUBROUTINE itr(mat, sp, si)
  INTEGER:: sp, si, i, j, k1, jpos
  DOUBLE PRECISION:: mat(sp,si), arow(sp), acol(si)

  acol=0
  arow=0

  DO 1 j=1, si
    DO 2 i=1, sp
      IF(mat(i,j) > 0) THEN
        acol(j)=acol(j)+mat(i,j)
      END IF
  2 CONTINUE
1 CONTINUE
  DO 3 i=1, sp
    DO 4 j=1, si
      arow(i)=arow(i)+mat(i,j)
  4 CONTINUE
3 CONTINUE
  mat=0
  DO 5 i=1, sp
    k1=0
      DO WHILE(k1 < arow(i))
        CALL prob(acol,si,jpos)
        mat(i,jpos)=mat(i,jpos)+1
        k1=k1+1
      END DO
5 CONTINUE
  RETURN
END

SUBROUTINE itc(mat, sp, si)
  INTEGER sp, si, i, j, k1, ipos
  DOUBLE PRECISION:: mat(sp,si), arow(sp), acol(si)

  acol=0
  arow=0

  DO 1 j=1, si
    DO 2 i=1, sp
      IF(mat(i,j) > 0) THEN
        acol(j)=acol(j)+mat(i,j)
      END IF
  2 CONTINUE
1 CONTINUE
  DO 3 i=1, sp
    DO 4 j=1, si
      arow(i)=arow(i)+mat(i,j)
  4 CONTINUE
3 CONTINUE
  mat=0
  DO 5 j=1, si
    k1=0
    DO WHILE(k1 < acol(j))
      CALL prob(arow,sp,ipos)
      mat(ipos,j)=mat(ipos,j)+1
      k1=k1+1
    END DO
5 CONTINUE
  RETURN
END

SUBROUTINE isr(mat, sp, si)
  INTEGER:: sp, si, i, j, k1, jpos
  DOUBLE PRECISION:: mat(sp,si), arow(sp), acol(si)
  DOUBLE PRECISION:: row(sp)

  row=0
  acol=0
  arow=0

  DO 1 j=1, si
    DO 2 i=1, sp
      IF(mat(i,j) > 0) THEN
        acol(j)=acol(j)+mat(i,j)
      END IF
  2 CONTINUE
1 CONTINUE
  DO 3 i=1, sp
    DO 4 j=1, si
      arow(i)=arow(i)+mat(i,j)
      row(i)=row(i)+1
  4 CONTINUE
3 CONTINUE
  mat=0
  DO 5 i=1, sp
    k1=0
    DO WHILE(k1 < row(i))
      CALL prob(acol,si,jpos)
      IF(mat(i,jpos) == 0) k1=k1+1
      mat(i,jpos)=mat(i,jpos)+1
    END DO
5 CONTINUE
  RETURN
END

SUBROUTINE isc(mat, sp, si)
  INTEGER:: sp, si, i, j, k1, ipos
  DOUBLE PRECISION:: mat(sp,si), arow(sp), acol(si)
  DOUBLE PRECISION:: col(si)

  col=0
  acol=0
  arow=0

  DO 1 j=1, si
    DO 2 i=1, sp
      IF(mat(i,j) > 0) THEN
        acol(j)=acol(j)+mat(i,j)
        col(j)=col(j)+1
      END IF
  2 CONTINUE
1 CONTINUE
  DO 3 i=1, sp
    DO 4 j=1, si
      arow(i)=arow(i)+mat(i,j)
  4 CONTINUE
3 CONTINUE
  mat=0
  DO 5 j=1, si
    k1=0
    DO WHILE(k1 < col(j))
      CALL prob(arow,sp,ipos)
      IF(mat(ipos,j) == 0) k1=k1+1
      mat(ipos,j)=mat(ipos,j)+1
    END DO
5 CONTINUE
  RETURN
END

SUBROUTINE ia(mat, sp, si)
  INTEGER:: sp, si, i, j, k1, ipos, jpos
  DOUBLE PRECISION:: mat(sp,si), arow(sp), acol(si)
  DOUBLE PRECISION:: item

  arow=0
  acol=0
  item=0

  DO 1 j=1, si
    DO 2 i=1, sp
      IF(mat(i,j) > 0) THEN
        acol(j)=acol(j)+mat(i,j)
        item=item+mat(i,j)
      END IF
  2 CONTINUE
1 CONTINUE
  DO 3 i=1, sp
    DO 4 j=1, si
      IF(mat(i,j) > 0) THEN
        arow(i)=arow(i)+mat(i,j)
      END IF
  4 CONTINUE
3 CONTINUE
  mat=0
  k1=0
  DO WHILE(k1 < item)
    CALL prob(arow,sp,ipos)
    CALL prob(acol,si,jpos)
    mat(ipos,jpos)=mat(ipos,jpos)+1
    k1=k1+1
  END DO
  RETURN
END

SUBROUTINE is(mat, sp, si)
  INTEGER:: sp, si, i, j, ipos, jpos, k1, k2, k3
  DOUBLE PRECISION:: mat(sp,si), arow(sp), acol(si)
  DOUBLE PRECISION:: srow(sp), scol(si), row(sp), col(si)
  DOUBLE PRECISION:: item, dens

  row=0
  col=0
  arow=0
  acol=0
  srow=0
  scol=0
  item=0
  dens=0

  DO 1 j=1, si
    DO 2 i=1, sp
      IF(mat(i,j) > 0) THEN
        acol(j)=acol(j)+mat(i,j)
        scol(j)=scol(j)+1
        dens=dens+mat(i,j)
        item=item+1
      END IF
  2 CONTINUE
1 CONTINUE
  DO 3 i=1, sp
    DO 4 j=1, si
      IF(mat(i,j) > 0) THEN
        arow(i)=arow(i)+mat(i,j)
        srow(i)=srow(i)+1
      END IF
  4 CONTINUE
3 CONTINUE
  mat=0
  row=0
  col=0
  k1=0
  k2=0
  k3=0
  DO WHILE(k3 < item)
  6 CALL prob(arow,sp,ipos)
    CALL prob(acol,si,jpos)
    k2=k2+1
    IF(k2 > 100*dens) RETURN
    IF(mat(ipos,jpos) == 0) THEN
      IF(row(ipos) >= srow(ipos)) GOTO 6
      IF(col(jpos) >= scol(jpos)) GOTO 6
      row(ipos)=row(ipos)+1
      col(jpos)=col(jpos)+1
      mat(ipos,jpos)=mat(ipos,jpos)+1
      k1=k1+1
      k3=k3+1
    ELSE
      mat(ipos,jpos)=mat(ipos,jpos)+1
      k1=k1+1
    END IF
  END DO
  RETURN
END

SUBROUTINE isa(mat, sp, si)
  INTEGER:: sp, si, i, j, ipos, jpos, k1, k2, k3
  DOUBLE PRECISION:: mat(sp,si), arow(sp), acol(si)
  DOUBLE PRECISION:: srow(sp), scol(si), row(sp), col(si)
  DOUBLE PRECISION:: item, dens

  row=0
  col=0
  arow=0
  acol=0
  srow=0
  scol=0
  item=0
  dens=0

  DO 1 j=1, si
    DO 2 i=1, sp
      IF(mat(i,j) > 0) THEN
        acol(j)=acol(j)+mat(i,j)
        scol(j)=scol(j)+1
        dens=dens+mat(i,j)
        item=item+1
      END IF
  2 CONTINUE
1 CONTINUE
  DO 3 i=1, sp
    DO 4 j=1, si
      IF(mat(i,j) > 0) THEN
        arow(i)=arow(i)+mat(i,j)
        srow(i)=srow(i)+1
      END IF
  4 CONTINUE
3 CONTINUE
  mat=0
  row=0
  col=0
  k1=0
  k2=0
  k3=0
  DO WHILE(k3 < item)
  6 CALL prob(arow,sp,ipos)
    CALL prob(acol,si,jpos)
    k2=k2+1
    IF(k2 > 100*dens) RETURN
    IF(mat(ipos,jpos) == 0) THEN
      IF(row(ipos) >= srow(ipos)) GOTO 6
      IF(col(jpos) >= scol(jpos)) GOTO 6
      row(ipos)=row(ipos)+1
      col(jpos)=col(jpos)+1
      mat(ipos,jpos)=mat(ipos,jpos)+1
      k1=k1+1
      k3=k3+1
    ELSE
      IF(k1 >= dens) RETURN
      mat(ipos,jpos)=mat(ipos,jpos)+1
      k1=k1+1
    END IF
  END DO
  RETURN
END

SUBROUTINE iff(mat, sp, si, replic)
  INTEGER:: sp, si, replic, i, j, ipos, jpos
  INTEGER:: k1, k2
  DOUBLE PRECISION:: mat(sp,si), arow(sp), acol(si)
  DOUBLE PRECISION:: row(sp), col(si), item

  row=0
  col=0
  arow=0
  acol=0
  item=0
  k1=0
  k2=0

  DO 1 j=1, si
    DO 2 i=1, sp
      IF(mat(i,j) > 0) THEN
        acol(j)=acol(j)+mat(i,j)
        item=item+mat(i,j)
      END IF
  2 CONTINUE
1 CONTINUE
  DO 3 i=1, sp
    DO 4 j=1, si
      IF(mat(i,j) > 0) THEN
        arow(i)=arow(i)+mat(i,j)
        mat(i,j)=1.0
        k1=k1+1
      END IF
  4 CONTINUE
3 CONTINUE
  DO 5 i=1, sp
    DO 6 j=1, si
      IF(mat(i,j) > 0) row(i)=row(i)+1
  6 CONTINUE
5 CONTINUE
  DO 7 j=1, si
    DO 8 i=1, sp
      IF(mat(i,j) > 0) col(j)=col(j)+1
  8 CONTINUE
7 CONTINUE
  CALL swap(mat,sp,si,replic)
  DO WHILE(k1 < item)
  9 CALL prob(arow,sp,ipos)
    CALL prob(acol,si,jpos)
    IF(mat(ipos,jpos) > 0) THEN
      k2=k2+1
      IF(k2 > 1000*item .OR. k1 == item-1) GOTO 10
      IF(row(ipos) >= arow(ipos)) GOTO 9
      IF(col(jpos) >= acol(jpos)) GOTO 9
      row(ipos)=row(ipos)+1
      col(jpos)=col(jpos)+1
      mat(ipos,jpos)=mat(ipos,jpos)+1
      k1=k1+1
    END IF
  END DO
  RETURN
10 CALL labelpr('Trials exceeded maximum',-1)
  RETURN
END

!! Helpers !!

SUBROUTINE rng(rn)
  DOUBLE PRECISION:: rn, frunif
  EXTERNAL:: rndstart, rndend, frunif

  CALL rndstart()
  rn=frunif()
  CALL rndend()

  RETURN
END

SUBROUTINE prob(in, sp, pos)
  INTEGER:: sp, pos
  DOUBLE PRECISION:: an, rn, in(sp), probf(sp)

  probf=0
  probf(1)=in(1)

  DO 100 i=2, sp
    probf(i)=in(i)+probf(i-1)
  100 CONTINUE
  CALL rng(rn)
  an=rn*probf(sp)
  DO 101 pos=1, sp
    IF(an < probf(pos)) RETURN
  101 CONTINUE
  pos=0
  RETURN
END

SUBROUTINE swap(mat, sp, si, rep)
  INTEGER:: si, sp, rep, jj, k1, k2, k3, k4, kswap
  DOUBLE PRECISION:: mat(sp,si), rn

  jj=0
  kswap=0

  DO WHILE(kswap < rep)
    jj=jj+1
    IF(jj > 100000) GOTO 200
    CALL rng(rn)
    k1=INT(sp*rn)+1
    CALL rng(rn)
    k2=INT(si*rn)+1
201 CALL rng(rn)
    k3=INT(sp*rn)+1
    IF(k3 ==k1) GOTO 201
202 CALL rng(rn)
    k4=INT(si*rn)+1
    IF(k4 == k2) GOTO 202
    IF((mat(k1,k2) /= 0 .AND. mat(k3,k4) /= 0) .OR. (mat(k3,k2) /= 0 .AND. mat(k1,k4) /= 0)) THEN
      IF((mat(k1,k2) == 0 .AND. mat(k3,k4) == 0) .OR. (mat(k3,k2) == 0 .AND. mat(k1,k4) == 0)) THEN
        kswap=kswap+1
        a=mat(k1,k2)
        mat(k1,k2)=mat(k1,k4)
        mat(k1,k4)=a
        a=mat(k3,k2)
        mat(k3,k2)=mat(k3,k4)
        mat(k3,k4)=a
      END IF
    END IF
  END DO
200 RETURN
END

