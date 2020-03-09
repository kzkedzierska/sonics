C written by amazing people from PyMC team: https://github.com/pymc-devs/pymc
C Chris Fonnesbeck, Anand Patil, David Huard, John Salvatier
C copied and compiled here, because PyMC is no longer maintained, and PyMC3
C lacks mvhyperg

      DOUBLE PRECISION FUNCTION gammln(xx)
C Returns the value ln[gamma(xx)] for xx > 0.

      DOUBLE PRECISION xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)

C Internal arithmetic will be done in double precision,
C a nicety that you can omit if five-figure accuracy is good enough.

      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     +24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     +-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*dlog(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
         y=y+1.d0
         ser=ser+cof(j)/y
      enddo
      gammln=tmp+dlog(stp*ser/x)
      return
      END


      DOUBLE PRECISION FUNCTION factln(n)
C USES gammln Returns ln(n!).

      INTEGER n
      DOUBLE PRECISION a(100),gammln, pass_val
      DOUBLE PRECISION infinity
      PARAMETER (infinity = 1.7976931348623157d308)

      SAVE a
C Initialize the table to negative values.
      DATA a/100*-1./
      pass_val = n + 1
      if (n.lt.0) then
c        write (*,*) 'negative factorial in factln'
        factln=-infinity
        return
      endif
C In range of the table.
      if (n.le.99) then
C If not already in the table, put it in.
        if (a(n+1).lt.0.) a(n+1)=gammln(pass_val)
        factln=a(n+1)
      else
C Out of range of the table.
        factln=gammln(pass_val)
      endif
      return
      END

      SUBROUTINE mvhyperg(x,color,k,like)

c Multivariate hypergeometric log-likelihood function
c Using the analogy of an urn filled with balls of different colors,
c the mv hypergeometric distribution describes the probability of
c drawing x(i) balls of a given color.
c
c x : (array) Number of draws for each color.
c color : (array) Number of balls of each color.

c Total number of draws = sum(x)
c Total number of balls in the urn = sum(color)

cf2py integer dimension(k),intent(in) :: x,color
cf2py integer intent(hide),depend(x) :: k=len(x)
cf2py double precision intent(out) :: like
cf2py threadsafe

      INTEGER x(k),color(k)
      INTEGER d,total,i,k
      DOUBLE PRECISION like
      DOUBLE PRECISION factln
      DOUBLE PRECISION infinity
      PARAMETER (infinity = 1.7976931348623157d308)

      total = 0
      d = 0
      like = 0.0
      do i=1,k
c Combinations of x balls of color i
        like = like + factln(color(i))-factln(x(i))
     +-factln(color(i)-x(i))
        if ((color(i) .LT. 0.0) .OR. (x(i) .LT. 0.0)) then
          like = -infinity
          RETURN
        endif
        d = d + x(i)
        total = total + color(i)
      enddo
      if (total .LE. 0.0) then
        like = -infinity
        RETURN
      endif
c Combinations of d draws from total
      like = like - (factln(total)-factln(d)-factln(total-d))
      return
      END
