!--------------------------------------------------------------------------------------
    !http://www.mpi-hd.mpg.de/astrophysik/HEA/internal/Numerical_Recipes/f6-2.pdf
!--------------------------------------------------------------------------------------
    !Cumulative Poisson Probability Function
    !Px(< k), for positive x and integer k >= 1, denotes the cumulative Poisson
    !probability function. It is defined as the probability that the number of Poisson
    !random events occurring will be between 0 and k - 1 inclusive, if the expected mean
    !number is x. It has the limiting values
    !    Px(< 1) = e-x
    !    Px(< inf) = 1
    !Its relation to the incomplete gamma function is simply
    !    Px(< k) = Q(k,x) = gammq (k,x)
    !Therefore, the probability of events >= k is given by 
    !    Px(>=k) = P(k,x) = 1 - gammq(k,x) = gammp(k,x)
!--------------------------------------------------------------------------------------
	program uvProb
	CHARACTER*10 zku,zkn
    INTEGER ku,kn
    DOUBLE PRECISION gammp,r,p	
	call getarg(1, zku) !the uv sample count
	call getarg(2, zkn) !the non-uv sample count
	read (zku, *) ku
	read (zkn, *) kn
	r = dble(ku + kn) / 2 !expected rate = avg of two counts
	p = gammp(dble(ku),r)
    write(*,1)p !the probability of having this many or more UV hits
1   format(F22.20) !20 digits after the decimal place i.e. 0.xxxxxxxxxxxxxxxxxxxx
    end
	
    !Returns the incomplete gamma function P(a, x).
    DOUBLE PRECISION FUNCTION gammp(a,x)
    DOUBLE PRECISION a,gammp,x
    DOUBLE PRECISION gammcf,gamser,gln
    if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
    else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
    endif
    END

    SUBROUTINE gser(gamser,a,x,gln)
    INTEGER ITMAX
    DOUBLE PRECISION a,gamser,gln,x,EPS
    PARAMETER (ITMAX=1000,EPS=3.e-7)
    INTEGER n
    DOUBLE PRECISION ap,del,sum,gammln
    gln=gammln(a)
    ap=a
    sum=1./a
    del=sum
    do n=1,ITMAX
    ap=ap+1.
    del=del*x/ap
    sum=sum+del
    if(abs(del).lt.abs(sum)*EPS)goto 1
    end do
1   gamser=sum*exp(-x+a*log(x)-gln)
    END

    SUBROUTINE gcf(gammcf,a,x,gln)
    INTEGER ITMAX
    DOUBLE PRECISION a,gammcf,gln,x,EPS,FPMIN
    PARAMETER (ITMAX=1000,EPS=3.e-7,FPMIN=1.e-30)
    INTEGER i
    DOUBLE PRECISION an,b,c,d,del,h,gammln
    gln=gammln(a)
    b=x+1.-a 
    c=1./FPMIN 
    d=1./b
    h=d
    do i=1,ITMAX 
    an=-i*(i-a)
    b=b+2.
    d=an*d+b
    if(abs(d).lt.FPMIN)d=FPMIN
    c=b+an/c
    if(abs(c).lt.FPMIN)c=FPMIN
    d=1./d
    del=d*c
    h=h*del
    if(abs(del-1.).lt.EPS)goto 1
    end do
1   gammcf=exp(-x+a*log(x)-gln)*h 
    END

    DOUBLE PRECISION FUNCTION gammln(xx)
    DOUBLE PRECISION xx
    INTEGER j
    DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
    SAVE cof,stp
    DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
     24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
     -.5395239384953d-5,2.5066282746310005d0/
    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do j=1,6
    y=y+1.d0
    ser=ser+cof(j)/y
    end do
    gammln=tmp+log(stp*ser/x)
    END

