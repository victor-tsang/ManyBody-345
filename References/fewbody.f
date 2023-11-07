
       program burrau

c..    Three body integrator from cartesian ic's
c..    currently set up to solve pythagorean problem.

c..    Uses the Bulirsch-Stoer method.

       implicit real*8(a-h,o-z), integer*4(i-n)

       parameter (nb=3)
       parameter(n=6*nb,pi=3.14159265358979323846)


       dimension y(n),dydx(n),yscal(n)
       dimension abod(nb),ebod(nb)

       common /path2/ rm(nb),vcx,vcy,vcz,dorig

       external derivs


c..    Overall accuracy criterion
       Accuracy=1.e-6
c..    Timestep accuracy for Bulirsch-Stoer
       Acci=1.e-14
c..    number of phase space dimensions
       nvar=n
c..    number of scattering trials
       nloop=1
c..    reintegration timescale (years)
       x2=100

c..    mass of body 1
       rM(1)=3.
c..    mass of body 2
       rM(2)=4.
c..    mass of body 3
       rM(3)=5.


c..    starting time for the integration
       x1=0.
   
         Acc=Acci
2110     continue
         
         
c..        initial conditions for the bodies:

c..        body 1 x
           y(1)=1.

c..        body 1 vx
           y(2)=0.

c..        body 1 y
           y(3)=3.

c..        body 1 vy
           y(4)=0.

c..        body 1 z
           y(5)=0.

c..        body 1 vz
           y(6)=0.

c..        and so on for the other bodies.

           y(7)=-2.
           y(8)=0.
           y(9)=-1.
           y(10)=0.
           y(11)=0.
           y(12)=0.
           
           y(13)=1.
           y(14)=0.
           y(15)=-1.
           y(16)=0.
           y(17)=0.
           y(18)=0.

2109     continue

c..      number of integration steps in current trial
         iflag=0


c..      Compute center of mass velocity

         rmtot=0.
         vcx=0.
         vcy=0.
         vcz=0.

         do i=1,nb
           ib=(i-1)*6
           rmtot=rmtot+rm(i)
           vcx=vcx+rm(i)*y(ib+2)
           vcy=vcy+rm(i)*y(ib+4)
           vcz=vcz+rm(i)*y(ib+6)
         end do

         vcx=vcx/rmtot
         vcy=vcy/rmtot
         vcz=vcz/rmtot

       
c..      Determine System Energy
         call EnergySum(y,Energy)

         Eorig=Energy
         x=x1
         iescape=0
         isscape=0
         ibelong1=1
         ibelong2=1
         istop=0

c..      Overall integration loop
2105     continue
           iflag=iflag+1

           if(iflag.eq.1) then
c..          Determine the initial timestep
             htry=0.001
           end if

c..        Use current timestep to take a Bulirsch-Stoer Integration Step
           call derivs(x,y,dydx)

           do iscale=1,nvar
             yscal(iscale)=abs(y(iscale))+abs(htry*dydx(iscale))+1.e-30
           end do


           call bsstep(y,dydx,nvar,x,htry,acc,yscal,hdid,hnext,derivs)  
           htry=hnext/2
     
c..        check for printout, or an unbound planet
           if(mod(iflag,1).eq.0) then
      
            if(x.gt.0.and.x.lt.70) then
            call printout(x,y)
            end if

c..          Check conservation
             call EnergySum(y,Energy)
             Efrac=Energy/Eorig
             check=abs(1.-Efrac)
             if(check.gt.Accuracy) then
               write(*,*) x,iflag,'accuracy stop'
               stop
             end if


c..          check if number of integration steps exceeded:
             if(iflag.gt.2000) then

               stop

             end if

c..          check if physical time exceeded:
             if(x.gt.x2) then

               stop

             end if

c..          Outcome checklist finished
           end if

         goto 2105
c..      End of overall integration loop

2106     continue

c..      Successful termination 


       end



      subroutine derivs(x,y,dydx)

      implicit real*8(a-h,o-z), integer*4(i-n)

      parameter (nb=3)
      parameter(n=6*nb,pi=3.14159265358979323846)


      common /path2/ rm(nb),vcx,vcy,vcz,dorig


      real*8 x,y(n),dydx(n)
      real*8 denom(nb,nb)

      

      do i=2,n,2
        dydx(i-1)=y(i)
      end do

      do i=1,nb
        do j=i+1,nb

          if (i.ne.j) then
            jb=(j-1)*6
            ib=(i-1)*6
            denom(i,j)=(y(jb+1)-y(ib+1))*(y(jb+1)-y(ib+1)) +
     +                 (y(jb+3)-y(ib+3))*(y(jb+3)-y(ib+3)) +
     +                 (y(jb+5)-y(ib+5))*(y(jb+5)-y(ib+5))

            denom(i,j)=denom(i,j)*sqrt(denom(i,j))
            denom(j,i)=denom(i,j)
          end if

        end do
      end do

2104  continue
 
      do i=1,nb
        ib=(i-1)*6
        do ic=1,3
          dydx(ib+(2*ic))=0.
          do j=1,nb
            jb=(j-1)*6
            if( i.ne.j) then
              dydx(ib+(2*ic))=dydx(ib+(2*ic)) -
     +        rm(j)*(y(ib+2*ic-1)-y(jb+2*ic-1))/denom(i,j)
            end if
          end do
        end do
      end do

 
      return
      end


c...................................................................

      subroutine EnergySum(y,energy)

c...................................................................

       implicit real*8(a-h,o-z), integer*4(i-n)

       parameter (nb=3)
       parameter(n=6*nb,pi=3.14159265358979323846)


      common /path2/ rm(nb),vcx,vcy,vcz,dorig
      real*8 y(6*nb)

      dimension denom(nb,nb)
      do i=1,nb
        do j=1,nb
          denom(i,j)=1.e+10
        end do
      end do
      do i=1,nb
        ib=(i-1)*6
        do j=1,nb
          jb=(j-1)*6
          if (i.ne.j) then
            denom(i,j) = (y(jb+1)-y(ib+1))**2
     +                 +(y(jb+3)-y(ib+3))**2
     +                 +(y(jb+5)-y(ib+5))**2
            denom(i,j)=sqrt(denom(i,j))
          end if
        end do
      end do

      Energy=0.
      do i=1,nb
        ib=(i-1)*6
        do j=1,nb
          if(i.ne.j) then
            Energy=Energy-0.5*rm(i)*rm(j)/denom(i,j)
          end if
        end do
        do ic=1,3
          Energy=Energy+0.5*rm(i)*y(ib+2*ic)**2
        end do
      end do

       return
       end


c..................................................................
c
      subroutine printout(x,y)
c
c..................................................................

       implicit real*8(a-h,o-z), integer*4(i-n)

       parameter (nb=3)
       parameter(n=6*nb,pi=3.14159265358979323846)


      common /path2/ rm(nb),vcx,vcy,vcz,dorig
      real*8 y(6*nb)

c..   write position data to files

      do ib=1,nb
        jb=(ib-1)*6
        ifb=10+ib
        write(ifb,2104) x,y(jb+1)-x*vcx,y(jb+3)-x*vcy,y(jb+5)-x*vcz  
      end do

2104       format(4e13.5)
      return
      end


c....................................................................
c
      subroutine bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs) 
c
c....................................................................


c..   take a bulirsch-stoer timestep (numerical recipes)

      implicit real*8(a-h,o-z), integer*4(i-n)
      integer*4 nv,nmax,kmaxx,imax
      real*8 eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv),safe1,
     +       safe2,redmax,redmin,tiny,scalmx
      parameter (nmax=50,kmaxx=8,imax=kmaxx+1,safe1=0.25,safe2=0.7,
     +       redmax=1.e-5,redmin=0.7,tiny=1.e-30,scalmx=0.1)

    
      integer*4 i,iq,k,kk,km,kmax,kopt,nseq(imax)
      real*8 eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest,
     +       xnew,a(imax),alf(kmaxx,kmaxx),err(kmaxx),yerr(nmax),
     +       ysav(nmax),yseq(nmax)
      logical first,reduct
      external derivs
      data first/.true./,epsold/-1./
      data nseq /2,4,6,8,10,12,14,16,18/
      save a,alf,epsold, first,kmax,kopt,nseq,xnew

      if(eps.ne.epsold)then
        hnext=-1.e29
        xnew=-1.e29
        eps1=safe1*eps
        a(1)=nseq(1)+1
        do k=1,kmaxx
          a(k+1)=a(k)+nseq(k+1)
        end do
        do iq=2,kmaxx
          do k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/
     +      ((a(iq+1)-a(1)+1.)*(2*k+1)))
          end do
        end do
        epsold=eps
        do kopt=2,kmaxx-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
        end do
1       kmax=kopt
      end if
      h=htry
      do i=1,nv
        ysav(i)=y(i)
      end do
      if(h.ne.hnext.or.x.ne.xnew) then
        first=.true.
        kopt=kmax
      end if
      reduct=.false.
2     do k=1,kmax
        xnew=x+h
        if(xnew.eq.x) pause 'step size underflow in bsstep'
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1) then
          errmax=TINY
          do i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
          end do
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/safe1)**(1./(2*km+1))
        end if
        if(k.ne.1.and.(k.ge.kopt-1.or.first)) then
          if(errmax.lt.1.) goto 4
          if(k.eq.kmax.or.k.eq.kopt+1) then
            red=safe2/err(km)
            goto 3
          else if(k.eq.kopt) then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1./err(km)
              goto 3
            end if
          else if(kopt.eq.kmax) then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*
     +          safe2/err(km)
               goto 3  
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          end if
        end if
      end do
3     red=min(red,redmin)
      red=max(red,redmax)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.e35
      do kk=1,km
        fact=max(err(kk),scalmx)
        work=fact*a(kk+1)
        if(work.lt.wrkmin) then
          scale=fact
          wrkmin=work
          kopt=kk+1
        end if
      end do
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),scalmx)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      end if
      return
      end

c....................................................................
c
      subroutine pzextr(iest,xest,yest,yz,dy,nv)
c
c....................................................................

c..   polynoical extrapolation (numerical recipes)
c..   (called by bsstep bulirsch-stoer integrator)

      implicit real*8(a-h,o-z), integer*4(i-n)
      integer*4 iest,nv,imax,nmax
      real*8 xest,dy(nv),yest(nv),yz(nv)
      parameter (imax=13,nmax=50)
      integer*4 j,k1
      real*8 delta,f1,f2,q,d(nmax),qcol(nmax,imax),x(imax)
      save qcol,x
      x(iest)=xest
      do j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
      end do
      if(iest.eq.1) then
        do j=1,nv
          qcol(j,1)=yest(j)
        end do
      else
        do j=1,nv
          d(j)=yest(j)
        end do
        do k1=1,iest-1
          delta=1./(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
          end do
        end do
        do j=1,nv
          qcol(j,iest)=dy(j)
        end do
      end if
      return
      end

c.......................................................................
c
      subroutine mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
c
c.......................................................................

c..   modified midpoint method to take one integration step (numerical recipes)
c..   (called by bsstep bulirsch-stoer integrator)

      integer*4 nstep,nvar,nmax
      real*8 htot,xs,dydx(nvar),y(nvar),yout(nvar)
      external derivs
      parameter (nmax=50)
      integer*4 i,n
      real*8 h,h2,swap,x,ym(nmax),yn(nmax)
      h=htot/nstep
      do i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
      end do
      x=xs+h
      call derivs(x,yn,yout)
      h2=2.*h
      do n=2,nstep
        do i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
        end do
        x=x+h
        call derivs(x,yn,yout)
      end do
      do i=1,nvar
        yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
      end do
      return
      end
