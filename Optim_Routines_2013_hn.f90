!==========================================================
!
!      Joint Latent class mixed model for continuous
!          Gaussian outcome
!
!        Cecile Proust, Helene Jacqmin-Gadda
!
!
!       Corresponding author :
!       Cecile Proust, INSERM U897, ISPED,
!       146  rue L\'eo Saignat,
!       33076 Bordeaux cedex, France.
!       Tel: (33) 5 57 57 45 79; Fax: (33) 5 56 24 00 81;
!       e-mail : cecile.proust@isped.u-bordeaux2.fr
!
!                                       21/12/2010

!
! Revu en 2013 pour Hessienne sans splines
!
!
!===========================================================
! - Version fortran 90
!----------- PROGRAM COMPOSITION --------------------------
!    - subroutine MARQ9
!    - subroutine DERIVA
!    - subroutine SEARPAS
!    - subroutine VALFPA
!    - subroutine MAXT
!    - subroutine DCHOLE
!    - subroutine DMFSD
!    - subroutine DSINV
!----------------------------------------------------------
!
!     INTERFACE TYPEc
!
!----------------------------------------------------------


!================================================
!  Module declaratif et Initialisation pour le MPI
!
!
!================================================
module ParametresPourParallelisation
  
!   include 'mpif.h'
  use mpi
  integer ::info
  integer ::CurrentProcessorID
  integer ::NombreDeProcesseurs
  
end module ParametresPourParallelisation
!=================================================
!
!
!=================================================


module type	
  
  interface verif1   
     subroutine marq98(b,m,ni,v,rl,ier,istop,ca,cb,dd,namefunc)
       integer,intent(in) :: m
       integer,intent(inout)::ni,ier,istop
       double precision,dimension(m*(m+3)/2),intent(out)::v
       double precision,intent(out)::rl
       double precision,dimension(m),intent(inout)::b	
       double precision,intent(inout)::ca,cb,dd 
       double precision,external::namefunc
     end subroutine marq98
     
     subroutine deriva(b,m,v,rl,namefunc)
       integer,intent(in)::m
       double precision,intent(inout)::rl
       double precision,dimension(m),intent(in)::b
       double precision,dimension((m*(m+3)/2)),intent(out)::v    
       double precision,external::namefunc   
     end subroutine deriva
     
     subroutine searpas(vw,step,b,bh,m,delta,fim,namefunc)
       integer,intent(in)::m      
       double precision,dimension(m),intent(in)::b
       double precision,dimension(m),intent(inout)::bh,delta
       double precision,intent(inout)::vw,fim,step  
       double precision,external::namefunc
     end subroutine searpas
     
     subroutine dmfsd(a,n,eps,ier)
       integer,intent(in)::n
       integer,intent(inout)::ier
       double precision,intent(inout)::eps 
       double precision,dimension(n*(n+1)/2),intent(inout)::A      
     end subroutine dmfsd
     
     subroutine valfpa(vw,fi,b,bk,m,delta,namefunc)
       integer,intent(in)::m  
       double precision,intent(in)::vw
       double precision,dimension(m),intent(in)::b,delta  
       double precision,dimension(m),intent(out)::bk 
       double precision,intent(out)::fi
       double precision,external::namefunc  
     end subroutine valfpa
     
     subroutine dmaxt(maxt,delta,m)
       integer,intent(in)::m
       double precision,dimension(m),intent(in)::delta 
       double precision,intent(out)::maxt
     end subroutine dmaxt
  end interface
  
  interface verif2
     subroutine dsinv(A,N,EPS,IER,DET)
       integer,intent(in)::n
       integer,intent(inout)::ier
       double precision,intent(inout)::eps      
       double precision,intent(inout),optional::det     
        double precision,dimension(n*(n+1)/2),intent(inout)::A  
      end subroutine dsinv
      
      subroutine dchole(a,k,nq,idpos)
        integer,intent(in)::k,nq
        integer,intent(inout)::idpos
        double precision,dimension(k*(k+3)/2),intent(inout)::a      
      end subroutine dchole
   end interface
   
 end module type
 
 !----------------------------------------------------------
 !
 !     MODULE PARAMETERS
 !
 ! Derniere mise a jour : 09/02/2011
 !-----------------------------------------------------------
 
 
 module parameters
   double precision,save:: epsa,epsb,epsd
   integer,dimension(:),allocatable,save::vectbHess
   integer,save::maxiter,nBnonHess
 end module parameters
 
 !-------------------------------------------------------------
 !    
 !          MODULE OPTIM avec MARQ98
 !
 !-------------------------------------------------------------
 
 
 module optim_parallele
   
   implicit none
   ! -Interface permettant la verification des type des arguments      
   interface verif1   
      module procedure marq98,deriva,searpas,dmfsd,valfpa
   end interface
   
   interface verif2
      module procedure dsinv,dchole,dmaxt
   end interface
   
   
   
 CONTAINS
   !-------------------------------------------------------------
   !                   MARQ98
   !-------------------------------------------------------------
   
   
   subroutine marq98(b,m,ni,v,rl,ier,istop,ca,cb,dd,namefunc,ng,nvc)
     
     !
     !  fu = matrice des derivees secondes et premieres
     !
     !  istop: raison de l'arret
     !  1: critere d'arret satisfait (prm=ca, vraisblce=cb, derivee=dd)
     !  2: nb max d'iterations atteints
     !  4: Erreur
     
     use parameters
     use ParametresPourParallelisation
     IMPLICIT NONE
     !   variables globales 
     integer,intent(in) :: m,ng,nvc
     integer,intent(inout)::ni,ier,istop
     double precision,dimension(m*(m+3)/2),intent(out)::v
     double precision,intent(out)::rl
     double precision,dimension(m),intent(inout)::b
     double precision,intent(inout)::ca,cb,dd
     !=============================
     !   variables locales 
     
     integer::nql,ii,nfmax,idpos,ncount,id,jd,m1,j,i,ij,i1,j1,ij1
     double precision,dimension(m*(m+3)/2)::fu,v1
	 double precision,dimension((m-nBnonHess)*(m-nBnonHess+1)/2)::fu1
     double precision,dimension(m)::delta,b1,bh
     double precision::da,dm,ga,tr
     double precision::GHG,det,step,eps,vw,fi,maxt, &
          z,rl1,th,ep
     double precision,external::namefunc
     
      !=================================  

	 
	 id=0
     jd=0
     z=0.d0
     th=1.d-5
     eps=1.d-7!1.d-6
     nfmax=m*(m+1)/2    
     ca=epsa+1.d0
     cb=epsb+1.d0
     rl1=-1.d+10    
     ni=0
     istop=0
     da=0.01d0
     dm=5.d0
     nql=1
     m1=m*(m+1)/2
     ep=1.d-20


     Main:Do 
        call deriva(b,m,v,rl,namefunc)

        if(CurrentProcessorID.eq.0) then ! le proc 0 seul
           write(*,*)
        !   write(*,*)'nbnonHess',nbnonHess,vectbhess
           write(*,*)'***     Iteration',ni,'      ***'
           write(*,*)'Log-likelihood=',rl  
           !write(*,*)'B=',(b(j),j=1,m)
           j=1
           write(*,*) 'B='
           do while(j.le.m)
              write(*,*)' ',b(j:min(j+4,m))
              j=j+5
           end do
        end if! proc 0 seul        
        rl1=rl      
        dd = 0.d0     
        fu=0.D0  !2d derivatives delta L/(delta x)2
        do i=1,m
           do j=i,m
              ij=(j-1)*j/2+i
              fu(ij)=v(ij)
           end do
        end do


      !  if(CurrentProcessorID.eq.0) write(*,*)'fu avant dsinv',fu(1:(m*(m+1)/2))

        call dsinv(fu,m,ep,ier,det)  

        if (ier.eq.-1) then
           dd=epsd+1.d0
        else
           GHG=0.d0
           do i=1,m
              do j=1,m
                 if(j.ge.i) then
                    ij=(j-1)*j/2+i
                 else
                    ij=(i-1)*i/2+j
                 end if
                 GHG = GHG + v(m1+i)*fu(ij)*V(m1+j)
              end do
           end do
           dd=GHG/dble(m)
        end if
   !     if(CurrentProcessorID.eq.0) write(*,*)'fu',fu(1:(m*(m+1)/2))

        if(ca.lt.epsa.and.cb.lt.epsb.and.dd.lt.epsd.and.ier.ne.-1) exit main
        !       if(ca.lt.epsa.and.cb.lt.epsb.and.dd.lt.epsd) exit main
        ! sous matrix pour Hessienne
        !     if(ca.lt.epsa.and.cb.lt.epsb.and.ier.eq.1.and.nBnonHess.ne.0) then

 !       if(CurrentProcessorID.eq.0) write(*,*) 'Avant if, ca =',ca,' cb =',cb,' dd =',dd
 !       if(CurrentProcessorID.eq.0) write(*,*) 'ier =',ier,'nBnonHess=',nBnonHess
        if(ca.lt.epsa.and.cb.lt.epsb.and.ier.eq.-1) then
        !   if(CurrentProcessorID.eq.0) write(*,*) 'entree if'

	i1=m-nbnonHess
        fu1=0.D0
        do i=1,i1
           do j=i,i1
              ij=(j-1)*j/2+i
              fu1(ij)=v(ij) ! derivees secondes dL:dx^2 rangees dans fu
           end do
        end do
	
	v1=0.d0
	do i=1,	i1*(i1+1)/2
		v1(i)=v(i) ! derivees secondes rangees dans v1
	end do	
	do i=1,	i1
		v1(i1*(i1+1)/2+i)=v(m1+i) ! derivees premieres rangees dans v1
	end do		
!           i1=0
!           do i=1,i 
!          !    if (vectbhess(i).ne.1.and.abs(b(i)).lt.1.d-2) THEN
!              if (vectbhess(i).ne.1) THEN
!                 i1=i1+1
!                 j1=i1-1
!                 do j=i,m
!              !      if (vectbhess(j).ne.1.and.abs(b(j)).lt.1.d-2) THEN
!                    if (vectbhess(j).ne.1) THEN
!                       j1=j1+1
!                       ij1=(j1-1)*j1/2+i1
!                       ij=(j-1)*j/2+i
!                       fu1(ij1)=v(ij)
!                       !  if(CurrentProcessorID.eq.0) write(*,*)'i,j',i,j,i1,j1
!                    end if
!                 end do
!              end if
!           end do

    !       if(CurrentProcessorID.eq.0) write(*,*)'i1,m-nbnonHess',i1,m-nbnonHess,m,nbnonHess,(m-nBnonHess)*(m-nBnonHess+1)/2

     !      if(CurrentProcessorID.eq.0) write(*,*)'fu1',fu1(1:(i1*(i1+1)/2))
           call dsinv(fu1,i1,ep,ier,det)  
           if (ier.eq.-1) then
              if(CurrentProcessorID.eq.0) write(*,*) 'ca =',ca,' cb =',cb,' dd =',dd
              dd=epsd+1.d0
              if(CurrentProcessorID.eq.0) write(*,*) 'With partial Hessian: ca =',ca,' cb =',cb,' dd =',dd
           else
              GHG=0.d0  

           do i=1,i1
              do j=1,i1
                 if(j.ge.i) then
                    ij=(j-1)*j/2+i
                 else
                    ij=(i-1)*i/2+j
                 end if
                 GHG = GHG + v1((i1)*(i1+1)/2+i)*fu1(ij)*V1((i1)*(i1+1)/2+j)
              end do
           end do
           dd=GHG/dble(i1)

              if (dd.lt.epsd) then 
              if(CurrentProcessorID.eq.0) write(*,*) 'ca =',ca,' cb =',cb,' dd =',epsd
              if(CurrentProcessorID.eq.0) write(*,*) 'With partial Hessian: ca =',ca,' cb =',cb,' dd =',dd
                 istop=5 
                 v=v1
                 goto 110
              end IF
           end IF
        end if
		
		
        tr=0.d0
        do i=1,m
           ii=i*(i+1)/2
           tr=tr+dabs(v(ii))
        end do
        tr=tr/dble(m)

        ncount=0
        ga=0.01d0
400     do i=1,nfmax+m
           fu(i)=v(i)
        end do
        do i=1,m
           ii=i*(i+1)/2
           fu(ii)=v(ii)+da*((1.d0-ga)*dabs(v(ii))+ga*tr)
        end do
        call dchole(fu,m,nql,idpos)

        if (idpos.ne.0) then
           ncount=ncount+1
           if (ncount.le.3.or.ga.ge.1.d0) then
              da=da*dm
           else
              ga=ga*dm
              if (ga.gt.1.d0) ga=1.d0
           endif
           goto 400
        else
           do i=1,m
              delta(i)=fu(nfmax+i)
              !               b1(i)=b(i)+delta(i)
           end do
           b1(:)=b(:)+delta(:)
           rl=namefunc(b1,m,id,z,jd,z)

           if (rl1.lt.rl) then
              if(da.lt.eps) then
                 da=eps
              else
                 da=da/(dm+2.d0)
              endif
              goto 800
           endif
        endif
        !      write(6,*) 'loglikelihood not improved '
        call dmaxt(maxt,delta,m)
        if(maxt.eq.0.D0) then
           vw=th
        else
           call dmaxt(maxt,delta,m)
            vw=th/maxt
         endif
         step=dlog(1.5d0)
         !        !      write(*,*) 'searpas'
         call searpas(vw,step,b,bh,m,delta,fi,namefunc)
         rl=-fi
         if(rl.eq.-1.D9) then
            istop=4
            !write(*,*)'searpas problem'
            goto 110
         end if
         
         do i=1,m
            delta(i)=vw*delta(i)
         end do
         da=(dm-3.d0)*da
         
800      cb=dabs(rl1-rl)
         ca=0.d0
         do i=1,m
            ca=ca+delta(i)*delta(i)
         end do
         if(CurrentProcessorID.eq.0) then    
            write(*,*) 'ca =',ca,' cb =',cb,' dd =',dd
         end if
         do i=1,m
            b(i)=b(i)+delta(i)
         end do
         ni=ni+1
         if (ni.ge.maxiter) then
            istop=2
            if(CurrentProcessorID.eq.0) then      
               write(*,*) 'maximum number of iteration reached'
            end if
            goto 110
         end if

      End do Main
     if(CurrentProcessorID.eq.0) then    
            write(*,*) 'ca =',ca,' cb =',cb,' dd =',dd
         end if
      v=0.D0	 
      v(1:m*(m+1)/2)=fu(1:m*(m+1)/2)
      istop=1
110   continue
      
      return    
    end subroutine marq98
    
    !------------------------------------------------------------
    !                          DERIVA
    !------------------------------------------------------------
    
    subroutine deriva(b,m,v,rl,namefunc)
      use ParametresPourParallelisation
      implicit none
      
      integer,intent(in)::m
      double precision,intent(inout)::rl
      double precision,dimension(m),intent(in)::b
      double precision,dimension((m*(m+3)/2)),intent(out)::v     
      double precision,dimension(m)::fcith
      integer ::i0,m1,ll,i,k,j
      double precision::thn,th,z,vl,temp,thi,thj,rlmax  
      double precision,external::namefunc
      !     partie declaratif pour mpi
      integer::m2
      double precision,dimension(m)::mytab 
      double precision,dimension(m*(m+3)/2)::vv
      integer::test_ok,max_ok
      !=================================    
      !    
      !     v:matrice d'information+score
      !     calcul de la derivee premiere par "central difference"
      !     calcul des derivees secondes par "forward difference"
      !	    v contient d'abord les derivees secondes delta L / delta x1 delta x2
      
      mytab(:)=0.D0 !initialisation du buffer  pour chaque procs
      fcith(:) =0.D0
      
      z=0.d0
      i0=0
      max_ok=0
      test_ok=0
      rl=namefunc(b,m,i0,z,i0,z)   

      ! write(*,*)'dans deriva',rl       
      if(rl.eq.-1.d9) then
         goto 123
      end if
      do i=CurrentProcessorID+1,m,NombreDeProcesseurs
         th=DMAX1(1.d-7, (1.d-4)*DABS(b(i)))
         fcith(i)=namefunc(b,m,i,th,i0,z)       
	if(fcith(i).eq.0)then
		print*,'fcith(i) pas bon',fcith(i),i, th,i0,z
		stop
	end if
         if(fcith(i).eq.-1.d9) then
	 print*,'fcith(i)',fcith(i)
            test_ok=1 
            go to 20013
         end if
      end do
      ! envoi du tableau a tt les procs 

      call MPI_ALLREDUCE(fcith,mytab,m, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      fcith=mytab
      k=0
      m1=m*(m+1)/2
      m2=m1+m
      ll=m1
      vv(:)=0.D0  !initialisation du buffer pour chaque procs 
      v(:)=0.D0
      call MPI_Barrier(MPI_COMM_WORLD,info)
      Main:do i=1,m
         ll=ll+1
         if(mod(ll,NombreDeProcesseurs).eq.CurrentProcessorID) then
            thn=-DMAX1(1.d-7, 1.d-4 * DABS(b(i)))
            temp=namefunc( b,m, i,thn,i0,z)  
       
            if(temp.eq.-1.d9) then
	       print*,'temp1',temp
               test_ok=1
               goto 20013
            end if
            vl=(fcith(i)-temp)/(2.d0*(-thn))
            v(ll)=vl

         elseif(mod(ll,NombreDeProcesseurs).ne.CurrentProcessorID)  then
            continue
         end if
         thi=DMAX1(1.d-7, 1.d-4 * DABS(b(i)))
         do j=1,i
            k=k+1
            if(mod(k,NombreDeProcesseurs).eq.CurrentProcessorID) then
               thj=DMAX1(1.d-7, 1.d-4 * DABS(b(j)))
               temp=namefunc(b,m,i,thi,j,thj)

               if(temp.eq.-1.d9) then
	          print*,'temp2',temp
                  test_ok=1
                  goto 20013
               end if
               v(k)=-(temp-fcith(j)-fcith(i)+rl)/(thi*thj)

            elseif(mod(k,NombreDeProcesseurs).ne.CurrentProcessorID)  then
               continue
            end if
         end do
      end do Main
20013 call MPI_ALLREDUCE(test_ok,max_ok,1, MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
      if(max_ok.ne.0) then
	 print*,'max_ok.ne.0 : ',max_ok,test_ok
         rl=-1.d9
         goto 123
      end if
      call MPI_ALLREDUCE(v,vv,m2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      v=vv         
123   continue 
      return
    end subroutine deriva
    
    
    !------------------------------------------------------------
    !                        SEARPAS
    !------------------------------------------------------------
    
    
    subroutine searpas(vw,step,b,bh,m,delta,fim,namefunc)
      !
      !  MINIMISATION UNIDIMENSIONNELLE
      !
      use ParametresPourParallelisation
      
      implicit none
      
      integer,intent(in)::m      
      double precision,dimension(m),intent(in)::b
      double precision,intent(inout)::vw
      double precision,dimension(m),intent(inout)::bh,delta      
      double precision,intent(inout)::fim,step   
      double precision::vlw,vlw1,vlw2,vlw3,vm,fi1,fi2,fi3    
      integer::i 
      double precision,external::namefunc
      
      vlw1=dlog(vw)
      vlw2=vlw1+step
      call valfpa(vlw1,fi1,b,bh,m,delta,namefunc)
      call valfpa(vlw2,fi2,b,bh,m,delta,namefunc)       
      
      if(fi2.ge.fi1) then
         vlw3=vlw2
         vlw2=vlw1
         fi3=fi2
         fi2=fi1
         step=-step
         
         vlw1=vlw2+step
         call valfpa(vlw1,fi1,b,bh,m,delta,namefunc)   
         if(fi1.gt.fi2) goto 50
      else 
         vlw=vlw1
         vlw1=vlw2
         vlw2=vlw
         fim=fi1
         fi1=fi2
         fi2=fim
      end if
      
      do i=1,40
         vlw3=vlw2
         vlw2=vlw1
         fi3=fi2
         fi2=fi1
         
         vlw1=vlw2+step
         call valfpa(vlw1,fi1,b,bh,m,delta,namefunc)
         if(fi1.gt.fi2) goto 50
         if(fi1.eq.fi2) then
            fim=fi2
            vm=vlw2 
            goto 100
         end if
      end do
      !
      !  PHASE 2 APPROXIMATION PAR QUADRIQUE
      !
50    continue
      !
      !  CALCUL MINIMUM QUADRIQUE
      !
      vm=vlw2-step*(fi1-fi3)/(2.d0*(fi1-2.d0*fi2+fi3))   
      call valfpa(vm,fim,b,bh,m,delta,namefunc)	
      if(fim.le.fi2) goto 100
      vm=vlw2
      fim=fi2
100   continue
      vw=dexp(vm)
      
      return
      
    end subroutine searpas
    
    
    
    !------------------------------------------------------------
    !                         DCHOLE
    !------------------------------------------------------------
    
    subroutine dchole(a,k,nq,idpos)
      
      implicit none
      
      integer,intent(in)::k,nq
      integer,intent(inout)::idpos
      double precision,dimension(k*(k+3)/2),intent(inout)::a
      
      integer::i,ii,i1,i2,i3,m,j,k2,jmk
      integer::ijm,irm,jji,jjj,l,jj,iil,jjl,il
      integer,dimension(k)::is	
      double precision ::term,xn,diag,p
      equivalence (term,xn)
      
      
      !      ss programme de resolution d'un systeme lineaire symetrique
      !
      !       k ordre du systeme /
      !       nq nombre de seconds membres
      !
      !       en sortie les seconds membres sont remplaces par les solutions
      !       correspondantes
      !
      
      i2=0
      ii=0
      idpos=0
      k2=k+nq
      !     calcul des elements de la matrice
      do i=1,k   
         ii=i*(i+1)/2
         !       elements diagonaux
         diag=a(ii)
         i1=ii-i
         if(i-1.ne.0) goto 1
         if(i-1.eq.0) goto 4
1        i2=i-1
         do l=1,i2
            m=i1+l
            p=a(m)
            p=p*p
            if(is(l).lt.0) goto 2
            if(is(l).ge.0) goto 3
2           p=-p
3           diag=diag-p
         end do
         
4        if(diag.lt.0) goto 5
         if(diag.eq.0) goto 50
         if(diag.gt.0) goto 6
5        is(i)=-1
         idpos=idpos+1
         diag=-dsqrt(-diag)
         a(ii)=-diag
         goto 7
6        is(i)=1
         diag=dsqrt(diag)
         a(ii)=diag
         !       elements non diagonaux
7        i3=i+1
         do j=i3,k2
            jj=j*(j-1)/2+i
            jmk=j-k-1
            if(jmk.le.0) goto 9
            if(jmk.gt.0) goto 8
8           jj=jj-jmk*(jmk+1)/2
9           term=a(jj)
            if(i-1.ne.0) goto 10
            if(i-1.eq.0) goto 13 
10          do l=1,i2
               iil=ii-l
               jjl=jj-l
               p=a(iil)*a(jjl)
               il=i-l
               if(is(il).lt.0) goto 11
               if(is(il).ge.0) goto 12
11             p=-p
12             term=term-p
            end do
13          a(jj)=term/diag
         end do
      end do
      
      !       calcul des solutions
      jj=ii-k+1
      do l=1,nq
         jj=jj+k
         i=k-1
14       jji=jj+i
         xn=a(jji)
         if(i-k+1.lt.0) goto 20
         if(i-k+1.ge.0) goto 22
20       j=k-1
21       jjj=jj+j
         ijm=i+1+j*(j+1)/2
         xn=xn-a(jjj)*a(ijm)
         if(j-i-1.le.0) goto 22
         if(j-i-1.gt.0) goto 30
30       j=j-1
         goto 21
22       irm=(i+1)*(i+2)/2
         a(jji)=xn/a(irm)
         if(i.le.0) cycle
         if(i.gt.0) goto 40
40       i=i-1
         go to 14
      end do
50    continue
      return
    end subroutine dchole
    
    
    subroutine dmfsd(a,n,eps,ier)
      !
      !   FACTORISATION DE CHOLESKY D'UNE MATRICE SDP
      !   MATRICE = TRANSPOSEE(T)*T
      !   ENTREE : TABLEAU A CONTENANT LA PARTIE SUPERIEURE STOCKEE COLONNE
      !            PAR COLONNE DE LA MATRICE A FACTORISER
      !   SORTIE : A CONTIENT LA PARTIE SUPPERIEURE DE LA MATRICE triangulaire T
      !
      !   SUBROUTINE APPELEE PAR DSINV
      !
      !   N : DIM. MATRICE
      !   EPS : SEUIL DE TOLERANCE
      !   IER = 0 PAS D'ERREUR
      !   IER = -1 ERREUR
      !   IER = K COMPRIS ENTRE 1 ET N, WARNING, LE CALCUL CONTINUE
      !
      implicit none
      
      integer,intent(in)::n
      integer,intent(inout)::ier
      double precision,intent(in)::eps 
      double precision,dimension(n*(n+1)/2),intent(inout)::A
      double precision :: dpiv,dsum,tol
      integer::i,k,l,kpiv,ind,lend,lanf,lind
      
      !
      !   TEST ON WRONG INPUT PARAMETER N
      !
      dpiv=0.d0
      if (n-1.lt.0) goto 12
      if (n-1.ge.0) ier=0
      !
      !   INITIALIZE DIAGONAL-LOOP
      !
      kpiv=0
      do k=1,n
         kpiv=kpiv+k
         ind=kpiv
         lend=k-1
         !
         !   CALCULATE TOLERANCE
         !
         tol=dabs(eps*sngl(A(kpiv)))
         !
         !   START FACTORIZATION-LOOP OVER K-TH ROW
         !
         do i=k,n
	    dsum=0.d0
     if (lend.lt.0) goto 2
     if (lend.eq.0) goto 4
     if (lend.gt.0) goto 2
     !
     !   START INNER LOOP
     !
2    do l=1,lend
        lanf=kpiv-l
        lind=ind-l
        dsum=dsum+A(lanf)*A(lind)
     end do
     
     !     
     !   END OF INNEF LOOP
     !
     !   TRANSFORM ELEMENT A(IND)
     ! 	
4    dsum=A(ind)-dsum
     if (i-k.ne.0) goto 10
     if (i-k.eq.0) goto 5
     !   TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE
     !	
     
     
5    if (sngl(dsum)-tol.le.0) goto 6
     if (sngl(dsum)-tol.gt.0) goto 9
6    if (dsum.le.0) then
	goto 12 
	end if
     if (dsum.gt.0) goto 7
7    if (ier.le.0) goto 8
     if (ier.gt.0) goto 9
8    ier=k-1
     !
     !   COMPUTE PIVOT ELEMENT
     !
9    dpiv=dsqrt(dsum)
     A(kpiv)=dpiv
     dpiv=1.D0/dpiv
     goto 11
     !
     !   CALCULATE TERMS IN ROW
     !
10   A(ind)=dsum*dpiv
11   ind=ind+i
  end do
end do

!
!   END OF DIAGONAL-LOOP
!
return
12 ier=-1
return

end subroutine dmfsd


!------------------------------------------------------------
!                            DSINV
!------------------------------------------------------------


subroutine dsinv(A,N,EPS,IER,DET)
  
  !
  !     INVERSION D'UNE MATRICE SYMETRIQUE DEFINIE POSITIVE :
  !
  !     MATRICE = TRANSPOSEE(T)*T
  !     INERSE(MATRICE) = INVERSE(T)*INVERSE(TRANSPOSEE(T))
  !
  !     A : TABLEAU CONTENANT LA PARTIE SUPERIEURE DE LA MATRICE A INVERSER
  !         STOCKEE COLONNE PAR COLONNE
  !     DIM. MATRICE A INVERSER = N
  !     DIM. TABLEAU A = N*(N+1)/2
  !
  !     EPS : SEUIL DE TOLERANCE AU-DESSOUS DUQUEL UN PIVOT EST CONSIDERE
  !           COMME NUL
  !
  !     IER : CODE D'ERREUR
  !         IER=0 PAS D'ERREUR
  !         IER=-1 ERREUR SUR LA DIM.N OU MATRICE PAS DEFINIE POSITIVE
  !         IER=1 PERTE DE SIGNIFICANCE, LE CALCUL CONTINUE
  !=========================================================================
  
  use ParametresPourParallelisation    
  implicit none
  
  integer,intent(in)::n
  integer,intent(inout)::ier
  double precision,intent(inout)::eps      
  double precision,intent(inout),optional::det     
  double precision,dimension(n*(n+1)/2),intent(inout)::A     
  double precision::din,work
  integer::ind,ipiv,i,j,k,l,min,kend,lhor,lver,lanf
  
  !=========================================================================
  !     FACTORIZE GIVEN MATRIX BY MEANS OF SUBROUTINE DMFSD
  !     A=TRANSPOSE(T) * T
  !=========================================================================
  
  call dmfsd(A,n,eps,ier)
  
  det=0.d0
  
  if (ier.lt.0) goto 9
  if (ier.ge.0) det=0.d0
  
  !========================================================================
  !     INVERT UPPER TRIANGULAR MATRIX T
  !     PREPARE INVERSION-LOOP
  !
  !
  ! calcul du log du determinant 
  !=======================================================================
  
  do i=1,n
     det=det+dlog(A(i*(i+1)/2))
  end do
  det=2*det
  ipiv=n*(n+1)/2
  ind=ipiv
  !
  !     INITIALIZE INVERSION-LOOP
  !
  do i=1,n
     din=1.d0/A(ipiv)
     A(ipiv)=din
     min=n
     kend=i-1
     lanf=n-kend
     if (kend.le.0) goto 5
     if (kend.gt.0) j=ind
     !
     !     INITIALIZE ROW-LOOP
     !
     do k=1,kend
        work=0.d0
        min=min-1
        lhor=ipiv
        lver=j
        !
        !     START INNER LOOP
        !
        do l=lanf,min 
           lver=lver+1
           lhor=lhor+l
           work=work+A(lver)*A(lhor)
        end do
        !
        !     END OF INNER LOOP
        !
        A(j)=-work*din
        j=j-min
     end do
     
     !
     !     END OF ROW-LOOP
     !
5    ipiv=ipiv-min 
     ind=ind-1
  end do
  
  !
  !     END OF INVERSION-LOOP
  !
  !     CALCULATE INVERSE(A) BY MEANS OF INVERSE(T)
  !     INVERSE(A) = INVERSE(T) * TRANSPOSE(INVERSE(T))
  !     INITIALIZE MULTIPLICATION-LOOP
  !
  do i=1,n
     ipiv=ipiv+i
     j=ipiv
     !
     !     INITIALIZE ROW-LOOP
     !
     do k=i,n
        work=0.d0
        lhor=j
        !
        !     START INNER LOOP
        !
        do l=k,n
           lver=lhor+k-i
           work=work+A(lhor)*A(lver)
           lhor=lhor+l
        end do
        !
        !     END OF INNER LOOP
        !       
        A(j)=work
        j=j+k
     end do
  end do
  
  !
  !     END OF ROW-AND MULTIPLICATION-LOOP
  !
9 return
end subroutine dsinv

!------------------------------------------------------------
!                          VALFPA
!------------------------------------------------------------

subroutine valfpa(vw,fi,b,bk,m,delta,namefunc)
  use ParametresPourParallelisation
  implicit none
  
  integer,intent(in)::m  
  double precision,dimension(m),intent(in)::b,delta  
  double precision,dimension(m),intent(out)::bk 
  double precision,intent(out)::fi 
  double precision::vw,z	
  integer::i0,i
  double precision,external::namefunc
  
  z=0.d0
  i0=1
  do i=1,m
     bk(i)=b(i)+dexp(vw)*delta(i)
  end do
  fi=-namefunc(bk,m,i0,z,i0,z)
  
  return
  
end subroutine valfpa

!------------------------------------------------------------
!                            MAXT
!------------------------------------------------------------


subroutine dmaxt(maxt,delta,m)
  use ParametresPourParallelisation     
  implicit none
  
  integer,intent(in)::m
  double precision,dimension(m),intent(in)::delta
  double precision,intent(out)::maxt
  integer::i 
  
  maxt=Dabs(delta(1))
  do i=2,m
     if(Dabs(delta(i)).gt.maxt)then
        maxt=Dabs(delta(i))
     end if
  end do
  
  return
end subroutine dmaxt

end module optim_parallele
