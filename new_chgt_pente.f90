!===========================================================
!
!      Program for estimating latent process mixed models for
!      multivariate longitudinal outcomes and times-to-events with multiple failure types 
!         using latent classes to link the outcomes together
!      times-to-events may be left-truncated & right-censored
!      Multivariate longitudinal outcomes are linked to the latent 
!      process through different link functions : parameterized nonlinear
!        continuous transformations or threshold models 
!
!      Cecile PROUST-LIMA-Florian
!
!      INSERM U897, BORDEAUX, FRANCE
!
!
!  The procedure is detailed in README.txt
!
!  The program depends on two modules:
!   Optim_routines.f90 and Gaussian_routines.f90
!  The algorithm of optimisation is a modified Marquardt
!  algorithm with three combined convergence criteria
!  (on parameter stability, log-likelihood stability and
!   first (g) and second (H) derivatives with g'H-1g<eps)
!
! This programs works with 0 or 3 mixture components only
! and no time-dependent covariates.
!==========================================================

!==========================================================
!      Hetmixsurv_v2 
!                                         8/01/2013
!==========================================================





!==========================================================
!      MODULES  
!                                         8/01/2013
!==========================================================




! {{{

!********************** interface peut être rajoutee : interface.f90
!

      module commun

      use ParametresPourParallelisation

      implicit none
      integer ,save::ns,ng,nv,nq,np,npp,idiag,ntr,iAR,nvdep,nrisq,npg,ntrspl,maxmestot,mesure,nprm,nxevt,tbeta,interv_censoring
      integer ,save::nea,nvarglobssG,nvarglobG,nvarglob,evt,nvarprob,nprmdep,nvdepssG,nvdepG,ncontdep
      integer ,save::ncont,nef,npm_ne,nvarxevt,nvc,nAR,idtrunc,nvdepsurv,nvardepsurv,paramsurv,vcclass,nlg
      integer ,save::nva01,nva02,nva12,troncature,lenb !pour smooth hazard
      integer,dimension(:),allocatable,save::typrisq,nprisq,risqcom,nrisqspec,nxevtspec,nz,idzi
      integer, save ::npmtot,quecont,npmtest,queord,spline
      integer,save :: nztr,idzitr,nevt,sumidxevt,evtact,nparestsurv,nparestvexpsurv
      integer,dimension(:),allocatable,save::numero,ind_survint,Devt
      integer,dimension(:,:),allocatable,save ::nmes
      integer,dimension(:),allocatable,save::ide,Navant_i
      integer,dimension(:),allocatable,save ::idglob,idg,idea,idtest,idprob,idevt,idxevt,idvdep,idmort,iddem
      double precision, dimension(:),allocatable,save::zi01,zi02,zi12
      double precision,save::min_save,max_save,stept,trn
      integer,dimension(:),allocatable,save ::idord
      real,dimension(:),allocatable,save ::Y,t
      real,dimension(:,:),allocatable,save ::Xdep
      real,dimension(:,:),allocatable,save ::X
      real,dimension(:),allocatable,save::Tsurv0,Tsurv,Tsurvint

      real,dimension(:),allocatable,save::t0,t1,t2,t3
	  real,dimension(:),allocatable,save::Tevt
	  real,dimension(:),allocatable,save::tt0
      real,dimension(:,:),allocatable,save::sumvect,tfinmoy
      integer,dimension(:),allocatable,save::c,idm,idd
      real,dimension(:),allocatable,save::pipre
      real,save::sumvectpi
      integer,save::contevt,contevt2,bparc2,numact2
      double precision,dimension(:,:),allocatable,save ::pig
      double precision,dimension(:),allocatable,save ::range,rangereal
      double precision, dimension(:),allocatable,save::bne
      double precision, dimension(:),allocatable,save::Tpred
      double precision,dimension(:),allocatable,save ::mm,mm1,mm2,im,im1,im2
      double precision,dimension(:,:),allocatable,save ::zitr,Y1all
      double precision,dimension(:),allocatable,save::Tmm,Tmm1,&
           Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,Tmm01,Tmm02,Tmm03,Tim0 &
           ,Tim01,Tim02,Tim03,Tmmt,Tmmt1,Tmmt2,Tmmt3,Timt,Timt1,&
           Timt2,Timt3
      double precision,dimension(:),allocatable,save::Tmm_est,Tmm1_est,&
                Tmm2_est,Tmm3_est,Tim_est,Tim1_est,Tim2_est,Tim3_est &
                ,Tmm_est2,Tmm1_est2,&
                Tmm2_est2,Tmm3_est2,Tim_est2,Tim1_est2,Tim2_est2,Tim3_est2

      end module commun


      module donnees_indiv
      implicit none
      double precision::jac
      double precision,dimension(:,:),allocatable::Ut1
      double precision, dimension(:),allocatable::mu,Y1,mu2
      double precision,dimension(:,:),allocatable::Z1
      integer::numpat,gcourant
      double precision,dimension(:,:),allocatable:: Sigmaq,Vw
      integer,parameter ::nf=1
      double precision,dimension(:),allocatable::b1

      end module donnees_indiv

! }}}
 




!----------------------------------------------------------
!                   MAIN PROGRAM
!----------------------------------------------------------

! {{{

      program hetmixbeta


      use commun
      use optim_parallele
      use donnees_indiv
      use parameters
      use lois_normales
      implicit none

      integer ::ier,nobs,id,istop,ni,npm,m,q,i,j,k,l,g,jj,kk,nrisq_curr,igl
      integer ::nsim,itemp,ndata,ind,ntrtot,nmestot,indice,m1,ONLYWeibull,nsim2
      integer ::kkk,contsortie,contfins,sex
      integer,dimension(:),allocatable ::classe,classe1
      integer,dimension(:,:),allocatable ::classement,classement1
      integer,dimension(:,:),allocatable ::claMCpl,claFCpl,claMCmo,claFCmo
      real ::nobsmoy,time,time_deb,time_fin
      real ::heure,minute,seconde
      integer,dimension(:),allocatable::nb_evt
!      integer,dimension(:),allocatable::valeur
      double precision ::binf,bsup,aic,bic,EPS,vrais,betai,icl
      double precision ::temp,ytemp,ca,cb,dd,pas,inf,sup,aa,bb,aa1,bb1,cc1,dd1
      double precision,dimension(:),allocatable ::mvcg,mvc,Bh1,se
      integer,dimension(:),allocatable :: pitest,pipost
      double precision,dimension(:),allocatable::valerr,vtest
      double precision,dimension(:),allocatable ::Vopt
      double precision,dimension(:,:),allocatable ::Ut,U,VC,V1,CORR
      double precision,dimension(:,:),allocatable ::PPI,ppitest,P_AI,pi,probapost
      double precision,dimension(:),allocatable::Ytemp2,bsex
      integer,dimension(:),allocatable::nmestotcur


!      double precision,dimension(:,:),allocatable::ymarg
      double precision,dimension(:,:),allocatable::mui1,ymarg1,ycond1,muicond1




      character(len=50) ::nomdata,nomoutput,nomoutput2
      double precision,dimension(:),allocatable::B,brisq_est,time_est,time_est01,time_est02,time_est12,time_est2,brisq,brisqtot
      double precision,dimension(:,:,:),allocatable::som,pond,stderr,temp_std,som_obs,IC_inf &
           ,IC_sup,pond2,som_ui,som_ui_obs
      double precision,dimension(:,:),allocatable::som_marg,som_marg_ui,som_obs_yi,IC_sup_yi,IC_inf_yi,stderr_yi,temp_std_yi
      double precision,dimension(:),allocatable::cl_tps
      integer,dimension(:,:),allocatable::ns_cl_tps
      double precision,dimension(:,:),allocatable::time_cl_tps,incidcum_est,risq_est,risqcum_est, &
           risqspec,risqcumspec
      integer::nb_cl_tps
      double precision,dimension(:),allocatable::splaa,splaa2
      double precision,dimension(:),allocatable::zitr2
!      double precision,dimension(:,:),allocatable::pred_ind
!      double precision,dimension(:,:,:),allocatable::Y_trans,mui
!      double precision,dimension(:,:),allocatable::alpha_i,u_i
      double precision,dimension(:),allocatable::test,transf
    double precision, external::funcpa
	!  integer, parameter :: QR_K = selected_real_kind (16)
     !  real (kind=QR_K) :: funcpa
      double precision,dimension(:,:),allocatable::gaussGL


      double precision,dimension(:),allocatable::valinit,uscoreYT
      double precision,dimension(:),allocatable::b3,b4,varV
      double precision,dimension(:,:),allocatable::varM
      double precision,dimension(:),allocatable::b3g,b4g,varVg,var_uYT
      double precision,dimension(:,:),allocatable::varMg,var_uscoreYT,var_YT
      double precision::Zscore,det
      integer:: cursor
      real::elapsed_time_op,time_op1,time_op2,time_max
      double precision::res_inst,res_cum,val_func,res_tot,uscoreTT,var_uscoreTT,test_scoreTT,test_scoreYT
      integer :: nidevt,nprob,nhazard,ncov,nlatent,ncontr,ntransf,nvarcov,num_frailty,scoretest


!--------- Lecture du fichier parametres ---------------------

      eps=1.d-10
      ntr=4
      trn=0.01d0
      

      contsortie=0
      contfins=0


      !==================================
      !      Initialisation de MPI:
      !==================================
      call MPI_INIT(info)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, NombreDeProcesseurs, info)
      call MPI_COMM_RANK(MPI_COMM_WORLD, CurrentProcessorID, info)
      !==================================
      










      call cpu_time(time_deb)


      open(2,file='step_verif.inf',status='old') ! changer entre les .inf

      read(2,*)     ! nom du fichier de donnees
      read(2,*)nomdata
      if(CurrentProcessorID.eq.0)then
         write(*,*)
         write(*,*)'Datafile: ',nomdata
      end IF
      read(2,*)     ! nom du fichier de resultats
      read(2,*)nomoutput
      if(CurrentProcessorID.eq.0)then
         write(*,*)'Output file: ',nomoutput
      end IF
      read(2,*)     ! nom du fichier pour predictions
      read(2,*)nomoutput2
      read(2,*)
      read(2,*)mesure,nprm
      if(CurrentProcessorID.eq.0)then
         write(*,*)'max number of measures:',mesure,'max number of parms',nprm
      end IF
      read(2,*)     ! nombre de sujets
      read(2,*)ns
      if(CurrentProcessorID.eq.0) write(*,*)'Number of subjects:',ns

      read(2,*)     ! nombre de tests
      read(2,*)nq
      if(CurrentProcessorID.eq.0) write(*,*)'Number of dependent variables:',nq

      allocate(idord(nq),range(2*nq),rangereal(2*nq))
      idord=0
      range=0
      rangereal=0


      read(2,*)     ! indicateur de nature de test : 0=continu et k=kmodalites
      read(2,*)(idord(k),k=1,nq)
      read(2,*)     !  Indicator if the outcome needs to be beta-transformed
      read(2,*)tbeta
      read(2,*)     !  Indicator if the estimation accounts for interval censoring
      read(2,*)interv_censoring
      read(2,*)     ! range des tests
      read(2,*)(range(k),k=1,2*nq)
      read(2,*)     ! range des tests
      read(2,*)(rangereal(k),k=1,2*nq)
      do k=1,nq
         if(idord(k).gt.1) then
            if (range(2+(k-1)*2)-range(1+(k-1)*2)+1.ne.idord(k)) then
               if(CurrentProcessorID.eq.0)then
                  write(*,*)'probleme range pour ordinaux',     &
                       ' : doit concorder avec nombre de modalites'
               end IF
               
               stop
            end if
         end if
      end do



      read(2,*)     ! degre du polynome du temps
      read(2,*)np
	npp=(2*np+1)
      if(CurrentProcessorID.eq.0)then
         write(*,*)'Degree of the polynomious time function: ',np
      end if

      allocate(idea(np+1))
      idea=0

      read(2,*)
      read(2,*)(idea(k),k=1,np+1)

      l=0
      Do k=1,np+1
         if (idea(k).eq.1) then
            l=l+1
	    if(k.gt.1) l=l+1
         end if
      end do
      nea=l

      read(2,*)     ! nombre de composantes
      read(2,*)ng

      allocate(pig(ns,ng)) ! P(c_i=g| X_i)
      pig=0.d0
      if(CurrentProcessorID.eq.0)then
         write(*,*)'Number of random-effects: ',nea
         write(*,*)'Number of latent classes: ',ng
      end if

      if (ng.lt.1) then
         
         if(CurrentProcessorID.eq.0)then
            write(*,*)'the number of components must be greater or equal to 1'
         end IF
         
         stop
      end if

      read(2,*)
      read(2,*)nevt
     ! nevt=3 !valeur fixee
      if (nevt.lt.0) then
         if(CurrentProcessorID.eq.0)then
            write(*,*)'number of events not negative'
         end IF
         stop
      end if

      read(2,*)

	  allocate(typrisq(nevt),risqcom(nevt),nrisqspec(nevt),nprisq(nevt),nxevtspec(nevt),nz(nevt),idzi(nevt))
      if (nevt.ne.0) then ! quelle fonction de base ? risqcom, noeuds et ecart        
         read(2,*) (typrisq(k),k=1,nevt),(risqcom(k),k=1,nevt),idtrunc,(nz(k),k=1,nevt),(idzi(k),k=1,nevt),paramsurv
         if (idtrunc.eq.1) then
            if(CurrentProcessorID.eq.0)then
            write(*,*)'accounts for the left-truncation given in the dataset'
            end if
         end if
         troncature=idtrunc

         nprisq=0
         ONLYWeibull=1
         do k=1,nevt
            if (typrisq(k).eq.1) then
               nprisq(k)=nz(k)-1
               ONLYWeibull=0
            else if (typrisq(k).eq.2) then
               nprisq(k)=2
            else
				write(*,*) "Risques de base autres que Weibull et step functions non geres"
            stop
         endif
         end do
      else 
         read(2,*)
      end if
	  !condition weibull, pour l'instant la partie smoothhazard ne gere que ce cas
      !do k=1,nevt

      !end do

! nombre de splines s'il y en a et indicateur:
! 1 : noeuds equi/ 0 : noeuds quantiles/2 noeuds en entree
      read(2,*)
      read(2,*)nztr, idzitr

      ntrspl=nztr+1+1

      allocate(zitr(-1:nztr+2,nq))

      zitr=0.d0

      do q=1,nq
         if (idord(q).eq.1) then
            if (idzitr.eq.2) then
               read(2,*)(zitr(j,q),j=2,nztr-1)
            end if
            zitr(1,q)=range(1+(q-1)*2)
            zitr(nztr,q)=range(2+(q-1)*2)
         else
            zitr(1,q)=rangereal(1+(q-1)*2)
            zitr(nztr,q)=rangereal(2+(q-1)*2)
         end if
      end do

      if (idzitr.ne.0.and.idzitr.ne.1.and.idzitr.ne.2) then

         if(CurrentProcessorID.eq.0) write(*,*)'0 : noeuds aux quantiles ou 1 : noeuds equidistants'
         if(CurrentProcessorID.eq.0) write(*,*)'ou 2 : noeuds interieurs entres manuellement'
         stop
      end if


      npmtest=0
      quecont=1
      queord=1
      spline=0

      do k=1,nq
         if (idord(k).gt.1) then
            quecont=0
            npmtest=npmtest+idord(k)-1
         else if (idord(k).eq.0) then
            queord=0
            npmtest=npmtest+ntr
         else if (idord(k).eq.-1) then
            quecont=0
            npmtest=npmtest+ntr
         else if (idord(k).eq.1) then
            queord=0
            npmtest=npmtest+ntrspl
            spline=1
         else
            if(CurrentProcessorID.eq.0) write(*,*)'problem: idord pas dans le range attendu'
            if(CurrentProcessorID.eq.0) write(*,*)'idord',idord(k),'k=',k
         end if
      end do

!C nombre de parms pour les fonctions de risq
      read(2,*)           ! fonction dep du temps dans la survie ?
      read(2,*)nvardepsurv,nvdepsurv
      if (nvardepsurv.gt.1) then
         if(CurrentProcessorID.eq.0) write(*,*)'possible values for nvdepsurv for the moment:0 and 1'
         stop
      end if

!C nombre de variables dependantes du temps dans le fichier de donnees + effet sur le temps (1=effet commun; 2=effet class-specific)

      read(2,*)
      read(2,*)nvdep

      allocate(idvdep(nvdep))
      read(2,*)
      read(2,*)(idvdep(j),j=1,nvdep)
!      if(CurrentProcessorID.eq.0) write(*,*)'nvdep',nvdep
!      if(CurrentProcessorID.eq.0) write(*,*)'idvdep',(idvdep(j),j=1,nvdep)
!C nombre de variables explicatives dans le fichier de donnees
      read(2,*)
      read(2,*)nv
      if(CurrentProcessorID.eq.0) write(*,*)'Number of independent covariates: ',nv

      allocate(idglob(nv),idg(nv),idtest(nv),idprob(nv),idevt(nv),idmort(nv),iddem(nv))
      if (nevt>1) then
         allocate(idxevt(nv*nevt))
      Else
         allocate(idxevt(nv))
      end IF
   


!C indicatrice d'effet global de la variable explicative
!C (i.e. sur la fonction cognitive globale)
      read(2,*)
      read(2,*)(idglob(k),k=1,nv)

      if(CurrentProcessorID.eq.0) write(*,*)'idglob',(idglob(k),k=1,nv)

      read(2,*)
      read(2,*)(idg(k),k=1,nv)
      if(CurrentProcessorID.eq.0) write(*,*)'idg',(idg(k),k=1,nv)

      itemp=0.d0
      do k=1,nv
         itemp=itemp+idg(k)
      end do
      if (ng.eq.1.and.itemp.gt.0) then
         if(CurrentProcessorID.eq.0) write(*,*) 'contradiction : no mixture but mixture effects specified'
         stop
      end if


      do k=1,nv
         if (idg(k).eq.1.and.idglob(k).eq.0) then
            if(CurrentProcessorID.eq.0) write(*,*)'Mixture effect only if global effect'
            if(CurrentProcessorID.eq.0) write(*,*)'covariate',k,'problem of definition'
            stop
         end if
      end do

      read(2,*)
      read(2,*)(idprob(k),k=1,nv)
      if(CurrentProcessorID.eq.0) write(*,*)'idprob',(idprob(k),k=1,nv)

!C  id evt : indicateur de la variable : evenement binaire
      read(2,*)
!cas un seul evt
      read(2,*)(idmort(k),k=1,nv)
      read(2,*)
!event mort et demence : deux indicateurs
      read(2,*)(iddem(k),k=1,nv)
!C verif qu'il y en a qu'un :
      evt=sum(idmort(:))+sum(iddem(:)) !variable de test
!idm puis idd puis variables exp : idm et idd forcement presents

!ne permettais pas la gestion de plusieurs evt : impossibilite de demence puis mort d'ou le commentaire

      if (evt.ne.2) then
         if(CurrentProcessorID.eq.0) write(*,*) 'Should include 2 outcome at max : death and dementia'
         stop

!on peut avoir les 2 evt
      ! else if(nevt*evt.ne.nevt) then  
!          if(CurrentProcessorID.eq.0) write(*,*)'inconsistency between the number of events and the indicator of event'
!          stop
      end if

!C variables explicatives dans la modelisation de l'evenement
      read(2,*)
      if (nevt.ge.1) then
         read(2,*)(idxevt(k),k=1,nv*nevt)
      Else
         read(2,*)
     end IF
     sumidxevt=0
      if(CurrentProcessorID.eq.0) write(*,*)'idxevt',idxevt
      do k=1,size(idxevt)
         if (idxevt(k).ge.1) then
            sumidxevt=sumidxevt+1
         endif
      enddo

      !fonctionnel pour allocation du vecteur vexpiact
      
      nrisq=0
      if (nevt.ne.0) then
         nrisqspec=0
         do k=1,nevt
            if (risqcom(k).eq.1) then
               nrisqspec(k)=nprisq(k)
            end if
            if (risqcom(k).eq.2) then
               nrisqspec(k)=nprisq(k)+ng-1
            end if
            if (risqcom(k).eq.0) then
               nrisqspec(k)=nprisq(k)*ng
            end if
            nrisq=nrisq+nrisqspec(k)
         end do
      end if
!C definir des contrastes plutôt
      read(2,*)
      read(2,*)(idtest(k),k=1,nv)
      if(CurrentProcessorID.eq.0) write(*,*)'idtest',(idtest(k),k=1,nv)

!c----------calculs du nombre de parms fixes ----------------
!C parms fixes pour les vexp dans les probas
      nvarprob=1
      Do k=1,nv
         If (idprob(k).eq.1) then
            nvarprob=nvarprob+1
         end if
      end do

!C nombre de parametres pour les variables dependantes du temps
      nprmdep=0
      nvdepssG=0
      nvdepG=0
      ncontdep=0
      do k=1,nvdep
         if (idvdep(k).eq.1) then
            nprmdep=nprmdep+1
            nvdepssG=nvdepssG+1
         elseif (idvdep(k).eq.2) then
            nprmdep=nprmdep+ng
            nvdepG=nvdepG+1
         elseif (idvdep(k).eq.3) then
            nprmdep=nprmdep+1
            nvdepssG=nvdepssG+1
            ncontdep=ncontdep+1
         elseif (idvdep(k).eq.4) then
            nprmdep=nprmdep+ng
            nvdepG=nvdepG+1
            ncontdep=ncontdep+1
         end if
      end do
!      if(CurrentProcessorID.eq.0) write(*,*)'nprmdep=',nprmdep
!C      if ((nvdepG+nvdepssG.ne.nprmdep) then
!C         if(CurrentProcessorID.eq.0) write(*,*) 'problem nombre parms nprmdep'
!         if(CurrentProcessorID.eq.0) write(*,*)'nvdepG',nvdepG,'nvdepssG',nvdepssG ,'ncontdep*(nq-1)',ncontdep*(nq-1)
!C         stop
!C      end if

!C parms fixes pour les vexp dans la modelisation de l'evenement. Attention
      nvarxevt=nvdepsurv
      nparestsurv=0
      nxevt=0  
	  
	  if (nevt.ne.0) then
      if (evt.ne.0) then
      nxevtspec=0
         do kk=1,NEVT
            Do k=1,nv
               IF (idxevt((kk-1)*nv+k).eq.1) THEN 
                  nvarxevt=nvarxevt+1
                  nxevtspec(kk)=nxevtspec(kk)+1
                  nparestvexpsurv=nparestvexpsurv+1
               end if
               If (idxevt((kk-1)*nv+k).eq.2) then
                  nvarxevt=nvarxevt+ng
                  nxevtspec(kk)=nxevtspec(kk)+1
                  nparestvexpsurv=nparestvexpsurv+1
               end if
            end do
            nxevt=nxevt+nxevtspec(kk)
         end do

!recuperation nva01
         nva01=nxevtspec(1)
         nva02=nxevtspec(2)
         nva12=nxevtspec(3)
!pour modele multi etat et smooth
        
         
!!!! ATTENTION: prevoir le cas ou des effets sont les mêmes sur plusieurs outcomes
      Else
         nvarxevt=0
      end if
end if

!C parms effets fixes dans partie latente : nvarglob
      nvarglob=0
      nvarglobssG=0
      nvarglobG=0

      Do k=1,nv
         If (idglob(k).eq.1.and.idg(k).eq.0) then
            nvarglob=nvarglob+1
            nvarglobssG=nvarglobssG+1
         else if (idglob(k).eq.2.and.idg(k).eq.0) then
            nvarglob=nvarglob+2*(np)+1
            nvarglobssG=nvarglobssG+2*(np)+1
         else if (idglob(k).eq.1.and.idg(k).eq.1) then
	    print*,'Does not account for class-specific changepoint + class-specific effects '
	    stop
            nvarglob=nvarglob+ng
            nvarglobG=nvarglobG+1
         else if (idglob(k).eq.2.and.idg(k).eq.1) then
	    print*,'Does not account for class-specific changepoint + class-specific effects '
	    stop
            nvarglob=nvarglob+(np+1)*ng
            nvarglobG=nvarglobG+np+1
         end if
      end do

      ncont=0.d0
      l=0
      do k=1,nv
         If (idtest(k).eq.1) then
            l=l+1
         else if (idtest(k).eq.2) then
            l=l+1+np
         end if
      end do
      ncont=l

      nef=ng*(2*np+1)+nvarprob*(ng-1)+nrisq+nvarxevt+ng+nprmdep+nvarglob+ ncont*(nq-1)  &
           +ncontdep*(nq-1)+npmtest

      if(CurrentProcessorID.eq.0) write(*,*)'nprmdep',nprmdep
      if(CurrentProcessorID.eq.0) write(*,*)'nvarprob',nvarprob
      if(CurrentProcessorID.eq.0) write(*,*)'nvarglob',nvarglob
      if(CurrentProcessorID.eq.0) write(*,*)'ncont',ncont
      if(CurrentProcessorID.eq.0) write(*,*)'ncontdep',ncontdep
      if(CurrentProcessorID.eq.0) write(*,*)'nef',nef
      if(CurrentProcessorID.eq.0) write(*,*)'npmtest',npmtest


      allocate(valinit(nprm),ide(nprm),bne(nprm))

      valinit=0.d0
      ide=0
      bne=0.d0
!C entrer les effets des variables sur l'evolution :
!C - d'abord les varexp dans le logit
!C - variables dans l'evt
!C - puis les variables du temps
!C - puis les variables explicatives par ordre d'entree (si mixture entrer les parms par mixture sur chaque variable (attention aux interactions avec le temps)
!C   dans le fichier de donnees. Pour les effets avec mixture :
!C   mettre par mixture puis par variable
     read(2,*)
  
      read(2,*)(valinit(j),j=1,nvarprob*(ng-1))

     
      read(2,*)
      read(2,*)(ide(j),j=1,nvarprob*(ng-1))
      if(CurrentProcessorID.eq.0) write(*,*)'b(prob)',(valinit(j),j=1,nvarprob*(ng-1))
	  
	  nprob=0
	  do j=1,nvarprob*(ng-1)
		if (ide(j).eq.1) then
			nprob=nprob+1
		end if
	  end do
	  
      read(2,*)
      read(2,*)(valinit(j),j=nvarprob*(ng-1)+1,nvarprob*(ng-1)+nrisq)
      read(2,*)
      read(2,*)(ide(j),j=nvarprob*(ng-1)+1,nvarprob*(ng-1)+nrisq)
      if(CurrentProcessorID.eq.0) write(*,*)'b(risq)',(valinit(j),j=nvarprob*(ng-1)+1,(nvarprob)*(ng-1)+nrisq)
	  
	  	  nhazard=0
	  do j=nvarprob*(ng-1)+1,nvarprob*(ng-1)+nrisq
		if (ide(j).eq.1) then
			nhazard=nhazard+1
		end if
	  end do

      read(2,*)
      read(2,*)(valinit(j),j=(nvarprob)*(ng-1)+nrisq+1,(nvarprob)*(ng-1)+nrisq+nvarxevt)
      read(2,*)
      read(2,*)(ide(j),j=(nvarprob)*(ng-1)+nrisq+1,(     &
 nvarprob)*(ng-1)+nrisq +nvarxevt)
      if(CurrentProcessorID.eq.0) write(*,*)'b(evt)',(valinit(j),j=(nvarprob)*(ng-1)+nrisq+1, &
 (nvarprob)*(ng-1)+nrisq+nvarxevt)
 
	  	  ncov=0
	  do j=(nvarprob)*(ng-1)+nrisq+1,(nvarprob)*(ng-1)+nrisq+nvarxevt
		if (ide(j).eq.1) then
			ncov=ncov+1
		end if
	  end do
      read(2,*)
      read(2,*)(valinit(j),j=(nvarprob)*(ng-1)+nrisq+nvarxevt+1,     & 
(nvarprob)*(ng-1)+nrisq+nvarxevt+ng)

      if(CurrentProcessorID.eq.0) write(*,*)'class-spec changepoints',(valinit(j),j=(nvarprob)*(ng-1)+nrisq & 
    +nvarxevt+1,(nvarprob)*(ng-1)+nrisq+nvarxevt+ng)
      read(2,*)
      read(2,*)(ide(j),j=(nvarprob)*(ng-1)+nrisq+nvarxevt+1,    &
  (nvarprob)*(ng-1)+nrisq+nvarxevt +ng)
      if(CurrentProcessorID.eq.0) write(*,*)'class-spec changepoints',(ide(j),j=(nvarprob)*(ng-1)+nrisq & 
    +nvarxevt+1,(nvarprob)*(ng-1)+nrisq+nvarxevt+ng)

      read(2,*)
      read(2,*)(valinit(j),j=(nvarprob)*(ng-1)+nrisq+nvarxevt+ng+1,     & 
(nvarprob)*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob)

      if(CurrentProcessorID.eq.0) write(*,*)'b(latent)',(valinit(j),j=(nvarprob)*(ng-1)+nrisq & 
    +nvarxevt+1,(nvarprob)*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1) & 
    +nprmdep +nvarglob)
      read(2,*)
      read(2,*)(ide(j),j=(nvarprob)*(ng-1)+nrisq+nvarxevt+ng+1,    &
  (nvarprob)*(ng-1)+nrisq+nvarxevt +ng+ng*(2*np+1)+nprmdep+nvarglob)
      if(CurrentProcessorID.eq.0) write(*,*)'b(latent)',(ide(j),j=(nvarprob)*(ng-1)+nrisq & 
    +nvarxevt+ng+1,(nvarprob)*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1) & 
    +nprmdep +nvarglob)


		  nlatent=0
		  do j=(nvarprob)*(ng-1)+nrisq+nvarxevt+ng+1,    &
	  (nvarprob)*(ng-1)+nrisq+nvarxevt +ng+ng*(2*np+1)+nprmdep+nvarglob
			if (ide(j).eq.1) then
				nlatent=nlatent+1
			end if
		  end do
	  
      IF (valinit((nvarprob)*(ng-1)+nrisq+nvarxevt+ng+1).ne.0.and.ide((nvarprob)*(ng-1)+nrisq &
     +nvarxevt+ng+1).ne.0) then
         if(CurrentProcessorID.eq.0) write(*,*)'la premiere intercept doit','  etre a 0 et non estimee'
!C         stop
      ENDIF

      read(2,*)
      read(2,*)(valinit(j),j=(nvarprob)*(ng-1)+nrisq+nvarxevt  &
    +ng+ng*(2*np+1)+nprmdep+nvarglob+1,(nvarprob)*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob & 
    +ncontdep*(nq-1)+ncont*(nq-1))
      read(2,*)
      read(2,*)(ide(j),j=(nvarprob)*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+1, & 
    (nvarprob)*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &  
    +ncontdep*(nq-1)+ncont*(nq-1))
      if(CurrentProcessorID.eq.0) write(*,*)'b(tests)',(valinit(j),j=(nvarprob)*(ng-1)+nrisq+nvarxevt &
     +ng+ng*(2*np+1)+nprmdep+nvarglob+1,(nvarprob)*(ng-1)+nrisq+nvarxevt &
      +ng+ng*(2*np+1)+nprmdep+nvarglob +ncontdep*(nq-1)+ncont*(nq-1))
	  
	  		  	  ncontr=0
	  do j=(nvarprob)*(ng-1)+nrisq+nvarxevt &
     +ng+ng*(2*np+1)+nprmdep+nvarglob+1,(nvarprob)*(ng-1)+nrisq+nvarxevt &
      +ng+ng*(2*np+1)+nprmdep+nvarglob +ncontdep*(nq-1)+ncont*(nq-1)
		if (ide(j).eq.1) then
			ncontr=ncontr+1
		end if
	  end do
	  
      read(2,*)
  !    if(CurrentProcessorID.eq.0) write(*,*)'lect',nvarprob*(ng-1)+nrisq+nvarxevt &
  !   +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1,nef

      read(2,*)(valinit(j),j=(nvarprob)*(ng-1)+nrisq+nvarxevt       & 
     +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1,nef)

	if(tbeta.eq.0)then
	  valinit((nvarprob)*(ng-1)+nrisq+nvarxevt       & 
     +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1)=0.d0
	  valinit((nvarprob)*(ng-1)+nrisq+nvarxevt       & 
     +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+2)=-dlog(2.d0)
	end if

      read(2,*)
      read(2,*)(ide(j),j=(nvarprob)*(ng-1)+nrisq+nvarxevt  & 
   +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1,nef)
      if(CurrentProcessorID.eq.0) write(*,*)'b(transfos)',(valinit(j),j=(nvarprob)*(ng-1) +nrisq+nvarxevt &
    +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1,nef)

	  		  	  ntransf=0
	  do j=(nvarprob)*(ng-1)+nrisq+nvarxevt  & 
   +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1,nef
		if (ide(j).eq.1) then
			ntransf=ntransf+1
		end if
	  end do
	  

      read(2,*)  ! indicateur de structure de variance
      read(2,*)idiag
      if(CurrentProcessorID.eq.0) write(*,*)'idiag',idiag

      if (idiag.eq.1) then
         nvc=nea
      else if(idiag.eq.0) then
         nvc=nea*(nea+1)/2
      end if

!c On entre les valeurs initiales pour les prms de
!c Variance covariance des effets aleatoires
!c on ajoute une contrainte sur la variance de l'intercept : elle est contrainte a zero

      allocate(mvc(nvc),mvcg(ng-1),vtest(nq))

      read(2,*)
      read(2,*)(mvc(j),j=1,nvc)
      read(2,*)
      read(2,*)(ide(j),j=nef+1,nef+nvc)
      if(CurrentProcessorID.eq.0) write(*,*)'mvc',(mvc(j),j=1,nvc)

	  	  		  	  nvarcov=0
	  do j=nef+1,nef+nvc
		if (ide(j).eq.1) then
			nvarcov=nvarcov+1
		end if
	  end do
	  
      read(2,*)     ! proportion par composante
      read(2,*)(mvcg(j),j=1,ng-1)
      read(2,*)
      read(2,*)(ide(j),j=nef+nvc+1,nef+nvc+ng-1)
	  
	  vcclass=0
	  do j=nef+nvc+1,nef+nvc+ng-1
		if(ide(j).eq.1)then
			vcclass=vcclass+1
		end if
	  end do
!      if(CurrentProcessorID.eq.0) write(*,*)'mvcg',(mvcg(j),j=1,ng-1)

      read(2,*)    ! ecart-type de l'effet aleatoire de chaque test
      read(2,*)(vtest(j),j=1,nq)

      read(2,*)
      read(2,*)(ide(j),j=nef+nvc+ng-1+1,nef+nvc+ng-1+nq)
	  
	  nlg=0
	  do j=nef+nvc+ng-1+1,nef+nvc+ng-1+nq
		if(ide(j).eq.1)then
			nlg=nlg+1
		end if
	  end do
	  
!C      DO k=1,nq
!C         if (idord(k).gt.1) then
!C         if (vtest(k).ne.1.or.ide(nef+nvc+ng).ne.0)
!C     &  then
!C            if(CurrentProcessorID.eq.0) write(*,*)'les tests ordinaux doivent avoir une ',
!C     &           '  variance spe du test = 1 et non estimee'
!C            stop
!C         ENDIF
!C         end if
!C      end do
      if (mvc(1).ne.1.and.ide(nef+1).ne.0)  then
         if(CurrentProcessorID.eq.0) write(*,*)'variance de intercept = 1 et non estimee'
!c         stop
      ENDIF

      read(2,*)  ! indicateur de brownien (1) et d'AR (2)
      read(2,*)iAR
      nAR=0
      if (iAR.eq.0) then
         nAR=0
      else if (iAR.eq.1) then
         nAR=1
      else if(iAR.eq.2) then
         nAR=2
      end if


      allocate(valerr(nq+nAR))
!c On entre les valeurs initiales pour les prms de brownien (sigma (se)) ou
!c AR (sigma (se) puis gamma) +erreurs ind (ecart-type) pour les nq tests
      read(2,*)
      read(2,*)(valerr(j),j=1,nAR+nq)
      read(2,*)
      read(2,*)(ide(j),j=nef+nvc+ng-1+nq+1,nef+nvc+ng-1+nq+nAR+nq)
      if(CurrentProcessorID.eq.0) write(*,*)'valerr(q)',(valerr(q),q=1,nAR+nq)
      read(2,*)
      read(2,*)nsim
      read(2,*)
      read(2,*)npg
      read(2,*)
      read(2,*)nb_cl_tps
      allocate(cl_tps(nb_cl_tps))


      read(2,*)
      read(2,*)(cl_tps(l),l=1,nb_cl_tps-1)
      read(2,*)
      read(2,*)maxiter
      read(2,*)
      read(2,*)epsa,epsb,epsd
      read(2,*)
      read(2,*),scoretest,num_frailty
      close(2)

	  nBnonHess=nvarcov+vcclass+nlg+nAR+nq

!c------------- creation du vecteur de prms -------------------
!c Transformation de cholesky pour G, matrice de variance
!c covariance des effets aleatoires
!c si idiag=1, on met dans le vecteur de parms seulement
!c les racines carrees des variances
!      if(CurrentProcessorID.eq.0) write(*,*)'Covariance parameters:',(mvc(k),k=1,nvc)
      if (idiag.eq.1) then
         DO j=1,nvc
            valinit(nef+j)=dsqrt(mvc(j))
         END DO
      end if

!c si idiag=0, on met dans le vecteur des parms, les parms
!c de la transformee de Cholesky
	ier=0
      if (idiag.eq.0) then
         CALL DMFSD(mvc,nea,EPS,IER,0)
         DO j=1,nvc
            valinit(nef+j)=mvc(j)
         END DO
!         if(CurrentProcessorID.eq.0) write(*,*)'apres cholesky',(mvc(j),j=1,nvc)
      end if
!      if(CurrentProcessorID.eq.0) write(*,*) 'Cholesky decomposition successful:',ier
!      if(CurrentProcessorID.eq.0) write(*,*)'parms var', (valinit(nef+j),j=1,nvc)

      do j=1,ng-1
         valinit(nef+nvc+j)=mvcg(j)
      end do
      do j=1,nq
         valinit(nef+nvc+ng-1+j)=vtest(j)
      end do
      do j=1,nAR+nq
         valinit(nef+nvc+ng-1+nq+j)=valerr(j)
      end do

      npmtot=nef+nvc+ng-1+nq*2+nAR

      if (npmtot.gt.nprm) then
         if(CurrentProcessorID.eq.0) write(*,*)'too many parameters compared to the max specified:',nprm,npmtot
      end if


!c creation du vecteur de prms pour marquardt :
      npm=0
      npm_ne=0
      DO j=1,npmtot
         if (ide(j).eq.1) then
            npm=npm+1
         else
            npm_ne=npm_ne+1
         end if
      END DO

      allocate(B(npm))

!c creation du vecteur de prms pour marquardt :
      npm=0
      npm_ne=0
      DO j=1,npmtot
         if (ide(j).eq.1) then
            npm=npm+1
            B(npm)=valinit(j)
         else
            npm_ne=npm_ne+1
            bne(npm_ne)=valinit(j)
         end if
      END DO

      deallocate(mvc,mvcg,vtest,valerr)

      if(CurrentProcessorID.eq.0) write(*,*)'npm=',npm
      if(CurrentProcessorID.eq.0) write(*,*)'npmtot=',npmtot
      if(CurrentProcessorID.eq.0) write(*,*)'npm_ne=',npm_ne
      if(CurrentProcessorID.eq.0) write(*,*)'b init=',(b(j),j=1,npm)
      if(CurrentProcessorID.eq.0) write(*,*)'b non est=',(bne(j),j=1,npm_ne)



!========== allocation

      allocate(numero(ns),ind_survint(ns),devt(ns),Tpred(ns),tsurv(ns),tsurv0(ns),tsurvint(ns))
      allocate(X(ns,nv),Y(ns*nq*mesure),t(ns*nq*mesure),nmes(ns,nq),Navant_i(ns),nmestotcur(ns))
      allocate(Xdep(ns*nq*mesure,nvdep))
      X=0.d0
      Y=0.d0
      t=0.d0
      Navant_i=0
      nmes=0
      Xdep=0.d0
      Devt=0

!c----------- Lecture du fichier de donnees -------------------
      open(3,file=nomdata,status='old')
      
      allocate(nb_evt(nevt+1))
      allocate(c(ns))
!vecteurs des donnees du multi etat
      allocate(idd(ns))
      allocate(idm(ns))
      allocate(t0(ns))
      allocate(t1(ns))
      allocate(t2(ns))
      allocate(t3(ns))
      nb_evt=0
      ind_survint=0
      maxmestot=0
      nmestot=0
      DO i = 1,ns
         nmestot=0
         read(3,*) numero(i)
         Do q=1,nq
            read(3,*) nmes(i,q)
            if (nmes(i,q).gt.mesure) then
               if(CurrentProcessorID.eq.0) write(*,*) 'Number of measures for the subject',i,  &
                    'for test ',q,' is greater than the maximum allowed',mesure  &
                    ,'vs.',  nmes(i,q)
               stop
            endif
            if (nmes(i,q).ne.0) then
               read(3,*)(Y(Navant_i(i)+nmestot+j),j=1,nmes(i,q))
               read(3,*)(t(Navant_i(i)+nmestot+j),j=1,nmes(i,q))
            endif
            do k=1,nvdep
               read(3,*)(Xdep(Navant_i(i)+nmestot+j,k),j=1,nmes(i,q))
            end do
            nmestot=nmestot+nmes(i,q)
         end do
         if (i.lt.ns) Navant_i(i+1)=Navant_i(i)+nmestot
         
         
         if (maxmestot.lt.nmestot) then
            maxmestot=nmestot
         end if
!lecture differents temps
         read(3,*) t0(i)!entree
         read(3,*) t1(i)!visite sain
         read(3,*) t2(i)!visite dement
!pour la censure par intervalle
         read(3,*) t3(i)!sortie


!Devt ne sera plus utilise : il ne contient que des 0, a supprimer???
!comme idevt??
         DO k = 1,nv
            read(3,*) X(i,k)
         END DO
         do k=1,nv
            if (idmort(k).eq.1) then
               idd(i)=X(i,k)
            endif
         enddo
         do k=1,nv
            if (iddem(k).eq.1) then
               idm(i)=X(i,k)
            endif
         enddo
         

!le vecteur nb_evt doit etre traite en regardant chaque parcours possible
if (nevt.eq.3)then
         if (idm(i).eq.0.and.idd(i).eq.0) THEN
            nb_evt(1)=nb_evt(1)+1
         elseif ((idm(i).eq.0).and.(idd(i).eq.1))THEN
            nb_evt(2)=nb_evt(2)+1
         elseif ((idm(i).eq.1).and.(idd(i).eq.0))THEN
            nb_evt(3)=nb_evt(3)+1
         elseif ((idm(i).eq.1).and.(idd(i).eq.1))THEN
            nb_evt(4)=nb_evt(4)+1
         end if
end if


!placement par categorie pour recuperer le parcours-->calcul de vraisemblance
	if(interv_censoring.eq.1)then
		 if((idm(i).eq.0).and.(idd(i).eq.0).and.(t1(i).eq.t3(i)))then
	!           censure a droite 01 et 02
		    c(i) = 1
		 endif
		 if((idm(i).eq.1).and.(idd(i).eq.0).and.(t1(i).lt.t2(i)))then
	!           censure par intervalle 01 et droite pour 12
		    c(i) = 2
		 endif
		 if((idm(i).eq.1).and.(idd(i).eq.0).and.(t1(i).eq.t2(i)))then
	!           observation 01 et censure a droite pour 12
		    c(i) = 3
		 endif
		 if((idd(i).eq.1).and.(idm(i).eq.1).and.(t1(i).lt.t2(i)))then
	!           censure par intervalle 01 et observation pour 12
		    c(i) = 4
		 endif
		 if((idd(i).eq.1).and.(idm(i).eq.1).and.(t1(i).eq.t2(i)))then
	!           observation 01 et observation pour 12
		    c(i) = 5
		 endif
		 if((idm(i).eq.0).and.(idd(i).eq.0).and.(t1(i).lt.t3(i)))then
	!           vivant
		    c(i) = 6
	!               write(*,*)i,t0(i),t1(i),t2(i),t3(i)
		 endif
		 if((idm(i).eq.0).and.(idd(i).eq.1))then
	!           mort
		    c(i) = 7
		 endif
	else
		 if((idm(i).eq.0).and.(idd(i).eq.0).and.(t1(i).eq.t3(i)))then
	!           censure a droite 01 et 02
		    c(i) = 1
		 endif
		 if((idm(i).eq.1).and.(idd(i).eq.0).and.(t1(i).lt.t2(i)))then
	!           censure par intervalle 01 et droite pour 12
		    t1(i)=(t1(i)+t2(i))/2
		    t2(i)=t1(i)
		    c(i) = 3
		 endif
		 if((idm(i).eq.1).and.(idd(i).eq.0).and.(t1(i).eq.t2(i)))then
	!           observation 01 et censure a droite pour 12
		    c(i) = 3
		 endif
		 if((idd(i).eq.1).and.(idm(i).eq.1).and.(t1(i).lt.t2(i)))then
	!           censure par intervalle 01 et observation pour 12
		    t1(i)=(t1(i)+t2(i))/2
		    t2(i)=t1(i)
		    c(i) = 5
		 endif
		 if((idd(i).eq.1).and.(idm(i).eq.1).and.(t1(i).eq.t2(i)))then
	!           observation 01 et observation pour 12
		    c(i) = 5
		 endif
		 if((idm(i).eq.0).and.(idd(i).eq.0).and.(t1(i).lt.t3(i)))then
	!           vivant
		    t1(i)=t3(i)
		    t2(i)=t3(i)
		    c(i) = 7
	!               write(*,*)i,t0(i),t1(i),t2(i),t3(i)
		 endif
		 if((idm(i).eq.0).and.(idd(i).eq.1))then
	!           mort
		    t1(i)=t3(i)
		    t2(i)=t3(i)
		    c(i) = 7
		 endif
		if(t1(i).ne.t2(i))then
			print*,'error for uncensored data'
			stop
		end if
	end if
      END DO

	if (tbeta.eq.0)then
	  range(1)=0.d0
	  range(2)=1.d0
	end if
!c-----------------------------------------------------------------------------
      close(3)
      if(CurrentProcessorID.eq.0) write(*,*),"Vivants Sains // Morts Sains // Dements vivants // Morts Dements // Total"
      if(CurrentProcessorID.eq.0) write(*,*),nb_evt,ns


!bonne gestion dans le modele multi etats 3 transitions
      if(CurrentProcessorID.eq.0) write(*,*)'Number of events:',nb_evt
      if (minval(nb_evt).eq.0) then 
         if(CurrentProcessorID.eq.0) write(*,*)'one of the competing event has no event recorded'
         stop 
      end IF 

!      if(CurrentProcessorID.eq.0) write(*,*)'npm=',npm,'npmtot=',npmtot,'npm_ne=',npm_ne
!      if(CurrentProcessorID.eq.0) write(*,*)'nombre de evenement',nb_evt


      do k=1,nv
         if(CurrentProcessorID.eq.0) write(*,*)'Mean for covariate ',k,':',sum(X(:,k))/dble(ns)
      end DO

!c---------- Ouverture fichier des resultats -------------------
      open(8,file=nomoutput)
      write(8,*)'        DATA             '
      write(8,*)'Datafile : ',nomdata
      write(8,*)'N= ',ns,'Q=',nq
      nobsmoy=0
      nobs=0
      DO i=1,ns
         do q=1,nq
            nobs=nobs+nmes(i,q)
         end do
      END DO
      nobsmoy=nobs/ns

      if (Navant_i(ns)+nmestot.ne.nobs) then 
         if(CurrentProcessorID.eq.0) write(*,*)'problem def Navant_i',Navant_i(ns),nobs
         stop
      end if
      write(8,*)'Average number of observations per subject :',nobsmoy
      write(8,*)
      write(8,*)'Indicator if the esitmation procedure accounts for interval censoring:', interv_censoring
      write(8,*)'Number of mixture components : ',ng
	  write(8,*)'Type of baseline risk function and hazard function : ', typrisq,risqcom
      write(8,*)'Degree of the polynomious time function :',np
      write(8,*)'Indicator of random effect on the time function',(idea(k),k=1,np+1)
      write(8,*)'Total number of covariates in the model : ',nv
      write(8,*)'Indicator that the covariate is in the latent part of the model', &
           (idglob(k),k=1,nv)
      write(8,*)'Indicator that the covariate is in the logit model',(idprob(k),k=1,nv)
!      write(8,*)'Indicator that the covariate is the binary outcome',(idevt(k),k=1,nv)
      write(8,*)'Indicator that the covariate is the death indicator',(idmort(k),k=1,nv)
      write(8,*)'Indicator that the covariate is the dementia indicator',(iddem(k),k=1,nv)
      write(8,*)'Indicator that the covariate is in the event model',(idxevt(k),k=1,nv*(nevt))
      write(8,*)'Indicator that the time-dependent covariates are included' ,(idvdep(k),k=1,nvdep)
      write(8,*)'Indicator that the covariate has contrasts on the tests',(idtest(k),k=1,nv)
      write(8,*)'Indicator of the covariance matrix structure:',idiag
      write(8,*)'Indicator of intra-correlation structure (0:No,1:BM,2:AR)',iAR
      write(8,*)'Number of parameters : ',npm
	  write(8,*) 'Initial values for parameters : ',(valinit(i),i=1,npmtot)
	   write(8,*)'Estimated parameters : ',(ide(i),i=1,npmtot)

!c ----------- Optimisation de Marquardt ------------------
     IF(CURRENTPROCESSORID.EQ.0) WRITE(*,*) 'npmtot',npmtot,'nef',nef,'nea',nea,'nvarglob',nvarglob  &
       ,'evt',evt,'nvarprob',nvarprob,'nprmdep',nprmdep,'ncont'  &
         ,ncont,'ncontdep',ncontdep  &
         ,'nef',nef,'npmtot',npmtot,'npm_ne',npm_ne,'nvarxevt'  &
         ,nvarxevt,'nvc',nvc,'nAR',nAR

!C min max
!###
if(nevt.ne.0)then
!#########################
!A time_est for the latent process
!#########################
        if (idtrunc.eq.0) then
           ndata=ns
       else
           ndata=ns*2
       end if
	   
       allocate(time_est(nsim),ytemp2(ndata))
        
       time_est=0.d0
         do i=1,ns
            Ytemp2(i)=t3(i)
            if (idtrunc.eq.1) then
               Ytemp2(ns+i)=t0(i)
            end if
         end do
!###         
!         C     appel de la procedure de tri d'un vecteur
!         ind=1
!         do while (ind.eq.1)
!            ind=0
!            do i=1,ndata-1
!               if (ytemp2(i).gt.ytemp2(i+1)) then
!                  temp=ytemp2(i)
!                  ytemp2(i)=ytemp2(i+1)
!                  ytemp2(i+1)=temp
!                  ind=1
!               end if
!            end do
!         end do
 
         !if(CurrentProcessorID.eq.0) write(*,*)'Tsurv',Tsurv
!         if(CurrentProcessorID.eq.0) write(*,*)'ytemp',ytemp2
!###
         min_save=0.d0
         if (idtrunc.eq.1) then 
            min_save=minval(ytemp2)
         end IF
         max_save=maxval(ytemp2)
         time_est(1)=min_save
         do j=2,nsim
            time_est(j)=time_est(j-1)+(max_save-min_save)/dble(nsim-1)
         end do
         time_est(nsim)=max_save
         
         deallocate(Ytemp2) 

!      !C         min_save=0.d0
       
         !if(CurrentProcessorID.eq.0) write(*,*)'Time of event: min=',min_save,'max=',max_save

!#########################
!A Nodes for step function for hazard
!#########################
	   allocate(zi01(nz(1)),zi02(nz(2)),zi12(nz(3)))
	   zi01=0.d0
	zi02=0.d0
	zi12=0.d0
      if (typrisq(1).eq.1) then
		   if (idzi(1).eq.1) then 
				call zi_noeuds(nz(1),idzi(1),size(t0),t3,zi01)

				if(CurrentProcessorID.eq.0) write(*,*)'Equidistant nodes for step function on time to healthy-dement event '
				if(CurrentProcessorID.eq.0) write(*,*)'zi01',zi01
				write(8,*)'zi01',zi01
			else
				allocate(Tevt(nb_evt(3)))
				numact2=1
				do i=1,ns
					 if(idm(i).eq.1.and.idd(i).eq.0)then
						Tevt(numact2)=(t1(i)+t2(i))/2
						numact2=numact2+1
					end if
				end do
				call zi_noeuds(nz(1),idzi(1),nb_evt(3),Tevt,zi01)
				if(CurrentProcessorID.eq.0) write(*,*)'Quantile nodes for step function on time to healthy-dement event '
				if(CurrentProcessorID.eq.0) write(*,*) 'zi01',zi01
				write(8,*)'zi01',zi01
				deallocate(Tevt)
			end if
	end if
		
	if(nevt.ge.2)then
		if (typrisq(2).eq.1) then
			   if (idzi(2).eq.1) then 
					call zi_noeuds(nz(2),idzi(2),size(t0),t3,zi02)

					if(CurrentProcessorID.eq.0) write(*,*)'Equidistant nodes for step function on time to healthy-dead event '
					if(CurrentProcessorID.eq.0) write(*,*)'zi02',zi02
					write(8,*)'zi02',zi02
				else
					allocate(Tevt(nb_evt(2)))
					numact2=1
					do i=1,ns
						 if(idm(i).eq.0.and.idd(i).eq.1)then
							Tevt(numact2)=t3(i)
							numact2=numact2+1
						end if
					end do
					call zi_noeuds(nz(2),idzi(2),nb_evt(2),Tevt,zi02)
					if(CurrentProcessorID.eq.0) write(*,*)'Quantile nodes for step function on time to healthy-dead event '
					if(CurrentProcessorID.eq.0) write(*,*) 'zi02',zi02
					write(8,*)'zi02',zi02
					deallocate(Tevt)
				end if
			end if
	end if
	
	if(nevt.ge.3)then
		 if (typrisq(3).eq.1) then
		 		if (idzi(3).eq.1) then 
					call zi_noeuds(nz(3),idzi(3),size(t0),t3,zi12)

					if(CurrentProcessorID.eq.0) write(*,*)'Equidistant nodes for step function on time to dement-dead event '
					if(CurrentProcessorID.eq.0) write(*,*)'zi12',zi12
					write(8,*)'zi12',zi12
				else
					allocate(Tevt(nb_evt(4)))
					numact2=1
					do i=1,ns
						 if(idm(i).eq.1.and.idd(i).eq.1)then
							Tevt(numact2)=t3(i)
							numact2=numact2+1
						end if
					end do
					call zi_noeuds(nz(3),idzi(3),nb_evt(4),Tevt,zi12)
					if(CurrentProcessorID.eq.0) write(*,*)'Quantile nodes for step function on time to dement-dead event '
					if(CurrentProcessorID.eq.0) write(*,*) 'zi12',zi12
					write(8,*)'zi12',zi12
					deallocate(Tevt)
				end if
			end if
		end if
  end if
	  
	  
      !C splines a rechercher:
       deallocate(nb_evt) 
  
      !C splines a rechercher:

      if (spline.eq.1) then


         allocate(mm(ns*nq*mesure),mm1(ns*nq*mesure),mm2(ns*nq*mesure), &
              im(ns*nq*mesure),im1(ns*nq*mesure),im2(ns*nq*mesure))
       

         call splines_tr()
!         do i=395,396
!            do j=1,nmes(i,1)
!               if(CurrentProcessorID.eq.0) write(*,*)'i',i,j,mm(Navant_i(i)+j),mm1(Navant_i(i)+j),mm2(Navant_i(i)+j) 
!               if(CurrentProcessorID.eq.0) write(*,*)'i',i,j,im(Navant_i(i)+j),im1(Navant_i(i)+j),im2(Navant_i(i)+j)
!            end do
!         end do




         write(8,*)'Nodes for splines transformations'
         do q=1,nq
            write(8,*)(zitr(j,q),j=-1,nztr+2)
         end do
         if(CurrentProcessorID.eq.0) write(*,*)'Nodes for splines transformations'
         do q=1,nq
            if(CurrentProcessorID.eq.0) write(*,*)(zitr(j,q),j=-1,nztr+2)
         end do
      end if
      
      
      ca=1.d-4
      cb=1.d-4
      dd=1.d-4

      allocate(vopt(npm*(npm+3)/2))



      allocate(Ut1(nea,nea),mu(maxmestot),mu2(maxmestot),Y1(maxmestot),Z1(maxmestot,nea),sigmaq(maxmestot,maxmestot) &
       ,Vw(maxmestot,maxmestot),b1(npmtot),Y1all(ns*nq,maxmestot))

!goto 1111
	   call marq98(b,npm,ni,Vopt,vrais,ier,istop,ca,cb,dd,funcpa,ng,nvc)
!1111 continue	


!       !C *** recreation du set de parm
      allocate(pipre(ng))
      if (ng.eq.1) then
         pipre=1.d0
      else
         sumvectpi=1.d0
         do g=1,ng-1
            sumvectpi=sumvectpi+exp(b(g)) !possible car b commence par les parametres de classe
         enddo

    
         
         pipre(ng)=1.0/sumvectpi
         do g=1,ng-1
            pipre(g)=exp(b(g))/sumvectpi
         enddo
       
      endif
      b1=0.d0
      l=0
      m=0

      do k=1,npmtot
         if (ide(k).eq.1) then
            l=l+1
            b1(k)=b(l)
         else
            m=m+1
            b1(k)=bne(m)
         end if
      end do
      
	  
      write(8,*)'convergence criteria',ca,cb,dd
      if(CurrentProcessorID.eq.0) write(*,*)
      if(CurrentProcessorID.eq.0) write(*,*)'       RESULTS'
      if(CurrentProcessorID.eq.0) write(*,*)
      
      if(CurrentProcessorID.eq.0) write(*,*)'loglikelihood', vrais
      if(CurrentProcessorID.eq.0) write(*,*)'number of iterations',ni
      if(CurrentProcessorID.eq.0) write(*,*)'Convergence criteria',istop
      write(8,*)
      write(8,*)'       RESULTS'
      write(8,*)
      write(8,*)'LogLikelihood l=',vrais,'/number of iterations:',ni
      aic=-2*vrais+2*npm
      bic=-2*vrais+log(dble(ns))*npm
      write(8,*)'-2l : ',-2*vrais
      write(8,*)'AIC : ',aic
      write(8,*)'BIC : ',bic
      
      if(CurrentProcessorID.eq.0) write(*,*)'AIC : ',aic
      if(CurrentProcessorID.eq.0) write(*,*)'BIC : ',bic
      if(CurrentProcessorID.eq.0) write(*,*)'convergence criteria',ca,cb,dd
      
      If (istop.eq.1) then
         write(8,*)'Convergence criteria OK, istop :',istop
      Elseif (istop.eq.2) then
         write(8,*)'Maximum number of iterations reached'
      Elseif (istop.eq.3) then
         write(8,*)'Convergence criteria OK but Hessian is singular'
      Elseif (istop.eq.4) then
         write(8,*)'LogLikelihood not improved after SEARPAS '
	 Elseif (istop.eq.5) then
         write(8,*)'Convergence criteria OK with partial Hessian '
      endif
      write(8,*)
      write(8,*)'Transformed parameter estimations : '
      write(8,*)'  esti,     se(esti),Z=esti/se ,95%IC: binf,bsup '
      write(8,*)
      id=1
      if(CurrentProcessorID.eq.0) write(*,*)
      if(CurrentProcessorID.eq.0) write(*,*)'Transformed parameter estimations :'

      allocate(se(npm))
      if (ng.gt.1) then
         write(8,*)"Class membership"
         contsortie=contsortie+nprob+1
	else
		contsortie=1
      endif
	  
      se=0.d0
	  DO i=1,npm
      if (istop.eq.1)then 
		    se(i)=dsqrt(Vopt(id))
      else if(istop.eq.5)then
		if (i.le.(npm-nBnonHess))then
		    se(i)=dsqrt(Vopt(id))
		 else
		    se(i)=1
		 end if
	else
	      se(i)=1
	end if

         if(CurrentProcessorID.eq.0) write(*,*)B(i)
         binf=B(i)-1.96*se(i)
         bsup=B(i)+1.96*se(i)
         if (i.eq.contsortie.and.contfins.eq.0) then
            write(8,*)"Basis Risks of event"
            contsortie=contsortie+nhazard
            contfins=1
         endif
         if (i.eq.contsortie.and.contfins.eq.1) then
            write(8,*)"Event covariates"
            contsortie=contsortie+ncov
            contfins=2
         endif
         if (i.eq.contsortie.and.contfins.eq.2) then
            write(8,*)"Changepoint times "
            contsortie=contsortie+ng!depend du nombre estime
            contfins=25
         endif
         if (i.eq.contsortie.and.contfins.eq.25) then
            write(8,*)"Latent process"
            contsortie=contsortie+nlatent!depend du nombre estime
            contfins=3
         endif
		 if (i.eq.contsortie.and.contfins.eq.3) then
            write(8,*)"Marker-specific parameters"
            contsortie=contsortie+ncontr !marche que pour beta
            contfins=4
         endif

         if (i.eq.contsortie.and.contfins.eq.4) then
            write(8,*)"Transformation parameters"
            contsortie=contsortie+ntransf !marche que pour beta
            contfins=5
         endif
         if (i.eq.contsortie.and.contfins.eq.5) then
            write(8,*)"Variance-covariance matrix"
            contfins=6
            contsortie=contsortie+nvarcov !depend du nombre estime
         endif
         if (i.eq.contsortie.and.contfins.eq.6) then
  write(8,*)"Class-spec proportional parameter VC matrix"
            contfins=7
            contsortie=contsortie+vcclass !depend du nombre estime
         endif
         if (i.eq.contsortie.and.contfins.eq.7) then
            write(8,*)"Others"
            contfins=7
         endif


         write (8,2100) B(i),se(i),B(i)/se(i),binf,bsup
         id=id+i+1
      end do

      deallocate(se)

2100  FORMAT(F10.6,2X,F10.6,2X,F8.4,2X,F8.4,2X,F8.4)
     
!c ------------- matrice de variance covariance -----------------------
      
      allocate(Ut(nea,nea),U(nea,nea),VC(nea,nea),CORR(nea,nea),bh1(nvc))
      allocate(ppi(ns,ng),ppitest(ns,ng),probapost(ns,ng),P_AI(ns,ng),pi(ns,ng))

      P_AI=0.d0
      Pi=0.d0
      bh1=0.d0
      do k=1,nvc
         bh1(k)=b1(nef+k)
      end do      
      U=0.d0
      Ut=0.d0
      If (idiag.eq.1) then
         do j=1,nea
            U(j,j)=bh1(j)
            Ut(j,j)=bh1(j)
         end do
      end if
      If (idiag.eq.0) then
         do j=1,nea
            do k=j,nea
               U(j,k)=bh1(j+k*(k-1)/2)
               Ut(k,j)=bh1(j+k*(k-1)/2)
            end do
         end do
      end if

      VC=Matmul(Ut,U)
      if(CurrentProcessorID.eq.0) write(*,*)
      if(CurrentProcessorID.eq.0) write(*,*)'Variance-covariance matrix :'
      do i=1,nea
         if(CurrentProcessorID.eq.0) write(*,*)(VC(i,j),j=1,i)
      end do
      write(8,*)
      write(8,*)'Variance-covariance matrix :'
      do i=1,nea
         write(8,*)(VC(i,j),j=1,i)
      end do
	
      do i=1,nea
	do j=1,nea
	  CORR(i,j)=VC(i,j)*((VC(i,i)*VC(j,j))**(-0.5d0))
	end do
      end do
      write(8,*)
      write(8,*)'Correlation matrix :'
      do i=1,nea
         write(8,*)(CORR(i,j),j=1,i)
      end do

      deallocate(Ut,U,VC,bh1,CORR)
      
      
     write(8,*)
     write(8,*)'Transformed parameter estimations :'
     DO i=1,npm
        write(8,*)B(i)
     end do
      
     write(8,*)

	      
     if (istop.ne.3) then


        open(345,file=nomoutput2)
  
        write(345,*)'Number of outcomes'
        write(345,*)nq
        write(345,*)'Indicator for transformation'
        write(345,*)idord
        write(345,*)'Range of outcomes for estimation'
        write(345,*)range
        write(345,*)'Real range of outcomes'
        write(345,*)rangereal
        write(345,*)'Degree of the polynomious time functions'
        write(345,*)np
        write(345,*)'Indicator of presence of random-effect'
        write(345,*)idea
        write(345,*)'Number of latent classes'
        write(345,*)ng
        write(345,*)'number of events'
        write(345,*)nevt
        if(nevt.gt.0) then
           write(345,*)'Indicators for the survival part'
           write(345,*)typrisq,risqcom,idtrunc,nz,idzi,paramsurv
	 if(maxval(nz).gt.0)then
           write(345,*) 'nodes for 01 : '
           write(345,*) zi01
	   write(345,*) 'nodes for 02 : '
           write(345,*)zi02
	   write(345,*) 'nodes for 12 : '
           write(345,*)zi12
        end if
	end if
        write(345,*)'Number of nodes and below list of nodes for Splines transformations'
        write(345,*)nztr
        do q=1,nq
           if (idord(q).eq.1) then
              write(345,*)(zitr(j,q),j=2,nztr-1)
	   end if
	end do
        write(345,*)'Number and indicator of inclusion for an intermediate survival time'
        write(345,*)nvardepsurv,nvdepsurv
        write(345,*)'Number of time-dependent covariates'
        write(345,*)nvdep
        write(345,*)'Indicator of time-dependent covariates inclusion'
        write(345,*)(idvdep(j),j=1,nvdep)
        write(345,*)'Number of time-independent covariates'
        write(345,*)nv
        write(345,*)'Indicator for inclusion in the latent process model'
        write(345,*)idglob
        write(345,*)'Indicator for inclusion of class-specific effects'
        write(345,*)idg
        write(345,*)'Indicator for inclusion in the class-membership model'
        write(345,*)idprob
        write(345,*)'Indicator of death variable'
        write(345,*)idmort
        write(345,*)'Indicator of dementia variable'
        write(345,*)iddem
        write(345,*)'Indicator for inclusion in the survival model'
        write(345,*)idxevt
        write(345,*)'Indicator for inclusion as an outcome-specific contrast'
        write(345,*)idtest
        write(345,*)'Indicator for diagonal or unstructured matrix of Var Cov'
        write(345,*)idiag
        write(345,*)'Indicator for autoregressive process'
	 write(345,*)iAR
	 write(345,*)'Parameter for smooth changepoint'
	 write(345,*) trn
	 write(345,*)'Indicators for parameters estimated'
	 write(345,*)(ide(j),j=1,npmtot)	
	 write(345,*)'Vector of parameters not estimated'
	 write(345,*)(bne(j),j=1,npm_ne)
	 write(345,*)'Vector of parameters estimated'
	 write(345,*)(b(i),i=1,npm)
	 write(345,*)'Matrix of Variance-Covariance'
	








		
 	  






         allocate(V1(npm,npm))

         V1=0.d0
         do j=1,npm
            do k=j,npm
               id=j+k*(k-1)/2
               V1(j,k)=Vopt(id)
               V1(k,j)=Vopt(id)
            end do
         end do
!         write(8,*)'correlation matrix for the parameters'
         do j=1,npm
!           write(8,*)(V1(j,k),k=1,npm)
            write(345,*)(V1(j,k),k=1,npm)
         end do

         close(345)
!###


!    WALD MULTIVARIE par variable avec interaction avec le temps



      if (istop.eq.1) then 
         
         if(CurrentProcessorID.eq.0) write(*,*)
         if(CurrentProcessorID.eq.0) write(*,*)'       MULTIVARIATE WALD TESTS'
         if(CurrentProcessorID.eq.0) write(*,*)
         if (np.gt.1) then
            if(CurrentProcessorID.eq.0) write(*,*)'===> Wald test for the overall effect of ', &
               'covariate interaction with the functions of time'
            allocate(b3(np),b4(np),varM(np,np),varv(np*(np+1)/2), &
            b3g(np*ng),b4g(np*ng),varMg(np*ng,np*ng) &
            ,varvg(np*ng*(np*ng+1)/2))
            indice=0
            do i=1,nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep
               if (ide(i).eq.1) then 
                  indice=indice+1
               end if
            end do
            jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep
            kk=0     
            do k=1,nv
               if (idglob(k).eq.2) then  
                  if (idg(k).eq.0) then 
                     b3=0.d0
                     kk=kk+1
                     indice=indice+ide(jj+kk)
                     ind=0
                     do j=1,np
                        ind=ind+ide(jj+kk+j)
                     end do
                     if(ind.eq.np) then  
                        Do m=1,np
                           b3(m)=b1(jj+kk+m)
                           do m1=1,np
                              varM(m,m1)=V1(indice+m,indice+m1)
                           End do
                        end do
                        if(CurrentProcessorID.eq.0) write(*,*)'  * Covariate #',k
                        if(CurrentProcessorID.eq.0) write(*,*)'    B=',(b3(j),j=1,np)
!                        if(CurrentProcessorID.eq.0) write(*,*)'variance'                        
                        do m=1,np
!                           if(CurrentProcessorID.eq.0) write(*,*)(varM(m,m1),m1=1,np)                     
                           do m1=m,np
                              varV(m+m1*(m1-1)/2)=VarM(m,m1)
                           end do
                        end do
                        
                        ier=0
                        det=0.d0
                        CALL DSINV(VarV,np,eps,ier,det,0)   
                        if (ier.eq.-1) then
                           if(CurrentProcessorID.eq.0) write(*,*)'dsinv does not work for Wald test'
                           stop
                        end if
                        VarM=0.d0
                        do m=1,np
                           do m1=1,np
                              if (m1.ge.m) then 
                                 VarM(m,m1)=VarV(m+m1*(m1-1)/2)
                              else
                                 VarM(m,m1)=VarV(m1+m*(m-1)/2)
                              end if
                           end do
                        end do    
                        b4=MATMUL(varM,b3)
                        Zscore=DOT_PRODUCT(b3,b4)
                  if(CurrentProcessorID.eq.0) write(*,*)'    Wald stat (chi-square with np df)',Zscore          
                     end if
                     indice=indice+ind
                     kk=kk+np                 
                     
                  elseif (idg(k).eq.1) then

                     b3g=0.d0
!                     if(CurrentProcessorID.eq.0) write(*,*)'---- Global Test over ng and np'
                     ind=0
                     do j=1,np*ng
                        ind=ind+ide(jj+kk+ng+j)
                     end do
                     if(ind.eq.np*ng) then               
                        indice=indice+sum(ide(jj+kk+1:jj+kk+ng))
                        varMg=0.d0
                        Do m=1,np*ng
                           b3g(m)=b1(jj+kk+ng+m)
                           do m1=1,np*ng
                              varMg(m,m1)=V1(indice+m,indice+m1)
                           End do
                        end do
                        
                        if(CurrentProcessorID.eq.0) write(*,*)'  * Covariate #',k
                        if(CurrentProcessorID.eq.0) write(*,*)'    B=',(b3g(j),j=1,np*ng)
!                        if(CurrentProcessorID.eq.0) write(*,*)'    variance'                        
                        do m=1,np*ng
!                           if(CurrentProcessorID.eq.0) write(*,*)'    ',(varMg(m,m1),m1=1,np*ng)                 
                           do m1=m,np*ng
                              varVg(m+m1*(m1-1)/2)=VarMg(m,m1)
                           end do
                        end do 
                        
                        ier=0
                        det=0.d0
                        CALL DSINV(VarVg,np*ng,eps,ier,det,0)   
                        if (ier.eq.-1) then
                           if(CurrentProcessorID.eq.0) write(*,*)'dsinv does not work for Wald test'
                           stop
                        end if
                        VarMg=0.d0
                        do m=1,np*ng
                           do m1=1,np*ng
                              if (m1.ge.m) then 
                                 VarMg(m,m1)=VarVg(m+m1*(m1-1)/2)
                              else
                                 VarMg(m,m1)=VarVg(m1+m*(m-1)/2)
                              end if
                           end do
                        end do    
                        b4g=MATMUL(varMg,b3g)
                        Zscore=DOT_PRODUCT(b3g,b4g)
               ! wortie
                        if(CurrentProcessorID.eq.0) write(*,*)'    Wald stat (chi-square ng*np df) =',Zscore  

                     end if
                     indice=indice+ind
                     kk=kk+(np+1)*ng  
                     
                  end if
 

               else if(idglob(k).eq.1.and.idg(k).eq.0) then
                  indice=indice+ide(jj+kk)
                  kk=kk+1
               Else if (idglob(k).eq.1.and.idg(k).eq.1) then
                  indice=indice+sum(ide(jj+kk+1:jj+kk+ng))
                  kk=kk+ng  
               end if
               
            end do
            if(CurrentProcessorID.eq.0) write(*,*)
               
            deallocate(b3,b4,varM,varv,b3g,b4g,varMg,varvg)

         end if


      end if

         if (ng.gt.1) then


            allocate(b3(ng),b4(ng),varM(ng,ng),varv(ng*(ng+1)/2))
 
            if(CurrentProcessorID.eq.0) write(*,*)'===> Wald test for the overall effect of', &
               ' class-specific covariate effects'

            jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep
            indice=0
            do i=1,jj
               if (ide(i).eq.1) then 
                  indice=indice+1
               end if
            end do
            kk=0     
            do k=1,nv
               if (idg(k).eq.1) then
                  cursor=0
!                  if(CurrentProcessorID.eq.0) write(*,*)'covariate',k
!                  if(CurrentProcessorID.eq.0) write(*,*)'indices',indice,jj+kk,cursor
                  do while (cursor.lt.np+1)
                     cursor=cursor+1
                     if (cursor.eq.1) then
                        if(CurrentProcessorID.eq.0) write(*,*)'  * For covariate #:',k
                     else
                        if(CurrentProcessorID.eq.0) write(*,*)'  * For covariate #:',k, &
                           'in interaction with t^',cursor-1
                     end if 
                     b3=0.d0
                     ind=0
                     do j=1,ng
                        ind=ind+ide(jj+kk+j)
                     end do
                     if(ind.eq.ng) then  
                        Do m=1,ng
                           b3(m)=b1(jj+kk+m)
                           do m1=1,ng
                              varM(m,m1)=V1(indice+m,indice+m1)
                           End do
                        end do                 
                        if(CurrentProcessorID.eq.0) write(*,*)'    B=',(b3(j),j=1,ng)
!                        if(CurrentProcessorID.eq.0) write(*,*)'    variance'                        
                        do m=1,ng
!                           if(CurrentProcessorID.eq.0) write(*,*)'    ',(varM(m,m1),m1=1,ng)                     
                           do m1=m,ng
                              varV(m+m1*(m1-1)/2)=VarM(m,m1)
                           end do
                        end do
                        ier=0
                        det=0.d0
                        CALL DSINV(VarV,ng,eps,ier,det,0)   
                        if (ier.eq.-1) then
                           if(CurrentProcessorID.eq.0) write(*,*)'dsinv does not work for Wald test'
                           stop
                        end if
                        VarM=0.d0
                        do m=1,ng
                           do m1=1,ng
                              if (m1.ge.m) then 
                                 VarM(m,m1)=VarV(m+m1*(m1-1)/2)
                              else
                                 VarM(m,m1)=VarV(m1+m*(m-1)/2)
                              end if
                           end do
                        end do    
                        b4=MATMUL(varM,b3)
                        Zscore=DOT_PRODUCT(b3,b4)
                        if(CurrentProcessorID.eq.0) write(*,*)'    Wald stat (chi-square with ng df) =',Zscore     
                             
                     end if
                     indice=indice+ind
                     kk=kk+ng    
                     if (idglob(k).eq.1) then
                        cursor=np+1
                     end if        
                     
                  end do
               elseif (idg(k).eq.0.and.idglob(k).eq.1) then
                  indice=indice+ide(jj+kk+1)
                  kk=kk+1
               elseif (idg(k).eq.0.and.idglob(k).eq.2) then
                  indice=indice+sum(ide(jj+kk+1:jj+kk+np+1))
                  kk=kk+np+1
 !                 if(CurrentProcessorID.eq.0) write(*,*)'indice',indice,sum(ide(jj+kk+1:jj+kk+np+1))
               end if
            end do
            deallocate(b3,b4,varM,varv)
         end if


         deallocate(V1)

      end if

	if(scoretest.eq.1)then
	allocate(uscoreYT(nea),var_uscoreYT(nea,nea),var_YT(nea,nea),var_uYT(nea*(nea-1)))
	UscoreTT=0.d0
	var_uscoreTT=0.d0
	UscoreYT=0.d0
	var_uscoreYT=0.d0
	var_YT=0.d0
	

        call postprobmult_beta(B,npm,PPI,ppitest,probapost,UscoreTT,var_uscoreTT,UscoreYT,var_uscoreYT,num_frailty)


	test_scoreTT=UscoreTT/var_uscoreTT*uscoreTT


!c     Vi en vecteur
	var_uYT=0.d0
         jj=0
         do j=1,nea
            do k=j,nea
               jj=j+k*(k-1)/2
               var_uYT(jj)=var_uscoreYT(j,k)
            end do
         end do

         CALL DSINV(var_uYT,nmestot,eps,ier,det,0)

         if (ier.eq.-1) then
          !  if(CurrentProcessorID.eq.0) write(*,*)'dsinv vi ier',ier,' i=',i,' nmes(i)=',nmestot
          !  if(CurrentProcessorID.eq.0) write(*,*)'b1', (b1(k),k=1,npm)

            if(CurrentProcessorID.eq.0) write(*,*)'var_uYT', (var_uYT(k),k=1,nea*(nea+1)/2)
	    write(*,*)'var_uYT'
            write(8,*)'var_uYT non inversible '!:stop
         end if

	var_YT=0.d0

       do m=1,nea
           do m1=1,nea
              if (m1.ge.m) then 
                 var_YT(m,m1)=var_uYT(m+m1*(m1-1)/2)
              else
                 var_YT(m,m1)=var_uYT(m1+m*(m-1)/2)
              end if
           end do
        end do 

	test_scoreYT=UscoreYT(1)*(UscoreYT(1)*var_YT(1,1)+UscoreYT(2)*var_YT(2,1)+UscoreYT(3)*var_YT(3,1))
	test_scoreYT=test_scoreYT+UscoreYT(2)*(UscoreYT(1)*var_YT(1,2)+UscoreYT(2)*var_YT(2,2)+UscoreYT(3)*var_YT(3,2))
	test_scoreYT=test_scoreYT+UscoreYT(3)*(UscoreYT(1)*var_YT(1,3)+UscoreYT(2)*var_YT(2,3)+UscoreYT(3)*var_YT(3,3))

         write(8,*)
	 write(8,*) 'Score test for independence between events (<2.70)'
	 write(8,*)
	 write(8,*) 'Test statistic :', test_scoreTT, 'Score :',uscoreTT,   'variance :', var_uscoreTT
	 write(8,*)
	 write(8,*) 'Score test for independence between Y and T (<3.84)'
	 write(8,*)
	 write(8,*)'Test statistic :', test_scoreYT,  'Score :',uscoreYT,  'variance :', var_uscoreYT
	 write(8,*)

	end if
!C------ estimation de la variance des probabilites des composantes

      if (ng.gt.1) then

               if(CurrentProcessorID.eq.0) write(*,*)
               if(CurrentProcessorID.eq.0) write(*,*)'       POSTERIOR CLASSIFICATION'
               if(CurrentProcessorID.eq.0) write(*,*)

1998     FORMAT(I5,1X,I3,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X, &
              F9.6,1X,F9.6,1X,F9.6)
         
         allocate(classe(ns),classe1(ns),classement(ng,4),classement1(ng,4),pipost(ng),pitest(ng))
         allocate(sumvect(ng,ng),tfinmoy(ng,ng))
         allocate(claMCpl(ng,4),claFCpl(ng,4),claMCmo(ng,4),claFCmo(ng,4))
         if(CurrentProcessorID.eq.0) write(*,*)'typrisq',typrisq

	!allocate(uscore(ns))

         if(CurrentProcessorID.eq.0) write(*,*)'postprob over'
         open(7,file='posterior_class_probs.txt')
         icl=0.d0
         sumvect=0.d0
         tfinmoy=0.d0
         DO i=1,ns
            classe(i)=1
            classe1(i)=1
            DO g=1,ng
               if (ppi(i,g).gt.ppi(i,classe(i))) then
                  classe(i)=g
               end if
               if (ppitest(i,g).gt.ppitest(i,classe1(i))) then
                  classe1(i)=g
               end if
            END DO
            if (g.ge.ng) then
               g=ng
            endif
            if (ppi(i,g).gt.0) then
               icl=icl+ppi(i,g)*log(ppi(i,g)) ! erreur g deviens > ng, correction inutile?
            end IF
         END DO
         pipost=0
         DO i=1,ns
            pipost(classe(i))=pipost(classe(i))+1
         END DO
         pitest=0.d0
         DO i=1,ns
            pitest(classe1(i))=pitest(classe1(i))+1
         END DO
         
         classement=0
         classement1=0
         claMCpl=0
         claMCmo=0
         claFCpl=0
         claFCmo=0
         if(CurrentProcessorID.eq.0) write(*,*)
         if(CurrentProcessorID.eq.0) write(*,*) 'BIC-ICL=',bic-2*icl
         if(CurrentProcessorID.eq.0) write(*,*) 'Number of subjects per component '
         if(CurrentProcessorID.eq.0) write(*,*) (pipost(g),g=1,ng)
         if(CurrentProcessorID.eq.0) write(*,*) (pitest(g),g=1,ng)
         write(8,*)
	 write(8,*) 'POSTERIOR CLASSIFICATION (P(c(i)=g|Y_i,T_i))'
	 write(8,*)
         write(8,*) 'Number of subjects per mixture component'
         write(8,*) (pipost(g),g=1,ng)

         do i=1,ns
            do g=1,ng
               sumvect(classe(i),g)=sumvect(classe(i),g)+PPI(i,g)
            enddo
!            write(7,1998)numero(i),(PPI(i,g),g=1,ng),classe(i)
            write(7,*)(PPI(i,g),g=1,ng),numero(i),classe(i)

            do g=1,ng
                if(classe(i).eq.g.and.idm(i).eq.1.and.idd(i).eq.0) then
                  classement(g,1)=classement(g,1)+1
		  if(nv.ge.4)then
                    if (X(i,3).eq.1.and.X(i,4).eq.1.d0) then
                       claFCpl(g,1)=claFCpl(g,1)+1
                    endif
                    if (X(i,3).eq.1.and.X(i,4).eq.0.d0) then
                       claMCpl(g,1)=claMCpl(g,1)+1
                    endif
                    if (X(i,3).eq.0.and.X(i,4).eq.1.d0) then
                       claFCmo(g,1)=claFCmo(g,1)+1
                    endif
                    if (X(i,3).eq.0.and.X(i,4).eq.0.d0) then
                       claMCmo(g,1)=claMCmo(g,1)+1
                    endif
		  end if
               end if
               if(classe(i).eq.g.and.idm(i).eq.0.and.idd(i).eq.1) then
                  classement(g,2)=classement(g,2)+1
		if(nv.ge.4)then
                  if (X(i,3).eq.1.and.X(i,4).eq.1.d0) then
                     claFCpl(g,2)=claFCpl(g,2)+1
                  endif
                  if (X(i,3).eq.1.and.X(i,4).eq.0.d0) then
                     claMCpl(g,2)=claMCpl(g,2)+1
                  endif
                  if (X(i,3).eq.0.and.X(i,4).eq.1.d0) then
                     claFCmo(g,2)=claFCmo(g,2)+1
                  endif
                  if (X(i,3).eq.0.and.X(i,4).eq.0.d0) then
                     claMCmo(g,2)=claMCmo(g,2)+1
                  endif
		end if
               end if
               if(classe(i).eq.g.and.idm(i).eq.1.and.idd(i).eq.1) then
                  classement(g,3)=classement(g,3)+1
		if(nv.ge.4)then
                  if (X(i,3).eq.1.and.X(i,4).eq.1.d0) then
                     claFCpl(g,3)=claFCpl(g,3)+1
                  endif
                  if (X(i,3).eq.1.and.X(i,4).eq.0.d0) then
                     claMCpl(g,3)=claMCpl(g,3)+1
                  endif
                  if (X(i,3).eq.0.and.X(i,4).eq.1.d0) then
                     claFCmo(g,3)=claFCmo(g,3)+1
                  endif
                  if (X(i,3).eq.0.and.X(i,4).eq.0.d0) then
                     claMCmo(g,3)=claMCmo(g,3)+1
                  endif
		end if
               end if
               if(classe(i).eq.g.and.idm(i).eq.0.and.idd(i).eq.0) then
                  classement(g,4)=classement(g,4)+1
		if(nv.ge.4)then
                  if (X(i,3).eq.1.and.X(i,4).eq.1.d0) then
                     claFCpl(g,4)=claFCpl(g,4)+1
                  endif
                  if (X(i,3).eq.1.and.X(i,4).eq.0.d0) then
                     claMCpl(g,4)=claMCpl(g,4)+1
                  endif
                  if (X(i,3).eq.0.and.X(i,4).eq.1.d0) then
                     claFCmo(g,4)=claFCmo(g,4)+1
                  endif
                  if (X(i,3).eq.0.and.X(i,4).eq.0.d0) then
                     claMCmo(g,4)=claMCmo(g,4)+1
		  end if
		end if
               end if
            end do
            do g=1,ng
               if(classe1(i).eq.g.and.idm(i).eq.1.and.idd(i).eq.0) then
                  classement1(g,1)=classement1(g,1)+1
               end if
               if(classe1(i).eq.g.and.idm(i).eq.0.and.idd(i).eq.1) then
                  classement1(g,2)=classement1(g,2)+1
               end if
               if(classe1(i).eq.g.and.idm(i).eq.1.and.idd(i).eq.1) then
                  classement1(g,3)=classement1(g,3)+1
               end if
               if(classe1(i).eq.g.and.idm(i).eq.0.and.idd(i).eq.0) then
                  classement1(g,4)=classement1(g,4)+1
               end if
            end do
         end do
       
         do i=1,ng
            do g=1,ng
               tfinmoy(i,g)=sumvect(i,g)/pipost(i)
            enddo
         enddo
         
         write(8,*)"Mean probabilities by mixture component"
         write(8,*)"   ",(g,g=1,ng)
!          write(8,*)tfinmoy
!          write(8,*)sumvect
!          write(8,*)pitest
         do g=1,ng
            write(8,*)g,tfinmoy(g,:)
         enddo
         

!tres detaille pour application PAQUID SmoothHazard, pas forcement utile car depend des vexp

         write(8,*)'Number of demented, demented and dead, dead, healthy subjects'            
         write(8,*)(classement(g,1),g=1,ng)
         write(8,*)(classement(g,2),g=1,ng)
         write(8,*)(classement(g,3),g=1,ng)
         write(8,*)(classement(g,4),g=1,ng)

         write(8,*)"   "
         write(8,*)'Number of demented, demented and dead, dead, healthy subjects (Sexe=F, CEP=+)'          
              
         write(8,*)(claFCpl(g,1),g=1,ng)
         write(8,*)(claFCpl(g,2),g=1,ng)
         write(8,*)(claFCpl(g,3),g=1,ng)
         write(8,*)(claFCpl(g,4),g=1,ng)


         write(8,*)"   "
         write(8,*)'Number of demented, demented and dead, dead, healthy subjects (Sexe=F, CEP=-)'
         write(8,*)(claFCmo(g,1),g=1,ng)
         write(8,*)(claFCmo(g,2),g=1,ng)
         write(8,*)(claFCmo(g,3),g=1,ng)
         write(8,*)(claFCmo(g,4),g=1,ng)

         write(8,*)"   "
         write(8,*)'Number of demented, demented and dead, dead, healthy subjects (Sexe=M, CEP=+)'
         write(8,*)(claMCpl(g,1),g=1,ng)
         write(8,*)(claMCpl(g,2),g=1,ng)
         write(8,*)(claMCpl(g,3),g=1,ng)
         write(8,*)(claMCpl(g,4),g=1,ng)


         write(8,*)"   "
         write(8,*)'Number of demented, demented and dead, dead, healthy subjects (Sexe=M, CEP=-)'
         write(8,*)(claMCmo(g,1),g=1,ng)
         write(8,*)(claMCmo(g,2),g=1,ng)
         write(8,*)(claMCmo(g,3),g=1,ng)
         write(8,*)(claMCmo(g,4),g=1,ng)


         close(7)
         write(8,*)"   "
         sumvect=0.d0
         tfinmoy=0.d0         

         open(7,file='posterior_class_probs_from_Y_only.txt')
         do i=1,ns
            write(7,1998)numero(i),classe1(i),(PPItest(i,g),g=1,ng)
            do g=1,ng
               sumvect(classe1(i),g)=sumvect(classe1(i),g)+PPItest(i,g)
            enddo
         end do
         do i=1,ng
            do g=1,ng
               tfinmoy(i,g)=sumvect(i,g)/pitest(i)
            enddo
         enddo
         
	 write(8,*) 'POSTERIOR CLASSIFICATION  (P(c(i)=g|Y_i))'
	 write(8,*)
         write(8,*) 'Number of subjects per mixture component'
         write(8,*) (pitest(g),g=1,ng)


         write(8,*)"Mean probabilities by mixture component"
         write(8,*)"   ",(g,g=1,ng)
         do g=1,ng
            write(8,*)g,tfinmoy(g,:)
         enddo

         if(CurrentProcessorID.eq.0) write(*,*)'Number of demented, demented and dead, dead, healthy subjects'
         if(CurrentProcessorID.eq.0) write(*,*)(classement1(g,1),g=1,ng)
         if(CurrentProcessorID.eq.0) write(*,*)(classement1(g,2),g=1,ng)
         if(CurrentProcessorID.eq.0) write(*,*)(classement1(g,3),g=1,ng)
         if(CurrentProcessorID.eq.0) write(*,*)(classement1(g,4),g=1,ng)

         if(CurrentProcessorID.eq.0) write(*,*)'PRIOR CLASSIFICATION'
         if(CurrentProcessorID.eq.0) write(*,*)'Number of demented, demented and dead, dead, healthy subjects'
         if(CurrentProcessorID.eq.0) write(*,*)(classement(g,1),g=1,ng)
         if(CurrentProcessorID.eq.0) write(*,*)(classement(g,2),g=1,ng)
         if(CurrentProcessorID.eq.0) write(*,*)(classement(g,3),g=1,ng)
         if(CurrentProcessorID.eq.0) write(*,*)(classement(g,4),g=1,ng)
         if(CurrentProcessorID.eq.0) write(*,*)

         write(8,*)'Number of demented, demented and dead, dead, healthy subjects '       
         write(8,*)(classement1(g,1),g=1,ng)
         write(8,*)(classement1(g,2),g=1,ng)
         write(8,*)(classement1(g,3),g=1,ng)
         write(8,*)(classement1(g,4),g=1,ng)
         write(8,*)


         close(7)
         write(8,*)"   "
         write(8,*)"Prior probabilities"
         write(8,*)(g,g=1,ng)         
         write(8,*)(pipre(g),g=1,ng)
         deallocate(classe,classe1,classement,classement1,pipost,pitest)
        
               if(CurrentProcessorID.eq.0) write(*,*)
         if(CurrentProcessorID.eq.0) write(*,*)'Class-membership probabilities saved in:'
         if(CurrentProcessorID.eq.0) write(*,*)'   prior_class_probs.txt, ', &
              'posterior_class_probs.txt ','and posterior_class_probs_from_Y_only.txt'
     

      end if

      i=1
      if(i.eq.1)then

!      if(CurrentProcessorID.eq.0) write(*,*)      
!      if(CurrentProcessorID.eq.0) write(*,*)'Computation of predictions may take time - be patient...'
!      allocate(ymarg(ns*nq*mesure,ng))
      allocate(mui1(ns*nq*mesure,ng),muicond1(ns*nq*mesure,ng),ymarg1(ns*nq*mesure,ng),ycond1(ns*nq*mesure,ng))

!      ymarg=0.d0
      mui1=0.d0
      ymarg1=0.d0
      ycond1=0.d0
      muicond1=0.d0
      !C           if(CurrentProcessorID.eq.0) write(*,*)'avant evol_pred'
      !C prediction marginale

      call evol_pred_marg_sim(b,npm,nsim,ymarg1,ycond1,mui1,muicond1)

      !C meme routine mais calcul direct pour outcomes ordinaux
!      if(CurrentProcessorID.eq.0) write(*,*)'apres evol_pred'
!      call evol_pred_marg_ordi(b,npm,nsim,ymarg,mui1)
      
      
      !if (quecont.eq.1.and.ng.eq.1) then
         
      !   allocate(alpha_i(ns,nq),u_i(ns,nea))
      !   allocate(pred_ind(ns,maxmestot))
      !   allocate(Y_trans(ns,nq,mesure),mui(ns,nq,mesure))

      !   u_i=0.d0
      !   alpha_i=0.d0
      !   pred_ind=0.d0
      !   Y_trans=0.d0
      !   mui=0.d0
         
         !C 2586         format(I5,1X,I3,1X,F8.5,1X,F8.5,1X,F8.5,1X,F8.5)
       !  call pred_indiv(b,npm,pred_ind,Y_trans,mui,u_i,alpha_i)
         
         !C              open(7,file='pred_latent.txt')
         !do i=1,ns
         !   jj=0
         !   do q=1,nq
         !      do j=1,nmes(i,q)
         !         jj=jj+1
                  !C                       write(7,2586)numero(i),q,t(i,q,j),Y_trans(i,q,j),
                  !c     &                      pred_ind(i,jj),mui(i,q,j)
         !      end do
         !   end do
         !end do
         !C              close(7)
         

        ! deallocate(alpha_i,u_i,pred_ind,Y_trans,mui)

      !end if
      
      if (ng.eq.1) then
         ppi=1.d0
      end if

      allocate(som(ng,nb_cl_tps,nq),pond(ng,nb_cl_tps,nq),stderr(ng,nb_cl_tps,nq),temp_std(ng,nb_cl_tps,nq) &
           ,som_obs(ng,nb_cl_tps,nq),som_ui(ng,nb_cl_tps,nq),som_ui_obs(ng,nb_cl_tps,nq),IC_inf(ng,nb_cl_tps,nq))
	allocate(IC_sup(ng,nb_cl_tps,nq),pond2(ng,nb_cl_tps,nq),ns_cl_tps(nb_cl_tps,nq),time_cl_tps(nb_cl_tps,nq))



2001  FORMAT(I4,1X,I4,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X, &
             F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6, &
             1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6, &
             1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6 &
             ,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X &
            ,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6, &
             1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6 &
             ,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6)


      
      open(9,file='mean_predicted_trajectories.txt')
      stderr=0.d0
      temp_std=0.d0
      som=0.d0
      som_obs=0.d0
      pond=0.d0
      pond2=0.d0
      ns_cl_tps=0
      time_cl_tps=0.d0
      IC_inf=0.d0
      IC_sup=0.d0
      som_ui=0.d0
      som_ui_obs=0.d0
	  
	allocate(som_marg(nb_cl_tps,nq),som_marg_ui(nb_cl_tps,nq),som_obs_yi(nb_cl_tps,nq))
	allocate(IC_sup_yi(nb_cl_tps,nq),IC_inf_yi(nb_cl_tps,nq),stderr_yi(nb_cl_tps,nq),temp_std_yi(nb_cl_tps,nq))
      som_marg=0.d0  
      temp_std_yi=0.d0
      som_marg_ui=0.d0
      som_obs_yi=0.d0   
      IC_inf_yi=0.d0
      IC_sup_yi=0.d0
      do l=1,nb_cl_tps-2
         nmestotcur=0
         do q=1,nq
            do i=1,ns
               do j=1,nmes(i,q)
                  if (l.eq.1.and.t(Navant_i(i)+nmestotcur(i)+j).lt.cl_tps(l)) then
                     ns_cl_tps(l,q)=ns_cl_tps(l,q)+1
                     time_cl_tps(l,q)=time_cl_tps(l,q)+t(Navant_i(i)+nmestotcur(i)+j)
                     do g=1,ng
                        som(g,l,q)=som(g,l,q)+ppi(i,g)*ymarg1(Navant_i(i)+nmestotcur(i)+j,g)
                        som_obs(g,l,q)=som_obs(g,l,q)+ppi(i,g)*y(Navant_i(i)+nmestotcur(i)+j)
			som_ui(g,l,q)=som_ui(g,l,q)+ppi(i,g)*ycond1(Navant_i(i)+nmestotcur(i)+j,g)
			som_marg(l,q)=som_marg(l,q)+pig(i,g)*ymarg1(Navant_i(i)+nmestotcur(i)+j,g)
			som_marg_ui(l,q)=som_marg_ui(l,q)+pig(i,g)*ycond1(Navant_i(i)+nmestotcur(i)+j,g)

                        pond(g,l,q)=pond(g,l,q)+ppi(i,g)
                        pond2(g,l,q)=pond2(g,l,q)+ppi(i,g)*ppi(i,g)
                        stderr(g,l,q)=stderr(g,l,q)+y(Navant_i(i)+nmestotcur(i)+j)*y(Navant_i(i)+nmestotcur(i)+j)
                        temp_std(g,l,q)=temp_std(g,l,q)+y(Navant_i(i)+nmestotcur(i)+j)
                     end do
		        stderr_yi(l,q)=stderr_yi(l,q)+y(Navant_i(i)+nmestotcur(i)+j)*y(Navant_i(i)+nmestotcur(i)+j)
                        temp_std_yi(l,q)=temp_std_yi(l,q)+y(Navant_i(i)+nmestotcur(i)+j)
                        som_obs_yi(l,q)=som_obs_yi(l,q)+y(Navant_i(i)+nmestotcur(i)+j)
                  else if (t(Navant_i(i)+nmestotcur(i)+j).ge.cl_tps(l).and.t(Navant_i(i)+nmestotcur(i)+j).lt.cl_tps(l+1)) then
                     ns_cl_tps(l+1,q)=ns_cl_tps(l+1,q)+1
                     time_cl_tps(l+1,q)=time_cl_tps(l+1,q)+t(Navant_i(i)+nmestotcur(i)+j)
                     do g=1,ng
                        som(g,l+1,q)=som(g,l+1,q)+ppi(i,g)*ymarg1(Navant_i(i)+nmestotcur(i)+j,g)
                        som_obs(g,l+1,q)=som_obs(g,l+1,q)+ppi(i,g)*y(Navant_i(i)+nmestotcur(i)+j)
                        som_ui(g,l+1,q)=som_ui(g,l+1,q)+ppi(i,g)*ycond1(Navant_i(i)+nmestotcur(i)+j,g)
			som_marg(l+1,q)=som_marg(l+1,q)+pig(i,g)*ymarg1(Navant_i(i)+nmestotcur(i)+j,g)
			som_marg_ui(l+1,q)=som_marg_ui(l+1,q)+pig(i,g)*ycond1(Navant_i(i)+nmestotcur(i)+j,g)

                        pond(g,l+1,q)=pond(g,l+1,q)+ppi(i,g)
                        pond2(g,l+1,q)=pond2(g,l+1,q)+ppi(i,g)*ppi(i,g)
                        stderr(g,l+1,q)=stderr(g,l+1,q)+y(Navant_i(i)+nmestotcur(i)+j)*y(Navant_i(i)+nmestotcur(i)+j)
                        temp_std(g,l+1,q)=temp_std(g,l+1,q)+y(Navant_i(i)+nmestotcur(i)+j)
                     end do
			stderr_yi(l+1,q)=stderr_yi(l+1,q)+y(Navant_i(i)+nmestotcur(i)+j)*y(Navant_i(i)+nmestotcur(i)+j)
                        temp_std_yi(l+1,q)=temp_std_yi(l+1,q)+y(Navant_i(i)+nmestotcur(i)+j)
                        som_obs_yi(l+1,q)=som_obs_yi(l+1,q)+y(Navant_i(i)+nmestotcur(i)+j)
                  else if (l.eq.nb_cl_tps-2.and.t(Navant_i(i)+nmestotcur(i)+j).ge.cl_tps(l+1)) then
                     ns_cl_tps(l+2,q)=ns_cl_tps(l+2,q)+1
                     time_cl_tps(l+2,q)=time_cl_tps(l+2,q)+t(Navant_i(i)+nmestotcur(i)+j)
                     do g=1,ng
                        som(g,l+2,q)=som(g,l+2,q)+ppi(i,g)*ymarg1(Navant_i(i)+nmestotcur(i)+j,g)
                        som_obs(g,l+2,q)=som_obs(g,l+2,q)+ppi(i,g)*y(Navant_i(i)+nmestotcur(i)+j)
                        som_ui(g,l+2,q)=som_ui(g,l+2,q)+ppi(i,g)*ycond1(Navant_i(i)+nmestotcur(i)+j,g)
			som_marg(l+2,q)=som_marg(l+2,q)+pig(i,g)*ymarg1(Navant_i(i)+nmestotcur(i)+j,g)
			som_marg_ui(l+2,q)=som_marg_ui(l+2,q)+pig(i,g)*ycond1(Navant_i(i)+nmestotcur(i)+j,g)

                        pond(g,l+2,q)=pond(g,l+2,q)+ppi(i,g)
                        pond2(g,l+2,q)=pond2(g,l+2,q)+ppi(i,g)*ppi(i,g)
                        stderr(g,l+2,q)=stderr(g,l+2,q)+y(Navant_i(i)+nmestotcur(i)+j)*y(Navant_i(i)+nmestotcur(i)+j)
                        temp_std(g,l+2,q)=temp_std(g,l+2,q)+y(Navant_i(i)+nmestotcur(i)+j)
                     end do
			stderr_yi(l+2,q)=stderr_yi(l+2,q)+y(Navant_i(i)+nmestotcur(i)+j)*y(Navant_i(i)+nmestotcur(i)+j)
                        temp_std_yi(l+2,q)=temp_std_yi(l+2,q)+y(Navant_i(i)+nmestotcur(i)+j)
                        som_obs_yi(l+2,q)=som_obs_yi(l+2,q)+y(Navant_i(i)+nmestotcur(i)+j)                     
                  end if
               end do
               nmestotcur(i)=nmestotcur(i)+nmes(i,q)
            end do
         end do
      end do

      do q=1,nq
         do l=1,nb_cl_tps
            do g=1,ng
               som(g,l,q)=som(g,l,q)/pond(g,l,q)
               som_obs(g,l,q)=som_obs(g,l,q)/pond(g,l,q)
               som_ui(g,l,q)=som_ui(g,l,q)/pond(g,l,q)
               stderr(g,l,q)=sqrt(pond2(g,l,q)*((stderr(g,l,q) -(temp_std(g,l,q) &
                    *temp_std(g,l,q))/dble(ns_cl_tps(l,q)))/ &
                    dble(ns_cl_tps(l,q)))/(pond(g,l,q)*pond(g,l,q)))
               IC_inf(g,l,q)=som_obs(g,l,q)-1.96*stderr(g,l,q)
               IC_sup(g,l,q)=som_obs(g,l,q)+1.96*stderr(g,l,q)
            end do
	       stderr_yi(l,q)=sqrt(stderr_yi(l,q)/dble(ns_cl_tps(l,q))-(temp_std_yi(l,q) &
                    *temp_std_yi(l,q))/(dble(ns_cl_tps(l,q))*dble(ns_cl_tps(l,q))))
	       som_marg(l,q)=som_marg(l,q)/dble(ns_cl_tps(l,q))
	       som_marg_ui(l,q)=som_marg_ui(l,q)/dble(ns_cl_tps(l,q))
	       som_obs_yi(l,q)=som_obs_yi(l,q)/dble(ns_cl_tps(l,q))
               IC_inf_yi(l,q)=som_obs_yi(l,q)-1.96*stderr_yi(l,q)
               IC_sup_yi(l,q)=som_obs_yi(l,q)+1.96*stderr_yi(l,q)
            write(9,2001)q,ns_cl_tps(l,q),time_cl_tps(l,q)/dble(ns_cl_tps(l,q)), &
                 (som(g,l,q),g=1,ng),(som_ui(g,l,q),g=1,ng),(som_obs(g,l,q),g=1,ng) &
              ,(stderr(g,l,q),g=1,ng),(IC_inf(g,l,q),g=1,ng) &
              ,(IC_sup(g,l,q),g=1,ng),som_marg(l,q),som_marg_ui(l,q),som_obs_yi(l,q)&
		,stderr_yi(l,q),IC_inf_yi(l,q),IC_sup_yi(l,q)
         end do
      end do

      if(CurrentProcessorID.eq.0) write(*,*)'Weighted mean of class-specific predicted outcome trajectories ', &
           'saved in:'
         if(CurrentProcessorID.eq.0) write(*,*)'    mean_predicted_trajectories.txt'
      close(9)
         

      deallocate(som,pond,stderr,temp_std,som_obs,IC_inf,IC_sup,pond2,ns_cl_tps,time_cl_tps)
      deallocate(IC_sup_yi,IC_inf_yi,som_obs_yi,som_marg_ui,som_marg)
!      deallocate(ymarg)
       deallocate(mui1,ymarg1)


      open(9,file='transformations.txt')
      
7800  format(I3,1X,F12.8,1X,F12.8)
      
      
      allocate(test(nsim),transf(nsim))
      write(9,*)'outcome  observed  tranformed'
      ntrtot=0
      do q=1,nq
         if (idord(q).eq.1) then

            allocate(splaa(-2:nztr),splaa2(nztr+3),zitr2(-1:nztr+2))


            k=0
            l=0
            splaa2=0.d0
            DO j=1,ntrspl
               splaa2(j)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &
                    +ncontdep*(nq-1)+ncont*(nq-1)+j+ntrtot)
            END DO
            
            
            
            zitr2=0.d0
            do j=-1,nztr+2
               zitr2(j)=zitr(j,q)
            end do
            
            !C               if(CurrentProcessorID.eq.0) write(*,*)
            !C               if(CurrentProcessorID.eq.0) write(*,*)'splaa',(splaa2(j),j=1,ntrspl)
            !C               if(CurrentProcessorID.eq.0) write(*,*)'splaa complet',splaa2
            


            
            test=0.d0
            transf=0.d0
            call   estim_splines_ssstd(nsim,nztr,zitr2,splaa2,test,transf)
            
            
            do j=1,nsim
               write(9,7800)q,test(j),transf(j)
            end do
            ntrtot=ntrtot+ntrspl
            
            deallocate(zitr2)
            deallocate(splaa,splaa2)
            
         else if (idord(q).le.0) then
            
            !C ------------ calcul transfo pour Beta ---------------------
            
            aa1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                 +ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot))/(1+exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1) &
                 +nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot)))
            bb1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &    
                 +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &    
                 +ncont*(nq-1)     &         
                 +2+ntrtot))/( 1.D0+ exp(b1(nvarprob*(ng-1)+nrisq     &   
                 +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+2+ntrtot)))
            bb1=aa1*(1.d0-aa1)*bb1
            cc1=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &   
                 +ncont*(nq-1)+3+ntrtot)
            dd1=abs(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1) &
                 +ncont*(nq-1)+4+ntrtot))
            aa=aa1*aa1*(1-aa1)/bb1-aa1
            bb=aa*(1-aa1)/aa1
            
            test=0.d0
            transf=0.d0
            
            pas=(rangereal((q-1)*2+2)-rangereal((q-1)*2+1))/dble(nsim-1)
            j=1
!            if(CurrentProcessorID.eq.0) write(*,*)'pas=',pas
            test(1)=rangereal((q-1)*2+1)
            do while(j.lt.nsim)
               j=j+1
               test(j)=test(j-1)+pas
            end do
            test(nsim)=rangereal((q-1)*2+2)
            
            do j=1,nsim
               ytemp=(test(j)-range(1+(q-1)*2))/(range(2+(q-1)*2)-range(1+(q-1)*2))
               transf(j)=(betai(aa,bb,ytemp,tbeta)-cc1)/dd1
               write(9,7800)q,test(j),transf(j)
            end do
            
            ntrtot=ntrtot+ntr
            
         else if (idord(q).gt.1) then
            
            sup=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                 +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
            inf=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                 +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
            
            
            write(9,7800)q,0.d0+range(1+(q-1)*2),inf
            write(9,7800)q,0.d0+range(1+(q-1)*2),sup
            
            

            do j=1,idord(q)-2
               
               
               sup=sup+(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                    +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+j+1)*b1(nvarprob*(ng-1)+nrisq &
                    +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+j+1))
               
               write(9,7800)q,dble(j+range(1+(q-1)*2)),inf
               write(9,7800)q,dble(j+range(1+(q-1)*2)),sup
               
               inf=sup
            end do
            write(9,7800)q,dble(idord(q)-1+range(1+(q-1)*2)),inf
            
            ntrtot=ntrtot+idord(q)-1
            
         end if
         
      end do


      deallocate(test,transf)
      
      close(9)
      if(CurrentProcessorID.eq.0) write(*,*)'Estimated transformations saved in:'
      if(CurrentProcessorID.eq.0) write(*,*)'   transformations.txt'
    

!Postfit suvie : routine dans programme a part
  
       if (evt.ne.0) then
         
         
         
! 5643     format(F9.4,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X, &
!               F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X, &
!               F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5,1X,F9.5)
                 
         open(78,file='survival.txt')

         allocate(gaussGL(2,npg))
         
         call gaussleg(gaussGL,npg)
         nrisq_curr=0
         nsim2=nsim*npg
         res_tot=0.d0
     
         allocate(risq_est(nevt,nsim*ng),risqcum_est(nevt,nsim*ng),brisq_est(nrisq)&
              ,incidcum_est(nevt,nsim*ng),time_est2(nsim2),risqspec(nevt,nsim2*ng), &
              risqcumspec(nevt,nsim2*ng),brisqtot(contevt),brisq(contevt2),time_est01(nsim),&
		time_est02(nsim),time_est12(nsim),bsex(nevt))  
  
	bsex=0.d0
	do i=1,nevt
		bsex(i)=0!b(ng-1+nhazard+(i-1)*(nv-2)+nv-2-1) ! si le sex est la derniere variable explicative
	end do
	sex=0 ! a changer si on veut calculer pour les femmes
         if(CurrentProcessorID.eq.0) write(*,*)'maxval(typrisq)',maxval(typrisq)
		 
		 if (typrisq(1).eq.1)then
			 stept=(100.d0-65.d0)/(nsim-1) !specifique application : moche!
			 time_est01(1)=65.d0
			 do i=2,nsim
				time_est01(i)=time_est01(i-1)+stept
			 enddo
		 else
		    stept=(100.d0-65.d0)/(nsim-1) !spécifique application : moche!
			time_est01(1)=65.d0
			 do i=2,nsim
				time_est01(i)=time_est01(i-1)+stept
			 enddo
		 end if
		 

		 if(nevt.ge.2)then
		 	if (typrisq(2).eq.1)then
				 stept=(100.d0-65.d0)/(nsim-1) !specifique application : moche!
				 time_est02(1)=65.d0
				 do i=2,nsim
					time_est02(i)=time_est02(i-1)+stept
				 enddo
			else
				stept=(100.d0-65.d0)/(nsim-1) !spécifique application : moche!
				time_est02(1)=65.d0
				 do i=2,nsim
					time_est02(i)=time_est02(i-1)+stept
				 enddo
			 end if
		end if
	
		
	 if(nevt.ge.3)then
	 	if (typrisq(3).eq.1)then
				 stept=(100.d0-65.d0)/(nsim-1) !specifique application : moche!
				 time_est12(1)=65.d0
				 do i=2,nsim
					time_est12(i)=time_est12(i-1)+stept
				 enddo
                else
		    stept=(100.d0-65.d0)/(nsim-1) !spécifique application : moche!
			time_est12(1)=65.d0
			 do i=2,nsim
				time_est12(i)=time_est12(i-1)+stept
			 enddo
		 end if
	end if
		
		
		! where is time_est2 used ?
         time_est2=0.d0
         jj=0
         do iGl=1,npg
            do i=1,nsim
               jj=jj+1
               time_est2(jj)=time_est(i)/2.d0+time_est(i)/2.d0 *gaussGL(1,iGL) 
!               if(CurrentProcessorID.eq.0) write(*,*)i,igl,jj,time_est2(jj),time_est(i),time_est(i)/2.d0
            end Do
         end DO


         risq_est=0.d0
         risqcum_est=0.d0
         risqspec=0.d0
         risqcumspec=0.d0
         
         brisq=0.d0

         if (evt.ne.0) then
            nrisq_curr=0
            numact2=1
            !numact permet le remplissage de brisq
            !on boucle sur le nombre d'evt
            do kk=1,nevt
               !boucle sur le nombre de classes
               if (risqcom(kk).eq.0) then
                  do g=1,ng
                     !boucle sur le nombre de parametres par evt
                     do k=1,nprisq(kk)
                        !positivation carre ou expo
                        if (paramsurv.eq.0) then 
                           !risque commun/specifique/proportionnel
                           if (risqcom(kk).eq.2) then
                              brisq(numact2)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)    
                              numact2=numact2+1
                           end if
                           if (risqcom(kk).eq.1) then
                              brisq(numact2)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)
                              numact2=numact2+1
                           end if
                           if (risqcom(kk).eq.0) then
                              brisq(numact2)=b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)&
                                   &+nrisq_curr+k)  !commente le carre pour tester!!
                              numact2=numact2+1
                           end if
                        else 
                           if (risqcom(kk).eq.2) then
                              brisq(numact2)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))    
                              numact2=numact2+1
                           end if
                           if (risqcom(kk).eq.1) then
                              brisq(numact2)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))
                              numact2=numact2+1
                           end if
                           if (risqcom(kk).eq.0) then
                              brisq(numact2)=exp(b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k))
                              numact2=numact2+1
                           end if
                        end IF !fin paramsurv==0
                     end do !fin boucle nb param par evt
                  end do !boucle sur les classes
               endif
               if (risqcom(kk).ge.1) then
               
                  do k=1,nprisq(kk)
                    !positivation carre ou expo
                     if (paramsurv.eq.0) then 
                       !risque commun/specifique/proportionnel
                        if (risqcom(kk).eq.2) then
                           brisq(numact2)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)    
                           numact2=numact2+1
                        end if
                        if (risqcom(kk).eq.1) then
                           brisq(numact2)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)
                           numact2=numact2+1
                        end if
                        if (risqcom(kk).eq.0) then
                           brisq(numact2)=b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)&
                                                    &+nrisq_curr+k)
                           numact2=numact2+1
                        end if
                     else 
                        if (risqcom(kk).eq.2) then
                           brisq(numact2)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))    
                           numact2=numact2+1
                        end if
                        if (risqcom(kk).eq.1) then
                           brisq(numact2)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))
                           numact2=numact2+1
                        end if
                        if (risqcom(kk).eq.0) then
                           brisq(numact2)=exp(b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k))
                           numact2=numact2+1
                        end if
                     end IF !fin paramsurv==0
                  end do !fin boucle nb param par evt
               endif
               nrisq_curr=nrisq_curr+nrisqspec(kk)
            end do
            nrisq_curr=0
            do kk=1,nevt
               if (risqcom(kk).eq.2) then
                  nrisq_curr=nrisq_curr+nrisqspec(kk)-ng+1
                  do g=1,ng-1
                     brisq(numact2)=b1(nvarprob*(ng-1)+nrisq_curr+g)
                     numact2=numact2+1
                  enddo
                  nrisq_curr=nrisq_curr+ng-1
               else
                  nrisq_curr=nrisq_curr+nrisqspec(kk)
               endif
            enddo
            
         endif

kk=1

            numact2=1 !progression remplissage vecteur           
            do i=1,nsim
               do g=1,ng
                  call qg_risqcum(time_est01(1),kk,g,time_est01(i),brisq,res_inst,res_cum,val_func,sex,bsex) !calculs risques
                  call qg_inccum(time_est01(1),time_est01(i),brisq,g,kk,gaussgl,res_tot,sex,bsex) !quadrature pour incidence
                  risq_est(kk,numact2)=res_inst
                  risqcum_est(kk,numact2)=res_cum !vecteurs risques/irsques cumules et incidences cumules
                  incidcum_est(kk,numact2)=res_tot
                  numact2=numact2+1
               enddo
			   if (nevt.gt.1) then
					if(CurrentProcessorID.eq.0) write(78,*)kk,time_est01(i),(risq_est(kk,(i-1)*ng+g),g=1,ng),    &
                       (risqcum_est(kk,(i-1)*ng+g),g=1,ng),(incidcum_est(kk,(i-1)*ng+g),g=1,ng)
				 ELSE
					     if(CurrentProcessorID.eq.0) write(78,*)kk,time_est(i),(risq_est(kk,(i-1)*ng+g),g=1,ng),    &
						(risqcum_est(kk,(i-1)*ng+g),g=1,ng),(incidcum_est(kk,(i-1)*ng+g),g=1,ng), &
						(1.d0-exp(-risqcum_est(kk,(i-1)*ng+g)),g=1,ng)
				end if
            enddo

	 if (nevt.ge.2) then				   

	 kk=2
            numact2=1 !progression remplissage vecteur           
            do i=1,nsim
               do g=1,ng
                  call qg_risqcum(time_est02(1),kk,g,time_est02(i),brisq,res_inst,res_cum,val_func,sex,bsex) !calculs risques
                  call qg_inccum(time_est02(1),time_est02(i),brisq,g,kk,gaussgl,res_tot,sex,bsex) !quadrature pour incidence
                  risq_est(kk,numact2)=res_inst
                  risqcum_est(kk,numact2)=res_cum !vecteurs risques/irsques cumules et incidences cumules
                  incidcum_est(kk,numact2)=res_tot
                  numact2=numact2+1
               enddo
			    if(CurrentProcessorID.eq.0) write(78,*)kk,time_est02(i),(risq_est(kk,(i-1)*ng+g),g=1,ng),    &
                       (risqcum_est(kk,(i-1)*ng+g),g=1,ng),(incidcum_est(kk,(i-1)*ng+g),g=1,ng)
            enddo
	end if
	
	 if (nevt.ge.3) then	 
	kk=3
            numact2=1 !progression remplissage vecteur           
            do i=1,nsim
               do g=1,ng
                  call qg_risqcum(time_est12(1),kk,g,time_est12(i),brisq,res_inst,res_cum,val_func,sex,bsex) !calculs risques
                  call qg_inccum(time_est12(1),time_est12(i),brisq,g,kk,gaussgl,res_tot,sex,bsex) !quadrature pour incidence
                  risq_est(kk,numact2)=res_inst
                  risqcum_est(kk,numact2)=res_cum !vecteurs risques/irsques cumules et incidences cumules
                  incidcum_est(kk,numact2)=res_tot
                  numact2=numact2+1
               enddo
			if(CurrentProcessorID.eq.0) write(78,*)kk,time_est12(i),(risq_est(kk,(i-1)*ng+g),g=1,ng),    &
                       (risqcum_est(kk,(i-1)*ng+g),g=1,ng),(incidcum_est(kk,(i-1)*ng+g),g=1,ng)
		enddo
	
	end if	   
            
         deallocate(risq_est,risqcum_est,risqspec,risqcumspec,brisq_est,time_est2,time_est01,time_est02,time_est12)
         deallocate(gaussGL)
 

        

         close(78)
         if(CurrentProcessorID.eq.0) write(*,*)'Estimated baseline risk and survival functions saved in:'
         if(CurrentProcessorID.eq.0) write(*,*)'    survival.txt'

      end if
      end if

! Postfit multi-state to compare with smoothHazard
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	call postfit(brisq,PPItest,ppi,t0,t1,t2,t3,sex,bsex)
	if(ng.gt.1) call descript_class(ppi)      
      
!============= desallouer 


      deallocate(valinit,ide,bne,bsex,pipre)

      deallocate(time_est)
      deallocate(ppi,ppitest,P_AI,pi)
      deallocate(B)
      deallocate(cl_tps)
      deallocate(numero,ind_survint,devt,Tpred,tsurv,tsurv0,tsurvint,X,Y,t)
      deallocate(idglob,idg,idtest,idprob,idevt,idxevt,idvdep,Xdep,zitr,zi01,zi02,zi12)
      deallocate(vopt,b1,Navant_i,nmestotcur)
     


      deallocate(Ut1,mu,Y1,Z1,sigmaq,Vw)      
      deallocate(idord,idea,range,rangereal)
      if (spline.eq.1) then 
         deallocate(mm,mm1,mm2,im,im1,im2)
      end if
      if (nevt.ne.0) then ! quelle fonction de base ? risqcom, noeuds et ecart
         deallocate(typrisq,risqcom,nprisq,nrisqspec,nxevtspec,idzi,nz)
      end if

!C----------------- Duree de l'optimisation -------------------

      call cpu_time(time_fin)
      time=time_fin-time_deb
      call time_hms(time,heure,minute,seconde)
      if(CurrentProcessorID.eq.0) write(*,*)
      if(CurrentProcessorID.eq.0) write(*,*)'CPU time :'
      if(CurrentProcessorID.eq.0) write(*,*)heure,'h ',minute,'m ',seconde,'s'
      if(CurrentProcessorID.eq.0) write(*,*)

      write(8,*)
      write(8,*)'CPU time :'
      write(8,*)heure,'h ',minute,'m ',seconde,'s'

      close(8)

      call MPI_Finalize(info)

      end program hetmixbeta


!C---------------- FIN DU PROGRAMME PRINCIPAL -----------------

! }}}



!-------------------------------------------------------------------
!
!!  FCT_RISQ ESTIME
!
!-------------------------------------------------------------------

! {{{

!      subroutine fct_risq_estime(k,brisq,time,nsim,g,risq,surv)
!
!      use commun

!      implicit none
! risq contient le risque instantane d'evenement
! surv contient le risq cumule d'evenement et non la fonction de survie directe


!     double precision, dimension(nrisq)::brisq
!      integer ::i,g,l,j,nsim,kk,ll,ii,k
!      double precision,dimension(nsim*ng):: risq,surv
!      double precision::som
!      double precision,dimension(nsim)::time
      
      
     
!      l=0
!      do i=1,nsim
!         if (typrisq(k).eq.2) then
!            if (paramsurv.eq.0) then 
!               surv(nsim*(g-1)+i)=exp(((time(i))*brisq(2))**brisq(1))
!               risq(nsim*(g-1)+i)=(brisq(2)**brisq(1))*brisq(1)*(time(i))**(brisq(1)-1)
!            ELSE
!               surv(nsim*(g-1)+i)=exp(((time(i))/brisq(2))**brisq(1))
!               risq(nsim*(g-1)+i)=(brisq(1)/brisq(2))*(time(i)/brisq(2))**(brisq(1)-1)
!            end if
!         end if

!         if (typrisq(k).eq.1) then
!            do j=1,nz(k)-1
!               som=0.d0
!               do l=1,j-1
!                  som=som+brisq(l)*(zi(l+1)-zi(l))
!               end do
!               if (time(i).ge.zi(j).and.time(i).le.zi(j+1)) then
!                  surv(nsim*(g-1)+i)=som+brisq(j)*(time(i)-zi(j))
!                  risq(nsim*(g-1)+i)=brisq(j)
!               end if
!            end do
!         end if

!         if (typrisq(k).eq.3) then
            !------------ survie et risq pour Tsurv ----------------
!            ll=0
!            if (time(i).eq.zi(nz)) then
!               ll=nz-1
!            end if
!            som=0.d0
!            do kk=2,nz
!               if ((time(i).ge.zi(kk-1)).and.(time(i).lt.zi(kk))) &
!                    then
!                  ll=kk-1
!               end if
!            end do
!            if (ll.gt.1) then
!               do ii=1,ll-1
!                  som=som+brisq(ii)
!               end do
!            end if

!            surv(nsim*(g-1)+i)=som+brisq(ll)*Tim3_est(i)+brisq(ll+1)*Tim2_est(i) &
!                 +brisq(ll+2)*Tim1_est(i)+brisq(ll+3)*Tim_est(i)
!            risq(nsim*(g-1)+i)=brisq(ll)*Tmm3_est(i)+brisq(ll+1)*Tmm2_est(i)     &
!                 +brisq(ll+2)*Tmm1_est(i)+brisq(ll+3)*Tmm_est(i)

!         end if


!      end do


!      end subroutine fct_risq_estime

! }}}



!-------------------------------------------------------------------
!
!!  FCT_RISQ ESTIME  2
!
!-------------------------------------------------------------------

! {{{

!      subroutine fct_risq_estime2(k,brisq,time,nsim,g,risq,surv)

!      use commun

!      implicit none
! risq contient le risque instantane d'evenement
! surv contient le risq cumule d'evenement et non la fonction de survie directe


!      double precision, dimension(nrisq)::brisq
!      integer ::i,g,l,j,nsim,kk,ll,ii,k
!      double precision,dimension(nsim*ng):: risq,surv
!      double precision::som
!      double precision,dimension(nsim)::time


!      l=0
!      do i=1,nsim
!         if (typrisq(k).eq.2) then
!            if (paramsurv.eq.0) then 
!               surv(nsim*(g-1)+i)=exp(((time(i))*brisq(2))**brisq(1))
!               risq(nsim*(g-1)+i)=(brisq(2)**brisq(1))*brisq(1)*(time(i))**(brisq(1)-1)
!            ELSE
!               surv(nsim*(g-1)+i)=exp(((time(i))/brisq(2))**brisq(1))
!               risq(nsim*(g-1)+i)=(brisq(1)/brisq(2))*(time(i)/brisq(2))**(brisq(1)-1)
!            end if
!         end if

!         if (typrisq(k).eq.1) then
!            do j=1,nz-1
!               som=0.d0
!               do l=1,j-1
!                  som=som+brisq(l)*(zi(l+1)-zi(l))
!               end do
!               if (time(i).ge.zi(j).and.time(i).le.zi(j+1)) then
!                  surv(nsim*(g-1)+i)=som+brisq(j)*(time(i)-zi(j))
!                  risq(nsim*(g-1)+i)=brisq(j)
!               end if
!            end do
!         end if

!         if (typrisq(k).eq.3) then
            !------------ survie et risq pour Tsurv ----------------
!            ll=0
!            if (time(i).eq.zi(nz)) then
!               ll=nz-1
!            end if
!            som=0.d0
!            do kk=2,nz
!               if ((time(i).ge.zi(kk-1)).and.(time(i).lt.zi(kk))) &
!                    then
!                  ll=kk-1
!               end if
!            end do
!            if (ll.gt.1) then
!               do ii=1,ll-1
!                  som=som+brisq(ii)
!               end do
!            end if

!            surv(nsim*(g-1)+i)=som+brisq(ll)*Tim3_est2(i)+brisq(ll+1)*Tim2_est2(i) &
!                 +brisq(ll+2)*Tim1_est2(i)+brisq(ll+3)*Tim_est2(i)
!            risq(nsim*(g-1)+i)=brisq(ll)*Tmm3_est2(i)+brisq(ll+1)*Tmm2_est2(i)     &
!                 +brisq(ll+2)*Tmm1_est2(i)+brisq(ll+3)*Tmm_est2(i)

!         end if


!      end do



!      end subroutine fct_risq_estime2

! }}}



!C------------------------------------------------------------
!C                        FUNCPA
!C------------------------------------------------------------



! {{{

      double precision function funcpa(b,npm,id,thi,jd,thj,ni)

        use commun
        use donnees_indiv
        use optim_parallele

        implicit none
        integer::npm,ier,nnn,ntrtot,nrisq_curr,numact,numprob,ni
        integer :: id,jd,i,j,k,l,jj,kk,q,m,g,nmoins
        integer::nmestot,kkk,jjj
        integer ::jj1,k1,q1,h,ll,ii,nxevtcurr
        double precision::thi,thj,Y5,expo,eps,temp,vrais,det,aa,bb,varexpsurv
        double precision::jacobien,ytemp,betai,beta_densite,retard,entretard
        double precision::vrais_survie,aa1,bb1,cc1,dd1,fevt,vraisobs,som &
             ,surv_glob,surv0_glob,survint_glob,valfin,retard2,testsurv,scorefin
        double precision::vrais2,SUTi0,STi0

        !C --------------   dimensions revues
        double precision,dimension(npm)::b
        double precision,dimension(ng)::vectvalfin,vectretard,scoreTT,scoreYT
        double precision,dimension(nvc)::bh1
        double precision,dimension(:),allocatable::brisq
        double precision,dimension(nvarxevt)::bvexp ! le vecteur des coef des var exp pour les fonctions de risque
        double precision,dimension(sumidxevt)::vexpiact
        double precision,dimension(ng)::pi
        double precision,dimension(nvarprob)::bprob,Xprob
        double precision,dimension(nxevt)::bevt,Xevt
        double precision,dimension(-2:ntrspl-3)::splaa
        double precision,dimension(nea,nea)::Ut
        double precision,dimension(ng,nevt)::surv0,surv,risq,survint
        double precision,dimension(4)::tacti


        double precision,dimension(1)::bevtint
        double precision,dimension(nvarglobssG)::b3
        double precision,dimension(2*np+1)::b0
        double precision,dimension(ncont+ncontdep)::sum
        double precision,dimension(nvdepssG)::bdepssG
        double precision,dimension(nvdepG)::bdepG
        double precision,dimension(nvarglobG)::b2
        double precision,dimension((ncont+ncontdep)*nq)::b4

        double precision, dimension(maxmestot)::Y2,Y3,Y4
        double precision,dimension(maxmestot,(2*np+1))::Z
        double precision,dimension(maxmestot,nvarglobssG)::XglobssG
        double precision,dimension(maxmestot,nvdepssG)::XdepssG
        double precision,dimension(maxmestot,(ncontdep+ncont)*nq)::Xcont
        double precision,dimension(maxmestot,nvarglobG)::XglobG
        double precision,dimension(maxmestot,nvdepG)::XdepG
        double precision,dimension(maxmestot,nea)::P
        double precision,dimension(maxmestot,maxmestot)::VC,VCAR

        double precision,dimension(maxmestot*(maxmestot+1)/2)::Vi


        eps=1.D-20
        b1=0.d0
        k=0
        l=0

        !placement dans b1 des parametres a estimer, boucle sur le nombre de paramtre total
        DO j=1,npmtot
           if (ide(j).eq.1) then
              k=k+1
              b1(j)=b(k)
              !les parametres marques sont incrementes de th
              if (id.eq.k) then
                 b1(j)=b1(j)+thi
                 !C               if(CurrentProcessorID.eq.0) write(*,*)'id',id,'thi',thi
              end if
              if (jd.eq.k) then
                 b1(j)=b1(j)+thj
                 !C               if(CurrentProcessorID.eq.0) write(*,*)'jd',jd,'thj',thj
              end if
           else
              l=l+1
              b1(j)=bne(l)
           end if
        END DO
        !fin boucle NPM

        ! test_valfin=0.d0

        if (k.ne.npm) then
           !      if(CurrentProcessorID.eq.0) write(*,*) 'problem with number of parameters'
           !     if(CurrentProcessorID.eq.0) write(*,*)'entre',npm,'calcule',k
           funcpa=-1.d9
           goto 589
        end if

        !c matrice de covariance des effets aleatoires

        bh1=0.d0
        Do j=1,nvc
           bh1(j)=b1(nef+j)
        End do


        Ut=0.d0
        Ut1=0.d0
        If (idiag.eq.1) then
           do j=1,nea
              do k=1,nea
                 if (j.eq.k) then
                    Ut(j,k)=bh1(j)
                    Ut1(j,k)=bh1(j)
                 end if
              end do
           end do
        else if (idiag.eq.0) then
           do j=1,nea
              do k=j,nea
                 Ut(k,j)=bh1(j+k*(k-1)/2)
                 Ut1(k,j)=bh1(j+k*(k-1)/2)
              end do
           end do
        end if
        contevt=0
        contevt2=0
        do g=1,nevt
           contevt=contevt+nprisq(g)
           if (risqcom(g).eq.0) then
              contevt2=contevt2+ng*nprisq(g)
           else if (risqcom(g).eq.1) then
              contevt2=contevt2+nprisq(g)
           else
              contevt2=contevt2+nprisq(g)+ng-1
           endif
        enddo
        allocate(brisq(contevt2))


	    ! Parametre des covariables dans latent
              b3=0.d0
              jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep
              l=0
              kk=0

              Do k=1,nv
                 If (idglob(k).eq.1.and.idg(k).eq.0) then
                    l=l+1
                    kk=kk+1
                    b3(l)=b1(jj+kk)  ! correspond a beta
                 Else if (idglob(k).eq.2.and.idg(k).eq.0) then
                    l=l+1
                    kk=kk+1 !kk ??
		    m=0
		    b3(l+m)=b1(jj+kk+m)!+(k-3)*(-1+2*np))

                    Do m=1,2*np
                       b3(l+m)=b1(jj+kk+m)!+(k-3)*(-1+2*np))
                    End do
                    l=l+np*2 ! 2 pentes *ng param
                    kk=kk+2*np
                 Else if (idglob(k).eq.1.and.idg(k).eq.1) then
                    kk=kk+ng
                 Else if (idglob(k).eq.2.and.idg(k).eq.1) then
                    kk=kk+ng*(np+1)
                 End if
		 if(idg(k).eq.1)then
		    print*,'verify program (kk) for idg =1'
		    stop
		 end if
              End do

        !C fonction de risq a calculer pour chaque i et chaque g


        !C------boucle sur les sujets i-----------------------------
        entretard=0.d0
        jacobien=0.d0
        vrais=0.d0
        vrais_survie=0.d0
        risq=0.d0
        surv=0.d0
        surv0=0.d0
        survint=0.d0
        brisq=0.d0

        !####################################################
        !si evt a eu lieu
        if (evt.ne.0) then
           nrisq_curr=0
           numact=1
           !numact permet le remplissage de brisq
           !on boucle sur le nombre d'evt
           do kk=1,nevt
              !boucle sur le nombre de classes
              if (risqcom(kk).eq.0) then
                 do g=1,ng
                    !boucle sur le nombre de parametres par evt
                    do k=1,nprisq(kk)
                       !positivation carre ou expo
                       if (paramsurv.eq.0) then 
                          !risque commun/specifique/proportionnel
                          if (risqcom(kk).eq.2) then
                             brisq(numact)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)    
                             numact=numact+1
                          end if
                          if (risqcom(kk).eq.1) then
                             brisq(numact)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)
                             numact=numact+1
                          end if
                          if (risqcom(kk).eq.0) then
                             brisq(numact)=b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)&
                             &+nrisq_curr+k)  !commente le carre pour tester!!
                             numact=numact+1
                          end if
                       else 
                          if (risqcom(kk).eq.2) then
                             brisq(numact)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))    
                             numact=numact+1
                          end if
                          if (risqcom(kk).eq.1) then
                             brisq(numact)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))
                             numact=numact+1
                          end if
                          if (risqcom(kk).eq.0) then
                             brisq(numact)=exp(b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k))
                             numact=numact+1
                          end if
                       end IF !fin paramsurv==0
                    end do !fin boucle nb param par evt
                 end do !boucle sur les classes
              endif
              if (risqcom(kk).ge.1) then
                 !auparavant les variables explicatrices etaient ajoutees apres la fonction, ce n'est plus le cas ici
                 do k=1,nprisq(kk)
                    !positivation carre ou expo
                    if (paramsurv.eq.0) then 
                       !risque commun/specifique/proportionnel
                       if (risqcom(kk).eq.2) then
                          brisq(numact)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)    
                          numact=numact+1
                       end if
                       if (risqcom(kk).eq.1) then
                          brisq(numact)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)
                          numact=numact+1
                       end if
                       if (risqcom(kk).eq.0) then
                          brisq(numact)=b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)&
                                                    &+nrisq_curr+k)
                          numact=numact+1
                       end if
                    else 
                       if (risqcom(kk).eq.2) then
                          brisq(numact)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))    
                          numact=numact+1
                       end if
                       if (risqcom(kk).eq.1) then
                          brisq(numact)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))
                          numact=numact+1
                       end if
                       if (risqcom(kk).eq.0) then
                          brisq(numact)=exp(b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k))
                          numact=numact+1
                       end if
                    end IF !fin paramsurv==0
                 end do !fin boucle nb param par evt
              endif

              nrisq_curr=nrisq_curr+nrisqspec(kk)
              
           end do

           nrisq_curr=0
           do kk=1,nevt
              !ajout des parametres de proportionnalie
              if (risqcom(kk).eq.2) then
                 nrisq_curr=nrisq_curr+nrisqspec(kk)-ng+1
                 do g=1,ng-1
                    brisq(numact)=b1(nvarprob*(ng-1)+nrisq_curr+g)
                    numact=numact+1
                 enddo
                 nrisq_curr=nrisq_curr+ng-1
              else
                 nrisq_curr=nrisq_curr+nrisqspec(kk)
              endif
           enddo
           
        endif
        
        !bvexp obtient les params des vexp pour la classe courante ranges par evt 01 , 02 , 12 voir remplissage de valinit
        numact=1 ! garde l'endroit actuel du vecteur

        bvexp=0.d0

        do l=1,nvarxevt
           bvexp(numact)=b1((ng-1)*nvarprob+nrisq+numact)
           numact=numact+1
        enddo

	      Y1all=0.d0
        !i : represente le sujet
        DO i=1,ns

           testsurv=0.d0
           vexpiact=0.d0
           evtact=1
           kkk=1
!variables explicatives pour i
           do evtact=1,nevt
              do k=1,nv
                 if (idxevt(evtact*nv-nv+k).ge.1) then
                    vexpiact(kkk)=X(i,k)
                    kkk=kkk+1
                 endif
              enddo
           enddo

           !C         if(CurrentProcessorID.eq.0) write(*,*)'i',i


           numpat=i

           tacti=0.d0
           tacti(1)=t0(i)
           tacti(2)=t1(i)
           tacti(3)=t2(i)
           tacti(4)=t3(i)

              retard=0.d0
              if (evt.ne.0) then
                 do g=1,ng
                    call vraisSmooth(i,contevt+nva01+nva02+nva12,brisq,bvexp,g,vexpiact,tacti,valfin,retard2,scoreTT,&
									&scoreYT,SUTi0,STi0,0,0,ni)
                    vectvalfin(g)=valfin ! vecteur pour un individu sa contribution en fonction de la classe
                    vectretard(g)=exp(-retard2) !vecteur entree retardee
                 enddo
              end if !fin condition evt != 0
              !###############################################

	      jj=0
              Do q=1,nq
                 Do j=1,nmes(i,q)
			jj=jj+1
		  end do
	      end do
              nmestot=jj


              !         do j=1,nmestot
              !            if(CurrentProcessorID.eq.0) write(*,*)'Z',(Z(jj,k),k=1,np+1)
              !         end do
              !         stop

              XdepssG=0.d0
              l=0
              Do k=1,nvdep
                 If (idvdep(k).eq.1) then
                    l=l+1
                    jj=0
                    do q=1,nq
                       do j=1,nmes(i,q)
                          jj=jj+1
                          XdepssG(jj,l)=Xdep(Navant_i(i)+jj,k)
                       end do
                    End do
                    if (jj.ne.nmestot) then
                       if(CurrentProcessorID.eq.0) write(*,*)'probleme XdepssG(jj,k)'
                       funcpa=-1.d9
                       goto 589
                    end if
                 end if
              end do
              if (l.ne.nvdepssG) then
                 if(CurrentProcessorID.eq.0) write(*,*)'problem def nvdepssG'
                 funcpa=-1.d9
                 goto 589
              end if

              !C creation du vecteur de variables explicatives dep du temps : XdepG avec mixture


              XdepG=0.d0
              l=0
              Do k=1,nvdep
                 If (idvdep(k).eq.2) then
                    l=l+1
                    jj=0
                    do q=1,nq
                       do j=1,nmes(i,q)
                          jj=jj+1
                          XdepG(jj,l)=Xdep(Navant_i(i)+jj,k)
                       end do
                    End do
                    if (jj.ne.nmestot) then
                       if(CurrentProcessorID.eq.0) write(*,*)'probleme XdepG(jj,k)'
                       funcpa=-1.d9
                       goto 589
                    end if
                 end if
              end do
              if (l.ne.nvdepG) then
                 if(CurrentProcessorID.eq.0) write(*,*)'problem def nvdepG'
                 funcpa=-1.d9
                 goto 589
              end if

              !C------ creation de XglobG(j,k) (avec mixture) ----------

              XglobG=0.d0

              l=0
              Do k=1,nv
                 If (idglob(k).eq.1.and.idg(k).eq.1) then
                    l=l+1
                    Do j=1,nmestot
                       XglobG(j,l)=X(i,k)
                    End do

                 Else if (idglob(k).eq.2.and.idg(k).eq.1) then
                    jj=0
                    l=l+1
                    Do q=1,nq
                       Do j=1,nmes(i,q)
                          jj=jj+1
                          Do m=0,np
                             if (m.ge.1) then
                                XglobG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)* (0.1**(m-1))
                             else
                                XglobG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)
                             end if
                          End do
                       End do
                    End do
                    l=l+np
                 End if
              End do


              if (nvarglobssG+nvarglobG*ng.ne.nvarglob) then
                 if(CurrentProcessorID.eq.0) write(*,*) 'pb : def nvarglobssG/G'
                 funcpa=-1.d9
                 goto 589
              end if


              !C------ creation de Xcont(q,k,j) -----------------
              !C creation des variables pour les contrastes

              Xcont=0.d0
              nmoins=0
              do q=1,nq
                 do j=1,nmes(i,q)
                    h=0
                    do k=1,nvdep
                       if (idvdep(k).eq.3.or.idvdep(k).eq.4) then
                          h=h+1
                          Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)=Xdep(Navant_i(i)+jj,k)
                       end if
                    end do
                    h=0
                    do k=1,nv
                       if (idtest(k).eq.1) then
                          h=h+1
                          Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)=X(i,k)
                       end if
                       if (idtest(k).eq.2) then
                          Do m=0,np
                             h=h+1
                             if (m.ge.1) then
                                Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)     & 
                                     =X(i,k)*(t(Navant_i(i)+nmoins+j)**m)*(0.1**(m-1))

                             else
                                Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)   &
                                     =X(i,k)*(t(Navant_i(i)+nmoins+j)**m)
                             end if
                          End do
                       end if
                    end do
                    if (h.ne.ncont) then
                       if(CurrentProcessorID.eq.0) write(*,*)'h=',h,'problem in cont 2'
                       funcpa=-1.d9
                       goto 589
                    end if
                 End do
                 nmoins=nmoins+nmes(i,q)
              end do
              if (nmoins.ne.nmestot) then
                 if(CurrentProcessorID.eq.0) write(*,*)'problem in cont'
                 funcpa=-1.d9
                 goto 589
              end if


              !ici marqueur
              VCAR=0.d0
              if (iAR.gt.0.and.QUECONT.eq.1) then
                 jj=0
                 Do q=1,nq
                    jj1=0
                    Do q1=1,nq
                       Do k=1,nmes(i,q)
                          do k1=1,nmes(i,q1)
                             if (iAR.eq.1) then
                                VCAR(jj+k,jj1+k1)=(b1(nef+nvc+ng-1+nq+1)**2)*min(t(Navant_i(i)+jj+k),t(Navant_i(i)+jj1+k1))
                             else if (iAR.eq.2) then
                                VCAR(jj+k,jj1+k1)=(b1(nef+nvc+ng-1+nq+1)**2) &
                                     *exp(-b1(nef+nvc+ng-1+nq+2)*abs(t(Navant_i(i)+jj+k)-t(Navant_i(i)+jj1+k1)))
                             end if
                          end do
                       end do
                       jj1=jj1+nmes(i,q1)
                    end do
                    jj=jj+nmes(i,q)
                 end do
              end if


              !C         if(CurrentProcessorID.eq.0) write(*,*)'i2',i


              !C------- Creation du vecteur Y1 --------------------------
              !C------------ et des matrices Kmat et Sigmaq ------


              jac=0.d0
              Y1=0.d0
              Sigmaq=0.d0
              Vw=0.d0
              jj=0.d0
              jj1=0.d0
              ntrtot=0
              Do q=1,nq
                 if (idord(q).gt.1) then
                    ntrtot=ntrtot+idord(q)-1
                    jj1=jj1+nmes(i,q)
                 else if (idord(q).le.0) then
                    aa1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &
                         +ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot))/     &
                         (1+exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)     &
                         +nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1     &
                         +ntrtot)))
                    bb1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &
                         +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &
                         +ncont*(nq-1)+2+ntrtot))/( 1.D0+ exp(b1(nvarprob*(ng-1)+nrisq &
                         +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+2+ntrtot)))


                    bb1=aa1*(1.d0-aa1)*bb1
                    !C je change la modelisation
                    !C               bb1=aa1*bb1/4.d0

                    cc1=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1) &
                         +ncont*(nq-1)+3+ntrtot)

                    dd1=abs(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1) &
                         +ncont*(nq-1)+4+ntrtot))


                    aa=aa1*aa1*(1-aa1)/bb1-aa1
                    bb=aa*(1-aa1)/aa1

                    if ((aa.ne.aa).or.(bb.ne.bb)) then
                       if(CurrentProcessorID.eq.0) write(*,*)i,'probleme Beta',aa,bb,b1(nvarprob*(ng-1)+nrisq &
                            +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot), &
                            b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                            +ncontdep*(nq-1)+ncont*(nq-1)+2+ntrtot)
                       funcpa=-1.d9
                       if(CurrentProcessorID.eq.0) write(*,*)'problem in Beta transform'
                       goto 589
                       !C                  stop
                    end if

                    Do k=1,nmes(i,q)
                       do l=1,nmes(i,q)
                          Vw(jj+k,jj+l)=b1(nef+nvc+ng-1+q)**2
                          if(k.eq.l) then
                             ytemp=(Y(Navant_i(i)+jj1+k)-range(1+(q-1)*2))/(range(2+(q-1)*2)-range(1+(q-1)*2))
                             Y1(jj+k)=(betai(aa,bb,ytemp,tbeta)-cc1)/dd1
                             jac = jac + log(abs(beta_densite(ytemp,aa,bb))/dd1)
                             jac=jac-log(abs(range(2+(q-1)*2)-range(1+(q-1)*2)))
                             if ((Y1(jj+k)+1).eq.Y1(jj+k)) then
                                if(CurrentProcessorID.eq.0) write(*,*)'Y1(jj+k).eq.NAN'
                                funcpa=-1.d9
                                goto 589
                             end if

                             !C     OK pour la fonction de repartition
                             Sigmaq(jj+k,jj+l)=b1(nef+nvc+ng-1+nq+nAR+q)**2
                          end if
                       end do
                    end do

                    ntrtot=ntrtot+ntr
                    jj=jj+nmes(i,q)
                    jj1=jj1+nmes(i,q)

                 elseif (idord(q).eq.1) then

                    bb=b1(nvarprob*(ng-1)+nrisq +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                         +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)

                    do kk=2,ntrspl
                       splaa(kk-3)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                            +ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot)*b1(nvarprob*(ng-1)+nrisq &
                            +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot)
                    end do

                    Do k=1,nmes(i,q)

                       do l=1,nmes(i,q)
                          Vw(jj+k,jj+l)=b1(nef+nvc+ng-1+q)**2

                          if(k.eq.l) then

                             !C     creation du vecteur de donnees transformees

                             !C     creation du vecteur de donnees transformees
                             ll=0
                             if (Y(Navant_i(i)+jj1+k).eq.zitr(nztr,q)) then
                                ll=nztr-1
                             end if
                             som=0.d0
                             do kk = 2,nztr

                                if ((Y(Navant_i(i)+jj1+k).ge.zitr(kk-1,q)).and.(Y(Navant_i(i)+jj1+k).lt.zitr(kk,q)))then
                                   ll=kk-1
                                end if
                             end do
                             if (ll.lt.1.or.ll.gt.nztr-1) then
                                if(CurrentProcessorID.eq.0) write(*,*) 'probleme dans funcpa splines'
                                if(CurrentProcessorID.eq.0) write(*,*) 'll=',ll,'Y=',Y(Navant_i(i)+jj1+k)
                                funcpa=-1.d9
                                goto 589
                             end if
                             if (ll.gt.1) then
                                do ii=2,ll
                                   som=som+splaa(ii-3)
                                end do
                             end if

                             Y1(jj+k)=bb+ som +splaa(ll-2)*im2(Navant_i(i)+jj1+k)+splaa(ll-1) &
                                  *im1(Navant_i(i)+jj1+k)+splaa(ll)*im(Navant_i(i)+jj1+k)

                             jac = jac + log(splaa(ll-2)*mm2(Navant_i(i)+jj1+k)+splaa(ll-1)  &
                                  *mm1(Navant_i(i)+jj1+k)+splaa(ll)*mm(Navant_i(i)+jj1+k))

                             !                        if(CurrentProcessorID.eq.0) write(*,*)i,jj+k,Y1(jj+k),jac
                             !                        if(CurrentProcessorID.eq.0) write(*,*)i,mm2(Navant_i(i)+jj+k),mm1(Navant_i(i)+jj+k),mm(Navant_i(i)+jj+k)
                             !                       if(CurrentProcessorID.eq.0) write(*,*)i,im2(Navant_i(i)+jj+k),im1(Navant_i(i)+jj+k),im(Navant_i(i)+jj+k)

                             !C                        end if
                             !C     OK pour la fonction de repartition
                             Sigmaq(jj+k,jj+l)=b1(nef+nvc+ng-1+nq+nAR+q)**2

                          end if
                       end do
                    end do

                    ntrtot=ntrtot+ntrspl
                    jj=jj+nmes(i,q)
                    jj1=jj1+nmes(i,q)

                 end if
		do k=1,nmes(i,q)
			Y1all(ns*(q-1)+i,k)=Y1(k)
		end do
              End do  !loop on nq ended

              !C         if(CurrentProcessorID.eq.0) write(*,*)'fin jac',jac,(Y1(j),j=1,nmestot)


              IF (QUECONT.eq.1) THEN
                 jacobien=jacobien+jac
                 If (Nmestot.ne.jj) THEN
                    IF(CURRENTPROCESSORID.EQ.0) WRITE(*,*)'problem nmestot <> som(nmes)'
                    funcpa=-1.d9
                    goto 589
                 ENDIF
              ENDIF



              !C------Vecteurs de parms sans mixtures -------------------
              bdepssG=0.d0
              jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)
              l=0
              kk=0
              Do k=1,nvdep
                 If (idvdep(k).eq.1) then
                    l=l+1
                    kk=kk+1
                    bdepssG(l)=b1(jj+kk)
                 end if
                 If (idvdep(k).eq.2) then
                    kk=kk+ng
                 end if
              end do
              if (l.ne.nvdepssG) then
                 if(CurrentProcessorID.eq.0) write(*,*) 'probleme longueur vecture bdepssG'
                 funcpa=-1.d9
                 goto 589
              end if


              !C     b4 : contrastes sur les tests
              b4=0.d0
              sum=0.d0
              do j=1,ncontdep
                 do q=1,nq-1
                    b4((ncont+ncontdep)*(q-1)+j)=b1(nvarprob*(ng-1)+nrisq   &
                         +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+(j-1)*(nq-1)+q)
                    sum(j)=sum(j)+b4((ncont+ncontdep)*(q-1)+j)
                 end do
                 b4((ncont+ncontdep)*(nq-1)+j)=-sum(j)
              end do
              do j=1,ncont
                 do q=1,nq-1
                    b4((ncont+ncontdep)*(q-1)+ncontdep+j)=b1(nvarprob*(ng-1) &
                         +nrisq +nvarxevt +ng+ng*(2*np+1)+nprmdep+nvarglob & 
                         + ncontdep*(nq-1)+(j-1)*(nq-1)+q)
                    sum(ncontdep+j)=sum(ncontdep+j)+b4((ncont+ncontdep)*(q-1)+ncontdep+j)
                 end do
                 b4((ncont+ncontdep)*(nq-1)+ncontdep+j)=-sum(ncontdep+j)
              end do


              !C------ vraisemblance -------------------------------------

              !###############################################
              !variable continue ou non
              if (QUECONT.eq.1) then
                 vrais=vrais-nmestot*dlog(dble(2*3.14159265))
              end if

              expo=0.d0

              !si nombre de classe=1
              if (ng.eq.1) then !suite dans un else

              !C------creation de la matrice des puissances du temps -----
              !partie : temps latents

              Z=0.d0
              jj=0

              Do q=1,nq
                 Do j=1,nmes(i,q)
                    jj=jj+1
		    jjj=1
                    l=0

		    Z(jj,jjj)=1.d0
		    jjj=jjj+1

                    Do k=2,np+1
			    temp=t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+1)-65.0d0)/10.d0
			    Z(jj,jjj)=(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+1)-65.0d0)/10.d0)**(k-1)
			    temp=(temp**2+trn)**0.5d0
			    Z(jj,jjj+1)=temp!((b1((nvarprob)*(ng-1)+nrisq+nvarxevt+1)-65.0d0)/10.d0-t(Navant_i(i)+jj))**(k-1)*temp
		            If (k.gt.2) then
				Z(jj,jjj)=Z(jj,jjj)*(0.1)**(k-2)
				Z(jj,jjj+1)=Z(jj,jjj+1)*(0.1)**(k-2)
				print*,'pas code pour np > 1'
				stop
		                !Z(jj,k)=Z(jj,k)*(0.1)**(k-2)!! a changer !!
		            end if

			    jjj=jjj+2
                    End do
                 End do
              End do


              Z1=0.d0
	      jj=0
              Do q=1,nq
                 Do j=1,nmes(i,q)
                    jj=jj+1
                    l=0
                    Do k=1,np+1
		       If(idea(k).eq.1) then
		          l=l+1
		          if(k.eq.1) Z1(jj,l)=Z(jj,k) 
			  if(k.eq.2)then
				Z1(jj,l)=Z(jj,k)/2.0d0*(1.0d0-(Z(jj,k)**2+trn)**0.5d0/Z(jj,k))
			  	Z1(jj,l+1)=Z(jj,k)/2.0d0*(1.0d0+(Z(jj,k)**2+trn)**0.5d0/Z(jj,k))
			  end if
		       end if
                    End do
                 End do
              End do

                 if (evt.ne.0) then 
                    if (nxevt.ne.0) then
                       bevt=0.d0
                       l=1
                       do kk=1,NEVT
                          do k=1,nv
                             if (idxevt((kk-1)*nv+k).eq.1.or.idxevt((kk-1)*nv+k).eq.2) then                           
                                bevt(l)=b1(nvarprob*(ng-1)+nrisq+l)
                                l=l+1
                             end if
                          end do
                       end do

                       IF (max(l-1,0).ne.nxevt) then
                          if(CurrentProcessorID.eq.0) write(*,*)'problem bevt l',l,'nxevt',nxevt
                          funcpa=-1.d9
                          goto 589
                       end if
                    end If !fin if nxevt!=0

                    bevtint=0.d0
                    if (nvdepsurv.ne.0) then
                       bevtint(1)=b1(nvarprob*(ng-1)+nrisq+nvarxevt)
                    end if

                 end if !fin if evt!=0
                 !###############################################

                 IF (QUECONT.eq.1) then

                    VC=0.d0
                    P=0.d0
                    P=MATMUL(Z1,Ut)
                    VC=0.d0
                    VC=MATMUL(P,transpose(P))+Sigmaq+Vw+VCAR

                    !c     Vi en vecteur
                    jj=0
                    do j=1,nmestot
                       do k=j,nmestot
                          jj=j+k*(k-1)/2
                          Vi(jj)=VC(j,k)
                       end do
                    end do

                    !C               if(CurrentProcessorID.eq.0) write(*,*)'avant dsinv'

                    ier=0
                    CALL DSINV(Vi,nmestot,eps,ier,det,0)
                    if (ier.eq.-1) then
                       !  if(CurrentProcessorID.eq.0) write(*,*)'dsinv vi ier',ier,' i=',i,' nmes(i)=',nmestot
                       !  if(CurrentProcessorID.eq.0) write(*,*)'b1', (b1(k),k=1,npm)
                       nnn=nmestot*(nmestot+1)/2
                       !  if(CurrentProcessorID.eq.0) write(*,*)'Vi', (Vi(k),k=1,nnn)
					   if(CurrentProcessorID.eq.0) write(*,*) 'Problem DSINV'
                       funcpa=-1.d9
                       goto 589
                    end if


                    !C               if(CurrentProcessorID.eq.0) write(*,*)'avpres dsinv'


                    !c     retransformation du vecteur Vi en matrice :
                    VC=0.d0
                    do j=1,nmestot
                       do k=1,nmestot
                          if (k.ge.j) then
                             VC(j,k)=Vi(j+k*(k-1)/2)
                          else
                             VC(j,k)=Vi(k+j*(j-1)/2)
                          end if
                       end do
                    end do


                    vrais=vrais-det

                 END IF 

                 do k=1,2*np+1
                    b0(k)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+k)
                 end do
                 mu=0.d0
                 Y3=0.d0

              XglobssG=0.d0

	      g=1
              l=0
              Do k=1,nv
                 If (idglob(k).eq.1.and.idg(k).eq.0) then
                    l=l+1
                    Do j=1,nmestot
                       XglobssG(j,l)=X(i,k)
                    End do

                 Else if (idglob(k).eq.2.and.idg(k).eq.0) then
                    jj=0
                    l=l+1
                    Do q=1,nq
                       Do j=1,nmes(i,q)
                          jj=jj+1

                          Do m=0,np
   			     temp=(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0
                             !C Attention : changer les Xglob avec interaction
                             if (m.ge.1) then
                                XglobssG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)-temp)!((t(Navant_i(i)+jj)-temp)**m)*(0.1**(m-1))
			        temp=dble(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)
			        temp=(temp**2+trn)**0.5d0!/temp
			        XglobssG(jj,l+m+1)=X(i,k)*temp!((t(Navant_i(i)+jj)-(b1((nvarprob)*&
						!(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)**m)*(0.1**(m-1))*temp
                             else
                                XglobssG(jj,l+m)=X(i,k)

                             end if
                          End do

                       End do
                    End do
                    l=l+2*np
                 End if
              End do


                 mu=MATMUL(Z,b0)+MATMUL(XglobssG,b3)+MATMUL(XdepssG,bdepssG)   &
                      +MATMUL(Xcont,b4)

                 !C            if(CurrentProcessorID.eq.0) write(*,*)'mu',(mu(j),j=1,nmestot)

                 IF (QUECONT.eq.1) THEN
                    Y3=Y1-mu
                    Y4=MATMUL(VC,Y3)
                    Y5=DOT_PRODUCT(Y3,Y4)
                    vrais=vrais-Y5
                 end if

                 IF (QUECONT.eq.0) then
                    temp=vraisobs()
                    !C               if(CurrentProcessorID.eq.0) write(*,*)'i',i,'temp',temp
                    vrais=vrais+2*log(temp)

                    !C               if(thi.eq.0.and.thj.eq.0) then
                    !C                  if(CurrentProcessorID.eq.0) write(*,*)'i',i,temp,log(temp),vrais
                    !C               end if

                 END IF

                 !###############################################
                 if (evt.ne.0) then
                    surv_glob=0.d0
                    survint_glob=0.d0
                    surv0_glob=0.d0
                    varexpsurv=0.d0
                    nxevtcurr=0
                    fevt=0.d0


                    vrais=vrais+2*(valfin)
                    entretard=entretard-retard2
                  
                    vrais_survie=vrais_survie+2*(valfin)
                  
                 end if !fin boucle evt!=0

              else !ng != 1


                 !calcul du vecteur de proba d'appartenir a une classe, inutile precedemment
                 Xprob=0.d0
                 Xprob(1)=1
                 l=0
                 !les variable exp pour la prob
                 do k=1,nv
                    if (idprob(k).eq.1) then
                       l=l+1
                       Xprob(1+l)=X(i,k)
                    end if
                 end do

                 if (l+1.ne.nvarprob) then
                    if(CurrentProcessorID.eq.0) write(*,*)'problem nvarprob funcpa'
                    funcpa=-1.d9
                    goto 589
                 end if

                 pi=0.d0
                 temp=0.d0
                 Do g=1,ng-1
                    do k=1,nvarprob
                       bprob(k)=b1((k-1)*(ng-1)+g)
                    end do
                    temp=temp+exp(DOT_PRODUCT(bprob,Xprob))
                    pi(g)=exp(DOT_PRODUCT(bprob,Xprob))
                    if ((pi(g)+1).eq.pi(g)) then
                       if(CurrentProcessorID.eq.0) write(*,*)'g',g,'pi(g)',pi(g)
                       if(CurrentProcessorID.eq.0) write(*,*)'bprob',(bprob(j),j=1,nvarprob)
                       if(CurrentProcessorID.eq.0) write(*,*)'Xprob',(Xprob(j),j=1,nvarprob)
                       funcpa=-1.d9
                       goto 589
                    end if
                 end do

                 pi(ng)=1/(1+temp)
                 retard=0.d0

                 do g=1,ng ! somme sur g

                    if (g.ne.ng) then
                       pi(g)=pi(g)*pi(ng)
                    end if
                
                    if ((pi(g)+1).eq.pi(g)) then
                       if(CurrentProcessorID.eq.0) write(*,*)'g',g,'pi(g)',pi(g)
                       if(CurrentProcessorID.eq.0) write(*,*)'bprob',(bprob(j),j=1,nvarprob)
                       if(CurrentProcessorID.eq.0) write(*,*)'Xprob',(Xprob(j),j=1,nvarprob)
                       funcpa=-1.d9
                       goto 589
                    end if




              !C------creation de la matrice des puissances du temps -----
              !partie : temps latents

              Z=0.d0
              jj=0

              Do q=1,nq
                 Do j=1,nmes(i,q)
                    jj=jj+1
		    jjj=1
                    l=0

		    Z(jj,jjj)=1.d0
		    jjj=jjj+1

                    Do k=2,np+1
			    Z(jj,jjj)=(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)**(k-1)
			    temp=dble(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)
			    !temp=(temp**2+trn)**0.5d0/temp
			    !Z(jj,jjj+1)=(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)**(k-1)*temp
			    temp=(temp**2+trn)**0.5d0
			    Z(jj,jjj+1)=temp!

		            If (k.gt.2) then
			!	Z(jj,jjj)=Z(jj,jjj)*(0.1)**(k-2)
		!		Z(jj,jjj+1)=Z(jj,jjj+1)*(0.1)**(k-2)*temp
		                !Z(jj,k)=Z(jj,k)*(0.1)**(k-2)!! a changer !!
				print*,'pas code pour np=2'
				stop
		            end if

			    jjj=jjj+2
                    End do
                 End do
              End do

              Z1=0.d0
	      jj=0
              Do q=1,nq
                 Do j=1,nmes(i,q)
                    jj=jj+1
                    l=0
                    Do k=1,np+1
		       If(idea(k).eq.1) then
		          l=l+1
		          if(k.eq.1) Z1(jj,l)=Z(jj,k) 
			  if(k.eq.2)then
				Z1(jj,l)=Z(jj,k)/2.0d0*(1.0d0-(Z(jj,k)**2+trn)**0.5d0/Z(jj,k))
			  	Z1(jj,l+1)=Z(jj,k)/2.0d0*(1.0d0+(Z(jj,k)**2+trn)**0.5d0/Z(jj,k))

			  end if
		       end if
                    End do
                 End do
              End do

              !C------ creation de XglobssG(j,k) (sans mixture) ----------

              XglobssG=0.d0

              l=0
              Do k=1,nv
                 If (idglob(k).eq.1.and.idg(k).eq.0) then
                    l=l+1
                    Do j=1,nmestot
                       XglobssG(j,l)=X(i,k)
                    End do
                 Else if (idglob(k).eq.2.and.idg(k).eq.0) then
                    jj=0
                    l=l+1
                    Do q=1,nq
                       Do j=1,nmes(i,q)
                          jj=jj+1

                          Do m=0,np
   			     temp=(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0
                             !C Attention : changer les Xglob avec interaction
                             if (m.ge.1) then
                                XglobssG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)-temp)!((t(Navant_i(i)+jj)-temp)**m)*(0.1**(m-1))
			        temp=dble(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)
			        temp=(temp**2+trn)**0.5d0!/temp
			        XglobssG(jj,l+m+1)=X(i,k)*temp!((t(Navant_i(i)+jj)-(b1((nvarprob)*&
						!(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)**m)*(0.1**(m-1))*temp
                             else
                                XglobssG(jj,l+m)=X(i,k)
                             end if

                          End do

                       End do
                    End do
                    l=l+2*np
                 End if
              End do

                    !C ---- matrice de variance covariance

                    VC=0.d0
                    P=0.d0
                    Ut1=0.d0
                    if (g.eq.ng) then
                       Ut1=Ut
                    else
                       Ut1=Ut*b1(nef+nvc+g)
                    end if

                    IF (QUECONT.eq.1) THEN
                       P=MATMUL(Z1,Ut1)
                       VC=MATMUL(P,transpose(P))+Sigmaq+Vw+VCAR

                       !c     Vi en vecteur
                       jj=0
                       do j=1,nmestot
                          do k=j,nmestot
                             jj=j+k*(k-1)/2
                             Vi(jj)=VC(j,k)
                          end do
                       end do

                       CALL DSINV(Vi,nmestot,eps,ier,det,0)

                       if (ier.eq.-1) then
                          ! if(CurrentProcessorID.eq.0) write(*,*)'dsinv vi ier',ier,' i=',i,' nmes(i)=',nmestot,'num',numero(i)
                          ! if(CurrentProcessorID.eq.0) write(*,*)'b1', (b1(k),k=1,npm)
                          nnn=nmestot*(nmestot+1)/2
                          ! if(CurrentProcessorID.eq.0) write(*,*)'Vi', (Vi(k),k=1,nnn)
						  print*,'probleme DSINV2'
                          funcpa=-1.d9
                          goto 589
                       end if

                       !c     retransformation du vecteur Vi en matrice :
                       VC=0.d0
                       do j=1,nmestot
                          do k=1,nmestot
                             if (k.ge.j) then
                                VC(j,k)=Vi(j+k*(k-1)/2)
                             else
                                VC(j,k)=Vi(k+j*(j-1)/2)
                             end if
                          end do
                       end do
                    ENDIF



                    !C----- parametres avec mixtures ---------------------------

                    b0=0.d0
                    Do j=1,(2*np+1)
                       b0(j)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+(j-1)*ng+g)
                    End do

                    !C creation parametres mixture vdep du temps
                    bdepG=0.d0
                    jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)

                    l=0
                    kk=0
                    Do k=1,nvdep
                       If (idvdep(k).eq.2) then
                          l=l+1
                          kk=kk+g
                          bdepG(l)=b1(jj+kk)
                       end if
                       If (idvdep(k).eq.1) then
                          kk=kk+1
                          kk=kk+1
                          kk=kk+1
                       end if
                    end do
                    if (l.ne.nvdepG) then
                       if(CurrentProcessorID.eq.0) write(*,*) 'probleme longueur vecteur bdepG'
                       funcpa=-1.d9
                       goto 589
                    end if

                    b2=0.d0
                    l=0
                    kk=0

                    Do k=1,nv
                       If (idglob(k).eq.1.and.idg(k).eq.1) then
                          l=l+1
                          kk=kk+g

                          b2(l)=b1(nvarprob*(ng-1)+nrisq+nvarxevt &
                               +ng+ng*(2*np+1)+nprmdep+kk) ! correspond a beta(g)

                       Else if (idglob(k).eq.2.and.idg(k).eq.1) then
                          l=l+1
                          kk=kk+g
                          Do m=1,np+1
                             b2(l+m-1)=b1(nvarprob*(ng-1)+nrisq &
                                  +nvarxevt +ng+ng*(2*np+1)+nprmdep+(m-1)*ng+kk)
                          End do
                          l=l+np
                          kk=kk+np*ng

                       Else if (idglob(k).eq.1.and.idg(k).eq.0) then
                          kk=kk+1

                       Else if (idglob(k).eq.2.and.idg(k).eq.0) then
                          kk=kk+(2*np+1)

                       End if

                    End do

                    !C     if(CurrentProcessorID.eq.0) write(*,*)'g=',g,'b2=',(b2(j),j=1,l)



                    !C------calcul de l'expo ---------------------------------
                    mu=0.d0
                    Y2=0.d0
                    Y3=0.d0

!                    Y2=MATMUL(Z,b0)+MATMUL(XglobG,b2)+MATMUL(XdepssG,bdepssG)
!                    Y2=Y2+MATMUL(XglobssG,b3)+MATMUL(XdepG,bdepG)

		! Changement de pente
                    Y2=MATMUL(Z,b0)+MATMUL(XglobG,b2)+MATMUL(XdepssG,bdepssG)
                    Y2=Y2+MATMUL(XglobssG,b3)+MATMUL(XdepG,bdepG)

                    mu=0.d0
                    mu=Y2+MATMUL(Xcont,b4)
                    !###############################################
                    IF (QUECONT.eq.0) THEN !existence de variables ordinales
                       gcourant=g
                       temp=vraisobs()
                       if (evt.eq.0) then !si pas d'evt 
                          expo=expo+pi(g)*temp

                       else !au contraire si on a evt

           
                          expo=expo+pi(g)*temp*exp(&
                               vectvalfin(g))
                          testsurv=testsurv+pi(g)*exp(vectvalfin(g))
                       end if !fin if devt(i)==0

                    end if

                    !variable continues uniquement

                    IF (QUECONT.eq.1) THEN
                       Y3=0.d0
                       Y3=Y1-mu
                       Y4=0.d0
                       Y4=MATMUL(VC,Y3)
                       Y5=0.d0
                       Y5=DOT_PRODUCT(Y3,Y4)

                       if (evt.ne.0) then 

                          expo=expo+pi(g)*exp(&
                               +(-det-Y5)/2.d0+vectvalfin(g))
                          testsurv=testsurv+pi(g)*exp(vectvalfin(g))

                       Else !si pas evt calcul de expo
                          expo=expo+pi(g)*exp((-det-Y5)/2.d0)
                       end if

                    end if !fin IF que continu

                    ! calcule retard en fonction de classe ou sortie de vraisSmooth
                    retard=retard+pi(g)*vectretard(g)

                 end do !fin de la boucle sur le nombre de classes
                 !###############################################
                 entretard=entretard+log(retard)
				 
                 if (expo.le.0) then !verification de la valeur de expo
                    !               if(CurrentProcessorID.eq.0) write(*,*)'probleme funcpa expo i',i,'expo',expo,'devt' &
                    !                   ,devt(i),(surv(i,g),g=1,ng),exp((-det-Y5)/2.d0),det,Y5
                    !		if(CurrentProcessorID.eq.0) write(*,*)'survint',(survint(i,g),g=1,ng)
                     if(CurrentProcessorID.eq.0) write(*,*)'problem in expo computation'
					 funcpa=-1.d9
				 print*,'probleme expo',i,c(i),expo

                    goto 589 
                 end if

                 vrais=vrais+2*log(expo)

                 vrais_survie=vrais_survie+2*log(testsurv)
              end if !fin du IF sur le nb de classes
           end do !fin de la boucle sur les individus

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           if (idtrunc.eq.0) then !troncature gauche inexistante
              entretard=0.d0
           end if

           IF (QUECONT.eq.1) THEN !calcul vraisemblance finale
              if (ng.ne.0) then
                 vrais_survie=vrais_survie/2.d0
              end if
              funcpa=vrais/2.d0+jacobien-entretard
           ELSE
              funcpa=vrais/2.d0-entretard
              if (ng.ne.0) then
                 vrais_survie=vrais_survie/2.d0
              end if
           end if
        
           if (funcpa.ne.funcpa) then 
              funcpa=-1.d9
              if(CurrentProcessorID.eq.0) write(*,*)'problem in funcpa computation'
           end If

           if (thi.eq.0.and.thj.eq.0) then
              ! if(CurrentProcessorID.eq.0) write(*,*)'funcpa',funcpa,entretard
              ! if(CurrentProcessorID.eq.0) write(*,*)'vrais_survie',vrais_survie
      
           end if

589        continue

           !      stop
           return

         end function FUNCPA

! }}}


!C ================================================================
!C
!C    VRAISOBS
!C
!C =================================================================

! {{{

      double precision function vraisobs()


      use lois_normales
      use donnees_indiv
      use commun


      implicit none
      double precision:: eps
      external ::vraistot
      integer::j
      INTEGER ::NDIM2, NF2, MINPTS, MAXPTS, RESTAR
      INTEGER ::NEVAL, IFAIL
      DOUBLE PRECISION ::EPSABS, EPSREL,funvls
      DOUBLE PRECISION,dimension(2) ::RESULT, ABSERR2
      DOUBLE PRECISION,dimension(1000) ::WORK
        double precision,dimension(nea)::Xea
!c     fin de declaration des variables de mvndst
!      integer::numpat1
        double precision, dimension(2,51)::Gauss


      ndim2=nea
      nf2=nf
      eps=1.d-20

      MINPTS=30
      MAXPTS=500                !  100 !500

      EPSABS=1.d-100
      EPSREL=1.d-100
      RESTAR=0
      result(1)=0.d0
      if (nea.gt.1) then
         !        if (nea.gt.1) then
         
         MINPTS=30
         MAXPTS=500               !  100 !500
         EPSABS=1.d-100
         EPSREL=1.d-100
         RESTAR=0
         call  hrmsym( NDIM2, NF2, MINPTS, MAXPTS, vraistot, EPSABS, &
              EPSREL, RESTAR, RESULT, ABSERR2, NEVAL, IFAIL, WORK)
         
         if (result(1).le.1.d-300) then
            result(1)=1.d-300
         end if
         vraisobs=result(1)
         
         return
         
      elseif (nea.eq.1) then 
         
         ! on definit les points
         call gaussher(gauss,npg)
         
         
         ! boucle pour faire l'integration multiple
         do j=1,npg
           
            Xea(1)=gauss(1,j)
            call vraistot(nea,Xea,nf2,funvls)

            result(1)=result(1)+funvls*gauss(2,j)
         end do
         
         if (result(1).le.1.d-300) then
            result(1)=1.d-300
         end if

!         if (numpat.lt.5) then
!            if(CurrentProcessorID.eq.0) write(*,*)numpat,'npg',npg,result(1)
!         end if
         vraisobs=result(1)
         return
         
      else

         Xea=0.d0
         call vraistot(nea,Xea,nf,funvls)
         
         if (funvls.le.1.d-300) then
            funvls=1.d-300
         end if
         
         vraisobs=funvls
         return
      end if
      end

! }}}

!C ===============================================================
!C
!C     Vraisemblance sachant les effets aleatoires
!C
!C ===============================================================


! {{{

      subroutine vraistot(NDIM2,Xea,NF2,FUNVLS)

      use donnees_indiv
      use commun
      use optim_parallele
      use lois_normales

      implicit none
      double precision ::vraisind,scal,eps,det,gamma1,gamma0,alnorm
      double precision :: funvls,sup,inf,result,vrais,aa1,bb1,aa,bb,cc1 &
           ,beta_densite,betai,ytemp,vraisind_old,ytemp2,dd1
      integer::ndim2,nf2,i,g,k,nmoins,q,j,nmoinsc,jj,ier,ind,ntrtot,qq
      logical ::upper
      double precision,dimension(nea)::Xea,Xea2
      double precision,dimension(2,51)::gauss
      double precision,dimension(nea)::ui
      double precision,dimension(maxmestot,maxmestot):: VC
      double precision, dimension(maxmestot)::mu1
      double precision,dimension(maxmestot*(maxmestot+1)/2)::Vi
      double precision,dimension(maxmestot)::Y2,P



      i=numpat
      g=gcourant
      ndim2=nea
      nf2=nf

      Xea2=0.d0
      do j=1,nea
         Xea2(j)=Xea(j)
      end do


      eps=1.D-20
!C     effets aleatoire REcorreles

      ui=MATMUL(Ut1,Xea2)
      mu1=mu+MATMUL(Z1,ui)
      vraisind=1
      nmoins=0
      nmoinsc=0
      ntrtot=0
      upper=.false.


      DO q=1,nq
         if (idord(q).gt.1) then



            if (ide(nef+nvc+ng-1+q).eq.0.and.nmes(i,q).ne.0) then


               do j=1,nmes(i,q)
                  ind=0
                  if (Y(Navant_i(i)+nmoins+j).eq.range(1+(q-1)*2)) then
                     gamma0=(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)-mu1(nmoins+j)  )/abs(b1(nef+nvc+ng-1+nq+nAR+q))
                     vraisind=vraisind*(alnorm(gamma0,upper))
                     ind=1
                  else
                     sup=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob  &
                       +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                     inf=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                         +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                     do k=1,idord(q)-2
                        sup=sup+(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                            +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1)*b1(nvarprob*(ng-1)+nrisq &
                            +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1))
                        if (Y(Navant_i(i)+nmoins+j).eq.dble(k+range(1+(q-1)*2))) then
                           gamma1=(sup-mu1(nmoins+j))/abs(b1(nef+nvc+ng-1+nq+nAR+q))
                           gamma0=(inf-mu1(nmoins+j) )/abs(b1(nef+nvc+ng-1+nq+nAR+q))
                           vraisind=vraisind*(alnorm(gamma1,upper)-alnorm(gamma0,upper))
                           ind=1
                        end if
                        inf=sup
                     end do
                     if (Y(Navant_i(i)+nmoins+j).eq.range(2+(q-1)*2)) then
                        gamma0=(sup-mu1(nmoins+j))/abs(b1(nef+nvc+ng-1+nq+nAR+q))
                        vraisind=vraisind*(1-alnorm(gamma0,upper))
                        ind=1
                     end if
                  end if
               end do

            else if (ide(nef+nvc+ng-1+q).ne.0.and.nmes(i,q).ne.0) then

!C si il y a des RE spec aux tests

               call gaussher(gauss,npg)

               result=0.d0

               do qq=1,npg

                  vrais=1.d0

                  do j=1,nmes(i,q)
                     ind=0
                     if (Y(Navant_i(i)+nmoins+j).eq.range(1+(q-1)*2)) then
                        gamma0=(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                            +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)-mu1(nmoins+j)- &
                        abs(b1(nef+nvc+ng-1+q))*gauss(1,qq))/abs(b1(nef+nvc+ng-1+nq+nAR+q)) 
                        vrais=vrais*(alnorm(gamma0,upper))
                        ind=1
                     else
                        sup=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                           +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                        inf=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob  &    
                      +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                        do k=1,idord(q)-2
                           sup=sup+(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                            +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1)*b1(nvarprob*(ng-1)+nrisq &
                        +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1))
                           if (Y(Navant_i(i)+nmoins+j).eq.dble(k+range(1+(q-1)*2)))then
                              gamma1=(sup-mu1(nmoins+j)-abs(b1(nef+nvc+ng-1+q))*gauss(1,qq)) &
                              /abs(b1(nef+nvc+ng-1+nq+nAR+q))
                              gamma0=(inf-mu1(nmoins+j)-abs(b1(nef+nvc+ng-1+q))*gauss(1,qq)) &
                                  /abs(b1(nef+nvc+ng-1+nq+nAR+q))
                              vrais=vrais*(alnorm(gamma1,upper)-alnorm(gamma0,upper))
                              ind=1
                           end if
                           inf=sup
                        end do
                        if (Y(Navant_i(i)+nmoins+j).eq.range(2+(q-1)*2)) then
                           gamma0=(sup-mu1(nmoins+j)-abs(b1(nef+nvc+ng-1+q))*gauss(1,qq))    &    
                     /abs(b1(nef+nvc+ng-1+nq+nAR+q))
                           vrais=vrais*(1-alnorm(gamma0,upper))
                           ind=1
                        end if
                     end if
                  end do




                  result=result+gauss(2,qq)*vrais

               end do

               vraisind=vraisind*result

            end if

            ntrtot=ntrtot+idord(q)-1
            nmoins=nmoins+nmes(i,q)


         else if (idord(q).eq.0) then
            ntrtot=ntrtot+ntr
            do j=1,nmes(i,q)
               Y2(nmoinsc+j)=Y1(nmoinsc+j)-mu1(nmoins+j)
            end do
            nmoinsc=nmoinsc+nmes(i,q)
            nmoins=nmoins+nmes(i,q)


         else if (idord(q).eq.1) then
            ntrtot=ntrtot+ntrspl
            do j=1,nmes(i,q)
               Y2(nmoinsc+j)=Y1(nmoinsc+j)-mu1(nmoins+j)
            end do
            nmoinsc=nmoinsc+nmes(i,q)
            nmoins=nmoins+nmes(i,q)
         else if (idord(q).eq.-1) then



               aa1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob    &        
      +ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot))/(1+exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)     &    
       +nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot)))
           bb1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont     &      
        *(nq-1)+2+ntrtot))/( 1.D0+ exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &   
        +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+2+ntrtot)))


           bb1=aa1*(1.d0-aa1)*bb1

           cc1=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &    
      +ncont*(nq-1)+3+ntrtot)

           dd1=abs(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &   
       +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+4+ntrtot))

           aa=aa1*aa1*(1-aa1)/bb1-aa1
           bb=aa*(1-aa1)/aa1


            if (ide(nef+nvc+ng-1+q).eq.0.and.nmes(i,q).ne.0) then


           do j=1,nmes(i,q)

              gamma0=(-(cc1/dd1)-mu1(nmoins+j)) /abs(b1(nef+nvc+ng-1+nq+nAR+q))
              gamma1=((1.d0-cc1)/dd1-mu1(nmoins+j)) /abs(b1(nef+nvc+ng-1+nq+nAR+q))

              if (Y(Navant_i(i)+nmoins+j).eq.range(1+(q-1)*2)) then

                 vraisind_old=vraisind
                 vraisind=vraisind*alnorm(gamma0,upper)

              else if (Y(Navant_i(i)+nmoins+j).eq.range(2+(q-1)*2)) then

                 vraisind_old=vraisind
                 vraisind=vraisind*(1.d0-alnorm(gamma1,upper))

              else if (Y(Navant_i(i)+nmoins+j).lt.range(2+(q-1)*2).and.Y(Navant_i(i)+nmoins+j).gt.range(1+(q-1)*2)) then

                 vraisind_old=vraisind
                 vraisind=vraisind/(sqrt(dble(2*3.14159265)))
                 vraisind=vraisind/b1(nef+nvc+ng-1+nq+nAR+q)


                  ytemp=(Y(Navant_i(i)+nmoins+j)-range(1+(q-1)*2))/(range(2+(q-1)*2)-range(1+(q-1)*2))
                  ytemp2=((betai(aa,bb,ytemp,tbeta)-cc1)/dd1)-mu1(nmoins+j)

                  vraisind=vraisind*exp(-(ytemp2*ytemp2/(2.d0*(b1(nef+nvc+ng-1+nq+nAR+q)*      &   
              b1(nef+nvc+ng-1+nq+nAR+q)))))

                  vraisind=vraisind*abs(beta_densite(ytemp,aa,bb))/abs(dd1*(range(2+(q-1)*2)     &  
                     -range(1+(q-1)*2)))

              else
                 if(CurrentProcessorID.eq.0) write(*,*)'should not come here ... problem'
                 if(CurrentProcessorID.eq.0) write(*,*)'Y',Y(Navant_i(i)+nmoins+j),range(1+(q-1)*2),range(2+(q-1)*2)
                 stop
               end if
            end do



            nmoinsc=nmoinsc+nmes(i,q)
            nmoins=nmoins+nmes(i,q)
            ntrtot=ntrtot+ntr

         else if (ide(nef+nvc+ng-1+q).ne.0.and.nmes(i,q).ne.0) then

!C si il y a des RE spec aux tests

               call gaussher(gauss,npg)

               result=0.d0

               do qq=1,npg

                  vrais=1.d0

                  do j=1,nmes(i,q)

                     gamma0=(-(cc1/dd1)-mu1(nmoins+j)-abs(b1(nef+nvc+ng-1+q))*gauss(1,qq))     &        
            /abs(b1(nef+nvc+ng-1+nq+nAR+q))
                     gamma1=((1.d0-cc1)/dd1-mu1(nmoins+j)-abs(b1(nef+nvc+ng-1+q))*gauss(1,qq))     &  
                  /abs(b1(nef+nvc+ng-1+nq+nAR+q))

                     if (Y(Navant_i(i)+nmoins+j).eq.range(1+(q-1)*2)) then

                        vrais=vrais*alnorm(gamma0,upper)
                     else if (Y(Navant_i(i)+nmoins+j).eq.range(2+(q-1)*2)) then

                        vrais=vrais*(1.d0-alnorm(gamma1,upper))
                     else if (Y(Navant_i(i)+nmoins+j).lt.range(2+(q-1)*2).and.Y(Navant_i(i)+nmoins+j).gt.range(1+(q-1)*2)) then


                        vrais=vrais/(sqrt(dble(2*3.14159265)))
                        vrais=vrais/b1(nef+nvc+ng-1+nq+nAR+q)

                        ytemp=(Y(Navant_i(i)+nmoins+j)-range(1+(q-1)*2))/(range(2+(q-1)*2)-range(1+(q-1)*2))
                        ytemp2=((betai(aa,bb,ytemp,tbeta)-cc1)/dd1)-mu1(nmoins+j)-     &         
              abs(b1(nef+nvc+ng-1+q))*gauss(1,qq)

                        vrais=vrais*exp(-(ytemp2*ytemp2/(2.d0*(b1(nef+nvc+ng-1+nq+nAR+q)*     &    
                   b1(nef+nvc+ng-1+nq+nAR+q)))))

                        vrais=vrais*abs(beta_densite(ytemp,aa,bb))/abs(dd1*(range(2+(q-1)*2)     &     
                  -range(1+(q-1)*2)))

                     else
                        if(CurrentProcessorID.eq.0) write(*,*)'ne devrait pas venir ici'
                        if(CurrentProcessorID.eq.0) write(*,*)'Y',Y(Navant_i(i)+nmoins+j),range(1+(q-1)*2),range(2+(q-1)*2)
                        stop
                     end if
                  end do

                  nmoinsc=nmoinsc+nmes(i,q)
                  nmoins=nmoins+nmes(i,q)
                  ntrtot=ntrtot+ntr




                  result=result+gauss(2,qq)*vrais

               end do

               vraisind=vraisind*result
            end if

         end if
      END DO

      IF (QUEORD.eq.0.and.nmoinsc.gt.0) THEN
         VC=Sigmaq+Vw

!c     Vi en vecteur
         jj=0
         do j=1,nmoinsc
            do k=j,nmoinsc
               jj=j+k*(k-1)/2
               Vi(jj)=VC(j,k)
            end do
         end do
         ier=0
         CALL DSINV(Vi,nmoinsc,eps,ier,det,0)
         if (ier.eq.-1) then
          !  if(CurrentProcessorID.eq.0) write(*,*)'dsinv vi ier',ier,' i=',i,' nmes(i)=', nmoinsc
            stop
         end if

!c     retransformation du vecteur Vi en matrice :
         VC=0.d0
         do j=1,nmoinsc
            do k=1,nmoinsc
               if (k.ge.j) then
                  VC(j,k)=Vi(j+k*(k-1)/2)
               else
                  VC(j,k)=Vi(k+j*(j-1)/2)
               end if
            end do
         end do

         P=MATMUL(VC,Y2)
         scal=DOT_PRODUCT(Y2,P)
         vraisind=vraisind*exp(jac)
         vraisind=vraisind/((sqrt(dble(2*3.14159265)))**nmoinsc)
         vraisind=vraisind/sqrt(exp(det))
         vraisind=vraisind*exp(-scal/2.d0)

      END IF

      FUNVLS=vraisind





      end subroutine vraistot

! }}}


! =========================================================
!
!   Fonction de risque en competitif
!
!
! =========================================================


! {{{

!C-----------------------------------------------------------
!C
!C     Calcul de la fonction de risque pour la
!C     vraisemblance : subroutine qui necessite le vecteur
!C     de parametres de risque et le sujet pour lequel, sont
!C     calculees les quantites : risq, surv et surv (entree
!C     retardee)
!C
!C-----------------------------------------------------------




! }}}


    
! =========================================================
!
!   SUbroutine calcul des noeuds de splines ou stepfunctions
!
!
! =========================================================

  
! {{{  

!C============== noeuds pour step functions ou splines =================


      subroutine zi_noeuds(nzz,idzii,s,tevtt,zii)

      use commun
      implicit none
      integer ::j,jj,i,ind,seuil,ndata
	  integer, intent(in) :: s
	  integer, intent(in) :: nzz,idzii
	  double precision , dimension (nzz), intent(out) :: zii
	  real , dimension (s), intent(in) :: tevtt


      double precision,dimension(2*s)::ytemp

      double precision::pas,temp
      if (nzz <2) then
         if(CurrentProcessorID.eq.0) write(*,*) 'pas assez de noeuds : nz pour l evenement considere=',nzz
         stop
      end if

   if (idzii.eq.0) then
		  do i=1,s
			 Ytemp(i)=tevtt(i)
		 end do
	!C     appel de la procedure de tri d'un vecteur
		  ind=1
		  do while (ind.eq.1)
			 ind=0
			 do i=1,s-1
				if (ytemp(i).gt.ytemp(i+1)) then
				   temp=ytemp(i)
				   ytemp(i)=ytemp(i+1)
				   ytemp(i+1)=temp
				   ind=1
				end if
			 end do
		  end do
	  
         if(CurrentProcessorID.eq.0) write(*,*)'noeuds aux quantiles'

!C     creation de la sequence de noeuds a partir des quantiles
         zii(1)=0.d0
         if(idtrunc.eq.1) then 
            zii(1)=minval(t0)
         endif
         pas=floor(dble(s)/dble(nzz-1))
         seuil=pas
         j=1
         jj=0
         do jj=1,s
            if (jj.eq.seuil.and.j.ne.(nzz-1)) then
               j=j+1
               zii(j)=ytemp(jj)
               seuil=seuil+pas
            end if
         end do

         if (j.ne.(nzz-1)) then
            if(CurrentProcessorID.eq.0) write(*,*)'probleme definition des quantiles'
            if(CurrentProcessorID.eq.0) write(*,*)'j+1=',j+1,'nz=',nzz
            stop
         end if
         zii(nzz)=maxval(t3)
     end if
	  
      if (idzii.eq.1) then

         if(CurrentProcessorID.eq.0) write(*,*) 'noeuds equidistants'
         zii(1)=minval(t0)
         zii(nzz)=maxval(t3)
         pas =(zii(nzz)-zii(1))/dble(nzz-1)
!!c        if(CurrentProcessorID.eq.0) write(*,*)'pas=',pas
         j=1
         do while(j.lt.nzz-1)
            j=j+1
            zii(j)=zii(j-1)+pas
         end do
         if (j.ne.(nzz-1)) then
            if(CurrentProcessorID.eq.0) write(*,*)'probleme definition des noeuds'
            if(CurrentProcessorID.eq.0) write(*,*)'j+1=',j+1,'nz-1=',nzz
            stop
         end if
      end if

      end subroutine zi_noeuds

! }}}



!C------------------------------------------------------------
!C                     TIME_HMS
!C------------------------------------------------------------

! {{{
!C------------------------------------------------------------
!C      Procedure pour transformer un temps en seconde
!C      en un temps decompte en heures minutes secondes
!C------------------------------------------------------------

      subroutine time_hms(time,hour,minut,second)



      implicit none
      real ::time,hour,minut,second,remain,mo1,mo2

      remain=0
      hour=0
      minut=0
      second=0
      mo1=3600
      mo2=60

      remain=modulo(time,mo1)
      hour=(time-remain)/mo1

      second=modulo(remain,mo2)
      minut=(remain-second)/mo2

      end subroutine time_hms



! }}}



    
! =========================================================
!
!   Subroutine pour la fonction Beta
!
!
! =========================================================



! {{{

!C -------------------------------------------------------
!C              Densite d'une beta
!C--------------------------------------------------------




      Function beta_densite(X,a,b)


      implicit none
      double precision :: beta,a,b,gammln
      double precision:: beta_densite
      double precision:: X


      beta=dexp(gammln(a+b)-gammln(a)-gammln(b))

if (a.eq.1.d0.and.b.eq.1.d0) beta=1.d0
      beta_densite=((X)**(a-1))*((1-X)**(b-1))*beta
      return
      end


!C Calcul de gamma (gammln(a))


      Function gammln(xx)


      implicit none
!C retourne la valeur ln(gamma(xx)) pour xx>0

      double precision :: gammln,xx

      integer j
      double precision ser,stp,tmp,x,y,cof(6)
      save cof,stp

      data cof,stp/76.18009172947146d0,-86.50532032941677d0,     & 
      24.01409824083091d0, -1.231739572450155d0,.1208650973866179d-2,     & 
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

      return
      end



!C ---------------- CDF beta -------------------------------------


      Function betai(a,b,x,indi)

      double precision :: betai,a,b,x
      double precision :: bt,betaCF,gammln
      integer :: indi
  if(indi.eq.1)then
      if(x.lt.0..or.x.gt.1.) stop 'bad argument x in betai'
!C      if(CurrentProcessorID.eq.0) write(*,*)'in betai',a,b
      if (x.eq.0..or.X.eq.1.) then
         bt=0.
      else
         bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
      end if

      if (x.lt.(a+1.)/(a+B+2.)) then
         betai=bt*betaCF(a,b,x)/A
         return
      else
         betai=1.-bt*betaCF(b,a,1.-x)/b
         return
      end if

   else
     betai=x
   end if
end


!C betaCF est utilise par betai


      Function betaCF(a,b,x)


      implicit none
      integer ::maxit
      double precision ::betaCF,a,b,x,eps,fpmin
      parameter (maxit=100,eps=3.e-7,fpmin=1.e-30)
      integer ::m,m2
      double precision ::aa,c,d,del,h,qab,qam,qap
      qab=a+b
      qap=a+1
      qam=a-1
      c=1.       ! first step of Lentz's method
      d=1.-qab*x/qap
      if (abs(d).lt.fpmin) d=fpmin
      d=1./d
      h=d
      do m=1,maxit
         m2=2*m
         aa=m*(b-m)*x/((qam+m2)*(a+m2))
         d=1.+aa*d        ! one step (the even one) of the recurrence
         if (abs(d).lt.fpmin) d=fpmin
         c=1.+aa/c
         if (abs(c).lt.fpmin) c=fpmin
         d=1./d
         h=h*d*c
         aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
         d=1.+aa*d        ! next step of the recurrence (the odd one)
         if (abs(d).lt.fpmin) d=fpmin
         c=1.+aa/c
         if (abs(c).lt.fpmin) c=fpmin
         d=1./d
         del=d*c
         h=h*del
         if (abs(del-1.).lt.eps) goto 1
      end do
!      if(CurrentProcessorID.eq.0) write(*,*) 'a or b too big, or maxit too small in betaCF'
!      if(CurrentProcessorID.eq.0) write(*,*) 'a',a,'b',b
      stop

 1    betaCF=h
      return
      end

!C----------------------------------------------------------------



!C********************************************************************
!C            calcul d'une inverse de Beta incomplete
!C********************************************************************




      double precision function xinbta(p,q,beta,alpha,ifault)





      implicit double precision (a-h,o-z)


      double precision :: beta,alpha
      integer :: ifault
      double precision ::iex


!c     algorithm as 109 appl. statist. (1977), vol.26, no.1
!c     (replacing algorithm as 64  appl. statist. (1973),
!c     vol.22, no.3)
!c
!c     Remark AS R83 and the correction in vol40(1) p.236 have been
!c     incorporated in this version.
!c
!c     Computes inverse of the incomplete beta function
!c     ratio for given positive values of the arguments
!c     p and q, alpha between zero and one.
!c     log of complete beta function, beta, is assumed to be known.
!c
!c     Auxiliary function required: BETAIN = algorithm AS63
!c
      logical indx
!c
!c     Define accuracy and initialise.
!c     SAE below is the most negative decimal exponent which does not
!c     cause an underflow; a value of -308 or thereabouts will often be
!c     OK in double precision.
!c  variable SAE in XINBTA changed from -37D.0 to -308D.0 to avoid
!c  infinite loop (only a problem on digital unix).

!c
!C      data acu/1.0d-14/
!C      data SAE/-308.D0/
!C      data zero/0.0d0/, one/1.0d0/, two/2.0d0/
!C      data three/3.0d0/, four/4.0d0/, five/5.0d0/, six/6.0d0/



      double precision ::SAE=-308.D0,zero=0.0d0,one=1.0d0,two=2.0d0     &  
   ,three=3.0d0,four=4.0d0,five=5.0d0,six=6.0d0

      double precision ::a,pp,qq,p,q,y,r,t,s,h,w,yprev,sq,prev,ACU
      double precision ::betain,xin,g,adj,tx,fpu

      fpu = 10.d0 ** sae
      xinbta = alpha
!c
!c     test for admissibility of parameters
!c
      ifault = 1
      if (p.le.zero .or. q.le.zero) return
      ifault = 2
      if (alpha.lt.zero .or. alpha.gt.one) return
      ifault = 0
      if (alpha.eq.zero .or. alpha.eq.one) return
!c
!c     change tail if necessary
!c

      if (alpha.le.0.5d0) goto 1
      a = one-alpha
      pp = q
      qq = p
      indx = .true.
      goto 2
    1 a = alpha
      pp = p
      qq = q
      indx = .false.
!c
!c     calculate the initial approximation
!c
    2 r = dsqrt(-dlog(a*a))
      y = r-(2.30753d0+0.27061d0*r)/(one+(0.99229d0+0.04481d0*r)*r)
      if(pp.gt.one .and. qq.gt.one) goto 5
      r = qq+qq
      t = one/(9.0d0*qq)
      t = r*(one-t+y*dsqrt(t))**3
      if(t.le.zero) goto 3
      t = (four*pp+r-two)/t
      if(t.le.one) goto 4
      xinbta = one-two/(t+one)
      goto 6
    3 xinbta = one-dexp((dlog((one-a)*qq)+beta)/qq)
      goto 6
    4 xinbta = dexp((dlog(a*pp)+beta)/pp)
      goto 6
    5 r = (y*y-three)/six
      s = one/(pp+pp-one)
      t = one/(qq+qq-one)
      h = two/(s+t)
      w = y*dsqrt(h+r)/h-(t-s)*(r+five/six-two/(three*h))
      xinbta = pp/(pp+qq*dexp(w+w))

!c
!c     solve for x by a modified newton-raphson method,
!c     using the function betain
!c
    6 r = one-pp
      t = one-qq
      yprev = zero
      sq = one
      prev = one
      if(xinbta.lt.0.0001d0) xinbta = 0.0001d0
      if(xinbta.gt.0.9999d0) xinbta = 0.9999d0
      IEX = MAX(-5.D0/PP**2 - 1.D0/A**.2 - 13.D0, SAE)
      ACU = 10.D0 ** IEX

      ACU=1.0D-30


    7 y = betain(xinbta,pp,qq,beta,ifault)

      if(ifault.eq.0) goto 8
      ifault = 3
      return
    8 continue
      xin = xinbta
      y = (y-a)*exp(beta+r*log(xin)+t*log(one-xin))
      if(y*yprev.le.zero) prev = max(sq, fpu)
      g = one
    9 adj = g*y
      sq = adj*adj
      if(sq.ge.prev) goto 10
      tx = xinbta-adj
      if(tx.ge.zero .and. tx.le.one) goto 11
   10 g = g/three
      goto 9
   11 if(prev.le.acu) goto 12
      if(y*y.le.acu) goto 12
      if(tx.eq.zero .or. tx.eq.one) goto 10
      if(tx.eq.xinbta) goto 12
      xinbta = tx
      yprev = y
      goto 7
   12 if (indx) xinbta = one-xinbta
      return
      end



!C****************************************************************



      double precision function betain(x, p, q, beta, ifault)
      implicit double precision (a-h, o-z)
!c
!c     algorithm as 63  appl. statist. (1973), vol.22, no.3
!c
!c     computes incomplete beta function ratio for arguments
!c     x between zero and one, p and q positive.
!c     log of complete beta function, beta, is assumed to be known
!c
      logical indx
      integer :: ifault,ns
      double precision :: beta
!c
!c     define accuracy and initialise
!c

      double precision ::zero=0.0d0,one=1.0d0,acu=0.1d-14

      double precision :: x,p,q,psq,cx,xx,pp,qq,term,ai,rx,temp




      betain=x

!c     test for admissibility of arguments
!c
      ifault=1
      if(p.le.zero .or. q.le.zero) return
      ifault=2
      if(x.lt.zero .or. x.gt.one) return
      ifault=0
      if(x.eq.zero .or. x.eq. one) return
!c
!c     change tail if necessary and determine s
!c
      psq=p+q
      cx=one-x
      if(p.ge.psq*x) goto 1
      xx=cx
      cx=x
      pp=q
      qq=p
      indx=.true.
      goto 2
    1 xx=x
      pp=p
      qq=q
      indx=.false.
    2 term=one
      ai=one
      betain=one
      ns=qq+cx*psq
!c
!c     user soper's reduction formulae.
!c
      rx=xx/cx
    3 temp=qq-ai
      if(ns.eq.0) rx=xx
    4 term=term*temp*rx/(pp+ai)
      betain=betain+term
      temp=abs(term)
      if(temp.le.acu .and. temp.le.acu*betain) goto 5
      ai=ai+one
      ns=ns-1
      if(ns.ge.0) goto 3
      temp=psq
      psq=psq+one
      goto 4
!cc     calculate result
!c
    5 betain=betain*exp(pp*log(xx)+(qq-one)*log(cx)-beta)/pp
      if(indx) betain=one-betain
      return
      end


!C****************************************************************

      DOUBLE PRECISION function beta_ln(z,w)

      implicit none
      double precision :: z,w
      double precision :: gammln

      beta_ln=gammln(z)+gammln(w)-gammln(z+w)
      return
      end function beta_ln

! }}}

!C------------------------------------------------------------
!C     
!                 POSTPROBMULT
!
!C------------------------------------------------------------


! {{{

!C-------------------------------------------------------------
!c
!c          Subroutine pour calculer les
!c      probabilites a posteriori de suivre chacune
!c      des composantes g pour chacun des sujets i
!c     sachant les informations de Y_i et D_i
!c     Modifie pour smooth hazard (en cours)
!c
!c-------------------------------------------------------------
      subroutine postprobmult_beta(b,npm,PPI,ppitest,probapost,scorefinalTT,var_scoreTT,scorefinalYT,var_scoreYT,n_frailty)!,P_AI,P_I)

!-------------------------------------------------------------
!PPI(i,g)=P(ci=g|Yi,Ti,Xi,theta)
!ppitest(i,g)=P(ci=g|Yi,Xi,theta)
!P_AI(i,g)=P(ci=g|Yi,u_i,Xi,theta)
!P_I(i,g)=P(ci=g|Xi,theta)
!U_i = E(Ui|Yi,Ti) bayesian estimates
!-------------------------------------------------------------
        use commun
        use donnees_indiv
        use optim_parallele

        implicit none
	integer,intent(in)::npm,n_frailty
        integer::ier,nnn,nmoins,ntrtot,nrisq_curr
        integer :: i,j,k,l,jj,kk,q,m,g,ii,ll
        integer::nmestot,nxevtcurr
        integer :: h,jj1,q1,k1,jjj

        double precision::Y5,Y5_2,eps,temp,det,det2,f,f1,fi_scoreYTf_a,ytemp,aa,bb,betai,fevt,aa1,bb1,vraisobs,cc1,dd1,som &
             ,surv_glob,survint_glob,varexpsurv,f_score,SS


        !C --------------   dimensions revues
        double precision,dimension(npm)::b
	    double precision,dimension(maxmestot,npg)::mu_ui
        double precision,dimension(nvc)::bh1
        double precision,dimension(nrisq)::brisq
        double precision,dimension(ng)::pi
        double precision,dimension(nvarprob)::bprob,Xprob
        double precision,dimension(nxevt)::bevt,Xevt
        double precision,dimension(-2:ntrspl-3)::splaa
        double precision,dimension(nea,nea)::Ut,Ut2
        double precision,dimension(ng,nevt)::surv0,surv,risq,survint
        double precision,dimension(ng)::vectretard,vectvalfin,scoreTT,scoreYT,SUTi0TT,STi0TT
        double precision,dimension(1)::bevtint
        double precision,dimension(nvarglobssG)::b3
        double precision,dimension(2*np+1)::b0
        double precision,dimension(ncontdep+ncont)::sum
        double precision,dimension(nvdepssG)::bdepssG
        double precision,dimension(nvdepG)::bdepG
        double precision,dimension(nvarglobG)::b2
        double precision,dimension((ncontdep+ncont)*nq)::b4
        double precision,dimension(4)::tacti
        double precision,dimension(sumidxevt)::vexpiact
        double precision::retard2,STi0,SUTi0,retard
        integer::kkk
        double precision,dimension(nvarxevt)::bvexp
        double precision,dimension(maxmestot)::Y2,Y3,Y4,Y4_2,VV,YY
        double precision,dimension(maxmestot,2*np+1)::Z
        double precision,dimension(maxmestot,nvarglobssG)::XglobssG
        double precision,dimension(maxmestot,nvdepssG)::XdepssG
        double precision,dimension(maxmestot,(ncontdep+ncont)*nq)::Xcont
        double precision,dimension(maxmestot,nvarglobG)::XglobG
        double precision,dimension(maxmestot,nvdepG)::XdepG
        double precision,dimension(maxmestot,nea)::P
        double precision,dimension(nea,maxmestot)::Uii,ZZ
        double precision,dimension(maxmestot,maxmestot)::VC,VCAR!,VC2

        double precision,dimension(maxmestot*(maxmestot+1)/2)::Vi
        double precision,dimension(ng)::fi,fi1,fi_scoreTT,fi_scoreYT
        double precision,dimension(nea)::f_scoreYT
        double precision,dimension(ng),intent(inout)::probapost
        double precision::valfin,scoregTT,scoregYT,scorefinTT,scorefinYT
	double precision,intent(out)::scorefinalTT,var_scoreTT
        double precision,dimension(nea*nq,ng)::ui_hat
	double precision,intent(out),dimension(nea)::scorefinalYT
	double precision,intent(out),dimension(nea,nea)::var_scoreYT
		
        !C NON REVU

        double precision,dimension(ns,ng)::PPI
	double precision,dimension(ns,ng),intent(inout)::ppitest!,U_i!,P_AI,P_I
        integer::numact

        if(CurrentProcessorID.eq.0) write(*,*)typrisq,"typrisq"

	ui_hat=0.d0
        eps=1.d-10
	scorefinalTT=0.d0
	var_scoreTT=0.d0
	scorefinalYT=0.d0
	var_scoreYT=0.d0
	
        DO i=1,ns
           DO g=1,ng
              PPI(i,g)=0.d0
              ppitest(i,g)=0.d0
           END DO
        END DO

        !C-------------------- Nombre d'ea -----------------------------

        !C-------------------- Nombre d'ea -----------------------------

        !      IF(CURRENTPROCESSORID.EQ.0) WRITE(*,*) 'npmtot',npmtot,'nef',nef


        k=0
        l=0
        DO j=1,npmtot
           if (ide(j).eq.1) then
              k=k+1
              b1(j)=b(k)
           else
              l=l+1
              b1(j)=bne(l)
           end if

        END DO

        if (k.ne.npm) then
           if(CurrentProcessorID.eq.0) write(*,*) 'problem with number of parameters'
           if(CurrentProcessorID.eq.0) write(*,*)'entre',npm,'calcule',k

           stop
        end if

        !C-------------------- Vecteur de probas -----------------------




        !C ------ Ut : ---------------------------------------------

        !C------- Creation de Ut la transformee de Cholesky ------

        Do j=1,nvc
           bh1(j)=b1(nef+j)
        End do

        !C     if(CurrentProcessorID.eq.0) write(*,*)'bh1',(bh1(k),k=1,nvc)

        Ut=0.d0
        If (idiag.eq.1) then
           do j=1,nea
              do k=1,nea
                 if (j.eq.k) then
                    Ut(j,k)=bh1(j)
                 end if
              end do
           end do
        else if (idiag.eq.0) then
           do j=1,nea
              do k=j,nea
                 Ut(k,j)=bh1(j+k*(k-1)/2)
              end do
           end do
        end if

        !C calcul de brisq a chaque composante et risq, surv et surv0 pour tous les i et tous les g

        !C------ boucle sur les sujets i-----------------------------
2222    FORMAT(I5,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X, &
             F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6)
        open(2,file='prior_class_probs.txt')


        !entretard=0.d0
        !jacobien=0.d0
        !vrais=0.d0
        !vrais_survie=0.d0
        !risq=0.d0
        !surv=0.d0
        !surv0=0.d0
        !survint=0.d0
        !brisq=0.d0

        !####################################################
        !   si evt a eu lieu
        if (evt.ne.0) then
           nrisq_curr=0
           numact=1
           !numact permet le remplissage de brisq
           !on boucle sur le nombre d'evt
           do kk=1,nevt
              !boucle sur le nombre de classes
              if (risqcom(kk).eq.0) then
                 do g=1,ng
                    !boucle sur le nombre de parametres par evt
                    do k=1,nprisq(kk)
                       !positivation carre ou expo
                       if (paramsurv.eq.0) then 
                          !risque commun/specifique/proportionnel
                          if (risqcom(kk).eq.2) then
                             brisq(numact)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)    
                             numact=numact+1
                          end if
                          if (risqcom(kk).eq.1) then
                             brisq(numact)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)
                             numact=numact+1
                          end if
                          if (risqcom(kk).eq.0) then
                             brisq(numact)=b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)&
                                  &+nrisq_curr+k)  !commente le carre pour tester!!
                             numact=numact+1
                          end if
                       else 
                          if (risqcom(kk).eq.2) then
                             brisq(numact)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))    
                             numact=numact+1
                          end if
                          if (risqcom(kk).eq.1) then
                             brisq(numact)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))
                             numact=numact+1
                          end if
                          if (risqcom(kk).eq.0) then
                             brisq(numact)=exp(b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k))
                             numact=numact+1
                          end if
                       end IF !fin paramsurv==0
                    end do !fin boucle nb param par evt

                    !brisq correctement rempli, sinon recouvrement indiçage par k des valeurs n'etaient pas remplies
                    !appel de la fonction de calcul des risques
            
                 end do !boucle sur les classes
              endif

              if (risqcom(kk).ge.1) then
                 !auparavant les variables explicatrices etaient ajoutees apres la fonction, ce n'est plus le cas ici
                 !brisq doit donc contenir les coef des varexp
                 do k=1,nprisq(kk)
                    !positivation carre ou expo
                    if (paramsurv.eq.0) then 
                       !risque commun/specifique/proportionnel
                       if (risqcom(kk).eq.2) then
                          brisq(numact)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)    
                          numact=numact+1
                       end if
                       if (risqcom(kk).eq.1) then
                          brisq(numact)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)
                          numact=numact+1
                       end if
                       if (risqcom(kk).eq.0) then
                          brisq(numact)=b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)&
                               &+nrisq_curr+k)
                          numact=numact+1
                       end if
                    else 
                       if (risqcom(kk).eq.2) then
                          brisq(numact)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))    
                          numact=numact+1
                       end if
                       if (risqcom(kk).eq.1) then
                          brisq(numact)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))
                          numact=numact+1
                       end if
                       if (risqcom(kk).eq.0) then
                          brisq(numact)=exp(b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k))
                          numact=numact+1
                       end if
                    end IF !fin paramsurv==0
                 end do !fin boucle nb param par evt
              endif

              !si weib calcul des valeurs de survie/risque
              ! if (risqcom(kk).eq.2.and.g.lt.ng.and.ng.gt.1) then               
              !                   risq(g,kk)=risq(g,kk)*exp(b1(nvarprob*(ng-1)+nrisq_curr+nprisq(kk)+g))

              !                   surv(g,kk)=surv(g,kk)*exp(b1(nvarprob*(ng-1)+nrisq_curr+nprisq(kk)+g))
              !                   survint(g,kk)=survint(g,kk)*exp(b1(nvarprob*(ng-1)+nrisq_curr+nprisq(kk)+g))
              !                   surv0(g,kk)=surv0(g,kk)*exp(b1(nvarprob*(ng-1)+nrisq_curr+nprisq(kk)+g))
              !             end if 
              nrisq_curr=nrisq_curr+nrisqspec(kk)
           end do

           nrisq_curr=0
           do kk=1,nevt
              if (risqcom(kk).eq.2) then
                 nrisq_curr=nrisq_curr+nrisqspec(kk)-ng+1
                 do g=1,ng-1
                    brisq(numact)=b1(nvarprob*(ng-1)+nrisq_curr+g)
                    numact=numact+1
                 enddo
                 nrisq_curr=nrisq_curr+ng-1
              else
                 nrisq_curr=nrisq_curr+nrisqspec(kk)
              endif
           enddo
        endif

        risq=0.d0
        surv=0.d0
        surv0=0.d0
        survint=0.d0

	bvexp=0.d0
	numact=1
        do l=1,nvarxevt
           bvexp(numact)=b1((ng-1)*nvarprob+nrisq+numact)
           numact=numact+1
        enddo
		
	Y1all=0.d0
        DO i=1,ns

           numpat=i

           ! creation vecteur de survie



           !C calcul de brisq a chaque composante et risq, surv et surv0 pour tous les i et tous les g

         
           !   if (evt.ne.0) then
           !               nrisq_curr=0
           !               do kk=1,nevt
           !                  do g=1,ng
           !                     brisq=0.d0
           !                     do k=1,nprisq(kk)
           !                        if (paramsurv.eq.0) then
           !                           if (risqcom(kk).eq.2) then
           !                              brisq(k)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)    
           !                           end if
           !                           if (risqcom(kk).eq.1) then
           !                              brisq(k)=b1(nvarprob*(ng-1)+nrisq_curr+k)*b1(nvarprob*(ng-1)+nrisq_curr+k)
           !                           end if
           !                           if (risqcom(kk).eq.0) then
           !                              brisq(k)=b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k)*&
           !                                   b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k)
           !                           end if
           !                        ELSE
           !                           if (risqcom(kk).eq.2) then
           !                              brisq(k)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))    
           !                           end if
           !                           if (risqcom(kk).eq.1) then
           !                              brisq(k)=exp(b1(nvarprob*(ng-1)+nrisq_curr+k))
           !                           end if
           !                           if (risqcom(kk).eq.0) then
           !                              brisq(k)=exp(b1(nvarprob*(ng-1)+nprisq(kk)*(g-1)+nrisq_curr+k))
           !                           end if
           !                        end if

           !                     end do


           vexpiact=0.d0
           evtact=1
           kkk=1
           do evtact=1,nevt
              do k=1,nv
                 if (idxevt(evtact*nv-nv+k).ge.1) then
                    vexpiact(kkk)=X(i,k)
                    kkk=kkk+1
                 endif
              enddo
           enddo


           tacti=0.d0
           tacti(1)=t0(i)
           tacti(2)=t1(i)
           tacti(3)=t2(i)
           tacti(4)=t3(i)

           if (evt.ne.0) then
              do g=1,ng
                 call vraisSmooth(i,contevt+nva01+nva02+nva12,brisq,bvexp,g,vexpiact,tacti,valfin,retard2,scorefinTT,scorefinYT,&
							&SUTi0,STi0,1,n_frailty,0)
                 vectvalfin(g)=valfin ! vecteur pour un individu sa contribution en fonction de la classe
                 vectretard(g)=exp(-retard2)
                 scoreTT(g)=scorefinTT
                 scoreYT(g)=scorefinYT
                 SUTi0TT(g)=SUTi0
                 STi0TT(g)=STi0
              enddo
           end if !fin condition evt != 0



           !C creation du vecteur de variables explicatives dep du temps : XdepssG sans mixture
	   nmestot=0
	   do q=1,nq
		do j=1,nmes(i,q)
		 nmestot=nmestot+1
		end do
	  end do

           XdepssG=0.d0
           l=0
           Do k=1,nvdep
              If (idvdep(k).eq.1) then
                 l=l+1
                 jj=0
                 do q=1,nq
                    do j=1,nmes(i,q)
                       jj=jj+1
                       XdepssG(jj,l)=Xdep(Navant_i(i)+jj,k)
                    end do
                 End do
                 if (jj.ne.nmestot) then
                    if(CurrentProcessorID.eq.0) write(*,*)'probleme XdepssG(jj,k)'
                    stop
                 end if
              end if
           end do



           XdepG=0.d0
           l=0
           Do k=1,nvdep
              If (idvdep(k).eq.2) then
                 l=l+1
                 jj=0
                 do q=1,nq
                    do j=1,nmes(i,q)
                       jj=jj+1
                       XdepG(jj,l)=Xdep(Navant_i(i)+jj,k)
                    end do
                 End do
                 if (jj.ne.nmestot) then
                    if(CurrentProcessorID.eq.0) write(*,*)'probleme XdepG(jj,k)'
                    stop
                 end if
              end if
           end do

           !C------ creation de XglobG(j,k) (avec mixture) ----------

           XglobG=0.d0

           l=0
           Do k=1,nv
              If (idglob(k).eq.1.and.idg(k).eq.1) then
                 l=l+1
                 Do j=1,nmestot
                    XglobG(j,l)=X(i,k)
                 End do

              Else if (idglob(k).eq.2.and.idg(k).eq.1) then
                 jj=0
                 l=l+1
                 Do q=1,nq
                    Do j=1,nmes(i,q)
                       jj=jj+1
                       Do m=0,np

                          if (m.ge.1) then
                             XglobG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)*(0.1**(m-1))

                          else
                             XglobG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)
                          end if


                       End do

                    End do
                 End do
                 l=l+np
              End if
           End do


           if (nvarglobssG+nvarglobG*ng.ne.nvarglob) then
              if(CurrentProcessorID.eq.0) write(*,*) 'pb : def nvarglobssG/G'
              stop
           end if


           !C------ creation de Xcont(q,k,j) -----------------

           !C creation des variables pour les contrastes

           Xcont=0.d0
           nmoins=0
           do q=1,nq
              do j=1,nmes(i,q)
                 h=0
                 do k=1,nvdep
                    if (idvdep(k).eq.3.or.idvdep(k).eq.4) then
                       h=h+1
                       Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)=Xdep(Navant_i(i)+nmoins+j,k)
                    end if
                 end do
                 h=0
                 do k=1,nv
                    if (idtest(k).eq.1) then
                       h=h+1
                       Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)=X(i,k)
                    end if
                    if (idtest(k).eq.2) then
                       Do m=0,np
                          h=h+1
                          if (m.ge.1) then
                             Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)     &
                                  =X(i,k)*(t(Navant_i(i)+nmoins+j)**m)*(0.1**(m-1))

                          else
                             Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)     & 
                                  =X(i,k)*(t(Navant_i(i)+nmoins+j)**m)
                          end if
                       End do
                    end if
                 end do
                 if (h.ne.ncont) then
                    if(CurrentProcessorID.eq.0) write(*,*)'h=',h,'problem in cont 2'
                    stop
                 end if
              End do
              nmoins=nmoins+nmes(i,q)
           end do

           if (nmoins.ne.nmestot) then
              if(CurrentProcessorID.eq.0) write(*,*)'problem in cont'
              stop
           end if


           !C creer Xevt:



!            if (evt.ne.0) then
!               Xevt=0.d0
!               l=0
!               do kk=1,NEVT
!                  do k=1,nv
!                     if (idxevt((kk-1)*nv+k).eq.1.or.idxevt((kk-1)*nv+k).eq.2) then
!                        l=l+1
!                        Xevt(l)=X(i,k)
!                     end if
!                  end do
!               end do
!            end if


           VCAR=0.d0
           if (iAR.gt.0) then
              jj=0
              Do q=1,nq
                 jj1=0
                 Do q1=1,nq

                    Do k=1,nmes(i,q)

                       do k1=1,nmes(i,q1)


                          if (iAR.eq.1) then

                             VCAR(jj+k,jj1+k1)=(b1(nef+nvc+ng-1+nq+1)**2)*min(t(Navant_i(i)+jj+k),t(Navant_i(i)+jj1+k1))

                          else if (iAR.eq.2) then

                             VCAR(jj+k,jj1+k1)=(b1(nef+nvc+ng-1+nq+1)**2) & 
                                  *exp(-b1(nef+nvc+ng-1+nq+2)*abs(t(Navant_i(i)+jj+k)-t(Navant_i(i)+jj1+k1)))

                          end if

                       end do

                    end do

                    jj1=jj1+nmes(i,q1)

                 end do

                 jj=jj+nmes(i,q)

              end do


           end if



           !C------- Creation du vecteur Y1 --------------------------
           !C------------ et des matrices Kmat et Sigmaq ------



           Y1=0.d0
           Sigmaq=0.d0
           Vw=0.d0
           jj=0.d0
           jj1=0
           ntrtot=0
           Do q=1,nq
              if (idord(q).gt.1) then
                 ntrtot=ntrtot+idord(q)-1
                 jj1=jj1+nmes(i,q)               
              else if (idord(q).eq.0) then
                 aa1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &     
                      +ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot))/(1+exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)     & 
                      +nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot)))
                 bb1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &
                      +ncont*(nq-1)+2+ntrtot))/(1.d0+exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &       
                      +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+2+ntrtot)) )
                 !C               bb1=bb1/4.d0
                 bb1=aa1*(1.d0-aa1)*bb1
                 cc1=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     & 
                      +ncont*(nq-1)+3+ntrtot)
                 dd1=abs(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &   
                      +ncont*(nq-1)+4+ntrtot))


                 aa=aa1*aa1*(1-aa1)/bb1-aa1
                 bb=aa*(1-aa1)/aa1
                 !C			   if(CurrentProcessorID.eq.0) write(*,*)'aa-dd',aa,bb,cc1,dd1

                 if (aa.ne.aa.or.bb.ne.bb) then
                    if(CurrentProcessorID.eq.0) write(*,*)'aa',aa,'bb',bb,'pb in funcpa'
                    if(CurrentProcessorID.eq.0) write(*,*)'b1',b1
                    stop
                 end If

                 Do k=1,nmes(i,q)
                    do l=1,nmes(i,q)
                       Vw(jj+k,jj+l)=b1(nef+nvc+ng-1+q)**2
                       if(k.eq.l) then
                          ytemp=(Y(Navant_i(i)+jj1+k)-range(1+(q-1)*2))/(range(2+(q-1)*2)-range(1+(q-1)*2))
                          Y1(jj+k)=(betai(aa,bb,ytemp,tbeta)-cc1)/dd1
                          !C     OK pour la fonction de repartition
                          Sigmaq(jj+k,jj+l)=b1(nef+nvc+ng-1+nq+nAR+q)**2
                       end if
                    end do
                 end do
                 ntrtot=ntrtot+ntr
                 jj=jj+nmes(i,q)
                 jj1=jj1+nmes(i,q)

              elseif (idord(q).eq.1) then

                 bb=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &        
                      +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)

                 do kk=2,ntrspl
                    splaa(kk-3)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &       
                         +ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot)*b1(nvarprob*(ng-1)+nrisq     &       
                         +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot)
                    !C     aa(kk-3)=exp(b1(kk+(q-1)*ntr))
                    !C     if(CurrentProcessorID.eq.0) write(*,*)'aa(kk-3)',aa(kk-3)
                 end do

                 Do k=1,nmes(i,q)

                    do l=1,nmes(i,q)
                       Vw(jj+k,jj+l)=b1(nef+nvc+ng-1+q)**2

                       if(k.eq.l) then

                          !C     creation du vecteur de donnees transformees

                          !C     creation du vecteur de donnees transformees
                          ll=0
                          if (Y(Navant_i(i)+jj1+k).eq.zitr(nztr,q)) then
                             ll=nztr-1
                          end if
                          som=0.d0
                          do kk = 2,nztr

                             if ((Y(Navant_i(i)+jj1+k).ge.zitr(kk-1,q)).and.(Y(Navant_i(i)+jj1+k).lt.zitr(kk,q)))then
                                ll=kk-1
                             end if
                          end do
                          if (ll.lt.1.or.ll.gt.nztr-1) then
                             if(CurrentProcessorID.eq.0) write(*,*) 'probleme dans funcpa splines'
                             if(CurrentProcessorID.eq.0) write(*,*) 'll=',ll,'Y=',Y(Navant_i(i)+jj1+k)
                             stop
                          end if
                          if (ll.gt.1) then
                             do ii=2,ll
                                som=som+splaa(ii-3)
                             end do
                          end if

                          Y1(jj+k)=bb+ som +splaa(ll-2)*im2(Navant_i(i)+jj1+k)+splaa(ll-1)   &
                               *im1(Navant_i(i)+jj1+k)+splaa(ll)*im(Navant_i(i)+jj1+k)
                          Sigmaq(jj+k,jj+l)=b1(nef+nvc+ng-1+nq+nAR+q)**2

                       end if
                    end do
                 end do

                 ntrtot=ntrtot+ntrspl
                 jj=jj+nmes(i,q)
                 jj1=jj1+nmes(i,q)

              else if (idord(q).eq.-1) then
                 if(CurrentProcessorID.eq.0) write(*,*)'a faire postprob pour tronque'
                 stop
              end if

		do k=1,nmes(i,q)
			Y1all(ns*(q-1)+i,k)=Y1(k)
		end do
           End do

           !C------ Vecteurs de parms sans mixtures -------------------

           bdepssG=0.d0
           jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)

           l=0
           kk=0
           Do k=1,nvdep
              If (idvdep(k).eq.1) then
                 l=l+1
                 kk=kk+1
                 bdepssG(l)=b1(jj+kk)
              end if
              If (idvdep(k).eq.2) then
                 kk=kk+ng
              end if
           end do
           if (l.ne.nvdepssG) then
              if(CurrentProcessorID.eq.0) write(*,*) 'probleme longueur vecture bdepssG'
              stop
           end if



      b3=0.d0
      b4=0.d0
      jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+(2*np+1)*ng+nprmdep
      l=0
      kk=0
      Do k=1,nv
         If (idglob(k).eq.1.and.idg(k).eq.0) then
            l=l+1
            kk=kk+1
            b3(l)=b1(jj+kk)  ! correspond à beta
         Else if (idglob(k).eq.2.and.idg(k).eq.0) then
            l=l+1
            kk=kk+1
            Do m=0,np
               b3(l+m)=b1(jj+kk+m)
               b3(l+m+1)=b1(jj+kk+m+1)
            End do
            l=l+np*2
            kk=kk+np*2
         Else if (idglob(k).eq.1.and.idg(k).eq.1) then
            kk=kk+ng
         Else if (idglob(k).eq.2.and.idg(k).eq.1) then
            kk=kk+ng*(np+1)
         End if
      End do

           !C        if(CurrentProcessorID.eq.0) write(*,*)'verif : kk=nvarglobssG+nvarglobG*ng',
           !C     & kk,nvarglobssG+nvarglobG*ng          !OK
           !C        if(CurrentProcessorID.eq.0) write(*,*)'verif : l=nvarglobssG',l,nvarglobssG      !OK



           !C     b4 : contrastes sur les tests

           !C     b4 : contrastes sur les tests
           b4=0.d0
           sum=0.d0
           do j=1,ncontdep
              do q=1,nq-1
                 b4((ncont+ncontdep)*(q-1)+j)=b1(nvarprob*(ng-1)+nrisq     &   
                      +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+(j-1)*(nq-1)+q)
                 sum(j)=sum(j)+b4((ncont+ncontdep)*(q-1)+j)
              end do
              b4((ncont+ncontdep)*(nq-1)+j)=-sum(j)
           end do
           do j=1,ncont
              do q=1,nq-1
                 b4((ncont+ncontdep)*(q-1)+ncontdep+j)=b1(nvarprob*(ng-1)     &  
                      +nrisq +nvarxevt +ng+ng*(2*np+1)+nprmdep+nvarglob     &      
                      + ncontdep*(nq-1)+(j-1)*(nq-1)+q)
                 sum(ncontdep+j)=sum(ncontdep+j)+b4((ncont+ncontdep)*(q-1)+ncontdep+j)
              end do
              b4((ncont+ncontdep)*(nq-1)+ncontdep+j)=-sum(ncontdep+j)
           end do

           !C------ contributions -------------------------------------

           f=0.d0
!           f_a=0.d0

           DO g=1,ng
              fi(g)=0.d0
!	      f_ai(g)=0.d0
           END DO

           Xprob=0.d0
           Xprob(1)=1.d0
           l=0
           do k=1,nv
              if (idprob(k).eq.1) then
                 l=l+1
                 Xprob(1+l)=X(i,k)
              end if
           end do

           if (l+1.ne.nvarprob) then
              if(CurrentProcessorID.eq.0) write(*,*)'problem nvarprob funcpa'
              stop
           end if

           pi=0.d0
           temp=0.d0
           Do g=1,ng-1

              do k=1,nvarprob
                 bprob(k)=b1((k-1)*(ng-1)+g)
              end do


              temp=temp+exp(DOT_PRODUCT(bprob,Xprob))

              pi(g)=exp(DOT_PRODUCT(bprob,Xprob))

              if ((pi(g)+1.d0).eq.pi(g)) then
                 if(CurrentProcessorID.eq.0) write(*,*)'g',g,'pi(g)',pi(g)
                 if(CurrentProcessorID.eq.0) write(*,*)'bprob',(bprob(j),j=1,nvarprob)
                 if(CurrentProcessorID.eq.0) write(*,*)'Xprob',(Xprob(j),j=1,nvarprob)
                 stop
              end if
           end do

           pi(ng)=1/(1+temp)

           Do g=1,ng


              if (g.ne.ng) then
                 pi(g)=pi(g)*pi(ng)
              end if

              !C     calcul de la survie a posteriori

              !C            if (evt.eq.1) then
              !C               Tpred(i)=pi(g)*exp(-surv(i,g))
              !C            end if

              !C bevt


!               if (evt.ne.0)then
!                  if(nxevt.ne.0) THEN
!                     bevt=0.d0
!                     m=0
!                     l=1
!                     do kk=1,nevt
!                        do k=1,nv
!                           if (idxevt((kk-1)*nv+k).eq.1) then
!                              bevt(l)=b1(nvarprob*(ng-1)+nrisq+m+1)
!                              l=l+1
!                              m=m+1
!                           end if
!                           if (idxevt((kk-1)*nv+k).eq.2) then
!                              bevt(l)=b1(nvarprob*(ng-1)+nrisq+m+g)
!                              l=l+1
!                              m=m+ng
!                           end if
!                        end do
!                     end do
!                  end IF

!                  bevtint=0.d0 
!                  if (nvdepsurv.ne.0) then
!                     bevtint(1)=b1(nvarprob*(ng-1)+nrisq+nvarxevt)
!                  end if
!               end if



           !C------ creation de la matrice des puissances du temps -----
             Z=0.d0
              jj=0

              Do q=1,nq
                 Do j=1,nmes(i,q)
                    jj=jj+1
		    jjj=1
                    l=0

		    Z(jj,jjj)=1.d0
		    jjj=jjj+1
		
                    Do k=2,np+1
			    Z(jj,jjj)=(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)**(k-1)
			    temp=dble(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)
			    !Z(jj,jjj+1)=(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)**(k-1)*temp
			    temp=(temp**2+trn)**0.5d0
			    Z(jj,jjj+1)=temp!

		            If (k.gt.2) then
			!	Z(jj,jjj)=Z(jj,jjj)*(0.1)**(k-2)
		!		Z(jj,jjj+1)=Z(jj,jjj+1)*(0.1)**(k-2)*temp
		                !Z(jj,k)=Z(jj,k)*(0.1)**(k-2)!! a changer !!
				print*,'pas code pour np=2'
				stop
		            end if

			    jjj=jjj+2
                    End do
                 End do
              End do
      	      nmestot=jj

              Z1=0.d0
	      jj=0
              Do q=1,nq
                 Do j=1,nmes(i,q)
                    jj=jj+1
                    l=0
                    Do k=1,np+1
		       If(idea(k).eq.1) then
		          l=l+1
		          if(k.eq.1) Z1(jj,l)=Z(jj,k) 
			  if(k.eq.2)then
				Z1(jj,l)=Z(jj,k)/2.0d0*(1.0d0-(Z(jj,k)**2+trn)**0.5d0/Z(jj,k))
			  	Z1(jj,l+1)=Z(jj,k)/2.0d0*(1.0d0+(Z(jj,k)**2+trn)**0.5d0/Z(jj,k))

			  end if
		       end if
                    End do
                 End do
              End do

           !C------ creation de XglobssG(j,k) (sans mixture) ----------

              XglobssG=0.d0

              l=0
              Do k=1,nv
                 If (idglob(k).eq.1.and.idg(k).eq.0) then
                    l=l+1
                    Do j=1,nmestot
                       XglobssG(j,l)=X(i,k)
                    End do
                 Else if (idglob(k).eq.2.and.idg(k).eq.0) then
                    jj=0
                    l=l+1
                    Do q=1,nq
                       Do j=1,nmes(i,q)
                          jj=jj+1

                          Do m=0,np
   			     temp=(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0
                             !C Attention : changer les Xglob avec interaction
                             if (m.ge.1) then
                                XglobssG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)-temp)!((t(Navant_i(i)+jj)-temp)**m)*(0.1**(m-1))
			        temp=dble(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)
			        temp=(temp**2+trn)**0.5d0
			        XglobssG(jj,l+m+1)=X(i,k)*temp!((t(Navant_i(i)+jj)-(b1((nvarprob)*&
						!(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)**m)*(0.1**(m-1))*temp
                             else
                                XglobssG(jj,l+m)=X(i,k)
                             end if

                          End do

                       End do
                    End do
                    l=l+2*np
                 End if
              End do

           !C     if(CurrentProcessorID.eq.0) write(*,*)'nvarglobssG (funcpa)',nvarglobssG

           !C     if(CurrentProcessorID.eq.0) write(*,*)'XglobssG(jj,k)'
           !C     do j=1,nmestot
           !C     if(CurrentProcessorID.eq.0) write(*,*)(XglobssG(j,k),k=1,nvarglobssG)
           !C     end do

           !C creation du vecteur de variables explicatives dep du temps : XdepG avec mixture

              !C------- creation de Vi ---------------------------------

              VC=0.d0
              P=0.d0
              Ut1=0.d0

              if (g.eq.ng) then
                 Ut1=Ut
              else
                 Ut1=Ut*(b1(nef+nvc+g))
              end if

              !C----- parametres avec mixtures ---------------------------

              b0=0.d0
              Do j=1,(2*np+1)
                 b0(j)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+(j-1)*ng+g)
              End do

              b2=0.d0
              l=0
              kk=0
              Do k=1,nv

                 If (idglob(k).eq.1.and.idg(k).eq.1) then
                    l=l+1
                    kk=kk+g

                    b2(l)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+kk) ! correspond a beta(g)

                 Else if (idglob(k).eq.2.and.idg(k).eq.1) then
                    l=l+1
                    kk=kk+g
                    Do m=1,np+1
                       b2(l+m-1)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+(m-1)*ng+kk)
                    End do
                    l=l+np
                    kk=kk+np*ng

                 Else if (idglob(k).eq.1.and.idg(k).eq.0) then
                    kk=kk+1

                 Else if (idglob(k).eq.2.and.idg(k).eq.0) then
                    kk=kk+(np+1)

                 End if

              End do

              ier=0

              !C------ calcul de l'expo ---------------------------------
              Y2=0.d0
              Y3=0.d0
              Y2=MATMUL(Z,b0)+MATMUL(XglobG,b2)+MATMUL(XdepssG,bdepssG)
              Y2=Y2+MATMUL(XglobssG,b3)+MATMUL(XdepG,bdepG)

              mu=0.d0

              mu=Y2+MATMUL(Xcont,b4)


              !C	        if(CurrentProcessorID.eq.0) write(*,*)'mu',(mu(j),j=1,nmestot)


              IF (QUECONT.eq.1) THEN

                 P=MATMUL(Z1,Ut1)
                 VC=MATMUL(P,transpose(P))+Sigmaq+Vw+VCAR
!		 VC2=Sigmaq+Vw+VCAR

                 !c     Vi en vecteur
		 Vi=0.d0
                 jj=0
                 do j=1,nmestot
                    do k=j,nmestot
                       jj=j+k*(k-1)/2
                       Vi(jj)=VC(j,k)
                    end do
                 end do
                 CALL DSINV(Vi,nmestot,eps,ier,det,0)
                 if (ier.eq.-1) then
                    !  if(CurrentProcessorID.eq.0) write(*,*)'dsinv vi ier',ier,' i=',i,' nmes(i)=',nmestot,'num',numero(i)
                    !  if(CurrentProcessorID.eq.0) write(*,*)'b1', (b1(k),k=1,npm)
                    nnn=nmestot*(nmestot+1)/2
                    !  if(CurrentProcessorID.eq.0) write(*,*)'Vi', (Vi(k),k=1,nnn)
                    stop
                 end if

                 !c     retransformation du vecteur Vi en matrice :
                 VC=0.d0
                 do j=1,nmestot
                    do k=1,nmestot
                       if (k.ge.j) then
                          VC(j,k)=Vi(j+k*(k-1)/2)
                       else
                          VC(j,k)=Vi(k+j*(j-1)/2)
                       end if
                    end do
                 end do


!                 jj=0
!                 do j=1,nmestot
!                    do k=j,nmestot
!                       jj=j+k*(k-1)/2
!                       Vi(jj)=VC2(j,k)
!                    end do
!                 end do
!                 CALL DSINV(Vi,nmestot,eps,ier,det2)
!                 if (ier.eq.-1) then
                    !  if(CurrentProcessorID.eq.0) write(*,*)'dsinv vi ier',ier,' i=',i,' nmes(i)=',nmestot,'num',numero(i)
                    !  if(CurrentProcessorID.eq.0) write(*,*)'b1', (b1(k),k=1,npm)
!                    nnn=nmestot*(nmestot+1)/2
                    !  if(CurrentProcessorID.eq.0) write(*,*)'Vi', (Vi(k),k=1,nnn)
!                    stop
!                 end if

                 !c     retransformation du vecteur Vi en matrice :
!                 VC2=0.d0
!                 do j=1,nmestot
!                    do k=1,nmestot
!                       if (k.ge.j) then
!                          VC2(j,k)=Vi(j+k*(k-1)/2)
!                       else
!                          VC2(j,k)=Vi(k+j*(j-1)/2)
!                       end if
!                    end do
!                 end do

                 !C 	           if(CurrentProcessorID.eq.0) write(*,*)'mu',(mu(j),j=1,nmestot)
		yy=Y1all(i,:)

		Ut2=0.d0
		vv=0.d0
		ZZ=0.d0
		Ut2=MATMUL(Ut1,transpose(Ut1))
		ZZ=MATMUL(Ut2,transpose(Z1))
	        do q=1,nq
		  VV=MATMUL(VC,(yy-mu))
		  ui_hat(:,g)=MATMUL(ZZ,VV)
	        end do


                 Y3=Y1-mu
                 !C 	           if(CurrentProcessorID.eq.0) write(*,*)'Y3',(Y3(j),j=1,nmestot)
                 Y4=MATMUL(VC,Y3)
                 Y5=DOT_PRODUCT(Y3,Y4)
!                 Y4_2=MATMUL(VC2,Y3)
!                 Y5_2=DOT_PRODUCT(Y3,Y4_2)

                 fi(g)=fi(g)- nmestot*log(dble(2*3.14159265))
                 fi(g)=fi(g) -det
                 fi(g)=fi(g) - Y5

                 fi(g)=fi(g)/(2.d0)
                 fi1(g)=exp(fi(g))

!                 f_ai(g)=f_ai(g)- nmestot*log(dble(2*3.14159265))
!                 f_ai(g)=f_ai(g) -det2
!                 f_ai(g)=f_ai(g) - Y5_2

!                 f_ai(g)=f_ai(g)/(2.d0)
              ELSE
                 gcourant=g
                 fi1(g)=vraisobs()
                 fi(g)=dlog(fi1(g))
		print*,'P(ci=g|u_i,Xi) not calculated'
              end if

              if (evt.ne.0) then

                 ! surv_glob=0.d0
!                  survint_glob=0.d0
!                  nxevtcurr=0
!                  fevt=0.d0
!                  do kk=1,nevt
!                     varexpsurv=0.d0
!                     if (nxevt.ne.0) then 
!                        varexpsurv=DOT_PRODUCT(Xevt((nxevtcurr+1):(nxevtcurr+nxevtspec(kk)))&
!                             ,bevt((nxevtcurr+1):(nxevtcurr+nxevtspec(kk))))
!                     end if
!                     Surv_glob=surv_glob+surv(g,kk)*exp(varexpsurv)
!                     survint_glob=survint_glob+survint(g,kk)*exp(varexpsurv)

!                     if (Devt(i).eq.kk) then
!                        fevt=log(risq(g,kk))+varexpsurv              
!                        if (ind_survint(i).eq.1) then
!                           fevt=fevt+ bevtint(1)
!                        end if
!                     end IF
!                     nxevtcurr=nxevtcurr+nxevtspec(kk)
!                  end do
!                  fi(g)=fi(g)+fevt-   &   
                 !                      (survint_glob+exp(bevtint(1))*(surv_glob-survint_glob))
		 fi_scoreTT(g)=exp(fi(g))*scoreTT(g)
		 fi_scoreYT(g)=exp(fi(g))*scoreYT(g)
                 fi(g)=fi(g)+vectvalfin(g)
              end if
              fi(g)=exp(fi(g))
           end do	
	   f1=DOT_PRODUCT(pi,fi1)
           f=DOT_PRODUCT(pi,fi)
	   SS=1.d0
	   if(idtrunc.eq.1) SS=DOT_PRODUCT(pi,STi0TT)

           if (f.eq.0.or.(f+1.d0).eq.f) then
              if(CurrentProcessorID.eq.0) write(*,*)'i',i,'f',f,'fi(g)',(fi(g),g=1,ng)
              if(CurrentProcessorID.eq.0) write(*,*)'pi(g)',(pi(g),g=1,ng)
			  print*,'fi pas bons a i=', i
			  print*, 'bvexp=',bvexp
			  print*, 'fi=', fi
              stop
           end if
	   scoregTT=0.d0
	   scoregYT=0.d0  
	   f_scoreYT=0.d0

     	  do g=1,ng
              PPI(i,g)=pi(g)*fi(g)/f
              PPItest(i,g)=pi(g)*fi1(g)/f1
	      !P_AI(i,g)=pi(g)*f_ai(g)/f_a
	      scoregTT=scoregTT+pi(g)/2.0d0*(fi_scoreTT(g)/f-SUTi0TT(g)/SS)
	      do jj=1,nea
		f_scoreYT(jj)=f_scoreYT(jj)+pi(g)*fi_scoreYT(g)*ui_hat(jj,g)
	      end do
	      retard=retard+vectretard(g)*pi(g)
           end do

!		   scorefinalTT_moy=scorefinalTT_moy+scoregTT


           ScorefinalYT=scorefinalYT+f_scoreYT(:)

	 do ii=1,nea
          do jj=1,nea
	   var_scoreYT(ii,jj)=f_scoreYT(ii)*f_scoreYT(jj)
	  end do
	end do


           ScorefinalTT=scorefinalTT+scoregTT
           var_scoreTT=var_scoreTT+scoregTT**2

           write(2,2222)numero(i),(pi(g),g=1,ng)
!	do g=1,ng
!	  P_I(i,g)=pi(g)
!	end do
        end do

	 do ii=1,nea
          do jj=1,nea
	   var_scoreYT(ii,jj)=var_scoreYT(ii,jj)-ScorefinalYT(ii)*ScorefinalYT(jj)/dble(ns)
	  end do
	end do

!	scorefinalYT_moy=scorefinalYT_moy/dble(ns)
	!var_scoreYT=var_scoreYT-scorefinalYT**2/dble(ns)

!	scorefinalTT_moy=scorefinalTT_moy/dble(ns)
	var_scoreTT=var_scoreTT-scorefinalTT**2/dble(ns)

        close(2)
	probapost=0
	do i=1,ns
	  temp=ppitest(i,1)
	  j=1
	  do g=2,ng
	      if(ppitest(i,g).gt.temp)then
			temp=ppitest(i,g)
			j=g
	      end if
	  end do
	  probapost(j)=probapost(j)+1
	end do



	open(99,file='ppitest.txt')
	open(98,file='ppi.txt')
	open(97,file='pipre.txt')
	do i=1,ns
		write(99,*)ppitest(i,:)
		write(98,*) ppi(i,:)
		write(97,*)pi
	end do





      end subroutine postprobmult_beta

! }}}




!C ******************** BGOS ********************************
! {{{

      SUBROUTINE BGOS(SX,ID,X1,X2,RO)



      implicit none
!C     ID=1:U(0,SX); ID DIFF DE 1 :N(0,SX)

      double precision ::RO,SX
      integer ::ID
      double precision ::F,V1,V2,S,DLS,RO2
      double precision ::X1,X2,UNIRAN
!C     if(CurrentProcessorID.eq.0) write(*,*)'dans bgos'


 5    CONTINUE

!C     if(CurrentProcessorID.eq.0) write(*,*)'avant rand :'

!C     X1=RAND()
!C     X2=RAND()

      X1=UNIRAN()
      X2=UNIRAN()

      IF(ID.NE.1) GO TO 10
      F=2.*SQRT(3.)
      X1=(X1-0.5)*F
      X2=(X2-0.5)*F
      GO TO 20
10    CONTINUE
      V1=2.*X1-1
      V2=2.*X2-1
      S=V1*V1+V2*V2
      IF(S.GE.1.) GO TO 5
      DLS=SQRT(-2.*LOG(S)/S)
      X1=V1*DLS
      X2=V2*DLS
20    CONTINUE
      RO2=RO*RO
      IF(ABS(RO).GT.1.E-10) X2=(X1+X2*SQRT(1./RO2-1.))*RO
      X1=X1*SX
      X2=X2*SX

!C      if(CurrentProcessorID.eq.0) write(*,*) 'X1 ',X1,' X2 ',X2
!C OK, X1 et X2 sont crees

!C      if(CurrentProcessorID.eq.0) write(*,*)'fin bgos'

      RETURN
      END subroutine bgos


!C ------------------- FIN SUBROUTINE BGOS -----------------


!C ------------------------------------------------------

      DOUBLE PRECISION FUNCTION UNIRAN()
!C
!C     Random number generator(RCARRY), adapted from F. James
!C     "A Review of Random Number Generators"
!C      Comp. Phys. Comm. 60(1990), pp. 329-344.
!C

      implicit none
      DOUBLE PRECISION SEEDS(24), TWOM24, CARRY, ONE
      PARAMETER ( ONE = 1, TWOM24 = ONE/16777216 )
      INTEGER I, J
      SAVE I, J, CARRY, SEEDS
      DATA I, J, CARRY / 24, 10, 0.0 /
      DATA SEEDS /     & 
           0.8804418, 0.2694365, 0.0367681, 0.4068699, 0.4554052, 0.2880635,     & 
           0.1463408, 0.2390333, 0.6407298, 0.1755283, 0.7132940, 0.4913043,     & 
           0.2979918, 0.1396858, 0.3589528, 0.5254809, 0.9857749, 0.4612127,     & 
           0.2196441, 0.7848351, 0.4096100, 0.9807353, 0.2689915, 0.5140357/
      UNIRAN = SEEDS(I) - SEEDS(J) - CARRY
      IF ( UNIRAN .LT. 0 ) THEN
         UNIRAN = UNIRAN + 1
         CARRY = TWOM24
      ELSE
         CARRY = 0
      ENDIF
      SEEDS(I) = UNIRAN
      I = 24 - MOD( 25-I, 24 )
      J = 24 - MOD( 25-J, 24 )
      END


! }}}



!C------------------------------------------------------------
!C     
!                CALCUL FONCTION REPARTITION NORMALE
!
!C------------------------------------------------------------

! {{{

!C==========================================================================
!C
!C
!C
!C     Fonctions pour calculer une fonction de repartition issue d'une
!C     loi normale
!C                                   ajout le 28 fevrier 2008
!C
!C===========================================================================




      double precision function alnorm ( x, upper )

!*****************************************************************************80
!
!! ALNORM computes the cumulative density of the standard normal distribution.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    David Hill
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    David Hill,
!    Algorithm AS 66:
!    The Normal Integral,
!    Applied Statistics,
!    Volume 22, Number 3, 1973, pages 424-427.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
!    over which the integration takes place.
!
!    Input, logical UPPER, determines whether the upper or lower
!    interval is to be integrated:
!    .TRUE.  => integrate from X to + Infinity;
!    .FALSE. => integrate from - Infinity to X.
!
!    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
!    distribution over the desired interval.
!
      implicit none

      real ( kind = 8 ), parameter :: a1 = 5.75885480458D+00
      real ( kind = 8 ), parameter :: a2 = 2.62433121679D+00
      real ( kind = 8 ), parameter :: a3 = 5.92885724438D+00
      real ( kind = 8 ), parameter :: b1 = -29.8213557807D+00
      real ( kind = 8 ), parameter :: b2 = 48.6959930692D+00
      real ( kind = 8 ), parameter :: c1 = -0.000000038052D+00
      real ( kind = 8 ), parameter :: c2 = 0.000398064794D+00
      real ( kind = 8 ), parameter :: c3 = -0.151679116635D+00
      real ( kind = 8 ), parameter :: c4 = 4.8385912808D+00
      real ( kind = 8 ), parameter :: c5 = 0.742380924027D+00
      real ( kind = 8 ), parameter :: c6 = 3.99019417011D+00
      real ( kind = 8 ), parameter :: con = 1.28D+00
      real ( kind = 8 ), parameter :: d1 = 1.00000615302D+00
      real ( kind = 8 ), parameter :: d2 = 1.98615381364D+00
      real ( kind = 8 ), parameter :: d3 = 5.29330324926D+00
      real ( kind = 8 ), parameter :: d4 = -15.1508972451D+00
      real ( kind = 8 ), parameter :: d5 = 30.789933034D+00
      real ( kind = 8 ), parameter :: ltone = 7.0D+00
      real ( kind = 8 ), parameter :: p = 0.398942280444D+00
      real ( kind = 8 ), parameter :: q = 0.39990348504D+00
      real ( kind = 8 ), parameter :: r = 0.398942280385D+00
      logical up
      logical upper
      real ( kind = 8 ), parameter :: utzero = 18.66D+00
      real ( kind = 8 ) x
      real ( kind = 8 ) y
      real ( kind = 8 ) z

      up = upper
      z = x

      if ( z < 0.0D+00 ) then
         up = .not. up
         z = - z
      end if

      if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

         if ( up ) then
            alnorm = 0.0D+00
         else
            alnorm = 1.0D+00
         end if

         return

      end if

      y = 0.5D+00 * z * z

      if ( z <= con ) then

         alnorm = 0.5D+00 - z * ( p - q * y     & 
   / ( y + a1 + b1/ ( y + a2 + b2 / ( y + a3 ))))

      else

         alnorm = r * exp ( - y )     &    
      / ( z + c1 + d1 / ( z + c2 + d2 / ( z + c3 + d3     &  
      / ( z + c4 + d4 / ( z + c5 + d5 / ( z + c6 ))))))

      end if

      if ( .not. up ) then
         alnorm = 1.0D+00 - alnorm
      end if

      return
      end function alnorm



      subroutine normp ( z, p, q, pdf )

!*****************************************************************************80
!
!! NORMP computes the cumulative density of the standard normal distribution.
!
!  Discussion:
!
!    This is algorithm 5666 from Hart, et al.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Alan Miller
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thacher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Z, divides the real ( kind = 8 ) line into two
!    semi-infinite intervals, over each of which the standard normal
!    distribution is to be integrated.
!
!    Output, real ( kind = 8 ) P, Q, the integrals of the standard normal
!    distribution over the intervals ( - Infinity, Z] and
!    [Z, + Infinity ), respectively.
!
!    Output, real ( kind = 8 ) PDF, the value of the standard normal distribution
!    at Z.
!
      implicit none

      real ( kind = 8 ) :: cutoff = 7.071D+00
      real ( kind = 8 ) expntl
      real ( kind = 8 ) p
      real ( kind = 8 ) :: p0 = 220.2068679123761D+00
      real ( kind = 8 ) :: p1 = 221.2135961699311D+00
      real ( kind = 8 ) :: p2 = 112.0792914978709D+00
      real ( kind = 8 ) :: p3 = 33.91286607838300D+00
      real ( kind = 8 ) :: p4 = 6.373962203531650D+00
      real ( kind = 8 ) :: p5 = 0.7003830644436881D+00
      real ( kind = 8 ) :: p6 = 0.03526249659989109D+00
      real ( kind = 8 ) pdf
      real ( kind = 8 ) q
      real ( kind = 8 ) :: q0 = 440.4137358247522D+00
      real ( kind = 8 ) :: q1 = 793.8265125199484D+00
      real ( kind = 8 ) :: q2 = 637.3336333788311D+00
      real ( kind = 8 ) :: q3 = 296.5642487796737D+00
      real ( kind = 8 ) :: q4 = 86.78073220294608D+00
      real ( kind = 8 ) :: q5 = 16.06417757920695D+00
      real ( kind = 8 ) :: q6 = 1.755667163182642D+00
      real ( kind = 8 ) :: q7 = 0.08838834764831844D+00
      real ( kind = 8 ) :: root2pi = 2.506628274631001D+00
      real ( kind = 8 ) z
      real ( kind = 8 ) zabs

      zabs = abs ( z )
!
!  37 < |Z|.
!
      if ( 37.0D+00 < zabs ) then

         pdf = 0.0D+00
         p = 0.0D+00
!
!     |Z| <= 37.
!
      else

         expntl = exp ( - 0.5D+00 * zabs * zabs )
         pdf = expntl / root2pi
!
!  |Z| < CUTOFF = 10 / sqrt(2).
!
         if (zabs.lt.cutoff) then

          p = expntl * ((((((  p6   * zabs + p5 ) * zabs + p4 ) * zabs     &   
    + p3 ) * zabs + p2 ) * zabs + p1 ) * zabs + p0 ) / (((((((     & 
      q7   * zabs + q6 ) * zabs+ q5 ) * zabs + q4 ) * zabs + q3 )     & 
      * zabs + q2 ) * zabs + q1 ) * zabs + q0 )
!
!  CUTOFF <= |Z|.
!
         else

            p = pdf / ( zabs + 1.0D+00 / ( zabs + 2.0D+00 / ( zabs + 3.0D+00 / (     &  
     zabs + 4.0D+00 / ( zabs + 0.65D+00 )))))

         end if

      end if

      if ( z < 0.0D+00 ) then
         q = 1.0D+00 - p
      else
         q = p
         p = 1.0D+00 - q
      end if

      return
      end subroutine normp
      subroutine nprob ( z, p, q, pdf )

!*****************************************************************************80
!
!! NPROB computes the cumulative density of the standard normal distribution.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    AG Adams
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39:
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, Number 2, May 1969, pages 197-198.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Z, divides the real ( kind = 8 ) line into
!    two semi-infinite intervals, over each of which the standard normal
!    distribution is to be integrated.
!
!    Output, real ( kind = 8 ) P, Q, the integrals of the standard normal
!    distribution over the intervals ( - Infinity, Z] and
!    [Z, + Infinity ), respectively.
!
!    Output, real ( kind = 8 ) PDF, the value of the standard normal
!    distribution at Z.
!
      implicit none

      real ( kind = 8 ), parameter :: a0 = 0.5D+00
      real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
      real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
      real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
      real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
      real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
      real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
      real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
      real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
      real ( kind = 8 ), parameter :: b1 = 0.000000038052D+00
      real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
      real ( kind = 8 ), parameter :: b3 = 0.000398064794D+00
      real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
      real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
      real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
      real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
      real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
      real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
      real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
      real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
      real ( kind = 8 ) p
      real ( kind = 8 ) pdf
      real ( kind = 8 ) q
      real ( kind = 8 ) y
      real ( kind = 8 ) z
      real ( kind = 8 ) zabs

      zabs = abs ( z )
!
!  |Z| between 0 and 1.28
!
      if ( abs ( z ) <= 1.28D+00 ) then

         y = a0 * z * z
         pdf = exp ( - y ) * b0

         q = a0 - zabs * ( a1 - a2 * y   / ( y + a3 - a4 / ( y + a5 + a6/ ( y + a7 ))))
!
!     |Z| between 1.28 and 12.7
!
      else if ( abs ( z ) <= 12.7D+00 ) then

         y = a0 * z * z
         pdf = exp ( - y ) * b0

         q = pdf / ( zabs - b1 + b2 / ( zabs + b3 + b4   / ( zabs - b5 + b6     &   
 / ( zabs + b7 - b8    / ( zabs + b9 + b10    / ( zabs + b11 ))))))
!
!     Z far out in tail.
!
      else

         q = 0.0D+00
         pdf = 0.0D+00

      end if

      if ( z < 0.0D+00 ) then
         p = q
         q = 1.0D+00 - p
      else
         p = 1.0D+00 - q
      end if

      return
      end subroutine nprob

! }}}





!C *************************************************************
!C
!C     evolution a posteriori (sim)
!C
!C *************************************************************
 ! {{{

      subroutine evol_pred_marg_sim(b,npm,nsim,ymarg,ycond,mui1,muicond1)
 !ymarg : y predits
 !mu : latent
      use commun
      use optim_parallele



      implicit none
      integer::npm,ier,nnn,ntrtot
      integer::nmestot,tronque
      integer :: i,j,k,l,jj,kk,q,m,g,nmoins,niter
      integer ::jj1,k1,q1,h,nsim,ll,ifault,ntronque,jjj

      double precision :: aa,bb,eps,X2,xinbta,SX,temp     &
           , aa1,bb1,beta,beta_ln,cc1,dd1,ytemp,diff,det,uniran,ytemp2



!C --------------   dimensions revues
      double precision,dimension(npm)::b
      double precision,dimension(npmtot)::b1
      double precision,dimension(nvc)::bh1
      double precision,dimension(ng)::pi
      double precision,dimension(nea*nq)::ui_hat
      double precision,dimension(nea,maxmestot)::ZZ
      double precision,dimension(nvarprob)::bprob,Xprob
      double precision,dimension(-2:ntrspl-3)::splaa
      double precision,dimension(nea,nea)::Ut,Ut1,Ut2

      double precision,dimension(nvarglobssG)::b3
      double precision,dimension(2*np+1)::b0
      double precision,dimension(ncontdep+ncont)::sum
      double precision,dimension(nvdepssG)::bdepssG
      double precision,dimension(nvdepG)::bdepG
      double precision,dimension(nvarglobG)::b2
      double precision,dimension((ncontdep+ncont)*nq)::b4

      double precision, dimension(maxmestot)::mu,usim,ysim,VV,ysimcond,esim,yy
      double precision,dimension(maxmestot,2*np+1)::Z
      double precision,dimension(maxmestot,nea)::Z1
      double precision,dimension(maxmestot,nvarglobssG)::XglobssG
      double precision,dimension(maxmestot,nvdepssG)::XdepssG
      double precision,dimension(maxmestot,(ncontdep+ncont)*nq)::Xcont
      double precision,dimension(maxmestot,nvarglobG)::XglobG
      double precision,dimension(maxmestot,nvdepG)::XdepG
      double precision,dimension(maxmestot,nea)::P
      double precision,dimension(maxmestot,maxmestot)::VC,VCAR,Vw,sigmaq
      double precision,dimension(maxmestot*(maxmestot+1)/2)::Vi
      double precision,dimension(-1:nztr+2)::zi_eval

      double precision,dimension(nsim,maxmestot)::ysim2,ysim2cond
      double precision::INV_ISPLINES


!C NON REVU

      double precision,dimension(ns*nq*mesure,ng)::ymarg,mui1,ycond,muicond1

      eps=1.D-20
      k=0
      l=0
      DO j=1,npmtot
         if (ide(j).eq.1) then
            k=k+1
            b1(j)=b(k)
         else
            l=l+1
            b1(j)=bne(l)
         end if
      END DO

!C	  if(CurrentProcessorID.eq.0) write(*,*)'dans evol_pred 1'
!C ------ Ut : ---------------------------------------------
!C------- Creation de Ut la transformee de Cholesky ------


      Do j=1,nvc
         bh1(j)=b1(nef+j)
      End do

!C     if(CurrentProcessorID.eq.0) write(*,*)'bh1',(bh1(k),k=1,nvc)

      Ut=0.d0
      Ut1=0.d0
      If (idiag.eq.1) then
         do j=1,nea
            do k=1,nea
               if (j.eq.k) then
                  Ut(j,k)=bh1(j)
                  Ut1(j,k)=bh1(j)
               end if
            end do
         end do
      else if (idiag.eq.0) then
         do j=1,nea
            do k=j,nea
               Ut(k,j)=bh1(j+k*(k-1)/2)
               Ut1(k,j)=bh1(j+k*(k-1)/2)
            end do
         end do
      end if

!C fonction de risq a calculer pour chaque i et chaque g


!C------ boucle sur les sujets i-----------------------------
      mui1=0.d0
      ymarg=0.d0
      ycond=0.d0

      DO i=1,ns

	nmestot=0
	jj=0
        do q=1,nq
          do j=1,nmes(i,q)
             jj=jj+1
          end do
       End do
	
	nmestot=jj

         XdepssG=0.d0
         l=0
         Do k=1,nvdep
            If (idvdep(k).eq.1) then
               l=l+1
               jj=0
               do q=1,nq
                  do j=1,nmes(i,q)
                     jj=jj+1
                     XdepssG(jj,l)=Xdep(Navant_i(i)+jj,k)
                  end do
               End do
               if (jj.ne.nmestot) then
                  if(CurrentProcessorID.eq.0) write(*,*)'probleme XdepssG(jj,k)'
                  stop
               end if
            end if
         end do


!C creation du vecteur de variables explicatives dep du temps : XdepG avec mixture


         XdepG=0.d0
         l=0
         Do k=1,nvdep
            If (idvdep(k).eq.2) then
               l=l+1
               jj=0
               do q=1,nq
                  do j=1,nmes(i,q)
                     jj=jj+1
                     XdepG(jj,l)=Xdep(Navant_i(i)+jj,k)
                  end do
               End do
               if (jj.ne.nmestot) then
                  if(CurrentProcessorID.eq.0) write(*,*)'probleme XdepG(jj,k)'
                  stop
               end if
            end if
         end do


!C------ creation de XglobG(j,k) (avec mixture) ----------

         XglobG=0.d0

         l=0
         Do k=1,nv
            If (idglob(k).eq.1.and.idg(k).eq.1) then
               l=l+1
               Do j=1,nmestot
                  XglobG(j,l)=X(i,k)
               End do

            Else if (idglob(k).eq.2.and.idg(k).eq.1) then
               jj=0
               l=l+1
               Do q=1,nq
                  Do j=1,nmes(i,q)
                     jj=jj+1
                     Do m=0,np
                        if (m.ge.1) then
                           XglobG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)* (0.1**(m-1))
                        else
                           XglobG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)
                        end if
                     End do
                  End do
               End do
               l=l+np
            End if
         End do


         if (nvarglobssG+nvarglobG*ng.ne.nvarglob) then
            if(CurrentProcessorID.eq.0) write(*,*) 'pb : def nvarglobssG/G'
            stop
         end if


!C------ creation de Xcont(q,k,j) -----------------
!C creation des variables pour les contrastes

        Xcont=0.d0
         nmoins=0
         do q=1,nq
            do j=1,nmes(i,q)
               h=0
               do k=1,nvdep
                  if (idvdep(k).eq.3.or.idvdep(k).eq.4) then
                     h=h+1
                     Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)=Xdep(Navant_i(i)+nmoins+j,k)
                  end if
               end do
               h=0
               do k=1,nv
                  if (idtest(k).eq.1) then
                     h=h+1
                     Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)=X(i,k)
                  end if
                  if (idtest(k).eq.2) then
                     Do m=0,np
                        h=h+1
                        if (m.ge.1) then
                           Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)     & 
                                =X(i,k)*(t(Navant_i(i)+nmoins+j)**m)*(0.1**(m-1))

                        else
                           Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)     & 
                                =X(i,k)*(t(Navant_i(i)+nmoins+j)**m)
                        end if
                     End do
                  end if
               end do
               if (h.ne.ncont) then
                  if(CurrentProcessorID.eq.0) write(*,*)'h=',h,'problem in cont 2'
                  stop
               end if
            End do
            nmoins=nmoins+nmes(i,q)
         end do
         if (nmoins.ne.nmestot) then
            if(CurrentProcessorID.eq.0) write(*,*)'problem in cont', nmoins,nmestot
            stop
         end if



         VCAR=0.d0
         if (iAR.gt.0.and.QUECONT.eq.1) then
            jj=0
            Do q=1,nq
               jj1=0
               Do q1=1,nq
                  Do k=1,nmes(i,q)
                     do k1=1,nmes(i,q1)
                        if (iAR.eq.1) then
                           VCAR(jj+k,jj1+k1)=(b1(nef+nvc+ng-1+nq+1)**2)*min(t(Navant_i(i)+jj+k),t(Navant_i(i)+jj1+k1))
                        else if (iAR.eq.2) then
                           VCAR(jj+k,jj1+k1)=(b1(nef+nvc+ng-1+nq+1)**2)*exp(-b1(nef+nvc+ng-1+nq+2)     &  
                        *abs(t(Navant_i(i)+jj+k)-t(Navant_i(i)+jj1+k1)))
                        end if
                     end do
                  end do
                  jj1=jj1+nmes(i,q1)
               end do
               jj=jj+nmes(i,q)
            end do
         end if


         Sigmaq=0.d0
         Vw=0.d0
         jj=0.d0
         ntrtot=0
         Do q=1,nq
            Do k=1,nmes(i,q)
               do l=1,nmes(i,q)
                  Vw(jj+k,jj+l)=b1(nef+nvc+ng-1+q)**2
                  if (k.eq.l) then
                     Sigmaq(jj+k,jj+l)=b1(nef+nvc+ng-1+nq+nAR+q)**2
                  end if
               end do
            end do
            jj=jj+nmes(i,q)
         End do


!C	     if(CurrentProcessorID.eq.0) write(*,*)'i',i,'nmes',nmestot
!C------Vecteurs de parms sans mixtures -------------------
         bdepssG=0.d0
	 jj=0
         jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)
         l=0
         kk=0
         Do k=1,nvdep
            If (idvdep(k).eq.1) then
               l=l+1
               kk=kk+1
               bdepssG(l)=b1(jj+kk)
            end if
            If (idvdep(k).eq.2) then
               kk=kk+ng
            end if
         end do
         if (l.ne.nvdepssG) then
            if(CurrentProcessorID.eq.0) write(*,*) 'probleme longueur vecture bdepssG'
            stop
         end if


         b3=0.d0
         b4=0.d0
	 jj=0

         jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep
         l=0
         kk=0
         Do k=1,nv
            If (idglob(k).eq.1.and.idg(k).eq.0) then
               l=l+1
               kk=kk+1
               b3(l)=b1(jj+kk)  ! correspond a beta
            Else if (idglob(k).eq.2.and.idg(k).eq.0) then
               l=l+1
               kk=kk+1
               Do m=0,np
                  b3(l+m)=b1(jj+kk+m)
               End do
               l=l+np
               kk=kk+np
            Else if (idglob(k).eq.1.and.idg(k).eq.1) then
               kk=kk+ng
            Else if (idglob(k).eq.2.and.idg(k).eq.1) then
               kk=kk+ng*(np+1)
            End if
         End do

!C     b4 : contrastes sur les tests
         b4=0.d0
         sum=0.d0
         do j=1,ncontdep
            do q=1,nq-1
               b4((ncont+ncontdep)*(q-1)+j)=b1(nvarprob*(ng-1)+nrisq     &   
        +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+(j-1)*(nq-1)+q)
               sum(j)=sum(j)+b4((ncont+ncontdep)*(q-1)+j)
            end do
            b4((ncont+ncontdep)*(nq-1)+j)=-sum(j)
         end do
         do j=1,ncont
            do q=1,nq-1
               b4((ncont+ncontdep)*(q-1)+ncontdep+j)=b1(nvarprob*(ng-1)     &     
         +nrisq +nvarxevt +ng+ng*(2*np+1)+nprmdep+nvarglob     &  
         + ncontdep*(nq-1)+(j-1)*(nq-1)+q)
               sum(ncontdep+j)=sum(ncontdep+j)+b4((ncont+ncontdep)*(q-1)+ncontdep+j)
            end do
            b4((ncont+ncontdep)*(nq-1)+ncontdep+j)=-sum(ncontdep+j)
         end do


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (ng.eq.1) then

            VC=0.d0
            P=0.d0
            P=MATMUL(Z1,Ut)


            do k=1,2*np+1
               b0(k)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+k)
            end do
            mu=0.d0
            mu=MATMUL(Z,b0)+MATMUL(XglobssG,b3)+MATMUL(XdepssG,bdepssG)+MATMUL(Xcont,b4)

			
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			
			!C ---- Bayesian estimates of ui
	    VC=0.d0
	    P=0.d0
            VC=MATMUL(P,transpose(P))+Sigmaq+Vw+VCAR

	    Vi=0.d0
               jj=0
               do j=1,nmestot
                  do k=j,nmestot
                     jj=j+k*(k-1)/2
                     Vi(jj)=VC(j,k)
                  end do
               end do

         CALL DSINV(Vi,nmestot,eps,ier,det,0)
         if (ier.eq.-1) then
          !  if(CurrentProcessorID.eq.0) write(*,*)'dsinv vi ier',ier,' i=',i,' nmes(i)=',nmestot
          !  if(CurrentProcessorID.eq.0) write(*,*)'b1', (b1(k),k=1,npm)
            nnn=nmestot*(nmestot+1)/2
          !  if(CurrentProcessorID.eq.0) write(*,*)'Vi', (Vi(k),k=1,nnn)
            stop
         end if

!c     retransformation du vecteur Vi en matrice :
         VC=0.d0
         do j=1,nmestot
            do k=1,nmestot
               if (k.ge.j) then
                  VC(j,k)=Vi(j+k*(k-1)/2)
               else
                  VC(j,k)=Vi(k+j*(j-1)/2)
               end if
            end do
         end do
		
		yy=Y1all(i,:)
		ui_hat=0.d0
		Ut2=0.d0
		vv=0.d0
		ZZ=0.d0
		Ut2=MATMUL(Ut,transpose(Ut))
		ZZ=MATMUL(Ut2,transpose(Z1))
	       do q=1,nq
		  VV=MATMUL(VC,(yy-mu))
		  ui_hat(:)=MATMUL(ZZ,VV)
	       end do

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            VC=0.d0
            VC=MATMUL(P,transpose(P))+Sigmaq+Vw+VCAR
!c     Vi en vecteur
           jj=0
            do j=1,nmestot
               do k=j,nmestot
                  jj=j+k*(k-1)/2
                  Vi(jj)=VC(j,k)
               end do
            end do
            ier=0
            CALL DMFSD(Vi,nmestot,EPS,IER,0)
            if (ier.eq.-1) then
               if(CurrentProcessorID.eq.0) write(*,*)'dmfsd vi ier',ier,' i=',i,' nmes(i)=',nmestot
               if(CurrentProcessorID.eq.0) write(*,*)'b1', (b1(k),k=1,npmtot)
               nnn=nmestot*(nmestot+1)/2
               if(CurrentProcessorID.eq.0) write(*,*)'Vi', (Vi(k),k=1,nnn)
               stop
            end if


!c     retransformation du vecteur Vi en matrice :
            VC=0.d0
            do j=1,nmestot
               do k=1,j
                  VC(j,k)=Vi(k+j*(j-1)/2)
               end do
            end do

!C     simulation des vecteurs normaux
               ysim2=0.d0
               ysim2cond=0.d0
			   
            do l=1,nsim
               ntronque=0
               tronque=1
               do while (tronque.eq.1)
                  tronque=0
             	  ysim=0.d0
           	    usim=0.d0
	   	    esim=0.d0
           	    do m=1,nmestot
           	       SX=1.d0
          	        call bgos(SX,0,usim(m),X2,0.d0)
			  call bgos(SX,0,esim(m),X2,0.d0)
          	     end do


!C     if (i.eq.8.and.l.eq.2) then
!C     if(CurrentProcessorID.eq.0) write(*,*)'nmestot',nmestot,'usim',(usim(m),m=1,nmestot)
!C     end if

               ysim=mu+MATMUL(VC,usim)
	       ysimcond=mu+MATMUL(Z1,ui_hat)+MATMUL(sigmaq,esim)

               ll=0

               if (l.eq.1) then

                  jj=0
                  do q=1,nq
                     do j=1,nmes(i,q)
                        jj=jj+1
                        mui1(Navant_i(i)+jj,1)=ysim(jj)
			muicond1(Navant_i(i)+jj,1)=ysimcond(jj)
                     end do
                  end do
               end if




               ntrtot=0
		ll=0
               Do q=1,nq
                  if (idord(q).le.0) then
                     aa1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &   
                 +ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot))/(1+exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &  
                  +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)     &  
                  +1+ntrtot)))
                     bb1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &  
                  +ncont*(nq-1)+2+ntrtot))/( 1.D0+ exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &      
              +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+2+ntrtot)))


                     bb1=aa1*(1.d0-aa1)*bb1

                     cc1=b1(nvarprob*(ng-1)+nrisq+nvarxevt +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &  
                  +ncont*(nq-1)+3+ntrtot)
                     dd1=abs(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &  
                  +ncont*(nq-1)+4+ntrtot))


                     aa=aa1*aa1*(1-aa1)/bb1-aa1
                     bb=aa*(1-aa1)/aa1
                     beta=beta_ln(aa,bb)

!!C					 if(CurrentProcessorID.eq.0) write(*,*)'aa-dd',aa,bb,cc1,dd1
!C
!C                     if(CurrentProcessorID.eq.0) write(*,*)'ntrtot',ntrtot,'ntr',ntr

                     do j=1,nmes(i,q)
                        ytemp=ysim(ll+j)*dd1+cc1
			ytemp2=ysimcond(ll+j)*dd1+cc1
                        if (ytemp.lt.0.or.ytemp.gt.1) then
!C                           tronque=1
!C                           ntronque=ntronque+1
!C                           if(CurrentProcessorID.eq.0) write(*,*)'ytemp en dehors de 0,1',ytemp
!C     &                          ,ysim(ll+j),dd1,cc1
                           if (ytemp.lt.0) then
                              ytemp=0.d0
                           end if
                           if (ytemp.gt.1) then
                              ytemp=1.d0
                           end if
!C                           stop
                        end if
						
						
			if (ytemp2.lt.0.or.ytemp2.gt.1) then
!C                              if(CurrentProcessorID.eq.0) write(*,*)'ytemp en dehors de 0,1',ytemp
!C     &                             ,ysim(ll+j),dd1,cc1
!C                              tronque=1
!C                              ntronque=ntronque+1
                              if (ytemp2.lt.0) then
                                 ytemp2=0.d0
                              end if
                              if (ytemp2.gt.1) then
                                 ytemp2=1.d0
                              end if
                           end if		  
							  
                        ysim2(l,ll+j)=xinbta(aa,bb,beta,ytemp,ifault)
			ysim2cond(l,ll+j)=xinbta(aa,bb,beta,ytemp2,ifault)
                     end do
                     ll=ll+nmes(i,q)
                     ntrtot=ntrtot+ntr
                  else if (idord(q).eq.1) then
                     zi_eval=0.d0
                     do j=-1,nztr+2
                        zi_eval(j)=zitr(j,q)
                     end do

                     bb=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &      
              +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                     splaa=0.d0
                     do kk=2,ntrspl
                        splaa(kk-3)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &   
                    +ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot)*b1(nvarprob*(ng-1)+nrisq     &        
               +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot)
                     end do

!C                     if(CurrentProcessorID.eq.0) write(*,*)
!c                     if(CurrentProcessorID.eq.0) write(*,*)'ijq',i,q,l
!c                     if(CurrentProcessorID.eq.0) write(*,*)


                     ntrtot=ntrtot+ntrspl
                     do j=1,nmes(i,q)
                        niter=0
                        diff=0.d0
                        ier=0
                        ysim2(l,ll+j)=INV_ISPLINES(ysim(ll+j),splaa,bb,zi_eval,ier,niter,diff)
                        if ((ier.eq.3).or.(ier.ne.1.and.diff.gt.1.d-3))then
                           if(CurrentProcessorID.eq.0) write(*,*)'problem INV_ISPLINE',ysim2(l,ll+j),ysim(ll+j),i,j,q
                           if(CurrentProcessorID.eq.0) write(*,*)'zi_eval apres',nztr,(zi_eval(k),k=-1,nztr+2)
                           if(CurrentProcessorID.eq.0) write(*,*)'i',i,j,q,l,ll,ysim(ll+j)
                           if(CurrentProcessorID.eq.0) write(*,*)splaa,bb,zi_eval,ier,niter,diff
                           stop
                        end if
			ysim2cond(l,ll+j)=INV_ISPLINES(ysimcond(ll+j),splaa,bb,zi_eval,ier,niter,diff)
                           if ((ier.ne.1.and.diff.gt.1.d-3).or.ier.eq.3) then
                              if(CurrentProcessorID.eq.0) write(*,*)'problem INV_ISPLINE',  ysim2cond(l,ll+j),ysimcond(ll+j),i,j,q,l
                              if(CurrentProcessorID.eq.0) write(*,*)'zi_eval apres',nztr,(zi_eval(k) ,k=-1,nztr+2)
                              if(CurrentProcessorID.eq.0) write(*,*)'i',i,j,q,l,ll,ysimcond(ll+j)
                              if(CurrentProcessorID.eq.0) write(*,*)'bb',bb
                              if(CurrentProcessorID.eq.0) write(*,*)'spl',splaa
                              stop
                           end if
                     end do
                     ll=ll+nmes(i,q)

                  else if (idord(q).eq.2) then
                     do j=1,nmes(i,q)
                        aa=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &   
                    +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                        if (ysim(ll+j).le.aa) then
                           ysim2(l,ll+j)=range(1+2*(q-1))
                        else
                           ysim2(l,ll+j)=range(2+2*(q-1))
                        end if
						 if (ysimcond(ll+j).le.aa) then
                           ysim2cond(l,ll+j)=range(1+2*(q-1))
                        else
                           ysim2cond(l,ll+j)=range(2+2*(q-1))
                        end if
                     end do
                     ll=ll+nmes(i,q)
                     ntrtot=ntrtot+idord(q)-1
                  else if  (idord(q).gt.2) then
                     do j=1,nmes(i,q)
                        ysim2(l,ll+j)=range(1+2*(q-1))
						ysim2cond(l,ll+j)=range(1+2*(q-1))
                        aa=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &     
                 +ncontdep*(nq-1) +ncont*(nq-1)+ntrtot+1)
                        do k=1,idord(q)-2
                           bb=aa+b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &   
                    +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1)*b1(nvarprob*(ng-1)+nrisq     &       
                +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1)
                           if (ysim(ll+j).gt.aa.and.ysim(ll+j).le.bb)then
                              ysim2(l,ll+j)=range(1+2*(q-1))+k
                           end if
						   if (ysimcond(ll+j).gt.aa.and.ysimcond(ll+j).le.bb)then
                              ysim2cond(l,ll+j)=range(1+2*(q-1))+k
                           end if
                           aa=bb
                        end do
                        if (ysim(ll+j).gt.bb) then
                           ysim2(l,ll+j)=range(2+2*(q-1))
                        end if
						 if (ysimcond(ll+j).gt.bb) then
                           ysim2cond(l,ll+j)=range(2+2*(q-1))
                        end if
                     end do
                     ll=ll+nmes(i,q)
                     ntrtot=ntrtot+idord(q)-1
                  end if
               end do
            end do
!C            if (ntronque.ne.0) then
!C               if(CurrentProcessorID.eq.0) write(*,*)'ntronque',i,ntronque
!C            end if
         end do
            ll=0
            do q=1,nq
               do j=1,nmes(i,q)
                  do l=1,nsim
                     ymarg(Navant_i(i)+ll+j,1)= ymarg(Navant_i(i)+ll+j,1)+ysim2(l,ll+j)
                     ycond(Navant_i(i)+ll+j,1)= ycond(Navant_i(i)+ll+j,1)+ysim2cond(l,ll+j)
                  end do
                  ymarg(Navant_i(i)+ll+j,1)=dble(ymarg(Navant_i(i)+ll+j,1)/dble(nsim))
                  ycond(Navant_i(i)+ll+j,1)=dble(ycond(Navant_i(i)+ll+j,1)/dble(nsim))
                  if (idord(q).le.0) then
                     ymarg(Navant_i(i)+ll+j,1)=range(1+2*(q-1))+ymarg(Navant_i(i)+ll+j,1)*(range(2+2*(q-1))-range(1+2*(q-1)))
	         	     ycond(Navant_i(i)+ll+j,1)=range(1+2*(q-1))+ycond(Navant_i(i)+ll+j,1)*(range(2+2*(q-1))-range(1+2*(q-1)))
                     if (ymarg(Navant_i(i)+ll+j,1).lt.range(1+2*(q-1))) then
                        ymarg(Navant_i(i)+ll+j,1)=range(1+2*(q-1))
                     end if
                     if (ycond(Navant_i(i)+ll+j,1).lt.range(1+2*(q-1))) then
                        ycond(Navant_i(i)+ll+j,1)=range(1+2*(q-1))
                     end if
                     if (ymarg(Navant_i(i)+ll+j,1).gt.range(2+2*(q-1))) then
                        ymarg(Navant_i(i)+ll+j,1)=range(2+2*(q-1))
                     end if
                     if (ycond(Navant_i(i)+ll+j,1).gt.range(2+2*(q-1))) then
                        ycond(Navant_i(i)+ll+j,1)=range(2+2*(q-1))
                     end if
                  end if
               end do
               ll=ll+nmes(i,q)
            end do



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   else !ng>1
!C------creation de la matrice des puissances du temps ----

            Xprob=0.d0
            Xprob(1)=1
            l=0
            do k=1,nv
               if (idprob(k).eq.1) then
                  l=l+1
                  Xprob(1+l)=X(i,k)
               end if
            end do
            if (l+1.ne.nvarprob) then
               if(CurrentProcessorID.eq.0) write(*,*)'problem nvarprob funcpa'
               stop
            end if

            pi=0.d0
            temp=0.d0
            Do g=1,ng-1
               do k=1,nvarprob
                  bprob(k)=b1((k-1)*(ng-1)+g)
               end do
               temp=temp+exp(DOT_PRODUCT(bprob,Xprob))
               pi(g)=exp(DOT_PRODUCT(bprob,Xprob))
               if ((pi(g)+1).eq.pi(g)) then
                  if(CurrentProcessorID.eq.0) write(*,*)'g',g,'pi(g)',pi(g)
                  if(CurrentProcessorID.eq.0) write(*,*)'bprob',(bprob(j),j=1,nvarprob)
                  if(CurrentProcessorID.eq.0) write(*,*)'Xprob',(Xprob(j),j=1,nvarprob)
                  stop
               end if
            end do

            pi(ng)=1/(1+temp)
			
            do g=1,ng

            Z=0.d0
              jj=0

              Do q=1,nq
                 Do j=1,nmes(i,q)
                    jj=jj+1
		    jjj=1
                    l=0

		    Z(jj,jjj)=1.d0
		    jjj=jjj+1

                    Do k=2,np+1
			    Z(jj,jjj)=(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)**(k-1)
			    temp=dble(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)
			    !temp=(temp**2+trn)**0.5d0/temp
			    !Z(jj,jjj+1)=(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)**(k-1)*temp
			    temp=(temp**2+trn)**0.5d0
			    Z(jj,jjj+1)=temp!

		            If (k.gt.2) then
			!	Z(jj,jjj)=Z(jj,jjj)*(0.1)**(k-2)
		!		Z(jj,jjj+1)=Z(jj,jjj+1)*(0.1)**(k-2)*temp
		                !Z(jj,k)=Z(jj,k)*(0.1)**(k-2)!! a changer !!
				print*,'pas code pour np=2'
				stop
		            end if

			    jjj=jjj+2
                    End do
                 End do
              End do
      	      nmestot=jj

              Z1=0.d0
	      jj=0
              Do q=1,nq
                 Do j=1,nmes(i,q)
                    jj=jj+1
                    l=0
                    Do k=1,np+1
		       If(idea(k).eq.1) then
		          l=l+1
		          if(k.eq.1) Z1(jj,l)=Z(jj,k) 
			  if(k.eq.2)then
				Z1(jj,l)=Z(jj,k)/2.0d0*(1.0d0-(Z(jj,k)**2+trn)**0.5d0/Z(jj,k))
			  	Z1(jj,l+1)=Z(jj,k)/2.0d0*(1.0d0+(Z(jj,k)**2+trn)**0.5d0/Z(jj,k))

			  end if
		       end if
                    End do
                 End do
              End do

!C------ creation de XglobssG(j,k) (sans mixture) ----------

              XglobssG=0.d0

              l=0
              Do k=1,nv
                 If (idglob(k).eq.1.and.idg(k).eq.0) then
                    l=l+1
                    Do j=1,nmestot
                       XglobssG(j,l)=X(i,k)
                    End do
                 Else if (idglob(k).eq.2.and.idg(k).eq.0) then
                    jj=0
                    l=l+1
                    Do q=1,nq
                       Do j=1,nmes(i,q)
                          jj=jj+1

                          Do m=0,np
   			     temp=(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0
                             !C Attention : changer les Xglob avec interaction
                             if (m.ge.1) then
                                XglobssG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)-temp)!((t(Navant_i(i)+jj)-temp)**m)*(0.1**(m-1))
			        temp=dble(t(Navant_i(i)+jj)-(b1((nvarprob)*(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)
			        temp=(temp**2+trn)**0.5d0!/temp
			        XglobssG(jj,l+m+1)=X(i,k)*temp!((t(Navant_i(i)+jj)-(b1((nvarprob)*&
						!(ng-1)+nrisq+nvarxevt+g)-65.0d0)/10.d0)**m)*(0.1**(m-1))*temp
                             else
                                XglobssG(jj,l+m)=X(i,k)
                             end if

                          End do

                       End do
                    End do
                    l=l+2*np
                 End if
              End do
               if (g.ne.ng) then
                  pi(g)=pi(g)*pi(ng)
               end if

               if ((pi(g)+1).eq.pi(g)) then
                  if(CurrentProcessorID.eq.0) write(*,*)'g',g,'pi(g)',pi(g)
                  if(CurrentProcessorID.eq.0) write(*,*)'bprob',(bprob(j),j=1,nvarprob)
                  if(CurrentProcessorID.eq.0) write(*,*)'Xprob',(Xprob(j),j=1,nvarprob)

                  stop
               end if

!C     if(CurrentProcessorID.eq.0) write(*,*)'i',i,'det',det
!C----- parametres avec mixtures ---------------------------

               b0=0.d0
               Do j=1,(2*np+1)
                  b0(j)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+(j-1)*ng+g)
               End do

!C     creation parametres mixture vdep du temps
               bdepG=0.d0
               jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)
               l=0
               kk=0
               Do k=1,nvdep
                  If (idvdep(k).eq.2) then
                     l=l+1
                     kk=kk+g
                     bdepG(l)=b1(jj+kk)
                  end if
                  If (idvdep(k).eq.1) then
                     kk=kk+1
                  end if
               end do
               if (l.ne.nvdepG) then
                  if(CurrentProcessorID.eq.0) write(*,*) 'probleme longueur vecture bdepG'
                  stop
               end if
               b2=0.d0
               l=0
               kk=0
               Do k=1,nv

                  If (idglob(k).eq.1.and.idg(k).eq.1) then
                     l=l+1
                     kk=kk+g

                     b2(l)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+kk) ! correspond a beta(g)

                  Else if (idglob(k).eq.2.and.idg(k).eq.1) then
                     l=l+1
                     kk=kk+g
                     Do m=1,np+1
                        b2(l+m-1)=b1(nvarprob*(ng-1)+nrisq+nvarxevt +ng+ng*(2*np+1)+nprmdep+(m-1)*ng+kk)
                     End do
                     l=l+np
                     kk=kk+np*ng

                  Else if (idglob(k).eq.1.and.idg(k).eq.0) then
                     kk=kk+1

                  Else if (idglob(k).eq.2.and.idg(k).eq.0) then
                     kk=kk+(np+1)

                  End if

               End do

!C     if(CurrentProcessorID.eq.0) write(*,*)'g=',g,'b2=',(b2(j),j=1,l)
!C------calcul de l'expo ---------------------------------
               mu=0.d0

               mu=MATMUL(Z,b0)+MATMUL(XglobG,b2)+MATMUL(XdepssG,bdepssG)     & 
             +MATMUL(XglobssG,b3)+MATMUL(XdepG,bdepG)+MATMUL(Xcont,b4)



!C ---- Bayesian estimates of ui
	    VC=0.d0
	    P=0.d0
               Ut1=0.d0
               if (g.eq.ng) then
                  Ut1=Ut
               else
                  Ut1=Ut*b1(nef+nvc+g)
               end if

               P=MATMUL(Z1,Ut1)
            VC=MATMUL(P,transpose(P))+Sigmaq+Vw+VCAR

	    Vi=0.d0
               jj=0
               do j=1,nmestot
                  do k=j,nmestot
                     jj=j+k*(k-1)/2
                     Vi(jj)=VC(j,k)
                  end do
               end do

         CALL DSINV(Vi,nmestot,eps,ier,det,0)
         if (ier.eq.-1) then
          !  if(CurrentProcessorID.eq.0) write(*,*)'dsinv vi ier',ier,' i=',i,' nmes(i)=',nmestot
          !  if(CurrentProcessorID.eq.0) write(*,*)'b1', (b1(k),k=1,npm)
            nnn=nmestot*(nmestot+1)/2
          !  if(CurrentProcessorID.eq.0) write(*,*)'Vi', (Vi(k),k=1,nnn)
            stop
         end if

!c     retransformation du vecteur Vi en matrice :
         VC=0.d0
         do j=1,nmestot
            do k=1,nmestot
               if (k.ge.j) then
                  VC(j,k)=Vi(j+k*(k-1)/2)
               else
                  VC(j,k)=Vi(k+j*(j-1)/2)
               end if
            end do
         end do
		
		yy=Y1all(i,:)
		ui_hat=0.d0
		Ut2=0.d0
		vv=0.d0
		ZZ=0.d0
		Ut2=MATMUL(Ut1,transpose(Ut1))
		ZZ=MATMUL(Ut2,transpose(Z1))
	       do q=1,nq
		  VV=MATMUL(VC,(yy-mu))
		  ui_hat(:)=MATMUL(ZZ,VV)
	       end do


!C Matrice de variance covariance

               VC=0.d0
               VC=MATMUL(P,transpose(P))+Sigmaq+Vw+VCAR

!c     Vi en vecteur
               jj=0
               do j=1,nmestot
                  do k=j,nmestot
                     jj=j+k*(k-1)/2
                     Vi(jj)=VC(j,k)
                  end do
               end do

               CALL DMFSD(Vi,nmestot,EPS,IER,0)
               if (ier.eq.-1) then
                  if(CurrentProcessorID.eq.0) write(*,*)'dmfsd vi ier',ier,' i=',i,' nmes(i)=',nmestot
                  if(CurrentProcessorID.eq.0) write(*,*)'b1', (b1(k),k=1,npmtot)
                  nnn=nmestot*(nmestot+1)/2
                  if(CurrentProcessorID.eq.0) write(*,*)'Vi', (Vi(k),k=1,nnn)
                  stop
               end if


!c     retransformation du vecteur Vi en matrice :
               VC=0.d0
               do j=1,nmestot
                  do k=1,j
                     VC(j,k)=Vi(k+j*(j-1)/2)
                  end do
               end do


!C     simulation des vecteurs normaux
               ysim2=0.d0
               ysim2cond=0.d0
               do l=1,nsim
                  ntronque=0
                  tronque=1
                  do while (tronque.eq.1)
                     tronque=0

                  ysim=0.d0
                  usim=0.d0
	         	  esim=0.d0
                     ysimcond=0.d0
                  do m=1,nmestot
                     SX=1.d0
                     call bgos(SX,0,usim(m),X2,0.d0)
                     call bgos(SX,0,esim(m),X2,0.d0)
                  end do
                  ysim=mu+MATMUL(VC,usim)
                  ysimcond=mu+MATMUL(Z1,ui_hat)+MATMUL(sigmaq,esim)
					
					ll=0
                  if (l.eq.1) then

                     jj=0
                     do q=1,nq
                        do j=1,nmes(i,q)
                           jj=jj+1
                           mui1(Navant_i(i)+jj,g)=ysim(jj)
                           muicond1(Navant_i(i)+jj,g)=ysimcond(jj)
                        end do
                     end do
                  end if



                  ntrtot=0
                  ll=0
                  Do q=1,nq
                     if (idord(q).le.0) then
                        aa1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &   
                    +ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot))/(1+exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &  
                     +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot)))
                        bb1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &    
                   +ncont*(nq-1)+2+ntrtot))/( 1.D0+ exp(b1(nvarprob*(ng-1)     &         
              +nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+2+ntrtot)))


                        bb1=aa1*(1.d0-aa1)*bb1

                        cc1=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &   
                    +ncont*(nq-1)+3+ntrtot)

                        dd1=abs(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &            
           +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+4+ntrtot))


                        aa=aa1*aa1*(1-aa1)/bb1-aa1
                        bb=aa*(1-aa1)/aa1
                        beta=beta_ln(aa,bb)

                        do j=1,nmes(i,q)
                           ytemp=ysim(ll+j)*dd1+cc1
                           ytemp2=ysimcond(ll+j)*dd1+cc1
                           if (ytemp.lt.0.or.ytemp.gt.1) then
!C                              if(CurrentProcessorID.eq.0) write(*,*)'ytemp en dehors de 0,1',ytemp
!C     &                             ,ysim(ll+j),dd1,cc1
!C                              tronque=1
!C                              ntronque=ntronque+1
                              if (ytemp.lt.0) then
                                 ytemp=0.d0
                              end if
                              if (ytemp.gt.1) then
                                 ytemp=1.d0
                              end if
!C     stop
                           end if

                           if (ytemp2.lt.0.or.ytemp2.gt.1) then
!C                              if(CurrentProcessorID.eq.0) write(*,*)'ytemp en dehors de 0,1',ytemp
!C     &                             ,ysim(ll+j),dd1,cc1
!C                              tronque=1
!C                              ntronque=ntronque+1
                              if (ytemp2.lt.0) then
                                 ytemp2=0.d0
                              end if
                              if (ytemp2.gt.1) then
                                 ytemp2=1.d0
                              end if
!C     stop
                           end if

                           ysim2(l,ll+j)=xinbta(aa,bb,beta,ytemp,ifault)
                           ysim2cond(l,ll+j)=xinbta(aa,bb,beta,ytemp2,ifault)
                        end do
                        ll=ll+nmes(i,q)
                        ntrtot=ntrtot+ntr
                     else if (idord(q).eq.1) then

                        zi_eval=0.d0
                        do j=-1,nztr+2
                           zi_eval(j)=zitr(j,q)
                        end do

                        bb=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &   
                    +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)


                        splaa=0.d0
                        do kk=2,ntrspl
                           splaa(kk-3)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &  
                        +ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot)*b1(nvarprob*(ng-1)+nrisq     &   
                       +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot)
                        end do

!C                        if(CurrentProcessorID.eq.0) write(*,*)
!C                        if(CurrentProcessorID.eq.0) write(*,*)'ijq',i,q,l
!C                        if(CurrentProcessorID.eq.0) write(*,*)




                        ntrtot=ntrtot+ntrspl
                        do j=1,nmes(i,q)
                           niter=0
                           diff=0.d0
                           ier=0
                           ysim2(l,ll+j)=INV_ISPLINES(ysim(ll+j),splaa,bb,zi_eval,ier,niter,diff)
						                              if ((ier.ne.1.and.diff.gt.1.d-3).or.ier.eq.3) then
                              if(CurrentProcessorID.eq.0) write(*,*)'problem INV_ISPLINE',  ysim2(l,ll+j),ysim(ll+j),i,j,q,l
                              if(CurrentProcessorID.eq.0) write(*,*)'zi_eval apres',nztr,(zi_eval(k) ,k=-1,nztr+2)
                              if(CurrentProcessorID.eq.0) write(*,*)'i',i,j,q,l,ll,ysim(ll+j)
                              if(CurrentProcessorID.eq.0) write(*,*)'bb',bb
                              if(CurrentProcessorID.eq.0) write(*,*)'spl',splaa
                              stop
                           end if
						   ysim2cond(l,ll+j)=INV_ISPLINES(ysimcond(ll+j),splaa,bb,zi_eval,ier,niter,diff)
                           if ((ier.ne.1.and.diff.gt.1.d-3).or.ier.eq.3) then
                              if(CurrentProcessorID.eq.0) write(*,*)'problem INV_ISPLINE',  ysim2cond(l,ll+j),ysimcond(ll+j),i,j,q,l
                              if(CurrentProcessorID.eq.0) write(*,*)'zi_eval apres',nztr,(zi_eval(k) ,k=-1,nztr+2)
                              if(CurrentProcessorID.eq.0) write(*,*)'i',i,j,q,l,ll,ysimcond(ll+j)
                              if(CurrentProcessorID.eq.0) write(*,*)'bb',bb
                              if(CurrentProcessorID.eq.0) write(*,*)'spl',splaa
                              stop
                           end if
                        end do
                        ll=ll+nmes(i,q)






                     else if (idord(q).eq.2) then
                        do j=1,nmes(i,q)
                           aa=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &    
                      +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                           if (ysim(ll+j).le.aa) then
                              ysim2(l,ll+j)=range(1+2*(q-1))
                           else
                              ysim2(l,ll+j)=range(2+2*(q-1))
                           end if
						   
						   if (ysimcond(ll+j).le.aa) then
                              ysim2cond(l,ll+j)=range(1+2*(q-1))
                           else
                              ysim2cond(l,ll+j)=range(2+2*(q-1))
                           end if
                        end do
                        ll=ll+nmes(i,q)
                        ntrtot=ntrtot+idord(q)-1
						
						
						
                     else if  (idord(q).gt.2) then
                        do j=1,nmes(i,q)
                           ysim2(l,ll+j)=range(1+2*(q-1))
						   ysim2cond(l,ll+j)=range(1+2*(q-1))
                              aa=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &  
                      +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
								do k=1,idord(q)-2
									bb=aa+b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     & 
									+ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1)*b1(nvarprob*(ng-1)+nrisq     &  
										+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1)
									if (ysim(ll+j).gt.aa.and.ysim(ll+j).le.bb)then
										ysim2(l,ll+j)=range(1+2*(q-1))+k
									end if
									if (ysimcond(ll+j).gt.aa.and.ysimcond(ll+j).le.bb)then
										ysim2cond(l,ll+j)=range(1+2*(q-1))+k
									end if
									aa=bb
								end do
								if (ysim(ll+j).gt.bb) then
									ysim2(l,ll+j)=range(2+2*(q-1))
								end if
									if (ysimcond(ll+j).gt.bb) then
									ysim2cond(l,ll+j)=range(2+2*(q-1))
								end if
                   end do
                        ll=ll+nmes(i,q)
                        ntrtot=ntrtot+idord(q)-1
                     end if   !   if idord
                  end do  !  l=1,nsim
               end do  !  while tronque

            if (ntronque.ne.0) then
               if(CurrentProcessorID.eq.0) write(*,*)'ntronque',i,ntronque
            end if
               end do

			   
               ll=0
               do q=1,nq
                  do j=1,nmes(i,q)
                     do l=1,nsim
                        ymarg(Navant_i(i)+ll+j,g)= ymarg(Navant_i(i)+ll+j,g)+ysim2(l,ll+j)
                        ycond(Navant_i(i)+ll+j,g)= ycond(Navant_i(i)+ll+j,g)+ysim2cond(l,ll+j)
                     end do
                     ymarg(Navant_i(i)+ll+j,g)=dble(ymarg(Navant_i(i)+ll+j,g)/dble(nsim))
                     ycond(Navant_i(i)+ll+j,g)=dble(ycond(Navant_i(i)+ll+j,g)/dble(nsim))
                     if (idord(q).le.0) then
                        ymarg(Navant_i(i)+ll+j,g)=range(1+2*(q-1))+ymarg(Navant_i(i)+ll+j,g)*(range(2+2*(q-1))-range(1+2*(q-1)))
                        ycond(Navant_i(i)+ll+j,g)=range(1+2*(q-1))+ycond(Navant_i(i)+ll+j,g)*(range(2+2*(q-1))-range(1+2*(q-1)))
                        if (ymarg(Navant_i(i)+ll+j,g).lt.range(1+2*(q-1))) then
                           ymarg(Navant_i(i)+ll+j,g)=range(1+2*(q-1))
                        end if
                        if (ymarg(Navant_i(i)+ll+j,g).gt.range(2+2*(q-1))) then
                           ymarg(Navant_i(i)+ll+j,g)=range(2+2*(q-1))
                        end if
                        if (ycond(Navant_i(i)+ll+j,g).lt.range(1+2*(q-1))) then
                           ycond(Navant_i(i)+ll+j,g)=range(1+2*(q-1))
                        end if
                        if (ycond(Navant_i(i)+ll+j,g).gt.range(2+2*(q-1))) then
                           ycond(Navant_i(i)+ll+j,g)=range(2+2*(q-1))
                        end if
                     end if
                  end do  !over j
                  ll=ll+nmes(i,q)
               end do  !over q
		   
			   !c sur G
            end do
!C ng.gt.1
         end if
      end do ! over i


      end subroutine evol_pred_marg_sim

! }}}

!C *************************************************************
!C
!C     evolution a posteriori (ordi)
!C
!C *************************************************************


! {{{

      subroutine evol_pred_marg_ordi(b,npm,nsim,ymarg,mui1)

      use commun
      use lois_normales
      use optim_parallele
      implicit none
      integer::npm,ier,nnn,ntrtot
      integer::nmestot,tronque
      integer :: i,j,k,l,jj,kk,q,m,g,nmoins,niter
      integer ::jj1,k1,q1,h,nsim,ll,ifault,ntronque

      double precision :: aa,bb,eps,X2,xinbta,SX,temp,inf, aa1,bb1,beta,beta_ln,cc1,dd1,ytemp,alnorm,diff



!C --------------   dimensions revues
      double precision,dimension(npm)::b
      double precision,dimension(npmtot)::b1
      double precision,dimension(nvc)::bh1
      double precision,dimension(ng)::pi
      double precision,dimension(nvarprob)::bprob,Xprob
      double precision,dimension(-2:ntrspl-3)::splaa
      double precision,dimension(nea,nea)::Ut,Ut1

      double precision,dimension(nvarglobssG)::b3
      double precision,dimension(2*np+1)::b0
      double precision,dimension(ncontdep+ncont)::sum
      double precision,dimension(nvdepssG)::bdepssG
      double precision,dimension(nvdepG)::bdepG
      double precision,dimension(nvarglobG)::b2
      double precision,dimension((ncontdep+ncont)*nq)::b4

      double precision, dimension(maxmestot)::mu,usim,ysim
      double precision,dimension(maxmestot,2*np+1)::Z
      double precision,dimension(maxmestot,nea)::Z1
      double precision,dimension(maxmestot,nvarglobssG)::XglobssG
      double precision,dimension(maxmestot,nvdepssG)::XdepssG
      double precision,dimension(maxmestot,(ncontdep+ncont)*nq)::Xcont
      double precision,dimension(maxmestot,nvarglobG)::XglobG
      double precision,dimension(maxmestot,nvdepG)::XdepG
      double precision,dimension(maxmestot,nea)::P
      double precision,dimension(maxmestot,maxmestot)::VC,VCAR,Vw,sigmaq
      double precision,dimension(maxmestot*(maxmestot+1)/2)::Vi
      double precision,dimension(-1:nztr+2)::zi_eval

      double precision,dimension(nsim,maxmestot)::ysim2
      double precision::INV_ISPLINES
      logical::upper


!C NON REVU

      double precision,dimension(ns*nq*mesure,ng)::ymarg,mui1





      eps=1.D-20

      upper=.false.

      k=0
      l=0
      DO j=1,npmtot
         if (ide(j).eq.1) then
            k=k+1
            b1(j)=b(k)
         else
            l=l+1
            b1(j)=bne(l)
         end if
      END DO

!C	  if(CurrentProcessorID.eq.0) write(*,*)'dans evol_pred 1'
!C ------ Ut : ---------------------------------------------
!C------- Creation de Ut la transformee de Cholesky ------


      Do j=1,nvc
         bh1(j)=b1(nef+j)
      End do

!C     if(CurrentProcessorID.eq.0) write(*,*)'bh1',(bh1(k),k=1,nvc)

      Ut=0.d0
      Ut1=0.d0
      If (idiag.eq.1) then
         do j=1,nea
            do k=1,nea
               if (j.eq.k) then
                  Ut(j,k)=bh1(j)
                  Ut1(j,k)=bh1(j)
               end if
            end do
         end do
      else if (idiag.eq.0) then
         do j=1,nea
            do k=j,nea
               Ut(k,j)=bh1(j+k*(k-1)/2)
               Ut1(k,j)=bh1(j+k*(k-1)/2)
            end do
         end do
      end if

!C fonction de risq a calculer pour chaque i et chaque g


!C------ boucle sur les sujets i-----------------------------
      mui1=0.d0
      ymarg=0.d0
      DO i=1,ns




!C------creation de la matrice des puissances du temps -----
         Z1=0.d0
         Z=0.d0
         jj=0
         Do q=1,nq
            Do j=1,nmes(i,q)
               jj=jj+1
               l=0
               Do k=1,np+1
                  Z(jj,k)=t(Navant_i(i)+jj)**(k-1)
                  If (k.gt.2) then
                     Z(jj,k)=Z(jj,k)*(0.1)**(k-2)
                  end if
                  If(idea(k).eq.1) then
                     l=l+1
                     Z1(jj,l)=Z(jj,k)
                  end if
               End do
            End do
         End do
         nmestot=jj


         XdepssG=0.d0
         l=0
         Do k=1,nvdep
            If (idvdep(k).eq.1) then
               l=l+1
               jj=0
               do q=1,nq
                  do j=1,nmes(i,q)
                     jj=jj+1
                     XdepssG(jj,l)=Xdep(Navant_i(i)+jj,k)
                  end do
               End do
               if (jj.ne.nmestot) then
                  if(CurrentProcessorID.eq.0) write(*,*)'probleme XdepssG(jj,k)'
                  stop
               end if
            end if
         end do


!C------ creation de XglobssG(j,k) (sans mixture) ----------

         XglobssG=0.d0

         l=0
         Do k=1,nv
            If (idglob(k).eq.1.and.idg(k).eq.0) then
               l=l+1
               Do j=1,nmestot
                  XglobssG(j,l)=X(i,k)
               End do

            Else if (idglob(k).eq.2.and.idg(k).eq.0) then
               jj=0
               l=l+1
               Do q=1,nq
                  Do j=1,nmes(i,q)
                     jj=jj+1
                     Do m=0,np



!C Attention : changer les Xglob avec interaction
                        if (m.ge.1) then
                           XglobssG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)*(0.1**(m-1))

                        else
                           XglobssG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)
                        end if

                     End do

                  End do
               End do
               l=l+np
            End if
         End do

!C creation du vecteur de variables explicatives dep du temps : XdepG avec mixture


         XdepG=0.d0
         l=0
         Do k=1,nvdep
            If (idvdep(k).eq.2) then
               l=l+1
               jj=0
               do q=1,nq
                  do j=1,nmes(i,q)
                     jj=jj+1
                     XdepG(jj,l)=Xdep(Navant_i(i)+jj,k)
                  end do
               End do
               if (jj.ne.nmestot) then
                  if(CurrentProcessorID.eq.0) write(*,*)'probleme XdepG(jj,k)'
                  stop
               end if
            end if
         end do


!C------ creation de XglobG(j,k) (avec mixture) ----------

         XglobG=0.d0

         l=0
         Do k=1,nv
            If (idglob(k).eq.1.and.idg(k).eq.1) then
               l=l+1
               Do j=1,nmestot
                  XglobG(j,l)=X(i,k)
               End do

            Else if (idglob(k).eq.2.and.idg(k).eq.1) then
               jj=0
               l=l+1
               Do q=1,nq
                  Do j=1,nmes(i,q)
                     jj=jj+1
                     Do m=0,np
                        if (m.ge.1) then
                           XglobG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)* (0.1**(m-1))
                        else
                           XglobG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)
                        end if
                     End do
                  End do
               End do
               l=l+np
            End if
         End do


         if (nvarglobssG+nvarglobG*ng.ne.nvarglob) then
            if(CurrentProcessorID.eq.0) write(*,*) 'pb : def nvarglobssG/G'
            stop
         end if


!C------ creation de Xcont(q,k,j) -----------------
!C creation des variables pour les contrastes

        Xcont=0.d0
         nmoins=0
         do q=1,nq
            do j=1,nmes(i,q)
               h=0
               do k=1,nvdep
                  if (idvdep(k).eq.3.or.idvdep(k).eq.4) then
                     h=h+1
                     Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)  =Xdep(Navant_i(i)+nmoins+j,k)
                  end if
               end do
               h=0
               do k=1,nv
                  if (idtest(k).eq.1) then
                     h=h+1
                     Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)=X(i,k)
                  end if
                  if (idtest(k).eq.2) then
                     Do m=0,np
                        h=h+1
                        if (m.ge.1) then
                           Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)=X(i,k)*(t(Navant_i(i)+nmoins+j)**m)*(0.1**(m-1))

                        else
                           Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)=X(i,k)*(t(Navant_i(i)+nmoins+j)**m)
                        end if
                     End do
                  end if
               end do
               if (h.ne.ncont) then
                  if(CurrentProcessorID.eq.0) write(*,*)'h=',h,'problem in cont 2'
                  stop
               end if
            End do
            nmoins=nmoins+nmes(i,q)
         end do
         if (nmoins.ne.nmestot) then
            if(CurrentProcessorID.eq.0) write(*,*)'problem in cont'
            stop
         end if



         VCAR=0.d0
         if (iAR.gt.0.and.QUECONT.eq.1) then
            jj=0
            Do q=1,nq
               jj1=0
               Do q1=1,nq
                  Do k=1,nmes(i,q)
                     do k1=1,nmes(i,q1)
                        if (iAR.eq.1) then
                           VCAR(jj+k,jj1+k1)=(b1(nef+nvc+ng-1+nq+1)**2)*min(t(Navant_i(i)+nmoins+jj+k),t(Navant_i(i)+nmoins+jj1+k1))
                        else if (iAR.eq.2) then
                           VCAR(jj+k,jj1+k1)=(b1(nef+nvc+ng-1+nq+1) **2)*exp(-b1(nef+nvc+ng-1+nq+2)     &   
                       *abs(t(Navant_i(i)+nmoins+jj+k)-t(Navant_i(i)+nmoins+jj1+k1)))
                        end if
                     end do
                  end do
                  jj1=jj1+nmes(i,q1)
               end do
               jj=jj+nmes(i,q)
            end do
         end if

         Sigmaq=0.d0
         Vw=0.d0
         jj=0.d0
         ntrtot=0
         Do q=1,nq
            Do k=1,nmes(i,q)
               do l=1,nmes(i,q)
                  Vw(jj+k,jj+l)=b1(nef+nvc+ng-1+q)**2
                  if (k.eq.l) then
                     Sigmaq(jj+k,jj+l)=b1(nef+nvc+ng-1+nq+nAR+q)**2
                  end if
               end do
            end do
            jj=jj+nmes(i,q)
         End do


!C	     if(CurrentProcessorID.eq.0) write(*,*)'i',i,'nmes',nmestot

!C------Vecteurs de parms sans mixtures -------------------
         bdepssG=0.d0
         jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)
         l=0
         kk=0
         Do k=1,nvdep
            If (idvdep(k).eq.1) then
               l=l+1
               kk=kk+1
               bdepssG(l)=b1(jj+kk)
            end if
            If (idvdep(k).eq.2) then
               kk=kk+ng
            end if
         end do
         if (l.ne.nvdepssG) then
            if(CurrentProcessorID.eq.0) write(*,*) 'probleme longueur vecture bdepssG'
            stop
         end if
         b3=0.d0
         b4=0.d0
         jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep
         l=0
         kk=0
         Do k=1,nv
            If (idglob(k).eq.1.and.idg(k).eq.0) then
               l=l+1
               kk=kk+1
               b3(l)=b1(jj+kk)  ! correspond a beta
            Else if (idglob(k).eq.2.and.idg(k).eq.0) then
               l=l+1
               kk=kk+1
               Do m=0,np
                  b3(l+m)=b1(jj+kk+m)
               End do
               l=l+np
               kk=kk+np
            Else if (idglob(k).eq.1.and.idg(k).eq.1) then
               kk=kk+ng
            Else if (idglob(k).eq.2.and.idg(k).eq.1) then
               kk=kk+ng*(np+1)
            End if
         End do

!C     b4 : contrastes sur les tests
         b4=0.d0
         sum=0.d0
         do j=1,ncontdep
            do q=1,nq-1
               b4((ncont+ncontdep)*(q-1)+j)=b1(nvarprob*(ng-1)+nrisq     &  
         +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+(j-1)*(nq-1)+q)
               sum(j)=sum(j)+b4((ncont+ncontdep)*(q-1)+j)
            end do
            b4((ncont+ncontdep)*(nq-1)+j)=-sum(j)
         end do
         do j=1,ncont
            do q=1,nq-1
               b4((ncont+ncontdep)*(q-1)+ncontdep+j)=b1(nvarprob*(ng-1)     &  
            +nrisq +nvarxevt +ng+ng*(2*np+1)+nprmdep+nvarglob     &       
    + ncontdep*(nq-1)+(j-1)*(nq-1)+q)
               sum(ncontdep+j)=sum(ncontdep+j)+b4((ncont+ncontdep)*(q-1)     &      
        +ncontdep+j)
            end do
            b4((ncont+ncontdep)*(nq-1)+ncontdep+j)=-sum(ncontdep+j)
         end do

         if (ng.eq.1) then

            VC=0.d0
            P=0.d0
            P=MATMUL(Z1,Ut)
            VC=0.d0
            VC=MATMUL(P,transpose(P))+Sigmaq+Vw+VCAR



!c     Vi en vecteur
           jj=0
            do j=1,nmestot
               do k=j,nmestot
                  jj=j+k*(k-1)/2
                  Vi(jj)=VC(j,k)
               end do
            end do

            CALL DMFSD(Vi,nmestot,EPS,IER,0)
            if (ier.eq.-1) then
               if(CurrentProcessorID.eq.0) write(*,*)'dmfsd vi ier',ier,' i=',i,' nmes(i)=',nmestot
               if(CurrentProcessorID.eq.0) write(*,*)'b1', (b1(k),k=1,npmtot)
               nnn=nmestot*(nmestot+1)/2
               if(CurrentProcessorID.eq.0) write(*,*)'Vi', (Vi(k),k=1,nnn)
               stop
            end if


!c     retransformation du vecteur Vi en matrice :
            VC=0.d0
            do j=1,nmestot
               do k=1,j
                  VC(j,k)=Vi(k+j*(j-1)/2)
               end do
            end do

            do k=1,2*np+1
               b0(k)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+k)
            end do
            mu=0.d0
            mu=MATMUL(Z,b0)+MATMUL(XglobssG,b3)+MATMUL(XdepssG,bdepssG)+MATMUL(Xcont,b4)


            jj=0
            do q=1,nq

               if (idord(q).le.0) then

                  cc1=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob &
                       +ncontdep*(nq-1)+ncont*(nq-1)+3+ntrtot)

                  dd1=abs(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &  
               +ncont*(nq-1)+4+ntrtot))
                  do j=1,nmes(i,q)
                     jj=jj+1

!==============================
! 2012/02/13
! attention, je pense que la ligne suivante etait pour comparer en lineaire l'approximation et la prediction directe ... peut-être est-ce a enlever.

                     mui1(Navant_i(i)+jj,1)=range(1+2*(q-1))+(mu(jj)*dd1+cc1)*(range(2+2*(q-1))-range(1+2*(q-1)))
                  end do
               end if
            end do

!C     simulation des vecteurs normaux



!C je ne touche pas au programme pour BETA, je rajoute ensuite le calcul pour ORDI sans integration numerique

            do l=1,nsim
               ntronque=0
               tronque=1
               do while (tronque.eq.1)
                  tronque=0
                  ysim=0.d0
                  usim=0.d0
                  do m=1,nmestot
                     SX=1.d0
                     call bgos(SX,0,usim(m),X2,0.d0)
                  end do


!C     if (i.eq.8.and.l.eq.2) then
!C     if(CurrentProcessorID.eq.0) write(*,*)'nmestot',nmestot,'usim',(usim(m),m=1,nmestot)
!C     end if

                  ysim=mu+MATMUL(VC,usim)

                  ll=0



                  ntrtot=0
                  Do q=1,nq
                     if (idord(q).le.0) then
                        aa1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &   
                    +ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot))/     &   
                    (1+exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)     &    
                   +nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot)))
                        bb1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &      
                 +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)     &    
                   +2+ntrtot))/( 1.D0+ exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &   
                    +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+2+ntrtot)))


                        bb1=aa1*(1.d0-aa1)*bb1

                        cc1=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &  
                     +ncont*(nq-1)+3+ntrtot)

                        dd1=abs(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &    
                   +ncont*(nq-1)+4+ntrtot))


                        aa=aa1*aa1*(1-aa1)/bb1-aa1
                        bb=aa*(1-aa1)/aa1
                        beta=beta_ln(aa,bb)

!C     if(CurrentProcessorID.eq.0) write(*,*)'aa-dd',aa,bb,cc1,dd1
!C
!C     if(CurrentProcessorID.eq.0) write(*,*)'ntrtot',ntrtot,'ntr',ntr
                        ntrtot=ntrtot+ntr
                        do j=1,nmes(i,q)
                           ytemp=ysim(ll+j)*dd1+cc1
                           if (ytemp.lt.0.or.ytemp.gt.1) then
!C     tronque=1
!C     ntronque=ntronque+1
!C     if(CurrentProcessorID.eq.0) write(*,*)'ytemp en dehors de 0,1',ytemp
!C     &                          ,ysim(ll+j),dd1,cc1
                              if (ytemp.lt.0) then
                                 ytemp=0.d0
                              end if
                              if (ytemp.gt.1) then
                                 ytemp=1.d0
                              end if
!C     stop
                           end if
                           ysim2(l,ll+j)=xinbta(aa,bb,beta,ytemp,ifault)
                        end do
                        ll=ll+nmes(i,q)
                     else if (idord(q).eq.1) then

                        zi_eval=0.d0
                        do j=-1,nztr+2
                           zi_eval(j)=zitr(j,q)
                        end do

                        bb=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &    
                   +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                        splaa=0.d0
                        do kk=2,ntrspl
                           splaa(kk-3)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     & 
                         +ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot)*b1(nvarprob*(ng-1)+nrisq     &      
                    +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot)
                        end do

!C     if(CurrentProcessorID.eq.0) write(*,*)
!c     if(CurrentProcessorID.eq.0) write(*,*)'ijq',i,q,l
!c     if(CurrentProcessorID.eq.0) write(*,*)


                        ntrtot=ntrtot+ntrspl
                        do j=1,nmes(i,q)
                           niter=0
                           diff=0.d0
                           ier=0
                           ysim2(l,ll+j)=INV_ISPLINES(ysim(ll+j),splaa,bb,zi_eval,ier,niter,diff)
                           if ((ier.eq.3).or.(ier.ne.1.and.diff.gt.1.d-3))  then
                              if(CurrentProcessorID.eq.0) write(*,*)'problem INV_ISPLINE',ysim2(l,ll+j),ysim(ll+j),i,j,q
                              if(CurrentProcessorID.eq.0) write(*,*)'zi_eval apres',nztr,(zi_eval(k),k=-1,nztr+2)
                              if(CurrentProcessorID.eq.0) write(*,*)'i',i,j,q,l,ll,ysim(ll+j)
                              if(CurrentProcessorID.eq.0) write(*,*)splaa,bb,zi_eval,ier,niter,diff
                              stop
                           end if
                        end do
                        ll=ll+nmes(i,q)


                   else if (idord(q).eq.2) then
                        do j=1,nmes(i,q)
                           aa=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &  
                        +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                           if (ysim(ll+j).le.aa) then
                              ysim2(l,ll+j)=range(1+2*(q-1))
                           else
                              ysim2(l,ll+j)=range(2+2*(q-1))
                           end if
                        end do
                        ll=ll+nmes(i,q)
                        ntrtot=ntrtot+idord(q)-1
                     else if  (idord(q).gt.2) then
                        do j=1,nmes(i,q)
                           ysim2(l,ll+j)=range(1+2*(q-1))
                           aa=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &   
                       +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                           do k=1,idord(q)-2
                              bb=aa+b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &   
                          +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1)*b1(nvarprob*(ng-1)+nrisq     &    
                         +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1)
                              if (ysim(ll+j).gt.aa.and.ysim(ll+j).le.bb) then
                                 ysim2(l,ll+j)=range(1+2*(q-1))+k
                              end if
                              aa=bb
                           end do
                           if (ysim(ll+j).gt.bb) then
                              ysim2(l,ll+j)=range(2+2*(q-1))
                           end if
                        end do
                        ll=ll+nmes(i,q)
                        ntrtot=ntrtot+idord(q)-1
                     end if
                  end do
               end do
!C     if (ntronque.ne.0) then
!C     if(CurrentProcessorID.eq.0) write(*,*)'ntronque',i,ntronque
!C     end if
            end do

            ll=0
            do q=1,nq
               do j=1,nmes(i,q)
                  do l=1,nsim
                     ymarg(Navant_i(i)+ll+j,1)= ymarg(Navant_i(i)+ll+j,1)+ysim2(l,ll+j)
                  end do
                  ymarg(Navant_i(i)+ll+j,1)=dble(ymarg(Navant_i(i)+ll+j,1)/dble(nsim))
                  if (idord(q).le.0) then
                     ymarg(Navant_i(i)+ll+j,1)=range(1+2*(q-1))+ymarg(Navant_i(i)+ll+j,1)*(range(2+2*(q-1))-range(1+2*(q-1)))
                     if (ymarg(Navant_i(i)+ll+j,1).lt.range(1+2*(q-1))) then
                        ymarg(Navant_i(i)+ll+j,1)=range(1+2*(q-1))
                     end if
                     if (ymarg(Navant_i(i)+ll+j,1).gt.range(2+2*(q-1))) then
                        ymarg(Navant_i(i)+ll+j,1)=range(2+2*(q-1))
                     end if
                  end if
               end do
               ll=ll+nmes(i,q)
            end do


            ntrtot=0
            jj=0
            Do q=1,nq
               if (idord(q).le.0) then
                  jj=jj+nmes(i,q)
                  ntrtot=ntrtot+ntr
               else if (idord(q).eq.1) then
                  jj=jj+nmes(i,q)
                  ntrtot=ntrtot+ntrspl
               else if (idord(q).eq.2) then
                  do j=1,nmes(i,q)
                     aa=(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &     
               +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)-mu(jj+j)  )/sqrt(VC(jj+j,jj+j))
                     ymarg(Navant_i(i)+jj+j,1)=range(2+2*(q-1)) -alnorm(aa,upper)
                  end do
                  jj=jj+nmes(i,q)
                  ntrtot=ntrtot+idord(q)-1
               else if  (idord(q).gt.2) then
                  do j=1,nmes(i,q)
                     inf=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &    
                +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                     aa=(inf-mu(jj+j))/sqrt(VC(jj+j,jj+j))
                     ymarg(Navant_i(i)+jj+j,1)=range(2+2*(q-1)) -alnorm(aa,upper)
                     do k=1,idord(q)-2
                        inf=inf+(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &    
                   +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1)*b1(nvarprob*(ng-1)+nrisq     &   
                    +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1))
                        aa=(inf-mu(jj+j))/sqrt(VC(jj+j,jj+j))
                        ymarg(Navant_i(i)+jj+j,1)=ymarg(Navant_i(i)+jj+j,1)-alnorm(aa ,upper)
                     end do
                  end do
                  jj=jj+nmes(i,q)
                  ntrtot=ntrtot+idord(q)-1
               end if

            end do

         else ! ng>1


            Xprob=0.d0
            Xprob(1)=1
            l=0
            do k=1,nv
               if (idprob(k).eq.1) then
                  l=l+1
                  Xprob(1+l)=X(i,k)
               end if
            end do
            if (l+1.ne.nvarprob) then
               if(CurrentProcessorID.eq.0) write(*,*)'problem nvarprob funcpa'
               stop
            end if

            pi=0.d0
            temp=0.d0
            Do g=1,ng-1
               do k=1,nvarprob
                  bprob(k)=b1((k-1)*(ng-1)+g)
               end do
               temp=temp+exp(DOT_PRODUCT(bprob,Xprob))
               pi(g)=exp(DOT_PRODUCT(bprob,Xprob))
               if ((pi(g)+1).eq.pi(g)) then
                  if(CurrentProcessorID.eq.0) write(*,*)'g',g,'pi(g)',pi(g)
                  if(CurrentProcessorID.eq.0) write(*,*)'bprob',(bprob(j),j=1,nvarprob)
                  if(CurrentProcessorID.eq.0) write(*,*)'Xprob',(Xprob(j),j=1,nvarprob)
                  stop
               end if
            end do

            pi(ng)=1/(1+temp)
            do g=1,ng

               if (g.ne.ng) then
                  pi(g)=pi(g)*pi(ng)
               end if

               if ((pi(g)+1).eq.pi(g)) then
                  if(CurrentProcessorID.eq.0) write(*,*)'g',g,'pi(g)',pi(g)
                  if(CurrentProcessorID.eq.0) write(*,*)'bprob',(bprob(j),j=1,nvarprob)
                  if(CurrentProcessorID.eq.0) write(*,*)'Xprob',(Xprob(j),j=1,nvarprob)

                  stop
               end if

!C ---- matrice de variance covariance


               VC=0.d0
               P=0.d0
               Ut1=0.d0
               if (g.eq.ng) then

                  Ut1=Ut

               else

                  Ut1=Ut*b1(nef+nvc+g)

               end if

               P=MATMUL(Z1,Ut1)
               VC=MATMUL(P,transpose(P))+Sigmaq+Vw+VCAR

!c     Vi en vecteur
               jj=0
               do j=1,nmestot
                  do k=j,nmestot
                     jj=j+k*(k-1)/2
                     Vi(jj)=VC(j,k)
                  end do
               end do

               CALL DMFSD(Vi,nmestot,EPS,IER,0)
               if (ier.eq.-1) then
                  if(CurrentProcessorID.eq.0) write(*,*)'dmfsd vi ier',ier,' i=',i,' nmes(i)=' ,nmestot
                  if(CurrentProcessorID.eq.0) write(*,*)'b1', (b1(k),k=1,npmtot)
                  nnn=nmestot*(nmestot+1)/2
                  if(CurrentProcessorID.eq.0) write(*,*)'Vi', (Vi(k),k=1,nnn)
                  stop
               end if


!c     retransformation du vecteur Vi en matrice :
               VC=0.d0
               do j=1,nmestot
                  do k=1,j
                     VC(j,k)=Vi(k+j*(j-1)/2)
                  end do
               end do

!C     if(CurrentProcessorID.eq.0) write(*,*)'i',i,'det',det

!C----- parametres avec mixtures ---------------------------

               b0=0.d0
               Do j=1,(2*np+1)
                  b0(j)=b1(nvarprob*(ng-1)+nrisq+nvarxevt +ng +(j-1)*ng+g)
               End do

!C     creation parametres mixture vdep du temps
               bdepG=0.d0
               jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)
               l=0
               kk=0
               Do k=1,nvdep
                  If (idvdep(k).eq.2) then
                     l=l+1
                     kk=kk+g
                     bdepG(l)=b1(jj+kk)
                  end if
                  If (idvdep(k).eq.1) then
                     kk=kk+1
                  end if
               end do
               if (l.ne.nvdepG) then
                  if(CurrentProcessorID.eq.0) write(*,*) 'probleme longueur vecture bdepG'
                  stop
               end if
               b2=0.d0
               l=0
               kk=0
               Do k=1,nv

                  If (idglob(k).eq.1.and.idg(k).eq.1) then
                     l=l+1
                     kk=kk+g

                     b2(l)=b1(nvarprob*(ng-1)+nrisq+nvarxevt        +ng+ng*(2*np+1)+nprmdep+kk) ! correspond a beta(g)

                  Else if (idglob(k).eq.2.and.idg(k).eq.1) then
                     l=l+1
                     kk=kk+g
                     Do m=1,np+1
                        b2(l+m-1)=b1(nvarprob*(ng-1)+nrisq     +nvarxevt +ng+ng*(2*np+1)+nprmdep+(m-1)*ng+kk)
                     End do
                     l=l+np
                     kk=kk+np*ng

                  Else if (idglob(k).eq.1.and.idg(k).eq.0) then
                     kk=kk+1

                  Else if (idglob(k).eq.2.and.idg(k).eq.0) then
                     kk=kk+(np+1)

                  End if

               End do

!C     if(CurrentProcessorID.eq.0) write(*,*)'g=',g,'b2=',(b2(j),j=1,l)



!C------calcul de l'expo ---------------------------------
               mu=0.d0

               mu=MATMUL(Z,b0)+MATMUL(XglobG,b2)+MATMUL(XdepssG,bdepssG)     &   
           +MATMUL(XglobssG,b3)+MATMUL(XdepG,bdepG)+MATMUL(Xcont,b4)





               jj=0
               do q=1,nq

               if (idord(q).le.0) then
                  cc1=b1(nvarprob*(ng-1)+nrisq+nvarxevt     &             
          +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+3+ntrtot)

                  dd1=abs(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &         
              +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+4+ntrtot))


                  do j=1,nmes(i,q)
                     jj=jj+1
                     mui1(Navant_i(i)+jj,g)=range(1+2*(q-1))+(mu(jj)*dd1+cc1)*(range(2+2*(q-1))-range(1+2*(q-1)))
                  end do
                  end if
               end do

!C     simulation des vecteurs normaux
               ysim2=0.d0
               do l=1,nsim
                  ntronque=0
                  tronque=1
                  do while (tronque.eq.1)
                     tronque=0

                  ysim=0.d0
                  usim=0.d0
                  do m=1,nmestot
                     SX=1.d0
                     call bgos(SX,0,usim(m),X2,0.d0)
                  end do

                  ysim=mu+MATMUL(VC,usim)

!C	              if(CurrentProcessorID.eq.0) write(*,*)'ysim',(ysim(j),j=1,nmestot)
                  ntrtot=0
                  ll=0
                  Do q=1,nq
                     if (idord(q).le.0) then
                        aa1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &   
                    +ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot))/(1+exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt     & 
                      +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot)))
                        bb1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &      
                 +ncont*(nq-1)+2+ntrtot))/( 1.D0+ exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &      
                 +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+2+ntrtot)))


                        bb1=aa1*(1.d0-aa1)*bb1

                        cc1=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     & 
                      +ncont*(nq-1)+3+ntrtot)

                        dd1=abs(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &  
                     +ncont*(nq-1)+4+ntrtot))


                        aa=aa1*aa1*(1-aa1)/bb1-aa1
                        bb=aa*(1-aa1)/aa1
                        beta=beta_ln(aa,bb)

                        do j=1,nmes(i,q)
                           ytemp=ysim(ll+j)*dd1+cc1
                           if (ytemp.lt.0.or.ytemp.gt.1) then
!C                              if(CurrentProcessorID.eq.0) write(*,*)'ytemp en dehors de 0,1',ytemp
!C     &                             ,ysim(ll+j),dd1,cc1
!C                              tronque=1
!C                              ntronque=ntronque+1
                              if (ytemp.lt.0) then
                                 ytemp=0.d0
                              end if
                              if (ytemp.gt.1) then
                                 ytemp=1.d0
                              end if
!C     stop
                           end if
                           ysim2(l,ll+j)=xinbta(aa,bb,beta,ytemp ,ifault)
                        end do
                        ll=ll+nmes(i,q)
                        ntrtot=ntrtot+ntr
                     else if (idord(q).eq.1) then

                    zi_eval=0.d0
                     do j=-1,nztr+2
                        zi_eval(j)=zitr(j,q)
                     end do

                     bb=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &    
                +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                     splaa=0.d0
                     do kk=2,ntrspl
                        splaa(kk-3)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &  
                     +ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot)*b1(nvarprob*(ng-1)+nrisq     &        
               +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot)
                     end do

!C                     if(CurrentProcessorID.eq.0) write(*,*)
!c                     if(CurrentProcessorID.eq.0) write(*,*)'ijq',i,q,l
!c                     if(CurrentProcessorID.eq.0) write(*,*)


                     ntrtot=ntrtot+ntrspl
                     do j=1,nmes(i,q)
                        niter=0
                        diff=0.d0
                        ier=0
                        ysim2(l,ll+j)=INV_ISPLINES(ysim(ll+j),splaa,bb,zi_eval,ier,niter,diff)
                        if ((ier.eq.3).or.(ier.ne.1.and.diff.gt.1.d-3))  then
                           if(CurrentProcessorID.eq.0) write(*,*)'problem INV_ISPLINE',ysim2(l,ll+j) ,ysim(ll+j),i,j,q
                           if(CurrentProcessorID.eq.0) write(*,*)'zi_eval apres',nztr,(zi_eval(k),k=-1,nztr+2)
                           if(CurrentProcessorID.eq.0) write(*,*)'i',i,j,q,l,ll,ysim(ll+j)
                           if(CurrentProcessorID.eq.0) write(*,*)splaa,bb,zi_eval,ier,niter,diff
                           stop
                        end if
                     end do
                     ll=ll+nmes(i,q)



                     else if (idord(q).eq.2) then
                        do j=1,nmes(i,q)
                           aa=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &  
                        +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                           if (ysim(ll+j).le.aa) then
                              ysim2(l,ll+j)=range(1+2*(q-1))
                           else
                              ysim2(l,ll+j)=range(2+2*(q-1))
                           end if
                        end do
                        ll=ll+nmes(i,q)
                        ntrtot=ntrtot+idord(q)-1
                     else if  (idord(q).gt.2) then
                        do j=1,nmes(i,q)
                           ysim2(l,ll+j)=range(1+2*(q-1))
                              aa=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &  
                           +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                           do k=1,idord(q)-2
                              bb=aa+b1(nvarprob*(ng-1)+nrisq +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &  
                           +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1)*b1(nvarprob *(ng-1)+nrisq     &
                             +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob +ncontdep*(nq-1)+ncont*(nq-1) +ntrtot+k+1)
                              if (ysim(ll+j).gt.aa.and.ysim(ll+j).le.bb)   then
                                 ysim2(l,ll+j)=range(1+2*(q-1))+k
                              end if
                              aa=bb
                           end do
                           if (ysim(ll+j).gt.bb) then
                              ysim2(l,ll+j)=range(2+2*(q-1))
                           end if
                        end do
                        ll=ll+nmes(i,q)
                        ntrtot=ntrtot+idord(q)-1
                     end if
                  end do
               end do



               if (ntronque.ne.0) then
                  if(CurrentProcessorID.eq.0) write(*,*)'ntronque',i,ntronque
               end if

               end do
               ll=0
               do q=1,nq
                  do j=1,nmes(i,q)
                     do l=1,nsim
                        ymarg(Navant_i(i)+ll+j,g)= ymarg(Navant_i(i)+ll+j,g)+ysim2(l,ll+j)
                     end do
                     ymarg(Navant_i(i)+ll+j,g)=dble(ymarg(Navant_i(i)+ll+j,g)/dble(nsim))
                     if (idord(q).le.0) then
                        ymarg(Navant_i(i)+ll+j,g)=range(1+2*(q-1))+ymarg(Navant_i(i)+ll+j,g) *(range(2+2*(q-1))-range(1+2*(q-1)))
                        if (ymarg(Navant_i(i)+ll+j,g).lt.range(1+2*(q-1))) then
                           ymarg(Navant_i(i)+ll+j,g)=range(1+2*(q-1))
                        end if
                        if (ymarg(Navant_i(i)+ll+j,g).gt.range(2+2*(q-1))) then
                           ymarg(Navant_i(i)+ll+j,g)=range(2+2*(q-1))
                        end if
                     end if
                  end do
                  ll=ll+nmes(i,q)
               end do


               ntrtot=0
               jj=0
               Do q=1,nq
                  if (idord(q).eq.0) then
                     jj=jj+nmes(i,q)
                     ntrtot=ntrtot+ntr
                  else if (idord(q).eq.1) then
                     jj=jj+nmes(i,q)
                     ntrtot=ntrtot+ntrspl
                  else if (idord(q).eq.-1) then
                     jj=jj+nmes(i,q)
                     ntrtot=ntrtot+ntr
                  else if (idord(q).eq.2) then
                     do j=1,nmes(i,q)
                        aa=(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &  
                     +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)-mu(jj+j)  )/sqrt(VC(jj+j,jj+j))
                        ymarg(Navant_i(i)+jj+j,g)=range(2+2*(q-1))-alnorm(aa,upper)
                     end do
                     jj=jj+nmes(i,q)
                     ntrtot=ntrtot+idord(q)-1
                  else if  (idord(q).gt.2) then
                     do j=1,nmes(i,q)
                        inf=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &  
                  +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)
                        aa=(inf-mu(jj+j))/sqrt(VC(jj+j,jj+j))
                        ymarg(Navant_i(i)+jj+j,g)=range(2+2*(q-1))-alnorm(aa ,upper)
                        do k=1,idord(q)-2
                           inf=inf+(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     & 
                         +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+k+1)*b1(nvarprob*(ng-1)+nrisq     &        
                  +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)     &  
                        +ntrtot+k+1))
                           aa=(inf-mu(jj+j))/sqrt(VC(jj+j,jj+j))
                           ymarg(Navant_i(i)+jj+j,g)=ymarg(Navant_i(i)+jj+j,g)-alnorm(aa,upper)
                        end do
                     end do
                     jj=jj+nmes(i,q)
                     ntrtot=ntrtot+idord(q)-1
                  end if
               end do

!c sur G
            end do

!C ng.gt.1
         end if

!C sur I
      end do


      end subroutine evol_pred_marg_ordi

! }}}

! ===============================================================
!
!
!  FONCTIONS TRANSFOS SPLINES
!
!
! ===============================================================


! {{{

      subroutine splines_tr ()

      use commun

      implicit none
      integer ::j,jj,q,i,k,n,l,ind,seuil

      integer,dimension(nq)::ntot
      double precision,dimension(mesure*ns)::ytemp
      double precision ::ht,htm,ht2,ht3,h,hh,h2,h3,h2n,hn,hht,pas,temp


      integer,dimension(ns)::nmestotcur

!      if(CurrentProcessorID.eq.0) write(*,*)'zitr',zitr(1,1),zitr(nztr,1)


!C debut de la subtroutine

      if (nztr <3) then
         if(CurrentProcessorID.eq.0) write(*,*) 'pas assez de noeuds : nz=',nztr
         stop
      end if

!C tri des donnees par test
!C boucle sur les tests
      ntot=0
      nmestotcur=0
      do q=1,nq
         Ytemp=0.d0
         jj=0
         do i=1,ns
            do j=1,nmes(i,q)
               jj=jj+1
               Ytemp(jj)=Y(Navant_i(i)+nmestotcur(i)+j)
            end do
         end do
         ntot(q)=jj


         if (idzitr.eq.0) then
!C     appel de la procedure de tri d'un vecteur
            ind=1
            do while (ind.eq.1)
               ind=0
               do i=1,ntot(q)-1
                  if (ytemp(i).gt.ytemp(i+1)) then
                     temp=ytemp(i)
                     ytemp(i)=ytemp(i+1)
                     ytemp(i+1)=temp
                     ind=1
                  end if
               end do
               !C     if(CurrentProcessorID.eq.0) write(*,*)'liste ',i,'=',(liste(i),i=1,n)
            end do
            
!            if(CurrentProcessorID.eq.0) write(*,*)'noeuds aux quantiles'
            
            !C     creation de la sequence de noeuds a partir des quantiles
            !C            zi(1,q)=ytemp(1)
            
            pas=floor(dble(ntot(q))/dble(nztr-1))
            seuil=pas
            j=1
            jj=0
            do jj=1,ntot(q)
               if (jj.eq.seuil.and.j.ne.(nztr-1)) then
                  j=j+1
                  zitr(j,q)=ytemp(jj)
                  seuil=seuil+pas
               end if
               
            end do
            
            if (j.ne.(nztr-1)) then
               if(CurrentProcessorID.eq.0) write(*,*)'probleme definition des quantiles'
               if(CurrentProcessorID.eq.0) write(*,*)'j+1=',j+1,'nz=',nztr
               stop
            end if
            
            zitr(-1,q)=zitr(1,q)
            zitr(0,q)=zitr(1,q)
            zitr(nztr+1,q)=zitr(nztr,q)
            zitr(nztr+2,q)=zitr(nztr,q)
         end if
         if (idzitr.eq.1) then
            
!            if(CurrentProcessorID.eq.0) write(*,*) 'noeuds equidistants'
            pas =(zitr(nztr,q)-zitr(1,q))/dble(nztr-1)
!            if(CurrentProcessorID.eq.0) write(*,*)'pas=',pas
            j=1
            do while(j.lt.nztr-1)
               j=j+1
               zitr(j,q)=zitr(j-1,q)+pas
            end do
            if (j.ne.(nztr-1)) then
               if(CurrentProcessorID.eq.0) write(*,*)'probleme definition des noeuds'
               if(CurrentProcessorID.eq.0) write(*,*)'j+1=',j+1,'nz-1=',nztr
               stop
            end if

            zitr(-1,q)=zitr(1,q)
            zitr(0,q)=zitr(1,q)
            zitr(nztr+1,q)=zitr(nztr,q)
            zitr(nztr+2,q)=zitr(nztr,q)

         end if
         if (idzitr.eq.2) then
            zitr(-1,q)=zitr(1,q)
            zitr(0,q)=zitr(1,q)
            zitr(nztr+1,q)=zitr(nztr,q)
            zitr(nztr+2,q)=zitr(nztr,q)
         end if


!C nombre de splines

         n=nztr+1


!         if(CurrentProcessorID.eq.0) write(*,*)'splines',(zitr(j,q),j=-1,ntrspl)
         
         l=0
         do i=1,ns
            do j=1,nmes(i,q)

!C ou se trouve la valeur de zi

               do k = 2,nztr
                  if ((Y(Navant_i(i)+nmestotcur(i)+j).ge.zitr(k-1,q)).and. (Y(Navant_i(i)+nmestotcur(i)+j).lt.zitr(k,q)))then
                     l=k-1
                  endif

               end do


               if (Y(Navant_i(i)+nmestotcur(i)+j).eq.zitr(nztr,q)) then
                  l=nztr-1
               end if

               ht2 = zitr(l+1,q)-Y(Navant_i(i)+nmestotcur(i)+j)
               htm= Y(Navant_i(i)+nmestotcur(i)+j)-zitr(l-1,q)
               ht = Y(Navant_i(i)+nmestotcur(i)+j)-zitr(l,q)
               ht3 = zitr(l+2,q)-Y(Navant_i(i)+nmestotcur(i)+j)
               hht = Y(Navant_i(i)+nmestotcur(i)+j)-zitr(l-2,q)
               h = zitr(l+1,q)-zitr(l,q)
               hh= zitr(l+1,q)-zitr(l-1,q)
               hn= zitr(l+1,q)-zitr(l-2,q)
               h2n=zitr(l+2,q)-zitr(l-1,q)
               h2= zitr(l+2,q)-zitr(l,q)
               h3= zitr(l+3,q)-zitr(l,q)

               if (Y(Navant_i(i)+nmestotcur(i)+j).ne.zitr(nztr,q)) then
                  mm2(Navant_i(i)+nmestotcur(i)+j) = (3.d0*ht2*ht2)/(hh*h*hn)
                  mm1(Navant_i(i)+nmestotcur(i)+j) = (3.d0*htm*ht2)/(h2n*hh*h) +(3.d0*ht*ht3)/(h2*h*h2n)
                  mm(Navant_i(i)+nmestotcur(i)+j)  = (3.d0*ht*ht)/(h3*h2*h)

               end if
               if (Y(Navant_i(i)+nmestotcur(i)+j).eq.zitr(nztr,q)) then

                  mm2(Navant_i(i)+nmestotcur(i)+j) = 0.d0
                  mm1(Navant_i(i)+nmestotcur(i)+j) = 0.d0
                  mm(Navant_i(i)+nmestotcur(i)+j)  = 3.d0/h

               end if

               if (mm2(Navant_i(i)+nmestotcur(i)+j).lt.0.or.mm1(Navant_i(i)+nmestotcur(i)+j)  &
                    .lt.0.or.mm(Navant_i(i)+nmestotcur(i)+j).lt.0) then
                  if(CurrentProcessorID.eq.0) write(*,*)'problem spline: -2=',mm2(Navant_i(i)+nmestotcur(i)+j),'-1='   &
                       , mm1(Navant_i(i)+nmestotcur(i)+j),'0=',mm(Navant_i(i)+nmestotcur(i)+j)
                  if(CurrentProcessorID.eq.0) write(*,*)'Y',Y(Navant_i(i)+nmestotcur(i)+j)
                  stop
               end if

!               if (i.eq.396) then
!                     if(CurrentProcessorID.eq.0) write(*,*)nmestotcur(i)+j,Y(Navant_i(i)+nmestotcur(i)+j),mm2(Navant_i(i)+nmestotcur(i)+j),mm1(Navant_i(i)+nmestotcur(i)+j),mm(Navant_i(i)+nmestotcur(i)+j)
!                  end if


               im2(Navant_i(i)+nmestotcur(i)+j)=hht*mm2(Navant_i(i)+nmestotcur(i)+j)/(3.d0)  &
                    + h2n*mm1(Navant_i(i)+nmestotcur(i)+j)/(3.d0) +h3*mm(Navant_i(i)+nmestotcur(i)+j)/(3.d0)
               im1(Navant_i(i)+nmestotcur(i)+j)=htm*mm1(Navant_i(i)+nmestotcur(i)+j)/(3.d0)  &
                    +h3*mm(Navant_i(i)+nmestotcur(i)+j)/(3.d0)

               im(Navant_i(i)+nmestotcur(i)+j)=ht*mm(Navant_i(i)+nmestotcur(i)+j)/(3.d0)


            end do

            nmestotcur(i)=nmestotcur(i)+nmes(i,q)
         end do


      end do

!      stop

      end subroutine splines_tr




!C ------------------------------------------------------------
!C
!C     INVERSE de I-splines
!C
!C     BY USING NEWTON-RAPHSON METHOD
!C ------------------------------------------------------------


!C ------------------------------------------------------------
!C
!C     INVERSE de I-splines
!C
!C     BY USING NEWTON-RAPHSON METHOD
!C ------------------------------------------------------------

      double precision FUNCTION INV_ISPLINES(X00,splaa,bb,zi_eval,istop,iter,eps)

      use commun
      implicit none
      double precision,dimension(-1:nztr+2)::zi_eval
      double precision::X0,X00,X1,fx0,f1x0,eps,bb,bb1
      double precision,dimension(-2:nztr)::splaa
      integer::iter,istop

      eps=1.d-6
      iter=1
      X0=1.d10
      call eval_splines(X0,fx0,f1x0,splaa,bb,zi_eval)

!C      if (fx0.eq.1.d-7.and.f1x0.eq.1.d-7) then
!C         INV_ISPLINES=fx0
!C         istop=3
!c         goto 1234
!C      end if


      if (X00.ge.fx0) then
         INV_ISPLINES=zi_eval(nztr)
         istop=1
         goto 1234
      end if
      X0=-1.d10
      call eval_splines(X0,fx0,f1x0,splaa,bb,zi_eval)
!C      if (fx0.eq.1.d-7.and.f1x0.eq.1.d-7) then
!C         INV_ISPLINES=fx0
!C         istop=3
!C         goto 1234
!C      end if
      if (X00.le.fx0) then
         INV_ISPLINES=zi_eval(1)
         istop=1
         goto 1234
      end if
      bb1=bb-X00
      X0=0
      call eval_splines(X0,fx0,f1x0,splaa,bb1,zi_eval)
!C      if (fx0.eq.1.d-7.and.f1x0.eq.1.d-7) then
!C         INV_ISPLINES=fx0
!C         istop=3
!C         goto 1234
!C      end if
      X1=X0-fx0/f1x0
      do while (ABS((X1-X0)/X0).GT.EPS.and.iter.lt.500)
         iter=iter+1
         X0=X1
         call eval_splines(X0,fx0,f1x0,splaa,bb1,zi_eval)
!C      if (fx0.eq.1.d-7.and.f1x0.eq.1.d-7) then
!C         INV_ISPLINES=fx0
!c         istop=3
!c         goto 1234
!C      end if
         X1=X0-fx0/f1x0
      end do
      INV_ISPLINES=zi_eval(1)+(zi_eval(nztr)-zi_eval(1))*exp(X1)/(1.d0+exp(X1))

      if (ABS((X1-X0)/X0).le.EPS) then
         istop=1
      else if (iter.ge.500) then
         istop=2
      else
         istop=3
      end if

      eps=ABS((X1-X0)/X0)

1234  continue

!C      if(CurrentProcessorID.eq.0) write(*,*)'INV_ISPLINES,X00,iter',INV_ISPLINES,X00,istop,iter


      return
      end function INV_ISPLINES



!C ------------------------------------------------------------
!C
!C     EVALUATION OF I-SPLINES and M-splines
!C ------------------------------------------------------------




      SUBROUTINE eval_splines(X00,Ispl,Mspl,splaa,bb,zi_eval)

      use commun

      implicit none
      double precision::X00,X0,Ispl,Mspl,som,bb
      double precision,dimension(-1:nztr+2)::zi_eval
      integer ::k,l,i
      double precision ::ht,htm,ht2,ht3,h,hh,h2,h3,h2n,hn,hht,mmeval,mm1eval,mm2eval,imeval,im1eval,im2eval

      double precision,dimension(-2:nztr)::splaa


!C ou se trouve la valeur de X


      X0=zi_eval(1)+(zi_eval(nztr)-zi_eval(1)) *(1.d0-1.d0/(1.d0+exp(X00)))

      do k = 2,nztr
         if ((X0.ge.zi_eval(k-1)).and. (X0.lt.zi_eval(k)))then
            l=k-1
         endif
      end do


      if (X0.eq.zi_eval(nztr)) then
         l=nztr-1
      end if

      ht2 = zi_eval(l+1)-X0
      htm= X0-zi_eval(l-1)
      ht = X0-zi_eval(l)
      ht3 = zi_eval(l+2)-X0
      hht = X0-zi_eval(l-2)
      h = zi_eval(l+1)-zi_eval(l)
      hh= zi_eval(l+1)-zi_eval(l-1)
      hn= zi_eval(l+1)-zi_eval(l-2)
      h2n=zi_eval(l+2)-zi_eval(l-1)
      h2= zi_eval(l+2)-zi_eval(l)
      h3= zi_eval(l+3)-zi_eval(l)

      if (h.eq.0.or.hh.eq.0.or.hn.eq.0.or.h2n.eq.0.or.h2.eq.0 .or.h3.eq.0)  then
!C         if(CurrentProcessorID.eq.0) write(*,*)'problem splines',h,hh,hn,h2n,h2,h3
!C         if(CurrentProcessorID.eq.0) write(*,*)'l',l,'nztr',nztr
!C         stop
      end if

      if (X0.ne.zi_eval(nztr)) then
         mm2eval = (3.d0*ht2*ht2)/(hh*h*hn)
         mm1eval = (3.d0*htm*ht2)/(h2n*hh*h)       +(3.d0*ht*ht3)/(h2*h*h2n)
         mmeval  = (3.d0*ht*ht)/(h3*h2*h)
      end if
      if (X0.eq.zi_eval(nztr)) then

         mm2eval = 0.d0
         mm1eval = 0.d0
         mmeval  = 3.d0/h

      end if

      if (mm2eval.lt.0.or.mm1eval.lt.0.or.mmeval.lt.0)  then
         if(CurrentProcessorID.eq.0) write(*,*)'problem spline: -2=',mm2eval,'-1=',   mm1eval,'0=',mmeval
         if(CurrentProcessorID.eq.0) write(*,*)'X0',X0,'l',l
         if(CurrentProcessorID.eq.0) write(*,*)'zi',(zi_eval(k),k=-1,nztr+2)
         if(CurrentProcessorID.eq.0) write(*,*)'splaa',splaa
         if(CurrentProcessorID.eq.0) write(*,*)'bb',bb
         if(CurrentProcessorID.eq.0) write(*,*)'X00',X00,'Ispl',Ispl,'Mspl',Mspl
!C         Mspl=1.d-7
!C         Ispl=1.d-7
!C         go to 587
         stop
      end if

      im2eval=hht*mm2eval/(3.d0)+ h2n*mm1eval/(3.d0) +h3*mmeval/(3.d0)
      im1eval=htm*mm1eval/(3.d0)    +h3*mmeval/(3.d0)
      imeval=ht*mmeval/(3.d0)

      som=0.d0
      if (l.gt.1) then
         do i=2,l
            som=som+splaa(i-3)
         end do
      end if

      Ispl=bb+ som +splaa(l-2)*im2eval  +splaa(l-1)*im1eval    +splaa(l)*imeval

      Mspl= (splaa(l-2)*mm2eval  +splaa(l-1)*mm1eval   +splaa(l)*mmeval)*  &
   (1.d0-1.d0/((1.d0+exp(X00))**2))* (zi_eval(nztr)-zi_eval(1))

! 587   continue

      end subroutine eval_splines



      subroutine estim_splines_ssstd(nsim,nz,zi,aa,test,transf)

      implicit none
      integer :: nsim,nz,lavant
      double precision,dimension(nsim) ::mm,mm1,mm2,im,im1,im2,transf,test
      double precision,dimension(-1:nz+2)::zi
      double precision,dimension(nz+3) ::aa,aa1,Xspl


      integer ::j,k,l,ntr
      double precision::min,max,pas
      double precision ::ht,htm,ht2,ht3,hht,h,hh,h2,h3,h2n,hn


!C nombre de parametres de noeuds + niveau
      ntr=nz+2

      test=0.d0
      aa1=0.d0
!C recuperation du vecteur de parametres : bb vecteur de niveau + aa1 : vecteur des splines
      aa1(1)=aa(1)
      do j=2,ntr
         aa1(j)=aa(j)*aa(j)
      end do

!C      if(CurrentProcessorID.eq.0) write(*,*)'aa',aa

!C matrice de transition pour delta-metho (carre des parms 2,..,ntr)


!C creation du vecteur de donnees

      min=zi(1)
      max=zi(nz)

!C     if(CurrentProcessorID.eq.0) write(*,*)'min=',min,'max=',max
      pas=(max-min)/dble(nsim-1)
      j=1
!      if(CurrentProcessorID.eq.0) write(*,*)'pas=',pas
      test(1)=min
      do while(j.lt.nsim)
         j=j+1
         test(j)=test(j-1)+pas
      end do
      test(nsim)=zi(nz)
      if (nsim.ne.j) then
!         if(CurrentProcessorID.eq.0) write(*,*)'probleme nsim'
         stop
      end if
!C     affectation des splines


      do j=1,nsim

!C ou se trouve la valeur
         l=0
         lavant=l
         if (test(j).eq.zi(nz)) then
            l=nz-1
         end if

         do k = 2,nz
            if ((test(j).ge.zi(k-1)).and. (test(j).lt.zi(k)))then
               l=k-1
            endif
         end do
         if (l.lt.1.or.l.gt.nz-1) then
!            if(CurrentProcessorID.eq.0) write(*,*)'probleme estim splines',l
!            if(CurrentProcessorID.eq.0) write(*,*)'j=',j,'test(j)',test(j)
            stop
         end if



               ht2 = zi(l+1)-test(j)
               htm= test(j)-zi(l-1)
               ht = test(j)-zi(l)
               ht3 = zi(l+2)-test(j)
               hht = test(j)-zi(l-2)
               h = zi(l+1)-zi(l)
               hh= zi(l+1)-zi(l-1)
               hn= zi(l+1)-zi(l-2)
               h2n=zi(l+2)-zi(l-1)
               h2= zi(l+2)-zi(l)
               h3= zi(l+3)-zi(l)

               if (test(j).ne.zi(nz)) then
                  mm2(j) = (3.d0*ht2*ht2)/(hh*h*hn)
                  mm1(j) = (3.d0*htm*ht2)/(h2n*hh*h) +(3.d0*ht*ht3)/(h2*h*h2n)
                  mm(j)  = (3.d0*ht*ht)/(h3*h2*h)
               end if
               if (test(j).eq.zi(nz)) then

                  mm2(j) = 0.d0
                  mm1(j) = 0.d0
                  mm(j)  = 3.d0/h

               end if

               im2(j)=hht*mm2(j)/(3.d0)+ h2n*mm1(j)/(3.d0) +h3*mm(j)/(3.d0)
               im1(j)=htm*mm1(j)/(3.d0)+h3*mm(j)/(3.d0)

               im(j)=ht*mm(j)/(3.d0)



!C-------- transformation et IC de la transformation :

            Xspl=0.d0
            Xspl(1)=1
            do k=2,l
               Xspl(k)=1
            end do
            Xspl(l+1)=im2(j)
            Xspl(l+2)=im1(j)
            Xspl(l+3)=im(j)


            transf(j)= dot_product(Xspl,aa1)

      end do



      end subroutine estim_splines_ssstd

! }}}

!C -----------------------------------------------------
!C
!C     ESTIMATION OF SUBJECT-SPECIFIC DEVIATIONS
!C                                        14/06/04
!C -----------------------------------------------------


! {{{

      subroutine pred_indiv(b,npm,pred_ind,Y_trans,mui,u_i,alpha_i)


      use commun
      use optim_parallele

      implicit none

      integer::npm,ier,nnn,ntrtot
      integer :: i,j,k,l,jj,kk,q,m,nmoins
      integer::nmestot
      integer ::jj1,k1,q1,h,ll,ii
      double precision::eps,aa,bb,det,ytemp,betai,aa1,bb1,cc1,dd1,som


!C --------------   dimensions revues
      double precision,dimension(npm)::b
      double precision,dimension(nvc)::bh1
      double precision,dimension(npmtot)::b1
      double precision,dimension(-2:ntrspl-3)::splaa
      double precision,dimension(nea,nea)::Ut

      double precision,dimension(nvarglobssG)::b3
      double precision,dimension(2*np+1)::b0
      double precision,dimension(ncont+ncontdep)::sum
      double precision,dimension(nvdepssG)::bdepssG
      double precision,dimension((ncontdep+ncont)*nq)::b4

      double precision, dimension(maxmestot)::mu,y1
      double precision,dimension(maxmestot,2*np+1)::Z
      double precision,dimension(maxmestot,nea)::Z1
      double precision,dimension(maxmestot,nvarglobssG)::XglobssG
      double precision,dimension(maxmestot,nvdepssG)::XdepssG
      double precision,dimension(maxmestot,(ncontdep+ncont)*nq)::Xcont
      double precision,dimension(maxmestot,nvarglobG)::XglobG
      double precision,dimension(maxmestot,nvdepG)::XdepG
      double precision,dimension(maxmestot,nea)::P
      double precision,dimension(maxmestot,maxmestot):: VC,VCAR,Sigmaq,Vw,Valea

      double precision,dimension(maxmestot*(maxmestot+1)/2)::Vi
      double precision, dimension(nq,maxmestot)::valea3
      double precision, dimension(nea,maxmestot)::valea2
      double precision,dimension(nea)::err4
      double precision,dimension(nq)::err5
      double precision,dimension(maxmestot)::err1,err2,err3



!C NON REVU

      double precision,dimension(ns,maxmestot)::pred_ind
      double precision,dimension(ns,nq,mesure)::Y_trans,mui
      double precision,dimension(ns,nea)::u_i
      double precision, dimension(ns,nq)::alpha_i



      eps=1.D-20



!C -----Vector of parameters
!

      k=0
      l=0
      DO j=1,npmtot
         if (ide(j).eq.1) then
            k=k+1
            b1(j)=b(k)
         else
            l=l+1
            b1(j)=bne(l)
         end if
      END DO
      if (k.ne.npm) then
         if(CurrentProcessorID.eq.0) write(*,*) 'problem with number of parameters'
         stop
      end if


!C------ Cholesky transformation
      Do j=1,nvc
         bh1(j)=b1(nef+j)
      End do
      Ut=0.d0
      If (idiag.eq.1) then
         do j=1,nea
            do k=1,nea
               if (j.eq.k) then
                  Ut(j,k)=bh1(j)
               end if
            end do
         end do
      else if (idiag.eq.0) then
         do j=1,nea
            do k=j,nea
               Ut(k,j)=bh1(j+k*(k-1)/2)
            end do
         end do
      end if


!C-----Loop across the subjects
      u_i=0.d0
      alpha_i=0.d0
      Y_trans=0.d0
      pred_ind=0.d0
      mui=0.d0
      DO i=1,ns

         Z1=0.d0
         Z=0.d0
         jj=0
         Do q=1,nq
            Do j=1,nmes(i,q)
               jj=jj+1
               l=0
               Do k=1,np+1
                  Z(jj,k)=t(Navant_i(i)+jj)**(k-1)
                  If (k.gt.2) then
                     Z(jj,k)=Z(jj,k)*(0.1)**(k-2)
                  end if
                  If(idea(k).eq.1) then
                     l=l+1
                     Z1(jj,l)=Z(jj,k)
                  end if
               End do
            End do
         End do
         nmestot=jj


         XdepssG=0.d0
         l=0
         Do k=1,nvdep
            If (idvdep(k).eq.1) then
               l=l+1
               jj=0
               do q=1,nq
                  do j=1,nmes(i,q)
                     jj=jj+1
                     XdepssG(jj,l)=Xdep(Navant_i(i)+jj,k)
                  end do
               End do
               if (jj.ne.nmestot) then
                  if(CurrentProcessorID.eq.0) write(*,*)'probleme XdepssG(jj,k)'
                  stop
               end if
            end if
         end do


!C------ creation de XglobssG(j,k) (sans mixture) ----------

         XglobssG=0.d0

         l=0
         Do k=1,nv
            If (idglob(k).eq.1.and.idg(k).eq.0) then
               l=l+1
               Do j=1,nmestot
                  XglobssG(j,l)=X(i,k)
               End do

            Else if (idglob(k).eq.2.and.idg(k).eq.0) then
               jj=0
               l=l+1
               Do q=1,nq
                  Do j=1,nmes(i,q)
                     jj=jj+1
                     Do m=0,np
!C Attention : changer les Xglob avec interaction
                        if (m.ge.1) then
                           XglobssG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)* (0.1**(m-1))
                        else
                           XglobssG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)
                        end if

                     End do

                  End do
               End do
               l=l+np
            End if
         End do

!C creation du vecteur de variables explicatives dep du temps : XdepG avec mixture


         XdepG=0.d0
         l=0
         Do k=1,nvdep
            If (idvdep(k).eq.2) then
               l=l+1
               jj=0
               do q=1,nq
                  do j=1,nmes(i,q)
                     jj=jj+1
                     XdepG(jj,l)=Xdep(Navant_i(i)+jj,k)
                  end do
               End do
               if (jj.ne.nmestot) then
                  if(CurrentProcessorID.eq.0) write(*,*)'probleme XdepG(jj,k)'
                  stop
               end if
            end if
         end do


!C------ creation de XglobG(j,k) (avec mixture) ----------

         XglobG=0.d0

         l=0
         Do k=1,nv
            If (idglob(k).eq.1.and.idg(k).eq.1) then
               l=l+1
               Do j=1,nmestot
                  XglobG(j,l)=X(i,k)
               End do

            Else if (idglob(k).eq.2.and.idg(k).eq.1) then
               jj=0
               l=l+1
               Do q=1,nq
                  Do j=1,nmes(i,q)
                     jj=jj+1
                     Do m=0,np
                        if (m.ge.1) then
                           XglobG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)* (0.1**(m-1))
                        else
                           XglobG(jj,l+m)=X(i,k)*(t(Navant_i(i)+jj)**m)
                        end if
                     End do
                  End do
               End do
               l=l+np
            End if
         End do


         if (nvarglobssG+nvarglobG*ng.ne.nvarglob) then
            if(CurrentProcessorID.eq.0) write(*,*) 'pb : def nvarglobssG/G'
            stop
         end if


!C------ creation de Xcont(q,k,j) -----------------
!C creation des variables pour les contrastes

        Xcont=0.d0
         nmoins=0
         do q=1,nq
            do j=1,nmes(i,q)
               h=0
               do k=1,nvdep
                  if (idvdep(k).eq.3.or.idvdep(k).eq.4) then
                     h=h+1
                     Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)    =Xdep(Navant_i(i)+nmoins+j,k)
                  end if
               end do
               h=0
               do k=1,nv
                  if (idtest(k).eq.1) then
                     h=h+1
                     Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)=X(i,k)
                  end if
                  if (idtest(k).eq.2) then
                     Do m=0,np
                        h=h+1
                        if (m.ge.1) then
                           Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h) =X(i,k)*(t(Navant_i(i)+nmoins+j)**m)*(0.1**(m-1))

                        else
                           Xcont(nmoins+j,(ncont+ncontdep)*(q-1)+h)=X(i,k)*(t(Navant_i(i)+nmoins+j)**m)
                        end if
                     End do
                  end if
               end do
               if (h.ne.ncont) then
                  if(CurrentProcessorID.eq.0) write(*,*)'h=',h,'problem in cont 2'
                  stop
               end if
            End do
            nmoins=nmoins+nmes(i,q)
         end do




         if (nmoins.ne.nmestot) then
            if(CurrentProcessorID.eq.0) write(*,*)'problem in cont'
            stop
         end if


         VCAR=0.d0
         if (iAR.gt.0.and.QUECONT.eq.1) then
            jj=0
            Do q=1,nq
               jj1=0
               Do q1=1,nq
                  Do k=1,nmes(i,q)
                     do k1=1,nmes(i,q1)
                        if (iAR.eq.1) then
                           VCAR(jj+k,jj1+k1)=(b1(nef+nvc+ng-1+nq+1)**2)*min(t(Navant_i(i)+jj+k),t(Navant_i(i)+jj1+k1))
                        else if (iAR.eq.2) then
                           VCAR(jj+k,jj1+k1)=(b1(nef+nvc+ng-1+nq+1)**2)*exp(-b1(nef+nvc+ng-1+nq+2)     & 
                         *abs(t(Navant_i(i)+jj+k)-t(Navant_i(i)+jj1+k1)))
                        end if
                     end do
                  end do
                  jj1=jj1+nmes(i,q1)
               end do
               jj=jj+nmes(i,q)
            end do
         end if





!C------- Creation du vecteur Y1 --------------------------
!C------------ et des matrices Kmat et Sigmaq ------

         Y1=0.d0
         Sigmaq=0.d0
         Vw=0.d0
         jj=0.d0
         ntrtot=0
         Do q=1,nq
            if (idord(q).le.0) then
               aa1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &     
         +ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot))/     &    
          (1+exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)     &  
            +nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+1+ntrtot)))
               bb1=exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &  
            +ncont*(nq-1)+2+ntrtot))/( 1.D0+ exp(b1(nvarprob*(ng-1)+nrisq+nvarxevt     &     
         +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1) +ncont*(nq-1) +2+ntrtot)))


               bb1=aa1*(1.d0-aa1)*bb1
!C     je change la modelisation
!C     bb1=aa1*bb1/4.d0

               cc1=b1(nvarprob*(ng-1)+nrisq+nvarxevt +ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &  
            +ncont*(nq-1)+3+ntrtot)

               dd1=abs(b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)     &     
         +ncont*(nq-1)+4+ntrtot))


               aa=aa1*aa1*(1-aa1)/bb1-aa1
               bb=aa*(1-aa1)/aa1
               Do k=1,nmes(i,q)
                  do l=1,nmes(i,q)
                     Vw(jj+k,jj+l)=b1(nef+nvc+ng-1+q)**2
                     if(k.eq.l) then
                        ytemp=(Y(Navant_i(i)+jj+k)-range(1+(q-1)*2))/(range(2+(q-1)*2)-range(1+(q-1)*2))
                        Y1(jj+k)=(betai(aa,bb,ytemp,tbeta)-cc1)/dd1
                        Valea3(q,jj+l)=b1(nef+nvc+q)**2
                        Y_trans(i,q,k)=Y1(jj+k)

                        if ((Y1(jj+k)+1).eq.Y1(jj+k)) then
                           if(CurrentProcessorID.eq.0) write(*,*)'Y1(jj+k).eq.NAN'
                           stop
                        end if

!C     OK pour la fonction de repartition
                        Sigmaq(jj+k,jj+l)=b1(nef+nvc+ng-1+nq+nAR+q)**2
                     end if
                  end do
               end do

               ntrtot=ntrtot+ntr
               jj=jj+nmes(i,q)

            else  if (idord(q).eq.1) then

              bb=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &     
         +ncontdep*(nq-1)+ncont*(nq-1)+ntrtot+1)

               do kk=2,ntrspl
                  splaa(kk-3)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob     &  
            +ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot) *b1(nvarprob*(ng-1)+nrisq     &  
            +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+ncontdep*(nq-1)+ncont*(nq-1)+kk+ntrtot)
               end do

               Do k=1,nmes(i,q)
                  do l=1,nmes(i,q)
                     Vw(jj+k,jj+l)=b1(nef+nvc+ng-1+q)**2
                     if(k.eq.l) then
                        ll=0
                        if (Y(Navant_i(i)+jj+k).eq.zitr(nztr,q)) then
                           ll=nztr-1
                        end if
                        som=0.d0
                        do kk = 2,nztr

                           if ((Y(Navant_i(i)+jj+k).ge.zitr(kk-1,q)).and.(Y(Navant_i(i)+jj+k).lt.zitr(kk,q)))then
                              ll=kk-1
                           end if
                        end do
                        if (ll.lt.1.or.ll.gt.nztr-1) then
                           if(CurrentProcessorID.eq.0) write(*,*) 'probleme dans funcpa splines'
                           if(CurrentProcessorID.eq.0) write(*,*) 'll=',ll,'Y=',Y(Navant_i(i)+jj+k)
                           stop
                        end if
                        if (ll.gt.1) then
                           do ii=2,ll
                              som=som+splaa(ii-3)
                           end do
                        end if

                        Y1(jj+k)=bb+ som +splaa(ll-2)*im2(Navant_i(i)+jj+k)+splaa(ll-1)  &
                             *im1(Navant_i(i)+jj+k)+splaa(ll)*im(Navant_i(i)+jj+k)

                        Valea3(q,jj+l)=b1(nef+nvc+q)**2
                        Y_trans(i,q,k)=Y1(jj+k)

                        if ((Y1(jj+k)+4).eq.Y1(jj+k)) then
                           if(CurrentProcessorID.eq.0) write(*,*)'Y1(jj+k).eq.NAN'
                           stop
                        end if

!C     OK pour la fonction de repartition
                        Sigmaq(jj+k,jj+l)=b1(nef+nvc+ng-1+nq+nAR+q)**2
                     end if
                  end do
               end do

               ntrtot=ntrtot+ntrspl
               jj=jj+nmes(i,q)



            else
               if(CurrentProcessorID.eq.0) write(*,*)'idord<>0 alors que pred_indiv que cont'
               stop
            end if

         End do


         Valea=0.d0
         VC=0.d0
         P=0.d0
         P=MATMUL(Z1,Ut)
         VC=MATMUL(P,transpose(P))+Sigmaq+Vw+VCAR
         Valea=MATMUL(P,transpose(P))+Vw+VCAR
         Valea2=MATMUL(Ut,transpose(P))

!c     Vi en vecteur
         jj=0
         do j=1,nmestot
            do k=j,nmestot
               jj=j+k*(k-1)/2
               Vi(jj)=VC(j,k)
            end do
         end do

         CALL DSINV(Vi,nmestot,eps,ier,det,0)
         if (ier.eq.-1) then
          !  if(CurrentProcessorID.eq.0) write(*,*)'dsinv vi ier',ier,' i=',i,' nmes(i)=',nmestot
          !  if(CurrentProcessorID.eq.0) write(*,*)'b1', (b1(k),k=1,npm)
            nnn=nmestot*(nmestot+1)/2
          !  if(CurrentProcessorID.eq.0) write(*,*)'Vi', (Vi(k),k=1,nnn)
            stop
         end if

!c     retransformation du vecteur Vi en matrice :
         VC=0.d0
         do j=1,nmestot
            do k=1,nmestot
               if (k.ge.j) then
                  VC(j,k)=Vi(j+k*(k-1)/2)
               else
                  VC(j,k)=Vi(k+j*(j-1)/2)
               end if
            end do
         end do



!C------Vecteurs de parms sans mixtures -------------------
         bdepssG=0.d0
         jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)
         l=0
         kk=0
         Do k=1,nvdep
            If (idvdep(k).eq.1) then
               l=l+1
               kk=kk+1
               bdepssG(l)=b1(jj+kk)
            end if
            If (idvdep(k).eq.2) then
               kk=kk+ng
            end if
         end do
         if (l.ne.nvdepssG) then
            if(CurrentProcessorID.eq.0) write(*,*) 'probleme longueur vecture bdepssG'
            stop
         end if
         b3=0.d0
         b4=0.d0
         jj=nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(2*np+1)+nprmdep
         l=0
         kk=0
         Do k=1,nv
            If (idglob(k).eq.1.and.idg(k).eq.0) then
               l=l+1
               kk=kk+1
               b3(l)=b1(jj+kk)  ! correspond a beta
            Else if (idglob(k).eq.2.and.idg(k).eq.0) then
               l=l+1
               kk=kk+1
               Do m=0,np
                  b3(l+m)=b1(jj+kk+m)
               End do
               l=l+np
               kk=kk+np
            Else if (idglob(k).eq.1.and.idg(k).eq.1) then
               kk=kk+ng
            Else if (idglob(k).eq.2.and.idg(k).eq.1) then
               kk=kk+ng*(np+1)
            End if
         End do

!C     b4 : contrastes sur les tests
         b4=0.d0
         sum=0.d0
         do j=1,ncontdep
            do q=1,nq-1
               b4((ncont+ncontdep)*(q-1)+j)=b1(nvarprob*(ng-1)+nrisq     &    
       +nvarxevt+ng+ng*(2*np+1)+nprmdep+nvarglob+(j-1)*(nq-1)+q)
               sum(j)=sum(j)+b4((ncont+ncontdep)*(q-1)+j)
            end do
            b4((ncont+ncontdep)*(nq-1)+j)=-sum(j)
         end do
         do j=1,ncont
            do q=1,nq-1
               b4((ncont+ncontdep)*(q-1)+ncontdep+j)=b1(nvarprob*(ng-1)     & 
             +nrisq +nvarxevt +ng+ng*(2*np+1)+nprmdep+nvarglob     &    
       + ncontdep*(nq-1)+(j-1)*(nq-1)+q)
               sum(ncontdep+j)=sum(ncontdep+j)+b4((ncont+ncontdep)*(q-1)+ncontdep+j)
            end do
            b4((ncont+ncontdep)*(nq-1)+ncontdep+j)=-sum(ncontdep+j)
         end do




         do k=1,2*np+1
            b0(k)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+k)
         end do
         mu=0.d0

         mu=MATMUL(Z,b0)+MATMUL(XglobssG,b3)+MATMUL(XdepssG,bdepssG) +MATMUL(Xcont,b4)


         jj=0
         do q=1,nq
            do j=1,nmes(i,q)
               jj=jj+1
               mui(i,q,j)=mu(jj)
            end do
         end do

!C------residual vector
         err1=0.d0
         err2=0.d0
         err3=0.d0
         err4=0.d0
         err5=0.d0
         err1=Y1-mu
         err2=MATMUL(VC,err1)
         err3=MATMUL(Valea,err2)
         err4=MATMUL(Valea2,err2)
         err5=MATMUL(Valea3,err2)
         do j=1,nmestot
            pred_ind(i,j)=err3(j)
         end do
         do j=1,nea
            u_i(i,j)=err4(j)
         end do
         do j=1,nq
            alpha_i(i,j)=err5(j)
         end do

      end do
      end subroutine pred_indiv

! }}}











! ===========================================================================
!
!
      !CREATION SPLINES POUR EVT
!
! ===========================================================================

! {{{

  
! }}}



! ===========================================================================
!
!
      !CREATION SPLINES POUR ESTIMATION SURVIE
!
! ===========================================================================



! {{{

 
! }}}





! ===========================================================================
!
!
      !CREATION SPLINES POUR ESTIMATION SURVIE 2 avec npg
!
! ===========================================================================



! {{{

        subroutine vraisSmooth(i,nlenb,b,bvexp,g,vexpi,tacti,valfin,tronc,score_gTT,score_gYT,SUTi0,STi0,postfit,trans,ni)
!	STi0= exp(-A01)
        use commun

        implicit none

!  double precision,intent(in)::thi,thj
        integer,intent(in)::i,nlenb,postfit,trans,ni
        integer::k,j,g,bparc,numact,numactbh,numactb,jb,delta_TiLi
        double precision,intent(out)::valfin,STi0,SUTi0
        double precision,intent(in),dimension(nva01+nva02+nva12)::vexpi
        double precision,intent(in),dimension(contevt2)::b
        double precision,intent(in),dimension(4)::tacti
        double precision,dimension(nvarxevt),intent(in)::bvexp
        double precision,dimension(nlenb)::bh
        double precision,dimension(ng,nevt):: risq,surv,surv0
        double precision,dimension(nprisq(1))::the01
	double precision,dimension(nprisq(2))::the02
	double precision,dimension(nprisq(3))::the12
        double precision,intent(out)::tronc,score_gTT,score_gYT
        double precision::res,res1,res2,vet01,vet02,vet12
        double precision::ri01,gl01,su01,ri02,gl02,su02,ri12,gl12,su12,temp,marting,martingYT

        bparc=0
        k=1
        bh=0.d0
        numact=1

        do k=1,nevt
           do j=1,nprisq(k)
              if (risqcom(k).ge.1) then !commun
                 bh(numact)=b(bparc+1)
                 numact=numact+1 !ok si risque commun
                 bparc=1+bparc
              endif
              if (risqcom(k).eq.0) then  !specific
                 bh(numact)=b(bparc+g)
                 numact=numact+1 !OK fonctionnel
                 bparc=bparc+ng
              endif
           enddo
        enddo

        bparc=0
!gestion vexp
        do k=1,nevt
           do j=1,nv
              if (idxevt(k*nv-nv+j).ge.1) then
                 if (idxevt(k*nv-nv+j).eq.1) then
                    bh(numact)=bvexp(bparc+1)
                    numact=numact+1
                    bparc=bparc+1
                 endif
                 if (idxevt(k*nv-nv+j).eq.2) then
                    bh(numact)=bvexp(bparc+g)
                    numact=numact+1
                    bparc=bparc+ng
                 endif
              endif
           enddo
        enddo

	do k=1,nprisq(1)
           the01(k)=bh(k)
        end do

        do k=1,nprisq(2)
           j = nprisq(1)+k
           the02(k)= bh(j) !très spécifique!
        end do

        do k=1,nprisq(3)
           j = nprisq(1) + nprisq(2) + k
           the12(k)=bh(j)
        end do

!i reserve au parcours des individus, j et k : indices de parcours local
! recuperation des parametres des fonctions de survie
!positivation au carre/ou expo deja faite dans funcpa



!---------- calcul de la vraisemblance ------------------
        res = 0.d0

        vet01 = 0.d0
        vet02 = 0.d0
        vet12 = 0.d0

!calcul des variables explicatrices
!nva : enregistres dans commun en fonction du .inf

        if(nva01.gt.0)then
           do j=1,nva01
              vet01 =vet01 + bh(nlenb-nva01-nva02-nva12+j)*dble(vexpi(j))
           end do
        endif
  
        if(nva02.gt.0)then
           do j=1,nva02
              vet02 =vet02 + bh(nlenb-nva02-nva12+j)*dble(vexpi(nva01+j))
           end do
        endif

        if(nva12.gt.0)then
           do j=1,nva12
              vet12 =vet12 + bh(nlenb-nva12+j)*dble(vexpi(nva01+nva02+j))
           end do
        endif

!if(ni.eq.50)then

!end if
!Fonctionnel pour n classes, proportionnalite par classe
      
        if (risqcom(1).eq.2.and.g.lt.ng.and.ng.gt.1) then
           vet01=vet01+b(size(b)-nevt*(ng-1)+g) !super spécifique à smoothhazard
        endif

        if (risqcom(2).eq.2.and.g.lt.ng.and.ng.gt.1) then
           vet02=vet02+b(size(b)-nevt*(ng-1)+g+(ng-1)) !super spécifique à smoothhazard
        endif

        if (risqcom(3).eq.2.and.g.lt.ng.and.ng.gt.1) then
           vet12=vet12+b(size(b)-nevt*(ng-1)+g+(ng-1)*2) !super spécifique à smoothhazard
        endif

       
!mis a l'esponentielle

        vet01 = dexp(vet01)
        vet02 = dexp(vet02)
        vet12 = dexp(vet12)

	STi0=0.d0
	SUTi0=0.d0
        if(troncature.eq.1)then
!calcul de la troncature gauche
           if(tacti(1).eq.0.d0)then
              tronc = 0.d0
           else   
              call fonct(tacti(1),typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
              call fonct(tacti(1),typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
              tronc = (gl01*vet01)+(gl02*vet02)
	      SUTi0=(su01**vet01)*(su02**vet02)*((-(gl01*vet01)-(gl02*vet02))**2-(gl01*vet01)-(gl02*vet02))
	      STi0=(su01**vet01)*(su02**vet02)
	   endif
        else
           tronc = 0.d0
        endif

        res1 = 0.d0

  !           censure a droite 01 et 02
        if(c(i).eq.1)then   
           call fonct(tacti(2),typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
           call fonct(tacti(4),typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
           res1 = -(gl01*vet01)-(gl02*vet02)

		   
		   
		   
        else if(c(i).eq.2)then  
              call fonct(tacti(4),typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
              call  qgauss1(tacti(2),tacti(3),size(the01),the01,size(the02),the02,size(the12),the12,res2,vet01,vet02,vet12)
              res1=dlog(res2*(su12**vet12))
        !observation 01 et censure a droite pour 12
		else if (c(i).eq.3)then  
              call fonct(tacti(2),typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
              call fonct(tacti(2),typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
              call fonct(tacti(2),typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
              res1 = -(gl01*vet01)-(gl02*vet02)+dlog(ri01*vet01)+(gl12*vet12)
  
              call fonct(tacti(4),typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
              res1 = res1 -(gl12*vet12)
		!censure par intervalle 01 et observation pour 12		 
        else if(c(i).eq.4)then  
              call fonct(tacti(4),typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
              call  qgauss1(tacti(2),tacti(3),size(the01),the01,size(the02),the02,size(the12),the12,res2,vet01,vet02,vet12)
				
			   res1=dlog(res2*(su12**vet12)*ri12*vet12) 
                   
                   ! res1=dlog(res2*(su12**vet12)*ri12*vet12)
        !observation 01 et observation pour 12
         else if(c(i).eq.5)then 
              call fonct(tacti(2),typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
              call fonct(tacti(2),typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
              call fonct(tacti(2),typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
              res1 = -(gl01*vet01)-(gl02*vet02)+dlog(ri01*vet01)+(gl12*vet12)
              call fonct(tacti(4),typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
                       res1 = res1 -(gl12*vet12) + dlog(ri12*vet12)
         !alive 
		 else if(c(i).eq.6)then
              call fonct(tacti(4),typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
              call fonct(tacti(4),typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
              call fonct(tacti(4),typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
              call  qgauss1(tacti(2),tacti(4),size(the01),the01,size(the02),the02,size(the12),the12,res2,vet01,vet02,vet12)
              res1 = (res2*(su12**vet12))+((su01**vet01)*(su02**vet02))
              res1 = dlog(res1)
		!passage 02?
         else if(c(i).eq.7)then
	      call fonct(tacti(4),typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
              call fonct(tacti(4),typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
              call fonct(tacti(4),typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)

              call  qgauss1(tacti(2),tacti(4),size(the01),the01,size(the02),the02,size(the12),the12,res2,vet01,vet02,vet12)
			  res1 = (res2*(su12**vet12)*ri12*vet12)+&
              ((su01**vet01)*(su02**vet02)*ri02*vet02)

			  !write(*,*) res1
              res1 = dlog(res1)

        endif
		
	if(postfit.eq.1)then

		score_gTT=0.d0
		score_gYT=0.d0
		if(interv_censoring.eq.0)then
		    if(trans.eq.1)then







		    else if(trans.eq.2)then




		    else if(trans.eq.3)then






		    end if

		      call fonct(tacti(4),typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
		      call fonct(tacti(4),typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)

		      score_gTT=(1-idm(i))*(su01**vet01)*(su02**vet02)*(ri02*vet02)**(idd(i))*&
					((idd(i)-gl01*vet01-gl02*vet02)**2-gl01*vet01-gl02*vet02)


		      score_gYT=(1-idm(i))*(su01**vet01)*(su02**vet02)*(ri02*vet02)**(idd(i))
		      if(trans.eq.1)then
		         score_gYT=score_gYT*(-gl01*vet01)
		      else if(trans.eq.2)then
		         score_gYT=score_gYT*(idd(i)-gl02*vet02)
		      else if(trans.eq.3)then
		         score_gYT=0.d0
		      end if

		      call fonct((tacti(2)+tacti(3))/2,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
		      call fonct((tacti(2)+tacti(3))/2,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
		      call fonct((tacti(2)+tacti(3))/2,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)

		      temp=idm(i)*(su01**vet01)*(su02**vet02)*(ri01*vet01)/(su12**vet12)
		      marting=(idd(i)+1-gl01*vet01-gl02*vet02+gl12*vet12)

		      if(trans.eq.1)then
		         martingYT=1-gl01*vet01
		      else if(trans.eq.2)then
		         martingYT=-gl02*vet02
		      else if(trans.eq.3)then
		         martingYT=gl12*vet12
		      end if

		      call fonct(tacti(4),typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)

		      temp=temp*(su12**vet12)*(ri12*vet12)**idd(i)
		      marting=(marting-gl12*vet12)**2-gl12*vet12

		      if(trans.eq.3)then
		         martingYT=martingYT+idd(i)-gl12*vet12
		      end if

		      call fonct((tacti(2)+tacti(3))/2,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
		      call fonct((tacti(2)+tacti(3))/2,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
		      call fonct((tacti(2)+tacti(3))/2,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)

		      marting=marting-gl01*vet01-gl02*vet02+gl12*vet12
		      score_gTT=(score_gTT+temp*marting)


		      score_gYT=(score_gYT+temp*martingYT)
		else
		      delta_TiLi=0
		      if(tacti(2).eq.tacti(4))delta_TiLi=1

		      call fonct(tacti(4),typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
		      call fonct(tacti(4),typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)

		      score_gTT=(1-idm(i))*(su01**vet01)*(su02**vet02)*(ri02*vet02)**(idd(i))*((idd(i)-&
				&gl01*vet01-gl02*vet02)**2-gl01*vet01-gl02*vet02)
		      call  qgauss_scoretest(tacti(2),tacti(3),tacti(4),size(the01),the01,size(the02),&
				&the02,size(the12),the12,res2,vet01,vet02,vet12,1,i,0)
		      score_gTT=score_gTT+(1-delta_TiLi)*res2

		      score_gYT=(1-idm(i))*(su01**vet01)*(su02**vet02)*(ri02*vet02)**(idd(i))
		      if(trans.eq.1)then
		         score_gYT=score_gYT*(-gl01*vet01)
		      else if(trans.eq.2)then
		         score_gYT=score_gYT*(idd(i)-gl01*vet01)
		      else if(trans.eq.3)then
		         score_gYT=0.d0
		      end if

		      call  qgauss_scoretest(tacti(2),tacti(3),tacti(4),size(the01),&
				&the01,size(the02),the02,size(the12),the12,res2,vet01,vet02,vet12,0,i,trans)
		      score_gYT=score_gYT+(1-delta_TiLi)*res2

		end if

	end if
        res =  res1

        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
           valfin=-1.d9
		   !print*, 'ne marche pas',res,c(i)
           goto 123
        end if
        
        valfin = res

123     continue

        return
  
      end subroutine vraisSmooth
!------------------------  FONCT  -------------------------
subroutine fonct(x,typerisq,nz,zi,dimP,p,risq,glam,surv)
!paramsurv pris en compte dans funcpa
  implicit none
        
  integer,intent(in) :: typerisq,dimP,nz
  double precision,intent(in),dimension(nz)::zi
  double precision,intent(in),dimension(dimP)::p
  double precision,intent(in)::x
  double precision,intent(out)::surv,risq,glam
  integer :: j,l
  double precision :: som


    if (x.le.0.d0) then
     surv = 1.d0
     glam = 0.d0
     risq = 0.d0
  endif
  
 if (typerisq.eq.2) then 
  surv = dexp(-(p(2)*x)**p(1))
  glam = (p(2)*x)**p(1)
  risq = p(1)*(p(2)**p(1))*(x**(p(1)-1.d0))
  
else if (typerisq.eq.1) then !Step function
		do j=1,nz-1
			som=p(1)*zi(1)
			if (x.ge.zi(j).and.x.le.zi(j+1)) then
				do l=1,j-1
					som=som+p(l)*(zi(l+1)-zi(l))			
				end do
				glam=som+p(j)*(x-zi(j))
				risq=p(j)
				surv=dexp(-glam)
			end if
		end do
end if


    return
end subroutine fonct



!================================  QGAUS1   ==========================
subroutine qgauss1(a,b,s1,the01,s2,the02,s3,the12,res,v01,v02,v12)
 use commun , only : typrisq,nz,zi01,zi02,zi12
  implicit none
integer, intent(in) :: s1,s2,s3
  double precision a,b,the01(s1),the02(s2),the12(s3)
  double precision dx,xm,xrr,w(5),xt(5),res,v01,v02,v12
  double precision xx,f1,su01,ri01,ri12,f2,su12,su02,ri02
  double precision gl01,gl12,gl02
   integer j


  save w,xt
  data w/0.2955242247d0,0.2692667193d0,0.2190863625d0,0.1494513491d0,0.0666713443d0/
  data xt/0.1488743389d0,0.4333953941d0,0.6794095682d0,0.8650633666d0,0.9739065285d0/
  
  
  xm = 0.5d0*(b+a)
  xrr = 0.5d0*(b-a)
  res = 0.d0
  if(a.eq.b)then
        call fonct(a,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
        call fonct(a,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
        call fonct(a,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
        res = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
  else
     do 11 j=1,5
        dx=xrr*xt(j)
        xx = xm+dx
        call fonct(xx,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
        call fonct(xx,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
        call fonct(xx,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
        f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
		       xx = xm-dx
        call fonct(xx,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
        call fonct(xx,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
        call fonct(xx,typrisq(3),nz(3),zi12, size(the12),the12,ri12,gl12,su12)
        f2 = ((su01**v01)*(su02**v02)*ri01*v01)/(su12**v12)
        res = res + w(j)*(f1+f2)
		


		
11      continue
  endif

  res = res*xrr

            
end subroutine qgauss1


subroutine qgauss_scoretest(a,b,c,s1,the01,s2,the02,s3,the12,res,v01,v02,v12,ind,ii,transition)
!ind =1 if score test between times-to-events, ind=0 if score test between Y and T
!transition : scire test independance entre Y et alpha_transition
 use commun , only : typrisq,nz,zi01,zi02,zi12,idd,idm
  implicit none
integer, intent(in) :: s1,s2,s3,ind,ii,transition
  double precision a,b,c,the01(s1),the02(s2),the12(s3)
  double precision dx,xm,xrr,w(5),xt(5),res,v01,v02,v12
  double precision xx,f1,su01,ri01,ri12,f2,su12,su02,ri02
  double precision gl01,gl12,gl02,temp
   integer j


  save w,xt
  data w/0.2955242247d0,0.2692667193d0,0.2190863625d0,0.1494513491d0,0.0666713443d0/
  data xt/0.1488743389d0,0.4333953941d0,0.6794095682d0,0.8650633666d0,0.9739065285d0/
  
  
  xm = 0.5d0*(b+a)
  xrr = 0.5d0*(b-a)
  res = 0.d0
  if(a.eq.b)then
        call fonct(a,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
        call fonct(a,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
        call fonct(a,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)

        res = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
	if(ind.eq.1) temp=(1-gl01*v01-gl02*v02+gl12*v12)
	if(ind.eq.0)then
	  if(transition.eq.1)temp=1-gl01*v01
	  if(transition.eq.2)temp=-gl02*v02
	  if(transition.eq.3)temp=gl12*v12
	end if

        call fonct(c,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
	res=res*(su12**v12)*(ri12*v12)**idd(ii)
	if(ind.eq.1) temp=(temp-gl12*v12)**2-gl12*v12
	if(ind.eq.0.and.transition.eq.3) temp=temp+idd(ii)-gl12*v12


	if(ind.eq.1)then
          call fonct(a,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
          call fonct(a,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
          call fonct(a,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
 	  temp=temp-gl01*v01-gl02*v02+gl12*v12
	end if
	res=res*temp
  else
     do 11 j=1,5
        dx=xrr*xt(j)
        xx = xm+dx
        call fonct(xx,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
        call fonct(xx,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
        call fonct(xx,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
        f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)

	if(ind.eq.1) temp=(1-gl01*v01-gl02*v02+gl12*v12+idd(ii))
	if(ind.eq.0)then
	  if(transition.eq.1)temp=1-gl01*v01
	  if(transition.eq.2)temp=-gl02*v02
	  if(transition.eq.3)temp=gl12*v12
	end if

        call fonct(c,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
	f1=f1*(su12**v12)*(ri12*v12)**idd(ii)
	if(ind.eq.1) temp=(temp-gl12*v12)**2-gl12*v12
	if(ind.eq.0.and.transition.eq.3) temp=temp+idd(ii)-gl12*v12


	if(ind.eq.1)then
          call fonct(xx,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
          call fonct(xx,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
          call fonct(xx,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
 	  temp=temp-gl01*v01-gl02*v02+gl12*v12
	end if
	f1=f1*temp

		       xx = xm-dx
        call fonct(xx,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
        call fonct(xx,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
        call fonct(xx,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
        f2 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)

	if(ind.eq.1) temp=(1-gl01*v01-gl02*v02+gl12*v12+idd(ii))
	if(ind.eq.0)then
	  if(transition.eq.1)temp=1-gl01*v01
	  if(transition.eq.2)temp=-gl02*v02
	  if(transition.eq.3)temp=gl12*v12
	end if

        call fonct(c,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
	f2=f2*(su12**v12)*(ri12*v12)**idd(ii)

	if(ind.eq.1) temp=(temp-gl12*v12)**2-gl12*v12
	if(ind.eq.0.and.transition.eq.3) temp=temp+idd(ii)-gl12*v12

	if(ind.eq.1)then 
          call fonct(xx,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
          call fonct(xx,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
          call fonct(xx,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
	  temp=temp-gl01*v01-gl02*v02+gl12*v12
	end if
	f2=f2*temp

        res = res + w(j)*(f1+f2)
		


		
11      continue
  endif

  res = res*xrr

            
end subroutine qgauss_scoretest
          
subroutine qg_risqcum(aa,kk,g,te,brisq,res_inst,res_cum,val_func,sex,bsex)

  use commun
  
  implicit none

!risque et risque cumule et incidence pour temps=te

 
  double precision,intent(in)::te,aa ! valeur temps
  integer,intent(in)::g,kk,sex
  double precision,intent(in),dimension(nevt)::bsex
!  double precision,intent(in),dimension(nz(1))::zi1
!  double precision,intent(in),dimension(nz(2))::zi2
!  double precision,intent(in),dimension(nz(3))::zi3
  double precision,intent(in),dimension(contevt2)::brisq
  double precision,dimension(nrisq)::bh !specifique smoothhazard
  double precision,intent(out)::res_inst,res_cum,val_func
  integer::bparc,numact,k,j,l
  double precision :: som,res
  double precision,dimension(nprisq(1))::the01
  double precision,dimension(nprisq(2))::the02
  double precision,dimension(nprisq(3))::the12
  double precision::ri01,gl01,su01,ri02,gl02,su02,ri12,gl12,su12
  bparc=0
  numact=1
  bh=0.d0
!recuperation des parametres de la classe g

  do k=1,nevt
     do j=1,nprisq(k)
        if (risqcom(k).ge.1) then
           bh(numact)=brisq(bparc+1)
           numact=numact+1
           bparc=1+bparc
        endif
        if (risqcom(k).eq.0) then ! class specific
           bh(numact)=brisq(bparc+g)
           numact=numact+1
           bparc=bparc+ng
        endif
     enddo
  enddo

!attention tres specifique
!suivant la transition calcul different, fonction versatile mais specifique

if (typrisq(kk).eq.2) then !Weibull
!    if(risqcom(kk).ne.0)then ! Common or proportional
	  if (kk.eq.1) then
		 res_inst=bh(1)*(bh(2)**bh(1))*(te**(bh(1)-1.d0))*dexp(bsex(1)*sex)
		 res_cum=(te*bh(2))**bh(1)*dexp(bsex(1)*sex)
	  endif

	  if (kk.eq.2) then
		 res_inst=bh(3)*(bh(4)**bh(3))*(te**(bh(3)-1.d0))*dexp(bsex(2)*sex)
		 res_cum=(te*bh(4))**bh(3)*dexp(bsex(2)*sex)
	  endif

	  if (kk.eq.3) then
		 res_inst=bh(5)*(bh(6)**bh(5))*(te**(bh(5)-1.d0))*dexp(bsex(3)*sex)
		 res_cum=(te*bh(6))**bh(5)*dexp(bsex(3)*sex)
	  endif
!     else ! class-specific
!	  if (kk.eq.1) then
!		 res_inst=bh(g)*(bh(ng+g)**bh(g))*(te**(bh(g)-1.d0))*dexp(bsex(1)*sex)
!		 res_cum=(te*bh(ng+g))**bh(g)*dexp(bsex(1)*sex)
!
!	  endif
!
!	  if (kk.eq.2) then
!		 res_inst=bh(2*ng+g)*(bh(3*ng+g)**bh(2*ng+g))*(te**(bh(2*ng+g)-1.d0))*dexp(bsex(2)*sex)
!		 res_cum=(te*bh(3*ng+g))**bh(2*ng+g)*dexp(bsex(2)*sex)
!	  endif
!
!	  if (kk.eq.3) then
!		 res_inst=bh(4*ng+g)*(bh(5*ng+g)**bh(4*ng+g))*(te**(bh(4*ng+g)-1.d0))*dexp(bsex(3)*sex)
!		 res_cum=(te*bh(5*ng+g))**bh(4*ng+g)*dexp(bsex(3)*sex)
!	  endif
!     end if

!     if(risqcom(kk).eq.2)then
!	  if(g.lt.ng.and.ng.gt.1) then
!			res_inst=res_inst*exp(brisq(size(brisq)-nevt*(ng-1)+(ng-1)*(kk-1)+g))
!			res_cum=res_cum*exp(brisq(size(brisq)-nevt*(ng-1)+(ng-1)*(kk-1)+g))
!	  endif
  
!calcules faits pour toutes les situations
!	  if ((kk.eq.1.or.kk.eq.2).and.g.lt.ng.and.ng.gt.1) then
!		 val_func=res_inst*dexp(-(te*bh(2))**bh(1)*dexp(bsex(1)*sex)*exp(brisq(size(brisq) & 
!			  -nevt*(ng-1)+(ng-1)*(0)+g))-(te*bh(4))**bh(3)*dexp(bsex(2)*sex)* & 
!			  dexp(brisq(size(brisq)-nevt*(ng-1)+(ng-1)+g)))
!		val_func=val_func/dexp(-(aa*bh(2))**bh(1)*dexp(bsex(1)*sex)*exp(brisq(size(brisq) & 
!			  -nevt*(ng-1)+(ng-1)*(0)+g))-(aa*bh(4))**bh(3)*dexp(bsex(2)*sex)* & 
!			  dexp(brisq(size(brisq)-nevt*(ng-1)+(ng-1)+g)))
!
!	  endif
!	  if ((kk.eq.1.or.kk.eq.2).and.g.eq.ng.and.ng.gt.1) then
!		 val_func=res_inst*dexp(-(te*bh(2))**bh(1)*dexp(bsex(1)*sex)-&
!			(te*bh(4))**bh(3)*dexp(bsex(2)*sex))
!		val_func=val_func/dexp(-(aa*bh(2))**bh(1)*dexp(bsex(1)*sex)-&
!			(aa*bh(4))**bh(3)*dexp(bsex(2)*sex))
!	  endif
!
!	  if (kk.eq.3.and.g.lt.ng.and.ng.gt.1) then
!		 val_func=dexp(-(te*bh(6))**bh(5)*dexp(bsex(3)*sex)*dexp(brisq(size(brisq) & 
!			  -nevt*(ng-1)+(ng-1)*2+g)))*res_inst
!		val_func=val_func/dexp(-(aa*bh(6))**bh(5)*dexp(bsex(3)*sex)*dexp(brisq(size(brisq) & 
!			  -nevt*(ng-1)+(ng-1)*2+g)))
!	  endif
!
!	  if (kk.eq.3.and.g.eq.ng.and.ng.gt.1) then
!		 val_func=dexp(-(te*bh(6))**bh(5)*dexp(bsex(3)*sex))*res_inst
!		val_func=val_func/dexp(-(aa*bh(6))**bh(5)*dexp(bsex(3)*sex))
!	  end if
 !    elseif(risqcom(kk).eq.1)then!common

	  if ((kk.eq.1.or.kk.eq.2)) then 
		 val_func=res_inst*dexp(-(te*bh(2))**bh(1)*dexp(bsex(1)*sex)-(te*bh(4))**bh(3)*dexp(bsex(2)*sex))
		val_func=val_func/dexp(-(aa*bh(2))**bh(1)*dexp(bsex(1)*sex)-(aa*bh(4))**bh(3)*dexp(bsex(2)*sex))
	  endif
	 
	  if (kk.eq.3) then
		 val_func=dexp(-(te*bh(6))**bh(5)*dexp(bsex(3)*sex))*res_inst
		val_func=val_func/dexp(-(aa*bh(6))**bh(5)*dexp(bsex(3)*sex))
	  endif


!     elseif(risqcom(kk).eq.0)then!class-spec
!	  if ((kk.eq.1.or.kk.eq.2)) then 
!		 val_func=res_inst*dexp(-(te*bh(ng+g))**bh(g)*dexp(bsex(1)*sex)&
!						-(te*bh(3*ng+g))**bh(2*ng+g)*dexp(bsex(2)*sex))
!
!		val_func=val_func/dexp(-(aa*bh(ng+g))**bh(g)*dexp(bsex(1)*sex)&
!
!						-(aa*bh(3*ng+g))**bh(2*ng+g)*dexp(bsex(2)*sex))
!	  endif
!	  if (kk.eq.3) then
!	!print*,'bh',bh,g
!		 val_func=dexp(-(te*bh(5*ng+g))**bh(4*ng+g)*dexp(bsex(3)*sex))*res_inst
!		val_func=val_func/dexp(-(aa*bh(5*ng+g))**bh(4*ng+g)*dexp(bsex(3)*sex))
!	  endif
 !   end if

else if (typrisq(kk).eq.1)then

	do k=1,nprisq(1)
           the01(k)=bh(k)
        end do

        do k=1,nprisq(2)
           j = nprisq(1)+k
           the02(k)= bh(j) !très spécifique!
        end do

        do k=1,nprisq(3)
           j = nprisq(1) + nprisq(2) + k
           the12(k)=bh(j)
        end do
		
	 if (kk.eq.1) then
		call fonct(te,typrisq(1),nz(1),zi01,size(the01),the01,res_inst,res_cum,su01)
		res_inst=res_inst*dexp(bsex(1)*sex)
		res_cum=res_cum*dexp(bsex(1)*sex)
		su01=su01**dexp(bsex(1)*sex)
	else if (kk.eq.2) then
		call fonct(te,typrisq(2),nz(2),zi02,size(the02),the02,res_inst,res_cum,su02)
		res_inst=res_inst*dexp(bsex(2)*sex)
		res_cum=res_cum*dexp(bsex(2)*sex)
		su02=su02**dexp(bsex(2)*sex)
	else if (kk.eq.3) then
		call fonct(te,typrisq(3),nz(3),zi12,size(the12),the12,res_inst,res_cum,su12)
		res_inst=res_inst*dexp(bsex(3)*sex)
		res_cum=res_cum*dexp(bsex(3)*sex)
		su12=su12**dexp(bsex(3)*sex)
	end if
		
	  if(risqcom(kk).eq.2.and.g.lt.ng.and.ng.gt.1) then
		res_inst=res_inst*exp(brisq(size(brisq)-nevt*(ng-1)+(ng-1)*(kk-1)+g))
		res_cum=res_cum*exp(brisq(size(brisq)-nevt*(ng-1)+(ng-1)*(kk-1)+g))
	endif
  
!calcules faits pour toutes les situations

		  if ((kk.eq.1.or.kk.eq.2).and.risqcom(kk).eq.2.and.g.lt.ng.and.ng.gt.1) then
		    call fonct(te,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
		    call fonct(te,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
			 val_func=res_inst*su01*exp(brisq(size(brisq) & 
				  -nevt*(ng-1)+(ng-1)*(0)+g))*su02* & 
				  exp(brisq(size(brisq)-nevt*(ng-1)+(ng-1)*(1)+g))
		    call fonct(aa,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
		    call fonct(aa,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
			 val_func=val_func/(su01*su02)
		  endif
		  if ((kk.eq.1.or.kk.eq.2).and.risqcom(kk).eq.2.and.g.eq.ng.and.ng.gt.1) then
		    call fonct(te,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
		    call fonct(te,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
			 val_func=res_inst*su01*su02
		    call fonct(aa,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
		    call fonct(aa,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
			 val_func=val_func/(su01*su02)
		  endif
		  if ((kk.eq.1.or.kk.eq.2).and.risqcom(kk).ne.2) then
		    call fonct(te,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
		    call fonct(te,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
			 val_func=res_inst*su01*su02
		    call fonct(aa,typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
		    call fonct(aa,typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
			 val_func=val_func/(su01*su02)
		  endif
		 
		  if (kk.eq.3.and.risqcom(kk).ne.2) then
			call fonct(te,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
			 val_func=su12*ri12
		    call fonct(aa,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
			 val_func=val_func/(su12)
		  endif
		  if (kk.eq.3.and.risqcom(kk).eq.2.and.g.lt.ng.and.ng.gt.1) then
		  call fonct(te,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
			 val_func=su12*exp(brisq(size(brisq) & 
				  -nevt*(ng-1)+(ng-1)*(2)+g))*ri12
		    call fonct(aa,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
			 val_func=val_func/(su12)
		  endif

		  if (kk.eq.3.and.risqcom(kk).eq.2.and.g.eq.ng.and.ng.gt.1) then
		  call fonct(te,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
			 val_func=su12*ri12
		    call fonct(aa,typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
			 val_func=val_func/(su12)
		  endif
end if

!si proportionnalite entre les classes

  return

end subroutine qg_risqcum


subroutine qg_inccum(a,te,brisq,g,kk,gaussgl,res_tot,sex,bsex)
  use commun
  
  implicit none

!quadrature gaussienne

 
  double precision,intent(in)::te,a! bornes integrale
  integer,intent(in)::g,kk,sex ! classe et evt
  double precision, intent(in),dimension(nevt)::bsex
  double precision,intent(in),dimension(contevt2)::brisq ! vecteur params
  double precision,intent(out)::res_tot
  double precision::res_inst,res_cum,val_func ! pour appel fonction
  integer::bparc,numact,k,j
  double precision::ftemp,xx,xr,xm,dx
  double precision,dimension(2,npg),intent(in)::gaussgl


  xm = 0.5d0*(te+a)
  xr = 0.5d0*(te-a)
  res_tot=0.d0 ! initialisation de la valeur a 0
  if(a.eq.te)then
     res_tot = 0.d0
  else
     do j=1,size(gaussgl(1,:))
        dx=xr*gaussgl(1,j)
        xx = xm+dx

        call qg_risqcum(a,kk,g,xx,brisq,res_inst,res_cum,val_func,sex,bsex)
        ftemp=val_func !on recupere que la valeur d'interêt correspondant a kk
        res_tot = res_tot + gaussgl(2,j)*ftemp !on somme l'integrale
     enddo
  endif
  res_tot = res_tot*xr

  return

end subroutine qg_inccum

subroutine postfit(bh,ppitest,ppi,ta,tb,tc,td,sex,bsex)
!Postfit survival : risque inst, risq cumul, incidence pour graphics :(Sexe=0,CEP=0)
	  use commun
  
	  implicit none
	integer::i,kk,g,nb_cl,numact,j,mmm

        !double precision,dimension(nb_cl)::class_X	
	double precision :: res_inst,res_cum,val_func,temp
        double precision,intent(in),dimension(ns,ng)::PPI,PPItest
        double precision,intent(in),dimension(nrisq)::bh
        real,intent(in),dimension(ns)::ta,tb,tc,td
	double precision, intent(in),dimension(nevt)::bsex
	integer,intent(in)::sex
	double precision,dimension(71)::temps !(100-65)/0.5+1  => 65,65.5,70...,100
	double precision,dimension(71,nevt)::alpha,Lambda,Incid

	open(9898,file='postfit_survival.txt')
	open(9899,file='SmoothHazard.txt')
	temps=0.d0
	alpha=0.d0
	Lambda=0.d0
	Incid=0.d0
	numact=0
	nb_cl=0

	  do j=1,size(temps)
	      temps(j)=0.5*(j-1)+65.d0
	  end do

      do kk=1,nevt
	do i=1,ns
	  if(X(i,3).eq.0.and.X(i,4).eq.0)then
  	    nb_cl=nb_cl+1
  	    do j=1,size(temps)
	        do g=1,ng
	           call qg_risqcum(temps(1),kk,g,temps(j),bh,res_inst,res_cum,val_func,sex,bsex)
	           alpha(j,kk)=alpha(j,kk)+ppitest(i,g)*res_inst
	           Lambda(j,kk)=Lambda(j,kk)+ppitest(i,g)*res_cum
	           Incid(j,kk)=Incid(j,kk)+ppitest(i,g)*val_func
	       end do
	     end do
	  end if
	end do
      end do

	do kk=1,nevt
    	    do i=1,size(temps)
		     write(9898,*) kk,temps(i),alpha(i,kk)/dble(nb_cl),&
				Lambda(i,kk)/dble(nb_cl),Incid(i,kk)/dble(nb_cl)
    	    end do
	end do

	DO i=1,ns
            mmm=1
	    temp=0.d0
            DO g=1,ng
               if (ppi(i,g).gt.temp) then
		  temp=ppi(i,g)
                  mmm=g
               end if
            END DO
		write(9899,*) X(i,1),X(i,2),ta(i),tb(i),tc(i),td(i),X(i,3),X(i,4)!,ppi(i,:)
	end do
end subroutine postfit

subroutine descript_class(ppi)
!Postfit survival
	  use commun
	implicit none
        double precision,intent(in),dimension(ns,ng)::ppi
 	integer::i,g
	double precision,dimension(ng,9) :: resultat
	integer,dimension(ns) :: class
	integer,dimension(ng) :: a,b
	double precision :: temp

!resultat : nb_indiv(1), mean age_entree (2),sexe(3),cep(4),nb_dements(5),nb deces sans demence(6), nb deces post demence, age deces post demence
write(8,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(8,*) 'g','nb_indiv','mean age entry','sexe','CEP','nb_demented','nb death without dem',&
'nb death with dem','mean age death with dem'
	a=0
	b=0
	resultat=0
	class=0
	do i=1,ns
	  temp=ppi(i,1)
	  class(i)=1
	  do g=2,ng
	    if(ppi(i,g).gt.temp) then
		class(i)=g
		temp=ppi(i,g)
	  end if
	  end do
	end do
	
	do i=1,ns
	do g=1,ng
	  if(class(i).eq.g)then
		resultat(g,1)=resultat(g,1)+1
		resultat(g,2)=resultat(g,2)+t0(i)
		resultat(g,3)=resultat(g,3)+X(i,3)
		resultat(g,4)=resultat(g,4)+X(i,4)
		resultat(g,5)=resultat(g,5)+X(i,1)
		if(X(i,1).eq.0.and.X(i,2).eq.1) resultat(g,6)=resultat(g,6)+X(i,2)
		if(X(i,1).eq.1.and.X(i,2).eq.1) then
			resultat(g,7)=resultat(g,7)+X(i,2)
			resultat(g,8)=resultat(g,8)+t3(i)
			a(g)=a(g)+1
		end if
		if(X(i,1).eq.1) then
			resultat(g,9)=resultat(g,9)+t2(i)
			b(g)=b(g)+1
		end if
	  end if
	end do
	end do

	do g=1,ng
	   resultat(g,2)=resultat(g,2)/dble(resultat(g,1))
	   resultat(g,8)=resultat(g,8)/dble(a(g))
	   resultat(g,9)=resultat(g,9)/dble(b(g))
	   write(8,*) g,real(resultat(g,:))
	end do

	write(8,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	do g=1,ng
	 do i=3,7
	   resultat(g,i)=resultat(g,i)/resultat(g,1)
	end do
	   !resultat(g,8)=resultat(g,8)/dble(a(g))
	   write(8,*) g,real(resultat(g,:))
	end do
	write(8,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
end subroutine descript_class


subroutine plot_postfit(n,b1,pi,nplot,time_plot,eq_plot)
	  use commun
	  implicit none

	integer::i,j,k,g,ntemps,bparc,numact
	integer,intent(in)::nplot,n
	double precision, dimension(2,nplot),intent(in):: time_plot
	integer, dimension(2,nplot),intent(in)::eq_plot
	real, dimension(ng),intent(in):: pi
	double precision, dimension(nprisq(1)):: the01
	double precision, dimension(nprisq(2)):: the02
	double precision, dimension(nprisq(3)):: the12
	double precision,dimension(nrisq)::bh
	double precision, dimension(n),intent(in):: b1
	double precision :: tinf,tsup,pas,ri01,gl01,su01,surv,survtot
	double precision :: ri02,gl02,su02,ri12,gl12,su12,temp
	double precision, dimension(:),allocatable::temps
	double precision, dimension(:,:),allocatable::traj
	double precision, dimension(2*np+1):: b0
	


	tinf=0.d0
	tsup= 3.5d0
	pas= 0.05d0
        ntemps=floor((Tsup-Tinf)/pas)+1

        allocate(temps(ntemps),traj(nplot,ntemps))
  
       temps(1)=Tinf
        if (ntemps.gt.1) then 
           do j=2,ntemps
              temps(j)=temps(j-1)+pas
           end do
           temps(ntemps)=Tsup
        end if

  	open(0103,file='postfit_trajectories.txt')

	do i=1,nplot
	  do g=1,ng

          bparc=0
          k=1
          bh=0.d0
          numact=1

          do k=1,nevt
             do j=1,nprisq(k)
                if (risqcom(k).ge.1) then !commun
                   bh(numact)=b1(bparc+1)
                   numact=numact+1 !ok si risque commun
                   bparc=1+bparc
                endif
                if (risqcom(k).eq.0) then  !specific
                   bh(numact)=b1(bparc+g)
                   numact=numact+1 !OK fonctionnel
                   bparc=bparc+ng
                endif
             enddo
          enddo
  
  	  do k=1,nprisq(1)
             the01(k)=bh(k)
          end do

          do k=1,nprisq(2)
             j = nprisq(1)+k
             the02(k)= bh(j) !très spécifique!
          end do

          do k=1,nprisq(3)
             j = nprisq(1) + nprisq(2) + k
             the12(k)=bh(j)
          end do
	  call fonct(time_plot(1,i),typrisq(1),nz(1),zi01,size(the01),the01,ri01,gl01,su01)
	  call fonct(time_plot(2,i),typrisq(2),nz(2),zi02,size(the02),the02,ri02,gl02,su02)
	  call fonct(time_plot(2,i),typrisq(3),nz(3),zi12,size(the12),the12,ri12,gl12,su12)
	  surv=su01*(ri01**eq_plot(1,i))*(su02*(ri02**eq_plot(2,i))+su12*(ri12**eq_plot(2,i)))
	  survtot=survtot + surv*pi(g)
	  do k=1,2*np+1
              b0(k)=b1(nvarprob*(ng-1)+nrisq+nvarxevt+ng+ng*(k-1)+g)
          end do
	  do j=1,ntemps
	     temp= pi(g)*(b0(1)+b0(2)*((temps(j)-65.d0)/10.d0)+b0(3)*((temps(j)-65.d0)/10.d0)**2/10.d0) !! a verifier !!
	     traj(i,j)=traj(i,j) +temp*surv
	  end do
      end do
    end do

      do j=1,ntemps
	do i=1,nplot
	  traj(i,j)=traj(i,j)/survtot
	end do
	write(0103,*)  temps(j), (traj(k,j),k=1,nplot)
      end do

	close(0103)

end subroutine plot_postfit
