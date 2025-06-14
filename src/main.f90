       module main_mod
       contains

      subroutine MinGCE(endoftime,sigmat,sigmah,psfr,pwind,&
     &     delay,time_wind)
      
      use io_mod
      use io_aux_mod
      use interpolation_mod
      implicit none
      integer nmax,elem
      PARAMETER (NMAX=15000)
      parameter (elem=33)
      character*4 nameout
      character*25 inputyield(5)
      character*10 siess
      CHARACTER(100) :: num1char,num2char,num3char,num4char,num5char
      CHARACTER(100) :: num6char,num7char
      
      real Q(elem), Norm
      integer k,j3,ss,ss2
      real XMU(200)
      real amu(115),snian
      real QQN(elem,NMAX,100)
      real stars(NMAX),remn(NMAX),gas(NMAX),zeta(NMAX)
      real W(23,35,15),WH(23,35,15)
      real massa(35),massac(13),massac2(4),massas(7)
      real MBa(4),WBa(4,3)
      real MSr(4),WSr(4,3)
      real MY(4),WY(4,3)
      real MLa(4),WLa(4,3)
      real MRb(4),WRb(4,3)
      real MZr(4),WZr(4,3)
      real MEu(4),WEu(4,3)      

   
      real mstars(1000),mstars1(1000)
      real alls(NMAX),all(NMAX),hot(NMAX),wind(NMAX),oldstars(NMAX)
      real threshold,tau
      real multi1(1000),multi2(1000)
      integer hottime
!     external ranf,assa
      integer i,ii,j,jj,t3,n,p,step,step3
      integer z
      integer tt,t
      real tdead(1000)
      real SFR,windist
      real BINMAX(1000),BINMAX1(1000),Hecore
      real SNIAnum(NMAX)
      
      real pSFR,superf,Pwind,Pinfall
      real mux,mu1,mu2,mmu,mumin,mumin2,mun
      integer endoftime,IMF,tautype,ninputyield
      integer lowmassive,MM
      integer ence,delay
      real temp1,temp2,boh,pi,spalla(NMAX)
c      real ampli,pi,sinsfr,azi,period cspiral

      
      real SigmaR,RD,SD,TauR
      real SigmaH,TauH,T1,kappa,Rm
      real Sigmasun,Tsolar,mingas
      real Raggio,sigmat
      integer time_wind
      real Ini(elem)

      
      real remntot, hecoretot, starstot, QIspecial(elem,1000)
      real difftot, hecores(1000)
      integer binmass
     
      real UM1,UM2,UM3,UM4
      real A,B,C,D
      real M1,M2,M3,M4

      COMMON UM1,UM2,UM3,UM4
      COMMON A,B,C,D
      COMMON M1,M2,M3,M4
      common ninputyield
      common Q,W,massa,massac,massac2,massas
      common MBa,WBa,MSr,WSr,MY,WY,MEu,WEu
      common MLa,WLa,MZr,WZr,WRb,MRb

ccc     Modification for implement output by Ivan
cc      real ivan(36,20,300),ti(300)
      real vdens
cc      integer ri,  fi
cc      character*3 namei(300)
      real winds(31)

      do ii=1,31               
         winds(ii)=1.
      enddo
      winds(9)=1.

      
      pi=3.141592654
      
      ss=0          !c Index for SNIa bins
      Norm=0.0      !c Normalization of the IMF (real, supposed to be 1, slightly lessc)
      snian=0.0     !c Number of SNIa in the mass steps
      hottime=0     ! Possible to insert a delay between the death of the stars and
c                    the chemical enrichment of the ISM

      
      IMF=1

      tautype=1

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      write (*,*)'Which yields for low itermediate mass stars?'
c      write (*,*) 'Van den Hoek & Groenewegen (1) or Karakas (2)'
c      read(*,*) lowmassive

      lowmassive=1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      write (*,*)'Which yields (He, C & O) for massive high metallicy stars?'
!      write (*,*) '(0) WW95 --  (1) MM02 -- (2) Maeder 1992'
!      read(*,*) MM
      MM=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      



      call second(temp1)         
      
      if (lowmassive.eq.1) then
         ninputyield=32
         inputyield(1)='WW95aMD.csv' ! VdH&G Z=0.0001 WW95 0.0            .Mg.original
         inputyield(2)='WW95bMD.csv' ! VdH&G Z=0.004  WW95 0.02 10^-4     .Mg.original
         inputyield(3)='WW95cMD.csv' ! VdH&G Z=0.008  WW95 0.02 10^-2     .Mg.original
         inputyield(4)='WW95dMD.csv' ! VdH&G Z=0.02   WW95 0.02 10^-1     .Mg.original
         inputyield(5)=  'WW95e.csv' ! VdH&G Z=0.04   WW95 0.02 10^0      .Mg.original
         if (mm.eq.1) inputyield(5)='WW95eMM.csv'  ! ! VdH&G Z=0.04   WW95 0.02 10^0 + Meynet&M02 per C&O 
         if (mm.eq.2) inputyield(5)='WW95eMD.csv'  ! ! VdH&G Z=0.04   WW95 0.02 10^0 + Maeder'92  per C&O
      else
         ninputyield=35
         inputyield(1)='WW95zeroMD.K0001.csv'! Karakas Z=0.0001 WW95 0.0                           
         inputyield(2)='WW95-4MD.K004.csv'   ! Karakas Z=0.004  WW95 0.02 10^-4                    
         inputyield(3)='WW95-2MD.K008.csv'   ! Karakas Z=0.008  WW95 0.02 10^-2                    
         inputyield(4)='WW95-1MD.K02.csv'    ! Karakas Z=0.02   WW95 0.02 10^-1                    
         inputyield(5)='WW9502.K02.csv'      ! Karakas Z=0.04   WW95 0.02 10^0                     


         if (mm.eq.1) inputyield(5)='WW9502MM.K02.csv'! Karakas Z=0.04   WW95 0.02 10^0 + Meynet&M02 per C&O
         if (mm.eq.2) inputyield(5)='WW9502MD.K02.csv'! Karakas Z=0.04   WW95 0.02 10^0 + Maeder'92  per C&O
      endif


!     Read the yields for all the considered elements but Eu, Ba, La, Sr, Y, Zr.
!     Note that the "z" can indicate the dependence to metallicities but also 
!     different type of yields    
!     
!


 101  format (24(e12.5))
 102  format (5(1x,e12.5))
 104  format (6e13.5)

      do z=1,5

         OPEN(17,file='./DATI/'//inputyield(z),STATUS='old')
         read (17,*)
!         write (*,*) ' massa        He4          C12         O16
!     $   remn'
         do n=1,Ninputyield
            read (17,101) massa(n),W(1,n,z),W(2,n,z),W(3,n,z),W(4,n,z),
     $           W(5,n,z),W(6,n,z), W(7,n,z),W(8,n,z),W(9,n,z),
     $           W(10,n,z),W(11,n,z),W(12,n,z),W(13,n,z),W(14,n,z), 
     $           W(15,n,z),W(16,n,z),W(17,n,z),W(18,n,z),W(19,n,z),
     $           W(20,n,z),W(21,n,z), W(22,n,z),W(23,n,z)


         enddo
         close(17)
      enddo

      OPEN(17,file='./DATI/Kobayashi-Iron.dat',STATUS='old')
      read (17,*)
      do n=15,32
         read (17,104) massa(n),W(9,n,5),W(9,n,4),W(9,n,3),W(9,n,2),
     $           W(9,n,1)
      enddo
      close(17)

      
      OPEN(17,file='./DATI/Kobayashi-IronHyper.dat',STATUS='old')
      read (17,*)
      do n=15,32
         read (17,104) massa(n),WH(9,n,5),WH(9,n,4),WH(9,n,3),
     $        WH(9,n,2), WH(9,n,1)
         do i=1,5
            W(9,n,i)=(WH(9,n,i)+W(9,n,i))*0.5
         enddo
      enddo
      close(17)

      

!! Yields used by cristina chiappini basically  for CNO !
!! z=0  (which is equal to the yields for z= 10^-8)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 100  format (15(e12.5))

      z=6 
      OPEN(17,file='./DATI/'//'Cris0.dat',STATUS='old')
      read(17,*)
      read(17,*)
      do n=1,13
         read (17,100) massac(n),W(1,n,z),W(2,n,z),W(3,n,z),W(4,n,z),
     $        W(5,n,z),W(6,n,z), W(7,n,z),W(8,n,z),W(9,n,z),W(10,n,z),
     $        W(11,n,z),W(12,n,z),W(13,n,z),W(14,n,z)
      enddo
      close(17)


      z=7 !z= 10^-8 
      OPEN(17,file='./DATI/'//'Cris8.dat',STATUS='old')
      read(17,*)
      read(17,*)
      do n=1,13
         read (17,100) massac(n),W(1,n,z),W(2,n,z),W(3,n,z),W(4,n,z),
     $        W(5,n,z),W(6,n,z), W(7,n,z),W(8,n,z),W(9,n,z),W(10,n,z),
     $        W(11,n,z),W(12,n,z),W(13,n,z),W(14,n,z)
      enddo
      close(17)


      z=8 !z= 10^-5 
      OPEN(17,file='./DATI/'//'Cris5.dat',STATUS='old')
      read(17,*)
      read(17,*)
      do n=1,13
         read (17,100) massac(n),W(1,n,z),W(2,n,z),W(3,n,z),W(4,n,z),
     $        W(5,n,z),W(6,n,z), W(7,n,z),W(8,n,z),W(9,n,z),W(10,n,z),
     $        W(11,n,z),W(12,n,z),W(13,n,z),W(14,n,z)
      enddo
      close(17)

      z=9 !z= 0.004
      OPEN(17,file='./DATI/'//'Cris004.dat',STATUS='old')
      read(17,*)
      read(17,*)
      do n=1,13
         read (17,100) massac(n),W(1,n,z),W(2,n,z),W(3,n,z),W(4,n,z),
     $        W(5,n,z),W(6,n,z), W(7,n,z),W(8,n,z),W(9,n,z),W(10,n,z),
     $        W(11,n,z),W(12,n,z),W(13,n,z),W(14,n,z)
      enddo
      close(17)

      z=10 !z= 0.02
      OPEN(17,file='./DATI/'//'Cris02.dat',STATUS='old')
      read(17,*)
      read(17,*)
      do n=1,13
         read (17,100) massac(n),W(1,n,z),W(2,n,z),W(3,n,z),W(4,n,z),
     $        W(5,n,z),W(6,n,z), W(7,n,z),W(8,n,z),W(9,n,z),W(10,n,z),
     $        W(11,n,z),W(12,n,z),W(13,n,z),W(14,n,z)
      enddo
      close(17)

 103  format (4(e12.5))
      OPEN(17,file='./DATI/Bariumnew.dat',STATUS='old')
      do n=1,4
         read (17,103) Mba(n),WBa(n,1),WBa(n,2),WBa(n,3)
         !write(*,103) MBa(n),WBa(n,1),WBa(n,2),WBa(n,3)
         WBa(n,3)=WBa(n,3)*640.!347.!
         WBa(n,2)=WBa(n,2)!*15.  !640.!347.!
      enddo
      close(17)
      OPEN(17,file='./DATI/Strontiumnew.dat',STATUS='old')
      do n=1,4
         read (17,103) MSr(n),WSr(n,1),WSr(n,2),WSr(n,3)
         !write(*,103) MSr(n),WSr(n,1),WSr(n,2),WSr(n,3)
         WSr(n,3)=WSr(n,3)*50.!66.
         WSr(n,2)=WSr(n,2)!*5.
      enddo
      close(17)
      OPEN(17,file='./DATI/Yttriumnew.dat',STATUS='old')
      do n=1,4
         read (17,103) MY(n),WY(n,1),WY(n,2),WY(n,3)
         
         WY(n,3)=WY(n,3)*150. !159.!150.
         WY(n,2)=WY(n,2)!*5.
         !write(*,103) MY(n),WY(n,1),WY(n,2),WY(n,3)
      enddo
      
      close(17)
      OPEN(17,file='./DATI/Lantanumnew.dat',STATUS='old')
      do n=1,4
         read (17,103) MLa(n),WLa(n,1),WLa(n,2),WLa(n,3)
      
         WLa(n,3)=WLa(n,3)*211.
         WLa(n,2)=WLa(n,2)!*15.
         !write(*,103) MLa(n),WLa(n,1),WLa(n,2),WLa(n,3)
      enddo
      
      close(17)
      OPEN(17,file='./DATI/Zirconiumnew.dat',STATUS='old')
      do n=1,4
         read (17,103) MZr(n),WZr(n,1),WZr(n,2),WZr(n,3)
         !write(*,103) MZr(n),WZr(n,1),WZr(n,2),WZr(n,3)
         WZr(n,3)=WZr(n,3)*346.
         WZr(n,2)=WZr(n,2)!*5.
         !write(*,103) MZr(n),WZr(n,1),WZr(n,2),WZr(n,3)
      enddo
      
      close(17)
      OPEN(17,file='./DATI/Rubidiumnew.dat',STATUS='old')
      do n=1,4
         read (17,103) MRb(n),WRb(n,1),WRb(n,2),WRb(n,3)
         !write(*,103) MRb(n),WRb(n,1),WRb(n,2),WRb(n,3)
         WRb(n,3)=WRb(n,3)*10.
         WRb(n,2)=WRb(n,2)!*5.
      enddo
      close(17)
      OPEN(17,file='./DATI/Europiumnew.dat',STATUS='old')
      do n=1,4
         read (17,*) MEu(n),WEu(n,1),WEu(n,2),WEu(n,3)

         WEu(n,3)=WEu(n,3)*1.e-20!/(-0.2)*10.
         WEu(n,2)=WEu(n,3)*1.e-20
         WEu(n,1)=WEu(n,3)*1.e-20
         !write(*,*) MEu(n),WEu(n,1),WEu(n,2),WEu(n,3)
      enddo
      close(17)


      
     

     

      call leggili
      call leggila
      call leggirb
      call leggizr
      call leggieu
      call leggiy
      call leggisr
      call leggi
      




      OPEN(17,FILE='./DATI/amu.dat',STATUS='old')
      READ(17,20) (AMU(K),K=1,115)

 20   FORMAT(5E12.5)
      close(17)
      
      if (imf.eq.1) call IMFScalo
      if (imf.eq.2) call IMFKroupa
      if (imf.eq.3) call IMFSalpeter
     

      binmass=115               ! number of grid points in the stellar masses
      ss2=binmass-1
      ss=0
      do jj=1,binmass-1         ! do in the grid of the masses
!!!!!
               
         mstars1(jj)=amu(binmass-jj) ! just change the order of the 
         mstars1(jj+1)=amu(binmass-jj+1) ! grid of the masses
!!!
!!!    As a Function of the IMF we define the weight of an interval in the grid of the masses
!!!    THE WEIGHT IS THE NUMBER OF STAR in that BIN for MSUN
!!!    so multiplied for the SFR gives the number of stars in such bin


         if (IMF.eq.1) 
     $        call MultiSCALO(mstars1(jj),Mstars1(jj+1),multi1(jj))
         
         if (IMF.eq.2) 
     $        call MultiKroupa(mstars1(jj),Mstars1(jj+1),multi1(jj))

         if (IMF.eq.3) 
     $        call MultiSALPETER(mstars1(jj),Mstars1(jj+1),multi1(jj))


         mstars(jj)=(mstars1(jj)+mstars1(jj+1))/2. !!! The stellar mass is defined as the average between the 2 

         binmax(jj)=0.0         !!! in general for a normal star there is not a mass for
!                                                        !!! the binary star!     
            
         
         If ((mstars(jj).ge.3.).and.(mstars(jj).le.16.)) then ! in this mass range we allow to form SNIa


            mmu=mstars(jj)      ! Total mass of the binary system
            mun=10              ! number of mass bin in which we split the mu
            
            mumin=0.8/mmu       ! Mu for minimum mass of the secondary star of the binary systems
            mumin2=1.-8./mmu    ! Constrain that the primary should be less than 8Msun
            


            if (mumin2.gt.mumin) mumin=mumin2 ! Choose the constrain
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   
!!!   We split (up to now) in the same way it has be done for the standard
!!!   model, not the best not the smartest 
!!!   
!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            mux=0.5-MUMIN
            XMU(1)=MUMIN
            XMU(2)=MUMIN+0.01*MUX
            XMU(3)=MUMIN+0.02*MUX
            XMU(4)=MUMIN+0.05*MUX
            XMU(5)=MUMIN+0.1*MUX
            XMU(6)=MUMIN+0.2*MUX
            XMU(7)=MUMIN+0.3*MUX
            XMU(8)=MUMIN+0.4*MUX
            XMU(9)=MUMIN+0.6*MUX
            XMU(10)=MUMIN+0.8*MUX
            XMU(11)=0.5
                  
            do j3=1,10          !!! do in the 10 bins of the grid for SNIa 
               ss=ss+1          !!! total number of iteration for calculation of SNIa rate
               ss2=ss+binmass-1 !!! total number of possible weights for the 
                                      !!! mass bin of single star and binary progenitor of SNIa
                     
                        
               mu1=XMU(j3)      ! Lower limit for the mass of the secondary star(/total mass) we integrate
               mu2=XMU(j3+1)    ! upper limit


               mstars(ss2)=mstars(jj)*
     $              (xmu(j3)+xmu(j3))/2 ! mass of the secondary star (average beteween the 2 limits)
                                                   ! for which we calculate the timescale for explotion

               multi1(ss2)=(8.*(mu2**3.-mu1**3.))*
     $              multi1(jj)*0.09 ! weight of this configuration!

               BINmax(ss2)= mstars(jj)-
     $              mstars(ss2) ! Mass of the primary star of the system
                        
               Norm=Norm+multi1(ss2)*mstars(JJ) ! normalization of the IMF for this weight

            enddo

            multi1(jj)=0.95*multi1(jj) ! the weight of the single star in the SNIa range is decreased because
                                          ! 0.05 is going in SNIa

         endif

         !If ((mstars(jj).ge.0.8).and.(mstars(jj).le.8.)) then ! addiction of novae: assumed a new bin with binmax negative
         If ((mstars(jj).ge.2.0).and.(mstars(jj).le.8.)) then ! addiction of novae: assumed a new bin with binmax negative
         !If ((mstars(jj).ge.4.).and.(mstars(jj).le.8.)) then ! addiction of novae: assumed a new bin with binmax negative

            ss=ss+1                
            ss2=ss+binmass-1
            mstars(ss2)=mstars(jj)
            BINmax(ss2)=-1
!            multi1(ss2)=0.015*multi1(jj)
            multi1(ss2)=0.015*multi1(jj)
         
         
         endif
! modification for for NSM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         If ((mstars(jj).ge.90.0).and.(mstars(jj).le.30.)) then 

            ss=ss+1                
            ss2=ss+binmass-1
            mstars(ss2)=mstars(jj)
            BINmax(ss2)=-mstars(JJ)
            multi1(ss2)=0.05*multi1(jj)
         
         
         endif
! modification for for NSM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



         
         Norm=Norm+multi1(jj)*mstars(JJ) ! total weigth of the sum of each weight, should be 1, instead is a bit less
                                               ! because of how it is normalized the Scalo-Zita IMF

               

      enddo

      do j3=1,ss2      
!   write(*,*) mstars(j3),binmax(j3),multi1(j3)
         mstars1(j3)=0. 
      enddo

!!!!!           
!  This part of the code is devoted to order as a function of lifetime the different mstar ...
!!!!
      do j3=1,ss2
         
         do ss=1,j3
            
            if (tau(mstars1(ss),tautype,binmax1(ss)).ge.    
     $           tau(mstars(j3),tautype,binmax(j3))) then
               
               do k=j3-1,ss,-1
                  
                  mstars1(k+1)=mstars1(k)
                  binmax1(k+1)=binmax1(k)
                  multi2(k+1)=multi2(k)
                  
               enddo

               mstars1(ss)=mstars(j3)
               binmax1(ss)=binmax(j3)
               multi2(ss)=multi1(j3)
               
               exit 
                      
            endif               
         enddo
              
      enddo


!!!!
!!!!!!!!!!!!!! IN ORDER e quick re-definition of the variable |
!!!!             this part of the code run once at the begin
!!!!             not a problem to do it ...

      do jj=1,ss2

         tdead(jj)=tau(mstars1(jj),tautype,binmax1(jj))
                        
!     write(*,*) mstars1(jj),tdead(jj),binmax1(jj),multi2(jj)
               
         mstars(jj)=Mstars1(jj)
         binmax(jj)=binmax1(jj)
         multi1(jj)=Multi2(jj)

      enddo



      !do ence=1,10000
      ence=1

      superf=20000.
      p=1                       ! Just one run for an homogenous model
      threshold=0.1
      step=1                    ! Timestep in Myr usually set in 1Myr 
      sigmasun=50.
      kappa=1.5
      Rm=8.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Parameters 
!!!!!!!
!!!!!!    Model ENCE
!!!!!!

!      CALL GET_COMMAND_ARGUMENT(1,num1char)   !first, read in the two values 
!      CALL GET_COMMAND_ARGUMENT(2,num2char)   !first, read in the two values 
!      CALL GET_COMMAND_ARGUMENT(3,num3char)   !first, read in the two values 
!      CALL GET_COMMAND_ARGUMENT(4,num4char)   !first, read in the two values 
!      CALL GET_COMMAND_ARGUMENT(5,num5char)   !first, read in the two values 
!      CALL GET_COMMAND_ARGUMENT(6,num6char)   !first, read in the two values 
!      CALL GET_COMMAND_ARGUMENT(7,num7char)   !first, read in the two values 
!
      

   

      ence=endoftime
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!  cleaning  variables
      jj=1
      do i=1,elem
         Ini(i)=0.
         do j=1,NMAX
            QQN(i,j,jj)=1.e-20  !sigmar/100
            
         enddo
      enddo
      Ini(31)=1.e-9
    

!!!   output files - open

      OPEN(17,file=
     $     './RISULTATI2/fis.encesmin.dat', STATUS='unknown')

      write(17,124) 'time', 'all','gas','stars','remn','hot','zeta',
     $     'SFR','nume','SFR2','SIaN','SIar'

      OPEN(13,file=
     $     './RISULTATI2/modencesmin.dat', STATUS='unknown')

      !write(13,*) endoftime,sigmat,sigmah,time_wind,pSFR,pwind,delay
     
      write(13,124) 'time','all','gas ','star','SFR ','run ','Hen ',
     $      'C12 ','O16 ','N14 ','C13 ','Ne  ','Mg  ','Si  ','Fe  ',
     $      'S14 ','C13S','S32 ','Ca  ','Remn','Zn  ','K   ','Sc  ',
     $      'Ti  ','V   ','Cr  ','Mn  ','Co  ','Ni  ','La  ','Ba  ',
     $      'Eu  ','Sr  ','Y   ','Zr  ','Rb   ','Li  ','H   ','He4 '

 124  format(39(3X,A4,6X))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!   set to zero the variables each run (which corresponds to a volume completely evolved)

      windist=0.
      do tt=1,NMAX
         hot(tt)=0.
         stars(tt)=0.
         remn(tt)=0.
         gas(tt)=0.
         all(tt)=0.
         zeta(tt)=0.
         wind(tt)=0.
         SNIAnum(tt)=0.
      enddo
           
      t=0 !!! SO WE START WITH 0


      do !!!! START OF THE LOOP IN THE TIME !!!


         step=1
         t=t+step
         
!!!   
!!!   Infall rate (note at each time we compute the total mass and the gas mass)
!!!   
!!!   

         t3=t
         step3= 1

            

         do  !!! UPDATE FOR ALL THE TIME THE AMOUNT OF TOTAL MASS AND GAS




!!!!    In this model just the infall for the thin disc 
!!!!

            

!            all(t3)=(sigmah)*(1./(1.-exp(-endoftime/delay)))*
!     $           (1.-1./exp(real(t3)/delay))

            all(t3)=all(t3-1)+sigmah*superf/(2.5*sigmat) 
     $           *exp(-(real(t3)-delay)**2/(2*sigmat**2))




          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   gas is the mass which is not in remnants or in stars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            gas(t3)=all(t3)-stars(t3)-remn(t3)
     $           -hot(t3)-wind(t3)


            if (t3.ge.endoftime) exit
            
            t3=t3+step3
            
         enddo

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         
         if (gas(t)/superf.gt.threshold) then
      
            SFR= psfr*(gas(t)/(superf*sigmah))**kappa !!! General equation for SF law (see portinari and chiosi 1999)
     $        *    (sigmah/sigmasun)**(kappa-1)      !!! General equation 
     $        *    (8./Rm)                           !!! Scaling law for different radii
     $        *        (superf/1000)*(sigmah) !!!! conversion from psi to SF law in mass in the considered volume in Myr


            
         else
            SFR=0.
         endif


!!!
!!!!
!!!!!!!!!!!!!        
        ! if (t.gt.1000) threshold=0.
         
!!!   
!!!   
!!!   Output on the screeen
!!!   
!         write(*,*) ence,psfr,pwind,real(t), SFR,  all(t+1) ,gas(t), 
!     $        stars(t),wind(t),zeta(t)
         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!     Threshold check (if any threshold)
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         
 ! if ((SFR.gt.threshold).or.(zeta(t).gt.1e-20)) then !!! if the "SFR" value  exceeds the 100Msun
         if ((gas(t)/superf).ge.threshold) then 

            DO JJ=1,ss2
               
               if ((tau(mstars(jj),tautype,binmax(jj))+real(t))
     $              .gt.13500) then
                  oldstars(t)=oldstars(t)+multi1(jj)*SFR
               endif

               call interp(mstars(jj),zeta(t),BINMAX(jj),
     $              Hecore)     ! FOR EACH BIN (SINGLE AND SNIA) we calculate the yields


               if (Binmax(jj).gt.0) then   ! if it is a SNIa we calculate the mass of the total stellar system
                  Mstars1(jj)=BINMAX(JJ)+mstars(jj)
               else
                  Mstars1(jj)=mstars(jj)
               endif

               hecores(jj)=hecore          ! we define as a function of JJ the heliumcore
    
               do i=1,elem-2

                  QIspecial(i,jj)=q(i)    ! and the yields for each element!!
                  q(i)=0.0                
               enddo

            ENDDO

                 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ! aggiorna tutte le variabili a causa delle nuove stelle a partire da 
!     !                      dal tempo attuale 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  t3=t            


                  starstot=SFR*Norm
                  difftot =SFR*Norm
                  remntot =0.0
                  hecoretot=0.0
                  snian=0.
                  jj=1

                  do    ! CYCLE IN THE TIME FROM T TO TENDOFTIME

                     step3=1
                   

                     do ii=1,(elem-2) !!!! ciclo sugli elementi  che aggiorna
                        
                           if (real(t3).ge.(real(t)+tdead(jj))) then

                              q(ii)=q(ii)+QIspecial(ii,jj)*   !!! We calculate the total mass for each element
     $                             multi1(jj)*SFR             !!! ejected by the mass bins when they enrich the ISM 
                              
                              if (ii.eq.1) then !             !!! JUST FOR Hecore - stars and diff we calculate once.
                          
                                hecoretot=hecoretot+hecores(jj)
     $                                *multi1(jj)*SFR
                                difftot=difftot-(mstars1(jj)-
     $                               hecores(jj)) *multi1(jj)*SFR
                                starstot=starstot-(mstars1(jj)* 
     $                               multi1(jj)*SFR)

                                if (binmax(jj).gt.0.0) then    !!! IN BINMAX greater than 0.0 then we have a SNIA!
                                   SNIAn=SNIAn+
     $                                 multi1(jj)*SFR     

                                endif
                                
                              endif                         
                              
                              
                              
                           endif
                           if (real(t3).ge.(real(t)+tdead(jj+1)))  then
                              
                              
                       
                              else
!!!   
!!!   Aggiunge lo yield dell'elemento i-mo da quando la stella scoppia in poi    
!!!   eventualmente implemententando una fase di gas caldo in cui non si mischia
!!!   

                                 if (II.ne.14) then

                                    
                                    QQN(ii,t3,p)=QQN(ii,t3,p)+Q(ii)
     $                                   -QQN(ii,t,p)*(difftot)/gas(t)


                                    if (II.eq.31) then ! astration per il LITIO!
                                       QQN(ii,t3,p)= QQN(ii,t3,p)
     $                                      -QQN(ii,t,p)*
     $                                      (mstars1(jj)+hecores(jj))* 
     $                                      multi1(jj)*SFR/gas(t)
                                    endif
!!!

                                 
!!!   sottrae la massa dell'elemento i-mo lockato nella stella (ma che verr espulso)
!!!   pari alla frazione di massa dell'elemento a t presente nel gas per
!!!   la massa della stella meno il remnant o piu' precisamente l'Helium core
!!!   in quanto bisogna considerare che parte viene sintetizzata e quindi trasformata in newly produced
                                 
                                    
!                                    QQN(ii,t3,p)=QQN(ii,t3,p)
!     $                                   
                                    
!!!   sottrae la massa dell'elemento i-mo lockato (o trasformato in qualcos'altro ...) per sempre nel remnant
!!!   pari alla frazione di massa dell'elemento a t presente nel gas per
!!!   la massa del remnant o meglio Helium core
                              
!                                    QQN(ii,t3,p)=QQN(ii,t3,p)
!     $                                   -QQN(ii,t,p)*Hecoretot/gas(t)




                                 endif !!! 14

                              endif


                              
                              
                              enddo ! ciclo per i vari elementi gestiti


                     

                     if (real(t3).ge.(real(t)+tdead(jj)))  then



                        if (real(t3).ge.(real(t)+tdead(jj+1)))  then
                           
                           jj=jj+1  !!! IF ALSO THE NEXT MASS BIN IS IN GONNA DIE IN THIS TIME STEP
                                    !!! IT JUST GOES TO THIS NEXT WITH UPDATE THE VARIABLES
 
                        else                                 ! IF NOT UPDATE THE VARIABLES
                        
                           stars(t3)=stars(t3)+starstot
                      
                           remn(t3)=remn(t3) +  q(14)


                           if (q(14).lt.0) exit                           

                           snianum(t3)=snianum(t3)+snian


                           if (t3.ge.endoftime) exit         ! CHECK IF IT HAS TO GO OUT
                           
                           t3=t3+step3                       ! UPDATE THE TIME
                          
                           
   
                           jj=jj+1                           ! GO TO THE NEXT MASS BIN
                        
                        endif
                     else                                    ! IF IT HAS NOT REACH THE NEXT MASS BIN

                           stars(t3)=stars(t3)+starstot      ! UPDATE THE VARIABLES
                      
                           remn(t3)=remn(t3) + q(14)
                           
                           if (q(14).lt.0) exit
                           
                           snianum(t3)=snianum(t3)+snian

                           
                           if (t3.ge.endoftime) exit         ! CHECK IF WE HAVE REACHED THE ENDOFTIME
                           
                           t3=t3+step3                       ! UPDATE THE TIME
                        





                  endif
                 
                 
                  enddo     ! CYCLE IN THE TIME FROM T TO TENDOFTIME



               


         endif                  !!! controllo legato alla threshold, se eventualmente presente
            


         t3=t                   !!! altro giro nel tempo per aggiornare altre variabili.
                                !!! e per tenere in conto del vento eventualmente presente
                                !!! POSSO METTERLO NEL CICLO PRECEDENTE MA NON ACCRESCE LA VELOCITA'!
 

         
         
         if ((t.gt.time_wind).and.(gas(t).gt.0)) then
            windist=pwind*SFR
         else
            windist=0.
         endif
         do 
            step3=1
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Sottrae dal gas i vari elementi a causa del vento se presente!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            wind(t3)=wind(t3)+windist

               
            do ii=1,31
               if ((ii.ne.14).and.(SFR.gt.0.)) then
                  QQN(ii,t3,p)=QQN(ii,t3,p) -
     $                 QQN(ii,t,p)*winds(ii)*windist/gas(t)


                  
                  if (QQN(ii,t3,p).le.1.e-20)  QQN(ii,t3,p)=1.e-20
               endif
            enddo

         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Calcola l'arricchimento 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            zeta(t3)=0.
            
            do i=2,(elem-2)          !!! da 2 poich elio non  metallo!!! (1  elio)
               

               if (i.ne.14) zeta(t3)=zeta(t3)+QQN(i,t3,p)
            

            enddo

            if (SFR.gt.0) then
           
                           
               zeta(t3)= zeta(t3)/(gas(t3))
            else
               zeta(t3)=1.e-20
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  QUESTA E' una ASSUNZIONE CHE PRENDO PER CALCOLARE
!!!!  L'ELIO E L'IDROGENO AD OGNI TIME STEP
!!!!  OVVERO CHE L'ELIO SIA UGUALE A QUELLO PRIMORDIALE + piu' quello che
!!!!  La nucleosintesi crea
!!!!  mentre l'idrogeno e' pari a quello primordiale meno quello trasformato
!!!!  in metalli o in elio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ! calcolo dell'elio totale ...
            QQN(elem,t3,p)=(gas(t3))*(0.241)+QQN(1,t3,p) 
!     ! calcolo dell'idrogeno rimasto         
            QQN(elem-1,t3,p)=(gas(t3))*(0.759-zeta(t3))-QQN(1,t3,p)
            
            if (t.gt.1) then
            
               QQN(31,t3,p)=QQN(31,t3,p)+(all(t)-all(t-1))*Ini(31) !infalled from primordial composition
           
!!!!  Spallation effect
               spalla(t3)=10**(
     $              -9.50+1.24*((alog10(QQN(9,t3,p)/QQN(elem-1,t3,p)))-
     $              (-2.75))+alog10(QQN(elem-1,t3,p)))
!write(*,*) spalla
               QQN(31,t3,p)=QQN(31,t3,p)+spalla(t3)-spalla(t3-1)
            else
               QQN(31,t3,p)=QQN(31,t3,p)+all(t)*Ini(31) !infalled from primordial composition
            endif
            
            !!!   
            
            if (t3.ge.endoftime) exit
            
            t3=t3+step3            
            
       
         enddo  

!!!   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   Scrive le abbondanze ad ogni time step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   
         write(13,140) real(t),all(t),gas(t),
     $        stars(t),SFR/real(step),oldstars(t),
     $        QQN(1,t,p), QQN(2,t,p), QQN(3,t,p), QQN(4,t,p),
     $        QQN(5,t,p), QQN(6,t,p), QQN(7,t,p), QQN(8,t,p),
     $        QQN(9,t,p), QQN(10,t,p), QQN(11,t,p), QQN(12,t,p),
     $        QQN(13,t,p), QQN(14,t,p), QQN(15,t,p), QQN(16,t,p),
     $        QQN(17,t,p), QQN(18,t,p), QQN(19,t,p), QQN(20,t,p),
     $        QQN(21,t,p), QQN(22,t,p), QQN(23,t,p), QQN(24,t,p),
     $        QQN(25,t,p), QQN(26,t,p), QQN(27,t,p), QQN(28,t,p),
     $        QQN(29,t,p), QQN(30,t,p), QQN(31,t,p), QQN(32,t,p),
     $        QQN(33,t,p)





         if (t.ge.endoftime) exit

      enddo                     !!! do per riprodurre nel tempo 
      




      t=1

      do 

         step=1
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   Scrive le caratteristiche "fisiche"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         write(17,122) real(t),all(t),gas(t),stars(t),remn(t),hot(t),
     $       zeta(t),SFR,real(p),20*(stars(t)-stars(t-step))/real(step),
     $        SNIanum(t),     
     $        (SNIAnum(t) -SNIAnum(t-step))/real(step)


         

         if (t.ge.endoftime) exit
         t=t+step
      enddo

      


      close(13)
      close(17)
      
 122  format (12(1x,e12.5))
 140  format (39(1x,e12.5))
 !     enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call second(temp2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     !
!     ! tempo che sta il modello a girare
!     !
!      write(*,*) temp2-temp1
!     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      STOP 
      END
      subroutine IMFScalo

      real UM1,UM2
      real MI,M1,M3,MS
      real A,B
      real zita
      real Azita,Bzita,CZITA,DZITA
      COMMON UM1,UM2,UM3,UM4
      COMMON A,B,C,D
      COMMON M1,M2,M3,M4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      SCALO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    
      UM1=-1.35                 !! esponente stelle piccole  
      UM2=-1.7                  !! esponente stelle massicce 
      MI=0.1                    !! Massa minima              
      M1=2.0                    !! Massa cambio pendenza     
      MS=100.0                  !! Massa massima                                          
      
      m3=1.                     !! SERVE PER IL CALCOLO DELLA NORMALIZZAZIONE "ZITA"

      ZITA=0.3
      Azita=MS**(1.+UM2)
      Bzita=M1**(1.+UM2)
      Czita=M1**(1.+UM1)
      Dzita=M3**(1.+UM1)        !mi senza zita o m3 se vuoi usare zita a 0.3!!!


      B=ZITA/(M1**(UM2-UM1)*(Czita-Dzita)/(1.+UM1)+
     $     (Azita-Bzita)/(1.+UM2))

      A=B*M1**(UM2-UM1)

      return
      end


      subroutine IMFSalpeter

      real UM1
      real MI, MS
      real A
      

      COMMON UM1,UM2,UM3,UM4
      COMMON A,B,C,D
      COMMON M1,M2,M3,M4
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     !   SALPETER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   
!!!   
      
      UM1=-1.35                 !! esponente salpeter                       
      MI=0.1                    !! Massa minima                  
      MS=80.0                   !! Massa massima                 
      
      
      Azita=MS**(1.+UM1)                                   
      Bzita=MI**(1.+UM1)                                   
      
      
      A=(1.+UM1)/(Azita-Bzita)                           
      write (*,*) A
 
      
      return 
      end



      subroutine MULTIKROUPA(Max,Min,Result)

      real Max,Min,Result
      real UM1,UM2,UM3,UM4
      real A,B,C,D
      real M1,M2,M3

      common UM1,UM2,UM3,UM4
      common A,B,C,D
      COMMON M1,M2,M3,M4


      if (MAX.le.M1) then 
         result=(MIN)**(UM1)-
     $        (MAX)**(UM1)
         result=result/(UM1)*A 
                                !! Frazione di stelle in numero in quell'intervallo di masse 
                                !  per unit di massa
                                !     
      else
         if (MAX.le.M2) then
                           
            result=(MIN)**(UM2)-
     $           (MAX)**(UM2)
            result=result/(UM2)*B 
                                !     ! Frazione di stelle in numero in quell'intervallo di masse 
                                !     per unit di massa
         else
            if (MAX.le.M3) then
               
               result=(MIN)**(UM3)-
     $              (MAX)**(UM3)
               result=result/(UM3)*C
            else
               result=(MIN)**(UM4)-
     $              (MAX)**(UM4)
               result=result/(UM4)*D
               
            endif
            
         endif
      endif
      if (MAX.gt.50) result=0.
      
      return
      end



      subroutine MULTISCALO(Max,Min,Result)

      real Max,Min,Result
      real UM1,UM2
      real A,B
      real M1

      common UM1,UM2,UM3,UM4
      common A,B,C,D
      COMMON M1,M2,M3,M4




      if (MAX.le.M1) then 
         
         result=(MIN)**(UM1)-
     $        (MAX)**(UM1)
         result=result/(UM1)*A 
!     !  Frazione di stelle in numero in quell'intervallo di masse 
!     !  per unit di massa
!     
      else
         result=(MIN)**(UM2)-
     $        (MAX)**(UM2)
         result=result/(UM2)*B 
!     ! Frazione di stelle in numero in quell'intervallo di masse 
!     ! per unit di massa


      endif

      return
      end
      
      subroutine MULTISALPETER(Max,Min,Result)

      real Max,Min,Result
      real UM1
      real A

      common UM1,UM2,UM3,UM4
      common A,B,C,D
      COMMON M1,M2,M3,M4

      result=(MIN)**(UM1)- 
     $     (MAX)**(UM1)     
      result=result/(UM1)*A        
      
      return 
      end


       end module main_mod
