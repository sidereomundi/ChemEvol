!f2py fixed
       module interpolation_mod
       contains
      subroutine interp(H2,zeta,BINMAX,Hecore)

      IMPLICIT none
      
      INTEGER NMAX,elem
      PARAMETER (NMAX=23)       ! number of elements (excluded s-r process)
      parameter (elem=33)
      real DD(2),mm(2),met(2)!,mbar(2),bar(2)
      real Q1(elem),Q2(elem)
      
      real cosn2(nmax)
      integer j,k,kk,i,z,zz
      real W(nmax,35,15),massa(35),massac(13),massac2(4),massas(7)
      real MBa(4),WBa(4,3)
      real MSr(4),WSr(4,3)
      real MY(4),WY(4,3)
      real MLa(4),WLa(4,3)
      real MRb(4),WRb(4,3)
      real MZr(4),WZr(4,3)
      real MEu(4),WEu(4,3)

      real Rando
      integer koba
      
      real Q(elem),qbar,qeu2,qla,qsrr,qzr,qy,qrb,aa(15)
      real qli
      real QIa(elem)
      real dy
      real mlow1,mlow2,mup1,mup2
      real costBa,CostLa,CostEu,CostSr,CostZr,CostY,CostRb
      real value1,value2,value3,value4
      
      real H,H2,zeta,BINMAX,ratio,Hecore
      integer ninputyield

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


      H=H2

!!!      Yields di SNIa 

      QIa(1)=  0.              !      1)He4         
      QIa(2)=  0.048           !      2)C12         
      QIa(3)=  0.143           !      3)O16         
      QIa(4)=  1.16e-6         !      4)N14         
      QIa(5)=  1.40e-6         !      5)C13         
      QIa(6)=  0.00202         !      6)Ne          
      QIa(7)=  0.0425  *1.2        !      7)Mg          
      QIa(8)=  0.154           !      8)Si          
      QIa(9)=  0.6             !      9)Fe          
      QIa(10)= 0.              !     10)neutron rich
      QIa(11)= 0.              !     11)C13S        
      QIa(12)= 0.0846          !     12)S32         
      QIa(13)= 0.0119          !     13)Ca          
      QIa(14)= 0.              !     14)Remn        
      QIa(15)= 6.345e-4        !     15)Zn          
      QIa(16)= 1.19E-3         !     16)K           
      QIa(17)= 1.28e-5         !     17)Sc          
      QIa(18)= 4.16e-4         !     18)Ti          
      QIa(19)= 7.49E-5         !     19)V           
      QIa(20)= 5.67E-3         !     20)Cr          
      QIa(21)= 6.21E-3         !     21)Mn          
      QIa(22)= 1.56E-3         !     22)Co          
      QIa(23)= 1.83e-2         !     23)Ni          
      QIa(24)= 0.              !     24)La          
      QIa(25)= 0.              !     25)Ba          
      QIa(26)= 0.              !     26)Eu          
      QIa(27)= 0.              !     27)Sr          
      QIa(28)= 0.              !     28)Y          
      QIa(29)= 0.              !     29)Zr          
      QIa(30)= 0.               !    30)Rb          
      QIa(31)= 0.               !    31)Li          
            
         


      !!! Setting delle variabili 

      k=33 !!! Fuori dal range (max  23)
      z=1000 !!! inizializzato a valore assurdo...

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inizializzazione delle variabili  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=1,elem
         Q(i)=0.
      enddo


    
!-----------------------------------------
!-----------------------------------------
      !Pone a 1 le costanti
      
      DO I=1,NMAX      
         COSN2(I)=1.000e+00
      end do
!-----------------------------------------------
!-----------------------------------------------
      !!! Inserisce le costanti determinate da Francois et al.(2004)
      
      
      COSN2(3)=1.               ! Oxygen
      COSN2(7)=6.!cris! 7.               ! Magnesium
      cosn2(8)=1.               ! Silicon
      cosn2(13)=1.              ! Calcium
      cosn2(23)=1.              ! Nickel
      COSN2(16)=2.              ! Potassium
      COSN2(17)=1.15            ! Scandium
      COSN2(18)=15.             ! Titanium
      COSN2(20)=3.              ! Chromium
      cosn2(22)=0.3             ! Cobalt
      cosn2(19)=4.5             ! Vanadium
      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Valuta che ratio bisogna applicare per rimediare all'errore
!!   di partenza sulle SNIa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (binmax.ge.0) then
      
      if (BINMAX.gt.0.) then 

         ratio=(H+BINMAX)/BINMAX
        
        
        H=BINMAX
      else
         ratio=0.0
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!! Determina in che bins di massa bisogna interpolare
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      if (H.gt.8) then


!!!! Valori di metallicit per WW95

      aa(1)=0.02*0.    !  0.     
      aa(2)=0.02*1.e-4 !  0.004  
      aa(3)=0.02*1.e-2 !  0.008  
      aa(4)=0.02*1.e-1 !  0.02   
      aa(5)=0.02*1.    !  0.04   

      else

!!!! Valori di metallicit per van den Hoek and Groenewegen
      aa(1)= 0.     
      aa(2)= 0.004  
      aa(3)= 0.008  
      aa(4)= 0.02   
      aa(5)= 0.04   

      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Calcolo in che bin di metallicit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=1,5
         If (aa(j).le.zeta) Z=J
      enddo

      zz=z+1
   
      If (z.eq.5) ZZ=Z


      met(1)=aa(z)
      met(2)=aa(ZZ)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!  Prove per valutare i valori da 8 a 11 ...
!!!!
!!! 
!!!
!!!!      if ((H.gt.8).and.(H.lt.11.)) then
!!!!         k=15
!!!
!!!      if (z.lt.5) then ! allora interpola almeno le metallicit, tendendo fissa la massa a 11
!!!
!!!      do j=1,23
!!!         D(1)=W(j,k,z)
!!!         D(2)=W(j,k,zz)
!!!         call polint(mm,D,2,H,Q1(j),dy)
!!!         Q(j)=Q1(j)*cosn2(j)
!!!      enddo
!!!
!!!      else   ! allora prende solo la metallicit massima e massa 11
!!!         do j=1,23
!!!            Q(j)=W(j,k,z)*cosn2(j)
!!!         enddo
!!!
!!!      endif
!!!
!!!      else ! caso in cui interpolo sia massa che metallicit
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Calcolo in che bin di massa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=1,ninputyield
         
         IF (massa(j).lt.H) K=J 
         
      enddo

      KK=K+1
      

      mm(1)=massa(k)
      mm(2)=massa(kk)

    
      !!!!
      !!!!  VALORI SOLARI!!! 
      !!!!
      
!      z=5
      !!!!!!!!!!!!!

      !!! Interpola i valori corretti per le stelle massicce 
      !!! con i valori di WW95B modificati da Francois et al 2004

!!!!!!!!!!!! Interpola su griglia massa valore metallicit pi basso 

      do j=1,23
         DD(1)=W(j,k,z)
         DD(2)=W(j,kk,z)
         call polint(mm,DD,2,H,Q1(j),dy)
         Q1(j)=Q1(j)*cosn2(j)
      enddo

      if (z.lt.5) then

!!!!!!!!!!!! Interpola su griglia massa valore metallicit pi alto

      do j=1,23
         DD(1)=W(j,k,zz)
         DD(2)=W(j,kk,zz)
         call polint(mm,DD,2,H,Q2(j),dy)
         Q2(j)=Q2(j)*cosn2(j)
      enddo


!!!!!!!!!!!! Interpola su griglia metallicit i due valori !!!!

      do j=1,23

         DD(1)=Q1(j)
         DD(2)=Q2(j)
         call polint(met,DD,2,zeta,Q(j),dy)

      enddo

      else 

!!!! Caso in cui la metallicit sia maggiore di quella massima considerata la prende 
!!!   sempre costante e pari alla massima considerata.


      do j=1,23             

         Q(j)=Q1(j)

      enddo
        

      endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cambiamenti pro Cristina new yields
! in pratica per quegli elementi che ci interessa
! ricalcoliamo lo yields con i suoi valori...
! direi Ossigeno Ferro Azoto e Carbonio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!! Possibili valori di metallicit 
      aa(6)=0.
      aa(7)=1.e-8
      aa(8)=1.e-5
      aa(9)=4.e-3
      aa(10)=2.e-2

      aa(11)=0.0 !!! non corretto sarebbe 0.0001
      aa(12)=0.001
      aa(13)=0.004
      aa(14)=0.008
      aa(15)=2.e-2
!
!!!!!!!

      !!!! Determina in che bins di massa bisogna interpolare

      k=0

      do j=1,13
         
         IF (massac(j).lt.H) K=J 
         
      enddo

      KK=K+1



      !!! Inserisce l'array corretto per le masse

      mm(1)=massac(k)
      mm(2)=massac(kk)

      z=100

      !!!! Determina in che bins di zeta bisogna interpolare

         do j=6,10
            If (aa(j).le.zeta) then
            Z=J
            ZZ=Z+1
            endif
         enddo

      
         met(1)=aa(z)
         met(2)=aa(zz)
         
       

      !!! Interpola i valori corretti per le stelle massicce 
      !!! con i valori di CRISTINA PER OSSIGENO CARBONIO E AZOTO

!!!!!!!!!!!! Interpola su griglia massa valore metallicit pi basso 

      do j=1,14 !carbonio ossigeno e azoto
         DD(1)=W(j,k,z)
         DD(2)=W(j,kk,z)
         call polint(mm,DD,2,H,Q1(j),dy)
         Q1(j)=Q1(j)
      enddo



!!!!!!!!!!!! Interpola su griglia massa valore metallicit pi alto
      if (z.lt.10) then 

      do j=1,14 !C12 O16 e N14
         DD(1)=W(j,k,zz)
         DD(2)=W(j,kk,zz)
         call polint(mm,DD,2,H,Q2(j),dy)
         Q2(j)=Q2(j)
      enddo

!!!!!!!!!!!! Interpola su griglia metallicit i due valori !!!!

      do j=1,14!He4,C12 O16 e N14

         DD(1)=Q1(j)
         DD(2)=Q2(j)
!         call polint(met,DD,2,zeta,Q(j),dy)
      enddo

      else
         do j=1,14              !He4,C12 O16 e N14
!            Q(j)=Q1(j)
         enddo

      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!
!!!!!!!       SAGB!!!!!!
!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!
!!!      if ((h.ge.7.5).and.(h.le.10.5)) then
!!!      
!!!      !!!! Determina in che bins di massa bisogna interpolare
!!!
!!!      k=0
!!!
!!!      do j=1,7
!!!         
!!!         IF (massas(j).lt.H) K=J 
!!!         
!!!      enddo
!!!
!!!      KK=K+1
!!!
!!!
!!!
!!!      !!! Inserisce l'array corretto per le masse
!!!
!!!      mm(1)=massas(k)
!!!      mm(2)=massas(kk)
!!!
!!!      z=100
!!!
!!!      !!!! Determina in che bins di zeta bisogna interpolare
!!!
!!!         do j=11,15
!!!            If (aa(j).le.zeta) then
!!!            Z=J
!!!            ZZ=Z+1
!!!            endif
!!!         enddo
!!!    
!!!         met(1)=aa(z)
!!!         met(2)=aa(zz)
!!!
!!!      !!! Interpola i valori corretti per le stelle SAGB
!!!      !!! con i valori di SIESS
!!!
!!!!!!!!!!!!!!! Interpola su griglia massa valore metallicit pi basso 
!!!
!!!      do j=1,14 
!!!         DD(1)=W(j,k,z)
!!!         DD(2)=W(j,kk,z)
!!!         call polint(mm,DD,2,H,Q1(j),dy)
!!!         Q1(j)=Q1(j)
!!!      enddo
!!!
!!!
!!!
!!!!!!!!!!!!!!! Interpola su griglia massa valore metallicit pi alto
!!!
!!!      do j=1,14 
!!!         DD(1)=W(j,k,zz)
!!!         DD(2)=W(j,kk,zz)
!!!         call polint(mm,DD,2,H,Q2(j),dy)
!!!         Q2(j)=Q2(j)
!!!      enddo
!!!
!!!      do j=1,14
!!!
!!!         DD(1)=Q1(j)
!!!         DD(2)=Q2(j)
!!!         call polint(met,DD,2,zeta,Q(j),dy)
!!!
!!!      enddo
!!!
!!!
!!!      endif
!!!
!!!
!!!
!!!      !!!!
!!!      !!! SAGB
!!!      !!!!
!!!
!!!      if ((h.ge.7).and.(H.le.12).and.(zeta.lt.1.e-5)) then
!!!         q(4)=q(4)*1.
!!!         q(5)=q(5)*1.
!!!      endif
!!!
!!!
!!!
!!!
!!!
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccCC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!     Modification to use properly the Kobayashi yields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      koba=1
      if (koba.eq.1) then


      
!!!! Valori di metallicit in Kobayashi yields

      aa(1)= 0.     
      aa(2)= 0.001  
      aa(3)= 0.004  
      aa(4)= 0.02   
      aa(5)= 0.05   

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Calcolo in che bin di metallicit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=1,5
         If (aa(j).le.zeta) Z=J
      enddo

      zz=z+1
   
      If (z.eq.5) ZZ=Z


      met(1)=aa(z)
      met(2)=aa(ZZ)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Calcolo in che bin di massa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=1,ninputyield
         
         IF (massa(j).lt.H) K=J 
         
      enddo

      KK=K+1
      

      mm(1)=massa(k)
      mm(2)=massa(kk)

    

      !!! Interpola i valori corretti per le stelle massicce 

!!!!!!!!!!!! Interpola su griglia massa valore metallicit pi basso 

      j=9
      DD(1)=W(j,k,z)
      DD(2)=W(j,kk,z)
      call polint(mm,DD,2,H,Q1(j),dy)
      Q1(j)=Q1(j)

      j=21
      DD(1)=W(j,k,z)
      DD(2)=W(j,kk,z)
      call polint(mm,DD,2,H,Q1(j),dy)
      Q1(j)=Q1(j)
      

      j=13
      DD(1)=W(j,k,z)
      DD(2)=W(j,kk,z)
      call polint(mm,DD,2,H,Q1(j),dy)
      Q1(j)=Q1(j)

      if (z.lt.5) then

!!!!!!!!!!!! Interpola su griglia massa valore metallicit pi alto

      j=9
      DD(1)=W(j,k,zz)
      DD(2)=W(j,kk,zz)
      call polint(mm,DD,2,H,Q2(j),dy)
      Q2(j)=Q2(j)
      
      j=21
      DD(1)=W(j,k,zz)
      DD(2)=W(j,kk,zz)
      call polint(mm,DD,2,H,Q2(j),dy)
      Q2(j)=Q2(j)
      
      j=13
      DD(1)=W(j,k,zz)
      DD(2)=W(j,kk,zz)
      call polint(mm,DD,2,H,Q2(j),dy)
      Q2(j)=Q2(j)
  

!!!!!!!!!!!! Interpola su griglia metallicit i due valori !!!!

      j=9
      
      DD(1)=Q1(j)
      DD(2)=Q2(j)
      call polint(met,DD,2,zeta,Q(j),dy)
      j=21
      DD(1)=Q1(j)
      DD(2)=Q2(j)
   !   call polint(met,DD,2,zeta,Q(j),dy)

      j=13
      DD(1)=Q1(j)
      DD(2)=Q2(j)
   !   call polint(met,DD,2,zeta,Q(j),dy)
     

      else 

!!!! Caso in cui la metallicit sia maggiore di quella massima considerata la prende 
!!!   sempre costante e pari alla massima considerata.


      j=9

      Q(j)=Q1(j)
      j=21
!     Q(j)=Q1(j)
      j=13
!     Q(j)=Q1(j)

      endif
      endif
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!     END Modification to use properly the Kobayashi yields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      !!! Azzera i valori per Ba, Eu e La che per alcuni valori delle masse
      !!! sono pari a zero

      qbar=1.e-30
      qeu2=1.e-30
      qla=1.e-30
      qsrr=1.e-30
      qy=1.e-30
      qzr=1.e-30
      qrb=1.e-30
      
      !!! Determina i  valori per Ba, Eu, La, Zr, Sr, Y
! modification for for NSM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! modification for for NSM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mlow1=10. !!!! switch off MRD!!!!

! modification for for NSM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! modification for for NSM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! modification for for NSM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      mup1=30.
      mlow2=8.
      mup2=0. ! set to zero to swith off the r-process contribution from 10 to 30
      
      value1=0.8e-6!*0.6
      value2= .3e-8
      value3= .3e-8
      value4= .3e-8

      costBa=1.00
      costSr=3.16 * 88./138.
      costLa=0.136* 1.
      costZr=2.53 * 90./138.
      costEu=0.117*151./138. 
      costY =1.625* 89./138./3
      costRb =3.16*  86./138.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! CONTRIBUTION by "R-PROCESS" 5% delle massive stars 15-40Msun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF(H.GE.mlow1.and.H.le.mup1) THEN !!! TUTTE LE STELLE MASSICCE POSSONO
         Rando=0.1 ! simply not extracted it is always present OMOG

         if (Rando.le.0.1) then 

            value1=value1!!!Modificato! ! because on average!!! OMOG


!            mbar(1)=mlow1
!            mbar(2)=mup1
!         
!            bar(1)=value1
!            bar(2)=value2
! Barium
!            call polint(mbar,bar,2,H,qbar,dy)

!            value1=(.99*rando/0.1+.01)*value1*2
  
!            value1=qbar

            qbar=value1*costBa
!            write(*,*) 'r-process!'
 
! Europium
 
            qeu2=value1*costEu

! Lanthanum

            qla=value1*costLa

! Strontium

            qsrr=value1*costSr          
!  Yttrium          

            qY=value1*costY          
!  Zirconium          

            qZr=value1*costZr

            qRb=value1*costRb
            
         endif
      endif
      
         Q(24)=qla 
         Q(25)=qbar
         Q(26)=qeu2
         Q(27)=qsrr
         Q(28)=qy
         Q(29)=qzr
         Q(30)=qrb


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! Yields di Urs per Yttrium Strontium and Barium
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         aa(1)= 1.4e-2
         aa(2)= 1.e-3 
         aa(3)= 1.e-5  
         
        if ((h.gt.15).and.(h.lt.80).and.(zeta.gt.1.e-30)) then
            
            k=4
            do j=1,3
               if ((H.ge.Mba(j)).and.(H.le.Mba(j+1))) k=j
            enddo
      
            if (k.eq.4) then
               if (h.ge.40) then
                  mm(1)=MBa(4)
                  mm(2)=MBa(4)
                  k=4
                  kk=4
               else
                  mm(1)=MBa(1)
                  mm(2)=MBa(1)
                  k=1
                  kk=1
               endif

            else
               k=k
               kk=k+1
               mm(1)=MBa(k)
               mm(2)=MBa(kk)
            endif

!!!! Part added in order to use the solar metallicity yields for 
!!!! stars above the solar metallicity.



             
            if (zeta.lt.1.e-5) then
      
               
               Dd(1)=WBa(k,3)
               Dd(2)=WBa(kk,3)
               call polint(mm,DD,2,H,qbar,dy)
               q(25)=q(25)+qbar

               
               dD(1)=WSr(k,3)
               dD(2)=WSr(kk,3)
               call polint(mm,DD,2,H,qbar,dy)
               q(27)=q(27)+qbar
               

               dD(1)=WY(k,3)
               dD(2)=WY(kk,3)
               call polint(mm,DD,2,H,qbar,dy)               
               q(28)=q(28)+qbar
               


               Dd(1)=WLa(k,3)
               Dd(2)=WLa(kk,3)
               call polint(mm,DD,2,H,qbar,dy)
               q(24)=q(24)+qbar

               
               dD(1)=WZr(k,3)
               dD(2)=WZr(kk,3)
               call polint(mm,DD,2,H,qbar,dy)
               q(29)=q(29)+qbar
               

               dD(1)=WRb(k,3)
               dD(2)=WRb(kk,3)
               call polint(mm,DD,2,H,qbar,dy)               
               q(30)=q(30)+qbar


               dD(1)=WEu(k,3)
               dD(2)=WEu(kk,3)
               call polint(mm,DD,2,H,qbar,dy)               
               q(26)=q(26)+qbar




            endif
               
                  ! FInd in which (of the 2: 10^-5--10^-3---10^-2 ) intervals of metallicity we interpolate 
            if ((zeta.ge.1.e-5).and.(zeta.lt.1.4e-2)) then 

                  z=100
                  do j=1,2
                     if ((zeta.gt.aa(j+1)).and.(zeta.le.aa(j))) z=j
                  enddo

                  

                  zz=z+1

                  met(1)=(aa(zz))
                  met(2)=(aa(z))
         
                  

                  ! Interpolate the mass lower grid of metallicity (z)


                  Dd(1)=WBa(k,zz)
                  Dd(2)=WBa(kk,zz)
                  call polint(mm,DD,2,H,qbar,dy)
                  q1(25)=qbar

                  dD(1)=WSr(k,zz)
                  dD(2)=WSr(kk,zz)
                  call polint(mm,DD,2,H,qbar,dy)
                  q1(27)=qbar
               
                  dD(1)=WY(k,zz)
                  dD(2)=WY(kk,zz)
                  call polint(mm,DD,2,H,qbar,dy)
                  q1(28)=qbar


                  Dd(1)=WLa(k,zz)
                  Dd(2)=WLa(kk,zz)
                  call polint(mm,DD,2,H,qbar,dy)
                  q1(24)=qbar
               
               
                  dD(1)=WZr(k,zz)
                  dD(2)=WZr(kk,zz)
                  call polint(mm,DD,2,H,qbar,dy)
                  q1(29)=qbar
               

                  dD(1)=WRb(k,zz)
                  dD(2)=WRb(kk,zz)
                  call polint(mm,DD,2,H,qbar,dy)               
                  q1(30)=qbar

                  dD(1)=WEu(k,zz)
                  dD(2)=WEu(kk,zz)
                  call polint(mm,DD,2,H,qbar,dy)               
                  q1(26)=qbar




                  !     Interpolate the mass lower grid of metallicity (z)

                  Dd(1)=WBa(k,z)
                  Dd(2)=WBa(kk,z)
                  call polint(mm,DD,2,H,qbar,dy)
                  q2(25)=qbar

                  dD(1)=WSr(k,z)
                  dD(2)=WSr(kk,z)
                  call polint(mm,DD,2,H,qbar,dy)
                  q2(27)=qbar
               
                  dD(1)=WY(k,z)
                  dD(2)=WY(kk,z)
                  call polint(mm,DD,2,H,qbar,dy)
                  q2(28)=qbar

                  Dd(1)=WLa(k,z)
                  Dd(2)=WLa(kk,z)
                  call polint(mm,DD,2,H,qbar,dy)
                  q2(24)=qbar
               
               
                  dD(1)=WZr(k,z)
                  dD(2)=WZr(kk,z)
                  call polint(mm,DD,2,H,qbar,dy)
                  q2(29)=qbar
               

                  dD(1)=WRb(k,z)
                  dD(2)=WRb(kk,z)
                  call polint(mm,DD,2,H,qbar,dy)               
                  q2(30)=qbar   
                  
                  dD(1)=WEu(k,z)
                  dD(2)=WEu(kk,z)
                  call polint(mm,DD,2,H,qbar,dy)               
                  q2(26)=qbar






                  !     Interpolate the between the 2 metallicity z and zz

                  Dd(1)=Q1(25)
                  Dd(2)=Q2(25)
                  call polint(met,DD,2,(zeta),qbar,dy)
                  q(25)=q(25)+qbar
 
                  


                  Dd(1)=Q1(27)
                  Dd(2)=Q2(27)
                  call polint(met,DD,2,(zeta),qbar,dy)
                  q(27)=q(27)+qbar

 
                  Dd(1)=Q1(28)
                  Dd(2)=Q2(28)
                  call polint(met,DD,2,(zeta),qbar,dy)
                  q(28)=q(28)+qbar

                  Dd(1)=Q1(24)
                  Dd(2)=Q2(24)
                  call polint(met,DD,2,(zeta),qbar,dy)
                  q(24)=q(24)+qbar

                  Dd(1)=Q1(29)
                  Dd(2)=Q2(29)
                  call polint(met,DD,2,(zeta),qbar,dy)
                  q(29)=q(29)+qbar

                  Dd(1)=Q1(30)
                  Dd(2)=Q2(30)
                  call polint(met,DD,2,(zeta),qbar,dy)
                  q(30)=q(30)+qbar

                  Dd(1)=Q1(26)
                  Dd(2)=Q2(26)
                  call polint(met,DD,2,(zeta),qbar,dy)
                  q(26)=q(26)+qbar

              
               else

                  if (zeta.gt.1.4e-2) then

                     Dd(1)=WBa(k,1)
                     Dd(2)=WBa(kk,1)
                     call polint(mm,DD,2,H,qbar,dy)
                     q(25)=q(25)+qbar
               
               
                     dD(1)=WSr(k,1)
                     dD(2)=WSr(kk,1)
                     call polint(mm,DD,2,H,qbar,dy)
               
                     q(27)=q(27)+qbar
               
                     dD(1)=WY(k,1)
                     dD(2)=WY(kk,1)
                     call polint(mm,DD,2,H,qbar,dy)
               
                     q(28)=q(28)+qbar



                     Dd(1)=WLa(k,1)
                     Dd(2)=WLa(kk,1)
                     call polint(mm,DD,2,H,qbar,dy)
                     q(24)=q(24)+qbar
                     
               
                     dD(1)=WZr(k,1)
                     dD(2)=WZr(kk,1)
                     call polint(mm,DD,2,H,qbar,dy)
                     q(29)=q(29)+qbar
               

                     dD(1)=WRb(k,1)
                     dD(2)=WRb(kk,1)
                     call polint(mm,DD,2,H,qbar,dy)               
                     q(30)=q(30)+qbar

                     dD(1)=WEu(k,1)
                     dD(2)=WEu(kk,1)
                     call polint(mm,DD,2,H,qbar,dy)               
                     q(26)=q(26)+qbar






                  endif


            endif
         endif
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c$$$  c$$$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c$$$  c$$$C     FINE r-process and s-process in massive  STARS
c$$$  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!         q(28)=1.e-30
!         q(25)=1.e-30
!         q(27)=1.e-30
!         q(30)=1.e-30
!
c$$$  c$$$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
c$$$  c$$$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c$$$  c$$$C     INIZIO s-process in STARS
c$$$  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         if (H.ge.1.3.and.H.le.3) then
            
            call bario(zeta,H,QBAr,qsrr,qy,qeu2,qzr,qla,qrb)
            

         endif


          if (H.ge.1..and.H.le.6.) then
            
            call litio(zeta,H,Qli)
            
            if (qli.gt.1.e-20) Q(31)=qli
            
         endif

         If (H.ge.1.3.and.H.le.3) then
            Q(24)=qla/2. 
            Q(25)=qbar/2. 
            Q(26)=qeu2/2.
            Q(27)=qsrr/2.
            Q(28)=qy/2.  
            Q(29)=qzr/2. 
            Q(30)=qrb/2. 
         endif


         
         if ((H.ge.12).and.(H.le.50)) then
            q(9)=0.07
         else
            q(9)=1.e-20
         endif

         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Se la label SNIa  pari a 1 si assume 
!!    che crei una SNIa e ne produca gli Yield
!!    alla fine della sua vita!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      
    
         Hecore=0.0
         do i=1,31
            if (H.lt.0.5) then 

               Q(i)=0.0
               if (i.eq.14) q(i)=H
            endif

            if (binmax.gt.0.0) then
               Q(i)=QIa(i)*ratio+q(i) ! l'output della primaria 
               HECORE=HECORE+Q(i) 
            else
               Q(i)=Q(i)
               Hecore=Hecore+Q(i)
            endif

         enddo
      else
! modification for for NSM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! modification for for NSM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (binmax.ge.-8) then

            Q(31)= 2.E-6*4.*1.

         else

            value1=0.8e-6*20
            costBa=1.00
            costSr=3.16 * 88./138.
            costLa=0.136* 1.
            costZr=2.53 * 90./138.
            costEu=0.117*151./138. 
            costY =1.625* 89./138./3
            costRb =3.16*  86./138


            
            Q(24)= value1*costLa
            Q(25)= value1*costBa
            Q(26)= value1*costEu
            Q(27)= value1*costSr
            Q(28)= value1*costY 
            Q(29)= value1*costZr
            Q(30)= value1*costRb
         endif
! modification for for NSM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! modification for for NSM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! modification for for NSM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      endif
         
      return


      end
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=1000)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
         dift=abs(x-xa(i))
         if (dift.lt.dif) then
            ns=i
            dif=dift
         endif
         c(i)=ya(i)
         d(i)=ya(i)
 11   continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
         do 12 i=1,n-m
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
            if(den.eq.0.) exit !pause 'failure in polint'
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
 12      continue
         if (2*ns.lt.n-m)then
            dy=c(ns+1)
         else
            dy=d(ns)
            ns=ns-1
         endif
         y=y+dy
 13   continue
      return
      END
      subroutine Litio (Zcerc2,massa,qli)
      implicit none
      include 'Letturanew2.inc'
      
      real yd(2),yd2(2)
      real qbav(2)
      real massav(2),zv(2)
      integer i,indiceM,indiceZ
      real massa,zcerc,zcerc2,y,dy,qli
      
      do i=1,15
         if (massa.ge.massali(i)) indiceM=i
      enddo
      
      massav(1)=massali(indiceM) ! 2 mass grid points between which interpolate
      massav(2)=massali(indiceM+1)
      
      yd(1)=YLi(indiceM,1) ! 2 yields grid points for lower mass between which interpolate
      yd(2)=YLi(indiceM+1,1)

      call polint(massav,yd,2,massa,qli,dy)

      return
      end
      subroutine Bario (Zcerc2,massa,qba,qsr,qy,qeu,qzr,qla,qrb) 
      implicit none
      include 'Letturanew2.inc'
      
      real yd(2),yd2(2)
      real qbav(2)
      real massav(2),zv(2)
      integer i,indiceM,indiceZ
      real massa,zcerc,zcerc2,y,dy,qba,qsr,qy,qeu,qzr,qla,qrb
      
      zcerc=zcerc2

      indicem=100
      indicez=100
      
      do i=1,5
         if (massa.ge.massaba(i)) indiceM=i
      enddo
      
      if (zcerc.lt.zbario(1)) zcerc=zbario(1) !!! Decidiamo di approssimare cosi' 
      if (zcerc.gt.zbario(9)) zcerc=zbario(9) !!! metallicita' estreme, potremmo anche prendere 0 come yield

      
      do i=1,9
         if (zcerc.ge.zbario(i)) indiceZ=i
      enddo
      
    

      zv(1)=zbario(indiceZ)    ! 2 metallicities grid points between which interpolate
      zv(2)=zbario(indiceZ+1)      !

      massav(1)=massaba(indiceM) ! 2 mass grid points between which interpolate
      massav(2)=massaba(indiceM+1)

      yd(1)=ba(indiceM,indiceZ) ! 2 yields grid points for lower mass between which interpolate
      yd(2)=ba(indiceM,indiceZ+1)

      yd2(1)=ba(indiceM+1,indiceZ) ! 2 yields grid points for higher mass between which interpolate
      yd2(2)=ba(indiceM+1,indiceZ+1)

!!! Interpolation between 2 metallicity for the lower mass

      call polint(zv,yd,2,zcerc,y,dy)
      qbav(1) = y

!!! Interpolation between 2 metallicity for the higher mass

      call polint(zv,yd2,2,zcerc,y,dy)
      qbav(2) = y

C
C E interpola sua volta fra le due masse la massa richiesta!
C            
             
      call polint(massav,qbav,2,massa,y,dy)
C
C che viene definita appunto come lo yield cercato
C
          
      qba=y

!!! Yttrium

      yd(1)=yt(indiceM,indiceZ) ! 2 yields grid points for lower mass between which interpolate
      yd(2)=yt(indiceM,indiceZ+1)

      yd2(1)=yt(indiceM+1,indiceZ) ! 2 yields grid points for higher mass between which interpolate
      yd2(2)=yt(indiceM+1,indiceZ+1)

!!! Interpolation between 2 metallicity for the lower mass

      call polint(zv,yd,2,zcerc,y,dy)
      qbav(1) = y

!!! Interpolation between 2 metallicity for the higher mass

      call polint(zv,yd2,2,zcerc,y,dy)
      qbav(2) = y

C
C E interpola sua volta fra le due masse la massa richiesta!
C            
             
      call polint(massav,qbav,2,massa,y,dy)
C
C che viene definita appunto come lo yield cercato
C
          
      qy=y


!!! STRONTIUM

      yd(1)=sr(indiceM,indiceZ) ! 2 yields grid points for lower mass between which interpolate
      yd(2)=sr(indiceM,indiceZ+1)

      yd2(1)=sr(indiceM+1,indiceZ) ! 2 yields grid points for higher mass between which interpolate
      yd2(2)=sr(indiceM+1,indiceZ+1)

!!! Interpolation between 2 metallicity for the lower mass

      call polint(zv,yd,2,zcerc,y,dy)
      qbav(1) = y

!!! Interpolation between 2 metallicity for the higher mass

      call polint(zv,yd2,2,zcerc,y,dy)
      qbav(2) = y

C
C E interpola sua volta fra le due masse la massa richiesta!
C            
             
      call polint(massav,qbav,2,massa,y,dy)
C
C che viene definita appunto come lo yield cercato
C
          
      qsr=y


!!! Europium???

      yd(1)=eu(indiceM,indiceZ) ! 2 yields grid points for lower mass between which interpolate
      yd(2)=eu(indiceM,indiceZ+1)

      yd2(1)=eu(indiceM+1,indiceZ) ! 2 yields grid points for higher mass between which interpolate
      yd2(2)=eu(indiceM+1,indiceZ+1)

!!! Interpolation between 2 metallicity for the lower mass

      call polint(zv,yd,2,zcerc,y,dy)
      qbav(1) = y

!!! Interpolation between 2 metallicity for the higher mass

      call polint(zv,yd2,2,zcerc,y,dy)
      qbav(2) = y

C
C E interpola sua volta fra le due masse la massa richiesta!
C            
             
      call polint(massav,qbav,2,massa,y,dy)
C
C che viene definita appunto come lo yield cercato
C
          
      qeu=y
!!! ZIRCONIUM

      yd(1)=zr(indiceM,indiceZ) ! 2 yields grid points for lower mass between which interpolate
      yd(2)=zr(indiceM,indiceZ+1)

      yd2(1)=zr(indiceM+1,indiceZ) ! 2 yields grid points for higher mass between which interpolate
      yd2(2)=zr(indiceM+1,indiceZ+1)

!!! Interpolation between 2 metallicity for the lower mass

      call polint(zv,yd,2,zcerc,y,dy)
      qbav(1) = y

!!! Interpolation between 2 metallicity for the higher mass

      call polint(zv,yd2,2,zcerc,y,dy)
      qbav(2) = y

C
C E interpola sua volta fra le due masse la massa richiesta!
C            
             
      call polint(massav,qbav,2,massa,y,dy)
C
C che viene definita appunto come lo yield cercato
C
          
      qzr=y

!!! Lanthanum

      yd(1)=la(indiceM,indiceZ) ! 2 yields grid points for lower mass between which interpolate
      yd(2)=la(indiceM,indiceZ+1)

      yd2(1)=la(indiceM+1,indiceZ) ! 2 yields grid points for higher mass between which interpolate
      yd2(2)=la(indiceM+1,indiceZ+1)

!!! Interpolation between 2 metallicity for the lower mass

      call polint(zv,yd,2,zcerc,y,dy)
      qbav(1) = y

!!! Interpolation between 2 metallicity for the higher mass

      call polint(zv,yd2,2,zcerc,y,dy)
      qbav(2) = y

C
C E interpola sua volta fra le due masse la massa richiesta!
C            
             
      call polint(massav,qbav,2,massa,y,dy)
C
C che viene definita appunto come lo yield cercato
C
          
      qla=y

!!! Rubidium

      yd(1)=rb(indiceM,indiceZ) ! 2 yields grid points for lower mass between which interpolate
      yd(2)=rb(indiceM,indiceZ+1)

      yd2(1)=rb(indiceM+1,indiceZ) ! 2 yields grid points for higher mass between which interpolate
      yd2(2)=rb(indiceM+1,indiceZ+1)

!!! Interpolation between 2 metallicity for the lower mass

      call polint(zv,yd,2,zcerc,y,dy)
      qbav(1) = y

!!! Interpolation between 2 metallicity for the higher mass

      call polint(zv,yd2,2,zcerc,y,dy)
      qbav(2) = y

C
C E interpola sua volta fra le due masse la massa richiesta!
C            
             
      call polint(massav,qbav,2,massa,y,dy)
C
C che viene definita appunto come lo yield cercato
C
          
      qrb=y





      
      
      return
      end
       end module interpolation_mod
