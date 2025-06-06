       module io_mod
       contains
      subroutine leggi
      implicit none
      include 'Letturanew2.inc'      

      real ba134,ba135,ba136,ba137,ba138,ba139,ba140
      integer i,j

      
C Dati in cui si  aggiunto un valore per z=0 pari a
c quello di z=0.00001 ps cambiamento piccolo ma apprezzabile
c una volta sostituite la AMU!!
c
      open (19,file='./YIELDSBA/CristalloBa2.dat',status='old')
      read(19,*)
      do j=1,5
         do i=1,9
            read(19,*) massaBa(j),zbario(i),ba134,ba135,ba136,ba137,
     $           ba138,ba139,ba140
            ba(j,i)=ba134+ba135+ba136+ba137+ba138+ba139+ba140
         enddo
      enddo
      
      close(19)
      
      return
      end


      subroutine leggiSr
      implicit none
      include 'Letturanew2.inc'      

      real Sr86, Sr87, Sr88
      integer i,j

      
      open (19,file='./YIELDSBA/CristalloSr.dat',status='old')
      read(19,*)
      do j=1,5
         do i=1,9
            read(19,*) massaBa(j),zbario(i),Sr86,Sr87,Sr88
            
            Sr(j,i)=Sr86+Sr87+Sr88

            
         enddo
      enddo
   
      close(19)
      
      return
      end


      subroutine leggiY
      implicit none
      include 'Letturanew2.inc'      

      real Y89
      integer i,j

      
      open (19,file='./YIELDSBA/CristalloY.dat',status='old')
      read(19,*)
      do j=1,5
         do i=1,9
            read(19,*) massaBa(j),zbario(i),Y89
            Yt(j,i)=Y89

            
         enddo
         
      enddo
   
      close(19)
      
      return
      end


      subroutine leggiEu
      implicit none
      include 'Letturanew2.inc'      

      real Eu51,Eu52,Eu53,Eu54,Eu55,Eu56,Eu57
      integer i,j

      
      open (19,file='./YIELDSBA/CristalloEu.dat',status='old')
      read(19,*)
      do j=1,5
         do i=1,9
            read(19,*) massaBa(j),zbario(i),Eu51,Eu52,Eu53,Eu54,Eu55,
     $           Eu56,Eu57
            Eu(j,i)=Eu51+Eu52+Eu53+Eu54+Eu55+Eu56+Eu57
 
   
         enddo
         
      enddo
    
      close(19)
      
      return
      end


      subroutine leggiZr
      implicit none
      include 'Letturanew2.inc'      

      real Zr90,Zr91,Zr92,Zr93,Zr94,Zr95,Zr96,Zr97
      integer i,j

      
      open (19,file='./YIELDSBA/CristalloZr.dat',status='old')
      read(19,*)
      do j=1,5
         do i=1,9
            read(19,*) massaBa(j),zbario(i),Zr90,Zr91,Zr92,Zr93,Zr94,
     $           Zr95,Zr96,Zr97

            Zr(j,i)=Zr90+Zr91+Zr92+Zr93+Zr94+Zr95+Zr96+Zr97
   
         enddo
         
      enddo
    
      close(19)
      
      return
      end




      subroutine leggiLa
      implicit none
      include 'Letturanew2.inc'      

      real La139
      integer i,j

      
      open (19,file='./YIELDSBA/CristalloLa.dat',status='old')
      read(19,*)
      do j=1,5
         do i=1,9
            read(19,*) massaBa(j),zbario(i),La139

            La(j,i)=La139


   
         enddo
         
      enddo
    
      close(19)
      
      return
      end

      subroutine leggiRb
      implicit none
      include 'Letturanew2.inc'      

      real Rb85,Rb86,Rb87,Rb88
      integer i,j

      
      open (19,file='./YIELDSBA/CristalloRb.dat',status='old')
      read(19,*)
      do j=1,5
         do i=1,9
            read(19,*) massaBa(j),zbario(i),Rb85,Rb86,Rb87,Rb88

            Rb(j,i)=Rb85+Rb86+Rb87+Rb88


            
         enddo
         
      enddo
      
      close(19)
      
      return
      end

      subroutine leggiLi
      implicit none
      include 'Letturanew2.inc'      
      integer i

      
      open (19,file='./YIELDSBA/KarakasLi.dat',status='old')
      read(19,*)
      do i=1,15
         
         read(19,*) massaLi(i),YLi(i,1),YLi(i,2),YLi(i,3),YLi(i,4)

      enddo
    
      close(19)
      
      return
      end

       end module io_mod
