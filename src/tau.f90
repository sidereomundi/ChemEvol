      FUNCTION TAU(H,tautype,binmax)
      real tau
      real H
      real binmax
      integer tautype
      if (tautype.eq.1) then
      IF(H.LE.1.3) THEN
         TAU=(1000*(10**(-0.6545*ALOG10(H)+1.0)))
      ELSE
         IF(H.LE.3.) THEN
            TAU=(1000*(10**(-3.7*ALOG10(H)+1.35)))
         ELSE
            IF(H.LE.7.) THEN
               TAU=(1000*(10**(-2.51*ALOG10(H)+0.77)))
            ELSE
               IF(H.LE.15) THEN
                  TAU=(1000*(10**(-1.78*ALOG10(H)+0.17)))
               ELSE
                  IF(H.LE.60) THEN
                     TAU=(1000*(10**(-0.86*ALOG10(H)-0.94)))
                  ELSE
                     TAU=(1000*(1.2*H**(-1.85)+3.E-3))
                  END IF
               END IF
            END IF
         END IF
      END IF
      else
         if(h.le.0.56) then
            tau=50
         else IF(H.LE.6.6) THEN
            TAU=10**((0.334-SQRT(1.79-0.2232*(7.764-ALOG10(H))))/0.1116)
         ELSE
         TAU=1.2*H**(-1.85)+3.E-3
      END IF

      tau=tau*1000.

      endif

      if (binmax.lt.0.) then
         if (binmax.gt.-8.) then
            tau=tau+500
         else
            tau=tau+4000.
         endif

      endif

      RETURN
      END
