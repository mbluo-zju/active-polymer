c===================================================================================
c PBC
c===================================================================================
      subroutine PBC(r1xyz,s_size)
          real*8 r1xyz(3),s_size(3)
          do k=1,3
              if(r1xyz(k).le.0d0) r1xyz(k)=r1xyz(k)+s_size(k)
              if(r1xyz(k).gt.s_size(k)) r1xyz(k)=r1xyz(k)-s_size(k)
          enddo
          return
      end



c===================================================================================
c within box
c===================================================================================
      subroutine within_box(nbead,rbead,s_size)
          integer*4 nbead
          real*8 rbead(3,nbead),s_size(3)
          do k=1,nbead
              do i=1,3
                  if(rbead(i,k).le.0.0d0) rbead(i,k)=rbead(i,k)+
     &                                               s_size(i)
                  if(rbead(i,k).gt.s_size(i)) rbead(i,k)=rbead(i,k)-
     &                                                   s_size(i)
              enddo
          enddo
          return
      end


c===================================================================================
c cell list 
c===================================================================================
      subroutine box(nbead,rbead,ncell,mxyz,cell,head,list)
          real*8 rbead(3,nbead),cell(3)
          integer*4 mxyz(3),head(ncell),list(nbead)
          
          mx=mxyz(1)
          mxy=mx*mxyz(2)
c          print*,list(7999)
          do 10 i=1,nbead
              list(i)=0
10        CONTINUE
          
          do 20 i=1,ncell
              head(i)=0
20        CONTINUE
          
          do 30 i=1,nbead
              
              icell=1+int(rbead(1,i)*cell(1))
     &               +int(rbead(2,i)*cell(2))*mx
     &               +int(rbead(3,i)*cell(3))*mxy
              
              j=head(icell)
              
              head(icell)=i
              list(i)=j
30        CONTINUE
          
          return
      end



c**********************************************
c check over lap between polymer beads
c**********************************************
      subroutine check_EV(r1xyz,rbead,n_bead_polymer,s_size,k_bead,
     &                   d_min,key_error)
          real*8 rbead(3,N_bead_polymer),s_size(3),r1xyz(3)
          real*8 DX,DY,DZ,d_min,d_min2,R
          real*8 x,y,z
          integer*4 key_error
          
          d_min2=d_min**2
          key_error=0
          x=r1xyz(1)
          y=r1xyz(2)
          z=r1xyz(3)
          
          do 103 i2=1,k_bead-1
              DX=x-rbead(1,i2)
              if(X.lt.d_min) then
                  if((DX.lt.0).AND.((s_size(1)+DX).lt.d_min)) then
                      DX=DX+s_size(1)
                  endif
              endif
              if(X.gt.(s_size(1)-d_min)) then
                  if((DX.gt.0).AND.((s_size(1)-DX).lt.d_min)) then
                      DX=DX-s_size(1)
                  endif
              endif
              if(ABS(DX).gt.d_min) goto 103
              DY=Y-rbead(2,I2)
              IF(Y.LT.d_min) THEN
                  IF((DY.LT.0).AND.((s_size(2)+DY).LT.d_min)) THEN
                      DY=DY+s_size(2)
                  ENDIF
              ENDIF
              IF(Y.GT.(s_size(2)-d_min)) THEN
                  IF((DY.GT.0).AND.((s_size(2)-DY).LT.d_min)) THEN
                      DY=DY-s_size(2)
                  ENDIF
              ENDIF
              IF(ABS(DY).GT.d_min) GOTO 103
              DZ=Z-rbead(3,I2)
              IF(Z.LT.d_min) THEN
                  IF((DZ.LT.0).AND.((s_size(3)+DZ).LT.d_min)) THEN
                      DZ=DZ+s_size(3)
                  ENDIF
              ENDIF
              IF(Z.GT.(s_size(3)-d_min)) THEN
                  IF((DZ.GT.0).AND.((s_size(3)-DZ).LT.d_min)) THEN
                      DZ=DZ-s_size(3)
                  ENDIF
              ENDIF
              IF(ABS(DZ).GT.d_min) GOTO 103
              R=DX**2+DY**2+DZ**2
              IF(R.lt.d_min2) then
                  key_error=1
                  exit
              endif
103       CONTINUE
          return
      end



c**********************************************
c create fortable
c**********************************************      
      subroutine create_fortable1(force_tab1,eps,rcut,sigma,nf)
          real*8 eps,force_tab1(nf),rcut,sigma
          real*8 r,dr,a,x
          dr=rcut/nf
          
          a=24.0d0*eps/(sigma*sigma)
	    r=0.0d0
          do k=1,nf
              r=r+dr
              x=sigma/r
              force_tab1(k)=a*(x**8)*(2.0d0*(x**6)-1.0d0)
          enddo
          return
      end

      
       subroutine create_fortable2(force_tab2,eps,rcut,sigma,nf 
     &	 ,Surf_type)
c  surface
          real*8 eps,force_tab2(nf),rcut,sigma
          real*8 r,dr,a,x,z,dz
	    integer*4 Surf_type

          dz= rcut/nf
	       if(Surf_type.eq.1) then
                a=(10.0d0/3.0d0)*(2.5**(2.0d0/3.0d0))*eps/(sigma*sigma)
	       elseif(Surf_type.eq.2) then
                a=(20.0d0/3.0d0)*eps/(sigma*sigma)
	       elseif(Surf_type.eq.3) then
                a=4.5d0*sqrt(0.4)*eps/(sigma*sigma)
	       elseif(Surf_type.eq.4) then
                a=4.5d0*sqrt(0.2)*eps/(sigma*sigma)
	       endif
	    z=0.0d0
	    do k=1,nf
 	      z=z+dz
	      x=sigma/z
	       if(Surf_type.eq.1) then
               force_tab2(k)=a*(x**6)*(5.0d0*(x**6)-2.0d0)*z
	       elseif(Surf_type.eq.2) then
               force_tab2(k)=a*(x**6)*(x**6-1.0d0)*z
	       elseif(Surf_type.eq.3) then
               force_tab2(k)=a*(x**5)*(0.4d0*(x**6)-1.0d0)*z
	       elseif(Surf_type.eq.4) then
               force_tab2(k)=a*(x**5)*(0.2d0*(x**6)-1.0d0)*z
	       endif
          enddo
       
          return
       end      



c**********************************************
c codes for Random Number end
c**********************************************

      SUBROUTINE ranecu(iseed1,iseed2,ranec1)
          REAL*8 ranec1
          INTEGER*4 iran,iseed1,iseed2,kkk
          kkk = iseed1 / 53668
          iseed1 = 40014 * (iseed1 - kkk * 53668) - kkk * 12211
          IF(iseed1 .LT. 0) iseed1 = iseed1 + 2147483563
          kkk = iseed2 / 52774
          iseed2 = 40692 * (iseed2 - kkk * 52774) - kkk * 3791
          IF(iseed2 .LT. 0) iseed2 = iseed2 + 2147483399
          iran = iseed1 - iseed2
          IF(iran .LT. 1) iran = iran + 2147483562
          ranec1 = iran * 4.656612873D-10
      RETURN
      END
      
      SUBROUTINE MP2(N,M,L,A,B,R)
          DOUBLE PRECISION A(N,M),B(M,L),R(N,L)
          DO 10 i=1,N
              DO 10 j=1,L
                  r(i,j)=0.0d0
                  DO 20 k=1,M
20                    r(i,j)=r(i,j)+A(i,k)*B(k,j)
10        continue
          RETURN
      END

      SUBROUTINE TFOR(theta,phi,X)
          real*8 theta,phi,X(3,3),CO,SI,COP,SIP
          co=cos(theta)
	      SI=sin(theta)
	      COP=cos(phi)
	      SIP=sin(phi)
          X(1,1)=CO
          X(1,2)=SI
          X(1,3)=0.0d0
          X(2,1)=si*COP
          X(2,2)=-co*COP
          X(2,3)=SIP
          X(3,1)=si*SIP
          X(3,2)=-co*SIP
          X(3,3)=-COP
          RETURN
      END
 

c**********************************************
c codes for cross product
c**********************************************
      subroutine cross(a,b,c)
          real*8 a(3),b(3),c(3)
          c(1)=a(2)*b(3)-a(3)*b(2)
          c(2)=a(3)*b(1)-a(1)*b(3)
          c(3)=a(1)*b(2)-a(2)*b(1)
          return
      end
