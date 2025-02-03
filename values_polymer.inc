c  conformation properties

c  g5data:
c       R2   R2z   Rg2  Rg2z  zc   b    A_2D
	subroutine values_polymer(rbeadcont,N_bead_polymer,
     & 	   lchain,nGn,g5data)
	
      implicit double precision(a-h,o-z)
      real*8 rbeadcont(3,N_bead_polymer)
		
	 real*8 g5data(nGn)
	 Real*8 cmrp(3),shape2D
	real*8 tyy(2,2),dtem(2)

c	write(*,*)n_bead_polymer,npoly,lchain
	    mhead=1
	    mtail=lchain
c  End-to-End distance          
	    rrx=(rbeadcont(1,mtail)-rbeadcont(1,mhead))**2
	    rry=(rbeadcont(2,mtail)-rbeadcont(2,mhead))**2
          rrz=(rbeadcont(3,mtail)-rbeadcont(3,mhead))**2
          g5data(1)=rrx+rry+rrz
          g5data(2)=rrz

c Centers of mass of polymer chains
		ccm_x=0.d0
          ccm_y=0.d0
          ccm_z=0.d0
          do j=mhead,mtail
            ccm_x=ccm_x+rbeadcont(1,j)
            ccm_y=ccm_y+rbeadcont(2,j)
            ccm_z=ccm_z+rbeadcont(3,j)
	    enddo
          ccm_x=ccm_x/lchain
          ccm_y=ccm_y/lchain
          ccm_z=ccm_z/lchain


  	    cmrp(1)=ccm_x
  	    cmrp(2)=ccm_y
  	    cmrp(3)=ccm_z


c Radius of gyration
	      diffx=0.0d0
	      diffy=0.0d0
	      diffz=0.0d0

          do j=mhead,mtail
            diffx=diffx+(rbeadcont(1,j)-ccm_x)**2
            diffy=diffy+(rbeadcont(2,j)-ccm_y)**2
            diffz=diffz+(rbeadcont(3,j)-ccm_z)**2
          enddo
		  
		  g5data(3)=(diffx+diffy+diffz)/lchain
            g5data(4)=diffz/lchain
 
	      g5data(5)=ccm_z
            
c bondlength
	     bond_length=0.0d0
	     bond_length2=0.0d0

		do j=1,lchain-1
		  kj=j
		  km=kj+1   
	      xx=rbeadcont(1,kj)-rbeadcont(1,km)
	      yy=rbeadcont(2,kj)-rbeadcont(2,km)
	      zz=rbeadcont(3,kj)-rbeadcont(3,km)
	      rr=xx*xx+yy*yy+zz*zz
	      bond_length=bond_length+dsqrt(rr)
c	      bond_length2=bond_length2+rr
      	enddo
		  g5data(6)=bond_length/(lchain-1)
c		  g5data(7)=bond_length2/(lchain-1)
		  

c 2D shape
	do 50 k=1,2
	do 50 kk=1,2
	    tyy(k,kk)=0.0               
50     continue      
	do 40 k=1,lchain
	tyy(1,1)=tyy(1,1)+(rbeadcont(1,k)-ccm_x)*(rbeadcont(1,k)-ccm_x)
	tyy(1,2)=tyy(1,2)+(rbeadcont(1,k)-ccm_x)*(rbeadcont(2,k)-ccm_y)
	tyy(2,2)=tyy(2,2)+(rbeadcont(2,k)-ccm_y)*(rbeadcont(2,k)-ccm_y)
40      continue
	tyy(2,1)=tyy(1,2)
	do 41 k1=1,2
	do 41 k2=1,2
	tyy(k1,k2)=tyy(k1,k2)/lchain
41      continue

c 2D shape
        delta=sqrt((tyy(1,1)-tyy(2,2))**2+4*tyy(1,2)**2)
        dtem(1)=0.5*((tyy(1,1)+tyy(2,2))+delta)
	  dtem(2)=0.5*((tyy(1,1)+tyy(2,2))-delta)
	  shape2D=(dtem(1)-dtem(2))**2/(tyy(1,1)+tyy(2,2))**2

		  g5data(7)=shape2D

	  return
 	  end
