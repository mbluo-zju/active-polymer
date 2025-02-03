c  Initial polymer conformation

c===================================================================================
c config0
c===================================================================================
      subroutine config0(n_bead_polymer,npoly, lchain,rbead,rbeadcont
     &                  ,s_size,dirbead,vbead,angvbead,species
     &                  ,pi,iseed1,iseed2)
          implicit double precision (a-h,o-z)
          real*8 ran0
          integer*4 n_bead_polymer,npoly, lchain,iseed1,iseed2
          real*8 rbead(3,n_bead_polymer),rbeadcont(3,n_bead_polymer)
          real*8 dirbead(3)
          real*8 vbead(3,n_bead_polymer),angvbead(3)
          real*8 r1xyz(3),theta,phi,TX(3,3),TS(3,3),TZ(3,3)
          real*8 s_size(3),d_min,pi
	    real*8 rtemp(3),rtempc(3)

		integer*4 species(n_bead_polymer),spec_tmp
      
	  nbead=n_bead_polymer

	  r_eq=1.0d0
	  d_length=r_eq

	  d_min=1.2
	 twopi=2.0d0*pi

      do 1000 i=1,npoly

c first monomer in the system
11       k_bead=(i-1)*lchain
          call ranecu(iseed1,iseed2,ran0)
          r1xyz(1)=ran0*s_size(1)
          call ranecu(iseed1,iseed2,ran0)
          r1xyz(2)=ran0*s_size(2)
          r1xyz(3)=1.0d0

	  call check_EV(r1xyz,rbead,Nbead,s_size,k_bead,d_min,key)
	  if (key.eq.1) goto 11

	     k_bead=k_bead+1
           rbead(1:3,k_bead)=r1xyz(1:3)
           rbeadcont(1:3,k_bead)=r1xyz(1:3)

           species(k_bead)=1

          numb=0
	  kbead_temp=k_bead
212     k_bead=kbead_temp
        do 25 j=2,lchain
12            if(j.eq.2) then
	         theta=0.5*pi
		     phi=0.5*pi

                  call TFOR(theta,phi,TX)
                  delx=d_length*TX(1,1)
                  dely=d_length*TX(2,1)
                  delz=d_length*TX(3,1)
                  r1xyz(1)=rbead(1,k_bead)+delx
                  r1xyz(2)=rbead(2,k_bead)+dely
                  r1xyz(3)=rbead(3,k_bead)+delz
              else
                  call ranecu(iseed1,iseed2,ran0)
                  theta=dacos(1.0d0-2.0d0*ran0)
                  call ranecu(iseed1,iseed2,ran0)
                  phi=twopi*ran0
                  
                  call TFOR(theta,phi,TX)
                  call MP2(3,3,3,TS,TX,TZ)

                  delx=d_length*TZ(1,1)
                  dely=d_length*TZ(2,1)
                  delz=d_length*TZ(3,1)
                  r1xyz(1)=rbead(1,k_bead)+delx
                  r1xyz(2)=rbead(2,k_bead)+dely
                  r1xyz(3)=rbead(3,k_bead)+delz
              endif
  
              call PBC(r1xyz,s_size)
	        call check_EV(r1xyz,rbead,Nbead,s_size,k_bead,d_min,key1)
              
	 	  if(r1xyz(3).lt.1.0.or.r1xyz(3).gt.s_size(3)-1.0) then
                key2=1
		  else
	          key2=0
	        endif  

          if(key1.eq.1.or.key2.eq.1) then
	  	  numb=numb+1
		    if(numb.eq.200) then
	          goto 212
	        elseif(numb.eq.3000) then
	          goto 11
	        endif
		  goto 12
           endif

	    k_bead=k_bead+1

          rbead(1,k_bead)=r1xyz(1)
          rbead(2,k_bead)=r1xyz(2)
          rbead(3,k_bead)=r1xyz(3)
	   
          rbeadcont(1,k_bead)=rbeadcont(1,k_bead-1)+delx
          rbeadcont(2,k_bead)=rbeadcont(2,k_bead-1)+dely
          rbeadcont(3,k_bead)=rbeadcont(3,k_bead-1)+delz

          species(k_bead)=1
		
		
		    do i1=1,3
                  do i2=1,3
                      TS(i2,i1)=TX(i2,i1)
                  enddo
              enddo
25        CONTINUE

c	  open(3,file='a.txt')
c		do k1=1,lchain
c	   write(3,'(i8,3f16.8)')k1,rbeadcont(:,k1)
c	enddo
c	 close(3)


c change ID
	    kstart=(i-1)*lchain+1
	    kend=i*lchain

	   do k1=1,lchain/2
	     k2=kstart+k1-1
	     rtemp(1:3)=rbead(1:3,k2)
	     rtempc(1:3)=rbeadcont(1:3,k2)
	     spec_tmp=species(k2)	
	     k3=kend+1-k1
	     rbead(1:3,k2)=rbead(1:3,k3)
	     rbeadcont(1:3,k2)=rbeadcont(1:3,k3)
	     species(k2)=	species(k3)		
	     rbead(1:3,k3)=rtemp(1:3)
	     rbeadcont(1:3,k3)=rtempc(1:3)
	     species(k3)=	spec_tmp		
	   enddo 

c	  open(3,file='b.txt')
c		do k1=1,lchain
c	   write(3,'(i8,3f16.8)')k1,rbeadcont(:,k1)
c	enddo
c	 close(3)

c  direction of the 1st bead
                  dirbead(1)=rbeadcont(1,1)-rbeadcont(1,2)
                  dirbead(2)=rbeadcont(2,1)-rbeadcont(2,2)
                  dirbead(3)=rbeadcont(3,1)-rbeadcont(3,2)
1000	   enddo
	   

          do k1=1,3
              angvbead(k1)=0
              do k2=1,n_bead_polymer
                  vbead(k1,k2)=0
              enddo
          enddo
          print*,'polymer finish'

	 
          return
      end