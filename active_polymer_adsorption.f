c   Program active_polymer_adsorption.f 

c   This codes are written for the adsorption of an active polymer on a homogeneous surface
c   citation: 
c             M. B. Luo and Y. F. Shen, Soft Matter, 2024, 20, 5113-5121. 

c     Two ways for initial conformations:
c                     (1) generate initial conformation; 
c                     (2) read in conformation from the specified file

c   parameters_input.txt   --   parameters for the simulation
c      change this file if you want to simulate for other parameters, like polymer length, active force

c   included files:
c     input_data.inc   ---   read parameters from "parameters_input.txt"
c     config0.inc    ---  generate initial conformation
c     polymer-surface_energy.inc  ---  calculate the polymer-surface interaction energy
c     small_codes.inc   ----  several small codes 
c     values_polymer.inc  ---  calculate polymer's conformation properties
c     conformation_output.inc  --- output polymer's conformation

      parameter(nbead_type=1)
      parameter(nforce=6000,NF=500)
      parameter(Intmax=2147483647)
      parameter(nGn=11,NEmax=50)
      
      implicit double precision (a-h,o-z)
      character*35 fname_input,fname_out,fname_data,fnamecfgin
      character*35 fname_tmp,fname_dfu
 	character*35,allocatable:: fname_cfg(:)
 	character*4 char7,char0

      real*8 deltat0,deltat,sigma,sigmar,d1,d2,vis,visr
      real*8 half_deltat,half_deltat2,lamda_deltat
      real*8 dtime,delta,mpoint0(3)
      real*8 mpoint(3),p0(3),pz(100000)
      integer*4 N_step_one_time 
      
      real*8,allocatable::rbead(:,:),vbead(:,:),afbead(:,:),
     &                    rbeadcont(:,:),Fnoise(:,:),
     &                   g3(:)
      integer*4,allocatable::list(:),head(:),species(:)

      real*8 dirbead(3),angvbead(3),fene_12(3)
      real*8 force_tab1(nforce),force_tab2(nforce)
      integer*4 nixyz(3,27) 

c parameters related to system size:
      real*8 s_size(3),cell(3),xxx10(10) 
      integer*4 mxyz(3) 

	integer*4 Surf_type

c size of beads,ENE,rcutoff(LJ force)
        real*8 ENE(nbead_type,nbead_type),r_cut(nbead_type,nbead_type)
     &		,bead_size(nbead_type),ENE_ps(NEmax),sigma_E

      logical have
      real*8 pi,twopi,e
	real*8 J_head

c  R2,  R2z, Rg2, Rg2z
	 real*8 g5data(nGn),g4(nGn),g5sum(nGn,NEmax),cmpr(3)
	 integer*4 Nconfg
	       
	 data nixyz/0,0,0, 1,0,0, 0,1,0, 0,0,1,
     &     1,1,0, -1,1,0, -1,0,1, 1,0,1, 0,1,1, 0,-1,1, 
     &     1,-1,1, -1,-1,1, -1,1,1, 1,1,1,
     &     -1,0,0, 0,-1,0, 0,0,-1, 
     &     1,-1,0, -1,-1,0, 1,0,-1, -1,0,-1, 0,1,-1, 0,-1,-1,
     &     1,1,-1, -1,-1,-1, -1,1,-1, 1,-1,-1/

      data sd/2.5d0/
      data rcut/6.0d0/
      data sigma_E/1.0d0/

       pi=2.0d0*dasin(1.0d0)
       twopi=2.0d0*pi

         call random_seed()
         call random_number(e)
         iseed1=37872631*e+129833
         call random_number(e)
         iseed2=6872631*e+729077
	
      
c===================================================================================
c read input parameters
c===================================================================================
       fname_tmp='temp.tmp'
       fname_input='parameters_input.txt'

       call input_data(fname_input,fname_out,kcfg,fnamecfgin
     &                ,s_size,nsamp,npoly,lchain,nbead_type
     &                ,n_bead_polymer,elsc,d_bend,thtat0,z_ads,z_min
     &                ,bead_size,itime_equil,itime_sta 
     &                ,itime_record,itime_snap
     &                ,T_solv,ENE,r_cut,r_cut_ps,force_sd,D_r,D_t
     &                ,J_head,deltat0,rm12,rm,NENEps0,NENEps
     &                ,ENE_ps,NEmax,Surf_type)

c      vis=3*pi*eta*bead_size(1)
c      visr=pi*eta*(bead_size(1)**3)
	  vis=T_solv/D_t
        visr=T_solv/D_r

        Nconfg=itime_sta/itime_record
	  Nsnap=itime_sta/itime_snap

        fname_data=Trim(fname_out)//'.dat'
     
	  allocate(fname_cfg(NENEps))
        do k=NENEps0,NENEps
	    if(k.lt.10) then
	      write(char0,'(I1)')k
	      char7='000'//Trim(char0)
	     elseif(k.lt.100) then
	       write(char0,'(I2)')k
	       char7='00'//Trim(char0)
	     elseif(k.lt.1000) then
	       write(char0,'(I3)')k
	       char7='0'//Trim(char0)
	     elseif(k.lt.10000) then
	       write(char7,'(I4)')k
	     endif

	    fname_cfg(k)=Trim(fname_out)//Trim(char7)//'.cfg'
	  enddo

       call create_fortable1(force_tab1,1.0d0,rcut,sigma_E,nforce)
       call create_fortable2(force_tab2,1.0d0,rcut,sigma_E,nforce
     &	 ,Surf_type)
       print*,'force tables are finished!'


c      do i=1,1000
c          print*,force_tab1(i)
c      enddo

        ALLOCATE(rbead(3,n_bead_polymer),rbeadcont(3,n_bead_polymer)
     &        ,vbead(3,n_bead_polymer),Fnoise(3,n_bead_polymer),
     &        afbead(3,n_bead_polymer),species(n_bead_polymer),
     &        list(n_bead_polymer),
     &        g3(itime_sta))

        do i=1,3
          mxyz(i)=int(s_size(i)/sd+0.1)
        enddo
        ncell=mxyz(1)*mxyz(2)*mxyz(3)
        do i=1,3
          cell(i)=DFLOAT(mxyz(i))/s_size(i)
        enddo
        ALLOCATE(head(ncell))
      
      
      write(*,*)'System size: ',int(s_size(1)),int(s_size(2))
     &                         ,int(s_size(3))
	print*,' Surface type    =   ',Surf_type
      write(*,*)'Solvent temperature: ',T_solv
      print*,'Chain number N = ',npoly
      print*,'Chain length n = ',lchain
      print*,'Self-driving force is',force_sd
	print*,'inertia of head J = ', J_head

c start
     
      inquire(file=fname_data,exist=have)
       if(have.eqv..true.) then
          open(2,file=fname_data)
          read(2,*) lchain1,xxx10,nxxx,nxxx,nxxx,iconf1
          close(2)
       else 
          iconf1=0
       endif
      
       if(iconf1.ne.0) then
          open(9,file=fname_data)
          read(9,*) nxxx
 	    do kE=NENEps0,NENEps
		  read(9,*)E_temp,g4
            do k1=1,nGn
	        g5sum(k1,KE)=g4(k1)*iconf1*Nconfg
	  	  enddo
	    enddo
          CLOSE(9)
	  else
	    do kE=NENEps0,NENEps
          do k1=1,nGn
	      g5sum(k1,KE)=0.0d0
	    enddo
	    enddo
       endif
      
c change random seed
        do km=1,iconf1*1037433
          call ranecu(iseed1,iseed2,e)
        enddo
          
      do 2000 iconf=iconf1+1,Nsamp
	  print*
	  print*, ' Sample number =',iconf

	 if(kcfg.eq.0) then
	   print*,'  New configuration ... '
         call config0(n_bead_polymer,npoly, lchain,rbead,rbeadcont
     &                  ,s_size,dirbead,vbead,angvbead,species
     &                  ,pi,iseed1,iseed2)
	 else
	   print*,'  Read configuration from ...  ',fnamecfgin
	   open(13,file=fnamecfgin)
	    read(13,*)KENE,ENE_ps0,Nsamp1,npoly,lchain,
     &	  n_bead_polymer,s_size,Nsnap1
	      do k1=1,iconf-1
	      do k2=1,Nsnap1
	         read(13,*)nxxx,k_snap,dirbead(1:3)
	         read(13,*)xxx
		   kxx=0
              do kp=1,npoly
              do kc=1,lchain 
               kxx=kxx+1
               read(13,*)kxx,kp1,kc1,rbeadcont(1,kxx),
     &           rbeadcont(2,kxx),rbeadcont(3,kxx),species(kxx)
              enddo
              enddo
	      enddo
	      enddo
c  polymer conformation
	         read(13,*)nxxx,k_snap,dirbead(1:3)
	         read(13,*)xxx
		   kxx=0
              do kp=1,npoly
              do kc=1,lchain 
               kxx=kxx+1
               read(13,*)kxx,kp1,kc1,rbeadcont(1,kxx),
     &           rbeadcont(2,kxx),rbeadcont(3,kxx),species(kxx)

	         do kxyz=1,3
                 rbead(kxyz,kxx)=rbeadcont(kxyz,kxx)

                if(rbead(kxyz,kxx).ge.s_size(kxyz)) then
                  rbead(kxyz,kxx)=mod(rbead(kxyz,kxx),s_size(kxyz)) 
                endif
                if(rbead(kxyz,kxx).lt.0.0d0) then
                  rbead(kxyz,kxx)=mod(rbead(kxyz,kxx),s_size(kxyz)) 
     &				+s_size(kxyz)
                endif
	         enddo  	         
              enddo
              enddo
          do k1=1,3
              angvbead(k1)=0.0
              do k2=1,n_bead_polymer
                  vbead(k1,k2)=0.0
              enddo
          enddo
	  close(13)
	 endif

cc   initial relax at small deltat  (otherwise bonds maybe broken)
	  PRINT*,'  initial equilibrium ...'
        deltat=deltat0*0.2
        half_deltat=0.5d0*deltat
        half_deltat2=deltat*half_deltat
        lamda_deltat=0.6*deltat
        sigma=sqrt(2.0d0*T_solv*vis/deltat)
        sigmar=sqrt(2.0d0*T_solv*visr/deltat)
        N_step_one_time=int(1.0/deltat+0.5)

	   ENEps=ENE_ps(1)	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         call box(n_bead_polymer,rbead,ncell,mxyz,cell,head,list)
         call bead_force(rbead,rbeadcont,n_bead_polymer,npoly,lchain
     &               ,species,ENE,r_cut,ENEps,r_cut_ps,
     &                key_error,afbead,sigma,sigma_E
     &               ,dirbead,bead_size,nbead_type,rm,rm12,elsc,fene_12
     &              ,mxyz,force_tab1,force_tab2,nforce,rcut,
     &               head,list,ncell,nixyz
     &               ,vis,iseed1,iseed2,twopi,force_sd
     &               ,s_size,vbead,z_min)

         do it=1,itime_equil/5
              do ikk=1,N_step_one_time
                  
              call relax(rbead,rbeadcont,vbead,afbead,n_bead_polymer
     &               ,npoly,lchain,species,ENE,r_cut,ENEps,r_cut_ps,
     &               deltat,half_deltat2,lamda_deltat,s_size,ncell
     &                ,mxyz,cell,key_error,sigma,sigma_E,
     &               dirbead,bead_size,nbead_type,rm,rm12,elsc,
     &                fene_12,force_tab1,force_tab2,nforce,rcut,nixyz
     &                ,J_head,vis,iseed1,iseed2,twopi
     &                ,force_sd,delta,angvbead,visr,sigmar
     &                ,half_deltat,z_min)
                  
             enddo
         enddo
	  PRINT*,'  initial equilibrium finished.'
ccccccccccc  finish   

        deltat=deltat0
        half_deltat=0.5d0*deltat
        half_deltat2=deltat*half_deltat
        lamda_deltat=0.6*deltat
        sigma=sqrt(2.0d0*T_solv*vis/deltat)
        sigmar=sqrt(2.0d0*T_solv*visr/deltat)

        N_step_one_time=int(1.0/deltat+0.5)

	 do 2222 kE=NENEps0,NENEps
	   ENEps=ENE_ps(kE)	
	   print*,kE,'   ENEps= ',ENEps
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         call box(n_bead_polymer,rbead,ncell,mxyz,cell,head,list)
         call bead_force(rbead,rbeadcont,n_bead_polymer,npoly,lchain
     &               ,species,ENE,r_cut,ENEps,r_cut_ps,
     &                key_error,afbead,sigma,sigma_E
     &               ,dirbead,bead_size,nbead_type,rm,rm12,elsc,fene_12
     &              ,mxyz,force_tab1,force_tab2,nforce,rcut,
     &               head,list,ncell,nixyz
     &               ,vis,iseed1,iseed2,twopi,force_sd
     &               ,s_size,vbead,z_min)

c  equilibrium stage
         do it=1,itime_equil
              do ikk=1,N_step_one_time
                  
              call relax(rbead,rbeadcont,vbead,afbead,n_bead_polymer
     &               ,npoly,lchain,species,ENE,r_cut,ENEps,r_cut_ps,
     &               deltat,half_deltat2,lamda_deltat,s_size,ncell
     &                ,mxyz,cell,key_error,sigma,sigma_E,
     &               dirbead,bead_size,nbead_type,rm,rm12,elsc,
     &                fene_12,force_tab1,force_tab2,nforce,rcut,nixyz
     &                ,J_head,vis,iseed1,iseed2,twopi
     &                ,force_sd,delta,angvbead,visr,sigmar
     &                ,half_deltat,z_min)
                  
             enddo
         enddo
      

c statistics           
          do it=1,itime_sta
              do itt=1,N_step_one_time
               call relax(rbead,rbeadcont,vbead,afbead,n_bead_polymer
     &               ,npoly,lchain,species,ENE,r_cut,ENEps,r_cut_ps,
     &                deltat,half_deltat2,lamda_deltat,s_size,ncell
     &                ,mxyz,cell,key_error,sigma,sigma_E,
     &               dirbead,bead_size,nbead_type,rm,rm12,elsc
     &                ,fene_12,force_tab1,force_tab2,nforce,rcut,nixyz
     &                ,J_head,vis,iseed1,iseed2,twopi
     &                ,force_sd,delta,angvbead,visr,sigmar
     &                ,half_deltat,z_min)
              
              enddo
              
	   if(mod(it,itime_record).eq.0) then
	     print*,'  statistics stage =  ',it

   	    call values_polymer(rbeadcont,N_bead_polymer,
     & 	 lchain,nGn,g5data)
           do k1=1,nGn-4
	       g5sum(k1,KE)=g5sum(k1,KE)+g5data(k1)
	     enddo

          call polymer_energy(n_bead_polymer,npoly,lchain
     &  	       ,Surf_type
     &          ,rbead,ENEps,r_cut_ps,sigma_E,Ene_p_s,z_ads,Num_ads)
c	  print*,Num_ads,Ene_p_s
	       g5data(nGn-3)=Num_ads
	       g5data(nGn-2)=Ene_p_s

 	       g5sum(nGn-3,kE)=g5sum(nGn-3,kE)+Num_ads
	       g5sum(nGn-2,kE)=g5sum(nGn-2,kE)+Num_ads*Num_ads

	       g5sum(nGn-1,kE)=g5sum(nGn-1,kE)+Ene_p_s
	       g5sum(nGn,kE)=g5sum(nGn,kE)+Ene_p_s**2

	   endif

 	   if(mod(it,itime_snap).eq.0) then
 	      k_snap=it/itime_snap
   	    call values_polymer(rbeadcont,N_bead_polymer,
     & 	 lchain,nGn,g5data)
          call polymer_energy(n_bead_polymer,npoly,lchain
     &  	       ,Surf_type
     &          ,rbead,ENEps,r_cut_ps,sigma_E,Ene_p_s,z_ads,Num_ads)
	       g5data(nGn-3)=Num_ads
	       g5data(nGn-2)=Ene_p_s

            call conf1_output(fname_cfg(kE),ENEps,Nsamp,
     &   	    npoly,lchain,n_bead_polymer,s_size,rbead,rbeadcont,
     &         species,iconf,kE,k_snap,Nsnap,g5data,nGn,dirbead)
	   endif
              
         enddo
2222     enddo

c to file Fname_data
	   open(20,file=Fname_data)
	   write(20,996)lchain,s_size,ENE(1,1),elsc,rm,force_sd,
     &       D_r,D_t,J_head,NENEps0,NENEps,Surf_type,itime_sta,iconf
	     do kE=NENEps0,NENEps
	     do k1=1,nGn
	        g4(k1)=g5sum(k1,KE)/(Nconfg*iconf)
	     enddo
              write(20,998)ENE_ps(kE),g4
	     enddo
     	   close(20)	


996	  format(i6,3f10.3,7f10.5,3i6,i9,i6)
997	  format(1x,i8,6f14.2,f17.8)
998     format(f10.4,2f14.4,2f14.5,3f14.6,2(f15.6,f18.5))

2000    enddo

      end
c Main program ends

c  subroutions



c===================================================================================
c relax
c===================================================================================
      subroutine relax(rbead,rbeadcont,vbead,afbead,n_bead_polymer
     &               ,npoly,lchain,species,ENE,r_cut,ENEps,r_cut_ps,
     &                deltat,half_deltat2,lamda_deltat,s_size,ncell
     &                ,mxyz,cell,key_error,sigma,sigma_E,
     &               dirbead,bead_size,nbead_type,rm,rm12,elsc,
     &                fene_12,force_tab1,force_tab2,nforce,rcut,nixyz
     &                ,J_head,vis,iseed1,iseed2,twopi
     &                ,force_sd,delta,angvbead,visr,sigmar
     &                ,half_deltat,z_min)
          implicit double precision (a-h,o-z)
          integer*4 n_bead_polymer,n_nano,ncell,mxyz(3)
          real*8 rbead(3,n_bead_polymer),rbeadcont(3,n_bead_polymer)
          real*8 vbead(3,n_bead_polymer),afbead(3,n_bead_polymer)
          real*8 afxyz1(3,n_bead_polymer),v_int(3,n_bead_polymer)
          real*8 rmove,s_size(3),cell(3),cell_nano(3),fene_12(3)
          real*8 sigma,dirbead(3),rm,rm12,elsc
          real*8 force_tab1(nforce),force_tab2(nforce),half_deltat
          real*8 d_min,vis,twopi,force_sd,delta,visr,sigmar
          real*8 angvbead(3),deltat,half_deltat2,lamda_deltat,J_head
          integer*4 nixyz(3,27)
          integer*4 head(ncell),list(n_bead_polymer)
          integer*4 iseed1,iseed2
          integer*4 key_error 
         real*8 ENE(nbead_type,nbead_type),r_cut(nbead_type,nbead_type)
		real*8 bead_size(nbead_type),r_cut2(nbead_type,nbead_type)
	    real*8 sigma_E
 	    integer*4 species(n_bead_polymer)

c  fixed monomers	
		do k3=1,npoly
	       k4=k3*lchain
	       vbead(1:3,k4)=0.0d0
	       afbead(1:3,k4)=0.0d0
	    enddo

          do kk=1,n_bead_polymer
              do i=1,3
                  rmove=deltat*vbead(i,kk)+half_deltat2*afbead(i,kk)
                  rbead(i,kk)=rbead(i,kk)+rmove
                  rbeadcont(i,kk)=rbeadcont(i,kk)+rmove
                  afxyz1(i,kk)=afbead(i,kk)
                  v_int(i,kk)=vbead(i,kk)
                  vbead(i,kk)=vbead(i,kk)+lamda_deltat*afbead(i,kk)
c                  print*,kk,rmove
              enddo
          enddo
          

c    
          key_error=0
          call relax_dir(dirbead,bead_size,nbead_type,
     &	   	  angvbead,visr,sigmar,twopi
     &              ,fene_12,deltat,n_bead_polymer,iseed1,iseed2
     &              ,J_head,rbeadcont,elsc,rm12,key_error)
          
          call within_box(n_bead_polymer,rbead,s_size)
          call box(n_bead_polymer,rbead,ncell,mxyz,cell,head,list)
          
          call bead_force(rbead,rbeadcont,n_bead_polymer,npoly,lchain
     &               ,species,ENE,r_cut,ENEps,r_cut_ps,
     &                key_error,afbead,sigma,sigma_E
     &               ,dirbead,bead_size,nbead_type,rm,rm12,elsc,fene_12
     &               ,mxyz,force_tab1,force_tab2,nforce,rcut,
     &                 head,list,ncell, nixyz
     &                     ,vis,iseed1,iseed2,twopi,force_sd
     &                     ,s_size,vbead,z_min)
          
          do kk=1,n_bead_polymer
              do i=1,3
                  vbead(i,kk)=v_int(i,kk)+
     &                        half_deltat*(afxyz1(i,kk)+afbead(i,kk))
              enddo
          enddo
          return
      end


c===================================================================================
c bead force
c===================================================================================
      subroutine bead_force(rbead,rbeadcont,n_bead_polymer,npoly,lchain
     &               ,species,ENE,r_cut,ENEps,r_cut_ps,
     &                key_error,afbead,sigma,sigma_E
     &                ,dirbead,bead_size,nbead_type,rm,rm12,elsc,fene_12
     &               ,mxyz,force_tab1,force_tab2,nforce,rcut,
     &                   head,list,ncell,nixyz
     &                     ,vis,iseed1,iseed2,twopi,force_sd
     &                     ,s_size,vbead,z_min )
          implicit double precision (a-h,o-z)
          real*8 ran1,ran2,twopi,force_sd
          integer*4 n_bead_polymer,iseed1,iseed2 
          real*8 rbead(3,n_bead_polymer),rbeadcont(3,n_bead_polymer)
          real*8 afbead(3,n_bead_polymer)
          real*8 vbead(3,n_bead_polymer)
         real*8 ENE(nbead_type,nbead_type),r_cut(nbead_type,nbead_type)
		real*8 bead_size(nbead_type),r_cut2(nbead_type,nbead_type)
          real*8 dirbead(3),s_size(3)
          real*8 xx,yy,zz,rr,rr_sqrt,rm,rm12,elsc,rmef
          real*8 rrr,force_tab1(nforce),force_tab2(nforce),dr_force
          real*8 fene_12(3),delta,vis,sigma,sigma_E,cell_nano(3)
          integer*4 mxyz(3),mx,my,mz,mxy,ncell
          integer*4 head(ncell),list(n_bead_polymer),nixyz(3,27)
	    integer*4 species(n_bead_polymer)
          integer*4 key_error

           dr_force=rcut/nforce
	    rm2=rm*rm

	  do k2=1,Nbead_type
	  do k1=1,Nbead_type
	    r_cut2(k1,k2)=r_cut(k1,k2)*r_cut(k1,k2)
	  enddo
	  enddo

          key_error=0
c initial
          do k=1,n_bead_polymer
              do i=1,3
                  afbead(i,k)=0
              enddo
          enddo
          
c1 FENE
          do k=1,n_bead_polymer-1
           if(mod(k1,lchain).ne.0) then
              k1=k
              k2=k+1
              if(k1.eq.1) then
                  xx=rbeadcont(1,k1)-0.5*bead_size(1)*dirbead(1)
     &               -rbeadcont(1,k2)
                  yy=rbeadcont(2,k1)-0.5*bead_size(1)*dirbead(2)
     &               -rbeadcont(2,k2)
                  zz=rbeadcont(3,k1)-0.5*bead_size(1)*dirbead(3)
     &               -rbeadcont(3,k2)
                  rmef=rm12*rm12
              else
                  xx=rbeadcont(1,k1)-rbeadcont(1,k2)
                  yy=rbeadcont(2,k1)-rbeadcont(2,k2)
                  zz=rbeadcont(3,k1)-rbeadcont(3,k2)
                  rmef=rm2
              endif

              rr=xx**2+yy**2+zz**2
	        rrr=elsc/((rr/rmef)-1.0d0)
             
              add=rrr*xx
              afbead(1,k1)=afbead(1,k1)+add
              afbead(1,k2)=afbead(1,k2)-add
              if(k1.eq.1) fene_12(1)=add
              add=rrr*yy
              afbead(2,k1)=afbead(2,k1)+add
              afbead(2,k2)=afbead(2,k2)-add
              if(k1.eq.1) fene_12(2)=add
              add=rrr*zz
              afbead(3,k1)=afbead(3,k1)+add
              afbead(3,k2)=afbead(3,k2)-add
              if(k1.eq.1) fene_12(3)=add
	     endif
          enddo

          
c2 LJ potential(polymer-polymer)
          mx=mxyz(1)
          my=mxyz(2)
          mz=mxyz(3)
          mxy=mx*my
          
          do iz=1,mz
              do iy=1,my
                  do ix=1,mx
                      
                      icell=(iz-1)*mxy+(iy-1)*mx+ix
                      ii=head(icell)
                      if(ii.gt.0) then
                          do 2000 kcell=1,14
                              i=ii
                              xxper=0.0d0
                              yyper=0.0d0
                              zzper=0.0d0
                              jx=ix+nixyz(1,kcell)
                              jy=iy+nixyz(2,kcell)
                              jz=iz+nixyz(3,kcell)
                              
                              if((ix.eq.mx).and.(jx.gt.ix)) then
                                  jx=1
                                  xxper=-s_size(1)
                              elseif((ix.eq.1).and.(jx.lt.ix)) then
                                  jx=mx
                                  xxper=s_size(1)
                              endif
                              if((iy.eq.my).and.(jy.gt.iy))then
                                  jy=1
                                  yyper=-s_size(2)
                              elseif((iy.eq.1).and.(jy.lt.iy))then
                                  jy=my
                                  yyper=s_size(2)
                              endif
                              if((iz.eq.mz).and.(jz.gt.iz))then
                                  jz=1
                                  zzper=-s_size(3)
                              elseif((iz.eq.1).and.(jz.lt.iz))then
                                  jz=mz
                                  zzper=s_size(3)
                              endif
                              jcell=jx+mx*(jy-1)+mxy*(jz-1)
                              j=head(jcell)
                              if(j.gt.0) then
800                               if(jcell.eq.icell) j=list(i)

                              do while(j.gt.0)
                                    xx=rbead(1,i)-rbead(1,j)+xxper
                                    yy=rbead(2,i)-rbead(2,j)+yyper
                                    zz=rbead(3,i)-rbead(3,j)+zzper
                                    rr=xx**2+yy**2+zz**2

								ksi=species(i)
								ksj=species(j)

                                  if(rr.lt.r_cut2(ksi,ksj)) then
                                       rr_sqrt=dsqrt(rr)
                                       kkk=rr_sqrt/dr_force
                                       c2=ENE(ksi,ksj)*force_tab1(kkk)
                                          add=xx*c2
                                          afbead(1,i)=afbead(1,i)+add
                                          afbead(1,j)=afbead(1,j)-add
                                          add=yy*c2
                                          afbead(2,i)=afbead(2,i)+add
                                          afbead(2,j)=afbead(2,j)-add
                                          add=zz*c2
                                          afbead(3,i)=afbead(3,i)+add
                                          afbead(3,j)=afbead(3,j)-add
                                  endif
                                  j=list(j)
                               enddo
                                  j=head(jcell)
                                  i=list(i)
                                  if(i.gt.0) goto 800
                              endif
2000                      continue
                      endif
                  enddo
              enddo
          enddo
 
         
c5  LJ force because of wall (along z ) 

	    do i=1,n_bead_polymer
	      ksi=species(i)

c lower surface at z = 0
	      rr_w=rbead(3,i)
            if(rr_w.lt.r_cut_ps) then
cc 1 - LJ
	        kkk=rr_w/dr_force
	        c2=ENEps*force_tab2(kkk)
              afbead(3,i)=afbead(3,i)+c2 
	      endif

c upper surface at z = s_size(3) [repulsive ]
 	      rr_w=s_size(3)-rbead(3,i)
            if(rr_w.lt.z_min) then
	        kkk=rr_w/dr_force
	        c2= force_tab2(kkk)
              afbead(3,i)=afbead(3,i)-c2 
	      endif
          enddo
 
c friction and thermal force
          do k=1,n_bead_polymer
              do i=1,3
                  call ranecu(iseed1,iseed2,ran1)
                  call ranecu(iseed1,iseed2,ran2)
                  ran2=twopi*ran2
                  afbead(i,k)=afbead(i,k)+
     &              sigma*dsqrt(2.0d0*abs(log(ran1)))*cos(ran2)
     &              -vis*vbead(i,k)
              enddo
          enddo

c self-driving force
          do i=1,3
              afbead(i,1)=afbead(i,1)+dirbead(i)*force_sd
          enddo

          return
       end


c===================================================================================
c direction relax
c===================================================================================
      subroutine relax_dir(dirbead,bead_size,nbead_type,
     &	   	  angvbead,visr,sigmar,twopi
     &               ,fene_12,deltat,n_bead_polymer,iseed1,iseed2
     &               ,J_head,rbeadcont,elsc,rm_12,key_error)
          real*8 ran0,a,b,sigmar,visr,vecr(3),deltat,rm_12
          real*8 twopi,mnoise(3),mbead(3),angvbead(3)
          real*8 bead_size(nbead_type),dirbead(3),fene_12(3)
          real*8 mfene(3),dangvbead1(3),angvbead1(3)
          real*8 ddirbead1(3),dirbead1(3),xx(3),rr,sqrt_rr
          real*8 rrr,dangvbead2(3),ddirbead2(3),pp,elsc,J_head
          real*8 rbeadcont(3,n_bead_polymer)
          integer*4 iseed1,iseed2

c thermal torque and resistance torque
          do i=1,3
              call ranecu(iseed1,iseed2,ran0)
              a=ran0
              if(a.lt.1.0d-5) a=1.0d-5
              call ranecu(iseed1,iseed2,ran0)
              b=twopi*ran0
              mnoise(i)=sigmar*dsqrt(2.0d0*abs(log(a)))*cos(b)
              mbead(i)=mnoise(i)-visr*angvbead(i)
          enddo
c fene torque
          do i=1,3
              vecr(i)=-0.5d0*bead_size(1)*dirbead(i)
          enddo
          call cross(vecr,fene_12,mfene)
          do i=1,3
              mbead(i)=mbead(i)+mfene(i)
              dangvbead1(i)=mbead(i)*deltat/J_head
c  J_head -- moment of inertia J
              angvbead1(i)=angvbead(i)+dangvbead1(i)
          enddo
          call cross(angvbead,dirbead,ddirbead1)
          pp=0
          do i=1,3
              ddirbead1(i)=ddirbead1(i)*deltat
              dirbead1(i)=dirbead(i)+ddirbead1(i)
              pp=pp+dirbead(i)**2
          enddo
          pp=dsqrt(pp)

          do i=1,3
              dirbead(i)=dirbead(i)/pp
          enddo
          rr=0
          do i=1,3
              xx(i)=rbeadcont(i,1)
     &              -0.5*bead_size(1)*dirbead1(i)
     &              -rbeadcont(i,2)
              rr=rr+xx(i)**2
          enddo

          sqrt_rr=dsqrt(rr)
          if(sqrt_rr.ge.rm_12) key_error=1
          rrr=elsc/(rr/(rm_12**2)-1.0d0)
          do i=1,3
              fene_12(i)=rrr*xx(i)
          enddo
          do i=1,3
              vecr(i)=-0.5d0*bead_size(1)*dirbead1(i)
          enddo
          call cross(vecr,fene_12,mfene)
          do i=1,3
              mbead(i)=mnoise(i)-visr*angvbead1(i)+mfene(i)
              dangvbead2(i)=mbead(i)*deltat/J_head
          enddo
          call cross(angvbead1,dirbead1,ddirbead2)
          do i=1,3
              ddirbead2(i)=ddirbead2(i)*deltat
          enddo
          pp=0
          do i=1,3
              angvbead(i)=angvbead(i)
     &                   +0.5*(dangvbead1(i)+dangvbead2(i))
              dirbead(i)=dirbead(i)
     &                  +0.5*(ddirbead1(i)+ddirbead2(i))
              pp=pp+dirbead(i)**2
          enddo
          pp=dsqrt(pp)
          do i=1,3
              dirbead(i)=dirbead(i)/pp
          enddo
          return
      end


 

	
	 include 'input_data.inc'	
	 include 'config0.inc'
	 include 'small_codes.inc'
	 include 'conformation_output.inc'
	 include 'values_polymer.inc'
	 include 'polymer-surface_energy.inc'