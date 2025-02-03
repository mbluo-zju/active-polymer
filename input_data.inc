c Read simulation parameters from file

c===================================================================================
c input data
c===================================================================================
      subroutine input_data(fname_input,fname_out,kcfg,fnamecfgin
     &                ,s_size,nsamp,npoly,lchain,nbead_type
     &                ,n_bead_polymer,elsc,d_bend,thtat0,z_ads,z_min
     &                ,bead_size,itime_equil,itime_sta
     &                ,itime_record,itime_snap
     &                ,T_solv,ENE,r_cut,r_cut_ps,force_sd,D_r,D_t
     &                ,J_head,deltat,rm12,rm,NENEps0,NENEps
     &                ,ENE_ps,NEmax,Surf_type)

          implicit double precision (a-h,o-z)
          character*35 fname_input,fname_out,fnamecfgin
          integer*4 nsamp
          integer*4 n_nano,npoly,lchain,nbead_type,Surf_type
          integer*4 iseed1,iseed2
          integer*4 itime_equil,itime_sta,itime_record
          real*8 elsc,rm,rm12,d_min,deltat
          real*8 d_bend,thtat0,force_sd,ENE_ps(NEmax),J_head
          real*8 s_size(3)
          real*8 ENE(nbead_type,nbead_type),r_cut(nbead_type,nbead_type)
     &		,bead_size(nbead_type)
  
          open(12,file=fname_input)
          read(12,*) fname_out
          read(12,*) kcfg,fnamecfgin

          read(12,*) nsamp
          read(12,*) s_size(1),s_size(2),s_size(3)
c size of the system
         

          read(12,*) npoly,lchain
c n_nano is the number of nanoparticles
c npoly is the number of polymer chains
c lchain is chain length
                
          n_bead_polymer=npoly*lchain
          read(12,*) bead_size(1)
          read(12,*) elsc
c coefficient of the FENE bond
          read(12,*) rm,rm12
c bond length max, bond between 1,2
          read(12,*) dbend,theta0
          
	    read(12,*)ENE(1,1),r_cut(1,1)          

          read(12,*) force_sd
c self-driving force
          read(12,*) D_r,D_t
c inertia 
          read(12,*) J_head
          
          read(12,*) T_solv
c kB*T

c time step
          read(12,*) deltat
          read(12,*) itime_equil,itime_sta,itime_record,itime_snap

c number of time steps needed for equilibration, static and record

	    read(12,*)Surf_type

	   read(12,*)r_cut_ps

	 if(Surf_type.eq.1) then
	   z_ads=1.55
	   z_min=1.15
	 elseif(Surf_type.eq.2) then
	   z_ads=1.324
	   z_min=1.0
	 elseif(Surf_type.eq.3) then
	   z_ads=1.221
	   z_min=0.8585
	 elseif(Surf_type.eq.4) then
	   z_ads=1.086
	   z_min=0.7647
	 endif

c  polymer - surface-particle interaction
         read(12,*)NENEps0,NENEps
         read(12,*)(ENE_ps(k4),k4=NENEps0,NENEps)
c ***********************************************

          CLOSE(12)
        return
      end
