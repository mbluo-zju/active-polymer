cccc   output conformation of polymer chains after removed some non-adsorbed chains
c
       subroutine conf1_output(fname_cfg,T_solv,Nsamp,npoly,lchain,
     &	  n_bead_polymer,s_size,rbead,rbeadcont,
     &         species,iconf,kT,k_snap,Nsnap,g5data,nGn,dirbead)

      implicit double precision (a-h,o-z)
      character*35 fname_cfg
      real*8 s_size(3),g5data(nGn),T_solv,dirbead(3)
      real*8 rbead(3,n_bead_polymer),rbeadcont(3,n_bead_polymer)
      integer*4 species(n_bead_polymer) 

c 2  configuration file
c
        open(22,file=fname_cfg,position='append') 
	   if(iconf.eq.1.and.k_snap.eq.1) then
         write(22,881)KT,T_solv,Nsamp,npoly,lchain,
     &	  n_bead_polymer,s_size,Nsnap
	   endif
	   write(22,'(2i9,3f14.6)')iconf,k_snap,dirbead(1:3)
	   write(22,998)g5data(1:9)
            k=0
          do kp=1,npoly
          do kc=1,lchain 
            k=k+1
            write(22,995)k,kp,kc,rbeadcont(1,k),
     &         rbeadcont(2,k),rbeadcont(3,k),species(k)
          enddo
          enddo
  
       close(22)
881    format(i7,f10.4,4i7,3f10.3,2i8)	  
995    format(3i7,3f17.8,i6)
998     format(2f14.4,2f14.5,3f14.6,f12.2,f18.7)
       return
       end
