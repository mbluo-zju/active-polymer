     
c  polymer-surface energy     
        subroutine polymer_energy(n_bead_polymer,npoly,lchain
     &  	       ,Surf_type
     &            ,rbead,ENEps,r_cut_ps,sigma,Energy,z_ads,Num_ads)
	   real*8 rbead(3,n_bead_polymer)
	   real*8 ENEps,r_cut_ps,Energy,z_ads,sigma,a,c2,x,x2
	   integer*4  Surf_type


	     Num_ads=0
	     Energy=0.0d0
 	   do i=1,n_bead_polymer
	     if(mod(i,lchain).eq.0) then
	       cycle
		 endif

c lower surface at z = 0
	      rr_w=rbead(3,i)
            if(rr_w.lt.r_cut_ps) then
cc 1 - LJ
	       if(Surf_type.eq.1) then
		      x=(sigma/rr_w)**2
		      x2=x*x
                a=(5.0d0/3.0d0)*(2.5**(2.0d0/3.0d0))*ENEps
                c2=a*x2*(x*x2-1.0d0)
	       elseif(Surf_type.eq.2) then
 		      x=(sigma/rr_w)**2
		      x2=x*x
                a=(5.0d0/3.0d0)*ENEps
                c2=a*x2*(0.4d0*x*x2-1.0d0)
	       elseif(Surf_type.eq.3) then
		      x=(sigma/rr_w)**3
		      x2=x*x
                a=(1.5d0)*sqrt(0.4d0)*ENEps
                c2=a*x*((2.0d0/15.0d0)*x2-1.0d0)
	       elseif(Surf_type.eq.4) then
		      x=(sigma/rr_w)**3
		      x2=x*x
                a=(1.5d0)*sqrt(0.2d0)*ENEps
                c2=a*x*((1.0d0/15.0d0)*x2-1.0d0)
	       endif

	          Energy=Energy+c2

	         if(rr_w.le.z_ads) then
                 Num_ads=Num_ads+1
               endif
	      endif
	    enddo
	   return
	  end   
