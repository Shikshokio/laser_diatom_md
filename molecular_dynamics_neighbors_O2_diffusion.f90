!This program is a classical molecular dynamics simulation for a gas of diatomic molecules, with a possible addition of rotational velocity "kicks" by one or several laser pulses. It is also possible to locally "heat" the translational motion of the molecules, or to start from an ensemble of coherently centrifuged molecules. 

!Units: the chosen units are convenient for 14^N_2 molecules. One unit of mass is the mass of a single 14^N_2 molecule (28 g/mol); unit of length is the equilibrium bond length of the nitrogen molecule (1.1 Angstrom); unit of time is the "historically" chosen I_{N_2}/hbar (=1.3 psec). In this weird unit of time, the rotational velocity change caused by the laser is convenient. In this unit of time, the quantum rotational revival period is equal to 2*pi.

!Parameters that can be changed: pressure in atm (p), temperature in K (Temp), number of molecules (small_nx+1) * (small_ny+1) * (small_nz+1)[these variables fix the volume, chosen according to the ideal gas law]. The time step (dt) and the total number of timesteps (Nt) can be also changed.

!Other diatomic molecules can be used. For this purpose several things should be changed: the molecule mass (M_mol), the bond length (r_e) [this determines the moment of inertia I_mol], the Lennard Jones constants for the force calculation [they are used in subroutines "forces_torques" and "pot_energy_tot"] 

!UNITS:
!mass - 0.028 kg/mol / 6.022e23
!length - 1.1 Angstrom
!time - 1.3 psec
!energy - 3.33e-22 J
!angular momentum - 4.33e-34 J*sec

program molecular_dynamics
implicit none
real(kind=8) :: p,Temp,T_heat,Lx,Ly,Lz,Lhalf_x,Lhalf_y,Lhalf_z,dt,dt_trans,dt_rot,t_max,& 
                E_trans_cent, E_rot_cent, E_pot_cent ,P_total(3),laser_range,laser_range_by2, &
            J_trans_cent(3),J_rot_cent(3),P_kick,pol(3), M_mol, re_mol, I_mol,V_max,u_max, &
           const,const2,pot_at_cutoff,sigma,sigmasq,cutoff,cutoffsq,skin,J_rot_sq_cent(3), &
	   xj,xk,yj,yk,correl_jk,x2_cent,y2_cent,z2_cent,xy_cent ,disp2_cent(3),disp4_cent(3)
integer :: small_nx,small_ny,small_nz,N,Nt,i,j,k,N_places=0,MaxPairsPerAtom,N_neigh, &
	maximum_box,N_cent,file_number=0
integer, allocatable :: Marker1(:),Marker2(:),List(:),position(:),number_density(:)
real(kind=8), allocatable :: R_ini(:,:), e_ini(:,:), V_ini(:,:), u_ini(:,:)
real(kind=8), allocatable :: R(:,:),V(:,:),e(:,:),u(:,:),u_half(:,:),F(:,:),DisplaceList(:,:), &
                     G_perp(:,:),lambda(:,:),Omega(:,:),zeros(:,:) ,Jbig(:,:),deltaR(:,:), &
		E_trans(:),E_rot(:),x2(:),y2(:),z2(:),xy(:),J_trans(:,:),J_rot(:,:),displacement(:,:)
logical, allocatable :: places1(:), places3(:,:), places_cent(:)
logical :: ListUpdateRequested=.TRUE.

character :: date*10, time1*10, time2*10, filename_V*100, filename_ang_mom*100,file_number_char*2
integer :: hour1, hour2, minut1, minut2 
integer :: tot_hour, tot_minut
real(kind=8) :: sec1, sec2, tot_sec
real(kind=8), parameter :: pi=3.14159265358979


! Recording the initial time:
call DATE_AND_TIME(date,time1)
read(time1(1:2),1) hour1
read(time1(3:4),1) minut1
1 format(i3)
read(time1(5:10),2) sec1
2 format(f7.3)

! Defining the parameters of the simulation and calculating the initial input:

!Here I define the dimensionless parameters, that define the SITE-SITE INTERACTION POTENTIAL / the forces. In the code I use the simplest Lennard-Jones potential. Here are the values of its 2 parameters for different atoms (from Allen's book, Computer simulation of liquids, p.21):
!------------ epsilon/k_B (Kelvin) -- sigma (Angstrom) --
!Nitrogen         37.3                    3.31           
!Oxigen           61.6                    2.95           (used in this code)
!Chlorine        173.5                    3.35
!Bromine         257.2                    3.54

 const=17.05 != (48*epsilon) * (tau^2)/(sigma^2*M^2) - M and r_e are for ^14 N_2 molecule; see the unit convention above
 sigma=2.6818 !=sigma/r_e
 const2=(sigma**2)*const/12. != 4*epsilon * (tau^2)/(r_e^2*M)
 sigmasq=sigma**2
 cutoff=sigma*2.5! cutoff distance of the Lennard-Jones force 
 cutoffsq=cutoff**2
 pot_at_cutoff=const2*((sigmasq/cutoffsq)**6-(sigmasq/cutoffsq)**3)! value of the potential at cutoff distance

 skin=13.! The skin parameter that can be optimized for minimal running time. For each molecule, a neighbors' list is created. The neighbors of a given molecule "live" in a sphere of radius cutoff+skin This list is valid for certain amount of time (until there is a chance that one molecule from outside the sphere penetrated through the "skin" inside the "cutoff" region). For each molecule, only the neighbors will be considered in the force calculation. 
! If "skin" is too small, update of the neighbors' list will be required to often, making the running time high. If "skin" is too large, the calculation of the neighbors' list will become time consuming, making the running time high again. An optimal value of "skin" should therefore exist.

 MaxPairsPerAtom=100 ! The maximal number of molecules that interact with a given molecule (should be increased when the density increases). This parameter is important for defining vector dimensions.

write(*,*) 'skin=',skin
!!!!!!!!!!!DEFINING THE PRESSURE, TEMPERATURE AND THE NUMBER OF PARTICLES !!!!!!!!!!!

Nt=300000!20000000! total number of timesteps
p=1. !pressure in atm
Temp=300. ! temperature in Kelvin
small_nx=20 !Initially the molecules are arranged in a lattice (with small random displacements) in a rectangular box.
small_ny=20 !"small_nx+1","small_ny+1","small_nz+1" are the number of molecules in 3 perpendicular directions
small_nz=20 !

!!!!!!!!!DEFINING MOLECULAR MASS, BOND LENGTH AND MOMENT OF INERTIA

M_mol=32./28.!A mass of the 16^O_2 molecule (in units of mass of 14^N_2 molecule)
re_mol=1.21/1.1 !Distance between the 2 oxygen "force centers" (in units of the N_2 molecule bond length0)
I_mol=0.25*M_mol*re_mol**2 !Moment of inertia
call initial_cond(p,Temp,small_nx,small_ny,small_nz, &
                   R_ini,e_ini,V_ini,u_ini,Lx,Ly,Lz,N,places_cent,N_cent)
!Here a subroutine is called that produces all the initial conditions: positions (R_ini), orientations (e_ini), velocities (V_ini), angular velocities (u_ini), lengths of the box edges (Lx,Ly,Lz), total number of molecules (N), list of the places of the centrifuged molecules (places_cent), number of the centrifuged molecules (N_cent).


Lhalf_x=Lx/2.
Lhalf_y=Ly/2.
Lhalf_z=Lz/2.
if (Lhalf_x<cutoff.or.Lhalf_y<cutoff.or.Lhalf_z<cutoff) then

	write(*,*) 'The dimensions of the box are too small.'
	write(*,*) 'Choose more molecules or lower pressure.'
	stop

end if

!Allocation of vectors and matrices after the total number of molecules "N" is known
allocate(R(N,3),deltaR(N,3),DisplaceList(N,3),V(N,3),e(N,3),u(N,3),u_half(N,3), &
          F(N,3),G_perp(N,3),Marker1(N),Marker2(N),List(N*MaxPairsPerAtom), &
          lambda(N,3),Omega(N,3),zeros(N,3),Jbig(N,3),places1(N),places3(N,3),position(N), &
	E_trans(N),E_rot(N),x2(N),y2(N),z2(N),xy(N),J_trans(N,3),J_rot(N,3),displacement(N,3))

!The following initializes the running variable matrices R,V,e,u,u_half by their subroutine values (here) or from a "dat" file (below)
R=R_ini
V=V_ini
V=V-spread(sum(V,1)/N,1,N)
e=e_ini
u=u_ini
u_half=u_ini
displacement=0.

! open(60,file='final_cond1e4_T300_10pct_centrifuge.dat')
! read(60,*)
! do j=1,N
!    read(60,*) R(j,:)
! end do
! read(60,*)
! do j=1,N
!    read(60,*) V(j,:)
! end do
! read(60,*)
! do j=1,N
!    read(60,*) e(j,:)
! end do
! read(60,*)
! do j=1,N
!    read(60,*) u(j,:)
! end do
! read(60,*)
! do j=1,N
!    read(60,*) u_half(j,:)
! end do
! do j=1,N
!    read(60,*) displacement(j,:)
! end do
! close(60)
! 
! open(1011,file='cent_num_and_place_1e4_T300_10pct_centrifuge.dat')
! read(1011,*) N_cent
! do i=1,N
! 	read(1011,*) places_cent(i)
! end do
! close(1011)


!open(10,file='inter_cond1000_T100.dat')
laser_range=1e10*1.4*2.*31.32*((1./p)*(Temp/300.))**(1./3.)
zeros=0
places1=abs(R_ini(:,3))<laser_range
places3=spread(places1,2,3)
do j=1,N
   if (places1(j)) N_places=N_places+1
end do
write(*,*) 'N_diff=',N_places, 'N_centrifuged=',N_cent
!open(7,file='diffusion1000_T100_out.dat')

! laser_range_by2=laser_range*2.
! position=nint(R_ini(:,3)/laser_range_by2)
! maximum_box=maxval(abs(position))
! write(*,*) 'The length of the number density vector is:',2*(maximum_box+1)+1
! allocate(number_density(-maximum_box:maximum_box))
! number_density=0
! do j=1,N
! 	number_density(position(j))=number_density(position(j))+1
! end do
! open(11,file='number2e4_T300_out_b.dat')
! do j=-maximum_box,maximum_box
! 	write(11,*) number_density(j)
! end do

!!!!DEFINITION OF THE POLARIZATION AND THE STRENGTH OF THE LASER PULSE!!!!!!!

pol=(/0.,0.,1./) !Polarization direction of the pulse
P_kick=15. !Kick strength of the pulse
!u=laser_pulse(u,e,places3,pol,P_kick)
!The function "laser_pulse" gives instantaneous new angular velocities to laser-kicked molecules, according to molecular old angular velocities (u), molecular orientations (e), a list of the molecules we want to excite, defined above (places3), polarization direction (pol) and laser kick strength (P_kick)

! T_heat=500.
! V=local_heating(V,places3,T_heat)
! The function "local_heating" instantaneously heats the translational velocities "V" of molecules at "places3" to temperature "T_heat"

write(*,*) 'V_ini_max=',maxval(abs(V)),'u_ini_max=', maxval(abs(u))
V_max=maxval(abs(V))
u_max=maxval(abs(u))

!!!DEFINITION OF THE TIME OF THE SIMULATION, AND OF THE TIMESTEP!!!!!!!!

dt_trans=(0.01)/V_max
dt_rot=(0.01)/u_max
dt=min(dt_trans,dt_rot)
!N_neigh=int((skin-cutoff)/2./V_max/dt)
!write(*,*) 'Neighbors check after',N_neigh,'steps'

t_max=Nt*dt !Total time of the simulation
write(*,*) 'dt=',dt

open(111,file='time_1e4_T300_10pct_centrifuge_check.dat') ! Time vector
write(111,*) 0.
do i=1000,Nt,1000
	write(111,*) i*dt
end do
close(111)

!"places_cent" is a logical vector that marks the molecules that were centrifuged
open(101,file='cent_num_and_place_1e4_T300_10pct_centrifuge_check.dat')
write(101,*) N_cent
do i=1,N
	write(101,*) places_cent(i)
end do
close(101)

!stop

! open(8,file='alignment1000_T100_out.dat')
! V_ini(:,3)=merge(e(:,3)**2,zeros(:,3),places1)
! write(8,*) sum(V_ini(:,3))/N_places 


! BELOW WE CALCULATE QUANTITIES AVERAGED OVER ALL THE MOLECULES
E_trans=0.5*M_mol*(V(:,1)**2+V(:,2)**2+V(:,3)**2) !Translational energy
Omega(:,1)=e(:,2)*u(:,3)-e(:,3)*u(:,2) !Angular velocity
Omega(:,2)=e(:,3)*u(:,1)-e(:,1)*u(:,3)
Omega(:,3)=e(:,1)*u(:,2)-e(:,2)*u(:,1)
E_rot=0.5*I_mol*(Omega(:,1)**2+Omega(:,2)**2+Omega(:,3)**2) !Rotational energy
!E_pot=pot_energy_tot(R,e)/N
open(2,file='energy1e4_T300_10pct_centrifuge_check.dat')
write(2,*) sum(E_trans)/N,sum(E_rot)/N !Average translational and rotational energy

open(12,file='disp2_cent1e4_T300_10pct_centrifuge_check.dat')
open(13,file='disp2_cent_std1e4_T300_10pct_centrifuge_check.dat')

! open(222,file='V_dist1e4_T300_10pct_centrifuge.dat') !Initial velocity distribution
! do j=1,N
! 	write(222,*) V(j,:)
! end do
! close(222)

! Quantities describing molecular orientation
! x2=e(:,1)**2
! y2=e(:,2)**2
! z2=e(:,3)**2
! xy=e(:,1)*e(:,2)
! open(33,file='orient1e4_T300_10pct_centrifuge.dat')
! write(33,*) sum(x2)/N,sum(y2)/N,sum(z2)/N,sum(xy)/N

!open(66,file='field1e4_T300_centrifuge.dat')
!write(66,*) sqrt((x2avg-y2avg)**2+4*xyavg**2)

! Correlation factor defined as 2/N/(N-1)*SUM_{i,j,(i<j)} cos(phi_i-phi_j)
! open(44,file='correl1e4_T300_10pct_centrifuge.dat')
!  correl_jk=0.
! do j=1,N-1
!    do k=j+1,N
! 	xj=e(j,1)
! 	xk=e(k,1)
! 	yj=e(j,2)
! 	yk=e(k,2)
! 	correl_jk=correl_jk+(xj*xk+yj*yk)/(sqrt((xj**2+yj**2)*(xk**2+yk**2)))
!    end do
! end do
!  correl_jk=2*correl_jk/(N*(N-1))
! write(44,*) correl_jk

!P_total=M_mol*sum(V,1)/N
! Jbig(:,1)=R(:,2)*V(:,3)-R(:,3)*V(:,2) ! R cross V
! Jbig(:,2)=R(:,3)*V(:,1)-R(:,1)*V(:,3)
! Jbig(:,3)=R(:,1)*V(:,2)-R(:,2)*V(:,1)
! J_trans=M_mol*Jbig !translational angular momentum
! J_rot=I_mol*Omega !rotational angular momentum
! open(3,file='ang_mom_trans1e4_T300_10pct_centrifuge.dat')
! write(3,*) sum(J_trans,1)/N
! open(4,file='ang_mom_rot1e4_T300_10pct_centrifuge.dat')
! write(4,*) sum(J_rot,1)/N
! open(55,file='ang_mom_sq_rot1e4_T300_10pct_centrifuge.dat')
! write(55,*) sum(J_rot**2,1)/N

! open(444,file='ang_mom_rot_dist1e4_T300_10pct_centrifuge.dat')
! do j=1,N
! 	write(444,*) J_rot(j,:)
! end do
! close(444)

! BELOW WE CALCULATE QUANTITIES AVERAGED ONLY OVER THE CENTRIFUGED MOLECULES
E_trans_cent=0.
E_rot_cent=0.
disp2_cent=0.
disp4_cent=0.
! x2_cent=0.
! y2_cent=0.
! z2_cent=0.
! xy_cent=0.
! J_trans_cent=0.
! J_rot_cent=0.
! J_rot_sq_cent=0.
do j=1,N

   if (places_cent(j)) then ! variables are accumulated only for centrifuged molecules

	E_trans_cent=E_trans_cent+E_trans(j)
	E_rot_cent=E_rot_cent+E_rot(j)
	disp2_cent=disp2_cent+displacement(j,:)**2
	disp4_cent=disp4_cent+displacement(j,:)**4
! 	x2_cent=x2_cent+x2(j)
! 	y2_cent=y2_cent+y2(j)
! 	z2_cent=z2_cent+z2(j)
! 	xy_cent=xy_cent+xy(j)
! 	J_trans_cent=J_trans_cent+J_trans(j,:)
! 	J_rot_cent=J_rot_cent+J_rot(j,:)
! 	J_rot_sq_cent=J_rot_sq_cent+J_rot(j,:)**2

   end if

end do
write(2,*) E_trans_cent/N_cent,E_rot_cent/N_cent
write(12,*) disp2_cent/N_cent
write(13,*) sqrt(disp4_cent/N_cent-(disp2_cent/N_cent)**2)
!write(33,*) x2_cent/N_cent,y2_cent/N_cent,z2_cent/N_cent,xy_cent/N_cent
! write(3,*) J_trans_cent/N_cent
! write(4,*) J_rot_cent/N_cent
! write(55,*) J_rot_sq_cent/N_cent

! These commands are explained inside the main loop
call UpdateList(cutoff+skin,R)
ListUpdateRequested=.FALSE.
call forces_torques(R,e,F,G_perp)

! HERE STARTS THE MAIN LOOP OF THE TIME PROPAGATION

do i=1,Nt

!!!PERIODIC BOUNDARY CONDITIONS (REFOLDING)!!!!!!!!	

   where (R(:,1)>Lhalf_x) R(:,1)=R(:,1)-Lx
   where (R(:,1)<-Lhalf_x) R(:,1)=R(:,1)+Lx
   where (R(:,2)>Lhalf_y) R(:,2)=R(:,2)-Ly
   where (R(:,2)<-Lhalf_y) R(:,2)=R(:,2)+Ly
   where (R(:,3)>Lhalf_z) R(:,3)=R(:,3)-Lz
   where (R(:,3)<-Lhalf_z) R(:,3)=R(:,3)+Lz

   !write(*,*) i
   V=V+F*0.5*dt/M_mol !half step in velocity of COM
   deltaR=V*dt
   R=R+deltaR !full step in position of COM
   lambda=spread(u_half(:,1)*e(:,1)+u_half(:,2)*e(:,2)+u_half(:,3)*e(:,3),2,3) !Lagrange multiplier
   u_half=u+G_perp*0.5*dt/I_mol-lambda*e !half step in rot velocity
   e=e+u_half*dt !full step in orientation

   displacement=displacement+deltaR

!   if (mod(i,N_neigh)==0) then
   if (ListUpdateRequested) then !"ListUpdateRequested" will become true if there is a pair of molecules that the sum of their displacements is higher than "skin". In this case there is a chance they moved one towards another. Then the neighbors' list has to be updated.
	call UpdateList(cutoff+skin,R) !The subroutine "UpdateList" is called and the neighbors' list is updated. This subroutine goes over all N*(N-1) pairs. In contrast, the force calculation subroutine goes over N*(average number of neighbors per molecule) pairs. This reduces the running time.
	ListUpdateRequested=.FALSE.
!	write(*,*) i
   end if

   call forces_torques(R,e,F,G_perp) ! Subroutine calculating the forces and the torques
   V=V+F*0.5*dt/M_mol !half step in velocity of COM
   lambda=spread(u_half(:,1)*e(:,1)+u_half(:,2)*e(:,2)+u_half(:,3)*e(:,3),2,3) !Lagrange multiplier
   u=u_half+G_perp*0.5*dt/I_mol-lambda*e !half step in rot velocity

   if (mod(i,1000)==0) then ! Averages are calculated only every 1000 time steps


! 	position=nint(R(:,3)/laser_range_by2)
! 	number_density=0
! 	do j=1,N
! 		number_density(position(j))=number_density(position(j))+1
! 	end do
! 	do j=-maximum_box,maximum_box
! 		write(11,*) number_density(j)
! 	end do

! 		if (mod(i,500000)==0) then ! Velocity distribution is saved every half million time steps
! 
! 			file_number=file_number+1
! 
! 			write(file_number_char,5151) file_number
! 			5151 format(i2)
! 
! 		filename_V='V_dist1e4_T300_10pct_centrifuge_'//file_number_char//'.dat'
! 		filename_V=trim(filename_V)
! 
! 			open(222,file=filename_V)
! 			do j=1,N
! 				write(222,*) V(j,:)
! 			end do
! 			close(222)
! 
! 		filename_ang_mom='ang_mom_rot_dist1e4_T300_10pct_centrifuge_'//file_number_char//'.dat'
! 		filename_ang_mom=trim(filename_ang_mom)
! 
! 			open(444,file=filename_ang_mom)
! 			do j=1,N
! 				write(444,*) J_rot(j,:)
! 			end do
! 			close(444)
! 
! 		end if

E_trans=0.5*M_mol*(V(:,1)**2+V(:,2)**2+V(:,3)**2)
Omega(:,1)=e(:,2)*u(:,3)-e(:,3)*u(:,2)
Omega(:,2)=e(:,3)*u(:,1)-e(:,1)*u(:,3)
Omega(:,3)=e(:,1)*u(:,2)-e(:,2)*u(:,1)
E_rot=0.5*I_mol*(Omega(:,1)**2+Omega(:,2)**2+Omega(:,3)**2)
!E_pot=pot_energy_tot(R,e)/N
write(2,*) sum(E_trans)/N,sum(E_rot)/N

! x2=e(:,1)**2
! y2=e(:,2)**2
! z2=e(:,3)**2
! xy=e(:,1)*e(:,2)
! write(33,*) sum(x2)/N,sum(y2)/N,sum(z2)/N,sum(xy)/N


!  correl_jk=0.
! do j=1,N-1
!    do k=j+1,N
! 	xj=e(j,1)
! 	xk=e(k,1)
! 	yj=e(j,2)
! 	yk=e(k,2)
! 	correl_jk=correl_jk+(xj*xk+yj*yk)/(sqrt((xj**2+yj**2)*(xk**2+yk**2)))
!    end do
! end do
!  correl_jk=2*correl_jk/(N*(N-1))
! write(44,*) correl_jk

!P_total=M_mol*sum(V,1)/N
! Jbig(:,1)=R(:,2)*V(:,3)-R(:,3)*V(:,2)
! Jbig(:,2)=R(:,3)*V(:,1)-R(:,1)*V(:,3)
! Jbig(:,3)=R(:,1)*V(:,2)-R(:,2)*V(:,1)
! J_trans=M_mol*Jbig
! J_rot=I_mol*Omega
! write(3,*) sum(J_trans,1)/N
! write(4,*) sum(J_rot,1)/N
! write(55,*) sum(J_rot**2,1)/N

! BELOW WE CALCULATE QUANTITIES AVERAGED ONLY OVER THE CENTRIFUGED MOLECULES

E_trans_cent=0.
E_rot_cent=0.
disp2_cent=0.
disp4_cent=0.
! x2_cent=0.
! y2_cent=0.
! z2_cent=0.
! xy_cent=0.
! J_trans_cent=0.
! J_rot_cent=0.
! J_rot_sq_cent=0.
do j=1,N

   if (places_cent(j)) then

	E_trans_cent=E_trans_cent+E_trans(j)
	E_rot_cent=E_rot_cent+E_rot(j)
	disp2_cent=disp2_cent+displacement(j,:)**2
	disp4_cent=disp4_cent+displacement(j,:)**4
! 	x2_cent=x2_cent+x2(j)
! 	y2_cent=y2_cent+y2(j)
! 	z2_cent=z2_cent+z2(j)
! 	xy_cent=xy_cent+xy(j)
! 	J_trans_cent=J_trans_cent+J_trans(j,:)
! 	J_rot_cent=J_rot_cent+J_rot(j,:)
! 	J_rot_sq_cent=J_rot_sq_cent+J_rot(j,:)**2

   end if

end do
write(2,*) E_trans_cent/N_cent,E_rot_cent/N_cent
write(12,*) disp2_cent/N_cent
write(13,*) sqrt(disp4_cent/N_cent-(disp2_cent/N_cent)**2)
!write(33,*) x2_cent/N_cent,y2_cent/N_cent,z2_cent/N_cent,xy_cent/N_cent
! write(3,*) J_trans_cent/N_cent
! write(4,*) J_rot_cent/N_cent
! write(55,*) J_rot_sq_cent/N_cent


!    	V_ini(:,3)=merge(e(:,3)**2,zeros(:,3),places1)
!    	write(8,*) sum(V_ini(:,3))/N_places 


   	!if (i==30000) then
   	!    pol=(/0.,0.,1./) !Polarization direction of the pulse
   	!    P_kick=10. !Kick strength of the pulse
   	!    u=laser_pulse(u,e,places3,pol,P_kick)
   	!end if

   	!if (i==31123) then
   	!    pol=(/1./sqrt(2.),0.,1./sqrt(2.)/) !Polarization direction of the pulse
   	!    P_kick=20. !Kick strength of the pulse
   	!    u=laser_pulse(u,e,places3,pol,P_kick)
   	!end if
   end if

   DisplaceList=DisplaceList+deltaR !Check the total displacement of the molecules sinse the last update of the neighbors' list
   ListUpdateRequested=MovedTooMuch(skin) !if a molecule exists for which the total displacement is higher than the "skin" parameter, "ListUpdateRequested" becomes .true. This is checked in the "MovedTooMuch" function

end do

! HERE ENDS THE MAIN LOOP OF THE TIME PROPAGATION

!close(1)
close(2)
!close(3)
!close(4)
!close(7)
!close(8)
!close(9)
!close(10)
!close(11)
!close(33)
!close(44)
!close(55)
close(12)
close(13)

!Writing the final conditions to a file, in order to continue the simulation from the final point, if needed
open(600,file='final_cond1e4_T300_10pct_centrifuge_check.dat')
write(600,4) p,Temp,N
4 format('p= ',f6.2,'Temp= ',f6.2,'N= ',i5)
do j=1,N
   write(600,*) R(j,:)
end do
write(600,*)
do j=1,N
   write(600,*) V(j,:)
end do
write(600,*)
do j=1,N
   write(600,*) e(j,:)
end do
write(600,*)
do j=1,N
   write(600,*) u(j,:)
end do
write(600,*)
do j=1,N
   write(600,*) u_half(j,:)
end do
do j=1,N
   write(600,*) displacement(j,:)
end do
close(600)

! Recording the final time and finding the total runnung time (works well only when the program is run during a single day):
call DATE_AND_TIME(date,time2)
read(time2(1:2),1) hour2
read(time2(3:4),1) minut2
read(time2(5:10),2) sec2
if (sec2<sec1) then
   minut2=minut2-1
   sec2=sec2+60.
   tot_sec=sec2-sec1
else
   tot_sec=sec2-sec1
end if
if (minut2<minut1) then
   hour2=hour2-1
   minut2=minut2+60
   tot_minut=minut2-minut1
else
   tot_minut=minut2-minut1
end if
tot_hour=hour2-hour1
write(6,3) tot_hour,tot_minut,tot_sec
3 format('This program took ',i3,' hours ',i2,' minutes and ',f6.3,' seconds.')

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
contains

!This subroutine generates the initial positions, orientations, translational and rotational
!velocities for the diatomic dumbbell molecules, that participate in the molecular dynamics
!simulation. The molecules are considered as an ideal gas with respect to
!the definition of their temperature, pressure, and density. The subroutine receives as an input !the gas pressure, temperature and the particle number. It calculates the lengths of the 3-D box, and the positions and velocities of the particles !(center of mass and internal) according to the temperature. It also chooses the molecules that are being centrifuged and gives them orientations and rotational velocities.
subroutine initial_cond(p,Temp,small_nx,small_ny,small_nz,&
                         R_ini,e_ini,V_ini,u_ini,Lx,Ly,Lz,N,places_cent,N_cent)
implicit none
real(kind=8), intent(in) :: p, Temp !pressure in atmospheres, temperature in Kelvin
integer, intent(in) :: small_nx,small_ny,small_nz !molecule number along each one of cube edges
real(kind=8), allocatable, intent(out) :: R_ini(:,:),e_ini(:,:),V_ini(:,:),u_ini(:,:) !initial positions !(R for center of mass, e for molecule orientation) and velocities (V for c.m., u for internal), where the angular velocity Omega=e cross u
real(kind=8), intent(out) :: Lx,Ly,Lz !lengths of the box edges
integer, intent(out) :: N,N_cent !number of molecules, number of centrifuged molecules
real(kind=8) :: small_l, sigma_V, sigma_u,theta_i,phi_i,vtheta_i,vphi_i,x,percent,v_by_V ! small_l - mean distance between molecules along any !dimension, sigma_V - translational temperature, sigma_u - rotational temperature
integer :: small_n,row, nx, ny, nz, i, x_int,j 
real(kind=8), allocatable :: n_mat(:,:), theta(:), phi(:), vtheta(:), vphi(:),u_length(:)
logical, allocatable, intent(out) :: places_cent(:) !which molecules on the list are centrifuged
real(kind=8), parameter :: pi=3.14159265358979

small_l=31.32*((1./p)*(Temp/300.))**(1./3.) !calculate mean distance between molecules from Temp and pressure, according to the ideal gas law

N=(small_nx+1)*(small_ny+1)*(small_nz+1) !total number of molecules

Lx=(small_nx+1)*small_l
Ly=(small_ny+1)*small_l
Lz=(small_nz+1)*small_l

write(*,*) 'N=',N,'l=',small_l
write(*,*) 'Lx=',Lx 
write(*,*) 'Ly=',Ly
write(*,*) 'Lz=',Lz

allocate(n_mat(N,3),R_ini(N,3),e_ini(N,3),V_ini(N,3),u_ini(N,3),theta(N),phi(N),vtheta(N),vphi(N), &
           u_length(N),places_cent(N))

!Here n_mat is created, the cubic lattice of the positions of the center of mass. The lattice parameter is 1.
row=0 
do nx=-small_nx,small_nx,2
    do ny=-small_ny,small_ny,2
        do nz=-small_nz,small_nz,2
            row=row+1
            n_mat(row,:)=(/ nx/2., ny/2., nz/2. /)
        end do
    end do
end do

call random_number(R_ini)
R_ini=small_l*n_mat+0.2*small_l*(2.*R_ini-1.) !each molecule is located at a cubic lattice edge !(the lattice parameter becomes small_l here) and is moved a bit by some random 3D vector with maximum length (in each dimension) of 0.2 times small_l

! call random_number(R_ini)
! R_ini=2.*R_ini-1.
! R_ini(:,1)=R_ini(:,1)*Lx/2
! R_ini(:,2)=R_ini(:,2)*Ly/2
! R_ini(:,3)=R_ini(:,3)*Lz/2

! open(600,file='check.dat')
! do j=1,N
!    write(600,*) R_ini(j,:)
! end do
! close(600)
! stop

call random_number(theta)
theta=2.*asin(sqrt(theta)) 

call random_number(phi)
phi=2.*pi*phi ! theta and phi of the molecular orientation are distributed such that the molecules are !randomly oriented
e_ini(:,1)=sin(theta)*cos(phi)
e_ini(:,2)=sin(theta)*sin(phi)
e_ini(:,3)=cos(theta)

sigma_V=3.528*sqrt(Temp/300.*1./M_mol)
call random_number(n_mat)
call random_number(u_ini)
V_ini=sigma_V*sqrt(-2.*log(n_mat))*cos(2.*pi*u_ini) !c.m. velocities are according to normal !distribution (Box Muller method)

sigma_u=2.*sigma_V/re_mol !For a definite temperature sigma_u depends on sigma_V and the length of the rotor
call random_number(n_mat)
vtheta=sigma_u*sqrt(-2.*log(n_mat(:,1)))*cos(2.*pi*n_mat(:,2)) !
vphi=sigma_u*sqrt(-2.*log(n_mat(:,1)))*sin(2.*pi*n_mat(:,2)) !
!vtheta and vphi are according to normal distribution
u_ini(:,1)=vtheta*cos(theta)*cos(phi)-vphi*sin(phi)
u_ini(:,2)=vtheta*cos(theta)*sin(phi)+vphi*cos(phi)
u_ini(:,3)=-vtheta*sin(theta)

!u_length=sqrt(u_ini(:,1)**2+u_ini(:,2)**2+u_ini(:,3)**2)
!places_cent=u_length<0.5*sigma_u !These velocities constitute 10% of the slowest molecules
percent=0.1
N_cent=nint(percent*N) !The number of centrifuged molecules is some percentage of the total number of molecules
places_cent=.false.
write(*,*) 'N_cent_initial=',N_cent

! do i=1,N_cent
! 
! 	call random_number(x)
! 	x_int=nint(N*x)
! 	places_cent(x_int)=.true.
! 
! end do
do i=1,N

	call random_number(x)
	if (x<percent) places_cent(i)=.true.

end do
!Below the real number of centrifuged molecules is calculated. This number is a little different, from the number defined above, for my random algorithm
N_cent=0
do i=1,N
	if(places_cent(i)==.true.) N_cent=N_cent+1
end do
write(*,*) 'N_cent=',N_cent

!Modifying initial conditions for the centrifuged molecules
v_by_V=6. !The ratio of the rotational and the (average) translational speed for the centrifuged molecules
sigma_u=v_by_V*sqrt(8./pi)*2.*sigma_V/re_mol
theta_i=pi/2.
phi_i=0.
vtheta_i=0.
vphi_i=sigma_u

do i=1,N

   if (places_cent(i)) then
	call random_number(phi_i)
	phi_i=2*pi*phi_i
	e_ini(i,1)=sin(theta_i)*cos(phi_i)
	e_ini(i,2)=sin(theta_i)*sin(phi_i)
	e_ini(i,3)=cos(theta_i)

	u_ini(i,1)=vtheta_i*cos(theta_i)*cos(phi_i)-vphi_i*sin(phi_i)
	u_ini(i,2)=vtheta_i*cos(theta_i)*sin(phi_i)+vphi_i*cos(phi_i)
	u_ini(i,3)=-vtheta_i*sin(theta_i)

   end if

end do

write(*,*) 'centrifuged ratio=',real(N_cent)/real(N)

deallocate(n_mat,theta,phi,vtheta,vphi)

end subroutine initial_cond

!-------------------------------------------------------------------------------------------------

!This subroutine calculates the forces and the torques acting on the
!molecules in a diatomic gas interacting via a Lennard-Jones potential
subroutine forces_torques(M,E,P,Q)
implicit none
real(kind=8), intent(in) :: M(N,3), E(N,3) !Input c.m. positions M and orientations E
real(kind=8), intent(inout) :: P(N,3), Q(N,3) !Output forces P and torques Q
integer :: i,j,ListIndex
real(kind=8) :: F1a(N,3), F1b(N,3), F2a(N,3), F2b(N,3)
real(kind=8) :: Rij(3),Rij1a(3),Rij1b(3),Rij2a(3),Rij2b(3),Force(3)
real(kind=8) :: Rsqij1a,Rsqij1b,Rsqij2a,Rsqij2b,rm2,rm6,rm12

F1a=0. !Force on atom 1 of molecule i from atom a of molecule j
F1b=0. !Force on atom 1 of molecule i from atom b of molecule j
F2a=0. !Force on atom 2 of molecule i from atom a of molecule j
F2b=0. !Force on atom 2 of molecule i from atom b of molecule j

do i=1,N !sum over all i molecules
   do ListIndex=Marker1(i),Marker2(i) !sum over neighbors of each molecule i
        j=List(ListIndex) !The "List" vector contains the lists of neighbors for all the molecules 1 through N
!For molecule "i" the indices of the neighbors "sit" in vector "List" from place number "Marker1(i)" to place number "Marker2(i)"
	Rij=M(i,:)-M(j,:)
!The following makes sure that the calculation of force is done with the closest "image" of the molecule
	if (abs(Rij(1))>Lhalf_x) then
		Rij(1)=Rij(1)-sign(Lx,Rij(1))
	end if
	if (abs(Rij(2))>Lhalf_y) then
		Rij(2)=Rij(2)-sign(Ly,Rij(2))
	end if
	if (abs(Rij(3))>Lhalf_z) then
		Rij(3)=Rij(3)-sign(Lz,Rij(3))
	end if

!Finding the positions of the 4 atoms of the 2 interacting molecules
	Rij1a=Rij+0.5*re_mol*(E(i,:)-E(j,:))
	Rij1b=Rij+0.5*re_mol*(E(i,:)+E(j,:))
	Rij2a=Rij-0.5*re_mol*(E(i,:)+E(j,:))
	Rij2b=Rij-0.5*re_mol*(E(i,:)-E(j,:))

	Rsqij1a=dot_product(Rij1a,Rij1a)
	Rsqij1b=dot_product(Rij1b,Rij1b)
	Rsqij2a=dot_product(Rij2a,Rij2a)
	Rsqij2b=dot_product(Rij2b,Rij2b)

	if (Rsqij1a<cutoffsq) then
		rm2=sigmasq/Rsqij1a
		rm6=rm2**3
		rm12=rm6**2
		Force=const*(rm12-0.5*rm6)*rm2*Rij1a !Force calculation
		F1a(i,:)=F1a(i,:)+Force
		F1a(j,:)=F1a(j,:)-Force !Using Newton's 3rd law (the fact that Fij=-Fji)
	end if
	if (Rsqij1b<cutoffsq) then
		rm2=sigmasq/Rsqij1b
		rm6=rm2**3
		rm12=rm6**2
		Force=const*(rm12-0.5*rm6)*rm2*Rij1b
		F1b(i,:)=F1b(i,:)+Force
		F2a(j,:)=F2a(j,:)-Force
	end if
	if (Rsqij2a<cutoffsq) then
		rm2=sigmasq/Rsqij2a
		rm6=rm2**3
		rm12=rm6**2
		Force=const*(rm12-0.5*rm6)*rm2*Rij2a
		F2a(i,:)=F2a(i,:)+Force
		F1b(j,:)=F1b(j,:)-Force
	end if
	if (Rsqij2b<cutoffsq) then
		rm2=sigmasq/Rsqij2b
		rm6=rm2**3
		rm12=rm6**2
		Force=const*(rm12-0.5*rm6)*rm2*Rij2b
		F2b(i,:)=F2b(i,:)+Force
		F2b(j,:)=F2b(j,:)-Force
	end if

   end do
end do


P=F1a+F1b+F2a+F2b !the total force on the molecule
Q=0.5*re_mol*(F1a+F1b-F2a-F2b) !The vector G. The torque is e cross G
F1a=spread(Q(:,1)*E(:,1)+Q(:,2)*E(:,2)+Q(:,3)*E(:,3),2,3)
Q=Q-F1a*E !the torque on the molecule, after substracting the part acting along the bond, or G_perpend

end subroutine forces_torques


!-------------------------------------------------------------------------------------------------
!This function checks if there exist two molecules that the sum of their total displacements is higher than the "skin" parameter
function MovedTooMuch(skin)
implicit none
real(kind=8), intent(in) :: skin
logical :: MovedTooMuch
real(kind=8) :: Displ1,Displ2,Displ
integer :: i

Displ1=0.!Largest displacement
Displ2=0.!Second largest displacement
do i=1,N
	Displ=dot_product(DisplaceList(i,:),DisplaceList(i,:))
	if (Displ >= Displ1) then
		Displ2=Displ1
		Displ1=Displ
 	else if (Displ >= Displ2) then
		Displ2=Displ
	end if
end do
MovedTooMuch=(Displ1+Displ2+2.*sqrt(Displ1)*sqrt(Displ2)>skin**2)

end function MovedTooMuch
!-----------------------------------------------------------------------------------------------------
!This subroutine finds (when required) the updated neighbors list for all the molecules. This is a heavy calculation with a running time proportional to N(N-1)/2. "M" is the matrix of the molecular c.m. positions. "dist" is the radius of the neighbors search, equal to the cutoff radius+skin.
subroutine UpdateList(dist,M)
implicit none
real(kind=8),intent(in) :: dist, M(N,3)
real(kind=8) :: distsq,Rij(3),Rsqij
integer :: i,j,ListIndex,advance(N)

distsq=dist**2

ListIndex=1
do i=1,N
   do j=i+1,N !No need to run over all pairs, because Fij=-Fji
	Rij=M(i,:)-M(j,:)

	if (abs(Rij(1))>Lhalf_x) then
		Rij(1)=Rij(1)-sign(Lx,Rij(1))
	end if
	if (abs(Rij(2))>Lhalf_y) then
		Rij(2)=Rij(2)-sign(Ly,Rij(2))
	end if
	if (abs(Rij(3))>Lhalf_z) then
		Rij(3)=Rij(3)-sign(Lz,Rij(3))
	end if

	Rsqij=dot_product(Rij,Rij)
	if (Rsqij<distsq) then
		advance(j)=1 !j is a neighbor of i
	else
		advance(j)=0 !j is not a neighbor of i
	end if
   end do
   Marker1(i)=ListIndex !this is where the list of neighbors for atom i begins
   do j=i+1,N
	if (ListIndex>N*MaxPairsPerAtom) goto 9999
	List(ListIndex)=j !write the number of neighbor in "List"
	ListIndex=ListIndex+advance(j) !add "1" to LIstIndex if neighbor is added
   end do
   Marker2(i)=ListIndex-1 !this is where the list of neighbors for atom i ends
end do
!write(*,*) 'List length is ',ListIndex-1

DisplaceList=0. !Reset the displacements after the neighbors list was updated
return

9999 continue
	print '(a)','Dimensions of the variable List are too small.'
	print '(a)','You can decrease the Skin parameter,'
	print '(a)','or increase the variable MaxPairsPerAtom.'
	stop

end subroutine UpdateList
!--------------------------------------------------------------------------------------------------
!This function calculates the total potential energy of the diatomic
!molecular gas, for the Lennard-Jones interaction potential. Simiar to force calculation, except the result is a scalar - the total potential energy
function pot_energy_tot(M,E)
implicit none
real(kind=8), intent(in) :: M(N,3), E(N,3)
real(kind=8) :: pot_energy_tot
real(kind=8) :: Rij(3),Rij1a(3),Rij1b(3),Rij2a(3),Rij2b(3)
real(kind=8) :: Rsqij1a,Rsqij1b,Rsqij2a,Rsqij2b,rm2,rm6,rm12,potential
integer :: i,j,ListIndex

pot_energy_tot=0.

do i=1,N
   do ListIndex=Marker1(i),Marker2(i)
        j=List(ListIndex)
	Rij=M(i,:)-M(j,:)

	if (abs(Rij(1))>Lhalf_x) then
		Rij(1)=Rij(1)-sign(Lx,Rij(1))
	end if
	if (abs(Rij(2))>Lhalf_y) then
		Rij(2)=Rij(2)-sign(Ly,Rij(2))
	end if
	if (abs(Rij(3))>Lhalf_z) then
		Rij(3)=Rij(3)-sign(Lz,Rij(3))
	end if

	Rij1a=Rij+0.5*re_mol*(E(i,:)-E(j,:))
	Rij1b=Rij+0.5*re_mol*(E(i,:)+E(j,:))
	Rij2a=Rij-0.5*re_mol*(E(i,:)+E(j,:))
	Rij2b=Rij-0.5*re_mol*(E(i,:)-E(j,:))

	Rsqij1a=dot_product(Rij1a,Rij1a)
	Rsqij1b=dot_product(Rij1b,Rij1b)
	Rsqij2a=dot_product(Rij2a,Rij2a)
	Rsqij2b=dot_product(Rij2b,Rij2b)

	if (Rsqij1a<cutoffsq) then
		rm2=sigmasq/Rsqij1a
		rm6=rm2**3
		rm12=rm6**2
		potential=const2*(rm12-rm6)-pot_at_cutoff
		pot_energy_tot=pot_energy_tot+potential
	end if
	if (Rsqij1b<cutoffsq) then
		rm2=sigmasq/Rsqij1b
		rm6=rm2**3
		rm12=rm6**2
		potential=const2*(rm12-rm6)-pot_at_cutoff
		pot_energy_tot=pot_energy_tot+potential
	end if
	if (Rsqij2a<cutoffsq) then
		rm2=sigmasq/Rsqij2a
		rm6=rm2**3
		rm12=rm6**2
		potential=const2*(rm12-rm6)-pot_at_cutoff
		pot_energy_tot=pot_energy_tot+potential
	end if
	if (Rsqij2b<cutoffsq) then
		rm2=sigmasq/Rsqij2b
		rm6=rm2**3
		rm12=rm6**2
		potential=const2*(rm12-rm6)-pot_at_cutoff
		pot_energy_tot=pot_energy_tot+potential
	end if

   end do
end do


end function pot_energy_tot


!-------------------------------------------------------------------------------------------------

! This function "kicks" the molecules in places defined by the variable "places3" by an ultrashort laser pulse. The pulse is characterized by polarization vector "pol" and a kick strength "P". The function gives back to the program the new value of the rotational velocity matrix "u".
function laser_pulse(u,e,places3,pol,P) 
implicit none
real(kind=8), intent(in) :: u(N,3), e(N,3), pol(3), P
logical, intent(in) :: places3(N,3)
real(kind=8) :: laser_pulse(N,3)
real(kind=8) :: pol_mat(N,3), cos_beta(N), cos_beta_mat(N,3), delta_u(N,3)

pol_mat=spread(pol,1,N)
 cos_beta=pol(1)*e(:,1)+pol(2)*e(:,2)+pol(3)*e(:,3) !cosine of the angle between the molecular direction and the pulse polarization direction
 cos_beta_mat=spread(cos_beta,2,3)
delta_u=2*P*cos_beta_mat*(pol_mat-cos_beta_mat*e)/I_mol
laser_pulse=merge(u+delta_u,u,places3)

end function laser_pulse

!--------------------------------------------------------------------------------------------------------
!This function replaces the translational velocities of molecules at "places3" by a thermal distribution with a temperature "T_heat"
function local_heating(V,places3,T_heat)
implicit none
real(kind=8), intent(in) :: V(N,3),T_heat
logical, intent(in) :: places3(N,3)
real(kind=8) :: local_heating(N,3)
real(kind=8) :: sigma_V,n_mat(N,3),m_mat(N,3),V_heat(N,3)
real(kind=8), parameter :: pi=3.14159265358979

sigma_V=3.528*sqrt(T_heat/300.*1./M_mol)
call random_number(n_mat)
call random_number(m_mat)
V_heat=sigma_V*sqrt(-2.*log(n_mat))*cos(2.*pi*m_mat) !c.m. velocities are according to normal !distribution

local_heating=merge(V_heat,V,places3)

end function local_heating

!--------------------------------------------------------------------------------------------------------

end program molecular_dynamics
