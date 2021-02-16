c  ***********************  EMC3D.f  ***********************************
c  This program is a modified version of "elas3d.f", which was developed
c  and distributed by Garboczi and others from NIST.
c
c  Adding fuctionalities are not major.
c   (1) File IO (input/output/stress field output files)
c   (2) Calculation of properties assuming isotropic & elastic media.
c
c  Following is the original documentation from NIST.
c
c-----------------------
c  BACKGROUND

c  This program solves the linear elastic equations in a
c  random linear elastic material, subject to an applied macroscopic strain,
c  using the finite element method.  Each pixel in the 3-D digital
c  image is a cubic tri-linear finite element,  having its own
c  elastic moduli tensor. Periodic boundary conditions are maintained.
c  In the comments below, (USER) means that this is a section of code that
c  the user might have to change for his particular problem. Therefore the
c  user is encouraged to search for this string.

c  PROBLEM AND VARIABLE DEFINITION

c  The problem being solved is the minimization of the energy
c  1/2 uAu + b u + C, where A is the Hessian matrix composed of the
c  stiffness matrices (dk) for each pixel/element, b is a constant vector
c  and C is a constant that are determined by the applied strain and
c  the periodic boundary conditions, and u is a vector of
c  all the displacements. The solution
c  method used is the conjugate gradient relaxation algorithm.
c  Other variables are:  gb is the gradient = Au+b, h and Ah are
c  auxiliary variables used in the conjugate gradient algorithm (in dembx),
c  dk(n,i,j) is the stiffness matrix of the n'th phase, cmod(n,i,j) is
c  the elastic moduli tensor of the n'th phase, pix is a vector that gives
c  the phase label of each pixel, ib is a matrix that gives the labels of
c  the 27 (counting itself) neighbors of a given node, prob is the volume
c  fractions of the various phases,
c  strxx, stryy, strzz, strxz, stryz, and strxy are the six Voigt
c  volume averaged total stresses, and
c  sxx, syy, szz, sxz, syz, and sxy are the six Voigt
c  volume averaged total strains.

c  DIMENSIONS

c  The vectors u,gb,b,h, and Ah are dimensioned to be the system size,
c  ns=nx*ny*nz, with three components, where the digital image of the
c  microstructure considered is a rectangular paralleliped, nx x ny x nz
c  in size.  The arrays pix and ib are are also dimensioned to the system size.
c  The array ib has 27 components, for the 27 neighbors of a node.
c  Note that the program is set up at present to have at most 100
c  different phases.  This can easily be changed, simply by changing
c  the dimensions of dk, prob, and cmod. The parameter nphase gives the
c  number of phases being considered in the problem.
c  All arrays are passed between subroutines using simple common statements.

c  STRONGLY SUGGESTED:  READ THE MANUAL BEFORE USING PROGRAM!!

c  (USER) Change these dimensions and in other subroutines at same time.
c  For example, search and replace all occurrences throughout the
c  program of "(64000000" by "(64000", to go from a
c  20 x 20 x 20 system to a 40 x 40 x 40 system.
	real u(64000000,3),gb(64000000,3),b(64000000,3)
      real h(64000000,3),Ah(64000000,3)
	real cmod(100,6,6),dk(100,8,3,8,3)
	real phasemod(100,2),prob(100), density(100)
	real fdensity, fporosity, BulkM, ShearM


	integer in(27),jn(27),kn(27),iwriteFlag
	integer*4 ib(64000000,27), nblank
	integer*2 pix(64000000)
	character*40 poreFileName

	common/list1/strxx,stryy,strzz,strxz,stryz,strxy
	common/list2/exx,eyy,ezz,exz,eyz,exy
	common/list3/ib
        common/list4/pix
	common/list5/dk,b,C
	common/list6/u
	common/list7/gb
        common/list8/cmod
        common/list9/h,Ah
	common/list10/sxx,syy,szz,sxz,syz,sxy

c (USER) Unit 9 is the microstructure input file,
c  unit 7 is the results output file.
	open (unit=8, file='emc3d.pam')

      open (unit=7,file='outputfile.out')


	read(8,'(a40)') poreFileName
	write(*,7901) poreFileName
7901  format(' Pore File Name : ', a40)

	read(8,*) nx, ny, nz
	write(*,7902) nx, ny, nz
7902  format(' (NX, NY, NZ) : ',3i5)

	read(8,*) gtest
	write(*,7903) gtest
7903  format(' Convergence criteria: ',e12.5)

	read(8,*) nblank
	write(*,7904) nblank
7904  format(' No. of Header line : ',i5)

	read(8,*) iwriteFlag
	write(*,7905) iwriteFlag
7905  format(' Write Stress Field (0=no, 1=yes) : ',i3)

	read(8,*) nphase
	write(*,7906) nphase
7906  format(' No. of Phases', i3)

	do 6060 n=1,nphase
		read(8,*) phasemod(n,1), phasemod(n,2), density(n)
		write(*,7907) n, phasemod(n,1), phasemod(n,2), density(n)
6060  continue
7907  format(i2,' phase (K,mu,density) : ', 3(e12.5,1x))

	read(8,*) exx,eyy,ezz,exz,eyz,exy
	print*, 'Applied Strain :'
	write(*,7908) exx,eyy,ezz,exz,eyz,exy
7908  format(3x,6(e12.5,1x))

	close(8)

	open (unit=9,file=poreFileName)

c (USER)  nx,ny,nz give the size of the lattice
c        nx=200
c        ny=200
c        nz=200
c ns=total number of sites
        ns=nx*ny*nz
      write(7,9010) nx,ny,nz,ns
9010  format('nx= ',i4,' ny= ',i4,' nz= ',i4,' ns= ',i8)
      gtest = gtest*ns

c  (USER) nphase is the number of phases being considered in the problem.
c  The values of pix(m) will run from 1 to nphase.
c	nphase=2

c  (USER) gtest is the stopping criterion, the number
c  to which the quantity gg=gb*gb is compared.
c  Usually gtest = abc*ns, so that when gg < gtest, the rms value
c  per pixel of gb is less than sqrt(abc).
c       gtest=1.e-7*ns

c  (USER)
c  The parameter phasemod(i,j) is the bulk (i,1) and shear (i,2) moduli of
c  the i'th phase. These can be input in terms of Young's moduli E(i,1) and
c  Poisson's ratio nu (i,2).  The program, in do loop 1144, then changes them
c  to bulk and shear moduli, using relations for isotropic elastic
c  moduli.  For anisotropic elastic material, one can directly input
c  the elastic moduli tensor cmod in subroutine femat, and skip this part.
c  If you wish to input in terms of bulk (1) and shear (2), then make sure
c  to comment out the do 1144 loop.
c	phasemod(1,1)=2.2
c	phasemod(1,2)=0.0
c	phasemod(2,1)=36.6
c	phasemod(2,2)=44.0

c---------------------------------**********************___________________
c  (USER) Program uses bulk modulus (1) and shear modulus (2), so transform
c  Young's modulis (1) and Poisson's ratio (2).
c        do 1144 i=1,nphase
c        save=phasemod(i,1)
c        phasemod(i,1)=phasemod(i,1)/3./(1.-2.*phasemod(i,2))
c        phasemod(i,2)=save/2./(1.+phasemod(i,2))
c1144    continue

c  Construct the neighbor table, ib(m,n)

c  First construct the 27 neighbor table in terms of delta i, delta j, and
c  delta k information (see Table 3 in manual)
      in(1)=0
      in(2)=1
      in(3)=1
      in(4)=1
      in(5)=0
      in(6)=-1
      in(7)=-1
      in(8)=-1

      jn(1)=1
      jn(2)=1
      jn(3)=0
      jn(4)=-1
      jn(5)=-1
      jn(6)=-1
      jn(7)=0
      jn(8)=1

      do 555 n=1,8
      kn(n)=0
      kn(n+8)=-1
      kn(n+16)=1
      in(n+8)=in(n)
      in(n+16)=in(n)
      jn(n+8)=jn(n)
      jn(n+16)=jn(n)
555   continue
      in(25)=0
      in(26)=0
      in(27)=0
      jn(25)=0
      jn(26)=0
      jn(27)=0
      kn(25)=-1
      kn(26)=1
      kn(27)=0

c  Now construct neighbor table according to 1-d labels
c  Matrix ib(m,n) gives the 1-d label of the n'th neighbor (n=1,27) of
c  the node labelled m.
      nxy=nx*ny
      do 1020 k=1,nz
      do 1020 j=1,ny
      do 1020 i=1,nx
      m=nxy*(k-1)+nx*(j-1)+i
      do 1004 n=1,27
      i1=i+in(n)
      j1=j+jn(n)
      k1=k+kn(n)
      if(i1.lt.1) i1=i1+nx
      if(i1.gt.nx) i1=i1-nx
      if(j1.lt.1) j1=j1+ny
      if(j1.gt.ny) j1=j1-ny
      if(k1.lt.1) k1=k1+nz
      if(k1.gt.nz) k1=k1-nz
      m1=nxy*(k1-1)+nx*(j1-1)+i1
      ib(m,n)=m1
1004  continue
1020  continue

c Compute the average stress and strain in each microstructure.
c (USER) npoints is the number of microstructures to use.

        npoints=1
        do 8000 micro=1,npoints
c  Read in a microstructure in subroutine ppixel, and set up pix(m)
c  with the appropriate phase assignments.
        call ppixel(nx,ny,nz,ns,nphase,nblank)
c Count and output the volume fractions of the different phases
        call assig(ns,nphase,prob)
        do 111 i=1,nphase
        write(7,9020) i,phasemod(i,1),phasemod(i,2)
9020    format(' Phase ',i3,' bulk = ',f12.6,' shear = ',f12.6)
111	continue

	do 8050 i=1,nphase
	write(7,9065) i,prob(i)
9065	format(' Volume fraction of phase ',i3,'  is ',f8.5)
8050	continue
	fporosity = prob(1)

c  (USER) Set applied strains
c  Actual shear strain applied in do 1050 loop is exy, exz, and eyz as
c  given in the statements below.  The engineering shear strain, by which
c  the shear modulus is usually defined, is twice these values.
c        exx=0.002
c        eyy=0.002
c        ezz=0.002
c        exz=0.002/2.
c        eyz=0.004/2.
c        exy=0.006/2.
	exz = exz/2.
	eyz = eyz/2.
	exy	= exy/2.
        write(7,*) 'Applied engineering strains'
        write(7,*) ' exx eyy ezz exz eyz exy'
        write(7,*) exx,eyy,ezz,2.*exz,2.*eyz,2.*exy

c Set up the elastic modulus variables, finite element stiffness matrices,
c the constant, C, and vector, b, required for computing the energy.
c  (USER) If anisotropic elastic moduli tensors are used, these need to be
c  input in subroutine femat.

	call femat(nx,ny,nz,ns,phasemod,nphase)

c Apply chosen strains as a homogeneous macroscopic strain
c as the initial condition.
	do 1050 k=1,nz
	do 1050 j=1,ny
        do 1050 i=1,nx
		m=nxy*(k-1)+nx*(j-1)+i
		x=float(i-1)
		y=float(j-1)
		z=float(k-1)
		u(m,1)=x*exx+y*exy+z*exz
                u(m,2)=x*exy+y*eyy+z*eyz
                u(m,3)=x*exz+y*eyz+z*ezz
1050	continue

c  RELAXATION LOOP
c  (USER) kmax is the maximum number of times dembx will be called, with
c  ldemb conjugate gradient steps performed during each call.  The total
c  number of conjugate gradient steps allowed for a given elastic
c  computation is kmax*ldemb.
        kmax=50
        ldemb=50
        ltot=0
c  Call energy to get initial energy and initial gradient
        call energy(nx,ny,nz,ns,utot)
c  gg is the norm squared of the gradient (gg=gb*gb)
 	gg=0.0
        do 100 m3=1,3
        do 100 m=1,ns
        gg=gg+gb(m,m3)*gb(m,m3)
100     continue
	write(7,*) 'Initial energy = ',utot,' gg = ',gg
c        call flush(7)

        do 5000 kkk=1,kmax
c  call dembx to go into the conjugate gradient solver
        call dembx(ns,Lstep,gg,dk,gtest,ldemb,kkk)
        ltot=ltot+Lstep
c  Call energy to compute energy after dembx call. If gg < gtest, this
c  will be the final energy.  If gg is still larger than gtest, then this
c  will give an intermediate energy with which to check how the
c  relaxation process is coming along.
        call energy(nx,ny,nz,ns,utot)
	write(7,*) 'Energy = ',utot,' gg = ',gg
	write(7,*) 'Number of conjugate steps  = ',ltot
c       call flush(7)
c  If relaxation process is finished, jump out of loop
        if(gg.le.gtest) goto 444
c  If relaxation process will continue, compute and output stresses
c  and strains as an additional aid to judge how the
c  relaxation procedure is progressing.
	call stress(nx,ny,nz,ns,0)
        write(7,*) ' stresses:  xx,yy,zz,xz,yz,xy'
	write(7,*) strxx,stryy,strzz,strxz,stryz,strxy
        write(7,*) ' strains:  xx,yy,zz,xz,yz,xy'
	write(7,*) sxx,syy,szz,sxz,syz,sxy
c        call flush(7)
5000    continue

444     call stress(nx,ny,nz,ns,iwriteFlag)

	write(7,*)
	write(7,*) ' stresses:  xx,yy,zz,xz,yz,xy'
	write(7,6072) strxx,stryy,strzz,strxz,stryz,strxy
      write(7,*) ' strains:  xx,yy,zz,xz,yz,xy'
	write(7,6072) sxx,syy,szz,sxz,syz,sxy
	write(7,*)

6072	format(6(e12.5, 2x))

	BulkM  = (strxx+stryy+strzz)/(sxx+syy+szz)/3.0
	ShearM = (strxz/sxz + stryz/syz + strxy/sxy)/3.0
	fdensity = 0.0;
	do 6061 n=1,nphase
		fdensity = fdensity + density(n)*prob(n)
6061	continue

6071  format(a10,f12.5)
	write(7,6071) 'K       : ', BulkM
	write(7,6071) 'G       : ', ShearM
	write(7,6071) 'density : ', fdensity
	write(7,6071) 'Porosity: ', fporosity
	write(7,6071) 'Vp      : ', sqrt((BulkM+4.0/3.0*ShearM)/fdensity)
	write(7,6071) 'Vs      : ', sqrt(ShearM/fdensity)
	write(7,6071) 'E       : ', 9.0*BulkM*ShearM/(3.0*BulkM+ShearM)
	write(7,6071) 'Poisson : ',(3.*BulkM-2.*ShearM)/
     &                           (6.*BulkM+2.*ShearM)

c	open(unit = 8, file = 'tResult.txt', status = 'unknown')
c	write(8, 6969) fporosity, BulkM, ShearM,
c     &   sqrt((BulkM+4.0/3.0*ShearM)/fdensity),sqrt(ShearM/fdensity)
c6969  format(5(e12.5, 2x))
c	close(8)

8000    continue

        end

c  Subroutine that sets up the elastic moduli variables,
c  the stiffness matrices,dk, the linear term in
c  displacements, b, and the constant term, C, that appear in the total energy
c  due to the periodic boundary conditions

      subroutine femat(nx,ny,nz,ns,phasemod,nphase)
      real dk(100,8,3,8,3),phasemod(100,2),dndx(8),dndy(8),dndz(8)
      real b(64000000,3),g(3,3,3),C,ck(6,6),cmu(6,6),cmod(100,6,6)
      real es(6,8,3),delta(8,3)
      integer is(8)
      integer*4 ib(64000000,27)
      integer*2 pix(64000000)

	common/list2/exx,eyy,ezz,exz,eyz,exy
	common/list3/ib
        common/list4/pix
	common/list5/dk,b,C
        common/list8/cmod

      nxy=nx*ny

c  (USER) NOTE:  complete elastic modulus matrix is used, so an anisotropic
c  matrix could be directly input at any point, since program is written
c  to use a general elastic moduli tensor, but is only explicitly
c  implemented for isotropic materials.

c  initialize stiffness matrices
      do 40 m=1,nphase
      do 40 l=1,3
      do 40 k=1,3
      do 40 j=1,8
      do 40 i=1,8
      dk(m,i,k,j,l)=0.0
40    continue

c set up elastic moduli matrices for each kind of element
c  ck and cmu are the bulk and shear modulus matrices, which need to be
c  weighted by the actual bulk and shear moduli

      ck(1,1)=1.0
      ck(1,2)=1.0
      ck(1,3)=1.0
      ck(1,4)=0.0
      ck(1,5)=0.0
      ck(1,6)=0.0
      ck(2,1)=1.0
      ck(2,2)=1.0
      ck(2,3)=1.0
      ck(2,4)=0.0
      ck(2,5)=0.0
      ck(2,6)=0.0
      ck(3,1)=1.0
      ck(3,2)=1.0
      ck(3,3)=1.0
      ck(3,4)=0.0
      ck(3,5)=0.0
      ck(3,6)=0.0
      ck(4,1)=0.0
      ck(4,2)=0.0
      ck(4,3)=0.0
      ck(4,4)=0.0
      ck(4,5)=0.0
      ck(4,6)=0.0
      ck(5,1)=0.0
      ck(5,2)=0.0
      ck(5,3)=0.0
      ck(5,4)=0.0
      ck(5,5)=0.0
      ck(5,6)=0.0
      ck(6,1)=0.0
      ck(6,2)=0.0
      ck(6,3)=0.0
      ck(6,4)=0.0
      ck(6,5)=0.0
      ck(6,6)=0.0

      cmu(1,1)=4.0/3.0
      cmu(1,2)=-2.0/3.0
      cmu(1,3)=-2.0/3.0
      cmu(1,4)=0.0
      cmu(1,5)=0.0
      cmu(1,6)=0.0
      cmu(2,1)=-2.0/3.0
      cmu(2,2)=4.0/3.0
      cmu(2,3)=-2.0/3.0
      cmu(2,4)=0.0
      cmu(2,5)=0.0
      cmu(2,6)=0.0
      cmu(3,1)=-2.0/3.0
      cmu(3,2)=-2.0/3.0
      cmu(3,3)=4.0/3.0
      cmu(3,4)=0.0
      cmu(3,5)=0.0
      cmu(3,6)=0.0
      cmu(4,1)=0.0
      cmu(4,2)=0.0
      cmu(4,3)=0.0
      cmu(4,4)=1.0
      cmu(4,5)=0.0
      cmu(4,6)=0.0
      cmu(5,1)=0.0
      cmu(5,2)=0.0
      cmu(5,3)=0.0
      cmu(5,4)=0.0
      cmu(5,5)=1.0
      cmu(5,6)=0.0
      cmu(6,1)=0.0
      cmu(6,2)=0.0
      cmu(6,3)=0.0
      cmu(6,4)=0.0
      cmu(6,5)=0.0
      cmu(6,6)=1.0

      do 31 k=1,nphase
      do 21 j=1,6
      do 11 i=1,6
      cmod(k,i,j)=phasemod(k,1)*ck(i,j)+phasemod(k,2)*cmu(i,j)
11    continue
21    continue
31    continue

c  set up Simpson's integration rule weight vector
      do 30 k=1,3
      do 30 j=1,3
      do 30 i=1,3
      nm=0
      if(i.eq.2) nm=nm+1
      if(j.eq.2) nm=nm+1
      if(k.eq.2) nm=nm+1
      g(i,j,k)=4.0**nm
30    continue

c  loop over the nphase kinds of pixels and Simpson's rule quadrature
c  points in order to compute the stiffness matrices.  Stiffness matrices
c  of trilinear finite elements are quadratic in x, y, and z, so that
c  Simpson's rule quadrature gives exact results.
      do 4000 ijk=1,nphase
      do 3000 k=1,3
      do 3000 j=1,3
      do 3000 i=1,3
      x=float(i-1)/2.0
      y=float(j-1)/2.0
      z=float(k-1)/2.0
c  dndx means the negative derivative, with respect to x, of the shape
c  matrix N (see manual, Sec. 2.2), dndy, and dndz are similar.
      dndx(1)=-(1.0-y)*(1.0-z)
      dndx(2)=(1.0-y)*(1.0-z)
      dndx(3)=y*(1.0-z)
      dndx(4)=-y*(1.0-z)
      dndx(5)=-(1.0-y)*z
      dndx(6)=(1.0-y)*z
      dndx(7)=y*z
      dndx(8)=-y*z
      dndy(1)=-(1.0-x)*(1.0-z)
      dndy(2)=-x*(1.0-z)
      dndy(3)=x*(1.0-z)
      dndy(4)=(1.0-x)*(1.0-z)
      dndy(5)=-(1.0-x)*z
      dndy(6)=-x*z
      dndy(7)=x*z
      dndy(8)=(1.0-x)*z
      dndz(1)=-(1.0-x)*(1.0-y)
      dndz(2)=-x*(1.0-y)
      dndz(3)=-x*y
      dndz(4)=-(1.0-x)*y
      dndz(5)=(1.0-x)*(1.0-y)
      dndz(6)=x*(1.0-y)
      dndz(7)=x*y
      dndz(8)=(1.0-x)*y
c  now build strain matrix
      do 2799 n1=1,6
      do 2799 n2=1,8
      do 2799 n3=1,3
      es(n1,n2,n3)=0.0
2799  continue
      do 2797 n=1,8
      es(1,n,1)=dndx(n)
      es(2,n,2)=dndy(n)
      es(3,n,3)=dndz(n)
      es(4,n,1)=dndz(n)
      es(4,n,3)=dndx(n)
      es(5,n,2)=dndz(n)
      es(5,n,3)=dndy(n)
      es(6,n,1)=dndy(n)
      es(6,n,2)=dndx(n)
2797  continue
c  Matrix multiply to determine value at (x,y,z), multiply by
c  proper weight, and sum into dk, the stiffness matrix
      do 900 mm=1,3
      do 900 nn=1,3
      do 900 ii=1,8
      do 900 jj=1,8
c  Define sum over strain matrices and elastic moduli matrix for
c  stiffness matrix
      sum=0.0
      do 890 kk=1,6
      do 890 ll=1,6
      sum=sum+es(kk,ii,mm)*cmod(ijk,kk,ll)*es(ll,jj,nn)
890   continue
      dk(ijk,ii,mm,jj,nn)=dk(ijk,ii,mm,jj,nn)+g(i,j,k)*sum/216.
900   continue
3000  continue
4000  continue

c  Set up vector for linear term, b, and constant term, C,
c  in the elastic energy.  This is done using the stiffness matrices,
c  and the periodic terms in the applied strain that come in at the
c  boundary pixels via the periodic boundary conditions and the
c  condition that an applied macroscopic strain exists (see Sec. 2.2
c  in the manual). It is easier to set b up this way than to analytically
c  write out all the terms involved.

c  Initialize b and C
      do 5000 m3=1,3
      do 5000 m=1,ns
      b(m,m3)=0.0
5000  continue
      C=0.0

c  For all cases, the correspondence between 1-8 finite element node
c  labels and 1-27 neighbor labels is (see Table 4 in manual):
c  1:ib(m,27), 2:ib(m,3),
c  3:ib(m,2),4:ib(m,1),
c  5:ib(m,26),6:ib(m,19)
c  7:ib(m,18),8:ib(m,17)
      is(1)=27
      is(2)=3
      is(3)=2
      is(4)=1
      is(5)=26
      is(6)=19
      is(7)=18
      is(8)=17

c  x=nx face
      do 2001 i3=1,3
      do 2001 i8=1,8
      delta(i8,i3)=0.0
      if(i8.eq.2.or.i8.eq.3.or.i8.eq.6.or.i8.eq.7) then
      delta(i8,1)=exx*nx
      delta(i8,2)=exy*nx
      delta(i8,3)=exz*nx
      end if
2001  continue
      do 2000 j=1,ny-1
      do 2000 k=1,nz-1
      m=nxy*(k-1)+j*nx
      do 1900 nn=1,3
      do 1900 mm=1,8
      sum=0.0
      do 1899 m3=1,3
      do 1899 m8=1,8
      sum=sum+delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)
      C=C+0.5*delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)*delta(mm,nn)
1899  continue
      b(ib(m,is(mm)),nn)=b(ib(m,is(mm)),nn)+sum
1900  continue
2000  continue
c  y=ny face
      do 2011 i3=1,3
      do 2011 i8=1,8
      delta(i8,i3)=0.0
      if(i8.eq.3.or.i8.eq.4.or.i8.eq.7.or.i8.eq.8) then
      delta(i8,1)=exy*ny
      delta(i8,2)=eyy*ny
      delta(i8,3)=eyz*ny
      end if
2011  continue
      do 2010 i=1,nx-1
      do 2010 k=1,nz-1
      m=nxy*(k-1)+nx*(ny-1)+i
      do 1901 nn=1,3
      do 1901 mm=1,8
      sum=0.0
      do 2099 m3=1,3
      do 2099 m8=1,8
      sum=sum+delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)
      C=C+0.5*delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)*delta(mm,nn)
2099  continue
      b(ib(m,is(mm)),nn)=b(ib(m,is(mm)),nn)+sum
1901  continue
2010  continue
c  z=nz face
      do 2021 i3=1,3
      do 2021 i8=1,8
      delta(i8,i3)=0.0
      if(i8.eq.5.or.i8.eq.6.or.i8.eq.7.or.i8.eq.8) then
      delta(i8,1)=exz*nz
      delta(i8,2)=eyz*nz
      delta(i8,3)=ezz*nz
      end if
2021  continue
      do 2020 i=1,nx-1
      do 2020 j=1,ny-1
      m=nxy*(nz-1)+nx*(j-1)+i
      do 1902 nn=1,3
      do 1902 mm=1,8
      sum=0.0
      do 2019 m3=1,3
      do 2019 m8=1,8
      sum=sum+delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)
      C=C+0.5*delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)*delta(mm,nn)
2019  continue
      b(ib(m,is(mm)),nn)=b(ib(m,is(mm)),nn)+sum
1902  continue
2020  continue
c  x=nx y=ny edge
      do 2031 i3=1,3
      do 2031 i8=1,8
      delta(i8,i3)=0.0
      if(i8.eq.2.or.i8.eq.6) then
      delta(i8,1)=exx*nx
      delta(i8,2)=exy*nx
      delta(i8,3)=exz*nx
      end if
      if(i8.eq.4.or.i8.eq.8) then
      delta(i8,1)=exy*ny
      delta(i8,2)=eyy*ny
      delta(i8,3)=eyz*ny
      end if
      if(i8.eq.3.or.i8.eq.7) then
      delta(i8,1)=exy*ny+exx*nx
      delta(i8,2)=eyy*ny+exy*nx
      delta(i8,3)=eyz*ny+exz*nx
      end if
2031  continue
      do 2030 k=1,nz-1
      m=nxy*k
      do 1903 nn=1,3
      do 1903 mm=1,8
      sum=0.0
      do 2029 m3=1,3
      do 2029 m8=1,8
      sum=sum+delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)
      C=C+0.5*delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)*delta(mm,nn)
2029  continue
      b(ib(m,is(mm)),nn)=b(ib(m,is(mm)),nn)+sum
1903  continue
2030  continue
c  x=nx z=nz edge
      do 2041 i3=1,3
      do 2041 i8=1,8
      delta(i8,i3)=0.0
      if(i8.eq.2.or.i8.eq.3) then
      delta(i8,1)=exx*nx
      delta(i8,2)=exy*nx
      delta(i8,3)=exz*nx
      end if
      if(i8.eq.5.or.i8.eq.8) then
      delta(i8,1)=exz*nz
      delta(i8,2)=eyz*nz
      delta(i8,3)=ezz*nz
      end if
      if(i8.eq.6.or.i8.eq.7) then
      delta(i8,1)=exz*nz+exx*nx
      delta(i8,2)=eyz*nz+exy*nx
      delta(i8,3)=ezz*nz+exz*nx
      end if
2041  continue
      do 2040 j=1,ny-1
      m=nxy*(nz-1)+nx*(j-1)+nx
      do 1904 nn=1,3
      do 1904 mm=1,8
      sum=0.0
      do 2039 m3=1,3
      do 2039 m8=1,8
      sum=sum+delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)
      C=C+0.5*delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)*delta(mm,nn)
2039  continue
      b(ib(m,is(mm)),nn)=b(ib(m,is(mm)),nn)+sum
1904  continue
2040  continue
c  y=ny z=nz edge
      do 2051 i3=1,3
      do 2051 i8=1,8
      delta(i8,i3)=0.0
      if(i8.eq.5.or.i8.eq.6) then
      delta(i8,1)=exz*nz
      delta(i8,2)=eyz*nz
      delta(i8,3)=ezz*nz
      end if
      if(i8.eq.3.or.i8.eq.4) then
      delta(i8,1)=exy*ny
      delta(i8,2)=eyy*ny
      delta(i8,3)=eyz*ny
      end if
      if(i8.eq.7.or.i8.eq.8) then
      delta(i8,1)=exy*ny+exz*nz
      delta(i8,2)=eyy*ny+eyz*nz
      delta(i8,3)=eyz*ny+ezz*nz
      end if
2051  continue
      do 2050 i=1,nx-1
      m=nxy*(nz-1)+nx*(ny-1)+i
      do 1905 nn=1,3
      do 1905 mm=1,8
      sum=0.0
      do 2049 m3=1,3
      do 2049 m8=1,8
      sum=sum+delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)
      C=C+0.5*delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)*delta(mm,nn)
2049  continue
      b(ib(m,is(mm)),nn)=b(ib(m,is(mm)),nn)+sum
1905  continue
2050  continue
c  x=nx y=ny z=nz corner
      do 2061 i3=1,3
      do 2061 i8=1,8
      delta(i8,i3)=0.0
      if(i8.eq.2) then
      delta(i8,1)=exx*nx
      delta(i8,2)=exy*nx
      delta(i8,3)=exz*nx
      end if
      if(i8.eq.4) then
      delta(i8,1)=exy*ny
      delta(i8,2)=eyy*ny
      delta(i8,3)=eyz*ny
      end if
      if(i8.eq.5) then
      delta(i8,1)=exz*nz
      delta(i8,2)=eyz*nz
      delta(i8,3)=ezz*nz
      end if
      if(i8.eq.8) then
      delta(i8,1)=exy*ny+exz*nz
      delta(i8,2)=eyy*ny+eyz*nz
      delta(i8,3)=eyz*ny+ezz*nz
      end if
      if(i8.eq.6) then
      delta(i8,1)=exx*nx+exz*nz
      delta(i8,2)=exy*nx+eyz*nz
      delta(i8,3)=exz*nx+ezz*nz
      end if
      if(i8.eq.3) then
      delta(i8,1)=exx*nx+exy*ny
      delta(i8,2)=exy*nx+eyy*ny
      delta(i8,3)=exz*nx+eyz*ny
      end if
      if(i8.eq.7) then
      delta(i8,1)=exx*nx+exy*ny+exz*nz
      delta(i8,2)=exy*nx+eyy*ny+eyz*nz
      delta(i8,3)=exz*nx+eyz*ny+ezz*nz
      end if
2061  continue
      m=nx*ny*nz
      do 1906 nn=1,3
      do 1906 mm=1,8
      sum=0.0
      do 2059 m3=1,3
      do 2059 m8=1,8
      sum=sum+delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)
      C=C+0.5*delta(m8,m3)*dk(pix(m),m8,m3,mm,nn)*delta(mm,nn)
2059  continue
      b(ib(m,is(mm)),nn)=b(ib(m,is(mm)),nn)+sum
1906  continue

      return
      end

c  Subroutine computes the total energy, utot, and the gradient, gb

      subroutine energy(nx,ny,nz,ns,utot)

	real u(64000000,3),gb(64000000,3)
	real b(64000000,3),C,utot
	real dk(100,8,3,8,3)
 	integer*4 ib(64000000,27)
 	integer*2 pix(64000000)

	common/list2/exx,eyy,ezz,exz,eyz,exy
	common/list3/ib
        common/list4/pix
	common/list5/dk,b,C
	common/list6/u
	common/list7/gb

        do 2090 m3=1,3
	do 2090 m=1,ns
	gb(m,m3)=0.0
2090	continue

c  Do global matrix multiply via small stiffness matrices, gb = A * u
c  The long statement below correctly brings in all the terms from
c  the global matrix A using only the small stiffness matrices.

        do 3000 j=1,3
        do 3000 n=1,3
	do 3000 m=1,ns
      gb(m,j)=gb(m,j)+u(ib(m,1),n)*( dk(pix(ib(m,27)),1,j,4,n)
     &+dk(pix(ib(m,7)),2,j,3,n)
     &+dk(pix(ib(m,25)),5,j,8,n)+dk(pix(ib(m,15)),6,j,7,n) )+
     &u(ib(m,2),n)*( dk(pix(ib(m,27)),1,j,3,n)
     &+dk(pix(ib(m,25)),5,j,7,n) )+
     &u(ib(m,3),n)*( dk(pix(ib(m,27)),1,j,2,n)+dk(pix(ib(m,5)),4,j,3,n)+
     &dk(pix(ib(m,13)),8,j,7,n)+dk(pix(ib(m,25)),5,j,6,n) )+
     &u(ib(m,4),n)*( dk(pix(ib(m,5)),4,j,2,n)
     &+dk(pix(ib(m,13)),8,j,6,n) )+
     &u(ib(m,5),n)*( dk(pix(ib(m,6)),3,j,2,n)+dk(pix(ib(m,5)),4,j,1,n)+
     &dk(pix(ib(m,14)),7,j,6,n)+dk(pix(ib(m,13)),8,j,5,n) )+
     &u(ib(m,6),n)*( dk(pix(ib(m,6)),3,j,1,n)
     &+dk(pix(ib(m,14)),7,j,5,n) )+
     &u(ib(m,7),n)*( dk(pix(ib(m,6)),3,j,4,n)+dk(pix(ib(m,7)),2,j,1,n)+
     &dk(pix(ib(m,14)),7,j,8,n)+dk(pix(ib(m,15)),6,j,5,n) )+
     &u(ib(m,8),n)*( dk(pix(ib(m,7)),2,j,4,n)
     &+dk(pix(ib(m,15)),6,j,8,n) )+
     &u(ib(m,9),n)*( dk(pix(ib(m,25)),5,j,4,n)
     &+dk(pix(ib(m,15)),6,j,3,n) )+
     &u(ib(m,10),n)*( dk(pix(ib(m,25)),5,j,3,n) )+
     &u(ib(m,11),n)*( dk(pix(ib(m,13)),8,j,3,n)
     &+dk(pix(ib(m,25)),5,j,2,n) )+
     &u(ib(m,12),n)*( dk(pix(ib(m,13)),8,j,2,n) )+
     &u(ib(m,13),n)*( dk(pix(ib(m,13)),8,j,1,n)
     &+dk(pix(ib(m,14)),7,j,2,n) )+
     &u(ib(m,14),n)*( dk(pix(ib(m,14)),7,j,1,n) )+
     &u(ib(m,15),n)*( dk(pix(ib(m,14)),7,j,4,n)
     &+dk(pix(ib(m,15)),6,j,1,n) )+
     &u(ib(m,16),n)*( dk(pix(ib(m,15)),6,j,4,n) )+
     &u(ib(m,17),n)*( dk(pix(ib(m,27)),1,j,8,n)
     &+dk(pix(ib(m,7)),2,j,7,n) )+
     &u(ib(m,18),n)*( dk(pix(ib(m,27)),1,j,7,n) )+
     &u(ib(m,19),n)*( dk(pix(ib(m,27)),1,j,6,n)
     &+dk(pix(ib(m,5)),4,j,7,n) )+
     &u(ib(m,20),n)*( dk(pix(ib(m,5)),4,j,6,n) )+
     &u(ib(m,21),n)*( dk(pix(ib(m,5)),4,j,5,n)
     &+dk(pix(ib(m,6)),3,j,6,n) )+
     &u(ib(m,22),n)*( dk(pix(ib(m,6)),3,j,5,n) )+
     &u(ib(m,23),n)*( dk(pix(ib(m,6)),3,j,8,n)
     &+dk(pix(ib(m,7)),2,j,5,n) )+
     &u(ib(m,24),n)*( dk(pix(ib(m,7)),2,j,8,n) )+
     &u(ib(m,25),n)*( dk(pix(ib(m,14)),7,j,3,n)
     &+dk(pix(ib(m,13)),8,j,4,n)+
     &dk(pix(ib(m,15)),6,j,2,n)+dk(pix(ib(m,25)),5,j,1,n) )+
     &u(ib(m,26),n)*( dk(pix(ib(m,6)),3,j,7,n)
     &+dk(pix(ib(m,5)),4,j,8,n)+
     &dk(pix(ib(m,27)),1,j,5,n)+dk(pix(ib(m,7)),2,j,6,n) )+
     &u(ib(m,27),n)*( dk(pix(ib(m,27)),1,j,1,n)
     &+dk(pix(ib(m,7)),2,j,2,n)+
     &dk(pix(ib(m,6)),3,j,3,n)+dk(pix(ib(m,5)),4,j,4,n)
     &+dk(pix(ib(m,25)),5,j,5,n)+
     &dk(pix(ib(m,15)),6,j,6,n)+dk(pix(ib(m,14)),7,j,7,n)+
     &dk(pix(ib(m,13)),8,j,8,n) )
3000	continue

	utot=C
        do 3100 m3=1,3
	do 3100 m=1,ns
	utot=utot+0.5*u(m,m3)*gb(m,m3)+b(m,m3)*u(m,m3)
	gb(m,m3)=gb(m,m3)+b(m,m3)
3100	continue

        return
        end

c  Subroutine that carries out the conjugate gradient relaxation process

      subroutine dembx(ns,Lstep,gg,dk,gtest,ldemb,kkk)
      real gb(64000000,3),u(64000000,3),dk(100,8,3,8,3)
      real h(64000000,3),Ah(64000000,3)
      real lambda,gamma
      integer*4 ib(64000000,27)
      integer*2 pix(64000000)

      common/list3/ib
      common/list4/pix
      common/list6/u
      common/list7/gb
      common/list9/h,Ah

c  Initialize the conjugate direction vector on first call to dembx only
c  For calls to dembx after the first, we want to continue using the
c  value of h determined in the previous call. Of course, if npoints is
c  greater than 1, this initialization step will be run for every new
c  microstructure used, as kkk is reset to 1 every time the counter micro
c  is increased.
      if(kkk.eq.1) then
      do 500 m3=1,3
      do 500 m=1,ns
      h(m,m3)=gb(m,m3)
500   continue
      end if
c  Lstep counts the number of conjugate gradient steps taken in
c  each call to dembx
      Lstep=0

      do 800 ijk=1,ldemb
      Lstep=Lstep+1

      do 290 m3=1,3
      do 290 m=1,ns
      Ah(m,m3)=0.0
290   continue
c  Do global matrix multiply via small stiffness matrices, Ah = A * h
c  The long statement below correctly brings in all the terms from
c  the global matrix A using only the small stiffness matrices dk.
        do 400 j=1,3
        do 400 n=1,3
	do 400 m=1,ns
      Ah(m,j)=Ah(m,j)+h(ib(m,1),n)*( dk(pix(ib(m,27)),1,j,4,n)
     &+dk(pix(ib(m,7)),2,j,3,n)
     &+dk(pix(ib(m,25)),5,j,8,n)+dk(pix(ib(m,15)),6,j,7,n) )+
     &h(ib(m,2),n)*( dk(pix(ib(m,27)),1,j,3,n)
     &+dk(pix(ib(m,25)),5,j,7,n) )+
     &h(ib(m,3),n)*( dk(pix(ib(m,27)),1,j,2,n)+dk(pix(ib(m,5)),4,j,3,n)+
     &dk(pix(ib(m,13)),8,j,7,n)+dk(pix(ib(m,25)),5,j,6,n) )+
     &h(ib(m,4),n)*( dk(pix(ib(m,5)),4,j,2,n)
     &+dk(pix(ib(m,13)),8,j,6,n) )+
     &h(ib(m,5),n)*( dk(pix(ib(m,6)),3,j,2,n)+dk(pix(ib(m,5)),4,j,1,n)+
     &dk(pix(ib(m,14)),7,j,6,n)+dk(pix(ib(m,13)),8,j,5,n) )+
     &h(ib(m,6),n)*( dk(pix(ib(m,6)),3,j,1,n)
     &+dk(pix(ib(m,14)),7,j,5,n) )+
     &h(ib(m,7),n)*( dk(pix(ib(m,6)),3,j,4,n)+dk(pix(ib(m,7)),2,j,1,n)+
     &dk(pix(ib(m,14)),7,j,8,n)+dk(pix(ib(m,15)),6,j,5,n) )+
     &h(ib(m,8),n)*( dk(pix(ib(m,7)),2,j,4,n)
     &+dk(pix(ib(m,15)),6,j,8,n) )+
     &h(ib(m,9),n)*( dk(pix(ib(m,25)),5,j,4,n)
     &+dk(pix(ib(m,15)),6,j,3,n) )+
     &h(ib(m,10),n)*( dk(pix(ib(m,25)),5,j,3,n) )+
     &h(ib(m,11),n)*( dk(pix(ib(m,13)),8,j,3,n)
     &+dk(pix(ib(m,25)),5,j,2,n) )+
     &h(ib(m,12),n)*( dk(pix(ib(m,13)),8,j,2,n) )+
     &h(ib(m,13),n)*( dk(pix(ib(m,13)),8,j,1,n)
     &+dk(pix(ib(m,14)),7,j,2,n) )+
     &h(ib(m,14),n)*( dk(pix(ib(m,14)),7,j,1,n) )+
     &h(ib(m,15),n)*( dk(pix(ib(m,14)),7,j,4,n)
     &+dk(pix(ib(m,15)),6,j,1,n) )+
     &h(ib(m,16),n)*( dk(pix(ib(m,15)),6,j,4,n) )+
     &h(ib(m,17),n)*( dk(pix(ib(m,27)),1,j,8,n)
     &+dk(pix(ib(m,7)),2,j,7,n) )+
     &h(ib(m,18),n)*( dk(pix(ib(m,27)),1,j,7,n) )+
     &h(ib(m,19),n)*( dk(pix(ib(m,27)),1,j,6,n)
     &+dk(pix(ib(m,5)),4,j,7,n) )+
     &h(ib(m,20),n)*( dk(pix(ib(m,5)),4,j,6,n) )+
     &h(ib(m,21),n)*( dk(pix(ib(m,5)),4,j,5,n)
     &+dk(pix(ib(m,6)),3,j,6,n) )+
     &h(ib(m,22),n)*( dk(pix(ib(m,6)),3,j,5,n) )+
     &h(ib(m,23),n)*( dk(pix(ib(m,6)),3,j,8,n)
     &+dk(pix(ib(m,7)),2,j,5,n) )+
     &h(ib(m,24),n)*( dk(pix(ib(m,7)),2,j,8,n) )+
     &h(ib(m,25),n)*( dk(pix(ib(m,14)),7,j,3,n)
     &+dk(pix(ib(m,13)),8,j,4,n)+
     &dk(pix(ib(m,15)),6,j,2,n)+dk(pix(ib(m,25)),5,j,1,n) )+
     &h(ib(m,26),n)*( dk(pix(ib(m,6)),3,j,7,n)
     &+dk(pix(ib(m,5)),4,j,8,n)+
     &dk(pix(ib(m,27)),1,j,5,n)+dk(pix(ib(m,7)),2,j,6,n) )+
     &h(ib(m,27),n)*( dk(pix(ib(m,27)),1,j,1,n)
     &+dk(pix(ib(m,7)),2,j,2,n)+
     &dk(pix(ib(m,6)),3,j,3,n)+dk(pix(ib(m,5)),4,j,4,n)
     &+dk(pix(ib(m,25)),5,j,5,n)+
     &dk(pix(ib(m,15)),6,j,6,n)+dk(pix(ib(m,14)),7,j,7,n)+
     &dk(pix(ib(m,13)),8,j,8,n) )
400    continue

      hAh=0.0
      do 530 m3=1,3
      do 530 m=1,ns
      hAh=hAh+h(m,m3)*Ah(m,m3)
530   continue

      lambda=gg/hAh
      do 540 m3=1,3
      do 540 m=1,ns
      u(m,m3)=u(m,m3)-lambda*h(m,m3)
      gb(m,m3)=gb(m,m3)-lambda*Ah(m,m3)
540   continue

      gglast=gg
      gg=0.0
      do 550 m3=1,3
      do 550 m=1,ns
      gg=gg+gb(m,m3)*gb(m,m3)
550   continue
      if(gg.lt.gtest) goto 1000

      gamma=gg/gglast
      do 570 m3=1,3
      do 570 m=1,ns
      h(m,m3)=gb(m,m3)+gamma*h(m,m3)
570   continue

800   continue

1000  continue

      return
      end

c  Subroutine that computes the six average stresses and six
c  average strains.

      subroutine stress(nx,ny,nz,ns,iwriteFlag)
      real u(64000000,3),gb(64000000,3),uu(8,3)
      real dndx(8),dndy(8),dndz(8),es(6,8,3),cmod(100,6,6)
      integer*4 ib(64000000,27)
      integer*2 pix(64000000)

	common/list1/strxx,stryy,strzz,strxz,stryz,strxy
	common/list2/exx,eyy,ezz,exz,eyz,exy
	common/list3/ib
        common/list4/pix
	common/list6/u
	common/list7/gb
        common/list8/cmod
	common/list10/sxx,syy,szz,sxz,syz,sxy

      nxy=nx*ny

c  set up single element strain matrix
c  dndx, dndy, and dndz are the components of the average strain
c  matrix in a pixel

      dndx(1)=-0.25
      dndx(2)=0.25
      dndx(3)=0.25
      dndx(4)=-0.25
      dndx(5)=-0.25
      dndx(6)=0.25
      dndx(7)=0.25
      dndx(8)=-0.25
      dndy(1)=-0.25
      dndy(2)=-0.25
      dndy(3)=0.25
      dndy(4)=0.25
      dndy(5)=-0.25
      dndy(6)=-0.25
      dndy(7)=0.25
      dndy(8)=0.25
      dndz(1)=-0.25
      dndz(2)=-0.25
      dndz(3)=-0.25
      dndz(4)=-0.25
      dndz(5)=0.25
      dndz(6)=0.25
      dndz(7)=0.25
      dndz(8)=0.25
c  Build averaged strain matrix, follows code in femat, but for average
c  strain over the pixel, not the strain at a point.
      do 2799 n1=1,6
      do 2799 n2=1,8
      do 2799 n3=1,3
      es(n1,n2,n3)=0.0
2799  continue
      do 2797 n=1,8
      es(1,n,1)=dndx(n)
      es(2,n,2)=dndy(n)
      es(3,n,3)=dndz(n)
      es(4,n,1)=dndz(n)
      es(4,n,3)=dndx(n)
      es(5,n,2)=dndz(n)
      es(5,n,3)=dndy(n)
      es(6,n,1)=dndy(n)
      es(6,n,2)=dndx(n)
2797  continue

c  Compute components of the average stress and strain tensors in each pixel
      strxx=0.0
      stryy=0.0
      strzz=0.0
      strxz=0.0
      stryz=0.0
      strxy=0.0
      sxx=0.0
      syy=0.0
      szz=0.0
      sxz=0.0
      syz=0.0
      sxy=0.0

	if(iwriteFlag .eq. 1) open(unit=3, file='stressField.dat')


      do 470 k=1,nz
      do 470 j=1,ny
      do 470 i=1,nx
      m=(k-1)*nxy+(j-1)*nx+i
c  load in elements of 8-vector using pd. bd. conds.
      do 9898 mm=1,3
      uu(1,mm)=u(m,mm)
      uu(2,mm)=u(ib(m,3),mm)
      uu(3,mm)=u(ib(m,2),mm)
      uu(4,mm)=u(ib(m,1),mm)
      uu(5,mm)=u(ib(m,26),mm)
      uu(6,mm)=u(ib(m,19),mm)
      uu(7,mm)=u(ib(m,18),mm)
      uu(8,mm)=u(ib(m,17),mm)
9898  continue
c  Correct for periodic boundary conditions, some displacements are wrong
c  for a pixel on a periodic boundary.  Since they come from an opposite
c  face, need to put in applied strain to correct them.
      if(i.eq.nx) then
      uu(2,1)=uu(2,1)+exx*nx
      uu(2,2)=uu(2,2)+exy*nx
      uu(2,3)=uu(2,3)+exz*nx
      uu(3,1)=uu(3,1)+exx*nx
      uu(3,2)=uu(3,2)+exy*nx
      uu(3,3)=uu(3,3)+exz*nx
      uu(6,1)=uu(6,1)+exx*nx
      uu(6,2)=uu(6,2)+exy*nx
      uu(6,3)=uu(6,3)+exz*nx
      uu(7,1)=uu(7,1)+exx*nx
      uu(7,2)=uu(7,2)+exy*nx
      uu(7,3)=uu(7,3)+exz*nx
      end if
      if(j.eq.ny) then
      uu(3,1)=uu(3,1)+exy*ny
      uu(3,2)=uu(3,2)+eyy*ny
      uu(3,3)=uu(3,3)+eyz*ny
      uu(4,1)=uu(4,1)+exy*ny
      uu(4,2)=uu(4,2)+eyy*ny
      uu(4,3)=uu(4,3)+eyz*ny
      uu(7,1)=uu(7,1)+exy*ny
      uu(7,2)=uu(7,2)+eyy*ny
      uu(7,3)=uu(7,3)+eyz*ny
      uu(8,1)=uu(8,1)+exy*ny
      uu(8,2)=uu(8,2)+eyy*ny
      uu(8,3)=uu(8,3)+eyz*ny
      end if
      if(k.eq.nz) then
      uu(5,1)=uu(5,1)+exz*nz
      uu(5,2)=uu(5,2)+eyz*nz
      uu(5,3)=uu(5,3)+ezz*nz
      uu(6,1)=uu(6,1)+exz*nz
      uu(6,2)=uu(6,2)+eyz*nz
      uu(6,3)=uu(6,3)+ezz*nz
      uu(7,1)=uu(7,1)+exz*nz
      uu(7,2)=uu(7,2)+eyz*nz
      uu(7,3)=uu(7,3)+ezz*nz
      uu(8,1)=uu(8,1)+exz*nz
      uu(8,2)=uu(8,2)+eyz*nz
      uu(8,3)=uu(8,3)+ezz*nz
      end if

c  local stresses and strains in a pixel
      str11=0.0
      str22=0.0
      str33=0.0
      str13=0.0
      str23=0.0
      str12=0.0
      s11=0.0
      s22=0.0
      s33=0.0
      s13=0.0
      s23=0.0
      s12=0.0
      do 465 n3=1,3
      do 465 n8=1,8
      s11=s11+es(1,n8,n3)*uu(n8,n3)
      s22=s22+es(2,n8,n3)*uu(n8,n3)
      s33=s33+es(3,n8,n3)*uu(n8,n3)
      s13=s13+es(4,n8,n3)*uu(n8,n3)
      s23=s23+es(5,n8,n3)*uu(n8,n3)
      s12=s12+es(6,n8,n3)*uu(n8,n3)
      do 465 n=1,6
      str11=str11+cmod(pix(m),1,n)*es(n,n8,n3)*uu(n8,n3)
      str22=str22+cmod(pix(m),2,n)*es(n,n8,n3)*uu(n8,n3)
      str33=str33+cmod(pix(m),3,n)*es(n,n8,n3)*uu(n8,n3)
      str13=str13+cmod(pix(m),4,n)*es(n,n8,n3)*uu(n8,n3)
      str23=str23+cmod(pix(m),5,n)*es(n,n8,n3)*uu(n8,n3)
      str12=str12+cmod(pix(m),6,n)*es(n,n8,n3)*uu(n8,n3)
465   continue
c  sum local strains and stresses into global values
      strxx=strxx+str11
      stryy=stryy+str22
      strzz=strzz+str33
      strxz=strxz+str13
      stryz=stryz+str23
      strxy=strxy+str12
      sxx=sxx+s11
      syy=syy+s22
      szz=szz+s33
      sxz=sxz+s13
      syz=syz+s23
      sxy=sxy+s12

	if(iwriteFlag .eq. 1) write(3,6070) str11, str22, str33, str13,
     &	str23, str12

6070	format(6(e12.5, 2x))



470   continue

c  Volume average of global stresses and strains
      strxx=strxx/float(ns)
      stryy=stryy/float(ns)
      strzz=strzz/float(ns)
      strxz=strxz/float(ns)
      stryz=stryz/float(ns)
      strxy=strxy/float(ns)
      sxx=sxx/float(ns)
      syy=syy/float(ns)
      szz=szz/float(ns)
      sxz=sxz/float(ns)
      syz=syz/float(ns)
      sxy=sxy/float(ns)

      return
      end

c  Subroutine that counts volume fractions

      subroutine assig(ns,nphase,prob)
      integer*2 pix(64000000)
      real prob(100)

      common/list4/pix

	do 90 i=1,nphase
	prob(i)=0.0
90      continue

        do 100 m=1,ns
	do 100 i=1,nphase
        if(pix(m).eq.i) then
	prob(i)=prob(i)+1
	end if
100     continue

	do 110 i=1,nphase
	prob(i)=prob(i)/float(ns)
110     continue

        return
        end

c  Subroutine that sets up microstructural image

      subroutine ppixel(nx,ny,nz,ns,nphase,nblank)
      integer*2 pix(64000000)
      common/list4/pix

c  (USER)  If you want to set up a test image inside the program, instead of
c  reading it in from a file, this should be done inside this subroutine.

      nxy=nx*ny

	do 291 k=1,nblank
	read(9,*)
291   continue

      do 200 k=1,nz
      do 200 j=1,ny
      do 200 i=1,nx
      m=nxy*(k-1)+nx*(j-1)+i
      read(9,*) pix(m)
	pix(m) = pix(m) + 1
200   continue

c  Check for wrong phase labels--less than 1 or greater than nphase
       do 500 m=1,ns
       if(pix(m).lt.1) then
        write(7,*) 'Phase label in pix < 1--error at ',m
       end if
       if(pix(m).gt.nphase) then
        write(7,*) 'Phase label in pix > nphase--error at ',m
       end if
500    continue

      return
      end
