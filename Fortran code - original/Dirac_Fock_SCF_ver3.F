CCCCCCCCCCCCCCCCCCC VERSION LOG CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								    C
CC    New  Features in version 3:  				    C 
C      								    C
C     MP2 Correlation energy is added.				    C
C								    C
CC    New Features in version 2:				    C
C								    C
C     Input structure is made more user friendly.		    C
C								    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CInput data format:
C----------------------------------------------------------
C
C Step 1: Open the input file in the name of 'DFINP.DAT'
C
C Step 2: Read the following items sequentially:
C
C   (a) title : character describing your intention
C
C   (b) atom : Symbol of atom (in two capital letters; for single letter
C              sumbol it should follow with a blank space
C
C   (c) nsym:  No. of symmetry to be used
C
C   (d) nbas(i) :  No. of basis for each symmetry
C
C   (e) nocorb(i) : No. of occupied orbitals for each symmetry
C
C   (f) kap(i)  :   Kappa values for each symmetry
C
C   (g) orbj(i) :  J- values for each symmetry
C
C   (h) method : Type of methods (UB, ET, WT and EE)
C
C   (i) alpha and beta : values for the above methods
C
C   (j) maxit : Max. no. of iterations it should check for convergence
C
C   (k) npower : Tolerence parameter for 10^(-npower)
C
C   (l) amass : Mass number
C
C   (m) z : atomic number
C
C   (n) Optional flags:
C
C       isgo   :  1  STOs -->  r^nk * exp[-alpha*r] 
C              :  2  GTOs -->  r^nk * exp[-alpha*r^2] 
C              :  3  QUAD -->  r^nk * exp[-alpha*r^4]
C              :  4  dQUAD ->  r^(2*nk) * exp[-alpha*r^4]
C              :  5  nMULT ->  r^nk * exp[-alpha*r^nk]
C    
C       jpnt   :  0  default, unless print out something
C    
C       ibrt   :  0  Coulomb
C              :  1  Magnetic part of Breit
C              :  2  Full Breit
C
C       ispd   :  0  Usual radial distribution
C              :  1  Larger radial distribution
C
C       idflt  :  0  Consider default values of step size, r_0 
C                      and no. of grids 
C              :  1  Read them through input file
C
C       ieng   :  0  Write all the orbitals in the output
C              :  1  Write only selective orbitals for active space
C                    and write the frozen orbital energies in fort.15
C
C       inuc   :  0  for point nuclei       
C              :  1  for Fermi charge distribution
C              :  2  for uniform charge distribution
C
C       ichrg  :  0  for neutral atom
C              :  1  for singly ion
C              :  2  for doubly ion and so on!
C
C       istp   :  0  no isotope shift calculation
C              :  1  specific isotope shift self-consistently METHOD 1
C              :  2  specific isotope shift self-consistently METHOD 2
C              :  3  first order correction due to SMS METHOD 1
C              :  4  first order correction due to SMS METHOD 2
C
C       ialpvr :  0  no variation in alpha (also, velocity of light c)
C              :  1  variation in alpha at second order
C              :  2  variation in alpha at fourth order
C
C       ihbd   :  0  Pure analytical
C              :  1  Hybridise the orbitals with GRASP orbitals
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  These options are not used and commented in subroutine "output"    CC
CC      intg   :  0  If do not want to write integrals separately      CC
CC             :  1  to write only one-body integrals in fort.10       CC
CC             :  2  both one-body and two-body in fort.10 and fort.11 CC
CC             :  3  only two-body in fort.11                          CC
CC  **NOTE: Follow the instruction below                               CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       intg   :  0  Not to append nuclear potential in wfn.dat
C              :  1  To append nuclear potential in wfn.dat
C        
C       ivn    :  0  V^N-1 potential for all the N-1 orbitals 
C              :  1  V^N potential for core and V^N-1 potential for virtuals
C              :  2  V^N potential for N orbitals
C
C       ival   :  0  V^N-1 potential
C              :  i  orbital index for V^N potential
c       
C       if ispd=1, read cla   (lies between 0.8-0.9)
C
C       if idflt=1, read h, r_0 and n
C       
C       if istp=1 or 2, read cstp
C       
C       if ialpvr=1  read x variation parameter
C
C       if ihbd=1, read no. of GRASP orbitals for each symm
C
C       if ieng=1, read minimum and maximum orbital indices for each 
C                  symmetry in order to freeze the inactive orbitals
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C---------------------------------------------------------------------C
C The output files are:                                               C
C                                                                     C
C     SCFOUT.OUT:  contains all information of orbital energies       C
C     wfn.dat  :   binary file containing DF wave functions           C
C     fort.15  :   frozen energies                                    C
C                                                                     C
C---------------------------------------------------------------------C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Note : FOR THE SAMPLE INPUT FILE SEE THE END OF THE FILE

C These are the input and output file numbers
#define STD12INT 2
#define STDENG   3
#define STDOUT   4
#define  STDIN   5
#define   IWFN   7
#define  STDXP   8
#define IGRASP   9

C These are the global dimensions
#define   MNS   15
#define MAXVAL  10
#define MNOCC   50
#define MNBAS  500
#define   MNB  900
#define    MN  1550

C This is the tiny number to check density consistency in each iteration
#define  TINY  1.0D-12
#define  ACCY  1.0d-18


      Program DFATOM

      Implicit Real*8(a-h,o-z)

      Character*2 Method,atom 
      character*80 title
      LOGICAL MRCC, IBASOPT

      Dimension Rho(MN)

      COMMON/splflags/MRCC,IBASOPT

C     SETTING SOME SPECIAL FLAGS; USE THEM ONLY IF REQUIRED
      MRCC    = .TRUE.    
      IBASOPT = .TRUE.  

      Call Readinp(title,atom,Method,cla,cstp,ival,dltx)
      Call DEFULTVAL(title,atom,METHOD,iele)
      Call CHECKNUC(atom,IELE,RRMS,NPROT)
      Call Setgrd
      Call NUCPARM(RRMS,dltx)
      Call NUCPOT(RHO,NPROT)
      Call Setqic
      Call SETBASIS(Method,cla)
      Call SETMATRIX(RHO)
      Call Scfiter(cstp,ival)

      Call Output(RHO)

      Stop
      End

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE READINP(title,atom,METHOD,cla,cstp,ival,dltx)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit real*8 (a-h,o-z)
      character*80 title
      character*2 atom,method

* New common blocks
      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     :      /UB/alpha1,beta1
     :      /WT/alpha2,beta2,gamma1,delta
     :      /ET/alpha3(MNS),beta3(MNS)
     :      /EE/expnts(MNB),trnopt(MNB)
     :      /conv/crit,maxit,npower
     :      /jvalue/orbj(MNS),kap(MNS)
     :     /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /propopt/istp,ialpvr
     :      /nucinf/amass,z
     :      /freez/minorb(MNS),maxorb(MNS)
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /basc/nbasis,nocc,mbasis
     :      /skipg/nbasg(MNS),nbasisg,mbasisg,nskipg(MNS)
     :      /splflags/MRCC,IBASOPT

      DATA ORBJ/0.5D0,0.5D0,1.5D0,1.5D0,2.5D0,2.5D0,3.5D0,3.5D0,4.5D0,
     :                            4.5D0,5.5D0,5.5D0,6.5D0,6.5D0,7.5D0/
      DATA KAP/-1,1,-2,2,-3,3,-4,4,-5,5,-6,6,-7,7,-8/

C  Inputfile name
      OPEN(STDIN,file='SCFINP.DAT',form='formatted',status='old')

C Outputfile name
      OPEN(STDOUT,file='SCFOUT.OUT',form='formatted',status='unknown')
      OPEN(IWFN,file='wfn.dat',form='unformatted',status='unknown')

      read(STDIN,'(a)')title
      read(STDIN,'(a)')atom
      read(STDIN,*)nsym,(nbas(i),i=1,nsym)
      read(STDIN,*)(nocorb(i),i=1,nsym)

        nbasis=0
        do isym=1,nsym
        nbasis=nbasis+nbas(isym)  
        enddo
        mbasis=2*nbasis

      read(STDIN,'(a)')method

      if(method.eq.'UB') read(STDIN,*)alpha1,beta1              ! UNIVERSAL BASIS
      if(method.eq.'WT') read(STDIN,*)alpha2,beta2,gamma1,delta ! WELL-TEMPERED BASIS
      if(method.eq.'ET') read(STDIN,*)(alpha3(i),i=1,nsym),     ! EVEN-TEMPERED BASIS
     :                                 (beta3(i),i=1,nsym)      

C     READING "EXPLICIT EXPONENTS" FROM A FILE EXPONENTS.DAT

      IF(IBASOPT.or.METHOD.eq.'EE')THEN
      OPEN(STDXP,file='EXPONENTS.DAT',form='formatted',status='unknown')
      ENDIF

      if(method.eq.'EE') then                                  

        read(STDXP,*)nsymee,(nbas(i),i=1,nsym)
        if (nsymee.ne.nsym)then
        print*,'CHECK THE NUMBER OF SYMMETRIES; TERMINATING NOW'
        endif

        nbasis=0
        do isym=1,nsym
        nbasis=nbasis+nbas(isym)  ! deliberatey overwriting nbas(isym),nbasis
        enddo
        mbasis=2*nbasis           ! deliberatey overwriting mbasis

        read(STDXP,*)(expnts(i),trnopt(i),i=1,nbasis)
C        read(STDXP,*)(expnts(i),i=1,nbasis)
      endif

      read(STDIN,*)maxit,npower,amass,z
      read(STDIN,*)

      read(STDIN,*)isgo,jpnt,ibrt,ispd,idflt,Ieng,inuc,ichrg
      read(STDIN,*)

      read(STDIN,*)istp,ialpvr,ihbd,intg,ivn,ival
      read(STDIN,*)

C     INITIALIZING SOME ARRAYS
      do isym=1,nsym
      nbasg(isym) =0
      minorb(isym)=0
      maxorb(isym)=0
      nvorb(isym) =0
      enddo     

      if(ispd.eq.1)read(STDIN,*)cla
      If(idflt.eq.1)read(STDIN,*)h,rnt,n

      if(istp.eq.1.or.istp.eq.2)then
      read(STDIN,*)cstp
      endif

      if(ialpvr.eq.1)then
      read(STDIN,*)dltx
      endif

      if(ihbd.eq.1)then
      read(STDIN,*)(nbasg(i),i=1,nsym)
      endif

      if(ivn.eq.1)then
      do i=1,nsym
      nbasg(i)=nocorb(i)
      enddo
      endif

      If(Ieng.eq.1)read(STDIN,*)(minorb(i),i=1,nsym),
     :(maxorb(i),i=1,nsym)

      if(ivn.eq.2)read(STDIN,*)(nvorb(i),i=1,nsym)

      nocc=0
      ichk=0
      do isym=1,nsym
      nocc=nocc+nocorb(isym)
      ichk=ichk+nocorb(isym)*(2*orbj(isym)+1)+nvorb(isym)
      enddo

C     CHECKING FOR THE CORRECTNESS OF THE INPUT DATA
      if(int(z).ne.(ichk+ichrg)) then
      print*,'Check for the correctness of input data; Terminating now'
      STOP
      endif

      print*,'                     *********                         '
      print*,'<<<<<<<<<<<<<<<<<<<<<<<<<O>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      print*,'                     *********                         '
      print*,'SCF RUN FOR DF WAVE FUNCTION'
      print*,'Atom= ',atom,' Method=',method,' Charge=',ichrg


      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DEFULTVAL(title,atom,METHOD,iele)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit real*8 (a-h,o-z)
      character*80 title
      character*2 atom,method

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     :      /UB/alpha1,beta1
     :      /WT/alpha2,beta2,gamma1,delta
     :      /ET/alpha3(MNS),beta3(MNS)
     :      /EE/expnts(MNB),trnopt(MNB)
     :      /conv/crit,maxit,npower
     :      /jvalue/orbj(MNS),kap(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /nucinf/amass,z
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /basc/nbasis,nocc,mbasis
     :      /freez/minorb(MNS),maxorb(MNS)

      If(idflt.eq.0)then
      n = 740
      h = 3.0d-02
      rnt=2.0d-06
      EndIf

      crit=(0.1d0)**npower

      iele=0
      ii=0

      do isym=1,nsym
      nn=nocorb(isym)

      if(nn.ne.0)then

      do ibas=1,nocorb(isym)

      iele=iele+2*orbj(isym)+1
      ii=ii+1

      enddo

      endif

      enddo

      If(ii.ne.nocc)then
      Print*,'Wrong in counting total occupied orbitals'
      Stop
      Endif

      write(STDOUT,9990)
      write(STDOUT,9999)
      write(STDOUT,9998)
      write(STDOUT,9998)
      write(STDOUT,9998)
      write(STDOUT,9997)
      write(STDOUT,9998)
      write(STDOUT,9996)
      write(STDOUT,9998)
      write(STDOUT,9995)
      write(STDOUT,9998)
      write(STDOUT,9998)
      write(STDOUT,9998)
      write(STDOUT,9999)

 9999 format(20x,'**************************************')
 9998 format(20x,'*                                    *')
 9997 format(20x,'*        Relativistic SCF Run        *')
 9996 format(20x,'*            written by              *')
 9995 format(20x,'* B. K. Sahoo. Start date 12-05-2008 *')
 9990 format(//)

      if(isgo.eq.1)write(STDOUT,1001)
      if(isgo.eq.2)write(STDOUT,1002)
      if(ispd.eq.1)write(STDOUT,1003)
      if(ibrt.eq.1)write(STDOUT,1004)
      if(ibrt.eq.2)write(STDOUT,1005)
      if(idflt.eq.1)write(STDOUT,1006)
      if(ieng.eq.1)write(STDOUT,1007)
      if(inuc.eq.0)write(STDOUT,1008)
      if(inuc.eq.1)write(STDOUT,1009)
      if(inuc.eq.2)write(STDOUT,1015)
      if(ichrg.eq.0)write(STDOUT,1010)
      if(ichrg.eq.1)write(STDOUT,1011)
      if(ichrg.eq.2)write(STDOUT,1012)
      if(ichrg.eq.3)write(STDOUT,1013)
      if(ichrg.gt.3)write(STDOUT,1014)
      write(STDOUT,9999)

  301 format('>>>>>>>>>>>>>>> INPUT DATA >>>>>>>>>>>>>>>>>>>>>>>>')
      write(STDOUT,9994)title
      write(STDOUT,301)
      write(STDOUT,302)nsym
      write(STDOUT,"((a20),x,15I4)")'NUM BASIS  --> ',(nbas(i),i=1,nsym)
      write(STDOUT,"((a20),x,15I4)")'OCC ORB  --> ',(nocorb(i),i=1,nsym)
      if(ieng.eq.1)write(STDOUT,"((a20),x,15I4)")
     :                             'FRZ CORE  --> ',(minorb(i),i=1,nsym)
      if(ieng.eq.1)write(STDOUT,"((a20),x,15I4)")
     :                            'HIGH VIRT  --> ',(maxorb(i),i=1,nsym)

      write(STDOUT,"(//,10x,(a20),2x,(a4),//)")'BASIS TYPE :', method
      if(method.eq.'UB')write(STDOUT,310)alpha1,beta1
      if(method.eq.'WT')write(STDOUT,311)alpha2,beta2,gamma,delta
      if(method.eq.'ET')then
      write(STDOUT,"('ALPHA  -->',3x,15f15.5)"),(alpha3(i),i=1,nsym)
      write(STDOUT,"('BETA   -->',3x,15f15.5)"),(beta3(i),i=1,nsym)
      endif
      if(method.eq.'EE')then
      write(STDOUT,"('EXPLICIT EXPONENTS -->',3x,300f30.10)")
     : (expnts(i),i=1,nbasis)
      endif

      write(STDOUT,303)ichrg
      write(STDOUT,304)maxit,crit
      write(STDOUT,305)jpnt
      write(STDOUT,306)inuc
      write(STDOUT,307)nocc,iele
      write(STDOUT,308)z,amass
      write(STDOUT,309)
      write(STDOUT,9999)

  302 format('NUMBER OF SYMMETRIES = ',i2)
  304 format('MAX ITER = ',i2,8x,'CONV. =',d10.4)
  303 format('CHARGE =',i2)
  305 format('IPRINT = ',i2)
  306 format('NUCLEAR OPTION = ',i2)
  307 format('Occupied orbitals ',i4,4x,'Occupied electrons =',i4)
  308 format('Atomic number ',d10.4,8x,'and mass = ',d10.4)
  309 format('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<',/)
  310 format('ALPHA1 = ',d10.4,8x,'BETA1 =',d10.4)
  311 format('ALPHA2 = ',d10.4,3x,'BETA2 =',d10.4,3x,'GAMMA = ',
     :   d10.4,3x,'DELTA =',d10.4)
 9994 format(//,27x,a,//)
 1001 format(10x,' Slater type of orbitals chosen  ')
 1002 format(10x,' Gaussian type of orbitals chosen')
 1003 format(10x,' Radial distributions are spreaded ')
 1004 format(10x,' Guant term of Breit interaction is taken ')
 1005 format(10x,' Full Breit interaction is taken ')
 1006 format(10x,' Default radial parameters are taken ')
 1007 format(10x,' Orbitals are frozen for active space ')
 1008 format(10x,' Point nucleus is considered          ')
 1009 format(10x,' Fermi charge distribution considered ')
 1015 format(10x,' Uniform charge distribution considered ')
 1010 format(10x,' Neutral atom has been considered     ')
 1011 format(10x,' Singly ion has been considered       ')
 1012 format(10x,' Doubly ion has been considered       ')
 1013 format(10x,' Triply ion has been considered       ')
 1014 format(10x,' Multiply ion has been considered     ')

      RETURN
      END


      SUBROUTINE CHECKNUC(atom,IELE,RRMS,NPROT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC This subroutine CHECKS whether nuclear info correct!           CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Implicit real*8 (a-h,o-z)

      Character*2 yes,y,atom,ATOMS(109)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /nucinf/amass,z

      data ATOMS/' H','HE','LI','BE',' B',' C',' N',' O',' F','NE',
     &'NA','MG','AL','SI',' P',' S','CL','AR',' K','CA','SC','TI',' V',
     &'CR','MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR',
     &'RB','SR',' Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN',
     &'SN','SB','TE',' I','XE','CS','BA','LA','CE','PR','ND','PM','SM',
     &'EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA',' W','RE',
     &'OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN','FR','RA',
     &'AC','TH','PA',' U','NP','PU','AM','CM','BK','CF','ES','FM','MD',
     &'NO','LR','DB','JL','RF','BH','HN','MT'/

      Do 1 Nl=1,109
      If(atom.Eq.ATOMS(Nl)) Goto 2
   1  Continue

      Write(STDOUT,100)atom,Nl
  100 format(//,' ATOM ',a2,i4,'  IS NOT CORRECT')
      Stop

    2 NPROT=Nl

      if(int(z).ne.nprot)then
      Print*,'Wrong proton number'
      stop
      endif

      Write(STDOUT,101)atom,NPROT
      Write(STDOUT,102)amass,z
  101 format(//,' ATOM =  ',a2,7x,' PROTON NUMBER = ',i4)
  102 format(//,' Mass = ',d10.4,7x,' ATOMIC NUMBER = ',d10.4)


      if(ivn.eq.2)then
      if(iele+ichrg+1.ne.nprot)then
      write(STDOUT,*)'>>>>>>>>> PROGRAMME TERMINATES DUE TO IO ERROR'
      write(STDOUT,*)' CHECK ATOMIC NUMBER, CHARGE AND OCCUPANCY'
      stop
      endif
      else
      if(iele+ichrg.ne.nprot)then
      write(STDOUT,*)'>>>>>>>>> PROGRAMME TERMINATES DUE TO IO ERROR'
      write(STDOUT,*)' CHECK ATOMIC NUMBER, CHARGE AND OCCUPANCY'
      stop
      endif
      endif

      If(inuc.eq.1)then

C Converting units
      ainfcm = 0.529177249d-08
      fmtoau = 1.0d-13/ainfcm

      rrms=fmtoau*(0.836d0*Amass**(1.d0/3.d0)+0.57d0)

      Else

      CP=2.2677d0*1.0d-5*Amass**(1.0d0/3.0d0)
      rrms=Dsqrt(3.0d0/5.0d0)*CP

      Endif

      If(inuc.eq.0)then
      Rnt= DEXP(-65.0D0/16.0D0)/Z
      Else
      Rnt=2.0d-06
      Endif

      WRITE(STDOUT,*)"Nuclear R0 = ",RNT, "Rrms =",RRMS

      RETURN
      END


      SUBROUTINE SETGRD
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  This subroutine defines the grids                           CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /nucinf/amass,z

      R(1) = 0.0d0
      Rp(1) = Rnt

      eph = dexp(h)
      ett = 1.d0
c ---------------------------------------------------------------------
c   set up the arrays r, rp, rpor
c ---------------------------------------------------------------------
      do i = 2,N
      ett = eph*ett
      ettm1 = ett-1.0d0
      r(i) = rnt*ettm1
      rp(i) = rnt*ett
      rpor(i) = ett/ettm1
      enddo

      WRITE(STDOUT,*)"Stepsize=",H," Total grids=",N
      WRITE(STDOUT,*)" Maximum radius R(N) =  ",R(N)

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE NUCPARM(RRMS,dltx)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

*********************************************************************
*  This subroutine determines nuclear parameters
*********************************************************************
      implicit real*8 (a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /nucinf/amass,z
     : /propopt/istp,ialpvr

C New common block
      Common/cons/zero,half,tenth,one,two,three,ten
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
      Common/npar/parm(2)

      c0=1.3703599972d002     
      pi = 3.14159265358979324d0 

      zero=0.0d00
      half=0.5d00
      one=1.0d00
      two=2.0d00
      three=3.0d00
      tenth=0.1d0
      ten=10.0d0

      if(ialpvr.eq.0)then
      c=c0
      else
      c=zero
      if(ialpvr.eq.1)c=c0/dsqrt(one+dltx)
      endif

C Converting units
      ainfcm = 0.529177249d-08
      fmtoau = 1.0d-13/ainfcm

C Setting proton parameters
      tfm=2.3d0
      t=tfm*fmtoau
      const= 4.0d0*log(3.0d0)

      parm(2)=t/const

      If(Amass.le.5.0d0)then

      parm(1)=(2.2291d0*Amass**(1.0d0/3.0d0)-0.90676d0)*1.0d-5

      else

      factor=rrms**2-(7.0d00/5.0d00)*(pi**2)*(parm(2)**2)
      parm(1)=dsqrt(5.0d00/3.0d00)*dsqrt(factor)

      endif

      RETURN
      END


      SUBROUTINE NUCPOT(RHO,NPROT)
c ********************************************************************
C It calculates nuclear potential (rho=rho_p+rho_n
c ********************************************************************
      implicit real*8 (a-h, o-z)

      Common/cons/zero,half,tenth,one,two,three,ten
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
      Common/npar/parm(2)
      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /nucinf/amass,z

      Dimension rho(MN)

         DO I = 1,N
            RHO(I) = zero
         ENDDO

      if(inuc.eq.0)then

c-----------------------------------------------------------------------
c                                 Point nucleus
c-----------------------------------------------------------------------

        DO 10 I = 1,N
          RHO(I) = dble(NPROT)
   10   CONTINUE

      endif

      if(inuc.eq.1)then

c-----------------------------------------------------------------------
c                                Fermi distribution nucleus
c-----------------------------------------------------------------------

        CP=PARM(1)
        AP=PARM(2)
        ABC = AP/CP
        TABC = TWO*ABC
        ABC2 = ABC*ABC
        THABC2 = THREE*ABC2
        ABC3 = ABC2*ABC
        CBA = CP/AP
        PI2 = PI*PI
        HPIAC2 = HALF*PI2*ABC2
        SIX = TWO*THREE
        H3 = HALF*THREE
        H3PHP = H3 + HPIAC2
        CALL ES(-CBA,S2MCBA,S3MCBA)
        SABC3 = SIX*ABC3
        DMSAS = -SABC3*S3MCBA
        EN = ONE + ABC2*PI2 + DMSAS
        ZBN = dble(NPROT)/EN

        DO 20 I = 1,N
          RI = R(I)
          RMC = RI - CP
          RMCBA = RMC/AP
          RBC = RI/CP
          IF (RBC.LE.ONE) THEN
          NNUC = I
            CALL ES(RMCBA,S2RCBA,S3RCBA)
            RHO(I) = ZBN* (DMSAS+SABC3*S3RCBA+
     A              RBC* (H3PHP-THABC2*S2RCBA-HALF*RBC*RBC))
          ELSE
            CALL ES(-RMCBA,S2RCBA,S3RCBA)
        RHO(I)=dble(NPROT)*(ONE+THABC2*(RBC*S2RCBA+TABC*S3RCBA)/EN)
          ENDIF
   20   CONTINUE

      WRITE(STDOUT,3030)NNUC
      IF((R(NNUC)-R(NNUC-1)).GT.TWO*TWO*AP)WRITE(STDOUT,3040)

       endif

      if(inuc.eq.2)then

c-----------------------------------------------------------------------
c                              Uniform distribution nucleus
c-----------------------------------------------------------------------

        CP=2.2677d0*1.0d-5*Amass**(1.0d0/3.0d0)

        DO 30 I = 1,N

        RI = R(I)
        RBC = RI/CP

        IF(RBC.LT.ONE)THEN
          NNUC = I
          RHO(I) = RBC*dble(NPROT)*HALF*(THREE-RBC**2)
        ELSE
          RHO(I) = dble(NPROT)
        ENDIF

   30   CONTINUE

      WRITE(STDOUT,3050)NNUC

      endif

 3030 FORMAT (4X,I4,' grid points inside radius C')
 3040 FORMAT (/' ***** WARNING: grid is not sufficiently dense at',
     A       ' nuclear half-charge radius *****')
 3050 FORMAT (4X,I4,' grid points inside radius RNUC')

      RETURN
      END


      SUBROUTINE ES(F,S2F,S3F)
c ********************************************************************
c Evaluate the sum of the series
c
c          k(f)=sum_(0)^(infinity)(-1)^n  exp(n*f)/ (n^k)
c
c       for k = 2, 3 to machine precision. this is a utility
c------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      Common/cons/zero,half,tenth,one,two,three,ten

      N = 0
      S2F = ZERO
      S3F = ZERO
      FASE = ONE
   10 CONTINUE
      N = N + 1
      EN = DBLE(N)
      OBN = ONE/EN
      FASE = -FASE
      ENF = EXP(EN*F)
      TERM2 = FASE*ENF*OBN*OBN
      TERM3 = TERM2*OBN
      S2LAST = S2F
      S2F = S2F + TERM2
      S3F = S3F + TERM3
      IF (ABS(S2F).NE.ABS(S2LAST)) GOTO 10

      RETURN
      END

c ********************************************************************
c SUBROUTINE setqic
c ********************************************************************
c This  subroutine sets up the coefficients for subroutines
c quad, rinti, yzk
c--------------------------------------------------------------------------
      subroutine setqic
      implicit real*8(a-h,o-z)

      logical FIRST
      DIMENSION B13(13,13)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c

C New common block
      Common/LIC13/A(13,13)

      data c5den/720.0d00/

C   Thirteen-point Lagrange interpolation coefficients for first
C   derivative
C
      DATA (B13(1,I),I=1,13)/-1486442880.0D00,5748019200.0D00,
     A     -15807052800.0D00,35126784000.0D00,-59276448000.0D00,
     B     75873853440.0D00,-73766246400.0D00,54195609600.0D00,
     C     -29638224000.0D00,11708928000.0D00,-3161410560.0D00,
     D     522547200.0D00,-39916800.0D00/
      DATA (B13(2,I),I=1,13)/-39916800.0D00,-967524480.0D00,
     A     2634508800.0D00,-4390848000.0D00,6586272000.0D00,
     B     -7903526400.0D00,7376624640.0D00,-5269017600.0D00,
     C     2822688000.0D00,-1097712000.0D00,292723200.0D00,
     D     -47900160.0D00,3628800.0D00/
      DATA (B13(3,I),I=1,13)/3628800.0D00,-87091200.0D00,
     A     -684478080.0D00,1596672000.0D00,-1796256000.0D00,
     B     1916006400.0D00,-1676505600.0D00,1149603840.0D00,
     C     -598752000.0D00,228096000.0D00,-59875200.0D00,9676800.0D00,
     D     -725760.0D00/
      DATA (B13(4,I),I=1,13)/-725760.0D00,13063680.0D00,-143700480.0D00,
     A     -476910720.0D00,1077753600.0D00,-862202880.0D00,
     B     670602240.0D00,-431101440.0D00,215550720.0D00,-79833600.0D00,
     C     20528640.0D00,-3265920.0D00,241920.0D00/
      DATA (B13(5,I),I=1,13)/241920.0D00,-3870720.0D00,31933440.0D00,
     A     -212889600.0D00,-303937920.0D00,766402560.0D00,
     B     -447068160.0D00,255467520.0D00,-119750400.0D00,42577920.0D00,
     C     -10644480.0D00,1658880.0D00,-120960.0D00/
      DATA (B13(6,I),I=1,13)/-120960.0D00,1814400.0D00,-13305600.0D00,
     A     66528000.0D00,-299376000.0D00,-148262400.0D00,558835200.0D00,
     B     -239500800.0D00,99792000.0D00,-33264000.0D00,7983360.0D00,
     C     -1209600.0D00,86400.0D00/
      DATA (B13(7,I),I=1,13)/86400.0D00,-1244160.0D00,8553600.0D00,
     A     -38016000.0D00,128304000.0D00,-410572800.0D00,0.0D00,
     B     410572800.0D00,-128304000.0D00,38016000.0D00,-8553600.0D00,
     C     1244160.0D00,-86400.0D00/
      DATA (B13(8,I),I=1,13)/-86400.0D00,1209600.0D00,-7983360.0D00,
     A     33264000.0D00,-99792000.0D00,239500800.0D00,-558835200.0D00,
     B     148262400.0D00,299376000.0D00,-66528000.0D00,13305600.0D00,
     C     -1814400.0D00,120960.0D00/
      DATA (B13(9,I),I=1,13)/120960.0D00,-1658880.0D00,10644480.0D00,
     A     -42577920.0D00,119750400.0D00,-255467520.0D00,447068160.0D00,
     B     -766402560.0D00,303937920.0D00,212889600.0D00,-31933440.0D00,
     C     3870720.0D00,-241920.0D00/
      DATA (B13(10,I),I=1,13)/-241920.0D00,3265920.0D00,-20528640.0D00,
     A     79833600.0D00,-215550720.0D00,431101440.0D00,-670602240.0D00,
     B     862202880.0D00,-1077753600.0D00,476910720.0D00,
     C     143700480.0D00,-13063680.0D00,725760.0D00/
      DATA (B13(11,I),I=1,13)/725760.0D00,-9676800.0D00,59875200.0D00,
     A     -228096000.0D00,598752000.0D00,-1149603840.0D00,
     B     1676505600.0D00,-1916006400.0D00,1796256000.0D00,
     C     -1596672000.0D00,684478080.0D00,87091200.0D00,-3628800.0D00/
      DATA (B13(12,I),I=1,13)/-3628800.0D00,47900160.0D00,
     A     -292723200.0D00,1097712000.0D00,-2822688000.0D00,
     B     5269017600.0D00,-7376624640.0D00,7903526400.0D00,
     C     -6586272000.0D00,4390848000.0D00,-2634508800.0D00,
     D     967524480.0D00,39916800.0D00/
      DATA (B13(13,I),I=1,13)/39916800.0D00,-522547200.0D00,
     A     3161410560.0D00,-11708928000.0D00,29638224000.0D00,
     B     -54195609600.0D00,73766246400.0D00,-75873853440.0D00,
     C     59276448000.0D00,-35126784000.0D00,15807052800.0D00,
     D     -5748019200.0D00,1486442880.0D00/
C
      DATA B13DEN/479001600.0D00/

      DATA FIRST/.TRUE./

c------------------------------------------------------------------------
c   Five-point newton-cotes coefficients for closed integration.
c   Expressed as  rational numbers
c-------------------------------------------------------------------------
      c5num(1,2)=251.0d00
      c5num(2,2)=646.0d00
      c5num(3,2)=-264.0d00
      c5num(4,2)=106.0d00
      c5num(5,2)=-19.0d00

      c5num(1,3)=232.0d00
      c5num(2,3)=992.0d00
      c5num(3,3)=192.0d00
      c5num(4,3)=32.0d00
      c5num(5,3)=-8.0d00

      c5num(1,4)=243.0d00
      c5num(2,4)=918.0d00
      c5num(3,4)=648.0d00
      c5num(4,4)=378.0d00
      c5num(5,4)=-27.0d00

      c5num(1,5)=224.0d00
      c5num(2,5)=1024.0d00
      c5num(3,5)=384.0d00
      c5num(4,5)=1024.0d00
      c5num(5,5)=224.0d00

      IF (FIRST) THEN
C
C   Thirteen-point coefficients for DPBDT
C
        FACTOR = 1.0d0/B13DEN
        DO 2 J = 1,13
          DO 1 I = 1,13
            A(I,J) = B13(I,J)*FACTOR
   1     CONTINUE
   2   CONTINUE
C
        FIRST = .FALSE.
C
      ENDIF

      factor = h/c5den
      do 5 j = 2,4
         do 4 i = 2,5
            cnc5c(i,j) = factor*c5num(i,j)
    4    continue
    5 continue

      c1 = factor*c5num(1,5)
      c2 = factor*c5num(2,5)
      c3 = factor*c5num(3,5)
      c4 = c1+c1

      do 7 j = 2,5
         do 6 i = 2,5
            c5num(i,j) = c5num(i,j)/c5den
    6    continue
    7 continue

      return
      end


c ********************************************************************
c FUNCTION idbf
c ********************************************************************
c     Calculates factorial
c-----------------------------------------------------------------------
      integer function idbf(m,i)
      implicit real*8(a-h,o-z)

      if(m.lt.-1) then

      write(*,*)'Factorial problem: m.lt.-1',m
      stop

      else

        if(m.eq.0.or.m.eq.-1)then

          if(m.eq.0)idbf=1

          if(m.eq.-1)then
          if(i.eq.2)then
           idbf=1
          else
           write(*,*)'Factorial problem: m.lt.0',m
           stop
          endif
          endif

        else

          idbf=1

          do j=m,1,-i
            idbf=j*idbf
          enddo

        endif

      endif

      return
      end


      SUBROUTINE SETBASIS(Method,cla)
c ********************************************************************
c SUBROUTINE gbasis
c --------------------------------------------------------------------
c     Generate the basis functions on a grid point. It generates
c     (at present) geometrical basis but can be easily transformed
c     to a general one and also to contracted gaussian type.
c --------------------------------------------------------------------

      implicit real*8(a-h,o-z)
      character*2 method

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     :      /UB/alpha1,beta1
     :      /WT/alpha2,beta2,gamma1,delta
     :      /ET/alpha3(MNS),beta3(MNS)
     :      /EE/expnts(MNB),trnopt(MNB)
     :      /jvalue/orbj(MNS),kap(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /basc/nbasis,nocc,mbasis
     :      /skipg/nbasg(MNS),nbasisg,mbasisg,nskipg(MNS)

C New common blocks
      Common/basis/gl(MN,MNBAS),gs(MN,MNBAS)
     :  /quan1/kappe(MNBAS),orbje(MNBAS),kappc(MNBAS),orbjc(MNBAS)
     :  /quan2/lorba(MNS),lorbb(MNS),iqe(MNBAS),iqc(MNBAS)
     :  /skip/nskipe(MNS),nskipc(MNS)
     :  /diffb/dgl(MN,MNBAS)
     :  /stpfac/dg1(MN,MNBAS),dg2(MN,MNBAS)
     :  /splflags/MRCC,IBASOPT

      dimension rnk(MN),NAK(15),nkk(15),alpha(MNBAS,MNS)
      dimension llb(15)

      data nak/-1,1,-2,2,-3,3,-4,4,-5,5,-6,6,-7,7,-8/
      data nkk/1,2,2,3,3,4,4,5,5,6,6,7,7,8,8/
      data llb/1,0,2,1,3,2,4,3,5,4,6,5,7,6,8/

      four = two*two
      onehalf = one + half
      sq2=dsqrt(two)
      spi=dsqrt(pi)

      nbasisg=0

      nskipe(1)=0

      do isym=1,nsym
      lorba(isym)=nkk(isym)-1
      lorbb(isym)=llb(isym)
      If(nak(isym).ne.kap(isym))then
      Print*,'Wrong in kappa values'
       stop
      EndIf
      nskipe(isym+1)=nskipe(isym)+nbas(isym)
      enddo

      ii=0
      do isym=1,nsym
      do ibas=1,nbas(isym)
      ii=ii+1
      kappe(ii)=kap(isym)
      orbje(ii)=orbj(isym)

      If(mod(isym,2).eq.0)then
      iqe(ii)=-1
      else
      iqe(ii)=1
      EndIf

      If(iqe(ii).ne.idint(2*(orbje(ii)-lorba(isym))))then
      Print*,'Wrong in assignment of a quantum numbers'
       stop
      EndIf

      If(kappe(ii).ne.-iqe(ii)*idint(Orbje(ii)+Half))then
      Print*,'Wrong in assignment of quantum numbers e'
       stop
      EndIf

      enddo
      enddo

      If(ii.ne.nbasis.or.2*ii.ne.mbasis)then
      Print*,'Wrong in counting total basis number'
      Stop
      Endif

      ii=0
      nskipc(1)=0
      do isym=1,nsym
      nskipc(isym+1)=nskipc(isym)+nocorb(isym)
      nn=nocorb(isym)
      if(nn.ne.0)then
      do ibas=1,nocorb(isym)
      ii=ii+1
      kappc(ii)=kap(isym)
      orbjc(ii)=orbj(isym)

      If(mod(isym,2).eq.0)then
      iqc(ii)=-1
      else
      iqc(ii)=1
      EndIf

      If(kappc(ii).ne.-iqc(ii)*idint(Orbjc(ii)+Half))then
      Print*,'Wrong in assignment of quantum numbers c'
       stop
      EndIf

      enddo
      endif

      enddo

      If(ii.ne.nocc)then
      Print*,'Wrong in counting total occupied orbitals'
      Stop
      Endif

      ii=0
      nskipg(1)=0
      do isym=1,nsym
      nskipg(isym+1)=nskipg(isym)+(nbas(isym)-nbasg(isym))

      if(ivn.eq.1)then
      if(nskipg(isym+1).ne.nskipe(isym+1)-nskipc(isym+1))then
      Print*,'Skiping problem for V^N or hybridization' 
      Stop
      endif
      endif

      nn=nbasg(isym)
      if(nn.ne.0)then
      do ibas=1,nbasg(isym)
      ii=ii+1
      enddo
      endif
      enddo

      if(ivn.eq.1)then
      if(ii.ne.nocc)print*,'V^N counting prob'
      stop
      endif

      nbasisg=nbasis-ii
      mbasisg=2*nbasisg

      write(STDOUT,200)nbasis
      write(STDOUT,201)

      do isym=1,nsym
      write(STDOUT,202)isym,nocorb(isym),nskipc(isym)
      enddo

      write(STDOUT,203)
      do isym=1,nsym
      write(STDOUT,204)isym,nbas(isym),nskipe(isym)
      enddo

      if(method.eq.'UB')iopt=1
      if(method.eq.'WT')iopt=2
      if(method.eq.'ET')iopt=3
      if(method.eq.'EE')iopt=4

      do ib=1,n
      rnk(ib)=zero 
      enddo

      do isym=1,nsym

      do ib=1,n
      rnk(ib) = r(ib)**(lorba(isym)+1)
      enddo

      do ic=1,nbas(isym)
      jc=ic+nskipe(isym)

c----------------------------------------------------------------
C Here the alpha is calculated using alpha(i) and beta(i) for the
C three different types of basis given by UB,WT and ET. For EE
C the exponents are used from outside.
c----------------------------------------------------------------

      goto (1,2,3,4),iopt

   1  alpha(jc,isym) = alpha1*(beta1**(ic-1))

      if(ispd.eq.1)alpha(jc,isym)=alpha1*
     :(beta1**(cla*ic-1)+one/beta1**(ic))

      goto 5

   2  alpha(jc,isym) = alpha2*(beta2**(ic-1))*(1+ gamma1 *
     :               (ic/nbas(isym))**delta)
      goto 5

   3  alpha(jc,isym) = alpha3(isym)*(beta3(isym)**(ic-1))

      if(ispd.eq.1)alpha(jc,isym)=alpha3(isym)*
     :(beta3(isym)**(cla*ic-1)+one/beta3(isym)**(ic))

      goto 5

   4  alpha(jc,isym) = expnts(jc)

   5  continue
      enddo

      if(isgo.eq.1)then

      dkap=dfloat(kap(isym))
      dnkk=dfloat(lorba(isym)+1)
      ifact=idbf(2*(lorba(isym)+1),1)
      power=(two*dnkk+one)

      do id=1,nbas(isym)
      jd=id+nskipe(isym)

      cnl=dsqrt(((two*alpha(jd,isym))**(power))/dfloat(ifact))

      fact1=two*(alpha(jd,isym)**2)*((dnkk+dkap)**2)/
     :(dnkk*(two*dnkk-one))
      fact2=alpha(jd,isym)*alpha(jd,isym)
      fact3=-two*(alpha(jd,isym)**2)*(dnkk+dkap)/dnkk
      cns=dsqrt(one/(fact1+fact2+fact3))

      gl(1,jd)=0
      gs(1,jd)=0

      do ie = 2,n

      expon=dexp(-alpha(jd,isym)*r(ie))*rnk(ie)

      gl(ie,jd)=cnl*expon
      gs(ie,jd)=cns*gl(ie,jd)*((dnkk+dkap)/r(ie)-
     :                   alpha(jd,isym))

      dgl(ie,jd)=gl(ie,jd)*((dnkk+dkap)/r(ie)-
     :                   alpha(jd,isym))

      enddo

      enddo

      endif

      if(isgo.eq.2)then

      dkap=dfloat(kap(isym))
      dnkk=dfloat(lorba(isym)+1)
      ifact=idbf(2*lorba(isym)+1,2)
      power=(two*dnkk+one)/two

      do id=1,nbas(isym)
      jd=id+nskipe(isym)

      cnl=dsqrt(two**(two*dnkk+onehalf)*(alpha(jd,isym)**power)/
     :(ifact*spi))

      fact1=four*alpha(jd,isym)*(dnkk+dkap)**2/(two*dnkk-one)
      fact2=(two*dnkk+one)*alpha(jd,isym)
      fact3=-four*alpha(jd,isym)*(dnkk+dkap)
      cns=dsqrt(one/(fact1+fact2+fact3))

      gl(1,jd)=0
      gs(1,jd)=0
      do ie = 2,n
      expon= dexp(-alpha(jd,isym)*r(ie)*r(ie))*rnk(ie)

      gl(ie,jd)=cnl*expon
      gs(ie,jd)=cns*gl(ie,jd)*((dnkk+dkap)/r(ie)-
     :                   two*alpha(jd,isym)*r(ie))

      dgl(ie,jd)=gl(ie,jd)*((dnkk+dkap)/r(ie)-
     :                   two*alpha(jd,isym)*r(ie))

      dg1(ie,jd)=(dnkk/r(ie))-(two*alpha(jd,isym)*r(ie))
      dg2(ie,jd)=-(((dnkk+dkap)/r(ie)**2)+two*alpha(jd,isym))/
     :((dnkk+dkap)/r(ie)-two*alpha(jd,isym)*r(ie))

      enddo

      enddo

      endif

      if(isgo.eq.3)then

      dkap=dfloat(kap(isym))
      dnkk=dfloat(lorba(isym)+1)
      power=(two*dnkk+one)/(two*four)
      power1=(two*dnkk-one)/(two*four)

      do id=1,nbas(isym)
      jd=id+nskipe(isym)

      cnl=two*(two*alpha(jd,isym))**power/dsqrt(gml(lorba(isym)+1))

      fact1=(dnkk+dkap)**2*gml(lorba(isym))/four
      fact2=-(dnkk+dkap)*gml(lorba(isym)+2)
      fact3=(2*lorba(isym)+9)*gml(lorba(isym)+2)/four

      cns=(two*alpha(jd,isym))**power1/(cnl*dsqrt(fact1+fact2+fact3))

      gl(1,jd)=0
      gs(1,jd)=0
      do ie = 2,n
      expon= rnk(ie)*dexp(-alpha(jd,isym)*r(ie)**4)

      gl(ie,jd)=cnl*expon
      gs(ie,jd)=cns*gl(ie,jd)*((dnkk+dkap)/r(ie)-
     :four*alpha(jd,isym)*r(ie)**3)

      dgl(ie,jd)=gl(ie,jd)*((dnkk+dkap)/r(ie)-
     :four*alpha(jd,isym)*r(ie)**3)

      enddo

      enddo

      endif

      if(isgo.eq.4)then

      dkap=dfloat(kap(isym))
      dnkk=dfloat(lorba(isym)+1)
      power=(four*dnkk+one)/(two*four)
      power1=(four*dnkk-one)/(two*four)

      do id=1,nbas(isym)
      jd=id+nskipe(isym)

      cnl=two*(two*alpha(jd,isym))**power/dsqrt(gmm(lorba(isym)+1))

      fact1=(two*dnkk+dkap)**2*gmn(lorba(isym)+1)/four
      fact2=-(two*dnkk+dkap)*gmn(lorba(isym)+2)
      fact3=(4*lorba(isym)+11)*gmn(lorba(isym)+2)/four

      cns=(two*alpha(jd,isym))**power1/(cnl*dsqrt(fact1+fact2+fact3))

      gl(1,jd)=0
      gs(1,jd)=0
      do ie = 2,n
      expon= rnk(ie)**2*dexp(-alpha(jd,isym)*r(ie)**4)

      gl(ie,jd)=cnl*expon
      gs(ie,jd)=cns*gl(ie,jd)*((two*dnkk+dkap)/r(ie)-
     :four*alpha(jd,isym)*r(ie)**3)

      dgl(ie,jd)=gl(ie,jd)*((two*dnkk+dkap)/r(ie)-
     :four*alpha(jd,isym)*r(ie)**3)

      enddo

      enddo

      endif

      if(isgo.eq.5)then

      dkap=dfloat(kap(isym))
      dnkk=dfloat(lorba(isym)+1)
      power=(two*dnkk+one)/dnkk
      power1=(two*dnkk-one)/(two*dnkk)

      do id=1,nbas(isym)
      jd=id+nskipe(isym)

      cnl=dsqrt(dnkk*(two*alpha(jd,isym))**power/gm2p1(lorba(isym)+1))

      fact1=(dnkk+dkap)**2*gm2n1(lorba(isym)+1)/dnkk
      fact2=-(dnkk+dkap)*gm3n1(lorba(isym)+1)
      fact3=dnkk*gm4n1(lorba(isym)+1)/four

      cns=(two*alpha(jd,isym))**power1/(cnl*dsqrt(fact1+fact2+fact3))

      gl(1,jd)=0
      gs(1,jd)=0
      do ie = 2,n
      expon= rnk(ie)*dexp(-alpha(jd,isym)*r(ie)**(lorba(isym)+1))

      gl(ie,jd)=cnl*expon
      gs(ie,jd)=cns*gl(ie,jd)*((dnkk+dkap)/r(ie)-
     :dnkk*alpha(jd,isym)*r(ie)**lorba(isym))

      dgl(ie,jd)=gl(ie,jd)*((dnkk+dkap)/r(ie)-
     :dnkk*alpha(jd,isym)*r(ie)**lorba(isym))

      enddo

      enddo

      endif

      enddo

      IF(IBASOPT.and.METHOD.ne.'EE')THEN
C     PRINTING OUT THE EXPONENTS TO A SEPERATE FILE
           do isym=1,nsym
              do j=1,nbas(isym)
              jck=j+nskipe(isym)
              trnopt(jck)=0.0
              expnts(jck)=alpha(jck,isym)
              enddo
           enddo

      write(STDXP,"(16I4)")nsym,(nbas(isym),isym=1,nsym)
      write(STDXP,"(F30.12,2x,F4.0)")(expnts(i),trnopt(i),i=1,nbasis)
      ENDIF

 200  format(/,2x,' TOTAL NUMBER OF BASIS FUNCTION=',i4,//)
 201  format(2x,'SYMMETRY   NO. OF OCC ORBITS   NSKIPC',/)
 202  format(2x,i4,12x,i4,12x,i4)
 203  format(2x,'SYMMETRY   NO. OF BASIS FUNCTION   NSKIPE',/)
 204  format(2x,i4,12x,i4,12x,i4)

      RETURN
      END


      SUBROUTINE SETMATRIX(RHO)
c ********************************************************************
C  This subroutine sets the diagonalization subroutine
c --------------------------------------------------------------------

      implicit real*8(a-h,o-z)

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     :      /jvalue/orbj(MNS),kap(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /basc/nbasis,nocc,mbasis
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS)
     :  /skip/nskipe(MNS),nskipc(MNS)
     :  /diffb/dgl(MN,MNBAS)

C New common blocks
      Common/dfm/emt(MNB,MNB),df_single(MNB,MNB),df_two(MNB,MNB)
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
c ---------------------------------------------------------------
c Local dimension
c ---------------------------------------------------------------
      dimension atemp(MNB,MNB),rho(MN)

      spi = dsqrt(pi)
      four=two*two
      sq2=dsqrt(two)

      do i=1,MNB
      do j=1,MNB
      emt(i,j)=zero
      df_single(i,j)=zero
      df_two(i,j)=zero
      atemp(i,j)=zero
      enddo
      enddo

      do i=1,n
      tmp1(i)=zero
      tmp2(i)=zero
      tmp3(i)=zero
      enddo

      do isym=1,nsym
      do ia=1,nbas(isym)
      ja=ia+nskipe(isym)
      do ib=1,nbas(isym)
      jb=ib+nskipe(isym)

      ka=ia+2*nskipe(isym)
      index1=ia+nbas(isym)+2*nskipe(isym)

      kb=ib+2*nskipe(isym)
      index2=ib+nbas(isym)+2*nskipe(isym)

      do la = 1,n
      tmp1(la)=gl(la,ja)*gl(la,jb)*rp(la)
      tmp2(la)=gs(la,ja)*gs(la,jb)*rp(la)
      tmp3(la)=gs(la,jb)*dgl(la,ja)*rp(la)
      enddo
      call valint(tmp1,result1)
      call valint(tmp2,result2)
      call valint(tmp3,result3)

      emt(ka,kb)=result1
      emt(index1,index2)=result2
      atemp(ja,jb)=c*result3

      enddo
      enddo

c ---------------------------------------
c     matrix df_single
c --------------------------------------
      do ia=1,nbas(isym)
      ka=ia+nskipe(isym)
      ja=ia+2*nskipe(isym)
      index1=ia+nbas(isym)+2*nskipe(isym)
      do ib= 1,nbas(isym)
      kb=ib+nskipe(isym)
      jb=ib+2*nskipe(isym)
      index2=ib+nbas(isym)+2*nskipe(isym)

      df_single(ja,index2)=df_single(ja,index2)+atemp(ka,kb)
      df_single(index2,ja)=df_single(ja,index2)

c -----------------------------------------
c     nuclear repulsion energy
c ----------------------------------------

      do la = 1,n
      tmp1(la)=gl(la,ka)*gl(la,kb)*rpor(la)*rho(la)
      tmp2(la)=gs(la,ka)*gs(la,kb)*rpor(la)*rho(la)
      enddo
      call valint(tmp1,result1)
      call valint(tmp2,result2)

      df_single(ja,jb)=df_single(ja,jb)-result1
      df_single(index1,index2)=df_single(index1,index2)
     :   -(result2+two*c*c*emt(index1,index2))

      enddo
      enddo
      enddo

c      if(jpnt.ne.0)then
c      call matout(df_single,MNB,mbasis,mbasis,mbasis,4hH   )
c      endif

      RETURN
      END
      

      SUBROUTINE VALINT(fx,res)
c ******************************************************************
c    This is the integration routine using numerical formula       *
c Simpson's formular truncate at h^4, but here it is at h^13       *
c ******************************************************************


      implicit real*8(a-h,o-z)
      dimension fx(MN)
      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten

      res=zero

      if((n+10).ge.MN)then
      print*,'Dimension of grid size not sufficient'
      stop
      endif

      Do i=n+1,n+10
      fx(i)=zero
      EndDo

      Do i=1,n+10,10
      res=res+(16067*(fx(i)+fx(i+10))+427368*fx(i+5)
     :+106300*(fx(i+1)+fx(i+9))-48525*(fx(i+2)+fx(i+8))
     :+272400*(fx(i+3)+fx(i+7))-260550*(fx(i+4)+fx(i+6)))
      EndDo

      res=5.0d0*h*(res/299376.0d0)

      return
      end

c***************************************************************
      subroutine eigen(a,b,nsub,valu,msub,anorm,nnmax)
      implicit real*8 (a-h,o-z)

c      parameter(nmax=1600)
c     eigenvalues and eigenvectors of a real symmetric matrix
c     a ----- secular determinant (input)
c             eigenvetors  column-wise (output)
c             a is destroyed
c     b ----- nsub x nsub array for intermidiate storage
c     nsub ----- order of a (input)
c                must be equal to or greater than 1
c     valu ----- eigenvalues in decreasing order (output)
c     msub -----  number of eigenvalues (input)
c                 if equal to -1, eigenvectors are not computed
c     anorm ----- matrix norm which is a sqrt of sum of a(i,j)**2
c                                        over all i & j (output)
c     nmax ----- row dimension of a & b
c

c      Dimension a(MNB,MNB), b(MNB,MNB), valu(MNB),VALI(MNB)
c      Dimension EIVEC(MNB,MNB),SCAL(MNB), VALR(MNB),D(MNB)
      real*8 MAT(MNB,MNB)
      integer indx(MNB)

      real*8 a(MNB,MNB), b(MNB,MNB), valu(MNB)
      real*8  anorm
      real*8  EIVEC(MNB,MNB)
      real*8  SCAL(MNB), VALR(MNB)
      real*8  D(MNB), VALI(MNB)



      integer out,BASIS, LOW, HIGH, CNT
      logical VEC, ORTHO, EVNORM

      VEC = .TRUE.                 ! Calculate eigenvectors
      ORTHO = .TRUE.               ! Use routine `ORTHES()' to trans-
                                   ! form `Atilde' into Hessenb. form
      EVNORM = .FALSE.              ! Normalize eigenvectors
      BASIS = 2                    ! Radix

      anorm = 0.0d0
      do i = 1, nsub
         do j = 1, nsub
            MAT(i,j) = dble(a(i,j))
            anorm = anorm + a(i,j)*a(i,j)
         enddo
      enddo
      anorm = dsqrt(anorm)


      out = IEIGEN (VEC,ORTHO,EVNORM,BASIS,nsub,MAT,
     &                        SCAL,D,EIVEC,VALR,VALI,CNT,LOW,HIGH)



      if (out .eq. 0) goto 10          ! Everything ok
      if (out .eq. 401) then
         print*,'eigen(): Order of matrix < 1'
         print*,'execution terminated'
         stop
      endif
      if (out .eq. 402) then
         print*,'eigen(): Matrix is zero'
         print*,'execution terminated'
         stop
      endif
      if (out .eq. 403) then
         print*,'eigen(): Max. no. of iterations in hqr2() exceeded'
         print*,'execution terminated'
      endif
 10   continue



      call indexx(nsub,VALR,indx)    ! sort eigenvalues
      sum = 0.0d0
      do i = 1, nsub
         valu(i) = VALR(indx(nsub+1-i))
         do j = 1, nsub
            a(i,j) = EIVEC(i,indx(nsub+1-j))
         enddo
         sum = sum + dabs(VALI(i))
      enddo

      if (sum .gt. 1.0e-10) then
         print*,'eigen(): Oops! Eigenvalues not real'
         print*,'||VALI|| = ', sum
         print*,'Ignore imaginary part'
      endif

      return
      end






!F 7.8 Eigenwerte und Eigenvektoren einer Matrix nach den Verfahren von Martin,
!Peters, Reinsch und Wilkinson


      INTEGER FUNCTION IEIGEN (VEC,ORTHO,EVNORM,BASIS,N,MAT,
     &                        SCAL,D,EIVEC,VALR,VALI,CNT,LOW,HIGH)

!*****************************************************************
!                                                                *
!     Dieses FUNCTION-Unterprogramm vom Typ INTEGER berechnet    *
!     alle Eigenwerte und Eigenvektoren einer beliebigen         *
!     reellen Matrix MAT.                                        *
!     Die Eigenwerte werden in den Feldern VALR(1:N)             *
!     (Realteil) und VALI(1:N) (Imaginärteil) und (bei gesetzter *
!     VEC-Flagge) die Eigenvektoren im Feld EIVEC(1:N,1:N)       *
!     abgespeichert.                                             *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     VEC:         Parameter für Eigenvektoren                   *
!                     = .TRUE.:  Bestimmung der Eigenvektoren    *
!                     = .FALSE.: Nur Berechnung der Eigenwerte   *
!     ORTHO        Flagge, die anzeigt, ob die Transformation    *
!                  von MAT auf Hessenbergform durch orthogonale  *
!                  Transformation per ORTHES (gesetzt) oder      *
!                  durch elementare Transformationen per         *
!                  ELMHES (gelöscht) vorgenommen werden soll.    *
!                  Die bei der orthogonalen Transformation       *
!                  verwendeten Householder-Matrizen haben den    *
!                  Vorteil, daß sie eine eventuelle Symmetrie    *
!                  der Ausgangsmatrix nicht zerstören.           *
!     EVNORM:      Flagge, die anzeigt, ob die Eigenvektoren     *
!                  normiert werden sollen (gesetzt) oder nicht   *
!                  (gelöscht)                                    *
!     BASIS:       die Basis der Gleitkommadarstellung des       *
!                  verwendeten Rechners (meistens 2 oder 16)     *
!     N:           Ordnung der Matrix MAT                        *
!     MAT:         ein (1:N,1:N)-Feld vom Typ DOUBLE PRECISION,  *
!                  das die Matrix enthält, deren Eigenwerte      *
!                  und Eigenvektoren berechnet werden sollen     *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     MAT:         der obere Teil dieses (1:N,1:N)-Feldes        *
!                  enthält die Eigenvektoren der                 *
!                  Fast-Dreiecksmatrix, die vom QR-Verfahren     *
!                  erzeugt wird.                                 *
!     SCAL:        ein (1:N)-Feld vom Typ DOUBLE PRECISION, das  *
!                  die Information über die durchgeführten       *
!                  Vertauschungen und die Skalierungsfaktoren    *
!                  enthält.                                      *
!     D:           (1:N)-Vektor mit Umformungsinformationen aus  *
!                  ORTHES                                        *
!     EIVEC:       ein (1:N,1:N)-Feld vom Typ DOUBLE PRECISION,  *
!                  das bei gesetzter VEC-Flagge die              *
!                  normalisierten Eigenvektoren der              *
!                  ursprünglichen vollen Matrix aufnimmt.        *
!                  Falls der i. Eigenwert reell ist, so ist      *
!                  die i. Spalte von EIVEC der dazugehörige      *
!                  reelle Eigenvektor. Falls die Eigenwerte i    *
!                  i+1 ein komplexes Paar bilden, geben i. und   *
!                  (i+1). Spalte Real- und Imaginärteil          *
!                  desjenigen Eigenvektors an, der zu dem        *
!                  Eigenwert mit positivem Imaginärteil gehört.  *
!     VALR,VALI:   zwei (1:N)-Felder vom Typ DOUBLE PRECISION,   *
!                  die Real- und Imaginärteil der Eigenwerte     *
!                  aufnehmen                                     *
!     CNT:         ein (1:N)-Feld vom Typ INTEGER, das die Zahl  *
!                  der Iterationsschritte für jeden Eigenwert    *
!                  aufnimmt. Falls zwei Eigenwerte als Paar      *
!                  gleichzeitig gefunden werden, wird die Zahl   *
!                  der Iterationsschritte mit einem positiven    *
!                  Vorzeichen für den ersten und einem negativen *
!                  Vorzeichen für den zweiten Eigenwert          *
!                  eingetragen.                                  *
!     LOW,HIGH:    Die Zeilen von 1 bis LOW-1 bzw. die Zeilen    *
!                  von HIGH+1 bis N enthalten die isolierten     *
!                  Eigenwerte, d. h. Eigenwerte zu Vielfachen    *
!                  von Einheitsvektoren.                         *
!                                                                *
!     RÜCKGABEWERT:                                              *
!     =============                                              *
!     0:      kein Fehler                                        *
!     401:    Die Ordnung N der Matrix MAT ist kleiner als 1.    *
!     402:    MAT ist die Nullmatrix.                            *
!     403:    Die in HQR2 definierte maximale Schrittzahl fuer   *
!             das QR-Verfahren wurde ueberschritten, ohne dass   *
!             alle Eigenwerte berechnet werden konnten.          *
!                                                                *
!     LOKALE NAMEN:                                              *
!     =============                                              *
!     ONE,TWO,HALF: Gleitkommakonstanten                         *
!     EPS:          Maschinengenauigkeit                         *
!     TEMP:         Hilfsvariable                                *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  Benötigte Unterprogramme: BALAN, ELMHES, ORTHES, ELMTRA,      *
!                            ORTTRA, HQR2, BALBAK, NORMAL, SWAP, *
!                            COMDIV, COMABS                      *
!                                                                *
!                                                                *
!  Quellen : 1. Martin, R. S. und Wilkinson, J. H.,              *
!               siehe [MART68].                                  *
!            2. Parlett, B. N. und Reinsch, C., siehe [PARL69].  *
!            3. Peters, G. und Wilkinson, J. H., siehe [PETE70]. *
!                                                                *
!*****************************************************************
!                                                                *
!  Autor     : Jürgen Dietel                                     *
!  Datum     : 15.7.1993                                         *
!  Quellcode : FORTRAN 77                                        *
!                                                                *
!*****************************************************************

      implicit real*8 (a-h,o-z)

      PARAMETER        (ONE = 1.0D0,TWO = 2.0D0,HALF = 0.5D0)
      integer BASIS, N, LOW, HIGH, CNT
      logical VEC,ORTHO,EVNORM
      real*8 MAT(MNB,MNB)
      Dimension EIVEC(MNB,MNB),SCAL(MNB), VALR(MNB)
      Dimension D(MNB), VALI(MNB)
      INTEGER RES, BALAN, BALBAK, ELMHES, ELMTRA, HQR2
      INTEGER ORTHES, ORTTRA

!     Berechnung der Maschinengenauigkeit EPS (d. h. der kleinsten
!     positiven Maschinenzahl, für die auf dem Rechner gilt:
!     1 + EPS > 1):

      TEMP = TWO
      EPS = ONE

   10 IF (ONE .LT. TEMP) THEN
         EPS = EPS * HALF
         TEMP = ONE + EPS
         GOTO 10
         ENDIF
      EPS = TWO * EPS
      RES = BALAN(N,MAT,SCAL,LOW,HIGH,BASIS)
      IF (RES .NE. 0) THEN
         IEIGEN = RES + 100
         RETURN
      ENDIF
      IF (ORTHO) THEN
         RES = ORTHES(N,LOW,HIGH,MAT,D,EPS)
      ELSE
         RES = ELMHES(N,LOW,HIGH,MAT,CNT)
      ENDIF
      IF (RES .NE. 0) THEN
         IEIGEN = RES + 200
         RETURN
      ENDIF
      IF (VEC) THEN
         IF (ORTHO) THEN
            RES = ORTTRA(N,LOW,HIGH,MAT,D,EIVEC)
         ELSE
            RES = ELMTRA(N,LOW,HIGH,MAT,CNT,EIVEC)
         ENDIF
         IF (RES .NE. 0) THEN
            IEIGEN = RES + 300
            RETURN
         ENDIF
      ENDIF
      RES = HQR2(VEC,N,LOW,HIGH,MAT,VALR,VALI,EIVEC,CNT,EPS)
      IF (RES .NE. 0) THEN
         IEIGEN = RES + 400
         RETURN
      ENDIF
      IF (VEC) THEN
         RES = BALBAK(N,LOW,HIGH,SCAL,EIVEC)
         IF (RES .NE. 0) THEN
            IEIGEN = RES + 500
            RETURN
         ENDIF
         IF (EVNORM) THEN
            RES = NORMAL(N,EIVEC,VALI)
            IF (RES .NE. 0) THEN
               IEIGEN = RES + 600
               RETURN
            ENDIF
         ENDIF
      ENDIF
      IEIGEN = 0
      END



      INTEGER FUNCTION BALAN (N,MAT,SCAL,LOW,HIGH,BASIS)

!*****************************************************************
!                                                                *
!     BALAN balanciert eine gegebene reelle Matrix in der 1-Norm *
!     aus.                                                       *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     N:        die Ordnung der gegebenen Matrix                 *
!     MAT:      ein (1:N,1:N)-Feld, das die Komponenten der      *
!               gegebenen Matrix enthält                         *
!     BASIS:    die Basis der Gleitkommadarstellung auf          *
!               der Maschine                                     *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     MAT:      die ausbalancierte Matrix                        *
!     LOW,HIGH: zwei Integerzahlen, für die MAT(I,J) gleich Null *
!               ist, falls gilt:                                 *
!               1. I>J und                                       *
!               2. J=1,...LOW-1 oder I=HIGH+1,...N               *
!     SCAL:     ein (1:N)-Feld, das die Informationen über       *
!               die durchgeführten Vertauschungen und die        *
!               Skalierungsfaktoren enthält.                     *
!                                                                *
!     RÜCKGABEWERT:                                              *
!     =============                                              *
!     0:        kein Fehler                                      *
!                                                                *
!     LOKALE NAMEN:                                              *
!     =============                                              *
!     ZERO,ONE,PT95: Gleitkommakonstanten                        *
!     I,J,K,L:       Zählvariablen                               *
!     B2:            Quadrat von Basis                           *
!     R,C,F,G,S:     Hilfsvariablen zur Berechnung von           *
!                    Zeilennormen, ihren Kehrwerten und          *
!                    Ähnlichem                                   *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  Benötigte Unterprogramme: SWAP                                *
!                                                                *
!                                                                *
!  Quellen : Parlett, B. N. und Reinsch, c., siehe [PARL69].     *
!                                                                *
!*****************************************************************
!                                                                *
!  Autor     : Jürgen Dietel                                     *
!  Datum     : 10.04.1987                                        *
!  Quellcode : FORTRAN 77                                        *
!                                                                *
!*****************************************************************


      implicit real*8 (a-h,o-z)

      integer BASIS, N, LOW, HIGH
      real*8 MAT(MNB, MNB)
      Dimension SCAL(MNB)


      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0,PT95 = 0.95D0)
      INTEGER          I,J,K,L,B2

!     die Norm von MAT(1:N,1:N) reduzieren durch exakte Ähnlich-
!     keitstransformationen, die in SCAL(1:N) abgespeichert werden

      B2 = BASIS*BASIS
      L = 1
      K = N

!     nach Zeilen mit einem isolierten Eigenwert suchen und sie
!     nach unten schieben

      do i = 1, N
         scal(i) = 0
      enddo

   10 DO 50 J=K,1,-1
         R = ZERO
         DO  I=1,K
            IF (I .NE. J) R = R+ABS(MAT(J,I))
         ENDDO
         IF (R .EQ. ZERO) THEN
            SCAL(K) = J
            IF (J .NE. K) THEN
               DO  I=1,K
                  CALL SWAP(MAT(I,J),MAT(I,K))
               ENDDO
               DO  I=L,N
                  CALL SWAP(MAT(J,I),MAT(K,I))
               ENDDO
            ENDIF
            K = K-1
            GOTO 10
         ENDIF
   50    CONTINUE

!     nach Spalten mit einem isolierten Eigenwert suchen und sie
!     nach links schieben

   60 DO 100 J=L,K
         C = ZERO
         DO  I=L,K
            IF (I .NE. J) C = C+ABS(MAT(I,J))
         ENDDO
         IF (C .EQ. ZERO) THEN
            SCAL(L) = J
            IF (J .NE. L) THEN
               DO  I=1,K
                  CALL SWAP(MAT(I,J),MAT(I,L))
               ENDDO
               DO  I=L,N
                  CALL SWAP(MAT(J,I),MAT(L,I))
               ENDDO
            ENDIF
            L = L+1
            GOTO 60
         ENDIF
  100    CONTINUE

!     nun die Teilmatrix in den Zeilen L bis K ausbalancieren

      LOW = L
      HIGH = K
      DO  I=L,K
         SCAL(I) = ONE
      ENDDO
  120 DO 180 I=L,K
         C = ZERO
         R = ZERO
         DO 130 J=L,K
            IF (J .NE. I) THEN
               C = C+ABS(MAT(J,I))
               R = R+ABS(MAT(I,J))
            ENDIF
  130       CONTINUE
         G = R/BASIS
         F = ONE
         S = C+R
  140    IF (C .LT. G) THEN
            F = F*BASIS
            C = C*B2
            GOTO 140
         ENDIF
         G = R*BASIS
  150    IF (C .GE. G) THEN
            F = F/BASIS
            C = C/B2
            GOTO 150
         ENDIF
         IF ((C+R)/F .LT. PT95*S) THEN
            G = ONE/F
            SCAL(I) = SCAL(I)*F
            DO  J=L,N
               MAT(I,J) = MAT(I,J)*G
            ENDDO
            DO  J=1,K
               MAT(J,I) = MAT(J,I)*F
            ENDDO
            GOTO 120
         ENDIF
  180    CONTINUE
      BALAN = 0
      END



      INTEGER FUNCTION BALBAK (N,LOW,HIGH,SCAL,EIVEC)

!*****************************************************************
!                                                                *
!     BALBAK führt eine Rücktransformation aller Rechtseigen-    *
!     vektoren einer ausbalancierten Matrix in die Eigenvektoren *
!     der Originalmatrix durch, von der die balancierte Matrix   *
!     durch Aufruf der Prozedur BALAN abgeleitet wurde.          *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     N:        die Ordnung der Eigenvektoren (Zahl der          *
!               Komponenten)                                     *
!     LOW,HIGH: zwei Integerzahlen, die von der Prozedur         *
!               BALAN stammen                                    *
!     SCAL:     Ausgabevektor der Prozedur BALAN                 *
!     EIVEC:    ein (1:N,1:N)-Feld, von dem jede Spalte einen    *
!               Eigenvektor (oder seinen Realteil oder seinen    *
!               Imaginärteil) der ausbalancierten Matrix         *
!               darstellt                                        *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     EIVEC:    die entsprechenden Eigenvektoren (oder Realteile *
!               oder Imaginärteile) der ursprünglichen Matrix    *
!                                                                *
!     RÜCKGABEWERT:                                              *
!     =============                                              *
!     0:        kein Fehler                                      *
!                                                                *
!     LOKALE NAMEN:                                              *
!     =============                                              *
!     I,J,K:    Hilfsvariablen zur Indexbildung                  *
!     S:        Skalierungswert                                  *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  Benötigte Unterprogramme: SWAP                                *
!                                                                *
!                                                                *
!  Quellen : Parlett, B. N. und Reinsch, C., siehe [PARL69].     *
!                                                                *
!*****************************************************************
!                                                                *
!  Autor     : Jürgen Dietel                                     *
!  Datum     : 10.04.1987                                        *
!  Quellcode : FORTRAN 77                                        *
!                                                                *
!*****************************************************************

      implicit real*8 (a-h,o-z)


      integer N, LOW, HIGH
      Dimension EIVEC(MNB,MNB),SCAL(MNB)
      INTEGER          I,J,K

      DO 20 I=LOW,HIGH
         S = SCAL(I)

!        Linkseigenvektoren werden zurücktransformiert, indem man
!        die vorige Anweisung ersetzt durch: 'S = 1.0D0/SCAL(I)'

         DO  J=1,N
            EIVEC(I,J) = EIVEC(I,J)*S
         ENDDO
   20    CONTINUE
      DO 40 I=LOW-1,1,-1
         K=SCAL(I)
         DO  J=1,N
            CALL SWAP(EIVEC(I,J),EIVEC(K,J))
         ENDDO
   40    CONTINUE
      DO 60 I=HIGH+1,N
         K=SCAL(I)
         DO  J=1,N
            CALL SWAP(EIVEC(I,J),EIVEC(K,J))
         ENDDO
   60    CONTINUE
      BALBAK = 0
      END



      INTEGER FUNCTION ELMHES (N,LOW,HIGH,MAT,PERM)

!*****************************************************************
!     Gegeben ist eine unsymmetrische Matrix a(1:N,1:N), die auf *
!     Hessenbergform durch nichtorthogonale                      *
!     Elementartransformationen überführt wird.                  *
!*****************************************************************
!                                                                *
!     ZWECK DES UNTERPROGRAMMS:                                  *
!     =========================                                  *
!     Gegeben ist eine unsymmetrische Matrix a(1:N,1:N). Dann    *
!     reduziert diese Prozedur die Teilmatrix der Ordnung        *
!     HIGH-LOW+1, die beim Element a(LOW,LOW) beginnt und beim   *
!     Element a(HIGH,HIGH) endet, auf Hessenbergform h durch     *
!     nichtorthogonale Elementartransformationen. Die Teilmatrix *
!     wird mit h überschrieben, wobei die Einzelheiten der       *
!     Transformationen in dem übrigbleibenden Dreieck unterhalb  *
!     von h und in dem Feld PERM abgespeichert werden.           *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     N:        Ordnung der vollen Matrix a                      *
!     LOW,HIGH: Ausgabeparameter einer Prozedur, die             *
!               a aufbereitet (siehe [PARL69]). Falls a nicht    *
!               derart aufbereitet ist, setze LOW:=1, HIGH:=N.   *
!     MAT:      die (1:N,1:N)-Matrix a, normalerweise in         *
!               aufbereiteter Form (siehe oben)                  *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     MAT:      ein (1:N,1:N)-Feld, das zu einem Teil aus der    *
!               abgeleiteten oberen Hessenbergmatrix besteht;    *
!               die Größe N(I,R+1), die bei der Reduktion        *
!               eine Rolle spielt, wird im (I,R)-Element         *
!               abgespeichert.                                   *
!     PERM:     ein Integerfeld, das die bei der Reduktion       *
!               ausgeführten Zeilen- und Spaltenvertauschungen   *
!               beschreibt                                       *
!                                                                *
!     RÜCKGABEWERT:                                              *
!     =============                                              *
!     0:        kein Fehler                                      *
!                                                                *
!     LOKALE NAMEN:                                              *
!     =============                                              *
!     ZERO,ONE: Gleitkommakonstanten                             *
!     I,J,M:    Zählvariablen                                    *
!     X,Y:      Hilfsvariablen zur Aufnahme von Matrixelementen  *
!               und Zwischenergebnissen                          *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  Benötigte Unterprogramme: SWAP                                *
!                                                                *
!                                                                *
!  Quellen : Martin, R. S. und Wilkinson, J. H., siehe [MART68]. *
!                                                                *
!*****************************************************************
!                                                                *
!  Autor     : Jürgen Dietel                                     *
!  Datum     : 10.04.1987                                        *
!  Quellcode : FORTRAN 77                                        *
!                                                                *
!*****************************************************************

      implicit real*8 (a-h,o-z)

      integer PERM(MNB)
      integer N, LOW, HIGH
      real*8 MAT(MNB,MNB)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0)

      DO 70 M=LOW+1,HIGH-1
         I = M
         X = ZERO
         DO 10 J=M,HIGH
            IF (ABS(MAT(J,M-1)) .GT. ABS(X)) THEN
               X = MAT(J,M-1)
               I=J
            ENDIF
   10       CONTINUE
         PERM(M) = I
         IF (I .NE. M) THEN

!           Zeilen und Spalten von MAT vertauschen

            DO  J=M-1,N
               CALL SWAP(MAT(I,J),MAT(M,J))
            ENDDO
            DO  J=1,HIGH
               CALL SWAP(MAT(J,I),MAT(J,M))
            ENDDO
         ENDIF
         IF (X .NE. ZERO) THEN
            DO 60 I=M+1,HIGH
               Y = MAT(I,M-1)
               IF (Y .NE. ZERO) THEN
                  Y = Y/X
                  MAT(I,M-1) = Y
                  DO  J=M,N
                     MAT(I,J) = MAT(I,J)-Y*MAT(M,J)
                  ENDDO
                  DO  J=1,HIGH
                     MAT(J,M) = MAT(J,M)+Y*MAT(J,I)
                  ENDDO
               ENDIF
   60          CONTINUE
         ENDIF
   70    CONTINUE
      ELMHES = 0
      END



      INTEGER FUNCTION ELMTRA (N,LOW,HIGH,MAT,PERM,H)

!*****************************************************************
!                                                                *
!     Diejenige Matrix in dem Feld H(1:N,1:N) abspeichern, die   *
!     sich aus den Informationen ergibt, die die Prozedur ELMHES *
!     hinterlassen hat im unteren Dreieck der Hessenbergmatrix   *
!     h, und zwar im Feld MAT(1:N,1:N) und im Integerfeld        *
!     PERM(1:N)                                                  *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     N:        Ordnung der Hessenbergmatrix h                   *
!     LOW,HIGH: Integerzahlen, die von der Prozedur BALAN        *
!               erzeugt wurden (falls sie verwandt wurde;        *
!               andernfalls setze LOW:=1, HIGH:=N.)              *
!     PERM:     ein von ELMHES erzeugtes (1:N)-Integerfeld       *
!     MAT:      ein (1:N,1:N)-Feld, das von ELMHES erzeugt wurde *
!               und die Hessenbergmatrix h und die               *
!               Multiplikatoren enthält, die benutzt wurden,     *
!               um es aus der allgemeinen Matrix a zu erzeugen   *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     H:        dasjenige (1:N,1:N)-Feld, das die Ähnlichkeits-  *
!               transformation von A in H definiert              *
!                                                                *
!     RÜCKGABEWERT:                                              *
!     =============                                              *
!     0:        kein Fehler                                      *
!                                                                *
!     LOKALE NAMEN:                                              *
!     =============                                              *
!     ZERO,ONE: Gleitkommakonstanten                             *
!     I,J,K:    Indexvariablen                                   *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  Benötigte Unterprogramme: keine                               *
!                                                                *
!                                                                *
!  Quellen : Peters, G. und Wilkinson, J. H., siehe [PETE70].    *
!                                                                *
!*****************************************************************
!                                                                *
!  Autor     : Jürgen Dietel                                     *
!  Datum     : 10.04.1987                                        *
!  Quellcode : FORTRAN 77                                        *
!                                                                *
!*****************************************************************

      implicit real*8 (a-h,o-z)

      integer PERM(MNB)
      integer N, LOW, HIGH
      Dimension H(MNB,MNB)
      real*8 MAT(MNB,MNB)
      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0)

      DO 20 I=1,N
         DO  J=1,N
            H(I,J) = ZERO
         ENDDO
         H(I,I) = ONE
   20    CONTINUE

      DO 50 I=HIGH-1,LOW+1,-1
         J=PERM(I)
         DO  K=I+1,HIGH
            H(K,I) = MAT(K,I-1)
         ENDDO
         IF (I .NE. J) THEN
            DO 40 K=I,HIGH
               H(I,K)=H(J,K)
               H(J,K)=ZERO
   40          CONTINUE
            H(J,I) = ONE
         ENDIF
   50    CONTINUE
      ELMTRA = 0
      END



      INTEGER FUNCTION ORTHES (N,LOW,HIGH,MAT,D,EPSM)

!*****************************************************************
! Diese Funktion transformiert die Matrix MAT durch              *
! Householdertransformationen auf obere Hessenberggestalt.       *
! Alle wichtigen Informationen über die Transformation werden    *
! dabei in dem ungenutzten unteren Dreieck der Hessenbergmatrix  *
! und im Vektor D abgespeichert.                                 *
!                                                                *
! Eingabeparameter:                                              *
! =================                                              *
! N        Ordnung der Matrix MAT                                *
! LOW  \   Die Zeilen von 1 bis LOW-1 bzw. die Zeilen von HIGH+1 *
! HIGH  >  bis N enthalten die isolierten Eigenwerte, d. h.      *
!      /   Eigenwerte zu Vielfachen von Einheitsvektoren.        *
! MAT      die umzuformende (1:N,1:N)-Matrix                     *
! EPSM     Maschinengenauigkeit                                  *
!                                                                *
! Ausgabeparameter:                                              *
! =================                                              *
! MAT      die gesuchte Hessenbergmatrix zusammen mit dem ersten *
!          Teil der Umformungsinformationen unterhalb der        *
!          Subdiagonalen                                         *
! D        (1:N)-Vektor mit den restlichen                       *
!          Umformungsinformationen                               *
!                                                                *
! Funktionswert:                                                 *
! ==============                                                 *
! Fehlercode. Hier kann nur der Wert 0 (kein Fehler) auftreten.  *
!                                                                *
! lokale Namen:                                                  *
! =============                                                  *
! I,J,M  Laufvariablen                                           *
! S      Euklidische Norm sigma des subdiagonalen Spaltenvektors *
!        v von MAT, der auf ein Vielfaches des Einheitsvektors   *
!        e1 = (1,0,...,0) gespiegelt werden soll                 *
!        (v = (v1,...,v(HIGH-M+1))                               *
! X      zunächst erstes Element von v, dann Summationsvariable  *
!        bei der eigentlichen Householder-Transformation         *
! Y      zunächst sigma^2, dann ||u||^2, wobei                   *
!        u := v +- sigma * e1                                    *
! EPS    Toleranzschwelle für die Pruefung, ob die               *
!        Transformation gültig ist                               *
! ZERO   doppeltgenaue Gleitkommakonstante 0                     *
!                                                                *
!*****************************************************************
! Quelle: Numerische Mathematik 12 (1968), Seiten 359 und 360    *
! Autor:  Jürgen Dietel, Rechenzentrum der RWTH Aachen           *
! Datum:  DO 15. 7. 1993                                         *
!*****************************************************************

      implicit real*8 (a-h,o-z)

      integer N, LOW, HIGH
      real*8 MAT(MNB,MNB)
      Dimension D(MNB)
      PARAMETER        (ZERO = 0.0D0)

      EPS = 128 * EPSM

      DO 80 M = LOW+1,HIGH-1

         Y = ZERO
         DO 10 I=HIGH,M,-1
            X    = MAT(I,M-1)
            D(I) = X
            Y    = Y + X * X
   10       CONTINUE
         IF (Y .LE. EPS) THEN
            S = ZERO
         ELSE

            IF (X .GE. ZERO) THEN
               S = -SQRT(Y)
            ELSE
               S = SQRT(Y)
            ENDIF
            Y    = Y - X * S
            D(M) = X - S

!           MAT von links mit  (E-(u * uT)/y)  multiplizieren

            DO 40 J=M,N
               X = ZERO
               DO 20 I=HIGH,M,-1
                  X = X + D(I) * MAT(I,J)
   20             CONTINUE
               X = X / Y
               DO 30 I=M,HIGH
                  MAT(I,J) = MAT(I,J) - X * D(I)
   30             CONTINUE
   40          CONTINUE

!           MAT von rechts mit  (E-(u * uT)/y)  multiplizieren

            DO 70 I=1,HIGH
               X = ZERO
               DO 50 J=HIGH,M,-1
                  X = X + D(J) * MAT(I,J)
   50             CONTINUE
               X = X / Y
               DO 60 J=M,HIGH
                  MAT(I,J) = MAT(I,J) - X * D(J)
   60             CONTINUE
   70          CONTINUE

         ENDIF

        MAT(M,M-1) = S
   80   CONTINUE

      ORTHES = 0

      END



      INTEGER FUNCTION ORTTRA (N,LOW,HIGH,MAT,D,V)

!*****************************************************************
! die Matrix V der Transformation berechnen, deren Informationen *
! sich bei der Householder-Reduktion der Matrix MAT auf obere    *
! Hessenbergform unterhalb der Hessenbergmatrix in MAT und im    *
! Vektor D angesammelt haben. Der Inhalt von D wird dabei        *
! zerstört.                                                      *
!                                                                *
! Eingabeparameter:                                              *
! =================                                              *
! N        Ordnung der Matrix MAT                                *
! LOW  \   Die Zeilen von 1 bis LOW-1 bzw. die Zeilen von HIGH+1 *
! HIGH  >  bis N enthalten die isolierten Eigenwerte, d. h.      *
!      /   Eigenwerte zu Vielfachen von Einheitsvektoren.        *
! MAT      (1:N,1:N)-Matrix, die durch ORTHES auf obere          *
!          Hessenbergform gebracht und dabei auch mit            *
!          Umformungsinformationen der orthogonalen Reduktion    *
!          versehen wurde                                        *
! D        (1:N)-Vektor mit den restlichen                       *
!          Umformungsinformationen der Transformation auf        *
!          Hessenbergform                                        *
!                                                                *
! Ausgabeparameter:                                              *
! =================                                              *
! D        durch das Unterprogramm zerstörter Eingabevektor      *
! V        (1:N,1:N)-Matrix, die die Ähnlichkeitstransformation  *
!          auf die obere Hessenbergmatrix in MAT definiert       *
!                                                                *
! Funktionswert:                                                 *
! ==============                                                 *
! Fehlercode. Hier kann nur der Wert 0 (kein Fehler) auftreten.  *
!                                                                *
! lokale Namen:                                                  *
! =============                                                  *
! I,J,M  Laufvariablen                                           *
! X      Summationsvariable bei der Householder-Transformation   *
! Y      sigma  bzw.  sigma * (v1 +- sigma)                      *
! ZERO   doppeltgenaue Gleitkommakonstante 0                     *
! ONE    doppeltgenaue Gleitkommakonstante 1                     *
!                                                                *
!*****************************************************************
! Quelle: Numerische Mathematik 16 (1970), Seite 191             *
! Autor:  Jürgen Dietel, Rechenzentrum der RWTH Aachen           *
! Datum:  DO 15. 7. 1993                                         *
!*****************************************************************

      implicit real*8 (a-h,o-z)

      integer N, LOW, HIGH
      logical VEC,ORTHO,EVNORM
      Dimension V(MNB,MNB),D(MNB)
      real*8 MAT(MNB,MNB)

      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0)

!     die Einheitsmatrix in V erzeugen

      DO 20 I=1,N
         DO 10 J=1,N
           V(I,J) = ZERO
   10      CONTINUE
         V(I,I) = ONE
   20    CONTINUE

!     die Transformationen, durch die MAT auf Hessenbergform
!     gebracht wurde, nun auch auf die Einheitsmatrix anwenden.
!     Auf diese Weise entsteht aus V die gesuchte
!     Transformationsmatrix.

      DO 70 M=HIGH-1,LOW+1,-1
         Y = MAT(M,M-1)

         IF (Y .NE. ZERO) THEN
            Y = Y * D(M)
            DO 30 I=M+1,HIGH
               D(I) = MAT(I,M-1)
   30          CONTINUE
            DO 60 J=M,HIGH
            X = ZERO
            DO 40 I=M,HIGH
               X = X + D(I) * V(I,J)
   40          CONTINUE
            X = X / Y
            DO 50 I=M,HIGH
               V(I,J) = V(I,J) + X * D(I)
   50          CONTINUE
   60       CONTINUE
         ENDIF
   70    CONTINUE

      ORTTRA = 0

      END



      INTEGER FUNCTION HQR2 (VEC,N,LOW,HIGH,H,VALR,VALI,
     &                       EIVEC,CNT,EPS)

!*****************************************************************
!     Findet die Eigenwerte und Eigenvektoren einer reellen      *
!     Matrix, die, auf obere Hessenbergform reduziert ist        *
!*****************************************************************
!                                                                *
!     Findet die Eigenwerte und (bei gesetzter VEC-Flagge)       *
!     Eigenvektoren einer reellen Matrix, die, auf obere         *
!     Hessenbergform reduziert, im Feld H(1:N,1:N) steht, wobei  *
!     das Produkt der bisher durchgeführten Transformationen im  *
!     Feld EIVEC(1:N,1:N) steht.                                 *
!     Die Real- und die Imaginärteile der Eigenwerte werden in   *
!     den Feldern VALR(1:N),VALI(1:N) und die Eigenvektoren im   *
!     Feld EIVEC(1:N,1:N) gebildet, wo nur ein komplexer Vektor, *
!     nämlich der zu der Wurzel mit positivem Imaginärteil       *
!     gehörige, für jedes komplexe Paar von Eigenwerten erzeugt  *
!     wird. LOW und HIGH sind zwei Integerzahlen, die beim Aus-  *
!     balancieren entstehen, wo Eigenwerte in den Positionen 1   *
!     bis LOW-1 und HIGH+1 bis N isoliert werden. Falls keine    *
!     Ausbalancierung durchgeführt wurde, setze LOW:=1, HIGH:=N. *
!     Das Unterprogramm bricht mit einer Fehlermeldung ab, falls *
!     irgendein Eigenwert mehr als MAXSTP Iterationsschritte     *
!     benötigt.                                                  *
!                                                                *
!     Eingabeparameter:                                          *
!     =================                                          *
!     VEC:      Parameter fuer Eigenvektoren                     *
!                  = .TRUE.:  Bestimmung der Eigenvektoren       *
!                  = .FALSE.: Nur Berechnung der Eigenwerte      *
!     N:        Ordnung der Hessenbergmatrix h                   *
!     LOW,HIGH: von BALAN erzeugte Integerzahlen, falls          *
!               BALAN benutzt wurde. Ansonsten setze             *
!               LOW:=1, HIGH:=N.                                 *
!     EPS:      die kleinste Zahl auf dem Computer, für          *
!               die gilt: 1 + EPS > 1.                           *
!     H:        ein (1:N,1:N)-Feld, das die Matrix h in ihren    *
!               relevanten Teilen enthält                        *
!     EIVEC:    ein (1:N,1:N)-Feld, das die Matrix enthält, die  *
!               die Ähnlichkeitstransformation von a in h        *
!               definiert. (Es wird von ELMTRA erzeugt.)         *
!               Falls h die ursprüngliche Matrix ist,            *
!               setze EIVEC := Einheitsmatrix.                   *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     H:           der obere Teil dieses (1:N,1:N)-Feldes        *
!                  enthält die Eigenvektoren der                 *
!                  Fast-Dreiecksmatrix, die vom QR-Verfahren     *
!                  erzeugt wird.                                 *
!     VALR,VALI:   zwei (1:N)-Felder, die Real- und Imaginärteil *
!                  der Eigenwerte aufnehmen                      *
!     CNT:         ein (N,1)-Integerfeld, das die Zahl der       *
!                  Iterationsschritte für jeden Eigenwert        *
!                  aufnimmt. Falls zwei Eigenwerte als Paar      *
!                  gleichzeitig gefunden werden, so wird         *
!                  die Zahl der Iterationsschritte mit einem     *
!                  positiven Vorzeichen für den ersten und       *
!                  einem negativen Vorzeichen für den            *
!                  zweiten Eigenwert eingetragen.                *
!     EIVEC:       ein (1:N,1:N)-Feld, das bei gesetzter         *
!                  VEC-Flagge die nichtnormalisierten            *
!                  Eigenvektoren der ursprünglichen vollen       *
!                  Matrix aufnimmt (falls h nicht die            *
!                  Ausgangsmatrix war). Falls der i. Eigenwert   *
!                  reell ist, so ist die i. Spalte von EIVE! der *
!                  dazugehörige reelle Eigenvektor. Falls die    *
!                  Eigenwerte i und i+1 ein komplexes Paar       *
!                  bilden, geben i. und (i+1). Spalte Real- und  *
!                  Imaginärteil desjenigen Eigenvektors an, der  *
!                  zu dem Eigenwert mit positivem Imaginärteil   *
!                  gehört.                                       *
!                                                                *
!     RÜCKGABEWERT:                                              *
!     =============                                              *
!     0:           kein Fehler                                   *
!     1:           die Parameter N, LOW oder HIGH haben          *
!                  unerlaubte Werte.                             *
!     2:           alle Eigenvektoren sind der Nullvektor.       *
!     3:           die maximale Schrittzahl ist überschritten.   *
!                                                                *
!     LOKALE NAMEN:                                              *
!     =============                                              *
!     ZERO,ONE,TWO,PT75,PT4375: wichtige Gleitkommakonstanten    *
!     MAXSTP:                   Konstante für die maximale       *
!                               Schrittzahl                      *
!     I,J,K,L,M,N,NA,EN:        Indexvariablen                   *
!     ITER:                     Schrittzähler                    *
!     P,Q,R,S,T,W,X,Y,Z,NORM,                                    *
!       RA,SA,VR,VI:            Hilfsvariablen für Gleitkomma-   *
!                               berechnungen                     *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  Benötigte Unterprogramme: COMDIV                              *
!                                                                *
!                                                                *
!  Quellen : Peters, G. und Wilkinson, J. H., siehe [PETE70].    *
!                                                                *
!*****************************************************************
!                                                                *
!  Autor     : Jürgen Dietel                                     *
!  Datum     : 14.7.1993                                         *
!  Quellcode : FORTRAN 77                                        *
!                                                                *
!*****************************************************************

      implicit real*8 (a-h,o-z)

      integer CNT(MNB)
      integer N, LOW, HIGH, MAXSTEP, EN
      logical VEC
      dimension EIVEC(MNB,MNB),H(MNB,MNB),VALR(MNB),VALI(MNB)
      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0,TWO = 2.0D0,
     &                  PT75 = 0.75D0,PT4375 = 0.4375D0,
     &                  MAXSTP = 100)
      real*8 NORM

!     Fehler 1: die Parameter N, LOW oder HIGH haben unerlaubte
!     Werte:

      IF (N .LT. 1 .OR. LOW .LT. 1 .OR. HIGH .GT. N) THEN
         HQR2 = 1
         RETURN
      ENDIF

!     Vorbesetzung für die bei der Ausbalancierung gefundenen
!     isolierten Eigenwerte:

      DO 10 I=1,N
         IF (I .LT. LOW .OR. I .GT. HIGH) THEN
            VALR(I) = H(I,I)
            VALI(I) = ZERO
            CNT(I) = 0
         ELSE
            CNT(I) = -1
         ENDIF
   10    CONTINUE

      EN = HIGH
      T = ZERO
   15    IF (EN .LT. LOW) GOTO 333
         ITER = 0
         NA = EN-1

!           nach einem einzelnen kleinen Subdiagonalelement
!           suchen:

   20       DO  L=EN,LOW+1,-1
               IF (ABS(H(L,L-1)) .LE. EPS*
     &             (ABS(H(L-1,L-1))+ABS(H(L,L)))) GOTO 40
            ENDDO
   40       X = H(EN,EN)
            IF (L .EQ. EN) THEN

!              eine Wurzel gefunden:

               VALR(EN) = X + T
               H(EN,EN) = VALR(EN)
               VALI(EN) = ZERO
               CNT(EN) = ITER
               EN = NA
               GOTO 15
            ENDIF

            Y = H(NA,NA)
            W = H(EN,NA)*H(NA,EN)
            IF (L .EQ. NA) THEN

!              zwei Wurzeln gefunden:

               P = (Y-X)/TWO
               Q = P*P+W
               Z = SQRT(ABS(Q))
               H(EN,EN) = X+T
               X = H(EN,EN)
               H(NA,NA) = Y+T
               CNT(EN) = -ITER
               CNT(NA) = ITER
               IF (Q .GE. ZERO) THEN

!                 ein reelles Paar gefunden:

                  IF (P .LT. ZERO) Z = -Z
                  Z = P+Z
                  VALR(NA) = X+Z
                  VALR(EN) = X-W/Z
                  VALI(NA) = ZERO
                  VALI(EN) = ZERO
                  X = H(EN,NA)
                  R = SQRT(X*X+Z*Z)
                  IF (VEC) THEN
                     P = X/R
                     Q = Z/R

!                    Zeilenmodifikation:

                     DO 50 J=NA,N
                        Z = H(NA,J)
                        H(NA,J) = Q*Z+P*H(EN,J)
                        H(EN,J) = Q*H(EN,J)-P*Z
   50                   CONTINUE

!                    Spaltenmodifikation:

                     DO 60 I=1,EN
                        Z = H(I,NA)
                        H(I,NA) = Q*Z+P*H(I,EN)
                        H(I,EN) = Q*H(I,EN)-P*Z
   60                   CONTINUE

!                    Akkumulation:

                     DO 70 I=LOW,HIGH
                        Z = EIVEC(I,NA)
                        EIVEC(I,NA) = Q*Z+P*EIVEC(I,EN)
                        EIVEC(I,EN) = Q*EIVEC(I,EN)-P*Z
   70                   CONTINUE
                  ENDIF
               ELSE

!                 komplexes Paar:

                  VALR(NA) = X+P
                  VALR(EN) = VALR(NA)
                  VALI(NA) = Z
                  VALI(EN) = -Z
               ENDIF
               EN = EN-2
               GOTO 15
            ENDIF

            IF (ITER .EQ. MAXSTP) THEN

!              Fehler 3: maximale Schrittzahl überschritten:

               CNT(EN) = MAXSTP+1
               HQR2 = 3
               RETURN
            ENDIF
            IF (MOD(ITER,10) .EQ. 0 .AND. ITER .NE. 0) THEN

!              einen ungewöhnlichen Shift durchführen:

               T = T+X
               DO  I=LOW,EN
                  H(I,I) = H(I,I)-X
               ENDDO
               S = ABS(H(EN,NA))+ABS(H(NA,EN-2))
               X = PT75*S
               Y = X
               W = -PT4375*S*S
            ENDIF
            ITER = ITER+1

!           nach zwei aufeinanderfolgenden kleinen
!           Subdiagonalelementen suchen:

            DO 90 M=EN-2,L,-1
               Z = H(M,M)
               R = X-Z
               S = Y-Z
               P = (R*S-W)/H(M+1,M)+H(M,M+1)
               Q = H(M+1,M+1)-Z-R-S
               R = H(M+2,M+1)
               S = ABS(P)+ABS(Q)+ABS(R)
               P = P/S
               Q = Q/S
               R = R/S
               IF (M .EQ. L) GOTO 100
               IF (ABS(H(M,M-1))*(ABS(Q)+ABS(R)) .LE. EPS*ABS(P)*
     &             (ABS(H(M-1,M-1))+ABS(Z)+ABS(H(M+1,M+1))))
     &             GOTO 100
   90          CONTINUE
  100       DO  I=M+2,EN
               H(I,I-2) = ZERO
            ENDDO
            DO  I=M+3,EN
               H(I,I-3) = ZERO
            ENDDO

!           ein doppelter QR-Schritt, der die Zeilen L bis EN und
!           die Spalten M bis EN des ganzen Feldes betrifft:

            DO 200 K=M,NA
               IF (K .NE. M) THEN
                  P = H(K,K-1)
                  Q = H(K+1,K-1)
                  IF (K .NE. NA) THEN
                     R = H(K+2,K-1)
                  ELSE
                     R = ZERO
                  ENDIF
                  X = ABS(P)+ABS(Q)+ABS(R)
                  IF (X .EQ. ZERO) GOTO 200
                  P = P/X
                  Q = Q/X
                  R = R/X
               ENDIF
               S = SQRT(P*P+Q*Q+R*R)
               IF (P .LT. ZERO) S = -S
               IF (K .NE. M) THEN
                  H(K,K-1) = -S*X
               ELSEIF (L .NE. M) THEN
                  H(K,K-1) = -H(K,K-1)
               ENDIF
               P = P+S
               X = P/S
               Y = Q/S
               Z = R/S
               Q = Q/P
               R = R/P

!              Zeilenmodifikation:

               DO 130 J=K,N
                  P = H(K,J)+Q*H(K+1,J)
                  IF (K .NE. NA) THEN
                     P = P+R*H(K+2,J)
                     H(K+2,J) = H(K+2,J)-P*Z
                  ENDIF
                  H(K+1,J) = H(K+1,J)-P*Y
                  H(K,J) = H(K,J)-P*X
  130             CONTINUE
               J = MIN(K+3,EN)

!              Spaltenmodifikation:

               DO 140 I=1,J
                  P = X*H(I,K)+Y*H(I,K+1)
                  IF (K .NE. NA) THEN
                     P = P+Z*H(I,K+2)
                     H(I,K+2) = H(I,K+2)-P*R
                  ENDIF
                  H(I,K+1) = H(I,K+1)-P*Q
                  H(I,K) = H(I,K)-P
  140             CONTINUE

               IF (VEC) THEN

!                 Transformationen akkumulieren:

                  DO 150 I=LOW,HIGH
                     P = X*EIVEC(I,K)+Y*EIVEC(I,K+1)
                     IF (K .NE. NA) THEN
                        P = P+Z*EIVEC(I,K+2)
                        EIVEC(I,K+2) = EIVEC(I,K+2)-P*R
                     ENDIF
                     EIVEC(I,K+1) = EIVEC(I,K+1)-P*Q
                     EIVEC(I,K) = EIVEC(I,K)-P
  150                CONTINUE
               ENDIF
  200          CONTINUE
            GOTO 20


  333 IF (.NOT. VEC) THEN
         HQR2 = 0
         RETURN
      ENDIF


!     alle Wurzeln gefunden, nun wird rücktransformiert:

!     1-Norm von H bestimmen:

      NORM = ZERO
      K=1
      DO  I=1,N
         DO  J=K,N
            NORM = NORM+ABS(H(I,J))
         ENDDO
         K = I
      ENDDO
      IF (NORM .EQ. ZERO) THEN
!        Fehler 2: 1-Norm von H ist gleich 0:
         HQR2 = 2
         RETURN
      ENDIF

!     Rücktransformation:

      DO 207 EN=N,1,-1
         P = VALR(EN)
         Q = VALI(EN)
         NA = EN - 1
         IF (Q .EQ. ZERO) THEN

!           reeller Vektor:

            M = EN
            H(EN,EN) = ONE
            DO 63 I=NA,1,-1
               W = H(I,I)-P
               R = H(I,EN)
               DO  J=M,NA
                  R = R+H(I,J)*H(J,EN)
               ENDDO
               IF (VALI(I) .LT. ZERO) THEN
                  Z = W
                  S = R
               ELSE
                  M = I
                  IF (VALI(I) .EQ. ZERO) THEN
                     IF (W .NE. ZERO) THEN
                        H(I,EN) = -R/W
                     ELSE
                        H(I,EN) = -R/(EPS*NORM)
                     ENDIF
                  ELSE

!                    löse das Gleichungssystem:
!                    [ W   X ] [ H(I,EN)   ]   [ -R ]
!                    [       ] [           ] = [    ]
!                    [ Y   Z ] [ H(I+1,EN) ]   [ -S ]

                     X = H(I,I+1)
                     Y = H(I+1,I)
                     Q = (VALR(I)-P)*(VALR(I) - P)+
     &                    VALI(I)*VALI(I)
                     T = (X*S-Z*R)/Q
                     H(I,EN) = T
                     IF (ABS(X) .GT. ABS(Z)) THEN
                        H(I+1,EN) = (-R-W*T)/X
                     ELSE
                        H(I+1,EN) = (-S-Y*T)/Z
                     ENDIF
                  ENDIF
               ENDIF
   63          CONTINUE
         ELSEIF (Q .LT. ZERO) THEN

!           komplexer Vektor, der zu LAMBDA = P - I * Q gehört:

            M = NA
            IF (ABS(H(EN,NA)) .GT. ABS(H(NA,EN))) THEN
               H(NA,NA) = -(H(EN,EN)-P)/H(EN,NA)
               H(NA,EN) = -Q/H(EN,NA)
            ELSE
               CALL COMDIV(-H(NA,EN),ZERO,H(NA,NA)-P,Q,
     &                     H(NA,NA),H(NA,EN))
            ENDIF
            H(EN,NA) = ONE
            H(EN,EN) = ZERO
            DO 190 I=NA-1,1,-1
               W = H(I,I)-P
               RA = H(I,EN)
               SA = ZERO
               DO 75 J=M,NA
                  RA = RA+H(I,J)*H(J,NA)
                  SA = SA+H(I,J)*H(J,EN)
   75             CONTINUE
               IF (VALI(I) .LT. ZERO) THEN
                  Z = W
                  R = RA
                  S = SA
               ELSE
                  M = I
                  IF (VALI(I) .EQ. ZERO) THEN
                      CALL COMDIV(-RA,-SA,W,Q,H(I,NA),H(I,EN))
                  ELSE

!                    löse die komplexen Gleichungen:
!           [ W+Q*I   X   ] [H(I,NA)+H(I,EN)*I    ]   [-RA-SA*I]
!           [             ] [                     ] = [        ]
!           [   Y   Z+Q*I ] [H(I+1,NA)+H(I+1,EN)*I]   [-R-S*I  ]

                     X = H(I,I+1)
                     Y = H(I+1,I)
                     VR = (VALR(I)-P)*(VALR(I)-P)+
     &                    VALI(I)*VALI(I)-Q*Q
                     VI = TWO*Q*(VALR(I)-P)
                     IF (VR .EQ. ZERO .AND. VI .EQ. ZERO) VR =
     &                   EPS*NORM*
     &                   (ABS(W)+ABS(Q)+ABS(X)+ABS(Y)+ABS(Z))
                     CALL COMDIV(X*R-Z*RA+Q*SA,X*S-Z*SA-Q*RA,
     &                           VR,VI,H(I,NA),H(I,EN))
                     IF (ABS(X) .GT. ABS(Z)+ABS(Q)) THEN
                        H(I+1,NA) = (-RA-W*H(I,NA)+Q*H(I,EN))/X
                        H(I+1,EN) = (-SA-W*H(I,EN)-Q*H(I,NA))/X
                     ELSE
                       CALL COMDIV(-R-Y*H(I,NA),-S-Y*H(I,EN),Z,Q,
     &                             H(I+1,NA),H(I+1,EN))
                     ENDIF
                  ENDIF
               ENDIF
  190          CONTINUE
         ENDIF
  207    CONTINUE

!     zu den isolierten Wurzeln gehörige Vektoren:

      DO 230 I=1,N
         IF (I .LT. LOW .OR. I .GT. HIGH) THEN
            DO  J=I+1,N
               EIVEC(I,J) = H(I,J)
            ENDDO
         ENDIF
  230    CONTINUE

!     mit der Transformationsmatrix multiplizieren, um die
!     Vektoren der ursprünglichen vollen Matrix zu erhalten:

      DO 300 J=N,LOW,-1
         IF (J .LE. HIGH) THEN
            M = J
         ELSE
            M = HIGH
         ENDIF
         L = J-1
         IF (VALI(J) .LT. ZERO) THEN
            DO 330 I=LOW,HIGH
               Y = ZERO
               Z = ZERO
               DO 320 K=LOW,M
                  Y = Y+EIVEC(I,K)*H(K,L)
                  Z = Z+EIVEC(I,K)*H(K,J)
  320             CONTINUE
               EIVEC(I,L) = Y
               EIVEC(I,J) = Z
  330          CONTINUE
         ELSE
            IF (VALI(J) .EQ. ZERO) THEN
               DO  I=LOW,HIGH
                  Z = ZERO
                  DO  K=LOW,M
                     Z =Z+EIVEC(I,K)*H(K,J)
                  ENDDO
                  EIVEC(I,J) = Z
               ENDDO
            ENDIF
         ENDIF
  300    CONTINUE

!     Rückgabe von 0: kein Fehler:

      HQR2 = 0
      END



      SUBROUTINE COMDIV (AR,AI,BR,BI,RESR,RESI)

!*****************************************************************
!                                                                *
!     komplexe Division: RESR+I*RESI := (AR+I*AI)/(BR+I*BI).     *
!     Diese Prozedur sollte nicht mit BR=BI=0 aufgerufen werden. *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     AR,AI:     Real- und Imaginärteil des Dividenden           *
!     BR,BI:     Real- und Imaginärteil des Divisors             *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     RESR,RESI: Real- und Imaginärteil des Quotienten           *
!                                                                *
!     LOKALE NAMEN:                                              *
!     =============                                              *
!     ZERO:              Gleitkommakonstante 0                   *
!     TEMP1,TEMP2,TEMP3: Hilfsvariablen zur Speicherung von      *
!                        Zwischenergebnissen                     *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  Benötigte Unterprogramme: keine                               *
!                                                                *
!                                                                *
!  Quellen : Martin, R. S. und Wilkinson, J. H., siehe [MART68]. *
!                                                                *
!*****************************************************************
!                                                                *
!  Autor     : Jürgen Dietel                                     *
!  Datum     : 10.04.1987                                        *
!  Quellcode : FORTRAN 77                                        *
!                                                                *
!*****************************************************************

      implicit real*8 (a-h,o-z)

      PARAMETER        (ZERO = 0.0D0)

      IF (BR .EQ. ZERO .AND. BI .EQ. ZERO) THEN
         RESR = ZERO
         RESI = ZERO
         RETURN
      ENDIF
      IF (ABS(BR) .GT. ABS(BI)) THEN
         TEMP1 = BI/BR
         TEMP2 = TEMP1*BI+BR
         TEMP3 = (AR+TEMP1*AI)/TEMP2
         RESI = (AI-TEMP1*AR)/TEMP2
         RESR = TEMP3
      ELSE
         TEMP1 = BR/BI
         TEMP2 = TEMP1*BR+BI
         TEMP3 = (TEMP1*AR+AI)/TEMP2
         RESI = (TEMP1*AI-AR)/TEMP2
         RESR = TEMP3
      ENDIF
      END



      DOUBLE PRECISION FUNCTION COMABS (AR,AI)

!*****************************************************************
!                                                                *
!     Berechnung des Betrages der komplexen Zahl AR+I*AI:        *
!     COMABS := SQRT(AR*AR+AI*AI)                                *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     AR,AI: Real- und Imaginärteil der komplexen Zahl, deren    *
!            Betrag zu berechnen ist                             *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     keine                                                      *
!                                                                *
!     RÜCKGABEWERT:                                              *
!     =============                                              *
!     Betrag des komplexen Parameters                            *
!                                                                *
!     LOKALE NAMEN:                                              *
!     =============                                              *
!     ZERO,ONE:    Konstanten                                    *
!     TEMP1,TEMP2: Hilfsvariablen zur Speicherung von            *
!                  Zwischenergebnissen                           *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  Benötigte Unterprogramme: SWAP                                *
!                                                                *
!                                                                *
!  Quellen : Martin, R. S. und Wilkinson, J. H., siehe [MART68]. *
!                                                                *
!*****************************************************************
!                                                                *
!  Autor     : Jürgen Dietel                                     *
!  Datum     : 10.04.1987                                        *
!  Quellcode : FORTRAN 77                                        *
!                                                                *
!*****************************************************************

      implicit real*8 (a-h,o-z)
      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0)

      TEMP1 = ABS(AR)
      TEMP2 = ABS(AI)
      IF (AR .EQ. ZERO .AND. AI .EQ. ZERO) THEN
         COMABS = ZERO
         RETURN
      ENDIF
      IF (TEMP2 .GT. TEMP1) CALL SWAP(TEMP1,TEMP2)
      IF (TEMP2 .EQ. ZERO) THEN
         COMABS = TEMP1
      ELSE
         COMABS = TEMP1*SQRT(ONE+(TEMP2/TEMP1)**2)
      ENDIF
      END



      INTEGER FUNCTION NORMAL (N,V,WI)

!*****************************************************************
!                                                                *
!     NORMAL normalisiert die Eigenvektoren in der Maximumnorm.  *
!                                                                *
!     EINGABEPARAMETER:                                          *
!     =================                                          *
!     N:      die Ordnung der Matrix V                           *
!     V:      ein (1:N,1:N)-Feld vom Typ DOUBLE PRECISION,       *
!             dessen Spalten die Eigenvektoren bilden            *
!             (siehe EIGEN)                                      *
!     WI:     ein Feld mit N Komponenten vom Typ                 *
!             DOUBLE PRECISION, dessen Komponenten die           *
!             Imaginärteile der Eigenwerte sind                  *
!                                                                *
!     AUSGABEPARAMETER:                                          *
!     =================                                          *
!     V:      Matrix der normalisierten Eigenvektoren            *
!                                                                *
!     LOKALE NAMEN:                                              *
!     =============                                              *
!     ZERO,ONE: Gleitkommakonstanten 0 und 1                     *
!     I,J:      Indexvariablen                                   *
!     MAXI:     Hilfsvariable zur Berechnung der reellen         *
!               Vektornorm                                       *
!     TR,TI:    Hilfsvariablen zur Berechnung der komplexen      *
!               Vektornorm                                       *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  Benötigte Unterprogramme: COMABS, COMDIV                      *
!                                                                *
!*****************************************************************
!                                                                *
!  Autor     : Jürgen Dietel                                     *
!  Datum     : 10.04.1987                                        *
!  Quellcode : FORTRAN 77                                        *
!                                                                *
!*****************************************************************


      implicit real*8 (a-h,o-z)

      Dimension V(MNB,MNB),WI(MNB)
      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0)
      real*8 MAXI

      J = 1
   10 IF (J .GT. N) GOTO 80
         IF (WI(J) .EQ. ZERO) THEN
             MAXI = V(1,J)
             DO  I=2,N
                IF (ABS(V(I,J)) .GT. ABS(MAXI)) MAXI = V(I,J)
             ENDDO
             IF (MAXI .NE. ZERO) THEN
                MAXI = ONE/MAXI
                DO  I=1,N
                   V(I,J) = V(I,J)*MAXI
                ENDDO
             ENDIF
             J = J+1
         ELSE
            TR = V(1,J)
            TI = V(1,J+1)
            DO 30 I=2,N
               IF (COMABS(V(I,J),V(I,J+1)) .GT. COMABS(TR,TI))
     &          THEN
                  TR = V(I,J)
                  TI = V(I,J+1)
               ENDIF
   30          CONTINUE
            IF (TR .NE. ZERO .OR. TI .NE. ZERO) THEN
               DO  I=1,N
                  CALL COMDIV(V(I,J),V(I,J+1),TR,TI,
     &                        V(I,J),V(I,J+1))
               ENDDO
            ENDIF
            J = J+2
         ENDIF
         GOTO 10
   80 NORMAL = 0
      END


      SUBROUTINE SWAP (X,Y)

!*****************************************************************
!                                                                *
!     SWAP vertauscht die Werte der beiden DOUBLE PRECISION-     *
!     Variablen X und Y.                                         *
!                                                                *
!*****************************************************************
      implicit real*8 (a-h,o-z)

      TEMP = X
      X    = Y
      Y    = TEMP
      END

      SUBROUTINE indexx(n,arr,indx)
!-------------------------------------------------------------------------
!  Is required to sort the eigenvalues of the eigen() routine.
!  Indexes an array `arr' in a way such that arr(index(j)) is in ascending
!  order for j = 1, 2, ..., n. The input `arr' remains unchanged.
!  A `Numerical Recipes' routine.
!-------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      Dimension arr(MNB)
      PARAMETER (M=7,NSTACK=50)
      INTEGER indx(MNB),istack(NSTACK)

      do 11 j=1,n
         indx(j)=j
 11   continue
      jstack=0
      l=1
      ir=n
 1    if(ir-l.lt.M)then
         do 13 j=l+1,ir
            indxt=indx(j)
            a=arr(indxt)
            do 12 i=j-1,1,-1
               if(arr(indx(i)).le.a)goto 2
               indx(i+1)=indx(i)
 12         continue
            i=0
 2          indx(i+1)=indxt
 13      continue
         if(jstack.eq.0)return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         itemp=indx(k)
         indx(k)=indx(l+1)
         indx(l+1)=itemp
         if(arr(indx(l+1)).gt.arr(indx(ir)))then
            itemp=indx(l+1)
            indx(l+1)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l)).gt.arr(indx(ir)))then
            itemp=indx(l)
            indx(l)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l+1)).gt.arr(indx(l)))then
            itemp=indx(l+1)
            indx(l+1)=indx(l)
            indx(l)=itemp
         endif
         i=l+1
         j=ir
         indxt=indx(l)
         a=arr(indxt)
 3       continue
         i=i+1
         if(arr(indx(i)).lt.a)goto 3
 4       continue
         j=j-1
         if(arr(indx(j)).gt.a)goto 4
         if(j.lt.i)goto 5
         itemp=indx(i)
         indx(i)=indx(j)
         indx(j)=itemp
         goto 3
 5       indx(l)=indx(j)
         indx(j)=indxt
         jstack=jstack+2
c         if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1
      END


c***************************************************************
c SUBROUTINE eigen
c***************************************************************
c     eigenvalues and eigenvectors of a real symmetric matrix
c     a ----- secular determinant (input)
c             eigenvetors  column-wise (output)
c             a is destroyed
c     b ----- nsub x nsub array for intermidiate storage
c     nsub ----- order of a (input)
c                must be equal to or greater than 1
c     valu ----- eigenvalues in decreasing order (output)
c     msub -----  number of eigenvalues (input)
c                 if equal to -1, eigenvectors are not computed
c     anorm ----- matrix norm which is a sqrt of sum of a(i,j)**2
c                                        over all i & j (output)
c     nmax ----- row dimension of a & b
c ------------------------------------------------------------------------
      subroutine eigen1 (a,b,nsub,valu,msub,anorm,nnmax)

      parameter(nmax=1600)
      implicit real*8 (a-h,o-z)
      integer index
      dimension a(nnmax,nnmax), b(nnmax,nnmax), valu(nnmax)
      dimension diag(nmax), superd(nmax), wvec(nmax),pvec(nmax),
     1          qvec(nmax), vall(nmax), q(nmax), u(nmax),
     2          index(nmax), factor(nmax), v(nmax), t(nmax,3)
      common /f202/ diag, superd, wvec, pvec, t
      equivalence (wvec(1),vall(1),factor(1),u(1)),
     1            (pvec(1),qvec(1),q(1),v(1))
      data e1/2.3d-15/
c ------------------------------------------------------------------------
c     Householder similarity transformation to co-diagonal form
c ------------------------------------------------------------------------
      lpout = 0
      n=nsub
      m=msub
      do 10 i = 1,n
      valu(i) = 0.0d+00
   10 continue
      if (n.lt.0) go to 220
      if (m.le.0.or.n.le.2) go to 50
      do 40 i=2,n
      do 40 j=2,n
      b(j,i)=0.0d0
      if (i.eq.j) b(j,i)=1.0d0
   40 continue
   50 do 200 i=1,n
      i1=i+1
      i2=i1+1
      if (i2.gt.n) go to 160
      sum=0.0d0
      do 70 j=i2,n
   70 sum=sum+a(j,i)**2
      if (sum.eq.0.0d0) go to 160
      j=i1
      temp=a(j,i)
      sum=dsqrt(sum+ temp  **2)
      a(j,i)=-dsign(sum, temp  )
      wvec(j )=dsqrt(     1.0d0+dabs( temp  )/sum)
      div=dsign(    wvec(j )*sum, temp  )
      do 85 j=i2,n
   85 wvec(j)=a(j,i)/div
      scalar=0.0d0
      do 95 j=i1,n
      pvec(j)=0.0d0
      do 90 k=i1,n
   90 pvec(j)=pvec(j)+a(k,j)*wvec(k)
      scalar=scalar+pvec(j)*wvec(j)
   95 continue
      scalar=scalar/2.0d0
      do 120 j=i1,n
      qvec(j)=pvec(j)-scalar*wvec(j)
      do 120 k=i1,j
      a(k,j)=a(k,j)-(wvec(k)*qvec(j)+wvec(j)*qvec(k))
      a(j,k)=a(k,j)
  120 continue
      if (m.le.0) go to 160
      do 150 k=2,n
      temp=0.0d0
      do 140 j=i1,n
  140 temp=temp+wvec(j)*b(j,k)
      do 150 j=i1,n
      b(j,k)=b(j,k)-wvec(j)*temp
  150 continue
  160 j=i
      diag(i)=a(j,i)
      if (i.ne.n) superd(i)=a(j+1,i)
  200 continue
c ------------------------------------------------------------------------
c     Givens eigenvalue iteration from sturm chain of co-diagonal minors
c ------------------------------------------------------------------------
  220 n=iabs (n)
      m=iabs (m)
c ------------------------------------------------------------------------
c     Calculate norm of matrix and initialize eigenvalue bounds
c ------------------------------------------------------------------------
      anorm2=diag(1)**2
      if (n.eq.1) go to 235
      do 230 l=2,n
      q(l-1)=superd(l-1)**2
      anorm2=diag(l)**2+q(l-1)+q(l-1)+anorm2
  230 continue
  235 anorm=dsqrt(anorm2)
      if (m.eq.0) return
      do 240 l=1,m
      valu(l)=anorm
      vall(l)=-anorm
  240 continue
      eps1=anorm*e1
      if (eps1.eq.0.0d0) return
  250 do 570 l=1,m
  260 tau=(valu(l)+vall(l))/2.0d0
      if (    (tau-vall(l)).le.e1*(dabs(tau)+e1)) go to 570
      match=0
      t2=0.0d0
      t1=1.0d0
      do 450 l1=1,n
      p=diag(l1)-tau
      if (t2.ne.0.0d0) go to 330
  300 t1=dsign(1.0d0,t1)
  330 if (t1.ne.0.0d0) go to 400
      t0=-dsign(1.0d0,t2)
      t2=0.0d0
      if (q(l1-1).ne.0.0d0) go to 410
      go to 300
  400 t0=p-q(l1-1)*t2/t1
      t2=1.0d0
  410 if (t0) 440, 420, 430
  420 t2=t1
      if (t2) 440, 430, 430
  430 match=match+1
  440 t1=t0
  450 continue
      do 530 l1=l,m
      if (l1.le.match) go to 500
      if (valu(l1).le.tau) go to 260
      valu(l1)=tau
      go to 530
  500 vall(l1)=tau
  530 continue
      go to 260
  570 continue
c ------------------------------------------------------------------------
c     Eigenvectors of co-diagonal symmetric matrix--inverse iteration
c ------------------------------------------------------------------------
      m=msub
      if (m.lt.0) return
      do 970 i=1,m
      if (i.eq.1) go to 725
      if (valu(i-1)-valu(i).lt.1.0e+4*eps1) go to 730
  725 i1=-1
  730 i1=i1+1
      do 760 l=1,n
      v(l)=eps1*(i1*l+1)
      t(l,2)=diag(l)-valu(i)
      if (l.lt.n) go to 740
      t(l,3)=0.0d0
      go to 760
  740 t(l,3)=superd(l)
      if (t(l,3).eq.0.0d0) t(l,3)=eps1
      t(l+1,1)=t(l,3)
  760 continue
      do 820 j=1,n
      t(j,1)=t(j,2)
      t(j,2)=t(j,3)
      t(j,3)=0.0d0
      go to 780
  770 t(j,1)=eps1
c ------------------------------------------------------------------------
c   Replace function name  in next      statement  in single precision
c ------------------------------------------------------------------------
  780 vtemp=dabs(t(j,1))
      if (vtemp.lt.eps1) go to 770
      if (j.eq.n) go to 820
      index(j)=0
      if (dabs(t(j+1,1)).le.vtemp) go to 810
      index(j)=1
      do 800 k=1,3
      vtemp=t(j,k)
      t(j,k)=t(j+1,k)
      t(j+1,k)=vtemp
  800 continue
  810 vtemp   =t(j+1,1)/t(j,1)
      factor(j)=vtemp
      t(j+1,2)=t(j+1,2)-  vtemp  *t(j,2)
      t(j+1,3)=t(j+1,3)-  vtemp  *t(j,3)
  820 continue
      iter=1
      if (i1.gt.0) go to 920
  860 do 870 l1=1,n
      l=n+1-l1
      v( l )=(v( l )-t(l,2)*v( l+1 )-t(l,3)*v( l+2 ))/t(l,1)
  870 continue
      go to (875,920), iter
  875 iter=2
      if (n.eq.1) go to 860
  880 do 910 l=2,n
      if (index(l-1).eq.0) go to 900
      vtemp=v(l-1)
      v(l-1)=v(l)
      v(l)=vtemp
  900 v(l)=v(l)-factor(l-1)*v(l-1)
  910 continue
      go to 860
  920 if (i1.eq.0) go to 945
      do 940 l1=1,i1
      k=i-l1
      vtemp=0.0d0
      do 930 j=1,n
  930 vtemp=vtemp+a(j,k)*v(j)
      do 940 j=1,n
  940 v(j)=v(j)-a(j,k)*vtemp
      go to (880,945), iter
  945 vnorm2=0.0d0
      do 950 l=1,n
  950 vnorm2=vnorm2+v( l )**2
      vnorm=dsqrt(vnorm2)
      do 960 j=1,n
  960 a(j,i)=v( j )/vnorm
  970 continue
      n=nsub
      if (n.le.2) return
      do 990 i=1,m
      do 980 k=2,n
      u(k)=0.0d0
      do 980 j=2,n
  980 u(k)=u(k)+b(j,k)*a(j,i)
      do 990 j=2,n
  990 a(j,i)=u(j)

      return
      end

c ********************************************************************
c SUBROUTINE mult
c ********************************************************************
c    Matrix  multiplication
c------------------------------------------------------------
      subroutine mult(aa,bb,cc,im,m)

      implicit double precision (a-h, o-z)
      dimension aa(im,im),bb(im,im),cc(im,im)
      do i=1,m
          do j=1,m
          cc(i,j)=0.0d0
            do k=1,m
            cc(i,j)=cc(i,j)+aa(i,k)*bb(k,j)
            end do
          end do
       end do

       return
       end


c ********************************************************************
C SUBROUTINE matout
c ********************************************************************
      subroutine matout(a,im,in,m,n,label)

      implicit double precision (a-h,o-z)
      character*6 label
      dimension a(im,im)
      ihigh=0
  10  low=ihigh+1
      ihigh=ihigh+5
      ihigh=min0(ihigh,n)
              write(STDOUT,20)label,(i,i=low,ihigh)
  20  format(///,3x,5h THE ,a6,6h ARRAY,/15x,5(10x,i3,6x)//)
      do i=1,m
              write(STDOUT,40) i,(a(i,j),j=low,ihigh)
      enddo
  40  format(i8,5x,8(1x,d12.6))
      if(n-ihigh)50,50,10

  50  return
      end


c ********************************************************************
c SUBROUTINE scfiter
c ********************************************************************
c     This performs the scf iteration.
c------------------------------------------------------------------
      subroutine scfiter(cstp,ival)

      implicit real*8(a-h,o-z)
      logical set
      character*3 nh_i

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     :      /jvalue/orbj(MNS),kap(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /basc/nbasis,nocc,mbasis
     :  /skip/nskipe(MNS),nskipc(MNS)
     :  /dfm/emt(MNB,MNB),df_single(MNB,MNB),df_two(MNB,MNB)
     :  /conv/crit,maxit,npower
     :  /basis/gl(MN,MNBAS),gs(MN,MNBAS)
     : /propopt/istp,ialpvr
     :      /freez/minorb(MNS),maxorb(MNS)

C New common blocks
      Common/wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     :      /splflags/MRCC,IBASOPT
      Common/totenergies/en1,en2,en3,engMP2
c------------------------------------------------------------------
c Local dimension
c------------------------------------------------------------------
      dimension can(MNB,MNB),eigv(MNB,MNB),eng(MNB),stemp(MNB,MNB)
      dimension pold(MNB,MNB),p1(MNB,MNB),p2(MNB,MNB),p3(MNB,MNB)
      dimension drs(MNB,MNB),df_mat(MNB,MNB),a(MNB,MNB)
      dimension st(MNB,MNB),b(MNB,MNB),scr(MNB,MNB)
      dimension nh_i(MNS),sa(MNB),vec(MNB,MNB),oe(MNB,MNB)

      data nh_i/' S ',' P-',' P ',' D-',' D ',' F-',' F ',
     :    ' G-',' G ',' H-',' H ',' I-',' I ',' K-', ' K '/
C     Note : No J in this sequence
c--------------------------------------------------------------------
c       initialize arrays and constants required for the computation.
c--------------------------------------------------------------------

      IF(IBASOPT) OPEN(STDENG,file='ENERGIES.DAT',form='formatted',
     :                                            status='unknown')

      set = .true.
      iter=0
      en=zero
      eold=zero
      delta=zero

      do ia=1,MNB
      sa(ia)=zero
      do ib=1,MNB
      can(ib,ia)=zero
      drs(ib,ia)=zero
      stemp(ib,ia)=zero
      p1(ia,ib)=zero
      p2(ia,ib)=zero
      p3(ia,ib)=zero
      enddo
      enddo

      write(STDOUT,199)maxit,npower
      write(STDOUT,200)
      write(STDOUT,201)iter,delta
      write(*,201)iter,delta

 100  iter=iter+1

      do isym=1,nsym
      ndim=nbas(isym)
      ndim2=2*ndim

      do ia=1,ndim2
      ja=ia+2*nskipe(isym)
      do ib=1,ndim2
      jb=ib+2*nskipe(isym)
      stemp(ia,ib)=emt(ja,jb)
      enddo
      enddo

c      if(jpnt.ne.0)call matout(stemp,MNB,ndim2,ndim2,ndim2,4hS   )

      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=stemp(i,j)
      enddo
      enddo

      call eigen1(drs,scr,ndim2,sa,ndim2,bnorm,MNB)

C      Alternatively, Use the following subroutine
c      call eigen(drs,scr,ndim2,sa,ndim2,bnorm,MNB)

c-------------------------------------------------------------------
c       Canonical transformation
c------------------------------------------------------------------
      do ia=1,ndim2
      do ib=1,ndim2
      if(dabs(sa(ib)).gt.TINY)then
      can(ia,ib)=drs(ia,ib)/dsqrt(dabs(sa(ib)))
      st(ib,ia)=can(ia,ib)
      else
      write(*,*)'There is dioginalisation problem',sa(ib)
       stop
      endif
      enddo
      enddo

c      if(jpnt.ne.0)call matout(st,MNB,ndim2,ndim2,ndim2,4hXT  )

c      if(jpnt.ne.0)
c     :          call matout(df_single,MNB,mbasis,mbasis,mbasis,4hH   )

      if(iter.gt.1)call dfpotiter(isym,cstp,ival)

      do ic=1,ndim2
      jc=ic+2*nskipe(isym)
      do id=1,ndim2
      jd=id+2*nskipe(isym)
      if(iter.ne.1)then
      a(ic,id)=df_single(jc,jd)+df_two(jc,jd)
      else
      a(ic,id)=df_single(jc,jd)
      endif

      enddo
      enddo

c      if(jpnt.ne.0)call matout(a,MNB,ndim2,ndim2,ndim2,4hFOCK )

c      if(jpnt.ne.0)
c     :call matout(df_two,MNB,mbasis,mbasis,mbasis,4hVHF )

      call mult(a,can,b,MNB,ndim2)
      call mult(st,b,df_mat,MNB,ndim2)

c      if(jpnt.ne.0)
c     :      call matout(df_mat,MNB,ndim2,ndim2,ndim2,6hFRIME )

      call eigen1(df_mat,scr,ndim2,eng,ndim2,bnorm,MNB)
C      call eigen(df_mat,scr,ndim2,eng,ndim2,bnorm,MNB)

      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=df_mat(i,j)
      enddo
      enddo

c      if(jpnt.ne.0)
c     : call matout(drs,MNB,ndim2,ndim2,ndim2,4hC   )

      call mult(can,drs,eigv,MNB,ndim2)

      loc = 1
      last=ndim-1
      dowhile(loc.le.last)
      samin = eng(loc)
      locus = loc

      do  ia = loc,ndim
      if(eng(ia).lt.samin)then
      samin = eng(ia)
      locus = ia
      endif
      enddo

      if(locus.gt.loc)then
      swap = eng(loc)
      eng(loc) = eng(locus)
      eng(locus) = swap
      do  ib = 1,ndim2
      swap = eigv(ib,loc)
      eigv(ib,loc) = eigv(ib,locus)
      eigv(ib,locus) = swap
      enddo
      endif
      loc = loc + 1
      enddo

      do ia=1,ndim2
      do ib=1,ndim2
      b(ia,ib)=zero
      enddo
      b(ia,ia)=eng(ia)
      oe(ia,isym)=eng(ia)
      enddo

c      if(jpnt.ne.0)call matout(b,MNB,ndim2,ndim2,ndim2,4hE   )

c      if(jpnt.ne.0)call matout(eigv,MNB,ndim2,ndim2,ndim2,4hCP  )

      do ia=1,ndim2
      ja=ia+2*nskipe(isym)

      do ib=1,ndim2
      jb=ib+2*nskipe(isym)

      pold(ja,jb)=p1(ja,jb)
      p1(ja,jb)=zero
      p2(ja,jb)=zero
      p3(ja,jb)=zero

      nn=nocorb(isym)
      if(nn.ne.0)then
      do ibas=1,nocorb(isym)
      twojp1=two*orbj(isym)+one
      p1(ja,jb)=p1(ja,jb)+twojp1*eigv(ia,ibas)*eigv(ib,ibas)
      enddo
      endif

      if(ieng.eq.1.and.minorb(isym).ne.0)then

      do ibas=1,minorb(isym)
      twojp1=two*orbj(isym)+one
      p2(ja,jb)=p2(ja,jb)+twojp1*eigv(ia,ibas)*eigv(ib,ibas)
      enddo

      do ibas=minorb(isym)+1,nocorb(isym)
      twojp1=two*orbj(isym)+one
      p3(ja,jb)=p3(ja,jb)+twojp1*eigv(ia,ibas)*eigv(ib,ibas)
      enddo

      endif

      if(ivn.gt.0)then
      if(nvorb(isym).eq.1)then
      ibas=nocorb(isym)+1
C      twojp1=two/(two*orbj(isym)+one)    ! Bijaya Apr 2, 2010
C      p1(ja,jb)=p1(ja,jb)+twojp1*eigv(ia,ibas)*eigv(ib,ibas)
C      p3(ja,jb)=p3(ja,jb)+twojp1*eigv(ia,ibas)*eigv(ib,ibas)
      p1(ja,jb)=p1(ja,jb)+eigv(ia,ibas)*eigv(ib,ibas)
      p3(ja,jb)=p3(ja,jb)+eigv(ia,ibas)*eigv(ib,ibas)
      endif
      endif

      vec(ja,jb)=eigv(ia,ib)

      enddo
      enddo

c-----------------------------------------------------------------
c   Initialisation
c-----------------------------------------------------------------
       do ia =1,nbas(isym)
       ka=ia+nskipe(isym)
       ja=ia+2*nskipe(isym)

       do ic=1,n
       pf(ic,ka)=zero
       qf(ic,ka)=zero
       enddo

       do ib =1,nbas(isym)
       kb=ib+nskipe(isym)
       jb=ib+2*nskipe(isym)

       do ic=1,n
       pf(ic,ka)=pf(ic,ka)+vec(jb,ja)*gl(ic,kb)
       qf(ic,ka)=qf(ic,ka)+vec(jb+nbas(isym),ja)*gs(ic,kb)
       enddo

       eng1(ka)=oe(ia,isym)

       enddo
       enddo

       enddo

c      if(jpnt.ne.0)call matout(p,MNB,mbasis,mbasis,mbasis,4hP   )

c-----------------------------------------------------------------
c     Calculate delta
c-----------------------------------------------------------------

      delta=0.0d0
      do i=1,mbasis
      do j=1,mbasis
      delta=delta+(p1(i,j)-pold(i,j))**2
      enddo
      enddo
      delta=dsqrt(delta/4.0d0)
      write(STDOUT,201)iter,delta
      write(*,201)iter,delta
      if(delta.lt.crit)go to 160
      if(iter.lt.maxit)go to 100
      write(STDOUT,202)
      stop
  160 continue

      val=0.0d0
      call scfeng(cstp,ival,val)

      if(istp.gt.2)write(STDOUT,*)'Dirac-Fock SMS constant',ival,
     :val*219474.6308*0.01644

      en1=zero
      en2=zero
      en3=zero
      do i=1,mbasis
      do j=1,mbasis
      en1=en1+half*p1(i,j)*(two*df_single(i,j)+df_two(i,j))
      en2=en2+half*p2(i,j)*(two*df_single(i,j)+df_two(i,j))
      en3=en3+half*p3(i,j)*(two*df_single(i,j)+df_two(i,j))
      enddo
      enddo

      write(STDOUT,203)iter
      write(STDOUT,204)en1
      IF(IBASOPT)write(STDENG,204)en1    ! needed for basis optimization

      if(ieng.gt.0)write(STDOUT,208)en2
      if(ieng.gt.0)write(STDOUT,209)en3
      write(STDOUT,205)

      If(ihbd.eq.0.and.ivn.ne.1)then
      do jsym=1,nsym
      do ibas=1,nbas(jsym)
      write(STDOUT,206)ibas,nh_i(jsym),oe(ibas,jsym),
     :                 oe(ibas+nbas(jsym),jsym)
      enddo
      enddo
      endif

 199  format(/,1x,'CONVERGENCE DATA',/,
     1        ' MAXIMUM NO. OF ITERATION=',i6,/,
     2        ' CONVERGENCE CRITERION   = 1.0D-',i2,//)
 200  format(/,1x,'CYCLE      DENSITY CONV',/)
 201  format(i3,1x,(2x,d20.7))
 202  format(///,4x,' SCF FAILS TO CONVERGES AT CYCLE', i4)
 203  format(///,4x,22hSCF CONVERGES AT CYCLE, i4)
 204  format(///,4x,22hSCF ELECTRONIC ENERGY=, F25.15)
 205  format(///,4x,32h ORBITAL ENERGIES (+ve and -ve) ,i4,///)
 206  format(/,4x,I2,A,2F25.15)
 208  format(///,4x,22hFRZ ELECTRONIC ENERGY=, F25.15)
 209  format(///,4x,22hREM ELECTRONIC ENERGY=, F25.15)

      RETURN
      END


      subroutine dfpotiter(isym,cstp,ival)
c ********************************************************************
c This evaluates the direct and exchange matrix elements numericaly.
c P(i,j) is the density matrix.
c------------------------------------------------------------------
      implicit real*8(a-h,o-z)

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     : /cons/zero,half,tenth,one,two,three,ten
     : /skip/nskipe(MNS),nskipc(MNS)
     : /dfm/emt(MNB,MNB),df_single(MNB,MNB),df_two(MNB,MNB)
     : /quan1/kappe(MNBAS),orbje(MNBAS),kappc(MNBAS),orbjc(MNBAS)
     : /quan2/lorba(MNS),lorbb(MNS),iqe(MNBAS),iqc(MNBAS)
     : /propopt/istp,ialpvr

C New common block15
c------------------------------------------------------------------
c Local dimension
c------------------------------------------------------------------
      dimension mulk(MNS+1)

      do ia = 1,MNB
      do ib = 1,MNB
      df_two(ib,ia) = zero
      enddo
      enddo

      do ia=1,nbas(isym)
      ka=ia+nskipe(isym)
      ja=ia+2*nskipe(isym)
      index1=ia+nbas(isym)+2*nskipe(isym)

      laa=lorba(isym)
      lab=lorbb(isym)
      orba=orbje(ka)
      iaa=iqe(ka)
      kapa=kappe(ka)

      do ib=1,nbas(isym)
      index2=ib+nbas(isym)+2*nskipe(isym)
      jb=ib+2*nskipe(isym)
      kb=ib+nskipe(isym)

      lba=lorba(isym)
      lbb=lorbb(isym)
      kapb=kappe(kb)
      orbb=orbje(kb)
      iab=iqe(kb)

      df_two(ja,jb)=zero
      df_two(ja,index2)=zero
      df_two(index1,jb)=zero
      df_two(index1,index2)=zero

      do jsym=1,nsym
      nn=nocorb(jsym)

      if(nn.ne.0)then

      do ic=1,nocorb(jsym)
      jc=ic+nskipc(jsym)
      kc=ic+nskipe(jsym)

      lca=lorba(jsym)
      lcb=lorbb(jsym)
      kapc=kappc(jc)
      orbc=orbjc(jc)
      iac=iqc(jc)

c-------------------------------------------------------------------
c     direct part
c-------------------------------------------------------------------

      fact=(two*orbc+one)

      call findk(orba,orbb,orbc,orbc,mulk,maxk)

      if(maxk.ne.0) then
      do ik=1,maxk
      k1=mulk(ik)

      if(k1.eq.0)then

      rkll=zero
      rkss=zero

c-------------------------------------------------------------------
c     Direct Coulomb part
c-------------------------------------------------------------------

      call slatcould(ka,kc,kb,kc,orba,orbc,orbb,orbc,
     :iaa,iac,iab,iac,k1,rkll,rkss)

      df_two(ja,jb)=df_two(ja,jb)+rkll*fact

      df_two(index1,index2)=df_two(index1,index2)+rkss*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     exchange part
c-------------------------------------------------------------------

      call findk(orba,orbc,orbc,orbb,mulk,maxk)

      if(maxk.ne.0)then
      do ik=1,maxk
      k1=mulk(ik)
      ak1=dfloat(k1)

      rkll=zero
      rkls=zero
      rksl=zero
      rkss=zero

      bkll=zero
      bkls=zero
      bksl=zero
      bkss=zero

      skll=zero
      skls=zero
      sksl=zero
      skss=zero

      fact=(two*orbc+one)*clrx(kapa,k1,kapc)**2

c-------------------------------------------------------------------
c     Exchange Coulomb part
c-------------------------------------------------------------------

      call slatcoule(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,k1,rkll,rkls,rksl,rkss)

c-------------------------------------------------------------------
c     Exchange Breit part
c-------------------------------------------------------------------

      if(ibrt.gt.0)then
      call slatbrt(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,k1,bkll,bkls,bksl,bkss)
      endif

c-------------------------------------------------------------------
c     Exchange isotope part
c-------------------------------------------------------------------

      if(k1.eq.1)then
      if(istp.eq.1.or.istp.eq.2)then
      call slatistp(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,laa,lab,lca,lcb,lca,lcb,lba,lbb,
     :k1,skll,skls,sksl,skss,cstp)
      endif
      endif

      df_two(ja,jb)=df_two(ja,jb)-(rkll+bkll+skll)*fact

      df_two(ja,index2)=df_two(ja,index2)-(rkls+bkls+skls)*fact

      df_two(index1,jb)=df_two(index1,jb)-(rksl+bksl+sksl)*fact

      df_two(index1,index2)=df_two(index1,index2)-
     :(rkss+bkss+skss)*fact

      enddo
      endif

      enddo
      endif
      enddo

      if(ivn.gt.0)then

      do jsym=1,nsym
      do ic=1,nbas(jsym)
      kc=ic+nskipe(jsym)
      orbc=orbje(kc)
      iac=iqe(kc)
      kapc=kappe(kc)

      if(kc.eq.ival)then

c-------------------------------------------------------------------
c     direct part
c-------------------------------------------------------------------

      sgn=(two*orbc+one)*(two*orba+one)

      call findk(orba,orbb,orbc,orbc,mulk,maxk)

      if(maxk.gt.0)then
      do ik=1,maxk
      k1=mulk(ik)
      ak1=dfloat(k1)

      rkll=zero
      rkss=zero

cC mj=orba
c      f1=dr(orba,ak1,orbb,-orba,zero,orbb)
c      f2=dr(orbc,ak1,orbc,-orbc,zero,orbc)
c      f3=dr(orba,ak1,orbb,half,zero,-half)
c      f4=dr(orbc,ak1,orbc,half,zero,-half)
c      fact=sgn*f1*f2*f3*f4*(-1)**(one+orba+orbc)

C mj=half
      f1=dr(orba,ak1,orbb,-half,zero,half)
      f2=dr(orbc,ak1,orbc,-half,zero,half)
      f3=dr(orba,ak1,orbb,half,zero,-half)
      f4=dr(orbc,ak1,orbc,half,zero,-half)
      fact=sgn*f1*f2*f3*f4

cC Average mj
c      if(k1.eq.0)then
c      fact=one
c      if(ka.eq.ival.or.kb.eq.ival)fact=fact/(2*orba+one)

      call slatcould(ka,kc,kb,kc,orba,orbc,orbb,orbc,
     :iaa,iac,iab,iac,k1,rkll,rkss)

      df_two(ja,jb)=df_two(ja,jb)+rkll*fact

      df_two(index1,index2)=df_two(index1,index2)+rkss*fact

c      endif

      enddo
      endif

c-------------------------------------------------------------------
c     exchange part
c-------------------------------------------------------------------

      call findk(orba,orbc,orbc,orbb,mulk,maxk)

      if(maxk.gt.0)then
      do ik=1,maxk
      k1=mulk(ik)
      ak1=dfloat(k1)

      rkll=zero
      rkls=zero
      rksl=zero
      rkss=zero

      bkll=zero
      bkls=zero
      bksl=zero
      bkss=zero

cC mj=orba
c      f1=dr(orba,ak1,orbc,-orba,orba-orbc,orbc)
c      f2=dr(orbc,ak1,orbb,-orbc,orbc-orbb,orbb)
c      f3=dr(orba,ak1,orbc,half,zero,-half)
c      f4=dr(orbc,ak1,orbb,half,zero,-half)
c      fact=sgn*f1*f2*f3*f4

C mj=half
      f1=dr(orba,ak1,orbc,-half,zero,half)
      f2=dr(orbc,ak1,orbb,-half,zero,half)
      f3=dr(orba,ak1,orbc,half,zero,-half)
      f4=dr(orbc,ak1,orbb,half,zero,-half)
      fact=sgn*f1*f2*f3*f4

cC Average mj
c      fact=clrx(kapa,k1,kapc)**2
c      if(ka.eq.ival.or.kb.eq.ival)fact=fact/(2*orba+one)

      call slatcoule(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,k1,rkll,rkls,rksl,rkss)

      if(ibrt.gt.0)then
      call slatbrt(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,k1,bkll,bkls,bksl,bkss)
      endif

      df_two(ja,jb)=df_two(ja,jb)-(rkll+bkll)*fact

      df_two(ja,index2)=df_two(ja,index2)-(rkls+bkls)*fact

      df_two(index1,jb)=df_two(index1,jb)-(rksl+bksl)*fact

      df_two(index1,index2)=df_two(index1,index2)-(rkss+bkss)*fact

      enddo
      endif

      endif

      enddo
      enddo

      endif

      enddo
      enddo

      return
      end


C This subroutine selects possible multipoles
      subroutine findk(orja,orjb,orjc,orjd,mulk,maxk)
      implicit real*8 (a-h,o-z)

      dimension mulk(MNS+1)

      Do ik=1,MNS+1
      mulk(ik)=0
      EndDo

      kmin1=idint(dabs(orja-orjb))
      kmax1=idint(orja+orjb)

      kmin2=idint(dabs(orjc-orjd))
      kmax2=idint(orjc+orjd)

      kcount=0
      maxk=0

      Do k1=max0(kmin1,kmin2),min0(kmax1,kmax2)

            kcount=kcount+1
            mulk(kcount)=k1

      Enddo

      maxk=kcount

      return
      end


      FUNCTION C3J(XJ1,XJ2,XJ3,XM1,XM2,XM3)
      implicit real*8(a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten

      DIMENSION FLOG(301)                                               TRO00020
      DATA(FLOG(I),I=2,31)/0.,.69314718,1.7917595,3.1780538,4.7874917,6.TRO00060
     15792511,8.5251613,10.604603,12.801827,15.104413,17.502307,19.98721TRO00070
     24,22.552163,25.191221,27.899271,30.671860,33.505072,36.395445,39.3TRO00080
     339884,42.335616,45.380139,48.471180,51.606674,54.784729,58.003604,TRO00090
     461.261702,64.557537,67.889743,71.257038,74.658235/                TRO00100
      DATA(FLOG(I),I=32,61)/78.092223,81.557959,85.054466,88.580827,92.1TRO00110
     136175,95.719694,99.330612,102.96820,106.63176,110.32064,114.03421,TRO00120
     2117.77188,121.53308,125.31727,129.12393,132.95257,136.80272,140.67TRO00130
     3392,144.56574,148.47776,152.40959,156.36083,160.33112,164.32011,16TRO00140
     48.32744,172.35279,176.39584,180.45629,184.53383,188.62817/        TRO00150
      DATA(FLOG(I),I=62,91)/192.73904,196.86618,201.00931,205.16820,209.TRO00160
     134258,213.53224,217.73693,221.95644,226.19054,230.43904,234.70172,TRO00170
     2238.97839,243.26885,247.57291,251.89040,256.22113,260.56494,264.92TRO00180
     3164,269.29110,273.67312,278.06757,282.47429,286.89313,291.32394,29TRO00190
     45.76659,300.22094,304.68685,309.16419,313.65283,318.15264/        TRO00200
      DATA(FLOG(I),I=92,121)/322.66349,327.18529,331.71788,336.26118,340TRO00210
     1.81505,345.37940,349.95411,354.53908,359.13420,363.73937,368.35449TRO00220
     2,372.97946,377.61419,382.25859,386.91255,391.57598,396.24881,400.9TRO00230
     33094,405.62230,410.32277,415.03230,419.75080,424.47819,429.21439,4TRO00240
     433.95932,438.71291,443.47508,448.24576,453.02489,457.81238/       TRO00250
      DATA(FLOG(I),I=122,151)/462.60817,467.41220,472.22438,477.04466,48TRO00260
     11.87298,486.70926,491.55345,496.40547,501.26529,506.13282,511.0080TRO00270
     22,515.89082,520.78117,525.67901,530.58428,535.49694,540.41692,545.TRO00280
     334417,550.27865,555.22029,560.16905,565.12488,570.08772,575.05753,TRO00290
     4580.03427,585.01787,590.00830,595.00552,600.00946,605.02010/      TRO00300
      DATA(FLOG(I),I=152,181)/610.03738,615.06126,620.09170,625.12866,63TRO00310
     10.17208,635.22193,640.27818,645.34077,650.40968,655.48486,660.5662TRO00320
     26,665.65385,670.74760,675.84747,680.95341,686.06541,691.18340,696.TRO00330
     330735,701.43726,706.57306,711.71472,716.86221,722.01551,727.17456,TRO00340
     4732.33934,737.50983,742.68598,747.86776,753.05516,758.24811/      TRO00350
      DATA(FLOG(I),I=182,211)/763.44661,768.65061,773.86010,779.07503,78TRO00360
     14.29539,789.52114,794.75224,799.98869,805.23044,810.47747,815.7297TRO00370
     23,820.98722,826.24991,831.51778,836.79078,842.06890,847.35209,852.TRO00380
     364036,857.93366,863.23199,868.53529,873.84356,879.15676,884.47488,TRO00390
     4889.79789,895.12577,900.45848,905.79603,911.13836,916.48547/      TRO00400
      DATA(FLOG(I),I=212,241)/921.83732,927.19391,932.55521,937.92118,94TRO00410
     13.29181,948.66710,954.04699,959.43148,964.82056,970.21419,975.6123TRO00420
     25,981.01503,986.42220,991.83385,997.24995,1002.6705,1008.0954,1013TRO00430
     3.5248,1018.9585,1024.3966,1029.8389,1035.2857,1040.7367,1046.1920,TRO00440
     41051.6516,1057.1155,1062.5836,1068.0558,1073.5323,1079.0129/      TRO00450
      DATA(FLOG(I),I=242,271)/1084.4977,1089.9866,1095.4797,1100.9768,11TRO00460
     106.4781,1111.9834,1117.4928,1123.0063,1128.5237,1134.0452,1139.570TRO00470
     26,1145.1001,1150.6335,1156.1708,1161.7120,1167.2573,1172.8063,1178TRO00480
     3.3593,1183.9161,1189.4768,1195.0413,1200.6097,1206.1818,1211.7577,TRO00490
     41217.3375,1222.9209,1228.5082,1234.0992,1239.6939,1245.2924/      TRO00500
      DATA(FLOG(I),I=272,301)/1250.8944,1256.5003,1262.1097,1267.7228,12TRO00510
     173.3396,1278.9600,1284.5840,1290.2117,1295.8429,1301.4777,1307.116TRO00520
     20,1312.7580,1318.4034,1324.0524,1329.7048,1335.3609,1341.0203,1346TRO00530
     3.6833,1352.3497,1358.0196,1363.6929,1369.3697,1375.0499,1380.7334,TRO00540
     41386.4204,1392.1107,1397.8045,1403.5016,1409.2020,1414.9058/      TRO00550
      DATA EPS1,EPS2/.1,-.2/                                            TRO00560
      XN=XJ2-XM2+EPS1                                                   TRO00600
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO00610
      N1=XN                                                             TRO00620
      XN=XJ3+XM3+EPS1                                                   TRO00630
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO00640
      N2=XN                                                             TRO00650
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO00660
      XN=XJ3-XM3+EPS1                                                   TRO00670
      N3=XN                                                             TRO00680
      XN=XJ1+XM1+EPS1                                                   TRO00690
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO00700
      N4=XN                                                             TRO00710
      XN=XJ2+XJ3-XJ1+EPS1                                               TRO00720
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO00730
      N5=XN                                                             TRO00740
      XN=XJ1+XJ3-XJ2+EPS1                                               TRO00750
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO00760
      N6=XN                                                             TRO00770
      XN=XJ1+XJ2-XJ3+EPS1                                               TRO00780
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO00790
      N7=XN                                                             TRO00800
      XN=XJ1-XM1+EPS1                                                   TRO00810
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO00820
      N8=XN                                                             TRO00830
      XN=XJ2+XM2+EPS1                                                   TRO00840
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO00850
      N9=XN                                                             TRO00860
      XN=XJ3-XJ1-XM2+EPS1                                               TRO00870
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO00880
      N10=XN                                                            TRO00890
      XN=XJ3-XJ2+XM1+EPS1                                               TRO00900
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO00910
      N11=XN                                                            TRO00920
      XN=XJ1+XJ2+XJ3+EPS1                                               TRO00930
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO00940
      N12=XN                                                            TRO00950
      N12=N12+1                                                         TRO00960
      K=N4*N8                                                           TRO01000
      IF(K.LT.0) GO TO 59                                               TRO01010
      K=N1*N9                                                           TRO01020
      IF(K.LT.0) GO TO 59                                               TRO01030
      K=N2*N3                                                           TRO01040
      IF(K.LT.0) GO TO 59                                               TRO01050
      IF(N5.LT.0) GO TO 59                                              TRO01060
      IF(N6.LT.0) GO TO 59                                              TRO01070
      IF(N7.LT.0) GO TO 59                                              TRO01080
      L=N1-N2+N3-N4+N8-N9                                               TRO01090
      IF(L.NE.0) GO TO 59                                               TRO01100
      K=N12-1                                                           TRO01110
      IF(K.GT.0) GO TO 60                                               TRO01120
      C3J=ONE                                                           TRO01130
      RETURN                                                            TRO01140
   59 C3J=ZERO                                                          TRO01150
      RETURN                                                            TRO01160
   60 K=0                                                               TRO01200
      L=-N10                                                            TRO01210
      IF(L.GT.K) K=L                                                    TRO01220
      L=-N11                                                            TRO01230
      IF(L.GT.K) K=L                                                    TRO01240
      L=N7                                                              TRO01250
      IF(N8.GT.L) L=N8                                                  TRO01260
      IF(N9.GT.L) L=N9                                                  TRO01270
      F=ONE                                                             TRO01280
      S=ONE                                                             TRO01290
      I=K+1                                                             TRO01300
   62 IF(I.GT.L) GO TO 80                                               TRO01310
      IM1=I-1                                                           TRO01320
      NN=(N7-IM1)*(N8-IM1)*(N9-IM1)                                     TRO01330
      ND=I*(N10+I)*(N11+I)                                              TRO01340
      F=-F*DFLOAT(NN)/DFLOAT(ND)                                        TRO01350
      S=S+F                                                             TRO01360
      I=I+1                                                             TRO01370
      GO TO 62                                                          TRO01380
   80 C2N=FLOG(N1+1)+FLOG(N2+1)+FLOG(N3+1)+FLOG(N4+1)+FLOG(N5+1)+FLOG(N6TRO01420
     1+1)+FLOG(N7+1)+FLOG(N8+1)+FLOG(N9+1)                              TRO01430
      C2N=.5*C2N                                                        TRO01440
      KM1=K-1                                                           TRO01450
      KP1=K+1                                                           TRO01460
      C2D=FLOG(KP1)+FLOG(N7-KM1)+FLOG(N8-KM1)+FLOG(N9-KM1)+FLOG(N10+KP1)TRO01470
     1+FLOG(N11+KP1)+.5*FLOG(N12+1)                                     TRO01480
      F=C2D-C2N                                                         TRO01520
      IF(F.GT.80.) GO TO 98                                             TRO01530
      F=C2N/C2D                                                         TRO01540
      IF((F.LT.1.01).AND.(F.GT.0.98)) GO TO 98                          TRO01550
      C3J=S*EXP(C2N-C2D)                                                TRO01560
      GO TO 106                                                         TRO01570
   98 IF(S) 100,59,102                                                  TRO01580
  100 S=DLOG(-S)                                                        TRO01590
      C3J=-DEXP(S+C2N-C2D)                                              TRO01600
      GO TO 106                                                         TRO01610
  102 S=DLOG(S)                                                         TRO01620
      C3J=DEXP(S+C2N-C2D)                                               TRO01630
  106 XN=XJ1-XJ2-XM3+EPS1                                               TRO01670
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO01680
      L=XN                                                              TRO01690
      L=L+K                                                             TRO01700
      K=L/2                                                             TRO01710
      K=2*K                                                             TRO01720

      IF(L.NE.K) C3J=-C3J                                               TRO01730

      RETURN                                                            TRO01740
      END                                                               TRO01750


      FUNCTION C6J(XJ1,XJ2,XJ3,XL1,XL2,XL3)
      implicit real*8(a-h,o-z)

      DIMENSION FLOG(301)                                               TRO01780

      DATA(FLOG(I),I=2,31)/0.,.69314718,1.7917595,3.1780538,4.7874917,6.TRO01820
     15792511,8.5251613,10.604603,12.801827,15.104413,17.502307,19.98721TRO01830
     24,22.552163,25.191221,27.899271,30.671860,33.505072,36.395445,39.3TRO01840
     339884,42.335616,45.380139,48.471180,51.606674,54.784729,58.003604,TRO01850
     461.261702,64.557537,67.889743,71.257038,74.658235/                TRO01860
      DATA(FLOG(I),I=32,61)/78.092223,81.557959,85.054466,88.580827,92.1TRO01870
     136175,95.719694,99.330612,102.96820,106.63176,110.32064,114.03421,TRO01880
     2117.77188,121.53308,125.31727,129.12393,132.95257,136.80272,140.67TRO01890
     3392,144.56574,148.47776,152.40959,156.36083,160.33112,164.32011,16TRO01900
     48.32744,172.35279,176.39584,180.45629,184.53383,188.62817/        TRO01910
      DATA(FLOG(I),I=62,91)/192.73904,196.86618,201.00931,205.16820,209.TRO01920
     134258,213.53224,217.73693,221.95644,226.19054,230.43904,234.70172,TRO01930
     2238.97839,243.26885,247.57291,251.89040,256.22113,260.56494,264.92TRO01940
     3164,269.29110,273.67312,278.06757,282.47429,286.89313,291.32394,29TRO01950
     45.76659,300.22094,304.68685,309.16419,313.65283,318.15264/        TRO01960
      DATA(FLOG(I),I=92,121)/322.66349,327.18529,331.71788,336.26118,340TRO01970
     1.81505,345.37940,349.95411,354.53908,359.13420,363.73937,368.35449TRO01980
     2,372.97946,377.61419,382.25859,386.91255,391.57598,396.24881,400.9TRO01990
     33094,405.62230,410.32277,415.03230,419.75080,424.47819,429.21439,4TRO02000
     433.95932,438.71291,443.47508,448.24576,453.02489,457.81238/       TRO02010
      DATA(FLOG(I),I=122,151)/462.60817,467.41220,472.22438,477.04466,48TRO02020
     11.87298,486.70926,491.55345,496.40547,501.26529,506.13282,511.0080TRO02030
     22,515.89082,520.78117,525.67901,530.58428,535.49694,540.41692,545.TRO02040
     334417,550.27865,555.22029,560.16905,565.12488,570.08772,575.05753,TRO02050
     4580.03427,585.01787,590.00830,595.00552,600.00946,605.02010/      TRO02060
      DATA(FLOG(I),I=152,181)/610.03738,615.06126,620.09170,625.12866,63TRO02070
     10.17208,635.22193,640.27818,645.34077,650.40968,655.48486,660.5662TRO02080
     26,665.65385,670.74760,675.84747,680.95341,686.06541,691.18340,696.TRO02090
     330735,701.43726,706.57306,711.71472,716.86221,722.01551,727.17456,TRO02100
     4732.33934,737.50983,742.68598,747.86776,753.05516,758.24811/      TRO02110
      DATA(FLOG(I),I=182,211)/763.44661,768.65061,773.86010,779.07503,78TRO02120
     14.29539,789.52114,794.75224,799.98869,805.23044,810.47747,815.7297TRO02130
     23,820.98722,826.24991,831.51778,836.79078,842.06890,847.35209,852.TRO02140
     364036,857.93366,863.23199,868.53529,873.84356,879.15676,884.47488,TRO02150
     4889.79789,895.12577,900.45848,905.79603,911.13836,916.48547/      TRO02160
      DATA(FLOG(I),I=212,241)/921.83732,927.19391,932.55521,937.92118,94TRO02170
     13.29181,948.66710,954.04699,959.43148,964.82056,970.21419,975.6123TRO02180
     25,981.01503,986.42220,991.83385,997.24995,1002.6705,1008.0954,1013TRO02190
     3.5248,1018.9585,1024.3966,1029.8389,1035.2857,1040.7367,1046.1920,TRO02200
     41051.6516,1057.1155,1062.5836,1068.0558,1073.5323,1079.0129/      TRO02210
      DATA(FLOG(I),I=242,271)/1084.4977,1089.9866,1095.4797,1100.9768,11TRO02220
     106.4781,1111.9834,1117.4928,1123.0063,1128.5237,1134.0452,1139.570TRO02230
     26,1145.1001,1150.6335,1156.1708,1161.7120,1167.2573,1172.8063,1178TRO02240
     3.3593,1183.9161,1189.4768,1195.0413,1200.6097,1206.1818,1211.7577,TRO02250
     41217.3375,1222.9209,1228.5082,1234.0992,1239.6939,1245.2924/      TRO02260
      DATA(FLOG(I),I=272,301)/1250.8944,1256.5003,1262.1097,1267.7228,12TRO02270
     173.3396,1278.9600,1284.5840,1290.2117,1295.8429,1301.4777,1307.116TRO02280
     20,1312.7580,1318.4034,1324.0524,1329.7048,1335.3609,1341.0203,1346TRO02290
     3.6833,1352.3497,1358.0196,1363.6929,1369.3697,1375.0499,1380.7334,TRO02300
     41386.4204,1392.1107,1397.8045,1403.5016,1409.2020,1414.9058/      TRO02310
      DATA EPS1,EPS2/.1,-.2/                                            TRO02320
      XN=-XJ1+XJ2+XJ3+EPS1                                              TRO02360
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02370
      N1=XN                                                             TRO02380
      XN=-XL1+XL2+XJ3+EPS1                                              TRO02390
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02400
      N2=XN                                                             TRO02410
      XN=-XL1+XJ2+XL3+EPS1                                              TRO02420
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02430
      N3=XN                                                             TRO02440
      XN=-XJ1+XL2+XL3+EPS1                                              TRO02450
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02460
      N4=XN                                                             TRO02470
      XN=XJ1-XJ2+XJ3+EPS1                                               TRO02480
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02490
      N5=XN                                                             TRO02500
      XN=XL1-XL2+XJ3+EPS1                                               TRO02510
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02520
      N6=XN                                                             TRO02530
      XN=XL1-XJ2+XL3+EPS1                                               TRO02540
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02550
      N7=XN                                                             TRO02560
      XN=XJ1-XL2+XL3+EPS1                                               TRO02570
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02580
      N8=XN                                                             TRO02590
      XN=XJ1+XJ2-XJ3+EPS1                                               TRO02600
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02610
      N9=XN                                                             TRO02620
      XN=XL1+XL2-XJ3+EPS1                                               TRO02630
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02640
      N10=XN                                                            TRO02650
      XN=XL1+XJ2-XL3+EPS1                                               TRO02660
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02670
      N11=XN                                                            TRO02680
      XN=XJ1+XL2-XL3+EPS1                                               TRO02690
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02700
      N12=XN                                                            TRO02710
      XN=-XJ1-XL1+XJ3+XL3+EPS1                                          TRO02720
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02730
      N13=XN                                                            TRO02740
      XN=-XJ2-XL2+XJ3+XL3+EPS1                                          TRO02750
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02760
      N14=XN                                                            TRO02770
      XN=XJ1+XL1+XJ2+XL2+EPS1                                           TRO02780
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02790
      N15=XN                                                            TRO02800
      N15=N15+1                                                         TRO02810
      XN=XJ1+XJ2+XJ3+EPS1                                               TRO02820
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02830
      N16=XN                                                            TRO02840
      N16=N16+1                                                         TRO02850
      XN=XL1+XL2+XJ3+EPS1                                               TRO02860
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02870
      N17=XN                                                            TRO02880
      N17=N17+1                                                         TRO02890
      XN=XL1+XJ2+XL3+EPS1                                               TRO02900
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02910
      N18=XN                                                            TRO02920
      N18=N18+1                                                         TRO02930
      XN=XJ1+XL2+XL3+EPS1                                               TRO02940
      IF (XN.LT.0.) XN=XN+EPS2                                          TRO02950
      N19=XN                                                            TRO02960
      N19=N19+1                                                         TRO02970
      IF(N9.LT.0) GO TO 50                                              TRO03010
      IF(N5.LT.0) GO TO 50                                              TRO03020
      IF(N1.LT.0) GO TO 50                                              TRO03030
      IF(N10.LT.0) GO TO 50                                             TRO03040
      IF(N6.LT.0) GO TO 50                                              TRO03050
      IF(N2.LT.0) GO TO 50                                              TRO03060
      IF(N11.LT.0) GO TO 50                                             TRO03070
      IF(N7.LT.0) GO TO 50                                              TRO03080
      IF(N3.LT.0) GO TO 50                                              TRO03090
      IF(N12.LT.0) GO TO 50                                             TRO03100
      IF(N8.LT.0) GO TO 50                                              TRO03110
      IF(N4.LT.0) GO TO 50                                              TRO03120
      K=N17+N18+N19-3                                                   TRO03130
      IF(K.GT.0) GO TO 54                                               TRO03140
      C6J=ONE                                                           TRO03150
      RETURN                                                            TRO03160
   50 C6J=ZERO                                                          TRO03170
      RETURN                                                            TRO03180
   54 K=0                                                               TRO03220
      L=-N13                                                            TRO03230
      IF(L.GT.K) K=L                                                    TRO03240
      L=-N14                                                            TRO03250
      IF(L.GT.K) K=L                                                    TRO03260
      L=N9                                                              TRO03270
      IF(N10.LT.L) L=N10                                                TRO03280
      IF(N11.LT.L) L=N11                                                TRO03290
      IF(N12.LT.L) L=N12                                                TRO03300
      IF(N15.LT.L) L=N15                                                TRO03310
      F=ONE                                                             TRO03320
      S=ONE                                                             TRO03330
      I=K+1                                                             TRO03340
   62 IF(I.GT.L) GO TO 80                                               TRO03350
      IM1=I-1                                                           TRO03360
      NN=(N9-IM1)*(N10-IM1)*(N11-IM1)*(N12-IM1)                         TRO03370
      ND=I*(N13+I)*(N14+I)*(N15-IM1)                                    TRO03380
      F=-F*DFLOAT(NN)/DFLOAT(ND)                                        TRO03390
      S=S+F                                                             TRO03400
      I=I+1                                                             TRO03410
      GO TO 62                                                          TRO03420
   80 C2N=FLOG(N1+1)+FLOG(N2+1)+FLOG(N3+1)+FLOG(N4+1)+FLOG(N5+1)+FLOG(N6TRO03460
     1+1)+FLOG(N7+1)+FLOG(N8+1)+FLOG(N9+1)+FLOG(N10+1)+FLOG(N11+1)+FLOG(TRO03470
     2N12+1)                                                            TRO03480
      C2N=HALF*C2N                                                      TRO03490
      C2D=FLOG(N16+1)+FLOG(N17+1)+FLOG(N18+1)+FLOG(N19+1)               TRO03500
      C2D=HALF*C2D                                                      TRO03510
      KM1=K-1                                                           TRO03520
      KP1=K+1                                                           TRO03530
      C2N=C2N+FLOG(N15-KM1)                                             TRO03540
      C2D=C2D+FLOG(KP1)+FLOG(N13+KP1)+FLOG(N14+KP1)+FLOG(N9-KM1)+FLOG(N1TRO03550
     10-KM1)+FLOG(N11-KM1)+FLOG(N12-KM1)                                TRO03560
      F=C2D-C2N                                                         TRO03600
      IF(F.GT.80.) GO TO 98                                             TRO03610
      F=C2N/C2D                                                         TRO03620
      IF((F.LT.1.01).AND.(F.GT.0.98)) GO TO 98                          TRO03630
      C6J=S*DEXP(C2N-C2D)                                               TRO03640
      GO TO 106                                                         TRO03650
   98 IF(S) 100,50,102                                                  TRO03660
  100 S=DLOG(-S)                                                        TRO03670
      C6J=-DEXP(S+C2N-C2D)                                              TRO03680
      GO TO 106                                                         TRO03690
  102 S=DLOG(S)                                                         TRO03700
      C6J=DEXP(S+C2N-C2D)                                               TRO03710
  106 L=N15+KM1                                                         TRO03750
      K=L/2                                                             TRO03760
      K=2*K                                                             TRO03770
      IF(L.NE.K) C6J=-C6J                                               TRO03780

      RETURN                                                            TRO03790
      END                                                               TRO03800


      FUNCTION C9J(XJ11,XJ12,XJ13,XJ21,XJ22,XJ23,XJ31,XJ32,XJ33)        TRO03820
      implicit real*8(a-h,o-z)

      DIMENSION FLOG(301)                                               TRO03830
      DATA(FLOG(I),I=2,31)/0.,.69314718,1.7917595,3.1780538,4.7874917,6.TRO03870
     15792511,8.5251613,10.604603,12.801827,15.104413,17.502307,19.98721TRO03880
     24,22.552163,25.191221,27.899271,30.671860,33.505072,36.395445,39.3TRO03890
     339884,42.335616,45.380139,48.471180,51.606674,54.784729,58.003604,TRO03900
     461.261702,64.557537,67.889743,71.257038,74.658235/                TRO03910
      DATA(FLOG(I),I=32,61)/78.092223,81.557959,85.054466,88.580827,92.1TRO03920
     136175,95.719694,99.330612,102.96820,106.63176,110.32064,114.03421,TRO03930
     2117.77188,121.53308,125.31727,129.12393,132.95257,136.80272,140.67TRO03940
     3392,144.56574,148.47776,152.40959,156.36083,160.33112,164.32011,16TRO03950
     48.32744,172.35279,176.39584,180.45629,184.53383,188.62817/        TRO03960
      DATA(FLOG(I),I=62,91)/192.73904,196.86618,201.00931,205.16820,209.TRO03970
     134258,213.53224,217.73693,221.95644,226.19054,230.43904,234.70172,TRO03980
     2238.97839,243.26885,247.57291,251.89040,256.22113,260.56494,264.92TRO03990
     3164,269.29110,273.67312,278.06757,282.47429,286.89313,291.32394,29TRO04000
     45.76659,300.22094,304.68685,309.16419,313.65283,318.15264/        TRO04010
      DATA(FLOG(I),I=92,121)/322.66349,327.18529,331.71788,336.26118,340TRO04020
     1.81505,345.37940,349.95411,354.53908,359.13420,363.73937,368.35449TRO04030
     2,372.97946,377.61419,382.25859,386.91255,391.57598,396.24881,400.9TRO04040
     33094,405.62230,410.32277,415.03230,419.75080,424.47819,429.21439,4TRO04050
     433.95932,438.71291,443.47508,448.24576,453.02489,457.81238/       TRO04060
      DATA(FLOG(I),I=122,151)/462.60817,467.41220,472.22438,477.04466,48TRO04070
     11.87298,486.70926,491.55345,496.40547,501.26529,506.13282,511.0080TRO04080
     22,515.89082,520.78117,525.67901,530.58428,535.49694,540.41692,545.TRO04090
     334417,550.27865,555.22029,560.16905,565.12488,570.08772,575.05753,TRO04100
     4580.03427,585.01787,590.00830,595.00552,600.00946,605.02010/      TRO04110
      DATA(FLOG(I),I=152,181)/610.03738,615.06126,620.09170,625.12866,63TRO04120
     10.17208,635.22193,640.27818,645.34077,650.40968,655.48486,660.5662TRO04130
     26,665.65385,670.74760,675.84747,680.95341,686.06541,691.18340,696.TRO04140
     330735,701.43726,706.57306,711.71472,716.86221,722.01551,727.17456,TRO04150
     4732.33934,737.50983,742.68598,747.86776,753.05516,758.24811/      TRO04160
      DATA(FLOG(I),I=182,211)/763.44661,768.65061,773.86010,779.07503,78TRO04170
     14.29539,789.52114,794.75224,799.98869,805.23044,810.47747,815.7297TRO04180
     23,820.98722,826.24991,831.51778,836.79078,842.06890,847.35209,852.TRO04190
     364036,857.93366,863.23199,868.53529,873.84356,879.15676,884.47488,TRO04200
     4889.79789,895.12577,900.45848,905.79603,911.13836,916.48547/      TRO04210
      DATA(FLOG(I),I=212,241)/921.83732,927.19391,932.55521,937.92118,94TRO04220
     13.29181,948.66710,954.04699,959.43148,964.82056,970.21419,975.6123TRO04230
     25,981.01503,986.42220,991.83385,997.24995,1002.6705,1008.0954,1013TRO04240
     3.5248,1018.9585,1024.3966,1029.8389,1035.2857,1040.7367,1046.1920,TRO04250
     41051.6516,1057.1155,1062.5836,1068.0558,1073.5323,1079.0129/      TRO04260
      DATA(FLOG(I),I=242,271)/1084.4977,1089.9866,1095.4797,1100.9768,11TRO04270
     106.4781,1111.9834,1117.4928,1123.0063,1128.5237,1134.0452,1139.570TRO04280
     26,1145.1001,1150.6335,1156.1708,1161.7120,1167.2573,1172.8063,1178TRO04290
     3.3593,1183.9161,1189.4768,1195.0413,1200.6097,1206.1818,1211.7577,TRO04300
     41217.3375,1222.9209,1228.5082,1234.0992,1239.6939,1245.2924/      TRO04310
      DATA(FLOG(I),I=272,301)/1250.8944,1256.5003,1262.1097,1267.7228,12TRO04320
     173.3396,1278.9600,1284.5840,1290.2117,1295.8429,1301.4777,1307.116TRO04330
     20,1312.7580,1318.4034,1324.0524,1329.7048,1335.3609,1341.0203,1346TRO04340
     3.6833,1352.3497,1358.0196,1363.6929,1369.3697,1375.0499,1380.7334,TRO04350
     41386.4204,1392.1107,1397.8045,1403.5016,1409.2020,1414.9058/      TRO04360
      DATA EPS1,EPS2/.1,-.2/                                            TRO04370
      XN=-XJ11+XJ21+XJ31+EPS1                                           TRO04410
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04420
      N1=XN                                                             TRO04430
      XN=-XJ32+XJ33+XJ31+EPS1                                           TRO04440
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04450
      N2=XN                                                             TRO04460
      XN=XJ11-XJ21+XJ31+EPS1                                            TRO04470
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04480
      N5=XN                                                             TRO04490
      XN=XJ32-XJ33+XJ31+EPS1                                            TRO04500
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04510
      N6=XN                                                             TRO04520
      XN=XJ11+XJ21-XJ31+EPS1                                            TRO04530
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04540
      N9=XN                                                             TRO04550
      XN=XJ32+XJ33-XJ31+EPS1                                            TRO04560
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04570
      N10=XN                                                            TRO04580
      XN=XJ11+XJ32+XJ21+XJ33+EPS1                                       TRO04590
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04600
      N15=XN                                                            TRO04610
      N15=N15+1                                                         TRO04620
      XN=XJ11+XJ21+XJ31+EPS1                                            TRO04630
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04640
      N16=XN                                                            TRO04650
      N16=N16+1                                                         TRO04660
      XN=XJ32+XJ33+XJ31+EPS1                                            TRO04670
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04680
      N17=XN                                                            TRO04690
      N17=N17+1                                                         TRO04700
      XN=-XJ12+XJ22+XJ32+EPS1                                           TRO04710
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04720
      N21=XN                                                            TRO04730
      XN=-XJ21+XJ22+XJ23+EPS1                                           TRO04740
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04750
      N23=XN                                                            TRO04760
      XN=XJ12-XJ22+XJ32+EPS1                                            TRO04770
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04780
      N25=XN                                                            TRO04790
      XN=XJ21-XJ22+XJ23+EPS1                                            TRO04800
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04810
      N27=XN                                                            TRO04820
      XN=XJ12+XJ22-XJ32+EPS1                                            TRO04830
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04840
      N29=XN                                                            TRO04850
      XN=XJ21+XJ22-XJ23+EPS1                                            TRO04860
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04870
      N31=XN                                                            TRO04880
      XN=-XJ12-XJ21+XJ32+XJ23+EPS1                                      TRO04890
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04900
      N33=XN                                                            TRO04910
      XN=XJ12+XJ22+XJ32+EPS1                                            TRO04920
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04930
      N36=XN                                                            TRO04940
      N36=N36+1                                                         TRO04950
      XN=XJ21+XJ22+XJ23+EPS1                                            TRO04960
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO04970
      N38=XN                                                            TRO04980
      N38=N38+1                                                         TRO04990
      XN=-XJ13+XJ23+XJ33+EPS1                                           TRO05000
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05010
      N41=XN                                                            TRO05020
      XN=-XJ13+XJ11+XJ12+EPS1                                           TRO05030
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05040
      N44=XN                                                            TRO05050
      XN=XJ13-XJ23+XJ33+EPS1                                            TRO05060
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05070
      N45=XN                                                            TRO05080
      XN=XJ13-XJ11+XJ12+EPS1                                            TRO05090
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05100
      N48=XN                                                            TRO05110
      XN=XJ13+XJ23-XJ33+EPS1                                            TRO05120
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05130
      N49=XN                                                            TRO05140
      XN=XJ13+XJ11-XJ12+EPS1                                            TRO05150
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05160
      N52=XN                                                            TRO05170
      XN=-XJ23-XJ11+XJ33+XJ12+EPS1                                      TRO05180
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05190
      N54=XN                                                            TRO05200
      XN=XJ13+XJ23+XJ33+EPS1                                            TRO05210
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05220
      N56=XN                                                            TRO05230
      N56=N56+1                                                         TRO05240
      XN=XJ13+XJ11+XJ12+EPS1                                            TRO05250
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05260
      N59=XN                                                            TRO05270
      N59=N59+1                                                         TRO05280
      IF(N9.LT.0) GO TO 50                                              TRO05320
      IF(N5.LT.0) GO TO 50                                              TRO05330
      IF(N1.LT.0) GO TO 50                                              TRO05340
      IF(N10.LT.0) GO TO 50                                             TRO05350
      IF(N6.LT.0) GO TO 50                                              TRO05360
      IF(N2.LT.0) GO TO 50                                              TRO05370
      IF(N29.LT.0) GO TO 50                                             TRO05380
      IF(N25.LT.0) GO TO 50                                             TRO05390
      IF(N21.LT.0) GO TO 50                                             TRO05400
      IF(N31.LT.0) GO TO 50                                             TRO05410
      IF(N27.LT.0) GO TO 50                                             TRO05420
      IF(N23.LT.0) GO TO 50                                             TRO05430
      IF(N49.LT.0) GO TO 50                                             TRO05440
      IF(N45.LT.0) GO TO 50                                             TRO05450
      IF(N41.LT.0) GO TO 50                                             TRO05460
      IF(N52.LT.0) GO TO 50                                             TRO05470
      IF(N48.LT.0) GO TO 50                                             TRO05480
      IF(N44.LT.0) GO TO 50                                             TRO05490
      K=N1+N2+N5+N6+N9+N10+N21+N23+N25+N27+N29+N31+N41+N44+N45+N48+N49+NTRO05500
     152                                                                TRO05510
      IF(K.GT.0) GO TO 54                                               TRO05520
      C9J=ONE                                                           TRO05530
      RETURN                                                            TRO05540
  50  C9J=ZERO                                                          TRO05550
      RETURN                                                            TRO05560
   54 XN=2.*(XJ21-XJ32)+EPS1                                            TRO05600
      IF(XN.LT.0.) XN=-(XN+EPS2)                                        TRO05610
      JMIN=XN                                                           TRO05620
      XN=2.*(XJ11-XJ33)+EPS1                                            TRO05630
      IF(XN.LT.0.) XN=-(XN+EPS2)                                        TRO05640
      N=XN                                                              TRO05650
      IF(N.GT.JMIN) JMIN=N                                              TRO05660
      XN=2.*(XJ12-XJ23)+EPS1                                            TRO05670
      IF(XN.LT.0.) XN=-(XN+EPS2)                                        TRO05680
      N=XN                                                              TRO05690
      IF(N.GT.JMIN) JMIN=N                                              TRO05700
      XN=2.*(XJ21+XJ32)+EPS1                                            TRO05710
      JMAX=XN                                                           TRO05720
      XN=2.*(XJ11+XJ33)+EPS1                                            TRO05730
      N=XN                                                              TRO05740
      IF(N.LT.JMAX) JMAX=N                                              TRO05750
      XN=2.*(XJ12+XJ23)+EPS1                                            TRO05760
      N=XN                                                              TRO05770
      IF(N.LT.JMAX) JMAX=N                                              TRO05780
      XJMIN=DFLOAT(JMIN)/2.                                             TRO05790
      XJMAX=DFLOAT(JMAX)/2.                                             TRO05800
      S=0.                                                              TRO05810
      XJ=XJMIN                                                          TRO05820
      XN=-XJ32+XJ21+XJ+EPS1                                             TRO05830
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05840
      N3=XN                                                             TRO05850
      XN=-XJ11+XJ33+XJ+EPS1                                             TRO05860
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05870
      N4=XN                                                             TRO05880
      XN=XJ32-XJ21+XJ+EPS1                                              TRO05890
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05900
      N7=XN                                                             TRO05910
      XN=XJ11-XJ33+XJ+EPS1                                              TRO05920
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05930
      N8=XN                                                             TRO05940
      XN=XJ32+XJ21-XJ+EPS1                                              TRO05950
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05960
      N11=XN                                                            TRO05970
      XN=XJ11+XJ33-XJ+EPS1                                              TRO05980
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO05990
      N12=XN                                                            TRO06000
      XN=-XJ11-XJ32+XJ31+XJ+EPS1                                        TRO06010
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06020
      N13=XN                                                            TRO06030
      XN=-XJ21-XJ33+XJ31+XJ+EPS1                                        TRO06040
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06050
      N14=XN                                                            TRO06060
      XN=XJ32+XJ21+XJ+EPS1                                              TRO06070
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06080
      N18=XN                                                            TRO06090
      N18=N18+1                                                         TRO06100
      XN=XJ11+XJ33+XJ+EPS1                                              TRO06110
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06120
      N19=XN                                                            TRO06130
      N19=N19+1                                                         TRO06140
      XN=-XJ21+XJ+XJ32+EPS1                                             TRO06150
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06160
      N22=XN                                                            TRO06170
      XN=-XJ12+XJ23+XJ+EPS1                                             TRO06180
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06190
      N24=XN                                                            TRO06200
      XN=XJ21-XJ+XJ32+EPS1                                              TRO06210
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06220
      N26=XN                                                            TRO06230
      XN=XJ12-XJ+XJ23+EPS1                                              TRO06240
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06250
      N28=XN                                                            TRO06260
      XN=XJ21+XJ-XJ32+EPS1                                              TRO06270
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06280
      N30=XN                                                            TRO06290
      XN=XJ12+XJ-XJ23+EPS1                                              TRO06300
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06310
      N32=XN                                                            TRO06320
      XN=-XJ22-XJ+XJ32+XJ23+EPS1                                        TRO06330
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06340
      N34=XN                                                            TRO06350
      XN=XJ12+XJ21+XJ22+XJ+EPS1                                         TRO06360
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06370
      N35=XN                                                            TRO06380
      N35=N35+1                                                         TRO06390
      XN=XJ21+XJ+XJ32+EPS1                                              TRO06400
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06410
      N37=XN                                                            TRO06420
      N37=N37+1                                                         TRO06430
      XN=XJ12+XJ+XJ23+EPS1                                              TRO06440
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06450
      N39=XN                                                            TRO06460
      N39=N39+1                                                         TRO06470
      XN=-XJ+XJ11+XJ33+EPS1                                             TRO06480
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06490
      N42=XN                                                            TRO06500
      XN=-XJ+XJ23+XJ12+EPS1                                             TRO06510
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06520
      N43=XN                                                            TRO06530
      XN=XJ-XJ11+XJ33+EPS1                                              TRO06540
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06550
      N46=XN                                                            TRO06560
      XN=XJ-XJ23+XJ12+EPS1                                              TRO06570
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06580
      N47=XN                                                            TRO06590
      XN=XJ+XJ11-XJ33+EPS1                                              TRO06600
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06610
      N50=XN                                                            TRO06620
      XN=XJ+XJ23-XJ12+EPS1                                              TRO06630
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06640
      N51=XN                                                            TRO06650
      XN=-XJ13-XJ+XJ33+XJ12+EPS1                                        TRO06660
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06670
      N53=XN                                                            TRO06680
      XN=XJ13+XJ+XJ23+XJ11+EPS1                                         TRO06690
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06700
      N55=XN                                                            TRO06710
      N55=N55+1                                                         TRO06720
      XN=XJ+XJ11+XJ33+EPS1                                              TRO06730
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06740
      N57=XN                                                            TRO06750
      N57=N57+1                                                         TRO06760
      XN=XJ+XJ23+XJ12+EPS1                                              TRO06770
      IF(XN.LT.0.) XN=XN+EPS2                                           TRO06780
      N58=XN                                                            TRO06790
      N58=N58+1                                                         TRO06800
      GO TO 10                                                          TRO06810
   52 IF(XJ.GT.XJMAX) GO TO 120                                         TRO06820
      N3=N3+1                                                           TRO06830
      N4=N4+1                                                           TRO06840
      N7=N7+1                                                           TRO06850
      N8=N8+1                                                           TRO06860
      N11=N11-1                                                         TRO06870
      N12=N12-1                                                         TRO06880
      N13=N13+1                                                         TRO06890
      N14=N14+1                                                         TRO06900
      N18=N18+1                                                         TRO06910
      N19=N19+1                                                         TRO06920
      N22=N22+1                                                         TRO06930
      N24=N24+1                                                         TRO06940
      N26=N26-1                                                         TRO06950
      N28=N28-1                                                         TRO06960
      N30=N30+1                                                         TRO06970
      N32=N32+1                                                         TRO06980
      N34=N34-1                                                         TRO06990
      N35=N35+1                                                         TRO07000
      N37=N37+1                                                         TRO07010
      N39=N39+1                                                         TRO07020
      N42=N42-1                                                         TRO07030
      N43=N43-1                                                         TRO07040
      N46=N46+1                                                         TRO07050
      N47=N47+1                                                         TRO07060
      N50=N50+1                                                         TRO07070
      N51=N51+1                                                         TRO07080
      N53=N53-1                                                         TRO07090
      N55=N55+1                                                         TRO07100
      N57=N57+1                                                         TRO07110
      N58=N58+1                                                         TRO07120
   10 K1=0.                                                             TRO07160
      L1=-N13                                                           TRO07170
      IF(L1.GT.K1) K1=L1                                                TRO07180
      L1=-N14                                                           TRO07190
      IF(L1.GT.K1) K1=L1                                                TRO07200
      L1=N9                                                             TRO07210
      IF(N10.LT.L1) L1=N10                                              TRO07220
      IF(N11.LT.L1) L1=N11                                              TRO07230
      IF(N12.LT.L1) L1=N12                                              TRO07240
      IF(N15.LT.L1) L1=N15                                              TRO07250
      F1=1.                                                             TRO07260
      S1=1.                                                             TRO07270
      I1=K1+1                                                           TRO07280
  62  IF(I1.GT.L1) GO TO 64                                             TRO07290
      I1M1=I1-1                                                         TRO07300
      NN1=(N9-I1M1)*(N10-I1M1)*(N11-I1M1)*(N12-I1M1)                    TRO07310
      ND1=I1*(N13+I1)*(N14+I1)*(N15-I1M1)                               TRO07320
      F1=-F1*DFLOAT(NN1)/DFLOAT(ND1)                                    TRO07330
      S1=S1+F1                                                          TRO07340
      I1=I1+1                                                           TRO07350
      GO TO 62                                                          TRO07360
   64 K2=0                                                              TRO07370
      L2=-N33                                                           TRO07380
      IF(L2.GT.K2) K2=L2                                                TRO07390
      L2=-N34                                                           TRO07400
      IF(L2.GT.K2) K2=L2                                                TRO07410
      L2=N29                                                            TRO07420
      IF(N30.LT.L2) L2=N30                                              TRO07430
      IF(N31.LT.L2) L2=N31                                              TRO07440
      IF(N32.LT.L2) L2=N32                                              TRO07450
      IF(N35.LT.L2) L2=N35                                              TRO07460
      F2=1.                                                             TRO07470
      S2=1.                                                             TRO07480
      I2=K2+1                                                           TRO07490
  70  IF(I2.GT.L2) GO TO 80                                             TRO07500
      I2M2=I2-1                                                         TRO07510
      NN2=(N29-I2M2)*(N30-I2M2)*(N31-I2M2)*(N32-I2M2)                   TRO07520
      ND2=I2*(N33+I2)*(N34+I2)*(N35-I2M2)                               TRO07530
      F2=-F2*DFLOAT(NN2)/DFLOAT(ND2)                                    TRO07540
      S2=S2+F2                                                          TRO07550
      I2=I2+1                                                           TRO07560
      GO TO 70                                                          TRO07570
   80 K3=0                                                              TRO07580
      L3=-N53                                                           TRO07590
      IF(L3.GT.K3) K3=L3                                                TRO07600
      L3=-N54                                                           TRO07610
      IF(L3.GT.K3) K3=L3                                                TRO07620
      L3=N49                                                            TRO07630
      IF(N50.LT.L3) L3=N50                                              TRO07640
      IF(N51.LT.L3) L3=N51                                              TRO07650
      IF(N52.LT.L3) L3=N52                                              TRO07660
      IF(N55.LT.L3) L3=N55                                              TRO07670
      F3=1.                                                             TRO07680
      S3=1.                                                             TRO07690
      I3=K3+1                                                           TRO07700
   84 IF(I3.GT.L3) GO TO 90                                             TRO07710
      I3M3=I3-1                                                         TRO07720
      NN3=(N49-I3M3)*(N50-I3M3)*(N51-I3M3)*(N52-I3M3)                   TRO07730
      ND3=I3*(N53+I3)*(N54+I3)*(N55-I3M3)                               TRO07740
      F3=-F3*DFLOAT(NN3)/DFLOAT(ND3)                                    TRO07750
      S3=S3+F3                                                          TRO07760
      I3=I3+1                                                           TRO07770
      GO TO 84                                                          TRO07780
   90 S2N=FLOG(N3+1)+FLOG(N4+1)+FLOG(N7+1)+FLOG(N8+1)+FLOG(N11+1)+FLOG(NTRO07820
     112+1)+FLOG(N22+1)+FLOG(N24+1)+FLOG(N26+1)+FLOG(N28+1)+FLOG(N30+1)+TRO07830
     2FLOG(N32+1)+FLOG(N42+1)+FLOG(N43+1)+FLOG(N46+1)+FLOG(N47+1)+FLOG(NTRO07840
     350+1)+FLOG(N51+1)                                                 TRO07850
      S2N=.5*S2N                                                        TRO07860
      S2D=FLOG(N18+1)+FLOG(N19+1)+FLOG(N37+1)+FLOG(N39+1)+FLOG(N57+1)+FLTRO07870
     1OG(N58+1)                                                         TRO07880
      S2D=.5*S2D                                                        TRO07890
      KM1=K1-1                                                          TRO07900
      KP1=K1+1                                                          TRO07910
      KM2=K2-1                                                          TRO07920
      KP2=K2+1                                                          TRO07930
      KM3=K3-1                                                          TRO07940
      KP3=K3+1                                                          TRO07950
      S2N=S2N+FLOG(N15-KM1)+FLOG(N35-KM2)+FLOG(N55-KM3)                 TRO07960
      S2D=S2D+FLOG(KP1)+FLOG(KP2)+FLOG(KP3)+FLOG(N9-KM1)+FLOG(N10-KM1)+FTRO07970
     1LOG(N11-KM1)+FLOG(N12-KM1)+FLOG(N13+KP1)+FLOG(N14+KP1)+FLOG(N29-KMTRO07980
     22)+FLOG(N30-KM2)+FLOG(N31-KM2)+FLOG(N32-KM2)+FLOG(N33+KP2)+FLOG(N3TRO07990
     34+KP2)+FLOG(N49-KM3)+FLOG(N50-KM3)+FLOG(N51-KM3)                  TRO08000
     4+FLOG(N52-KM3)+FLOG(N53+KP3)+FLOG(N54+KP3)                        TRO08010
      F=S2D-S2N                                                         TRO08020
      IF (F.GT.80) GO TO 100                                            TRO08030
      F=S2N/S2D                                                         TRO08040
      IF((F.LT.1.01).AND.(F.GT.0.98)) GO TO 100                         TRO08050
      F=S1*S2*S3*DEXP(S2N-S2D)                                          TRO08060
      GO TO 110                                                         TRO08070
  100 F=S1*S2*S3                                                        TRO08080
      IF(F) 102,112,104                                                 TRO08090
  102 F=DLOG(-F)                                                        TRO08100
      F=-DEXP(F+S2N-S2D)                                                TRO08110
      GO TO 110                                                         TRO08120
  104 F=DLOG(F)                                                         TRO08130
      F=DEXP(F+S2N-S2D)                                                 TRO08140
  110 F=F*(2.*XJ+1.)                                                    TRO08150
      L1=K1+K2+K3                                                       TRO08190
      K1=L1/2                                                           TRO08200
      K1=2*K1                                                           TRO08210
      IF(L1.NE.K1) F=-F                                                 TRO08220
      S=S+F                                                             TRO08230
  112 XJ=XJ+1.                                                          TRO08240
      GO TO 52                                                          TRO08250
  120 C2N=FLOG(N1+1)+FLOG(N2+1)+FLOG(N5+1)+FLOG(N6+1)+FLOG(N9+1)+FLOG(N1TRO08290
     10+1)+FLOG(N21+1)+FLOG(N23+1)+FLOG(N25+1)+FLOG(N27+1)+FLOG(N29+1)+FTRO08300
     2LOG(N31+1)+FLOG(N41+1)+FLOG(N44+1)+FLOG(N45+1)+FLOG(N48+1)+FLOG(N4TRO08310
     39+1)+FLOG(N52+1)                                                  TRO08320
      C2N=.5*C2N                                                        TRO08330
      C2D=FLOG(N16+1)+FLOG(N17+1)+FLOG(N36+1)+FLOG(N38+1)+FLOG(N56+1)+FLTRO08340
     1OG(N59+1)                                                         TRO08350
      C2D=.5*C2D                                                        TRO08360
      F=C2D-C2N                                                         TRO08370
      IF(F.GT.80.) GO TO 122                                            TRO08380
      F=C2N/C2D                                                         TRO08390
      IF((F.LT.1.01).AND.(F.GT.0.98)) GO TO 122                         TRO08400
      C9J=S*DEXP(C2N-C2D)                                               TRO08410
      GO TO 130                                                         TRO08420
  122 IF(S) 124,50,126                                                  TRO08430
  124 S=DLOG(-S)                                                        TRO08440
      C9J=-DEXP(S+C2N-C2D)                                              TRO08450
      GO TO 130                                                         TRO08460
  126 S=DLOG(S)                                                         TRO08470
      C9J=DEXP(S+C2N-C2D)                                               TRO08480
  130 K=N9+N16+N36+N56-1                                                TRO08520
      L=K/2                                                             TRO08530
      L=2*L                                                             TRO08540
      IF (L.NE.K) C9J=-C9J                                              TRO08550

      RETURN                                                            TRO08560
      END                                                               TRO08570


      subroutine slatcould(ka,kb,kc,kd,orba,orbb,orbc,orbd,
     :iaa,iab,iac,iad,k1,rkll,rkss)
      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /basis/gl(MN,MNBAS),gs(MN,MNBAS)
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)

      do i=1,MN
         tmp1(i)=zero
         tmp2(i)=zero
         tmp3(i)=zero
         tmp4(i)=zero
      enddo

      rkll=zero
      rkss=zero
      fact=zero

      fact=selectc(iaa,iac,orba,orbc,k1)*
     : selectc(iab,iad,orbb,orbd,k1)

      if(fact.ne.zero)then

      call yzk(k1,kb,kd)

      do ir = 2,n
         tmp2(ir)=gl(ir,ka)*gl(ir,kc)*rpor(ir)*tmp1(ir)
      enddo
      call valint(tmp2,rkll)

      do ir = 2,n
         tmp3(ir)=gs(ir,ka)*gs(ir,kc)*rpor(ir)*tmp1(ir)
      enddo

      call valint(tmp3,rkss)

      endif

      return
      end


      SUBROUTINE YZK(K,I,J)
*CCC
*CCC      $Id: YZK.f,v 1.4 1994/09/15 13:34:04 pn Exp $
*CCC
*CCC      $Log: YZK.f,v $
*CCC      Revision 1.4  1994/09/15 13:34:04  pn
*CCC      Version 2c  31 May 1994
*CCC      Modified version of 2b.
*CCC
*CCC Revision 1.3  1994/03/31  11:57:06  pn
*CCC Version 2b 6 July 1990
*CCC
*CCC
C***********************************************************************
C***********************************************************************
C-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)

      DIMENSION RHOP(MN),RTTK(MN),WK(MN),ZK(MN)

C-----------------------------------------------------------------------

C Coefficients for 10 point integration formula as used in VALINT.
      D1=0.26834148361926139703d0*h
      D2=7.13764630431297097963d0*h
      D3=1.77535941424830313719d0*h
      D4=-.81043570626903960237d0*h
      D5=4.54946288279621612954d0*h
      D6=-4.35155122655122655122d0*h

      if((n+10).ge.MN)then
      print*,'Dimension of grid size not sufficient'
      stop
      endif

C Initializing all the variables
      Do II=1,MN
      WK(II) = ZERO
      ZK(II) = ZERO
      TMP1(II) = ZERO
      TMP4(II) = ZERO
      RHOP(II) = ZERO
      ENDDO

C Calculates r(i)**k, where K is the multipole
        DO 1 II = 2,N
          RTTK(II) = R(II)**K
   1    CONTINUE

C Calculates Integrand for Z(i)**k at different grids
      DO 2 II = 2,N
        RHOP(II) = RP(II)*(PF(II,I)*PF(II,J)+QF(II,I)*QF(II,J))
   2  CONTINUE

C Defines temporray functions as r(i)**k*Integrand
        DO 3  II = 2,N
          TMP4(II) = RTTK(II)*RHOP(II)
   3    CONTINUE

C First 10 points of Z(i)**k are taken same as given below
      DO 4  II = 2,10

      ZK(II) = D1*(TMP4(1)+TMP4(11))+D2*TMP4(6)
     :  +D3*(TMP4(2)+TMP4(10))+D4*(TMP4(3)+TMP4(9))
     :  +D5*(TMP4(4)+TMP4(8))+D6*(TMP4(5)+TMP4(7))

   4  CONTINUE

C Z(i)**k has been calculated using 10 integration points.
      DO 5  II = 11,N+10

      ZK(II) = ZK(II-10)+D1*(TMP4(II-10)+TMP4(II))+D2*TMP4(II-5)
     :  +D3*(TMP4(II-1)+TMP4(II-9))+D4*(TMP4(II-2)+TMP4(II-8))
     :  +D5*(TMP4(II-3)+TMP4(II-7))+D6*(TMP4(II-4)+TMP4(II-6))

   5  CONTINUE

C Z(i)**k has been defined as Z(i)**k*(r</r>)**k
        DO 6 II = 2,N
          ZK(II) = ZK(II)/RTTK(II)
   6    CONTINUE

C R(i)**(K+1) has been defined
        DO 7 II = 2,N
          RTTK(II) = RTTK(II)*R(II)
   7    CONTINUE

C Calculates Integrand for W(i)**k/r(i)**(k+1) at different grids
        DO 8 II = 2,N
          TMP4(II) = RHOP(II)/RTTK(II)
   8    CONTINUE

C Last 10 points of W(i)**k are taken same as given below
      DO 9 II = N,N+10

      WK(II) = D1*(TMP4(N)+TMP4(N+10))+D2*TMP4(N+5)
     :  +D3*(TMP4(N+1)+TMP4(N+9))+D4*(TMP4(N+2)+TMP4(N+8))
     :  +D5*(TMP4(N+3)+TMP4(N+7))+D6*(TMP4(N+4)+TMP4(N+6))

   9  CONTINUE

C W(i)**k has been calculated using 10 integration points.
      DO 10 II = N,2,-1

      WK(II) = WK(II+10)+D1*(TMP4(II)+TMP4(II+10))+D2*TMP4(II+5)
     :  +D3*(TMP4(II+1)+TMP4(II+9))+D4*(TMP4(II+2)+TMP4(II+8))
     :  +D5*(TMP4(II+3)+TMP4(II+7))+D6*(TMP4(II+4)+TMP4(II+6))

   10 CONTINUE

C W(i)**k has been defined as Z(i)**k*(r</r>)**(k+1)
        DO 11 II = 2,N
          WK(II) = WK(II)*RTTK(II)
   11   CONTINUE

C YZK functions are determined as Z(i)**k+W(i)**k
      DO 12 II = 2,N
        TMP1(II) = ZK(II) + WK(II)
   12 CONTINUE


      RETURN
      END


      function selectc(ia,ib,orba,orbb,k1)
      implicit real*8 (a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten

      selectc=zero

      dk1=dfloat(k1)
      da=dfloat(ia)
      db=dfloat(ib)

      selectc=(one-(da*db*(-1)**(orba+orbb+dk1)))/two

      return
      end

      function selectb(ia,ib,orba,orbb,k1)
      implicit real*8 (a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten

      selectb=zero

      dk1=dfloat(k1)
      da=dfloat(ia)
      db=dfloat(ib)

      selectb=(one+(da*db*(-1)**(orba+orbb+dk1)))/two

      return
      end


      subroutine slatcoule(ka,kb,kc,kd,orba,orbb,orbc,orbd,
     :iaa,iab,iac,iad,k1,rkll,rkls,rksl,rkss)
      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /basis/gl(MN,MNBAS),gs(MN,MNBAS)
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)

C New common block
      Common/atmp/tmp5(MN),tmp6(MN)

      do i=1,MN
         tmp1(i)=zero
         tmp2(i)=zero
         tmp3(i)=zero
         tmp4(i)=zero
         tmp5(i)=zero
         tmp6(i)=zero
      enddo

      rkll=zero
      rkls=zero
      rksl=zero
      rkss=zero
      fact=zero

      fact=selectc(iaa,iac,orba,orbc,k1)*
     : selectc(iab,iad,orbb,orbd,k1)

      if(fact.ne.zero)then

      call yzkdfe(k1,kb,kd)

      do ir = 2,n
         tmp2(ir)=gl(ir,ka)*pf(ir,kc)*rpor(ir)*tmp1(ir)
      enddo
      call valint(tmp2,rkll)

      do ir = 2,n
         tmp3(ir)=gl(ir,ka)*pf(ir,kc)*rpor(ir)*tmp4(ir)
      enddo
      call valint(tmp3,rkls)

      do i=1,MN
         tmp2(i)=zero
         tmp3(i)=zero
      enddo

      do ir = 2,n
         tmp2(ir)=gs(ir,ka)*qf(ir,kc)*rpor(ir)*tmp1(ir)
      enddo
      call valint(tmp2,rksl)

      do ir = 2,n
         tmp3(ir)=gs(ir,ka)*qf(ir,kc)*rpor(ir)*tmp4(ir)
      enddo
      call valint(tmp3,rkss)

      endif

      return
      end

      SUBROUTINE YZKDFE(K,I,J)
*CCC
*CCC      $Id: YZK.f,v 1.4 1994/09/15 13:34:04 pn Exp $
*CCC
*CCC      $Log: YZK.f,v $
*CCC      Revision 1.4  1994/09/15 13:34:04  pn
*CCC      Version 2c  31 May 1994
*CCC      Modified version of 2b.
*CCC
*CCC Revision 1.3  1994/03/31  11:57:06  pn
*CCC Version 2b 6 July 1990
*CCC
*CCC
C***********************************************************************
C***********************************************************************
C-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /basis/gl(MN,MNBAS),gs(MN,MNBAS)
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     :      /atmp/tmp5(MN),tmp6(MN)

      DIMENSION RHOP1(MN),RTTK(MN),WK1(MN),ZK1(MN)
      DIMENSION RHOP2(MN),WK2(MN),ZK2(MN)

C-----------------------------------------------------------------------

C Coefficients for 10 point integration formula as used in VALINT.
      D1=.26834148361926139703d0*h
      D2=7.13764630431297097963d0*h
      D3=1.77535941424830313719d0*h
      D4=-.81043570626903960237d0*h
      D5=4.54946288279621612954d0*h
      D6=-4.35155122655122655122d0*h

      if((n+10).ge.MN)then
      print*,'Dimension of grid size not sufficient'
      stop
      endif

C Initializing all the variables
      Do II=1,MN
      WK1(II) = ZERO
      ZK1(II) = ZERO
      WK2(II) = ZERO
      ZK2(II) = ZERO
      TMP1(II) = ZERO
      TMP4(II) = ZERO
      TMP5(II) = ZERO
      TMP6(II) = ZERO
      RHOP1(II) = ZERO
      RHOP2(II) = ZERO 
      ENDDO

C Calculates r(i)**k, where K is the multipole
        DO 1 II = 2,N
          RTTK(II) = R(II)**K
   1    CONTINUE

C Calculates Integrand for Z(i)**k at different grids
      DO 2 II = 2,N
        RHOP1(II) = RP(II)*PF(II,I)*GL(II,J)
        RHOP2(II) = RP(II)*QF(II,I)*GS(II,J)
   2  CONTINUE

C Defines temporray functions as r(i)**k*Integrand
        DO 3  II = 2,N
          TMP5(II) = RTTK(II)*RHOP1(II)
          TMP6(II) = RTTK(II)*RHOP2(II)
   3    CONTINUE

C First 10 points of Z(i)**k are taken same as given below
      DO 4  II = 2,10

      ZK1(II) = D1*(TMP5(1)+TMP5(11))+D2*TMP5(6)
     :  +D3*(TMP5(2)+TMP5(10))+D4*(TMP5(3)+TMP5(9))
     :  +D5*(TMP5(4)+TMP5(8))+D6*(TMP5(5)+TMP5(7))

      ZK2(II) = D1*(TMP6(1)+TMP6(11))+D2*TMP6(6)
     :  +D3*(TMP6(2)+TMP6(10))+D4*(TMP6(3)+TMP6(9))
     :  +D5*(TMP6(4)+TMP6(8))+D6*(TMP6(5)+TMP6(7))

   4  CONTINUE

C Z(i)**k has been calculated using 10 integration points.
      DO 5  II = 11,N+10

      ZK1(II) = ZK1(II-10)+D1*(TMP5(II-10)+TMP5(II))+D2*TMP5(II-5)
     :  +D3*(TMP5(II-1)+TMP5(II-9))+D4*(TMP5(II-2)+TMP5(II-8))
     :  +D5*(TMP5(II-3)+TMP5(II-7))+D6*(TMP5(II-4)+TMP5(II-6))

      ZK2(II) = ZK2(II-10)+D1*(TMP6(II-10)+TMP6(II))+D2*TMP6(II-5)
     :  +D3*(TMP6(II-1)+TMP6(II-9))+D4*(TMP6(II-2)+TMP6(II-8))
     :  +D5*(TMP6(II-3)+TMP6(II-7))+D6*(TMP6(II-4)+TMP6(II-6))

   5  CONTINUE

C Z(i)**k has been defined as Z(i)**k*(r</r>)**k
        DO 6 II = 2,N
          ZK1(II) = ZK1(II)/RTTK(II)
          ZK2(II) = ZK2(II)/RTTK(II)
   6    CONTINUE

C R(i)**(K+1) has been defined
        DO 7 II = 2,N
          RTTK(II) = RTTK(II)*R(II)
   7    CONTINUE

C Calculates Integrand for W(i)**k/r(i)**(k+1) at different grids
        DO 8 II = 2,N
          TMP5(II) = RHOP1(II)/RTTK(II)
          TMP6(II) = RHOP2(II)/RTTK(II)
   8    CONTINUE

C Last 10 points of W(i)**k are taken same as given below
      DO 9 II = N,N+10

      WK1(II) = D1*(TMP5(N)+TMP5(N+10))+D2*TMP5(N+5)
     :  +D3*(TMP5(N+1)+TMP5(N+9))+D4*(TMP5(N+2)+TMP5(N+8))
     :  +D5*(TMP5(N+3)+TMP5(N+7))+D6*(TMP5(N+4)+TMP5(N+6))

      WK2(II) = D1*(TMP6(N)+TMP6(N+10))+D2*TMP6(N+5)
     :  +D3*(TMP6(N+1)+TMP6(N+9))+D4*(TMP6(N+2)+TMP6(N+8))
     :  +D5*(TMP6(N+3)+TMP6(N+7))+D6*(TMP6(N+4)+TMP6(N+6))

   9  CONTINUE

C W(i)**k has been calculated using 10 integration points.
      DO 10 II = N,2,-1

      WK1(II) = WK1(II+10)+D1*(TMP5(II)+TMP5(II+10))+D2*TMP5(II+5)
     :  +D3*(TMP5(II+1)+TMP5(II+9))+D4*(TMP5(II+2)+TMP5(II+8))
     :  +D5*(TMP5(II+3)+TMP5(II+7))+D6*(TMP5(II+4)+TMP5(II+6))

      WK2(II) = WK2(II+10)+D1*(TMP6(II)+TMP6(II+10))+D2*TMP6(II+5)
     :  +D3*(TMP6(II+1)+TMP6(II+9))+D4*(TMP6(II+2)+TMP6(II+8))
     :  +D5*(TMP6(II+3)+TMP6(II+7))+D6*(TMP6(II+4)+TMP6(II+6))

   10 CONTINUE

C W(i)**k has been defined as Z(i)**k*(r</r>)**(k+1)
        DO 11 II = 2,N
          WK1(II) = WK1(II)*RTTK(II)
          WK2(II) = WK2(II)*RTTK(II)
   11   CONTINUE

C YZK functions are determined as Z(i)**k+W(i)**k
      DO 12 II = 2,N
        TMP1(II) = ZK1(II) + WK1(II)
        TMP4(II) = ZK2(II) + WK2(II)
   12 CONTINUE

      RETURN
      END


      subroutine slatbrt(ka,kb,kc,kd,orba,orbb,orbc,orbd,
     :iaa,iab,iac,iad,k1,bkll,bkls,bksl,bkss)
      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /basis/gl(MN,MNBAS),gs(MN,MNBAS)
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /atmp/tmp5(MN),tmp6(MN)
     :  /quan1/kappe(MNBAS),orbje(MNBAS),kappc(MNBAS),orbjc(MNBAS)

      DIMENSION TMP(MN)

      bkll=zero
      bkls=zero
      bksl=zero
      bkss=zero
      fact=zero

      kapa=-idint((orba+half)*iaa)
      kapb=-idint((orbb+half)*iab)
      kapc=-idint((orbc+half)*iac)
      kapd=-idint((orbd+half)*iad)

      if(kapa.ne.kappe(ka).or.kapb.ne.kappe(kb).or.
     :kapc.ne.kappe(kc).or.kapd.ne.kappe(kd))then
      Print*,'Quantum numbers are wrongly assigned'
      Stop
      Endif
      
      If(k1.gt.0)then

      Do L1=k1+1,k1+3

      L=L1-2

      fact=selectb(iaa,iac,orba,orbc,L)*
     : selectb(iab,iad,orbb,orbd,L)

      if(fact.gt.zero)then

      do i=1,MN
         tmp(i)=zero
         tmp1(i)=zero
         tmp2(i)=zero
         tmp3(i)=zero
         tmp4(i)=zero
         tmp5(i)=zero
         tmp6(i)=zero
      enddo

         res=zero

      call yzkbrt(L,kb,kd)

      do ie = 2,n
         tmp(ie)=gl(ie,ka)*qf(ie,kc)*rpor(ie)*tmp1(ie)
      enddo
      call valint(tmp,res)

      bkll=bkll+res*GK(2,k1,L,kapa,kapb,kapc,kapd)

      do ie = 2,n
         tmp(ie)=gl(ie,ka)*qf(ie,kc)*rpor(ie)*tmp2(ie)
      enddo
      call valint(tmp,res)

      bkls=bkls+res*GK(1,k1,L,kapa,kapb,kapc,kapd)

      do ie = 2,n
         tmp(ie)=gs(ie,ka)*pf(ie,kc)*rpor(ie)*tmp1(ie)
      enddo
      call valint(tmp,res)

      bksl=bksl+res*GK(4,k1,L,kapa,kapb,kapc,kapd)

      do ie = 2,n
         tmp(ie)=gs(ie,ka)*pf(ie,kc)*rpor(ie)*tmp2(ie)
      enddo
      call valint(tmp,res)

      bkss=bkss+res*GK(3,k1,L,kapa,kapb,kapc,kapd)

      If(ibrt.gt.1)then

      do ie = 2,n
         tmp(ie)=gl(ie,ka)*qf(ie,kc)*rpor(ie)*tmp3(ie)
      enddo
      call valint(tmp,res)

      bkll=bkll+res*RA(2,k1,L,kapa,kapb,kapc,kapd)

      do ie = 2,n
         tmp(ie)=gl(ie,ka)*qf(ie,kc)*rpor(ie)*tmp4(ie)
      enddo
      call valint(tmp,res)

      bkls=bkls+res*RA(1,k1,L,kapa,kapb,kapc,kapd)

      do ie = 2,n
         tmp(ie)=gs(ie,ka)*pf(ie,kc)*rpor(ie)*tmp3(ie)
      enddo
      call valint(tmp,res)

      bksl=bksl+res*RA(4,k1,L,kapa,kapb,kapc,kapd)

      do ie = 2,n
         tmp(ie)=gs(ie,ka)*pf(ie,kc)*rpor(ie)*tmp4(ie)
      enddo
      call valint(tmp,res)

      bkss=bkss+res*RA(3,k1,L,kapa,kapb,kapc,kapd)

      do ie = 2,n
         tmp(ie)=gl(ie,ka)*qf(ie,kc)*rpor(ie)*tmp5(ie)
      enddo
      call valint(tmp,res)

      bkll=bkll+res*RB(2,k1,L,kapa,kapb,kapc,kapd)

      do ie = 2,n
         tmp(ie)=gl(ie,ka)*qf(ie,kc)*rpor(ie)*tmp6(ie)
      enddo
      call valint(tmp,res)

      bkls=bkls+res*RB(1,k1,L,kapa,kapb,kapc,kapd)

      do ie = 2,n
         tmp(ie)=gs(ie,ka)*pf(ie,kc)*rpor(ie)*tmp5(ie)
      enddo
      call valint(tmp,res)

      bksl=bksl+res*RB(4,k1,L,kapa,kapb,kapc,kapd)

      do ie = 2,n
         tmp(ie)=gs(ie,ka)*pf(ie,kc)*rpor(ie)*tmp6(ie)
      enddo
      call valint(tmp,res)

      bkss=bkss+res*RB(3,k1,L,kapa,kapb,kapc,kapd)

      EndIf

      endif

      EndDo

      EndIf

      return
      end


      SUBROUTINE YZKBRT(K,I,J)
*CCC
*CCC      $Id: YZK.f,v 1.4 1994/09/15 13:34:04 pn Exp $
*CCC
*CCC      $Log: YZK.f,v $
*CCC      Revision 1.4  1994/09/15 13:34:04  pn
*CCC      Version 2c  31 May 1994
*CCC      Modified version of 2b.
*CCC
*CCC Revision 1.3  1994/03/31  11:57:06  pn
*CCC Version 2b 6 July 1990
*CCC
*CCC
C***********************************************************************
C***********************************************************************
C-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /basis/gl(MN,MNBAS),gs(MN,MNBAS)
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     :      /atmp/tmp5(MN),tmp6(MN)

      DIMENSION RHOP1(MN),RTTK(MN),WK1(MN),ZK1(MN)
      DIMENSION RHOP2(MN),WK2(MN),ZK2(MN)

C-----------------------------------------------------------------------

C Coefficients for 10 point integration formula as used in VALINT.
      D1=.26834148361926139703d0*h
      D2=7.13764630431297097963d0*h
      D3=1.77535941424830313719d0*h
      D4=-.81043570626903960237d0*h
      D5=4.54946288279621612954d0*h
      D6=-4.35155122655122655122d0*h

      if((n+10).ge.MN)then
      print*,'Dimension of grid size not sufficient'
      stop
      endif

C Initializing all the variables
      Do II=1,MN
      WK1(II) = ZERO
      ZK1(II) = ZERO
      WK2(II) = ZERO
      ZK2(II) = ZERO
      TMP1(II) = ZERO
      TMP2(II) = ZERO
      TMP3(II) = ZERO
      TMP4(II) = ZERO
      TMP5(II) = ZERO
      TMP6(II) = ZERO
      RHOP1(II) = ZERO
      RHOP2(II) = ZERO
      ENDDO

C Calculates r(i)**k, where K is the multipole
        DO 1 II = 2,N
          RTTK(II) = R(II)**K
   1    CONTINUE

C Calculates Integrand for Z(i)**k at different grids
      DO 2 II = 2,N
        RHOP1(II) = RP(II)*QF(II,I)*GL(II,J)
        RHOP2(II) = RP(II)*PF(II,I)*GS(II,J)
   2  CONTINUE

C Defines temporray functions as r(i)**k*Integrand
        DO 3  II = 2,N
          TMP5(II) = RTTK(II)*RHOP1(II)
          TMP6(II) = RTTK(II)*RHOP2(II)
   3    CONTINUE

C First 10 points of Z(i)**k are taken same as given below
      DO 4  II = 2,10

      ZK1(II) = D1*(TMP5(1)+TMP5(11))+D2*TMP5(6)
     :  +D3*(TMP5(2)+TMP5(10))+D4*(TMP5(3)+TMP5(9))
     :  +D5*(TMP5(4)+TMP5(8))+D6*(TMP5(5)+TMP5(7))

      ZK2(II) = D1*(TMP6(1)+TMP6(11))+D2*TMP6(6)
     :  +D3*(TMP6(2)+TMP6(10))+D4*(TMP6(3)+TMP6(9))
     :  +D5*(TMP6(4)+TMP6(8))+D6*(TMP6(5)+TMP6(7))

   4  CONTINUE

C Z(i)**k has been calculated using 10 integration points.
      DO 5  II = 11,N+10

      ZK1(II) = ZK1(II-10)+D1*(TMP5(II-10)+TMP5(II))+D2*TMP5(II-5)
     :  +D3*(TMP5(II-1)+TMP5(II-9))+D4*(TMP5(II-2)+TMP5(II-8))
     :  +D5*(TMP5(II-3)+TMP5(II-7))+D6*(TMP5(II-4)+TMP5(II-6))

      ZK2(II) = ZK2(II-10)+D1*(TMP6(II-10)+TMP6(II))+D2*TMP6(II-5)
     :  +D3*(TMP6(II-1)+TMP6(II-9))+D4*(TMP6(II-2)+TMP6(II-8))
     :  +D5*(TMP6(II-3)+TMP6(II-7))+D6*(TMP6(II-4)+TMP6(II-6))

   5  CONTINUE

C Z(i)**k has been defined as Z(i)**k*(r</r>)**k
        DO 6 II = 2,N
          ZK1(II) = ZK1(II)/RTTK(II)
          ZK2(II) = ZK2(II)/RTTK(II)
   6    CONTINUE

C R(i)**(K+1) has been defined
        DO 7 II = 2,N
          RTTK(II) = RTTK(II)*R(II)
   7    CONTINUE

C Calculates Integrand for W(i)**k/r(i)**(k+1) at different grids
        DO 8 II = 2,N
          TMP5(II) = RHOP1(II)/RTTK(II)
          TMP6(II) = RHOP2(II)/RTTK(II)
   8    CONTINUE

C Last 10 points of W(i)**k are taken same as given below
      DO 9 II = N,N+10

      WK1(II) = D1*(TMP5(N)+TMP5(N+10))+D2*TMP5(N+5)
     :  +D3*(TMP5(N+1)+TMP5(N+9))+D4*(TMP5(N+2)+TMP5(N+8))
     :  +D5*(TMP5(N+3)+TMP5(N+7))+D6*(TMP5(N+4)+TMP5(N+6))

      WK2(II) = D1*(TMP6(N)+TMP6(N+10))+D2*TMP6(N+5)
     :  +D3*(TMP6(N+1)+TMP6(N+9))+D4*(TMP6(N+2)+TMP6(N+8))
     :  +D5*(TMP6(N+3)+TMP6(N+7))+D6*(TMP6(N+4)+TMP6(N+6))

   9  CONTINUE

C W(i)**k has been calculated using 10 integration points.
      DO 10 II = N,2,-1

      WK1(II) = WK1(II+10)+D1*(TMP5(II)+TMP5(II+10))+D2*TMP5(II+5)
     :  +D3*(TMP5(II+1)+TMP5(II+9))+D4*(TMP5(II+2)+TMP5(II+8))
     :  +D5*(TMP5(II+3)+TMP5(II+7))+D6*(TMP5(II+4)+TMP5(II+6))

      WK2(II) = WK2(II+10)+D1*(TMP6(II)+TMP6(II+10))+D2*TMP6(II+5)
     :  +D3*(TMP6(II+1)+TMP6(II+9))+D4*(TMP6(II+2)+TMP6(II+8))
     :  +D5*(TMP6(II+3)+TMP6(II+7))+D6*(TMP6(II+4)+TMP6(II+6))

   10 CONTINUE

C W(i)**k has been defined as Z(i)**k*(r</r>)**(k+1)
        DO 11 II = 2,N
          WK1(II) = WK1(II)*RTTK(II)
          WK2(II) = WK2(II)*RTTK(II)
   11   CONTINUE

C YZK functions are determined as Z(i)**k+W(i)**k
      DO 12 II = 2,N
        TMP1(II) = ZK1(II) + WK1(II)
        TMP2(II) = ZK2(II) + WK2(II)
        TMP3(II) = ZK1(II)
        TMP4(II) = ZK2(II)
        TMP5(II) = WK1(II)
        TMP6(II) = WK2(II)
   12 CONTINUE

      RETURN
      END


      FUNCTION GK(I,K,L,KAPA,KAPB,KAPC,KAPD)
      implicit real*8 (a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg

      GK=ZERO
      AK=DFLOAT(K)

      KAP1=KAPC-KAPA
      KAP2=KAPD-KAPB
      KAP3=KAPA+KAPC
      KAP4=KAPB+KAPD

      AKAP1=DFLOAT(KAP1)
      AKAP2=DFLOAT(KAP2)
      AKAP3=DFLOAT(KAP3)
      AKAP4=DFLOAT(KAP4)

      IF(L.EQ.K-1)THEN

      GK=ZERO
      FACT=ZERO

      if(ibrt.eq.1)then
      FACT=ONE/(AK*(TWO*AK-ONE))
      else
      FACT=(AK+ONE)/(AK*(TWO*AK-ONE)*(TWO*AK+ONE))
      endif

      IF(I.EQ.1)GK=FACT*(AKAP1+AK)*(AKAP2+AK)
      IF(I.EQ.2)GK=FACT*(AKAP1+AK)*(AKAP2-AK)
      IF(I.EQ.3)GK=FACT*(AKAP1-AK)*(AKAP2+AK)
      IF(I.EQ.4)GK=FACT*(AKAP1-AK)*(AKAP2-AK)

      ENDIF

      IF(L.EQ.K)THEN

      GK=ZERO
      GK=-AKAP3*AKAP4/(AK*(AK+ONE))

      ENDIF

      IF(L.EQ.K+1)THEN

      GK=ZERO
      FACT=ZERO

      if(ibrt.eq.1)then
      FACT=ONE/((AK+ONE)*(TWO*AK+THREE))
      else
      FACT=AK/((AK+ONE)*(TWO*AK+ONE)*(TWO*AK+THREE))
      endif

      IF(I.EQ.1)GK=FACT*(AKAP1-AK-ONE)*(AKAP2-AK-ONE)
      IF(I.EQ.2)GK=FACT*(AKAP1-AK-ONE)*(AKAP2+AK+ONE)
      IF(I.EQ.3)GK=FACT*(AKAP1+AK+ONE)*(AKAP2-AK-ONE)
      IF(I.EQ.4)GK=FACT*(AKAP1+AK+ONE)*(AKAP2+AK+ONE)

      ENDIF

      RETURN
      END


      FUNCTION RA(I,K,L,KAPA,KAPB,KAPC,KAPD)
      implicit real*8 (a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten

      RA=ZERO
      AK=DFLOAT(K)

      KAP1=KAPC-KAPA
      KAP2=KAPD-KAPB

      AKAP1=DFLOAT(KAP1)
      AKAP2=DFLOAT(KAP2)

      FACT=ZERO

      IF(L.EQ.K-1)FACT=-ONE/(TWO*(TWO*AK+ONE))
      IF(L.EQ.K+1)FACT=ONE/(TWO*(TWO*AK+ONE))

      IF(I.EQ.1)RA=FACT*(AKAP1-AK-ONE)*(AKAP2+AK)
      IF(I.EQ.2)RA=FACT*(AKAP1-AK-ONE)*(AKAP2-AK)
      IF(I.EQ.3)RA=FACT*(AKAP1+AK+ONE)*(AKAP2+AK)
      IF(I.EQ.4)RA=FACT*(AKAP1+AK+ONE)*(AKAP2-AK)

      RETURN
      END


      FUNCTION RB(I,K,L,KAPA,KAPB,KAPC,KAPD)
      implicit real*8 (a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten

      RB=ZERO
      AK=DFLOAT(K)

      KAP1=KAPC-KAPA
      KAP2=KAPD-KAPB

      AKAP1=DFLOAT(KAP1)
      AKAP2=DFLOAT(KAP2)

      FACT=ZERO

      IF(L.EQ.K-1)FACT=-ONE/(TWO*(TWO*AK+ONE))
      IF(L.EQ.K+1)FACT=ONE/(TWO*(TWO*AK+ONE))

      IF(I.EQ.1)RB=FACT*(AKAP1+AK)*(AKAP2-AK-ONE)
      IF(I.EQ.2)RB=FACT*(AKAP1+AK)*(AKAP2+AK+ONE)
      IF(I.EQ.3)RB=FACT*(AKAP1-AK)*(AKAP2-AK-ONE)
      IF(I.EQ.4)RB=FACT*(AKAP1-AK)*(AKAP2+AK+ONE)

      RETURN
      END

      subroutine slatistp(ka,kb,kc,kd,orba,orbb,orbc,orbd,
     :iaa,iab,iac,iad,laa,lab,lba,lbb,lca,lcb,lda,ldb,
     :k1,skll,skls,sksl,skss,cstp)
      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /basis/gl(MN,MNBAS),gs(MN,MNBAS)
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     : /propopt/istp,ialpvr
     :  /quan1/kappe(MNBAS),orbje(MNBAS),kappc(MNBAS),orbjc(MNBAS)
     :  /stpfac/dg1(MN,MNBAS),dg2(MN,MNBAS)
     :  /atmp/tmp5(MN),tmp6(MN)
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c

      do i=1,MN
         tmp1(i)=zero
         tmp2(i)=zero
         tmp3(i)=zero
         tmp4(i)=zero
         tmp5(i)=zero
         tmp6(i)=zero
      enddo

      skll=zero
      skls=zero
      sksl=zero
      skss=zero
      fact=zero

      kapa=-idint((orba+half)*iaa)
      kapb=-idint((orbb+half)*iab)
      kapc=-idint((orbc+half)*iac)
      kapd=-idint((orbd+half)*iad)

      if(kapa.ne.kappe(ka).or.kapb.ne.kappe(kb).or.
     :kapc.ne.kappe(kc).or.kapd.ne.kappe(kd))then
      Print*,'Quantum numbers are wrongly assigned'
      Stop
      Endif

      if(k1.ne.1)then
      print*,'Wrong in k1 value as 1',k1
      stop
      endif

      fact=selectc(iaa,iac,orba,orbc,k1)*
     : selectc(iab,iad,orbb,orbd,k1)

      if(fact.ne.zero)then

      if(istp.eq.1)then

      do ie = 2,n
         tmp1(ie)=gl(ie,ka)*qf(ie,kc)*rp(ie)
      enddo
      call valint(tmp1,res1)

      do ie = 2,n
         tmp2(ie)=gs(ie,ka)*pf(ie,kc)*rp(ie)
      enddo
      call valint(tmp2,res2)

      do ie = 2,n
         tmp3(ie)=qf(ie,kb)*gl(ie,kd)*rp(ie)
      enddo
      call valint(tmp3,res3)

      do ie = 2,n
         tmp4(ie)=pf(ie,kb)*gs(ie,kd)*rp(ie)
      enddo
      call valint(tmp4,res4)

      skll=-(kapa-kapc-1)*res1*(kapb-kapd+1)*res3*cstp
      skls=-(kapa-kapc-1)*res1*(kapb-kapd-1)*res4*cstp
      sksl=-(kapa-kapc+1)*res2*(kapb-kapd+1)*res3*cstp
      skss=-(kapa-kapc+1)*res2*(kapb-kapd-1)*res4*cstp

      else

      call dpbdt(kc)

      do ie = 2,n
         tmp1(ie)=(gl(ie,ka)*tmp5(ie)/h)+(ieta(laa,lca)*
     :                        gl(ie,ka)*pf(ie,kc)*rpor(ie))
      enddo
      call valint(tmp1,res1)

      do ie = 2,n
         tmp2(ie)=(gs(ie,ka)*tmp6(ie)/h)+(izeta(lab,lcb)*
     :                        gs(ie,ka)*qf(ie,kc)*rpor(ie))
      enddo
      call valint(tmp2,res2)

      do ie = 2,n
         tmp3(ie)=(gl(ie,kd)*pf(ie,kb)*dg1(ie,kd)*rp(ie))+
     :             (ieta(lba,lda)*gl(ie,kd)*pf(ie,kb)*rpor(ie))
      enddo
      call valint(tmp3,res3)

      do ie = 2,n
         tmp4(ie)=(gs(ie,kd)*qf(ie,kb)*(dg1(ie,kd)+dg2(ie,kd))*
     :rp(ie))+(izeta(lbb,ldb)*gs(ie,kd)*qf(ie,kb)*rpor(ie))
      enddo
      call valint(tmp4,res4)

      skll=-res1*res3*cstp/c**2
      skls=-res1*res4*cstp/c**2
      sksl=-res2*res3*cstp/c**2
      skss=-res2*res4*cstp/c**2

      endif

      endif
      
      return
      end

      subroutine dfistp(ka,kb,kc,kd,orba,orbb,orbc,orbd,iaa,
     :iab,iac,iad,laa,lab,lba,lbb,lca,lcb,lda,ldb,k1,stpval)
      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /basis/gl(MN,MNBAS),gs(MN,MNBAS)
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     : /propopt/istp,ialpvr
     :  /quan1/kappe(MNBAS),orbje(MNBAS),kappc(MNBAS),orbjc(MNBAS)
     :  /stpfac/dg1(MN,MNBAS),dg2(MN,MNBAS)
     :  /atmp/tmp5(MN),tmp6(MN)
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c

      do i=1,MN
         tmp1(i)=zero
         tmp2(i)=zero
         tmp3(i)=zero
         tmp4(i)=zero
         tmp5(i)=zero
         tmp6(i)=zero
      enddo

      stpval=zero
      fact=zero

      kapa=-idint((orba+half)*iaa)
      kapb=-idint((orbb+half)*iab)
      kapc=-idint((orbc+half)*iac)
      kapd=-idint((orbd+half)*iad)

      if(kapa.ne.kappe(ka).or.kapb.ne.kappe(kb).or.
     :kapc.ne.kappe(kc).or.kapd.ne.kappe(kd))then
      Print*,'Quantum numbers are wrongly assigned'
      Stop
      Endif

      if(k1.ne.1)then
      print*,'Wrong in k1 value as 1',k1
      stop
      endif

      fact=selectc(iaa,iac,orba,orbc,k1)*
     : selectc(iab,iad,orbb,orbd,k1)

      if(fact.ne.zero)then

      if(istp.eq.3)then

      do ie = 2,n
         tmp1(ie)=pf(ie,ka)*qf(ie,kc)*rp(ie)
      enddo
      call valint(tmp1,res1)

      do ie = 2,n
         tmp2(ie)=qf(ie,ka)*pf(ie,kc)*rp(ie)
      enddo
      call valint(tmp2,res2)

      do ie = 2,n
         tmp3(ie)=qf(ie,kb)*pf(ie,kd)*rp(ie)
      enddo
      call valint(tmp3,res3)

      do ie = 2,n
         tmp4(ie)=pf(ie,kb)*qf(ie,kd)*rp(ie)
      enddo
      call valint(tmp4,res4)

      stpval=-((kapa-kapc-1)*res1+(kapa-kapc+1)*res2)*
     :        ((kapb-kapd-1)*res4+(kapb-kapd+1)*res3)*c**2

      else

      call dpbdt(kc)

      do ie = 2,n
         tmp1(ie)=(pf(ie,ka)*tmp5(ie)/h)+(ieta(laa,lca)*
     :                        pf(ie,ka)*pf(ie,kc)*rpor(ie))
      enddo
      call valint(tmp1,res1)

      do ie = 2,n
         tmp2(ie)=(qf(ie,ka)*tmp6(ie)/h)+(izeta(lab,lcb)*
     :                        qf(ie,ka)*qf(ie,kc)*rpor(ie))
      enddo
      call valint(tmp2,res2)

      call dpbdt(kd)

      do ie = 2,n
         tmp3(ie)=(pf(ie,kb)*tmp5(ie)/h)+(ieta(lba,lda)*
     :                        pf(ie,kb)*pf(ie,kd)*rpor(ie))
      enddo
      call valint(tmp3,res3)

      do ie = 2,n
         tmp4(ie)=(qf(ie,kb)*tmp6(ie)/h)+(izeta(lbb,ldb)*
     :                        qf(ie,kb)*qf(ie,kd)*rpor(ie))
      enddo
      call valint(tmp4,res4)

      stpval=-(res1+res2)*(res3+res4)

      endif

      endif
      
      return
      end


c ********************************************************************
c FUNCTION clrx
c ********************************************************************
c   The value of clrx is the 3-j symbol:
c
c                    ( ja        k        jb  )
c                    ( 1/2       0       -1/2 )
c
c The  k's are kappa angular quantum numbers. the formula is taken from D.M.
c Brink and G R Satchler, <Angular Momentum>, second edition (Oxford:
c Clarendon Press, 1968), p 138.   The logarithms of the first  mfact
c factorials must be available in  common/facts/ for this program to function
c correctly. note that  n!  is stored in fact(n+1).
c----------------------------------------------------------------------
      real*8 function clrx(kappaa,k,kappab)

      implicit real*8(a-h, o-z)
      dimension gam(500)
      Common/cons/zero,half,tenth,one,two,three,ten

       call factt(gam)
c----------------------------------------------------------------------
c               Determine the absolute values of the kappas
c----------------------------------------------------------------------
      ka =iabs (kappaa)
      kb =iabs (kappab)
c----------------------------------------------------------------------
c                   Perform the triangularity check
c----------------------------------------------------------------------
      if ((iabs(ka-kb) .le. k) .and. (ka+kb-1 .ge. k)) then
c----------------------------------------------------------------------
c   Triangularity satisfied; compute the 3j coefficient
c   begin with the logarithm of the square of the leading term
c----------------------------------------------------------------------
         exptrm = -dlog (dble (ka*kb))
c----------------------------------------------------------------------
c    Compute the logarithm of the square root of the leading term and the
c    factorial part that doesn't depend on the parity of ka+kb+k (the delta
c    factor)
c----------------------------------------------------------------------
         kapkb = ka+kb
         kabkp = kapkb+k
         kamkb = ka-kb
         kbmka = kb-ka
         exptrm = half*(exptrm+gam(kapkb-k  )+gam(kamkb+k+1)
     :                        +gam(kbmka+k+1)-gam(kabkp  +1) )
c----------------------------------------------------------------------
c   The remainder depends on the parity of ka+kb+k
c----------------------------------------------------------------------
         if (mod (kabkp,2) .eq. 0) then
c----------------------------------------------------------------------
c   Computation for even parity case
c   Includes the phase factor: a minus sign if necessary
c----------------------------------------------------------------------
            if (mod (3*kabkp/2,2) .eq. 0) then
               clrx =  one
            else
               clrx = -one
            endif
c----------------------------------------------------------------------
c   Include the contribution from the factorials
c----------------------------------------------------------------------
            exptrm = exptrm+gam((kabkp  +2)/2)-gam((kapkb-k  )/2)
     :                     -gam((kamkb+k+2)/2)-gam((kbmka+k+2)/2)

         else
c----------------------------------------------------------------------
c                   computation for odd parity case
c         include the phase factor: a minus sign if necessary
c----------------------------------------------------------------------
            if (mod ((3*kabkp-1)/2,2) .eq. 0) then
               clrx =  one
            else
               clrx = -one
            endif
c----------------------------------------------------------------------
c   include the contribution from the factorials
c----------------------------------------------------------------------
            exptrm = exptrm+gam((kabkp  +1)/2)-gam((kapkb-k+1)/2)
     :                     -gam((kamkb+k+1)/2)-gam((kbmka+k+1)/2)

         endif
c----------------------------------------------------------------------
c   final assembly
c----------------------------------------------------------------------
         clrx = clrx*dexp (exptrm)

      else
c----------------------------------------------------------------------
c   triangularity violated; set the coefficient to zero
c----------------------------------------------------------------------
         clrx = zero

      endif

      return
      end
c ********************************************************************


c ********************************************************************
c SUBROUTINE factt
c ********************************************************************
c Calculates the logs  of factorials required by the racah coefficient
c routine  dracah. written by n.s. scott.
c---------------------------------------------------------------------
      subroutine factt(gam)

      implicit real*8(a-h,o-z)
      dimension gam(500)
      common /cons/zero,half,tenth,one,two,three,ten
      data thirty/3.0d 01/
*
      gam(1) = one
      gam(2) = one
      x = two
      do i = 3,30
         gam(i) = gam(i-1)*x
         x = x+one
      enddo
      do  i = 1,30
         gam(i) = dlog(gam(i))
      enddo
      x = thirty
      do i = 31,500
         gam(i) = gam(i-1)+dlog(x)
         x = x+one
      enddo
      return
      end


      subroutine scfeng(cstp,ival,valistp)
c ********************************************************************
c This evaluates the direct and exchange matrix elements numericaly.
c P(i,j) is the density matrix.
c------------------------------------------------------------------
      implicit real*8(a-h,o-z)

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /basc/nbasis,nocc,mbasis
     :  /skip/nskipe(MNS),nskipc(MNS)
     :  /dfm/emt(MNB,MNB),df_single(MNB,MNB),df_two(MNB,MNB)
     :  /quan1/kappe(MNBAS),orbje(MNBAS),kappc(MNBAS),orbjc(MNBAS)
     :  /quan2/lorba(MNS),lorbb(MNS),iqe(MNBAS),iqc(MNBAS)
     : /propopt/istp,ialpvr

C New common blocks

c------------------------------------------------------------------
c Local dimension
c------------------------------------------------------------------
      dimension mulk(MNS+1)

      do ia = 1,MNB
      do ib = 1,MNB
      df_two(ib,ia) = zero
      enddo
      enddo

      valistp=zero

      do isym=1,nsym
      do ia=1,nbas(isym)
      ka=ia+nskipe(isym)
      ja=ia+2*nskipe(isym)
      index1=ia+nbas(isym)+2*nskipe(isym)

      laa=lorba(isym)
      lab=lorbb(isym)
      orba=orbje(ka)
      iaa=iqe(ka)
      kapa=kappe(ka)

      do ib=1,nbas(isym)
      index2=ib+nbas(isym)+2*nskipe(isym)
      jb=ib+2*nskipe(isym)
      kb=ib+nskipe(isym)

      lba=lorba(isym)
      lbb=lorbb(isym)
      kapb=kappe(kb)
      orbb=orbje(kb)
      iab=iqe(kb)

      df_two(ja,jb)=zero
      df_two(ja,index2)=zero
      df_two(index1,jb)=zero
      df_two(index1,index2)=zero

      do jsym=1,nsym
      nn=nocorb(jsym)

      if(nn.ne.0)then

      do ic=1,nocorb(jsym)
      jc=ic+nskipc(jsym)
      kc=ic+nskipe(jsym)

      lca=lorba(jsym)
      lcb=lorbb(jsym)
      kapc=kappc(jc)
      orbc=orbjc(jc)
      iac=iqc(jc)

c-------------------------------------------------------------------
c     direct part
c-------------------------------------------------------------------

      fact=(two*orbc+one)

      call findk(orba,orbb,orbc,orbc,mulk,maxk)

      if(maxk.ne.0) then
      do ik=1,maxk
      k1=mulk(ik)

      if(k1.eq.0)then

      rkll=zero
      rkss=zero

c-------------------------------------------------------------------
c     Direct Coulomb part
c-------------------------------------------------------------------

      call slatcould(ka,kc,kb,kc,orba,orbc,orbb,orbc,
     :iaa,iac,iab,iac,k1,rkll,rkss)

      df_two(ja,jb)=df_two(ja,jb)+rkll*fact

      df_two(index1,index2)=df_two(index1,index2)+rkss*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     exchange part
c-------------------------------------------------------------------

      call findk(orba,orbc,orbc,orbb,mulk,maxk)

      if(maxk.ne.0)then
      do ik=1,maxk
      k1=mulk(ik)
      ak1=dfloat(k1)

      rkll=zero
      rkls=zero
      rksl=zero
      rkss=zero

      bkll=zero
      bkls=zero
      bksl=zero
      bkss=zero

      skll=zero
      skls=zero
      sksl=zero
      skss=zero
      val=zero

      fact=(two*orbc+one)*clrx(kapa,k1,kapc)**2

c-------------------------------------------------------------------
c     Exchange Coulomb part
c-------------------------------------------------------------------

      call slatcoule(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,k1,rkll,rkls,rksl,rkss)

c-------------------------------------------------------------------
c     Exchange Breit part
c-------------------------------------------------------------------

      if(ibrt.ne.0)then
      call slatbrt(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,k1,bkll,bkls,bksl,bkss)
      endif

c-------------------------------------------------------------------
c     Exchange isotope part
c-------------------------------------------------------------------

      if(k1.eq.1)then

      if(istp.eq.1.or.istp.eq.2)then
      call slatistp(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,laa,lab,lca,lcb,lca,lcb,lba,lbb,
     :k1,skll,skls,sksl,skss,cstp)
      endif

      if(ka.eq.ival.and.kb.eq.ival)then
      if(istp.eq.3.or.istp.eq.4)then
      call dfistp(ka,kc,kc,kb,orba,orbc,orbc,orbb,iaa,iac,
     :iac,iab,laa,lab,lca,lcb,lca,lcb,lba,lbb,k1,val)
      valistp=valistp-val*fact
      endif
      endif

      endif

      df_two(ja,jb)=df_two(ja,jb)-(rkll+bkll+skll)*fact

      df_two(ja,index2)=df_two(ja,index2)-(rkls+bkls+skls)*fact

      df_two(index1,jb)=df_two(index1,jb)-(rksl+bksl+sksl)*fact

      df_two(index1,index2)=df_two(index1,index2)-
     :(rkss+bkss+skss)*fact

      enddo
      endif

      enddo
      endif
      enddo

      if(ivn.gt.0)then

      do ksym=1,nsym
      do ic=1,nbas(ksym)
      kc=ic+nskipe(ksym)
      orbc=orbje(kc)
      iac=iqe(kc)
      kapc=kappe(kc)

      if(kc.eq.ival)then

c-------------------------------------------------------------------
c     direct part
c-------------------------------------------------------------------

      sgn=(two*orbc+one)*(two*orba+one)

      call findk(orba,orbb,orbc,orbc,mulk,maxk)

      if(maxk.ne.0)then
      do ik=1,maxk
      k1=mulk(ik)
      ak1=dfloat(k1)

      rkll=zero
      rkss=zero

cC mj=orba
c      f1=dr(orba,ak1,orbb,-orba,zero,orbb)
c      f2=dr(orbc,ak1,orbc,-orbc,zero,orbc)
c      f3=dr(orba,ak1,orbb,half,zero,-half)
c      f4=dr(orbc,ak1,orbc,half,zero,-half)
c      fact=sgn*f1*f2*f3*f4*(-1)**(one+orba+orbc)

C mj=half
      f1=dr(orba,ak1,orbb,-half,zero,half)
      f2=dr(orbc,ak1,orbc,-half,zero,half)
      f3=dr(orba,ak1,orbb,half,zero,-half)
      f4=dr(orbc,ak1,orbc,half,zero,-half)
      fact=sgn*f1*f2*f3*f4

cC Average mj
c      if(k1.eq.0)then
c      fact=one
c      if(ka.eq.ival.or.kb.eq.ival)fact=fact/(2*orba+one)

      call slatcould(ka,kc,kb,kc,orba,orbc,orbb,orbc,
     :iaa,iac,iab,iac,k1,rkll,rkss)

      df_two(ja,jb)=df_two(ja,jb)+rkll*fact

      df_two(index1,index2)=df_two(index1,index2)+rkss*fact

c      endif

      enddo
      endif

c-------------------------------------------------------------------
c     exchange part
c-------------------------------------------------------------------

      call findk(orba,orbc,orbc,orbb,mulk,maxk)

      if(maxk.ne.0)then
      do ik=1,maxk
      k1=mulk(ik)
      ak1=dfloat(k1)

      rkll=zero
      rkls=zero
      rksl=zero
      rkss=zero

      bkll=zero
      bkls=zero
      bksl=zero
      bkss=zero

cC mj=orba
c      f1=dr(orba,ak1,orbc,-orba,orba-orbc,orbc)
c      f2=dr(orbc,ak1,orbb,-orbc,orbc-orbb,orbb)
c      f3=dr(orba,ak1,orbc,half,zero,-half)
c      f4=dr(orbc,ak1,orbb,half,zero,-half)
c      fact=sgn*f1*f2*f3*f4

C mj=half
      f1=dr(orba,ak1,orbc,-half,zero,half)
      f2=dr(orbc,ak1,orbb,-half,zero,half)
      f3=dr(orba,ak1,orbc,half,zero,-half)
      f4=dr(orbc,ak1,orbb,half,zero,-half)
      fact=sgn*f1*f2*f3*f4

cC Average mj
c      fact=clrx(kapa,k1,kapc)**2
c      if(ka.eq.ival.or.kb.eq.ival)fact=fact/(2*orba+one)

      call slatcoule(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,k1,rkll,rkls,rksl,rkss)

      if(ibrt.gt.0)then
      call slatbrt(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,k1,bkll,bkls,bksl,bkss)
      endif

      df_two(ja,jb)=df_two(ja,jb)-(rkll+bkll)*fact

      df_two(ja,index2)=df_two(ja,index2)-(rkls+bkls)*fact

      df_two(index1,jb)=df_two(index1,jb)-(rksl+bksl)*fact

      df_two(index1,index2)=df_two(index1,index2)-(rkss+bkss)*fact

      enddo
      endif

      endif

      enddo
      enddo

      endif

      enddo
      enddo

      enddo

      return
      end


      subroutine output(rho)
      implicit real*8(a-h,o-z)

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /basc/nbasis,nocc,mbasis
     :      /skip/nskipe(MNS),nskipc(MNS)
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     :      /freez/minorb(MNS),maxorb(MNS)

      common/jvalue/orbj(MNS),kap(MNS) 
     :      /nucinf/amass,z
     :      /splflags/MRCC,IBASOPT

      Common/totenergies/en1,en2,en3,engMP2

       LOGICAL IORTHO

       dimension icntb(MNS),icntc(MNS),rho(MN)

        IORTHO = .FALSE.

       IF (IBASOPT) CALL CORRENG
c      IF(intg.eq.1.or.intg.eq.2)CALL DHFMAT(RHO)
c      IF(intg.eq.2.or.intg.eq.3)CALL TWOINT


      IF(ihbd.eq.1.or.ivn.eq.1)THEN

C Diagonalizes only for V^N-1
      Call SETMATHYBVN1(RHO)
      Call scfhybvn1

cC Alternatively, diagonalizes for all and keeps only V^N-1 
c      Call SETMATHYBVN(RHO)
c      Call scfhybvn

      ENDIF

      write(STDOUT,*)" "
      write(STDOUT,208)

      write(IWFN)h,n
      write(IWFN)(r(i),i=1,n),(rp(i),i=1,n),(rpor(i),i=1,n)

c     write only the relevant orbitals i.e skip writting the
c     negavtive energy orbitals

      do isym=1,nsym

       icntb(isym)=0
       icntc(isym)=0

      do ia=1,nbas(isym)
       ka=ia+nskipe(isym)

       if(ieng.eq.0)goto 8889

       if(maxorb(isym).ge.0.and.ia.gt.maxorb(isym)) goto 8888

        icntc(isym)=icntc(isym)+1

         if(minorb(isym).ge.0.and.ia.le.minorb(isym)) goto 8888

 8889 continue
         
          icntb(isym)=icntb(isym)+1

      write(15,*)ia,eng1(ka)
      write(IWFN)eng1(ka)
      write(IWFN)(pf(ii,ka),ii=1,n),(qf(ii,ka),ii=1,n)

 8888 continue

      IF(IORTHO)THEN
      do ib=1,ia
      kb=ib+nskipe(isym)

      r0=rint(ia,ib,0)
      r1=zero
      r2=zero
      r3=zero

      if(ia.eq.ib)r1=rint(ia,ib,1)
      if(ia.eq.ib)r2=rint(ia,ib,2)
      if(ia.eq.ib)r3=rint(ia,ib,-1)

      write(STDOUT,'(2i4,2x,4(2x,D15.9))')ka,kb,r0,r1,r2,r3
      enddo

C!      It is commented out by HSN on Apr. 02, 2010

      ENDIF ! IORTHO

      enddo
      enddo

      write(*,"((a20),x,I4)")'No. of sym:', nsym

      if(ieng.eq.0)write(*,"((a20),x,15I4)")
     :                              'No. of basis:',(icntb(i),i=1,nsym)
      if(ieng.eq.1)write(*,"((a20),x,15I4)")
     :                              'No. of basis:',(icntc(i),i=1,nsym)

      write(*,"((a20),x,15I4)")'No. of occ. orbs:',(nocorb(i),i=1,nsym)

      if(ieng.eq.1)write(*,"((a20),x,15I4)")
     :                         'No. of frz. orbs:',(minorb(i),i=1,nsym)
      if(ivn.eq.2)write(*,"((a20),x,I4)")'No. of valence orb:',1


      write(15,"((a20),x,I4)")'No. of sym:', nsym

      if(ieng.eq.0)write(15,"((a20),x,15I4)")
     :                              'No. of basis:',(icntb(i),i=1,nsym)
      if(ieng.eq.1)write(15,"((a20),x,15I4)")
     :                              'No. of basis:',(icntc(i),i=1,nsym)

      write(15,"((a20),x,15I4)")'No. of occ. orbs:',(nocorb(i),i=1,nsym)

      if(ieng.eq.1)write(15,"((a20),x,15I4)")
     :                          'No. of frz. orbs:',(minorb(i),i=1,nsym)

      if(ivn.eq.2)write(15,"((a20),x,I4)")'No. of valence orb:',1


C     To append the nuclear density to wfn.dat file
      if(intg.eq.1)write(IWFN)(rho(ii),ii=1,n)

      IF(IBASOPT) THEN
C     To print out all the Bound Orbital Energies to a seperate file
      write(STDENG,*)''
      do i=1,nbasis
      if(eng1(i).le.0.0D0)write(STDENG,"(I4,100F20.10)")i,eng1(i)
      enddo

C     Print the MP2 energy at the end
            
      write(STDENG,*)''
      WRITE(STDENG,"(4x,(a22),F25.15)")' TOTAL MP2 ENERGY = ', engMP2

      ENDIF

C**********************************************************************C
      IF (MRCC) THEN

C     INPUT GENERATION FOR THE INTERFACE FILE ONETWOINT.F

C     Note : For Frozen core calculation with MRCC do not freeze core in SCF step 
C     Just Freeze the desired core in the interface step which generates fort.55

      open(STD12INT,file='input.int',form="formatted",status="unknown")

      if(ieng.eq.0)write(STD12INT,"(16I4)") nsym,(icntb(i),i=1,nsym)
      if(ieng.eq.1)write(STD12INT,"(16I4)") nsym,(icntc(i),i=1,nsym)
      write(STD12INT,"(4x,15I4)") (nocorb(i),i=1,nsym)
      write(STD12INT,"(4x,15I4)") (minorb(i),i=1,nsym)    

      numofjm = 0
      do isym=1,nsym
      if(ieng.eq.0)numofjm = numofjm + icntb(isym)*(int(2*orbj(isym)+1))
      if(ieng.eq.1)numofjm = numofjm + (icntc(isym)-minorb(isym))
     :                                     *(int(2*orbj(isym)+1))
      enddo

      write(STD12INT,"(4x,3I4)")ibrt,0,0   !iprop=0 and igag=0 default options
      write(STD12INT,"(4x,2I4)")numofjm,int(z)
      close(STD12INT)
     
      ENDIF   ! MRCC
C**********************************************************************C

 208  format(2x,'THESE ARE THE ORTHOGONALITY OF THE  ',
     :  1x,'ORBITALS',//,2x,'Sl. No.','    <A|B>      ',
     :  '    <A|R|B>     ','    <A|R*R|B>   ',
     :  '     <A|(1/R)|B>',//)

 210  format(///,2x,'THESE ARE THE DIPOLE MOMENTS OF THE  ',
     :  1x,'ORBITALS',//,2x,'Sl. No.',5x,'DIPOLE MOMENTS'/)

      return
      end


      function rint(i,j,k)
      implicit real*8(a-h,o-z)

      common/cons/zero,half,tenth,one,two,three,ten
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)

      rint=zero

      Do ic=1,n
      tmp1(ic)=zero
      EndDo

      tmp1(1)=zero

      do 1 l=2,n
      tmp1(l)=(r(l)**k)*(pf(l,i)*pf(l,j)+qf(l,i)*qf(l,j))*rp(l)
    1 continue

      call valint(tmp1,result)

      rint=result

      return
      end


      SUBROUTINE DHFMAT(RHO)
c ********************************************************************
C  This subroutine determines matrix elements of kinetic, rest mass energy,
c  nuclear potential and Dirac-Fock matrix elements
c   <a| c alpha.p + (beta-1) c^2 + V_nuc(r) + u_DF (r) | b >
c --------------------------------------------------------------------

      implicit real*8(a-h,o-z)

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /basc/nbasis,nocc,mbasis
     :  /skip/nskipe(MNS),nskipc(MNS)
     :  /quan1/kappe(MNBAS),orbje(MNBAS),kappc(MNBAS),orbjc(MNBAS)
     :  /quan2/lorba(MNS),lorbb(MNS),iqe(MNBAS),iqc(MNBAS)
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     :      /atmp/tmp5(MN),tmp6(MN)

C New common block
      Common/dhfm/dhf(MNB,MNB)

c ---------------------------------------------------------------
c Local dimension
c ---------------------------------------------------------------
      dimension tmp(MN),rho(MN),mulk(MNS+1)

      do i=1,MNB
      do j=1,MNB
      dhf(i,j)=zero
      enddo
      enddo

      do isym=1,nsym

      do ib=1,nbas(isym)
      kb=ib+nskipe(isym)

      orbb=orbje(kb)
      iab=iqe(kb)
      kapb=kappe(kb)

      do i=1,MN
      tmp5(i)=zero
      tmp6(i)=zero
      enddo
      call dpbdt(kb)

      do ia=1,nbas(isym)
      ka=ia+nskipe(isym)

      orba=orbje(ka)
      iaa=iqe(ka)
      kapa=kappe(ka)

      piece1=zero
      piece2=zero
      piece3=zero
      piece4=zero
      result=zero
      result1=zero
      result2=zero

      do i=1,MN
      tmp(i)=zero
      tmp1(i)=zero
      tmp2(i)=zero
      tmp3(i)=zero
      tmp4(i)=zero
      enddo

      result=zero

      piece1=zero
      do ir=2,n
      tmp(ir)=(pf(ir,ka)*pf(ir,kb)+qf(ir,ka)*qf(ir,kb))*
     :rho(ir)*rpor(ir)
      enddo
      call valint(tmp,result)

      piece1=-result

      piece2=zero
      do ir=2,n
      tmp(ir)=(pf(ir,ka)*qf(ir,kb)+pf(ir,kb)*qf(ir,ka))*rpor(ir)
      enddo
      call valint(tmp,result)

      piece2=c*dble(kapb)*result

      piece3=zero
      do ir=2,n
      tmp(ir)=qf(ir,ka)*tmp5(ir)-pf(ir,ka)*tmp6(ir)
      enddo
      call valint(tmp,result)

      piece3=c*result/h

      do ir=2,n
      tmp(ir)=qf(ir,ka)*qf(ir,kb)*rp(ir)
      enddo
      call valint(tmp,result)

      piece4=-two*c*c*result

      do i=1,MN
      tmp(i)=zero
      tmp1(i)=zero
      enddo

      result1=zero
      result2=zero

      do jsym=1,nsym
      nn=nocorb(jsym)

      if(nn.ne.0)then

      do ic=1,nocorb(jsym)
      jc=ic+nskipc(jsym)
      kc=ic+nskipe(jsym)

      kapc=kappc(jc)
      orbc=orbjc(jc)
      iac=iqc(jc)

c-------------------------------------------------------------------
c     direct part
c-------------------------------------------------------------------

      fact=zero

      call findk(orba,orbb,orbc,orbc,mulk,maxk)

      if(maxk.ne.0)then
      do ik=1,maxk
      k1=mulk(ik)

      if(k1.eq.0)then

      fact1=zero

      fact1=selectc(iaa,iab,orba,orbb,k1)*
     : selectc(iac,iac,orbc,orbc,k1)

      if(fact1.ne.zero)then

      call yzk(k1,kc,kc)

      resd=zero
      do ir = 2,n
         tmp(ir)=(pf(ir,ka)*pf(ir,kb)+qf(ir,ka)*qf(ir,kb))*
     :rpor(ir)*tmp1(ir)
      enddo
      call valint(tmp,resd)

      fact=(two*orbc+one)
      result1=result1+resd*fact*fact1

      endif

      endif

      enddo
      endif

c-------------------------------------------------------------------
c     exchange part
c-------------------------------------------------------------------

      call findk(orba,orbc,orbc,orbb,mulk,maxk)

      if(maxk.ne.0)then
      do ik=1,maxk
      k1=mulk(ik)

      fact1=zero

      fact1=selectc(iaa,iac,orba,orbc,k1)*
     : selectc(iac,iab,orbc,orbb,k1)

      if(fact1.ne.zero)then

      call yzk(k1,kc,kb)

      rese=zero
      do ir = 2,n
         tmp(ir)=(pf(ir,ka)*pf(ir,kc)+qf(ir,ka)*qf(ir,kc))*
     :rpor(ir)*tmp1(ir)
      enddo
      call valint(tmp,rese)

      fact=(two*orbc+one)*clrx(kapa,k1,kapc)*clrx(kapb,k1,kapc)
      result2=result2-rese*fact*fact1

      endif

      enddo
      endif

      enddo
      endif
      enddo
       
      dhf(ka,kb)=piece1+piece2+piece3+piece4+result1+result2

      write(10,*)ka,kb,piece1,piece2,piece3,piece4
c      write(10,*)ka,kb,(piece1+piece2+piece3+piece4),
c     :(result1+result2),dhf(ka,kb)

      enddo
      enddo

      enddo

      RETURN
      END

      SUBROUTINE DPBDT(J)
*CCC
*CCC      $Id: DPBDT.f,v 1.4 1994/09/15 13:34:04 pn Exp $
*CCC
*CCC      $Log: DPBDT.f,v $
*CCC      Revision 1.4  1994/09/15 13:34:04  pn
*CCC      Version 2c  31 May 1994
*CCC      Modified version of 2b.
*CCC
*CCC Revision 1.3  1994/03/31  11:57:06  pn
*CCC Version 2b 6 July 1990
*CCC
*CCC
C***********************************************************************

      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :  /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     :      /atmp/tmp5(MN),tmp6(MN)
      COMMON /LIC13/A(13,13)

C
      EQUIVALENCE (A1,A(7,1)), (A2,A(7,2)), (A3,A(7,3)), (A4,A(7,4)),
     A            (A5,A(7,5)), (A6,A(7,6))
C-----------------------------------------------------------------------

      DO I=1,MN
      TMP5(I)=ZERO
      TMP6(I)=ZERO
      ENDDO

C   Compute derivative in three separate regions
C
C   First, points 1 to 6
C
      DO 20 I = 1,6
        HDPBDT = ZERO
        HDQBDT = ZERO
        DO 10 K = 1,13
          AIK = A(I,K)
          HDPBDT = HDPBDT + AIK*PF(K,J)
          HDQBDT = HDQBDT + AIK*QF(K,J)
   10   CONTINUE
        TMP5(I) = HDPBDT
        TMP6(I) = HDQBDT
   20 CONTINUE
C
C   Next, points 7 to N-6
C
C   Special treatment for this region because of the symmetry of
C   the differentiation formula
C
      DO 30 I = 7,N - 6
       TMP5(I)= A1* (PF(I-6,J)-PF(I+6,J)) + A2* (PF(I-5,J)-PF(I+5,J)) +
     A          A3* (PF(I-4,J)-PF(I+4,J)) + A4* (PF(I-3,J)-PF(I+3,J)) +
     B          A5* (PF(I-2,J)-PF(I+2,J)) + A6* (PF(I-1,J)-PF(I+1,J))
       TMP6(I)= A1* (QF(I-6,J)-QF(I+6,J)) + A2* (QF(I-5,J)-QF(I+5,J)) +
     A          A3* (QF(I-4,J)-QF(I+4,J)) + A4* (QF(I-3,J)-QF(I+3,J)) +
     B          A5* (QF(I-2,J)-QF(I+2,J)) + A6* (QF(I-1,J)-QF(I+1,J))
   30 CONTINUE
C
C   Last, points N-5 to N
C
      DO 50 I = N - 5,N
        IROW = I - N + 13
        HDPBDT = ZERO
        HDQBDT = ZERO
        DO 40 K = 1,13
          AIK = A(IROW,K)
          LOC = N - 13 + K
          HDPBDT = HDPBDT + AIK*PF(LOC,J)
          HDQBDT = HDQBDT + AIK*QF(LOC,J)
   40   CONTINUE
        TMP5(I) = HDPBDT
        TMP6(I) = HDQBDT
   50 CONTINUE

      RETURN
      END

      SUBROUTINE CORRENG
c ********************************************************************
C  This subroutine computes MBPT(2) correlation energy
c --------------------------------------------------------------------

      implicit real*8(a-h,o-z)

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /basc/nbasis,nocc,mbasis
     :  /skip/nskipe(MNS),nskipc(MNS)
     :  /quan1/kappe(MNBAS),orbje(MNBAS),kappc(MNBAS),orbjc(MNBAS)
     :  /quan2/lorba(MNS),lorbb(MNS),iqe(MNBAS),iqc(MNBAS)
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     :      /atmp/tmp5(MN),tmp6(MN)
     :      /symmetry/mtbl(MNS,MNS)

      Common/totenergies/en1,en2,en3,engMP2
c ---------------------------------------------------------------
c Local dimension
c ---------------------------------------------------------------
      dimension imulk(MNS+1),jmulk(MNS+1)

      corengd=zero
      corenge=zero

      do i=1,MN
      tmp1(i)=zero
      tmp2(i)=zero
      enddo

      call parityproduct

      do 1 isym=1,nsym
      lasym=lorba(isym)
      do 2 ia=1,nbas(isym)
      ka=ia+nskipe(isym)

      orba=orbje(ka)
      iaa=iqe(ka)
      kapa=kappe(ka)

      if(ia.gt.nocorb(isym))then

      do 3 jsym=1,nsym
      lbsym=lorba(jsym)
      do 4 ib=1,nbas(jsym)
      kb=ib+nskipe(jsym)

      orbb=orbje(kb)
      iab=iqe(kb)
      kapb=kappe(kb)

      if(ib.gt.nocorb(jsym))then

      do 5 ksym=1,nsym
      lcsym=lorba(ksym)
      do 6 ic=1,nbas(ksym)
      kc=ic+nskipe(ksym)

      orbc=orbje(kc)
      iac=iqe(kc)
      kapc=kappe(kc)

      if(ic.le.nocorb(ksym))then

      do 7 lsym=1,nsym
      ldsym=lorba(lsym)
      do 8 id=1,nbas(lsym)
      kd=id+nskipe(lsym)

      orbd=orbje(kd)
      iad=iqe(kd)
      kapd=kappe(kd)

      if(id.le.nocorb(lsym))then

C Direct term
      lacsym=mtbl(lasym,lcsym)
      lbdsym=mtbl(lbsym,ldsym)

      if(lacsym.eq.lbdsym)then
      call findk(orba,orbc,orbb,orbd,imulk,imaxk)
      if(imaxk.ne.0)then
      do ik=1,imaxk
      k1=imulk(ik)
      ak1=dfloat(k1)

      fact=zero

      fact=selectc(iaa,iac,orba,orbc,k1)*
     : selectc(iab,iad,orbb,orbd,k1)

      if(fact.ne.zero)then

      call yzk(k1,kb,kd)

      result=zero

      do ir = 2,n
         tmp2(ir)=(pf(ir,ka)*pf(ir,kc)+qf(ir,ka)*qf(ir,kc))*
     :rpor(ir)*tmp1(ir)
      enddo

      call valint(tmp2,result)

      rk1=result

      call yzk(k1,kd,kb)

      result=zero

      do ir = 2,n
         tmp2(ir)=(pf(ir,kc)*pf(ir,ka)+qf(ir,kc)*qf(ir,ka))*
     :rpor(ir)*tmp1(ir)
      enddo

      call valint(tmp2,result)

      rk2=result

      f1=dr(orba,ak1,orbc,half,zero,-half)
      f2=dr(orbb,ak1,orbd,half,zero,-half)
      f3=dr(orbc,ak1,orba,half,zero,-half)
      f4=dr(orbd,ak1,orbb,half,zero,-half)

      denom=one/(eng1(kc)+eng1(kd)-eng1(ka)-eng1(kb))
      fac=(two*orba+one)*(two*orbc+one)*(two*orbb+one)*
     :(two*orbd+one)*denom*f1*f2*f3*f4

      corengd=corengd+half*fact*fac*rk1*rk2/(2*ak1+one)

C Exchange term
      ladsym=mtbl(lasym,ldsym)
      lbcsym=mtbl(lbsym,lcsym)

      if(ladsym.eq.lbcsym)then
      call findk(orba,orbd,orbb,orbc,jmulk,jmaxk)
      if(jmaxk.ne.0)then
      do jk=1,jmaxk
      k2=jmulk(jk)
      ak2=dfloat(k2)

      fact1=zero

      fact1=selectc(iaa,iad,orba,orbd,k2)*
     : selectc(iab,iac,orbb,orbc,k2)

      if(fact1.ne.zero)then

      call yzk(k1,kd,ka)

      result=zero

      do ir = 2,n
         tmp2(ir)=(pf(ir,kc)*pf(ir,kb)+qf(ir,kc)*qf(ir,kb))*
     :rpor(ir)*tmp1(ir)
      enddo

      call valint(tmp2,result)

      rk2=result

      f3=dr(orbc,ak2,orbb,half,zero,-half)
      f4=dr(orbd,ak2,orba,half,zero,-half)

      isgn=(-1)**(2*orba+ak1+ak2)

      fac=(two*orba+one)*(two*orbc+one)*(two*orbb+one)*
     :(two*orbd+one)*denom*isgn*f1*f2*f3*f4

      call sixj(ak2,orba,orbd,ak1,orbb,orbc,s6j)
      corenge=corenge-half*fact*fact1*fac*rk1*rk2*s6j

      endif

      enddo
      endif
      endif

      endif

      enddo
      endif
      endif

       endif
    8 continue
    7 continue
       endif
    6 continue
    5 continue
       endif
    4 continue
    3 continue
       endif
    2 continue
    1 continue

      write(STDOUT,*)''
      write(STDOUT,*)
     :'****************************************************************'
      write(STDOUT,*)''

      write(STDOUT,"(x,(a37),x,F25.15)")
     :                ' MP2 corr. energy direct part : ',corengd
      write(STDOUT,*)''

      write(STDOUT,"(x,(a37),x,F25.15)")
     :                ' MP2 corr. energy xchange part : ',corenge
      write(STDOUT,*)''

      write(STDOUT,"(x,(a37),x,F25.15)")
     :                ' MP2 correlation energy : ',corengd+corenge

      write(STDOUT,*)''

      engMP2 = en1+corengd+corenge

      write(STDOUT,"(x,(a37),x,F25.15)")
     :                ' Total MP2 energy : ',engMP2

      write(STDOUT,*)''
      write(STDOUT,*)
     :'****************************************************************'
      write(STDOUT,*)''


      RETURN
      END

      SUBROUTINE TWOINT
c ********************************************************************
C  This subroutine computes reduced two-body matrix elements 
c     <ab || (1/r_12) || cd>
c --------------------------------------------------------------------

      implicit real*8(a-h,o-z)

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /basc/nbasis,nocc,mbasis
     :  /skip/nskipe(MNS),nskipc(MNS)
     :  /quan1/kappe(MNBAS),orbje(MNBAS),kappc(MNBAS),orbjc(MNBAS)
     :  /quan2/lorba(MNS),lorbb(MNS),iqe(MNBAS),iqc(MNBAS)
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     :      /atmp/tmp5(MN),tmp6(MN)
     :      /symmetry/mtbl(MNS,MNS)

c ---------------------------------------------------------------
c Local dimension
c ---------------------------------------------------------------
      dimension imulk(MNS+1)

      do i=1,MN
      tmp1(i)=zero
      tmp2(i)=zero
      enddo

      call parityproduct

      twointg=zero

      do 1 isym=1,nsym
      lasym=lorba(isym)
      do 2 ia=1,nbas(isym)
      ka=ia+nskipe(isym)

      orba=orbje(ka)
      iaa=iqe(ka)
      kapa=kappe(ka)

      do 3 jsym=1,nsym
      lbsym=lorba(jsym)
      do 4 ib=1,nbas(jsym)
      kb=ib+nskipe(jsym)

      orbb=orbje(kb)
      iab=iqe(kb)
      kapb=kappe(kb)

      do 5 ksym=1,nsym
      lcsym=lorba(ksym)
      do 6 ic=1,nbas(ksym)
      kc=ic+nskipe(ksym)

      orbc=orbje(kc)
      iac=iqe(kc)
      kapc=kappe(kc)

      do 7 lsym=1,nsym
      ldsym=lorba(lsym)
      do 8 id=1,nbas(lsym)
      kd=id+nskipe(lsym)

      orbd=orbje(kd)
      iad=iqe(kd)
      kapd=kappe(kd)

      lacsym=mtbl(lasym,lcsym)
      lbdsym=mtbl(lbsym,ldsym)

      if(lacsym.eq.lbdsym)then
      call findk(orba,orbc,orbb,orbd,imulk,maxk)

      if(maxk.ne.0)then
      do ik=1,maxk
      k1=imulk(ik)
      ak1=dfloat(k1)

      fact=zero

      fact=selectc(iaa,iac,orba,orbc,k1)*
     : selectc(iab,iad,orbb,orbd,k1)

      if(fact.ne.zero)then

      call yzk(k1,kb,kd)

      result=zero

      do ir = 2,n
         tmp2(ir)=(pf(ir,ka)*pf(ir,kc)+qf(ir,ka)*qf(ir,kc))*
     :rpor(ir)*tmp1(ir)
      enddo

      call valint(tmp2,result)

      isign=(-1)**(orba+orbb+one+ak1)
      fact1=dsqrt((two*orba+one)*(two*orbc+one))
      fact2=dsqrt((two*orbb+one)*(two*orbd+one))
      f1=dr(orba,ak1,orbc,half,zero,-half)
      f2=dr(orbb,ak1,orbd,half,zero,-half)

      twointg=fact*fact1*fact2*f1*f2*isign*result

      write(11,*)ka,kb,kc,kd,k1,twointg

      endif

      enddo
      endif
      endif

    8 continue
    7 continue
    6 continue
    5 continue
    4 continue
    3 continue
    2 continue
    1 continue

      RETURN
      END


      subroutine parityproduct
      implicit real*8 (a-h,o-z)

C New common block
      common/symmetry/mtbl(MNS,MNS)

c     1 stands for even and 2 stands for odd parity

* i and j value depends (2*jmax+1)/2, where jmax is the j value of
* the maximum symmetry has been considered. For example, here we have
* taken up to i_{11/2} symmetry and hence k = (2* 11/2 + 1)/2 = 6

      do i=1,11
      do j=1,11

      if(mod(i+j,2).eq.0)then
      mtbl(i,j)=1
      else
      mtbl(i,j)=2
      endif

      enddo
      enddo

      return
      end


*-----------------------------------------------------------------
*    programe to calculate the 3j symbols
*-----------------------------------------------------------------
      function dr(cj1,cj2,j3,cm1,cm2,m3)
      implicit real*8(a-h,o-z)

      real*8 m3,j3
      ddr=0.d0
      if(dabs(cm1+cm2+m3).gt.0.1d0)go to 80
      if(dabs(cm1).gt.cj1+0.1d0)go to 80
      if(dabs(cm2).gt.cj2+0.1d0)go to 80
      if(dabs(cm1+cm2).gt.j3+0.1d0)go to 80
      if((cj1+cj2).lt.j3-0.1d0)go to 80
      if(dabs(cj1-cj2).gt.j3+0.1d0)go to 80
      xnum=factl(j3+cj1-cj2,j3+cj2+cj1+1)*factl(j3-cj1+cj2,cj1-cm1)*
     1factl(cj1+cj2-j3,cj1+cm1)*factl(j3+m3,cj2-cm2)*
     2factl(j3-m3,cj2+cm2)
      xnum=dsqrt(xnum)
      a=-cj1+cm1
      b=cj2-cj1-m3
      c=cj2+j3+cm1
      d=-cj1+cj2+j3
      e=j3-m3
      numin=0
      if(b.gt.0.d0)numin=b+0.1d0
      numax=c+0.1d0
      if(c.gt.d+0.1d0)numax=d+0.1d0
      f=numax
      if(f.gt.e+0.1d0)numax=e+0.1d0
      cj=cj1+cm2-m3
      kmi=numin+1
      kma=numax+1
      do 1 k=kmi,kma
      ck=k-1
      term=factl(1.d0,d-ck)*factl(c-ck,e-ck)*factl(ck-a,ck-b)*
     :factl(1.d0,ck)
      ci=dmod(cj+ck,2.d0)
      if(dabs(ci).gt.0.5d0)term=-term
      ddr=ddr+term
   1  continue
      ddr=ddr*xnum
  80  dr=ddr
      return
      end

*************************************************************************
C Calculates factorial
      function factl(cj1,cj2)
      implicit real*8(a-h,o-z)

      n=idint(dabs(cj1-cj2))
      r=1.d0
      if(n.eq.0)go to 4
      k=1
      if(cj1.lt.cj2)k=-1
      cjx=dmax1(cj1,cj2)
      cjy=dmin1(cj1,cj2)
      hr=cjy
      do 1 i=1,n
      hr=hr+1.d0
      r=r*hr
  1   continue
      if(k.lt.0)r=1.d0/r
  4   factl=r
      return
      end


      subroutine slathybd(ka,kb,kc,kd,orba,orbb,orbc,orbd,
     :iaa,iab,iac,iad,k1,rkll,rkss)
      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)

      do i=1,MN
         tmp1(i)=zero
         tmp2(i)=zero
         tmp3(i)=zero
         tmp4(i)=zero
      enddo

      rkll=zero
      rkss=zero
      fact=zero

      fact=selectc(iaa,iac,orba,orbc,k1)*
     : selectc(iab,iad,orbb,orbd,k1)

      if(fact.ne.zero)then

      call yzk(k1,kb,kd)

      do ir = 2,n
         tmp2(ir)=pf(ir,ka)*pf(ir,kc)*rpor(ir)*tmp1(ir)
      enddo
      call valint(tmp2,rkll)

      do ir = 2,n
         tmp3(ir)=qf(ir,ka)*qf(ir,kc)*rpor(ir)*tmp1(ir)
      enddo

      call valint(tmp3,rkss)

      endif

      return
      end


      subroutine slathybe(ka,kb,kc,kd,orba,orbb,orbc,orbd,
     :iaa,iab,iac,iad,k1,rkll,rkls,rksl,rkss)
      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)

      do i=1,MN
         tmp1(i)=zero
         tmp2(i)=zero
         tmp3(i)=zero
         tmp4(i)=zero
      enddo

      rkll=zero
      rkls=zero
      rksl=zero
      rkss=zero
      fact=zero

      fact=selectc(iaa,iac,orba,orbc,k1)*
     : selectc(iab,iad,orbb,orbd,k1)

      if(fact.ne.zero)then

      call yzkhybe(k1,kb,kd)

      do ir = 2,n
         tmp2(ir)=pf(ir,ka)*pf(ir,kc)*rpor(ir)*tmp1(ir)
      enddo
      call valint(tmp2,rkll)

      do ir = 2,n
         tmp3(ir)=pf(ir,ka)*pf(ir,kc)*rpor(ir)*tmp4(ir)
      enddo
      call valint(tmp3,rkls)

      do i=1,MN
         tmp2(i)=zero
         tmp3(i)=zero
      enddo

      do ir = 2,n
         tmp2(ir)=qf(ir,ka)*qf(ir,kc)*rpor(ir)*tmp1(ir)
      enddo
      call valint(tmp2,rksl)

      do ir = 2,n
         tmp3(ir)=qf(ir,ka)*qf(ir,kc)*rpor(ir)*tmp4(ir)
      enddo
      call valint(tmp3,rkss)

      endif

      return
      end

      SUBROUTINE YZKHYBE(K,I,J)
*CCC
*CCC      $Id: YZK.f,v 1.4 1994/09/15 13:34:04 pn Exp $
*CCC
*CCC      $Log: YZK.f,v $
*CCC      Revision 1.4  1994/09/15 13:34:04  pn
*CCC      Version 2c  31 May 1994
*CCC      Modified version of 2b.
*CCC
*CCC Revision 1.3  1994/03/31  11:57:06  pn
*CCC Version 2b 6 July 1990
*CCC
*CCC
C***********************************************************************
C***********************************************************************
C-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)

      DIMENSION RHOP1(MN),RTTK(MN),WK1(MN),ZK1(MN)
      DIMENSION RHOP2(MN),WK2(MN),ZK2(MN)

C-----------------------------------------------------------------------

C Coefficients for 10 point integration formula as used in VALINT.
      D1=.26834148361926139703d0*h
      D2=7.13764630431297097963d0*h
      D3=1.77535941424830313719d0*h
      D4=-.81043570626903960237d0*h
      D5=4.54946288279621612954d0*h
      D6=-4.35155122655122655122d0*h

      if((n+10).ge.MN)then
      print*,'Dimension of grid size not sufficient'
      stop
      endif

C Initializing all the variables
      Do II=1,MN
      WK1(II) = ZERO
      ZK1(II) = ZERO
      WK2(II) = ZERO
      ZK2(II) = ZERO
      TMP1(II) = ZERO
      TMP2(II) = ZERO
      TMP3(II) = ZERO
      TMP4(II) = ZERO
      RHOP1(II) = ZERO
      RHOP2(II) = ZERO
      ENDDO

C Calculates r(i)**k, where K is the multipole
        DO 1 II = 2,N
          RTTK(II) = R(II)**K
   1    CONTINUE

C Calculates Integrand for Z(i)**k at different grids
      DO 2 II = 2,N
        RHOP1(II) = RP(II)*PF(II,I)*PF(II,J)
        RHOP2(II) = RP(II)*QF(II,I)*QF(II,J)
   2  CONTINUE

C Defines temporray functions as r(i)**k*Integrand
        DO 3  II = 2,N
          TMP2(II) = RTTK(II)*RHOP1(II)
          TMP3(II) = RTTK(II)*RHOP2(II)
   3    CONTINUE

C First 10 points of Z(i)**k are taken same as given below
      DO 4  II = 2,10

      ZK1(II) = D1*(TMP2(1)+TMP2(11))+D2*TMP2(6)
     :  +D3*(TMP2(2)+TMP2(10))+D4*(TMP2(3)+TMP2(9))
     :  +D5*(TMP2(4)+TMP2(8))+D6*(TMP2(5)+TMP2(7))

      ZK2(II) = D1*(TMP3(1)+TMP3(11))+D2*TMP3(6)
     :  +D3*(TMP3(2)+TMP3(10))+D4*(TMP3(3)+TMP3(9))
     :  +D5*(TMP3(4)+TMP3(8))+D6*(TMP3(5)+TMP3(7))

   4  CONTINUE

C Z(i)**k has been calculated using 10 integration points.
      DO 5  II = 11,N+10

      ZK1(II) = ZK1(II-10)+D1*(TMP2(II-10)+TMP2(II))+D2*TMP2(II-5)
     :  +D3*(TMP2(II-1)+TMP2(II-9))+D4*(TMP2(II-2)+TMP2(II-8))
     :  +D5*(TMP2(II-3)+TMP2(II-7))+D6*(TMP2(II-4)+TMP2(II-6))

      ZK2(II) = ZK2(II-10)+D1*(TMP3(II-10)+TMP3(II))+D2*TMP3(II-5)
     :  +D3*(TMP3(II-1)+TMP3(II-9))+D4*(TMP3(II-2)+TMP3(II-8))
     :  +D5*(TMP3(II-3)+TMP3(II-7))+D6*(TMP3(II-4)+TMP3(II-6))

   5  CONTINUE

C Z(i)**k has been defined as Z(i)**k*(r</r>)**k
        DO 6 II = 2,N
          ZK1(II) = ZK1(II)/RTTK(II)
          ZK2(II) = ZK2(II)/RTTK(II)
   6    CONTINUE

C R(i)**(K+1) has been defined
        DO 7 II = 2,N
          RTTK(II) = RTTK(II)*R(II)
   7    CONTINUE

C Calculates Integrand for W(i)**k/r(i)**(k+1) at different grids
        DO 8 II = 2,N
          TMP2(II) = RHOP1(II)/RTTK(II)
          TMP3(II) = RHOP2(II)/RTTK(II)
   8    CONTINUE

C Last 10 points of W(i)**k are taken same as given below
      DO 9 II = N,N+10

      WK1(II) = D1*(TMP2(N)+TMP2(N+10))+D2*TMP2(N+5)
     :  +D3*(TMP2(N+1)+TMP2(N+9))+D4*(TMP2(N+2)+TMP2(N+8))
     :  +D5*(TMP2(N+3)+TMP2(N+7))+D6*(TMP2(N+4)+TMP2(N+6))

      WK2(II) = D1*(TMP3(N)+TMP3(N+10))+D2*TMP3(N+5)
     :  +D3*(TMP3(N+1)+TMP3(N+9))+D4*(TMP3(N+2)+TMP3(N+8))
     :  +D5*(TMP3(N+3)+TMP3(N+7))+D6*(TMP3(N+4)+TMP3(N+6))

   9  CONTINUE

C W(i)**k has been calculated using 10 integration points.
      DO 10 II = N,2,-1

      WK1(II) = WK1(II+10)+D1*(TMP2(II)+TMP2(II+10))+D2*TMP2(II+5)
     :  +D3*(TMP2(II+1)+TMP2(II+9))+D4*(TMP2(II+2)+TMP2(II+8))
     :  +D5*(TMP2(II+3)+TMP2(II+7))+D6*(TMP2(II+4)+TMP2(II+6))

      WK2(II) = WK2(II+10)+D1*(TMP3(II)+TMP3(II+10))+D2*TMP3(II+5)
     :  +D3*(TMP3(II+1)+TMP3(II+9))+D4*(TMP3(II+2)+TMP3(II+8))
     :  +D5*(TMP3(II+3)+TMP3(II+7))+D6*(TMP3(II+4)+TMP3(II+6))

   10 CONTINUE

C W(i)**k has been defined as Z(i)**k*(r</r>)**(k+1)
        DO 11 II = 2,N
          WK1(II) = WK1(II)*RTTK(II)
          WK2(II) = WK2(II)*RTTK(II)
   11   CONTINUE

C YZK functions are determined as Z(i)**k+W(i)**k
      DO 12 II = 2,N
        TMP1(II) = ZK1(II) + WK1(II)
        TMP4(II) = ZK2(II) + WK2(II)
   12 CONTINUE

      RETURN
      END



      SUBROUTINE SETMATHYBVN1(RHO)
c ********************************************************************
C  This subroutine sets the diagonalization subroutine
c --------------------------------------------------------------------

      implicit real*8(a-h,o-z)

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     :      /jvalue/orbj(MNS),kap(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /basc/nbasis,nocc,mbasis
     :  /skip/nskipe(MNS),nskipc(MNS)
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     :  /quan1/kappe(MNBAS),orbje(MNBAS),kappc(MNBAS),orbjc(MNBAS)
     :  /quan2/lorba(MNS),lorbb(MNS),iqe(MNBAS),iqc(MNBAS)
     :      /atmp/tmp5(MN),tmp6(MN)
     :  /dfm/emt(MNB,MNB),df_single(MNB,MNB),df_two(MNB,MNB)
     :      /skipg/nbasg(MNS),nbasisg,mbasisg,nskipg(MNS)
     :      /npar/parm(2)

c ---------------------------------------------------------------
c Local dimension
c ---------------------------------------------------------------

      dimension rho(MN),mulk(MNS+1),atemp(MNB,MNB)
      DIMENSION PARM1(10),RG(MN),RPG(MN),RPORG(MN)
      DIMENSION EORBG(MNBAS),PZG(MNBAS),MFG(MNBAS)
      DIMENSION PFG(MN,MNBAS),QFG(MN,MNBAS)

      if(ihbd.eq.1)then

      OPEN(IGRASP,FILE='GRASP09.DAT',FORM='UNFORMATTED',STATUS='OLD')

      REWIND(IGRASP)
      READ(IGRASP)HG,NG
      READ(IGRASP)(rg(i),i = 1,ng),(rpg(i),i = 1,ng),(rporg(i),i=1,ng)
      do i=1,nbasis-nbasisg
      READ(IGRASP)EORBG(I)
      READ(IGRASP)(PFG(J,I),J=1,NG),(QFG(J,I),J=1,NG)
      ENDDO
      REWIND(IGRASP)

      IF(NG-N.NE.0)THEN
      Print*,'Mismatch in no. of grid points',NG,N
      STOP
      ENDIF

      IF(DABS(HG-H).GE.ACCY)THEN
      Print*,'Mismatch in stepsize',HG,H
      STOP
      ENDIF

      DO I=1,N
      IF(RG(I).NE.R(I))THEN
      Print*,'Mismatch in radial distribution',I,RG(I),R(I)
      STOP
      ENDIF
      ENDDO

      DO I=1,N
      IF(RPG(I).NE.RP(I))THEN
      Print*,'Mismatch in radial distribution',I,RPG(I),RP(I)
      STOP
      ENDIF
      ENDDO

      DO I=1,N
      IF(RPORG(I).NE.RPOR(I))THEN
      Print*,'Mismatch in radial distribution',I,RPORG(I),RPOR(I)
      STOP
      ENDIF
      ENDDO

      do isym=1,nsym
      nn=nbasg(isym)

      if(nn.ne.0)then
      do ia=1,nbasg(isym)
      ka=ia+nskipe(isym)

      ifact1=int(EORBG(ka)/dabs(EORBG(ka)))
      ifact2=int(eng1(ka)/dabs(eng1(ka)))
c
      IF(ifact1.NE.ifact2)THEN
      Print*,'Mismatch in sign of energies',ka,isym,ifact1,ifact2
      STOP
      ENDIF

      eng1(ka)=EORBG(ka)

      do ir=1,n

c      ifact1=int(pfg(ir,ka)/dabs(pfg(ir,ka)))
c      ifact2=int(pf(ir,ka)/dabs(pf(ir,ka)))
c
c      IF(ifact1.NE.ifact2)THEN
c      Print*,'Mismatch in sign of large component',ir,ka,ifact1,ifact2
c      STOP
c      ENDIF

c      ifact1=int(qfg(ir,ka)/dabs(qfg(ir,ka)))
c      ifact2=int(qf(ir,ka)/dabs(qf(ir,ka)))
c
c      IF(ifact1.NE.ifact2)THEN
c      Print*,'Mismatch in sign of small component',ir,ka,ifact1,ifact2
c      STOP
c      ENDIF

      pf(ir,ka)=pfg(ir,ka)
      qf(ir,ka)=qfg(ir,ka)

      enddo

      enddo
      endif
      enddo

c SCHIMIDT ORTHOGONALIZATION

         do isym=1,nsym
            do ia=nbasg(isym)+1,nbas(isym)
               ka=ia+nskipe(isym)

               do ib=1,ia-1,1
                 kb=ib+nskipe(isym)

                 r0=rint(ka,kb,0)

                 do ir=1,n
                   pf(ir,ka)=pf(ir,ka)-r0*pf(ir,kb)
                   qf(ir,ka)=qf(ir,ka)-r0*qf(ir,kb)
                 enddo

               enddo

                 r1=rint(ka,ka,0)

             if(dabs(one-r1).gt.ACCY)then
                do ir=1,n
                 pf(ir,ka)=pf(ir,ka)/r1
                 qf(ir,ka)=qf(ir,ka)/r1
                enddo
             endif

            enddo
         enddo
C Orthogonalization is over

      endif

      do i=1,MNB
      do j=1,MNB
      df_single(i,j)=zero
      df_two(i,j)=zero
      emt(i,j)=zero
      atemp(i,j)=zero
      enddo
      enddo

      do isym=1,nsym

      do ia=1+nbasg(isym),nbas(isym)
      ka=ia+nskipe(isym)
      ja=ia+2*nskipg(isym)-nbasg(isym)
      index1=ia+nbas(isym)-2*nbasg(isym)+2*nskipg(isym)

      kapa=kappe(ka)

      call dpbdt(ka)

      do ib=1+nbasg(isym),nbas(isym)
      kb=ib+nskipe(isym)
      jb=ib+2*nskipg(isym)-nbasg(isym)
      index2=ib+nbas(isym)-2*nbasg(isym)+2*nskipg(isym)

      kapb=kappe(kb)

      do la = 1,n
      tmp1(la)=pf(la,ka)*pf(la,kb)*rp(la)
      tmp2(la)=qf(la,ka)*qf(la,kb)*rp(la)
      tmp3(la)=qf(la,kb)*(kapa*pf(la,ka)*rpor(la)+tmp5(la)/h)
      enddo

      call valint(tmp1,result1)
      call valint(tmp2,result2)
      call valint(tmp3,result3)

      emt(ja,jb)=result1
      emt(index1,index2)=result2
      atemp(ka,kb)=c*result3

      enddo
      enddo

c ---------------------------------------
c     matrix df_single
c --------------------------------------
      do ia=1+nbasg(isym),nbas(isym)
      ka=ia+nskipe(isym)
      ja=ia+2*nskipg(isym)-nbasg(isym)
      index1=ia+nbas(isym)-2*nbasg(isym)+2*nskipg(isym)

      orba=orbje(ka)
      iaa=iqe(ka)
      kapa=kappe(ka)

      do ib=1+nbasg(isym),nbas(isym)
      kb=ib+nskipe(isym)
      jb=ib+2*nskipg(isym)-nbasg(isym)
      index2=ib+nbas(isym)-2*nbasg(isym)+2*nskipg(isym)

      orbb=orbje(kb)
      iab=iqe(kb)
      kapb=kappe(kb)

      df_single(ja,index2)=df_single(ja,index2)+atemp(ka,kb)
      df_single(index2,ja)=df_single(ja,index2)

c -----------------------------------------
c     nuclear repulsion energy
c ----------------------------------------

      do la = 1,n
      tmp1(la)=pf(la,ka)*pf(la,kb)*rpor(la)*rho(la)
      tmp2(la)=qf(la,ka)*qf(la,kb)*rpor(la)*rho(la)
      enddo

      call valint(tmp1,result1)
      call valint(tmp2,result2)

      df_single(ja,jb)=df_single(ja,jb)-result1
      df_single(index1,index2)=df_single(index1,index2)
     :   -(result2+two*c*c*emt(index1,index2))

      do jsym=1,nsym
      nn=nocorb(jsym)

      if(nn.ne.0)then

      do ic=1,nocorb(jsym)
      jc=ic+nskipc(jsym)
      kc=ic+nskipe(jsym)

      kapc=kappc(jc)
      orbc=orbjc(jc)
      iac=iqc(jc)

c-------------------------------------------------------------------
c     direct part
c-------------------------------------------------------------------

      fact=(two*orbc+one)

      call findk(orba,orbb,orbc,orbc,mulk,maxk)

      if(maxk.ne.0) then
      do ik=1,maxk
      k1=mulk(ik)

      if(k1.eq.0)then

      rkll=zero
      rkss=zero

      call slathybd(ka,kc,kb,kc,orba,orbc,orbb,orbc,
     :iaa,iac,iab,iac,k1,rkll,rkss)

      df_two(ja,jb)=df_two(ja,jb)+rkll*fact

      df_two(index1,index2)=df_two(index1,index2)+rkss*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     exchange part
c-------------------------------------------------------------------

      call findk(orba,orbc,orbc,orbb,mulk,maxk)

      if(maxk.ne.0)then
      do ik=1,maxk
      k1=mulk(ik)
      ak1=dfloat(k1)

      fact=(two*orbc+one)*clrx(kapa,k1,kapc)**2

      rkll=zero
      rkls=zero
      rksl=zero
      rkss=zero

      bkll=zero
      bkls=zero
      bksl=zero
      bkss=zero

      call slathybe(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,k1,rkll,rkls,rksl,rkss)

      if(ibrt.gt.0)then
      call slatvnb(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,k1,bkll,bkls,bksl,bkss)
      endif

      df_two(ja,jb)=df_two(ja,jb)-(rkll+bkll)*fact

      df_two(ja,index2)=df_two(ja,index2)-(rkls+bkls)*fact

      df_two(index1,jb)=df_two(index1,jb)-(rksl+bksl)*fact

      df_two(index1,index2)=df_two(index1,index2)-(rkss+bkss)*fact

      enddo
      endif

      enddo
      endif
      enddo

      enddo
      enddo

      enddo

      RETURN
      END


      subroutine scfhybvn1
      implicit real*8(a-h,o-z)

      logical set
      character*3 nh_i

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     :      /jvalue/orbj(MNS),kap(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /basc/nbasis,nocc,mbasis
     :  /skip/nskipe(MNS),nskipc(MNS)
     :  /conv/crit,maxit,npower
     :/wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     :  /dfm/emt(MNB,MNB),df_single(MNB,MNB),df_two(MNB,MNB)
     :   /skipg/nbasg(MNS),nbasisg,mbasisg,nskipg(MNS)

c------------------------------------------------------------------
c Local dimension
c------------------------------------------------------------------
      dimension eng(MNB),scr(MNB,MNB),vec(MNB,MNB)
      dimension nh_i(MNS),a(MNB,MNB),oe(MNB,MNB)
      dimension pfg(MN,MNBAS),qfg(MN,MNBAS),stemp(MNB,MNB)
      dimension drs(MNB,MNB),can(MNB,MNB),eigv(MNB,MNB)
      dimension st(MNB,MNB),b(MNB,MNB)
      dimension df_mat(MNB,MNB),sa(MNB)

      data nh_i/' S ',' P-',' P ',' D-',' D ',' F-',' F ',
     :   ' G-',' G ',' H-',' H ',' I-',' I ',' K-',' K '/
C     NOTE: No J in this sequence
c--------------------------------------------------------------------
c       initialize arrays and constants required for the computation.
c--------------------------------------------------------------------

      set = .true.

      do ia=1,MNB
      sa(ia)=zero
      do ib=1,MNB
      can(ib,ia)=zero
      drs(ib,ia)=zero
      stemp(ib,ia)=zero
      enddo
      enddo

      do isym=1,nsym
      ndim=nbas(isym)-nbasg(isym)
      ndim2=2*ndim

      do ia=1,ndim2
      ja=ia+2*nskipg(isym)
      do ib=1,ndim2
      jb=ib+2*nskipg(isym)
      stemp(ia,ib)=emt(ja,jb)
      enddo
      enddo

c      if(jpnt.ne.0)call matout(stemp,MNB,ndim2,ndim2,ndim2,4hS   )

      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=stemp(i,j)
      enddo
      enddo

      call eigen1(drs,scr,ndim2,sa,ndim2,bnorm,MNB)
C      call eigen(drs,scr,ndim2,sa,ndim2,bnorm,MNB)
c-------------------------------------------------------------------
c       Canonical transformation
c------------------------------------------------------------------
      do ia=1,ndim2
      do ib=1,ndim2
      if(dabs(sa(ib)).gt.TINY)then
      can(ia,ib)=drs(ia,ib)/dsqrt(dabs(sa(ib)))
      st(ib,ia)=can(ia,ib)
      else
      write(*,*)'There is dioginalisation problem',sa(ib)
       stop
      endif
      enddo
      enddo

c      if(jpnt.ne.0)call matout(st,MNB,ndim2,ndim2,ndim2,4hXT  )

c      if(jpnt.ne.0)
c     :        call matout(df_single,MNB,mbasis,mbasis,mbasis,4hH   )

      do ic=1,ndim2
      jc=ic+2*nskipg(isym)
      do id=1,ndim2
      jd=id+2*nskipg(isym)
      a(ic,id)=df_single(jc,jd)+df_two(jc,jd)
      enddo
      enddo

c      if(jpnt.ne.0)call matout(a,MNB,ndim2,ndim2,ndim2,4hFOCK )

c      if(jpnt.ne.0)
c     :call matout(df_two,MNB,mbasis,mbasis,mbasis,4hVHF )

      call mult(a,can,b,MNB,ndim2)
      call mult(st,b,df_mat,MNB,ndim2)

c      if(jpnt.ne.0)
c     :      call matout(df_mat,MNB,ndim2,ndim2,ndim2,6hFRIME )

      call eigen1(df_mat,scr,ndim2,eng,ndim2,bnorm,MNB)
C      call eigen(df_mat,scr,ndim2,eng,ndim2,bnorm,MNB)

      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=df_mat(i,j)
      enddo
      enddo

c      if(jpnt.ne.0)
c     : call matout(drs,MNB,ndim2,ndim2,ndim2,4hC   )

      call mult(can,drs,eigv,MNB,ndim2)

      loc = 1
      last=ndim-1
      dowhile(loc.le.last)
      samin = eng(loc)
      locus = loc

      do  ia = loc,ndim
      if(eng(ia).lt.samin)then
      samin = eng(ia)
      locus = ia
      endif
      enddo

      if(locus.gt.loc)then
      swap = eng(loc)
      eng(loc) = eng(locus)
      eng(locus) = swap
      do  ib = 1,ndim2
      swap = eigv(ib,loc)
      eigv(ib,loc) = eigv(ib,locus)
      eigv(ib,locus) = swap
      enddo
      endif
      loc = loc + 1
      enddo

      do ia=1,ndim2
      do ib=1,ndim2
      b(ia,ib)=zero
      enddo
      b(ia,ia)=eng(ia)
      oe(ia,isym)=eng(ia)
      enddo

c      if(jpnt.ne.0)call matout(b,MNB,ndim2,ndim2,ndim2,4hE   )

c      if(jpnt.ne.0)call matout(eigv,MNB,ndim2,ndim2,ndim2,4hCP  )

      do ia=1,ndim2
      ja=ia+2*nskipg(isym)

      do ib=1,ndim2
      jb=ib+2*nskipg(isym)

      vec(ja,jb)=eigv(ia,ib)

      enddo
      enddo

c-----------------------------------------------------------------
c   Initialisation
c-----------------------------------------------------------------

       do ia=1+nbasg(isym),nbas(isym)
       ka=ia+nskipe(isym)
       ja=ia+2*nskipg(isym)-nbasg(isym)

       do ic=1,n
       pfg(ic,ka)=zero
       qfg(ic,ka)=zero
       enddo

       do ib=1+nbasg(isym),nbas(isym)
       kb=ib+nskipe(isym)
       jb=ib+2*nskipg(isym)-nbasg(isym)

       do ic=1,n
       pfg(ic,ka)=pfg(ic,ka)+vec(jb,ja)*pf(ic,kb)
       qfg(ic,ka)=qfg(ic,ka)+
     :vec(jb+nbas(isym)-nbasg(isym),ja)*qf(ic,kb)
       enddo

       enddo

       eng1(ka)=oe(ia-nbasg(isym),isym)

       enddo

       do ia=1+nbasg(isym),nbas(isym)
       ka=ia+nskipe(isym)

       do ic=1,n
       pf(ic,ka)=pfg(ic,ka)
       qf(ic,ka)=qfg(ic,ka)
       enddo

       enddo

       do ibas=1,nbas(isym)
       ka=ibas+nskipe(isym)
       write(STDOUT,206)ibas,nh_i(isym),eng1(ka)
       enddo

      enddo


 206  format(/,4x,I2,A,D20.11)

      RETURN
      END

      SUBROUTINE SETMATHYBVN(RHO)
c ********************************************************************
C  This subroutine sets the diagonalization subroutine
c --------------------------------------------------------------------

      implicit real*8(a-h,o-z)

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     :      /jvalue/orbj(MNS),kap(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /basc/nbasis,nocc,mbasis
     :  /skip/nskipe(MNS),nskipc(MNS)
     :  /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     :  /quan1/kappe(MNBAS),orbje(MNBAS),kappc(MNBAS),orbjc(MNBAS)
     :  /quan2/lorba(MNS),lorbb(MNS),iqe(MNBAS),iqc(MNBAS)
     :      /atmp/tmp5(MN),tmp6(MN)
     :  /dfm/emt(MNB,MNB),df_single(MNB,MNB),df_two(MNB,MNB)
     :      /skipg/nbasg(MNS),nbasisg,mbasisg,nskipg(MNS)
     :      /npar/parm(2)

c ---------------------------------------------------------------
c Local dimension
c ---------------------------------------------------------------

      dimension rho(MN),mulk(MNS+1),atemp(MNB,MNB)
      DIMENSION PARM1(10),RG(MN),RPG(MN),RPORG(MN)
      DIMENSION EORBG(MNBAS),PZG(MNBAS),MFG(MNBAS)
      DIMENSION PFG(MN,MNBAS),QFG(MN,MNBAS)

      if(ihbd.eq.1)then

      OPEN(IGRASP,FILE='GRASP09.DAT',FORM='UNFORMATTED',STATUS='OLD')

      REWIND(IGRASP)
      READ(IGRASP)HG,NG
      READ(IGRASP)(rg(i),i = 1,ng),(rpg(i),i = 1,ng),(rporg(i),i=1,ng)
      do i=1,nbasis
      READ(IGRASP)EORBG(I)
      READ(IGRASP)(PFG(J,I),J=1,NG),(QFG(J,I),J=1,NG)
      ENDDO
      REWIND(IGRASP)

      IF(NG-N.NE.0)THEN
      Print*,'Mismatch in no. of grid points',NG,N
      STOP
      ENDIF

      IF(DABS(HG-H).GE.ACCY)THEN
      Print*,'Mismatch in stepsize',HG,H
      STOP
      ENDIF

      DO I=1,N
      IF(RG(I).NE.R(I))THEN
      Print*,'Mismatch in radial distribution',I,RG(I),R(I)
      STOP
      ENDIF
      ENDDO

      DO I=1,N
      IF(RPG(I).NE.RP(I))THEN
      Print*,'Mismatch in radial distribution',I,RPG(I),RP(I)
      STOP
      ENDIF
      ENDDO

      DO I=1,N
      IF(RPORG(I).NE.RPOR(I))THEN
      Print*,'Mismatch in radial distribution',I,RPORG(I),RPOR(I)
      STOP
      ENDIF
      ENDDO

      do isym=1,nsym
      nn=nbasg(isym)

      if(nn.ne.0)then
      do ia=1,nbasg(isym)
      ka=ia+nskipe(isym)

      ifact1=int(EORBG(ka)/dabs(EORBG(ka)))
      ifact2=int(eng1(ka)/dabs(eng1(ka)))
c
      IF(ifact1.NE.ifact2)THEN
      Print*,'Mismatch in sign of energies',ka,isym,ifact1,ifact2
      STOP
      ENDIF

      eng1(ka)=EORBG(ka)

      do ir=1,n

c      ifact1=int(pfg(ir,ka)/dabs(pfg(ir,ka)))
c      ifact2=int(pf(ir,ka)/dabs(pf(ir,ka)))
c
c      IF(ifact1.NE.ifact2)THEN
c      Print*,'Mismatch in sign of large component',ir,ka,ifact1,ifact2
c      STOP
c      ENDIF

c      ifact1=int(qfg(ir,ka)/dabs(qfg(ir,ka)))
c      ifact2=int(qf(ir,ka)/dabs(qf(ir,ka)))
c
c      IF(ifact1.NE.ifact2)THEN
c      Print*,'Mismatch in sign of small component',ir,ka,ifact1,ifact2
c      STOP
c      ENDIF

      pf(ir,ka)=pfg(ir,ka)
      qf(ir,ka)=qfg(ir,ka)

      enddo

      enddo
      endif
      enddo

c SCHIMIDT ORTHOGONALIZATION

         do isym=1,nsym
            do ia=nbasg(isym)+1,nbas(isym)
               ka=ia+nskipe(isym)

               do ib=1,ia-1,1
                 kb=ib+nskipe(isym)

                 r0=rint(ka,kb,0)

                 do ir=1,n
                   pf(ir,ka)=pf(ir,ka)-r0*pf(ir,kb)
                   qf(ir,ka)=qf(ir,ka)-r0*qf(ir,kb)
                 enddo

               enddo

                 r1=rint(ka,ka,0)

             if(dabs(one-r1).gt.ACCY)then
                do ir=1,n
                 pf(ir,ka)=pf(ir,ka)/r1
                 qf(ir,ka)=qf(ir,ka)/r1
                enddo
             endif

            enddo
         enddo
C Orthogonalization is over

      endif

      do i=1,MNB
      do j=1,MNB
      df_single(i,j)=zero
      df_two(i,j)=zero
      emt(i,j)=zero
      atemp(i,j)=zero
      enddo
      enddo

      do isym=1,nsym

      do ia=1,nbas(isym)
      ka=ia+nskipe(isym)
      ja=ia+2*nskipe(isym)
      index1=ia+nbas(isym)+2*nskipe(isym)

      kapa=kappe(ka)

      call dpbdt(ka)

      do ib=1,nbas(isym)
      kb=ib+nskipe(isym)
      jb=ib+2*nskipe(isym)
      index2=ib+nbas(isym)+2*nskipe(isym)

      kapb=kappe(kb)

      do la = 1,n
      tmp1(la)=pf(la,ka)*pf(la,kb)*rp(la)
      tmp2(la)=qf(la,ka)*qf(la,kb)*rp(la)
      tmp3(la)=qf(la,kb)*(kapa*pf(la,ka)*rpor(la)+tmp5(la)/h)
      enddo

      call valint(tmp1,result1)
      call valint(tmp2,result2)
      call valint(tmp3,result3)

      emt(ja,jb)=result1
      emt(index1,index2)=result2
      atemp(ka,kb)=c*result3

      enddo
      enddo

c ---------------------------------------
c     matrix df_single
c --------------------------------------
      do ia=1,nbas(isym)
      ka=ia+nskipe(isym)
      ja=ia+2*nskipe(isym)
      index1=ia+nbas(isym)+2*nskipe(isym)

      orba=orbje(ka)
      iaa=iqe(ka)
      kapa=kappe(ka)

      do ib=1,nbas(isym)
      kb=ib+nskipe(isym)
      jb=ib+2*nskipe(isym)
      index2=ib+nbas(isym)+2*nskipe(isym)

      orbb=orbje(kb)
      iab=iqe(kb)
      kapb=kappe(kb)

      df_single(ja,index2)=df_single(ja,index2)+atemp(ka,kb)
      df_single(index2,ja)=df_single(ja,index2)

c -----------------------------------------
c     nuclear repulsion energy
c ----------------------------------------

      do la = 1,n
      tmp1(la)=pf(la,ka)*pf(la,kb)*rpor(la)*rho(la)
      tmp2(la)=qf(la,ka)*qf(la,kb)*rpor(la)*rho(la)
      enddo

      call valint(tmp1,result1)
      call valint(tmp2,result2)

      df_single(ja,jb)=df_single(ja,jb)-result1
      df_single(index1,index2)=df_single(index1,index2)
     :   -(result2+two*c*c*emt(index1,index2))

      do jsym=1,nsym
      nn=nocorb(jsym)

      if(nn.ne.0)then

      do ic=1,nocorb(jsym)
      jc=ic+nskipc(jsym)
      kc=ic+nskipe(jsym)

      kapc=kappc(jc)
      orbc=orbjc(jc)
      iac=iqc(jc)

c-------------------------------------------------------------------
c     direct part
c-------------------------------------------------------------------

      rkll=zero
      rkss=zero

      fact=(two*orbc+one)

      call findk(orba,orbb,orbc,orbc,mulk,maxk)

      if(maxk.ne.0) then
      do ik=1,maxk
      k1=mulk(ik)

      if(k1.eq.0)then

      call slathybd(ka,kc,kb,kc,orba,orbc,orbb,orbc,
     :iaa,iac,iab,iac,k1,rkll,rkss)

      df_two(ja,jb)=df_two(ja,jb)+rkll*fact

      df_two(index1,index2)=df_two(index1,index2)+rkss*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     exchange part
c-------------------------------------------------------------------

      rkll=zero
      rkls=zero
      rksl=zero
      rkss=zero

      bkll=zero
      bkls=zero
      bksl=zero
      bkss=zero

      call findk(orba,orbc,orbc,orbb,mulk,maxk)

      if(maxk.ne.0)then
      do ik=1,maxk
      k1=mulk(ik)
      ak1=dfloat(k1)

      fact=(two*orbc+one)*clrx(kapa,k1,kapc)**2

      call slathybe(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,k1,rkll,rkls,rksl,rkss)

      if(ibrt.gt.0)then
      call slatvnb(ka,kc,kc,kb,orba,orbc,orbc,orbb,
     :iaa,iac,iac,iab,k1,bkll,bkls,bksl,bkss)
      endif

      df_two(ja,jb)=df_two(ja,jb)-(rkll+bkll)*fact

      df_two(ja,index2)=df_two(ja,index2)-(rkls+bkls)*fact

      df_two(index1,jb)=df_two(index1,jb)-(rksl+bksl)*fact

      df_two(index1,index2)=df_two(index1,index2)-(rkss+bkss)*fact

      enddo
      endif

      enddo
      endif
      enddo

      enddo
      enddo

      enddo

      RETURN
      END


      subroutine scfhybvn
      implicit real*8(a-h,o-z)

      logical set
      character*3 nh_i

      common/datt/nsym,nbas(MNS),nocorb(MNS),nvorb(MNS)
     :      /jvalue/orbj(MNS),kap(MNS)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /basc/nbasis,nocc,mbasis
     :  /skip/nskipe(MNS),nskipc(MNS)
     :  /conv/crit,maxit,npower
     :/wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     :  /dfm/emt(MNB,MNB),df_single(MNB,MNB),df_two(MNB,MNB)
     :   /skipg/nbasg(MNS),nbasisg,mbasisg,nskipg(MNS)

c------------------------------------------------------------------
c Local dimension
c------------------------------------------------------------------
      dimension eng(MNB),scr(MNB,MNB),vec(MNB,MNB)
      dimension nh_i(MNS),a(MNB,MNB),oe(MNB,MNB)
      dimension pfg(MN,MNBAS),qfg(MN,MNBAS),stemp(MNB,MNB)
      dimension drs(MNB,MNB),can(MNB,MNB),eigv(MNB,MNB)
      dimension st(MNB,MNB),b(MNB,MNB)
      dimension df_mat(MNB,MNB),sa(MNB)

      data nh_i/' S ',' P-',' P ',' D-',' D ',' F-',' F ',
     :   ' G-',' G ',' H-',' H ',' I-',' I ',' K-',' K '/

C     Note : No J in this sequence
c--------------------------------------------------------------------
c       initialize arrays and constants required for the computation.
c--------------------------------------------------------------------

      set = .true.

      do ia=1,MNB
      sa(ia)=zero
      do ib=1,MNB
      can(ib,ia)=zero
      drs(ib,ia)=zero
      stemp(ib,ia)=zero
      enddo
      enddo

      do isym=1,nsym
      ndim=nbas(isym)
      ndim2=2*ndim

      do ia=1,ndim2
      ja=ia+2*nskipe(isym)
      do ib=1,ndim2
      jb=ib+2*nskipe(isym)
      stemp(ia,ib)=emt(ja,jb)
      enddo
      enddo

c      if(jpnt.ne.0)call matout(stemp,MNB,ndim2,ndim2,ndim2,4hS   )

      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=stemp(i,j)
      enddo
      enddo

      call eigen1(drs,scr,ndim2,sa,ndim2,bnorm,MNB)
C      call eigen(drs,scr,ndim2,sa,ndim2,bnorm,MNB)
c-------------------------------------------------------------------
c       Canonical transformation
c------------------------------------------------------------------
      do ia=1,ndim2
      do ib=1,ndim2
      if(dabs(sa(ib)).gt.TINY)then
      can(ia,ib)=drs(ia,ib)/dsqrt(dabs(sa(ib)))
      st(ib,ia)=can(ia,ib)
      else
      write(*,*)'There is dioginalisation problem',sa(ib)
       stop
      endif
      enddo
      enddo

c      if(jpnt.ne.0)call matout(st,MNB,ndim2,ndim2,ndim2,4hXT  )

c      if(jpnt.ne.0)
c     :        call matout(df_single,MNB,mbasis,mbasis,mbasis,4hH   )

      do ic=1,ndim2
      jc=ic+2*nskipe(isym)
      do id=1,ndim2
      jd=id+2*nskipe(isym)
      a(ic,id)=df_single(jc,jd)+df_two(jc,jd)
      enddo
      enddo

c      if(jpnt.ne.0)call matout(a,MNB,ndim2,ndim2,ndim2,4hFOCK )

c      if(jpnt.ne.0)
c     :call matout(df_two,MNB,mbasis,mbasis,mbasis,4hVHF )

      call mult(a,can,b,MNB,ndim2)
      call mult(st,b,df_mat,MNB,ndim2)

c      if(jpnt.ne.0)
c     :      call matout(df_mat,MNB,ndim2,ndim2,ndim2,6hFRIME )

      call eigen1(df_mat,scr,ndim2,eng,ndim2,bnorm,MNB)
C      call eigen(df_mat,scr,ndim2,eng,ndim2,bnorm,MNB)

      do i=1,ndim2
      do j=1,ndim2
      drs(i,j)=df_mat(i,j)
      enddo
      enddo

c      if(jpnt.ne.0)
c     : call matout(drs,MNB,ndim2,ndim2,ndim2,4hC   )

      call mult(can,drs,eigv,MNB,ndim2)

      loc = 1
      last=ndim-1
      dowhile(loc.le.last)
      samin = eng(loc)
      locus = loc

      do  ia = loc,ndim
      if(eng(ia).lt.samin)then
      samin = eng(ia)
      locus = ia
      endif
      enddo

      if(locus.gt.loc)then
      swap = eng(loc)
      eng(loc) = eng(locus)
      eng(locus) = swap
      do  ib = 1,ndim2
      swap = eigv(ib,loc)
      eigv(ib,loc) = eigv(ib,locus)
      eigv(ib,locus) = swap
      enddo
      endif
      loc = loc + 1
      enddo

      do ia=1,ndim2
      do ib=1,ndim2
      b(ia,ib)=zero
      enddo
      b(ia,ia)=eng(ia)
      oe(ia,isym)=eng(ia)
      enddo

c      if(jpnt.ne.0)call matout(b,MNB,ndim2,ndim2,ndim2,4hE   )

c      if(jpnt.ne.0)call matout(eigv,MNB,ndim2,ndim2,ndim2,4hCP  )

      do ia=1,ndim2
      ja=ia+2*nskipe(isym)

      do ib=1,ndim2
      jb=ib+2*nskipe(isym)

      vec(ja,jb)=eigv(ia,ib)

      enddo
      enddo

c-----------------------------------------------------------------
c   Initialisation
c-----------------------------------------------------------------

       do ia=1+nbasg(isym),nbas(isym)
       ka=ia+nskipe(isym)
       ja=ia+2*nskipe(isym)

       do ic=1,n
       pfg(ic,ka)=zero
       qfg(ic,ka)=zero
       enddo

       do ib=1+nbasg(isym),nbas(isym)
       kb=ib+nskipe(isym)
       jb=ib+2*nskipe(isym)

       do ic=1,n
       pfg(ic,ka)=pfg(ic,ka)+vec(jb,ja)*pf(ic,kb)
       qfg(ic,ka)=qfg(ic,ka)+vec(jb+nbas(isym),ja)*qf(ic,kb)
       enddo

       enddo

       eng1(ka)=oe(ia,isym)

       enddo

       do ia=1+nbasg(isym),nbas(isym)
       ka=ia+nskipe(isym)

       do ic=1,n
       pf(ic,ka)=pfg(ic,ka)
       qf(ic,ka)=qfg(ic,ka)
       enddo

       enddo

       do ibas=1,nbas(isym)
       ka=ibas+nskipe(isym)
       write(STDOUT,206)ibas,nh_i(isym),eng1(ka)
       enddo

      enddo

 206  format(/,4x,I2,A,D20.11)

      RETURN
      END

      subroutine slatvnb(ka,kb,kc,kd,orba,orbb,orbc,orbd,
     :iaa,iab,iac,iad,k1,bkll,bkls,bksl,bkss)
      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     : /iptn/isgo,jpnt,ibrt,ispd,ieng,inuc,idflt,ichrg,ivn,ihbd,intg
     :      /atmp/tmp5(MN),tmp6(MN)
     :  /quan1/kappe(MNBAS),orbje(MNBAS),kappc(MNBAS),orbjc(MNBAS)

      DIMENSION TMP(MN)

      do i=1,MN
         tmp(i)=zero
         tmp1(i)=zero
         tmp2(i)=zero
         tmp3(i)=zero
         tmp4(i)=zero
         tmp5(i)=zero
         tmp6(i)=zero
      enddo

      bkll=zero
      bkls=zero
      bksl=zero
      bkss=zero
      fact=zero

      kapa=-idint((orba+half)*iaa)
      kapb=-idint((orbb+half)*iab)
      kapc=-idint((orbc+half)*iac)
      kapd=-idint((orbd+half)*iad)

      if(kapa.ne.kappe(ka).or.kapb.ne.kappe(kb).or.
     :kapc.ne.kappe(kc).or.kapd.ne.kappe(kd))then
      Print*,'Quantum numbers are wrongly assigned'
      Stop
      Endif

      if(k1.gt.0)then

      Do L1=k1+1,k1+3

      L=L1-2

      fact=selectb(iaa,iac,orba,orbc,L)*
     : selectb(iab,iad,orbb,orbd,L)

      if(fact.ne.zero)then

      call yzkvnbrt(L,kb,kd)

      do ir = 2,n
         tmp(ir)=pf(ir,ka)*qf(ir,kc)*rpor(ir)*tmp1(ir)
      enddo
      call valint(tmp,result)

      bkll=bkll+result*GK(2,k1,L,kapa,kapb,kapc,kapd)

      do ir = 2,n
         tmp(ir)=pf(ir,ka)*qf(ir,kc)*rpor(ir)*tmp2(ir)
      enddo
      call valint(tmp,result)

      bkls=bkls+result*GK(1,k1,L,kapa,kapb,kapc,kapd)

      do ir = 2,n
         tmp(ir)=qf(ir,ka)*pf(ir,kc)*rpor(ir)*tmp1(ir)
      enddo
      call valint(tmp,result)

      bksl=bksl+result*GK(4,k1,L,kapa,kapb,kapc,kapd)

      do ir = 2,n
         tmp(ir)=qf(ir,ka)*pf(ir,kc)*rpor(ir)*tmp2(ir)
      enddo
      call valint(tmp,result)

      bkss=bkss+result*GK(3,k1,L,kapa,kapb,kapc,kapd)

      If(ibrt.gt.1)then

      do ir = 2,n
         tmp(ir)=pf(ir,ka)*qf(ir,kc)*rpor(ir)*tmp3(ir)
      enddo
      call valint(tmp,result)

      bkll=bkll+result*RA(2,k1,L,kapa,kapb,kapc,kapd)

      do ir = 2,n
         tmp(ir)=pf(ir,ka)*qf(ir,kc)*rpor(ir)*tmp4(ir)
      enddo
      call valint(tmp,result)

      bkls=bkls+result*RA(1,k1,L,kapa,kapb,kapc,kapd)

      do ir = 2,n
         tmp(ir)=qf(ir,ka)*pf(ir,kc)*rpor(ir)*tmp3(ir)
      enddo
      call valint(tmp,result)

      bksl=bksl+result*RA(4,k1,L,kapa,kapb,kapc,kapd)

      do ir = 2,n
         tmp(ir)=qf(ir,ka)*pf(ir,kc)*rpor(ir)*tmp4(ir)
      enddo
      call valint(tmp,result)

      bkss=bkss+result*RA(3,k1,L,kapa,kapb,kapc,kapd)

      do ir = 2,n
         tmp(ir)=pf(ir,ka)*qf(ir,kc)*rpor(ir)*tmp5(ir)
      enddo
      call valint(tmp,result)

      bkll=bkll+result*RB(2,k1,L,kapa,kapb,kapc,kapd)

      do ir = 2,n
         tmp(ir)=pf(ir,ka)*qf(ir,kc)*rpor(ir)*tmp6(ir)
      enddo
      call valint(tmp,result)

      bkls=bkls+result*RB(1,k1,L,kapa,kapb,kapc,kapd)

      do ir = 2,n
         tmp(ir)=qf(ir,ka)*pf(ir,kc)*rpor(ir)*tmp5(ir)
      enddo
      call valint(tmp,result)

      bksl=bksl+result*RB(4,k1,L,kapa,kapb,kapc,kapd)

      do ir = 2,n
         tmp(ir)=qf(ir,ka)*pf(ir,kc)*rpor(ir)*tmp6(ir)
      enddo
      call valint(tmp,result)

      bkss=bkss+result*RB(3,k1,L,kapa,kapb,kapc,kapd)

      EndIf

      endif

      EndDo

      endif

      return
      end


      SUBROUTINE YZKVNBRT(K,I,J)
*CCC
*CCC      $Id: YZK.f,v 1.4 1994/09/15 13:34:04 pn Exp $
*CCC
*CCC      $Log: YZK.f,v $
*CCC      Revision 1.4  1994/09/15 13:34:04  pn
*CCC      Version 2c  31 May 1994
*CCC      Modified version of 2b.
*CCC
*CCC Revision 1.3  1994/03/31  11:57:06  pn
*CCC Version 2b 6 July 1990
*CCC
*CCC
C***********************************************************************
C***********************************************************************
C-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)

      Common/grid/r(MN),rp(MN),rpor(MN),rnt,h,n
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /dumm/tmp1(MN),tmp2(MN),tmp3(MN),tmp4(MN)
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :      /wave/pf(MN,MNBAS),qf(MN,MNBAS),eng1(MNB)
     :      /atmp/tmp5(MN),tmp6(MN)

      DIMENSION RHOP1(MN),RTTK(MN),WK1(MN),ZK1(MN)
      DIMENSION RHOP2(MN),WK2(MN),ZK2(MN)

C-----------------------------------------------------------------------

C Coefficients for 10 point integration formula as used in VALINT.
      D1=.26834148361926139703d0*h
      D2=7.13764630431297097963d0*h
      D3=1.77535941424830313719d0*h
      D4=-.81043570626903960237d0*h
      D5=4.54946288279621612954d0*h
      D6=-4.35155122655122655122d0*h

      if((n+10).ge.MN)then
      print*,'Dimension of grid size not sufficient'
      stop
      endif

C Initializing all the variables
      Do II=1,MN
      WK1(II) = ZERO
      ZK1(II) = ZERO
      WK2(II) = ZERO
      ZK2(II) = ZERO
      TMP1(II) = ZERO
      TMP2(II) = ZERO
      TMP3(II) = ZERO
      TMP4(II) = ZERO
      TMP5(II) = ZERO
      TMP6(II) = ZERO
      RHOP1(II) = ZERO
      RHOP2(II) = ZERO
      ENDDO

C Calculates r(i)**k, where K is the multipole
        DO 1 II = 2,N
          RTTK(II) = R(II)**K
   1    CONTINUE

C Calculates Integrand for Z(i)**k at different grids
      DO 2 II = 2,N
        RHOP1(II) = RP(II)*QF(II,I)*PF(II,J)
        RHOP2(II) = RP(II)*PF(II,I)*QF(II,J)
   2  CONTINUE

C Defines temporray functions as r(i)**k*Integrand
        DO 3  II = 2,N
          TMP5(II) = RTTK(II)*RHOP1(II)
          TMP6(II) = RTTK(II)*RHOP2(II)
   3    CONTINUE

C First 10 points of Z(i)**k are taken same as given below
      DO 4  II = 2,10

      ZK1(II) = D1*(TMP5(1)+TMP5(11))+D2*TMP5(6)
     :  +D3*(TMP5(2)+TMP5(10))+D4*(TMP5(3)+TMP5(9))
     :  +D5*(TMP5(4)+TMP5(8))+D6*(TMP5(5)+TMP5(7))

      ZK2(II) = D1*(TMP6(1)+TMP6(11))+D2*TMP6(6)
     :  +D3*(TMP6(2)+TMP6(10))+D4*(TMP6(3)+TMP6(9))
     :  +D5*(TMP6(4)+TMP6(8))+D6*(TMP6(5)+TMP6(7))

   4  CONTINUE

C Z(i)**k has been calculated using 10 integration points.
      DO 5  II = 11,N+10

      ZK1(II) = ZK1(II-10)+D1*(TMP5(II-10)+TMP5(II))+D2*TMP5(II-5)
     :  +D3*(TMP5(II-1)+TMP5(II-9))+D4*(TMP5(II-2)+TMP5(II-8))
     :  +D5*(TMP5(II-3)+TMP5(II-7))+D6*(TMP5(II-4)+TMP5(II-6))

      ZK2(II) = ZK2(II-10)+D1*(TMP6(II-10)+TMP6(II))+D2*TMP6(II-5)
     :  +D3*(TMP6(II-1)+TMP6(II-9))+D4*(TMP6(II-2)+TMP6(II-8))
     :  +D5*(TMP6(II-3)+TMP6(II-7))+D6*(TMP6(II-4)+TMP6(II-6))

   5  CONTINUE

C Z(i)**k has been defined as Z(i)**k*(r</r>)**k
        DO 6 II = 2,N
          ZK1(II) = ZK1(II)/RTTK(II)
          ZK2(II) = ZK2(II)/RTTK(II)
   6    CONTINUE

C R(i)**(K+1) has been defined
        DO 7 II = 2,N
          RTTK(II) = RTTK(II)*R(II)
   7    CONTINUE

C Calculates Integrand for W(i)**k/r(i)**(k+1) at different grids
        DO 8 II = 2,N
          TMP5(II) = RHOP1(II)/RTTK(II)
          TMP6(II) = RHOP2(II)/RTTK(II)
   8    CONTINUE

C Last 10 points of W(i)**k are taken same as given below
      DO 9 II = N,N+10

      WK1(II) = D1*(TMP5(N)+TMP5(N+10))+D2*TMP5(N+5)
     :  +D3*(TMP5(N+1)+TMP5(N+9))+D4*(TMP5(N+2)+TMP5(N+8))
     :  +D5*(TMP5(N+3)+TMP5(N+7))+D6*(TMP5(N+4)+TMP5(N+6))

      WK2(II) = D1*(TMP6(N)+TMP6(N+10))+D2*TMP6(N+5)
     :  +D3*(TMP6(N+1)+TMP6(N+9))+D4*(TMP6(N+2)+TMP6(N+8))
     :  +D5*(TMP6(N+3)+TMP6(N+7))+D6*(TMP6(N+4)+TMP6(N+6))

   9  CONTINUE

C W(i)**k has been calculated using 10 integration points.
      DO 10 II = N,2,-1

      WK1(II) = WK1(II+10)+D1*(TMP5(II)+TMP5(II+10))+D2*TMP5(II+5)
     :  +D3*(TMP5(II+1)+TMP5(II+9))+D4*(TMP5(II+2)+TMP5(II+8))
     :  +D5*(TMP5(II+3)+TMP5(II+7))+D6*(TMP5(II+4)+TMP5(II+6))

      WK2(II) = WK2(II+10)+D1*(TMP6(II)+TMP6(II+10))+D2*TMP6(II+5)
     :  +D3*(TMP6(II+1)+TMP6(II+9))+D4*(TMP6(II+2)+TMP6(II+8))
     :  +D5*(TMP6(II+3)+TMP6(II+7))+D6*(TMP6(II+4)+TMP6(II+6))

   10 CONTINUE

C W(i)**k has been defined as Z(i)**k*(r</r>)**(k+1)
        DO 11 II = 2,N
          WK1(II) = WK1(II)*RTTK(II)
          WK2(II) = WK2(II)*RTTK(II)
   11   CONTINUE

C YZK functions are determined as Z(i)**k+W(i)**k
      DO 12 II = 2,N
        TMP1(II) = ZK1(II) + WK1(II)
        TMP2(II) = ZK2(II) + WK2(II)
        TMP3(II) = ZK1(II)
        TMP4(II) = ZK2(II)
        TMP5(II) = WK1(II)
        TMP6(II) = WK2(II)
   12 CONTINUE

      RETURN
      END

*************************************************************************
      function gml(nk)
      implicit real*8(a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c


      gml=zero

      if(mod(nk,2).eq.0)then

      lc=0

      pdt=one

  1   continue

      pdt=pdt*(2*lc+1)/4.0d0

      if(lc.eq.nk)goto 2

      lc=lc+2

        goto 1

  2   continue

      gml=4.0d0*pdt*3.6256099082219083119d0

      else

      lc=1

      pdt=one

  3   continue

      pdt=pdt*(2*lc+1)/4.0d0

      if(lc.eq.nk)goto 4

      lc=lc+2

        goto 3

  4   continue

      gml=4.0d0*pdt*pi*dsqrt(two)/(3.0d0*3.6256099082219083119d0)

      endif

      return
      end
*************************************************************************

*************************************************************************
      function gmm(nk)
      implicit real*8(a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c


      gmm=zero

      lc=0

      pdt=one

  1   continue

      pdt=pdt*(4*lc+1)/4.0d0

      if(lc.eq.nk)goto 2

      lc=lc+1

        goto 1

  2   continue

      gmm=4.0d0*pdt*3.6256099082219083119d0

      return
      end
*************************************************************************

*************************************************************************
      function gmn(nk)
      implicit real*8(a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten
     :      /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c


      gmn=zero

      lc=1

      pdt=one

  1   continue

      pdt=pdt*(4*lc+1)/4.0d0

      if(lc.eq.nk)goto 2

      lc=lc+1

        goto 1

  2   continue

      gmn=4.0d0*pdt*pi*dsqrt(two)/(5.0d0*3.6256099082219083119d0)

      return
      end
*************************************************************************

*************************************************************************
      function gm1(nk)
      implicit real*8(a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten

      gm1=zero

      if(nk.eq.1)gm1=one
      if(nk.eq.2)gm1=1.7724538509055160273d0
      if(nk.eq.3)gm1=2.6789385347077476337d0
      if(nk.eq.4)gm1=3.6256099082219083119d0
      if(nk.eq.5)gm1=4.5908437119988030532d0
      if(nk.eq.6)gm1=5.5663160017802352043d0
      if(nk.eq.7)gm1=6.5480629402478244377d0

      return
      end
*************************************************************************

*************************************************************************
      function gm12(nk)
      implicit real*8(a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c

      gm12=zero

      if(nk.eq.0)gm12=gm1(1)
      if(nk.eq.1)gm12=gm1(2)
      if(nk.eq.2)gm12=two*pi/(gm1(3)*dsqrt(three))
      if(nk.eq.3)gm12=dsqrt(two)*pi/gm1(4)
      if(nk.eq.4)gm12=dsqrt(two)*pi*dsqrt(5.0d0+dsqrt(5.0d0))/
     :(dsqrt(5.0d0)*gm1(5))
      if(nk.eq.5)gm12=pi**(three/two)*two**(two*two/three)/
     :(dsqrt(three)*gm1(3)**2)
      if(nk.eq.6)gm12=5.5663160017802352043
      if(nk.eq.7)gm12=pi*two**(three/(two*two))*
     :dsqrt((dsqrt(two))+one)/gm1(8)

      return
      end
*************************************************************************

*************************************************************************
      function gm2p1(nk)
      implicit real*8(a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten


      gm2p1=zero
      gm2p1=(two+one/dfloat(nk))*(one+one/dfloat(nk))*gm1(nk)

      return
      end
*************************************************************************

*************************************************************************
      function gm2n1(nk)
      implicit real*8(a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten


      gm2n1=zero
      gm2n1=(one+dfloat((nk-1)/nk))*gm12(nk-1)

      return
      end
*************************************************************************

*************************************************************************
      function gm3n1(nk)
      implicit real*8(a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten


      gm3n1=zero
      gm3n1=(two+dfloat((nk-1)/nk))*(one+dfloat((nk-1)/nk))**gm12(nk-1)

      return
      end
*************************************************************************

*************************************************************************
      function gm4n1(nk)
      implicit real*8(a-h,o-z)

      Common/cons/zero,half,tenth,one,two,three,ten


      gm4n1=zero
      gm4n1=(three+dfloat((nk-1)/nk))*(two+dfloat((nk-1)/nk))*
     :(one+dfloat((nk-1)/nk))*gm12(nk-1)

      return
      end
*************************************************************************

*************************************************************************
      function ieta(lba,laa)

      ieta=0

      if(lba.eq.laa-1)ieta=laa
      if(lba.eq.laa+1)ieta=-laa-1

      return
      end
*************************************************************************

*************************************************************************
      function izeta(lbb,lab)

      izeta=0

      if(lbb.eq.lab-1)izeta=lab
      if(lbb.eq.lab+1)izeta=-lab-1

      return
      end

        subroutine sixj(rj1,rj2,rj3,rl1,rl2,rl3,s6j)
        implicit real*8 (a-h,o-z)
        common/factor/fct(0:50),mfd
!comment function s6j(j1,j2,j3,l1,l2,l3) calculates the 6j symbol
!        of its arguments.  it returns an error message if the
!        sum of angular momenta involved in a triangular
!        relationship is not an integer;
!for all j1,j2,j3,l1,l2,l3 such that
!        fixp(2*j1) and fixp(2*j2) and fixp(2*j3) and
!        fixp(2*l1) and fixp(2*l2) and fixp(2*l3)
! let s6j(j1,j2,j3,l1,l2,l3) =
!    begin scalar ws6j,delprod,zmin,zmax,tria,trib,tric,trid,s;

        if(mfd.eq.33)goto 704
        mfd=33
        fct(0)=1.
        fct(1)=1.
        do 20 i=2,33
        fct(i)=fct(i-1)*i
20      continue

704    tria=tri6j(rj1,rj2,rj3)
       trib=tri6j(rl1,rj2,rl3)
       tric=tri6j(rj1,rl2,rl3)
       trid=tri6j(rl1,rl2,rj3)
c       type *,tria,trib,tric,trid
       if(tria.eq.-1..or.trib.eq.-1..or.tric.eq.-1..or.trid.eq.-1.)then
       print *,'inconsistent arguments to 6j symbol'
       print *,tria,trib,tric,trid
       print *,rj1,rj2,rj3
       print *,rl1,rl2,rl3
       stop
        endif
       if(tria.eq.0..or.trib.eq.0..or.tric.eq.0..or.trid.eq.0.)then
c         type *,'triangle condition'
          s6j=0
        return
        endif
!    comment apply formula in edmonds [(6.3.7) p.99]

       delprod=tridel(rj1,rj2,rj3)*tridel(rj1,rl2,rl3)*
     *         tridel(rl1,rj2,rl3)*tridel(rl1,rl2,rj3)
c       type *,'delprod=',delprod
c       type *,'args to min',(rj1+rj2+rl1+rl2),(rj2+rj3+rl2+rl3),
c     *                       (rj3+rj1+rl3+rl1)
       izmin=dmax1((rj1+rj2+rj3),(rj1+rl2+rl3),(rl1+rj2+rl3)
     *           ,(rl1+rl2+rj3))
       izmax=dmin1((rj1+rj2+rl1+rl2),(rj2+rj3+rl2+rl3)
     *           ,(rj3+rj1+rl3+rl1))
        ws6j=0
c       type *,'izmin,izmax=',izmin,izmax
        do 10 iz=izmin,izmax
        z=iz
      ws6j=ws6j+((-1)**iz)*fct(int(z+1.))/(fct(int(z-rj1-rj2-rj3))*
     *           fct(int(z-rj1-rl2-rl3))*
     *           fct(int(z-rl1-rj2-rl3))*
     *           fct(int(z-rl1-rl2-rj3))*
     *           fct(int(rj1+rj2+rl1+rl2-z))*
     *           fct(int(rj2+rj3+rl2+rl3-z))*
     *           fct(int(rj3+rj1+rl3+rl1-z)))
10      continue
       s6j=delprod*ws6j
       return
       end


        function tridel(rm1,rm2,rm3)
        implicit real*8 (a-h,o-z)
        common/factor/fct(0:50),mfd

c       type *,'+++',rm1,rm2,rm3
!comment function tridel evaluates the delta symbol defined on
!        p.99 of edmonds;
!for all m1,m2,m3 such that fixp(2*m1) and fixp(2*m2) and fixp(2*m3)
c       type *,m1+m2-m3,m1-m2+m3,-m1+m2+m3,m1+m2+m3+1

        i1=rm1+rm2-rm3
        i2=rm1-rm2+rm3
        i3=-rm1+rm2+rm3
        i4=rm1+rm2+rm3+1

        tridel=sqrt(fct(i1)*fct(i2)*fct(i3)/fct(i4))

c       type *,'***',rm1,rm2,rm3

        return
        end
        function tri6j(rj1,rj2,rj3)
        implicit real*8 (a-h,o-z)
        sumj=rj1+rj2+rj3
        isumj=int(sumj)
        if(abs(sumj-isumj).gt.1e-3)then !if j1+j2+j3 is not an integer
          tri6j=-1                      !then tri6j=-1
        goto 999
        endif
        if((2.*abs(rj1-rj2).le.2.*rj3).and.
     *        (2.*(rj1+rj2).ge.2.*rj3))then
          tri6j=1                               !triangle condition is satisfied
        goto 999
        endif

        tri6j=0                         !j1+j2+j3 is an integer
c       type *,rj1,rj2,rj3
999     return
        end

* ------------------------------------------------------------------

*************************************************************************
CC       **************************************************************    
CC       SAMPLE INPUT FILE FOR MG WITH UNIVERSAL BASIS
CC
CC       $CAT SCFINP.DAT
CC       Mg DF SCF Calculation
CC       MG
CC       7   35  28   28   25  25   20  20
CC            3   1    1   0    0   0   0
CC       UB
CC       0.00625   2.63
CC       50      7    24.305  12
CC       maxit,npower,amass,   z
CC       2     0     0    0    0     1    1   0
CC       isgo,jpnt,ibrt,ispd,idflt,Ieng,inuc,ichrg
CC       0       0     0    1   0   0
CC       istp,ialpvr,ihbd,intg,ivn,ival
CC       1  1  1  0  0  0  0
CC       7  5  5  3  3  2  2
CC       **************************************************************    
