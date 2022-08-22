C ======================================================================
C User Subroutine UEL for Abaqus: two elements for the split scheme
C  operator phase and displacement problem:
C Type 1: C3D12 displacement triangle element
C Type 2: C3D24 displacement rectangular element
C Type 3: C3P4 phase-field triangle element
C Type 4: C3P8 phase-field rectangular element
C ======================================================================
C
C Copyright© Gergely Molnar, CNRS, INSA of Lyon, LaMCoS, UMR5259
C
C gergely.molnar@insa-lyon.fr
C http://www.molnar-research.com
C
C This software is a computer program whose purpose is to calculate
C the fracture resistance of elasto-plastic solids in Abaqus under
C static and dynamic conditions.
C
C This software is governed by the CeCILL-B license under French law and
C abiding by the rules of distribution of free software. You can  use, 
C modify and/or redistribute the software under the terms of the CeCILL-B
C license as circulated by CEA, CNRS and INRIA at the following URL
C "http://www.cecill.info". 
C
C As a counterpart to the access to the source code and  rights to copy,
C modify and redistribute granted by the license, users are provided only
C with a limited warranty  and the software's author,  the holder of the
C economic rights,  and the successive licensors  have only  limited
C liability. 
C
C In this respect, the user's attention is drawn to the risks associated
C with loading,  using,  modifying and/or developing or reproducing the
C software by the user in light of its specific status of free software,
C that may mean  that it is complicated to manipulate,  and  that  also
C therefore means  that it is reserved for developers  and  experienced
C professionals having in-depth computer knowledge. Users are therefore
C encouraged to load and test the software's suitability as regards their
C requirements in conditions enabling the security of their systems and/or 
C data to be ensured and,  more generally, to use and operate it in the 
C same conditions as regards security. 
C
C The fact that you are presently reading this means that you have had
C knowledge of the CeCILL-B
C
C Disclaimer
C The restart option has not been tested.
C
C ======================================================================
C Material properties to be given through the input file (*.inp):
C
C For Type 1 element (stress-strain):
C PROPS(1) = Young's modulus (E)
C PROPS(2) = Poisson's ratio (nu)
C PROPS(3) = Yield stress (sig_y)
C PROPS(4) = Hardening modulus (H)
C PROPS(5) = Critical plastic strain (eps_pl_crit)
C PROPS(6) = Density (rho)
C PROPS(7) = Aniso energy degradation switch (if 1 - yes, 0 - no)
C PROPS(8) = Plasticuty switch (if 1 - yes, 0 - no)
C PROPS(9) = Length scale parameter (lc)
C PROPS(10) = Crack surface energy (gc)
C
C For Type 2 element (phase field):
C PROPS(1) = Length scale parameter (lc)
C PROPS(2) = Crack surface energy (gc)
C PROPS(3) = Elastic switch  (if 1 - yes, 0 - no)
C
C ---- Used variables ---------------------
C N_ELEM - number of elements used in the model divided
C            by 3 - (N_phase+N_stress+N_UMAT)/3 (to be changed for each model)
Cs
C NSTVTT - solution dependent variables for the displacement element
C            (displacements, strains, stresses,energies, phase, etc.)
C NSTVTO - solution dependent variables for the phase-field element
C            (phase, energy history)
C NSTV - overall solution dependent variables (NSTVTO+NSTVTT+4), where
C           the additional 4 variables are the: time and iteration number
C
C     ==================================================================
C     Comments on solution dependent variables
C     ==================================================================
C
C     Stress/strain element
C      SVARS(1-320): SDV(1-40)x8(or 1)
C                               SDV(1) - X translation
C                               SDV(2) - Y translation
C                               SDV(3) - Z translation
C                               SDV(4) - X normal strain
C                               SDV(5) - Y normal strain
C                               SDV(6) - Z normal strain
C                               SDV(7) - XY engineering shear strain
C                               SDV(8) - XZ engineering shear strain
C                               SDV(9) - YZ engineering shear strain
C                               SDV(10-15) - Elastic strains (x, y, z, xy, xz, yz)
C                               SDV(16-21) - Plastic strains
C                               SDV(22) - Eq. plastic strain
C                               SDV(23-28) - Stresses 
C                               SDV(29) - Hydrostatic stress
C                               SDV(30) - von Mises stress
C                               SDV(31) - plastic energy
C                               SDV(32) - tensile elastic energy
C                               SDV(33) - potential strain energy
C                               SDV(34) - phase-field
C                               SDV(35-40) - El. strain at the beginning of the step
C      SVARS(321-344): RHS(t_{n-1}) - 24(or 12) components, previous internal
C                                    force vector for HHT (dynamic)
C
C     Phase-field element
C      SVARS(1-16): SDV(1-2)x8(or 1)
C                               SDV(1) - phase-field
C                               SDV(2) - history energy
C      SVARS(17-24): RHS(t_{n-1}), 8(or 4) components previous internal
C                                    force vector for HHT (dynamic)
C
C ======================================================================
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
C     ==================================================================
      INCLUDE 'ABA_PARAM.INC'
C     ==================================================================
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,THREE=3.D0,
     1 TOLER=1.0D-8,FOUR=4.D0,RP25 = 0.25D0,HALF=0.5D0,SIX=6.D0,
     2 TEN=10.D0,PHCALCMAX=0.95D0,DEPSCR=0.1D0,DENGMAX=1,DTMIN=1.0D-7,
     3 N_ELEM=6400,NSTVTT=40,NSTVTO=2,NSTV=46)
C     ==================================================================
C     Initialization for all the element types
C     ==================================================================
      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,1),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
     
       INTEGER I,J,L,K,K1,K2,K3,K4,IX,IY,NALL,CNT

       REAL*8 AINTW(8),XII(8,3),XI(3),dNdxi(NNODE,3),
     1 VJACOB(3,3),dNdx(NNODE,3),VJABOBINV(3,3),AN(8),BP(3,NDOFEL),
     2 DP(3),SDV(NSTV),BB(6,NDOFEL),CMAT(6,6),EPS(6),STRESS(6),
     3 VNI(3,NDOFEL),ULOC(3),PHASENOD(NNODE),AMASS(NDOFEL,NDOFEL),
     4 EIGV(3),ALPHAI(3),CMATG(6,6),VECTI(3),EPSZ(6),
     5 EPSC(6),ASTIFF(NDOFEL,NDOFEL),ALOC(3),VLOC(3),RHSINI(NDOFEL),
     6 RHSK(NDOFEL),DEPS(6),EELAS(6),EPLAS(6),SFULL(6),FLOW(6),
     7 CMATP(6,6),EPSP(6),DEPLAS(6),DDSDDEEQ(6,6),
     8 FLOWG(6),PHMAX(8),PNEWDTIP(8)

       REAL*8 DTM,HIST,CLPAR,GCPAR,EMOD,ENU,PARK,ENG,ENGK,ENGEG,
     1 ENGKG,ENGDG,ENGPG,ENGD,PLSWT,ANISOSWT,EQPLAS,YSHARD,YIELDS,
     2 DLAMB,VNEVEZ,ALPHA,ENGPL,REITER
C
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,8)
C
C     History variables
        ENGKG=ZERO
        ENGEG=ZERO
        ENGDG=ZERO
        ENGPG=ZERO
        NALL=NSTVTO+NSTVTT  
C     ==================================================================
C     ******************************************************************
C     Constructing elemet TYPE 1 and 2 (stress/strain - displacement)
C     ******************************************************************
C     ==================================================================
C      
       IF ((JTYPE.EQ.ONE).OR.(JTYPE.EQ.TWO)) THEN
C     ==================================================================
C     Time an iteration variables
C     ==================================================================
       IF (TIME(2).EQ.ZERO) THEN
        TIMEZ=-999.D0
        TIMEZLOC=-99.D0
        DO K2=1,NSTV
         DO K3=1,8
           USRVAR(JELEM,K2,K3)=ZERO
         END DO
        END DO        
       ELSE
        TIMEZ=USRVAR(JELEM,NALL+1,1)
        TIMEZLOC=TIME(2)-DTIME
       ENDIF
       DTZERO=USRVAR(JELEM,NALL+2,1)
       IF (TIMEZ.LT.TIMEZLOC) THEN
        USRVAR(JELEM,NALL+1,1)=TIMEZLOC
        USRVAR(JELEM,NALL+2,1)=DTIME
        USRVAR(JELEM,NALL+3,1)=ZERO
        USRVAR(JELEM,NALL+4,1)=ZERO
       ELSE
        IF (DTZERO.GT.DTIME*(ONE+TOLER)) THEN
C -----   New correcting iteration   -----
         USRVAR(JELEM,NALL+2,1)=DTIME
         USRVAR(JELEM,NALL+3,1)=USRVAR(JELEM,NALL+3,1)+ONE
         USRVAR(JELEM,NALL+4,1)=ZERO
        ELSE
C -----   New local step   -----
         USRVAR(JELEM,NALL+4,1)=USRVAR(JELEM,NALL+4,1)+ONE
        ENDIF
       ENDIF      
       REITER=USRVAR(JELEM,NALL+3,1)
       STEPITER=USRVAR(JELEM,NALL+4,1)
C     ==================================================================
C     Additional plasticity control
C     ==================================================================
       IF ((REITER.EQ.ZERO).AND.(STEPITER.EQ.ZERO)) THEN
        PLSWTGLOB=ONE
        USRVAR(JELEM,NALL+1,2)=PLSWTGLOB
       ELSE
        PLSWTGLOB=USRVAR(JELEM,NALL+1,2)
       ENDIF
C       
C     ==================================================================
C     Material parameters
C     ==================================================================
       EMOD = PROPS(1)
       ENU = PROPS(2)
       YIELDS = PROPS(3)
       HMOD = PROPS(4)
       EQPLASCR = PROPS(5)
       DENS = PROPS(6)
       ANISOSWT = PROPS(7)
       PLSWT = PROPS(8)
       CLPAR = PROPS(9)
       GCPAR = PROPS(10)
       PARK = TOLER
       ELAMEL=EMOD*ENU/((ONE+ENU)*(ONE-TWO*ENU))
       ELAMEG=EMOD/(TWO*(ONE+ENU))
C     ==================================================================
C     Initial preparations
C     ==================================================================
       DO K1 = 1, NDOFEL                      
        DO KRHS = 1, NRHS
         RHS(K1,KRHS) = ZERO
        END DO
        RHSK(K1) = ZERO
        RHSINI(K1) = ZERO
        DO K2 = 1, NDOFEL
         AMATRX(K2,K1) = ZERO
         AMASS(K2,K1) = ZERO
         ASTIFF(K2,K1) = ZERO
        END DO
       END DO
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
       IF (JTYPE.EQ.ONE) THEN
        XII(1,1) = ONE/FOUR
        XII(1,2) = ONE/FOUR
        XII(1,3) = ONE/FOUR
        INNODE=1
        DO I=1,INNODE
         AINTW(I) = ONE/SIX
        END DO
       ELSEIF (JTYPE.EQ.TWO) THEN
        XII(1,1) = MONE/THREE**HALF
        XII(1,2) = MONE/THREE**HALF
        XII(1,3) = MONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = MONE/THREE**HALF
        XII(2,3) = MONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(3,3) = MONE/THREE**HALF
        XII(4,1) = MONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        XII(4,3) = MONE/THREE**HALF
        XII(5,1) = MONE/THREE**HALF
        XII(5,2) = MONE/THREE**HALF
        XII(5,3) = ONE/THREE**HALF
        XII(6,1) = ONE/THREE**HALF
        XII(6,2) = MONE/THREE**HALF
        XII(6,3) = ONE/THREE**HALF
        XII(7,1) = ONE/THREE**HALF
        XII(7,2) = ONE/THREE**HALF
        XII(7,3) = ONE/THREE**HALF
        XII(8,1) = MONE/THREE**HALF
        XII(8,2) = ONE/THREE**HALF
        XII(8,3) = ONE/THREE**HALF
        INNODE = 8
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
       ENDIF
C
C     ==================================================================
C     Determining maximum phase value in the element
C     ==================================================================
        PHELEMAX=ZERO
        DO K1=1,8
         PHMAX(K1)=ZERO
         PNEWDTIP(K1)=TEN
        END DO
        DO K1=1,INNODE
         IF ((STEPITER.EQ.ZERO).AND.(REITER.EQ.ZERO)) THEN
          PHMAX(K1)=USRVAR(JELEM,NSTVTT+1,K1)
         ELSE
          PHMAX(K1)=USRVAR(JELEM,34,K1)
         ENDIF        
        END DO
        PHELEMAX=MAXVAL(PHMAX)
C
C     ==================================================================
C     Calculating properties at each integration point
C     ==================================================================
       DO INPT=1,INNODE
C     Initializing solution dependent variables (phase,history)
        DO I=1,NSTVTT
          SDV(I)=SVARS(NSTVTT*(INPT-1)+I)
        END DO
C
C     Local coordinates of the integration point
        XI(1) = XII(INPT,1)
        XI(2) = XII(INPT,2)
        XI(3) = XII(INPT,3) 
C     Shape functions and local derivatives
        IF (JTYPE.EQ.ONE) THEN
         CALL SHAPEFUNT(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.TWO) THEN
         CALL SHAPEFUN(AN,dNdxi,XI)
        ENDIF
C     Shape functions
        IZ=ZERO
        DO I = 1,NNODE
         IX=IZ+ONE
         IY=IX+ONE
         IZ=IY+ONE
         VNI(1,IX)=AN(I)
         VNI(2,IX)=ZERO
         VNI(3,IX)=ZERO
         VNI(1,IY)=ZERO
         VNI(2,IY)=AN(I)
         VNI(3,IY)=ZERO
         VNI(1,IZ)=ZERO
         VNI(2,IZ)=ZERO
         VNI(3,IZ)=AN(I)
        END DO
C     Jacobian
        DO I = 1,3
         DO J = 1,3
          VJACOB(I,J) = ZERO
          DO K = 1,NNODE
           VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
          END DO
         END DO
        END DO
C        
        DTM = ZERO
        DTM = VJACOB(1,1)*VJACOB(2,2)*VJACOB(3,3)+VJACOB(1,2)*
     1   VJACOB(2,3)*VJACOB(3,1)+VJACOB(1,3)*VJACOB(2,1)*
     2   VJACOB(3,2)-VJACOB(3,1)*VJACOB(2,2)*VJACOB(1,3)-
     3   VJACOB(3,2)*VJACOB(2,3)*VJACOB(1,1)-VJACOB(3,3)*
     4   VJACOB(2,1)*VJACOB(1,2)
        IF (DTM.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DTM
         CALL XIT
        ENDIF
C     Inverse of Jacobian
        VJABOBINV(1,1)=(VJACOB(2,2)*VJACOB(3,3)-VJACOB(2,3)*
     1   VJACOB(3,2))/DTM
        VJABOBINV(1,2)=-(VJACOB(1,2)*VJACOB(3,3)-VJACOB(3,2)*
     1   VJACOB(1,3))/DTM
        VJABOBINV(1,3)=(VJACOB(1,2)*VJACOB(2,3)-VJACOB(1,3)*
     1   VJACOB(2,2))/DTM
        VJABOBINV(2,1)=-(VJACOB(2,1)*VJACOB(3,3)-VJACOB(2,3)*
     1   VJACOB(3,1))/DTM
        VJABOBINV(2,2)=(VJACOB(1,1)*VJACOB(3,3)-VJACOB(1,3)*
     1   VJACOB(3,1))/DTM
        VJABOBINV(2,3)=-(VJACOB(1,1)*VJACOB(2,3)-VJACOB(1,3)*
     1   VJACOB(2,1))/DTM
        VJABOBINV(3,1)=(VJACOB(2,1)*VJACOB(3,2)-VJACOB(2,2)*
     1   VJACOB(3,1))/DTM
        VJABOBINV(3,2)=-(VJACOB(1,1)*VJACOB(3,2)-VJACOB(1,2)*
     1   VJACOB(3,1))/DTM
        VJABOBINV(3,3)=(VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*
     1   VJACOB(2,1))/DTM
C        
C     Derivatives of shape functions respect to global ccordinates
        DO K = 1,NNODE
         DO I = 1,3
          dNdx(K,I) = ZERO
          DO J = 1,3
           dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
          END DO
         END DO
        END DO
C
C     Calculating B matrix (B=LN)
       IZ=ZERO
       DO INODE=1,NNODE
        IX=IZ+ONE
        IY=IX+ONE
        IZ=IY+ONE
        BB(1,IX)= dNdx(INODE,1)
        BB(2,IX)= ZERO
        BB(3,IX)= ZERO
        BB(4,IX)= dNdx(INODE,2)
        BB(5,IX)= dNdx(INODE,3)
        BB(6,IX)= ZERO
        BB(1,IY)= ZERO
        BB(2,IY)= dNdx(INODE,2)
        BB(3,IY)= ZERO
        BB(4,IY)= dNdx(INODE,1)
        BB(5,IY)= ZERO
        BB(6,IY)= dNdx(INODE,3)
        BB(1,IZ)= ZERO
        BB(2,IZ)= ZERO
        BB(3,IZ)= dNdx(INODE,3)
        BB(4,IZ)= ZERO
        BB(5,IZ)= dNdx(INODE,1)
        BB(6,IZ)= dNdx(INODE,2)
       END DO
C       
C     ==================================================================
C     Nodal displacements
C     ==================================================================
        DO J=1,3
         ULOC(J)=ZERO
        END DO
        DO J=1,3
         DO I=1,NDOFEL
          ULOC(J)=ULOC(J)+VNI(J,I)*U(I)
         END DO
        END DO  
        DO J=1,3
         SDV(J)=ULOC(J)
        END DO
C   
C     ==================================================================
C     Nodal velocities
C     ==================================================================
        DO J=1,3
         VLOC(J)=ZERO
        END DO
        DO J=1,3
         DO I=1,NDOFEL
          VLOC(J)=VLOC(J)+VNI(J,I)*V(I)
         END DO
        END DO  
C
C     ==================================================================
C     Nodal accelerations
C     ==================================================================
        DO J=1,3
         ALOC(J)=ZERO
        END DO
        DO J=1,3
         DO I=1,NDOFEL
          ALOC(J)=ALOC(J)+VNI(J,I)*A(I)
         END DO
        END DO  
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
        IF ((STEPITER.EQ.ZERO).AND.(REITER.EQ.ZERO)) THEN
         PHASE=USRVAR(JELEM,NSTVTT+1,INPT)
        ELSE
         PHASE=USRVAR(JELEM,34,INPT)
        ENDIF
        IF (PHASE.GT.ONE) THEN
         PHASE=ONE
        ELSEIF (PHASE.LT.ZERO) THEN
          PHASE=ZERO
        ENDIF
C
        SDV(34)=PHASE
C     ==================================================================
C     Calculating strain
C     ==================================================================
        DO J=1,6
         EPS(J)=ZERO
         DEPS(J)=ZERO
        END DO
        DO I=1,6
         DO J=1,NDOFEL
          EPS(I)=EPS(I)+BB(I,J)*U(J)    
         END DO
        END DO
        DO J=1,6
         DEPS(J)=EPS(J)-SDV(J+3)
        END DO
        DO J=1,6
         SDV(J+3)=EPS(J)
        END DO
C     ==================================================================
C     Recovering elastic and plastic strains from previous step
C     ==================================================================
       DO K1=1,6 
        EELAS(K1)=SDV(K1+9)
        EPLAS(K1)=SDV(K1+15)
        EPSZ(K1)=SDV(K1+9)
       END DO
       EQPLAS=SDV(22)
C
       DO K1=1,6 
        EELAS(K1)=EELAS(K1)+DEPS(K1)
       END DO
C
C     ==================================================================
C     Calculating eigenvalue decomposition
C     ==================================================================
C    Only updating the stiffness matrix in the fist iteration step
       DO K1=1,6
        EPSC(K1)=ZERO
       END DO
C
       IF ((STEPITER.LE.(3-REITER)).AND.(REITER.LT.4)) THEN
C ------- Normal case: the stiffness is refreshed in every iteration
C         based on the actual strain state
         DO K1=1,6
          EPSC(K1)=EELAS(K1)
          SDV(34+K1)=EPSC(K1)
         END DO
C     
       ELSEIF ((REITER.GE.4).AND.(STEPITER.EQ.ZERO)) THEN
C ------- If maximum trials are exceeded, the strain state from
C         the original (converged) step is taken
        DO K1=1,6
         EPSC(K1)=EPSZ(K1)
         SDV(34+K1)=EPSC(K1)
        END DO
       ELSE
C ------- If the stiffness is not refreshed, the strain is recovered
C         from an older NR iteration
        DO K1=1,6
         EPSC(K1)=USRVAR(JELEM,K1+34,INPT)
         SDV(34+K1)=EPSC(K1)
        END DO
       ENDIF
C       
       VALMDEPS=MAXVAL(ABS(DEPS))
       IF (VALMDEPS.GT.(0.1)) THEN
        USRVAR(JELEM,NALL+1,2)=ZERO
        PLSWTGLOB=ZERO
       ENDIF
C
C     ==================================================================
C     Calculating degrated elastic stiffness matrix
C     ==================================================================
       DO I=1,6
        DO J=1,6
        CMATG(I,J)=ZERO
        END DO
       END DO
       IF (PHASE.LT.TOLER) THEN
        CMATG(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMATG(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMATG(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMATG(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMATG(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMATG(1,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMATG(3,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMATG(2,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMATG(3,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMATG(4,4)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMATG(5,5)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMATG(6,6)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
         DO I=1,6
          DO J=1,6
           CMATG(I,J)=CMATG(I,J)*((ONE-PHASE)**TWO+PARK)
          END DO
         END DO
        ELSE
         CALL HMAT(EPSC,CMATG,ENU,EMOD,ANISOSWT,PLSWT,PHASE,PARK)
        ENDIF
C
       DO I=1,6
        DO J=1,6
         CMATP(I,J)=CMATG(I,J)
        END DO
       END DO
C
C     ==================================================================
C
C     ------------------ Testing the yield criterion -------------------
C
C     ==================================================================
C   Expending the 2D stress to 3D (as in plane-strain sig_z~=0)
C   Strain vectors
      DO K1=1,6
       EPSP(K1)=EELAS(K1)
       DEPLAS(K1)=ZERO
      END DO
C
C     Initial stress prediction
        DO K1=1,6
         SFULL(K1)=ZERO
        END DO
         DO I=1,6
          DO J=1,6
           SFULL(I)=SFULL(I)+CMATP(I,J)*EPSP(J)
          END DO
         END DO
C   Hardening yield strength
C
       YSHARD=(YIELDS+EQPLAS*HMOD)*((ONE-PHASE)**TWO+PARK)
C     Stress invariants
C       
       SMISES=(SFULL(1)-SFULL(2))*(SFULL(1)-SFULL(2)) + 
     1 (SFULL(2)-SFULL(3))*(SFULL(2)-SFULL(3)) + 
     2 (SFULL(3)-SFULL(1))*(SFULL(3)-SFULL(1)) 
       DO K1=4,6
        SMISES=SMISES+SIX*SFULL(K1)*SFULL(K1)
       END DO
       SMISES=SQRT(SMISES/TWO)
       FFUN=SMISES-YSHARD
C
       DLAMB=ZERO
       IF ((PLSWT.GT.TOLER).AND.(PLSWTGLOB.GT.TOLER)) THEN
C ==========	Testing yield criterion =============
C
        IF ((FFUN.GT.TOLER).AND.(PHELEMAX.LT.PHCALCMAX)) THEN
         FLOW(1)=TWO*SFULL(1)-SFULL(2)-SFULL(3)
         FLOW(2)=TWO*SFULL(2)-SFULL(1)-SFULL(3)
         FLOW(3)=TWO*SFULL(3)-SFULL(2)-SFULL(1)
         DO K1=4,6
          FLOW(K1)=SIX*SFULL(K1)
         END DO
         DO K1=1,6
          FLOW(K1)=FLOW(K1)*HALF/SMISES
         END DO
         DO K1=1,6
          FLOWG(K1)=FLOW(K1)
         END DO
C
C   Calculating the relationship between EQPLAS and DLAMB
        VNEVEZ=ZERO
        DO K1=1,6
         DO K2=1,6
          VNEVEZ=VNEVEZ+FLOW(K2)*CMATP(K2,K1)*FLOWG(K1)
         END DO
        END DO
        VNEVEZ=VNEVEZ+(HMOD*((ONE-PHASE)**TWO+TOLER))
        DLAMB=FFUN/VNEVEZ
C
        DO K1=1,6
         DEPLAS(K1)=DLAMB*FLOWG(K1)
        END DO
        EQPLAS=EQPLAS+DLAMB
C
C   Updating stresses and strains
        DO K2=1,6
         DO K1=1,6
          SFULL(K2)=SFULL(K2)-CMATP(K2,K1)*DEPLAS(K1)
         END DO
         EELAS(K2)=EELAS(K2)-DEPLAS(K2)
         EPLAS(K2)=EPLAS(K2)+DEPLAS(K2)
        END DO
C       
C   Tangent matrix      
        DO K2=1,6
         DO K4=1,6
          DDSDDEEQ(K2,K4)= ZERO
         END DO
        END DO
C
        DO K1=1,6
         DO K2=1,6
          DO K3=1,6
           DO K4=1,6
            DDSDDEEQ(K2,K4)= DDSDDEEQ(K2,K4)+CMATP(K2,K1)*FLOWG(K1)*
     1        FLOW(K3)*CMATP(K3, K4)
           END DO
          END DO
         END DO
        END DO
C       
       SMISES=(SFULL(1)-SFULL(2))*(SFULL(1)-SFULL(2)) + 
     1 (SFULL(2)-SFULL(3))*(SFULL(2)-SFULL(3)) + 
     2 (SFULL(3)-SFULL(1))*(SFULL(3)-SFULL(1)) 
       DO K1=4,6
        SMISES=SMISES+SIX*SFULL(K1)*SFULL(K1)
       END DO
       SMISES=SQRT(SMISES/TWO)        
C             
        DO K2=1,6
         DO K1=1,6
          CMATP(K2,K1)=CMATP(K2,K1)-DDSDDEEQ(K2,K1)/VNEVEZ
         END DO
        END DO
         
        ENDIF
       ENDIF
C
C   Converting every result back to 2D
C
      DO K1=1,6
       SDV(K1+9)=EELAS(K1)
       SDV(K1+15)=EPLAS(K1)
      END DO
      SDV(22)=EQPLAS
C
C     Stresses
        DO J=1,6        
         STRESS(J)=SFULL(J)
        END DO
        DO J=1,6
         SDV(J+22)=SFULL(J)
        END DO
        HYDRO=(SFULL(1)+SFULL(2)+SFULL(3))/THREE
        SDV(29)=HYDRO
        SDV(30)=SMISES
C
C      Materials stiffness matrix
        DO K2=1,6
         DO K1=1,6
          CMAT(K1,K2)=CMATP(K1,K2)
         END DO
        END DO        
C
C     ==================================================================
C     Calculating elastic ENERGY
C     ==================================================================
        ENG=ZERO
        ENGP=ZERO
        ENGK=ZERO
C
C   Kinetic enery        
        DO K2=1,3
         ENGK=ENGK+VLOC(K2)**TWO*DENS*HALF
        END DO
C
C   Elastic strain enery        
        DO I=1,6
         EPSC(I)=EELAS(I)
        END DO
        CALL EIGOWN(EPSC,EIGV,ALPHA,ALPHAI,VECTI)
        IF (ANISOSWT.LT.TOLER) THEN
         DO K1=1,3
          ALPHAI(K1)=ONE
         END DO
         ALPHA=ONE
        ENDIF
        IF (PLSWT.GT.TOLER) THEN
         DO K1=1,3
          ALPHAI(K1)=ONE
         END DO
        ENDIF
C      
        ENGP=(ELAMEL*(ALPHA*(EIGV(1)+EIGV(2)+EIGV(3)))**TWO)/
     1  TWO+ELAMEG*((EIGV(1)*ALPHAI(1))**TWO+(EIGV(2)*
     2  ALPHAI(2))**TWO+(EIGV(3)*ALPHAI(3))**TWO)
        ENGN=(ELAMEL*((ONE-ALPHA)*(EIGV(1)+EIGV(2)+EIGV(3)))**
     1  TWO)/TWO+ELAMEG*((EIGV(1)*(ONE-ALPHAI(1)))**TWO+(EIGV(2)*
     2  (ONE-ALPHAI(2)))**TWO+(EIGV(3)*(ONE-ALPHAI(3)))**TWO)
C
C   Plastic energy
C        
        PLPEN=ZERO
        IF (EQPLASCR.EQ.ZERO) THEN
         PLPEN=ZERO        
        ELSE
         ENGMAXEL=(YIELDS+HMOD*(DEPSCR+ONE)*EQPLASCR)**TWO/
     1    SIX/ELAMEG
         ENGMAXPL=YIELDS*(DEPSCR+ONE)*EQPLASCR+HMOD*
     1    ((DEPSCR+ONE)*EQPLASCR)**TWO/TWO
         ENGMAXCR=GCPAR/TWO/CLPAR
C         
         PLPEN=ENGMAXCR-ENGMAXEL-ENGMAXPL
         PLPEN=(PLPEN+ABS(PLPEN))/TWO
         PLPEN=TWO*PLPEN/(DEPSCR*EQPLASCR)**TWO
        ENDIF
C
        ENGPL=ZERO
        IF (EQPLAS.GT.EQPLASCR) THEN
         ENGPL=HALF*(EQPLAS-EQPLASCR)**TWO*PLPEN
        ELSE
         ENGPL=ZERO
        ENDIF
        ENGPL=ENGPL+EQPLAS*(YIELDS+HALF*EQPLAS*HMOD)
C
C
C     ---------   Energy increment and time step control -------------
C
        DENG=ENGPL+ENGP-SDV(31)-SDV(32)
C        
        DENGMAXV=GCPAR/TWO/CLPAR*DENGMAX
        IF ((LFLAGS(1).EQ.1).OR.(LFLAGS(1).EQ.11)) THEN
C --------- If automatic time integration is used
         IF ((DENG.GT.DENGMAXV).AND.(REITER.LT.4)) THEN
C --------- If energy criterion is violated (dE>dE_max)
          IF (DTIME.GT.DTMIN*(ONE+TOLER)) THEN
C --------- If time-step still can be reduced
           PNEWDTIP(INPT)=DENGMAXV/DENG/TEN
           DTIMENEXT=PNEWDTIP(INPT)*DTIME
           IF (DTIMENEXT.LT.(DTMIN*(ONE+TOLER))) THEN
            PNEWDTIP(INPT)=PNEWDTIP(INPT)*DTMIN/DTIMENEXT
           ENDIF
          ELSE
C --------- If time-step is already too small
           PNEWDTIP(INPT)=ONE 
          ENDIF
         ENDIF
        ENDIF
C
        IF ((STEPITER.EQ.ZERO).AND.(REITER.EQ.ZERO)) THEN
        ELSE
         SDV(31)=ENGPL
         SDV(32)=ENGP
        ENDIF
        ENGKG=ENGKG+ENGK*DTM
        ENGEG=ENGEG+(ENGP*((ONE-PHASE)**TWO+PARK)+ENGN)*DTM
        ENGPG=ENGPG+ENGPL*((ONE-PHASE)**TWO+PARK)*DTM
        SDV(33)=((ENGPL+ENGP)*((ONE-PHASE)**TWO+PARK)+ENGN)
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
C
        DO K=1,NDOFEL
         DO L=1,NDOFEL
          DO I=1,6
           DO J=1,6
            ASTIFF(K,L)=ASTIFF(K,L)+AINTW(INPT)*BB(I,K)*CMAT(I,J)*
     1       BB(J,L)*DTM
           END DO
          END DO
         END DO
        END DO
C
C     ==================================================================
C     Calculating element mass matrix
C     ==================================================================
        DO K=1,NDOFEL
         DO L=1,NDOFEL
          DO I=1,3
            AMASS(K,K)=AMASS(K,K)+AINTW(INPT)*VNI(I,K)*VNI(I,L)
     1       *DTM*DENS
          END DO
         END DO
        END DO
C       
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
        IF ((LFLAGS(1).EQ.11).OR.(LFLAGS(1).EQ.12)) THEN
         DO K1=1,NDOFEL
          DO K2=1,3
           RHSK(K1)=RHSK(K1)-AINTW(INPT)*VNI(K2,K1)*A(K1)*
     1      DTM*DENS
          END DO
         END DO
        ENDIF
C
        DO K1=1,NDOFEL
         DO K4=1,6
           RHS(K1,1)=RHS(K1,1)-AINTW(INPT)*BB(K4,K1)*STRESS(K4)*DTM
         END DO
        END DO
C       
C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================
        DO I=1,NSTVTT
         SVARS(NSTVTT*(INPT-1)+I)=SDV(I)
         IF (LFLAGS(3).EQ.5) THEN
         ELSE
          USRVAR(JELEM,I,INPT)=SVARS(NSTVTT*(INPT-1)+I)
         ENDIF
        END DO
       END DO
C        
C     ==================================================================
C     Jacobien for the element
C     ==================================================================
        IF ((LFLAGS(1).EQ.1).OR.(LFLAGS(1).EQ.2)) THEN
C       Simple static calculation
         DO K=1,NDOFEL
          DO L=1,NDOFEL
             AMATRX(K,L)=AMATRX(K,L)+ASTIFF(K,L)
          END DO
         END DO
C
         ELSEIF ((LFLAGS(1).EQ.11).OR.(LFLAGS(1).EQ.12)) THEN
C         Dynamic calculation
C
          PARALPHA=PARAMS(1)
          PARBETA=PARAMS(2)
          DO I=1,NDOFEL
           RHSINI(I)=SVARS(I+NSTVTT*INNODE)
           SVARS(I+NSTVTT*INNODE)=RHS(I,1)
          END DO
C
          IF (LFLAGS(3).EQ.4) THEN
           DO I=1,NDOFEL
            RHS(I,1)=ZERO
           END DO
          ELSE
           DO I=1,NDOFEL
            RHS(I,1)=RHS(I,1)*(ONE+PARALPHA)-RHSINI(I)*PARALPHA+
     1        RHSK(I)
           END DO
          ENDIF  
C          
          IF ((LFLAGS(3).EQ.4).OR.(LFLAGS(3).EQ.6)) THEN
           DO K=1,NDOFEL
            DO L=1,NDOFEL
              AMATRX(K,L)=AMATRX(K,L)+AMASS(K,L)
            END DO
           END DO         
C           
          ELSEIF (LFLAGS(3).EQ.1) THEN
           DADU=ONE/(PARBETA*DTIME**TWO)
           DO K=1,NDOFEL
            DO L=1,NDOFEL
              AMATRX(K,L)=AMATRX(K,L)+AMASS(K,L)*DADU+(ONE+PARALPHA)*
     1         ASTIFF(K,L)
            END DO
           END DO 
          ENDIF
        ENDIF       
C     New time increment
       PNEWDTE=MINVAL(PNEWDTIP)
       IF (PNEWDTE.LT.(ONE+TOLER)) THEN
        PNEWDT=PNEWDTE
       ENDIF
C      
C     ==================================================================
C     ******************************************************************
C     Constructing elemet TYPE 3 and 4 (damage phase-field)
C     ******************************************************************
C     ==================================================================
      ELSEIF ((JTYPE.EQ.THREE).OR.(JTYPE.EQ.FOUR)) THEN
       REITER=USRVAR(JELEM-N_ELEM,NALL+3,1)
       STEPITER=USRVAR(JELEM-N_ELEM,NALL+4,1)
C     ==================================================================
C     Material parameters
C     ==================================================================
       CLPAR=PROPS(1)
       GCPAR =PROPS(2)
       ELSWT = PROPS(3)
C     ==================================================================
C     Initial preparations
C     ==================================================================
       DO K1 = 1, NDOFEL                      
        DO KRHS = 1, NRHS
         RHS(K1,KRHS) = ZERO
        END DO
        DO K2 = 1, NDOFEL
         AMATRX(K2,K1) = ZERO
        END DO
       END DO
C     ==================================================================
C     Local coordinates and weights
C     ==================================================================
       IF (JTYPE.EQ.THREE) THEN
        XII(1,1) = ONE/FOUR
        XII(1,2) = ONE/FOUR
        XII(1,3) = ONE/FOUR
        INNODE=1.0
        AINTW(1) = ONE/SIX
       ELSEIF (JTYPE.EQ.FOUR) THEN
        XII(1,1) = MONE/THREE**HALF
        XII(1,2) = MONE/THREE**HALF
        XII(1,3) = MONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = MONE/THREE**HALF
        XII(2,3) = MONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(3,3) = MONE/THREE**HALF
        XII(4,1) = MONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        XII(4,3) = MONE/THREE**HALF
        XII(5,1) = MONE/THREE**HALF
        XII(5,2) = MONE/THREE**HALF
        XII(5,3) = ONE/THREE**HALF
        XII(6,1) = ONE/THREE**HALF
        XII(6,2) = MONE/THREE**HALF
        XII(6,3) = ONE/THREE**HALF
        XII(7,1) = ONE/THREE**HALF
        XII(7,2) = ONE/THREE**HALF
        XII(7,3) = ONE/THREE**HALF
        XII(8,1) = MONE/THREE**HALF
        XII(8,2) = ONE/THREE**HALF
        XII(8,3) = ONE/THREE**HALF
        INNODE = 8.0
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
       ENDIF
C
C     ==================================================================
C     Calculating properties at each integration point
C     ==================================================================
       DO INPT=1,INNODE
C     Initializing solution dependent variables (phase,history)
        DO I=1,NSTVTO
          SDV(I)=SVARS(NSTVTO*(INPT-1)+I)
        END DO
C
C     Local coordinates of the integration point
        XI(1) = XII(INPT,1)
        XI(2) = XII(INPT,2) 
        XI(3) = XII(INPT,3) 
C     Shape functions and local derivatives
        IF (JTYPE.EQ.THREE) THEN
         CALL SHAPEFUNT(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.FOUR) THEN
         CALL SHAPEFUN(AN,dNdxi,XI)
        ENDIF
C     Jacobian
        DO I = 1,3
         DO J = 1,3
          VJACOB(I,J) = ZERO
          DO K = 1,NNODE
           VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
          END DO
         END DO
        END DO
C        
        DTM = ZERO
        DTM = VJACOB(1,1)*VJACOB(2,2)*VJACOB(3,3)+VJACOB(1,2)*
     1   VJACOB(2,3)*VJACOB(3,1)+VJACOB(1,3)*VJACOB(2,1)*
     2   VJACOB(3,2)-VJACOB(3,1)*VJACOB(2,2)*VJACOB(1,3)-
     3   VJACOB(3,2)*VJACOB(2,3)*VJACOB(1,1)-VJACOB(3,3)*
     4   VJACOB(2,1)*VJACOB(1,2)
C     
        IF (DTM.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DTM
         CALL XIT	
        ENDIF
C     Inverse of Jacobian
        VJABOBINV(1,1)=(VJACOB(2,2)*VJACOB(3,3)-VJACOB(2,3)*
     1   VJACOB(3,2))/DTM
        VJABOBINV(1,2)=-(VJACOB(1,2)*VJACOB(3,3)-VJACOB(3,2)*
     1   VJACOB(1,3))/DTM
        VJABOBINV(1,3)=(VJACOB(1,2)*VJACOB(2,3)-VJACOB(1,3)*
     1   VJACOB(2,2))/DTM
        VJABOBINV(2,1)=-(VJACOB(2,1)*VJACOB(3,3)-VJACOB(2,3)*
     1   VJACOB(3,1))/DTM
        VJABOBINV(2,2)=(VJACOB(1,1)*VJACOB(3,3)-VJACOB(1,3)*
     1   VJACOB(3,1))/DTM
        VJABOBINV(2,3)=-(VJACOB(1,1)*VJACOB(2,3)-VJACOB(1,3)*
     1   VJACOB(2,1))/DTM
        VJABOBINV(3,1)=(VJACOB(2,1)*VJACOB(3,2)-VJACOB(2,2)*
     1   VJACOB(3,1))/DTM
        VJABOBINV(3,2)=-(VJACOB(1,1)*VJACOB(3,2)-VJACOB(1,2)*
     1   VJACOB(3,1))/DTM
        VJABOBINV(3,3)=(VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*
     1   VJACOB(2,1))/DTM
C        
C     Derivatives of shape functions respect to global ccordinates
        DO K = 1,NNODE
         DO I = 1,3
          dNdx(K,I) = ZERO
          DO J = 1,3
           dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
          END DO
         END DO
        END DO
C
C     Calculating B matrix (B=LN)
       DO INODE=1,NNODE
        BP(1,INODE)=dNdx(INODE,1)
        BP(2,INODE)=dNdx(INODE,2)
        BP(3,INODE)=dNdx(INODE,3)
       END DO
C
C     ==================================================================
C     Nodal phase-field
C     ==================================================================
        PHASE=ZERO
        DPHASE=ZERO
        DO I=1,NDOFEL
         PHASE=PHASE+AN(I)*U(I)
        END DO
        DO I=1,NDOFEL
         DPHASE=DPHASE+AN(I)*DU(I,1)
        END DO
        SDV(1)=PHASE
C
C     Gradient
        DO I=1,3
         DP(I)=ZERO
        END DO
        DO I=1,3
         DO J=1,NNODE
          DP(I)=DP(I)+BP(I,J)*U(J)
         END DO
        END DO
C
C     ==================================================================
C     Calculating elastic ENERGY history
C     ==================================================================
        IF ((STEPITER.EQ.ZERO).AND.(REITER.EQ.ZERO)) THEN
         ENGN=USRVAR(JELEM-N_ELEM,31,INPT)+USRVAR(JELEM-N_ELEM,32,INPT)
         IF (ELSWT.GT.ZERO) THEN
           ENGN=ENGN-GCPAR/TWO/CLPAR
         ENDIF
        ELSE
         ENGN=USRVAR(JELEM-N_ELEM,NSTVTT+2,INPT)
        ENDIF
C       
        HISTN=USRVAR(JELEM-N_ELEM,NSTVTT+2,INPT)
        IF (ENGN.GT.HISTN) THEN
         HIST=ENGN
        ELSE
         HIST=HISTN
        ENDIF
        SDV(2)=HIST
C        
C     ==================================================================
C     Calculating fracture energy for history output
C     ==================================================================
C
        ENGD=ZERO
        ENGD=PHASE**TWO/TWO/CLPAR*DTM*GCPAR
C
        DO J=1,3
         ENGD=ENGD+DP(J)*DP(J)*CLPAR*DTM/TWO*GCPAR
        END DO
        ENGDG=ENGDG+ENGD
C
C     ==================================================================
C     Calculating element stiffness matrix
C     ==================================================================
        DO I=1,NNODE
         DO K=1,NNODE
          DO J=1,3
           AMATRX(I,K)=AMATRX(I,K)+BP(J,I)*BP(J,K)*DTM*
     1      GCPAR*CLPAR*AINTW(INPT)
          END DO
          AMATRX(I,K)=AMATRX(I,K)+AN(I)*AN(K)*DTM*
     1     AINTW(INPT)*(GCPAR/CLPAR+TWO*HIST)
         END DO
        END DO
C        
C     ==================================================================
C     Internal forces (residual vector)
C     ==================================================================
        DO I=1,NDOFEL
         DO J=1,3
           RHS(I,1)=RHS(I,1)-BP(J,I)*DP(J)*GCPAR*CLPAR*
     1      AINTW(INPT)*DTM
         END DO
         RHS(I,1)=RHS(I,1)-AN(I)*AINTW(INPT)*DTM*
     1    ((GCPAR/CLPAR+TWO*HIST)*PHASE-TWO*HIST)
        END DO
C
C     ==================================================================
C     Uploading solution dep. variables
C     ==================================================================
        DO I=1,NSTVTO
         SVARS(NSTVTO*(INPT-1)+I)=SDV(I)
         IF (LFLAGS(3).EQ.5) THEN
         ELSE
          USRVAR(JELEM-N_ELEM,I+NSTVTT,INPT)=SVARS(NSTVTO*(INPT-1)+I)
         ENDIF
        END DO
       END DO
C       
       DO I=1,NDOFEL
        RHSINI(I)=SVARS(I+NSTVTO*INNODE)
        SVARS(I+NSTVTO*INNODE)=RHS(I,1)
       END DO
       IF ((LFLAGS(1).EQ.11).OR.(LFLAGS(1).EQ.12)) THEN
        PARALPHA=PARAMS(1)
        DO I=1,NDOFEL
         RHS(I,1)=RHS(I,1)*(ONE+PARALPHA)+RHSINI(I)*PARALPHA
        END DO
        DO I=1,NNODE
         DO K=1,NNODE
           AMATRX(I,K)=AMATRX(I,K)*(ONE+PARALPHA)
         END DO
        END DO
       ENDIF
      ENDIF
        ENERGY(1)=ENGKG
        ENERGY(2)=ENGEG
        ENERGY(4)=ENGPG
        ENERGY(7)=ENGDG
      RETURN
      END
C
C --------------------------------------------------------------      
C Shape functions for tetrahedral elements
C --------------------------------------------------------------      
C
      SUBROUTINE SHAPEFUNT(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(4,3)
      Real*8 XI(3)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)

C     Values of shape functions as a function of local coord.
      AN(1) = ONE-XI(1)-XI(2)-XI(3)
      AN(2) = XI(1)
      AN(3) = XI(2)
      AN(4) = XI(3)
C
C     Derivatives of shape functions respect to local coordinates
      DO I=1,4
        DO J=1,3
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE
      dNdxi(1,2) =  MONE
      dNdxi(1,3) =  MONE
C
      dNdxi(2,1) =  ONE
      dNdxi(2,2) =  ZERO
      dNdxi(2,3) =  ZERO
C
      dNdxi(3,1) =  ZERO
      dNdxi(3,2) =  ONE
      dNdxi(3,3) =  ZERO
C
      dNdxi(4,1) =  ZONE
      dNdxi(4,2) =  ZONE
      dNdxi(4,3) =  ONE
C
      RETURN
      END
C
C --------------------------------------------------------------      
C Shape functions for brick elements
C --------------------------------------------------------------      
C
      SUBROUTINE SHAPEFUN(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      REAL*8 AN(8),dNdxi(8,3)
      REAL*8 XI(3)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0,EIGHT=8.D0)

C     Values of shape functions as a function of local coord.
      AN(1) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE-XI(3))
      AN(2) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE-XI(3))
      AN(3) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE-XI(3))
      AN(4) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE-XI(3))
      AN(5) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE+XI(3))
      AN(6) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE+XI(3))
      AN(7) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE+XI(3))
      AN(8) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE+XI(3))
      
C     Derivatives of shape functions respect to local coordinates
      DO I=1,8
        DO J=1,3
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
      dNdxi(1,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
      dNdxi(1,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
      dNdxi(2,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
      dNdxi(2,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
      dNdxi(2,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
      dNdxi(3,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
      dNdxi(3,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
      dNdxi(3,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
      dNdxi(4,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
      dNdxi(4,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
      dNdxi(4,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
      dNdxi(5,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
      dNdxi(5,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
      dNdxi(5,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
      dNdxi(6,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
      dNdxi(6,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
      dNdxi(6,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
      dNdxi(7,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
      dNdxi(7,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
      dNdxi(7,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
      dNdxi(8,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
      dNdxi(8,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
      dNdxi(8,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
      
      RETURN
      END
C
C Eigenstrains from Voigt notation
C
      SUBROUTINE EIGOWN(EPS,EIGV,ALPHA,ALPHAI,VECTI)
       INCLUDE 'ABA_PARAM.INC'
       PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,
     1  TS=27.D0,THREE=3.D0,HALF=0.5D0,TOLER=1.0D-12,FOUR=4.D0,
     2  CNTN=100,TOLERE=1.0D-12)
       INTEGER I, J, K
       REAL*8 EPS(6), EIGV(3), ALPHAI(3), VECTI(3)
       REAL*8 PC, QC, ALPHA, DISC, PI, CNT
C
       PI=FOUR*ATAN(ONE)
C Scaling the strain vector
       VMAXE=MAXVAL(ABS(EPS))
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPS(K1)/VMAXE
        END DO
       ENDIF
C              
C   Calculating eigenvalues
       VECTI(1)=EPS(1)+EPS(2)+EPS(3)
       VECTI(2)=EPS(1)*EPS(2)+EPS(2)*EPS(3)+EPS(1)*EPS(3)-
     1 EPS(4)**TWO/FOUR-EPS(5)**TWO/FOUR-EPS(6)**TWO/FOUR
       VECTI(3)=EPS(1)*EPS(2)*EPS(3)+EPS(4)*EPS(5)*EPS(6)/
     1 FOUR-EPS(1)*EPS(6)**TWO/FOUR-EPS(2)*EPS(5)**TWO/FOUR-
     2 EPS(3)*EPS(4)**TWO/FOUR
C
C   Depressed coefficients    
       PC=VECTI(2)-VECTI(1)**TWO/THREE
       QC=VECTI(1)*VECTI(2)/THREE-TWO*VECTI(1)**THREE/TS-VECTI(3)
       DISC=MONE*(FOUR*PC**THREE+TS*QC**TWO)
C
       DO I=1,3
        EIGV(I)=ZERO
       END DO
       CNT=ZERO
       IF (ABS(DISC).LT.TOLER) THEN
        IF ((ABS(QC).LT.TOLER).AND.(ABS(PC).LT.TOLER)) THEN
         EIGV(1)=VECTI(1)/THREE
         EIGV(2)=VECTI(1)/THREE
         EIGV(3)=VECTI(1)/THREE
        ELSE
         EIGV(1)=-THREE*QC/TWO/PC+VECTI(1)/THREE
         EIGV(2)=-THREE*QC/TWO/PC+VECTI(1)/THREE
         EIGV(3)=THREE*QC/PC+VECTI(1)/THREE
         IF (EIGV(1).GT.EIGV(3)) THEN
          EONE=EIGV(1)
          EIGV(1)=EIGV(3)
          EIGV(3)=EONE
         ENDIF
        ENDIF
       ELSE
        DO I=1,3
         EIGV(I)=VECTI(1)/THREE+TWO*(MONE*PC/THREE)**HALF*
     1   COS(ONE/THREE*ACOS(MONE*QC/TWO*(TS/(MONE*PC**THREE))**
     2   HALF)+TWO*I*PI/THREE)
        END DO
       ENDIF
C       
       ALPHA=ZERO
       IF ((EIGV(1)+EIGV(2)+EIGV(3)).GT.TOLER) THEN
        ALPHA=ONE
       ENDIF
       DO K1=1,3
        ALPHAI(K1)=ZERO
        IF (EIGV(K1).GT.TOLER) THEN
         ALPHAI(K1)=ONE
        ENDIF
       END DO
C
C    Rescaling eigenvalues       
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPS(K1)*VMAXE
        END DO
        DO K1=1,3
         EIGV(K1)=EIGV(K1)*VMAXE
        END DO
         VECTI(1)=EIGV(1)+EIGV(2)+EIGV(3)
         VECTI(2)=EIGV(1)*EIGV(2)+EIGV(1)*EIGV(3)+EIGV(3)*EIGV(2)
         VECTI(3)=EIGV(1)*EIGV(2)*EIGV(3)
       ENDIF
C      
       RETURN
       END
C        
       SUBROUTINE HMAT(EPSI,CMAT,ENU,EMOD,ANISOSWT,PLSWT,PHASE,PARK)
       INCLUDE 'ABA_PARAM.INC'
       PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,FOUR=4.D0,
     1  TOLER=1.0D-12,SIX=6.D0,FT=50.D0,THREE=3.D0,HALF=0.5D0,
     2  TS=27.D0,CNTM=1000,TEN=10.D0,TOLERE=1.0D-12,TOLD=1.0D-7)
       INTEGER I, J, K, NDIM
       REAL*8 EPS(6),VECTI(3),CLMAT(3,3),EIGV(3),ALPHAI(3),
     1  DLDI(3,3),DLDIDI(3,3,3),DIDE(3,6),DIDEDE(3,6,6),
     2  DLDEDE(3,6,6),CMAT(6,6),DLDE(3,6),EPSI(6),
     3  EIGVEC(3,3),EIGVAL(3,3),EPSNEW(3,3)
       REAL*8 VMAXE, ALPHA, ANISOSWT, DENOM, ENU, EMOD,
     1 PHASE, PC, QC, CNT, DISC, PLSWT
C
C Rescaling the strain vector
       VMAXE=MAXVAL(ABS(EPSI))
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPSI(K1)/VMAXE
        END DO
       ELSE
        DO K1=1,6
         EPS(K1)=EPSI(K1)
        END DO
       ENDIF
C       
C   Calculating invariants
       VECTI(1)=EPS(1)+EPS(2)+EPS(3)
       VECTI(2)=EPS(1)*EPS(2)+EPS(2)*EPS(3)+EPS(1)*EPS(3)-
     1 EPS(4)**TWO/FOUR-EPS(5)**TWO/FOUR-EPS(6)**TWO/FOUR
       VECTI(3)=EPS(1)*EPS(2)*EPS(3)+EPS(4)*EPS(5)*EPS(6)/
     1 FOUR-EPS(1)*EPS(6)**TWO/FOUR-EPS(2)*EPS(5)**TWO/FOUR-
     2 EPS(3)*EPS(4)**TWO/FOUR
       PC=VECTI(2)-VECTI(1)**TWO/THREE
       QC=VECTI(1)*VECTI(2)/THREE-TWO*VECTI(1)**THREE/TS-VECTI(3)
       DISC=MONE*(FOUR*PC**THREE+TS*QC**TWO)
C
C   Calculating eigenvalues
       CALL EIGOWN(EPS,EIGV,ALPHA,ALPHAI,VECTI)
       IF (PLSWT.GT.TOLER) THEN
        DO K1=1,3
         ALPHAI(K1)=ONE
        END DO
       ENDIF
C
C	Initialising CMAT
        DO I=1,6
         DO J=1,6
          CMAT(I,J)=ZERO
         END DO
        END DO
C
C ************ Starting spectral decomposition ************************        
       IF ((MINVAL(EIGV).GE.-TOLER).OR.(ANISOSWT.LT.TOLER)) THEN
C All eigenvalues are in tension, therefore the stiffness matrix can be
C     degraded without the decomposition
        CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(1,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(4,4)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(5,5)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(6,6)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        DO I=1,6
         DO J=1,6
          CMAT(I,J)=CMAT(I,J)*((ONE-PHASE)**TWO+PARK)
         END DO
        END DO
C
       ELSEIF ((MAXVAL(EIGV).LT.-TOLER).AND.(PLSWT.LT.HALF)) THEN
C All eigenvalues are in compression, therefore the stiffness
C     matrix does not need to be degraded 
        CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
        CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(1,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(2,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(3,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
        CMAT(4,4)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(5,5)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        CMAT(6,6)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
C
       ELSE
C   Calculating materials stiffness matrix
        DO I=1,3
         DO J=1,3
          CLMAT(I,J)=ZERO
         END DO
        END DO
        GAMMA=ENU/(ONE-TWO*ENU)
        CLMAT(1,1)=((ONE-ALPHAI(1)*PHASE)**TWO+PARK)+GAMMA*
     1  ((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(2,2)=((ONE-ALPHAI(2)*PHASE)**TWO+PARK)+GAMMA*
     1  ((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(3,3)=((ONE-ALPHAI(3)*PHASE)**TWO+PARK)+GAMMA*
     1  ((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(1,2)=GAMMA*((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(1,3)=GAMMA*((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(2,3)=GAMMA*((ONE-ALPHA*PHASE)**TWO+PARK)
        CLMAT(2,1)=CLMAT(1,2)
        CLMAT(3,1)=CLMAT(1,3)
        CLMAT(3,2)=CLMAT(2,3)
        DO I=1,3
         DO J=1,3
          CLMAT(I,J)=CLMAT(I,J)*EMOD/(ONE+ENU)
         END DO
        END DO
C
C Adding small permutation if two eigenvalues are close or equal 
        IF (ABS(DISC).LT.TOLD) THEN
         CALL JACOBYEIG(EPS,EIGVEC,EIGVAL)
         IF (ABS(EIGV(2)).LT.TOLER) THEN
          EIGV(2)=-0.1
         ELSE
          EIGV(2)=EIGV(2)*1.1D0 
         ENDIF
C
C Reconstructing strain matrix from the modified eigenstrains         
         DO I=1,3
          DO J=1,3
           EPSNEW(I,J)=ZERO
           DO K=1,3
            EPSNEW(I,J)=EPSNEW(I,J)+EIGVEC(I,K)*EIGV(K)*EIGVEC(J,K)
           END DO
          END DO
         END DO
         EPS(1)=EPSNEW(1,1)
         EPS(2)=EPSNEW(2,2)
         EPS(3)=EPSNEW(3,3)
         EPS(4)=EPSNEW(1,2)*TWO
         EPS(5)=EPSNEW(1,3)*TWO
         EPS(6)=EPSNEW(2,3)*TWO
         VECTI(1)=EIGV(1)+EIGV(2)+EIGV(3)
         VECTI(2)=EIGV(1)*EIGV(2)+EIGV(1)*EIGV(3)+EIGV(3)*EIGV(2)
         VECTI(3)=EIGV(1)*EIGV(2)*EIGV(3)
         PC=VECTI(2)-VECTI(1)**TWO/THREE
         QC=VECTI(1)*VECTI(2)/THREE-TWO*VECTI(1)**THREE/TS-VECTI(3)
         DISC=MONE*(FOUR*PC**THREE+TS*QC**TWO)       
        ENDIF
C        
C   Calculating derivatives of eigenvalues respect to invariants
        DO I=1,3
         DO J=1,3
          DLDI(I,J)=ZERO
         END DO
        END DO
        DO K=1,3
         DENOM=THREE*EIGV(K)**TWO-VECTI(1)*TWO*EIGV(K)+
     1   VECTI(2)
         IF (ABS(DENOM).LT.TOLER) THEN
          WRITE(7,*) 'EPS: ',EPS
          WRITE(7,*) 'EIGV: ',EIGV
          WRITE(7,*) 'VECTI: ',VECTI
          WRITE(7,*) 'DENOM: ',DENOM
          WRITE(7,*) 'Denominator is close to 0.'
C          CALL XIT
         ENDIF
         DLDI(K,1)=EIGV(K)**TWO/DENOM
         DLDI(K,2)=MONE*EIGV(K)/DENOM
         DLDI(K,3)=ONE/DENOM
        END DO 
C
        DO I=1,3
         DO J=1,3
          DO K=1,3
           DLDIDI(K,I,J)=ZERO
          END DO
         END DO
        END DO
        DO K=1,3
         DENOM=THREE*EIGV(K)**TWO-VECTI(1)*TWO*EIGV(K)+
     1   VECTI(2)
         DLDIDI(K,1,1)=TWO*EIGV(K)*DLDI(K,1)/DENOM-EIGV(K)**
     1   TWO*(SIX*EIGV(K)*DLDI(K,1)-TWO*VECTI(1)*DLDI(K,1)-
     2   TWO*EIGV(K))/DENOM**TWO
         DLDIDI(K,1,2)=TWO*EIGV(K)*DLDI(K,2)/DENOM-EIGV(K)**
     1   TWO*(SIX*EIGV(K)*DLDI(K,2)-TWO*VECTI(1)*DLDI(K,2)+
     2   ONE)/DENOM**TWO
         DLDIDI(K,1,3)=TWO*EIGV(K)*DLDI(K,3)/DENOM-EIGV(K)**
     1   TWO*(SIX*EIGV(K)*DLDI(K,3)-TWO*VECTI(1)*DLDI(K,3))/
     2   DENOM**TWO
         DLDIDI(K,2,1)=EIGV(K)*(SIX*EIGV(K)*DLDI(K,1)-TWO*
     1   VECTI(1)*DLDI(K,1)-TWO*EIGV(K))/DENOM**TWO-
     2   DLDI(K,1)/DENOM
         DLDIDI(K,2,2)=EIGV(K)*(SIX*EIGV(K)*DLDI(K,2)-TWO*
     1   VECTI(1)*DLDI(K,2)+ONE)/DENOM**TWO-DLDI(K,2)/DENOM
         DLDIDI(K,2,3)=EIGV(K)*(SIX*EIGV(K)*DLDI(K,3)-TWO*
     1   VECTI(1)*DLDI(K,3))/DENOM**TWO-DLDI(K,3)/DENOM
         DLDIDI(K,3,1)=MONE*(SIX*EIGV(K)*DLDI(K,1)-TWO*
     1   VECTI(1)*DLDI(K,1)-TWO*EIGV(K))/DENOM**TWO
         DLDIDI(K,3,2)=MONE*(SIX*EIGV(K)*DLDI(K,2)-TWO*
     1   VECTI(1)*DLDI(K,2)+ONE)/DENOM**TWO
         DLDIDI(K,3,3)=MONE*(SIX*EIGV(K)*DLDI(K,3)-TWO*
     1   VECTI(1)*DLDI(K,3))/DENOM**TWO
        END DO
C
C   Calculating derivatives of invariants respect to matrix components
        DO I=1,3
         DO J=1,6
          DIDE(I,J)=ZERO
         END DO
        END DO
        DIDE(1,1)=ONE
        DIDE(1,2)=ONE
        DIDE(1,3)=ONE
        DIDE(2,1)=EPS(2)+EPS(3)
        DIDE(2,2)=EPS(1)+EPS(3)
        DIDE(2,3)=EPS(2)+EPS(1)
        DIDE(2,4)=MONE*EPS(4)/TWO
        DIDE(2,5)=MONE*EPS(5)/TWO
        DIDE(2,6)=MONE*EPS(6)/TWO
        DIDE(3,1)=EPS(2)*EPS(3)-EPS(6)**TWO/FOUR
        DIDE(3,2)=EPS(1)*EPS(3)-EPS(5)**TWO/FOUR
        DIDE(3,3)=EPS(2)*EPS(1)-EPS(4)**TWO/FOUR
        DIDE(3,4)=EPS(5)*EPS(6)/FOUR-EPS(4)*EPS(3)/TWO
        DIDE(3,5)=EPS(4)*EPS(6)/FOUR-EPS(5)*EPS(2)/TWO
        DIDE(3,6)=EPS(5)*EPS(4)/FOUR-EPS(6)*EPS(1)/TWO
C
        DO I=1,6
         DO J=1,6
          DO K=1,3
           DIDEDE(K,I,J)=ZERO
          END DO
         END DO
        END DO
        DIDEDE(2,4,4)=MONE*HALF
        DIDEDE(2,5,5)=MONE*HALF
        DIDEDE(2,6,6)=MONE*HALF
        DIDEDE(2,1,2)=ONE
        DIDEDE(2,2,1)=ONE
        DIDEDE(2,1,3)=ONE
        DIDEDE(2,3,1)=ONE
        DIDEDE(2,2,3)=ONE
        DIDEDE(2,3,2)=ONE
        DIDEDE(3,1,2)=EPS(3)
        DIDEDE(3,1,3)=EPS(2)
        DIDEDE(3,1,6)=MONE*EPS(6)/TWO
        DIDEDE(3,2,3)=EPS(1)
        DIDEDE(3,2,5)=MONE*EPS(5)/TWO
        DIDEDE(3,3,4)=MONE*EPS(4)/TWO
        DIDEDE(3,4,4)=MONE*EPS(3)/TWO
        DIDEDE(3,4,5)=EPS(6)/FOUR
        DIDEDE(3,4,6)=EPS(5)/FOUR
        DIDEDE(3,5,5)=MONE*EPS(2)/TWO
        DIDEDE(3,5,6)=EPS(4)/FOUR
        DIDEDE(3,6,6)=MONE*EPS(1)/TWO
        DO I=1,6
         DO J=1,I-1
          DIDEDE(3,I,J)=DIDEDE(3,J,I)
         END DO
        END DO
C
C    Calculating derivatives of eigenvalues respect to matrix components
        DO I=1,3
         DO J=1,6
          DLDE(I,J)=ZERO
         END DO
        END DO
        DO K=1,3
         DO J=1,6
          DO I=1,3
           DLDE(K,J)=DLDE(K,J)+DLDI(K,I)*DIDE(I,J)
          END DO
         END DO
        END DO
C
        DO I=1,6
         DO J=1,6
          DO K=1,3
           DLDEDE(K,I,J)=ZERO
          END DO
         END DO
        END DO
        DO N=1,3
         DO I=1,6
          DO J=1,6
           DO M=1,3
            DLDEDE(N,I,J)=DLDEDE(N,I,J)+DLDI(N,M)*DIDEDE(M,I,J)
             DO P=1,3
              DLDEDE(N,I,J)=DLDEDE(N,I,J)+DLDIDI(N,M,P)*
     1         DIDE(P,J)*DIDE(M,I)
             END DO
            END DO
           END DO
          END DO
         END DO
C
C   Stiffness matrix
        DO I=1,6
         DO J=1,6
          DO N=1,3
           DO M=1,3
            CMAT(I,J)=CMAT(I,J)+DLDE(M,J)*CLMAT(N,M)*DLDE(N,I)
     1       +EIGV(M)*CLMAT(M,N)*DLDEDE(N,I,J)
           END DO
          END DO
         END DO
        END DO
       ENDIF
       RETURN
       END
C
       SUBROUTINE JACOBYEIG(EPS,EIGVEC,EIGVAL)
       INCLUDE 'ABA_PARAM.INC'
       PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,FOUR=4.D0,
     1  TOLER=1.0D-12,SIX=6.D0,FT=50.D0,THREE=3.D0,HALF=0.5D0,
     2  TS=27.D0,CNTM=1000,TEN=10.D0,TOLERE=1.0D-12)
       INTEGER I, J, K, N, CNT
       REAL*8 EPS(6),EIGVEC(3,3),EIGVAL(3,3)
       REAL*8 B2, BAR, BETA, COEFF, S, C, CS, SC, VAG
C
       DO I=1,3
        DO J=1,3
         EIGVEC(I,J) = ZERO
         EIGVAL(I,J) = ZERO
        END DO
       END DO
C       
       DO I=1,3
         EIGVEC(I,I) = ONE
       END DO
C       
       EIGVAL(1,1) = EPS(1)
       EIGVAL(2,2) = EPS(2)
       EIGVAL(3,3) = EPS(3)
       EIGVAL(1,2) = EPS(4)/TWO
       EIGVAL(1,3) = EPS(5)/TWO
       EIGVAL(2,3) = EPS(6)/TWO
       EIGVAL(2,1) = EIGVAL(1,2)
       EIGVAL(3,1) = EIGVAL(1,3)
       EIGVAL(3,2) = EIGVAL(2,3)
C     
       B2 = ZERO
       DO I=1,3
         DO J=1,3
           IF (I.NE.J) THEN
               B2 = B2 + EIGVAL(I,J)**TWO
           ENDIF
         END DO
       END DO
       BAR = B2/(THREE*THREE)/TWO
C      
       CNT=ONE
       DO WHILE ((B2.GT.TOLER).AND.(CNT.LT.CNTM))
          DO I=1,2
            DO J=I+1,3
              IF (EIGVAL(J,I)**TWO.LE.BAR) THEN
              ELSE
               B2 = B2 - TWO*EIGVAL(J,I)**TWO
               BAR = HALF*B2/(THREE*THREE)
               BETA = (EIGVAL(J,J)-EIGVAL(I,I))/(TWO*EIGVAL(J,I))
               COEFF = HALF*BETA/SQRT(ONE+BETA**TWO)
               S = SQRT(MAX(HALF+COEFF,ZERO))
               C = SQRT(MAX(HALF-COEFF,ZERO))
               DO K=1,3
                CS =  C*EIGVAL(I,K)+S*EIGVAL(J,K)
                SC = -S*EIGVAL(I,K)+C*EIGVAL(J,K)
                EIGVAL(I,K) = CS
                EIGVAL(J,K) = SC
               END DO
               DO K=1,3
                CS =  C*EIGVAL(K,I)+S*EIGVAL(K,J)
                SC = -S*EIGVAL(K,I)+C*EIGVAL(K,J)
                EIGVAL(K,I) = CS
                EIGVAL(K,J) = SC
                CS =  C*EIGVEC(K,I)+S*EIGVEC(K,J)
                SC = -S*EIGVEC(K,I)+C*EIGVEC(K,J)
                EIGVEC(K,I) = CS
                EIGVEC(K,J) = SC
               END DO
              ENDIF
            END DO
          END DO
       CNT=CNT+ONE
       END DO
C
C ---------- Sorting eigenvalues and vectors ---------------
        VAG=ZERO
        IF (EIGVAL(1,1).GT.EIGVAL(2,2)) THEN
            VAG=EIGVAL(1,1)
            EIGVAL(1,1)=EIGVAL(2,2)
            EIGVAL(2,2)=VAG
            DO I=1,3
                VAG=EIGVEC(I,1)
                EIGVEC(I,1)=EIGVEC(I,2)
                EIGVEC(I,2)=VAG
            END DO
        ENDIF
        IF (EIGVAL(1,1).GT.EIGVAL(3,3)) THEN
            VAG=EIGVAL(1,1)
            EIGVAL(1,1)=EIGVAL(3,3)
            EIGVAL(3,3)=VAG
            DO I=1,3
                VAG=EIGVEC(I,1)
                EIGVEC(I,1)=EIGVEC(I,3)
                EIGVEC(I,3)=VAG
            END DO
        ENDIF 
        IF (EIGVAL(2,2).GT.EIGVAL(3,3)) THEN
            VAG=EIGVAL(2,2)
            EIGVAL(2,2)=EIGVAL(3,3)
            EIGVAL(3,3)=VAG
            DO I=1,3
                VAG=EIGVEC(I,2)
                EIGVEC(I,2)=EIGVEC(I,3)
                EIGVEC(I,3)=VAG
            END DO
        ENDIF        
       RETURN
       END    
C      
C Subroutine UMAT  : 
C Dummy material
C
C ==============================================================
C !!! NOTE: N_ELEM has to be changed according to the UEL !!!!!
C ==============================================================
C
       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
       CHARACTER*80 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C 
       PARAMETER (ONE=1.0,TWO=2.0,THREE=3.0,SIX=6.0, HALF =0.5,
     1 N_ELEM=6400,NSTV=46) 
       DATA NEWTON,TOLER/40,1.D-6/ 
C       
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,8)
C 
C ----------------------------------------------------------- 
C          Material properties
C ----------------------------------------------------------- 
C          PROPS(1) - Young's modulus 
C          PROPS(2) - Poisson ratio 
C ----------------------------------------------------------- 
C
C   Elastic properties
C
       EMOD=PROPS(1)
       ENU=PROPS(2)
       EG=EMOD/(TWO*(ONE+ENU))
       EG2=EG*TWO
       ELAM=EG2*ENU/(ONE-TWO*ENU)
C
C   Stiffness tensor
C
       DO K1=1, NTENS
        DO K2=1, NTENS
         DDSDDE(K2, K1)=0.0
        END DO
       END DO
C
       DO K1=1, NDI
        DO K2=1, NDI
         DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
       END DO
C
       DO K1=NDI+1, NTENS
        DDSDDE(K1, K1)=EG
       END DO
C
C     Calculate Stresses
C
       DO K1=1, NTENS
        DO K2=1, NTENS
         STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
        END DO
       END DO
C
       NELEMAN=NOEL-TWO*N_ELEM
       IF (NPT.EQ.3) THEN
        NPT=4
       ELSEIF (NPT.EQ.4) THEN
        NPT=3
       ENDIF
       IF (NPT.EQ.7) THEN
        NPT=8
       ELSEIF (NPT.EQ.8) THEN
        NPT=7
       ENDIF
       
       DO I=1,NSTATV
        STATEV(I)=USRVAR(NELEMAN,I,NPT)
       END DO
       
       RETURN
       END      
      