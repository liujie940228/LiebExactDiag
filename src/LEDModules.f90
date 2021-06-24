!     Last change:  RAR  18 Sep 1999   10:54 am
!!********************************************************************
!!
!! TMSE2D - Transfer matrix method for the Anderson
!! model with diagonal disorder in two dimensions
!!
!!********************************************************************

!!********************************************************************
!!
!! $Header: /home/cvs/phsht/AML/src/AMLModules.f90,v 1.3 2007/09/25 20:45:08 phrfar Exp $
!!
!!********************************************************************

!!**************************************************************************
!!$Log: AMLModules.f90,v $
!!Revision 1.3  2007/09/25 20:45:08  phrfar
!!added IChEVal, IChEVec
!!
!!Revision 1.2  2007/09/25 10:30:33  phsht
!!rewrote to make "look" nicer and removed an ERROR for NEvals and Energy; now makefile has proper flags to check for these things automatically
!!
!!Revision 1.1  2007/09/20 16:53:39  phrfar
!!previous files from a project of Rudo's. To be used as templates.
!!
!!Revision 1.6  2007/03/21 17:58:45  phsht
!!deleted the SB stuff, included the 3D TMM, included square/stripe
!!flags, 2D seems to work, 3D still segfaults, more work needed
!!
!!Revision 1.5  2007/03/21 12:37:20  phsht
!!added variables and flags for 2D/3D, square and stripes as well as kappa and magnetic flux, now ready to start changing the main code
!!
!!Revision 1.4  2006/08/01 11:00:21  phsht
!!this version now generates .evl data according to NPrint
!!
!!Revision 1.3  2006/06/26 08:43:13  phsht
!!changes that correct some typos/errors and now seems to work with
!!GFortran
!!as well
!!
!!Revision 1.2  2006/05/08 19:23:40  phsht
!!this version writes .psi files correctly
!!
!!Revision 1.1  2006/05/08 08:42:17  phsht
!!1st installement of f90 tmse2dCOE files converted from tmse2dSB
!!
!!Revision 1.5  2004/03/17 11:38:51  phsht
!!corrected version
!!
!!Revision 1.4  2004/01/06 15:21:34  phsht
!!TMM/SB_CONVERGED included
!!
!!Revision 1.3  2003/10/17 14:20:18  phsht
!!increased MAXWidth
!!
!!Revision 1.2  2003/07/23 13:30:13  phsht
!!reorthonormalized symplectic basis
!!
!!Revision 1.2  2003/07/16 12:40:58  phsht
!!this version seems to run with restarting options in IKeepFlag
!!
!!Revision 1.1  2003/07/16 10:43:08  phsht
!!Initial revision
!!
!!**************************************************************************

!!--------------------------------------------------------------------
MODULE CConstants
  CHARACTER*18, PARAMETER :: RStr= "$Revision: 1.3 $ "
  CHARACTER*30, PARAMETER :: DStr= "$Date: 2007/09/25 20:45:08 $ "
  CHARACTER*16, PARAMETER :: AStr= "$Author: phrfar $ "
END MODULE CConstants

!! MAXGamma needs to be equal to MAXWidth, as we need to find ALL
!! Lyapunov exponents, so do not change!
MODULE IConstants
  INTEGER, PARAMETER :: MAXWidth= 1000, MAXGamma= MAXWidth, MAXIter=2147483646
  INTEGER, PARAMETER :: MAXKeepFlag= 3, MAXWriteFlag= 4, MAXFluxFlag= 3, MAXRNGFlag=2
  INTEGER, PARAMETER :: MAXSortFlag=1, MAXBCFlag=2, MAXLevelFlag=0, MAXConvFlag=3
  INTEGER, PARAMETER :: MAXStripeFlag=3, MINDimenFlag=2,MAXDimenFlag=3
  INTEGER, PARAMETER :: MAXFiles= 5, MINIter=3
END MODULE IConstants

!!--------------------------------------------------------------------
MODULE IPara
  USE MyNumbers
  INTEGER(KIND=IKIND) :: ISeed, NSeed, Seed
  INTEGER(KIND=IKIND) :: Dim, Nx
  INTEGER(KIND=IKIND) :: Width0, Width1, dWidth, IKeepFlag, IWriteFlag, ISortFlag
  INTEGER(KIND=IKIND) :: IFluxFlag, IBCFlag, IRNGFlag, ILevelFlag, IConvFlag
  INTEGER(KIND=IKIND) :: IStripeFlag, IDimenFlag
END MODULE IPara

!!--------------------------------------------------------------------
MODULE DPara
  USE MyNumbers
  REAL(KIND=RKIND) :: DiagDis0,DiagDis1,dDiagDis
  REAL(KIND=RKIND) :: DiagDis
  REAL(KIND=RKIND) :: RimDiagDis
  REAL(KIND=RKIND) :: Kappa, MagFlux
  REAL(KIND=RKIND) :: MyEpsilon

  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE :: PsiA, PsiB
END MODULE DPara

!!--------------------------------------------------------------------
!!      Input- and Outputchannels
MODULE IChannels
  USE MyNumbers
  INTEGER(KIND=IKIND), PARAMETER :: IChInp= 40, IChOut= 41, IChOutGam= 42, &
       IChOutPsi= 43, IChOutRHO= 44, IChOutAvgRHO= 45, &
       IChOutAvgRHO1= 46, IChOutAvgRHOL= 47, &
       ICHtmp= 48, IChOutHAV= 49, &
       IChEVal=50, IChEVec=51
END MODULE IChannels






