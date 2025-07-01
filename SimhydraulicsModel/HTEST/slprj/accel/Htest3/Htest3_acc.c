#include "__cf_Htest3.h"
#include <math.h>
#include "Htest3_acc.h"
#include "Htest3_acc_private.h"
#include <stdio.h>
#include "simstruc.h"
#include "fixedpoint.h"
#define CodeFormat S-Function
#define AccDefine1 Accelerator_S-Function
static void mdlOutputs ( SimStruct * S , int_T tid ) { n3qi1whofz * _rtB ;
loikxjbxjg * _rtP ; ew10rzwqr2 * _rtDW ; _rtDW = ( ( ew10rzwqr2 * )
ssGetRootDWork ( S ) ) ; _rtP = ( ( loikxjbxjg * ) ssGetDefaultParam ( S ) )
; _rtB = ( ( n3qi1whofz * ) _ssGetBlockIO ( S ) ) ; if ( ssIsSampleHit ( S ,
1 , 0 ) ) { _rtB -> lpdpr4nb21 = _rtP -> P_0 ; } ssCallAccelRunBlock ( S , 0
, 1 , SS_CALL_MDL_OUTPUTS ) ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB ->
cdngl0iimw = _rtP -> P_1 ; } ssCallAccelRunBlock ( S , 0 , 4 ,
SS_CALL_MDL_OUTPUTS ) ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB ->
gkassldj01 = _rtP -> P_2 ; } ssCallAccelRunBlock ( S , 0 , 7 ,
SS_CALL_MDL_OUTPUTS ) ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB ->
mmklnhvead = _rtP -> P_3 ; } ssCallAccelRunBlock ( S , 0 , 10 ,
SS_CALL_MDL_OUTPUTS ) ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB ->
hty41mqnn3 = _rtP -> P_4 ; } ssCallAccelRunBlock ( S , 0 , 13 ,
SS_CALL_MDL_OUTPUTS ) ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB ->
af0a5ugqp2 = _rtP -> P_5 ; } ssCallAccelRunBlock ( S , 0 , 16 ,
SS_CALL_MDL_OUTPUTS ) ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB ->
pt0tcb1ftk = _rtP -> P_6 ; } ssCallAccelRunBlock ( S , 0 , 19 ,
SS_CALL_MDL_OUTPUTS ) ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB ->
gdktlb2pco = _rtP -> P_7 ; } ssCallAccelRunBlock ( S , 0 , 22 ,
SS_CALL_MDL_OUTPUTS ) ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB ->
d4t5u1rtmg = _rtP -> P_8 ; } ssCallAccelRunBlock ( S , 0 , 25 ,
SS_CALL_MDL_OUTPUTS ) ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB ->
p1xyblousu = _rtP -> P_9 ; } ssCallAccelRunBlock ( S , 0 , 28 ,
SS_CALL_MDL_OUTPUTS ) ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB ->
j02uxeb3pc = _rtP -> P_10 ; } ssCallAccelRunBlock ( S , 0 , 31 ,
SS_CALL_MDL_OUTPUTS ) ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB ->
ftiunazdag = _rtP -> P_11 ; } ssCallAccelRunBlock ( S , 0 , 34 ,
SS_CALL_MDL_OUTPUTS ) ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB ->
fl5l50wmyi = _rtP -> P_12 ; } ssCallAccelRunBlock ( S , 0 , 37 ,
SS_CALL_MDL_OUTPUTS ) ; _rtB -> kcfrx03ybf = muDoubleScalarSin ( _rtP -> P_15
* ssGetTaskTime ( S , 0 ) + _rtP -> P_16 ) * _rtP -> P_13 + _rtP -> P_14 ;
ssCallAccelRunBlock ( S , 0 , 40 , SS_CALL_MDL_OUTPUTS ) ;
ssCallAccelRunBlock ( S , 0 , 41 , SS_CALL_MDL_OUTPUTS ) ;
ssCallAccelRunBlock ( S , 0 , 42 , SS_CALL_MDL_OUTPUTS ) ; if ( ssIsSampleHit
( S , 1 , 0 ) ) { ssCallAccelRunBlock ( S , 0 , 43 , SS_CALL_MDL_OUTPUTS ) ;
ssCallAccelRunBlock ( S , 0 , 45 , SS_CALL_MDL_OUTPUTS ) ; }
ssCallAccelRunBlock ( S , 0 , 48 , SS_CALL_MDL_OUTPUTS ) ; UNUSED_PARAMETER (
tid ) ; }
#define MDL_UPDATE
static void mdlUpdate ( SimStruct * S , int_T tid ) { ssCallAccelRunBlock ( S
, 0 , 41 , SS_CALL_MDL_UPDATE ) ; UNUSED_PARAMETER ( tid ) ; }
#define MDL_DERIVATIVES
static void mdlDerivatives ( SimStruct * S ) { ssCallAccelRunBlock ( S , 0 ,
41 , SS_CALL_MDL_DERIVATIVES ) ; }
#define MDL_FORCINGFUNCTION
static void mdlForcingFunction ( SimStruct * S ) { ssCallAccelRunBlock ( S ,
0 , 41 , SS_CALL_MDL_FORCINGFUNCTION ) ; }
#define MDL_MASSMATRIX
static void mdlMassMatrix ( SimStruct * S ) { ssCallAccelRunBlock ( S , 0 ,
41 , SS_CALL_MDL_MASSMATRIX ) ; } static void mdlInitializeSizes ( SimStruct
* S ) { ssSetChecksumVal ( S , 0 , 2450949343U ) ; ssSetChecksumVal ( S , 1 ,
3490669517U ) ; ssSetChecksumVal ( S , 2 , 185591479U ) ; ssSetChecksumVal (
S , 3 , 327628992U ) ; { mxArray * slVerStructMat = NULL ; mxArray * slStrMat
= mxCreateString ( "simulink" ) ; char slVerChar [ 10 ] ; int status =
mexCallMATLAB ( 1 , & slVerStructMat , 1 , & slStrMat , "ver" ) ; if ( status
== 0 ) { mxArray * slVerMat = mxGetField ( slVerStructMat , 0 , "Version" ) ;
if ( slVerMat == NULL ) { status = 1 ; } else { status = mxGetString (
slVerMat , slVerChar , 10 ) ; } } mxDestroyArray ( slStrMat ) ;
mxDestroyArray ( slVerStructMat ) ; if ( ( status == 1 ) || ( strcmp (
slVerChar , "8.3" ) != 0 ) ) { return ; } } ssSetOptions ( S ,
SS_OPTION_EXCEPTION_FREE_CODE ) ; if ( ssGetSizeofDWork ( S ) != sizeof (
ew10rzwqr2 ) ) { ssSetErrorStatus ( S ,
"Unexpected error: Internal DWork sizes do "
"not match for accelerator mex file." ) ; } if ( ssGetSizeofGlobalBlockIO ( S
) != sizeof ( n3qi1whofz ) ) { ssSetErrorStatus ( S ,
"Unexpected error: Internal BlockIO sizes do "
"not match for accelerator mex file." ) ; } { int ssSizeofParams ;
ssGetSizeofParams ( S , & ssSizeofParams ) ; if ( ssSizeofParams != sizeof (
loikxjbxjg ) ) { static char msg [ 256 ] ; sprintf ( msg ,
"Unexpected error: Internal Parameters sizes do "
"not match for accelerator mex file." ) ; } } _ssSetDefaultParam ( S , (
real_T * ) & o2iu0a2jke ) ; if ( ssGetSizeofDWork ( S ) == sizeof (
ew10rzwqr2 ) ) { { ( ( ew10rzwqr2 * ) ssGetRootDWork ( S ) ) -> anikra2yw1 =
0 ; } } } static void mdlInitializeSampleTimes ( SimStruct * S ) { } static
void mdlTerminate ( SimStruct * S ) { }
#include "simulink.c"
