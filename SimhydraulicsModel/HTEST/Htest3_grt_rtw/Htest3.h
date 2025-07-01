/*
 * Htest3.h
 *
 * Code generation for model "Htest3".
 *
 * Model version              : 1.45
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Fri Oct 03 13:20:02 2014
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#ifndef RTW_HEADER_Htest3_h_
#define RTW_HEADER_Htest3_h_
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#ifndef Htest3_COMMON_INCLUDES_
# define Htest3_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_logging.h"
#include "nesl_rtw.h"
#include "Htest3_f0298a86_gateway.h"
#endif                                 /* Htest3_COMMON_INCLUDES_ */

#include "Htest3_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "math.h"
#include "rt_matrixlib.h"
#include "rt_nonfinite.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetBlkStateChangeFlag
# define rtmGetBlkStateChangeFlag(rtm) ((rtm)->ModelData.blkStateChange)
#endif

#ifndef rtmSetBlkStateChangeFlag
# define rtmSetBlkStateChangeFlag(rtm, val) ((rtm)->ModelData.blkStateChange = (val))
#endif

#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->ModelData.contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->ModelData.contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->ModelData.contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->ModelData.contStates = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->ModelData.derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->ModelData.derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->ModelData.intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->ModelData.intgData = (val))
#endif

#ifndef rtmGetMassMatrixIr
# define rtmGetMassMatrixIr(rtm)       ((rtm)->ModelData.massMatrixIr)
#endif

#ifndef rtmSetMassMatrixIr
# define rtmSetMassMatrixIr(rtm, val)  ((rtm)->ModelData.massMatrixIr = (val))
#endif

#ifndef rtmGetMassMatrixJc
# define rtmGetMassMatrixJc(rtm)       ((rtm)->ModelData.massMatrixJc)
#endif

#ifndef rtmSetMassMatrixJc
# define rtmSetMassMatrixJc(rtm, val)  ((rtm)->ModelData.massMatrixJc = (val))
#endif

#ifndef rtmGetMassMatrixNzMax
# define rtmGetMassMatrixNzMax(rtm)    ((rtm)->ModelData.massMatrixNzMax)
#endif

#ifndef rtmSetMassMatrixNzMax
# define rtmSetMassMatrixNzMax(rtm, val) ((rtm)->ModelData.massMatrixNzMax = (val))
#endif

#ifndef rtmGetMassMatrixPr
# define rtmGetMassMatrixPr(rtm)       ((rtm)->ModelData.massMatrixPr)
#endif

#ifndef rtmSetMassMatrixPr
# define rtmSetMassMatrixPr(rtm, val)  ((rtm)->ModelData.massMatrixPr = (val))
#endif

#ifndef rtmGetMassMatrixType
# define rtmGetMassMatrixType(rtm)     ((rtm)->ModelData.massMatrixType)
#endif

#ifndef rtmSetMassMatrixType
# define rtmSetMassMatrixType(rtm, val) ((rtm)->ModelData.massMatrixType = (val))
#endif

#ifndef rtmGetOdeDELTA
# define rtmGetOdeDELTA(rtm)           ((rtm)->ModelData.odeDELTA)
#endif

#ifndef rtmSetOdeDELTA
# define rtmSetOdeDELTA(rtm, val)      ((rtm)->ModelData.odeDELTA = (val))
#endif

#ifndef rtmGetOdeDFDX
# define rtmGetOdeDFDX(rtm)            ((rtm)->ModelData.odeDFDX)
#endif

#ifndef rtmSetOdeDFDX
# define rtmSetOdeDFDX(rtm, val)       ((rtm)->ModelData.odeDFDX = (val))
#endif

#ifndef rtmGetOdeE
# define rtmGetOdeE(rtm)               ((rtm)->ModelData.odeE)
#endif

#ifndef rtmSetOdeE
# define rtmSetOdeE(rtm, val)          ((rtm)->ModelData.odeE = (val))
#endif

#ifndef rtmGetOdeEDOT
# define rtmGetOdeEDOT(rtm)            ((rtm)->ModelData.odeEDOT)
#endif

#ifndef rtmSetOdeEDOT
# define rtmSetOdeEDOT(rtm, val)       ((rtm)->ModelData.odeEDOT = (val))
#endif

#ifndef rtmGetOdeF0
# define rtmGetOdeF0(rtm)              ((rtm)->ModelData.odeF0)
#endif

#ifndef rtmSetOdeF0
# define rtmSetOdeF0(rtm, val)         ((rtm)->ModelData.odeF0 = (val))
#endif

#ifndef rtmGetOdeF1
# define rtmGetOdeF1(rtm)              ((rtm)->ModelData.odeF1)
#endif

#ifndef rtmSetOdeF1
# define rtmSetOdeF1(rtm, val)         ((rtm)->ModelData.odeF1 = (val))
#endif

#ifndef rtmGetOdeFAC
# define rtmGetOdeFAC(rtm)             ((rtm)->ModelData.odeFAC)
#endif

#ifndef rtmSetOdeFAC
# define rtmSetOdeFAC(rtm, val)        ((rtm)->ModelData.odeFAC = (val))
#endif

#ifndef rtmGetOdeFMXDOT
# define rtmGetOdeFMXDOT(rtm)          ((rtm)->ModelData.odeFMXDOT)
#endif

#ifndef rtmSetOdeFMXDOT
# define rtmSetOdeFMXDOT(rtm, val)     ((rtm)->ModelData.odeFMXDOT = (val))
#endif

#ifndef rtmGetOdeMASSMATRIX_M
# define rtmGetOdeMASSMATRIX_M(rtm)    ((rtm)->ModelData.odeMASSMATRIX_M)
#endif

#ifndef rtmSetOdeMASSMATRIX_M
# define rtmSetOdeMASSMATRIX_M(rtm, val) ((rtm)->ModelData.odeMASSMATRIX_M = (val))
#endif

#ifndef rtmGetOdeMASSMATRIX_M1
# define rtmGetOdeMASSMATRIX_M1(rtm)   ((rtm)->ModelData.odeMASSMATRIX_M1)
#endif

#ifndef rtmSetOdeMASSMATRIX_M1
# define rtmSetOdeMASSMATRIX_M1(rtm, val) ((rtm)->ModelData.odeMASSMATRIX_M1 = (val))
#endif

#ifndef rtmGetOdePIVOTS
# define rtmGetOdePIVOTS(rtm)          ((rtm)->ModelData.odePIVOTS)
#endif

#ifndef rtmSetOdePIVOTS
# define rtmSetOdePIVOTS(rtm, val)     ((rtm)->ModelData.odePIVOTS = (val))
#endif

#ifndef rtmGetOdeW
# define rtmGetOdeW(rtm)               ((rtm)->ModelData.odeW)
#endif

#ifndef rtmSetOdeW
# define rtmSetOdeW(rtm, val)          ((rtm)->ModelData.odeW = (val))
#endif

#ifndef rtmGetOdeX0
# define rtmGetOdeX0(rtm)              ((rtm)->ModelData.odeX0)
#endif

#ifndef rtmSetOdeX0
# define rtmSetOdeX0(rtm, val)         ((rtm)->ModelData.odeX0 = (val))
#endif

#ifndef rtmGetOdeX1START
# define rtmGetOdeX1START(rtm)         ((rtm)->ModelData.odeX1START)
#endif

#ifndef rtmSetOdeX1START
# define rtmSetOdeX1START(rtm, val)    ((rtm)->ModelData.odeX1START = (val))
#endif

#ifndef rtmGetOdeXDOT
# define rtmGetOdeXDOT(rtm)            ((rtm)->ModelData.odeXDOT)
#endif

#ifndef rtmSetOdeXDOT
# define rtmSetOdeXDOT(rtm, val)       ((rtm)->ModelData.odeXDOT = (val))
#endif

#ifndef rtmGetOdeXTMP
# define rtmGetOdeXTMP(rtm)            ((rtm)->ModelData.odeXTMP)
#endif

#ifndef rtmSetOdeXTMP
# define rtmSetOdeXTMP(rtm, val)       ((rtm)->ModelData.odeXTMP = (val))
#endif

#ifndef rtmGetOdeZTMP
# define rtmGetOdeZTMP(rtm)            ((rtm)->ModelData.odeZTMP)
#endif

#ifndef rtmSetOdeZTMP
# define rtmSetOdeZTMP(rtm, val)       ((rtm)->ModelData.odeZTMP = (val))
#endif

#ifndef rtmGetRTWLogInfo
# define rtmGetRTWLogInfo(rtm)         ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->ModelData.zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->ModelData.zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->ModelData.derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->ModelData.derivs = (val))
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

/* Block signals (auto storage) */
typedef struct {
  real_T dv0[1656];
  real_T dv1[1656];
  real_T EXEC_INPUT_1[4];              /* '<S20>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_l[4];            /* '<S26>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_m[4];            /* '<S24>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_f[4];            /* '<S25>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_h[4];            /* '<S29>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_g[4];            /* '<S21>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_g2[4];           /* '<S28>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_e[4];            /* '<S30>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_p[4];            /* '<S31>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_px[4];           /* '<S32>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_i[4];            /* '<S33>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_n[4];            /* '<S22>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_f2[4];           /* '<S23>/EXEC_INPUT_1' */
  real_T SineWave;                     /* '<Root>/Sine Wave' */
  real_T EXEC_INPUT_1_n4[4];           /* '<S27>/EXEC_INPUT_1' */
  real_T EXEC_STATE_1[1600];           /* '<S34>/EXEC_STATE_1' */
  real_T EXEC_OUTPUT_3[2];             /* '<S34>/EXEC_OUTPUT_3' */
} B_Htest3_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T EXEC_INPUT_1_DiscStates;      /* '<S20>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_DiscStates_k;    /* '<S26>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_DiscStates_g;    /* '<S24>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_DiscStates_o;    /* '<S25>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_DiscStates_kh;   /* '<S29>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_DiscStates_e;    /* '<S21>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_DiscStates_j;    /* '<S28>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_DiscStates_jf;   /* '<S30>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_DiscStates_i;    /* '<S31>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_DiscStates_gh;   /* '<S32>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_DiscStates_p;    /* '<S33>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_DiscStates_ir;   /* '<S22>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_DiscStates_pn;   /* '<S23>/EXEC_INPUT_1' */
  real_T EXEC_INPUT_1_DiscStates_a;    /* '<S27>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator;        /* '<S20>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData;          /* '<S20>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr;          /* '<S20>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger;           /* '<S20>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx;    /* '<S20>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator_p;      /* '<S26>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData_n;        /* '<S26>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr_e;        /* '<S26>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger_n;         /* '<S26>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx_d;  /* '<S26>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator_c;      /* '<S24>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData_l;        /* '<S24>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr_f;        /* '<S24>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger_k;         /* '<S24>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx_c;  /* '<S24>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator_a;      /* '<S25>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData_h;        /* '<S25>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr_l;        /* '<S25>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger_f;         /* '<S25>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx_j;  /* '<S25>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator_g;      /* '<S29>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData_o;        /* '<S29>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr_n;        /* '<S29>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger_i;         /* '<S29>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx_h;  /* '<S29>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator_e;      /* '<S21>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData_b;        /* '<S21>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr_g;        /* '<S21>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger_n2;        /* '<S21>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx_b;  /* '<S21>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator_l;      /* '<S28>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData_nt;       /* '<S28>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr_ei;       /* '<S28>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger_h;         /* '<S28>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx_c0; /* '<S28>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator_lc;     /* '<S30>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData_e;        /* '<S30>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr_ex;       /* '<S30>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger_nj;        /* '<S30>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx_a;  /* '<S30>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator_n;      /* '<S31>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData_j;        /* '<S31>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr_a;        /* '<S31>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger_m;         /* '<S31>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx_f;  /* '<S31>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator_b;      /* '<S32>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData_d;        /* '<S32>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr_af;       /* '<S32>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger_h0;        /* '<S32>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx_p;  /* '<S32>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator_i;      /* '<S33>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData_k;        /* '<S33>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr_li;       /* '<S33>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger_p;         /* '<S33>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx_dr; /* '<S33>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator_ei;     /* '<S22>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData_p;        /* '<S22>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr_k;        /* '<S22>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger_md;        /* '<S22>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx_l;  /* '<S22>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator_k;      /* '<S23>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData_ey;       /* '<S23>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr_j;        /* '<S23>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger_j;         /* '<S23>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx_py; /* '<S23>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Simulator_lq;     /* '<S27>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SimData_dt;       /* '<S27>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_DiagMgr_c;        /* '<S27>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_Logger_c;         /* '<S27>/EXEC_INPUT_1' */
  void* EXEC_INPUT_1_SampleTimeIdx_ce; /* '<S27>/EXEC_INPUT_1' */
  void* EXEC_STATE_1_Simulator;        /* '<S34>/EXEC_STATE_1' */
  void* EXEC_STATE_1_SimData;          /* '<S34>/EXEC_STATE_1' */
  void* EXEC_STATE_1_DiagMgr;          /* '<S34>/EXEC_STATE_1' */
  void* EXEC_STATE_1_Logger;           /* '<S34>/EXEC_STATE_1' */
  void* EXEC_STATE_1_SampleTimeIdx;    /* '<S34>/EXEC_STATE_1' */
  void* EXEC_OUTPUT_3_Simulator;       /* '<S34>/EXEC_OUTPUT_3' */
  void* EXEC_OUTPUT_3_SimData;         /* '<S34>/EXEC_OUTPUT_3' */
  void* EXEC_OUTPUT_3_DiagMgr;         /* '<S34>/EXEC_OUTPUT_3' */
  void* EXEC_OUTPUT_3_Logger;          /* '<S34>/EXEC_OUTPUT_3' */
  void* EXEC_OUTPUT_3_SampleTimeIdx;   /* '<S34>/EXEC_OUTPUT_3' */
  struct {
    void *LoggedData;
  } Extension_PWORK;                   /* '<Root>/Extension' */

  struct {
    void *LoggedData;
  } Extension1_PWORK;                  /* '<Root>/Extension1' */

  void* EXEC_SINK_2_Simulator;         /* '<S34>/EXEC_SINK_2' */
  void* EXEC_SINK_2_SimData;           /* '<S34>/EXEC_SINK_2' */
  void* EXEC_SINK_2_DiagMgr;           /* '<S34>/EXEC_SINK_2' */
  void* EXEC_SINK_2_Logger;            /* '<S34>/EXEC_SINK_2' */
  void* EXEC_SINK_2_SampleTimeIdx;     /* '<S34>/EXEC_SINK_2' */
  int_T EXEC_STATE_1_Modes[1525];      /* '<S34>/EXEC_STATE_1' */
  int32_T EXEC_STATE_1_MASS_MATRIX_PR; /* '<S34>/EXEC_STATE_1' */
} DW_Htest3_T;

/* Continuous states (auto storage) */
typedef struct {
  real_T Htest3Double_Acting_Hydraulic_C[75];/* '<S34>/EXEC_STATE_1' */
} X_Htest3_T;

/* State derivatives (auto storage) */
typedef struct {
  real_T Htest3Double_Acting_Hydraulic_C[75];/* '<S34>/EXEC_STATE_1' */
} XDot_Htest3_T;

/* State disabled  */
typedef struct {
  boolean_T Htest3Double_Acting_Hydraulic_C[75];/* '<S34>/EXEC_STATE_1' */
} XDis_Htest3_T;

/* Mass Matrix (global) */
typedef struct {
  int_T ir[26];
  int_T jc[75+1];
  real_T pr[26];
} MassMatrix_Htest3_T;

#ifndef ODE14X_INTG
#define ODE14X_INTG

/* ODE14X Integration Data */
typedef struct {
  real_T *x0;
  real_T *f0;
  real_T *x1start;
  real_T *f1;
  real_T *Delta;
  real_T *E;
  real_T *fac;
  real_T *DFDX;
  real_T *W;
  int_T *pivots;
  real_T *xtmp;
  real_T *ztmp;
  real_T *M;
  real_T *M1;
  real_T *Edot;
  real_T *xdot;
  real_T *fminusMxdot;
  boolean_T isFirstStep;
} ODE14X_IntgData;

#endif

/* Parameters (auto storage) */
struct P_Htest3_T_ {
  real_T Constant_Value;               /* Expression: 188
                                        * Referenced by: '<Root>/Constant'
                                        */
  real_T Constant1_Value;              /* Expression: 1
                                        * Referenced by: '<Root>/Constant1'
                                        */
  real_T Constant10_Value;             /* Expression: 0
                                        * Referenced by: '<Root>/Constant10'
                                        */
  real_T Constant11_Value;             /* Expression: 1
                                        * Referenced by: '<Root>/Constant11'
                                        */
  real_T Constant12_Value;             /* Expression: 1
                                        * Referenced by: '<Root>/Constant12'
                                        */
  real_T Constant2_Value;              /* Expression: 1
                                        * Referenced by: '<Root>/Constant2'
                                        */
  real_T Constant3_Value;              /* Expression: 0
                                        * Referenced by: '<Root>/Constant3'
                                        */
  real_T Constant4_Value;              /* Expression: 1
                                        * Referenced by: '<Root>/Constant4'
                                        */
  real_T Constant5_Value;              /* Expression: 1
                                        * Referenced by: '<Root>/Constant5'
                                        */
  real_T Constant6_Value;              /* Expression: 1
                                        * Referenced by: '<Root>/Constant6'
                                        */
  real_T Constant7_Value;              /* Expression: 1
                                        * Referenced by: '<Root>/Constant7'
                                        */
  real_T Constant8_Value;              /* Expression: 1
                                        * Referenced by: '<Root>/Constant8'
                                        */
  real_T Constant9_Value;              /* Expression: 1
                                        * Referenced by: '<Root>/Constant9'
                                        */
  real_T SineWave_Amp;                 /* Expression: 3e-3
                                        * Referenced by: '<Root>/Sine Wave'
                                        */
  real_T SineWave_Bias;                /* Expression: 0
                                        * Referenced by: '<Root>/Sine Wave'
                                        */
  real_T SineWave_Freq;                /* Expression: 2
                                        * Referenced by: '<Root>/Sine Wave'
                                        */
  real_T SineWave_Phase;               /* Expression: 0
                                        * Referenced by: '<Root>/Sine Wave'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_Htest3_T {
  const char_T *errorStatus;
  RTWLogInfo *rtwLogInfo;
  RTWSolverInfo solverInfo;

  /*
   * ModelData:
   * The following substructure contains information regarding
   * the data used in the model.
   */
  struct {
    X_Htest3_T *contStates;
    real_T *derivs;
    boolean_T *contStateDisabled;
    boolean_T zCCacheNeedsReset;
    boolean_T derivCacheNeedsReset;
    boolean_T blkStateChange;
    int_T massMatrixType;
    int_T massMatrixNzMax;
    int_T *massMatrixIr;
    int_T *massMatrixJc;
    real_T *massMatrixPr;
    real_T odeX0[75];
    real_T odeF0[75];
    real_T odeX1START[75];
    real_T odeF1[75];
    real_T odeDELTA[75];
    real_T odeE[4*75];
    real_T odeFAC[75];
    real_T odeDFDX[75*75];
    real_T odeW[75*75];
    int_T odePIVOTS[75];
    real_T odeXTMP[75];
    real_T odeZTMP[75];
    real_T odeMASSMATRIX_M[26];
    real_T odeMASSMATRIX_M1[26];
    real_T odeEDOT[4*75];
    real_T odeXDOT[75];
    real_T odeFMXDOT[75];
    ODE14X_IntgData intgData;
  } ModelData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T tFinal;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block parameters (auto storage) */
extern P_Htest3_T Htest3_P;

/* Block signals (auto storage) */
extern B_Htest3_T Htest3_B;

/* Continuous states (auto storage) */
extern X_Htest3_T Htest3_X;

/* global MassMatrix */
extern MassMatrix_Htest3_T Htest3_MassMatrix;

/* Block states (auto storage) */
extern DW_Htest3_T Htest3_DW;

/* Model entry point functions */
extern void Htest3_initialize(void);
extern void Htest3_step(void);
extern void Htest3_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Htest3_T *const Htest3_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'Htest3'
 * '<S1>'   : 'Htest3/PS-Simulink Converter'
 * '<S2>'   : 'Htest3/PS-Simulink Converter1'
 * '<S3>'   : 'Htest3/Simulink-PS Converter'
 * '<S4>'   : 'Htest3/Simulink-PS Converter1'
 * '<S5>'   : 'Htest3/Simulink-PS Converter10'
 * '<S6>'   : 'Htest3/Simulink-PS Converter11'
 * '<S7>'   : 'Htest3/Simulink-PS Converter12'
 * '<S8>'   : 'Htest3/Simulink-PS Converter13'
 * '<S9>'   : 'Htest3/Simulink-PS Converter2'
 * '<S10>'  : 'Htest3/Simulink-PS Converter3'
 * '<S11>'  : 'Htest3/Simulink-PS Converter4'
 * '<S12>'  : 'Htest3/Simulink-PS Converter5'
 * '<S13>'  : 'Htest3/Simulink-PS Converter6'
 * '<S14>'  : 'Htest3/Simulink-PS Converter7'
 * '<S15>'  : 'Htest3/Simulink-PS Converter8'
 * '<S16>'  : 'Htest3/Simulink-PS Converter9'
 * '<S17>'  : 'Htest3/Solver Configuration'
 * '<S18>'  : 'Htest3/PS-Simulink Converter/EVAL_KEY'
 * '<S19>'  : 'Htest3/PS-Simulink Converter1/EVAL_KEY'
 * '<S20>'  : 'Htest3/Simulink-PS Converter/EVAL_KEY'
 * '<S21>'  : 'Htest3/Simulink-PS Converter1/EVAL_KEY'
 * '<S22>'  : 'Htest3/Simulink-PS Converter10/EVAL_KEY'
 * '<S23>'  : 'Htest3/Simulink-PS Converter11/EVAL_KEY'
 * '<S24>'  : 'Htest3/Simulink-PS Converter12/EVAL_KEY'
 * '<S25>'  : 'Htest3/Simulink-PS Converter13/EVAL_KEY'
 * '<S26>'  : 'Htest3/Simulink-PS Converter2/EVAL_KEY'
 * '<S27>'  : 'Htest3/Simulink-PS Converter3/EVAL_KEY'
 * '<S28>'  : 'Htest3/Simulink-PS Converter4/EVAL_KEY'
 * '<S29>'  : 'Htest3/Simulink-PS Converter5/EVAL_KEY'
 * '<S30>'  : 'Htest3/Simulink-PS Converter6/EVAL_KEY'
 * '<S31>'  : 'Htest3/Simulink-PS Converter7/EVAL_KEY'
 * '<S32>'  : 'Htest3/Simulink-PS Converter8/EVAL_KEY'
 * '<S33>'  : 'Htest3/Simulink-PS Converter9/EVAL_KEY'
 * '<S34>'  : 'Htest3/Solver Configuration/EVAL_KEY'
 */
#endif                                 /* RTW_HEADER_Htest3_h_ */
