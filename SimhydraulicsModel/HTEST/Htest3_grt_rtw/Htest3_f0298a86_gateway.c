/* Simscape target specific file.
 * This file is generated for the Simscape network associated with the solver block 'Htest3/Solver Configuration'.
 */

#include "nesl_rtw.h"
#include "Htest3_f0298a86_1.h"

void Htest3_f0298a86_gateway(void)
{
  NeModelParameters modelparams = { (NeSolverType) 0, 0.001, 0, 0.001, 0, 0, 0,
    0, (SscLoggingSetting) 0, };

  NeSolverParameters solverparams = { 0, 0, 1, 0, 0, 0.001, 1e-06, 1e-09, 0, 0,
    1e-09, 0, (NeAdvancerChoice) 0, 0.001, 0, 12, 2, (NeLinearAlgebraChoice) 0,
    1024, 1, 0.001, };

  const NeInputParameters* inputparameters = NULL;
  const NeOutputParameters* outputparameters = NULL;
  const NeLinearAlgebra* linear_algebra_ptr = ((solverparams.mLinearAlgebra ==
    NE_FULL_LA) ? get_rtw_linear_algebra() : neu_get_csparse_linear_algebra());
  NeDae* dae[1];
  size_t numInputs = 0;
  size_t numOutputs = 0;

  {
    static const NeInputParameters inputparameters_init[] = { { 0, 0, 0, 0.001,
        1, 0, "Htest3/Simulink-PS\nConverter", }, { 0, 0, 0, 0.001, 1, 0,
        "Htest3/Simulink-PS\nConverter1", }, { 0, 0, 0, 0.001, 1, 0,
        "Htest3/Simulink-PS\nConverter10", }, { 0, 0, 0, 0.001, 1, 0,
        "Htest3/Simulink-PS\nConverter11", }, { 0, 0, 0, 0.001, 1, 0,
        "Htest3/Simulink-PS\nConverter12", }, { 0, 0, 0, 0.001, 1, 0,
        "Htest3/Simulink-PS\nConverter13", }, { 0, 0, 0, 0.001, 1, 0,
        "Htest3/Simulink-PS\nConverter2", }, { 0, 0, 0, 0.001, 1, 0,
        "Htest3/Simulink-PS\nConverter3", }, { 0, 0, 0, 0.001, 1, 0,
        "Htest3/Simulink-PS\nConverter4", }, { 0, 0, 0, 0.001, 1, 0,
        "Htest3/Simulink-PS\nConverter5", }, { 0, 0, 0, 0.001, 1, 0,
        "Htest3/Simulink-PS\nConverter6", }, { 0, 0, 0, 0.001, 1, 0,
        "Htest3/Simulink-PS\nConverter7", }, { 0, 0, 0, 0.001, 1, 0,
        "Htest3/Simulink-PS\nConverter8", }, { 0, 0, 0, 0.001, 1, 0,
        "Htest3/Simulink-PS\nConverter9", }, };

    inputparameters = inputparameters_init;
    numInputs = sizeof(inputparameters_init)/sizeof(inputparameters_init[0]);
  }

  {
    static const NeOutputParameters outputparameters_init[] = { { 0, 0,
        "Htest3/Solver\nConfiguration", }, };

    outputparameters = outputparameters_init;
    numOutputs = sizeof(outputparameters_init)/sizeof(outputparameters_init[0]);
  }

  Htest3_f0298a86_1_dae(&dae[0],
                        &modelparams,
                        &solverparams,
                        linear_algebra_ptr);
  nesl_register_simulator_group("Htest3/Solver Configuration",
    1,
    dae,
    solverparams,
    modelparams,
    numInputs,
    inputparameters,
    numOutputs,
    outputparameters);
}
