/* Simscape target specific file.
 * This file is generated for the Simscape network associated with the solver block 'Htest3/Solver Configuration'.
 */

#include "ne_std.h"
#include "ne_default_allocator.h"
#include "ne_dae_fwd.h"
#include "ne_profiler_fwd.h"
#include "ne_dae_construct.h"
#include "Htest3_f0298a86_1_ds.h"

void Htest3_f0298a86_1_dae( NeDae **dae, const NeModelParameters *modelParams,
  const NeSolverParameters *solverParams,
  const NeLinearAlgebra *linear_algebra_ptr)
{
  NeAllocator *ne_allocator;
  ne_allocator = ne_default_allocator();
  ne_dae_create( dae,
                Htest3_f0298a86_1_dae_ds( ne_allocator, NULL),
                *solverParams,
                *modelParams,
                linear_algebra_ptr,
                NULL,
                NULL,
                NULL,
                ne_allocator);
}
