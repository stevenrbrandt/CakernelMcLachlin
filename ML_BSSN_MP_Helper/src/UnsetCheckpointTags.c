#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <assert.h>

static void
unset_checkpoint_tag (char const * restrict gn);

void
ML_BSSN_MP_UnsetCheckpointTags (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  unset_checkpoint_tag ("ADMBase::metric");
  unset_checkpoint_tag ("ADMBase::curv");
  unset_checkpoint_tag ("ADMBase::lapse");
  unset_checkpoint_tag ("ADMBase::shift");
  unset_checkpoint_tag ("ADMBase::dtlapse");
  unset_checkpoint_tag ("ADMBase::dtshift");
}

static void
unset_checkpoint_tag (char const * restrict const gn)
{
  assert (gn);
  
  int const gi = CCTK_GroupIndex (gn);
  assert (gi >= 0);
  
  int const table = CCTK_GroupTagsTableI (gi);
  assert (table >= 0);
  
  int const ierr = Util_TableSetString (table, "no", "Checkpoint");
  assert (! ierr);
}
