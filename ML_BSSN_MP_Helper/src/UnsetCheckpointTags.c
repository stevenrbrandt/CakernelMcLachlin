#include <cctk.h>
#include <util_Table.h>

#include <assert.h>

static void
unset_checkpoint_tag (char const * restrict gn);

int
ML_BSSN_MP_UnsetCheckpointTags (void)
{
  unset_checkpoint_tag ("ADMBase::metric");
  unset_checkpoint_tag ("ADMBase::curv");
  unset_checkpoint_tag ("ADMBase::lapse");
  unset_checkpoint_tag ("ADMBase::shift");
  unset_checkpoint_tag ("ADMBase::dtlapse");
  unset_checkpoint_tag ("ADMBase::dtshift");
  
  unset_checkpoint_tag ("ML_BSSN_MP::ML_cons_detg");
  unset_checkpoint_tag ("ML_BSSN_MP::ML_cons_Gamma");
  unset_checkpoint_tag ("ML_BSSN_MP::ML_cons_traceA");
  unset_checkpoint_tag ("ML_BSSN_MP::ML_Ham");
  unset_checkpoint_tag ("ML_BSSN_MP::ML_mom");
  unset_checkpoint_tag ("ML_BSSN_MP::ML_curvrhs");
  unset_checkpoint_tag ("ML_BSSN_MP::ML_dtlapserhs");
  unset_checkpoint_tag ("ML_BSSN_MP::ML_dtshiftrhs");
  unset_checkpoint_tag ("ML_BSSN_MP::ML_Gammarhs");
  unset_checkpoint_tag ("ML_BSSN_MP::ML_lapserhs");
  unset_checkpoint_tag ("ML_BSSN_MP::ML_log_confacrhs");
  unset_checkpoint_tag ("ML_BSSN_MP::ML_metricrhs");
  unset_checkpoint_tag ("ML_BSSN_MP::ML_shiftrhs");
  unset_checkpoint_tag ("ML_BSSN_MP::ML_trace_curvrhs");
  
  return 0;
}

static void
unset_checkpoint_tag (char const * restrict const gn)
{
  assert (gn);
  
  int const gi = CCTK_GroupIndex (gn);
  assert (gi >= 0);
  
  int const table = CCTK_GroupTagsTableI (gi);
  assert (table >= 0);
  
  int ierr;
  ierr = Util_TableSetString (table, "no", "Checkpoint");
  assert (! ierr);
  
  ierr = Util_TableSetString (table, "none", "Prolongation");
  assert (! ierr);
  
  ierr = Util_TableSetString (table, "no", "Persistent");
  assert (! ierr);
}
