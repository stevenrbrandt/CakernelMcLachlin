#include <cctk.h>
#include <util_Table.h>

#include <assert.h>

static void
unset_checkpoint_tag (int prolongate, char const * restrict gn);

int
ML_BSSNUp_UnsetCheckpointTags (void)
{
  unset_checkpoint_tag (1, "ADMBase::metric");
  unset_checkpoint_tag (1, "ADMBase::curv");
  unset_checkpoint_tag (1, "ADMBase::lapse");
  unset_checkpoint_tag (1, "ADMBase::shift");
  unset_checkpoint_tag (1, "ADMBase::dtlapse");
  unset_checkpoint_tag (1, "ADMBase::dtshift");
  
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_cons_detg");
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_cons_Gamma");
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_cons_traceA");
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_Ham");
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_mom");
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_curvrhs");
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_dtlapserhs");
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_dtshiftrhs");
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_Gammarhs");
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_lapserhs");
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_log_confacrhs");
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_metricrhs");
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_shiftrhs");
  unset_checkpoint_tag (0, "ML_BSSNUp::ML_trace_curvrhs");
  
  return 0;
}

static void
unset_checkpoint_tag (int const prolongate, char const * restrict const gn)
{
  assert (gn);
  
  int const gi = CCTK_GroupIndex (gn);
  assert (gi >= 0);
  
  int const table = CCTK_GroupTagsTableI (gi);
  assert (table >= 0);
  
  int ierr;
  ierr = Util_TableSetString (table, "no", "Checkpoint");
  assert (! ierr);
  
  ierr = Util_TableSetString (table, "no", "Persistent");
  assert (! ierr);
  
  if (! prolongate) {
    ierr = Util_TableSetString (table, "none", "Prolongation");
    assert (! ierr);
  }
}
