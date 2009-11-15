#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <assert.h>

static void
set_group_tags (int const checkpoint, int const prolongate,
                char const * restrict const gn);

int
ML_BSSN6_SetGroupTags (void)
{
  DECLARE_CCTK_PARAMETERS;
  
  set_group_tags (0, 1, "ADMBase::metric");
  set_group_tags (0, 1, "ADMBase::curv");
  set_group_tags (0, 1, "ADMBase::lapse");
  set_group_tags (0, 1, "ADMBase::shift");
  set_group_tags (0, 1, "ADMBase::dtlapse");
  set_group_tags (0, 1, "ADMBase::dtshift");
  
  set_group_tags (0, 0, "ML_BSSN6::ML_cons_detg");
  set_group_tags (0, 0, "ML_BSSN6::ML_cons_Gamma");
  set_group_tags (0, 0, "ML_BSSN6::ML_cons_traceA");
  set_group_tags (0, 0, "ML_BSSN6::ML_Ham");
  set_group_tags (0, 0, "ML_BSSN6::ML_mom");
  set_group_tags (0, 0, "ML_BSSN6::ML_curvrhs");
  
  int const checkpoint = rhs_timelevels > 1;
  set_group_tags (checkpoint, 0, "ML_BSSN6::ML_dtlapserhs");
  set_group_tags (checkpoint, 0, "ML_BSSN6::ML_dtshiftrhs");
  set_group_tags (checkpoint, 0, "ML_BSSN6::ML_Gammarhs");
  set_group_tags (checkpoint, 0, "ML_BSSN6::ML_lapserhs");
  set_group_tags (checkpoint, 0, "ML_BSSN6::ML_log_confacrhs");
  set_group_tags (checkpoint, 0, "ML_BSSN6::ML_metricrhs");
  set_group_tags (checkpoint, 0, "ML_BSSN6::ML_shiftrhs");
  set_group_tags (checkpoint, 0, "ML_BSSN6::ML_trace_curvrhs");
  
  return 0;
}

static void
set_group_tags (int const checkpoint, int const prolongate,
                char const * restrict const gn)
{
  assert (gn);
  
  int const gi = CCTK_GroupIndex (gn);
  assert (gi >= 0);
  
  int const table = CCTK_GroupTagsTableI (gi);
  assert (table >= 0);
  
  if (! checkpoint) {
    int ierr;
    ierr = Util_TableSetString (table, "no", "Checkpoint");
    assert (! ierr);
    
    ierr = Util_TableSetString (table, "no", "Persistent");
    assert (! ierr);
  }
  
  if (! prolongate) {
    int ierr;
    ierr = Util_TableSetString (table, "none", "Prolongation");
    assert (! ierr);
  }
}
