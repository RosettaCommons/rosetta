// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/evaluation.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_evaluation_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_evaluation_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace evaluation { extern BooleanOptionKey const evaluation; }
namespace evaluation { extern FileVectorOptionKey const rmsd_target; }
namespace evaluation { extern StringVectorOptionKey const rmsd_column; }
namespace evaluation { extern FileVectorOptionKey const rmsd_select; }
namespace evaluation { extern FileVectorOptionKey const align_rmsd_target; }
namespace evaluation { extern FileVectorOptionKey const structural_similarity; }
namespace evaluation { extern BooleanOptionKey const contact_map; }
namespace evaluation { extern StringVectorOptionKey const jscore_evaluator; }
namespace evaluation { extern StringVectorOptionKey const align_rmsd_column; }
namespace evaluation { extern FileVectorOptionKey const align_rmsd_fns; }
namespace evaluation { extern StringOptionKey const align_rmsd_format; }
namespace evaluation { extern StringOptionKey const predicted_burial_fn; }
namespace evaluation { extern FileOptionKey const pool; }
namespace evaluation { extern FileVectorOptionKey const rmsd; }
namespace evaluation { extern FileVectorOptionKey const chirmsd; }
namespace evaluation { extern BooleanOptionKey const gdtmm; }
namespace evaluation { extern BooleanOptionKey const gdttm; }
namespace evaluation { extern BooleanOptionKey const score_with_rmsd; }
namespace evaluation { extern FileVectorOptionKey const constraints; }
namespace evaluation { extern FileVectorOptionKey const constraints_column; }
namespace evaluation { extern FileVectorOptionKey const combined_constraints; }
namespace evaluation { extern FileVectorOptionKey const combined_constraints_column; }
namespace evaluation { extern IntegerOptionKey const combine_statistics; }
namespace evaluation { extern StringVectorOptionKey const chemical_shifts; }
namespace evaluation { extern StringOptionKey const sparta_dir; }
namespace evaluation { extern StringVectorOptionKey const cam_shifts; }
namespace evaluation { extern StringVectorOptionKey const pales; }
namespace evaluation { extern FileVectorOptionKey const extra_score; }
namespace evaluation { extern FileVectorOptionKey const extra_score_patch; }
namespace evaluation { extern StringVectorOptionKey const extra_score_column; }
namespace evaluation { extern FileVectorOptionKey const extra_score_select; }
namespace evaluation { extern FileVectorOptionKey const rdc_select; }
namespace evaluation { extern FileVectorOptionKey const rdc_target; }
namespace evaluation { extern BooleanOptionKey const symmetric_rmsd; }
namespace evaluation { extern StringVectorOptionKey const rdc_column; }
namespace evaluation { extern StringVectorOptionKey const rdc; }
namespace evaluation { extern StringOptionKey const built_in_rdc; }
namespace evaluation { extern BooleanOptionKey const jump_nr; }
namespace evaluation { extern IntegerVectorOptionKey const score_exclude_res; }
namespace evaluation { extern IntegerOptionKey const score_sscore_short_helix; }
namespace evaluation { extern IntegerOptionKey const score_sscore_maxloop; }
namespace evaluation { extern BooleanOptionKey const rpf; }
namespace evaluation { extern IntegerOptionKey const window_size; }
namespace evaluation { extern StringOptionKey const I_sc; }
namespace evaluation { extern BooleanOptionKey const Irms; }
namespace evaluation { extern BooleanOptionKey const Ca_Irms; }
namespace evaluation { extern BooleanOptionKey const Fnat; }
namespace evaluation { extern BooleanOptionKey const Lrmsd; }
namespace evaluation { extern BooleanOptionKey const Fnonnat; }
namespace evaluation { extern BooleanOptionKey const DockMetrics; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
