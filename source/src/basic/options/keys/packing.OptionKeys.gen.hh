// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/packing.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_packing_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_packing_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace packing { extern BooleanOptionKey const packing; }
namespace packing { extern BooleanOptionKey const repack_only; }
namespace packing { extern BooleanOptionKey const prevent_repacking; }
namespace packing { extern RealOptionKey const cenrot_cutoff; }
namespace packing { extern BooleanOptionKey const ignore_ligand_chi; }
namespace packing { extern IntegerOptionKey const ndruns; }
namespace packing { extern BooleanOptionKey const soft_rep_design; }
namespace packing { extern BooleanOptionKey const use_electrostatic_repulsion; }
namespace packing { extern BooleanOptionKey const dump_rotamer_sets; }
namespace packing { extern RealOptionKey const dunbrack_prob_buried; }
namespace packing { extern RealOptionKey const dunbrack_prob_nonburied; }
namespace packing { extern RealOptionKey const dunbrack_prob_nonburied_semirotameric; }
namespace packing { extern BooleanOptionKey const no_optH; }
namespace packing { extern BooleanOptionKey const optH_MCA; }
namespace packing { extern BooleanOptionKey const pack_missing_sidechains; }
namespace packing { extern BooleanOptionKey const preserve_c_beta; }
namespace packing { extern BooleanOptionKey const flip_HNQ; }
namespace packing { extern IntegerVectorOptionKey const fix_his_tautomer; }
namespace packing { extern BooleanOptionKey const print_pymol_selection; }
namespace packing { namespace ex1 { extern BooleanOptionKey const ex1; } }
namespace packing { namespace ex1 { extern IntegerOptionKey const level; } }
namespace packing { namespace ex1 { extern BooleanOptionKey const operate; } }
namespace packing { namespace ex2 { extern BooleanOptionKey const ex2; } }
namespace packing { namespace ex2 { extern IntegerOptionKey const level; } }
namespace packing { namespace ex2 { extern BooleanOptionKey const operate; } }
namespace packing { namespace ex3 { extern BooleanOptionKey const ex3; } }
namespace packing { namespace ex3 { extern IntegerOptionKey const level; } }
namespace packing { namespace ex3 { extern BooleanOptionKey const operate; } }
namespace packing { namespace ex4 { extern BooleanOptionKey const ex4; } }
namespace packing { namespace ex4 { extern IntegerOptionKey const level; } }
namespace packing { namespace ex4 { extern BooleanOptionKey const operate; } }
namespace packing { namespace ex1aro { extern BooleanOptionKey const ex1aro; } }
namespace packing { namespace ex1aro { extern IntegerOptionKey const level; } }
namespace packing { namespace ex2aro { extern BooleanOptionKey const ex2aro; } }
namespace packing { namespace ex2aro { extern IntegerOptionKey const level; } }
namespace packing { namespace ex1aro_exposed { extern BooleanOptionKey const ex1aro_exposed; } }
namespace packing { namespace ex1aro_exposed { extern IntegerOptionKey const level; } }
namespace packing { namespace ex2aro_exposed { extern BooleanOptionKey const ex2aro_exposed; } }
namespace packing { namespace ex2aro_exposed { extern IntegerOptionKey const level; } }
namespace packing { namespace exdna { extern BooleanOptionKey const exdna; } }
namespace packing { namespace exdna { extern IntegerOptionKey const level; } }
namespace packing { extern IntegerOptionKey const extrachi_cutoff; }
namespace packing { extern FileVectorOptionKey const resfile; }
namespace packing { extern RealOptionKey const outeriterations_scaling; }
namespace packing { extern RealOptionKey const inneriterations_scaling; }
namespace packing { extern BooleanOptionKey const explicit_h2o; }
namespace packing { extern StringVectorOptionKey const adducts; }
namespace packing { extern BooleanOptionKey const solvate; }
namespace packing { extern BooleanOptionKey const use_input_sc; }
namespace packing { extern FileVectorOptionKey const unboundrot; }
namespace packing { extern RealOptionKey const max_rotbump_energy; }
namespace packing { extern BooleanOptionKey const lazy_ig; }
namespace packing { extern BooleanOptionKey const double_lazy_ig; }
namespace packing { extern IntegerOptionKey const double_lazy_ig_mem_limit; }
namespace packing { extern IntegerOptionKey const linmem_ig; }
namespace packing { extern IntegerOptionKey const multi_cool_annealer; }
namespace packing { extern RealVectorOptionKey const minpack_temp_schedule; }
namespace packing { extern IntegerOptionKey const minpack_inner_iteration_scale; }
namespace packing { extern BooleanOptionKey const minpack_disable_bumpcheck; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
