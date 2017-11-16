// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/abinitio/util.cc
/// @author Andy Watkins
// keep these headers first for compilation with Visual Studio C++

// Package Headers
#include <core/io/silent/SilentStruct.hh>
#include <protocols/abinitio/AbrelaxMover.hh>

// Project Headers
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
// AUTO-REMOVED #include <core/scoring/constraints/NamedAtomPairConstraint.hh>

#include <core/scoring/constraints/BoundConstraint.hh>
#include <protocols/constraints_additional/AdditionalConstraintCreators.hh>
#include <core/scoring/func/FuncFactory.hh>

//archive headers
#include <protocols/abinitio/IterativeAbrelax.hh>
// AUTO-REMOVED #include <numeric/random/random.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/loopfcst.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>
#include <basic/options/keys/chunk.OptionKeys.gen.hh>
#include <basic/options/keys/broker.OptionKeys.gen.hh>

// C/C++ headers
#include <iostream>

#include <utility/vector1.hh>


static basic::Tracer tr( "main" );

namespace protocols  {
namespace abinitio {

void register_common_options() {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	option.add_relevant(in::file::native);
	option.add_relevant(in::file::silent);
	option.add_relevant(in::file::frag3);
	option.add_relevant(in::file::frag9);
	option.add_relevant(in::file::fasta);
	option.add_relevant(in::file::native_exclude_res);
	option.add_relevant(in::file::tags);

	option.add_relevant(out::file::silent);
	option.add_relevant(out::nstruct);

	option.add_relevant(run::proc_id);
	option.add_relevant(run::nproc);
	option.add_relevant(run::condor);

	option.add_relevant(OptionKeys::abinitio::fastrelax);
	option.add_relevant(OptionKeys::abinitio::relax);
	option.add_relevant(OptionKeys::abinitio::multifastrelax);
	option.add_relevant(OptionKeys::abinitio::relax_with_jumps);
	option.add_relevant(OptionKeys::abinitio::use_filters);
	option.add_relevant(OptionKeys::abinitio::detect_disulfide_before_relax);
	option.add_relevant(OptionKeys::abinitio::debug);
	option.add_relevant(OptionKeys::abinitio::number_3mer_frags);
	option.add_relevant(OptionKeys::abinitio::number_9mer_frags);
	option.add_relevant(OptionKeys::abinitio::process_store);
	option.add_relevant(OptionKeys::abinitio::fix_residues_to_native);
	option.add_relevant(OptionKeys::abinitio::return_full_atom);
	option.add_relevant(OptionKeys::abinitio::rerun);
	option.add_relevant(OptionKeys::abinitio::jdist_rerun);
	option.add_relevant(OptionKeys::abinitio::cyclic_peptide);

	// starting conditions
	option.add_relevant(OptionKeys::abinitio::start_native);
	option.add_relevant(OptionKeys::abinitio::perturb);
	option.add_relevant(OptionKeys::abinitio::close_loops);

	// evaluation
	option.add_relevant(OptionKeys::abinitio::rmsd_residues);
	option.add_relevant(OptionKeys::abinitio::bGDT);
	option.add_relevant(OptionKeys::run::no_prof_info_in_silentout);

	// use fragments from native structure
	option.add_relevant(OptionKeys::abinitio::steal_3mers);
	option.add_relevant(OptionKeys::abinitio::steal_9mers);
	option.add_relevant(OptionKeys::abinitio::dump_frags);
	option.add_relevant(OptionKeys::abinitio::no_write_failures);

	option.add_relevant(loopfcst::use_general_protocol);
	option.add_relevant(loopfcst::coord_cst_weight);
	option.add_relevant(loopfcst::coord_cst_all_atom);
	option.add_relevant(loopfcst::coord_cst_weight_array);
	option.add_relevant(loopfcst::dump_coord_cst_weight_array);

	option.add_relevant(OptionKeys::in::file::pca);
	option.add_relevant(OptionKeys::out::sf);

	// jumping
	option.add_relevant(jumps::fix_jumps);
	option.add_relevant(jumps::jump_lib);
	option.add_relevant(jumps::fix_chainbreak);
	option.add_relevant(jumps::pairing_file);
	option.add_relevant(jumps::sheets);
	option.add_relevant(jumps::random_sheets);
	option.add_relevant(jumps::evaluate);
	option.add_relevant(jumps::extra_frags_for_ss);
	option.add_relevant(jumps::loop_definition_from_file);
	option.add_relevant(jumps::no_chainbreak_in_relax);
	option.add_relevant(jumps::residue_pair_jump_file);
	option.add_relevant(jumps::topology_file);

	//loop closure
	option.add_relevant(OptionKeys::loops::loop_file);
	option.add_relevant(OptionKeys::loops::alternative_closure_protocol);
	option.add_relevant(OptionKeys::loops::short_frag_cycles);
	option.add_relevant(OptionKeys::loops::scored_frag_cycles);
	option.add_relevant(OptionKeys::loops::debug_loop_closure);
	option.add_relevant(OptionKeys::loops::non_ideal_loop_closing);
	option.add_relevant(OptionKeys::loops::chainbreak_max_accept);
	option.add_relevant(OptionKeys::loops::extended);

	// constraints
	option.add_relevant(constraints::cst_file);
	option.add_relevant(constraints::forest_file);
	option.add_relevant(constraints::compute_total_dist_cst);
	option.add_relevant(constraints::no_linearize_bounded);
	option.add_relevant(constraints::dump_cst_set);
	option.add_relevant(constraints::no_cst_in_relax);
	option.add_relevant(constraints::evaluate_max_seq_sep);
	option.add_relevant(constraints::cull_with_native);
	option.add_relevant(constraints::named);
	option.add_relevant(constraints::viol);
	option.add_relevant(constraints::viol_level);
	option.add_relevant(constraints::viol_type);

	// homologs
	option.add_relevant(templates::config);
	option.add_relevant(templates::pairings);

	//large default number means all frags are used if this option is not specified
	option.add_relevant(templates::min_nr_large_frags);
	option.add_relevant(templates::min_nr_small_frags);

	option.add_relevant(templates::nr_large_copies);
	option.add_relevant(templates::nr_small_copies);
	option.add_relevant(templates::vary_frag_size);
	option.add_relevant(templates::fix_aligned_residues);
	option.add_relevant(templates::fix_margin);
	option.add_relevant(templates::fix_frag_file);
	option.add_relevant(templates::no_pick_fragments);
	option.add_relevant(templates::pick_multiple_sizes);
	option.add_relevant(templates::strand_constraint);

	option.add_relevant(frags::nr_large_copies);
	option.add_relevant(frags::annotate);

}

}
}
