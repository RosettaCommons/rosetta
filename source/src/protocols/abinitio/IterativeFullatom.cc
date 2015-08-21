// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file IterativeAbrelax
/// @brief iterative protocol starting with abinitio and getting progressively more concerned with full-atom relaxed structures
/// @details
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/abinitio/IterativeFullatom.hh>
#include <protocols/jd2/archive/ArchiveManager.hh>

// Package Headers

// Project Headers
#include <core/types.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>

// Option Headers
#include <basic/Tracer.hh>
#include <basic/MemTracer.hh>

//// C++ headers
#include <cstdlib>
#include <string>

// Utility headers
#include <basic/options/option_macros.hh>

#include <core/scoring/ScoreFunction.hh>
#include <protocols/noesy_assign/NoesyModule.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.iterative" );
using basic::mem_tr;

using core::Real;
using namespace core;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_1GRP_KEY( Real, iterative, perturb_fa_resampling )
OPT_1GRP_KEY( Real, iterative, fapool_noesy_cst_weight )
OPT_1GRP_KEY( Real, iterative, fapool_chemicalshift_weight )
OPT_1GRP_KEY( Real, iterative, fapool_first_noesy_cycle_nr )

bool protocols::abinitio::IterativeFullatom::options_registered_( false );

void protocols::abinitio::IterativeFullatom::register_options() {
	IterativeBase::register_options();
	if ( !options_registered_ ) {
		NEW_OPT( iterative::perturb_fa_resampling, "perturb resample_stage2 start structures by this amount", 2.0 );
		NEW_OPT( iterative::fapool_noesy_cst_weight, "weight to apply to fullatom pool for noesy-autoassigned constraints", 5);
		NEW_OPT( iterative::fapool_chemicalshift_weight, "weight to apply to chemical shifts in centroid pool rescoring", 5 );
		NEW_OPT( iterative::fapool_first_noesy_cycle_nr, "start noesy assignment with this cycle selector", 6.0);
		options_registered_ = true;
	}
}

namespace protocols {
namespace abinitio {
using namespace jd2::archive;

IterativeFullatom::IterativeFullatom()
: IterativeBase( "fullatom_pool" )
{
	perturb_start_structures_ = option[ iterative::perturb_fa_resampling ];
}

void IterativeFullatom::initialize() {
	Parent::initialize();
	// --- setup scorefxn
	set_stage( LAST_CENTROID_START );
	set_finish_stage( FINISHED );
	test_for_stage_end();
	mem_tr << "before setup fa-score function" << std::endl;

	core::scoring::ScoreFunctionOP scorefxn =
		core::scoring::ScoreFunctionFactory::create_score_function( fa_score(), fa_score_patch() );

	// if local evaluation done by cmdline_cst evaluator
	set_scorefxn( scorefxn );

	//Base class sets chainbreak scores as convenience if not in patches..
	// this will go wrong in fullatom mode, since chainbreaks not returned ... remove here
	if ( !evaluate_local() ) {
		remove_evaluation( "linear_chainbreak");
		remove_evaluation( "overlap_chainbreak" );

		if ( noesy_assign::NoesyModule::cmdline_options_activated() ) {
			set_weight( "noesy_autoassign_cst", option[ iterative::fapool_noesy_cst_weight ]()*overall_cstfilter_weight());
		}
		if ( option[ iterative::fapool_chemicalshift_weight ].user() ) {
			set_weight( chemshift_column(), option[ iterative::fapool_chemicalshift_weight ]() );
		}
	} else { //evaluate local
		set_weight( "score", 1.0 );
		set_weight( "atom_pair_constraint", 0 ); //this is now done via FILTER mechanism of ConstraintClaimer only !
		if ( std::abs( scorefxn->get_weight( scoring::atom_pair_constraint ) - overall_cstfilter_weight() ) > 0.1 ) {
			set_overall_cstfilter_weight( scorefxn->get_weight( scoring::atom_pair_constraint ) );
			setup_filter_cst( overall_cstfilter_weight() );
		}
		scorefxn->set_weight( scoring::atom_pair_constraint, 0 );

		set_weight( "rdc", scorefxn->get_weight( scoring::rdc ) );
		scorefxn->set_weight( scoring::rdc, 0 );

		set_scorefxn( scorefxn );
	}
	set_noesy_assign_float_cycle( option[ iterative::fapool_first_noesy_cycle_nr ]() );
	if ( super_quick_relax_of_centroids() ) {
		set_weight( "score_fa", evaluate_local() ? 0.0 : 1.0 );
		set_weight( "prefa_centroid_score", 0.0 );
		set_weight( "prefa_clean_score3", 0.0 );
	}
}

bool IterativeFullatom::ready_for_batch() const {
	return ( stage() > LAST_CENTROID_START );
}

/// @details generate new batch...
/// type of batch depends on stage_. we switch to next stage based on some convergence criteria:
/// right now it is how many decoys were accepted from last batch.. if this number drops sufficiently ---> next stage...
///    (maybe need to put a safeguard in here: ratio small but at least XXX decoys proposed since last batch... )
///
core::Size IterativeFullatom::generate_batch( jd2::archive::Batch& batch, core::Size /* repeat_id */ ) {
	batch.set_intermediate_structs( false );
	cluster();
	gen_noe_assignments( batch );
	gen_resample_fragments( batch );
	gen_evaluation_output( batch, true /*fullatom*/ );
	if ( stage() == RIGID_CORE_RESAMPLING ) {
		gen_resample_core( batch, false /*flex*/ );
	}
	return 0;
}


/// ============================================================================
/// -----------           methods to make new batches               ------------
/// -----          each method should only "append" to broker and flags    -----
/// ============================================================================


void IterativeFullatom::gen_resample_core( Batch& batch, bool flex ) {
	//copy pool as input decoys
	batch.set_has_silent_in();
	io::silent::SilentFileData sfd;
	for ( SilentStructs::const_iterator it = decoys().begin(); it != decoys().end(); ++it ) {
		sfd.add_structure( *it ); //only add OP to sfd
	}
	sfd.write_all( batch.silent_in() );

	//take rigid-core definition that has most individual loops --- most jumps
	Size most_jumps( 0 ), nr_jumps( 0 );
	for ( Size i = 2; i<=4; ++i ) {
		if ( core( i ).num_loop() > nr_jumps ) {
			nr_jumps = core( i ).num_loop();
			most_jumps = i;
		}
	}

	utility::io::ozstream broker( batch.broker_file(), std::ios::app );
	broker << "\nUSE_INPUT_POSE" << std::endl;

	if ( most_jumps ) {
		core( most_jumps ).write_loops_to_file( batch.dir()+"core.rigid", "RIGID" );
		broker << "\nCLAIMER RigidChunkClaimer \n"
			<< "REGION_FILE "<< batch.dir() << "core.rigid\n"
			<< ( flex ? "KEEP_FLEXIBLE\n" : "" )
			<< "END_CLAIMER\n\n" << std::endl;
	}

	broker << "\nCLAIMER StartStructClaimer\n"
		<< "PERTURB " << perturb_start_structures_ << "\n"
		<< "END_CLAIMER\n\n" << std::endl;

	broker.close();

	utility::io::ozstream flags( batch.flag_file(), std::ios::app );
	flags << "-abinitio::skip_stages 1 2" << std::endl;

	add_fullatom_flags( batch );
	batch.nstruct() = std::max( 1, int( 1.0*batch.nstruct() / ( 1.0*decoys().size() ) ) );
}

} //abinitio
} //protocols
