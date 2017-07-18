// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief      Mover partners apart and relax them separately
/// @details Run quick relax on separated partners; this emulates unbound docking
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_docking_membrane_QuickRelaxPartnersSeparately_cc
#define INCLUDED_protocols_docking_membrane_QuickRelaxPartnersSeparately_cc

// Unit Headers
#include <protocols/docking/membrane/QuickRelaxPartnersSeparately.hh>
#include <protocols/docking/membrane/QuickRelaxPartnersSeparatelyCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/MPQuickRelaxMover.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/util.hh>
#include <protocols/docking/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <protocols/jd2/util.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <protocols/docking/metrics.hh>
#include <protocols/membrane/util.hh>
#include <numeric/random/random.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.membrane.QuickRelaxPartnersSeparately" );

namespace protocols {
namespace docking {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::moves;
using namespace protocols::docking;
using namespace protocols::membrane;
using namespace protocols::membrane::geometry;

/////////////////////
/// Constructors  ///
/////////////////////


/// @brief Default Constructor
/// @details Gets the jump from docking partners
QuickRelaxPartnersSeparately::QuickRelaxPartnersSeparately() :
	protocols::moves::Mover(),
	topo_( new SpanningTopology() ),
	topo_up_( new SpanningTopology() ),
	topo_down_( new SpanningTopology() )
{}

/// @brief Copy Constructor
QuickRelaxPartnersSeparately::QuickRelaxPartnersSeparately( QuickRelaxPartnersSeparately const & src ) :
	protocols::moves::Mover( src ),
	native_( src.native_ ),
	partners_( src.partners_ ),
	jump_( src.jump_ ),
	jumps_( src.jumps_ ),
	topo_( src.topo_ ),
	topo_up_( src.topo_up_ ),
	topo_down_( src.topo_down_ ),
	sfxn_( src.sfxn_ )
{}

/// @brief Destructor
QuickRelaxPartnersSeparately::~QuickRelaxPartnersSeparately() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
QuickRelaxPartnersSeparately::clone() const {
	return ( protocols::moves::MoverOP( new QuickRelaxPartnersSeparately( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
QuickRelaxPartnersSeparately::fresh_instance() const {
	return protocols::moves::MoverOP( new QuickRelaxPartnersSeparately() );
}

/// @brief Parse Rosetta Scripts Options for this Mover
void
QuickRelaxPartnersSeparately::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	// Read in docking partners
	if ( tag->hasOption( "partners" ) ) {
		partners_ = tag->getOption< std::string >( "partners" );
	}

	// TODO: to implement

}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
QuickRelaxPartnersSeparatelyCreator::create_mover() const {
	return protocols::moves::MoverOP( new QuickRelaxPartnersSeparately() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
QuickRelaxPartnersSeparatelyCreator::keyname() const {
	return QuickRelaxPartnersSeparatelyCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
QuickRelaxPartnersSeparatelyCreator::mover_name() {
	return "QuickRelaxPartnersSeparately";
}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (QuickRelaxPartnersSeparately)
std::string
QuickRelaxPartnersSeparately::get_name() const {
	return "QuickRelaxPartnersSeparately";
}

/// @brief Moving partners apart and relax them separately
void
QuickRelaxPartnersSeparately::apply( Pose & pose ) {

	using namespace core::conformation::membrane;
	using namespace protocols::membrane::geometry;
	using namespace protocols::membrane;
	using namespace protocols::rigid;
	using namespace protocols::docking;
	using namespace protocols::simple_moves;

	TR << "Moving partners apart and quick relaxing them separately ..." << std::endl;

	// initializations
	register_options();
	init_from_cmd();
	finalize_setup( pose );

	////////////////////// MOVE PARTNERS APART //////////////////////////

	TR << "Moving partners apart..." << std::endl;

	// get membrane axis from docking metrics function
	core::Vector slide_axis = membrane_axis( pose, jump_ );
	bool vary_stepsize = false;

	// compute embedding for partners (compute structure-based embedding with split topologies)
	EmbeddingDefOP emb_up( compute_structure_based_embedding( pose, *topo_up_ ) );
	EmbeddingDefOP emb_down( compute_structure_based_embedding( pose, *topo_down_ ) );

	// get distance between points
	core::Real dist1 = ( emb_down->center() - emb_up->center() ).length();
	TR << "distance between partners: " << dist1 << std::endl;

	slide_axis.negate();
	slide_axis.normalize();
	TR << "slide axis: " << slide_axis.to_string() << std::endl;

	// move apart
	TR << "Moving apart" << std::endl;
	RigidBodyTransMoverOP mover( new rigid::RigidBodyTransMover( slide_axis, jump_, vary_stepsize ) );
	mover->step_size( 100 );
	mover->apply( pose );

	///////////////////////// RUN RELAX ////////////////////////////

	// run quick relax
	MPQuickRelaxMoverOP relax( new MPQuickRelaxMover() );
	relax->add_membrane_again( false );
	relax->membrane_from_topology( true );
	relax->apply( pose );

	/////////////////// SLIDE BACK TOGETHER //////////////////////////////

	slide_axis.negate();
	TR << "Moving back together" << std::endl;
	TR << "slide axis: " << slide_axis.to_string() << std::endl;

	// move apart
	RigidBodyTransMoverOP together( new rigid::RigidBodyTransMover( slide_axis, jump_, vary_stepsize ) );
	together->step_size( 100 );
	together->apply( pose );

	////// SUPERIMPOSE THE FIRST PARTNER OF THE POSE WITH THE NATIVE ///////

	TR << "foldtree before superimposition: " << std::endl;
	pose.fold_tree().show( TR );

	// get partner 1: ALL CHAINS MUST BE IN PDB CONSECUTIVELY!!!
	utility::vector1< std::string > partners( utility::string_split( partners_, '_' ) );
	utility::vector1< std::string > partner1( utility::split( partners[1] ) );

	// get residue range for superposition: get start residue
	Size start(0);
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( start == 0 &&
				partner1[1] == utility::to_string( pose.pdb_info()->chain( i ) ) ) {
			start = i;
		}
	}

	// get end residue
	Size end(0);
	for ( Size j = pose.size(); j >= 1; --j ) {
		if ( end == 0 &&
				partner1[partner1.size()] == utility::to_string( pose.pdb_info()->chain( j ) ) ) {
			end = j;
		}
	}
	TR << "superimposing from: " << start << " to " << end << std::endl;

	// superimpose partner 1 with starting pose for easier rmsd calculation
	SuperimposeMoverOP super( new SuperimposeMover( native_, start, end, start, end, true ) );
	super->apply( pose );

	TR << "dumping pose..." << std::endl;
	pose.dump_pdb("superimposed_pose.pdb");

	// calculate and store the rmsds in the score file
	protocols::jd2::add_string_real_pair_to_current_job( "Lrms", protocols::docking::calc_Lrmsd( pose, native_, jumps_ ) );
	protocols::jd2::add_string_real_pair_to_current_job( "P1rms", protocols::docking::calc_P1rmsd( pose, native_, jumps_ ) );

}// apply

////////////////////////////////////////////////////////////////////////////////

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options with JD2
void
QuickRelaxPartnersSeparately::register_options() {

	using namespace basic::options;

	option.add_relevant( OptionKeys::in::file::native );
	option.add_relevant( OptionKeys::docking::partners );

}

////////////////////////////////////////////////////////////////////////////////

/// @brief Initialize Mover options from the comandline
void
QuickRelaxPartnersSeparately::init_from_cmd() {

	using namespace basic::options;

	// docking partners
	if ( option[ OptionKeys::docking::partners ].user() ) {
		partners_ = option[ OptionKeys::docking::partners ]();
	}

	// native
	if ( option[ OptionKeys::in::file::native ].user() ) {
		core::import_pose::pose_from_file( native_, option[ OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
	}

}// init from cmd

////////////////////////////////////////////////////////////////////////////////

/// @brief Finalize setup
void QuickRelaxPartnersSeparately::finalize_setup( Pose & pose ) {

	using namespace utility;
	using namespace basic::options;
	TR << "Finalizing setup... " << std::endl;

	// setup scorefunction
	sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );

	// call AddMembraneMover on native for RMSD calculation
	AddMembraneMoverOP addmem( new AddMembraneMover() );
	addmem->apply( native_ );

	///////////////////// FOLDTREE STUFF ///////////////////////

	// get foldtree from partners (setup_foldtree) and movable jump
	// we are using this function to add the membrane add the anchor point
	// closest to COM
	setup_foldtree( pose, partners_, jumps_);
	jump_ = jumps_[1];
	TR << "jump_ from foldtree: " << jump_ << std::endl;

	// get anchor point for membrane from jump; I think this is residue closest to
	// COM of the upstream partner
	int anchor = pose.fold_tree().upstream_jump_residue( jump_ );
	TR << "anchor point: " << anchor << std::endl;

	// Add Membrane to pose, appends MEM as jump1
	AddMembraneMoverOP add_memb( new AddMembraneMover( anchor, 0 ) );
	add_memb->apply( pose );

	// reorder foldtree to have membrane at root
	core::kinematics::FoldTree ft = pose.fold_tree();
	ft.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
	pose.fold_tree( ft );
	TR << "reordered foltree: " << std::endl;
	pose.fold_tree().show( TR );
	TR << "jump: " << jump_ << std::endl;

	///////////////////// TOPOLOGY ///////////////////////

	// SpanningTopology objects
	topo_ = pose.conformation().membrane_info()->spanning_topology();
	topo_->show();
	TR << "jump number: " << jump_ << std::endl;

	// splitting topology by jump into upstream and downstream topology
	split_topology_by_jump_noshift( pose, jump_, topo_, topo_up_, topo_down_ );

} // finalize setup


} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_QuickRelaxPartnersSeparately_cc
