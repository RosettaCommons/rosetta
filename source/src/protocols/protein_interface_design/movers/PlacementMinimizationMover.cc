// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/PlacementMinimizationMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/PlacementMinimizationMover.hh>
#include <protocols/protein_interface_design/movers/PlacementMinimizationMoverCreator.hh>

// Package headers
#include <protocols/protein_interface_design/movers/BuildAlaPose.hh>
#include <protocols/protein_interface_design/movers/PlaceUtils.hh>
#include <basic/datacache/DataMap.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <basic/options/option.hh>
// Unit Headers

#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AA.hh>
#include <numeric/xyzVector.hh>

// Unit Headers
#include <protocols/protein_interface_design/design_utils.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

// C++ headers
#include <map>
#include <boost/foreach.hpp>

#include <core/scoring/constraints/Constraint.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>
#include <protocols/simple_moves/DesignRepackMover.hh>

using namespace protocols::protein_interface_design;

static thread_local basic::Tracer TR( "protocols.protein_interface_design.movers.PlacementMinimizationMover" );

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace protocols::moves;
using namespace core;

std::string
PlacementMinimizationMoverCreator::keyname() const
{
	return PlacementMinimizationMoverCreator::mover_name();
}

protocols::moves::MoverOP
PlacementMinimizationMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PlacementMinimizationMover );
}

std::string
PlacementMinimizationMoverCreator::mover_name()
{
	return "PlacementMinimization";
}

protocols::moves::MoverOP
PlacementMinimizationMover::clone() const {
	return( protocols::moves::MoverOP( new PlacementMinimizationMover( *this ) ) );
}

void
PlacementMinimizationMover::refresh_bbstub_constraints( core::pose::Pose & pose )
{
	using namespace protocols::hotspot_hashing;
	using namespace core::pack::task;
	using namespace core::kinematics;

//setting up a deflt foldtree for the constraints to work with
	FoldTree const orig_foldtree( pose.fold_tree() );
	FoldTree new_ft;
	new_ft.clear();
	new_ft.add_edge( pose.conformation().chain_begin( 1 ), pose.conformation().chain_end( 1 ), Edge::PEPTIDE );
	new_ft.add_edge( pose.conformation().chain_begin( 2 ), pose.conformation().chain_end( 2 ), Edge::PEPTIDE );
	new_ft.add_edge( pose.conformation().chain_end( 1 ), pose.conformation().chain_begin( 2 ), 1 );
	new_ft.reorder(1);
	pose.fold_tree( new_ft );
	remove_hotspot_constraints_from_pose( pose );
	core::Size fixed_res(1);
	if( host_chain_ == 1 ) fixed_res = pose.total_residue();
	core::id::AtomID const fixed_atom_id = core::id::AtomID( pose.residue(fixed_res).atom_index("CA"), fixed_res );
	HotspotStubSetOP all_stubs( new HotspotStubSet );
	BOOST_FOREACH( StubSetStubPos const stubset_pos_pair, stub_sets_ )
		all_stubs->add_stub_set( *stubset_pos_pair.first );

	PackerTaskOP restricted_packer_task;
	if( task_factory() )
		restricted_packer_task = task_factory()->create_task_and_apply_taskoperations( pose );
	else
		restricted_packer_task = TaskFactory::create_packer_task( pose );
	core::pack::task::PackerTaskOP stub_task = all_stubs->prepare_hashing_packer_task_( pose, host_chain_ );
	core::Size const host_chain_begin( pose.conformation().chain_begin( host_chain_ ) );
	core::Size const host_chain_end  ( pose.conformation().chain_end  ( host_chain_ ) );

	for( core::Size resi( host_chain_begin ); resi<=host_chain_end; ++resi ){
		using namespace core::chemical;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		if( std::find( prevent_repacking_.begin(), prevent_repacking_.end(), resi ) != prevent_repacking_.end() ||
			  !restricted_packer_task->nonconst_residue_task( resi ).being_packed() )
			stub_task->nonconst_residue_task( resi ).prevent_repacking();
		if( (pose.residue( resi ).aa() == aa_gly && !option[hotspot::allow_gly]() ) || ( pose.residue( resi ).aa() == aa_pro && !option[hotspot::allow_proline ] ))
			stub_task->nonconst_residue_task( resi ).prevent_repacking();
	}

	runtime_assert( cb_force_ > -0.000001 );
	all_stubs->add_hotspot_constraints_to_pose( pose, fixed_atom_id, stub_task, all_stubs, cb_force_, 0/*worst allowed stub bonus*/, false/*apply self energies*/, 10.0/*bump cutoff*/, true/*apply ambiguous constraints*/ );

	core::scoring::constraints::ConstraintCOPs stub_constraints( pose.add_constraints( all_stubs->constraints() ) );
 	core::Size const constraint_num( stub_constraints.size() );
 	TR<<"adding "<<constraint_num<<" stub constraints to pose"<<std::endl;
	pose.fold_tree( orig_foldtree );
}

void
PlacementMinimizationMover::apply( core::pose::Pose & pose )
{
  core::pose::Pose const saved_pose( pose ); // the pose should not actually be changed within this function

	using namespace protocols::hotspot_hashing;
	using namespace core::scoring;

	ScoreFunctionOP bbcst_scorefxn( new ScoreFunction );
	bbcst_scorefxn->reset();
	bbcst_scorefxn->set_weight( backbone_stub_constraint, 1.0 );
	(*bbcst_scorefxn)( pose );

	// Remove old hotspot constraints from pose
	remove_hotspot_constraints_from_pose( pose );
	refresh_bbstub_constraints( pose );
	// Switch to Ala unless we are doing place scaffold as a replacement for docking
	BuildAlaPose toAla( host_chain_ == 1/*partner1*/, host_chain_ == 2 /*partner2*/ );
	utility::vector1< core::Size > no_repack;
	if( !prevent_repacking().empty() ) no_repack = prevent_repacking();
	if( !no_repack.empty() ){
		std::sort( no_repack.begin(), no_repack.end() );
		std::unique( no_repack.begin(), no_repack.end() );
		toAla.prevent_repacking( no_repack );
	}
	toAla.task_factory( task_factory() );
	TR<<"switching interface to alanine\n";
	toAla.apply( pose );

	setup_packer_and_movemap( pose ); // for min_sc below
	core::scoring::ScoreFunctionCOP stub_scorefxn( make_stub_scorefxn() );
	//for minimization (rb and sc of previous placed stubs)
	utility::vector1< bool > const no_min( pose.total_residue(), false );

	utility::vector1< core::Size > no_targets;
	MinimizeInterface( pose, stub_scorefxn, no_min/*bb*/, curr_min_sc_, min_rb(), optimize_foldtree(), no_targets, true/*simultaneous optimization*/);
	TR.flush();
}

void
PlacementMinimizationMover::cb_force( core::Real const cf )
{
	cb_force_ = cf;
}

void
PlacementMinimizationMover::host_chain( core::Size const hc ){
	host_chain_ = hc;
}

void
PlacementMinimizationMover::stub_sets( utility::vector1< StubSetStubPos > const & stub_sets ){
	stub_sets_ = stub_sets;
}

void
PlacementMinimizationMover::parse_my_tag( TagCOP const tag,
		basic::datacache::DataMap &data,
		protocols::filters::Filters_map const &,
		Movers_map const &,
		core::pose::Pose const & pose )
{
	host_chain( tag->getOption<core::Size>( "host_chain", 2 ) );
	cb_force( tag->getOption< core::Real >( "cb_force", 0.5 ) );
	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );
/// PlaceSim calls this parse_my_tag with its own tag, and there, cb_force is
/// set within a child tag
	BOOST_FOREACH( TagCOP const btag, branch_tags ){
		if( btag->hasOption( "cb_force" ) )
			cb_force( btag->getOption< core::Real >( "cb_force" ) );
	}
	runtime_assert( cb_force_ >= -0.000001 );
	design_partner1_ = host_chain_ == 1 ? true : false;
	design_partner2_ = host_chain_ == 2 ? true : false;
	repack_partner1_ = true;
	repack_partner2_ = true;
	optimize_foldtree_ = tag->getOption< bool >( "optimize_foldtree", false );
	min_rb( tag->getOption< bool >( "minimize_rb", 1 ));
	stub_sets( parse_stub_sets( tag, pose, host_chain_, data ) );
	runtime_assert( stub_sets_.size() );
	TR<<"optimize_foldtree set to: "<<optimize_foldtree_<<'\n';
	TR<<"cb_force is set to "<<cb_force_<<std::endl;
}

protocols::moves::MoverOP
PlacementMinimizationMover::fresh_instance() const {
   return protocols::moves::MoverOP( new PlacementMinimizationMover );
}

PlacementMinimizationMover::~PlacementMinimizationMover(){}

PlacementMinimizationMover::PlacementMinimizationMover() :
	simple_moves::DesignRepackMover( PlacementMinimizationMoverCreator::mover_name() ),
	host_chain_( 2 ),
	cb_force_( 0.0 )
{}

std::string
PlacementMinimizationMover::get_name() const {
	return PlacementMinimizationMoverCreator::mover_name();
}

} //movers
} //protein_interface_design
} //protocols

