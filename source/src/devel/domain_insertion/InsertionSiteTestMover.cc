// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/domain_insertion/InsertionSiteTestMover.cc
/// @brief  cc file for InsertionSiteTestMover
/// @author Florian Richter, flosopher@gmail.com, february 2013


// Unit headers
#include <devel/domain_insertion/InsertionSiteTestMover.hh>
#include <devel/domain_insertion/InsertionSiteTestMoverCreator.hh>

// package headers

// Project headers
#include <basic/Tracer.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

#include <devel/enzdes/EnzdesRemodelProtocol.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>


namespace devel {
namespace domain_insertion {

static basic::Tracer tr( "devel.domain_insertion.InsertionSiteTestMover" );

InsertionSiteTestMover::InsertionSiteTestMover()
:
	sfxn_(NULL), flex_window_(2),
	test_insert_ss_("LLLLLLLL"), insert_allowed_score_increase_(10.0),
	enz_flexbb_prot_( new protocols::enzdes::EnzdesFlexBBProtocol() )
	//mostly arbitrary numbers, but can be modified through RS tag
{
	insert_test_pos_.clear();
}

InsertionSiteTestMover::InsertionSiteTestMover( InsertionSiteTestMover const & other )
: parent( other ),
	sfxn_(other.sfxn_), insert_test_pos_(other.insert_test_pos_),
	flex_window_(other.flex_window_), test_insert_ss_(other.test_insert_ss_),
	 insert_allowed_score_increase_(other.insert_allowed_score_increase_), enz_flexbb_prot_( new protocols::enzdes::EnzdesFlexBBProtocol() )
{}


InsertionSiteTestMover::~InsertionSiteTestMover(){}

protocols::moves::MoverOP
InsertionSiteTestMover::clone() const{
	return new InsertionSiteTestMover( *this );
}


std::string
InsertionSiteTestMover::get_name() const
{
	return "InsertionSiteTestMover";
}


/// @details
void
InsertionSiteTestMover::apply( core::pose::Pose & pose )
{

	core::pose::Pose input_pose = pose;
	core::Real input_score = (*sfxn_)(input_pose);

	tr << "instest apply called, insert_test_pos_ has size " << insert_test_pos_.size() << std::endl;

	for( Size i = 1; i <= insert_test_pos_.size(); ++i){

		Size insert_pos = insert_test_pos_[ i ];
		tr << "Starting insertion test at position " << insert_pos << "..." << std::endl;

		enz_flexbb_prot_->add_flexible_region( insert_pos - 1, insert_pos + 2, pose, true );

		devel::enzdes::EnzdesRemodelMoverOP enzremodel_mover(make_enzremodel_mover( pose, insert_pos ) );

		enzremodel_mover->apply( pose );

		core::Real diff_score = (*sfxn_)(pose) - input_score;
		tr << "Test insertion at pos " << insert_pos << " produced a score diff of " << diff_score << std::endl;
		pose.dump_pdb("test_insertion_pos"+utility::to_string(insert_pos)+".pdb" );
	}

} //apply function

void
InsertionSiteTestMover::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose
)
{
	insert_test_pos_.clear();

	sfxn_ = protocols::rosetta_scripts::parse_score_function( tag, "sfxn", data )->clone();
	enz_flexbb_prot_->set_scorefxn( sfxn_ );

	if( tag->hasOption("test_insert_ss") ){
		test_insert_ss_ = tag->getOption<std::string>( "test_insert_ss");
	}

	if( tag->hasOption("allowed_sc_increase") ){
		insert_allowed_score_increase_ = tag->getOption<Real>( "allowed_sc_increase");
	}

	if( tag->hasOption("seqpos")){
		insert_test_pos_.push_back( tag->getOption<Size>("seqpos") );
	}

	utility::vector0< utility::tag::TagCOP > const subtags( tag->getTags() );
	for( utility::vector0< utility::tag::TagCOP >::const_iterator it= subtags.begin(); it!=subtags.end(); ++it ) {

		utility::tag::TagCOP const subtag = *it;
		if( subtag->getName() == "span" ) {
			core::Size const begin( core::pose::get_resnum( subtag, pose, "begin_" ) );
			core::Size const end( core::pose::get_resnum( subtag, pose, "end_" ) );
			runtime_assert( end > begin );
			runtime_assert( begin>=1);
			//runtime_assert( end<=reference_pose_->total_residue() );
			for( core::Size i=begin; i<=end; ++i ) insert_test_pos_.push_back( i );
		}
	}

	tr << "InsertionSiteTestMover parse_my_tag: trying insert string " << test_insert_ss_ << " at positions: ";
	for( Size i = 1; i <= insert_test_pos_.size(); ++i ) tr << insert_test_pos_[i] << ", ";
	tr << std::endl;
}


core::pack::task::PackerTaskCOP
InsertionSiteTestMover::make_insert_task(
	core::pose::Pose const & pose,
	Size insert_pos
) const
{
	core::pack::task::PackerTaskOP task = new core::pack::task::PackerTask_( pose );

	utility::vector1< bool > repack_pos( pose.total_residue(), false );
	repack_pos[ insert_pos ] = true;
	repack_pos[ insert_pos - 1] = true;
	repack_pos[ insert_pos + 1] = true;
	repack_pos[ insert_pos + 2] = true;

	enz_flexbb_prot_->get_tenA_neighbor_residues( pose, repack_pos );

	for( Size i = 1; i <= pose.total_residue(); ++i ){
		if( repack_pos[ i ] ) task->nonconst_residue_task( i ).restrict_to_repacking();
		else  task->nonconst_residue_task( i ).prevent_repacking();
	}
	return task;
}

devel::enzdes::EnzdesRemodelMoverOP
InsertionSiteTestMover::make_enzremodel_mover(
	core::pose::Pose & pose,
	Size insert_pos
) const
{

	core::pack::task::PackerTaskCOP insert_task( make_insert_task( pose, insert_pos ));

	devel::enzdes::EnzdesRemodelMoverOP enzremodel_mover = new devel::enzdes::EnzdesRemodelMover( enz_flexbb_prot_, insert_task, enz_flexbb_prot_->enz_flexible_region( 1 ) );

	std::string ss_edges;
	for( core::Size i = 1; i <= flex_window_; ++i ) ss_edges = ss_edges + "L";
	utility::vector1< std::string > ss_strings;
	ss_strings.push_back( ss_edges + test_insert_ss_ + ss_edges );
	enzremodel_mover->set_user_provided_ss( ss_strings );
	enzremodel_mover->set_max_allowed_score_increase(  insert_allowed_score_increase_ );

	return enzremodel_mover;
}


std::string
InsertionSiteTestMoverCreator::keyname() const
{
	return InsertionSiteTestMoverCreator::mover_name();
}

protocols::moves::MoverOP
InsertionSiteTestMoverCreator::create_mover() const {
	return new InsertionSiteTestMover;
}

std::string
InsertionSiteTestMoverCreator::mover_name()
{
	return "InsertionSiteTestMover";
}



}
}
