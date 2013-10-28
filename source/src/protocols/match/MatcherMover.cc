// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/MatcherMover.cc
/// @brief  Implementation of mover wrapper for matcher
/// @author Florian Richter (floric@u.washington.edu), june 2010


// Unit headers
#include <protocols/match/MatcherMover.hh>
#include <protocols/match/MatcherMoverCreator.hh>

//package headers
#include <protocols/match/Matcher.hh>
#include <protocols/match/MatcherTask.hh>
#include <protocols/match/output/ProcessorFactory.hh>
#include <protocols/match/output/MatchProcessor.hh>
#include <protocols/match/output/PDBWriter.hh>

//project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/match.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif


namespace protocols {
namespace match {

static basic::Tracer tr( "protocols.match.MatcherMover" );


std::string
MatcherMoverCreator::keyname() const
{
	return MatcherMoverCreator::mover_name();
}

protocols::moves::MoverOP
MatcherMoverCreator::create_mover() const {
	return new MatcherMover;
}

std::string
MatcherMoverCreator::mover_name()
{
	return "MatcherMover";
}

MatcherMover::MatcherMover( bool incorporate_matches_into_pose ):
	Mover( "MatcherMover" ),
	incorporate_matches_into_pose_( incorporate_matches_into_pose ),
	ligres_(NULL)
{
	//we need this for the output to be correct
	basic::options::option[basic::options::OptionKeys::run::preserve_header ].value(true);
}

MatcherMover::~MatcherMover(){}

MatcherMover::MatcherMover( MatcherMover const & rval ) :
	//utility::pointer::ReferenceCount(),
	Mover( rval ),
	incorporate_matches_into_pose_( rval.incorporate_matches_into_pose_ ),
	ligres_( rval.ligres_ ),
	match_positions_( rval.match_positions_ )
{}

/// @brief clone this object
MatcherMover::MoverOP MatcherMover::clone() const
{
	return new MatcherMover( *this );
}

/// @brief create this type of object
MatcherMover::MoverOP MatcherMover::fresh_instance() const
{
	return new MatcherMover();
}

///
void
MatcherMover::apply( core::pose::Pose & pose )
{

	protocols::match::MatcherTaskOP mtask = new protocols::match::MatcherTask;

	core::pose::Pose ligpose;
	core::pose::Pose save_pose = pose;
	if( !ligres_ ) ligres_ = core::conformation::ResidueFactory::create_residue(
		core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map(
			basic::options::option[ basic::options::OptionKeys::match::lig_name ] ) );

	if( !ligres_->type().is_ligand() ) std::cerr << "WARNING: downstream residue " << ligres_->type().name3() << " set in the matcher mover does not seem to be a ligand residue, matcher will likely not behave properly." << std::endl;

	ligpose.append_residue_by_jump( *ligres_, 1 );

	//we might have to remove the downstream pose from the input
	if( incorporate_matches_into_pose_ ){
		for( core::Size i = pose.total_residue(); i >  0; --i ){
			if( pose.residue_type( i ).name3() == ligres_->type().name3() )	pose.conformation().delete_residue_slow( i );
		}
	}

	//if ( option[ OptionKeys::match::ligand_rotamer_index ].user() ) {
	//	set_ligpose_rotamer( ligpose );
	//}
	Size cent, nbr1, nbr2;
	ligres_->select_orient_atoms( cent, nbr1, nbr2 );

	tr << "Selecting orientation atoms:";
	tr << " " << ligres_->atom_name( cent );
	tr << " " << ligres_->atom_name( nbr1 );
	tr << " " << ligres_->atom_name( nbr2 ) << std::endl;
	//pose.dump_pdb("pose_going_into_matcher.pdb");
	mtask->set_upstream_pose( pose );

	utility::vector1< core::id::AtomID > oats( 3 );
	oats[ 1 ] = core::id::AtomID( nbr2, 1 ); oats[ 2 ] = core::id::AtomID( nbr1, 1 ); oats[ 3 ] = core::id::AtomID( cent, 1 );

	mtask->set_downstream_pose( ligpose,  oats );
	if( match_positions_.size() != 0 ){
		mtask->set_ignore_cmdline_for_build_points( true );
		mtask->set_original_scaffold_build_points( match_positions_ );
	}

	mtask->initialize_from_command_line();

	if( incorporate_matches_into_pose_ ) mtask->output_writer_name("PoseMatchOutputWriter");

	time_t matcher_start_time = time(NULL);
	protocols::match::MatcherOP matcher = new protocols::match::Matcher;
	matcher->initialize_from_task( *mtask );


	protocols::match::output::MatchProcessorOP processor = protocols::match::output::ProcessorFactory::create_processor( matcher, mtask );

	if( matcher->find_hits() ){
		matcher->process_matches( *processor );
	}

	time_t matcher_end_time = time(NULL);
	std::string success_str( processor->match_processing_successful() ? "successful." : "not sucessful." );
	tr << "Matcher ran for " << (long)(matcher_end_time - matcher_start_time) << " seconds and was " << success_str << std::endl;

	if( processor->match_processing_successful() ) this->set_last_move_status( protocols::moves::MS_SUCCESS );
	else{
		this->set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
		pose = save_pose;
		return;
	}

	if( incorporate_matches_into_pose_ ){
		protocols::match::output::PoseMatchOutputWriterOP outputter( static_cast< protocols::match::output::PoseMatchOutputWriter * >( processor->output_writer().get() ) );
		outputter->insert_match_into_pose( pose );
		// should make another mover to set these constraints  flo and nobu
		//protocols::enzdes::AddOrRemoveMatchCsts addcsts = protocols::enzdes::AddOrRemoveMatchCsts();
		//addcsts.set_cst_action(protocols::enzdes::ADD_NEW);
		//addcsts.apply(pose);
	}
	//pose.dump_pdb("pose_coming_from_matcher.pdb");

} //MatcherMover::apply function

std::string
MatcherMover::get_name() const
{
	return "MatcherMover";
}

void
MatcherMover::set_ligres(
	core::conformation::ResidueCOP ligres )
{
	ligres_ = new core::conformation::Residue(*ligres);
}

void
MatcherMover::set_match_positions(
	utility::vector1< core::Size > const & match_positions )
{
	match_positions_ = match_positions;
}

void
MatcherMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	incorporate_matches_into_pose_ = tag->getOption<bool>( "incorporate_matches_into_pose", 1 );
}

}
}
