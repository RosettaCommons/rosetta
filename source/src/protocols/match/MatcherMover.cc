// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/pack/rotamers/SingleLigandRotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/match.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <fstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


namespace protocols {
namespace match {

static THREAD_LOCAL basic::Tracer tr( "protocols.match.MatcherMover" );

// XRW TEMP std::string
// XRW TEMP MatcherMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return MatcherMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP MatcherMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new MatcherMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP MatcherMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "MatcherMover";
// XRW TEMP }

MatcherMover::MatcherMover( bool incorporate_matches_into_pose ):
	protocols::rosetta_scripts::MultiplePoseMover(),
	incorporate_matches_into_pose_( incorporate_matches_into_pose ),
	return_single_random_match_( false ),
	ligres_(/* NULL */),
	selectors_()
{
	//we need this for the output to be correct
	basic::options::option[ basic::options::OptionKeys::run::preserve_header ].value(true);
}

MatcherMover::~MatcherMover() = default;

MatcherMover::MatcherMover( MatcherMover const & ) = default;

/// @brief clone this object
MatcherMover::MoverOP MatcherMover::clone() const
{
	return MatcherMover::MoverOP( new MatcherMover( *this ) );
}

/// @brief create this type of object
MatcherMover::MoverOP MatcherMover::fresh_instance() const
{
	return MatcherMover::MoverOP( new MatcherMover() );
}

void
MatcherMover::set_return_single_random_match( bool const single_random )
{
	return_single_random_match_ = single_random;
}

void
MatcherMover::apply( core::pose::Pose & pose )
{
	if ( !pose.pdb_info() ) {
		pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( pose ) ) );
	}
	protocols::rosetta_scripts::MultiplePoseMover::apply( pose );
}

bool
MatcherMover::process_pose( core::pose::Pose & pose, utility::vector1 < core::pose::PoseOP > & poselist )
{
	protocols::match::MatcherTaskOP mtask( new protocols::match::MatcherTask );

	core::pose::Pose ligpose;
	core::pose::Pose save_pose = pose;
	if ( !ligres_ ) {
		ligres_ = core::conformation::ResidueFactory::create_residue(
			pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t )->name_map(
			basic::options::option[ basic::options::OptionKeys::match::lig_name ] ) );
	}

	if ( !ligres_->type().is_ligand() ) tr.Error << "downstream residue " << ligres_->type().name3() << " set in the matcher mover does not seem to be a ligand residue, matcher will likely not behave properly." << std::endl;

	ligpose.append_residue_by_jump( *ligres_, 1 );

	if ( basic::options::option[ basic::options::OptionKeys::match::ligand_rotamer_index ].user() ) {
		set_ligpose_rotamer( ligpose );
	}

	//we might have to remove the downstream pose from the input
	if ( incorporate_matches_into_pose_ ) {
		for ( core::Size i = pose.size(); i >  0; --i ) {
			if ( pose.residue_type( i ).name3() == ligres_->type().name3() ) pose.conformation().delete_residue_slow( i );
		}
	}

	Size cent, nbr1, nbr2;
	ligres_->select_orient_atoms( cent, nbr1, nbr2 );

	tr << "Selecting orientation atoms:";
	tr << " " << ligres_->atom_name( cent );
	tr << " " << ligres_->atom_name( nbr1 );
	tr << " " << ligres_->atom_name( nbr2 ) << std::endl;
	core::scoring::get_score_function()->setup_for_scoring( pose );
	mtask->set_upstream_pose( pose );

	utility::vector1< core::id::AtomID > oats( 3 );
	oats[ 1 ] = core::id::AtomID( nbr2, 1 ); oats[ 2 ] = core::id::AtomID( nbr1, 1 ); oats[ 3 ] = core::id::AtomID( cent, 1 );

	core::scoring::get_score_function()->setup_for_scoring( ligpose );
	mtask->set_downstream_pose( ligpose,  oats );
	if ( match_positions_.size() != 0 ) {
		mtask->set_ignore_cmdline_for_build_points( true );
		mtask->set_original_scaffold_build_points( match_positions_ );
	}


	// if selectors are set, use them to generate a match.pos file
	if ( selectors_.size() ) {
		setup_seqpos_from_selectors( *mtask, pose );
	}
	mtask->initialize_from_command_line();

	if ( incorporate_matches_into_pose_ ) mtask->output_writer_name("PoseMatchOutputWriter");

	time_t const matcher_start_time = time(nullptr);
	protocols::match::MatcherOP matcher( new protocols::match::Matcher );
	matcher->initialize_from_task( *mtask );

	protocols::match::output::MatchProcessorOP processor = protocols::match::output::ProcessorFactory::create_processor( matcher, mtask );

	time_t find_hits_end_time = 0;
	time_t processing_time = 0;
	if ( matcher->find_hits() ) {
		find_hits_end_time = time(nullptr);
		time_t process_start_time( time(nullptr) );
		matcher->process_matches( *processor );
		processing_time = (long) (time(nullptr) - process_start_time);

	} else {
		find_hits_end_time = time(nullptr);
	}
	long find_hits_time = (long)(find_hits_end_time - matcher_start_time );
	time_t matcher_end_time = time(nullptr);
	std::string const success_str( processor->match_processing_successful() ? "successful." : "not sucessful." );
	tr << "Matcher ran for " << (long)(matcher_end_time - matcher_start_time)
		<< " seconds, where finding hits took " << find_hits_time
		<< " seconds and processing the matches took " << processing_time << " seconds. Matching was "
		<< success_str << std::endl;

	if ( !processor->match_processing_successful() ) {
		set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
		pose = save_pose;
		return false;
	}
	this->set_last_move_status( protocols::moves::MS_SUCCESS );

	if ( incorporate_matches_into_pose_ ) {
		protocols::match::output::PoseMatchOutputWriterOP outputter(
			utility::pointer::static_pointer_cast< protocols::match::output::PoseMatchOutputWriter >( processor->output_writer() ) );
		core::pose::Pose const origpose = pose;
		if ( return_single_random_match_ ) {
			outputter->insert_match_into_pose( pose );
			return true;
		}
		core::Size const num_match_groups = outputter->match_groups_ushits().size();
		for ( core::Size mgroup=1; mgroup<=num_match_groups; ++mgroup ) {
			if ( mgroup == 1 ) {
				outputter->insert_match_into_pose( pose, mgroup );
				tr << "Incorporated match " << mgroup << " into pose." << std::endl;
			} else {
				core::pose::PoseOP matchedpose = origpose.clone();
				if ( origpose.pdb_info() ) {
					matchedpose->pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( *origpose.pdb_info() ) ) );
				}
				outputter->insert_match_into_pose( *matchedpose, mgroup );
				poselist.push_back( matchedpose );
				tr << "Incorporated match " << mgroup << " into a stored pose." << std::endl;
			}
		}
	}
	return true;
} //MatcherMover::process_pose function

// XRW TEMP std::string
// XRW TEMP MatcherMover::get_name() const
// XRW TEMP {
// XRW TEMP  return "MatcherMover";
// XRW TEMP }

void
MatcherMover::setup_seqpos_from_selectors( protocols::match::MatcherTask & mtask, core::pose::Pose const & pose ) const
{
	mtask.set_ignore_cmdline_for_build_points( true );
	mtask.use_different_build_points_for_each_geometric_constraint( selectors_.size() );

	core::Size cst_idx = 1;
	for ( auto s=selectors_.begin(); s!=selectors_.end(); ++s, ++cst_idx ) {
		core::select::residue_selector::ResidueSubset const subset = (*s)->apply( pose );
		utility::vector1< core::Size > residues;
		for ( core::Size resid=1; resid<=subset.size(); ++resid ) {
			if ( subset[ resid ] ) {
				residues.push_back( resid );
			}
		}
		tr << "Residues for geomcst " << cst_idx << ": " << residues << std::endl;
		mtask.set_original_scaffold_build_points_for_geometric_constraint( cst_idx, residues );
	}
}

void
MatcherMover::set_ligres(
	core::conformation::ResidueCOP ligres )
{
	ligres_ = core::conformation::ResidueCOP( core::conformation::ResidueOP( new core::conformation::Residue(*ligres) ) );
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
	basic::datacache::DataMap & data,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	incorporate_matches_into_pose_ = tag->getOption<bool>( "incorporate_matches_into_pose", 1 );
	if ( tag->hasOption( "residues_for_geomcsts" ) ) {
		std::string const selector_str = tag->getOption< std::string >( "residues_for_geomcsts" );
		utility::vector1< std::string > const selector_strs = utility::string_split( selector_str, ',' );
		for ( auto const & selector_str : selector_strs ) {
			core::select::residue_selector::ResidueSelectorCOP selector;
			try {
				selector = data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_str );
			} catch ( utility::excn::EXCN_Msg_Exception & e ) {
				std::stringstream error_msg;
				error_msg << "Failed to find ResidueSelector named '" << selector_str << "' from the Datamap from MatcherMover::parse_my_tag.\n";
				error_msg << e.msg();
				throw utility::excn::EXCN_RosettaScriptsOption( error_msg.str() );
			}
			debug_assert( selector );
			selectors_.push_back( selector );
		}
		tr << "Obtained residue selectors for " << selectors_.size() << " geomcsts." << std::endl;
	}

	// Exactly one of three things needs to be specified
	// -match:scaffold_active_site_residues
	// -match:scaffold_active_site_residues_for_geomcsts
	// residue selectors as XML "residues_for_geomcsts" option
	std::stringstream msg;
	msg << "MatcherMover: Exactly one of the following three things needs to be specified: " << std::endl;
	msg << "1) -match:scaffold_active_site_residues <filename> command line option" << std::endl;
	msg << "2) -match:scaffold_active_site_residues_for_geomcsts <filename> command line option" << std::endl;
	msg << "3) residues_for_geomcsts XML option to MatcherMover" << std::endl;

	bool bad_options = false;
	if ( basic::options::option[ basic::options::OptionKeys::match::scaffold_active_site_residues ].user() ) {
		if ( basic::options::option[ basic::options::OptionKeys::match::scaffold_active_site_residues_for_geomcsts ].user()
				|| selectors_.size() ) {
			bad_options = true;
			msg << "You have specified -match:scaffold_active_site_residues and ";
			if ( selectors_.empty() ) {
				msg << "-match:scaffold_active_site_residues_for_geomcsts";
			} else {
				msg << "residues_for_geomcsts (XML option)";
			}
		}
	} else if ( basic::options::option[ basic::options::OptionKeys::match::scaffold_active_site_residues_for_geomcsts ].user() ) {
		if ( selectors_.size() ) {
			bad_options = true;
			msg << "You have specified -match:scaffold_active_site_residues_for_geomcsts and residues_for_geomcsts (XML option)";
		}
	} else if ( selectors_.empty() ) {
		bad_options = true;
		msg << "You have specified none of these";
	}
	if ( bad_options ) {
		msg << "." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
}

std::string MatcherMover::get_name() const {
	return mover_name();
}

std::string MatcherMover::mover_name() {
	return "MatcherMover";
}

void MatcherMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	//residues_for_geomcsts: comma-separated list of residue selectors
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "incorporate_matches_into_pose", xsct_rosetta_bool, "Incorporate the identified matches into the input pose", "true" )
		+ XMLSchemaAttribute( "residues_for_geomcsts", xs_string, "A comma-separated list of residue selectors specifying residues to be used in geometric constraints" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Wrapper mover for the matcher. Constraints must be specified on the command line.", attlist );
}

std::string MatcherMoverCreator::keyname() const {
	return MatcherMover::mover_name();
}

protocols::moves::MoverOP
MatcherMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MatcherMover );
}

void MatcherMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MatcherMover::provide_xml_schema( xsd );
}


void
set_ligpose_rotamer( core::pose::Pose & ligpose )
{
	// Retrieve the rotamer library for this ligand;
	// check that the requested ligand-rotamer-index is in-bounds.
	// Relplace-residue on the ligpose with the rotamer from the library.

	using namespace core;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::pack::dunbrack;
	using namespace core::pack::rotamers;

	runtime_assert( ligpose.size() == 1 ); // we're expecting a one-residue pose.

	core::Size const lig_rotamer_index =
		basic::options::option[ basic::options::OptionKeys::match::ligand_rotamer_index ];

	if ( lig_rotamer_index < 1 ) {
		utility_exit_with_message( "Illegal rotamer index given in command line flag match::ligand_rotamer_index ("
			+ utility::to_string( lig_rotamer_index ) + ").  Must be greater than 0." );
	}

	SingleResidueRotamerLibraryFactory const & rotlib( *SingleResidueRotamerLibraryFactory::get_instance() );
	SingleResidueRotamerLibraryCOP res_rotlib( rotlib.get( ligpose.residue_type( 1 ) ) );

	if ( res_rotlib != nullptr ) {
		SingleLigandRotamerLibraryCOP lig_rotlib(
			utility::pointer::dynamic_pointer_cast< SingleLigandRotamerLibrary const > ( res_rotlib ));

		if ( lig_rotlib == nullptr ) {
			utility_exit_with_message( "Failed to retrieve a ligand rotamer library for "
				+ ligpose.residue_type(1).name() + " after finding the flag match::ligand_rotamer_index <int> on the command line");
		}

		core::pack::rotamers::RotamerVector rot_vector;
		lig_rotlib->build_base_rotamers( ligpose.residue_type( 1 ), rot_vector );
		Size const nligrots = rot_vector.size();

		if ( lig_rotamer_index > nligrots ) {
			utility_exit_with_message( "Illegal rotamer index given in command line flag match::ligand_rotamer_index ("
				+ utility::to_string( lig_rotamer_index ) + "). Index exceeds the number"
				" of ligand rotamers ( " + utility::to_string( nligrots ) + ")" );
		}

		ResidueCOP ligrot = rot_vector[ lig_rotamer_index ];
		ligpose.replace_residue( 1, *ligrot, false );
	} else {
		utility_exit_with_message( "Failed to find ligand rotamer library for " +
			ligpose.residue(1).name() + " after finding the flag -match::ligand_rotamer_index on the command line." );
	}
}

}
}
