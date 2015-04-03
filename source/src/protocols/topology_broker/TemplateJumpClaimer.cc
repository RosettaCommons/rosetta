// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange

// Unit Headers
#include <protocols/topology_broker/TemplateJumpClaimer.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/Exceptions.hh>

#ifdef WIN32
#include <core/fragment/Frame.hh>
#endif

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <protocols/abinitio/TemplateJumpSetup.hh>
#include <protocols/abinitio/PairingStatistics.hh>
#include <core/scoring/dssp/PairingsList.hh>
#include <protocols/jumping/RandomSheetBuilder.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>

//// C++ headers
#include <fstream>

#include <utility/vector1.hh>
#include <boost/algorithm/string/erase.hpp>


// option key includes

static thread_local basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace core;

TemplateJumpClaimer::TemplateJumpClaimer() : templates_( /* NULL */ )
{
	set_keep_jumps_from_input_pose( false );
}

TemplateJumpClaimer::TemplateJumpClaimer( std::string config_file,  weights::AbinitioMoverWeightOP weight )
	: FragmentJumpClaimer( NULL, "template-jumps", weight )
{
	set_keep_jumps_from_input_pose( false );
	read_config_file( config_file );
}

void TemplateJumpClaimer::read_config_file( std::string const& file ) {
	templates_ = abinitio::TemplatesOP( new abinitio::Templates( file, NULL /* native_pose_*/ ) );
	//	templates_->target_sequence() = sequence_; // a hack until class SequenceMapping works better
	// want to pick fragments from templates... make sure they are not initialized yet
	//runtime_assert( !fragset_large_ );

	if( !templates_->is_good() ){
		utility_exit_with_message("ERRORS occured during template setup. check BAD_SEQUENCES file!");
	}
	set_jump_def( templates_->create_jump_def( NULL ) );
	/*NULL for ss_def_ ---- it will be determined from templates.. and is basically not used, since broker builds fold-tree himself;*/

}

void TemplateJumpClaimer::read_topol_file( std::string const& file ) {
		utility::io::izstream is( file );
		if ( !is.good() ) {
			utility_exit_with_message(" did not find topology_file: " + file );
		}
		abinitio::PairingStatisticsOP ps( new abinitio::PairingStatistics );
		is >> *ps;
		tr.Info << *ps << std::endl;
		core::scoring::dssp::PairingList helix_pairings; //empty for now
		core::fragment::SecondaryStructureOP ss_def( new core::fragment::SecondaryStructure );
		ss_def->extend( 1500 ); //HACK number of residues
		set_jump_def( jumping::BaseJumpSetupOP( new abinitio::TemplateJumpSetup( NULL, ss_def, ps, helix_pairings ) ) );
}

bool TemplateJumpClaimer::read_tag( std::string tag, std::istream& is ) {
	if ( tag == "file" ) {
			std::string file;
			is >> file;
			if ( jump_def() ) EXCN_Input( "Define either template config-file or topology_file for " + type() );
			read_config_file( file );
	} else if ( tag == "PAIRING_FILE" ) {
		// get pairings file
		std::string file;
		is >> file;
		pairings_.clear();
		read_pairing_list( file, pairings_ );
		tr.Info << "found " << pairings_.size() << " pairings in file " << file << std::endl;
	} else if ( tag == "SHEETS" || tag == "RANDOM_SHEETS" ) {
		//interpret rest of line as sheet topology
		std::string line;
		getline( is, line );
		std::istringstream is_line( line );
		copy( std::istream_iterator< core::Size >( is_line ), std::istream_iterator< core::Size >(), std::back_inserter( sheets_ ) );

		if ( tag == "RANDOM_SHEETS" ) bRandomSheet_ = true;
		else bRandomSheet_ = false;
	} else if ( tag == "SS_INFO" ) {
		std::string file;
		is >> file;
		ss_def_ = core::fragment::SecondaryStructureOP( new core::fragment::SecondaryStructure );
		ss_def_->read_psipred_ss2( file );
	} else if ( tag == "mover_weight" ) {
		read_mover_weight( is );
	} else if ( tag == "topol_file" || tag == "TOPOL_FILE" ) {
		std::string file;
		is >> file;
		if ( jump_def() ) EXCN_Input( "Define either template config-file or topology_file for " + type() );
		read_topol_file( file );
	} else return FragmentJumpClaimer::read_tag( tag, is );
	return true;
}

void TemplateJumpClaimer::init_after_reading() {
	if ( !jump_def() && basic::options::option[ basic::options::OptionKeys::templates::config ].user() && !templates_ ) {  //read from cmd-flag if not set explicitly
		read_config_file( basic::options::option[ basic::options::OptionKeys::templates::config ]() );
	}
	if ( pairings_.size() ) {
		if ( !ss_def_ ) EXCN_Input( "If you use PAIRING_FILE you also have to specify SS_INFO " );
		if ( !sheets_.size() ) EXCN_Input( "If you use PAIRING_FILE you also have to specify SHEET or RANDOM_SHEET " );
		if ( bRandomSheet_ ) set_jump_def( jumping::BaseJumpSetupOP( new jumping::RandomSheetBuilder( ss_def_, pairings_, sheets_ ) ) );
		else set_jump_def( jumping::BaseJumpSetupOP( new jumping::SheetBuilder( ss_def_, pairings_, sheets_ ) ) );
	}
}


} //topology_broker
} //protocols
