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
/// @detailed responsibilities:
/// @author Oliver Lange

// Unit Headers
#include <protocols/topology_broker/TemplateFragmentClaimer.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/Exceptions.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameList.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <core/fragment/BBTorsionSRFD.hh>
// AUTO-REMOVED #include <protocols/jumping/JumpSetup.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>

// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>

#include <core/fragment/FragData.hh>
#include <utility/vector1.hh>


//// C++ headers

// option key includes

static thread_local basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace core;

TemplateFragmentClaimer::TemplateFragmentClaimer() : frag_size_( 9 ), config_file_( "NoFile" )
{}

TemplateFragmentClaimer::TemplateFragmentClaimer( std::string config_file, core::Size fragsize, weights::AbinitioMoverWeightOP weight )
	: FragmentClaimer( NULL, "template-frags", weight ),
		frag_size_( fragsize )
{
	read_config_file( config_file );
}

void TemplateFragmentClaimer::read_config_file( std::string const& file ) {
	templates_ = abinitio::TemplatesOP( new abinitio::Templates( file /*native_pose_*/ ) );
	//templates_->target_sequence() = sequence_; // a hack until class SequenceMapping works better
	// want to pick fragments from templates... make sure they are not initialized yet
	//	runtime_assert( !fragset_large_ );

	if( !templates_->is_good() ){
		utility_exit_with_message("ERROR occured during template setup. check BAD_SEQUENCES file!");
	}

	fragment::FragSetOP fragset_large( new fragment::ConstantLengthFragSet( frag_size_ ) );
	Size nr = templates_->pick_frags( *fragset_large, core::fragment::FragDataCOP( new fragment::FragData( core::fragment::SingleResidueFragDataOP( new fragment::BBTorsionSRFD ), frag_size_ ) ) );
	tr.Info << nr << " " << fragset_large->max_frag_length() << "mer fragments picked from homolog structures" << std::endl;
	set_mover( simple_moves::FragmentMoverOP( new simple_moves::ClassicFragmentMover( fragset_large ) ) );
}

bool TemplateFragmentClaimer::read_tag( std::string tag, std::istream& is ) {
	if ( tag == "FILE" ) {
		is >> config_file_;
	} else if ( tag == "FRAG_SIZE" ) {
		is >> frag_size_;
	} else if ( tag == "MOVER_WEIGHT" ) {
		read_mover_weight( is );
	} else if ( tag == "CMD_FLAG" ) {
		if ( basic::options::option[ basic::options::OptionKeys::templates::config ].user() ) {  //read from cmd-flag if not set explicitly
			config_file_ = basic::options::option[ basic::options::OptionKeys::templates::config ]();
		}
	} else return FragmentClaimer::read_tag( tag, is );
	return true;
}

void TemplateFragmentClaimer::init_after_reading() {
	if ( config_file_ != "NoFile" ) {
		read_config_file( config_file_ );
	}
	if ( !templates_ ) {
		throw EXCN_Input( "TemplateFragmentClaimer not initialized properly, use CMD_FLAG and -templates:config or FILE configfile" );
	}
}


} //topology_broker
} //protocols
