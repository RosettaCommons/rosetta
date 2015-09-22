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
#include <protocols/topology_broker/ConstraintClaimer.hh>

// Package Headers
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>

//#include <protocols/topology_broker/DofClaim.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loops.tmpl.hh>
#include <protocols/loops/LoopsFileIO.hh>
// #include <core/kinematics/MoveMap.hh>
// #include <core/fragment/FragSet.hh>
// #include <protocols/simple_moves/FragmentMover.hh>


#include <core/id/Exceptions.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/excn/Exceptions.hh>
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
//// C++ headers

#include <core/id/SequenceMapping.hh>
#include <utility/vector1.hh>


// option key includes

static THREAD_LOCAL basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {


using namespace core;
using namespace scoring::constraints;
ConstraintClaimer::ConstraintClaimer() :
	filename_( "NO_FILE"),
	tag_( "NO_TAG" ),
	constraints_( /* NULL */ ),
	bCentroid_( true ),
	bFullatom_( false ),
	bCmdFlag_( false ),
	combine_ratio_( 1 ),
	drop_random_rate_( 0.0 ),
	skip_redundant_( false ),
	skip_redundant_width_( 1 ),
	filter_weight_( 0.0 ),
	filter_name_( "" )
{}

ConstraintClaimer::ConstraintClaimer( std::string filename, std::string tag ) :
	filename_( filename ),
	tag_( tag ),
	constraints_( /* NULL */ ),
	bCentroid_( true ),
	bFullatom_( false ),
	bCmdFlag_( false ),
	combine_ratio_( 1 ),
	drop_random_rate_( 0 ),
	skip_redundant_( false ),
	skip_redundant_width_( 1 ),
	filter_weight_( 0.0 ),
	filter_name_( "" )
{}

ConstraintClaimer::ConstraintClaimer( bool CmdFlag, bool centroid, bool fullatom )
: filename_( "" ),
	tag_( "" ),
	constraints_( /* NULL */ ),
	bCentroid_( centroid ),
	bFullatom_( fullatom ),
	bCmdFlag_( CmdFlag ),
	combine_ratio_( 1 ),
	drop_random_rate_( 0 ),
	skip_redundant_( false ),
	skip_redundant_width_( 1 ),
	filter_weight_( 0.0 ),
	filter_name_( "" )
{
	runtime_assert( CmdFlag );
}

void ConstraintClaimer::generate_claims( claims::DofClaims& /*new_claims*/ ) {
}

void ConstraintClaimer::new_decoy() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	tr.Debug << "ConstraintClaimer::new_decoy: cst-modus: " << ( bFullatom_ ? " fullatom " : "no fullatom" )
		<<  ( bCentroid_ ? " centroid " : " no centroid " )
		<< std::endl;
	if ( bCmdFlag_ && option[ constraints::combine ].user() ) {
		combine_ratio_ = option[ constraints::combine ]();
	}

	if ( bCmdFlag_ && ( option[ constraints::skip_redundant ]() || option[ constraints::skip_redundant_width ].user() ) ) {
		skip_redundant_ = true;
		skip_redundant_width_ = option[ constraints::skip_redundant_width ];
	}


	if ( bCmdFlag_ && option[ constraints::cst_file ].user() ) {
		// reads and sets constraints -- this might be different each time we call this function
		std::string old_filename = filename_;
		filename_ = core::scoring::constraints::get_cst_file_option();
		if ( old_filename != filename_ ) constraints_ = NULL;
	}
	if ( bCmdFlag_ && option[ constraints::cst_fa_file ].user() && bFullatom_ && !bCentroid_ ) {
		std::string old_filename = filename_;
		filename_ = core::scoring::constraints::get_cst_fa_file_option();
		if ( old_filename != filename_ ) constraints_ = NULL;
	}
}


void ConstraintClaimer::add_constraints( core::pose::Pose& pose ) const {
	using namespace basic::options;
	bool fullatom( pose.is_fullatom() );
	if ( fullatom && !bFullatom_ ) return;
	if ( !fullatom && !bCentroid_ ) return;
	tr.Debug << "add constraints "<< tag_ << std::endl;
	if ( constraints_ ) {
		tr.Debug << " constraint set is currently for a " <<( constraint_ref_pose_.is_fullatom() ? "fullatom" : "centroid") << " pose "
			<< "\n will now remap them to a " <<  (fullatom ? "fullatom" : "centroid") << " pose" << std::endl;
	}
	std::string const new_sequence ( pose.annotated_sequence( true ) );
	if ( bCmdFlag_ && option[ OptionKeys::constraints::combine_exclude_region ].user() && combine_exclude_res_.size() == 0 && sequence_ != new_sequence ) {
		std::string const file( option[ OptionKeys::constraints::combine_exclude_region ]() );
		std::ifstream is( file.c_str() );

		if ( !is.good() ) {
			utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + file + "'" );
		}

		loops::PoseNumberedLoopFileReader reader;
		reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
		loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file(is, file, false );
		loops::Loops rigid_core = loops::Loops( loops );

		combine_exclude_res_.resize( pose.total_residue(), false );
		rigid_core.transfer_to_residue_vector( combine_exclude_res_, true );
	}
	if ( !constraints_ || sequence_ != new_sequence ) {
		tr.Info << " read constraints from " << filename_ << "\n for pose " << new_sequence << "..." << std::endl;
		constraints_ = ConstraintIO::get_instance()->read_constraints( filename_, ConstraintSetOP( new ConstraintSet ), pose );
		sequence_ = new_sequence;
	} else {
		ConstraintSetOP new_cst(NULL);
		try {
			new_cst = constraints_->remapped_clone( constraint_ref_pose_, pose );
		} catch( core::id::EXCN_AtomNotFound& excn ) {
			tr.Error << "[ERROR] failed attempt to add constraints to the "
				<< (fullatom ? "fullatom" : "centroid") << " pose" << std::endl;
			tr.Error << excn << std::endl;
			if ( tr.Debug.visible() ) {
				pose.dump_pdb("new_pose_failed_constraints.pdb");
				constraint_ref_pose_.dump_pdb("cst_ref_pose.pdb");
			}
			constraints_->show_definition( tr.Error, constraint_ref_pose_ );
			tr.Error << std::endl;
			tr.Error << " try to recover by reading in original constraints from " << filename_ << "\n for pose " << new_sequence << "..." << std::endl;
			new_cst = ConstraintIO::get_instance()->read_constraints( filename_, ConstraintSetOP( new ConstraintSet ), pose );
		}

		constraints_ = new_cst;
	}
	constraint_ref_pose_ = pose;

	scoring::constraints::ConstraintCOPs added_constraints = constraints_->get_all_constraints();
	if ( skip_redundant_ ) scoring::constraints::skip_redundant_constraints( added_constraints, pose.total_residue(), skip_redundant_width_ );
	if ( drop_random_rate_ > 0.0 ) scoring::constraints::drop_constraints( added_constraints, drop_random_rate_ );

	kinematics::ShortestPathInFoldTree sp( pose.fold_tree() );
	scoring::constraints::choose_effective_sequence_separation( sp, added_constraints );
	scoring::constraints::combine_constraints( added_constraints, combine_ratio_, combine_exclude_res_, sp ); // if combine_ratio_ > 1 this will randomly combine constraints into multipletts with OR logic
	pose.add_constraints( added_constraints );
	if ( tr.Trace.visible() ) {
		pose.constraint_set()->show_definition( tr.Trace, pose );
	}

}

bool ConstraintClaimer::read_tag( std::string tag, std::istream& is ) {
	if ( tag == "file" || tag == "FILE" || tag == "CST_FILE" ) {
		is >> filename_;
	} else if ( tag == "NO_CENTROID" ) {
		bCentroid_ = false;
	} else if ( tag == "CENTROID" ) {
		bCentroid_ = true;
	} else if ( tag == "FULLATOM" ) {
		bFullatom_ = true;
	} else if ( tag == "SKIP_REDUNDANT" ) {
		skip_redundant_ = true;
		is >> skip_redundant_width_;
		if ( !is.good() ) skip_redundant_width_ = 1;
	} else if ( tag == "CMD_FLAG" ) {
		bCmdFlag_ = true;
	} else if ( tag == "COMBINE_RATIO" ) {
		is >> combine_ratio_;
	} else if ( tag == "DROP_RANDOM_RATE" ) {
		is >> drop_random_rate_;
	} else if ( tag == "FILTER_WEIGHT" ) {
		is >> filter_weight_;
	} else if ( tag == "FILTER_NAME" ) {
		is >> filter_name_;
	} else return Parent::read_tag( tag, is );
	return true;
}

void ConstraintClaimer::set_cst_file( std::string const& file ) {
	filename_ = file;
	constraints_ = NULL;
	fa_constraints_ = NULL;
}

void ConstraintClaimer::set_fullatom( bool setting ) {
	bFullatom_ = setting;
	constraints_ = NULL;
	fa_constraints_ = NULL;
}

void ConstraintClaimer::set_centroid( bool setting ) {
	bCentroid_ = setting;
	constraints_ = NULL;
	fa_constraints_ = NULL;
}

void ConstraintClaimer::set_skip_redundant( core::Size setting ) {
	skip_redundant_ = setting > 0;
	skip_redundant_width_=setting;
	constraints_ = NULL;
	fa_constraints_ = NULL;
}

void ConstraintClaimer::set_combine_ratio( core::Size setting ) {
	combine_ratio_ = setting;
	constraints_ = NULL;
	fa_constraints_ = NULL;
}

void ConstraintClaimer::set_filter_weight( core::Real setting ) {
	filter_weight_ = setting;
}


} //topology_broker
} //protocols
