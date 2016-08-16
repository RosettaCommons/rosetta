// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/VLB.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/VLB.hh>
#include <protocols/protein_interface_design/movers/VLBCreator.hh>

// Project headers
#include <protocols/moves/Mover.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/BuildInstruction.hh>
#include <protocols/forge/build/Bridge.hh>
#include <protocols/forge/build/ConnectRight.hh>
#include <protocols/forge/build/GrowLeft.hh>
#include <protocols/forge/build/GrowRight.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/build/SegmentRebuild.hh>
#include <protocols/forge/build/SegmentSwap.hh>

#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.hh>

#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <string>
#include <basic/options/option.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh> // used for internal tracking of ntrials (fragment caching)

#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;
using namespace protocols::forge::build;
using namespace protocols::forge::components;

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.VLB" );

std::string
VLBCreator::keyname() const
{
	return VLBCreator::mover_name();
}

protocols::moves::MoverOP
VLBCreator::create_mover() const {
	return protocols::moves::MoverOP( new VLB );
}

std::string
VLBCreator::mover_name()
{
	return "VLB";
}

VLB::VLB() :
	protocols::moves::Mover( VLBCreator::mover_name() )
{
	manager_ = protocols::forge::build::BuildManagerOP( new BuildManager );
	scorefxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
} // default ctor//design_ = true;

VLB::VLB( BuildManagerCOP manager, ScoreFunctionCOP scorefxn ) :
	protocols::moves::Mover( VLBCreator::mover_name() )
{
	manager_ = protocols::forge::build::BuildManagerOP( new BuildManager( *manager ) );
	scorefxn_ = scorefxn->clone();
}

void
VLB::apply( pose::Pose & pose ) {

	// internally handles ntrials, so that fragments can be cached!
	VarLengthBuild vlb( *manager_ );
	vlb.scorefunction( scorefxn_ );
	core::Size ntrials(1);
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( option[ jd2::ntrials ].user() ) ntrials = option[ jd2::ntrials ]() ;
	while ( ntrials > 0 ) {
		vlb.apply( pose );
		if ( vlb.get_last_move_status() == MS_SUCCESS ) break;
		else ntrials--;
	}

	if ( vlb.get_last_move_status() == FAIL_DO_NOT_RETRY ) set_last_move_status( FAIL_RETRY );
	else if ( vlb.get_last_move_status() == MS_SUCCESS ) set_last_move_status( vlb.get_last_move_status() );
	else {
		TR << "VLB mover status was not properly set! Aborting.";
		set_last_move_status(FAIL_DO_NOT_RETRY); // we should never get to this point!
	}

}

std::string
VLB::get_name() const {
	return VLBCreator::mover_name();
}

protocols::moves::MoverOP VLB::clone() const {
	return( protocols::moves::MoverOP( new VLB( *this ) ));
}

protocols::moves::MoverOP VLB::fresh_instance() const {
	return protocols::moves::MoverOP( new VLB );
}

VLB::VLB( VLB const & init ) :
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( init ) {
	manager_ = protocols::forge::build::BuildManagerOP( new BuildManager( *init.manager_ ) );
	scorefxn_ = init.scorefxn_->clone();
}

VLB & VLB::operator=( VLB const & init ) {
	manager_ = protocols::forge::build::BuildManagerOP( new BuildManager( *init.manager_ ) );
	scorefxn_ = init.scorefxn_->clone();
	return *this;
}


VLB::~VLB() {}

void
VLB::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	Movers_map const &,
	core::pose::Pose const & pose
)
{
	using namespace protocols::rosetta_scripts;
	//RelativeConnectRight( rel_seq_pos, right, connect_pose ); /// version of ConnectRight instruction that depends upon results from another BuildInstruction

	std::string const scorefxn( tag->getOption<std::string>( "scorefxn", "score4L" ) ); // VarLengthBuild uses remodel_cen by default. score4L is the same, but with Rg=2.0
	scorefxn_ = data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn )->clone();

	BuildInstructionOP instruction;
	utility::vector0< TagCOP > const & tags( tag->getTags() );
	for ( utility::vector0< TagCOP >::const_iterator it=tags.begin(); it!=tags.end(); ++it ) {
		TagCOP const tag = *it;
		if ( tag->getName() == "Bridge" ) {
			// connect two contiguous but disjoint sections of a Pose into one continuous section
			string const res1( tag->getOption< std::string >( "left" ) );
			Size const left = core::pose::parse_resnum( res1, pose );
			string const res2( tag->getOption< std::string >( "right" ) );
			Size const right = core::pose::parse_resnum( res2, pose );

			string const ss( tag->getOption< std::string >( "ss", "" ) );
			string const aa( tag->getOption< std::string >( "aa", "" ) );

			Interval const ival( left, right);
			instruction = BuildInstructionOP( new Bridge( ival, ss, aa ) );
		}
		if ( tag->getName() == "ConnectRight" ) {
			//instruction to jump-connect one Pose onto the right side of another
			string const res1( tag->getOption< std::string >( "left" ) );
			Size const left = core::pose::parse_resnum( res1, pose );
			string const res2( tag->getOption< std::string >( "right" ) );
			Size const right = core::pose::parse_resnum( res2, pose );

			string const pose_fname( tag->getOption< std::string >( "pdb", "" ) );
			runtime_assert( pose_fname != "" );
			pose::Pose connect_pose;
			core::import_pose::pose_from_file( connect_pose, pose_fname , core::import_pose::PDB_file);
			instruction = BuildInstructionOP( new ConnectRight( left, right, connect_pose ) );
		}
		if ( tag->getName() == "GrowLeft" ) {
			/// Use this for n-side insertions, but typically not n-terminal
			///  extensions unless necessary.  It does not automatically cover the
			///  additional residue on the right endpoint that needs to move during
			///  n-terminal extensions due to invalid phi torsion.  For that case,
			///  use the SegmentRebuild class replacing the n-terminal residue with
			///  desired length+1.

			string const res1( tag->getOption< std::string >( "pos" ) );
			Size const pos = core::pose::parse_resnum( res1, pose );

			string const ss( tag->getOption< std::string >( "ss", "" ) );
			string const aa( tag->getOption< std::string >( "aa", "" ) );

			instruction = BuildInstructionOP( new GrowLeft( pos, ss, aa ) );
		}
		if ( tag->getName() == "GrowRight" ) {
			/// instruction to create a c-side extension
			string const res1( tag->getOption< std::string >( "pos" ) );
			Size const pos = core::pose::parse_resnum( res1, pose );

			string const ss( tag->getOption< std::string >( "ss", "" ) );
			string const aa( tag->getOption< std::string >( "aa", "" ) );

			instruction = BuildInstructionOP( new GrowRight( pos, ss, aa ) );
		}


		if ( tag->getName() == "SegmentInsert" ) {
			/// interval: The interval between which the insert will span.
			///  To perform a pure insertion without replacing any residues
			///  within a region, use an interval with a zero as the left endpoint, e.g.
			///  [0, insert_after_this_residue].  If inserting before the first residue
			///  of the Pose then interval = [0,0].  If inserting after the last residue
			///  of the Pose then interval = [0, last_residue].
			/// ss: The secondary structure specifying the flanking regions,
			///  with a character '^' specifying where the insert is to be placed.
			/// insert: The Pose to insert.
			/// keep_known_bb_torsions_at_junctions: Attempt to keep the omega
			///  at original_interval().left-1, the phi at original_interval().left, and
			///  the psi+omega at original_interval().right present from the original Pose
			///  in the modified Pose.  This should be false for pure insertions.
			/// connection_scheme: Connect insertion on its N-side, C-side,
			///  or decide randomly between the two (default RANDOM). Random is only random on parsing, not per ntrial.

			string const res1( tag->getOption< std::string >( "left" ) );
			Size const left = core::pose::parse_resnum( res1, pose );
			string const res2( tag->getOption< std::string >( "right" ) );
			Size const right = core::pose::parse_resnum( res2, pose );
			Interval const ival( left, right );

			string const ss( tag->getOption< std::string >( "ss", "L^L" ) );
			//string const aa( tag->getOption< std::string >( "aa", "" ) );
			bool const keep_bb_torsions( tag->getOption< bool >( "keep_bb_torsions", "0" ));
			string const pose_fname( tag->getOption< std::string >( "pdb", "" ) );

			runtime_assert( pose_fname != "" );
			pose::Pose insert_pose;
			core::import_pose::pose_from_file( insert_pose, pose_fname , core::import_pose::PDB_file);

			string const side( tag->getOption< std::string >( "side", "" ) );
			SegmentInsertConnectionScheme::Enum connect_side;
			if ( side == "N" ) connect_side = SegmentInsertConnectionScheme::N;
			else if ( side == "C" ) connect_side = SegmentInsertConnectionScheme::C;
			else connect_side = SegmentInsertConnectionScheme::RANDOM_SIDE;

			instruction = BuildInstructionOP( new SegmentInsert( ival, ss, insert_pose, keep_bb_torsions, connect_side ) );
		}


		if ( tag->getName() == "SegmentRebuild" ) {
			/// @brief instruction to rebuild a segment
			string const res1( tag->getOption< std::string >( "left" ) );
			Size const left = core::pose::parse_resnum( res1, pose );
			string const res2( tag->getOption< std::string >( "right" ) );
			Size const right = core::pose::parse_resnum( res2, pose );
			Interval const ival( left, right);

			string const ss( tag->getOption< std::string >( "ss", "" ) );
			string const aa( tag->getOption< std::string >( "aa", "" ) );

			instruction = BuildInstructionOP( new SegmentRebuild( ival, ss, aa ) );
		}

		if ( tag->getName() == "SegmentSwap" ) {
			/// instruction to swap a segment with an external segment
			/// interval: swap out this range of residues
			/// move_map: fixed backbone residues in this movemap will be used for new jumps
			/// swap_in: swap in this pose
			string const res1( tag->getOption< std::string >( "left" ) );
			Size const left = core::pose::parse_resnum( res1, pose );
			string const res2( tag->getOption< std::string >( "right" ) );
			Size const right = core::pose::parse_resnum( res2, pose );

			Interval const ival( left, right);

			string const pose_fname( tag->getOption< std::string >( "pdb", "" ) );
			runtime_assert( pose_fname != "" );
			pose::Pose swap_pose;
			core::import_pose::pose_from_file( swap_pose, pose_fname , core::import_pose::PDB_file);

			core::kinematics::MoveMap movemap; // empty = place jump anywhere
			SegmentSwap( ival, movemap, swap_pose );
		}

		manager_->add( instruction );
	}

	runtime_assert( manager_->compatibility_check() );
	TR<<"defined VLB mover with " << manager_->size() << " instructions." << std::endl;
}

} //movers
} //protein_interface_design
} //protocols
