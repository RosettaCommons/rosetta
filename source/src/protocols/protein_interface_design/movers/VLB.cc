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
#include <core/pose/ResidueIndexDescription.hh>
#include <core/types.hh>
#include <string>
#include <basic/options/option.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh> // used for internal tracking of ntrials (fragment caching)

#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;
using namespace protocols::forge::build;
using namespace protocols::forge::components;

static basic::Tracer TR( "protocols.protein_interface_design.movers.VLB" );




VLB::VLB() :
	protocols::moves::Mover( VLB::mover_name() )
{
	scorefxn_ = utility::pointer::make_shared< core::scoring::ScoreFunction >();
} // default ctor//design_ = true;

void
VLB::apply( pose::Pose & pose ) {

	protocols::forge::build::BuildManagerOP manager = make_manager( pose );

	// internally handles ntrials, so that fragments can be cached!
	VarLengthBuild vlb( *manager );
	vlb.scorefunction( scorefxn_->clone() );
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


protocols::moves::MoverOP VLB::clone() const {
	return( utility::pointer::make_shared< VLB >( *this ));
}

protocols::moves::MoverOP VLB::fresh_instance() const {
	return utility::pointer::make_shared< VLB >();
}


void
VLB::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data
)
{
	using namespace protocols::rosetta_scripts;
	//RelativeConnectRight( rel_seq_pos, right, connect_pose ); /// version of ConnectRight instruction that depends upon results from another BuildInstruction

	std::string const scorefxn( tag->getOption<std::string>( "scorefxn", "score4L" ) ); // VarLengthBuild uses remodel_cen by default. score4L is the same, but with Rg=2.0
	scorefxn_ = data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn )->clone();

	utility::vector0< TagCOP > const & tags( tag->getTags() );
	for ( auto tag : tags ) {
		VLBInstruction instr;

		if ( tag->getName() == "Bridge" ) {
			// connect two contiguous but disjoint sections of a Pose into one continuous section
			instr.type = VLBInstructionType::Bridge;

			string const res1( tag->getOption< std::string >( "left" ) );
			instr.res1 = core::pose::parse_resnum( res1 );
			string const res2( tag->getOption< std::string >( "right" ) );
			instr.res2 = core::pose::parse_resnum( res2 );

			instr.ss = tag->getOption< std::string >( "ss", "" );
			instr.aa = tag->getOption< std::string >( "aa", "" );

		} else if ( tag->getName() == "ConnectRight" ) {
			//instruction to jump-connect one Pose onto the right side of another
			instr.type = VLBInstructionType::ConnectRight;

			string const res1( tag->getOption< std::string >( "left" ) );
			instr.res1 = core::pose::parse_resnum( res1 );
			string const res2( tag->getOption< std::string >( "right" ) );
			instr.res2 = core::pose::parse_resnum( res2 );

			instr.filename = tag->getOption< std::string >( "pdb", "" );
			runtime_assert( instr.filename != "" );

		} else if ( tag->getName() == "GrowLeft" ) {
			/// Use this for n-side insertions, but typically not n-terminal
			///  extensions unless necessary.  It does not automatically cover the
			///  additional residue on the right endpoint that needs to move during
			///  n-terminal extensions due to invalid phi torsion.  For that case,
			///  use the SegmentRebuild class replacing the n-terminal residue with
			///  desired length+1.
			instr.type = VLBInstructionType::GrowLeft;

			string const res1( tag->getOption< std::string >( "pos" ) );
			instr.res1 = core::pose::parse_resnum( res1 );

			instr.ss = tag->getOption< std::string >( "ss", "" );
			instr.aa = tag->getOption< std::string >( "aa", "" );

		} else if ( tag->getName() == "GrowRight" ) {
			/// instruction to create a c-side extension
			instr.type = VLBInstructionType::GrowRight;

			string const res1( tag->getOption< std::string >( "pos" ) );
			instr.res1 = core::pose::parse_resnum( res1 );

			instr.ss = tag->getOption< std::string >( "ss", "" );
			instr.aa = tag->getOption< std::string >( "aa", "" );

		} if ( tag->getName() == "SegmentInsert" ) {
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
			instr.type = VLBInstructionType::SegmentInsert;

			string const res1( tag->getOption< std::string >( "left" ) );
			instr.res1 = core::pose::parse_resnum( res1 );
			string const res2( tag->getOption< std::string >( "right" ) );
			instr.res2 = core::pose::parse_resnum( res2 );

			instr.ss = tag->getOption< std::string >( "ss", "L^L" );
			//instr.aa = tag->getOption< std::string >( "aa", "" );
			instr.keep_bb = tag->getOption< bool >( "keep_bb_torsions", false );
			instr.filename = tag->getOption< std::string >( "pdb", "" );
			runtime_assert( instr.filename != "" );

			instr.side = tag->getOption< std::string >( "side", "" );

		} else if ( tag->getName() == "SegmentRebuild" ) {
			/// @brief instruction to rebuild a segment
			instr.type = VLBInstructionType::SegmentRebuild;

			string const res1( tag->getOption< std::string >( "left" ) );
			instr.res1 = core::pose::parse_resnum( res1 );
			string const res2( tag->getOption< std::string >( "right" ) );
			instr.res2 = core::pose::parse_resnum( res2 );

			instr.ss = tag->getOption< std::string >( "ss", "" );
			instr.aa = tag->getOption< std::string >( "aa", "" );

		} else if ( tag->getName() == "SegmentSwap" ) {
			/// instruction to swap a segment with an external segment
			/// interval: swap out this range of residues
			/// move_map: fixed backbone residues in this movemap will be used for new jumps
			/// swap_in: swap in this pose
			instr.type = VLBInstructionType::SegmentSwap;

			string const res1( tag->getOption< std::string >( "left" ) );
			instr.res1 = core::pose::parse_resnum( res1 );
			string const res2( tag->getOption< std::string >( "right" ) );
			instr.res2 = core::pose::parse_resnum( res2 );

			instr.filename = tag->getOption< std::string >( "pdb", "" );
			runtime_assert( instr.filename != "" );
		}

		instruction_list_.push_back( instr );
	}

	TR<<"defined VLB mover with " << instruction_list_.size() << " instructions." << std::endl;
}

protocols::forge::build::BuildManagerOP
VLB::make_manager( core::pose::Pose const & pose ) const {
	auto manager = utility::pointer::make_shared< BuildManager >();
	for ( VLBInstruction const & instr: instruction_list_ ) {
		BuildInstructionOP instruction;
		switch ( instr.type ) {
		case VLBInstructionType::Bridge :
			{
			Interval const ival( instr.res1->resolve_index(pose), instr.res2->resolve_index(pose) );
			instruction = utility::pointer::make_shared< Bridge >( ival, instr.ss, instr.aa );
			break;
		}
		case VLBInstructionType::ConnectRight :
			{
			pose::Pose connect_pose;
			core::import_pose::pose_from_file( connect_pose, instr.filename , core::import_pose::PDB_file);
			instruction = utility::pointer::make_shared< ConnectRight >( instr.res1->resolve_index(pose), instr.res2->resolve_index(pose), connect_pose );
			break;
		}
		case VLBInstructionType::GrowLeft :
			{
			instruction = utility::pointer::make_shared< GrowLeft >( instr.res1->resolve_index(pose), instr.ss, instr.aa );
			break;
		}
		case VLBInstructionType::GrowRight :
			{
			instruction = utility::pointer::make_shared< GrowRight >( instr.res1->resolve_index(pose), instr.ss, instr.aa );
			break;
		}
		case VLBInstructionType::SegmentInsert :
			{
			pose::Pose insert_pose;
			core::import_pose::pose_from_file( insert_pose, instr.filename, core::import_pose::PDB_file);

			Interval const ival( instr.res1->resolve_index(pose), instr.res2->resolve_index(pose) );
			SegmentInsertConnectionScheme::Enum connect_side;
			if ( instr.side == "N" ) connect_side = SegmentInsertConnectionScheme::N;
			else if ( instr.side == "C" ) connect_side = SegmentInsertConnectionScheme::C;
			else connect_side = SegmentInsertConnectionScheme::RANDOM_SIDE;

			instruction = utility::pointer::make_shared< SegmentInsert >( ival, instr.ss, insert_pose, instr.keep_bb, connect_side );
			break;
		}
		case VLBInstructionType::SegmentRebuild :
			{
			Interval const ival( instr.res1->resolve_index(pose), instr.res2->resolve_index(pose) );
			instruction = utility::pointer::make_shared< SegmentRebuild >( ival, instr.ss, instr.aa );
			break;
		}
		case VLBInstructionType::SegmentSwap :
			{
			Interval const ival( instr.res1->resolve_index(pose), instr.res2->resolve_index(pose) );

			pose::Pose swap_pose;
			core::import_pose::pose_from_file( swap_pose, instr.filename, core::import_pose::PDB_file);
			core::kinematics::MoveMap movemap; // empty = place jump anywhere
			instruction = utility::pointer::make_shared< SegmentSwap >( ival, movemap, swap_pose );
			break;
		}
		default :
			utility_exit_with_message("Unknown instruction type in VLB");
		}
		manager->add( instruction );
	}

	runtime_assert( manager->compatibility_check() );

	return manager;
}

std::string VLB::get_name() const {
	return mover_name();
}

std::string VLB::mover_name() {
	return "VLB";
}

std::string mangled_name_for_VLB( std::string const & foo ) {
	return "VLB_subtag_" + foo + "_type";
}

void VLB::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "scorefxn", xs_string, "Scorefunction to be used", "score4L" );

	AttributeList bridge_subtag_attlist;
	bridge_subtag_attlist + XMLSchemaAttribute::required_attribute( "left", xsct_refpose_enabled_residue_number, "Left residue" )
		+ XMLSchemaAttribute::required_attribute( "right", xsct_refpose_enabled_residue_number, "Right residue" )
		+ XMLSchemaAttribute::required_attribute( "ss", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "aa", xs_string, "XRW TO DO" );

	AttributeList connectright_subtag_attlist;
	connectright_subtag_attlist + XMLSchemaAttribute::required_attribute( "left", xsct_refpose_enabled_residue_number, "Left residue" )
		+ XMLSchemaAttribute::required_attribute( "right", xsct_refpose_enabled_residue_number, "Right residue" )
		+ XMLSchemaAttribute::required_attribute( "pdb", xs_string, "PDB file name to be read in" );

	AttributeList grow_subtag_attlist;
	grow_subtag_attlist + XMLSchemaAttribute::required_attribute( "pos", xsct_refpose_enabled_residue_number, "The single residue from which to build" )
		+ XMLSchemaAttribute::required_attribute( "ss", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "aa", xs_string, "XRW TO DO" );

	AttributeList segmentinsert_subtag_attlist;
	segmentinsert_subtag_attlist + XMLSchemaAttribute::required_attribute( "left", xsct_refpose_enabled_residue_number, "Left residue" )
		+ XMLSchemaAttribute::required_attribute( "right", xsct_refpose_enabled_residue_number, "Right residue" )
		+ XMLSchemaAttribute::required_attribute( "ss", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "keep_bb_torsions", xsct_rosetta_bool, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::required_attribute( "pdb", xs_string, "PDB file name to be read in" )
		+ XMLSchemaAttribute::required_attribute( "side", xs_string, "XRW TO DO" );

	AttributeList segmentrebuild_subtag_attlist;
	segmentrebuild_subtag_attlist+ XMLSchemaAttribute::required_attribute( "left", xsct_refpose_enabled_residue_number, "Left residue" )
		+ XMLSchemaAttribute::required_attribute( "right", xsct_refpose_enabled_residue_number, "Right residue" )
		+ XMLSchemaAttribute::required_attribute( "ss", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "aa", xs_string, "XRW TO DO" );

	AttributeList segmentswap_subtag_attlist;
	segmentswap_subtag_attlist+ XMLSchemaAttribute::required_attribute( "left", xsct_refpose_enabled_residue_number, "Left residue" )
		+ XMLSchemaAttribute::required_attribute( "right", xsct_refpose_enabled_residue_number, "Right residue" )
		+ XMLSchemaAttribute::required_attribute( "pdb", xs_string, "PDB file name to be read in" );

	utility::tag::XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement( "Bridge", bridge_subtag_attlist, "XRW TODO"/*, 0 minoccurs*/ )
		.add_simple_subelement( "ConnectRight", connectright_subtag_attlist, "XRW TODO"/*, 0 minoccurs*/ )
		.add_simple_subelement( "GrowLeft", grow_subtag_attlist, "XRW TODO"/*, 0 minoccurs*/ )
		.add_simple_subelement( "GrowRight", grow_subtag_attlist, "XRW TODO"/*, 0 minoccurs*/ )
		.add_simple_subelement( "SegmentInsert", segmentinsert_subtag_attlist, "XRW TODO"/*, 0 minoccurs*/ )
		.add_simple_subelement( "SegmentRebuild", segmentrebuild_subtag_attlist, "XRW TODO"/*, 0 minoccurs*/ )
		.add_simple_subelement( "SegmentSwap", segmentswap_subtag_attlist, "XRW TODO"/*, 0 minoccurs*/ )
		.complex_type_naming_func( & mangled_name_for_VLB );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, ssl );
}

std::string VLBCreator::keyname() const {
	return VLB::mover_name();
}

protocols::moves::MoverOP
VLBCreator::create_mover() const {
	return utility::pointer::make_shared< VLB >();
}

void VLBCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	VLB::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
