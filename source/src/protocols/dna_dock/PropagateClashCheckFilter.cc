// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/dna_dock/PropagateClashCheckFilter.cc
/// @brief  Checks for fa_rep clashes between docked pose and DNA backbone.
/// @author Carl Walkey (cwalkey@uw.edu)

// Unit Headers
#include <protocols/dna_dock/PropagateClashCheckFilter.hh>
#include <protocols/dna_dock/PropagateClashCheckFilterCreator.hh>

// Basic Headers
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/out.OptionKeys.gen.hh>

// Core Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
//#include <core/pose/PDBInfo.hh>
#include <core/pose/subpose_manipulation_util.hh>

#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
//#include <core/pose/symmetry/util.hh>
//#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
//#include <core/chemical/ResidueConnection.hh>
//#include <core/conformation/Conformation.hh>
//#include <core/conformation/symmetry/SymmetryInfo.hh>
//#include <core/scoring/sasa.hh>
//#include <core/chemical/ChemicalManager.hh>
//#include <core/id/AtomID_Map.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>



// Utility headers
//#include <utility/vector1.fwd.hh>

// Protocols Headers
#include <protocols/rosetta_scripts/util.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
//#include <core/kinematics/FoldTree.hh>
//#include <core/conformation/Residue.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
//#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

//#include <utility/vector0.hh>
//#include <utility/vector1.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

//// C++ headers
static basic::Tracer TR( "protocols.matdes.ClashCheckFilter" );

namespace protocols {
namespace dna_dock {

// @brief default constructor
PropagateClashCheckFilter::PropagateClashCheckFilter():
	fa_rep_thresh_( 10.0 ),
	default_bridge_type_( "VAL" ),
	omega_( 36.0 ),
	rise_( 3.38 ),
	num_repeats_( 22 ),
	prop_dir_( 0 )
{}

// @brief constructor with arguments
PropagateClashCheckFilter::PropagateClashCheckFilter(
	core::Real const fa_rep_thresh,
	std::string default_bridge_type,
	core::Real omega,
	core::Real rise,
	core::Size num_repeats,
	core::Size prop_dir ):
	fa_rep_thresh_( fa_rep_thresh ),
	default_bridge_type_( default_bridge_type ),
	omega_( omega ),
	rise_( rise ),
	num_repeats_( num_repeats ),
	prop_dir_( prop_dir )
{}

// @brief copy constructor
PropagateClashCheckFilter::PropagateClashCheckFilter( PropagateClashCheckFilter const & )= default;

// @brief destructor
PropagateClashCheckFilter::~PropagateClashCheckFilter() = default;

protocols::filters::FilterOP
PropagateClashCheckFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new PropagateClashCheckFilter() );
}

protocols::filters::FilterOP
PropagateClashCheckFilter::clone() const{
	return protocols::filters::FilterOP( new PropagateClashCheckFilter( *this ) );
}

// @brief getters
core::Real PropagateClashCheckFilter::get_fa_rep_thresh() const { return fa_rep_thresh_; }
std::string PropagateClashCheckFilter::get_default_bridge_type() const { return default_bridge_type_; }
core::Real PropagateClashCheckFilter::get_omega() const { return omega_; }
core::Real PropagateClashCheckFilter::get_rise() const { return rise_; }
core::Size PropagateClashCheckFilter::get_num_repeats() const { return num_repeats_; }

// @brief setters
void PropagateClashCheckFilter::set_fa_rep_thresh( core::Real fa_rep_thresh ) { fa_rep_thresh_ = fa_rep_thresh; }
void PropagateClashCheckFilter::set_default_bridge_type( std::string default_bridge_type ) { default_bridge_type_ = default_bridge_type; }
void PropagateClashCheckFilter::set_omega( core::Real omega ) { omega_ = omega; }
void PropagateClashCheckFilter::set_rise( core::Real rise ) { rise_ = rise; }
void PropagateClashCheckFilter::set_num_repeats( core::Size num_repeats ) { num_repeats_ = num_repeats; }

// @brief returns true if the set of residues defined by the TaskOperations have a no clashes. False otherwise.
bool PropagateClashCheckFilter::apply( Pose const & pose ) const
{

	TR << "fa_rep_thresh_: " << std::to_string(fa_rep_thresh_) << std::endl;
	TR << "default_bridge_type_: " << default_bridge_type_ << std::endl;
	TR << "omega_: " << std::to_string(omega_) << std::endl;
	TR << "rise_: " << std::to_string(rise_) << std::endl;
	TR << "num_repeats_: " << std::to_string(num_repeats_) << std::endl;
	TR << "prop_dir_: " << std::to_string(prop_dir_) << std::endl;

	core::pose::Pose in_pose(pose);

	core::scoring::ScoreFunctionOP fa_sfxn (new core::scoring::ScoreFunction() );
	protocols::simple_moves::SwitchResidueTypeSetMoverOP to_fa(new protocols::simple_moves::SwitchResidueTypeSetMover("fa_standard"));
	fa_sfxn->set_weight( core::scoring::fa_rep, 1.0 );
	to_fa->apply(in_pose);

	core::Size last_subunit_pos=0;
	// Create sub-pose of single subunit unit
	for ( core::Size i=1; i<=in_pose.size(); i++ ) {
		std::string res_name = in_pose.residue(i).type().name3();
		if ( res_name == default_bridge_type_ ) {
			last_subunit_pos=i;
		}
	}

	if ( last_subunit_pos > 0 ) {

		core::pose::Pose in_pose_sub = core::pose::Pose(in_pose, 2, last_subunit_pos);

		core::Real omega, rise;
		if ( prop_dir_ == 1 ) {
			omega = omega_;
			rise = rise_;
		} else {
			omega = -omega_;
			rise = -rise_;
		}

		numeric::xyzMatrix< core::Real > rot_mat = numeric::z_rotation_matrix_degrees(omega); // Create rotation matrix around z axis
		numeric::xyzVector< core::Real > trans(0, 0, rise); // Create translation matrix

		fa_sfxn->score(in_pose_sub);
		core::Real fa_rep_in_pose_sub = in_pose_sub.energies().total_energies()[ core::scoring::fa_rep ];

		core::pose::Pose in_pose_rot = core::pose::Pose(in_pose_sub);
		core::pose::Pose trans_pose = core::pose::Pose(in_pose_sub);

		// Adding rotated subunits to dock
		fa_sfxn->score(in_pose_rot);
		core::Real fa_rep_in_pose_rot_s = in_pose_rot.energies().total_energies()[ core::scoring::fa_rep ];


		for ( core::Size i=1; i<=num_repeats_; i++ ) {
			trans_pose.apply_transform_Rx_plus_v(rot_mat, trans);
			append_pose_to_pose(in_pose_rot, trans_pose, true); // add trans_pose to in_pose as a new chain
			if ( i==1 ) {
				fa_sfxn->score(in_pose_rot);
				fa_rep_in_pose_rot_s = in_pose_rot.energies().total_energies()[ core::scoring::fa_rep ];
				fa_rep_in_pose_rot_s = fa_rep_in_pose_rot_s - fa_rep_in_pose_sub;
			}
		}

		fa_sfxn->score(in_pose_rot);
		core::Real fa_rep_in_pose_rot = in_pose_rot.energies().total_energies()[ core::scoring::fa_rep ];

		TR << "fa_rep_in_pose_sub: " << std::to_string(fa_rep_in_pose_sub) << std::endl;
		TR << "fa_rep_in_pose_rot_s: " << std::to_string(fa_rep_in_pose_rot_s) << std::endl;
		TR << "fa_rep_in_pose_sub + fa_rep_in_pose_rot_s*num_repeats_: " << std::to_string(fa_rep_in_pose_sub + fa_rep_in_pose_rot_s*num_repeats_) << std::endl;
		TR << "fa_rep_in_pose_rot: " << std::to_string(fa_rep_in_pose_rot) << std::endl;
		TR << "fa_rep_thresh_:" << std::to_string(fa_rep_thresh_) << std::endl;

		if ( fa_rep_in_pose_rot - (fa_rep_in_pose_sub + fa_rep_in_pose_rot_s*num_repeats_) < fa_rep_thresh_*(num_repeats_+1) ) {
			TR << "Passed rotation clash check" << std::endl;
			return true;
		} else {
			TR << "Passed rotation clash check" << std::endl;
			return false;
		}

	} else {
		TR << "Error! Could not resolve subunit" << std::endl;
		return false;
	}
}

/// @brief parse xml
void
PropagateClashCheckFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	fa_rep_thresh_ = tag->getOption< core::Real >( "fa_rep_thresh", 10.0 );
	default_bridge_type_ = tag->getOption< std::string >( "default_bridge_type", "VAL" );
	omega_ = tag->getOption< core::Real >( "omega", 36.0 );
	rise_ = tag->getOption< core::Real >( "rise", 3.38 );
	num_repeats_ = tag->getOption< core::Size >( "num_repeats", 22 );
	prop_dir_ = tag->getOption< core::Size >( "prop_dir", 0 );
}

protocols::filters::FilterOP
PropagateClashCheckFilterCreator::create_filter() const { return protocols::filters::FilterOP( new PropagateClashCheckFilter ); }

std::string
PropagateClashCheckFilterCreator::keyname() const { return "PropagateClashCheck"; }

void
PropagateClashCheckFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	using namespace utility::tag;
	AttributeList attlist; // XRW TO DO: check
	attlist + XMLSchemaAttribute::attribute_w_default( "fa_rep_thresh" , xsct_real , "Threshold fa_rep clash score." , "10.0" )
		+ XMLSchemaAttribute::attribute_w_default( "default_bridge_type" , xs_string , "Default type for bridge" , "VAL" )
		+ XMLSchemaAttribute::attribute_w_default( "omega" , xsct_real , "Rotation around z axis (per subunit)" , "36.0" )
		+ XMLSchemaAttribute::attribute_w_default( "rise" , xsct_real , "Rise along z axis (per subunit)" , "3.38" )
		+ XMLSchemaAttribute::attribute_w_default( "num_repeats" , xsct_real , "Number of subunit repeats" , "22" )
		+ XMLSchemaAttribute::attribute_w_default( "prop_dir" , xsct_real , "Propagation direction (0: forward, 1: reverse)" , "0" ) ;


	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist ) ;

	protocols::filters::xsd_type_definition_w_attributes( xsd, keyname(), "Checks for clashes between adjacent rotated layers of a subunit", attlist );
}

} // dna_dock
} // protocols

