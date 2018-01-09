// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/dna_dock/DNAClashCheckFilter.cc
/// @brief  Checks for fa_rep clashes between docked pose and DNA backbone.
/// @author Carl Walkey (cwalkey@uw.edu)

// Unit Headers
#include <protocols/dna_dock/DNAClashCheckFilter.hh>
#include <protocols/dna_dock/DNAClashCheckFilterCreator.hh>

// Basic Headers
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/out.OptionKeys.gen.hh>

// Core Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
//#include <core/pose/PDBInfo.hh>
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
#include <core/pose/subpose_manipulation_util.hh>

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
DNAClashCheckFilter::DNAClashCheckFilter():
	fa_rep_thresh_( 50.0 ),
	dna_a_path_( "" ),
	dna_b_path_( "" )
{}

// @brief constructor with arguments
DNAClashCheckFilter::DNAClashCheckFilter( core::Real const fa_rep_thresh, std::string dna_a_path, std::string dna_b_path ):
	fa_rep_thresh_( fa_rep_thresh ),
	dna_a_path_( dna_a_path ),
	dna_b_path_( dna_b_path )
{}

// @brief copy constructor
DNAClashCheckFilter::DNAClashCheckFilter( DNAClashCheckFilter const & )= default;

// @brief destructor
DNAClashCheckFilter::~DNAClashCheckFilter() = default;

protocols::filters::FilterOP
DNAClashCheckFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new DNAClashCheckFilter() );
}

protocols::filters::FilterOP
DNAClashCheckFilter::clone() const{
	return protocols::filters::FilterOP( new DNAClashCheckFilter( *this ) );
}

// @brief getters
core::Real DNAClashCheckFilter::get_fa_rep_thresh() const { return fa_rep_thresh_; }
std::string DNAClashCheckFilter::get_dna_a_path() const { return dna_a_path_; }
std::string DNAClashCheckFilter::get_dna_b_path() const { return dna_b_path_; }

// @brief setters
void DNAClashCheckFilter::set_fa_rep_thresh( core::Real fa_rep_thresh ) { fa_rep_thresh_ = fa_rep_thresh; }
void DNAClashCheckFilter::set_dna_a_path( std::string dna_a_path ) { dna_a_path_ = dna_a_path; }
void DNAClashCheckFilter::set_dna_b_path( std::string dna_b_path ) { dna_b_path_ = dna_b_path; }

// @brief returns true if the set of residues defined by the TaskOperations have a no clashes. False otherwise.
bool DNAClashCheckFilter::apply( Pose const & pose ) const
{
	core::pose::Pose in_pose(pose);

	core::pose::Pose dna_a_pose;
	core::import_pose::pose_from_file( dna_a_pose, dna_a_path_, core::import_pose::PDB_file ); // Get pose for DNA chain A

	core::pose::Pose dna_b_pose;
	core::import_pose::pose_from_file( dna_b_pose, dna_b_path_, core::import_pose::PDB_file ); // Get pose for DNA chain B

	core::scoring::ScoreFunctionOP fa_sfxn (new core::scoring::ScoreFunction() );
	protocols::simple_moves::SwitchResidueTypeSetMoverOP to_fa(new protocols::simple_moves::SwitchResidueTypeSetMover("fa_standard"));
	fa_sfxn->set_weight( core::scoring::fa_rep, 1.0 );

	to_fa->apply(in_pose);
	fa_sfxn->score(in_pose);
	core::Real fa_rep = in_pose.energies().total_energies()[ core::scoring::fa_rep ]; // Get fa_rep score of pose without DNA

	// Add DNA chains to pose
	core::pose::append_pose_to_pose(in_pose, dna_a_pose);
	core::pose::append_pose_to_pose(in_pose, dna_b_pose);

	fa_sfxn->score(in_pose); // Re-score with DNA
	core::Real fa_rep_dna = in_pose.energies().total_energies()[ core::scoring::fa_rep ];

	TR << "fa_rep: " << std::to_string(fa_rep) << std::endl;
	TR << "fa_rep_dna: " << std::to_string(fa_rep_dna) << std::endl;
	TR << "fa_rep_threshold: " << std::to_string(fa_rep_thresh_) << std::endl;

	// Check if clashes are under threshold
	if ( (fa_rep_dna-fa_rep) < fa_rep_thresh_ ) {
		TR << "Passed DNA clash check" << std::endl;
		return true;
	} else {
		TR << "Failed DNA clash check" << std::endl;
		return false;
	}
}

/// @brief parse xml
void
DNAClashCheckFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	dna_a_path_ = tag->getOption< std::string >( "dna_a_path", "" );
	dna_b_path_ = tag->getOption< std::string >( "dna_b_path", "" );
	fa_rep_thresh_ = tag->getOption< core::Real >( "fa_rep_thresh", 50.0);
}

protocols::filters::FilterOP
DNAClashCheckFilterCreator::create_filter() const { return protocols::filters::FilterOP( new DNAClashCheckFilter ); }

std::string
DNAClashCheckFilterCreator::keyname() const { return "DNAClashCheck"; }

void
DNAClashCheckFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	using namespace utility::tag;
	AttributeList attlist; // XRW TO DO: check
	attlist + XMLSchemaAttribute::attribute_w_default( "fa_rep_thresh" , xsct_real , "Threshold fa_rep clash score." , "50.0" )
		+ XMLSchemaAttribute::attribute_w_default( "dna_a_path" , xs_string , "File path to DNA A chain" , "" )
		+ XMLSchemaAttribute::attribute_w_default( "dna_b_path" , xs_string , "File path to DNa B chain" , "" ) ;

	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist ) ;

	protocols::filters::xsd_type_definition_w_attributes( xsd, keyname(), "Checks for clashes between docking residues and DNA backbone", attlist );
}

} // dna_dock
} // protocols

