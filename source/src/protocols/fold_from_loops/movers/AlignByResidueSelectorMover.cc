// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/AlignByResidueSelectorMover.cc
/// @brief Aligns two poses through the selected residues
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#include <protocols/fold_from_loops/movers/AlignByResidueSelectorMover.hh>
#include <protocols/fold_from_loops/movers/AlignByResidueSelectorMoverCreator.hh>
#include <protocols/fold_from_loops/utils/utils.hh>
#include <protocols/fold_from_loops/selectors/ConstraintResidueSelector.hh>

// Protocol headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/superimpose.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/selection.hh>
#include <core/kinematics/util.hh>
#include <core/import_pose/import_pose.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.fold_from_loops.AlignByResidueSelectorMover" );

namespace protocols {
namespace fold_from_loops {
namespace movers {

AlignByResidueSelectorMover::AlignByResidueSelectorMover():
	protocols::moves::Mover( mover_name() ),
	reference_select_( new core::select::residue_selector::TrueResidueSelector ),
	query_select_( new core::select::residue_selector::TrueResidueSelector ),
	reference_pose_( new core::pose::Pose )
{}

AlignByResidueSelectorMover::~AlignByResidueSelectorMover()= default;

void
AlignByResidueSelectorMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & reference_pose )
{
	/// Call the SavePoseMover, then it can be used from here... or from in:native
	if ( tag->hasOption("reference_name") ) {
		reference_pose_ = rosetta_scripts::saved_reference_pose(tag, data_map );
		TR<<"Loaded reference pose: "<<tag->getOption< std::string >( "reference_name" )<< " with " << reference_pose_->size() << " residues" << std::endl;
	} else {
		reference_pose_ = core::pose::PoseOP( new core::pose::Pose( reference_pose ) );
		if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
			core::import_pose::pose_from_file( *reference_pose_, basic::options::option[ basic::options::OptionKeys::in::file::native ] , core::import_pose::PDB_file);
		}
	}
	TR.Trace << TR.Green << tag->hasOption("reference_name") << " taken from " << reference_pose_ << TR.Reset << std::endl;
	reference_selector( core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "reference_selector" ), data_map ) );
	query_selector( core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "query_selector" ), data_map ) );
}

protocols::moves::MoverOP
AlignByResidueSelectorMover::clone() const
{
	return protocols::moves::MoverOP( new AlignByResidueSelectorMover( *this ) );
}

protocols::moves::MoverOP
AlignByResidueSelectorMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new AlignByResidueSelectorMover );
}

utility::vector1< numeric::xyzVector< core::Real > >
Ca_coord( core::pose::Pose const & pose, core::select::residue_selector::ResidueSubset const & positions ) {
	utility::vector1< numeric::xyzVector< core::Real > > coords;

	coords.clear();
	for ( core::Size i=1; i<=positions.size(); ++i ) {
		if ( positions[i] ) {
			TR.Trace << "    Add coordinates residue " << i << std::endl;
			coords.push_back( pose.residue( i ).xyz( "CA" ) );
		}
	}
	return coords;
}

void
AlignByResidueSelectorMover::apply( core::pose::Pose & pose )
{
	using namespace core::select::residue_selector;
	using namespace protocols::toolbox;

	core::pose::Pose native;
	native.detached_copy( *reference_pose_ );

	ResidueSubset query_subset     = query_select_->apply( pose );
	ResidueSubset reference_subset = reference_select_->apply( native );

	TR.Trace << "Evaluate selectors" << std::endl;
	TR.Trace << "query:     " << represent_residue_selector(query_subset) << std::endl;
	TR.Trace << "reference: " << represent_residue_selector(reference_subset) << std::endl;
	runtime_assert_msg( count_selected( query_subset ) > 0, "No residues are selected for the query pose" );
	runtime_assert_msg( count_selected( reference_subset ) > 0, "No residues are selected for the reference pose" );
	runtime_assert_msg( count_selected( query_subset ) == count_selected( reference_subset ),
		"Selections must have the same number of residues.");

	TR.Trace << "Get CA coordinates for query" << std::endl;
	utility::vector1< numeric::xyzVector< core::Real > > init_coords( Ca_coord( pose, query_subset ) );
	TR.Trace << "Get CA coordinates for reference" << std::endl;
	utility::vector1< numeric::xyzVector< core::Real > >  ref_coords( Ca_coord( native, reference_subset ) );

	TR.Trace << "Calculating the Matrices" << std::endl;
	numeric::xyzMatrix< core::Real > rotation;
	numeric::xyzVector< core::Real > to_init_center, to_fit_center;
	superposition_transform( init_coords, ref_coords, rotation, to_init_center, to_fit_center );
	TR.Trace << "Apply superimposition" << std::endl;
	apply_superposition_transform( pose, rotation, to_init_center, to_fit_center );
}

std::string AlignByResidueSelectorMover::get_name() const {
	return mover_name();
}

std::string AlignByResidueSelectorMover::mover_name() {
	return "AlignByResidueSelectorMover";
}

void AlignByResidueSelectorMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_saved_reference_pose( attlist );
	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "reference_selector", "Selector specifying residues to take into account in the reference pose" );
	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "query_selector", "Selector specifying residues to take into account in the query pose" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Aligns two poses through the selected residues", attlist );
}

std::string AlignByResidueSelectorMoverCreator::keyname() const {
	return AlignByResidueSelectorMover::mover_name();
}

protocols::moves::MoverOP
AlignByResidueSelectorMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AlignByResidueSelectorMover );
}

void AlignByResidueSelectorMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AlignByResidueSelectorMover::provide_xml_schema( xsd );
}

}
} //protocols
} //fold_from_loops
