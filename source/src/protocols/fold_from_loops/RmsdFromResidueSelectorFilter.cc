// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fold_from_loops/RmsdFromResidueSelectorFilter.cc
/// @brief  Evaluate RMSD between two poses allowing to select the regions to compare in each pose through ResidueSelector
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#include <protocols/fold_from_loops/RmsdFromResidueSelectorFilter.hh>
#include <protocols/fold_from_loops/RmsdFromResidueSelectorFilterCreator.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
// #include <core/conformation/Residue.hh>
#include <core/pose/datacache/cacheable_observers.hh>
// #include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.hh> //Movers_map
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
// #include <core/pose/PDBInfo.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
// #include <numeric/model_quality/rms.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>

// #include <core/id/AtomID.hh>
// #include <core/id/AtomID_Map.hh>
// #include <core/id/NamedAtomID.hh>
// #include <ObjexxFCL/FArray1D.hh>
// #include <ObjexxFCL/FArray2D.hh>

#include <algorithm>
#include <list>
#include <map>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// #include <utility/io/izstream.hh>
// #include <protocols/toolbox/superimpose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace fold_from_loops {

static basic::Tracer TR( "protocols.fold_from_loops.RmsdFromResidueSelectorFilter" );

RmsdFromResidueSelectorFilter::RmsdFromResidueSelectorFilter() :
	protocols::filters::Filter( class_name() ),
	reference_select_( new core::select::residue_selector::TrueResidueSelector ),
	query_select_( new core::select::residue_selector::TrueResidueSelector ),
	threshold_( default_rmsd_threshold() ),
	reference_pose_( /* NULL */ ),
	sensitivity_( 0.001 ),
	CA_only_ ( default_ca_selection() ),
	gdt_( default_gdt_selection() )
{
}

RmsdFromResidueSelectorFilter::~RmsdFromResidueSelectorFilter() {}


core::Real
RmsdFromResidueSelectorFilter::compute( core::pose::Pose const & pose ) const
{

	using namespace core::select::residue_selector;

	core::pose::Pose copy_pose = pose;
	core::pose::Pose native = *reference_pose_;
	core::Real rmsd( 0.0 );


	ResidueSubset query_subset = query_select_->apply( pose );
	ResidueSubset reference_subset = reference_select_->apply( native );

	runtime_assert_msg( count_selected( query_subset ) > 0, "No residues are selected for the query pose" );
	runtime_assert_msg( count_selected( reference_subset ) > 0, "No residues are selected for the reference pose" );

	TR << "Reference selects " << count_selected( reference_subset ) << " residues." << std::endl;
	TR << "Query selects " << count_selected( query_subset ) << " residues." << std::endl;

	grafting::simple_movers::DeleteRegionMover deleter;
	ResidueSelectorCOP query_not( new NotResidueSelector( query_select_ ) ); // We need to turn it the other way arround.
	deleter.set_residue_selector( query_not );
	deleter.apply( copy_pose );

	ResidueSelectorCOP reference_not( new NotResidueSelector( reference_select_ ) ); // We need to turn it the other way arround.
	deleter.set_residue_selector( reference_not );
	deleter.apply( native );

	if ( not gdt_ ) {
		protocols::simple_filters::SelectRmsdEvaluator evaluator( native, "fromResidueSelector", CA_only_ );
		rmsd = evaluator.apply( copy_pose );
	} else {
		protocols::simple_filters::SelectGdtEvaluator evaluator( native, "fromResidueSelector" );
		rmsd = evaluator.apply( copy_pose );
	}

	if ( rmsd < sensitivity_ ) rmsd = 0;

	return rmsd;
}

bool
RmsdFromResidueSelectorFilter::apply( core::pose::Pose const & pose ) const {

	core::Real const rmsd( compute( pose ));
	TR << "RMSD over selected residues: " << rmsd ;
	if ( rmsd <= threshold_ ) {
		TR<<" passing."<<std::endl;
		return( true );
	} else {
		TR<<" failing." << std::endl;
		return( false );
	}
}

void
RmsdFromResidueSelectorFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	out<<"RMSD: " << rmsd<<'\n';
}

core::Real
RmsdFromResidueSelectorFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	return( (core::Real) rmsd );
}

void
RmsdFromResidueSelectorFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & reference_pose )
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

	threshold( tag->getOption<core::Real>( "threshold", default_rmsd_threshold() ) );
	CA_only( tag->getOption<bool>( "CA_only", default_ca_selection() ) );
	GDT( tag->getOption<bool>( "use_gdt", default_gdt_selection() ) );

	reference_selector( core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "reference_selector" ), data_map ) );
	query_selector( core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "query_selector" ), data_map ) );

}

void RmsdFromResidueSelectorFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "CA_only", xsct_rosetta_bool, "When selected, use only CA RMSD", std::to_string( default_ca_selection() ) )
		+ XMLSchemaAttribute::attribute_w_default( "use_gdt", xsct_rosetta_bool, "When selected, use GDTm algorithm", std::to_string( default_gdt_selection() ) )
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "Threshold in RMSD above which the filter fails", std::to_string( default_rmsd_threshold() ) );
	rosetta_scripts::attributes_for_saved_reference_pose( attlist );
	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "reference_selector", "Selector specifying residues to take into account in the reference pose" );
	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "query_selector", "Selector specifying residues to take into account in the query pose" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Filter based on the C-alpha RMSD to a reference structure. ResiduesSelectors can be applied to both reference and query poses.", attlist );
}

std::string RmsdFromResidueSelectorFilterCreator::keyname() const {
	return RmsdFromResidueSelectorFilter::class_name();
}

protocols::filters::FilterOP
RmsdFromResidueSelectorFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new RmsdFromResidueSelectorFilter );
}

void RmsdFromResidueSelectorFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RmsdFromResidueSelectorFilter::provide_xml_schema( xsd );
}

} // fold_from_loops
} // protocols
