// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/simple_metrics/ConstraintsMetric.cc
/// @brief A simple metric that writes out all the constraints in a pose or sub-region of a pose,
/// in a format matching the (non-enzdes) constraints file format.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Unit headers
#include <protocols/analysis/simple_metrics/ConstraintsMetric.hh>
#include <protocols/analysis/simple_metrics/ConstraintsMetricCreator.hh>

// Core headers
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraint.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationCollection.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

// C++ STL headers
#include <sstream>

static basic::Tracer TR( "protocols.analysis.simple_metrics.ConstraintsMetric" );


namespace protocols {
namespace analysis {
namespace simple_metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ConstraintsMetric::ConstraintsMetric():
	core::simple_metrics::StringMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
ConstraintsMetric::~ConstraintsMetric(){}

core::simple_metrics::SimpleMetricOP
ConstraintsMetric::clone() const {
	return utility::pointer::make_shared< ConstraintsMetric >( *this );
}

std::string
ConstraintsMetric::name() const {
	return name_static();
}

std::string
ConstraintsMetric::name_static() {
	return "ConstraintsMetric";

}
std::string
ConstraintsMetric::metric() const {
	return "Constraints";
}

void
ConstraintsMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap  )
{
	SimpleMetric::parse_base_tag( tag );

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector( core::select::residue_selector::parse_residue_selector( tag, datamap, "residue_selector" ) );
	}
}

void
ConstraintsMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attributes_for_parse_residue_selector( attlist, "residue_selector", "Selector specifying subset of pose.  Constraints "
		"that involve selected residues (including those that involve some residues that are not selected) will "
		"be written out.  If a residue selector is not provided, all constraints will be written out."
	);

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric that summarizes all constraints in a pose or sub-region of a pose defined by a residue "
		"selector, and dumps the summary in a format matching the (non-enzdes) constraints file format.",
		attlist
	);
}

std::string
ConstraintsMetric::calculate(
	const core::pose::Pose & pose
) const {
	core::select::residue_selector::ResidueSubset const residue_selection(
		residue_selector_ == nullptr ?
		utility::vector1<bool>( pose.total_residue(), true )
		:
		residue_selector_->apply(pose)
	);
	utility::vector1< core::scoring::constraints::ConstraintCOP > cstset( pose.constraint_set()->get_all_constraints() );
	std::ostringstream ss;
	ss << "\n";
	for ( auto const & cst : cstset ) {
		if ( residue_selector_ == nullptr || selected_residue_included_in_cst( *cst, residue_selection ) ) {
			cst->show_def( ss, pose );
		}
	}
	return ss.str();
}

/// @brief This simple metric is unpublished.  It returns Vikram K. Mulligan as its author.
void
ConstraintsMetric::provide_citation_info( basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"ConstraintsMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
		"Vikram K. Mulligan",
		"Systems Biology, Center for Computational Biology, Flatiron Institute",
		"vmulligan@flatironinstitute.org",
		"Wrote the ConstraintsMetric."
		)
	);
}

/// @brief Set the residue selector.  Use nullptr to clear the residue selector.
/// @details Does not clone the input residue selector, but stores the const owning pointer directly.
void
ConstraintsMetric::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	residue_selector_ = selector_in;
}

/// @brief Of the residues that a particular constraint acts on, do any lie in the selected residues?
bool
ConstraintsMetric::selected_residue_included_in_cst(
	core::scoring::constraints::Constraint const & cst,
	core::select::residue_selector::ResidueSubset const & residue_selection
) const {
	//Get the resiues that this constraint acts upon:
	utility::vector1< core::Size > const resnums( cst.residues() );
	//Compare to selection:
	if ( cst.type() == "CoordinateConstraint" ) {
		//Special case: only consider first atom.
		if ( residue_selection[resnums[1]] ) return true;
	} else {
		for ( core::Size const resnum : resnums ) {
			runtime_assert( residue_selection.size() >= resnum );
			if ( residue_selection[resnum] ) return true;
		}
	}
	return false;
}

void
ConstraintsMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ConstraintsMetric::provide_xml_schema( xsd );
}

std::string
ConstraintsMetricCreator::keyname() const {
	return ConstraintsMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
ConstraintsMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< ConstraintsMetric >();
}

} //simple_metrics
} //analysis
} //protocols

#ifdef    SERIALIZATION



template< class Archive >
void
protocols::analysis::simple_metrics::ConstraintsMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::StringMetric>( this ) );
	arc( CEREAL_NVP( residue_selector_ ) );
}

template< class Archive >
void
protocols::analysis::simple_metrics::ConstraintsMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::StringMetric >( this ) );
	arc( residue_selector_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::analysis::simple_metrics::ConstraintsMetric );
CEREAL_REGISTER_TYPE( protocols::analysis::simple_metrics::ConstraintsMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_analysis_simple_metrics_ConstraintsMetric )
#endif // SERIALIZATION





