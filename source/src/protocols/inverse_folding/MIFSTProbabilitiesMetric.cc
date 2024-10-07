// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/inverse_folding/MIFSTProbabilitiesMetric.cc
/// @brief A PerResidueProbabilitiesMetric that stores amino acid probabilities predicted by the MIF-ST model.
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)

// Unit headers
#include <protocols/inverse_folding/MIFSTProbabilitiesMetric.hh>
#include <protocols/inverse_folding/MIFSTProbabilitiesMetricCreator.hh>

// protocols headers
#include <protocols/inverse_folding/MIFST.hh>

// Core headers
#include <core/simple_metrics/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.inverse_folding.MIFSTProbabilitiesMetric" );


namespace protocols {
namespace inverse_folding {

using namespace core::select;
using namespace core::select::residue_selector;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
MIFSTProbabilitiesMetric::MIFSTProbabilitiesMetric():
	core::simple_metrics::PerResidueProbabilitiesMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
MIFSTProbabilitiesMetric::~MIFSTProbabilitiesMetric(){}

core::simple_metrics::SimpleMetricOP
MIFSTProbabilitiesMetric::clone() const {
	return utility::pointer::make_shared< MIFSTProbabilitiesMetric >( *this );
}

std::string
MIFSTProbabilitiesMetric::name() const {
	return name_static();
}

std::string
MIFSTProbabilitiesMetric::name_static() {
	return "MIFSTProbabilitiesMetric";

}
std::string
MIFSTProbabilitiesMetric::metric() const {

	return "MIFSTProbs";
}

/// @brief Get the residue selector.
/// @details If this returns nullptr, it means that no residue selector is being used.
core::select::residue_selector::ResidueSelectorCOP
MIFSTProbabilitiesMetric::residue_selector() const {
	return residue_selector_;
}

void
MIFSTProbabilitiesMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{

	SimpleMetric::parse_base_tag( tag );

	bool multirun = tag->getOption< bool >( "multirun", true );
	set_multirun( multirun );

	bool use_gpu = tag->getOption< bool >( "use_gpu", false );
	set_use_gpu( use_gpu );

	set_residue_selector(
		core::select::residue_selector::parse_residue_selector( tag, datamap, "residue_selector" )
	);
	if ( tag->hasOption("feature_selector") ) {
		set_feature_selector( core::select::residue_selector::parse_residue_selector( tag, datamap, "feature_selector"));
	}

#ifndef USE_TORCH
	utility_exit_with_message( "To use the MIFSTProbabilitiesMetric you need to compile Rosetta with extras=pytorch!" );
#endif //USE_TORCH
}

void
MIFSTProbabilitiesMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attributes_for_parse_residue_selector( attlist, "residue_selector", "A residue selector specifying which residue or residues to predict for." );
	attributes_for_parse_residue_selector( attlist, "feature_selector", "A residue selector specifying which parts of the posed are used as features in prediction.");

	attlist + XMLSchemaAttribute::attribute_w_default( "multirun", xsct_rosetta_bool, "Whether to run multirun the network (one inference pass for all selected residues", "true");

	attlist + XMLSchemaAttribute::attribute_w_default( "use_gpu", xsct_rosetta_bool, "Whether to run the network on the GPU (if one is available)", "false");

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A PerResidueProbabilitiesMetric that stores amino acid probabilities predicted by the MIF-ST model. ", attlist);
}

std::map< core::Size, std::map< core::chemical::AA, core::Real > >
MIFSTProbabilitiesMetric::calculate(const core::pose::Pose & pose) const {

#ifndef USE_TORCH
	(void)pose; // avert error: unused parameter
	utility_exit_with_message( "To use the MIFSTProbabilitiesMetric you need to compile Rosetta with extras=pytorch!" );
#else

    core::select::residue_selector::ResidueSubset selection( pose.total_residue(), true );
    if ( residue_selector_ != nullptr ) {
        selection = residue_selector_->apply(pose);
    }
    core::select::residue_selector::ResidueSubset feature_selection( pose.total_residue(), true );
    if ( selector_two_ != nullptr ) {
        feature_selection = selector_two_->apply(pose);
    }

    std::map< core::Size, std::map< core::chemical::AA, core::Real >> values;
    values = MIFST::get_instance()->sample( pose, selection, feature_selection, multirun_, use_gpu_);

    return values;

#endif //USE_TORCH
}
/// @brief Set the residue selector that we'll be using.
/// @details Passing nullptr results in no residue selector being used.
void
MIFSTProbabilitiesMetric::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	residue_selector_ = selector_in;
}
/// @brief set the optional residue selector for feature selection
void
MIFSTProbabilitiesMetric::set_feature_selector(core::select::residue_selector::ResidueSelectorCOP selector) {
	selector_two_ = std::move(selector);
}

///@brief set the multirun option
void
MIFSTProbabilitiesMetric::set_multirun( bool const multirun ) {
	multirun_ = multirun;
}

///@brief set the gpu option
void
MIFSTProbabilitiesMetric::set_use_gpu( bool const use_gpu ) {
	use_gpu_ = use_gpu;
}

/// @brief This simple metric is unpublished.  It returns Moritz Ertelt as its author.
void
MIFSTProbabilitiesMetric::provide_citation_info( basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"MIFSTProbabilitiesMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
		"Moritz Ertelt",
		"University of Leipzig",
		"moritz.ertelt@gmail.com",
		"Wrote the MIFSTProbabilitiesMetric."
		)
	);
#ifdef USE_TORCH
    citations.add( MIFST::get_MIFST_neural_net_citation() );
#endif //USE_TORCH
}

void
MIFSTProbabilitiesMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	MIFSTProbabilitiesMetric::provide_xml_schema( xsd );
}

std::string
MIFSTProbabilitiesMetricCreator::keyname() const {
	return MIFSTProbabilitiesMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
MIFSTProbabilitiesMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< MIFSTProbabilitiesMetric >();
}

} //inverse_folding
} //protocols


#ifdef    SERIALIZATION



template< class Archive >
void
protocols::inverse_folding::MIFSTProbabilitiesMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueProbabilitiesMetric>( this ) );
	arc( CEREAL_NVP ( residue_selector_ ) );
	arc( CEREAL_NVP (selector_two_));
	arc( CEREAL_NVP ( multirun_ ) );
	arc( CEREAL_NVP ( use_gpu_ ) );

}

template< class Archive >
void
protocols::inverse_folding::MIFSTProbabilitiesMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueProbabilitiesMetric >( this ) );
	core::select::residue_selector::ResidueSelectorOP residue_selector_local_;
	arc( residue_selector_local_);
	residue_selector_ = residue_selector_local_;
	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector );
	selector_two_ = local_selector;
	arc( multirun_ );
	arc( use_gpu_ );

}

SAVE_AND_LOAD_SERIALIZABLE( protocols::inverse_folding::MIFSTProbabilitiesMetric );
CEREAL_REGISTER_TYPE( protocols::inverse_folding::MIFSTProbabilitiesMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_inverse_folding_MIFSTProbabilitiesMetric )
#endif // SERIALIZATION


