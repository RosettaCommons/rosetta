// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SelectedResiduesPyMOLMetric.cc
/// @brief A utility metric to output a string of selected residues from a residue selector as a pymol selection.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/metrics/SelectedResiduesPyMOLMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/chains_util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.metrics.SelectedResiduesPyMOLMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SelectedResiduesPyMOLMetric::SelectedResiduesPyMOLMetric():
	core::simple_metrics::StringMetric()
{}

SelectedResiduesPyMOLMetric::SelectedResiduesPyMOLMetric( select::residue_selector::ResidueSelectorCOP selector):
	core::simple_metrics::StringMetric()
{
	set_residue_selector( selector );
}
////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SelectedResiduesPyMOLMetric::~SelectedResiduesPyMOLMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
SelectedResiduesPyMOLMetric::SelectedResiduesPyMOLMetric( SelectedResiduesPyMOLMetric const & src ):
	core::simple_metrics::StringMetric( src )
{
	selector_ = src.selector_;
}

core::simple_metrics::SimpleMetricOP
SelectedResiduesPyMOLMetric::clone() const {
	return utility::pointer::make_shared< SelectedResiduesPyMOLMetric >( *this );

}

std::string
SelectedResiduesPyMOLMetric::name() const {
	return name_static();
}

std::string
SelectedResiduesPyMOLMetric::name_static() {
	return "SelectedResiduesPyMOLMetric";

}
std::string
SelectedResiduesPyMOLMetric::metric() const {
	return "pymol_selection";
}

void
SelectedResiduesPyMOLMetric::set_residue_selector(select::residue_selector::ResidueSelectorCOP selector ){
	selector_ = selector;
}

void
SelectedResiduesPyMOLMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &  datamap )
{
	SimpleMetric::parse_base_tag( tag );

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(select::residue_selector::parse_residue_selector( tag, datamap ));
	} else {
		utility_exit_with_message("SelectedResiduesMetric: Residue Selector is required.");
	}

	set_use_rosetta_numbering( tag->getOption<bool>( "use_rosetta_numbering", use_rosetta_numbering() ) );
}

void
SelectedResiduesPyMOLMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector",
		"Required.  Output those residues selected. " );

	attlist + XMLSchemaAttribute::attribute_w_default( "use_rosetta_numbering", xsct_rosetta_bool, "If true, the output uses Rosetta numbers.  If false, it uses PDB numbers.  False by default.", "false" );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"A utility metric that outputs the residue selection as a PyMOL selection. ", attlist);
}

std::string
SelectedResiduesPyMOLMetric::calculate(const pose::Pose & pose) const {
	if ( ! selector_ ) {
		utility_exit_with_message("SelectedResiduesMetricPyMOL: Residue Selector required!");
	} else if ( pose.pdb_info() == nullptr ) {
		utility_exit_with_message("SelectedResiduesMetricPyMOL Requires an intact PDBInfo!");
	} else if ( pose.pdb_info()->obsolete() ) {
		utility_exit_with_message("SelectedResiduesMetricPyMOL Requires an up-to-date PDBInfo!");
	}

	utility::vector1< bool > subset = selector_->apply( pose );
	utility::vector1< core::Size > selected = core::select::get_residues_from_subset( subset );

	//We don't have prefix/suffix, but we have custom_type. Use it here.
	std::string selection;
	if ( get_custom_type() != "" ) {
		selection = "select "+get_custom_type()+", ";
	} else {
		selection = "select rosetta_sele, ";
	}


	std::map< std::string, utility::vector1<std::string >> chain_residues;
	utility::vector1< std::string > chains; //So we can iterate over only the chains, then the residues ala python.
	core::pose::PDBInfoCOP pdb_info( pose.pdb_info() );

	for ( core::Size i : selected ) {
		std::string chain = core::pose::get_chain_from_chain_id( pose.chain(i), pose );
		std::string num = ( pdb_info != nullptr && !use_rosetta_numbering()? utility::to_string( pdb_info->number(i) ) : utility::to_string(i) );
		std::string icode = ( pdb_info != nullptr && !use_rosetta_numbering()? utility::to_string( pdb_info->icode(i) ) : utility::to_string(' ') );

		std::string res = "";
		if ( icode == utility::to_string(' ') ) {
			res = num;
		} else {
			res = num+icode;
		}

		if ( ! chain_residues.count( chain ) ) {
			utility::vector1< std::string > residues;
			chain_residues[chain] = residues;
			chain_residues[chain].push_back( res );
			chains.push_back( chain );
		} else {
			chain_residues[chain].push_back( res );
		}
	}

	for ( core::Size i = 1; i <= chains.size(); ++i ) {

		std::string chain = chains[i];
		std::string subselection = "";

		if ( i == 1 ) {
			subselection = "(chain "+chain+" and resid ";
		} else {
			subselection = subselection + " or "+"(chain "+chain+" and resid ";
		}
		for ( core::Size x = 1; x <= chain_residues[ chain ].size(); ++x ) {
			std::string res = chain_residues[chain][x];
			subselection += res;

			if ( x != chain_residues[chain].size() ) {
				subselection += ",";
			}
		}
		subselection+=")";
		//TR << chain << " " << subselection << std::endl;
		selection += subselection;
	}
	return selection;
}

/// @brief Provide the citation.
void
SelectedResiduesPyMOLMetric::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	basic::citation_manager::UnpublishedModuleInfoOP authors (
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		name(), basic::citation_manager::CitedModuleType::SimpleMetric,
		"Jared Adolf-Bryfogle",
		"Scripps Research Institute",
		"jadolfbr@gmail.com",
		"Wrote the SelectedResiduesPyMOLMetric."
		)
	);
	authors->add_author(
		"Vikram K. Mulligan",
		"Systems Biology, Center for Computational Biology, Flatiron Institute",
		"vmulligan@flatironinstitute.org",
		"Made the SelectedResiduesPyMOLMetric compatible with poses that lack PDBInfo."
	);

	citations.add( authors );
	citations.add( selector_ );
}

void
SelectedResiduesPyMOLMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SelectedResiduesPyMOLMetric::provide_xml_schema( xsd );
}

std::string
SelectedResiduesPyMOLMetricCreator::keyname() const {
	return SelectedResiduesPyMOLMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
SelectedResiduesPyMOLMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< SelectedResiduesPyMOLMetric >();

}

} //core
} //simple_metrics
} //metrics






