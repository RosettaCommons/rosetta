// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/protein_interface_design/filters/DisulfideFilter.hh
/// @brief Filters for interfaces which could form a disulfide bond between
/// docking partners.
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date Created 4/30/2009
#include <protocols/protein_interface_design/filters/DisulfideFilter.hh>
#include <protocols/protein_interface_design/filters/DisulfideFilterCreator.hh>


// Project Headers
#include <core/types.hh>


//parsing
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/util.hh>
#include <protocols/protein_interface_design/movers/DisulfideMover.hh>

// Utility Headers

// Unit Headers

// C++ headers

#include <core/scoring/disulfides/CentroidDisulfidePotential.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

using namespace core;
using namespace core::scoring;
using core::pose::Pose;
using namespace std;
using std::pair;
using utility::vector1;
using utility::vector1_int;
using protocols::protein_interface_design::movers::DisulfideMover;

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.filters.DisulfideFilter" );

const core::scoring::disulfides::CentroidDisulfidePotential DisulfideFilter::potential_;

/// @brief default ctor
DisulfideFilter::DisulfideFilter() :
	parent( "DisulfideFilter" ),
	targets_(),
	rb_jump_(1)
{}

/// @brief copy ctor
DisulfideFilter::DisulfideFilter(DisulfideFilter const& ) = default;

/// @brief Constructor with a single target residue
DisulfideFilter::DisulfideFilter( core::Size targetResidue ) :
	parent( "DisulfideFilter" ),
	targets_( new core::select::residue_selector::ResidueIndexSelector( utility::to_string(targetResidue) ) ),
	rb_jump_(1)
{}

/// @brief Constructor with multiple target residues
DisulfideFilter::DisulfideFilter( core::select::residue_selector::ResidueSelectorCOP targetResidues ) :
	parent( "DisulfideFilter" ),
	targets_(targetResidues),
	rb_jump_(1)
{}

DisulfideFilter::~DisulfideFilter() {}

/// @return Whether a disulfide bond is possible between any of the targets
bool DisulfideFilter::apply(Pose const & pose ) const
{
	vector1< pair<Size,Size> > disulfides;
	runtime_assert( targets_ );
	vector1< core::Size > targets = core::select::get_residues_from_subset( targets_->apply( pose ) );
	DisulfideMover::disulfide_list(pose, targets, rb_jump_, disulfides);
	if ( disulfides.empty() ) {
		TR << "Failing."<<std::endl;
		return false;
	} else {
		TR << "Passing."<<std::endl;
		return true;
	}
}

void DisulfideFilter::report( ostream & out, Pose const & pose ) const
{
	vector1< pair<Size,Size> > disulfides;
	runtime_assert( targets_ );
	vector1< core::Size > targets = core::select::get_residues_from_subset( targets_->apply( pose ) );
	DisulfideMover::disulfide_list(pose, targets, rb_jump_, disulfides);

	out << disulfides.size() << " disulfides possible: ";
	for ( auto const & disulf : disulfides ) {
		out << disulf.first << '-' << disulf.second << ", ";
	}
	out << "\n";
}
/// @return The number of disulfides possible
Real DisulfideFilter::report_sm( Pose const & pose ) const
{
	vector1< pair<Size,Size> > disulfides;
	runtime_assert( targets_ );
	vector1< core::Size > targets = core::select::get_residues_from_subset( targets_->apply( pose ) );
	DisulfideMover::disulfide_list(pose, targets, rb_jump_, disulfides);

	return disulfides.size();
}

/**
* @details Parameters recognized:
*  - targets. A list of possible target residues, seperated by commas.
*/
void DisulfideFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{

	// Set target to the residue specified by "target_pdb_num" or "target_res_num"
	if ( tag->hasOption("targets") ) {
		targets_ = core::pose::get_resnum_selector(tag, "targets");
	}
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP DisulfideFilterCreator::create_filter() const { return protocols::filters::FilterOP( new DisulfideFilter ); }

// XRW TEMP std::string
// XRW TEMP DisulfideFilterCreator::keyname() const { return "DisulfideFilter"; }

std::string DisulfideFilter::name() const {
	return class_name();
}

std::string DisulfideFilter::class_name() {
	return "DisulfideFilter";
}

void DisulfideFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "targets", xsct_refpose_enabled_residue_number_cslist, "Residues to evaluate for possible cross-interface disulfides" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string DisulfideFilterCreator::keyname() const {
	return DisulfideFilter::class_name();
}

protocols::filters::FilterOP
DisulfideFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new DisulfideFilter );
}

void DisulfideFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DisulfideFilter::provide_xml_schema( xsd );
}



} // filters
} // protein_interface_design
} // devel
