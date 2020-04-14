// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/LeastNativeLike9merFilter
/// @brief filter structures by IntraRepeatContacts
/// @details
/// @author TJ Brunette

// Unit Headers
#include <protocols/simple_filters/LeastNativeLike9merFilter.hh>
#include <protocols/simple_filters/LeastNativeLike9merFilterCreator.hh>

#include <numeric/xyzVector.hh>

// Project Headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>


#include <core/scoring/dssp/Dssp.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>
#include <protocols/indexed_structure_store/SSHashedFragmentStore.hh>
#include <protocols/indexed_structure_store/FragmentStore.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>


#include <protocols/jd2/util.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <numeric/util.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


//// C++ headers
static basic::Tracer TR("protocols.filters.LeastNativeLike9merFilter");

namespace protocols {
namespace simple_filters {

using core::Size;

// @brief default constructor
LeastNativeLike9merFilter::LeastNativeLike9merFilter():
	Filter( "worst9mer" ),
	report_mean_median_( false )
{}

// @brief destructor
LeastNativeLike9merFilter::~LeastNativeLike9merFilter() = default;

void
LeastNativeLike9merFilter::set_residue_selector( core::select::residue_selector::ResidueSelector const & selector ) {
	residue_selector_ = selector.clone();
}


/// @brief
LeastNativeLike9merFilter::Real
LeastNativeLike9merFilter::report_sm( const Pose & pose ) const
{
	return compute( pose );
}

/// @brief
void
LeastNativeLike9merFilter::report( std::ostream & out, Pose const & pose ) const
{
	out << "worst9mer: " <<  compute( pose ) << std::endl;
}


/// @brief
LeastNativeLike9merFilter::Real
LeastNativeLike9merFilter::compute( const Pose & pose ) const
{
	using namespace protocols::indexed_structure_store;
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	std::string dssp_string = dssp.get_dssp_secstruct();
	core::Size fragment_length = SSHashedFragmentStore_->get_fragment_length();
	utility::vector1< Real > rmsds;
	core::Size startRes = 1;
	core::Size endRes = pose.size()-fragment_length+1;
	core::conformation::symmetry::SymmetryInfoCOP symminfo(nullptr);
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >( pose.conformation()).Symmetry_Info();
		endRes = symminfo->num_independent_residues()-fragment_length+1;
	}
	if ( ignore_terminal_res_ ) {
		startRes = startRes+1;
		endRes = endRes-1;
	}
	core::select::residue_selector::ResidueSubset subset(pose.size(), true);
	if ( residue_selector_ ) subset = residue_selector_->apply(pose);

	for ( core::Size resid=startRes; resid<=endRes; ++resid ) {

		bool all_within_subset = true;
		for ( core::Size jj=0; jj<fragment_length; ++jj ) { if ( ! subset[resid+jj] ) all_within_subset=false; }
		if ( ! all_within_subset ) continue;

		std::string frag_ss = dssp_string.substr(resid-1,fragment_length);
		if ( only_helices_ && frag_ss == "HHHHHHHHH" ) {
			Real tmpRmsd = SSHashedFragmentStore_->lookback(pose,resid,frag_ss,false);
			if ( tmpRmsd>rmsd_lookup_thresh_ ) {
				TR << "position:" << resid << " rmsd:" << tmpRmsd <<std::endl;
			}
			rmsds.push_back(tmpRmsd);
		}
		if ( !only_helices_ ) {
			if ( frag_ss != "HHHHHHHHH" ) {
				Real tmpRmsd = SSHashedFragmentStore_->lookback_account_for_dssp_inaccuracy(pose,resid,frag_ss,true,rmsd_lookup_thresh_);
				if ( tmpRmsd>rmsd_lookup_thresh_ ) {
					TR << "position:" << resid << " rmsd:" << tmpRmsd <<std::endl;
				}
				rmsds.push_back(tmpRmsd);
			}
		}
	}

	if ( rmsds.size() == 0 ) rmsds.push_back( 0 );
	if ( report_mean_median_ ) {
		Real mean = numeric::mean( rmsds );
		Real median = numeric::median( rmsds );
		write_mean_median( pose, mean, median );
	}

	Real worstRmsd = numeric::max( rmsds );
	return(worstRmsd);
}


// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool LeastNativeLike9merFilter::apply(const Pose & pose ) const
{
	Real value = compute( pose );
	TR << "value" << value << "filtered_value_" << filtered_value_ << std::endl;
	if ( value <= filtered_value_ ) {
		TR << "Successfully passed filt " << value << std::endl;
		return true;
	} else {
		TR << "Filter failed current threshold=" << value << "/" << filtered_value_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
LeastNativeLike9merFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data
)
{
	using namespace protocols::indexed_structure_store;
	// set threshold
	filtered_value_ = tag->getOption<Real>( "threshold", 99.0 );
	if ( filtered_value_!=99.0 ) {
		rmsd_lookup_thresh_ = filtered_value_;
	} else {
		rmsd_lookup_thresh_ = tag->getOption<Real>("rmsd_lookup_threshold",0.40);
	}
	only_helices_ = tag->getOption<bool>( "only_helices", false ); //lower threshold to remove kinked helices
	SSHashedFragmentStore_ = SSHashedFragmentStore::get_instance();
	SSHashedFragmentStore_->set_threshold_distance(rmsd_lookup_thresh_);
	TR << "Structures which has the best fragment RMSD at the worst position greater than " << filtered_value_ << " will be filtered." << std::endl;
	report_mean_median_ = tag->getOption<bool>( "report_mean_median", report_mean_median_ );
	ignore_terminal_res_ = tag->getOption<bool>("ignore_terminal_residue",true);
	if ( tag->hasOption( "residue_selector" ) ) {
		core::select::residue_selector::ResidueSelectorCOP selector( core::select::residue_selector::parse_residue_selector( tag, data, "residue_selector" ) );
		if ( selector != nullptr ) {
			set_residue_selector( *selector );
		}
	}
}


void
LeastNativeLike9merFilter::write_mean_median( Pose const & , core::Real mean, core::Real median ) const {

	TR << "Mean: " << mean << " Median: " << median << std::endl;
	{
		std::string column_header = this->get_user_defined_name() + "_mean";
		protocols::jd2::add_string_real_pair_to_current_job( column_header, mean );
	}
	{
		std::string column_header = this->get_user_defined_name() + "_median";
		protocols::jd2::add_string_real_pair_to_current_job( column_header, median );
	}
}



std::string LeastNativeLike9merFilter::name() const {
	return class_name();
}

std::string LeastNativeLike9merFilter::class_name() {
	return "worst9mer";
}

void LeastNativeLike9merFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("threshold", xsct_real, "filter threshold", "99.0")
		+ XMLSchemaAttribute::attribute_w_default("rmsd_lookup_threshold", xsct_real, "if the filtered threshold is not 99, the filtered threshold is equal to the rmsd_lookup_threshold; else it is the default value of 0.40", "0.40")
		+ XMLSchemaAttribute::attribute_w_default("only_helices", xsct_rosetta_bool, "look at only helices?", "false")
		+ XMLSchemaAttribute::attribute_w_default("report_mean_median", xsct_rosetta_bool, "Calculate the mean and median rmsd as well?", "false")
		+ XMLSchemaAttribute( "residue_selector", xs_string, "Only look at worst fragment for residues within residue selector" )
		+ XMLSchemaAttribute::attribute_w_default("ignore_terminal_resiude", xsct_rosetta_bool, "ignore the terminal residue?", "true");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Filters out structures that don't have good fragment RMSD at the worst position.", attlist );
}

std::string LeastNativeLike9merFilterCreator::keyname() const {
	return LeastNativeLike9merFilter::class_name();
}

protocols::filters::FilterOP
LeastNativeLike9merFilterCreator::create_filter() const {
	return utility::pointer::make_shared< LeastNativeLike9merFilter >();
}

void LeastNativeLike9merFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LeastNativeLike9merFilter::provide_xml_schema( xsd );
}

} // filters
} // protocols
