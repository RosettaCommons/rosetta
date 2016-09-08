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


#include <core/sequence/ABEGOManager.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>
#include <core/indexed_structure_store/ABEGOHashedFragmentStore.hh>
#include <core/indexed_structure_store/FragmentStore.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoringManager.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//// C++ headers
static THREAD_LOCAL basic::Tracer tr("protocols.filters.LeastNativeLike9merFilter");

namespace protocols {
namespace simple_filters {

// @brief default constructor
LeastNativeLike9merFilter::LeastNativeLike9merFilter():
	Filter( "worst9mer" )
{}

// @brief destructor
LeastNativeLike9merFilter::~LeastNativeLike9merFilter() = default;


/// @brief
LeastNativeLike9merFilter::Real
LeastNativeLike9merFilter::report_sm( const Pose & pose ) const
{
	return  compute( pose );
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
	using namespace core::indexed_structure_store;
	core::sequence::ABEGOManager AM;
	utility::vector1< std::string > abegoSeq = AM.get_symbols( pose,1 );//1 stands for class of ABEGO strings
	Size fragment_length = ABEGOHashedFragmentStore_->get_fragment_length();
	Real worstRmsd = 0;
	Size startRes = 1;
	Size endRes = pose.size()-fragment_length+1;
	core::conformation::symmetry::SymmetryInfoCOP symminfo(nullptr);
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >( pose.conformation()).Symmetry_Info();
		endRes = symminfo->num_independent_residues()-fragment_length+1;
	}
	if ( ignore_terminal_res_ ) {
		startRes = startRes+1;
		endRes = endRes-1;
	}
	for ( Size resid=startRes; resid<=endRes; ++resid ) {
		std::string fragAbegoStr = "";
		for ( Size ii=0; ii<fragment_length; ++ii ) {
			fragAbegoStr += abegoSeq[resid+ii];
		}
		if ( only_helices_ && fragAbegoStr == "AAAAAAAAA" ) {
			Real tmpRmsd = ABEGOHashedFragmentStore_->lookback(pose,resid,fragAbegoStr);
			if ( tmpRmsd>worstRmsd ) {
				worstRmsd = tmpRmsd;
			}
		}
		if ( !only_helices_ ) {
			if ( fragAbegoStr != "AAAAAAAAA" ) {
				Real tmpRmsd = ABEGOHashedFragmentStore_->lookback(pose,resid,fragAbegoStr);
				if ( tmpRmsd>worstRmsd ) {
					worstRmsd = tmpRmsd;
				}
			}
		}
	}
	return(worstRmsd);
}


// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool LeastNativeLike9merFilter::apply(const Pose & pose ) const
{
	Real value = compute( pose );
	tr << "value" << value << "filtered_value_" << filtered_value_ << std::endl;
	if ( value <= filtered_value_ ) {
		tr << "Successfully passed filt " << value << std::endl;
		return true;
	} else {
		tr << "Filter failed current threshold=" << value << "/" << filtered_value_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
LeastNativeLike9merFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	using namespace core::indexed_structure_store;
	// set threshold
	filtered_value_ = tag->getOption<Real>( "threshold", 99.0 );
	if ( filtered_value_!=99.0 ) {
		rmsd_lookup_thresh_ = filtered_value_;
	} else {
		rmsd_lookup_thresh_ = tag->getOption<Real>("rmsd_lookup_threshold",0.40);
	}
	only_helices_ = tag->getOption<bool>( "only_helices", false ); //lower threshold to remove kinked helices
	ABEGOHashedFragmentStore_ = ABEGOHashedFragmentStore::get_instance();
	ABEGOHashedFragmentStore_->set_threshold_distance(rmsd_lookup_thresh_);
	tr << "Structures which has the best fragment RMSD at the worst position greater than " << filtered_value_ << " will be filtered." << std::endl;
	ignore_terminal_res_ = tag->getOption<bool>("ignore_terminal_residue",true);

}

filters::FilterOP
LeastNativeLike9merFilterCreator::create_filter() const { return protocols::filters::FilterOP(new LeastNativeLike9merFilter); }

std::string
LeastNativeLike9merFilterCreator::keyname() const { return "worst9mer"; }
} // filters
} // protocols
