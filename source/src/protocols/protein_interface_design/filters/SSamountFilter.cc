// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/SSamountFilter.cc
/// @brief SS fraction filtering
/// @author Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/filters/SSamountFilter.hh>
#include <protocols/protein_interface_design/filters/SSamountFilterCreator.hh>
#include <protocols/filters/Filter.hh>


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <core/pose/PDBInfo.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <numeric/model_quality/rms.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>

//Include DSSP
#include <core/scoring/dssp/Dssp.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <algorithm>
#include <list>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

//Include Rosetta protocols
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/toolbox/superimpose.hh>

//Include Boost
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

//Include ObjexxFCL
#include <ObjexxFCL/FArray.all.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

using namespace ObjexxFCL;

SSamountFilter::SSamountFilter() :
	protocols::filters::Filter( "SSamount" ),
	upper_threshold_(1.0),
	lower_threshold_(0.5),
	target_chain_(0),
	b_target_chain_(false),
	b_discard_lonely_SS_(false)
{
}

SSamountFilter::SSamountFilter(
	core::Real const upper_threshold,
	core::Real const lower_threshold,
	core::Size target_chain,
	bool b_target_chain,
	bool b_discard_lonely_SS) : protocols::filters::Filter( "SSamount" )
{
	upper_threshold_=upper_threshold;
	lower_threshold_=lower_threshold;
	target_chain_=target_chain;
	b_target_chain_=b_target_chain;
	b_discard_lonely_SS_=b_discard_lonely_SS;
}

SSamountFilter::~SSamountFilter() {}

protocols::filters::FilterOP
SSamountFilter::clone() const {
	return protocols::filters::FilterOP( new SSamountFilter( *this ) );
}

static thread_local basic::Tracer TR( "protocols.protein_interface_design.filters.SSamountFilter" );
core::Real
SSamountFilter::compute( core::pose::Pose const & pose ) const
{
	core::Real ssfraction( 0.0 );

	core::pose::Pose target_pose = pose;

	if ( TR.visible() ) TR << boost::format("apply( pose=<%s residues> )") % target_pose.n_residue() << std::endl;
	//Apply to a particular chain

	if ( b_target_chain_ ) {
		TR <<  "Target chain: " << target_chain_ << "in a pose with #Chains: " << pose.conformation().num_chains() << std::endl;
		if ( target_chain_ > pose.conformation().num_chains() ) {
			utility_exit_with_message(" FragmentLookupFilter invalid chain" );
		}
		target_pose = *pose.split_by_chain(target_chain_);
		if ( TR.visible() ) TR << boost::format("apply mod by chain! (Now pose=<%s residues> )") % target_pose.n_residue() << std::endl;
	}
	//Use reduced DSSP notation to calculate SS (notation is 'H','E','L',' ')
	core::scoring::dssp::Dssp dssp_profile(target_pose);
	std::string s_dssp_profile(dssp_profile.get_dssp_secstruct());
	TR << s_dssp_profile << std::endl;

	//This is to count contigous SS elements
	bool did_start_SS=false;
	bool did_end_SS=false;
	bool did_start_second_SS=false;
	//Count SS
	core::Size ssCount=0;
	for ( core::Size i=0; i< s_dssp_profile.length(); i++ ) {
		if ( (s_dssp_profile[i] == 'H') || (s_dssp_profile[i] == 'E') ) {
			ssCount+=1;
			if ( !did_start_SS ) {
				did_start_SS=true;
			} else if ( did_end_SS ) {
				//Improve to actually count
				did_start_second_SS=true;
			}
		} else {
			if ( did_start_SS ) {
				did_end_SS=true;
			}
		}
	}

	//Calculate number of residues participating in SS vs total lenght
	ssfraction=((core::Real)ssCount)/((core::Real)target_pose.n_residue());

	//hack to avoid lonely SS
	if ( b_discard_lonely_SS_ && !did_start_second_SS ) {
		ssfraction = -1.0;
	}

	return ssfraction;
}

//Shared Filter Stuff
bool SSamountFilter::apply( core::pose::Pose const & pose ) const {

	TR << "Calculating SS content." << std::endl;
	core::Real const ssfraction( compute( pose ));
	if ( (ssfraction >= lower_threshold_) && (ssfraction <= upper_threshold_) ) {
		TR << " SS amount/content filter Pass. SS fraction: "<< ssfraction << std::endl;
		return( true );
	}
	TR << "SS amount/content filter Fail. SS fraction: "<< ssfraction << std::endl;
	return( false );


}

//Common Filter Stuff
void SSamountFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const ssfraction( compute( pose ));
	out << "SS fraction: " << ssfraction << std::endl;
	//TR<< "SS fraction: " << ssfraction << std::endl;
}

core::Real SSamountFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const ssfraction( compute( pose ));
	return( (core::Real) ssfraction );
}

void SSamountFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose )
{

	/// @brief
	//Next line is to make the function pass the compiler, FixMe
	pose.n_residue();


	if ( tag->hasOption("reference_name") ) {
		TR.Warning << "reference_name not implemented yet for SSamountFilter, FixMe if you can!" << std::endl;
		utility_exit_with_message("This option hasn't been implemented yet");
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag, data_map );
	}//else if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ){
	//                core::import_pose::pose_from_pdb( *reference_pose_, basic::options::option[ basic::options::OptionKeys::in::file::native ] );
	//}else{
	//        utility_exit_with_message("Not reference structure defined! Use [reference_name] or [-in::file::native] ");
	//}

	/// @brief
	target_chain_ = 0;
	upper_threshold_ = 0.0;
	lower_threshold_ = 0.0;
	b_discard_lonely_SS_ = false;

	if ( tag->hasOption("chain") ) {
		target_chain_ = tag->getOption<core::Size>( "chain" );
		b_target_chain_ = true;
	} else {
		TR.Warning << "No chain specified, using them-all for calculation" << std::endl;
	}

	if ( tag->hasOption("upper_threshold") ) {
		upper_threshold_ = tag->getOption<core::Real>( "upper_threshold");
	} else {
		upper_threshold_=1.0;
		TR.Warning << "Threshold option not specified, using default" << upper_threshold_ << std::endl;
	}
	if ( tag->hasOption("lower_threshold") ) {
		lower_threshold_ = tag->getOption<core::Real>( "lower_threshold");
	} else {
		lower_threshold_=0.5;
		TR.Warning << "Threshold option not specified, using default" << lower_threshold_ << std::endl;
	}
	if ( tag->hasOption("discard_lonely_SS") ) {
		core::Size tmp_val= tag->getOption<core::Size>("discard_lonely_SS");
		if ( tmp_val == 1 ) {
			b_discard_lonely_SS_=true;
			TR.Info << "I'll return -1.0 for matches with a single SS element" << std::endl;
		} else if ( tmp_val == 0 ) {
			b_discard_lonely_SS_=false;
		} else {
			utility_exit_with_message("Invalid discard_lonely_SS value, valid options are 0 and 1 (disabled/enabled)" );
		}
	} else {
		TR.Warning << "discard_lonely_SS option not specified, using default: " << b_discard_lonely_SS_ << std::endl;
	}
}

protocols::filters::FilterOP SSamountFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new SSamountFilter );
}

std::string SSamountFilterCreator::keyname() const {
	return "SSamount";
}


} // filters
} // protein_interface_design
} // devel


