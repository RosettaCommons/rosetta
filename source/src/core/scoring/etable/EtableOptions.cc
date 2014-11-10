// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

// Unit headers
#include <core/scoring/etable/EtableOptions.hh>

// AUTO-REMOVED #include <core/scoring/types.hh>
#include <basic/options/option.hh>

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace etable {

EtableOptions::EtableOptions() :
	max_dis( 6.0 ),
	bins_per_A2( 20 ),
	Wradius( 1.0 ),
	lj_switch_dis2sigma( 0.6 ),
	no_lk_polar_desolvation( false ),
	lj_hbond_OH_donor_dis(3.0),
	lj_hbond_hdis(1.95),
	enlarge_h_lj_wdepth(false)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	max_dis = option[ score::fa_max_dis ];
	if( option[ score::no_smooth_etables ] && !option[ score::fa_max_dis ].user() ) {
		basic::T("core.scoring.etable") << "no_smooth_etables requested and fa_max_dis not specified: using 5.5 as default" << std::endl;
		max_dis = 5.5;
	}
	if ( option[ score::no_lk_polar_desolvation ] ) {
		no_lk_polar_desolvation = false;
	}

	lj_hbond_OH_donor_dis = option[ corrections::score::lj_hbond_OH_donor_dis ];
	lj_hbond_hdis = option[ corrections::score::lj_hbond_hdis ];

}

// Another constructor
EtableOptions::EtableOptions( EtableOptions const & src ) :
	ReferenceCount( src )
{
	*this = src;
}

EtableOptions::~EtableOptions(){}

EtableOptions const &
EtableOptions::operator=( EtableOptions const & src )
{
	max_dis = src.max_dis;
	bins_per_A2 = src.bins_per_A2;
	Wradius = src.Wradius;
	lj_switch_dis2sigma = src.lj_switch_dis2sigma;
	no_lk_polar_desolvation = src.no_lk_polar_desolvation;
	lj_hbond_OH_donor_dis = src.lj_hbond_OH_donor_dis;
	lj_hbond_hdis = src.lj_hbond_hdis;
	enlarge_h_lj_wdepth = src.enlarge_h_lj_wdepth;
	return *this;
}

bool
operator==( EtableOptions const & a, EtableOptions const & b )
{
	return (( a.max_dis == b.max_dis ) &&
					( a.bins_per_A2 == b.bins_per_A2 ) &&
					( a.Wradius == b.Wradius ) &&
					( a.lj_switch_dis2sigma == b.lj_switch_dis2sigma ) &&
					( a.no_lk_polar_desolvation == b.no_lk_polar_desolvation ) &&
					( a.lj_hbond_OH_donor_dis == b.lj_hbond_OH_donor_dis ) &&
					( a.lj_hbond_hdis == b.lj_hbond_hdis ) &&
					( a.enlarge_h_lj_wdepth == b.enlarge_h_lj_wdepth ) ) ;
}


////////////////////////////
std::ostream &
operator<< ( std::ostream & out, const EtableOptions & options ){
	options.show( out );
	return out;
}


///
void
EtableOptions::show( std::ostream & out ) const
{
	out <<"EtableOptions::max_dis: " << max_dis << std::endl;
	out <<"EtableOptions::bins_per_A2: " << bins_per_A2 << std::endl;
	out <<"EtableOptions::Wradius: " << Wradius << std::endl;
	out <<"EtableOptions::lj_switch_dis2sigma: " << lj_switch_dis2sigma << std::endl;
	out <<"EtableOptions::no_lk_polar_desolvation: " << no_lk_polar_desolvation << std::endl;
	out <<"EtableOptions::lj_hbond_OH_donor_dis: " << lj_hbond_OH_donor_dis << std::endl;
	out <<"EtableOptions::lj_hbond_hdis: " << lj_hbond_hdis << std::endl;
	out <<"EtableOptions::enlarge_h_lj_wdepth: " << enlarge_h_lj_wdepth << std::endl;
}

void
EtableOptions::parse_my_tag(
	utility::tag::TagCOP tag
) {
	if( tag->hasOption( "lj_hbond_OH_donor_dis" )) {
		lj_hbond_OH_donor_dis = tag->getOption<core::Real>( "lj_hbond_OH_donor_dis" );
	}

	if( tag->hasOption( "lj_hbond_hdis" )) {
		lj_hbond_hdis = tag->getOption<core::Real>( "lj_hbond_hdis" );
	}
}

} // etable
} // scoring
} // core


