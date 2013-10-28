// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/FavorNativeResiduePreCycle.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/movers/FavorNativeResiduePreCycle.hh>
#include <protocols/protein_interface_design/movers/FavorNativeResiduePreCycleCreator.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
// AUTO-REMOVED #include <boost/foreach.hpp>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#define foreach BOOST_FOREACH


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.FavorNativeResiduePreCycle" );

std::string FavorNativeResiduePreCycleCreator::keyname() const
{
        return FavorNativeResiduePreCycleCreator::mover_name();
}

protocols::moves::MoverOP
FavorNativeResiduePreCycleCreator::create_mover() const {
        return new FavorNativeResiduePreCycle;
}

std::string
FavorNativeResiduePreCycleCreator::mover_name() {
        return "FavorNativeResidue";
}


std::string
FavorNativeResiduePreCycle::get_name() const {
	return "FavorNativeResidue";
}

FavorNativeResiduePreCycle::~FavorNativeResiduePreCycle() {}

void
FavorNativeResiduePreCycle::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	using namespace utility::pointer;

	bonus_ = tag->getOption<core::Real>( "bonus", 1.5 );
	for( std::map< std::string, ReferenceCountOP >::const_iterator it = (data)[ "scorefxns" ].begin(); it!=(data)[ "scorefxns" ].end(); ++it ){
		ScoreFunctionOP scorefxn( *data.get< ScoreFunction * >( "scorefxns", it->first ) );
		if( scorefxn->get_weight( res_type_constraint ) == 0.0 ){
			scorefxn->set_weight( res_type_constraint, bonus_ );
			TR<<"Setting res_type_constraint weight in scorefxn "<<it->first<<" to "<<bonus_<<'\n';
		}
	}
/*
for( std::map< std::string, ReferenceCountOP >::const_iterator it=(data)[ "scorefxns_hshash" ].begin(); it!=(data)[ "scorefxns" ].end(); ++it ){ // scorefxns where the user defined hs_hash but not fnr
		ScoreFunctionOP scorefxn( *data.get< ScoreFunction * >( "scorefxns", it->first ) );
		scorefxn->set_weight( res_type_constraint, bonus_ );
		TR<<"Setting res_type_constraint weight in scorefxn "<<it->first<<" to "<<bonus_<<'\n';
	}
*/
	TR<<"applying favor native residue to pose with weight: "<<bonus_<<std::endl;
}

} //movers
} //protein_interface_design
} //protocols

