// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/SetupAPBSMover.hh
/// @brief Setup for Poisson-Boltzmann energy term use in scorefunction
/// This mover assumes you have a path to the APBS (Adaptive Poisson-Boltzmann Solver) program
/// in the system PATH.  The executable name must be "apbs", case-sensitive.
/// @author Sachko Honda (honda@apl.washington.edu)

// basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <protocols/pb_potential/SetupPoissonBoltzmannPotential.hh>
#include <protocols/pb_potential/SetupPoissonBoltzmannPotentialCreator.hh>
#include <protocols/simple_moves/ddG.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/PoissonBoltzmannEnergy.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/datacache/DataMap.hh>

// command line options
#include <basic/options/option.hh>
#include <basic/options/keys/pb_potential.OptionKeys.gen.hh>

#include <string>
#include <fstream> // for ifstream
#include <cstdio> // for remove()

namespace protocols{
namespace pb_potential{

typedef SetupPoissonBoltzmannPotential SetupPB;
typedef SetupPoissonBoltzmannPotentialCreator SetupPBCreator;
typedef core::Size Size;

static thread_local basic::Tracer TR( "protocols.pb_potential.SetupPoissonBoltzmannPotential" );

SetupPBCreator::SetupPoissonBoltzmannPotentialCreator()
{}
SetupPBCreator::~SetupPoissonBoltzmannPotentialCreator()
{}
protocols::moves::MoverOP
SetupPBCreator::create_mover() const
{
  return new SetupPoissonBoltzmannPotential;
}
std::string
SetupPBCreator::keyname() const
{
  return SetupPBCreator::mover_name();
}
std::string
SetupPBCreator::mover_name()
{
  return "SetupPoissonBoltzmannPotential";
}

const std::string SetupPB::DEFAULT_APBS_PATH = "apbs";

SetupPB::SetupPoissonBoltzmannPotential()
{}

SetupPB::~SetupPoissonBoltzmannPotential() {}

void
SetupPB::apply(core::pose::Pose & pose ) {

	using namespace core;
	using namespace scoring;
	using namespace methods;

	// Register the empty cache holder if not done so yet.
	if( !pose.data().has( pose::datacache::CacheableDataType::PB_LIFETIME_CACHE ) ){
		PoissonBoltzmannEnergy::PBLifetimeCacheOP new_cache( new PoissonBoltzmannEnergy::PBLifetimeCache() );
		pose.data().set( pose::datacache::CacheableDataType::PB_LIFETIME_CACHE, new_cache );
	}

	// Cache the "which chain" info
	PoissonBoltzmannEnergy::PBLifetimeCacheOP cached_data =
		static_cast< PoissonBoltzmannEnergy::PBLifetimeCacheOP >(pose.data().get_ptr<	PoissonBoltzmannEnergy::PBLifetimeCache>(pose::datacache::CacheableDataType::PB_LIFETIME_CACHE ));
	pose.data().set( pose::datacache::CacheableDataType::PB_LIFETIME_CACHE, cached_data );


	// Let's clean up the previous run's mess.
	remove("*.dx");
  remove("*.in");
	remove("*.pqr");

	// Prescore to cache bound/unbound poses.  This is necessary for filters.
	// Bound, unbound
 	ddg_->apply(pose);

}

std::string
SetupPB::get_name() const {
  return "SetupPoissonBoltzmannPotential";
}
protocols::moves::MoverOP
SetupPB::clone() const {
  return new SetupPoissonBoltzmannPotential( *this );
}

void
SetupPB::parse_my_tag( utility::tag::TagCOP tag,
			    basic::datacache::DataMap & data_map,
			    protocols::filters::Filters_map const & filters_map,
			    protocols::moves::Movers_map const & movers_map,
			    core::pose::Pose const & pose ) {
	// This param is required when the app is NOT linked against the apbs libraries.
	// Validate only when it is a requirement, but register in any way.
	std::string apbs_path;  // path to the apbs executable.
	if( tag->hasOption( "apbs_path" ) ) {
		apbs_path = tag->getOption<std::string>("apbs_path");
	}
#ifdef LINK_APBS_LIBS
	std::ifstream apbsstream( apbs_path.c_str() );
	if( !apbsstream.good() ){
		TR << "APBS not found.  Check the path: " << apbs_path << std::endl;
		TR.flush();
		runtime_assert(false);
	}
	apbsstream.close();
#endif
	basic::options::option[basic::options::OptionKeys::pb_potential::apbs_path](apbs_path);

	utility::vector1<Size> charged_chains;
  if( tag->hasOption( "charged_chains" )) {
    // comma delimited list of residue numbers in string (e.g. "1,2" for chain 1 and chain 2).
    utility::vector1<std::string> temp = utility::string_split( tag->getOption< std::string >( "charged_chains" ), ',');
    for( core::Size i=1; i<=temp.size(); ++i ) {
      charged_chains.push_back(atoi(temp[i].c_str()));
    }
	}
	else{
		TR << "No user defined charged chains.  Default to : 1" << std::endl;
		charged_chains.push_back(1);
	}
	basic::options::option[basic::options::OptionKeys::pb_potential::charged_chains](charged_chains);

	utility::vector1<Size> revamp_near_chain;
	if( tag->hasOption("revamp_near_chain") ) {
		 utility::vector1<std::string> temp = utility::string_split( tag->getOption< std::string >( "revamp_near_chain" ), ',');
    for( core::Size i=1; i<=temp.size(); ++i ) {
      revamp_near_chain.push_back(atoi(temp[i].c_str()));
    }
		basic::options::option[basic::options::OptionKeys::pb_potential::revamp_near_chain]( revamp_near_chain);
	}

	core::Real potential_cap;
	if( tag->hasOption("potential_cap") ) {
		potential_cap = tag->getOption<core::Real>( "potential_cap" );
		basic::options::option[basic::options::OptionKeys::pb_potential::potential_cap]( potential_cap );
	}

	bool sidechain_only;
	if( tag->hasOption("sidechain_only") ) {
		sidechain_only = tag->getOption<bool>( "sidechain_only" );
		basic::options::option[basic::options::OptionKeys::pb_potential::sidechain_only]( sidechain_only );
	}

	core::Real epsilon;
	if( tag->hasOption("epsilon") ) {
		epsilon = tag->getOption<core::Real>( "epsilon" );
		basic::options::option[basic::options::OptionKeys::pb_potential::epsilon]( epsilon );
	}
	bool calcenergy;
	if( tag->hasOption("calcenergy") ) {
		calcenergy = tag->getOption<bool>( "calcenergy" );
		basic::options::option[basic::options::OptionKeys::pb_potential::calcenergy]( calcenergy );
	}
	int apbs_debug;
	if( tag->hasOption("apbs_debug") ) {
		apbs_debug = tag->getOption<int>( "apbs_debug" );
		basic::options::option[basic::options::OptionKeys::pb_potential::apbs_debug]( apbs_debug );
	}
	//-------------------------------------------------------------------------
	// Initialize DDG for pre-scoring, which compute bound & unbound energies.
	//-------------------------------------------------------------------------
	ddg_ = new protocols::simple_moves::ddG();

	// Must turn this ON to enable caculation of bound/unbound states.
	// HIGHLY ILLEGAL CONST_CAST TO MODIFY THE INPUT TAG!
	// INSTEAD, THIS CODE SHOULD REQUIRE THAT "repack 1" IS SET
	// IN THE INPUT, AND FAIL IF IT IS NOT SET.
	( utility::pointer::const_pointer_cast< utility::tag::Tag > (tag) )->setOption<bool>("repack",1);

	std::string scorefxn_name = tag->getOption<std::string>("scorefxn");
	if( scorefxn_name != "" ) {
		core::scoring::ScoreFunction * scorefxn = data_map.get<core::scoring::ScoreFunction*>("scorefxns", scorefxn_name);
		TR << "Scorefxn weigths: " << scorefxn->serialize_weights() << std::endl;
		if( scorefxn->get_weight(core::scoring::PB_elec) == 0. ){
			TR.Error << "PB_elec term is required.  Not found in the scorefxn.  Terminating the program..." << std::endl;
			TR.Error.flush();
runtime_assert(false);
		}
	}
	ddg_->parse_my_tag( tag, data_map, filters_map, movers_map, pose );

}
protocols::moves::MoverOP
SetupPB::fresh_instance() const {

  return new SetupPoissonBoltzmannPotential();
}
}
}
