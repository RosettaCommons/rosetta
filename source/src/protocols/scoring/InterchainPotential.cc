// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  Statistically derived rotamer pair potentials
/// @details For docking (or between chains) only those residues at the interface
///      and between the two interfaces need to be evaluated
/// @author Monica Berrondo


// Unit headers
#include <protocols/scoring/InterchainPotential.hh>
#include <protocols/scoring/InterfaceInfo.hh>

#include <core/scoring/AtomVDW.hh>
#include <core/scoring/EnvPairPotential.hh>
#include <core/scoring/ScoringManager.hh>

// Package headers

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <basic/database/open.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <utility/io/izstream.hh>

// just for debugging
//#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.scoring.InterchainPotential" );

// Singleton instance and mutex static data members
namespace utility {

using protocols::scoring::InterchainPotential;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< InterchainPotential >::singleton_mutex_{};
template <> std::atomic< InterchainPotential * > utility::SingletonBase< InterchainPotential >::instance_( 0 );
#else
template <> InterchainPotential * utility::SingletonBase< InterchainPotential >::instance_( 0 );
#endif

}


namespace protocols {
namespace scoring {

InterchainPotential *
InterchainPotential::create_singleton_instance()
{
	return new InterchainPotential;
}

InterchainPotential::InterchainPotential() :
	atom_vdw_( core::scoring::ScoringManager::get_instance()->get_AtomVDW( core::chemical::CENTROID ) ) // need to make the table choice configurable
{

	using core::Size;
	// load the data
	Size const max_aa( 20 ); // just the standard aa's for now
	Size const env_log_table_size( 4 );

	std::string tag,line;
	core::chemical::AA aa;

	{ // interchain_env_log
		interchain_env_log_.dimension( max_aa, env_log_table_size );

		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/InterchainPotential/interchain_env_log.txt" );
		while ( getline( stream, line ) ) {
			std::istringstream l(line);
			l >> tag >> aa;
			for ( Size i=1; i<= env_log_table_size; ++i ) {
				l >> interchain_env_log_(aa,i);
			}
			if ( l.fail() || tag != "INT_CHAIN_ENV_LOG:"  ) utility_exit_with_message("bad format for scoring/score_functions/InterchainPotential/interchain_env_log.txt");
		}
	}
	{ // interchain_pair_log
		interchain_pair_log_.dimension( max_aa, max_aa );

		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/InterchainPotential/interchain_pair_log.txt" );
		while ( getline( stream, line ) ) {
			std::istringstream l(line);
			l >> tag >> aa;
			for ( Size i=1; i<= max_aa; ++i ) {
				l >> interchain_pair_log_(aa,i);
			}
			if ( l.fail() || tag != "INT_CHAIN_PAIR_LOG:"  ) utility_exit_with_message("bad format for scoring/score_functions/InterchainPotential/interchain_pair_log.txt");
		}
	}
}

void
InterchainPotential::compute_interface(
	core::pose::Pose & pose
) const
{
	InterfaceInfo & interface( nonconst_interface_from_pose( pose ) );

	/// initialize the cenlist info:
	/// only if they have not been calculated since the last score
	if ( !interface.calculated() ) {

		// initialize values
		interface.initialize();

		// compute interpolated number of neighbors at various distance cutoffs
		interface.calculate( pose );
	}

	interface.calculated() = true;
}

void
InterchainPotential::finalize( core::pose::Pose & pose ) const
{
	using core::scoring::CenListInfo;
	CenListInfo & cenlist( core::scoring::EnvPairPotential::nonconst_cenlist_from_pose( pose ));
	cenlist.calculated() = false;

	InterfaceInfo & interface( nonconst_interface_from_pose( pose ) );
	interface.calculated() = false;
}

////////////////////////////////////////////////////////////////////////////////////
void
InterchainPotential::evaluate_env_score(
	core::pose::Pose const & pose,
	core::conformation::Residue const & rsd,
	core::Real & env_score
) const
{
	int env;
	using core::Real;
	Real const fcen10( core::scoring::EnvPairPotential::cenlist_from_pose( pose ).fcen10(rsd.seqpos()) );

	InterfaceInfo const & interface( interface_from_pose( pose ) );

	//reset env_score
	env_score = 0.0;

	if ( rsd.is_protein() == false ) return;

	bool is_interface = interface.is_interface( rsd );

	if ( is_interface ) {
		if ( fcen10 > 16 ) {
			env = 1;
		} else {
			env = 2;
		}
	} else {
		if ( fcen10 > 16 ) {
			env = 3;
		} else {
			env = 4;
		}
	}

	env_score = interchain_env_log_( rsd.aa(), env );
	//std::cout << " residue: " << rsd.seqpos() << " res(i) " << rsd.aa()
	// << " env " << env << " score " << env_score << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////////
void
InterchainPotential::evaluate_contact_score(
	core::pose::Pose const & pose,
	core::Real & contact_score
) const
{
	protocols::scoring::InterfaceInfo const & interface( interface_from_pose( pose ) );

	//reset contact score
	contact_score = 0;

	//calculate contact score for each jump, sum total is contact score
	using core::Size;
	using core::Real;
	for ( Size i = 1; i <= interface.num_jump(); i++ ) {
		int interface_residues = interface.interface_nres(i);

		Real contact_score_jump = ( 20 - interface_residues ) * 0.5;
		if ( interface_residues == 0 ) contact_score_jump += 2.0;
		if ( interface_residues == 1 ) contact_score_jump += 1.0;
		if ( interface_residues == 2 ) contact_score_jump += 0.5;

		contact_score += contact_score_jump;
	}
	if ( contact_score < -10 ) contact_score = -10;
}

///////////////////////////////////////////////////////////////////////////////////////////////

void
InterchainPotential::evaluate_pair_and_vdw_score(
	core::pose::Pose const & pose,
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::Real & pair_contribution,
	core::Real & vdw_contribution
) const
{
	InterfaceInfo const & interface( interface_from_pose( pose ) );

	pair_contribution = 0.0;
	vdw_contribution = 0.0;

	if ( !rsd1.is_protein() || !rsd2.is_protein() ) return;

	if ( interface.is_pair( rsd1, rsd2) == false ) return;

	//fpd to match pre-48179 only compute when centroids are within 12.05A
	core::conformation::Atom const & cen1 ( rsd1.atom( rsd1.nbr_atom() ) ), cen2 (rsd2.atom( rsd2.nbr_atom() ) );
	core::Real const cendist = cen1.xyz().distance_squared( cen2.xyz() );
	if ( cendist > 12.05*12.05 ) return;

	pair_contribution = interchain_pair_log_( rsd1.aa(), rsd2.aa() );

	using namespace core;

	// calculation for vdw?
	// atoms between two chains are guaranteed to not be bonded
	// no countpair!
	for ( Size i=1, i_end = rsd1.natoms(); i<= i_end; ++i ) {

		Vector const & i_xyz( rsd1.xyz(i) );
		Size const i_type( rsd1.atom_type_index(i) );

		// if virtual atom, continue
		if ( rsd1.is_virtual( i ) ) {
			continue;
		}

		// get VDW radii
		utility::vector1< Real > const & i_atom_vdw( atom_vdw_( i_type ) );

		// iterate over interacting residue
		for ( Size j=1, j_end = rsd2.natoms(); j<= j_end; ++j ) {

			// if virtual atom, continue
			if ( rsd2.is_virtual( j ) ) {
				continue;
			}

			Real const bump_dis( i_atom_vdw[ rsd2.atom_type_index(j) ] );
			Real const clash( bump_dis - i_xyz.distance_squared( rsd2.xyz(j) ) );
			if ( clash > 0.0 ) {
				vdw_contribution += ( clash * clash ) / bump_dis;
			}
		}
	}
}

/// @details Pose must already contain a Interface object or this method will fail
InterfaceInfo const &
InterchainPotential::interface_from_pose( core::pose::Pose const & pose ) const
{
	//using core::pose::datacache::CacheableDataType::INTERFACE_INFO;
	return *( utility::pointer::static_pointer_cast< protocols::scoring::InterfaceInfo const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::INTERFACE_INFO ) ) );
}

/// @details Either returns a non-const reference to the Interface object that already exists
/// in the pose, or creates a new Interface object, places it in the pose, and then returns
/// a non-const reference to it
InterfaceInfo &
InterchainPotential::nonconst_interface_from_pose( core::pose::Pose & pose ) const
{
	//using core::pose::datacache::CacheableDataType::INTERFACE_INFO;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::INTERFACE_INFO ) ) {
		return *( utility::pointer::static_pointer_cast< protocols::scoring::InterfaceInfo > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::INTERFACE_INFO ) ) );
	}
	// else
	InterfaceInfoOP interface( new InterfaceInfo );
	pose.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, interface );
	return *interface;
}

}
}
