// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MembraneData.cc
///
/// @brief  Membrane Scorefunction Statistics Class
/// @details Access to centroid rotamer pair potential and membrane
///    environemnt potential statistics. Loads data in database tables
///    in construction. All immutable data - no setters please and const!.
///    Last Modified: 5/13/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MembraneData_cc
#define INCLUDED_core_scoring_membrane_MembraneData_cc

// Unit Headers
#include <core/scoring/membrane/MembraneData.hh>
#include <core/scoring/EnvPairPotential.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

#include <basic/database/open.hh>
#include <utility/io/izstream.hh>

#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>


static THREAD_LOCAL basic::Tracer TR( "core.scoring.membrane.MembraneData" );

namespace core {
namespace scoring {
namespace membrane {

/// @brief Default Constructor
MembraneData::MembraneData() :
	core::scoring::EnvPairPotential(),

	//cems transition regions between environment bins
	//cems transition is from +/- sqrt(36+pad6) +/- sqrt(100+pad10) etc
	cen_dist5_pad_( 0.5 ),
	cen_dist6_pad_( 0.6 ),
	cen_dist7_pad_( 0.65 ),
	cen_dist10_pad_( 1.0 ),
	cen_dist12_pad_( 1.2 ),

	cen_dist5_pad_plus_( cen_dist5_pad_  + 25.0 ),
	cen_dist6_pad_plus_( cen_dist6_pad_ + 36.0 ),
	cen_dist7_pad_plus_( cen_dist7_pad_  + 56.25 ),
	cen_dist10_pad_plus_( cen_dist10_pad_ + 100.0 ),
	cen_dist12_pad_plus_( cen_dist12_pad_ + 144.0 ),

	cen_dist5_pad_minus_( cen_dist5_pad_  - 25.0 ),
	cen_dist7_pad_minus_( cen_dist7_pad_  - 56.25 ),
	cen_dist10_pad_minus_( cen_dist10_pad_ - 100.0 ),
	cen_dist12_pad_minus_( cen_dist12_pad_ - 144.0 ),

	cen_dist5_pad_hinv_( 0.5 / cen_dist5_pad_ ),
	cen_dist6_pad_hinv_( 0.5 / cen_dist6_pad_ ),
	cen_dist7_pad_hinv_( 0.5 / cen_dist7_pad_ ),
	cen_dist10_pad_hinv_( 0.5 / cen_dist10_pad_ ),
	cen_dist12_pad_hinv_( 0.5 / cen_dist12_pad_ )
{
	// Load mp env info
	load_menv_info();
}

/// @brief Destructor
MembraneData::~MembraneData() {}

/// @brief Finalize Setup of MP Base Potential Class
void
MembraneData::finalize( pose::Pose & pose ) const {
	CenListInfo & cenlist( nonconst_cenlist_from_pose( pose ));
	cenlist.calculated() = false;
}

/// @brief Access Cenlist from Pose
/// @details Pose must already contain a cenlist object or this method will fail.
CenListInfo const &
MembraneData::get_cenlist_from_pose( pose::Pose const & pose ) const
{
	using namespace core::pose::datacache;
	return *( utility::pointer::static_pointer_cast< core::scoring::CenListInfo const > ( pose.data().get_const_ptr( CacheableDataType::CEN_LIST_INFO ) ));

}

/// @brief Database IO Helper Methods for Membrane
void
MembraneData::load_menv_info() {

	// Initialize Database Info
	Size const max_aa( 20 ); // restrict to canonical AAs for now
	Size const env_log_table_cen6_bins( 15 );
	Size const env_log_table_cen10_bins( 40 );
	Size const pair_log_table_size( 5 );
	Size const cbeta_den_table_size( 45 );
	Size const max_mem_layers( 3 );
	Size const min_mem_layers( 2 );

	std::string tag,line;
	chemical::AA aa;

	{ // mem_env_cen6:
		mem_env_log6_.dimension( max_aa, max_mem_layers,env_log_table_cen6_bins );
		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/MembranePotential/CEN6_mem_env_log.txt" );
		for ( Size i = 1; i <= max_aa; ++i ) {
			getline( stream, line );
			std::istringstream l( line );
			l >> tag >> aa;
			if ( l.fail() || tag != "MEM_ENV_LOG_CEN6:"  ) {
				utility_exit_with_message( "bad format for scoring/score_functions/MembranePotential/CEN6_mem_env_log.txt (cen6)" );
			}
			for ( Size j = 1; j <= max_mem_layers; ++j ) {
				for ( Size k = 1; k <= env_log_table_cen6_bins; ++k ) {
					l >> mem_env_log6_( aa, j, k );
				}
			}
		}
	}
	{ // mem_env_cen10:
		mem_env_log10_.dimension( max_aa, max_mem_layers,env_log_table_cen10_bins );

		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/MembranePotential/CEN10_mem_env_log.txt" );
		for ( Size i = 1; i <= max_aa; ++i ) {
			getline( stream, line );
			std::istringstream l( line );
			l >> tag >> aa;
			if ( l.fail() || tag != "MEM_ENV_LOG_CEN10:"  ) {
				utility_exit_with_message( "bad format for scoring/score_functions/MembranePotential/CEN10_mem_env_log.txt (cen10)" );
			}
			for ( Size j = 1; j <= max_mem_layers; ++j ) {
				for ( Size k = 1; k <= env_log_table_cen10_bins; ++k ) {
					l >> mem_env_log10_( aa, j, k );
				}
			}
		}
	}

	{ // cbeta_den_6/12
		mem_cbeta_den6_.dimension( cbeta_den_table_size );
		mem_cbeta_den12_.dimension( cbeta_den_table_size );
		mem_cbeta_2TM_den6_.dimension( cbeta_den_table_size );
		mem_cbeta_2TM_den12_.dimension( cbeta_den_table_size );
		mem_cbeta_4TM_den6_.dimension( cbeta_den_table_size );
		mem_cbeta_4TM_den12_.dimension( cbeta_den_table_size );


		// Read in etables into private vars
		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/MembranePotential/memcbeta_den.txt" );

		{ // den6
			getline( stream, line );
			{
				std::istringstream l(line);
				l >> tag;
				for ( Size i = 1; i <= cbeta_den_table_size; ++i ) {
					l >>mem_cbeta_den6_( i );
				}
				if ( l.fail() || tag != "MEMCBETA_DEN6:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (DEN6)");
			}
			getline( stream, line );
			{
				std::istringstream l( line );
				l >> tag;
				for ( Size i = 1; i<= cbeta_den_table_size; ++i ) {
					l >> mem_cbeta_2TM_den6_( i );
				}
				if ( l.fail() || tag != "MEMCBETA_2TM_DEN6:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (2TM_DEN6)");
			}
			getline( stream, line );
			{
				std::istringstream l( line );
				l >> tag;
				for ( Size i = 1; i <= cbeta_den_table_size; ++i ) {
					l >> mem_cbeta_4TM_den6_( i );
				}
				if ( l.fail() || tag != "MEMCBETA_4TM_DEN6:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (4TM_DEN6)");
			}
		} // end den6

		{ // den12
			{
				getline( stream, line );
				std::istringstream l( line );
				l >> tag;
				for ( Size i = 1; i <= cbeta_den_table_size; ++i ) {
					l >> mem_cbeta_den12_( i );
				}
				if ( l.fail() || tag != "MEMCBETA_DEN12:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (DEN12)");
			}
			getline( stream, line );
			{
				std::istringstream l(line);
				l >> tag;
				for ( Size i = 1; i <= cbeta_den_table_size; ++i ) {
					l >> mem_cbeta_2TM_den12_( i );
				}
				if ( l.fail() || tag != "MEMCBETA_2TM_DEN12:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (2TM_DEN12)");
			}
			getline( stream, line );
			{
				std::istringstream l(line);
				l >> tag;
				for ( Size i = 1; i <= cbeta_den_table_size; ++i ) {
					l >> mem_cbeta_4TM_den12_( i );
				}
				if ( l.fail() || tag != "MEMCBETA_4TM_DEN12:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (4TM_DEN12)");
			}

		}
	} // end den 12

	{ // pair_log
		mem_pair_log_.dimension( min_mem_layers, pair_log_table_size, max_aa, max_aa );

		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/MembranePotential/mem_pair_log.txt" );
		for ( Size i = 1; i <= min_mem_layers; ++i ) {
			for ( Size j = 1; j <= pair_log_table_size; ++j ) {
				for ( Size k = 1; k <= max_aa; ++k ) {
					getline( stream, line );
					std::istringstream l( line );
					Size ii,jj;
					l >> tag >> ii >> jj >> aa;
					debug_assert( Size( aa ) == k );
					for ( Size n=1; n <= max_aa; ++n ) {
						l >> mem_pair_log_(i, j, aa, n );
					}
					if ( l.fail() || ii != i || jj != j || tag != "MEM_PAIR_LOG:"  ) utility_exit_with_message("scoring/score_functions/MembranePotential/mem_pair_log.txt");
				}

			}
		}
	} // end pair log
}

/// @brief Membrane Base Potential Statistics
ObjexxFCL::FArray3D< Real > MembraneData::mem_env_log6() const { return mem_env_log6_; }
ObjexxFCL::FArray3D< Real > MembraneData::mem_env_log10() const { return mem_env_log10_; }

ObjexxFCL::FArray1D< Real > MembraneData::mem_cbeta_den6() const { return mem_cbeta_den6_; }
ObjexxFCL::FArray1D< Real > MembraneData::mem_cbeta_den12() const { return mem_cbeta_den12_; }
ObjexxFCL::FArray1D< Real > MembraneData::mem_cbeta_2TM_den6() const { return mem_cbeta_2TM_den6_; }
ObjexxFCL::FArray1D< Real > MembraneData::mem_cbeta_2TM_den12() const { return mem_cbeta_2TM_den12_; }
ObjexxFCL::FArray1D< Real > MembraneData::mem_cbeta_4TM_den6() const { return mem_cbeta_4TM_den6_; }
ObjexxFCL::FArray1D< Real > MembraneData::mem_cbeta_4TM_den12() const { return mem_cbeta_4TM_den12_; }

ObjexxFCL::FArray4D< Real > MembraneData::mem_pair_log() const { return mem_pair_log_; }

/// @brief Env Pair Potential Statistics
Real MembraneData::cen_dist5_pad() const { return cen_dist5_pad_; }
Real MembraneData::cen_dist6_pad() const { return cen_dist6_pad_; }
Real MembraneData::cen_dist7_pad() const { return cen_dist7_pad_; }
Real MembraneData::cen_dist10_pad() const { return cen_dist10_pad_; }
Real MembraneData::cen_dist12_pad() const { return cen_dist12_pad_; }

Real MembraneData::cen_dist5_pad_plus() const { return cen_dist5_pad_plus_; }
Real MembraneData::cen_dist6_pad_plus() const { return cen_dist6_pad_plus_; }
Real MembraneData::cen_dist7_pad_plus() const { return cen_dist7_pad_plus_; }
Real MembraneData::cen_dist10_pad_plus() const { return cen_dist10_pad_plus_; }
Real MembraneData::cen_dist12_pad_plus() const { return cen_dist12_pad_plus_; }

Real MembraneData::cen_dist5_pad_minus() const { return cen_dist5_pad_minus_; }
Real MembraneData::cen_dist7_pad_minus() const { return cen_dist7_pad_minus_; }
Real MembraneData::cen_dist10_pad_minus() const { return cen_dist10_pad_minus_; }
Real MembraneData::cen_dist12_pad_minus() const { return cen_dist12_pad_minus_; }

Real MembraneData::cen_dist5_pad_hinv() const { return cen_dist5_pad_hinv_; }
Real MembraneData::cen_dist6_pad_hinv() const { return cen_dist6_pad_hinv_; }
Real MembraneData::cen_dist7_pad_hinv() const { return cen_dist7_pad_hinv_; }
Real MembraneData::cen_dist10_pad_hinv() const { return cen_dist10_pad_hinv_; }
Real MembraneData::cen_dist12_pad_hinv() const { return cen_dist12_pad_hinv_; }

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MembraneData_cc

