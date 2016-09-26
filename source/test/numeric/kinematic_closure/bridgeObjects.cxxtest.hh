// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   bridgeObjects.cxxtest.hh
/// @brief  Unit tests for the bridgeObjects kinematic closure function.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Headers {{{1

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <numeric/types.hh>
#include <numeric/kinematic_closure/dixon.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
#include <numeric/linear_algebra/GeneralizedEigenSolver.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/vector1.hh>
#include <boost/foreach.hpp>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>

static THREAD_LOCAL basic::Tracer TR("numeric.kinematic_closure.bridgeObjectsTests.cxxtest");

class bridgeObjectsTests : public CxxTest::TestSuite {
public:

	bridgeObjectsTests() {}
	
	void setUp() { 
		core_init();
		peptide_pose_ = core::import_pose::pose_from_file( "numeric/kinematic_closure/test_peptide.pdb", false, core::import_pose::PDB_file);
	}

	void tearDown() {}
	
	/// @brief Are two angles within a threshhold of one another?
	///
	bool within_thresh( core::Real const &val1, core::Real const &val2, core::Real const &thresh ) {
		if ( std::abs( val1 - val2 ) <= thresh ) return true;
		core::Real val1prime = val1;
		core::Real val2prime = val2;
		if ( val1<val2 ) {
			val1prime += 360.0;
		} else {
			val2prime += 360.0;
		}
		return ( std::abs( val1prime - val2prime ) <= thresh );
	}

	/// @brief Test that mirror-image inputs yield mirror-image solutions.
	///
	void test_symm_solutions() {
		TR << "Starting bridgeObjectsTests::test_symm_solutions." << std::endl;
	
		using namespace numeric::kinematic_closure;
		
		utility::vector1< utility::vector1< numeric::Real > > atoms1, atoms2;
		atoms1.reserve( peptide_pose_->total_residue() * 3 );
		atoms2.reserve( peptide_pose_->total_residue() * 3 );
		
		for(numeric::Size ir=1, irmax=peptide_pose_->total_residue(); ir<=irmax; ++ir) {
			numeric::xyzVector < numeric::Real > xyzN( peptide_pose_->xyz( core::id::NamedAtomID("N", ir) ) );
			numeric::xyzVector < numeric::Real > xyzCA( peptide_pose_->xyz( core::id::NamedAtomID("CA", ir) ) );
			numeric::xyzVector < numeric::Real > xyzC( peptide_pose_->xyz( core::id::NamedAtomID("C", ir) ) );
			utility::vector1 < numeric::Real > Nvec(3), CAvec(3), Cvec(3);
			Nvec[1] = xyzN.x(); Nvec[2] = xyzN.y(); Nvec[3] = xyzN.z();
			CAvec[1] = xyzCA.x(); CAvec[2] = xyzCA.y(); CAvec[3] = xyzCA.z();
			Cvec[1] = xyzC.x(); Cvec[2] = xyzC.y(); Cvec[3] = xyzC.z();
			atoms1.push_back( Nvec );
			atoms1.push_back( CAvec );
			atoms1.push_back( Cvec );
		}
		
		atoms2 = atoms1;
		for(numeric::Size i=1, imax=atoms2.size(); i<=imax; ++i) {
			atoms2[i][1] *= -1.0;
		}
		

		utility::vector1< numeric::Real > dt1, dt2, da1, da2, db1, db2;
		utility::vector1<utility::vector1<numeric::Real> > q1(3); //Used by numeric::kinematic_closure::chainTORS
		utility::vector1<numeric::Real> r1(3); //Used by numeric::kinematic_closure::chainTORS
		utility::vector1<utility::vector1<numeric::Real> > q2(3); //Used by numeric::kinematic_closure::chainTORS
		utility::vector1<numeric::Real> r2(3); //Used by numeric::kinematic_closure::chainTORS
		
		chainTORS( atoms1.size(), atoms1, dt1, da1, db1, r1, q1 );
		chainTORS( atoms2.size(), atoms2, dt2, da2, db2, r2, q2 );
		
		utility::vector1< numeric::Size > pivots;
		pivots.resize(3);
		pivots[1] = 5;
		pivots[2] = 15;
		pivots[3] = 20;
		utility::vector1 < numeric::Size > order;
		order.resize(3); order[1]=1; order[2]=2; order[3]=3;
		int nsol1, nsol2;
		utility::vector1<utility::vector1<numeric::Real> > t_ang1, t_ang2;
		utility::vector1<utility::vector1<numeric::Real> > b_ang1, b_ang2;
		utility::vector1<utility::vector1<numeric::Real> > b_len1, b_len2;
		
		bridgeObjects(atoms1, dt1, da1, db1, pivots, order, t_ang1, b_ang1, b_len1, nsol1);		
		bridgeObjects(atoms2, dt2, da2, db2, pivots, order, t_ang2, b_ang2, b_len2, nsol2);
		
		TS_ASSERT_EQUALS( nsol1, nsol2 );
		TR << "Found " << nsol1 << " solutions for left-handed, and " << nsol2 << " solutions for right-handed." << std::endl;
		
		for( numeric::Size i=1; i<=static_cast<numeric::Size>(nsol1); ++i) {
			TR << "SOLUTION " << i << " TORSIONS:" << std::endl;			
			TS_ASSERT_EQUALS( t_ang1[i].size(), t_ang2[i].size() );
			for(numeric::Size j=1, jmax=t_ang1[i].size(); j<=jmax; ++j) {
				TR << t_ang1[i][j] << "\t" << t_ang2[i][j] << std::endl;
				TS_ASSERT( within_thresh( t_ang1[i][j], -1.0*t_ang2[i][j], 0.01 ) );
			}
		}		
		
		//TODO!!!
		TR.flush();
	}

private:

	core::pose::PoseOP peptide_pose_;
	
};

