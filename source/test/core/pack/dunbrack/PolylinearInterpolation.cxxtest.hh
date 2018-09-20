// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/dunbrack/PolylinearInterpolation.cxxtest.hh
/// @brief  Test the polylinear interpolation function used by the Dunbrack libraries.
/// @author Andy Watkins
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/pack/dunbrack/DunbrackEnergy.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/optimization/atom_tree_minimize.hh>
#include <core/optimization/AtomTreeMultifunc.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/types.hh>

#include <utility/fixedsizearray1.hh>

// Protocols Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh> //For convenience to build pose

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <numeric/conversions.hh>
#include <numeric/angle.functions.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("PolylinearInterpolationTests");

#define DERIV_DELTA 1e-9
#define DERIV_DELTA_LARGE 1e-6

//Uncomment the following line to do a rigourous loop though all the combinations of dihedrals for lysine (takes several minutes):
//#define LYSINE_DERIV_TEST_FULL

class PolylinearInterpolationTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-out:level 400" );
	}

	void tearDown(){
	}


	void test_polylinear_interpolation_lysine_derivs() {
		core::pose::Pose pose;
		{
			protocols::cyclic_peptide::PeptideStubMover make_pep;
			make_pep.set_reset_mode(true);
			make_pep.add_residue( "Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "" );
			make_pep.add_residue( "Append", "LYS", 2, false, "", 0, 1, "" );
			make_pep.add_residue( "Append", "GLY:CtermProteinFull", 3, false, "", 0, 2, "" );
			make_pep.apply(pose);
		}

		core::scoring::ScoreFunctionOP sfxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
		sfxn->set_weight(core::scoring::fa_dun, 1.0);
		(*sfxn)(pose);

		utility::vector1< core::Real > angle_vals;
		angle_vals.push_back( 0.823 );
		angle_vals.push_back( 36.135 );
		angle_vals.push_back( 179.541 );
		angle_vals.push_back( -178.442 );
		angle_vals.push_back( -63.51 );
		angle_vals.push_back( -0.31 );

		core::kinematics::MoveMap movemap;
		movemap.set_bb(true);
		movemap.set_chi(true);
		core::optimization::MinimizerMap minmap;
		minmap.setup(pose, movemap);

		pose.energies().set_use_nblist(pose, minmap.domain_map(), false);
		sfxn->setup_for_minimizing(pose, minmap);

		core::optimization::Multivec start_vars( minmap.nangles() );
		core::optimization::Multivec dE_dvars( minmap.nangles() );

		for ( core::Size phi(1); phi<=angle_vals.size(); ++phi ) {
			pose.set_phi(2, angle_vals[phi]);
			for ( core::Size psi(1); psi<=angle_vals.size(); ++psi ) {
				pose.set_psi(2, angle_vals[psi]);
				for ( core::Size chi1(1); chi1 <= angle_vals.size(); ++chi1 ) {
					pose.set_chi( 1, 2, angle_vals[chi1] );
#ifdef LYSINE_DERIV_TEST_FULL
					for( core::Size chi2(1); chi2 <= angle_vals.size(); ++chi2 ) {
#else
					{
						core::Size const chi2(chi1);
#endif
						pose.set_chi( 2, 2, angle_vals[chi2] );

#ifdef LYSINE_DERIV_TEST_FULL
						for( core::Size chi3(1); chi3 <= angle_vals.size(); ++chi3 ) {
#else
						{
							core::Size const chi3(chi2);
#endif
							pose.set_chi( 3, 2, angle_vals[chi3] );

#ifdef LYSINE_DERIV_TEST_FULL
							for( core::Size chi4(1); chi4 <= angle_vals.size(); ++chi4 ) {
#else
							{
								core::Size const chi4( chi3 );
#endif
								pose.set_chi( 4, 2, angle_vals[chi4] );
								pose.update_residue_neighbors();

								sfxn->setup_for_derivatives(pose);
								core::optimization::AtomTreeMultifunc func( pose, minmap, *sfxn, false, false );
								minmap.copy_dofs_from_pose( pose, start_vars );
								func.dfunc( start_vars, dE_dvars );
								core::optimization::SimpleDerivCheckResult result( core::optimization::simple_numeric_deriv_check( func, start_vars, dE_dvars, false, false, 5, 1.0e-5) );
								bool failure(false);
								for ( core::Size ii(1), iimax( result.nangles()); ii<=iimax; ++ii ) {
									TS_ASSERT_DELTA( result.step_data( ii, 1 ).ana_deriv(), result.step_data( ii, 1 ).num_deriv(), 1e-4 );
									if ( std::abs( result.step_data( ii, 1 ).ana_deriv() - result.step_data( ii, 1 ).num_deriv() ) > 1e-4 ) {
										failure=true;
										break;
									}
								}
								if ( failure ) {
									TR << "00000 FAILURE!\tphi: " << pose.phi(2) << "\tpsi: " << pose.psi(2) << "\tchi1: " << pose.chi(1,2) << "\tchi2: " << pose.chi(2,2) << "\tchi3: " << pose.chi(3,2) << "\tchi4: " << pose.chi(4,2) << std::endl;
								} else {
									TR << "11111 SUCCESS!\tphi: " << pose.phi(2) << "\tpsi: " << pose.psi(2) << "\tchi1: " << pose.chi(1,2) << "\tchi2: " << pose.chi(2,2) << "\tchi3: " << pose.chi(3,2) << "\tchi4: " << pose.chi(4,2) << std::endl;
								}
							}
						}
					}
				}
			}
		}
	}



	void test_1D_wraparound() {
		using namespace core::pack::dunbrack;

		utility::fixedsizearray1< double, 2 > vals;
		vals[1] = -170; vals[2] = 170;
		utility::fixedsizearray1< double, 1 > binrange;
		binrange[1] = 20;
		double val_out, val_out_plus_dx;
		utility::fixedsizearray1< double, 1 > dval_dbb_out, dval_dbb_out_plus_dx;

		utility::fixedsizearray1< double, 1 > bbd, bbd_plus_dx, bbd_plus_dy;
		for ( int ii=0; ii<=100; ii+=5 ) {
			bbd[1] = double(ii)/100;
			bbd_plus_dx = bbd;  bbd_plus_dx[1] += DERIV_DELTA_LARGE;
			interpolate_polylinear_by_value( vals, bbd, binrange, true, val_out, dval_dbb_out );
			TR << ii << " " << vals << " " << bbd << " " << val_out << " " << dval_dbb_out << std::endl;
			core::Real const num( numeric::principal_angle_degrees( 170 * bbd[1] + 190 * (1.0 - bbd[1] ) ) );
			TS_ASSERT_DELTA( num, val_out, 1e-3);

			interpolate_polylinear_by_value( vals, bbd_plus_dx, binrange, true, val_out_plus_dx, dval_dbb_out_plus_dx );
			TR << ii << " " << vals << " " << bbd_plus_dx << " " << val_out_plus_dx << " " << dval_dbb_out_plus_dx << std::endl;

			TS_ASSERT_DELTA( dval_dbb_out[1]*binrange[1], (val_out_plus_dx - val_out)/DERIV_DELTA_LARGE, 1e-3 );
		}
	}

	void test_1D() {
		using namespace core::pack::dunbrack;

		utility::fixedsizearray1< double, 2 > vals;
		vals[1] = 45; vals[2] = 60;
		utility::fixedsizearray1< double, 1 > binrange;
		binrange[1] = 15;
		double val_out, val_out_plus_dx;
		utility::fixedsizearray1< double, 1 > dval_dbb_out, dval_dbb_out_plus_dx;

		utility::fixedsizearray1< double, 1 > bbd, bbd_plus_dx, bbd_plus_dy;
		for ( int ii=0; ii<=100; ii+=5 ) {
			bbd[1] = double(ii)/100;
			bbd_plus_dx = bbd;  bbd_plus_dx[1] += DERIV_DELTA_LARGE;
			interpolate_polylinear_by_value( vals, bbd, binrange, true, val_out, dval_dbb_out );
			TR << ii << " " << vals << " " << bbd << " " << val_out << " " << dval_dbb_out << std::endl;
			core::Real const num( numeric::principal_angle_degrees( 60 * bbd[1] + 45 * (1.0 - bbd[1] ) ) );
			TS_ASSERT_DELTA( num, val_out, 1e-3);

			interpolate_polylinear_by_value( vals, bbd_plus_dx, binrange, true, val_out_plus_dx, dval_dbb_out_plus_dx );
			TR << ii << " " << vals << " " << bbd_plus_dx << " " << val_out_plus_dx << " " << dval_dbb_out_plus_dx << std::endl;

			TS_ASSERT_DELTA( dval_dbb_out[1]*binrange[1], (val_out_plus_dx - val_out)/DERIV_DELTA_LARGE, 1e-3 );
		}
	}

	/// @brief Test the derivatives with polylinear interpolation.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_polylinear_interpolation_2D_derivs() {
		using namespace core::pack::dunbrack;

		utility::fixedsizearray1< double, 4 > vals;
		vals[1] = 1; vals[2] = 2; vals[3] = 10; vals[4] = 12;
		utility::fixedsizearray1< double, 2 > binrange;
		binrange[1] = 10; binrange[2] = 10;
		double val_out, val_out_plus_dx, val_out_plus_dy;
		utility::fixedsizearray1< double, 2 > dval_dbb_out, dval_dbb_out_plus_dx, dval_dbb_out_plus_dy;

		utility::fixedsizearray1< double, 2 > bbd, bbd_plus_dx, bbd_plus_dy;
		for ( int jj=0; jj<=100; jj+=5 ) {
			for ( int ii=0; ii<=100; ii+=5 ) {
				bbd[1] = double(ii)/100.0; bbd[2] = double(jj)/100.0;
				bbd_plus_dx = bbd;  bbd_plus_dx[1] += DERIV_DELTA;
				bbd_plus_dy = bbd;  bbd_plus_dy[2] += DERIV_DELTA;
				interpolate_polylinear_by_value( vals, bbd, binrange, true, val_out, dval_dbb_out );
				TR << ii << " " << jj << " " << val_out << " " << dval_dbb_out << std::endl;
				core::Real const a( 10*bbd[1] + 1*(1.0-bbd[1])  );
				core::Real const b( 12*bbd[1] + 2*(1.0-bbd[1])  );
				core::Real const num( b*bbd[2] + a*(1.0-bbd[2]) );
				TS_ASSERT_DELTA( num, val_out, 1e-3);

				interpolate_polylinear_by_value( vals, bbd_plus_dx, binrange, true, val_out_plus_dx, dval_dbb_out_plus_dx );
				interpolate_polylinear_by_value( vals, bbd_plus_dy, binrange, true, val_out_plus_dy, dval_dbb_out_plus_dy );

				TS_ASSERT_DELTA( dval_dbb_out[1]*binrange[1], (val_out_plus_dx - val_out)/DERIV_DELTA, 1e-3 );
				TS_ASSERT_DELTA( dval_dbb_out[2]*binrange[2], (val_out_plus_dy - val_out)/DERIV_DELTA, 1e-3 );

			}
		}
	}

	/// @brief Test the derivatives with polylinear interpolation.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_polylinear_interpolation_3D_derivs() {
		using namespace core::pack::dunbrack;

		utility::fixedsizearray1< double, 8 > vals;
		vals[1] = 1; vals[2] = 2; vals[3] = 10; vals[4] = 12;
		vals[5] = 5; vals[6] = 18; vals[7] = 29; vals[8] = 22;
		utility::fixedsizearray1< double, 3 > binrange;
		binrange[1] = 10; binrange[2] = 14; binrange[3] = 17;
		double val_out, val_out_plus_dx, val_out_plus_dy, val_out_plus_dz;
		utility::fixedsizearray1< double, 3 > dval_dbb_out, dval_dbb_out_plus_dx, dval_dbb_out_plus_dy, dval_dbb_out_plus_dz;

		utility::fixedsizearray1< double, 3 > bbd, bbd_plus_dx, bbd_plus_dy, bbd_plus_dz;
		for ( int kk=0; kk<=100; kk+=5 ) {
			for ( int jj=0; jj<=100; jj+=5 ) {
				for ( int ii=0; ii<=100; ii+=5 ) {
					bbd[1] = double(ii)/100; bbd[2] = double(jj)/100; bbd[3] = double(kk)/100;
					bbd_plus_dx = bbd;  bbd_plus_dx[1] += DERIV_DELTA;
					bbd_plus_dy = bbd;  bbd_plus_dy[2] += DERIV_DELTA;
					bbd_plus_dz = bbd;  bbd_plus_dz[3] += DERIV_DELTA;
					interpolate_polylinear_by_value( vals, bbd, binrange, true, val_out, dval_dbb_out );
					TR << ii << " " << jj << " " << val_out << " " << dval_dbb_out << std::endl;

					core::Real const a( 1*(1.0-bbd[1]) + 5*bbd[1]);
					core::Real const b( 2*(1.0-bbd[1]) + 18*bbd[1]);
					core::Real const c( 10*(1.0-bbd[1]) + 29*bbd[1]);
					core::Real const d( 12*(1.0-bbd[1]) + 22*bbd[1]);

					core::Real const aa( a*(1.0-bbd[2]) + c*bbd[2]);
					core::Real const bb( b*(1.0-bbd[2]) + d*bbd[2]);

					core::Real const num( aa*(1.0-bbd[3]) + bb*bbd[3]);

					TS_ASSERT_DELTA( val_out, num, 1e-3 );

					interpolate_polylinear_by_value( vals, bbd_plus_dx, binrange, true, val_out_plus_dx, dval_dbb_out_plus_dx );
					interpolate_polylinear_by_value( vals, bbd_plus_dy, binrange, true, val_out_plus_dy, dval_dbb_out_plus_dy );
					interpolate_polylinear_by_value( vals, bbd_plus_dz, binrange, true, val_out_plus_dz, dval_dbb_out_plus_dz );

					TS_ASSERT_DELTA( dval_dbb_out[1]*binrange[1], (val_out_plus_dx - val_out)/DERIV_DELTA, 1e-3 );
					TS_ASSERT_DELTA( dval_dbb_out[2]*binrange[2], (val_out_plus_dy - val_out)/DERIV_DELTA, 1e-3 );
					TS_ASSERT_DELTA( dval_dbb_out[3]*binrange[3], (val_out_plus_dz - val_out)/DERIV_DELTA, 1e-3 );
				}
			}
		}
	}


	/// @brief Test the derivatives with polylinear interpolation near the -180/180 point.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void test_polylinear_interpolation_3D_derivs_wraparound() {
		using namespace core::pack::dunbrack;

		utility::fixedsizearray1< double, 8 > vals;
		vals[1] = -175; vals[2] = 175; vals[3] = 174; vals[4] = -172;
		vals[5] = -178; vals[6] = -179; vals[7] = 172; vals[8] = 173;
		utility::fixedsizearray1< double, 3 > binrange;
		binrange[1] = 10; binrange[2] = 14; binrange[3] = 17;
		double val_out, val_out_plus_dx, val_out_plus_dy, val_out_plus_dz;
		utility::fixedsizearray1< double, 3 > dval_dbb_out, dval_dbb_out_plus_dx, dval_dbb_out_plus_dy, dval_dbb_out_plus_dz;

		utility::fixedsizearray1< double, 3 > bbd, bbd_plus_dx, bbd_plus_dy, bbd_plus_dz;
		for ( int kk=0; kk<=100; kk+=5 ) {
			for ( int jj=0; jj<=100; jj+=5 ) {
				for ( int ii=0; ii<=100; ii+=5 ) {
					bbd[1] = double(ii)/100; bbd[2] = double(jj)/100; bbd[3] = double(kk)/100;
					bbd_plus_dx = bbd;  bbd_plus_dx[1] += DERIV_DELTA;
					bbd_plus_dy = bbd;  bbd_plus_dy[2] += DERIV_DELTA;
					bbd_plus_dz = bbd;  bbd_plus_dz[3] += DERIV_DELTA;
					interpolate_polylinear_by_value( vals, bbd, binrange, true, val_out, dval_dbb_out );
					TR << ii << " " << jj << " " << val_out << " " << dval_dbb_out << std::endl;

					core::Real const a( 185*(1.0-bbd[1]) + 182*bbd[1]);
					core::Real const b( 175*(1.0-bbd[1]) + 181*bbd[1]);
					core::Real const c( 174*(1.0-bbd[1]) + 172*bbd[1]);
					core::Real const d( 188*(1.0-bbd[1]) + 173*bbd[1]);

					core::Real const aa( a*(1.0-bbd[2]) + c*bbd[2]);
					core::Real const bb( b*(1.0-bbd[2]) + d*bbd[2]);

					core::Real const num( numeric::principal_angle_degrees(aa*(1.0-bbd[3]) + bb*bbd[3]));

					TS_ASSERT_DELTA( val_out, num, 1e-3 );

					interpolate_polylinear_by_value( vals, bbd_plus_dx, binrange, true, val_out_plus_dx, dval_dbb_out_plus_dx );
					interpolate_polylinear_by_value( vals, bbd_plus_dy, binrange, true, val_out_plus_dy, dval_dbb_out_plus_dy );
					interpolate_polylinear_by_value( vals, bbd_plus_dz, binrange, true, val_out_plus_dz, dval_dbb_out_plus_dz );

					TS_ASSERT_DELTA( dval_dbb_out[1]*binrange[1], (numeric::principal_angle_degrees(val_out_plus_dx - val_out))/DERIV_DELTA, 1e-3 );
					TS_ASSERT_DELTA( dval_dbb_out[2]*binrange[2], (numeric::principal_angle_degrees(val_out_plus_dy - val_out))/DERIV_DELTA, 1e-3 );
					TS_ASSERT_DELTA( dval_dbb_out[3]*binrange[3], (numeric::principal_angle_degrees(val_out_plus_dz - val_out))/DERIV_DELTA, 1e-3 );
				}
			}
		}
	}

};
