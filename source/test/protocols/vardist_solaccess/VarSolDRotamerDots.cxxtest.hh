// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/vardist_solaccess/VarSolDRotamerDots.cxxtest.hh
/// @brief  test suite for Variable-distance hydrogen-bonds to solvent calculator
/// @author Kevin Houlihan (khouli@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/vardist_solaccess/VarSolDRotamerDots.hh>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Core headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/interaction_graph/RotamerDots.hh>


// Basic headers
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>

#include <cmath>


static basic::Tracer TR("VarSolDRotamerDots.cxxtest");

using namespace core;
using namespace core::chemical;
using namespace core::conformation;

class VarSolDRotamerDotsTest : public CxxTest::TestSuite
{
public:
	typedef core::pose::PoseOP PoseOP;
	typedef protocols::vardist_solaccess::VarSolDRotamerDots VarSolDRotamerDots;

public:
	void setUp() {
		core_init();
	}

	void test_self_overlap() {

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueOP res10( new core::conformation::Residue( trpcage.residue( 10 ) ) );
		protocols::vardist_solaccess::VarSolDistSasaCalculatorOP vsasa_calc( new protocols::vardist_solaccess::VarSolDistSasaCalculator );
		//protocols::vardist_solaccess::VarSolDistSasaCalculator vsasa_calc;
		VarSolDRotamerDots dots(res10, vsasa_calc);
		dots.increment_self_overlap();
		//ResidueKinWriter writer;
		//std::ofstream fout( "test.kin" );
		//writer.write_kin_header( fout, *res10 );
		//writer.write_rsd_coords( fout, *res10 );
		//dots.write_dot_kinemage( fout );

		core::Size res10natoms = res10->natoms();
		utility::vector1< core::Real > collision_radii( res10natoms );
		utility::vector1< core::Vector > res10coords( res10natoms );
		for ( core::Size ii = 1; ii <= res10natoms ; ++ii ) {
			collision_radii[ ii ] = dots.get_atom_collision_radius(ii);
			res10coords[ ii ] = res10->xyz( ii );
		}

		core::Real unscaled_fudge = 0;
		for ( core::Size ii = 1; ii <= res10natoms ; ++ii ) {
			for ( core::Size jj = 1; jj <= dots.nshells_for_atom(ii); ++jj ) {
				core::Real jjradius = dots.shell_radius_for_atom(ii,jj);
				core::Real fudge = jjradius + unscaled_fudge;
				for ( core::Size kk = 1; kk <= dots.ndots(); ++kk ) {
					core::Vector kkcoord = dots.get_atom_coords_xyz(ii) + jjradius * core::pack::interaction_graph::RotamerDots::dot_coord( kk );
					bool dot_covered = dots.get_dot_covered(ii,jj,kk);
					bool actually_covered = false;
					core::Real max_penetration_depth = -1.0;
					core::Real closest_contact = -1.0;
					for ( core::Size ll = 1; ll <= res10natoms; ++ll ) {
						if ( ll == ii ) continue; // don't consider self coverage
						core::Real ll_kk_d = kkcoord.distance( res10coords[ ll ] );
						if ( ll_kk_d < collision_radii[ ll ] ) {
							actually_covered = true;
							core::Real penetration_depth = collision_radii[ ll ] - ll_kk_d;
							if ( max_penetration_depth < 0 || penetration_depth > max_penetration_depth ) {
								max_penetration_depth = penetration_depth;
							}
						} else if ( closest_contact < 0 || ll_kk_d - collision_radii[ ll ] < closest_contact ) {
							closest_contact = ll_kk_d - collision_radii[ ll ];
						}
					}
					TS_ASSERT( dot_covered == actually_covered || ( max_penetration_depth > 0 && max_penetration_depth < fudge ) || ( closest_contact < fudge ) );
				}
			}
		}
	}

	void test_residue_pair_overlap() {

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueOP res10( new core::conformation::Residue( trpcage.residue( 10 ) ) );
		ResidueOP res11( new core::conformation::Residue( trpcage.residue( 11 ) ) );
		protocols::vardist_solaccess::VarSolDistSasaCalculatorOP vsasa_calc( new protocols::vardist_solaccess::VarSolDistSasaCalculator );
		//protocols::vardist_solaccess::VarSolDistSasaCalculator vsasa_calc;
		VarSolDRotamerDots dots10(res10, vsasa_calc);
		VarSolDRotamerDots dots11(res11, vsasa_calc);
		//dots10.increment_self_overlap();
		//dots11.increment_self_overlap();
		dots10.intersect_residues( dots11 );

		core::Size res10natoms = res10->natoms();
		core::Size res11natoms = res11->natoms();

		utility::vector1< core::Real > r10_collision_radii( res10natoms );
		utility::vector1< core::Real > r11_collision_radii( res11natoms );

		utility::vector1< core::Vector > res10coords( res10natoms );
		utility::vector1< core::Vector > res11coords( res11natoms );

		for ( core::Size ii = 1; ii <= res10natoms ; ++ii ) {
			r10_collision_radii[ ii ] = dots10.get_atom_collision_radius(ii);
			res10coords[ ii ] = res10->xyz( ii );
		}
		for ( core::Size ii = 1; ii <= res11natoms ; ++ii ) {
			r11_collision_radii[ ii ] = dots11.get_atom_collision_radius(ii);
			res11coords[ ii ] = res11->xyz( ii );
		}

		core::Real unscaled_fudge = 0;

		for ( core::Size ii = 1; ii <= res10natoms ; ++ii ) {
			for ( core::Size jj = 1; jj <= dots10.nshells_for_atom(ii); ++jj ) {
				core::Real jjradius = dots10.shell_radius_for_atom(ii,jj);
				core::Real fudge = jjradius + unscaled_fudge;
				for ( core::Size kk = 1; kk <= dots10.ndots(); ++kk ) {
					core::Vector kkcoord = dots10.get_atom_coords_xyz(ii) + jjradius * core::pack::interaction_graph::RotamerDots::dot_coord( kk );
					bool dot_covered = dots10.get_dot_covered(ii,jj,kk);
					bool actually_covered = false;
					core::Real max_penetration_depth = -1.0;
					core::Real closest_contact = -1.0;
					for ( core::Size ll = 1; ll <= res11natoms; ++ll ) {
						if ( ll == ii ) continue; // don't consider self coverage
						core::Real ll_kk_d = kkcoord.distance( res11coords[ ll ] );
						if ( ll_kk_d < r11_collision_radii[ ll ] ) {
							actually_covered = true;
							core::Real penetration_depth = r11_collision_radii[ ll ] - ll_kk_d;
							if ( max_penetration_depth < 0 || penetration_depth > max_penetration_depth ) {
								max_penetration_depth = penetration_depth;
							}
						} else if ( closest_contact < 0 || ll_kk_d - r11_collision_radii[ ll ] < closest_contact ) {
							closest_contact = ll_kk_d - r11_collision_radii[ ll ];
						}
					}
					TS_ASSERT( dot_covered == actually_covered || ( max_penetration_depth > 0 && max_penetration_depth < fudge ) || ( closest_contact < fudge ) );
				}
			}
		}


		for ( core::Size ii = 1; ii <= res11natoms ; ++ii ) {
			for ( core::Size jj = 1; jj <= dots11.nshells_for_atom(ii); ++jj ) {
				core::Real jjradius = dots11.shell_radius_for_atom(ii,jj);
				core::Real fudge = jjradius + unscaled_fudge;
				for ( core::Size kk = 1; kk <= dots11.ndots(); ++kk ) {
					core::Vector kkcoord = dots11.get_atom_coords_xyz(ii) + jjradius * core::pack::interaction_graph::RotamerDots::dot_coord( kk );
					bool dot_covered = dots11.get_dot_covered(ii,jj,kk);
					bool actually_covered = false;
					core::Real max_penetration_depth = -1.0;
					core::Real closest_contact = -1.0;
					for ( core::Size ll = 1; ll <= res10natoms; ++ll ) {
						if ( ll == ii ) continue; // don't consider self coverage
						core::Real ll_kk_d = kkcoord.distance( res10coords[ ll ] );
						if ( ll_kk_d < r10_collision_radii[ ll ] ) {
							actually_covered = true;
							core::Real penetration_depth = r10_collision_radii[ ll ] - ll_kk_d;
							if ( max_penetration_depth < 0 || penetration_depth > max_penetration_depth ) {
								max_penetration_depth = penetration_depth;
							}
						} else if ( closest_contact < 0 || ll_kk_d - r10_collision_radii[ ll ] < closest_contact ) {
							closest_contact = ll_kk_d - r10_collision_radii[ ll ];
						}
					}
					TS_ASSERT( dot_covered == actually_covered || ( max_penetration_depth > 0 && max_penetration_depth < fudge ) || ( closest_contact < fudge ) );
				}
			}
		}


		//ResidueKinWriter writer;
		//std::ofstream fout( "test2.kin" );
		//writer.write_kin_header( fout, *res10 );
		//writer.write_rsd_coords( fout, *res10 );
		//writer.write_rsd_coords( fout, *res11 );
		//dots10.write_dot_kinemage( fout );
		//dots11.write_dot_kinemage( fout );

	}

	void radii_monotonically_increasing() {
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueOP res1( new core::conformation::Residue( trpcage.residue( 1 ) ) );
		protocols::vardist_solaccess::VarSolDistSasaCalculatorOP vsasa_calc( new protocols::vardist_solaccess::VarSolDistSasaCalculator );
		VarSolDRotamerDots dots1(res1, vsasa_calc);
		core::Size res1natoms = res1->natoms();
		utility::vector1< core::Real > r1_collision_radii( res1natoms );
		utility::vector1< core::Vector > res1coords( res1natoms );

		for ( core::Size ii = 1; ii <= res1natoms ; ++ii ) {
			if ( dots1.nshells_for_atom(ii) < 2 ) continue;
			for ( core::Size jj = 1; jj <= dots1.nshells_for_atom(ii) - 1; ++jj ) {
				core::Real jjradius = dots1.shell_radius_for_atom(ii,jj);
				core::Real kkradius = dots1.shell_radius_for_atom(ii, jj + 1);
				TS_ASSERT( (kkradius - jjradius) > 0);
			}
		}
	}

	// void test_wobble_decreases_collision_radius() {
	//  core::pose::Pose trpcage = create_trpcage_ideal_pose();
	//  ResidueOP res10 = new core::conformation::Residue( trpcage.residue( 10 ) );
	//  ResidueOP res11 = new core::conformation::Residue( trpcage.residue( 11 ) );
	//  core::Real probe_radius(1.4);
	//  core::Real wobble_0p0(0.0);
	//  core::Real wobble_0p3(0.3);
	//  VarSolDRotamerDots dots10_wb0p0( res10, true, probe_radius, wobble_0p0 );
	//  VarSolDRotamerDots dots10_wb0p3( res10, true, probe_radius, wobble_0p3 );
	//
	//  core::Size res10natoms = res10->natoms();
	//
	//  TR << "wb0p0 wobble = " << dots10_wb0p0.get_wobble() << std::endl;
	//  TR << "wb0p3 wobble = " << dots10_wb0p3.get_wobble() << std::endl;
	//  for ( core::Size ii = 1; ii <= res10natoms ; ++ii ) {
	//   core::Real wb0p0_c_radius = dots10_wb0p0.get_atom_collision_radius(ii);
	//   core::Real wb0p3_c_radius = dots10_wb0p3.get_atom_collision_radius(ii);
	//   TR << "wb0p0 collision radius = " << wb0p0_c_radius << std::endl;
	//   TR << "wb0p3 collision radius = " << wb0p3_c_radius << std::endl;
	//   TS_ASSERT(
	//     wb0p3_c_radius < wb0p0_c_radius ||
	//     wb0p0_c_radius == 0.0
	//     );
	//  }
	// }
	//
	// void test_wobble_increases_hbond_radius() {
	//  core::pose::Pose trpcage = create_trpcage_ideal_pose();
	//  ResidueOP res10 = new core::conformation::Residue( trpcage.residue( 10 ) );
	//  core::Real probe_radius(1.4);
	//  core::Real wobble_0p0(0.0);
	//  core::Real wobble_0p3(0.3);
	//  VarSolDRotamerDots dots10_wb0p0( res10, true, probe_radius, wobble_0p0 );
	//  VarSolDRotamerDots dots10_wb0p3( res10, true, probe_radius, wobble_0p3 );
	//
	//  core::Size res10natoms = res10->natoms();
	//
	//  for ( core::Size ii = 1; ii <= res10natoms ; ++ii ) {
	//   core::Real wb0p0_i_radius = dots10_wb0p0.get_atom_interaction_radius(ii);
	//   core::Real wb0p3_i_radius = dots10_wb0p3.get_atom_interaction_radius(ii);
	//   TR << "wb0p0 interaction radius = " << wb0p0_i_radius << std::endl;
	//   TR << "wb0p3 interaction radius = " << wb0p3_i_radius << std::endl;
	//   TS_ASSERT(
	//     wb0p3_i_radius > wb0p0_i_radius ||
	//     wb0p3_i_radius == 0.0
	//     );
	//  }
	// }

	// void test_probe_radius_increases_collision_radius() {
	//  core::pose::Pose trpcage = create_trpcage_ideal_pose();
	//  ResidueOP res10 = new core::conformation::Residue( trpcage.residue( 10 ) );
	//  core::Real probe_radius_pr1p4(1.4);
	//  core::Real probe_radius_pr1p2(1.2);
	//  core::Real wobble_0p0(0.0);
	//  VarSolDRotamerDots dots10_pr1p4( res10, true, probe_radius_pr1p4, wobble_0p0 );
	//  VarSolDRotamerDots dots10_pr1p2( res10, true, probe_radius_pr1p2, wobble_0p0 );
	//
	//  core::Size res10natoms = res10->natoms();
	//
	//  for ( core::Size ii = 1; ii <= res10natoms ; ++ii ) {
	//   TS_ASSERT(
	//     dots10_pr1p4.get_atom_collision_radius(ii) > dots10_pr1p2.get_atom_collision_radius(ii)
	//     );
	//   TS_ASSERT(
	//     dots10_pr1p2.get_atom_collision_radius(ii) > 0
	//     );
	//  }
	// }

	void test_shells_evenly_spaced() {
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		ResidueOP res1( new core::conformation::Residue( trpcage.residue( 1 ) ) );
		protocols::vardist_solaccess::VarSolDistSasaCalculatorOP vsasa_calc( new protocols::vardist_solaccess::VarSolDistSasaCalculator );
		//protocols::vardist_solaccess::VarSolDistSasaCalculator vsasa_calc;
		VarSolDRotamerDots dots1(res1, vsasa_calc);
		core::Size res1natoms = res1->natoms();

		core::Real fudge = 1e-5;
		for ( core::Size ii = 1; ii <= res1natoms ; ++ii ) {
			if ( dots1.nshells_for_atom(ii) < 3 ) continue;
			for ( core::Size jj = 1; jj <= dots1.nshells_for_atom(ii) - 2; ++jj ) {
				core::Real jjradius = dots1.shell_radius_for_atom(ii,jj);
				core::Real kkradius = dots1.shell_radius_for_atom(ii, jj + 1);
				core::Real llradius = dots1.shell_radius_for_atom(ii, jj + 2);
				TS_ASSERT_DELTA( (llradius - kkradius), (kkradius - jjradius), fudge );
			}
		}

	}

	void test_radii_setting() {
		TS_ASSERT( true );
		// It's not clear where I was going with this unit test.
		//
		// core::pose::Pose trpcage = create_trpcage_ideal_pose();
		// ResidueOP res1( new core::conformation::Residue( trpcage.residue( 1 ) ) );
		// protocols::vardist_solaccess::VarSolDistSasaCalculatorOP vsasa_calc( new protocols::vardist_solaccess::VarSolDistSasaCalculator );
		//protocols::vardist_solaccess::VarSolDistSasaCalculator vsasa_calc1;
		//protocols::vardist_solaccess::VarSolDistSasaCalculator vsasa_calc2;
	}

};
