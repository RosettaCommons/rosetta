// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/modifications/ValenceHandler.cxxtest.hh
/// @brief  test suite for core/chemical/modifications/ValenceHandler.hh
/// @author Steven Combs


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>


#include <core/chemical/modifications/ValenceHandler.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>
#include <core/chemical/gasteiger/GasteigerAtomTyper.hh>
#include <utility/vector1.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>
///#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("core.chemical.modifications.ValenceHandler.cxxtest");

class ValenceHandlerTests : public CxxTest::TestSuite {

public:

	//Shared data elements go here.
	core::chemical::gasteiger::GasteigerAtomTypeSetCOP atom_type_set_;
	core::chemical::ResidueTypeSetCOP const_residue_set_;

	ValenceHandlerTests() {}
	virtual ~ValenceHandlerTests() {}
	void tearDown() {
	}

	void setUp(){
		core_init();
		utility::vector1< std::string > params_files;

		atom_type_set_ = core::chemical::ChemicalManager::get_instance()->gasteiger_atom_type_set();

		const_residue_set_ = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	}

	// Return loaded residue type, but without adding it to a ResidueTypeSet
	core::chemical::MutableResidueTypeOP load_residue_type( std::string name ) {
		std::string const filename( "core/chemical/modifications/"+name+".params" );
		core::chemical::MutableResidueTypeOP residue_type = core::chemical::read_topology_file( filename, const_residue_set_ );

		//this code relies heavily on gasteiger atom types. Make sure to assign them before trying to go through other functions
		core::chemical::gasteiger::assign_gasteiger_atom_types( *residue_type, atom_type_set_, /*keep_existing=*/ false );
		return residue_type;
	}

	void test_linear_coordinates() {
		// Pythagorean Quadruple: 2,3,6,7
		numeric::xyzVector<core::Real> a(1,2,3);
		numeric::xyzVector<core::Real> b(3,5,9);

		numeric::xyzVector<core::Real> x( core::chemical::modifications::linear_coordinates(a,b,7) );

		TS_ASSERT_DELTA( x.distance(a), 7.0, 1e-4 );
		TS_ASSERT_DELTA( numeric::angle_degrees(x,a,b), 180.0, 1e-4);
		//TS_ASSERT_DELTA( c.x(),  5.0, 1e-4);
		//TS_ASSERT_DELTA( c.y(),  8.0, 1e-4);
		//TS_ASSERT_DELTA( c.z(), 15.0, 1e-4);
	}

	void test_angle_coordinates() {
		numeric::xyzVector<core::Real> a(1,2,3);
		numeric::xyzVector<core::Real> b(3,5,9);
		numeric::xyzVector<core::Real> c(2,1,4);

		core::Real length( 1.5 );
		core::Real angle_d( 145.0 );
		core::Real angle_d2( 113.0 );
		core::Real angle( angle_d / 180.0 * numeric::constants::d::pi );
		core::Real angle2( angle_d2 / 180.0 * numeric::constants::d::pi );

		numeric::xyzVector<core::Real> x( core::chemical::modifications::angle_coordinates(
			a, b, c,
			length, angle, angle2,
			false, numeric::xyzVector<core::Real>( 0.0) ) );

		TS_ASSERT_DELTA( x.distance( a ), length, 1e-4);
		TS_ASSERT_DELTA( numeric::angle_degrees(x,a,b), angle_d, 1e-4);
		TS_ASSERT_DELTA( numeric::angle_degrees(x,a,c), angle_d2, 1e-4);
	}

	void test_trigonal_coordinates() {
		numeric::xyzVector<core::Real> a(1,2,3);
		numeric::xyzVector<core::Real> b(3,5,9);
		numeric::xyzVector<core::Real> c(2,1,4);

		core::Real length( 3.4 );
		numeric::xyzVector<core::Real> x( core::chemical::modifications::triganol_coordinates(
			a, b, c, length ) );

		TS_ASSERT_DELTA( x.distance( a ), length, 1e-4);
		core::Real angle( numeric::angle_degrees(x,a,b) );
		// It should be split equally between the two vectors.
		TS_ASSERT_DELTA( numeric::angle_degrees(x,a,c), angle, 1e-2);
		// Also, it should be planar, so all the angle should add up to 360
		TS_ASSERT_DELTA( numeric::angle_degrees(b,a,c) + 2*angle, 360, 1e-2);
	}

	void test_tetrahedral_coordinates() {
		using numeric::constants::d::pi;
		numeric::xyzVector<core::Real> a(1,3.1,3);
		numeric::xyzVector<core::Real> b(1+1.5*sin(60.0/180.0*pi),2,3-1.5*cos(60.0/180.0*pi));
		numeric::xyzVector<core::Real> c(1-1.5*sin(60.0/180.0*pi),2,3-1.5*cos(60.0/180.0*pi));
		numeric::xyzVector<core::Real> d(1,2,3+1.5);

		core::Real length( 1.4 );
		numeric::xyzVector<core::Real> x( core::chemical::modifications::tetrahedral_coordinates(
			a, b, c, d, length ) );

		TS_ASSERT_DELTA( x.distance( a ), length, 1e-4);
		// Since the b,c,d are symmetric around the x-a line, the angles should all match
		core::Real angle( numeric::angle_degrees(x,a,b) );
		TS_ASSERT_DELTA( numeric::angle_degrees(x,a,c), angle, 1e-2);
		TS_ASSERT_DELTA( numeric::angle_degrees(x,a,d), angle, 1e-2);

		TR << "Tet: " << x.x() << " , " << x.y() << " , " << x.z() << std::endl;
	}

	void test_one_bond_sp_missing_zero(){
		//N1 is supposed to be an sp nitrogen with one bond and should have no free valences open for bonds!
		core::chemical::MutableResidueTypeOP restype( load_residue_type("1p1m0") );
		TR << restype->name() << std::endl;
		core::chemical::VD atom = restype->atom_vertex( "N1" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *restype, atom ) );

		TS_ASSERT_EQUALS( coordinates.size(), 0 );
	}

	void test_one_bond_sp_missing_one(){
		//C4 is supposed to be an sp carbon with one bond and should have one free valences open for bonds!
		core::chemical::MutableResidueTypeOP one_bond_sp_( load_residue_type("1p1m1") );
		TR << one_bond_sp_->name() << std::endl;
		core::chemical::VD C4 = one_bond_sp_->atom_vertex( "C4" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *one_bond_sp_, C4 ) );

		TS_ASSERT_EQUALS( coordinates.size(), 1 );
		for ( core::Size i(1); i <= coordinates.size(); ++i ) {
			core::Real sp_angle( numeric::angle_degrees( coordinates[i], one_bond_sp_->atom(C4).ideal_xyz(), one_bond_sp_->atom("C2").ideal_xyz() ) );
			TS_ASSERT_DELTA( sp_angle, 180.0, 1.0 );
		}

		using namespace ObjexxFCL::format;
		for ( core::chemical::VD atm: one_bond_sp_->all_atoms() ) {
			core::PointPosition coords = one_bond_sp_->atom(atm).ideal_xyz();
			std::string outline( "HETATM" + I( 5,  1 ) + "  C   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
		}
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			core::PointPosition coords = coordinates[i];
			std::string outline( "HETATM" + I( 5,  1 ) + " HH   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
			//TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}

		TR << "coordinates: " << std::endl;
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}
	}

	void test_one_bond_sp2_missing_zero(){
		//O1 is supposed to be an sp2 oxygen with one bond and should have no free valences open for bonds!
		core::chemical::MutableResidueTypeOP restype( load_residue_type("1p2m0") );
		TR << restype->name() << std::endl;
		core::chemical::VD atom = restype->atom_vertex( "O1" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *restype, atom ) );

		TS_ASSERT_EQUALS( coordinates.size(), 0 );
	}

	void test_one_bond_sp2_missing_one(){
		//N1 is supposed to be an sp2 nitrogen with one bond and should have one free valences open for bonds!
		core::chemical::MutableResidueTypeOP restype( load_residue_type("1p2m1") );
		TR << restype->name() << std::endl;
		core::chemical::VD atom = restype->atom_vertex( "N1" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *restype, atom ) );

		TS_ASSERT_EQUALS( coordinates.size(), 1 );
		for ( core::Size i(1); i <= coordinates.size(); ++i ) {
			core::Real angle( numeric::angle_degrees( coordinates[i], restype->atom(atom).ideal_xyz(), restype->atom("C1").ideal_xyz() ) );
			TS_ASSERT_DELTA( angle, 120.0, 2.5 );
			core::Real dihedral( numeric::dihedral_degrees( coordinates[i], restype->atom(atom).ideal_xyz(), restype->atom("C1").ideal_xyz(), restype->atom("C2").ideal_xyz()  ) );
			if ( dihedral < -90 ) {
				TS_ASSERT_DELTA( dihedral, -180.0, 1.0 );
			} else if ( dihedral > 90 ) {
				TS_ASSERT_DELTA( dihedral, 180.0, 1.0 );
			} else {
				TS_ASSERT_DELTA( dihedral, 0.0, 1.0 );
			}
		}
	}

	void test_one_bond_sp2_missing_two(){
		//C4 is supposed to be an sp2 carbon with one bond and should have two free valences open for bonds!
		core::chemical::MutableResidueTypeOP one_bond_sp2_( load_residue_type("1p2m2") );
		TR << one_bond_sp2_->name() << std::endl;
		core::chemical::VD C4 = one_bond_sp2_->atom_vertex( "C4" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *one_bond_sp2_, C4 ) );

		TS_ASSERT_EQUALS( coordinates.size(), 2 );
		for ( core::Size i(1); i <= coordinates.size(); ++i ) {
			core::Real sp2_angle( numeric::angle_degrees( coordinates[i], one_bond_sp2_->atom(C4).ideal_xyz(), one_bond_sp2_->atom("C2").ideal_xyz() ) );
			TS_ASSERT_DELTA( sp2_angle, 120.0, 2.5 );
			core::Real sp2_dihedral( numeric::dihedral_degrees( coordinates[i], one_bond_sp2_->atom(C4).ideal_xyz(), one_bond_sp2_->atom("C2").ideal_xyz(), one_bond_sp2_->atom("C1").ideal_xyz()  ) );
			if ( sp2_dihedral < -90 ) {
				TS_ASSERT_DELTA( sp2_dihedral, -180.0, 1.0 );
			} else if ( sp2_dihedral > 90 ) {
				TS_ASSERT_DELTA( sp2_dihedral, 180.0, 1.0 );
			} else {
				TS_ASSERT_DELTA( sp2_dihedral, 0.0, 1.0 );
			}
		}

		using namespace ObjexxFCL::format;
		for ( core::chemical::VD atm: one_bond_sp2_->all_atoms() ) {
			core::PointPosition coords = one_bond_sp2_->atom(atm).ideal_xyz();
			std::string outline( "HETATM" + I( 5,  1 ) + "  C   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
		}
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			core::PointPosition coords = coordinates[i];
			std::string outline( "HETATM" + I( 5,  1 ) + " HH   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
			//TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}

		TR << "coordinates: " << std::endl;
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}
	}

	void test_one_bond_sp3_missing_zero(){
		//CL1 is supposed to be an sp3 chlorine with one bond and should have no free valences open for bonds!
		core::chemical::MutableResidueTypeOP restype( load_residue_type("1p3m0") );
		TR << restype->name() << std::endl;
		core::chemical::VD atom = restype->atom_vertex( "CL1" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *restype, atom ) );

		TS_ASSERT_EQUALS( coordinates.size(), 0 );
	}

	void test_one_bond_sp3_missing_one(){
		//O1 is supposed to be an sp3 oxygen with one bond and should have one valence open for bonds!
		core::chemical::MutableResidueTypeOP one_bond_sp3_( load_residue_type("1p3m1") );
		TR << one_bond_sp3_->name() << std::endl;
		core::chemical::VD O1 = one_bond_sp3_->atom_vertex( "O1" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *one_bond_sp3_, O1 ) );

		TS_ASSERT_EQUALS( coordinates.size(), 1 );
		for ( core::Size i(1); i <= coordinates.size(); ++i ) {
			core::Real sp3_angle( numeric::angle_degrees( coordinates[i], one_bond_sp3_->atom(O1).ideal_xyz(), one_bond_sp3_->atom("C2").ideal_xyz() ) );
			TS_ASSERT_DELTA( sp3_angle, 109.5, 1.0 );
		}

		using namespace ObjexxFCL::format;
		for ( core::chemical::VD atm: one_bond_sp3_->all_atoms() ) {
			core::PointPosition coords = one_bond_sp3_->atom(atm).ideal_xyz();
			std::string outline( "HETATM" + I( 5,  1 ) + "  C   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
		}
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			core::PointPosition coords = coordinates[i];
			std::string outline( "HETATM" + I( 5,  1 ) + " HH   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
			//TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}

		TR << "coordinates: " << std::endl;
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}
	}

	void test_one_bond_sp3_missing_two(){
		//N1 is supposed to be an sp3 nitrogen with one bond and should have two free valences open for bonds!
		core::chemical::MutableResidueTypeOP restype( load_residue_type("1p3m2") );
		TR << restype->name() << std::endl;
		core::chemical::VD atom = restype->atom_vertex( "N1" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *restype, atom ) );

		TS_ASSERT_EQUALS( coordinates.size(), 2 );
		for ( core::Size i(1); i <= coordinates.size(); ++i ) {
			core::Real angle( numeric::angle_degrees( coordinates[i], restype->atom(atom).ideal_xyz(), restype->atom("C1").ideal_xyz() ) );
			TS_ASSERT_DELTA( angle, 109.5, 1.0 );
		}

		using namespace ObjexxFCL::format;
		for ( core::chemical::VD atm: restype->all_atoms() ) {
			core::PointPosition coords = restype->atom(atm).ideal_xyz();
			std::string outline( "HETATM" + I( 5,  1 ) + "  C   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
		}
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			core::PointPosition coords = coordinates[i];
			std::string outline( "HETATM" + I( 5,  1 ) + " HH   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
			//TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}

		TR << "coordinates: " << std::endl;
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}
	}

	void test_one_bond_sp3_missing_three(){
		//C2 is supposed to be an sp3 carbon with one bond and should have three free valences open for bonds!
		core::chemical::MutableResidueTypeOP restype( load_residue_type("1p3m3") );
		TR << restype->name() << std::endl;
		core::chemical::VD atom = restype->atom_vertex( "C2" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *restype, atom ) );

		TS_ASSERT_EQUALS( coordinates.size(), 3 );
		for ( core::Size i(1); i <= coordinates.size(); ++i ) {
			core::Real angle( numeric::angle_degrees( coordinates[i], restype->atom(atom).ideal_xyz(), restype->atom("C1").ideal_xyz() ) );
			TS_ASSERT_DELTA( angle, 109.5, 1.0 );
		}
	}

	void test_two_bond_sp_missing_zero(){
		//C2 is supposed to be an sp carbon with two bonds and should have no free valences open for bonds!
		core::chemical::MutableResidueTypeOP restype( load_residue_type("2p1m0") );
		TR << restype->name() << std::endl;
		core::chemical::VD atom = restype->atom_vertex( "C2" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *restype, atom ) );

		TS_ASSERT_EQUALS( coordinates.size(), 0 );
	}

	void test_two_bond_sp2_missing_zero(){
		//N1 is supposed to be an sp2 nitrogen with two bonds and should have no free valences open for bonds!
		core::chemical::MutableResidueTypeOP restype( load_residue_type("2p2m0") );
		TR << restype->name() << std::endl;
		core::chemical::VD atom = restype->atom_vertex( "N1" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *restype, atom ) );

		TS_ASSERT_EQUALS( coordinates.size(), 0 );
	}

	void test_two_bond_sp2_missing_one(){
		//C4 is supposed to be an sp2 carbon with two bonds and should have one free valence open for bonds!
		core::chemical::MutableResidueTypeOP two_bond_sp2_( load_residue_type("2p2m1") );
		TR << two_bond_sp2_->name() << std::endl;
		core::chemical::VD C4 = two_bond_sp2_->atom_vertex( "C4" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *two_bond_sp2_, C4 ) );

		TS_ASSERT_EQUALS( coordinates.size(), 1 );
		for ( core::Size i(1); i <= coordinates.size(); ++i ) {
			core::Real sp2_angle( numeric::angle_degrees( coordinates[i], two_bond_sp2_->atom(C4).ideal_xyz(), two_bond_sp2_->atom("C2").ideal_xyz() ) );
			TS_ASSERT_DELTA( sp2_angle, 120.0, 1.0 );
			core::Real sp2_dihedral( numeric::dihedral_degrees( coordinates[i], two_bond_sp2_->atom(C4).ideal_xyz(), two_bond_sp2_->atom("C2").ideal_xyz(), two_bond_sp2_->atom("C1").ideal_xyz()  ) );
			if ( sp2_dihedral < -90 ) {
				TS_ASSERT_DELTA( sp2_dihedral, -180.0, 1.0 );
			} else if ( sp2_dihedral > 90 ) {
				TS_ASSERT_DELTA( sp2_dihedral, 180.0, 1.0 );
			} else {
				TS_ASSERT_DELTA( sp2_dihedral, 0.0, 1.0 );
			}
		}

		using namespace ObjexxFCL::format;
		for ( core::chemical::VD atm: two_bond_sp2_->all_atoms() ) {
			core::PointPosition coords = two_bond_sp2_->atom(atm).ideal_xyz();
			std::string outline( "HETATM" + I( 5,  1 ) + "  C   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
		}
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			core::PointPosition coords = coordinates[i];
			std::string outline( "HETATM" + I( 5,  1 ) + " HH   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
			//TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}

		TR << "coordinates: " << std::endl;
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}
	}

	void test_two_bond_sp3_missing_zero(){
		//O1 is supposed to be an sp3 oxygen with two bonds and should have no free valences open for bonds!
		core::chemical::MutableResidueTypeOP restype( load_residue_type("2p3m0") );
		TR << restype->name() << std::endl;
		core::chemical::VD atom = restype->atom_vertex( "O1" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *restype, atom ) );

		TS_ASSERT_EQUALS( coordinates.size(), 0 );
	}

	void test_two_bond_sp3_missing_one(){
		//N1 is supposed to be an sp3 nitrogen with two bonds and should have one free valences open for bonds!
		core::chemical::MutableResidueTypeOP restype( load_residue_type("2p3m1") );
		TR << restype->name() << std::endl;
		core::chemical::VD atom = restype->atom_vertex( "N1" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *restype, atom ) );

		TS_ASSERT_EQUALS( coordinates.size(), 1 );
		for ( core::Size i(1); i <= coordinates.size(); ++i ) {
			core::Real angle( numeric::angle_degrees( coordinates[i], restype->atom(atom).ideal_xyz(), restype->atom("C2").ideal_xyz() ) );
			TS_ASSERT_DELTA( angle, 109.5, 1.0 );
			angle = numeric::angle_degrees( coordinates[i], restype->atom(atom).ideal_xyz(), restype->atom("C1").ideal_xyz() );
			TS_ASSERT_DELTA( angle, 109.5, 1.0 );
		}
	}

	void test_two_bond_sp3_missing_two(){
		//C4 is supposed to be an sp3 carbon with two bonds and should have two free valences open for bonds!
		core::chemical::MutableResidueTypeOP two_bond_sp3_( load_residue_type("2p3m2") );
		TR << two_bond_sp3_->name() << std::endl;
		core::chemical::VD C3 = two_bond_sp3_->atom_vertex( "C3" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *two_bond_sp3_, C3 ) );

		TS_ASSERT_EQUALS( coordinates.size(), 2 );
		for ( core::Size i(1); i <= coordinates.size(); ++i ) {
			core::Real sp3_angle( numeric::angle_degrees( coordinates[i], two_bond_sp3_->atom(C3).ideal_xyz(), two_bond_sp3_->atom("C2").ideal_xyz() ) );
			TS_ASSERT_DELTA( sp3_angle, 109.5, 1.0 );
		}

		using namespace ObjexxFCL::format;
		for ( core::chemical::VD atm: two_bond_sp3_->all_atoms() ) {
			core::PointPosition coords = two_bond_sp3_->atom(atm).ideal_xyz();
			std::string outline( "HETATM" + I( 5,  1 ) + "  C   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
		}
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			core::PointPosition coords = coordinates[i];
			std::string outline( "HETATM" + I( 5,  1 ) + " HH   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
			//TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}

		TR << "coordinates: " << std::endl;
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}
	}

	void test_three_bond_sp3_missing_zero(){
		//N1 is supposed to be an sp3 nitrogen with three bonds and should have no free valences open for bonds!
		core::chemical::MutableResidueTypeOP restype( load_residue_type("3p3m0") );
		TR << restype->name() << std::endl;
		core::chemical::VD atom = restype->atom_vertex( "N1" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *restype, atom ) );

		TS_ASSERT_EQUALS( coordinates.size(), 0 );
	}

	void test_three_bond_sp3_missing_one(){
		//C3 is supposed to be an sp3 carbon with three bonds and should have one free valence open for bonds!
		core::chemical::MutableResidueTypeOP three_bond_sp3_( load_residue_type("3p3m1") );
		TR << three_bond_sp3_->name() << std::endl;
		core::chemical::VD C3 = three_bond_sp3_->atom_vertex( "C3" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *three_bond_sp3_, C3 ) );

		TS_ASSERT_EQUALS( coordinates.size(), 1 );
		for ( core::Size i(1); i <= coordinates.size(); ++i ) {
			core::Real sp3_angle( numeric::angle_degrees( coordinates[i], three_bond_sp3_->atom(C3).ideal_xyz(), three_bond_sp3_->atom("C2").ideal_xyz() ) );
			TS_ASSERT_DELTA( sp3_angle, 109.5, 1.0 );
		}

		using namespace ObjexxFCL::format;
		for ( core::chemical::VD atm: three_bond_sp3_->all_atoms() ) {
			core::PointPosition coords = three_bond_sp3_->atom(atm).ideal_xyz();
			std::string outline( "HETATM" + I( 5,  1 ) + "  C   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
		}
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			core::PointPosition coords = coordinates[i];
			std::string outline( "HETATM" + I( 5,  1 ) + " HH   AAA A"
				+ I( 4, 1 ) + "    "
				+ F( 8, 3, coords.x() ) + F( 8, 3, coords.y() ) + F( 8, 3, coords.z() )
				+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 ));
			TR <<outline <<std::endl;
			//TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}

		TR << "coordinates: " << std::endl;
		for ( core::Size i=1; i<= coordinates.size(); ++i ) {
			TR << coordinates[i].x() << " " << coordinates[i].y() << " " << coordinates[i].z() << std::endl;
		}

	}

	void test_four_bond_sp3_missing_zero(){
		//C1 is supposed to be an sp3 carbon with four bonds and should have no free valences open for bonds!
		core::chemical::MutableResidueTypeOP restype( load_residue_type("3p3m0") );
		TR << restype->name() << std::endl;
		core::chemical::VD atom = restype->atom_vertex( "C1" );
		utility::vector1<numeric::xyzVector<core::Real> > coordinates( core::chemical::modifications::determine_coordinates( *restype, atom ) );

		TS_ASSERT_EQUALS( coordinates.size(), 0 );
	}

};


