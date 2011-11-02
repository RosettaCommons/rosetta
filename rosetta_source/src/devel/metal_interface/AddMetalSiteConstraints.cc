// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file    devel/metal_interface/AddMetalSiteConstraints.cc
/// @brief   Adds metal binding site geometry constraints to pose.
/// @details I typically work with zinc binding sites having tetrahedral coordination geometry, in which case there are 4 distance constraints, 4 angle constraints, 6 tetrahedral constraints, and 4 dihedral constraints.  The code is flexibile enough to accommodate 2, 3, or 4-residue metal binding sites.
/// @author Bryan Der

// Unit Headers
#include <devel/metal_interface/AddMetalSiteConstraints.hh>
#include <devel/metal_interface/MetalSiteResidue.hh>

//Constraint Headers
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/CircularHarmonicFunc.hh>
#include <core/conformation/Residue.hh>

// Project Headers
#include <core/pose/Pose.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/xyz.io.hh> // to print a point
#include <numeric/conversions.hh> //degrees-radians
#include <numeric/xyz.functions.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>
// AUTO-REMOVED #include <utility/file/FileName.hh>
#include <core/id/AtomID.hh>

// C++ headers
#include <string>

#include <numeric/xyzVector.io.hh>


typedef numeric::xyzVector<core::Real> point;

static basic::Tracer TR("devel.metal_interface.AddMetalSiteConstraints");

using namespace core;

namespace devel {
namespace metal_interface {

///@details Adds zinc coordination constraints to a pose.  Zinc site should be parsed with devel/metal_interface/ParseMetalSite, and the resulting vector of MetalSiteResidue objects is needed to initialize this class.
AddMetalSiteConstraints::AddMetalSiteConstraints( utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr )
	: msr_(msr)
{
}

AddMetalSiteConstraints::~AddMetalSiteConstraints()
{
}

///@details Adds distance, tetrahedral angle, angle, and dihedral constraints to pose metal site
void
AddMetalSiteConstraints::add_constraints( pose::Pose & pose ) {

	using namespace scoring::constraints;
	using namespace conformation;
	using namespace id;
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	Real const S_dist ( 2.33 );
	Real const S_dist_dev ( 0.15 );

	Real const N_dist ( 2.05 );
	Real const N_dist_dev ( 0.15 );

	Real const O_dist ( 2.05 );
	Real const O_dist_dev ( 0.15 );

	Real const S_angle ( radians(109.5) );
	Real const S_angle_dev ( radians(20.0) );

	Real const N_angle ( radians(125.0) );
	Real const N_angle_dev ( radians(20.0) );

	Real const O_angle ( radians(120.0) );
	Real const O_angle_dev ( radians(20.0) );

	Real const tetrahedral ( radians(109.5) );
	Real const tetr_dev ( radians(20.0) );

	Real const NE2_dihedral ( radians(180.0) );
	Real const NE2_dev ( radians(20.0) );

	Real const ND1_dihedral ( radians(0.0) );
	Real const ND1_dev ( radians(20.0) );



	FuncOP const S_dist_func( new HarmonicFunc( S_dist, S_dist_dev ) );
	FuncOP const N_dist_func( new HarmonicFunc( N_dist, N_dist_dev ) );
	FuncOP const O_dist_func( new HarmonicFunc( O_dist, O_dist_dev ) );

	FuncOP const S_angle_func( new HarmonicFunc( S_angle, S_angle_dev ) );
	FuncOP const N_angle_func( new HarmonicFunc( N_angle, N_angle_dev ) );
	FuncOP const O_angle_func( new HarmonicFunc( O_angle, O_angle_dev ) );

	FuncOP const tetr_func( new HarmonicFunc( tetrahedral, tetr_dev ) );

	FuncOP const NE2_func( new CircularHarmonicFunc( NE2_dihedral, NE2_dev ) );
	FuncOP const ND1_func( new CircularHarmonicFunc( ND1_dihedral, ND1_dev ) );
	//	FuncOP const O_func( new MinMultiHarmonicFunc( O_dihedrals, O_dihed_devs ) );

	//Print out what we've got in msr_ as a check
	for (Size ii(1); ii <=msr_.size(); ii++) {
		TR << "msr_" << ii << " " << msr_[ii]->get_seqpos() << " " << msr_[ii]->get_ligand_atom_xyz() << " " << msr_[ii]->get_ligand_atom_name() << std::endl;
	}

	//vectors of atom_ids and pre_atom_ids
	utility::vector1< AtomID > atom_ids;
	utility::vector1< AtomID > pre_atom_ids; // zinc is 1, 4 ligands are numbered 2 through 5
	AtomID const zinc_id( msr_[1]->get_ligand_atom_id() );
	atom_ids.push_back( zinc_id );
	pre_atom_ids.push_back( zinc_id );  // not needed, just to keep numbering consistent

	for (Size ii(2); ii <=msr_.size(); ii++) { //start with 2 because zinc is 1

		atom_ids.push_back( msr_[ii]->get_ligand_atom_id() );
		TR << "Ligand " << ii << " atom_id: " << msr_[ii]->get_ligand_atom_id() << std::endl;

		pre_atom_ids.push_back( msr_[ii]->get_pre_ligand_atom_id() );
		TR << "Ligand " << ii << " pre_atom_id: " << msr_[ii]->get_pre_ligand_atom_id() << std::endl;
	}


	// Distances and Angles
	for( core::Size i(2); i <= msr_.size(); ++i ) { // zinc is 1
		std::string resn( msr_[i]->get_resname() );
		TR << "RES NAME " << resn << std::endl;
		if ( resn == "CYS" ) { // CYS
			TR << "Adding Cys atom_pair_constraint " << atom_ids[i] << " " << zinc_id << std::endl;
			pose.add_constraint( new AtomPairConstraint( atom_ids[i], zinc_id, S_dist_func ) );
			TR << "Adding Cys angle_constraint " << pre_atom_ids[i] << " " << atom_ids[i] << " " << zinc_id << std::endl;
			pose.add_constraint( new AngleConstraint( pre_atom_ids[i], atom_ids[i], zinc_id, S_angle_func ) );
		}
		else if ( resn == "HIS") { // HIS
			TR << "Adding His atom_pair_constraint " << atom_ids[i] << " " << zinc_id << std::endl;
			pose.add_constraint( new AtomPairConstraint( atom_ids[i], zinc_id, N_dist_func ) );
			TR << "Adding His angle_constraint " << pre_atom_ids[i] << " " << atom_ids[i] << " " << zinc_id << std::endl;
			pose.add_constraint( new AngleConstraint( pre_atom_ids[i], atom_ids[i], zinc_id, N_angle_func ) );
		}
		else if ( resn == "ASP" || resn == "GLU" ) { // ASP/GLU
			TR << "Adding Asp/Glu atom_pair_constraint " << atom_ids[i] << " " << zinc_id << std::endl;
			pose.add_constraint( new AtomPairConstraint( atom_ids[i], zinc_id, O_dist_func ) );
			TR << "Adding Asp/Glu angle_constraint " << pre_atom_ids[i] << " " << atom_ids[i] << " " << zinc_id << std::endl;
			pose.add_constraint( new AngleConstraint( pre_atom_ids[i], atom_ids[i], zinc_id, O_angle_func ) );
		}
	}


	// Tetrahedral Angles (not totally generalized code, but it was the easy and quick solution)
	if( msr_.size() == 5 ) {
		TR << "Adding tetrahedral angle constraints" << std::endl;
		pose.add_constraint( new AngleConstraint( msr_[2]->get_ligand_atom_id(), zinc_id, msr_[3]->get_ligand_atom_id(), tetr_func ) );
		pose.add_constraint( new AngleConstraint( msr_[2]->get_ligand_atom_id(), zinc_id, msr_[4]->get_ligand_atom_id(), tetr_func ) );
		pose.add_constraint( new AngleConstraint( msr_[2]->get_ligand_atom_id(), zinc_id, msr_[5]->get_ligand_atom_id(), tetr_func ) );
		pose.add_constraint( new AngleConstraint( msr_[3]->get_ligand_atom_id(), zinc_id, msr_[4]->get_ligand_atom_id(), tetr_func ) );
		pose.add_constraint( new AngleConstraint( msr_[3]->get_ligand_atom_id(), zinc_id, msr_[5]->get_ligand_atom_id(), tetr_func ) );
		pose.add_constraint( new AngleConstraint( msr_[4]->get_ligand_atom_id(), zinc_id, msr_[5]->get_ligand_atom_id(), tetr_func ) );
	}
	else if( msr_.size() == 4) {
		TR << "Adding tetrahedral angle constraints" << std::endl;
		pose.add_constraint( new AngleConstraint( msr_[2]->get_ligand_atom_id(), zinc_id, msr_[3]->get_ligand_atom_id(), tetr_func ) );
		pose.add_constraint( new AngleConstraint( msr_[2]->get_ligand_atom_id(), zinc_id, msr_[4]->get_ligand_atom_id(), tetr_func ) );
		pose.add_constraint( new AngleConstraint( msr_[3]->get_ligand_atom_id(), zinc_id, msr_[4]->get_ligand_atom_id(), tetr_func ) );
	}
	else if( msr_.size() == 3) {
		TR << "Adding tetrahedral angle constraints" << std::endl;
		pose.add_constraint( new AngleConstraint( msr_[2]->get_ligand_atom_id(), zinc_id, msr_[3]->get_ligand_atom_id(), tetr_func ) );
	}

	// Dihedrals
	for( core::Size i(2); i <= msr_.size(); ++i ) { // zinc is 1
		TR << "Adding dihedral constraints for i= " << i << std::endl;
		if( msr_[i]->get_resname() == "CYS" ) { continue; }

		core::id::AtomID id1( msr_[i]->get_pre_pre_ligand_atom_id() );
		core::id::AtomID id2( msr_[i]->get_pre_ligand_atom_id() );
		core::id::AtomID id3( msr_[i]->get_ligand_atom_id() );
		core::id::AtomID id4( zinc_id );

		TR << "id1 " << id1 << std::endl;
		TR << "id2 " << id2 << std::endl;
		TR << "id3 " << id3 << std::endl;
		TR << "id4 " << id4 << std::endl;

		point p1( pose.residue(id1.rsd()).atom(id1.atomno()).xyz() );
		point p2( pose.residue(id2.rsd()).atom(id2.atomno()).xyz() );
		point p3( pose.residue(id3.rsd()).atom(id3.atomno()).xyz() );
		point p4( pose.residue(id4.rsd()).atom(id4.atomno()).xyz() );

		TR << "Check_point " << id1 << " " << p1 << std::endl;
		TR << "Check_point " << id2 << " " << p2 << std::endl;
		TR << "Check_point " << id3 << " " << p3 << std::endl;
		TR << "Check_point " << id4 << " " << p4 << std::endl;

		FuncOP dihedral_func;

		if ( msr_[i]->get_ligand_atom_name() == " ND1" ) {
			TR << "Adding dihedral constraint for ND1 His" << std::endl;
			dihedral_func = new CircularHarmonicFunc( ND1_dihedral, ND1_dev );
		}
		else if ( msr_[i]->get_ligand_atom_name() == " NE2" ) {
			TR << "Adding dihedral constraint for NE2 His" << std::endl;
			dihedral_func = new CircularHarmonicFunc( NE2_dihedral, NE2_dev );
		}

		else if ( msr_[i]->get_resname() == "ASP" || msr_[i]->get_resname() == "GLU" ) {
			TR << "Adding dihedral constraint for Asp or Glu" << std::endl;
			Real dihed = numeric::dihedral_degrees( p1, p2, p3, p4);
			TR << "Computed_DIHEDRAL " << i << " " << dihed << std::endl;
			Real dihed_diff_180 = fabs(fabs(dihed) - 180.0);
			Real dihed_diff_0 = fabs(dihed);

			if(dihed_diff_180 < dihed_diff_0) {
				TR << "Adding dihedral constraint for Oxy_anti" << std::endl;
				dihedral_func = new CircularHarmonicFunc( NE2_dihedral, NE2_dev );
			}

			else if(dihed_diff_0 < dihed_diff_180) {
				TR << "Adding dihedral constraint for Oxy_syn" << std::endl;
				dihedral_func = new CircularHarmonicFunc( ND1_dihedral, ND1_dev );
			}
		}

		//this will be performed foreach His/Asp/Glu dihedral
		scoring::constraints::ConstraintOP dihed_constraint = new scoring::constraints::DihedralConstraint( id1, id2, id3, id4, dihedral_func );
		scoring::constraints::DihedralConstraint dihed_constraint2( id1, id2, id3, id4, dihedral_func );
		Real dihed_constraint_score = dihed_constraint2.score( p1, p2, p3, p4 );
		TR << "DIHEDRAL_SCORE= " << dihed_constraint_score << " " << id1 << " " << id2 << " " << id3 <<  " " << id4 << std::endl;
		TR << "Adding dihed constraint " << std::endl;
		pose.add_constraint( dihed_constraint);
		TR << "Added dihed constraint " << std::endl;
		//pose.add_constraint( new DihedralConstraint( id1, id2, id3, id4, dihedral_func ) );

	}



	//pose.constraint_set()->show(std::cout); //show and show_definition do not compile, for some reason
	//pose.constraint_set()->show_definition(std::cout, pose);
	pose.constraint_set()->show_numbers(std::cout);

	return;
}


} //metal_interface
} //devel
