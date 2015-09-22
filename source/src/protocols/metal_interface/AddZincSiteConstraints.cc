// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/metal_interface/AddZincSiteConstraints.cc
/// @brief   Adds metal binding site geometry constraints to pose.
/// @details I typically work with zinc binding sites having tetrahedral coordination geometry, in which case there are 4 distance constraints, 4 angle constraints, 6 tetrahedral constraints, and 4 dihedral constraints.  The code is flexibile enough to accommodate 2, 3, or 4-residue metal binding sites.
/// @author Bryan Der

// Unit Headers
#include <protocols/metal_interface/AddZincSiteConstraints.hh>
#include <protocols/metal_interface/MetalSiteResidue.hh>

//Constraint Headers
#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/MinMultiHarmonicFunc.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh> //degrees-radians
#include <numeric/xyz.functions.hh>
#include <core/id/AtomID.hh>

// C++ headers
#include <sstream>
#include <utility/io/ozstream.hh> // used to create a resfile
#include <utility/file/FileName.hh>
#include <numeric/xyzVector.io.hh>


typedef numeric::xyzVector<core::Real> point;

static THREAD_LOCAL basic::Tracer TR( "protocols.metal_interface.AddZincSiteConstraints" );
static THREAD_LOCAL basic::Tracer TR_PYMOL( "TR_PYMOL" );

using namespace core;
using core::id::AtomID;

namespace protocols {
namespace metal_interface {

/// @details Adds zinc coordination constraints to a pose.  Zinc site should be parsed with protocols/metal_interface/ZincSiteFinder, and the resulting vector of MetalSiteResidue objects is needed to initialize this class.
AddZincSiteConstraints::AddZincSiteConstraints( utility::vector1< protocols::metal_interface::MetalSiteResidueOP > msr )
: msr_(msr)
{
}

AddZincSiteConstraints::~AddZincSiteConstraints()
{
}

/// @details Adds distance, tetrahedral angle, angle, and dihedral constraints to pose metal site
void
AddZincSiteConstraints::add_constraints( pose::Pose & pose ) {

	pdbname_ = pose.pdb_info()->name();

	distance_constraints_.clear();
	angle_constraints_.clear();
	dihedral_constraints_.clear();
	tetrahedral_constraints_.clear();


	using namespace scoring::constraints;
	using namespace conformation;
	using namespace id;
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	Real const S_dist ( 2.33 );
	Real const S_dist_dev ( 0.12 );

	Real const N_dist ( 2.10 );
	Real const N_dist_dev ( 0.12 );

	Real const O_dist ( 2.05 );
	Real const O_dist_dev ( 0.12 );

	Real const S_angle ( radians(105.0) );
	Real const S_angle_dev ( radians(10.0) );

	Real const N_angle ( radians(125.0) );
	Real const N_angle_dev ( radians(10.0) );

	Real const O_angle ( radians(120.0) );
	Real const O_angle_dev ( radians(10.0) );

	Real const tetrahedral ( radians(109.0) );
	Real const tetr_dev ( radians(15.0) );

	Real const NE2_dihedral ( radians(180.0) );
	Real const NE2_dev ( radians(10.0) );

	Real const ND1_dihedral ( radians(0.0) );
	Real const ND1_dev ( radians(10.0) );

	utility::vector1< Real > SG_dihedrals;
	SG_dihedrals.clear();
	SG_dihedrals.push_back( radians(-180.0) );
	SG_dihedrals.push_back( radians( 180.0) );
	SG_dihedrals.push_back( radians( -90.0) );
	SG_dihedrals.push_back( radians(  90.0) );

	utility::vector1< Real > SG_devs;
	SG_devs.clear();
	SG_devs.push_back( radians(15.0) );
	SG_devs.push_back( radians(15.0) );
	SG_devs.push_back( radians(15.0) );
	SG_devs.push_back( radians(15.0) );


	core::scoring::func::FuncOP const S_dist_func( new core::scoring::func::HarmonicFunc( S_dist, S_dist_dev ) );
	core::scoring::func::FuncOP const N_dist_func( new core::scoring::func::HarmonicFunc( N_dist, N_dist_dev ) );
	core::scoring::func::FuncOP const O_dist_func( new core::scoring::func::HarmonicFunc( O_dist, O_dist_dev ) );

	core::scoring::func::FuncOP const S_angle_func( new core::scoring::func::HarmonicFunc( S_angle, S_angle_dev ) );
	core::scoring::func::FuncOP const N_angle_func( new core::scoring::func::HarmonicFunc( N_angle, N_angle_dev ) );
	core::scoring::func::FuncOP const O_angle_func( new core::scoring::func::HarmonicFunc( O_angle, O_angle_dev ) );

	core::scoring::func::FuncOP const tetr_func( new core::scoring::func::HarmonicFunc( tetrahedral, tetr_dev ) );

	//these dihedral funcs are initialized later
	//FuncOP const NE2_func( new CircularHarmonicFunc( NE2_dihedral, NE2_dev ) );
	//FuncOP const ND1_func( new CircularHarmonicFunc( ND1_dihedral, ND1_dev ) );
	//FuncOP const SG_func( new MinMultiHarmonicFunc( SG_dihedrals, SG_devs ) );


	//Print out what we've got in msr_ as a check
	for ( Size ii(1); ii <=msr_.size(); ii++ ) {
		TR << "msr_" << ii << " " << msr_[ii]->get_seqpos() << " " << msr_[ii]->get_ligand_atom_xyz() << " " << msr_[ii]->get_ligand_atom_name() << std::endl;
	}

	//vectors of atom_ids and pre_atom_ids
	utility::vector1< AtomID > atom_ids;
	utility::vector1< AtomID > pre_atom_ids; // zinc is 1, 4 ligands are numbered 2 through 5
	AtomID const zinc_id( msr_[1]->get_ligand_atom_id() );
	atom_ids.push_back( zinc_id );
	pre_atom_ids.push_back( zinc_id );  // not needed, just to keep numbering consistent

	for ( Size ii(2); ii <=msr_.size(); ii++ ) { //start with 2 because zinc is 1

		atom_ids.push_back( msr_[ii]->get_ligand_atom_id() );
		TR << "Ligand " << ii << " atom_id: " << msr_[ii]->get_ligand_atom_id() << std::endl;

		pre_atom_ids.push_back( msr_[ii]->get_pre_ligand_atom_id() );
		TR << "Ligand " << ii << " pre_atom_id: " << msr_[ii]->get_pre_ligand_atom_id() << std::endl;
	}


	// Distances and Angles
	for ( core::Size i(2); i <= msr_.size(); ++i ) { // zinc is 1
		std::string resn( msr_[i]->get_resname() );
		TR << "RES NAME " << resn << std::endl;
		if ( resn == "CYS" ) { // CYS
			TR << "Adding Cys atom_pair_constraint " << atom_ids[i] << " " << zinc_id << std::endl;
			distance_constraints_.push_back( core::scoring::constraints::AtomPairConstraintCOP( core::scoring::constraints::AtomPairConstraintOP( new AtomPairConstraint( atom_ids[i], zinc_id, S_dist_func ) ) ) );
			pose.add_constraint( distance_constraints_[ distance_constraints_.size() ] );

			TR << "Adding Cys angle_constraint " << pre_atom_ids[i] << " " << atom_ids[i] << " " << zinc_id << std::endl;
			angle_constraints_.push_back( core::scoring::constraints::AngleConstraintCOP( core::scoring::constraints::AngleConstraintOP( new AngleConstraint( pre_atom_ids[i], atom_ids[i], zinc_id, S_angle_func ) ) ) );
			pose.add_constraint( angle_constraints_[ angle_constraints_.size() ] );
		} else if ( resn == "HIS" ) { // HIS
			TR << "Adding His atom_pair_constraint " << atom_ids[i] << " " << zinc_id << std::endl;
			distance_constraints_.push_back( core::scoring::constraints::AtomPairConstraintCOP( core::scoring::constraints::AtomPairConstraintOP( new AtomPairConstraint( atom_ids[i], zinc_id, N_dist_func ) ) ) );
			pose.add_constraint( distance_constraints_[ distance_constraints_.size() ] );

			TR << "Adding His angle_constraint " << pre_atom_ids[i] << " " << atom_ids[i] << " " << zinc_id << std::endl;
			angle_constraints_.push_back( core::scoring::constraints::AngleConstraintCOP( core::scoring::constraints::AngleConstraintOP( new AngleConstraint( pre_atom_ids[i], atom_ids[i], zinc_id, N_angle_func ) ) ) );
			pose.add_constraint( angle_constraints_[ angle_constraints_.size() ] );
		} else if ( resn == "ASP" || resn == "GLU" ) { // ASP/GLU
			TR << "Adding Asp/Glu atom_pair_constraint " << atom_ids[i] << " " << zinc_id << std::endl;
			distance_constraints_.push_back( core::scoring::constraints::AtomPairConstraintCOP( core::scoring::constraints::AtomPairConstraintOP( new AtomPairConstraint( atom_ids[i], zinc_id, O_dist_func ) ) ) );
			pose.add_constraint( distance_constraints_[ distance_constraints_.size() ] );

			TR << "Adding Asp/Glu angle_constraint " << pre_atom_ids[i] << " " << atom_ids[i] << " " << zinc_id << std::endl;
			angle_constraints_.push_back( core::scoring::constraints::AngleConstraintCOP( core::scoring::constraints::AngleConstraintOP( new AngleConstraint( pre_atom_ids[i], atom_ids[i], zinc_id, O_angle_func ) ) ) );
			pose.add_constraint( angle_constraints_[ angle_constraints_.size() ] );
		}
	}


	// Tetrahedral Angles (not totally generalized code, but it was the easy and quick solution)
	TR << "Adding tetrahedral angle constraints" << std::endl;
	if ( msr_.size() == 5 ) {
		AngleConstraintCOP tetr1_constraint( AngleConstraintOP( new AngleConstraint( msr_[2]->get_ligand_atom_id(), zinc_id, msr_[3]->get_ligand_atom_id(), tetr_func ) ) );
		AngleConstraintCOP tetr2_constraint( AngleConstraintOP( new AngleConstraint( msr_[2]->get_ligand_atom_id(), zinc_id, msr_[4]->get_ligand_atom_id(), tetr_func ) ) );
		AngleConstraintCOP tetr3_constraint( AngleConstraintOP( new AngleConstraint( msr_[2]->get_ligand_atom_id(), zinc_id, msr_[5]->get_ligand_atom_id(), tetr_func ) ) );
		AngleConstraintCOP tetr4_constraint( AngleConstraintOP( new AngleConstraint( msr_[3]->get_ligand_atom_id(), zinc_id, msr_[4]->get_ligand_atom_id(), tetr_func ) ) );
		AngleConstraintCOP tetr5_constraint( AngleConstraintOP( new AngleConstraint( msr_[3]->get_ligand_atom_id(), zinc_id, msr_[5]->get_ligand_atom_id(), tetr_func ) ) );
		AngleConstraintCOP tetr6_constraint( AngleConstraintOP( new AngleConstraint( msr_[4]->get_ligand_atom_id(), zinc_id, msr_[5]->get_ligand_atom_id(), tetr_func ) ) );
		tetrahedral_constraints_.push_back( tetr1_constraint );
		tetrahedral_constraints_.push_back( tetr2_constraint );
		tetrahedral_constraints_.push_back( tetr3_constraint );
		tetrahedral_constraints_.push_back( tetr4_constraint );
		tetrahedral_constraints_.push_back( tetr5_constraint );
		tetrahedral_constraints_.push_back( tetr6_constraint );
		pose.add_constraint( tetr1_constraint );
		pose.add_constraint( tetr2_constraint );
		pose.add_constraint( tetr3_constraint );
		pose.add_constraint( tetr4_constraint );
		pose.add_constraint( tetr5_constraint );
		pose.add_constraint( tetr6_constraint );
	} else if ( msr_.size() == 4 ) {
		AngleConstraintCOP tetr1_constraint( AngleConstraintOP( new AngleConstraint( msr_[2]->get_ligand_atom_id(), zinc_id, msr_[3]->get_ligand_atom_id(), tetr_func ) ) );
		AngleConstraintCOP tetr2_constraint( AngleConstraintOP( new AngleConstraint( msr_[2]->get_ligand_atom_id(), zinc_id, msr_[4]->get_ligand_atom_id(), tetr_func ) ) );
		AngleConstraintCOP tetr3_constraint( AngleConstraintOP( new AngleConstraint( msr_[3]->get_ligand_atom_id(), zinc_id, msr_[4]->get_ligand_atom_id(), tetr_func ) ) );
		tetrahedral_constraints_.push_back( tetr1_constraint );
		tetrahedral_constraints_.push_back( tetr2_constraint );
		tetrahedral_constraints_.push_back( tetr3_constraint );
		pose.add_constraint( tetr1_constraint );
		pose.add_constraint( tetr2_constraint );
		pose.add_constraint( tetr3_constraint );
	} else if ( msr_.size() == 3 ) {
		AngleConstraintCOP tetr1_constraint( AngleConstraintOP( new AngleConstraint( msr_[2]->get_ligand_atom_id(), zinc_id, msr_[3]->get_ligand_atom_id(), tetr_func ) ) );
		tetrahedral_constraints_.push_back( tetr1_constraint );
		pose.add_constraint( tetr1_constraint );
	}

	// Dihedrals
	for ( core::Size i(2); i <= msr_.size(); ++i ) { // zinc is 1
		//if( msr_[i]->get_resname() == "CYS" ) { continue; }

		AtomID id1( msr_[i]->get_pre_pre_ligand_atom_id() );
		AtomID id2( msr_[i]->get_pre_ligand_atom_id() );
		AtomID id3( msr_[i]->get_ligand_atom_id() );
		AtomID id4( zinc_id );

		TR << "id1 " << id1 << std::endl;
		TR << "id2 " << id2 << std::endl;
		TR << "id3 " << id3 << std::endl;
		TR << "id4 " << id4 << std::endl;

		point p1( pose.residue(id1.rsd()).atom(id1.atomno()).xyz() );
		point p2( pose.residue(id2.rsd()).atom(id2.atomno()).xyz() );
		point p3( pose.residue(id3.rsd()).atom(id3.atomno()).xyz() );
		point p4( pose.residue(id4.rsd()).atom(id4.atomno()).xyz() );

		// TR << "Check_point " << id1 << " " << p1 << std::endl;
		// TR << "Check_point " << id2 << " " << p2 << std::endl;
		// TR << "Check_point " << id3 << " " << p3 << std::endl;
		// TR << "Check_point " << id4 << " " << p4 << std::endl;

		core::scoring::func::FuncOP dihedral_func;

		if ( msr_[i]->get_ligand_atom_name() == " ND1" ) {
			TR << "Adding dihedral constraint for ND1 His" << std::endl;
			dihedral_func = core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( ND1_dihedral, ND1_dev ) );
			dihedral_constraints_.push_back( core::scoring::constraints::DihedralConstraintCOP( core::scoring::constraints::DihedralConstraintOP( new DihedralConstraint( id1, id2, id3, id4, dihedral_func ) ) ) );
			pose.add_constraint( dihedral_constraints_[ dihedral_constraints_.size() ] );
		} else if ( msr_[i]->get_ligand_atom_name() == " NE2" ) {
			TR << "Adding dihedral constraint for NE2 His" << std::endl;
			dihedral_func = core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( NE2_dihedral, NE2_dev ) );
			dihedral_constraints_.push_back( core::scoring::constraints::DihedralConstraintCOP( core::scoring::constraints::DihedralConstraintOP( new DihedralConstraint( id1, id2, id3, id4, dihedral_func ) ) ) );
			pose.add_constraint( dihedral_constraints_[ dihedral_constraints_.size() ] );
		} else if ( msr_[i]->get_resname() == "ASP" || msr_[i]->get_resname() == "GLU" ) {
			TR << "Adding dihedral constraint for Asp or Glu" << std::endl;
			Real dihed = numeric::dihedral_degrees( p1, p2, p3, p4);
			TR << "Computed_DIHEDRAL " << i << " " << dihed << std::endl;
			Real dihed_diff_180 = fabs(fabs(dihed) - 180.0);
			Real dihed_diff_0 = fabs(dihed);

			if ( dihed_diff_180 < dihed_diff_0 ) {
				TR << "Adding dihedral constraint for Oxy_anti" << std::endl;
				dihedral_func = core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( NE2_dihedral, NE2_dev ) );
				dihedral_constraints_.push_back( core::scoring::constraints::DihedralConstraintCOP( core::scoring::constraints::DihedralConstraintOP( new DihedralConstraint( id1, id2, id3, id4, dihedral_func ) ) ) );
				pose.add_constraint( dihedral_constraints_[ dihedral_constraints_.size() ] );
			} else if ( dihed_diff_0 < dihed_diff_180 ) {
				TR << "Adding dihedral constraint for Oxy_syn" << std::endl;
				dihedral_func = core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( ND1_dihedral, ND1_dev ) );
				dihedral_constraints_.push_back( core::scoring::constraints::DihedralConstraintCOP( core::scoring::constraints::DihedralConstraintOP( new DihedralConstraint( id1, id2, id3, id4, dihedral_func ) ) ) );
				pose.add_constraint( dihedral_constraints_[ dihedral_constraints_.size() ] );
			}
		} else if ( msr_[i]->get_ligand_atom_name() == " SG " ) {
			TR << "Adding dihedral constraint for SG Cys" << std::endl;
			dihedral_func = core::scoring::func::FuncOP( new core::scoring::func::MinMultiHarmonicFunc( SG_dihedrals, SG_devs ) );
			dihedral_constraints_.push_back( core::scoring::constraints::DihedralConstraintCOP( core::scoring::constraints::DihedralConstraintOP( new DihedralConstraint( id1, id2, id3, id4, dihedral_func ) ) ) );
			pose.add_constraint( dihedral_constraints_[ dihedral_constraints_.size() ] );
		}


	}


	//pose.constraint_set()->show(std::cout); //show and show_definition do not compile, for some reason
	//pose.constraint_set()->show_definition(std::cout, pose);
	pose.constraint_set()->show_numbers(std::cout);

	return;
}


void
AddZincSiteConstraints::evaluate_constraints( core::pose::Pose const & pose ) {
	if ( distance_constraints_.size() == 0 ) {
		TR << "The constraints have not been added, no info to give." << std::endl;
		return;
	}

	Real distance_total_score( 0.0 );
	Real angle_total_score( 0.0 );
	Real dihedral_total_score( 0.0 );
	Real tetrahedral_total_score( 0.0 );

	//Distances
	for ( Size i(1); i <= distance_constraints_.size(); ++i ) {
		AtomID lig_atom_id( distance_constraints_[i]->atom(1) );
		AtomID zinc_id    ( distance_constraints_[i]->atom(2) );

		TR << "CONSTRAINT: Distance: "
			<< pose.residue(lig_atom_id.rsd()).name3() << " "
			<< lig_atom_id.rsd()  << " "
			<< pose.residue(lig_atom_id.rsd()).atom_name(lig_atom_id.atomno()) << " "
			<< distance_constraints_[i]->dist( pose ) << " "
			<< distance_constraints_[i]->score( pose ) << std::endl;

		distance_total_score += distance_constraints_[i]->score( pose );
	}

	//Angles
	for ( Size i(1); i <= angle_constraints_.size(); ++i ) {
		AtomID pre_lig_atom_id( angle_constraints_[i]->atom(1) );
		AtomID lig_atom_id    ( angle_constraints_[i]->atom(2) );
		AtomID zinc_id        ( angle_constraints_[i]->atom(3) );

		point xyz1 = pose.residue( pre_lig_atom_id.rsd() ).atom( pre_lig_atom_id.atomno() ).xyz();
		point xyz2 = pose.residue( lig_atom_id.rsd() ).atom( lig_atom_id.atomno() ).xyz();
		point xyz3 = pose.residue( zinc_id.rsd() ).atom( zinc_id.atomno() ).xyz();

		//AngleConstraint does not have a function to return the current value of the angle
		Real angle = numeric::conversions::degrees( angle_of( xyz1, xyz2, xyz3 ) );

		TR << "CONSTRAINT: Angle   : "
			<< pose.residue(lig_atom_id.rsd()).name3() << " "
			<< lig_atom_id.rsd()  << " "
			<< pose.residue(lig_atom_id.rsd()).atom_name(lig_atom_id.atomno()) << " "
			<< angle << " "
			<< angle_constraints_[i]->score( xyz1, xyz2, xyz3 ) << std::endl;

		angle_total_score += angle_constraints_[i]->score( xyz1, xyz2, xyz3 );
	}

	//Dihedrals
	for ( Size i(1); i <= dihedral_constraints_.size(); ++i ) {

		AtomID pre_pre_lig_atom_id( dihedral_constraints_[i]->atom(1) );
		AtomID pre_lig_atom_id    ( dihedral_constraints_[i]->atom(2) );
		AtomID lig_atom_id        ( dihedral_constraints_[i]->atom(3) );
		AtomID zinc_id            ( dihedral_constraints_[i]->atom(4) );

		//get points from the current pose, not the pose from ZincSiteFinder when the constraints were first added
		point xyz1 = pose.residue( pre_pre_lig_atom_id.rsd() ).atom( pre_pre_lig_atom_id.atomno() ).xyz();
		point xyz2 = pose.residue( pre_lig_atom_id.rsd() ).atom( pre_lig_atom_id.atomno() ).xyz();
		point xyz3 = pose.residue( lig_atom_id.rsd() ).atom( lig_atom_id.atomno() ).xyz();
		point xyz4 = pose.residue( zinc_id.rsd() ).atom( zinc_id.atomno() ).xyz();

		//DihedralConstraint does not have a function to return the current value of the dihedral
		Real dihedral = numeric::dihedral_degrees( xyz1, xyz2, xyz3, xyz4 );

		TR << "CONSTRAINT: Dihedral: "
			<< pose.residue(lig_atom_id.rsd()).name3() << " "
			<< lig_atom_id.rsd() << " "
			<< pose.residue(lig_atom_id.rsd()).atom_name(lig_atom_id.atomno()) << " "
			<< dihedral << " "
			<< dihedral_constraints_[i]->score( xyz1, xyz2, xyz3, xyz4 ) << std::endl;

		dihedral_total_score += dihedral_constraints_[i]->score( xyz1, xyz2, xyz3, xyz4 );
	}


	//Tetrahedral
	for ( Size i(1); i <= tetrahedral_constraints_.size(); ++i ) {

		AtomID lig_atom_id_1( tetrahedral_constraints_[i]->atom(1) );
		AtomID       zinc_id( tetrahedral_constraints_[i]->atom(2) );
		AtomID lig_atom_id_2( tetrahedral_constraints_[i]->atom(3) );

		point xyz1 = pose.residue( lig_atom_id_1.rsd() ).atom( lig_atom_id_1.atomno() ).xyz();
		point xyz2 = pose.residue(       zinc_id.rsd() ).atom(       zinc_id.atomno() ).xyz();
		point xyz3 = pose.residue( lig_atom_id_2.rsd() ).atom( lig_atom_id_2.atomno() ).xyz();

		//AngleConstraint does not have a function to return the current value of the angle
		Real tetr_angle = numeric::conversions::degrees( angle_of( xyz1, xyz2, xyz3 ) );

		TR << "CONSTRAINT: Tetrahedral   : "
			<< lig_atom_id_1 << " and "
			<< lig_atom_id_2 << " "
			<< tetr_angle << " "
			<< tetrahedral_constraints_[i]->score( xyz1, xyz2, xyz3 ) << std::endl;

		tetrahedral_total_score += tetrahedral_constraints_[i]->score( xyz1, xyz2, xyz3 );
	}

	Real distance_score    = distance_total_score / distance_constraints_.size();
	Real angle_score       = angle_total_score / distance_constraints_.size();
	Real dihedral_score    = dihedral_total_score / distance_constraints_.size();
	Real tetrahedral_score = tetrahedral_total_score / distance_constraints_.size();

	TR << pdbname_ << "  Total: Distance   : " << distance_score    << std::endl;
	TR << pdbname_ << "  Total: Angle      : " << angle_score       << std::endl;
	TR << pdbname_ << "  Total: Dihedral   : " << dihedral_score    << std::endl;
	TR << pdbname_ << "  Total: Tetrahedral: " << tetrahedral_score << std::endl;

	TR << pdbname_ << "         TOTAL_SCORE: " << distance_score + angle_score + dihedral_score + tetrahedral_score << std::endl;

	return;
}


//assumes only one pdb opened in pymol
void
AddZincSiteConstraints::view_constraints_in_pymol( core::pose::Pose const & pose ) {
	if ( distance_constraints_.size() == 0 ) {
		TR_PYMOL << "The constraints have not been added, no info to give." << std::endl;
		return;
	}

	//Distances
	for ( Size i(1); i <= distance_constraints_.size(); ++i ) {
		AtomID lig_atom_id( distance_constraints_[i]->atom(1) );
		AtomID zinc_id    ( distance_constraints_[i]->atom(2) );
		std::stringstream ss;
		ss << "dist_" << lig_atom_id.rsd();
		TR_PYMOL << "select " << ss.str() << ", resi "
			<< lig_atom_id.rsd() << " and name " << pose.residue(lig_atom_id.rsd()).atom_name(lig_atom_id.atomno()) << " + resi "
			<< zinc_id.rsd() << " and name " << pose.residue(zinc_id.rsd()).atom_name(zinc_id.atomno()) << std::endl;
	}

	//Angles
	for ( Size i(1); i <= angle_constraints_.size(); ++i ) {
		AtomID pre_lig_atom_id( angle_constraints_[i]->atom(1) );
		AtomID lig_atom_id    ( angle_constraints_[i]->atom(2) );
		AtomID zinc_id        ( angle_constraints_[i]->atom(3) );
		std::stringstream ss;
		ss << "angle_" << lig_atom_id.rsd();
		TR_PYMOL << "select " << ss.str() << ", resi "
			<< pre_lig_atom_id.rsd() << " and name " << pose.residue(pre_lig_atom_id.rsd()).atom_name(pre_lig_atom_id.atomno()) << " + resi "
			<< lig_atom_id.rsd() << " and name " << pose.residue(lig_atom_id.rsd()).atom_name(lig_atom_id.atomno()) << " + resi "
			<< zinc_id.rsd() << " and name " << pose.residue(zinc_id.rsd()).atom_name(zinc_id.atomno()) << std::endl;

	}

	//Dihedrals
	for ( Size i(1); i <= dihedral_constraints_.size(); ++i ) {
		AtomID pre_pre_lig_atom_id( dihedral_constraints_[i]->atom(1) );
		AtomID pre_lig_atom_id    ( dihedral_constraints_[i]->atom(2) );
		AtomID lig_atom_id        ( dihedral_constraints_[i]->atom(3) );
		AtomID zinc_id            ( dihedral_constraints_[i]->atom(4) );
		std::stringstream ss;
		ss << "dihed_" << lig_atom_id.rsd();
		TR_PYMOL << "select " << ss.str() << ", resi "
			<< pre_pre_lig_atom_id.rsd() << " and name " << pose.residue(pre_pre_lig_atom_id.rsd()).atom_name(pre_pre_lig_atom_id.atomno()) << " + resi "
			<< pre_lig_atom_id.rsd() << " and name " << pose.residue(pre_lig_atom_id.rsd()).atom_name(pre_lig_atom_id.atomno()) << " + resi "
			<< lig_atom_id.rsd() << " and name " << pose.residue(lig_atom_id.rsd()).atom_name(lig_atom_id.atomno()) << " + resi "
			<< zinc_id.rsd() << " and name " << pose.residue(zinc_id.rsd()).atom_name(zinc_id.atomno()) << std::endl;
	}


	//Tetrahedral
	for ( Size i(1); i <= tetrahedral_constraints_.size(); ++i ) {
		AtomID lig_atom_id_1( tetrahedral_constraints_[i]->atom(1) );
		AtomID       zinc_id( tetrahedral_constraints_[i]->atom(2) );
		AtomID lig_atom_id_2( tetrahedral_constraints_[i]->atom(3) );
		std::stringstream ss;
		ss << "tetr_" << lig_atom_id_1.rsd() << "_" << lig_atom_id_2.rsd();
		TR_PYMOL << "select " << ss.str() << ", resi "
			<< lig_atom_id_1.rsd() << " and name " << pose.residue(lig_atom_id_1.rsd()).atom_name(lig_atom_id_1.atomno()) << " + resi "
			<< zinc_id.rsd() << " and name " << pose.residue(zinc_id.rsd()).atom_name(zinc_id.atomno()) << " + resi "
			<< lig_atom_id_2.rsd() << " and name " << pose.residue(lig_atom_id_2.rsd()).atom_name(lig_atom_id_2.atomno()) << std::endl;

	}


	return;

}


void
AddZincSiteConstraints::output_constraints_file( core::pose::Pose const & pose ) {
	if ( distance_constraints_.size() == 0 ) {
		TR_PYMOL << "The constraints have not been added, no info to give." << std::endl;
		return;
	}

	utility::file::FileName filename( pose.pdb_info()->name() );

	utility::io::ozstream OutputConstraints;
	std::string zinc_cst_file_name = filename.base() + ".zinc_cst";
	OutputConstraints.open(zinc_cst_file_name);
	TR << "Zinc constraints will be written to a file called " << zinc_cst_file_name << std::endl;

	//Distances
	for ( Size i(1); i <= distance_constraints_.size(); ++i ) {
		AtomID lig_atom_id( distance_constraints_[i]->atom(1) );
		AtomID zinc_id    ( distance_constraints_[i]->atom(2) );
		const scoring::func::HarmonicFunc & dist_func = dynamic_cast< const scoring::func::HarmonicFunc& >( distance_constraints_[i]->get_func() );
		OutputConstraints << "AtomPair  "
			<< pose.residue(lig_atom_id.rsd()).atom_name(lig_atom_id.atomno()) << " " << lig_atom_id.rsd() << pose.pdb_info()->chain(lig_atom_id.rsd()) << " "
			<< pose.residue(zinc_id.rsd()).atom_name(zinc_id.atomno()) << " " << zinc_id.rsd() << pose.pdb_info()->chain(zinc_id.rsd()) << " "
			<< " \t\t\tHARMONIC " << dist_func.x0() << "  " << dist_func.sd() << std::endl;

	}

	//Angles
	for ( Size i(1); i <= angle_constraints_.size(); ++i ) {
		AtomID pre_lig_atom_id( angle_constraints_[i]->atom(1) );
		AtomID lig_atom_id    ( angle_constraints_[i]->atom(2) );
		AtomID zinc_id        ( angle_constraints_[i]->atom(3) );
		const scoring::func::HarmonicFunc & angle_func = dynamic_cast< const scoring::func::HarmonicFunc& >( angle_constraints_[i]->get_func() );
		OutputConstraints << "Angle     "
			<< pose.residue(pre_lig_atom_id.rsd()).atom_name(pre_lig_atom_id.atomno()) << " " << pre_lig_atom_id.rsd() << pose.pdb_info()->chain(pre_lig_atom_id.rsd()) << " "
			<< pose.residue(lig_atom_id.rsd()).atom_name(lig_atom_id.atomno()) << " " << lig_atom_id.rsd() << pose.pdb_info()->chain(lig_atom_id.rsd()) << " "
			<< pose.residue(zinc_id.rsd()).atom_name(zinc_id.atomno()) << " " << zinc_id.rsd() << pose.pdb_info()->chain(zinc_id.rsd()) << " "
			<< " \t\tHARMONIC " << angle_func.x0() << "  " << angle_func.sd() << std::endl; // radians or degrees

	}

	//Dihedrals
	for ( Size i(1); i <= dihedral_constraints_.size(); ++i ) {
		AtomID pre_pre_lig_atom_id( dihedral_constraints_[i]->atom(1) );
		AtomID pre_lig_atom_id    ( dihedral_constraints_[i]->atom(2) );
		AtomID lig_atom_id        ( dihedral_constraints_[i]->atom(3) );
		AtomID zinc_id            ( dihedral_constraints_[i]->atom(4) );
		const scoring::func::CircularHarmonicFunc & dihed_func = dynamic_cast< const scoring::func::CircularHarmonicFunc& >( dihedral_constraints_[i]->get_func() );
		OutputConstraints << "Dihedral  "
			<< pose.residue(pre_pre_lig_atom_id.rsd()).atom_name(pre_pre_lig_atom_id.atomno()) << " " << pre_pre_lig_atom_id.rsd() << pose.pdb_info()->chain(pre_pre_lig_atom_id.rsd()) << " "
			<< pose.residue(pre_lig_atom_id.rsd()).atom_name(pre_lig_atom_id.atomno()) << " " << pre_lig_atom_id.rsd() << pose.pdb_info()->chain(pre_lig_atom_id.rsd()) << " "
			<< pose.residue(lig_atom_id.rsd()).atom_name(lig_atom_id.atomno()) << " " << lig_atom_id.rsd() << pose.pdb_info()->chain(lig_atom_id.rsd()) << " "
			<< pose.residue(zinc_id.rsd()).atom_name(zinc_id.atomno()) << " " << zinc_id.rsd() << pose.pdb_info()->chain(zinc_id.rsd()) << " "
			<< " HARMONIC " << dihed_func.x0() << "  " << dihed_func.sd() << std::endl; // radians or degrees
	}


	//Tetrahedral
	for ( Size i(1); i <= tetrahedral_constraints_.size(); ++i ) {
		AtomID lig_atom_id_1( tetrahedral_constraints_[i]->atom(1) );
		AtomID       zinc_id( tetrahedral_constraints_[i]->atom(2) );
		AtomID lig_atom_id_2( tetrahedral_constraints_[i]->atom(3) );
		const scoring::func::HarmonicFunc & tetr_func = dynamic_cast< const scoring::func::HarmonicFunc& >( tetrahedral_constraints_[i]->get_func() );
		OutputConstraints << "Angle     "
			<< pose.residue(lig_atom_id_1.rsd()).atom_name(lig_atom_id_1.atomno()) << " " << lig_atom_id_1.rsd() << pose.pdb_info()->chain(lig_atom_id_1.rsd()) << " "
			<< pose.residue(zinc_id.rsd()).atom_name(zinc_id.atomno()) << " " << zinc_id.rsd() << pose.pdb_info()->chain(zinc_id.rsd()) << " "
			<< pose.residue(lig_atom_id_2.rsd()).atom_name(lig_atom_id_2.atomno()) << " " << lig_atom_id_2.rsd() << pose.pdb_info()->chain(lig_atom_id_2.rsd()) << " "
			<< " \t\tHARMONIC " << tetr_func.x0() << "  " << tetr_func.sd() << std::endl; // radians or degrees
	}


	return;

}


} //metal_interface
} //protocols
