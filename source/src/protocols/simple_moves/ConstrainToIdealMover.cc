// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/ConstrainToIdealMover.cc
/// @brief Adds to your pose constraints suitable for idealization of bond lengths and angles
/// @author Steven Lewis (smlewi@gmail.com); Rhiju Das


// Unit Headers
#include <protocols/simple_moves/ConstrainToIdealMover.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>

#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/rna/util.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>

#include <core/conformation/Residue.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>

#include <utility/vector1.hh>


/////////////////////////////////////////////////////////////////////////////////////////////////
//
// NOTES
//
//  -- This object is supposed to install constraints into the Pose, and update the MoveMap
//      to move bond lengths and angles of atoms that are free to move. Its actually not
//      clear to me if we need both MoveMap and AtomLevelDomainMap given -- either one would be OK?
//
//  -- This object's name should not be ConstrainToIdealMover.
//
//  -- Does this reproduce functionality already developed in cartbond framework (which was
//      set up later by Frank Dimaio but is more general?). If so let's get rid of this one.
//
//  -- note also Steven Lewis comments in ConstrainToIdealMover.hh function.
//
//  -- setup of dihedrals will be important for torsions that are not covered by torsion potential.
//      There's something bare bones below for polar hydrogens, out of necessity.
//
//      -- rd, 2014
/////////////////////////////////////////////////////////////////////////////////////////////////

using basic::T;
using basic::Error;
using basic::Warning;
using namespace basic::options;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.ConstrainToIdealMover" );

namespace protocols {
namespace simple_moves {

ConstrainToIdealMover::ConstrainToIdealMover() :
	protocols::moves::Mover( "ConstrainToIdealMover" ),
	atom_level_domain_map_( /*0*/ ), //requires pose to initialize
	bond_length_sd_( 0.05 ), // might be better to have in a database for each atom
	bond_length_sd_polar_hydrogen_( 0.20 ), // might be better to have in a database for each atom
	bond_angle_sd_( numeric::conversions::radians( 5.0 ) ), // might be better to have in a database for each atom
	bond_angle_sd_polar_hydrogen_( numeric::conversions::radians( option[ OptionKeys::score::bond_angle_sd_polar_hydrogen ]() ) ), // might be better to have in a database for each atom
	bond_torsion_sd_( numeric::conversions::radians( 5.0 ) ), // might be better to have in a database for each atom
	bond_torsion_sd_polar_hydrogen_( numeric::conversions::radians( option[ OptionKeys::score::bond_torsion_sd_polar_hydrogen ]() ) ), // might be better to have in a database for each atom
	score_type_( core::scoring::bond_geometry ),
	just_rna_backbone_( false ),
	just_polar_hydrogens_( false ),
	legacy_dof_allow_move_( false ), // old setting -- did not move all atoms that should have been moving.
	verbose_( false )
{}

ConstrainToIdealMover::~ConstrainToIdealMover(){}

// Is this really in use?
// ConstrainToIdealMover::ConstrainToIdealMover(ConstrainToIdealMover const & rhs) :
//  protocols::moves::Mover(rhs),

// {
//  *this = rhs;
//  return;
// }

// Is this really in use?
ConstrainToIdealMover & ConstrainToIdealMover::operator = ( ConstrainToIdealMover const & rhs ){

	//abort self-assignment
	if ( this == &rhs ) return *this;

	atom_level_domain_map_     = rhs.get_atom_level_domain_map()->clone();
	just_rna_backbone_ = rhs.just_rna_backbone_;
	just_polar_hydrogens_ = rhs.just_polar_hydrogens_;

	return *this;
}

std::string
ConstrainToIdealMover::get_name() const {
	return "ConstrainToIdealMover"; //ConstrainToIdealMoverCreator::mover_name();
}

protocols::moves::MoverOP ConstrainToIdealMover::fresh_instance() const { return protocols::moves::MoverOP( new ConstrainToIdealMover ); }
protocols::moves::MoverOP ConstrainToIdealMover::clone() const { return protocols::moves::MoverOP( new ConstrainToIdealMover( *this ) ); }

/// @brief setter for AtomLevelDomainMap; shallow copy
void ConstrainToIdealMover::set_atom_level_domain_map( protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map ) { atom_level_domain_map_ = atom_level_domain_map; }

/// @brief getter for AtomLevelDomainMap
protocols::toolbox::AtomLevelDomainMapCOP ConstrainToIdealMover::get_atom_level_domain_map() const { return atom_level_domain_map_; }

/// @details This code will modify your input pose by adding constraints which will trend bond lengths and angles towards ideal.
void ConstrainToIdealMover::apply( core::pose::Pose & pose ){
	//handle mm - if we don't have a freshly externally supplied one, throw the old one out.
	//this is not a clear() operation in case someone is still sharing the old pointer
	core::kinematics::MoveMap mm;
	apply( pose, mm );
	if ( !atom_level_domain_map_ ) atom_level_domain_map_ = protocols::toolbox::AtomLevelDomainMapOP( new protocols::toolbox::AtomLevelDomainMap( pose ) );
} //apply


/// @details This code will modify your input pose by adding constraints which will trend bond lengths and angles towards ideal.  If you input a movemap via set_movemap, that movemap will be modified to free the same set of bond lengths and angles (needs some testing).
void ConstrainToIdealMover::apply( core::pose::Pose & pose, core::kinematics::MoveMap & mm ){
	core::pose::Pose pose_reference;
	create_pose_reference( pose, pose_reference );
	vary_bond_geometry( pose, mm, pose_reference );
} //apply

//////////////////////////////////START MOVED CODE/////////////////////////////
//Most of the code for this mover was moved from Rhiju's RNA_Minimizer.  It was
//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
///////////////////////////////////////////////////////////////////////////////
bool
ConstrainToIdealMover::check_in_bonded_list(
	core::id::AtomID const & atom_id1,
	core::id::AtomID const & atom_id2,
	utility::vector1< std::pair< core::id::AtomID, core::id::AtomID > > & bonded_atom_list ) const
{

	for ( core::Size n = 1; n <= bonded_atom_list.size(); n++ ) {
		if ( atom_id1 == bonded_atom_list[ n ].first && atom_id2 == bonded_atom_list[ n ].second ) return true;
		if ( atom_id2 == bonded_atom_list[ n ].first && atom_id1 == bonded_atom_list[ n ].second ) return true;
	}
	return false;
}

bool
ConstrainToIdealMover::check_in_bond_angle_list(
	core::id::AtomID const & atom_id1,
	core::id::AtomID const & atom_id2,
	core::id::AtomID const & atom_id3,
	utility::vector1< std::pair< core::id::AtomID, std::pair< core::id::AtomID, core::id::AtomID > > > & bond_angle_list ) const
{

	for ( core::Size n = 1; n <= bond_angle_list.size(); n++ ) {
		if ( atom_id1 == bond_angle_list[ n ].first ) {
			if ( atom_id2 == bond_angle_list[ n ].second.first && atom_id3 == bond_angle_list[ n ].second.second ) return true;
			if ( atom_id3 == bond_angle_list[ n ].second.first && atom_id2 == bond_angle_list[ n ].second.second ) return true;
		}
	}
	return false;
}

void
ConstrainToIdealMover::add_bond_length_constraint(
	core::id::AtomID const & atom_id1,
	core::id::AtomID const & atom_id2,
	utility::vector1< std::pair< core::id::AtomID, core::id::AtomID > > & bonded_atom_list,
	core::pose::Pose const & pose,
	core::pose::Pose const & pose_reference,
	core::scoring::constraints::ConstraintSetOP & cst_set ) const
{

	using namespace core::scoring;
	using namespace core::scoring::constraints;

	std::string const & atom_name1 = pose.residue( atom_id1.rsd() ).atom_name( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue( atom_id2.rsd() ).atom_name( atom_id2.atomno() );

	// if ( !pose_reference.residue( atom_id1.rsd() ).has( atom_name1 ) ) return;
	// if ( !pose_reference.residue( atom_id2.rsd() ).has( atom_name2 ) ) return;

	runtime_assert( pose_reference.residue( atom_id1.rsd() ).has( atom_name1 ) );
	runtime_assert( pose_reference.residue( atom_id2.rsd() ).has( atom_name2 ) );

	if ( pose.residue( atom_id1.rsd() ).is_virtual( atom_id1.atomno() ) ) return;
	if ( pose.residue( atom_id2.rsd() ).is_virtual( atom_id2.atomno() ) ) return;

	if ( !check_in_bonded_list( atom_id1, atom_id2,  bonded_atom_list ) ) {
		bonded_atom_list.push_back( std::make_pair( atom_id1, atom_id2 ) );

		core::Real const bond_length = ( pose_reference.residue( atom_id1.rsd() ).xyz( atom_name1 ) -
			pose_reference.residue( atom_id2.rsd() ).xyz( atom_name2 ) ).length();

		core::Real bond_length_sd( bond_length_sd_ );
		if ( pose.residue_type( atom_id1.rsd() ).atom_type( atom_id1.atomno() ).is_polar_hydrogen() ||
				pose.residue_type( atom_id2.rsd() ).atom_type( atom_id2.atomno() ).is_polar_hydrogen() ) {
			bond_length_sd = bond_length_sd_polar_hydrogen_;
		}

		core::scoring::func::FuncOP dist_harm_func_( new core::scoring::func::HarmonicFunc( bond_length, bond_length_sd ) );

		cst_set->add_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint( atom_id1,
			atom_id2,
			dist_harm_func_,
			score_type_ ) ) ) );
		if ( verbose_ ) {
			TR << "PUTTING CONSTRAINT ON DISTANCE: " <<
				atom_id2.rsd() << " " << atom_name1 << "; "  <<
				atom_id1.rsd() << " " << atom_name2 << " ideal length "  <<
				bond_length <<  " compared to current length " <<
				( pose.xyz( atom_id1 ) - pose.xyz( atom_id2 ) ).length() <<
				std::endl;
		}
	}

}

void
ConstrainToIdealMover::add_bond_angle_constraint(
	core::id::AtomID const & atom_id1,
	core::id::AtomID const & atom_id2,
	core::id::AtomID const & atom_id3,
	utility::vector1< std::pair < core::id::AtomID, std::pair< core::id::AtomID, core::id::AtomID > > > & bond_angle_list,
	core::pose::Pose const & pose,
	core::pose::Pose const & pose_reference,
	core::scoring::constraints::ConstraintSetOP & cst_set ) const
{

	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace numeric::conversions;

	if ( atom_id2 == atom_id3 ) return;

	std::string const & atom_name1 = pose.residue( atom_id1.rsd() ).atom_name( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue( atom_id2.rsd() ).atom_name( atom_id2.atomno() );
	std::string const & atom_name3 = pose.residue( atom_id3.rsd() ).atom_name( atom_id3.atomno() );

	// if ( !pose_reference.residue( atom_id1.rsd() ).has( atom_name1 ) ) return;
	// if ( !pose_reference.residue( atom_id2.rsd() ).has( atom_name2 ) ) return;
	// if ( !pose_reference.residue( atom_id3.rsd() ).has( atom_name3 ) ) return;

	runtime_assert( pose_reference.residue( atom_id1.rsd() ).has( atom_name1 ) );
	runtime_assert( pose_reference.residue( atom_id2.rsd() ).has( atom_name2 ) );
	runtime_assert( pose_reference.residue( atom_id3.rsd() ).has( atom_name3 ) );

	if ( pose.residue( atom_id1.rsd() ).is_virtual( atom_id1.atomno() ) ) return;
	if ( pose.residue( atom_id2.rsd() ).is_virtual( atom_id2.atomno() ) ) return;
	if ( pose.residue( atom_id3.rsd() ).is_virtual( atom_id3.atomno() ) ) return;

	if ( !check_in_bond_angle_list( atom_id1, atom_id2, atom_id3, bond_angle_list ) ) {
		bond_angle_list.push_back( std::make_pair( atom_id1, std::make_pair( atom_id2, atom_id3 ) ) );

		core::Real const bond_angle = angle_radians(
			pose_reference.residue( atom_id2.rsd() ).xyz( atom_name2 ) ,
			pose_reference.residue( atom_id1.rsd() ).xyz( atom_name1 ) ,
			pose_reference.residue( atom_id3.rsd() ).xyz( atom_name3 )
		);

		runtime_assert( bond_angle >= 0.001 );

		core::Real bond_angle_sd( bond_angle_sd_ );
		if ( pose.residue_type( atom_id2.rsd() ).atom_type( atom_id2.atomno() ).is_polar_hydrogen() ||
				pose.residue_type( atom_id3.rsd() ).atom_type( atom_id3.atomno() ).is_polar_hydrogen() ) {
			bond_angle_sd = bond_angle_sd_polar_hydrogen_;
		}

		core::scoring::func::FuncOP angle_harm_func_( new core::scoring::func::HarmonicFunc( bond_angle, bond_angle_sd ) );
		cst_set->add_constraint(ConstraintCOP( ConstraintOP( new AngleConstraint(
			atom_id2, atom_id1, atom_id3, angle_harm_func_, score_type_ ) ) ) );

		if ( verbose_ ) {
			TR << "PUTTING CONSTRAINT ON ANGLE: " <<
				atom_id2.rsd() << " " << pose_reference.residue( atom_id2.rsd() ).atom_name( atom_id2.atomno() ) << "; "  <<
				atom_id1.rsd() << " " << pose_reference.residue( atom_id1.rsd() ).atom_name( atom_id1.atomno() ) << "; "  <<
				atom_id3.rsd() << " " << pose_reference.residue( atom_id3.rsd() ).atom_name( atom_id3.atomno() ) << " ==> "  << degrees( bond_angle ) << " " << degrees( bond_angle_sd ) <<
				std::endl;
		}

	}

}

// This does not check for redundancy yet -- would be easy to put in.
void
ConstrainToIdealMover::add_bond_dihedral_constraint(
	core::id::AtomID const & atom_id1,
	core::id::AtomID const & atom_id2,
	core::id::AtomID const & atom_id3,
	core::id::AtomID const & atom_id4,
	core::pose::Pose const & pose,
	core::pose::Pose const & pose_reference,
	core::scoring::constraints::ConstraintSetOP & cst_set ) const
{

	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace numeric::conversions;

	std::string const & atom_name1 = pose.residue( atom_id1.rsd() ).atom_name( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue( atom_id2.rsd() ).atom_name( atom_id2.atomno() );
	std::string const & atom_name3 = pose.residue( atom_id3.rsd() ).atom_name( atom_id3.atomno() );
	std::string const & atom_name4 = pose.residue( atom_id4.rsd() ).atom_name( atom_id4.atomno() );

	// if ( !pose_reference.residue( atom_id1.rsd() ).has( atom_name1 ) ) return;
	// if ( !pose_reference.residue( atom_id2.rsd() ).has( atom_name2 ) ) return;
	// if ( !pose_reference.residue( atom_id3.rsd() ).has( atom_name3 ) ) return;
	// if ( !pose_reference.residue( atom_id4.rsd() ).has( atom_name4 ) ) return;

	runtime_assert( pose_reference.residue( atom_id1.rsd() ).has( atom_name1 ) );
	runtime_assert( pose_reference.residue( atom_id2.rsd() ).has( atom_name2 ) );
	runtime_assert( pose_reference.residue( atom_id3.rsd() ).has( atom_name3 ) );
	runtime_assert( pose_reference.residue( atom_id4.rsd() ).has( atom_name4 ) );

	if ( pose.residue( atom_id1.rsd() ).is_virtual( atom_id1.atomno() ) ) return;
	if ( pose.residue( atom_id2.rsd() ).is_virtual( atom_id2.atomno() ) ) return;
	if ( pose.residue( atom_id3.rsd() ).is_virtual( atom_id3.atomno() ) ) return;
	if ( pose.residue( atom_id4.rsd() ).is_virtual( atom_id4.atomno() ) ) return;

	// if ( !check_in_bond_angle_list( atom_id1, atom_id2, atom_id3, bond_angle_list ) ) {
	//  bond_angle_list.push_back( std::make_pair( atom_id1, std::make_pair( atom_id2, atom_id3 ) ) );

	core::Real const bond_torsion = dihedral_radians(
		pose_reference.residue( atom_id1.rsd() ).xyz( atom_name1 ) ,
		pose_reference.residue( atom_id2.rsd() ).xyz( atom_name2 ) ,
		pose_reference.residue( atom_id3.rsd() ).xyz( atom_name3 ) ,
		pose_reference.residue( atom_id4.rsd() ).xyz( atom_name4 )
	);

	core::Real bond_torsion_sd( bond_torsion_sd_ );
	if ( pose.residue_type( atom_id1.rsd() ).atom_type( atom_id1.atomno() ).is_polar_hydrogen() ||
			pose.residue_type( atom_id4.rsd() ).atom_type( atom_id4.atomno() ).is_polar_hydrogen() ) {
		bond_torsion_sd = bond_torsion_sd_polar_hydrogen_;
	}

	core::scoring::func::FuncOP angle_harm_func_( new core::scoring::func::CircularHarmonicFunc( bond_torsion, bond_torsion_sd ) );
	cst_set->add_constraint( DihedralConstraintOP( new DihedralConstraint(
		atom_id1, atom_id2, atom_id3, atom_id4, angle_harm_func_, score_type_ ) ) );

	if ( verbose_ ) {
		TR << "PUTTING CONSTRAINT ON TORSION: " <<
			atom_id2.rsd() << " " << pose_reference.residue( atom_id1.rsd() ).atom_name( atom_id1.atomno() ) << "; "  <<
			atom_id1.rsd() << " " << pose_reference.residue( atom_id2.rsd() ).atom_name( atom_id2.atomno() ) << "; "  <<
			atom_id1.rsd() << " " << pose_reference.residue( atom_id3.rsd() ).atom_name( atom_id3.atomno() ) << "; "  <<
			atom_id3.rsd() << " " << pose_reference.residue( atom_id4.rsd() ).atom_name( atom_id4.atomno() ) << " ==> "  << degrees( bond_torsion ) << " " << degrees( bond_torsion_sd ) <<
			std::endl;
	}


}

//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
bool
ConstrainToIdealMover::check_if_really_connected(
	core::pose::Pose const & pose,
	core::id::AtomID const & atom_id1,
	core::id::AtomID const & atom_id2 ) const
{
	if (  atom_id1.rsd() == atom_id2.rsd() ) return true;

	core::kinematics::tree::AtomCOP atom1 ( pose.atom_tree().atom( atom_id1 ).get_self_ptr() );
	core::kinematics::tree::AtomCOP atom2 ( pose.atom_tree().atom( atom_id2 ).get_self_ptr() );

	if ( atom1->parent() == atom2 ) return true;
	if ( atom2->parent() == atom1 ) return true;

	return false;
}

//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
// AMW TODO: adapt version from erraser_minimizer plus figure out how to incorporate the atom_level_domain_map bit...
bool
ConstrainToIdealMover::i_want_this_atom_to_move( core::conformation::Residue const & residue2, core::Size const & k ) const
{

	if ( just_rna_backbone_ ) { //default
		if ( !residue2.is_RNA() ) return false;
		if ( k > residue2.first_sidechain_atom() &&
				k != core::chemical::rna::first_base_atom_index( residue2.type() ) ) return false;
	}

	if ( just_polar_hydrogens_ ) {
		if ( !residue2.Hpos_polar().has_value( k ) )  return false;
	}

	core::id::AtomID id( k, residue2.seqpos() );
	if ( !atom_level_domain_map_->get( id ) ) return false;

	if ( residue2.is_virtual( k ) ) return false;

	return true;

}

//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
bool
ConstrainToIdealMover::i_want_this_atom_to_move( core::pose::Pose const & pose, core::id::AtomID const & atom_id ) const
{
	core::conformation::Residue const & residue( pose.residue( atom_id.rsd() ) );
	core::Size const & k( atom_id.atomno() );
	return i_want_this_atom_to_move( residue, k );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
core::Size
check_if_proton_chi_atom( core::pose::Pose const & pose, core::Size const i, core::Size const j ){
	core::chemical::ResidueType rsd_type = pose.residue_type( i );
	for ( core::Size n = 1; n <= rsd_type.n_proton_chi(); n++ ) {
		core::Size chino = rsd_type.proton_chi_2_chi( n );
		core::Size const & proton_chi_atom = rsd_type.chi_atoms( chino )[4];
		if ( proton_chi_atom == j ) return n;
	}
	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Following has not (yet) been carefully debugged. --Rhiju
//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
void
ConstrainToIdealMover::vary_bond_geometry(
	core::pose::Pose & pose,
	core::kinematics::MoveMap & mm,
	core::pose::Pose const & pose_reference ) const {

	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace numeric::conversions;

	ConstraintSetOP cst_set = pose.constraint_set()->clone();
	pose.constraint_set( cst_set );

	core::Size const nres( pose.total_residue() );

	//std::map< AtomID, utility::vector1< AtomID > > lists_of_angle_bonded_atoms;

	for  ( core::Size i = 1; i <= nres; i++ )  {

		if ( !atom_level_domain_map_->get( i ) ) continue;

		core::conformation::Residue const & residue( pose.residue( i )  );

		for ( core::Size j = 1; j <= residue.natoms(); j++ ) {

			//   TR << "checking: res " << i << " atom " << j << " " << residue.atom_name( j ) << " " << i_want_this_atom_to_move( residue, j ) << std::endl;
			if ( !i_want_this_atom_to_move( residue, j ) ) continue;

			core::kinematics::tree::AtomCOP current_atom ( pose.atom_tree().atom( AtomID(j,i) ).get_self_ptr() );
			if ( current_atom->is_jump() ) continue;

			///////////////////
			core::kinematics::tree::AtomCOP input_stub_atom1( current_atom->input_stub_atom1() );
			if ( !input_stub_atom1 ) continue;
			if ( legacy_dof_allow_move_ && !i_want_this_atom_to_move( pose, input_stub_atom1->id() ) ) continue;
			mm.set( DOF_ID( AtomID( j, i ), D ), true );
			if ( input_stub_atom1->is_jump() ) continue;

			core::kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );

			///////////////////
			if ( !input_stub_atom2 ) continue;
			if ( input_stub_atom2 == current_atom ) continue;
			if ( legacy_dof_allow_move_ && !i_want_this_atom_to_move( pose, input_stub_atom2->id() ) ) continue;
			mm.set( DOF_ID( AtomID( j, i ), THETA ), true );
			if ( input_stub_atom2->is_jump() ) continue;

			///////////////////
			core::kinematics::tree::AtomCOP input_stub_atom3( current_atom->input_stub_atom3() );
			if ( !input_stub_atom3 ) continue;
			if ( legacy_dof_allow_move_ && !i_want_this_atom_to_move( pose, input_stub_atom3->id() ) ) continue;
			if ( input_stub_atom3 == current_atom ) continue;
			mm.set( DOF_ID( AtomID( j, i ), PHI ), true );

		}
	}

	utility::vector1< std::pair< AtomID, AtomID > >  bond_list;
	utility::vector1< std::pair< AtomID, std::pair< AtomID, AtomID > > > bond_angle_list;

	for  ( core::Size i = 1; i <= nres; i++ )  {

		//Go through all bonds in pose...
		if ( !atom_level_domain_map_->get( i ) ) continue;

		core::conformation::Residue const & residue( pose.residue( i )  );

		for ( core::Size j = 1; j <= residue.natoms(); j++ ) {

			if ( legacy_dof_allow_move_ && !i_want_this_atom_to_move( residue, j ) ) continue;

			AtomID atom_id1( j, i );

			utility::vector1< AtomID >  nbrs( pose.conformation().bonded_neighbor_all_res( atom_id1 ) );

			// Bond lengths.
			for ( core::Size n = 1; n <= nbrs.size(); n++ ) {

				AtomID const & atom_id2( nbrs[ n ] );

				core::conformation::Residue const & residue2( pose.residue( atom_id2.rsd() ) ) ;
				core::Size const & k( atom_id2.atomno() ) ;

				if ( ! check_if_really_connected( pose, atom_id1, atom_id2 ) ) continue;

				if ( !i_want_this_atom_to_move( residue2, k ) ) continue;
				add_bond_length_constraint( atom_id1, atom_id2,
					bond_list,
					pose, pose_reference, cst_set );

			}

			// Bond angles
			for ( core::Size m = 1; m <= nbrs.size(); m++ ) {
				AtomID const & atom_id2( nbrs[ m ] );

				core::conformation::Residue const & residue2( pose.residue( atom_id2.rsd() ) ) ;
				if ( ! check_if_really_connected( pose, atom_id1, atom_id2 ) ) continue;

				core::Size const & k( atom_id2.atomno() ) ;

				for ( core::Size n = 1; n <= nbrs.size(); n++ ) {

					if ( m == n ) continue;
					AtomID const & atom_id3( nbrs[ n ] );

					core::conformation::Residue const & residue3( pose.residue( atom_id3.rsd() ) ) ;
					if ( ! check_if_really_connected( pose, atom_id1, atom_id3 ) ) continue;

					core::Size const & q( atom_id3.atomno() ) ;

					if ( !i_want_this_atom_to_move( residue2, k ) && !i_want_this_atom_to_move( residue3, q ) ) continue;
					if ( legacy_dof_allow_move_ &&  !( i_want_this_atom_to_move( residue2, k ) &&  i_want_this_atom_to_move( residue3, q ) ) ) continue;

					add_bond_angle_constraint( atom_id1, atom_id2, atom_id3, bond_angle_list,
						pose, pose_reference, cst_set );

				}
			}
		}
	}

	// following is not general -- need improper dihedrals for hydrogens that are not fixed.
	// could be generalized pretty easily -- look at all possible quartets and decide if they
	// are moving or not (like in movemap) -- just need to be smart about not double-counting
	// in rings. This may already be taken care of in cart_bonded...
	if ( just_polar_hydrogens_ ) {

		for  ( core::Size i = 1; i <= nres; i++ )  {

			if ( !atom_level_domain_map_->get( i ) ) continue;

			core::conformation::Residue const & residue( pose.residue( i )  );

			for ( core::Size j = 1; j <= residue.natoms(); j++ ) {

				//   TR << "checking: res " << i << " atom " << j << " " << residue.atom_name( j ) << " " << i_want_this_atom_to_move( residue, j ) << std::endl;

				if ( !i_want_this_atom_to_move( residue, j ) ) continue;

				core::kinematics::tree::AtomCOP current_atom ( pose.atom_tree().atom( AtomID( j, i ) ).get_self_ptr() );
				if ( current_atom->is_jump() ) continue;

				///////////////////
				core::kinematics::tree::AtomCOP input_stub_atom1( current_atom->input_stub_atom1() );
				if ( !input_stub_atom1 ) continue;
				if ( legacy_dof_allow_move_ && !i_want_this_atom_to_move( pose, input_stub_atom1->id() ) ) continue;
				if ( input_stub_atom1->is_jump() ) continue;

				core::kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );
				if ( !input_stub_atom2 ) continue;
				if ( input_stub_atom2 == current_atom ) continue;
				if ( legacy_dof_allow_move_ && !i_want_this_atom_to_move( pose, input_stub_atom2->id() ) ) continue;
				if ( input_stub_atom2->is_jump() ) continue;

				///////////////////
				core::kinematics::tree::AtomCOP input_stub_atom3( current_atom->input_stub_atom3() );
				if ( !input_stub_atom3 ) continue;
				if ( legacy_dof_allow_move_ && !i_want_this_atom_to_move( pose, input_stub_atom3->id() ) ) continue;
				if ( input_stub_atom3 == current_atom ) continue;

				// 'sister' atom
				// TR << pose.residue( i ).name() << " " << pose.residue( i ).atom_name( j ) << " " << ( input_stub_atom3->input_stub_atom1() == current_atom->input_stub_atom1() ) << " "
				//   << atom_id_to_named_atom_id( input_stub_atom3->input_stub_atom1()->id(), pose  ) << " "
				//   <<  atom_id_to_named_atom_id( current_atom->input_stub_atom1()->id(), pose  )
				//   << " " << check_if_proton_chi_atom( pose, i, j ) << std::endl;
				// if ( input_stub_atom3->input_stub_atom1() == current_atom->input_stub_atom1() ) continue;

				if ( check_if_proton_chi_atom( pose, i, j ) > 0 ) continue; // no phi constraints.

				add_bond_dihedral_constraint( current_atom->id(), input_stub_atom1->id(), input_stub_atom2->id(), input_stub_atom3->id(),
					pose, pose_reference, cst_set );


			}
		}

	}

	pose.constraint_set( cst_set );

}


///////////////////////////////////////////////////
//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
void
apply_ideal_coordinates_for_alternative_pucker( core::pose::Pose const & pose, core::pose::Pose & pose_reference )
{

	using namespace core::chemical::rna;
	using namespace core::pose::rna;
	RNA_FittedTorsionInfo const rna_fitted_torsion_info;

	core::Real const DELTA_CUTOFF( rna_fitted_torsion_info.delta_cutoff() );
	bool const use_phenix_geo = option[ OptionKeys::rna::corrected_geo ]();

	for ( core::Size n = 1; n <= pose.total_residue(); n++ ) {
		if ( !pose.residue_type( n ).is_RNA() ) continue;
		core::Real const delta = pose.residue( n ).mainchain_torsion( DELTA );
		if ( use_phenix_geo ) {
			apply_pucker( pose_reference, n,
				( delta < DELTA_CUTOFF ) ? NORTH : SOUTH,
				false /*skip_same_state*/,
				true /*idealize_coord*/ );
		} else {
			if ( delta > DELTA_CUTOFF ) core::pose::rna::apply_ideal_c2endo_sugar_coords( pose_reference, n );
		}
	}

	// pose_reference.dump_pdb( "REFERENCE.pdb" );
}

///////////////////////////////////////////////////
//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
void
ConstrainToIdealMover::create_pose_reference(
	core::pose::Pose const & pose,
	core::pose::Pose & pose_reference )
{
	using namespace core::chemical;
	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	make_pose_from_sequence( pose_reference, pose.annotated_sequence(), *rsd_set );
	apply_ideal_coordinates_for_alternative_pucker( pose, pose_reference );
}

void
ConstrainToIdealMover::set_score_type( core::scoring::ScoreType const setting ){ score_type_ = setting; }

///////////////////////////////////////////////////////////////////////////////////
// Let additional degrees of freedom vary -- but apply constraints to stay near
// ideal bond lengths and angles! This came from RNA_Minimizer.
void
setup_vary_rna_bond_geometry( core::kinematics::MoveMap & mm,
	core::pose::Pose & pose, toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
	core::scoring::ScoreType score_type /* = rna_bond_geometry*/ ) {
	ConstrainToIdealMover CTIMover;
	CTIMover.set_atom_level_domain_map( atom_level_domain_map );
	CTIMover.set_just_rna_backbone( true );
	CTIMover.set_score_type( score_type );
	CTIMover.apply( pose, mm );
}


///////////////////////////////////////////////////////////////////////////////////
// Let additional degrees of freedom vary -- but apply constraints to stay near
// ideal bond lengths and angles! This came from RNA_Minimizer.
void
setup_vary_polar_hydrogen_geometry( core::kinematics::MoveMap & mm,
	core::pose::Pose & pose, toolbox::AtomLevelDomainMapCOP atom_level_domain_map ) {
	ConstrainToIdealMover CTIMover;
	CTIMover.set_atom_level_domain_map( atom_level_domain_map );
	CTIMover.set_just_polar_hydrogens( true );
	CTIMover.set_score_type( core::scoring::bond_geometry );
	CTIMover.apply( pose, mm );
}


} // namespace moves
} // namespace protocols
