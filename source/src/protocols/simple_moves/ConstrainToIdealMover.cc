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
//#include <protocols/simple_moves/ConstrainToIdealMoverCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <protocols/toolbox/AllowInsert.hh> //need to move AllowInsert to toolbox XRW2

#include <core/conformation/Residue.hh>

#include <core/chemical/ChemicalManager.hh>
//#include <core/chemical/util.hh>

//#include <core/id/AtomID_Map.hh>
//TRANSCLUDED in MoveMap #include <core/id/AtomID.hh>
//#include <core/id/NamedAtomID.hh>
//TRANSCLUDED in MoveMap #include <core/id/DOF_ID.hh>
//#include <core/id/TorsionID.hh>

#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/MoveMap.hh>

//#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

// Utility Headers
#include <basic/Tracer.hh>
//#include <core/types.hh>
//#include <utility/exit.hh>

#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>

// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.simple_moves.ConstrainToIdealMover" );

namespace protocols {
namespace simple_moves {

/*std::string
ConstrainToIdealMoverCreator::keyname() const
{
	return ConstrainToIdealMoverCreator::mover_name();
}

protocols::moves::MoverOP
ConstrainToIdealMoverCreator::create_mover() const {
	return new ConstrainToIdealMover;
}

std::string
ConstrainToIdealMoverCreator::mover_name()
{
	return "ConstrainToIdealMover";
}
*/

///@brief parse XML (specifically in the context of the parser/scripting scheme)
// void
// ConstrainToIdealMover::parse_my_tag(
// 	TagCOP const tag,
// 	basic::datacache::DataMap & datamap,
// 	Filters_map const & filters,
// 	protocols::moves::Movers_map const & movers,
// 	Pose const & pose
// )
// {
// 	if ( tag->getName() != "ConstrainToIdealMover" ) {
// 		TR << " received incompatible Tag " << tag << std::endl;
// 		assert(false);
// 		return;
// 	}
//}

ConstrainToIdealMover::ConstrainToIdealMover()
	: protocols::moves::Mover("ConstrainToIdealMover"),
		allow_insert_(NULL), //requires pose to initialize
		mm_(new core::kinematics::MoveMap()),
		supplied_movemap_(false)
{}

ConstrainToIdealMover::~ConstrainToIdealMover(){}

ConstrainToIdealMover::ConstrainToIdealMover(ConstrainToIdealMover const & rhs) : protocols::moves::Mover(rhs) {
	*this = rhs;
	return;
}

ConstrainToIdealMover & ConstrainToIdealMover::operator=( ConstrainToIdealMover const & rhs ){

	//abort self-assignment
	if (this == &rhs) return *this;

	//	allow_insert_ = rhs.get_AllowInsert(); //would prefer to clone, but that class author offered no way to do it
	allow_insert_ = rhs.get_AllowInsert()->clone(); //OK stephen, I did it. --Rhiju
	mm_ = rhs.get_movemap()->clone();
	supplied_movemap_ = rhs.supplied_movemap_;
	return *this;
}

std::string
ConstrainToIdealMover::get_name() const {
	return "ConstrainToIdealMover"; //ConstrainToIdealMoverCreator::mover_name();
}

protocols::moves::MoverOP ConstrainToIdealMover::fresh_instance() const { return new ConstrainToIdealMover; }
protocols::moves::MoverOP ConstrainToIdealMover::clone() const { return new ConstrainToIdealMover( *this ); }

///@details supply a MoveMap to be further modified by this mover; otherwise it will modify a blank one; notice it sets a boolean flag to mark the source of the movemap
void ConstrainToIdealMover::set_movemap( core::kinematics::MoveMapOP movemap ) {
	mm_ = movemap;
	supplied_movemap_ = true;
}

///@brief get the MoveMap last created during apply(); it has the bond and angles free for minimization
core::kinematics::MoveMapOP ConstrainToIdealMover::get_movemap() const {return mm_;}

///@brief setter for AllowInsert; shallow copy
void ConstrainToIdealMover::set_AllowInsert(protocols::toolbox::AllowInsertOP allow_insert) {allow_insert_ = allow_insert;}

///@brief getter for AllowInsert
protocols::toolbox::AllowInsertOP ConstrainToIdealMover::get_AllowInsert() const {return allow_insert_;}

///@details This code will modify your input pose by adding constraints which will trend bond lengths and angles towards ideal.  If you input a movemap via set_movemap, that movemap will be modified to free the same set of bond lengths and angles (needs some testing) and will be available from get_movemap.
void ConstrainToIdealMover::apply( core::pose::Pose & pose ){
  using core::kinematics::MoveMapOP;
  using core::kinematics::MoveMap;

	//handle mm - if we don't have a freshly externally supplied one, throw the old one out.
	//this is not a clear() operation in case someone is still sharing the old pointer
	if(!supplied_movemap_ || !mm_) mm_ = new core::kinematics::MoveMap();
	if(!allow_insert_) allow_insert_ = new protocols::toolbox::AllowInsert(pose);

	core::pose::Pose pose_reference;
	create_pose_reference(pose, pose_reference);

	vary_bond_geometry(pose, pose_reference);

	//the movemap in mm_ is now modified and cannot be reused
	supplied_movemap_ = false;
	return;
}//apply

//////////////////////////////////START MOVED CODE/////////////////////////////
//Most of the code for this mover was moved from Rhiju's RNA_Minimizer.  It was
//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
///////////////////////////////////////////////////////////////////////////////

//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
bool
ConstrainToIdealMover::check_in_bonded_list(
	core::id::AtomID const & atom_id1,
	core::id::AtomID const & atom_id2,
	utility::vector1< std::pair< core::id::AtomID, core::id::AtomID > > & bonded_atom_list ) {

	for (core::Size n = 1; n <= bonded_atom_list.size(); n++ ) {
		if( atom_id1 == bonded_atom_list[ n ].first && atom_id2 == bonded_atom_list[ n ].second ) return true;
		if( atom_id2 == bonded_atom_list[ n ].first && atom_id1 == bonded_atom_list[ n ].second ) return true;
	}
	return false;
}

//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
bool
ConstrainToIdealMover::check_in_bond_angle_list(
	core::id::AtomID const & atom_id1,
	core::id::AtomID const & atom_id2,
	core::id::AtomID const & atom_id3,
	utility::vector1< std::pair< core::id::AtomID, std::pair< core::id::AtomID, core::id::AtomID > > > & bond_angle_list ) {

	for (core::Size n = 1; n <= bond_angle_list.size(); n++ ) {
		if( atom_id1 == bond_angle_list[ n ].first ) {
			if( atom_id2 == bond_angle_list[ n ].second.first && atom_id3 == bond_angle_list[ n ].second.second ) return true;
			if( atom_id3 == bond_angle_list[ n ].second.first && atom_id2 == bond_angle_list[ n ].second.second ) return true;
		}
	}
	return false;
}

//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
void
ConstrainToIdealMover::add_bond_constraint(
	core::id::AtomID const & atom_id1,
	core::id::AtomID const & atom_id2,
	utility::vector1< std::pair< core::id::AtomID, core::id::AtomID > > & bonded_atom_list,
	core::pose::Pose const & pose,
	core::pose::Pose const & pose_reference,
	core::scoring::constraints::ConstraintSetOP & cst_set )
{

	using namespace core::scoring;
	using namespace core::scoring::constraints;

	std::string const & atom_name1 = pose.residue( atom_id1.rsd() ).atom_name( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue( atom_id2.rsd() ).atom_name( atom_id2.atomno() );

	if ( !pose_reference.residue( atom_id1.rsd() ).has( atom_name1 ) ) return;
	if ( !pose_reference.residue( atom_id2.rsd() ).has( atom_name2 ) ) return;

	if ( !check_in_bonded_list( atom_id1, atom_id2,  bonded_atom_list ) ) {
		bonded_atom_list.push_back( std::make_pair( atom_id1, atom_id2 ) );

		core::Real const bond_length_sd_( 0.05 );

		core::Real const bond_length = ( pose_reference.residue( atom_id1.rsd() ).xyz( atom_name1 ) -
															 pose_reference.residue( atom_id2.rsd() ).xyz( atom_name2 ) ).length();


		core::scoring::func::FuncOP dist_harm_func_( new core::scoring::func::HarmonicFunc( bond_length, bond_length_sd_ ));

		cst_set->add_constraint( new AtomPairConstraint( atom_id1 ,
																										 atom_id2,
																										 dist_harm_func_,
																										 rna_bond_geometry ) );
		if ( false ) {
			TR << "PUTTING CONSTRAINT ON DISTANCE: " <<
				atom_id2.rsd() << " " << atom_name1 << "; "  <<
				atom_id1.rsd() << " " << atom_name2 << " "  <<
				bond_length <<
				std::endl;
		}
	}

}


//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
void
ConstrainToIdealMover::add_bond_angle_constraint(
	core::id::AtomID const & atom_id1,
	core::id::AtomID const & atom_id2,
	core::id::AtomID const & atom_id3,
	utility::vector1< std::pair < core::id::AtomID, std::pair< core::id::AtomID, core::id::AtomID > > > & bond_angle_list,
	core::pose::Pose const & pose,
	core::pose::Pose const & pose_reference,
	core::scoring::constraints::ConstraintSetOP & cst_set )
{

	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace numeric::conversions;

	if (atom_id2 == atom_id3) return;

	std::string const & atom_name1 = pose.residue( atom_id1.rsd() ).atom_name( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue( atom_id2.rsd() ).atom_name( atom_id2.atomno() );
	std::string const & atom_name3 = pose.residue( atom_id3.rsd() ).atom_name( atom_id3.atomno() );

	if ( !pose_reference.residue( atom_id1.rsd() ).has( atom_name1 ) ) return;
	if ( !pose_reference.residue( atom_id2.rsd() ).has( atom_name2 ) ) return;
	if ( !pose_reference.residue( atom_id3.rsd() ).has( atom_name3 ) ) return;

	if ( !check_in_bond_angle_list( atom_id1, atom_id2, atom_id3, bond_angle_list ) ) {
		bond_angle_list.push_back( std::make_pair( atom_id1, std::make_pair( atom_id2, atom_id3 ) ) );


		core::Real const bond_angle_sd_( radians( 5.0 ) );

		core::Real const bond_angle = angle_radians(
																								 pose_reference.residue( atom_id2.rsd() ).xyz( atom_name2 )  ,
																								 pose_reference.residue( atom_id1.rsd() ).xyz( atom_name1 )  ,
																								 pose_reference.residue( atom_id3.rsd() ).xyz( atom_name3 )
																					 );

		if (bond_angle < 0.001 ) TR << "WHAT THE HELL????????? " << std::endl;

		core::scoring::func::FuncOP angle_harm_func_( new core::scoring::func::HarmonicFunc( bond_angle, bond_angle_sd_ ));
		cst_set->add_constraint( new AngleConstraint(
																								 atom_id2 , atom_id1, atom_id3, angle_harm_func_,	rna_bond_geometry ) );

		if ( false ) {
			TR << "PUTTING CONSTRAINT ON ANGLE: " <<
				atom_id2.rsd() << " " << pose_reference.residue( atom_id2.rsd() ).atom_name( atom_id2.atomno() ) << "; "  <<
				atom_id1.rsd() << " " << pose_reference.residue( atom_id1.rsd() ).atom_name( atom_id1.atomno() ) << "; "  <<
				atom_id3.rsd() << " " << pose_reference.residue( atom_id3.rsd() ).atom_name( atom_id3.atomno() ) << " ==> "  << degrees( bond_angle ) << " " << degrees( bond_angle_sd_ ) <<
				std::endl;
		}

	}

}

//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
bool
ConstrainToIdealMover::check_if_really_connected(
  core::pose::Pose const & pose,
	core::id::AtomID const & atom_id1,
	core::id::AtomID const & atom_id2)
{
	if (  atom_id1.rsd() == atom_id2.rsd() ) return true;


	core::kinematics::tree::AtomCOP atom1 ( & pose.atom_tree().atom( atom_id1 ) );
	core::kinematics::tree::AtomCOP atom2 ( & pose.atom_tree().atom( atom_id2 ) );

	if ( atom1->parent() == atom2 ) return true;
	if ( atom2->parent() == atom1 ) return true;

	return false;
}

//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
bool
ConstrainToIdealMover::i_want_this_atom_to_move( core::conformation::Residue const & residue2, core::Size const & k )
{

	//if (k > residue2.first_sidechain_atom() &&
		//k != core::chemical::rna::first_base_atom_index( residue2 ) ) return false;
	core::id::AtomID id( k, residue2.seqpos() );
	if ( !allow_insert_->get( id ) ) return false;

	if ( residue2.is_virtual( k ) ) {
		//		std::cout << "Is this virtual? " << residue2.atom_name( k ) << std::endl;
		return false;
	}

	return true;

}

//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
bool
ConstrainToIdealMover::i_want_this_atom_to_move( core::pose::Pose const & pose, core::id::AtomID const & atom_id )
{
	core::conformation::Residue const & residue( pose.residue( atom_id.rsd() ) );
	core::Size const & k( atom_id.atomno() );
	
	if ( !allow_insert_->get( atom_id ) ) return false;
	
	if ( residue.is_virtual( k ) ) return false;
	
	return true;
	//return i_want_this_atom_to_move( pose.residue( atom_id.rsd() ) ,
																	 //atom_id.atomno() );

}
//////////////////////////////////////////////////////////////////////////////
// Following has not (yet) been carefully debugged. --Rhiju
//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
void
ConstrainToIdealMover::vary_bond_geometry(
	//core::kinematics::MoveMap & mm,
	core::pose::Pose & pose,
	core::pose::Pose const & pose_reference ) {

	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace numeric::conversions;

	ConstraintSetOP cst_set = pose.constraint_set()->clone();
	pose.constraint_set( cst_set );

	core::Size const nres( pose.total_residue() );

	std::map< AtomID, utility::vector1< AtomID > > lists_of_angle_bonded_atoms;

	for  (core::Size i = 1; i <= nres; i++ )  {

		if ( !allow_insert_->get( i ) ) continue;

		core::conformation::Residue const & residue( pose.residue( i )  );

		for (core::Size j = 1; j <= residue.natoms(); j++ ) {

			if ( !i_want_this_atom_to_move( residue, j ) ) continue;

			core::kinematics::tree::AtomCOP current_atom ( & pose.atom_tree().atom( AtomID(j,i) ) );
			if ( current_atom->is_jump() ) continue;

			///////////////////
			core::kinematics::tree::AtomCOP input_stub_atom1( current_atom->input_stub_atom1() );
			if ( !input_stub_atom1 ) continue;
			if ( !i_want_this_atom_to_move( pose, input_stub_atom1->id() ) ) continue;
			mm_->set( DOF_ID( AtomID( j, i ), D ), true );
			if ( input_stub_atom1->is_jump() ) continue;

			core::kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );

			///////////////////
			if ( !input_stub_atom2 ) continue;
			if ( input_stub_atom2 == current_atom ) continue;
			if ( !i_want_this_atom_to_move( pose, input_stub_atom2->id() ) ) continue;
			mm_->set( DOF_ID( AtomID( j, i ), THETA ), true );
			if ( input_stub_atom2->is_jump() ) continue;

			///////////////////
			core::kinematics::tree::AtomCOP input_stub_atom3( current_atom->input_stub_atom3() );
			if ( !input_stub_atom3 ) continue;
			if ( !i_want_this_atom_to_move( pose, input_stub_atom3->id() ) ) continue;
			if ( input_stub_atom3 == current_atom ) continue;
			mm_->set( DOF_ID( AtomID( j, i ), PHI ), true );

		}
	}

	utility::vector1< std::pair< AtomID, AtomID > >  bond_list;
	utility::vector1< std::pair< AtomID, std::pair< AtomID, AtomID > > > bond_angle_list;

 	for  (core::Size i = 1; i <= nres; i++ )  {

		//Go through all bonds in pose...
		if ( !allow_insert_->get( i ) ) continue;

		core::conformation::Residue const & residue( pose.residue( i )  );

		for (core::Size j = 1; j <= residue.natoms(); j++ ) {

			if ( !i_want_this_atom_to_move( residue, j ) ) continue;

			AtomID atom_id1( j, i );

			utility::vector1< AtomID >  nbrs(	pose.conformation().bonded_neighbor_all_res( atom_id1 ) );

			// Bond lengths.
			for ( core::Size n = 1; n <= nbrs.size(); n++ ) {

				AtomID const & atom_id2( nbrs[ n ] );

				core::conformation::Residue const & residue2( pose.residue( atom_id2.rsd() ) ) ;
				core::Size const & k( atom_id2.atomno() ) ;

				if ( ! check_if_really_connected( pose, atom_id1, atom_id2) ) continue;

				if ( i_want_this_atom_to_move( residue2, k ) )  {
					add_bond_constraint( atom_id1, atom_id2,
															 bond_list,
															 pose, pose_reference, cst_set );
				}

			}

			// Bond angles
			for ( core::Size m = 1; m <= nbrs.size(); m++ ) {
				AtomID const & atom_id2( nbrs[ m ] );

				core::conformation::Residue const & residue2( pose.residue( atom_id2.rsd() ) ) ;
				if ( ! check_if_really_connected( pose, atom_id1, atom_id2) ) continue;

				core::Size const & k( atom_id2.atomno() ) ;

				for ( core::Size n = 1; n <= nbrs.size(); n++ ) {
					AtomID const & atom_id3( nbrs[ n ] );

					core::conformation::Residue const & residue3( pose.residue( atom_id3.rsd() ) ) ;
					if ( ! check_if_really_connected( pose, atom_id1, atom_id3) ) continue;

					core::Size const & q( atom_id3.atomno() ) ;

					if ( i_want_this_atom_to_move( residue2, k ) &&
							 i_want_this_atom_to_move( residue3, q ) )  {
						add_bond_angle_constraint( atom_id1, atom_id2, atom_id3, bond_angle_list,
																			 pose, pose_reference, cst_set );
					}

				}
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
	RNA_FittedTorsionInfo const rna_fitted_torsion_info;

	core::Real const DELTA_CUTOFF( rna_fitted_torsion_info.delta_cutoff() );

	for (core::Size n = 1; n <= pose.total_residue(); n++ ) {
		core::Real const delta = pose.residue( n ).mainchain_torsion( DELTA );
		if ( delta > DELTA_CUTOFF ) {
			core::pose::rna::apply_ideal_c2endo_sugar_coords( pose_reference, n );
		}
	}

}

///////////////////////////////////////////////////
//copied from src/protocols/rna/RNA_Minimizer.cc, SVN #44771
void
ConstrainToIdealMover::create_pose_reference(
  core::pose::Pose const & pose,
	core::pose::Pose & pose_reference )
{
	using namespace core::chemical;
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	make_pose_from_sequence( pose_reference, pose.sequence(),	*rsd_set );
	apply_ideal_coordinates_for_alternative_pucker( pose, pose_reference );
}

} // namespace moves
} // namespace protocols
