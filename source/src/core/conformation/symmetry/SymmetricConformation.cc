// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/symmetry/SymmetricConformation.hh
/// @brief  symmetry conformation container.
//     Contains overloaded functions needed to
//     make changes in conformation symmetric
/// @author Phil Bradley, Ingemar Andre

// Unit headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryTransform.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <core/id/TorsionID.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>

#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.hh>
#include <core/conformation/ResidueFactory.hh>
#include <utility/graph/UpperEdgeGraph.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/conformation/util.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/HomogeneousTransform.srlz.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace symmetry {

static basic::Tracer TR( "core.conformation.symmetry.Conformation" );

/// @brief  Default CTOR
SymmetricConformation::SymmetricConformation():
	Conformation(),
	symm_info_()
{
	clear_Tsymm( );
}

/// @brief  Default CTOR
SymmetricConformation::SymmetricConformation(Conformation const & conf, SymmetryInfo const & symm_info):
	Conformation(conf),
	symm_info_()
{
	symm_info_ = symm_info.clone();
	clear_Tsymm( );  // force recompute
}


Conformation &
SymmetricConformation::operator=( Conformation const & src )
{
	auto const * sym_conf = dynamic_cast< SymmetricConformation const * > ( &src );
	if ( sym_conf ) {
		SymmetricConformation::operator=( *sym_conf );
	} else {
		//Conformation::operator=( src );
		throw CREATE_EXCEPTION(utility::excn::Exception,  "SymmetricConformation::operator= was handed a non-SymmetricConformation" );
	}

	return *this;
}

void
SymmetricConformation::detached_copy( Conformation const & src ) {
	Conformation::detached_copy( src );
	auto const * sym_conf = dynamic_cast< SymmetricConformation const * > ( &src );
	if ( sym_conf ) {

		// Copy the private data members of the symmetric conformation.
		symm_info_ = sym_conf->symm_info_->clone();
		Tsymm_ = sym_conf->Tsymm_;

	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "SymmetricConformation detached_copy was handed a non-SymmetricConformation" );
	}
}


/// @details make a copy of this conformation( allocate actual memory for it )
ConformationOP
SymmetricConformation::clone() const
{
	SymmetricConformationOP copy( new SymmetricConformation(*this, *symm_info_) );
	return copy;
}

bool
SymmetricConformation::same_type_as_me( Conformation const & other, bool const recurse  /* = true */ ) const
{
	if ( ! dynamic_cast< SymmetricConformation const * > ( &other) ) {
		return false;
	}
	if ( recurse ) {
		return other.same_type_as_me( *this, false );
	} else {
		return true;
	}
}

SymmetryInfoCOP
SymmetricConformation::Symmetry_Info() const
{
	return symm_info_;
}

SymmetryInfoOP
SymmetricConformation::Symmetry_Info()
{
	return symm_info_;
}

SymmetricConformation::~SymmetricConformation() { clear(); }

/// DOF
void
SymmetricConformation::set_dof( DOF_ID const & id, Real const setting )
{
	using DOF_IDs = SymmetryInfo::DOF_IDs;

	core::Size parent_rsd;
	if ( symm_info_->torsion_changes_move_other_monomers() ) clear_Tsymm( );

	if ( !symm_info_->dof_is_independent( id, *this ) ) {
		TR.Debug << "SymmetricConformation:: directly setting a dependent DOF!, try to set its parent" << std::endl;
		if ( id.type() >= id::RB1 && id.type() <= id::RB6 ) {
			int parent_jump = symm_info_->jump_follows( fold_tree().get_jump_that_builds_residue( id.rsd() ) );
			parent_rsd = fold_tree().downstream_jump_residue( parent_jump );
		} else {
			parent_rsd = symm_info_->bb_follows( id.rsd() ) ;
		}
	} else {
		parent_rsd = id.rsd();
	}

	if ( id.type() >= id::RB1 && id.type() <= id::RB6 ) {
		clear_Tsymm( );
	}

	id::DOF_ID parent_id ( id::AtomID( id.atomno(), parent_rsd ), id.type() );

	// set this DOF using base-class method
	Conformation::set_dof( parent_id, setting );
	{
		DOF_IDs const & dofs( symm_info_->dependent_dofs( parent_id, *this ) );
		for ( auto const & dof : dofs ) {
			Conformation::set_dof( dof, setting );
		}
	}
}

/// @brief set the secondary structure of a sequence position
/// @details Sets secondary structure character of a sequence position.
///  Will resize the secondary structure array if the requested sequence
///  position is larger than the length of the array.
void
SymmetricConformation::set_secstruct( Size const seqpos, char const setting )
{
	Conformation::set_secstruct( seqpos, setting );
	if ( symm_info_->bb_is_independent( seqpos ) ) {
		for ( unsigned long pos : symm_info_->bb_clones( seqpos ) ) {
			Conformation::set_secstruct( pos, setting );
		}
	} else {
		TR.Debug << "SymmetricConformation:: Setting secstruct for dependent residue!, try to set its parent" << std::endl;
	}
}

/// set the fold_tree, update symmetryinfo
void
SymmetricConformation::fold_tree( FoldTree const & fold_tree_in )
{
	Conformation::fold_tree( fold_tree_in );
	FoldTree asymm_f = get_asymm_unit_fold_tree(*this);
	symm_info_->update_nmonomer_jumps( asymm_f.num_jump()-1 );  // always an extra jump
}

/// TORSIONS
void
SymmetricConformation::set_torsion( id::TorsionID const & id, Real const setting )
{
	using TorsionIDs = SymmetryInfo::TorsionIDs;

	TR.Trace << "SymmetricConformation: set_torsion: " << id << ' ' << setting << std::endl;
	if ( symm_info_->torsion_changes_move_other_monomers() ) clear_Tsymm( );

	core::Size parent_rsd;

	if ( !symm_info_->torsion_is_independent( id ) ) {
		TR.Debug << "SymmetricConformation:: directly setting a dependent TORSION!, try to set its parent" << std::endl;
		parent_rsd = symm_info_->bb_follows( id.rsd() ) ;
	} else {
		parent_rsd = id.rsd();
	}

	id::TorsionID parent_id ( parent_rsd, id.type(), id.torsion() );

	// set this DOF using base-class method
	Conformation::set_torsion( parent_id, setting );
	{
		TorsionIDs const & tors( symm_info_->dependent_torsions( parent_id ) );
		for ( auto const & tor : tors ) {
			Conformation::set_torsion( tor, setting );
		}
	}
}

/// TORSION ANGLES
void
SymmetricConformation::set_torsion_angle(
	AtomID const & atom1,
	AtomID const & atom2,
	AtomID const & atom3,
	AtomID const & atom4,
	Real const setting,
	bool const quiet
)
{
	AtomID parent_atom1, parent_atom2, parent_atom3, parent_atom4;

	TR.Trace << "SymmetricConformation: set_torsion_angle: " << atom1 << "+ to " << setting << std::endl;
	if ( symm_info_->torsion_changes_move_other_monomers() ) clear_Tsymm( );

	// assume if 1st atom is dependent, all are
	if ( !symm_info_->bb_is_independent( atom1.rsd() ) ) {
		TR.Debug << "SymmetricConformation:: directly setting a dependent ANGLE!, try to set its parent" << std::endl;
		parent_atom1 = id::AtomID( atom1.atomno(), symm_info_->bb_follows( atom1.rsd() ) );
		parent_atom2 = id::AtomID( atom2.atomno(), symm_info_->bb_follows( atom2.rsd() ) );
		parent_atom3 = id::AtomID( atom3.atomno(), symm_info_->bb_follows( atom3.rsd() ) );
		parent_atom4 = id::AtomID( atom4.atomno(), symm_info_->bb_follows( atom4.rsd() ) );
	} else {
		parent_atom1=atom1;
		parent_atom2=atom2;
		parent_atom3=atom3;
		parent_atom4=atom4;
	}

	// set this DOF using base-class method
	Conformation::set_torsion_angle( parent_atom1, parent_atom2, parent_atom3, parent_atom4, setting, quiet );

	{
		core::Size nclones = symm_info_->num_bb_clones( );
		SymmetryInfo::Clones clones1 = symm_info_->bb_clones( parent_atom1.rsd() );
		SymmetryInfo::Clones clones2 = symm_info_->bb_clones( parent_atom2.rsd() );
		SymmetryInfo::Clones clones3 = symm_info_->bb_clones( parent_atom3.rsd() );
		SymmetryInfo::Clones clones4 = symm_info_->bb_clones( parent_atom4.rsd() );
		AtomID dep_atom1, dep_atom2, dep_atom3, dep_atom4;
		for ( int i=1; i<=(int)nclones; ++i ) {
			dep_atom1 = id::AtomID( atom1.atomno(), clones1[i] );
			dep_atom2 = id::AtomID( atom2.atomno(), clones2[i] );
			dep_atom3 = id::AtomID( atom3.atomno(), clones3[i] );
			dep_atom4 = id::AtomID( atom4.atomno(), clones4[i] );
			Conformation::set_torsion_angle( dep_atom1, dep_atom2, dep_atom3, dep_atom4, setting, quiet );
		}
	}
}


/// BOND ANGLES
void
SymmetricConformation::set_bond_angle(
	AtomID const & atom1,
	AtomID const & atom2,
	AtomID const & atom3,
	Real const setting
)
{
	AtomID parent_atom1, parent_atom2, parent_atom3;

	TR.Trace << "SymmetricConformation: set_bond_angle: " << atom1 << "+ to " << setting << std::endl;
	if ( symm_info_->torsion_changes_move_other_monomers() ) clear_Tsymm( );

	// assume if 1st atom is dependent, all are
	if ( !symm_info_->bb_is_independent( atom1.rsd() ) ) {
		TR.Debug << "SymmetricConformation:: directly setting a dependent ANGLE!, try to set its parent" << std::endl;
		parent_atom1 = id::AtomID( atom1.atomno(), symm_info_->bb_follows( atom1.rsd() ) );
		parent_atom2 = id::AtomID( atom2.atomno(), symm_info_->bb_follows( atom2.rsd() ) );
		parent_atom3 = id::AtomID( atom3.atomno(), symm_info_->bb_follows( atom3.rsd() ) );
	} else {
		parent_atom1=atom1;
		parent_atom2=atom2;
		parent_atom3=atom3;
	}

	// set this DOF using base-class method
	Conformation::set_bond_angle( parent_atom1, parent_atom2, parent_atom3, setting );

	{
		core::Size nclones = symm_info_->num_bb_clones( );
		SymmetryInfo::Clones clones1 = symm_info_->bb_clones( parent_atom1.rsd() );
		SymmetryInfo::Clones clones2 = symm_info_->bb_clones( parent_atom2.rsd() );
		SymmetryInfo::Clones clones3 = symm_info_->bb_clones( parent_atom3.rsd() );
		AtomID dep_atom1, dep_atom2, dep_atom3;
		for ( int i=1; i<=(int)nclones; ++i ) {
			dep_atom1 = id::AtomID( atom1.atomno(), clones1[i] );
			dep_atom2 = id::AtomID( atom2.atomno(), clones2[i] );
			dep_atom3 = id::AtomID( atom3.atomno(), clones3[i] );
			Conformation::set_bond_angle( dep_atom1, dep_atom2, dep_atom3, setting );
		}
	}
}

/// BOND LENGTHS
void
SymmetricConformation::set_bond_length(
	AtomID const & atom1,
	AtomID const & atom2,
	Real const setting
)
{
	AtomID parent_atom1, parent_atom2;

	TR.Trace << "SymmetricConformation: set_bond_length: " << atom1 << "+ to " << setting << std::endl;
	if ( symm_info_->torsion_changes_move_other_monomers() ) clear_Tsymm( );

	// assume if 1st atom is dependent, all are
	if ( !symm_info_->bb_is_independent( atom1.rsd() ) ) {
		TR.Debug << "SymmetricConformation:: directly setting a dependent BONDLENGTH!, try to set its parent" << std::endl;
		parent_atom1 = id::AtomID( atom1.atomno(), symm_info_->bb_follows( atom1.rsd() ) );
		parent_atom2 = id::AtomID( atom2.atomno(), symm_info_->bb_follows( atom2.rsd() ) );
	} else {
		parent_atom1=atom1;
		parent_atom2=atom2;
	}

	// set this DOF using base-class method
	Conformation::set_bond_length( parent_atom1, parent_atom2, setting );

	{
		core::Size nclones = symm_info_->num_bb_clones( );
		SymmetryInfo::Clones clones1 = symm_info_->bb_clones( parent_atom1.rsd() );
		SymmetryInfo::Clones clones2 = symm_info_->bb_clones( parent_atom2.rsd() );
		for ( int i=1; i<=(int)nclones; ++i ) {
			AtomID dep_atom1( atom1.atomno(), clones1[i] );
			AtomID dep_atom2( atom2.atomno(), clones2[i] );
			Conformation::set_bond_length( dep_atom1, dep_atom2, setting );
		}
	}
}


/// JUMPS
void
SymmetricConformation::set_jump( int const jump_number, kinematics::Jump const & new_jump )
{
	TR.Trace << "SymmetricConformation: set_jump jump_number: " << jump_number << std::endl;

	// clear cached transforms
	clear_Tsymm( );

	id::AtomID const id( Conformation::jump_atom_id( jump_number ) );
	Conformation::set_jump( id, new_jump );

	if ( !symm_info_->jump_is_independent( jump_number ) ) {
		TR.Warning << "SymmetricConformation:: directly setting a dependent ATOM!" << std::endl;
		TR.Warning << "the jump " << jump_number << " is controlled by " << symm_info_->jump_follows( jump_number ) << std::endl;
	} else {
		for ( unsigned long pos : symm_info_->jump_clones( jump_number ) ) {
			id::AtomID const id_clone( Conformation::jump_atom_id( pos ) );
			Conformation::set_jump( id_clone, new_jump );
		}
	}
}

void
SymmetricConformation::set_jump( id::AtomID const & id, kinematics::Jump const & new_jump )
{
	TR.Trace << "SymmetricConformation: set_jump id:" << id << std::endl;

	// clear cached transforms
	clear_Tsymm( );

	Conformation::set_jump( id, new_jump );

	int const jump_number ( fold_tree().get_jump_that_builds_residue( id.rsd() ) );
	if ( !symm_info_->jump_is_independent( jump_number ) ) {
		TR.Warning << "SymmetricConformation:: directly setting a dependent ATOM!" << std::endl;
	} else {
		for ( unsigned long pos : symm_info_->jump_clones( jump_number ) ) {
			id::AtomID const id_clone( Conformation::jump_atom_id( pos ) );
			Conformation::set_jump( id_clone, new_jump );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @details Returns a mask of residues over which scoring is restricted.
/// Only these residues will be used in constructing the neighbor list
utility::vector1<bool>
SymmetricConformation::get_residue_mask() const {
	return symm_info_->independent_residues();
}

core::Real
SymmetricConformation::get_residue_weight(core::Size resid1, core::Size resid2) const {
	return symm_info_->score_multiply(resid1,resid2);
}

/// @details symmetry-safe replace residue
void
SymmetricConformation::replace_residue( Size const seqpos, Residue const & new_rsd_in, bool const orient_backbone ) {
	TR.Debug << "SymmetricConformation: replace_residue: " << seqpos << std::endl;

	core::Size parent_rsd;
	Residue new_rsd = new_rsd_in;

	if ( !symm_info_->bb_is_independent( seqpos ) ) {
		TR.Debug << "SymmetricConformation::replace_residue(2) directly setting a dependent TORSION!, try to set its parent" << std::endl;

		parent_rsd = symm_info_->bb_follows( seqpos ) ;

		// if we're not orienting backbone, we need to transform the coord frame to that of the ind. subunit
		if ( !orient_backbone ) {
			// make the new res
			for ( int i=1; i<=(int)new_rsd.natoms(); ++i ) {
				new_rsd.set_xyz(i , apply_transformation( new_rsd_in.xyz(i), seqpos, parent_rsd ) );
			}
		}
	} else {
		parent_rsd = seqpos;
	}

	// set this residue using base-class method
	Conformation::replace_residue( parent_rsd, new_rsd, orient_backbone );

	// now the copies
	for ( auto pos=symm_info_->bb_clones( parent_rsd ).begin(),
			epos=symm_info_->bb_clones( parent_rsd ).end(); pos != epos; ++pos ) {
		// if we're not orienting backbone, update symm copies using local coord frame for each
		//fpd --> we already oriented backbone for master, that should be good enough
		if ( !orient_backbone ) {
			// make the new res
			Residue new_new_rsd = new_rsd;
			for ( int i=1; i<=(int)new_new_rsd.natoms(); ++i ) {
				new_new_rsd.set_xyz(i , apply_transformation( new_rsd.xyz(i), parent_rsd, *pos ) );
			}
			Conformation::replace_residue( *pos, new_new_rsd, false );
		} else {
			Conformation::replace_residue( *pos, new_rsd, orient_backbone );
		}
	}

	//fpd may have implicitly changed a jump
	clear_Tsymm( );
}

/// @details symmetry-safe replace residue
void
SymmetricConformation::replace_residue( Size const seqpos,
	Residue const & new_rsd,
	utility::vector1< std::pair< std::string, std::string > > const & atom_pairs )  {
	core::Size parent_rsd;

	if ( !symm_info_->bb_is_independent( seqpos ) ) {
		TR.Debug << "SymmetricConformation::replace_residue(3) setting a dependent TORSION!, try to set its parent" << std::endl;
		parent_rsd = symm_info_->bb_follows( seqpos ) ;
	} else {
		parent_rsd = seqpos;
	}

	// set this residue using base-class method
	Conformation::replace_residue( parent_rsd, new_rsd, atom_pairs );

	// now the copies
	for ( auto pos=symm_info_->bb_clones( parent_rsd ).begin(),
			epos=symm_info_->bb_clones( parent_rsd ).end(); pos != epos; ++pos ) {
		//fpd Same logic as above; we already oriented backbone for master, that should be good enough
		// make the new res
		Residue new_new_rsd = new_rsd;
		for ( int i=1; i<=(int)new_new_rsd.natoms(); ++i ) {
			new_new_rsd.set_xyz(i , apply_transformation( new_rsd.xyz(i), parent_rsd, *pos ) );
		}
		Conformation::replace_residue( *pos, new_new_rsd, false );
	}

	//fpd may have implicitly changed a jump (this should really check for downstream jump residues)
	clear_Tsymm( );
}


core::Size
SymmetricConformation::get_upstream_vrt( Size seqpos ) const {
	if ( this->residue( seqpos ).aa() == core::chemical::aa_vrt ) {
		return seqpos;
	}

	// find peptide segment that contains this res (?)
	core::kinematics::Edge const &e = fold_tree().get_residue_edge( seqpos );
	core::Size curr_res = e.start();

	while ( this->residue( curr_res ).aa() != core::chemical::aa_vrt ) {
		core::kinematics::Edge const &e_i = fold_tree().get_residue_edge( curr_res );
		curr_res = e_i.start();
	}

	return curr_res;
}

// Get the transformation controlling resid i
numeric::HomogeneousTransform< core::Real >
SymmetricConformation::get_transformation( core::Size resid, bool input_is_subunit_id /*=false*/ )  const {
	if ( Tsymm_.size() == 0 ) {
		recalculate_transforms( );
	}

	std::string compid = "A";
	if ( symm_info_->get_num_components() >= 2 ) {
		compid = symm_info_->get_component_of_residue( resid );
	}

	core::Size subunit = input_is_subunit_id ? resid : symm_info_->subunit_index( resid );
	return Tsymm_[ compid ][ subunit ].ht();
}

//  Remap coordinate X from resid i to j
PointPosition
SymmetricConformation::apply_transformation(
	PointPosition Xin,
	core::Size residfrom,
	core::Size residto,
	bool rotationonly ) {
	if ( Tsymm_.size() == 0 ) {
		recalculate_transforms( );
	}

	std::string compid = "A";
	if ( symm_info_->get_num_components() >= 2 ) {
		debug_assert( symm_info_->get_component_of_residue(residfrom) == symm_info_->get_component_of_residue(residto) );
		compid = symm_info_->get_component_of_residue( residfrom );
	}

	runtime_assert(Tsymm_.find(compid) != Tsymm_.end());
	utility::vector1< SymmetryTransform > const &T_i = (Tsymm_.find(compid))->second;
	numeric::HomogeneousTransform< core::Real > Tsymm_from = T_i[ symm_info_->subunit_index( residfrom ) ].ht();
	numeric::HomogeneousTransform< core::Real > Tsymm_to = T_i[ symm_info_->subunit_index( residto ) ].ht();

	if ( rotationonly ) {
		Tsymm_from.set_identity_transform();
		Tsymm_to.set_identity_transform();
	}

	PointPosition Xout;
	Xout = (Tsymm_from.inverse() * Xin);
	Xout = (Tsymm_to * Xout);
	return Xout;
}

//  Remap coordinate X from resid i to j
PointPosition
SymmetricConformation::apply_transformation_norecompute(
	PointPosition Xin,
	core::Size residfrom,
	core::Size residto,
	bool rotationonly ) const {
	std::string compid = "A";
	if ( symm_info_->get_num_components() >= 2 ) {
		debug_assert( symm_info_->get_component_of_residue(residfrom) == symm_info_->get_component_of_residue(residto) );
		compid = symm_info_->get_component_of_residue( residfrom );
	}

	runtime_assert(Tsymm_.find(compid) != Tsymm_.end());
	utility::vector1< SymmetryTransform > const &T_i = (Tsymm_.find(compid))->second;
	numeric::HomogeneousTransform< core::Real > Tsymm_from = T_i[ symm_info_->subunit_index( residfrom ) ].ht();
	numeric::HomogeneousTransform< core::Real > Tsymm_to = T_i[ symm_info_->subunit_index( residto ) ].ht();

	if ( rotationonly ) {
		Tsymm_from.set_identity_transform();
		Tsymm_to.set_identity_transform();
	}

	PointPosition Xout;
	Xout = (Tsymm_from.inverse() * Xin);
	Xout = (Tsymm_to * Xout);
	return Xout;
}

//  Symmetric set_xyz
void
SymmetricConformation::set_xyz(
	AtomID const & id,
	PointPosition const & position ) {
	TR.Debug << "SymmetricConformation::set_xyz: " << id << std::endl;

	AtomID parent_id;
	PointPosition parent_pos;

	// pass set_xyz on symmetric vrts directly to the base class
	// this is potentially dangerous but may be useful
	if ( id.rsd() > symm_info_->num_total_residues_without_pseudo() ) {
		Conformation::set_xyz( id, position );
		TR<< "***WARN*** XYZ set vrt!" << std::endl;
		clear_Tsymm( );
		return;
	}

	if ( !symm_info_->bb_is_independent( id.rsd() ) ) {
		TR.Debug << "SymmetricConformation::set_xyz setting a dependent XYZ; remapping to its parent" << std::endl;
		parent_id = AtomID( id.atomno(), symm_info_->bb_follows( id.rsd() ) );
		parent_pos = apply_transformation( position, id.rsd(), parent_id.rsd() );
	} else {
		parent_id = id;
		parent_pos = position;
	}

	// set parent XYZ using base-class method, followed by all copies
	Conformation::set_xyz( parent_id, parent_pos );
	for ( auto pos=symm_info_->bb_clones( parent_id.rsd() ).begin(),
			epos=symm_info_->bb_clones( parent_id.rsd() ).end(); pos != epos; ++pos ) {
		AtomID id_i = AtomID( id.atomno(), *pos );
		PointPosition pos_i = apply_transformation( parent_pos, parent_id.rsd(), *pos );
		Conformation::set_xyz( id_i, pos_i );
	}
}

// Symmetric batch_set_xyz
void
SymmetricConformation::batch_set_xyz(
	utility::vector1<AtomID> const & ids,
	utility::vector1<PointPosition> const & positions
) {
	runtime_assert( ids.size() == positions.size() );
	TR.Debug << "SymmetricConformation::batch_set_xyz" << std::endl;

	utility::vector1<AtomID> ids_with_symm;
	utility::vector1<PointPosition> positions_with_symm;

	// pass set_xyz on symmetric vrts directly to the base class
	// this is potentially dangerous but may be useful
	for ( int i=1; i<=(int)ids.size(); ++i ) {
		AtomID id=ids[i], parent_id;
		PointPosition position=positions[i], parent_pos;

		if ( id.rsd() > symm_info_->num_total_residues_without_pseudo() ) {
			TR<< "***WARN*** XYZ batch set vrt!" << std::endl;
			ids_with_symm.push_back( id );
			positions_with_symm.push_back( position );
			clear_Tsymm( );
		} else {
			if ( !symm_info_->bb_is_independent( id.rsd() ) ) {
				TR.Debug << "SymmetricConformation::set_xyz setting a dependent XYZ; remapping to its parent" << std::endl;
				parent_id = AtomID( id.atomno(), symm_info_->bb_follows( id.rsd() ) );
				parent_pos = apply_transformation( position, id.rsd(), parent_id.rsd() );
			} else {
				parent_id = id;
				parent_pos = position;
			}

			// set parent XYZ using base-class method, followed by all copies
			ids_with_symm.push_back( parent_id );
			positions_with_symm.push_back( parent_pos );
			for ( auto pos=symm_info_->bb_clones( parent_id.rsd() ).begin(),
					epos=symm_info_->bb_clones( parent_id.rsd() ).end(); pos != epos; ++pos ) {
				AtomID id_i = AtomID( id.atomno(), *pos );
				PointPosition pos_i = apply_transformation( parent_pos, parent_id.rsd(), *pos );
				ids_with_symm.push_back( id_i );
				positions_with_symm.push_back( pos_i );
			}
		}
	}

	Conformation::batch_set_xyz( ids_with_symm, positions_with_symm );
}

void
SymmetricConformation::apply_transform_Rx_plus_v(
	numeric::xyzMatrix< Real > const & R,
	Vector const & v
) {
	utility::vector1<AtomID> ids;
	utility::vector1<Vector> xyzs;
	for ( Size i = 1; i <= size(); ++i ) {
		for ( Size j = 1; j <= residue_type(i).natoms(); ++j ) {
			AtomID id( j, i );
			ids.push_back(id);
			xyzs.push_back( R * xyz(id) + v );
		}
	}

	//fpd Explicitly call base class
	Conformation::batch_set_xyz( ids, xyzs );

	//fpd then invalidate Tsymm
	clear_Tsymm( );
}


// @brief invalidate current Tsymm settings
void
SymmetricConformation::clear_Tsymm( ) const {
	Tsymm_.clear();
}

// @brief invert one transform about Z
void
SymmetricConformation::invert_Tsymm( std::string const & sub, core::Size subunit ) const {
	Tsymm_[sub][subunit].set_mirror_z(true);
}


// @brief recalculate the Tsymm_ transforms using the current pose
void
SymmetricConformation::recalculate_transforms( ) const {
	using namespace numeric;

	// clear current xforms
	clear_Tsymm( );

	// rebuild based on virtuals
	core::Size nsubunits = symm_info_->subunits();
	//core::Size nresMonomer = symm_info_->num_independent_residues();

	FoldTree const & f( fold_tree() );

	ResidueOP vrt_res_op( nullptr ); // maybe unused

	// multicomp support
	Size ncomps = symm_info_->get_num_components();

	// check to see if all subunits have an upstream jump to vrt res
	// do this because in some setups only one subunit has an upstream jump from a vrt res, would be
	// bad to use for that guy and not others
	bool each_subunit_has_vrt_control( true );
	{ // scope
		for ( Size icomp=1; icomp <= ncomps && each_subunit_has_vrt_control; ++icomp ) {
			std::string comptag = (ncomps==1)? "A" : symm_info_->get_component(icomp);
			for ( Size isub=1; isub <= nsubunits && each_subunit_has_vrt_control; ++isub ) {
				Size vrt_ctrl(0);
				for ( Size j=1; j<= f.num_jump(); ++j ) {
					Size downstream = f.downstream_jump_residue(j), upstream = f.upstream_jump_residue(j);
					if ( this->residue( upstream ).aa() == chemical::aa_vrt
							&& downstream <= symm_info_->num_total_residues_without_pseudo()
							&& symm_info_->subunit_index( downstream ) == isub
							&& ( ncomps==1 || symm_info_->get_component_of_residue( downstream ) == comptag ) ) {
						vrt_ctrl = f.upstream_jump_residue(j);
						//TR.Trace << "Found control virtual: subunit= " << isub << " vrt_ctrl= " << vrt_ctrl << std::endl;
						break;
					}
				}
				if ( !vrt_ctrl ) {
					each_subunit_has_vrt_control = false;
				}
			}
		}
	}

	for ( Size icomp=1; icomp <= ncomps; ++icomp ) {
		std::string comptag = (ncomps==1)? "A" : symm_info_->get_component(icomp);
		Tsymm_[comptag] = utility::vector1< SymmetryTransform >(nsubunits);

		for ( Size isub=1; isub <= nsubunits; ++isub ) {
			Size vrt_ctrl(0);
			if ( each_subunit_has_vrt_control ) {
				for ( Size j=1; j<= f.num_jump(); ++j ) {
					Size downstream = f.downstream_jump_residue(j), upstream = f.upstream_jump_residue(j);
					if ( this->residue( upstream ).aa() == chemical::aa_vrt
							&& downstream <= symm_info_->num_total_residues_without_pseudo()
							&& symm_info_->subunit_index( downstream ) == isub
							&& ( ncomps==1 || symm_info_->get_component_of_residue( downstream ) == comptag ) ) {
						vrt_ctrl = f.upstream_jump_residue(j);
						//TR.Trace << "Found control virtual: subunit= " << isub << " vrt_ctrl= " << vrt_ctrl << std::endl;
						break;
					}
				}
			}

			Residue const * vrt_res_cap = nullptr;
			if ( vrt_ctrl ) {
				vrt_res_cap = &( residue( vrt_ctrl ) ); /// the standard logic, if virtual rsds are present
			} else {
				/// construct a pseudo virtual using the "first residue in each monomer"
				Size first_independent_res( 0 );
				for ( Size i=1; i<= this->size(); ++i ) {
					if ( symm_info_->chi_is_independent( i ) ) {
						first_independent_res = i;
						break;
					}
				}
				runtime_assert( first_independent_res );

				Size coordframe_pos( 0 );
				if ( symm_info_->subunit_index( first_independent_res ) == isub ) {
					// we are the independent monomer; note that indep monomer not necessarily numbered 1...
					coordframe_pos = first_independent_res;
				} else {
					for ( Size i=1; i<= this->size(); ++i ) {
						if ( symm_info_->chi_follows(i) == first_independent_res &&
								symm_info_->subunit_index(i) == isub ) {
							coordframe_pos = i;
							break;
						}
					}
				}
				if ( !coordframe_pos ) {
					std::cout << "unable to get coordframe_pos: " << isub << ' ' << first_independent_res << std::endl;
					for ( Size i=1; i<= this->size(); ++i ) {
						std::cout << "details: " << i << ' ' << symm_info_->chi_follows(i) << ' ' <<
							symm_info_->subunit_index(i) << std::endl;
					}
				}
				runtime_assert( coordframe_pos );
				runtime_assert( symm_info_->subunit_index( coordframe_pos ) == isub );

				if ( !vrt_res_op ) {
					// slow, create residueop
					// Uses the typeset for the conformation as a whole.
					// Is there a reason we'd want to use the typeset for the coordframe_rsd instead, if they differ?
					vrt_res_op = conformation::ResidueFactory::create_residue( *conformation::virtual_type_for_conf( *this ) );
				}

				Residue const & coordframe_rsd( residue( coordframe_pos ) );

				/// create some reasonable coords based on coordframe_rsd
				if ( coordframe_rsd.is_protein() ) {
					Vector const orig( coordframe_rsd.xyz("CA") );
					Vector const ihat( ( coordframe_rsd.xyz("N") - orig ).normalized() );
					Vector jhat( coordframe_rsd.xyz("C") - orig );
					jhat -= ihat.dot( jhat )*ihat;
					try {
						jhat.normalize();
					} catch( utility::excn::Exception & excn ) {
						throw CREATE_EXCEPTION(utility::excn::Exception, "Error in SymmetricConformation::recalculate transforms(): " + excn.msg() + "  Error occurred when trying to normalize coordinate vector for atom \"C\", residue " + std::to_string( coordframe_pos ) + "." );
					}
					runtime_assert( std::abs( ihat.dot( jhat ) ) < 1e-6 ); // should make this an assert
					vrt_res_op->set_xyz("ORIG", orig );
					vrt_res_op->set_xyz("X", orig + ihat );
					vrt_res_op->set_xyz("Y", orig + jhat );
				} else {
					/// could pretty easily add dna,rna,etc here
					utility_exit_with_message("not setup for non-protein residues as anchor-rsds for virtual-rsd-less symm poses");
				}
				vrt_res_cap = vrt_res_op.get();
			}
			runtime_assert( vrt_res_cap ); // ensure success

			// find the VRT that builds this res
			//core::Size vrt_ctrl = get_upstream_vrt( nresMonomer*(i-1) + 1 );
			Residue const &parent_vrt_res = *vrt_res_cap;

			xyzVector< core::Real > const &orig = parent_vrt_res.atom("ORIG").xyz();
			xyzVector< core::Real > X = parent_vrt_res.atom("X").xyz() - orig;
			xyzVector< core::Real > Y = parent_vrt_res.atom("Y").xyz() - orig;
			xyzVector< core::Real > Z = X.cross( Y );

			Tsymm_[comptag][isub] = SymmetryTransform( orig-Y,orig-Z, orig );
		}
	}
}


//fpd
void
SymmetricConformation::append_residue_by_jump(
	conformation::Residue const & new_rsd,
	Size const anchor_pos,
	std::string const& anchor_atom, // could be zero
	std::string const& root_atom, // ditto
	bool const start_new_chain // default false
)
{
	if ( start_new_chain ) {
		TR.Warning << "SymmetricConformation::append_residue_by_jump ignores start_new_chain" << std::endl;
	}

	core::Size nmonomer_jumps = symm_info_->get_njumps_subunit();
	core::Size nres_monomer = symm_info_->get_nres_subunit();
	core::Size nsubunits = symm_info_->subunits();
	core::Size asymm_anchor = ((anchor_pos-1)%nres_monomer) + 1;

	// add to the end of each subunit
	// transform to the coordinate frame of the scoring subunit
	// go from last->first so we don't have to worry about offsets
	for ( int i=nsubunits; i>=1; --i ) {
		core::Size seqpos = i*nres_monomer;
		core::Size my_anchor_pos = (i-1)*nres_monomer+asymm_anchor;
		Residue new_new_rsd = new_rsd;
		if ( !symm_info_->bb_is_independent( seqpos ) ) {
			// transform coords
			for ( int j=1; j<=(int)new_new_rsd.natoms(); ++j ) {
				new_new_rsd.set_xyz(j , apply_transformation( new_rsd.xyz(j), nres_monomer, seqpos ) );
			}
		}
		insert_residue_by_jump( new_new_rsd, seqpos+1, my_anchor_pos, anchor_atom, root_atom, start_new_chain );
	}
	// update symminfo
	symm_info_->resize_asu( nres_monomer + 1 );
	symm_info_->update_nmonomer_jumps( nmonomer_jumps + 1 );

	// update ft
	FoldTree f_in = core::conformation::symmetry::get_asymm_unit_fold_tree( *this );
	core::conformation::symmetry::symmetrize_fold_tree( *this, f_in );
	fold_tree( f_in );
}

void
SymmetricConformation::append_polymer_residue_after_seqpos(
	conformation::Residue const & new_rsd,
	Size const seqpos,
	bool const build_ideal_geometry // default false
)
{
	//core::Size nmonomer_jumps = symm_info_->get_njumps_subunit();
	core::Size nres_monomer = symm_info_->get_nres_subunit();
	core::Size nsubunits = symm_info_->subunits();
	core::Size asymm_anchor = ((seqpos-1)%nres_monomer) + 1;

	// add to the end of each subunit
	// transform to the coordinate frame of the scoring subunit
	// go from last->first so we don't have to worry about offsets
	for ( int i=nsubunits; i>=1; --i ) {
		core::Size final_residue_in_current_monomer = i*nres_monomer;
		core::Size symmetric_seqpos = (i-1)*nres_monomer+asymm_anchor;
		Residue new_new_rsd = new_rsd;
		if ( !symm_info_->bb_is_independent( final_residue_in_current_monomer ) ) {
			// transform coords
			for ( int j=1; j<=(int)new_new_rsd.natoms(); ++j ) {
				new_new_rsd.set_xyz(j , apply_transformation( new_rsd.xyz(j), nres_monomer, final_residue_in_current_monomer ) );
			}
		}
		Conformation::append_polymer_residue_after_seqpos( new_new_rsd, symmetric_seqpos, build_ideal_geometry );
	}
	// update symminfo
	symm_info_->resize_asu( nres_monomer + 1 );

	// update ft
	FoldTree f_in = core::conformation::symmetry::get_asymm_unit_fold_tree( *this );
	core::conformation::symmetry::symmetrize_fold_tree( *this, f_in );
	fold_tree( f_in );
}

void
SymmetricConformation::safely_append_polymer_residue_after_seqpos(
	conformation::Residue const & new_rsd,
	Size const seqpos,
	bool const build_ideal_geometry // default false
)
{
	//core::Size nmonomer_jumps = symm_info_->get_njumps_subunit();
	core::Size nres_monomer = symm_info_->get_nres_subunit();
	core::Size nsubunits = symm_info_->subunits();
	core::Size asymm_anchor = ((seqpos-1)%nres_monomer) + 1;

	// Remove N-terminal types of all the subunits
	for ( int i=nsubunits; i>=1; --i ) {
		core::Size symmetric_seqpos = (i-1)*nres_monomer+asymm_anchor;
		core::conformation::remove_upper_terminus_type_from_conformation_residue( *this, symmetric_seqpos );
	}

	append_polymer_residue_after_seqpos( new_rsd, seqpos, build_ideal_geometry );

}

void
SymmetricConformation::prepend_polymer_residue_before_seqpos(
	conformation::Residue const & new_rsd,
	Size const seqpos,
	bool const build_ideal_geometry // default false
)
{
	core::Size nres_monomer = symm_info_->get_nres_subunit();
	core::Size nsubunits = symm_info_->subunits();
	core::Size asymm_anchor = ((seqpos-1)%nres_monomer) + 1;

	// add to the end of each subunit
	// transform to the coordinate frame of the scoring subunit
	// go from last->first so we don't have to worry about offsets
	for ( int i=nsubunits; i>=1; --i ) {
		core::Size final_residue_in_current_monomer = i*nres_monomer;
		core::Size symmetric_seqpos = (i-1)*nres_monomer+asymm_anchor;

		ResidueOP new_new_rsd = new_rsd.clone();
		if ( !symm_info_->bb_is_independent( final_residue_in_current_monomer ) ) {
			// transform coords
			for ( int j=1; j<=(int)new_new_rsd->natoms(); ++j ) {
				new_new_rsd->set_xyz(j , apply_transformation( new_rsd.xyz(j), nres_monomer, final_residue_in_current_monomer ) );
			}
		}
		if ( build_ideal_geometry ) {
			Residue const & anchor_rsd( const_residue_( symmetric_seqpos ) );
			debug_assert( !anchor_rsd.is_lower_terminus() );
			orient_residue_for_ideal_bond(*new_new_rsd, new_new_rsd->upper_connect(), anchor_rsd, anchor_rsd.lower_connect(), *this );
		}
		// build_ideal_geometry always set to false because we handle it in this function
		Conformation::prepend_polymer_residue_before_seqpos( *new_new_rsd, symmetric_seqpos, false );
	}

	// update symminfo
	symm_info_->resize_asu( nres_monomer + 1 );
	// update ft
	FoldTree f_in = core::conformation::symmetry::get_asymm_unit_fold_tree( *this );
	core::conformation::symmetry::symmetrize_fold_tree( *this, f_in );
	fold_tree( f_in );

	if ( build_ideal_geometry ) {
		for ( int i=nsubunits; i>=1; --i ) {
			core::Size symmetric_seqpos = (i-1)*nres_monomer+asymm_anchor;
			rebuild_polymer_bond_dependent_atoms( symmetric_seqpos );
		}
	}


}

void
SymmetricConformation::safely_prepend_polymer_residue_before_seqpos(
	conformation::Residue const & new_rsd,
	Size const seqpos,
	bool const build_ideal_geometry // default false
)
{
	//core::Size nmonomer_jumps = symm_info_->get_njumps_subunit();
	core::Size nres_monomer = symm_info_->get_nres_subunit();
	core::Size nsubunits = symm_info_->subunits();
	core::Size asymm_anchor = ((seqpos-1)%nres_monomer) + 1;

	// Remove N-terminal types of all the subunits
	for ( int i=nsubunits; i>=1; --i ) {
		core::Size symmetric_seqpos = (i-1)*nres_monomer+asymm_anchor;
		core::conformation::remove_lower_terminus_type_from_conformation_residue( *this, symmetric_seqpos );
	}

	prepend_polymer_residue_before_seqpos( new_rsd, seqpos, build_ideal_geometry );

}

void
SymmetricConformation::insert_conformation_by_jump(
	Conformation const & new_conf,
	Size const insert_seqpos,
	Size const,
	Size const anchor_pos,
	Size const anchor_jump_number,
	std::string const & anchor_atom,
	std::string const & root_atom
) {
	if ( anchor_jump_number != 0 ) {
		TR.Warning << "SymmetricConformation::insert_conformation_by_jump ignores anchor_jump_number" << std::endl;
	}

	core::Size nmonomer_jumps = symm_info_->get_njumps_subunit();
	core::Size nres_monomer = symm_info_->get_nres_subunit();
	core::Size nsubunits = symm_info_->subunits();
	core::Size asymm_insert = ((insert_seqpos-2)%nres_monomer) + 2;  //?
	core::Size asymm_anchor = ((anchor_pos-1)%nres_monomer) + 1;

	// insert at the end of each subunit
	// transform to the coordinate frame of the scoring subunit
	// go from last->first so we don't have to worry about offsets
	for ( int i=nsubunits; i>=1; --i ) {
		core::Size insert_i = (i-1)*nres_monomer+asymm_insert;
		core::Size anchor_i = (i-1)*nres_monomer+asymm_anchor;

		ConformationOP new_new_conf = new_conf.clone();
		if ( !symm_info_->bb_is_independent( anchor_i ) ) {
			// transform coords
			for ( core::Size j=1; j<=new_new_conf->size(); ++j ) {
				for ( core::Size k=1; k<=new_new_conf->residue(j).natoms(); ++k ) {
					new_new_conf->set_xyz(core::id::AtomID(k,j) , apply_transformation( new_conf.residue(j).xyz(k), nres_monomer, anchor_i ) );
				}
			}
		}
		Conformation::insert_conformation_by_jump( *new_new_conf, insert_i, nmonomer_jumps+i, anchor_i, 0, anchor_atom, root_atom );
	}

	// update symminfo
	symm_info_->resize_asu( nres_monomer + new_conf.size() );
	symm_info_->update_nmonomer_jumps( nmonomer_jumps + new_conf.fold_tree().num_jump() + 1 );

	// update ft
	FoldTree f_in = core::conformation::symmetry::get_asymm_unit_fold_tree( *this );
	core::conformation::symmetry::symmetrize_fold_tree( *this, f_in );
	fold_tree( f_in );
}

/**
* @brief Detect existing disulfides from the protein structure.
* @details For full atom confomations, looks at SG-SG distance. If the SG-SG
*  are about 2.02 A apart, calls it a disulfide bond. For centroid and other
*  conformations, the less accurate CB-CB distance is used instead. In this
*  case a CB-CB distance of 3.72 A is optimal.
*/
void
SymmetricConformation::detect_disulfides( utility::vector1< std::pair<Size, Size> > const & disulf_bonds )
{
	using namespace utility::graph;
	using namespace basic::options;
	if ( ! option[ OptionKeys::in::detect_disulf ].user() || // if the option is not specified
			( option[ OptionKeys::in::detect_disulf ].user() && option[ OptionKeys::in::detect_disulf ]() )  // option specified and is true
			) {

		// gather all cys, construct mapping from resid to cys index
		utility::vector1< Size > resid_2_cysid( size(), 0 );
		Size num_cys( 0 );
		for ( Size ii = 1; ii <= size(); ++ii ) {
			if ( residue(ii).type().is_sidechain_thiol() || residue(ii).type().is_disulfide_bonded() ) {
				++num_cys;
				resid_2_cysid[ ii ] = num_cys;
			}
		}
		if ( num_cys == 0 ) return;

		// construct reverse mapping from cys index to resid
		utility::vector1< Size > cysid_2_resid( num_cys );
		for ( Size ii = 1; ii <= size(); ++ii ) {
			if ( resid_2_cysid[ ii ] != 0 ) {
				cysid_2_resid[ resid_2_cysid[ ii ]] = ii;
			}
		}

		// If all the cys are fullatom, use stricter criteria
		bool fullatom(true);
		for ( Size ii = 1; ii <= num_cys; ++ii ) {
			if ( residue_type(cysid_2_resid[ii]).mode()
					!= core::chemical::FULL_ATOM_t ) {
				fullatom = false;
				break;
			}
		}
		// SG-SG distance for fullatom, CB-CB distance otherwise
		Real const typical_disulfide_distance = fullatom? 2.02 : 3.72;
		Real const tolerance = option[OptionKeys::in::detect_disulf_tolerance].user()
			? option[OptionKeys::in::detect_disulf_tolerance]()
			: ( fullatom? 0.5 : 1.0 );

		// Create point graph
		PointGraphOP pg( new PointGraph );
		pg->set_num_vertices( num_cys );
		Distance maxrad( 0.0 );
		Distance maxd( typical_disulfide_distance + tolerance );
		for ( Size ii = 1; ii <= num_cys; ++ii ) {
			Residue const & ii_res = residue( cysid_2_resid[ ii ] );
			pg->get_vertex(ii).data().xyz() = ii_res.atoms()[ ii_res.nbr_atom() ].xyz();
			if ( ii_res.nbr_radius() > maxrad ) maxrad = ii_res.nbr_radius();
		}
		Distance neighbor_cutoff = maxrad + maxd;
		find_neighbors( pg, neighbor_cutoff );

		// Iterate across neighbors of cys residues; examine SG-SG distances.
		// Note that since graph only stores upper neighbors iterating through
		// the (upper) edge list will automatically prevent double counting.
		std::set< Size > processed_cys; // track cys that have already been processed
		// Process given pairs first if any are available

		for ( auto const & pair: disulf_bonds ) {

			Size ii_resid = pair.first;
			Size best_neighbor = pair.second;

			TR << "Adding disulfide between residues " << ii_resid << " " << best_neighbor << std::endl;

			add_disulfide_bond( ii_resid, best_neighbor );

			// mark both cys as processed
			processed_cys.insert( ii_resid );
			processed_cys.insert( best_neighbor );
		}

		// Now do everything that remains!
		for ( Size ii = 1; ii <= num_cys; ++ii ) {
			Size const ii_resid = cysid_2_resid[ ii ];
			if ( ii_resid > symm_info_->num_independent_residues() ) continue;
			//Size const ii_n_conn = residue( ii_resid ).type().n_possible_residue_connections();
			Residue const & ii_res( residue( ii_resid ) );
			//if ii already processed, continue
			if ( processed_cys.find( ii_resid) != processed_cys.end() ) {
				continue;
			}

			//Determine which atom makes the disulfide bond
			Size ii_sg_atomno(0);
			if ( ii_res.type().get_disulfide_atom_name() == "NONE" ) {
				TR.Error << "Can't find an atom to disulfide bond from at residue "<< ii_resid <<std::endl;
				utility_exit();
			} else {
				ii_sg_atomno = ii_res.type().atom_index( ii_res.type().get_disulfide_atom_name() );
			}
			Size ii_distance_atom_id = fullatom ? ii_sg_atomno : ii_res.atom_index( "CB" );


			Distance best_match( 0.0 );
			Size best_neighbor( 0 );
			//Size best_neighbor_cysid( 0 );

			for ( PointGraph::UpperEdgeListConstIter
					ii_iter = pg->get_vertex( ii ).upper_edge_list_begin(),
					ii_end_iter = pg->get_vertex( ii ).upper_edge_list_end();
					ii_iter != ii_end_iter; ++ii_iter ) {
				Size const jj = ii_iter->upper_vertex();

				Size const jj_resid = cysid_2_resid[ jj ];

				//TR << "looking for valid distance to res" << jj_resid << std::endl;

				Residue const & jj_res( residue( jj_resid ) );

				//if jj already processed, continue
				if ( processed_cys.find( jj_resid) != processed_cys.end() ) continue;

				Size jj_sg_atomno(0);
				if ( jj_res.type().get_disulfide_atom_name() == "NONE" ) {
					TR.Error << "Can't find an atom to disulfide bond from at residue "<< jj_resid <<std::endl;
					utility_exit();
				} else {
					jj_sg_atomno = jj_res.type().atom_index( jj_res.type().get_disulfide_atom_name() );
				}

				Size jj_distance_atom_id = fullatom ? jj_sg_atomno : jj_res.atom_index( "CB" );

				//TR << "distance between " << ii_distance_atom << " and " << jj_distance_atom << std::endl;
				Distance dist = ii_res.atom( ii_distance_atom_id ).xyz().distance( jj_res.atom( jj_distance_atom_id ).xyz() );
				if ( best_neighbor == 0 || dist < best_match ) {
					best_neighbor = jj_resid;
					best_match = dist;
				}
			}

			if ( best_neighbor == 0 || best_match >= typical_disulfide_distance + tolerance ) {
				// handle case where old disulfide doesn't exist anymore and
				// needs to be cleared
				if ( processed_cys.find( ii_resid ) == processed_cys.end() && ii_res.has_variant_type( chemical::DISULFIDE ) ) {
					TR << "Reverting out-of-date disulfide to thiol type at resid " << ii_resid << std::endl;

					bool const successful_revert = conformation::change_cys_state( ii_resid, "CYS", *this ) && !residue(ii_resid).has_variant_type( chemical::DISULFIDE );
					if ( !successful_revert ) {
						TR.Error << "unable to revert disulfide to thiol type for removal of disulfide at resid " << ii_resid << std::endl;
					}
				}

				// mark cys as processed
				processed_cys.insert( ii_resid );

				continue;

			} else { // found disulfide bond

				TR << "Found disulfide between residues " << ii_resid << " " << best_neighbor << std::endl;

				add_disulfide_bond( ii_resid, best_neighbor );

				// mark both cys as processed
				processed_cys.insert( ii_resid );
				processed_cys.insert( best_neighbor );
			}
		}
	} // use option "-detect_disulf" to control this
}

/// @brief Declare that a chemical bond exists between two residues
/// @details This updates all symmetry copies, so that each one has a chemical
/// bond between the residues in question.
/// @author Frank DiMaio
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
void
SymmetricConformation::declare_chemical_bond(
	Size const seqpos1,
	std::string const & atom_name1,
	Size const seqpos2,
	std::string const & atom_name2
)
{
	numeric::xyzVector< core::Real > const xform1( residue(seqpos1).xyz(1) ); //Transform for residue 1
	numeric::xyzVector< core::Real > const xform2( residue(seqpos2).xyz(1) ); //Transform for residue 2

	//We need to find all pairs (seqpos1', seqpos2') of residues that have the same transform when transformed into the old pair's coordinate frame.
	core::Size const nclones( symm_info_->num_bb_clones() + 1 );
	core::Size const masterpos1( symm_info_->bb_is_independent( seqpos1 ) ? seqpos1 : symm_info_->bb_follows( seqpos1 ) );
	core::Size const masterpos2( symm_info_->bb_is_independent( seqpos2 ) ? seqpos2 : symm_info_->bb_follows( seqpos2 ) );
	SymmetryInfo::Clones clones1( symm_info_->bb_clones(masterpos1) );
	clones1.push_back(masterpos1);
	SymmetryInfo::Clones clones2( symm_info_->bb_clones(masterpos2) );
	clones2.push_back(masterpos2);

	//Now we have lists of all the positions equivalent to seqpos1 (including seqpos1), and all the positions equivalent to seqpos2 (including seqpos2).
	//Now we iterate over both lists:
	for ( core::Size i=1; i<=nclones; ++i ) {
		bool found(false);
		for ( core::Size j=1; j<=nclones; ++j ) {
			numeric::xyzVector < core::Real > const xform3( residue(clones1[i]).xyz(1) - residue(clones2[j]).xyz(1) ); //Difference of transforms of residue 1prime and 2prime.
			numeric::xyzVector < core::Real > const xform4( apply_transformation( xform1, seqpos1, clones1[i] ) - apply_transformation( xform2, seqpos1, clones1[i] ) ); //Difference of (transform of residue 1 in residue 1prime's frame) and (transform of residue2 in residue 1prime's frame).
			core::Real const dist( ( xform3 - xform4 ).length() ); //Take the difference of the two inter-residue transforms.
			if ( dist < 1e-4 ) {
				TR << "Add symmetric chemical bond " << clones1[i] << " to " << clones2[j] << std::endl;
				Conformation::declare_chemical_bond( clones1[i], atom_name1, clones2[j], atom_name2 );
				found=true;
				break;
			}
		}
		if ( !found ) {
			TR << "Error in declare_chemical_bond: unable to find corresponding partner!" << std::endl;
		}
	}

}

void
SymmetricConformation::update_noncanonical_connection(
	Size const lower_seqpos,
	Size const lr_conn_id,
	Size const upper_seqpos,
	Size const ur_conn_id)
{
	if ( lower_seqpos < 1 || lower_seqpos > size() ) return;
	if ( upper_seqpos < 1 || upper_seqpos > size() ) return;

	//We need to find all pairs (seqpos1', seqpos2') of residues that have the same transform when transformed into the old pair's coordinate frame.

	if ( ur_conn_id == 0 || ur_conn_id > residue(upper_seqpos).connect_map_size() ) {
		TR.Error << "Trying to update a noncanonical connection from " << lower_seqpos << " conn " << lr_conn_id
			<< " to " << upper_seqpos << " conn " << ur_conn_id << std::endl;
		TR.Error << "but the latter residue doesn't have that many connections... Something is out of date." << std::endl;
		return;
	}

	std::vector < std::pair < Size, Size > > res_pairs = symm_info_->map_symmetric_res_pairs( lower_seqpos, upper_seqpos, *this);
	res_pairs.push_back(std::make_pair(lower_seqpos, upper_seqpos));
	for ( auto res_pair : res_pairs ) {
		core::Size lower = res_pair.first;
		core::Size upper = res_pair.second;
		set_noncanonical_connection( lower, lr_conn_id, upper, ur_conn_id );
	}
}


}
} // conformation
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::symmetry::SymmetricConformation::save( Archive & arc ) const {
	arc( cereal::base_class< core::conformation::Conformation >( this ) );
	arc( CEREAL_NVP( symm_info_ ) ); // SymmetryInfoOP
	arc( CEREAL_NVP( Tsymm_ ) ); // std::map<char, utility::vector1<numeric::HomogeneousTransform<core::Real> > >
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::symmetry::SymmetricConformation::load( Archive & arc ) {
	arc( cereal::base_class< core::conformation::Conformation >( this ) );
	arc( symm_info_ ); // SymmetryInfoOP
	arc( Tsymm_ ); // std::map<char, utility::vector1<numeric::HomogeneousTransform<core::Real> > >
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::symmetry::SymmetricConformation );
CEREAL_REGISTER_TYPE( core::conformation::symmetry::SymmetricConformation )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_symmetry_SymmetricConformation )
#endif // SERIALIZATION
