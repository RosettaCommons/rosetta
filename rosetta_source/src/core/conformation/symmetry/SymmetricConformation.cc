// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/conformation/symmetry/SymmetricConformation.hh
/// @brief  symmetry conformation container.
//					Contains overloaded functions needed to
//					make changes in conformation symmetric
/// @author Phil Bradley, Ingemar Andre

// Unit headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/id/TorsionID.hh>
#include <numeric/xyz.functions.hh>



namespace core {
namespace conformation {
namespace symmetry {

static basic::Tracer TR("core.conformation.symmetry.Conformation");

/// @brief  Default CTOR
SymmetricConformation::SymmetricConformation():
		Conformation(),
		symm_info_()
{
	Tsymm_.clear();
}

	/// @brief  Default CTOR
SymmetricConformation::SymmetricConformation( Conformation const & conf, SymmetryInfo const & symm_info ):
		Conformation( conf ),
		symm_info_( symm_info.clone() )
{
	Tsymm_.clear();
}

/// @brief copy constructor
SymmetricConformation::SymmetricConformation( SymmetricConformation const & src ):
	Conformation( src ),
	symm_info_( src.symm_info_->clone() ),
	Tsymm_( src.Tsymm_ )
{}

/// @brief operator=
Conformation &
SymmetricConformation::operator=( SymmetricConformation const & src )
{
	// will this work?
	Conformation::operator=( src );
	symm_info_ = src.symm_info_->clone();
	Tsymm_ = src.Tsymm_;
	return *this;
}


///@details make a copy of this conformation( allocate actual memory for it )
ConformationOP
SymmetricConformation::clone() const
{
  return new SymmetricConformation( *this );
}

bool
SymmetricConformation::same_type_as_me( Conformation const & other, bool recurse  /* = true */ ) const
{
   SymmetricConformation const * symm_other = dynamic_cast< SymmetricConformation const * > ( &other);
   if ( ! symm_other ) {
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
	typedef SymmetryInfo::DOF_IDs DOF_IDs;

	core::Size parent_rsd;

	if ( !symm_info_->dof_is_independent( id, *this ) ) {
		TR.Debug << "SymmetricConformation:: directly setting a dependent DOF!, try to set its parent" << std::endl;
		if (id.type() >= id::RB1 && id.type() <= id::RB6) {
			int parent_jump = symm_info_->jump_follows( fold_tree().get_jump_that_builds_residue( id.rsd() ) );
			parent_rsd = fold_tree().downstream_jump_residue( parent_jump );

			// clear cached transforms  ... is this necessary??
			Tsymm_.clear();
		} else {
			parent_rsd = symm_info_->bb_follows( id.rsd() ) ;
		}
	} else {
		parent_rsd = id.rsd();
	}

	id::DOF_ID parent_id ( id::AtomID( id.atomno(), parent_rsd ), id.type() );

	// set this DOF using base-class method
	Conformation::set_dof( parent_id, setting );
	{
		DOF_IDs const & dofs( symm_info_->dependent_dofs( parent_id, *this ) );
		for ( DOF_IDs::const_iterator dof =dofs.begin(), dofe= dofs.end(); dof != dofe; ++dof ) {
			Conformation::set_dof( *dof, setting );
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
		for ( SymmetryInfo::Clones::const_iterator pos=symm_info_->bb_clones( seqpos ).begin(),
			epos=symm_info_->bb_clones( seqpos ).end(); pos != epos; ++pos ) {
			Conformation::set_secstruct( *pos, setting );
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
	typedef SymmetryInfo::TorsionIDs TorsionIDs;

	TR.Trace << "SymmetricConformation: set_torsion: " << id << ' ' << setting << std::endl;

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
		for ( TorsionIDs::const_iterator tor =tors.begin(), tore= tors.end(); tor != tore; ++tor ) {
			Conformation::set_torsion( *tor, setting );
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
	Real const setting
)
{
	AtomID parent_atom1, parent_atom2, parent_atom3, parent_atom4;

	TR.Trace << "SymmetricConformation: set_torsion_angle: " << atom1 << "+ to " << setting << std::endl;

	// assume if 1st atom is dependent, all are
	if ( !symm_info_->bb_is_independent( atom1.rsd() ) ) {
		TR.Debug << "SymmetricConformation:: directly setting a dependent ANGLE!, try to set its parent" << std::endl;
		parent_atom1 = id::AtomID( atom1.atomno(), symm_info_->bb_follows( atom1.rsd() ) );
		parent_atom2 = id::AtomID( atom2.atomno(), symm_info_->bb_follows( atom2.rsd() ) );
		parent_atom3 = id::AtomID( atom3.atomno(), symm_info_->bb_follows( atom3.rsd() ) );
		parent_atom4 = id::AtomID( atom3.atomno(), symm_info_->bb_follows( atom3.rsd() ) );
	} else {
		parent_atom1=atom1;
		parent_atom2=atom2;
		parent_atom3=atom3;
		parent_atom4=atom4;
	}

	// set this DOF using base-class method
	Conformation::set_torsion_angle( parent_atom1, parent_atom2, parent_atom3, parent_atom4, setting );

	{
		core::Size nclones = symm_info_->num_bb_clones( );
		SymmetryInfo::Clones clones1 = symm_info_->bb_clones( parent_atom1.rsd() );
		SymmetryInfo::Clones clones2 = symm_info_->bb_clones( parent_atom2.rsd() );
		SymmetryInfo::Clones clones3 = symm_info_->bb_clones( parent_atom3.rsd() );
		SymmetryInfo::Clones clones4 = symm_info_->bb_clones( parent_atom4.rsd() );
		AtomID dep_atom1, dep_atom2, dep_atom3, dep_atom4;
		for (int i=1; i<=(int)nclones; ++i) {
			dep_atom1 = id::AtomID( atom1.atomno(), clones1[i] );
			dep_atom2 = id::AtomID( atom2.atomno(), clones2[i] );
			dep_atom3 = id::AtomID( atom3.atomno(), clones3[i] );
			dep_atom4 = id::AtomID( atom4.atomno(), clones4[i] );
			Conformation::set_torsion_angle( dep_atom1, dep_atom2, dep_atom3, dep_atom4, setting );
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
		for (int i=1; i<=(int)nclones; ++i) {
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
		AtomID dep_atom1, dep_atom2, dep_atom3;
		for (int i=1; i<=(int)nclones; ++i) {
			dep_atom1 = id::AtomID( atom1.atomno(), clones1[i] );
			dep_atom2 = id::AtomID( atom2.atomno(), clones2[i] );
			Conformation::set_bond_length( dep_atom1, dep_atom2, setting );
		}
	}
}


	/// JUMPS
void
SymmetricConformation::set_jump( int const jump_number, kinematics::Jump const & new_jump )
{
	typedef SymmetryInfo::AtomIDs AtomIDs;
	typedef SymmetryInfo::Clones Clones;

	TR.Trace << "SymmetricConformation: set_jump jump_number: " << jump_number << std::endl;

	// clear cached transforms
	Tsymm_.clear();

	id::AtomID const id( Conformation::jump_atom_id( jump_number ) );
	Conformation::set_jump( id, new_jump );

	if ( !symm_info_->jump_is_independent( jump_number ) ) {
		TR.Warning << "SymmetricConformation:: directly setting a dependent ATOM!" << std::endl;
		TR.Warning << "the jump " << jump_number << " is controlled by " << symm_info_->jump_follows( jump_number ) << std::endl;
	} else {
		for ( Clones::const_iterator pos= symm_info_->jump_clones( jump_number ).begin(),
		      epos=symm_info_->jump_clones( jump_number ).end(); pos != epos; ++pos ) {
			id::AtomID const id_clone( Conformation::jump_atom_id( *pos ) );
			Conformation::set_jump( id_clone, new_jump );
		}
	}
}

void
SymmetricConformation::set_jump_now( int const jump_number, kinematics::Jump const & new_jump )
{
	typedef SymmetryInfo::AtomIDs AtomIDs;
	typedef SymmetryInfo::Clones Clones;

	TR.Trace << "SymmetricConformation: set_jump jump_number: " << jump_number << std::endl;

	// clear cached transforms
	Tsymm_.clear();

	id::AtomID const id( Conformation::jump_atom_id( jump_number ) );
	Conformation::set_jump_now( jump_number, new_jump );

	if ( !symm_info_->jump_is_independent( jump_number ) ) {
		TR.Warning << "SymmetricConformation:: directly setting a dependent ATOM!" << std::endl;
	} else {
		for ( Clones::const_iterator pos= symm_info_->jump_clones( jump_number ).begin(),
		      epos=symm_info_->jump_clones( jump_number ).end(); pos != epos; ++pos ) {
			id::AtomID const id_clone( Conformation::jump_atom_id( *pos ) );
			Conformation::set_jump( id_clone, new_jump );
		}
	}
}

// This doesn't work with
void
SymmetricConformation::set_jump( id::AtomID const & id, kinematics::Jump const & new_jump )
{
	typedef SymmetryInfo::AtomIDs AtomIDs;
	typedef SymmetryInfo::Clones Clones;

	TR.Trace << "SymmetricConformation: set_jump id:" << id << std::endl;

	// clear cached transforms
	Tsymm_.clear();

	Conformation::set_jump( id, new_jump );

	int const jump_number ( fold_tree().get_jump_that_builds_residue( id.rsd() ) );
 	if ( !symm_info_->jump_is_independent( jump_number ) ) {
   TR.Warning << "SymmetricConformation:: directly setting a dependent ATOM!" << std::endl;
	} else {
		for ( Clones::const_iterator pos= symm_info_->jump_clones( jump_number ).begin(),
		      epos=symm_info_->jump_clones( jump_number ).end(); pos != epos; ++pos ) {
			id::AtomID const id_clone( Conformation::jump_atom_id( *pos ) );
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
		if (!orient_backbone) {
			// make the new res
			for (int i=1; i<=(int)new_rsd.natoms(); ++i) {
				new_rsd.set_xyz(i , apply_transformation( new_rsd_in.xyz(i), seqpos, parent_rsd ) );
			}
		}
	} else {
		parent_rsd = seqpos;
	}

	// set this residue using base-class method
	Conformation::replace_residue( parent_rsd, new_rsd, orient_backbone );

	// now the copies
	for ( SymmetryInfo::Clones::const_iterator pos=symm_info_->bb_clones( parent_rsd ).begin(),
	      epos=symm_info_->bb_clones( parent_rsd ).end(); pos != epos; ++pos ) {
		// if we're not orienting backbone, update symm copies using local coord frame for each
		//fpd --> we already oriented backbone for master, that should be good enough
		if (!orient_backbone) {
			// make the new res
			Residue new_new_rsd = new_rsd;
			for (int i=1; i<=(int)new_new_rsd.natoms(); ++i) {
				new_new_rsd.set_xyz(i , apply_transformation( new_rsd.xyz(i), parent_rsd, *pos ) );
			}
			Conformation::replace_residue( *pos, new_new_rsd, false );
		} else {
			Conformation::replace_residue( *pos, new_rsd, orient_backbone );
		}
	}
}

/// @details symmetry-safe replace residue
void
SymmetricConformation::replace_residue( Size const seqpos,
                 Residue const & new_rsd,
                 utility::vector1< std::pair< std::string, std::string > > const & atom_pairs )  {
	TR.Debug << "SymmetricConformation: replace_residue: " << seqpos << std::endl;

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
	for ( SymmetryInfo::Clones::const_iterator pos=symm_info_->bb_clones( parent_rsd ).begin(),
	      epos=symm_info_->bb_clones( parent_rsd ).end(); pos != epos; ++pos ) {
		//fpd Same logic as above; we already oriented backbone for master, that should be good enough
		// make the new res
		Residue new_new_rsd = new_rsd;
		for (int i=1; i<=(int)new_new_rsd.natoms(); ++i) {
			new_new_rsd.set_xyz(i , apply_transformation( new_rsd.xyz(i), parent_rsd, *pos ) );
		}
		Conformation::replace_residue( *pos, new_new_rsd, false );
	}
}


core::Size
SymmetricConformation::get_upstream_vrt( Size seqpos ) const {
	if (this->residue( seqpos ).aa() == core::chemical::aa_vrt)
		return seqpos;

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
SymmetricConformation::get_transformation( core::Size resid ) {
	if (Tsymm_.size() != symm_info_->subunits())
		recalculate_transforms( );

	return Tsymm_[ symm_info_->subunit_index( resid ) ];
}

//  Remap coordinate X from resid i to j
PointPosition
SymmetricConformation::apply_transformation(
		PointPosition Xin,
		core::Size residfrom,
		core::Size residto ) {
	if (Tsymm_.size() != symm_info_->subunits())
		recalculate_transforms( );

	numeric::HomogeneousTransform< core::Real > const &Tsymm_from
		= Tsymm_[ symm_info_->subunit_index( residfrom ) ];
	numeric::HomogeneousTransform< core::Real > const &Tsymm_to
		= Tsymm_[ symm_info_->subunit_index( residto ) ];

	PointPosition Xout;
	Xout = (Tsymm_from.inverse() * Xin);
	Xout = (Tsymm_to * Xout);
	// TR.Debug << "XFORM " << residfrom << " to " << residto << std::endl;
	// TR.Debug << "                 Xin =  (" << Xin[0] << "," << Xin[1] << "," << Xin[2] << ")" << std::endl;
	// TR.Debug << "Tto * Tfrom^-1 * Xin =  (" << Xout[0] << "," << Xout[1] << "," << Xout[2] << ")" << std::endl;
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
	for ( SymmetryInfo::Clones::const_iterator pos=symm_info_->bb_clones( parent_id.rsd() ).begin(),
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
	for (int i=1; i<=(int)ids.size(); ++i) {
		AtomID id=ids[i], parent_id;
		PointPosition position=positions[i], parent_pos;

		if ( id.rsd() > symm_info_->num_total_residues_without_pseudo() ) {
			ids_with_symm.push_back( id );
			positions_with_symm.push_back( position );
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
			for ( SymmetryInfo::Clones::const_iterator pos=symm_info_->bb_clones( parent_id.rsd() ).begin(),
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

// recalculate the Tsymm_ transforms using the current pose
void
SymmetricConformation::recalculate_transforms( ) {
	using namespace numeric;

	// clear current xforms
	Tsymm_.clear();

	// rebuild based on virtuals
	core::Size nsubunits = symm_info_->subunits();
	core::Size nresMonomer = symm_info_->num_independent_residues();

	for (int i=1; i<=(int)nsubunits; ++i) {
		// find the VRT that builds this res
		core::Size vrt_ctrl = get_upstream_vrt( nresMonomer*(i-1) + 1 );
		Residue const &parent_vrt_res = residue( vrt_ctrl );

		xyzVector< core::Real > const &orig = parent_vrt_res.atom("ORIG").xyz();
		xyzVector< core::Real > X = parent_vrt_res.atom("X").xyz() - orig;
		xyzVector< core::Real > Y = parent_vrt_res.atom("Y").xyz() - orig;
		xyzVector< core::Real > Z = X.cross( Y );

		//TR.Debug << "[TRANSFORM " << i << "/" << vrt_ctrl << "]" << std::endl
		//         << "     X =  (" << X[0] << "," << X[1] << "," << X[2] << ")" << std::endl
		// 				 << "     Y =  (" << Y[0] << "," << Y[1] << "," << Y[2] << ")" << std::endl
		// 				 << "     Z =  (" << Z[0] << "," << Z[1] << "," << Z[2] << ")" << std::endl
		// 				 << "  orig =  (" << orig[0] << "," << orig[1] << "," << orig[2] << ")" << std::endl;

		// we can't assume the VRT is "proper" (orthogonal and normal)
		Tsymm_.push_back( HomogeneousTransform< core::Real>( orig-Y,orig-Z, orig ) );
	}
}

}
} // conformation
} // core


