// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/symmetry/MirrorSymmetricConformation.cc
/// @brief  Symmetry conformation container.
/// @details Contains overloaded functions needed to make changes in conformation symmetric and
/// to handle mirror symmetry operations.
/// @author Frank DiMaio
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit headers
#include <core/conformation/symmetry/MirrorSymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
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
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
//#include <core/chemical/AtomType.hh>
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

static THREAD_LOCAL basic::Tracer TR( "core.conformation.symmetry.Conformation" );

/// @brief  Default CTOR
MirrorSymmetricConformation::MirrorSymmetricConformation():
	SymmetricConformation()
{
	jump_is_mirrored_.clear();
	res_is_mirrored_.clear();
}

/// @brief  Default CTOR
MirrorSymmetricConformation::MirrorSymmetricConformation(
	Conformation const & conf,
	SymmetryInfo const & symm_info
):
	SymmetricConformation(conf, symm_info)
{
	update_njumps_nres();
	//update_residue_identities();
}

MirrorSymmetricConformation::MirrorSymmetricConformation(
	MirrorSymmetricConformation const & conf
):
	SymmetricConformation(conf, *conf.Symmetry_Info())
{
	res_is_mirrored_=conf.res_is_mirrored_;
	jump_is_mirrored_=conf.jump_is_mirrored_;
}

Conformation &
MirrorSymmetricConformation::operator=( Conformation const & src )
{
	MirrorSymmetricConformation const * mirr_sym_conf = dynamic_cast< MirrorSymmetricConformation const * > ( &src );
	if ( mirr_sym_conf ) {
		MirrorSymmetricConformation::operator=( *mirr_sym_conf );
	} else {
		throw utility::excn::EXCN_Msg_Exception( "MirrorSymmetricConformation::operator= was handed a non-MirrorSymmetricConformation" );
	}

	return *this;
}

void
MirrorSymmetricConformation::detached_copy( Conformation const & src ) {
	SymmetricConformation::detached_copy( src );
	MirrorSymmetricConformation const * mirr_sym_conf = dynamic_cast< MirrorSymmetricConformation const * > ( &src );
	if ( mirr_sym_conf ) {
		// Copy the private data members of the symmetric conformation.
		res_is_mirrored_ = mirr_sym_conf->res_is_mirrored_;
		jump_is_mirrored_ = mirr_sym_conf->jump_is_mirrored_;
	} else {
		throw utility::excn::EXCN_Msg_Exception( "MirrorSymmetricConformation detached_copy was handed a non-MirrorSymmetricConformation" );
	}
}


/// @details make a copy of this conformation( allocate actual memory for it )
ConformationOP
MirrorSymmetricConformation::clone() const
{
	MirrorSymmetricConformationOP copy( new MirrorSymmetricConformation(*this) );
	return copy;
}

bool
MirrorSymmetricConformation::same_type_as_me( Conformation const & other, bool recurse  /* = true */ ) const
{
	if ( ! dynamic_cast< MirrorSymmetricConformation const * > ( &other) ) {
		return false;
	}
	if ( recurse ) {
		return other.same_type_as_me( *this, false );
	} else {
		return true;
	}
}

MirrorSymmetricConformation::~MirrorSymmetricConformation() { clear(); }

/// DOF
void
MirrorSymmetricConformation::set_dof( DOF_ID const & id, Real const setting )
{
	typedef SymmetryInfo::DOF_IDs DOF_IDs;

	core::Size parent_rsd;
	//bool parent_mirrored=false;
	if ( Symmetry_Info()->torsion_changes_move_other_monomers() ) {
		SymmetricConformation::clear_Tsymm( );
	}

	// if parent is mirrored, we need to apply mirror op to jump
	if ( !Symmetry_Info()->dof_is_independent( id, *this ) ) {
		TR.Debug << "MirrorSymmetricConformation:: directly setting a dependent DOF!, trying to set its parent" << std::endl;
		if ( id.type() >= id::RB1 && id.type() <= id::RB6 ) {
			TR << "MirrorSymmetricConformation:: directly setting a dependent jump!  Failing." << std::endl;
			return;
		} else {
			parent_rsd = Symmetry_Info()->bb_follows( id.rsd() ) ;
		}
	} else {
		parent_rsd = id.rsd();
	}

	if ( id.type() >= id::RB1 && id.type() <= id::RB6 ) {
		SymmetricConformation::clear_Tsymm( );
	}

	id::DOF_ID parent_id ( id::AtomID( id.atomno(), parent_rsd ), id.type() );

	// set this DOF using base-class method
	Conformation::set_dof( parent_id, setting );

	DOF_IDs const & dofs( Symmetry_Info()->dependent_dofs( parent_id, *this ) );
	for ( DOF_IDs::const_iterator dof =dofs.begin(), dofe= dofs.end(); dof != dofe; ++dof ) {
		Conformation::set_dof( *dof, setting );
	}
}

/// TORSIONS
void
MirrorSymmetricConformation::set_torsion( id::TorsionID const & id, Real const setting )
{
	typedef SymmetryInfo::TorsionIDs TorsionIDs;

	TR.Trace << "MirrorSymmetricConformation: set_torsion: " << id << ' ' << setting << std::endl;
	if ( Symmetry_Info()->torsion_changes_move_other_monomers() ) {
		SymmetricConformation::clear_Tsymm( );
	}

	core::Size parent_rsd;
	//bool parent_mirrored=false;

	if ( !Symmetry_Info()->torsion_is_independent( id ) ) {
		TR.Debug << "MirrorSymmetricConformation:: directly setting a dependent TORSION!, trying to set its parent" << std::endl;
		//parent_mirrored = res_is_mirrored_[id.rsd()];
		parent_rsd = Symmetry_Info()->bb_follows( id.rsd() ) ;
	} else {
		parent_rsd = id.rsd();
	}

	id::TorsionID parent_id ( parent_rsd, id.type(), id.torsion() );

	Conformation::set_torsion( parent_id, setting );
	{
		TorsionIDs const & tors( Symmetry_Info()->dependent_torsions( parent_id ) );
		for ( TorsionIDs::const_iterator tor =tors.begin(), tore= tors.end(); tor != tore; ++tor ) {
			Conformation::set_torsion( *tor, setting );
		}
	}
}

/// TORSION ANGLES
void
MirrorSymmetricConformation::set_torsion_angle(
	AtomID const & atom1,
	AtomID const & atom2,
	AtomID const & atom3,
	AtomID const & atom4,
	Real const setting,
	bool const quiet
)
{
	AtomID parent_atom1, parent_atom2, parent_atom3, parent_atom4;

	TR.Trace << "MirrorSymmetricConformation: set_torsion_angle: " << atom1 << "+ to " << setting << std::endl;
	if ( Symmetry_Info()->torsion_changes_move_other_monomers() ) {
		SymmetricConformation::clear_Tsymm( );
	}

	// assume if 1st atom is dependent, all are
	//bool parent_mirrored=false;
	if ( !Symmetry_Info()->bb_is_independent( atom1.rsd() ) ) {
		TR.Debug << "MirrorSymmetricConformation:: directly setting a dependent ANGLE!, try to set its parent" << std::endl;
		//parent_mirrored = res_is_mirrored_[atom1.rsd()];
		parent_atom1 = id::AtomID( atom1.atomno(), Symmetry_Info()->bb_follows( atom1.rsd() ) );
		parent_atom2 = id::AtomID( atom2.atomno(), Symmetry_Info()->bb_follows( atom2.rsd() ) );
		parent_atom3 = id::AtomID( atom3.atomno(), Symmetry_Info()->bb_follows( atom3.rsd() ) );
		parent_atom4 = id::AtomID( atom4.atomno(), Symmetry_Info()->bb_follows( atom4.rsd() ) );
	} else {
		parent_atom1=atom1;
		parent_atom2=atom2;
		parent_atom3=atom3;
		parent_atom4=atom4;
	}

	Conformation::set_torsion_angle( parent_atom1, parent_atom2, parent_atom3, parent_atom4, setting, quiet );

	{
		core::Size nclones = Symmetry_Info()->num_bb_clones( );
		SymmetryInfo::Clones clones1 = Symmetry_Info()->bb_clones( parent_atom1.rsd() );
		SymmetryInfo::Clones clones2 = Symmetry_Info()->bb_clones( parent_atom2.rsd() );
		SymmetryInfo::Clones clones3 = Symmetry_Info()->bb_clones( parent_atom3.rsd() );
		SymmetryInfo::Clones clones4 = Symmetry_Info()->bb_clones( parent_atom4.rsd() );
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



/// JUMPS
void
MirrorSymmetricConformation::set_jump( int const jump_number, kinematics::Jump const & new_jump )
{
	//typedef SymmetryInfo::AtomIDs AtomIDs;
	typedef SymmetryInfo::Clones Clones;

	TR.Debug << "MirrorSymmetricConformation: set_jump jump_number: " << jump_number << std::endl;

	// clear cached transforms
	SymmetricConformation::clear_Tsymm( );

	id::AtomID const id( Conformation::jump_atom_id( jump_number ) );
	Conformation::set_jump( id, new_jump );

	if ( !Symmetry_Info()->jump_is_independent( jump_number ) ) {
		TR.Warning << "MirrorSymmetricConformation:: directly setting a dependent jump!" << std::endl;
		TR.Warning << "the jump " << jump_number << " is controlled by " << Symmetry_Info()->jump_follows( jump_number ) << std::endl;
	} else {
		for ( Clones::const_iterator pos= Symmetry_Info()->jump_clones( jump_number ).begin(),
				epos=Symmetry_Info()->jump_clones( jump_number ).end(); pos != epos; ++pos ) {
			id::AtomID const id_clone( Conformation::jump_atom_id( *pos ) );

			//fpd the logic is a bit convoluted here
			//    to make the inverse jump we need not just the rotation & translation
			//    but also the upstream stub
			//fpd this code could be made cleaner if jump stored the upstream stub positions=
			//    but that is probably inefficient?
			kinematics::Jump new_jump_copy = new_jump;
			new_jump_copy.set_invert(
				jump_is_mirrored_[*pos].first,
				jump_is_mirrored_[*pos].second
			);
			Conformation::set_jump( id_clone, new_jump_copy );
		}
	}
}

void
MirrorSymmetricConformation::set_jump( id::AtomID const & id, kinematics::Jump const & new_jump )
{
	typedef SymmetryInfo::Clones Clones;

	TR.Debug << "MirrorSymmetricConformation: set_jump id:" << id << std::endl;

	// clear cached transforms
	SymmetricConformation::clear_Tsymm( );

	Conformation::set_jump( id, new_jump );

	int const jump_number ( Conformation::fold_tree().get_jump_that_builds_residue( id.rsd() ) );
	if ( !Symmetry_Info()->jump_is_independent( jump_number ) ) {
		TR.Warning << "MirrorSymmetricConformation:: directly setting a dependent jump!" << std::endl;
	} else {
		for ( Clones::const_iterator pos= Symmetry_Info()->jump_clones( jump_number ).begin(),
				epos=Symmetry_Info()->jump_clones( jump_number ).end(); pos != epos; ++pos ) {
			id::AtomID const id_clone( Conformation::jump_atom_id( *pos ) );

			//fpd the logic is a bit convoluted here
			//    to make the inverse jump we need not just the rotation & translation
			//    but also the upstream stub
			//fpd this code could be made cleaner if jump stored the upstream stub positions=
			//    but that is probably inefficient?
			kinematics::Jump new_jump_copy = new_jump;
			new_jump_copy.set_invert(
				jump_is_mirrored_[*pos].first,
				jump_is_mirrored_[*pos].second
			);
			Conformation::set_jump( id_clone, new_jump_copy );
		}
	}
}

core::Real
MirrorSymmetricConformation::get_residue_weight(core::Size resid1, core::Size resid2) const {
	return Symmetry_Info()->score_multiply(resid1,resid2);
}

/// @details mirror-symmetry-safe replace residue
void
MirrorSymmetricConformation::replace_residue( Size const seqpos, Residue const & new_rsd_in, bool const orient_backbone )
{
	TR.Debug << "MirrorSymmetricConformation: replace_residue: " << seqpos << std::endl;

	core::Size parent_rsd;
	ResidueOP new_rsd( new_rsd_in.clone() );

	if ( !Symmetry_Info()->bb_is_independent( seqpos ) ) {
		TR.Debug << "MirrorSymmetricConformation::replace_residue(2) directly setting a dependent TORSION!, try to set its parent" << std::endl;

		parent_rsd = Symmetry_Info()->bb_follows( seqpos ) ;

		// unlike non-mirror symmetries, call this function _even if_ we are orienting the backbone
		//  --> this ensures residues are properly mirrored D<->L
		for ( int i=1; i<=(int)new_rsd->natoms(); ++i ) {
			new_rsd->set_xyz(i , apply_transformation( new_rsd_in.xyz(i), seqpos, parent_rsd ) );
		}
		if ( res_is_mirrored_[seqpos] ) {
			flip_chirality( new_rsd );
		}
	} else {
		parent_rsd = seqpos;
	}

	// set this residue using base-class method
	Conformation::replace_residue( parent_rsd, *new_rsd, orient_backbone );

	// now the copies
	for ( SymmetryInfo::Clones::const_iterator pos=Symmetry_Info()->bb_clones( parent_rsd ).begin(),
			epos=Symmetry_Info()->bb_clones( parent_rsd ).end(); pos != epos; ++pos ) {
		ResidueOP new_new_rsd( new_rsd->clone() );
		for ( int i=1; i<=(int)new_new_rsd->natoms(); ++i ) {
			Vector newpos = apply_transformation( new_rsd->xyz(i), parent_rsd, *pos );
			new_new_rsd->set_xyz(i , newpos );
		}

		if ( res_is_mirrored_[*pos] ) {
			flip_chirality( new_new_rsd );
		}

		//fpd: replace_residue works a bit different than other Conformation functions
		//     the residue is created rooted on a jump atom
		//     then it is stitched on the base conformation
		Conformation::replace_residue( *pos, *new_new_rsd, false );
	}

	//fpd may have implicitly changed a jump
	SymmetricConformation::clear_Tsymm( );
}

/// @details symmetry-safe replace residue
void
MirrorSymmetricConformation::replace_residue( Size const seqpos,
	Residue const & new_rsd_in,
	utility::vector1< std::pair< std::string, std::string > > const & atom_pairs )
{
	core::Size parent_rsd;

	ResidueOP new_rsd( new_rsd_in.clone() ); // need a copy since we may be modifying this

	if ( !Symmetry_Info()->bb_is_independent( seqpos ) ) {
		TR.Debug << "MirrorSymmetricConformation::replace_residue(3) setting a dependent TORSION!, try to set its parent" << std::endl;
		parent_rsd = Symmetry_Info()->bb_follows( seqpos ) ;

		// unlike non-mirror symmetries, we need to call this function to ensures residues are properly mirrored D<->L
		for ( int i=1; i<=(int)new_rsd->natoms(); ++i ) {
			new_rsd->set_xyz(i , apply_transformation( new_rsd_in.xyz(i), seqpos, parent_rsd ) );
		}
		if ( res_is_mirrored_[seqpos] ) {
			flip_chirality( new_rsd );
		}
	} else {
		parent_rsd = seqpos;
	}

	// set this residue using base-class method
	Conformation::replace_residue( parent_rsd, *new_rsd, atom_pairs );

	// now the copies
	for ( SymmetryInfo::Clones::const_iterator pos=Symmetry_Info()->bb_clones( parent_rsd ).begin(),
			epos=Symmetry_Info()->bb_clones( parent_rsd ).end(); pos != epos; ++pos ) {
		ResidueOP new_new_rsd( new_rsd->clone() );
		for ( int i=1; i<=(int)new_new_rsd->natoms(); ++i ) {
			new_new_rsd->set_xyz(i , apply_transformation( new_rsd->xyz(i), parent_rsd, *pos ) );
		}
		if ( res_is_mirrored_[*pos] ) {
			flip_chirality( new_new_rsd );
		}
		Conformation::replace_residue( *pos, *new_new_rsd, false );
	}

	//fpd may have implicitly changed a jump
	SymmetricConformation::clear_Tsymm( );
}

void
MirrorSymmetricConformation::fold_tree( FoldTree const & fold_tree_in ) {
	SymmetricConformation::fold_tree( fold_tree_in );
	update_njumps_nres();
}

// recalculate the Tsymm_ transforms using the current pose
void
MirrorSymmetricConformation::recalculate_transforms( ) {
	using namespace numeric;

	SymmetricConformation::recalculate_transforms( );

	// now flip based on mirroring
	core::Size nres_per_sub = Symmetry_Info()->get_nres_subunit();
	core::Size nsubunits = Symmetry_Info()->subunits();
	Size ncomps = Symmetry_Info()->get_num_components();

	//fpd this should be multicomp friendly ... though I'm not sure if this is necessary?
	for ( Size icomp=1; icomp <= ncomps; ++icomp ) {
		char comptag = (ncomps==1)? 'A' : Symmetry_Info()->get_component(icomp);
		for ( Size isub=1; isub <= nsubunits; ++isub ) {
			Size substart=(isub-1)*nres_per_sub+1, substop=isub*nres_per_sub;
			Size ires=substart;

			// if this is too inefficient we can store this in symminfo
			if ( ncomps!=1 ) {
				for ( ires=substart; ires<=substop && Symmetry_Info()->get_component(ires) != comptag; ++ires ) ;
				assert(Symmetry_Info()->get_component(ires) == comptag);
			}

			if ( res_is_mirrored_[ires] ) {
				SymmetricConformation::invert_Tsymm( comptag, isub );
			}
		}
	}
}


//fpd
void
MirrorSymmetricConformation::append_residue_by_jump(
	conformation::Residue const & new_rsd,
	Size const anchor_pos,
	std::string const& anchor_atom, // could be zero
	std::string const& root_atom, // ditto
	bool const start_new_chain // default false
)
{
	SymmetricConformation::append_residue_by_jump(new_rsd, anchor_pos, anchor_atom, root_atom, start_new_chain);

	update_njumps_nres();
	update_residue_identities();
}

void
MirrorSymmetricConformation::insert_conformation_by_jump(
	Conformation const & new_conf,
	Size const insert_seqpos,
	Size const,
	Size const anchor_pos,
	Size const anchor_jump_number,
	std::string const & anchor_atom,
	std::string const & root_atom
) {
	if ( anchor_jump_number != 0 ) {
		TR.Warning << "MirrorSymmetricConformation::insert_conformation_by_jump ignores anchor_jump_number" << std::endl;
	}

	core::Size nmonomer_jumps = Symmetry_Info()->get_njumps_subunit();
	core::Size nres_monomer = Symmetry_Info()->get_nres_subunit();
	core::Size nsubunits = Symmetry_Info()->subunits();
	core::Size asymm_insert = ((insert_seqpos-2)%nres_monomer) + 2;  //?
	core::Size asymm_anchor = ((anchor_pos-1)%nres_monomer) + 1;

	// insert at the end of each subunit
	// transform to the coordinate frame of the scoring subunit
	// go from last->first so we don't have to worry about offsets
	for ( int i=nsubunits; i>=1; --i ) {
		core::Size insert_i = (i-1)*nres_monomer+asymm_insert;
		core::Size anchor_i = (i-1)*nres_monomer+asymm_anchor;

		ConformationOP new_new_conf = new_conf.clone();
		if ( !Symmetry_Info()->bb_is_independent( anchor_i ) ) {
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
	Symmetry_Info()->resize_asu( nres_monomer + new_conf.size() );
	Symmetry_Info()->update_nmonomer_jumps( nmonomer_jumps + new_conf.fold_tree().num_jump() + 1 );

	// update ft
	FoldTree f_in = core::conformation::symmetry::get_asymm_unit_fold_tree( *this );
	core::conformation::symmetry::symmetrize_fold_tree( *this, f_in );
	fold_tree( f_in );
}

//
// @brief Detect existing disulfides from the protein structure.
// @details For full atom confomations, looks at SG-SG distance. If the SG-SG
//  are about 2.02 A apart, calls it a disulfide bond. For centroid and other
//  conformations, the less accurate CB-CB distance is used instead. In this
//  case a CB-CB distance of 3.72 A is optimal.
void
MirrorSymmetricConformation::detect_disulfides( utility::vector1< Size > const & disulf_one, utility::vector1< Size > const & disulf_two )
{
	SymmetricConformation::detect_disulfides(disulf_one, disulf_two);

	// TODO: Vikram
	// coordinates will be correct at this point
	// however, you need to handle the correct variant types for inter-subunit disulfides
	// this can be done after the fact (calling the base class implementation) or on-the-fly
}

/// @brief Helper function to flip the chirality of a residue type.
/// @details Assumes that the coordinates are already correct for the flipped type (i.e. this function does
/// not alter atom positions -- only the ResidueType identity).
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
MirrorSymmetricConformation::flip_chirality( ResidueOP & new_rsd ) {
	if ( new_rsd->type().is_achiral_backbone() ) return; //Do nothing if this residue is achiral.

	// Now copy the replacement residue to the input residue:
	new_rsd = new_rsd->clone_flipping_chirality();
}

/// @brief helper functions called when #jumps or #residues changes, updating internal data
///   just re-traverses the fold tree, if this needs to be made faster (why would it though?)
///   this could use logic similar to SymmetryInfo::update_nmonomer_jumps()
void
MirrorSymmetricConformation::update_njumps_nres(  ) {
	jump_is_mirrored_.clear(); jump_is_mirrored_.resize( Conformation::fold_tree().num_jump() );
	res_is_mirrored_.clear(); res_is_mirrored_.resize( Symmetry_Info()->num_total_residues(), false);
	utility::vector1<bool> subunit_is_mirrored( Symmetry_Info()->subunits(), false );

	calculate_inverting_virtuals( Conformation::fold_tree(), *this, *Symmetry_Info(), subunit_is_mirrored, jump_is_mirrored_ );

	Size nsubs = Symmetry_Info()->subunits(), nres_monomer=Symmetry_Info()->num_independent_residues();
	TR << "MIRRORED subunits: ";
	for ( Size i=1; i<=nsubs; ++i ) {
		if ( subunit_is_mirrored[i] ) TR << i << " ";
	}
	TR << std::endl;
	TR << "MIRRORED jumps: ";
	for ( Size i=1; i<=jump_is_mirrored_.size(); ++i ) {
		if ( jump_is_mirrored_[i].first ) TR << i << "u ";
		if ( jump_is_mirrored_[i].second ) TR << i << "d ";
	}
	TR << std::endl;

	// from subunit map, make residue map
	for ( Size i=1; i<=nsubs; ++i ) {
		for ( Size j=1; j<=nres_monomer; ++j ) {
			res_is_mirrored_[(i-1)*nres_monomer+j] = subunit_is_mirrored[i];
		}
	}

	synch_mirror_jumps_with_atomtree();
}

void
MirrorSymmetricConformation::synch_mirror_jumps_with_atomtree(  ) {
	core::Size njump = Conformation::fold_tree().num_jump();
	for ( core::Size i=1; i<=njump; ++i ) {
		kinematics::Jump jump_copy = Conformation::jump(i);
		jump_copy.set_invert(
			jump_is_mirrored_[i].first,
			jump_is_mirrored_[i].second
		);
		Conformation::set_jump(i, jump_copy);
	}

	// to be safe
	SymmetricConformation::clear_Tsymm();
}



/// @brief Updates residue identities in symmetric subunits, ensuring that they are mirrored relative to the ASU in mirrored subunits
/// and identical to the ASU in non-mirrored subunits.
/// @details Assumes that the residue identities and variants (aside from D/L variants) already match.  That is, if I have ASN at position
/// 5 in my asymmetric unit, I either have ASN or DASN at the equivalent position in each symmetry copy.  Safe to call repeatedly.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
MirrorSymmetricConformation::update_residue_identities( ) {
	TR << "Inverting restypes" << std::endl;

	for ( Size i=1; i<=Symmetry_Info()->num_total_residues_without_pseudo(); ++i ) {
		if ( Symmetry_Info()->bb_is_independent(i) ) {
			core::conformation::Residue new_res = residue(i);
			this->replace_residue( i, new_res, false );
		}
	}
}

/// @brief Is this residue mirrored relative to the asymmetric unit?
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
MirrorSymmetricConformation::res_is_mirrored (
	core::Size const seqpos
) const {
	runtime_assert_string_msg(
		seqpos > 0 && seqpos <= res_is_mirrored_.size(),
		"Error in core::conformation::symmetry::MirrorSymmetricConformation::res_is_mirrored(): The sequence position is out of range."
	);
	return res_is_mirrored_[seqpos];
}

}
} // conformation
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::symmetry::MirrorSymmetricConformation::save( Archive & arc ) const {
	arc( cereal::base_class< core::conformation::Conformation >( this ) );
	arc( CEREAL_NVP( jump_is_mirrored_ ) ); // utility::vector1< bool >
	arc( CEREAL_NVP( res_is_mirrored_ ) ); // utility::vector1< bool >
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::symmetry::MirrorSymmetricConformation::load( Archive & arc ) {
	arc( cereal::base_class< core::conformation::Conformation >( this ) );
	arc( jump_is_mirrored_ );// utility::vector1< bool >
	arc( res_is_mirrored_ );// utility::vector1< bool >
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::symmetry::MirrorSymmetricConformation );
CEREAL_REGISTER_TYPE( core::conformation::symmetry::MirrorSymmetricConformation )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_symmetry_MirrorSymmetricConformation )
#endif // SERIALIZATION
