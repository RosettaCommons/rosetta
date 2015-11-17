// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/pose/carbohydrates/util.cc
/// @brief   Utility function definitions for carbohydrate-containing poses.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit headers
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/Pose.hh>

// Package headers
#include <core/io/carbohydrates/pose_io.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>

// Project headers
#include <core/id/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/types.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/angle.functions.hh>

// External headers
#include <boost/lexical_cast.hpp>

// C++ header
#include <list>


// Construct tracer.
static THREAD_LOCAL basic::Tracer TR( "core.pose.carbohydrates.util" );


namespace core {
namespace pose {
namespace carbohydrates {

using namespace std;
using namespace core;

// Helper Functions ///////////////////////////////////////////////////////////
// Use a saccharide residue's connections to find the residue from which it follows or branches.
/// @return  The sequence position of the residue before this one (n-1) or the residue in the parent chain from which
/// the branch occurs or zero if N/A, i.e., if this is the lower terminus.
core::uint
find_seqpos_of_saccharides_parent_residue( conformation::Residue const & residue ) {
	debug_assert( residue.is_carbohydrate() );

	if ( ! residue.is_lower_terminus() ) {
		uint const id_of_connection_to_parent(
			residue.type().residue_connection_id_for_atom( residue.carbohydrate_info()->anomeric_carbon_index() ) );
		return residue.residue_connection_partner( id_of_connection_to_parent );
	} else /* residue is lower terminus */ {
		if ( TR.Debug.visible() ) {
			TR.Debug << "This residue is a lower terminus! Returning 0." << endl;
		}
		return 0;
	}
}


// Return pointers to the two residues of the glycosidic bond.
/// @return  Pointers to the residue at <sequence_position> and its parent or else the same pointer twice if undefined.
std::pair< conformation::ResidueCOP, conformation::ResidueCOP >
get_glycosidic_bond_residues( Pose const & pose, uint const sequence_position )
{
	using namespace conformation;

	// Get the 1st residue of interest.
	ResidueCOP res_n( pose.residue( sequence_position ).get_self_ptr() );

	if ( res_n->is_lower_terminus() ) {
		if ( TR.Info.visible() ) {
			TR.Info << "Glycosidic torsions are undefined for the first polysaccharide residue of a chain unless part "
				"of a branch." << endl;
		}
		return make_pair( res_n, res_n );
	}

	// Get the 2nd residue of interest.
	// (res_n_minus_1 is a misnomer for the lower termini of branches.)
	ResidueCOP res_n_minus_1 = pose.residue( find_seqpos_of_saccharides_parent_residue( *res_n ) ).get_self_ptr();

	return make_pair( res_n, res_n_minus_1 );
}


// Return the AtomIDs of the four phi torsion reference atoms.
/// @details For aldopyranoses, phi is defined as O5(n)-C1(n)-OX(n-1)-CX(n-1),
/// where X is the position of the glycosidic linkage.\n
/// For aldofuranoses, phi is defined as O4(n)-C1(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketopyranoses, phi is defined as O6(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketofuranoses, phi is defined as O5(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// Et cetera...\n
utility::vector1< id::AtomID >
get_reference_atoms_for_phi( Pose const & pose, uint const sequence_position )
{
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1< AtomID > ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair< ResidueCOP, ResidueCOP > const residues( get_glycosidic_bond_residues( pose, sequence_position ) );

	if ( residues.first->seqpos() == residues.second->seqpos() ) {  // This occurs when there is no parent residue.
		return ids;
	}
	ids.resize( 4 );  // A torsion has 4 reference atoms.

	// Set the atom names of the four reference atoms.
	// Reference 1 is O(cyclic) for cyclic saccharides, C? for linear saccharides.
	// Because the cyclic O is not connected by the atom tree,
	// we actually will return the AtomID for the virtual atom that superimposes with the real atom.
	AtomID ref1;
	if ( residues.first->carbohydrate_info()->is_cyclic() ) {
		ref1 = AtomID( residues.first->carbohydrate_info()->virtual_cyclic_oxygen_index(), residues.first->seqpos() );
	} else /* is linear */ {
		;  // TODO: Figure out how linear polysaccharides are handled by IUPAC.
	}
	ids[ 1 ] = ref1;

	// Reference 2 is always the anomeric carbon.
	AtomID const ref2( residues.first->carbohydrate_info()->anomeric_carbon_index(), residues.first->seqpos() );
	ids[ 2 ] = ref2;

	// Reference 3 is OX(n-1) for polysaccharides.
	AtomID const ref3( residues.second->connect_atom( *residues.first ), residues.second->seqpos() );
	ids[ 3 ] = ref3;

	// Reference 4 is CX(n-1) for polysaccharides.
	AtomID const ref4( residues.second->type().atom_base( ref3.atomno() ), residues.second->seqpos() );
	//AtomID ref4( residues.second->first_adjacent_heavy_atom( ref3.atomno() ), residues.second->seqpos() );
	ids[ 4 ] = ref4;

	if ( TR.Debug.visible() ) {
		TR.Debug << "Reference atoms for phi: " << ref1 << ", " << ref2 << ", " << ref3 << ", " << ref4 << endl;
	}

	return ids;
}


// Return the AtomIDs of the four psi torsion reference atoms.
/// @details For saccharides, psi is defined as: C(anomeric)(n)-OX(n-1)-CX(n-1)-CX-1(n-1),\n
/// where X is the position of the glycosidic linkage.
utility::vector1< id::AtomID >
get_reference_atoms_for_psi( Pose const & pose, uint const sequence_position )
{
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1< AtomID > ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair< ResidueCOP, ResidueCOP > const residues( get_glycosidic_bond_residues( pose, sequence_position ) );

	if ( residues.first->seqpos() == residues.second->seqpos() ) {  // This occurs when there is no parent residue.
		return ids;
	}

	ids.resize( 4 );  // A torsion has 4 reference atoms.

	// Set the atom names of the four reference atoms.
	// Reference 1 is always the anomeric carbon.
	AtomID const ref1( residues.first->carbohydrate_info()->anomeric_carbon_index(), residues.first->seqpos() );
	ids[ 1 ] = ref1;

	// Reference 2 is OX(n-1) for polysaccharides.
	AtomID const ref2( residues.second->connect_atom( *residues.first ), residues.second->seqpos() );
	ids[ 2 ] = ref2;

	// Reference 3 is CX(n-1) for polysaccharides.
	AtomID const ref3( residues.second->type().atom_base( ref2.atomno() ), residues.second->seqpos() );
	//AtomID ref3( residues.second->first_adjacent_heavy_atom( ref2.atomno() ), residues.second->seqpos() );
	ids[ 3 ] = ref3;

	// Reference 4 is CX-1(n-1) for polysaccharides.
	AtomID const ref4( residues.second->type().atom_base( ref3.atomno() ), residues.second->seqpos() );
	ids[ 4 ] = ref4;

	if ( TR.Debug.visible() ) {
		TR.Debug << "Reference atoms for psi: " << ref1 << ", " << ref2 << ", " << ref3 << ", " << ref4 << endl;
	}

	return ids;
}


// Return the AtomIDs of the four omega torsion reference atoms.
/// For carbohydrates glycosylated at an exocyclic position,
/// omega of residue n is defined as OX(n-1)-CX(n-1)-CX-1(n-1)-CX-2(n-1),
/// where X is the position of the glycosidic linkage.
utility::vector1< id::AtomID >
get_reference_atoms_for_1st_omega( Pose const & pose, uint const sequence_position )
{
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1< AtomID > ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair< ResidueCOP, ResidueCOP > const residues( get_glycosidic_bond_residues( pose, sequence_position ) );

	if ( residues.first->seqpos() == residues.second->seqpos() ) {  // This occurs when there is no parent residue.
		return ids;
	}
	if ( residues.second->is_carbohydrate() && ( ! residues.second->carbohydrate_info()->has_exocyclic_linkage() ) ) {
		TR.Warning << "Omega is undefined for this residue, because the glycosidic linkage is not exocyclic." << endl;
		return ids;
	}

	ids.resize( 4 );  // A torsion has 4 reference atoms.

	// Set the atom names of the four reference atoms.
	// Reference 1 is OX(n-1) for polysaccharides.
	AtomID const ref1( residues.second->connect_atom( *residues.first ), residues.second->seqpos() );
	ids[ 1 ] = ref1;

	// Reference 2 is CX(n-1) for polysaccharides.
	AtomID const ref2( residues.second->type().atom_base( ref1.atomno() ), residues.second->seqpos() );
	//AtomID ref2( residues.second->first_adjacent_heavy_atom( ref1.atomno() ), residues.second->seqpos() );
	ids[ 2 ] = ref2;

	// Reference 3 is CX-1(n-1) for polysaccharides.
	AtomID const ref3( residues.second->type().atom_base( ref2.atomno() ), residues.second->seqpos() );
	ids[ 3 ] = ref3;

	// Reference 4 is CX-2(n-1) for polysaccharides.
	AtomID const ref4( residues.second->type().atom_base( ref3.atomno() ), residues.second->seqpos() );
	ids[ 4 ] = ref4;

	if ( TR.Debug.visible() ) {
		TR.Debug << "Reference atoms for omega: " << ref1 << ", " << ref2 << ", " << ref3 << ", " << ref4 << endl;
	}

	return ids;
}


// Return the AtomIDs of the four reference atoms for the requested torsion.
utility::vector1< id::AtomID >
get_reference_atoms( uint const torsion_id, Pose const & pose, uint const sequence_position )
{
	using namespace id;
	using namespace utility;

	vector1< AtomID > ref_atoms;
	switch ( torsion_id ) {
	case phi_torsion :
		ref_atoms = get_reference_atoms_for_phi( pose, sequence_position );
		break;
	case psi_torsion :
		ref_atoms = get_reference_atoms_for_psi( pose, sequence_position );
		break;
	case omega_torsion :
		ref_atoms = get_reference_atoms_for_1st_omega( pose, sequence_position );
		break;
	default :
		utility_exit_with_message( "An invalid torsion angle was requested." );
	}
	return ref_atoms;
}


// Virtual Atom Alignment /////////////////////////////////////////////////////
// Set coordinates of virtual atoms (used as angle reference points) within a saccharide residue of the given
// conformation.
/// @details  This method aligns virtual atom VOX, where X is the position of the cyclic oxygen, OY and HOY, where Y is
/// the position of the anomeric carbon, (provided the residue is not the reducing end, where OY and HOY would be real
/// atoms), and HOZ, where Z is the mainchain glycosidic bond location.  OY and HOY are aligned with the last two "main
/// chain" atoms of the parent residue.  This ensures that torsion angles with duplicate names, e.g., chi1 and phi for
/// internal linked aldoses, will always return the same values.  The same concept applies for HOZ, which aligns with
/// the anomeric carbon of the downstream residue.  This method should be called after any coordinate change for a sac-
/// charide residue and after loading a saccharide residue from a file or sequence for the first time.
/// @note     Do I need to worry about aligning the virtual atoms left over from modified sugar patches?  Can such
/// virtual atoms be deleted?
void
align_virtual_atoms_in_carbohydrate_residue( conformation::Conformation & conf, uint const sequence_position ) {
	using namespace id;
	using namespace conformation;

	TR.Debug << " Aligning virtual atoms on residue " << sequence_position << "..." << endl;

	ResidueCOP res( conf.residue( sequence_position ).get_self_ptr() );

	// Find and align VOX, if applicable.
	if ( res->carbohydrate_info()->is_cyclic() ) {
		if ( TR.Debug.visible() ) {
			TR.Debug << "  Aligning VOX..." << endl;
		}
		uint const x( res->carbohydrate_info()->cyclic_oxygen() );
		uint const OX( res->atom_index( res->carbohydrate_info()->cyclic_oxygen_name() ) );
		uint const VOX( res->atom_index( "VO" + string( 1, x + '0' ) ) );

		conf.set_xyz( AtomID( VOX, sequence_position ), conf.xyz( AtomID( OX, sequence_position ) ) );
		if ( TR.Debug.visible() ) {
			TR.Debug << "  VOX aligned." << endl;
		}
	}

	// Find and align OY and HOY, if applicable.
	if ( ! res->is_lower_terminus() ) {
		if ( TR.Debug.visible() ) {
			TR.Debug << "  Aligning OY and HOY..." << endl;
		}
		uint const y( res->carbohydrate_info()->anomeric_carbon() );
		uint const OY( res->atom_index( "O" + string( 1, y + '0') ) );
		uint const HOY( res->atom_index( "HO" + string( 1, y + '0' ) ) );

		uint const parent_res_seqpos( find_seqpos_of_saccharides_parent_residue( *res ) );
		ResidueCOP parent_res( conf.residue( parent_res_seqpos ).get_self_ptr() );
		uint const OY_ref( parent_res->connect_atom( *res ) );
		uint const HOY_ref( parent_res->first_adjacent_heavy_atom( OY_ref ) );

		conf.set_xyz( AtomID( HOY, sequence_position ), conf.xyz( AtomID( HOY_ref, parent_res_seqpos ) ) );
		if ( TR.Debug.visible() ) {
			TR.Debug << "  HOY aligned with atom " << parent_res->atom_name( HOY_ref ) <<
				" of residue " << parent_res_seqpos << endl;

			TR.Debug << "   Updating torsions..." << endl;
		}
		ResidueCOP dummy( conf.residue( sequence_position ).get_self_ptr() );  // to trigger private method commented below
		//conf.update_residue_torsions( res->seqpos(), false );
		if ( TR.Debug.visible() ) {
			TR.Debug << "   Torsions updated." << endl;
		}

		conf.set_xyz( AtomID( OY, sequence_position ), conf.xyz( AtomID( OY_ref, parent_res_seqpos ) ) );
		if ( TR.Debug.visible() ) {
			TR.Debug << "  OY aligned with atom " << parent_res->atom_name( OY_ref ) <<
				" of residue " << parent_res_seqpos << endl;
		}
	}

	// Find and align HOZ(s), if applicable.
	if ( ! res->is_upper_terminus() ) {
		if ( TR.Debug.visible() ) {
			TR.Debug << "  Aligning HOZ..." << endl;
		}
		uint const z( res->carbohydrate_info()->mainchain_glycosidic_bond_acceptor() );
		uint const HOZ( res->atom_index( "HO" + string( 1, z + '0' ) ) );

		uint const downstream_res_seqpos( sequence_position + 1 );
		ResidueCOP downstream_res( conf.residue( downstream_res_seqpos ).get_self_ptr() );
		uint const HOZ_ref( downstream_res->atom_index( downstream_res->carbohydrate_info()->anomeric_carbon_name() ) );

		conf.set_xyz( AtomID( HOZ, sequence_position ), conf.xyz( AtomID( HOZ_ref, downstream_res_seqpos ) ) );
		if ( TR.Debug.visible() ) {
			TR.Debug << "  HOZ aligned." << endl;
		}
	}
	Size const n_branches( res->carbohydrate_info()->n_branches() );
	for ( uint branch_num( 1 ); branch_num <= n_branches; ++branch_num ) {
		uint const z( res->carbohydrate_info()->branch_point( branch_num ) );
		uint const OZ( res->atom_index( "O" + string( 1, z + '0' ) ) );
		uint const HOZ( res->atom_index( "HO" + string( 1, z + '0' ) ) );

		uint const branch_connection_id( res->type().residue_connection_id_for_atom( OZ ) );
		uint const branch_res_seqpos( res->residue_connection_partner( branch_connection_id ) );
		ResidueCOP branch_res( conf.residue( branch_res_seqpos ).get_self_ptr() );
		uint const HOZ_ref( branch_res->atom_index( branch_res->carbohydrate_info()->anomeric_carbon_name() ) );

		conf.set_xyz( AtomID( HOZ, sequence_position ), conf.xyz( AtomID( HOZ_ref, branch_res_seqpos ) ) );
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << " All virtual atoms aligned." << endl;
	}
}


// Set coordinates of virtual atoms (used as angle reference points) within a saccharide residue of the given pose.
void
align_virtual_atoms_in_carbohydrate_residue( Pose & pose, uint const sequence_position ) {
	align_virtual_atoms_in_carbohydrate_residue( pose.conformation(), sequence_position );
}


// TorsionID Queries //////////////////////////////////////////////////////////
// Is this is the phi torsion angle of a glycosidic linkage?
/// @details  Carbohydrate linkages are defined as the torsion angles leading back to the previous residue.  Much
/// of Rosetta code relies on TorsionIDs and assumes TorsionID( n, BB, 1 ) is phi.  For a sugar, phi (of the next
/// residue) is the last torsion, and the number of main-chain torsions varies per saccharide residue.
bool
is_glycosidic_phi_torsion( Pose const & pose, id::TorsionID const & torsion_id )
{
	using namespace id;

	if ( torsion_id.type() == BB || torsion_id.type() == BRANCH ) {  // Phi for a branched saccharide has a BRANCH type.
		conformation::Residue const & residue( pose.residue( torsion_id.rsd() ) );
		uint next_rsd_num( 0 );  // We will need to see if the "next" residue is a saccharide.

		switch( torsion_id.type() ) {
		case  BB :
			if ( ! residue.is_upper_terminus() ) {
				// If this is a main-chain torsion, we need the next residue on the main chain.
				next_rsd_num = torsion_id.rsd() + 1;
				if ( pose.residue( next_rsd_num ).is_carbohydrate() ) {
					return ( torsion_id.torsion() == residue.n_mainchain_atoms() );  // The last BB is phi.
				}
			}
			break;
		case BRANCH :
			{
			Size const n_mainchain_connections( residue.n_polymeric_residue_connections() );
			next_rsd_num = residue.residue_connection_partner( n_mainchain_connections + torsion_id.torsion() );
			if ( pose.residue( next_rsd_num ).is_carbohydrate() ) {
				return true;  // If it's a branch to a sugar, it must be the phi from the branching residue.
			}
		}
			break;
		default :
			break;
		}

	}
	return false;
}

// Is this is the psi torsion angle of a glycosidic linkage?
/// @details  Carbohydrates linkages are defined as the torsion angles leading back to the previous residue.  Much
/// of Rosetta code relies on TorsionIDs and assumes TorsionID( n, BB, 2 ) is psi.  For a sugar, psi (of the next
/// residue) is the penultimate torsion, and the number of main-chain torsions varies per saccharide residue.
bool
is_glycosidic_psi_torsion( Pose const & pose, id::TorsionID const & torsion_id )
{
	using namespace id;

	if ( torsion_id.type() == BB || torsion_id.type() == CHI ) {  // Psi for a branched saccharide has a CHI type.
		conformation::Residue const & residue( pose.residue( torsion_id.rsd() ) );
		uint next_rsd_num( 0 );  // We will need to see if the "next" residue is a saccharide.

		switch( torsion_id.type() ) {
		case BB :
			if ( ! residue.is_upper_terminus() ) {
				// If this is a main-chain torsion, we need the next residue on the main chain.
				next_rsd_num = torsion_id.rsd() + 1;
				if ( pose.residue( next_rsd_num ).is_carbohydrate() ) {
					return ( torsion_id.torsion() == residue.n_mainchain_atoms() - 1 );
				}
			}
			break;
		case CHI :
			{
			Size const n_mainchain_connections( residue.n_polymeric_residue_connections() );
			// A psi angle will always have the third atom of its definition be a connect atom.
			uint const third_atom( residue.chi_atoms( torsion_id.torsion() )[ 3 ] );
			Size const n_branches( residue.n_non_polymeric_residue_connections() );
			for ( uint branch_num( 1 ); branch_num <= n_branches; ++branch_num ) {
				next_rsd_num = residue.residue_connection_partner( n_mainchain_connections + branch_num );
				conformation::Residue const & next_rsd( pose.residue( next_rsd_num ) );
				if ( next_rsd.is_carbohydrate() ) {
					return ( residue.connect_atom( next_rsd ) == third_atom );
				}
			}
		}
			break;
		default :
			break;
		}
	}
	return false;
}

// Is this is an omega torsion angle of a glycosidic linkage?
/// @details  Carbohydrates linkages are defined as the torsion angles leading back to the previous residue.  Much
/// of Rosetta code relies on TorsionIDs and assumes TorsionID( n, BB, 3 ) is omega.  For a sugar, omega (of the next
/// residue) is the 3rd-to-last torsion, and the number of main-chain torsions varies per saccharide residue.
/// @remarks  Currently, this only works properly if there is a single omega and not yet for exocyclic branches either.
bool
is_glycosidic_omega_torsion( Pose const & pose, id::TorsionID const & torsion_id )
{
	using namespace id;

	if ( torsion_id.type() == BB || torsion_id.type() == CHI ) {  // Omega for a branched saccharide has a CHI type.
		conformation::Residue const & residue( pose.residue( torsion_id.rsd() ) );
		uint next_rsd_num( 0 );  // We will need to see if the "next" residue is a saccharide.

		switch( torsion_id.type() ) {
		case BB :
			if ( ! residue.is_upper_terminus() ) {
				// If this is a main-chain torsion, we need the next residue on the main chain.
				next_rsd_num = torsion_id.rsd() + 1;
				if ( pose.residue( next_rsd_num ).is_carbohydrate() ) {
					chemical::carbohydrates::CarbohydrateInfoCOP info( residue.carbohydrate_info() );
					if ( info->has_exocyclic_linkage() ) {
						return ( torsion_id.torsion() == residue.n_mainchain_atoms() - 2 );
					}
				}
			}
			break;
		case CHI :
			{
			Size const n_mainchain_connections( residue.n_polymeric_residue_connections() );
			// An omega angle will always have the fourth atom of its definition be a connect atom.
			uint const fourth_atom( residue.chi_atoms( torsion_id.torsion() )[ 4 ] );
			Size const n_branches( residue.n_non_polymeric_residue_connections() );
			for ( uint branch_num( 1 ); branch_num <= n_branches; ++branch_num ) {
				next_rsd_num = residue.residue_connection_partner( n_mainchain_connections + branch_num );
				conformation::Residue const & next_rsd( pose.residue( next_rsd_num ) );
				if ( next_rsd.is_carbohydrate() ) {
					return ( residue.connect_atom( next_rsd ) == fourth_atom );
				}
			}
		}
			break;
		default :
			break;
		}
	}
	return false;
}


// Torsion Access /////////////////////////////////////////////////////////////
// Getters ////////////////////////////////////////////////////////////////////
// Return the requested torsion angle between a saccharide residue of the given pose and the previous residue.
/// @details This method is used in place of Residue::mainchain_torsion() since the main-chain torsions of saccharide
/// residues only make sense in the context of two residues.  Moreover, the reference atoms as defined by the IUPAC are
/// different from the ones that Rosetta uses by default for mainchain torsions for sugars.
/// @param   <torsion_id> is an integer representing the specific torsion angle requested, as defined in
/// core/id/types.hh:\n
///     phi_torsion = 1\n
///     psi_torsion = 2\n
///     omega_torsion = 3
/// @note    I would rather the torsion id were an enum, but as it was already defined, I'm leaving it as a constant
/// for now.
core::Angle
get_glycosidic_torsion( uint const torsion_id, Pose const & pose, uint const sequence_position )
{
	using namespace id;
	using namespace numeric;
	using namespace utility;

	vector1< AtomID > const ref_atoms = get_reference_atoms( torsion_id, pose, sequence_position );

	if ( ref_atoms.size() == 0 ) {
		// This occurs when there is no parent residue or when the glycosidic bond is not exocyclic (omega only).
		TR.Warning << "Returning zero." << endl;
		return 0.0;
	}

	Angle const angle_in_radians(
		pose.conformation().torsion_angle( ref_atoms[ 1 ], ref_atoms[ 2 ], ref_atoms[ 3 ], ref_atoms[ 4 ]) );
	return principal_angle_degrees( conversions::degrees( angle_in_radians ) );
}


// Setters ////////////////////////////////////////////////////////////////////
// Set the requested torsion angle between a saccharide residue of the given pose and the previous residue.
/// @details This method is used in place of Conformation::set_torsion() since the reference atoms as defined by the
/// IUPAC are different from the ones that Rosetta uses by default for main-chain torsions for sugars.
/// @param   <torsion_id> is an integer representing the specific torsion angle requested, as defined in
/// core/id/types.hh:\n
///     phi_torsion = 1\n
///     psi_torsion = 2\n
///     omega_torsion = 3\n
/// <setting> is in degrees.
/// @note    I would rather the torsion id were an enum, but as it was already defined, I'm leaving it as a constant
/// for now.
void
set_glycosidic_torsion( uint const torsion_id, Pose & pose, uint const sequence_position, core::Angle const setting )
{
	using namespace id;
	using namespace numeric;
	using namespace utility;

	vector1< AtomID > const ref_atoms = get_reference_atoms( torsion_id, pose, sequence_position );

	if ( ref_atoms.size() == 0 ) {
		// This occurs when there is no parent residue or when the glycosidic bond is not exocyclic (omega only).
		return;
	}

	Angle const setting_in_radians( conversions::radians( setting ) );
	pose.conformation().set_torsion_angle(
		ref_atoms[ 1 ], ref_atoms[ 2 ], ref_atoms[ 3 ], ref_atoms[ 4 ], setting_in_radians );
	align_virtual_atoms_in_carbohydrate_residue( pose, sequence_position );
}


// Glycosylation //////////////////////////////////////////////////////////////
// Glycosylate the Pose at the given sequence position and atom using an IUPAC sequence.
/// @details   Format for <iupac_sequence>:\n
/// Prefixes apply to the residue to which they are attached, below indicated by residue n.\n
/// Residues are listed from N to 1, where N is the total number of residues in the saccharide.\n
/// The sequence is parsed by reading to the next hyphen, so hyphens are crucial.\n
/// Linkage indication: "(a->x)-" specifies the linkage of residue n, where a is the anomeric carbon number of residue
/// (n+1) and x is the oxygen number of residue n.  The first residue listed in the annotated sequence (residue N)
/// need not have the linkage prefix.  A ->4) ResidueType will automatically be assigned by default if not specified.\n
/// Anomer indication: The strings "alpha-" or "beta-" are supplied next, which determines the stereochemistry of the
/// anomeric carbon of the residue to which it is prefixed.  An alpha ResidueType will automatically be assigned by
/// default.\n
/// Stereochemical indication: "L-" or "D-" specifies whether residue n is an L- or D-sugar.  The default is "D-".\n
/// 3-Letter code: A three letter code (in sentence case) MUST be supplied next.  This specifies the "base sugar name",
/// e.g., Glc is for glucose.  (A list of all recognized 3-letter codes for sugars can be found in
/// database/chemical/carbohydrates/codes_to_roots.map.)\n
/// 1-Letter suffix: If no suffix follows, residue n will be linear.  If a letter is present, it indicates the ring
/// size, where "f" is furanose, "p" is puranose, and "s" is septanose.\n
/// Branches are indicated using nested brackets and are best explained by example:\n
/// beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpNAc- is:\n
/// beta-D-Galp-(1->4)-D-GlcpNAc-\n
///                       |\n
///     alpha-L-Fucp-(1->3)
/// @remarks   The order of residues in the created pose is in the opposite direction as the annotated sequence of
/// saccharide residues, as sugars are named with the 1st residue as the "main chain", with all other residues named as
/// substituents and written as prefixes.  In other words, sugars are usually drawn and named with residue 1 to the
/// right.
/// At present time, param files only exist for a limited number of sugars! ~ Labonte
/// Also, I've not written this to handle creation of glycolipids yet.
void
glycosylate_pose(
	Pose & pose,
	uint const sequence_position,
	std::string const & atom_name,
	std::string const & iupac_sequence )
{
	using namespace utility;
	using namespace chemical;
	using namespace conformation;

	conformation::Residue const & residue( pose.residue( sequence_position ) );
	ResidueTypeSet const & residue_set( residue.type().residue_type_set() );

	// Get list of carbohydrate ResidueTypes from which to construct the Pose.
	ResidueTypeCOPs residue_types( residue_types_from_saccharide_sequence( iupac_sequence, residue_set ) );

	Size const n_types( residue_types.size() );
	if ( ! n_types ) {
		TR.Warning << "No saccharide residues in sequence to append to pose!" << endl;
		return;
	}

	// We need to reorder the ResidueTypes, such that each chain is complete before a new branch is added.
	reorder_saccharide_residue_types( residue_types );

	// Prepare the glycosylation site for branching.
	// TODO: Add code to attach to lipids.
	if ( residue.is_carbohydrate() ) {
		;  // TODO: Add complicated check here to figure out which VariantType to assign.
	} else {  // It's a typical branch point, in that it has a single type of branch point variant.
		add_variant_type_to_pose_residue( pose, SC_BRANCH_POINT, sequence_position );
	}


	// Now we can extend the Pose.
	// Keep track of branch points as we go.
	list< pair< uint, string > > branch_points;

	// Begin with the first sugar.
	ResidueType const & first_sugar_type( *residue_types.front() );
	ResidueOP first_sugar( ResidueFactory::create_residue( first_sugar_type ) );
	string const upper_atom( first_sugar->carbohydrate_info()->anomeric_carbon_name() );
	pose.append_residue_by_atoms( *first_sugar, true, upper_atom, sequence_position, atom_name, true );

	// Build any other sugars.
	append_pose_with_glycan_residues( pose, ResidueTypeCOPs( residue_types.begin() + 1, residue_types.end() ) );

	// Let the Conformation know that it (now) contains sugars.
	pose.conformation().contains_carbohydrate_residues( true );

	// Finally, change the PDB information.
	// TODO: Be more intelligent about this by assigning new chain letters, etc.
	PDBInfoOP info( new PDBInfo( pose ) );
	info->name( pose.sequence() );  // Use the sequence as the default name.
	pose.pdb_info( info );

	TR << "Glycosylated pose with " << iupac_sequence << '-' << atom_name <<
		pose.residue( sequence_position ).name3() << sequence_position << endl;
}

// Glycosylate the Pose at the given sequence position using an IUPAC sequence.
/// @details  This is a wrapper function for standard AA cases, i.e., glycosylation at Asn, Thr, or Ser.
void
glycosylate_pose( Pose & pose, uint const sequence_position, std::string const & iupac_sequence )
{
	std::string const & glycosylation_site( pose.residue( sequence_position ).name3() );
	if ( glycosylation_site == "ASN" ) {
		glycosylate_pose( pose, sequence_position, "ND2", iupac_sequence );
	} else if ( glycosylation_site == "SER" ) {
		glycosylate_pose( pose, sequence_position, "OG", iupac_sequence );
	} else if ( glycosylation_site == "THR" ) {
		glycosylate_pose( pose, sequence_position, "OG1", iupac_sequence );
		// TODO: Add Trp, after creating an appropriate patch file.
	} else {
		utility_exit_with_message( glycosylation_site + " is not a common site of glycosylation or else it is "
			"ambiguous; Rosetta cannot determine attachment atom.  Use glycosylate_pose( Pose & pose, uint const "
			"sequence_position, std::string const & atom_name, std::string const & iupac_sequence ) instead." );
	}
}


// Glycosylate the Pose at the given sequence position and atom using a .GWS or IUPAC sequence file.
void
glycosylate_pose_by_file(
	Pose & pose,
	uint const sequence_position,
	std::string const & atom_name,
	std::string const & filename )
{
	using namespace std;
	using namespace io::carbohydrates;

	string const & sequence_file( find_glycan_sequence_file( filename ) );
	if ( sequence_file.find( ".gws" ) != string::npos ) {
		utility_exit_with_message( "Rosetta cannot yet read .gws format sequence files." );  // TODO: Change this.
		//string const & gws_sequence( read_glycan_sequence_file( sequence_file ) );
	} else {  // Assume any other file type contains an IUPAC sequence.
		string const & iupac_sequence( read_glycan_sequence_file( sequence_file ) );
		glycosylate_pose( pose, sequence_position, atom_name, iupac_sequence );
	}
}

// Glycosylate the Pose at the given sequence position using a .GWS or IUPAC sequence file.
/// @details  This is a wrapper function for standard AA cases, i.e., glycosylation at Asn, Thr, or Ser.
void
glycosylate_pose_by_file( Pose & pose, uint const sequence_position, std::string const & filename )
{
	std::string const & glycosylation_site( pose.residue( sequence_position ).name3() );
	if ( glycosylation_site == "ASN" ) {
		glycosylate_pose_by_file( pose, sequence_position, "ND2", filename );
	} else if ( glycosylation_site == "SER" ) {
		glycosylate_pose_by_file( pose, sequence_position, "OG", filename );
	} else if ( glycosylation_site == "THR" ) {
		glycosylate_pose_by_file( pose, sequence_position, "OG1", filename );
		// TODO: Add Trp, after creating an appropriate patch file.
	} else {
		utility_exit_with_message( glycosylation_site + " is not a common site of glycosylation or else it is "
			"ambiguous; Rosetta cannot determine attachment atom.  Use glycosylate_pose_by_file( Pose & pose, uint "
			"const sequence_position, std::string const & atom_name, std::string const & filename ) instead." );
	}
}

}  // namespace carbohydrates
}  // namespace pose
}  // namespace core
