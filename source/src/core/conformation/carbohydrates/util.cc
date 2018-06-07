// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/carbohydrates/util.cc
/// @brief Utility functions that DO NOT require a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Labonte <JWLabonte@jhu.edu>


// Unit Header
#include <core/conformation/carbohydrates/util.hh>

// Package Headers
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/carbohydrates/database_io.hh>

// Project Headers
#include <core/id/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/types.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

// Utility Headers
#include <utility/string_constants.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

// Numeric Headers
#include <numeric/conversions.hh>
#include <numeric/angle.functions.hh>
#include <numeric/random/random.hh>

//Options
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// External Headers
#include <boost/lexical_cast.hpp>
#include <basic/basic.hh>

// C++ Header
#include <list>
#include <utility>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.conformation.carbohydrates.util" );


namespace core {
namespace conformation {
namespace carbohydrates {

// Helper Functions ///////////////////////////////////////////////////////////
// Use a saccharide residue's connections to find the residue from which it follows or branches.
/// @return  The sequence position of the residue before this one (n-1) or the residue in the parent chain from which
/// the branch occurs or zero if N/A, i.e., if this is the lower terminus.
core::uint
find_seqpos_of_saccharides_parent_residue( conformation::Residue const & residue )
{
	debug_assert( residue.is_carbohydrate() );


	if ( ! residue.is_lower_terminus() ) {
		uint const id_of_connection_to_parent(
			residue.type().residue_connection_id_for_atom( residue.carbohydrate_info()->anomeric_carbon_index() ) );
		return residue.residue_connection_partner( id_of_connection_to_parent );
	} else  {
		TR.Debug << "This residue is a lower terminus! Returning 0." << std::endl;
		return 0;
	}


	//JAB - this fails during pose loading, even though the residue types should be finalized.
	// Not sure exactly why this would fail.  So, for now, we use the original code.
	core::Size anomeric_carbon = residue.carbohydrate_info()->anomeric_carbon_index();
	uint const id_of_connection_to_parent(
		residue.type().residue_connection_id_for_atom( anomeric_carbon ) );


	return residue.residue_connection_partner( id_of_connection_to_parent );

}


core::uint
find_seqpos_of_saccharides_mainchain_child( conformation::Residue const & residue )
{
	core::uint linkage_position = residue.carbohydrate_info()->mainchain_glycosidic_bond_acceptor();
	return find_seqpos_of_saccharides_child_residue_at( residue, linkage_position);
}


// TODO: What if this is is a sialic acid as reducing end?
// Use a saccharide residue's connections to find the residue following it from a given linkage position.
/// @param   <linkage_position>: an integer n of (1->n) of polysaccharide nomenclature, where n specifies the attachment
/// point on the parent monosaccharide residue; e.g., 4 specifies O4
/// @return  The sequence position of the residue downstream of this one attached to the given linkage positions.  This
/// is n+1, if the linkage position is the same as the main chain connectivity, or zero if N/A, i.e., if this is the
/// upper terminus.
core::uint
find_seqpos_of_saccharides_child_residue_at( conformation::Residue const & residue, core::uint linkage_position )
{
	using namespace std;
	using namespace chemical::carbohydrates;

	debug_assert( residue.is_carbohydrate() );
	debug_assert( linkage_position <= CarbohydrateInfo::MAX_C_SIZE_LIMIT );

	CarbohydrateInfoCOP info( residue.carbohydrate_info() );

	if ( residue.is_upper_terminus() ) {
		TR.Debug << "Residue " << residue.seqpos() << " is an upper terminus! Returning 0." << endl;
		return 0;
	} else if ( linkage_position > info->n_carbons() ) {
		TR.Debug << "Residue " << residue.seqpos() << " does not have a position " << linkage_position << '!' << endl;
		return 0;
	} else {
		string atom_name( "O" + string( 1, ( '0' + linkage_position ) ) );  // to get "O1", "O2", etc.

		if ( ! residue.has( atom_name ) ) {
			TR.Warning << "Residue " << residue.seqpos() << " does not have an " << atom_name << endl;
			return 0;
		}

		TR.Debug << "Finding seqpos of child residue of " << residue.seqpos() << " at " << atom_name;

		uint index( residue.atom_index( atom_name ) );
		if ( info->cyclic_oxygen_index() == index ) {  // This means that this is an exocyclic linkage.
			atom_name = "O" + string( 1, ( '0' + linkage_position + 1 ) );  // Try the next oxygen instead.
			index = residue.atom_index( atom_name );
			TR.Debug << "; cyclic oxygen, finding " << atom_name << " instead.";
		}

		TR.Debug << endl;

		chemical::ResidueType const & type( residue.type() );
		if ( type.atom_forms_residue_connection( index ) ) {
			uint const id_of_connection_to_child( type.residue_connection_id_for_atom( index ) );
			return residue.residue_connection_partner( id_of_connection_to_child );
		} else {
			return 0;
		}
	}
}


core::uint
get_linkage_position_of_saccharide_residue( conformation::Residue const & rsd, conformation::Residue const & parent_rsd )
{
	using namespace std;

	if ( rsd.seqpos() == parent_rsd.seqpos() ) {  // This occurs when there is no parent residue.
		TR.Debug << "This residue is a lower terminus! Returning 0." << endl;
		return 0;
	}
	if ( ! parent_rsd.is_carbohydrate() ) {
		TR.Debug << "Parent is a non-saccharide. Returning 0." << endl;
		return 0;
	}
	if ( rsd.seqpos() == parent_rsd.seqpos() + 1 ) {
		// We have a main chain connection. ( JAB - can we trust this numbering )?

		return parent_rsd.carbohydrate_info()->mainchain_glycosidic_bond_acceptor();
	}

	// If we get this far, we have a branch and need more information.

	//JAB - this looks OK to me - Jason - please check this over!
	// This will be replaced with a graph-search of the ResidueType once I know how to do that.
	core::Size parent_resnum = parent_rsd.seqpos();
	core::Size parent_connect = 0;
	for ( core::Size con = 1; con <= rsd.n_possible_residue_connections(); ++con ) {
		if ( rsd.connected_residue_at_resconn( con ) == parent_resnum ) {

			parent_connect = rsd.residue_connection_conn_id( con );
			break;
		}
	}
	core::Size connect_atom = parent_rsd.residue_connect_atom_index( parent_connect );

	core::Size c_atom_index = 0;

	// Get all bonded neighbors of this atom from the ResidueType.
	// This is important.  Our upper neighbor Carbon is not present, and the only other bonded atom is the H that would
	// be a monosaccharide.  So, we are good to go here as long as other carbons wouldn't be coming off of the O in any
	// monosaccharid,e (which I don't think is possible to still have a glycan linkage?) ~JAB
	utility::vector1< core::Size > bonded_atoms = parent_rsd.bonded_neighbor( connect_atom );
	for ( core::Size i  = 1; i <= bonded_atoms.size(); ++i ) {
		c_atom_index = bonded_atoms[ i ];
		std::string element = parent_rsd.atom_type( c_atom_index ).element();
		if ( element == "C" ) {
			break;
		}
	}

	std::string connect_atom_name = parent_rsd.atom_name( c_atom_index );
	return utility::string2Size(utility::to_string( connect_atom_name[ 2 ] ) );
}


// Get whether the glycosidic linkage between the residue and previous residue (parent residue) has an exocyclic carbon.
/// @details  Does not currently work for aa->glycan.  Returns false if previous residue is not carbohydrate.
/// @author   Jared Adolf-Bryfogle (jadolfbr@gmail.com)
bool
has_exocyclic_glycosidic_linkage( conformation::Conformation const & conf, uint const seqpos )
{
	conformation::Residue const & rsd = conf.residue( seqpos );
	core::Size lower_resnum = find_seqpos_of_saccharides_parent_residue( rsd );
	//Lowest of saccharide chain.  Return false.
	if ( lower_resnum == 0 ) {
		TR.Debug << "has_exocyclic_glycosidic_linkage: This residue is a lower terminus! Returning false." << std::endl;
		return false;
	}
	//TR << "lower resnum: " << lower_resnum << std::endl;

	conformation::Residue const & prev_rsd = conf.residue( lower_resnum );
	return has_exocyclic_glycosidic_linkage( rsd, prev_rsd );
}

// Get whether the glycosidic linkage between the residue and previous residue (parent residue) has an exocyclic carbon.
/// @details  Does not currently work for aa->glycan.  Returns false if previous residue is not carbohydrate.
/// @author   Jared Adolf-Bryfogle (jadolfbr@gmail.com)
bool
has_exocyclic_glycosidic_linkage( conformation::Residue const & rsd, conformation::Residue const & parent_rsd )
{
	// What does this mean for ASN-glycan connections?? Technically, it won't be an exocyclic atom -
	// but it WILL have omega and omega2 if ASN - so be careful here!
	if ( ! parent_rsd.is_carbohydrate() ) {
		TR.Debug << "has_exocyclic_glycosidic_linkage: Previous residue is not a carbohydrate! Returning false. " << std::endl;
		return false;
	}

	core::Size const n_carbons = parent_rsd.carbohydrate_info()->n_carbons();
	core::Size linkage_position = get_linkage_position_of_saccharide_residue( rsd, parent_rsd );
	core::Size last_carbon = parent_rsd.carbohydrate_info()->last_carbon_in_ring();

	if ( ( n_carbons == linkage_position ) && ( last_carbon != linkage_position ) ) {
		return true;
	} else {
		return false;
	}
}

// Return pointers to the two residues of the glycosidic bond.
/// @return  Pointers to the residue at <sequence_position> and its parent or else the same pointer twice if undefined.
std::pair< conformation::ResidueCOP, conformation::ResidueCOP >
get_glycosidic_bond_residues( Conformation const & conf, uint const sequence_position )
{
	using namespace std;
	using namespace conformation;

	// Get the 1st residue of interest.
	ResidueCOP res_n( conf.residue( sequence_position ).get_self_ptr() );

	if ( res_n->is_lower_terminus() ) {
		if ( TR.Debug.visible() ) {
			TR.Debug << "Glycosidic torsions are undefined for the first polysaccharide residue of a chain unless part "
				"of a branch." << endl;
		}
		return make_pair( res_n, res_n );
	}

	// Get the 2nd residue of interest.
	// (res_n_minus_1 is a misnomer for the lower termini of branches.)
	ResidueCOP res_n_minus_1;
	if ( conf.glycan_tree_set() ) {
		res_n_minus_1 = conf.residue( conf.glycan_tree_set()->get_parent( sequence_position ) ).get_self_ptr();
	} else {
		res_n_minus_1 = conf.residue( find_seqpos_of_saccharides_parent_residue( *res_n ) ).get_self_ptr();
	}

	return make_pair( res_n, res_n_minus_1 );
}

/// Return the AtomIDs of the four phi torsion reference atoms.
/// @details For aldopyranoses, phi is defined as O5(n)-C1(n)-OX(n-1)-CX(n-1),
/// where X is the position of the glycosidic linkage.\n
/// For aldofuranoses, phi is defined as O4(n)-C1(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketopyranoses, phi is defined as O6(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketofuranoses, phi is defined as O5(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// Et cetera...\n
utility::vector1< id::AtomID >
get_reference_atoms_for_phi( Conformation const & conf, uint const sequence_position )
{
	using namespace std;
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1< AtomID > ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair< ResidueCOP, ResidueCOP > const residues( get_glycosidic_bond_residues( conf, sequence_position ) );

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
		utility_exit_with_message("Linear Polysacharides are not yet handled!");  // TODO: Figure out how linear polysaccharides are handled by IUPAC.
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
	ids[ 4 ] = ref4;

	if ( TR.Trace.visible() ) {
		TR.Trace << "Reference atoms for phi: " << ref1 << ", " << ref2 << ", " << ref3 << ", " << ref4 << endl;
	}

	return ids;
}

// Return the AtomIDs of the four psi torsion reference atoms.
/// @details For saccharides, psi is defined as: C(anomeric)(n)-OX(n-1)-CX(n-1)-CX-1(n-1),\n
/// where X is the position of the glycosidic linkage.
utility::vector1< id::AtomID >
get_reference_atoms_for_psi( Conformation const & conf, uint const sequence_position )
{
	using namespace std;
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1< AtomID > ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair< ResidueCOP, ResidueCOP > const residues( get_glycosidic_bond_residues( conf, sequence_position ) );

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
get_reference_atoms_for_1st_omega( Conformation const & conf, uint const sequence_position )
{
	using namespace std;
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1< AtomID > ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair< ResidueCOP, ResidueCOP > const residues( get_glycosidic_bond_residues( conf, sequence_position ) );

	if ( residues.first->seqpos() == residues.second->seqpos() ) {  // This occurs when there is no parent residue.
		return ids;
	}


	//if ( residues.second->is_carbohydrate() && ( ! has_exocyclic_glycosidic_linkage( *residues.first, *residues.second ) ) ) {
	bool exocylic_linkage;
	if ( conf.glycan_tree_set() ) {
		exocylic_linkage = conf.glycan_tree_set()->has_exocyclic_glycosidic_linkage( residues.first->seqpos() );
	} else {
		exocylic_linkage = has_exocyclic_glycosidic_linkage(conf, residues.first->seqpos() );
	}

	if ( residues.second->is_carbohydrate() && ( ! exocylic_linkage) ) {
		TR.Debug << "Omega is undefined for residue " << sequence_position <<
			" because the glycosidic linkage is not exocyclic." << endl;
		return ids;
	}

	ids.resize( 4 );  // A torsion has 4 reference atoms.

	// Set the atom names of the four reference atoms.
	// Reference 1 is OX(n-1) for polysaccharides.
	AtomID const ref1( residues.second->connect_atom( *residues.first ), residues.second->seqpos() );
	ids[ 1 ] = ref1;

	// Reference 2 is CX(n-1) for polysaccharides.
	AtomID const ref2( residues.second->type().atom_base( ref1.atomno() ), residues.second->seqpos() );
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

// Return the AtomIDs of the four omega2 torsion reference atoms.
/// For carbohydrates glycosylated at positions with more than 1 exocyclic
/// position, or for linkages to longer amino acid residues,
/// omega2 of residue n is defined as CX(n-1)-CX-1(n-1)-CX-2(n-1)-CX-3(n-1),
/// where X is the position of the glycosidic linkage.
utility::vector1< id::AtomID >
get_reference_atoms_for_2nd_omega( Conformation const & conf, uint const sequence_position )
{
	using namespace std;
	using namespace id;
	using namespace utility;
	using namespace conformation;

	vector1< AtomID > ids;

	// Get the two residues.  (The first is the "current" residue; the second is the parent.)
	pair< ResidueCOP, ResidueCOP > const residues( get_glycosidic_bond_residues( conf, sequence_position ) );

	if ( residues.first->seqpos() == residues.second->seqpos() ) {  // This occurs when there is no parent residue.
		return ids;
	}

	bool exocyclic_linkage;
	if ( conf.glycan_tree_set() ) {
		exocyclic_linkage = conf.glycan_tree_set()->has_exocyclic_glycosidic_linkage( residues.first->seqpos() );
	} else {
		exocyclic_linkage = has_exocyclic_glycosidic_linkage(conf, residues.first->seqpos() );
	}

	if ( residues.second->is_carbohydrate() && ( ! exocyclic_linkage ) ) {
		TR.Debug << "Omega2 is undefined for residue " << sequence_position <<
			" because the glycosidic linkage is not exocyclic." << endl;
		return ids;
	}

	// Set the atom names of the four reference atoms.
	// Reference 0 is OX(n-1) for polysaccharides.
	AtomID const ref0( residues.second->connect_atom( *residues.first ), residues.second->seqpos() );

	// Reference 1 is CX(n-1) for polysaccharides.
	AtomID const ref1( residues.second->type().atom_base( ref0.atomno() ), residues.second->seqpos() );


	// Reference 2 is CX-1(n-1) for polysaccharides.
	AtomID const ref2( residues.second->type().atom_base( ref1.atomno() ), residues.second->seqpos() );

	//If Ref2 is a ring atom and we do not have an omega2 angle
	if ( residues.second->is_carbohydrate() && residues.second->type().is_ring_atom( 1, ref2.atomno() ) ) {
		return ids;
	} else if ( (! residues.second->is_carbohydrate() ) && (residues.second->aa() == core::chemical::aa_ser || residues.second->aa() == core::chemical::aa_thr ) ) {
		//O-linked Glycosylation has 3 dihedrals.  We should address this better.
		return ids;
	}

	ids.resize( 4 );  // A torsion has 4 reference atoms.
	ids[ 1 ] = ref1;
	ids[ 2 ] = ref2;

	// Reference 3 is CX-2(n-1) for polysaccharides.
	AtomID const ref3( residues.second->type().atom_base( ref2.atomno() ), residues.second->seqpos() );
	ids[ 3 ] = ref3;

	// Reference 4 is CX-3(n-1) for polysaccharides.
	AtomID const ref4( residues.second->type().atom_base( ref3.atomno() ), residues.second->seqpos() );
	ids[ 4 ] = ref4;

	if ( TR.Debug.visible() ) {
		TR.Debug << "Reference atoms for omega2: " << ref1 << ", " << ref2 << ", " << ref3 << ", " << ref4 << endl;
	}

	return ids;
}


// Return the AtomIDs of the four reference atoms for the requested torsion.
// Works for AA->glycan connection as well
utility::vector1< id::AtomID >
get_reference_atoms( uint const torsion_id, Conformation const & conf, uint const sequence_position )
{
	using namespace id;
	using namespace utility;

	vector1< AtomID > ref_atoms;
	switch ( torsion_id ) {
	case phi_torsion :
		ref_atoms = get_reference_atoms_for_phi( conf, sequence_position );
		break;
	case psi_torsion :
		ref_atoms = get_reference_atoms_for_psi( conf, sequence_position );
		break;
	case omega_torsion :
		ref_atoms = get_reference_atoms_for_1st_omega( conf, sequence_position );
		break;
	case 4 :
		ref_atoms = get_reference_atoms_for_2nd_omega( conf, sequence_position );
		break;
	default :
		utility_exit_with_message( "An invalid torsion angle was requested." );
	}
	return ref_atoms;
}


/// @brief   Get the main chain or branch TorsionID defined by these four AtomIDs.
/// @params  <atoms>: A vector1 of AtomIDs with a size of 4.
/// @return  A TorsionID guaranteed to have a BB, CHI, or BRANCH TorsionType.
/// @note    Sometimes it is helpful to know which TorsionID corresponds to which four atoms.
/// Since some torsions have multiple IDs, (e.g., for a aldohexopyranose, BB 2 and NU 2 are the same,) one cannot
/// write such a function that returns a single TorsionID, hence the specificity of this function to return either a
/// BB or BRANCH TorsionType or a CHI torsion defining a connection to a branch.\n
/// This code is probably general enough for non-sugar residues, but I'm rushing it in for CAPRI Round 41 and keeping it
/// here for now, just to be safe.
/// @author  Labonte <JWLabonte@jhu.edu>
core::id::TorsionID
get_non_NU_TorsionID_from_AtomIDs( Conformation const & conf, utility::vector1< core::id::AtomID > const & atoms )
{
	using namespace std;
	using namespace utility;
	using namespace id;

	// First, we assure that we have four atoms, since four atoms define a torsion angle.
	if ( atoms.size() != 4 ) { return TorsionID::BOGUS_TORSION_ID(); }

	// Second, if passed a bogus AtomID, return a bogus TorsionID.
	if ( ! atoms[ 1 ].valid() ) { return TorsionID::BOGUS_TORSION_ID(); }

	// Third, check to see if the order of atoms needs to be reversed.
	// Then, get the desired atoms.
	bool const reversed( atoms[ 4 ] < atoms[ 1 ] );

	AtomID const & upstream_bond_atom_parent( ( ! reversed ) ? atoms[ 1 ] : atoms[ 4 ] );
	AtomID const & upstream_bond_atom( ( ! reversed ) ? atoms[ 2 ] : atoms[ 3 ] );
	AtomID const & downstream_bond_atom( ( ! reversed ) ? atoms[ 3 ] : atoms[ 2 ] );
	AtomID const & downstream_bond_atom_child( ( ! reversed ) ? atoms[ 4 ] : atoms[ 1 ] );

	// Prepare variables for the three parameters of a TorsionID.
	// (The torsion always "belongs" to the same residue as the second atom.)
	uint const resnum( upstream_bond_atom.rsd() );
	Residue const & parent( conf.residue( resnum ) );

	TorsionType type;
	uint torsionnum( 0 );

	if ( upstream_bond_atom_parent.rsd() != resnum ) {
		// In this case, we must be at the first BB torsion.
		// e.g., the phi of an AA has LOWER as the first atom of the definition.

#ifndef NDEBUG
		// Make sure all is right with the world.
		vector1< uint > const & mainchain_atoms( parent.mainchain_atoms() );
		if ( mainchain_atoms[ 2 ] != upstream_bond_atom.atomno() ||
				mainchain_atoms[ 3 ] != downstream_bond_atom.atomno() ||
				mainchain_atoms[ 4 ] != downstream_bond_atom_child.atomno() ) {
			utility_exit_with_message( "get_non_NU_TorsionID_from_AtomIDs():"
				"Four atoms with the pattern res(n-1) atom, res(n) atom, res(n) atom, res(n) atom "
				"passed to this function, which can only be a 1st BB torsion, yet the 2nd, 3rd, and 4th atoms "
				"do not match the atoms stored in the ResidueType for BB1!" );
		}
#endif

		type = BB;
		torsionnum = 1;
	} else if ( upstream_bond_atom_parent.rsd() == downstream_bond_atom.rsd() ) {
		// If both atoms across the bond are on the same residue, this is NOT a BRANCH torsion.
		// Search the residue's main-chain atoms to find a match, followed by searching the chi atoms.

		bool torsionnum_found( false );

		vector1< uint > const & mainchain_atoms( parent.mainchain_atoms() );
		Size const n_mainchain_atoms( mainchain_atoms.size() );
		// Start at 2; BB1 covered above.
		for ( uint i( 2 ); i < n_mainchain_atoms; ++i ) {
			if ( mainchain_atoms[ i - 1 ] == upstream_bond_atom_parent.atomno() &&
					mainchain_atoms[ i ] == upstream_bond_atom.atomno() &&
					mainchain_atoms[ i + 1 ] == downstream_bond_atom.atomno() ) {
				type = BB;
				torsionnum = i;
				torsionnum_found = true;
				break;
			}
		}

		if ( ! torsionnum_found ) {
			vector1< vector1< uint > > const & chi_atoms( parent.chi_atoms() );
			Size const n_chis( chi_atoms.size() );
			for ( uint i( 1 ); i <= n_chis; ++i ) {
				if ( chi_atoms[ i ][ 1 ] == upstream_bond_atom_parent.atomno() &&
						chi_atoms[ i ][ 2 ] == upstream_bond_atom.atomno() &&
						chi_atoms[ i ][ 3 ] == downstream_bond_atom.atomno() ) {
					type = CHI;
					torsionnum = i;
					break;
				}
			}
		}
	} else /* the two atoms across the bond are not on the same residue */ {
		// If both atoms across the bond are NOT on the same residue,
		// this could be the final BB torsion or else a BRANCH torsion.

		// See if it is the final BB torsion.
		vector1< uint > const & mainchain_atoms( parent.mainchain_atoms() );
		if ( mainchain_atoms[ mainchain_atoms.size() - 1 ] == upstream_bond_atom_parent.atomno() &&
				mainchain_atoms[ mainchain_atoms.size() ] == upstream_bond_atom.atomno() ) {
			type = BB;
			torsionnum = parent.mainchain_atoms().size();
		} else {  // It must be a branch.
			type = BRANCH;

			// Now figure out the BRANCH number from the ResidueConnection.
			Size const n_polymer_connections( parent.n_polymeric_residue_connections() );
			for ( uint i( n_polymer_connections + 1 );
					i <= parent.n_possible_residue_connections(); ++i ) {
				if ( parent.connected_residue_at_resconn( i ) == downstream_bond_atom.rsd() ) {
					torsionnum = i - n_polymer_connections;
					break;
				}
			}
		}
	}

	if ( ! torsionnum ) {
		TR.Error << "get_non_NU_TorsionID_from_AtomIDs() could not determine the BB id for the";
		TR.Error << "TorsionID defined by these atoms:";
		TR.Error << upstream_bond_atom_parent << "; " << upstream_bond_atom << "; ";
		TR.Error << downstream_bond_atom << "; " << downstream_bond_atom_child << endl;
	}

	return TorsionID( resnum, type, torsionnum );
}

/// @brief   Get a list of the TorsionIDs of all glycosidic torsions for the residue at this position.
/// @author  Labonte <JWLabonte@jhu.edu>
utility::vector1< core::id::TorsionID >
get_glycosidic_TorsionIDs( core::conformation::Conformation const & conf, uint const seq_pos )
{
	using namespace core::id;
	using namespace core::conformation::carbohydrates;

	debug_assert( conf.residue( seq_pos ).is_carbohydrate() );

	utility::vector1< TorsionID > torsions;
	for ( core::uint i( 1 ); i <= 4; ++i ) {  // 4 covers phi, psi, and up to two omegas.
		TorsionID const torsion(
			get_non_NU_TorsionID_from_AtomIDs( conf, get_reference_atoms( i, conf, seq_pos ) ) );
		if ( torsion == TorsionID::BOGUS_TORSION_ID() ) { break; }
		torsions.push_back( torsion );
	}

	return torsions;
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
align_virtual_atoms_in_carbohydrate_residue( conformation::Conformation & conf, uint const sequence_position )
{
	using namespace std;
	using namespace id;
	using namespace conformation;

	TR.Debug << " Aligning virtual atoms on residue " << sequence_position << "..." << endl;

	ResidueCOP res( conf.residue( sequence_position ).get_self_ptr() );

	core::Real offset_value = basic::options::option[ basic::options::OptionKeys::in::glycan_virtual_offset ].value();
	numeric::xyzVector< core::Real > offset(offset_value,offset_value,offset_value);
	// Find and align VOX, if applicable.
	if ( res->carbohydrate_info()->is_cyclic() ) {
		if ( TR.Debug.visible() ) {
			TR.Debug << "  Aligning VOX..." << endl;
		}
		uint const x( res->carbohydrate_info()->cyclic_oxygen() );
		uint const OX( res->atom_index( res->carbohydrate_info()->cyclic_oxygen_name() ) );
		uint const VOX( res->atom_index( "VO" + string( 1, x + '0' ) ) );

		numeric::xyzVector< core::Real > coord = conf.xyz( AtomID( OX, sequence_position ) );
		//coord+=offset; don't offset VOX for now
		conf.set_xyz( AtomID( VOX, sequence_position ), coord );
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

		numeric::xyzVector< core::Real > coord = conf.xyz( AtomID( HOY_ref, parent_res_seqpos ) );
		coord+=offset;
		conf.set_xyz( AtomID( HOY, sequence_position ), coord );
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

		coord = conf.xyz( AtomID( OY_ref, parent_res_seqpos ) );
		coord+=offset;
		conf.set_xyz( AtomID( OY, sequence_position ), coord );
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

		numeric::xyzVector< core::Real > coord = conf.xyz( AtomID( HOZ_ref, downstream_res_seqpos ) );
		coord+=offset;
		conf.set_xyz( AtomID( HOZ, sequence_position ), coord );
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
		if ( branch_res_seqpos != 0 ) {
			ResidueCOP branch_res( conf.residue( branch_res_seqpos ).get_self_ptr() );
			uint const HOZ_ref( branch_res->atom_index( branch_res->carbohydrate_info()->anomeric_carbon_name() ) );
			numeric::xyzVector< core::Real > coord = conf.xyz( AtomID( HOZ_ref, branch_res_seqpos ) );
			coord+=offset;
			conf.set_xyz( AtomID( HOZ, sequence_position ), coord );
		} else {
			TR << "Branch connection " << branch_connection_id << " on residue " << res->seqpos() << " is not chemically connected." << std::endl;
		}
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << " All virtual atoms aligned." << endl;
	}
}


// TorsionID Queries //////////////////////////////////////////////////////////
// Is this is the phi torsion angle of a glycosidic linkage?
/// @details  Carbohydrate linkages are defined as the torsion angles leading back to the previous residue.  Much
/// of Rosetta code relies on TorsionIDs and assumes TorsionID( n, BB, 1 ) is phi.  For a sugar, phi (of the next
/// residue) is the last torsion, and the number of main-chain torsions varies per saccharide residue.
bool
is_glycosidic_phi_torsion( Conformation const & conf, id::TorsionID const & torsion_id )
{
	using namespace id;

	if ( torsion_id.type() == BB || torsion_id.type() == BRANCH ) {  // Phi for a branched saccharide has a BRANCH type.
		conformation::Residue const & residue( conf.residue( torsion_id.rsd() ) );
		uint next_rsd_num( 0 );  // We will need to see if the "next" residue is a saccharide.

		switch( torsion_id.type() ) {
		case  BB :
			if ( ! residue.is_upper_terminus() ) {
				// If this is a main-chain torsion, we need the next residue on the main chain.
				next_rsd_num = torsion_id.rsd() + 1;
				if ( conf.residue( next_rsd_num ).is_carbohydrate() ) {
					return ( torsion_id.torsion() == residue.n_mainchain_atoms() );  // The last BB is phi.
				}
			}
			break;
		case BRANCH :
			if ( torsion_id.torsion() <= residue.n_non_polymeric_residue_connections() ) {
				Size const n_mainchain_connections( residue.n_polymeric_residue_connections() );
				next_rsd_num = residue.residue_connection_partner( n_mainchain_connections + torsion_id.torsion() );
				if ( conf.residue( next_rsd_num ).is_carbohydrate() ) {
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
/// @note  This function will return false if the TorsionID is a CHI and the torsion is already covered by an equivalent
/// BB.  In other words, TorsionIDs for CHIs only return true for branch connections, where the ONLY way to access that
/// torsion is through CHI.
bool
is_glycosidic_psi_torsion( Conformation const & conf, id::TorsionID const & torsion_id )
{
	using namespace id;

	if ( torsion_id.type() == BB || torsion_id.type() == CHI ) {  // Psi for a branched saccharide has a CHI type.
		conformation::Residue const & residue( conf.residue( torsion_id.rsd() ) );
		uint next_rsd_num( 0 );  // We will need to see if the "next" residue is a saccharide.

		switch( torsion_id.type() ) {
		case BB :
			if ( ! residue.is_upper_terminus() ) {
				// If this is a main-chain torsion, we need the next residue on the main chain.
				next_rsd_num = torsion_id.rsd() + 1;
				if ( conf.residue( next_rsd_num ).is_carbohydrate() ) {
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
				conformation::Residue const & next_rsd( conf.residue( next_rsd_num ) );
				if ( ( next_rsd.is_carbohydrate() ) && ( residue.connect_atom( next_rsd ) == third_atom ) ) {
					return true;
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
/// @note     This function currently will not recognize additional omegas beyond omega1.\n
/// Also, this function will return false if the TorsionID is a CHI and the torsion is already covered by an equivalent
/// BB.  In other words, TorsionIDs for CHIs only return true for branch connections, where the ONLY way to access that
/// torsion is through CHI.
bool
is_glycosidic_omega_torsion( Conformation const & conf, id::TorsionID const & torsion_id )
{
	using namespace id;

	if ( torsion_id.type() == BB || torsion_id.type() == CHI ) {  // Omega for a branched saccharide has a CHI type.
		conformation::Residue const & residue( conf.residue( torsion_id.rsd() ) );
		uint next_rsd_num( 0 );  // We will need to see if the "next" residue is a saccharide.

		switch( torsion_id.type() ) {
		case BB :
			if ( ! residue.is_upper_terminus() ) {
				// If this is a main-chain torsion, we need the next residue on the main chain.
				next_rsd_num = torsion_id.rsd() + 1;
				if ( conf.residue( next_rsd_num ).is_carbohydrate() ) {
					chemical::carbohydrates::CarbohydrateInfoCOP info( residue.carbohydrate_info() );
					if ( info->has_exocyclic_linkage_to_child_mainchain() ) {
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
				conformation::Residue const & next_rsd( conf.residue( next_rsd_num ) );
				if ( ( next_rsd.is_carbohydrate() ) && ( residue.connect_atom( next_rsd ) == fourth_atom ) ) {
					return true;
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

// Return the sequence position of the immediate downstream (child) residue affected by this torsion.
/// @note  This method is designed for saccharide residues, because it makes assumptions about CHI and BRANCH numbering.
core::uint
get_downstream_residue_that_this_torsion_moves( Conformation const & conf, id::TorsionID const & torsion_id )
{
	using namespace id;

	conformation::Residue const & rsd( conf.residue( torsion_id.rsd() ) );
	uint next_rsd( 0 );
	if ( torsion_id.type() == BB ) {
		debug_assert( torsion_id.torsion() <= rsd.n_mainchain_atoms() );
		// If this is a main-chain torsion, we can be confident that the next residue on this chain MUST be n+1.
		next_rsd = torsion_id.rsd() + 1;
	} else if ( torsion_id.type() == BRANCH ) {
		debug_assert( torsion_id.torsion() <= rsd.n_non_polymeric_residue_connections() );
		Size const n_mainchain_connections( rsd.n_polymeric_residue_connections() );
		next_rsd = rsd.residue_connection_partner( n_mainchain_connections + torsion_id.torsion() );
	} else if ( torsion_id.type() == CHI ) {
		next_rsd = find_seqpos_of_saccharides_child_residue_at( rsd, torsion_id.torsion() );
	}
	if ( ! next_rsd ) {
		TR.Debug << "Torsion " << torsion_id << " does not move any downstream residues!  Returning 0." << std::endl;
	}
	return next_rsd;
}


// Torsion Access /////////////////////////////////////////////////////////////
// Getters ////////////////////////////////////////////////////////////////////

///@brief Get the number of glycosidic torsions for this residue.  Up to 4 (omega2).
Size
get_n_glycosidic_torsions_in_res( Conformation const & conf, uint const sequence_position )
{
	core::Size n_torsions( 0 );

	//std::cout << "seq position: " << sequence_position << std::endl;
	for ( core::Size torsion_id = 1; torsion_id <= 4; ++torsion_id ) {
		//std::cout << "Torsion: " << torsion_id << std::endl;
		utility::vector1< id::AtomID > const ref_atoms = get_reference_atoms( torsion_id, conf, sequence_position );
		//std::cout << "size: " << ref_atoms.size() << std::endl;
		if ( ref_atoms.size() != 0 ) {
			//std::cout << "adding torsion: "<< torsion_id << std::endl;
			n_torsions+=1;
		}
	}

	//std::cout << "n_torsions: " << n_torsions << std::endl;
	return n_torsions;
}


// Return the requested torsion angle between a saccharide residue of the given pose and the previous residue.
/// @details This method is used in place of Residue::mainchain_torsion() since the main-chain torsions of saccharide
/// residues only make sense in the context of two residues.  Moreover, the reference atoms as defined by the IUPAC are
/// different from the ones that Rosetta uses by default for mainchain torsions for sugars.
/// @param   <named_torsion> is an integer representing the specific torsion angle requested, as defined in
/// core/id/types.hh:\n
///     phi_torsion = 1\n
///     psi_torsion = 2\n
///     omega_torsion = 3
/// @note    I would rather named_torsion were an enum, but as it was already defined, I'm leaving it as a constant
/// for now.
core::Angle
get_glycosidic_torsion( uint const named_torsion, Conformation const & conf, uint const sequence_position )
{
	using namespace id;
	using namespace numeric;
	using namespace utility;

	utility::vector1< id::AtomID > const ref_atoms = get_reference_atoms( named_torsion, conf, sequence_position );

	if ( ref_atoms.size() == 0 ) {
		// This occurs when there is no parent residue or when the glycosidic bond is not exocyclic (omega only).
		TR.Debug << "Returning zero." << std::endl;
		return 0.0;
	}

	Angle const angle_in_radians(
		conf.torsion_angle( ref_atoms[ 1 ], ref_atoms[ 2 ], ref_atoms[ 3 ], ref_atoms[ 4 ]) );
	return principal_angle_degrees( conversions::degrees( angle_in_radians ) );
}


// Setters ////////////////////////////////////////////////////////////////////
// Set the requested torsion angle between a saccharide residue of the given pose and the previous residue.
/// @details This method is used in place of Conformation::set_torsion() since the reference atoms as defined by the
/// IUPAC are different from the ones that Rosetta uses by default for main-chain torsions for sugars.
/// @param   <named_torsion> is an integer representing the specific torsion angle requested, as defined in
/// core/id/types.hh:\n
///     phi_torsion = 1\n
///     psi_torsion = 2\n
///     omega_torsion = 3\n
/// <setting> is in degrees.
/// @note    I would rather named_torsion were an enum, but as it was already defined, I'm leaving it as a constant
/// for now.
void
set_glycosidic_torsion(
	uint const named_torsion,
	Conformation & conf,
	uint const sequence_position,
	core::Angle const setting )
{
	using namespace id;
	using namespace numeric;
	using namespace utility;

	vector1< AtomID > const ref_atoms = get_reference_atoms( named_torsion, conf, sequence_position );

	if ( ref_atoms.size() == 0 ) {
		// This occurs when there is no parent residue or when the glycosidic bond is not exocyclic (omega only).
		return;
	}

	Angle const setting_in_radians( conversions::radians( setting ) );
	conf.set_torsion_angle(
		ref_atoms[ 1 ], ref_atoms[ 2 ], ref_atoms[ 3 ], ref_atoms[ 4 ], setting_in_radians );
	//align_virtual_atoms_in_carbohydrate_residue( conf, sequence_position );
}


/// @brief  Recursive function to get branches of a set of residues, etc.
///  list_of_residues and tips are arrays are non-const references and modified by this function.
///
///  Children Residues:  Residue nums of parent residue connected that we are interested in finding connected branchs.
///  List Of  Residues:  All the residue nums of the branching from children residues
///  Tips:  All 'ends' of all the branches found using this function.
///
///  See Also: get_carbohydrate_residues_and_tips_of_branch
///            trim_carbohydrate_branch_from_X
void
get_branching_residues( conformation::Conformation const & conf,
	Size parent_residue,
	utility::vector1< Size > & children_residues,
	utility::vector1< Size > & list_of_residues,
	utility::vector1< Size > & tips,
	std::set< Size > const & ancestors )
{
	for ( core::Size i =1; i <= children_residues.size(); ++i ) {
		Size res = children_residues[ i ];
		if ( ancestors.count( res ) == 1 ) {
			TR.Debug << "In get_branching_residues, skipping residue " << res << " because it leads to a cycle." << std::endl;
			continue;
		}
		utility::vector1< Size > children;
		fill_upstream_children_res_and_tips( conf, res, parent_residue, children, list_of_residues, tips );

		if ( children.size() != 0 ) {
			std::set< Size > childs_ancestors( ancestors );
			childs_ancestors.insert( res );
			get_branching_residues( conf, res, children, list_of_residues, tips, childs_ancestors);
		}
	}
}

/// @brief  Find all children residues, list of residues, and any found tips from a given residue not including parent
///
///  Children Residues:  Filled in list of children residues found if not tips.
///  List Of  Residues:  All the residue nums found.
///  Tips:  All 'ends' of of children found.
///
///  See Also: get_carbohydrate_residues_and_tips_of_branch
///            trim_carbohydrate_branch_from_X
void
fill_upstream_children_res_and_tips( conformation::Conformation const & conf,
	Size res,
	Size parent_residue,
	utility::vector1< Size > & children_residues,
	utility::vector1< Size > & list_of_residues,
	utility::vector1< Size > & tips )
{
	Size connections = conf.residue( res ).n_possible_residue_connections(); //Want the index to match here.
	for ( core::Size connection = 1; connection <= connections; ++connection ) {

		Size connecting_res = conf.residue( res ).connected_residue_at_resconn( connection );

		if ( connecting_res != 0 && connecting_res != parent_residue && conf.residue( connecting_res ).is_carbohydrate() ) {
			list_of_residues.push_back( connecting_res );
			if ( conf.residue( connecting_res ).n_current_residue_connections() == 1 ) {
				tips.push_back( connecting_res );
			} else {
				children_residues.push_back( connecting_res );
			}
		}
	}
}

core::Size
get_glycan_tree_size( conformation::Conformation const & conf, core::Size const first_glycan_resnum ){
	utility::vector1< core::Size > glycan_resnums = get_carbohydrate_residues_of_branch(conf, first_glycan_resnum);
	return glycan_resnums.size() + 1;
}

core::Size
get_largest_glycan_tree_size( conformation::Conformation const & conf ){

	utility::vector1< core::Size > tree_sizes;
	utility::vector1< bool > glycan_start_points = get_glycan_start_points( conf );
	for ( core::Size glycan_start = 1; glycan_start <= glycan_start_points.size(); ++glycan_start ) {
		if ( glycan_start_points[ glycan_start ] ) {
			core::Size glycan_length = get_glycan_tree_size( conf, glycan_start );
			tree_sizes.push_back( glycan_length );
		}
	}
	return utility::max( tree_sizes );
}

core::Size
get_distance_to_start( conformation::Conformation const & conf, core::Size const position){
	core::Size res_distance = 0;
	core::Size parent = find_seqpos_of_saccharides_parent_residue(conf.residue(position));
	if ( parent == 0 ) {
		return 0;
	}
	while ( ( parent != 0 ) && ( conf.residue(parent).is_carbohydrate() ) ) {
		res_distance+=1;
		parent = find_seqpos_of_saccharides_parent_residue(conf.residue(parent));
	}
	return res_distance;
}

utility::vector1< bool >
get_glycan_start_points(conformation::Conformation const & conf){

	utility::vector1< bool > glycan_start_points(conf.size(), false);

	for ( core::Size i = 1; i <= conf.size(); ++i ) {

		//Branch point - definitely a start the glycan
		if ( conf.residue( i ).is_branch_point() && ! conf.residue(i).is_carbohydrate() ) {
			Size glycan_plus_one = get_glycan_connecting_protein_branch_point(conf, i);
			if ( glycan_plus_one != 0 ) {
				glycan_start_points[ glycan_plus_one ] = true;
			}
		} else if ( conf.residue( i ).is_carbohydrate() ) {
			//Indicates a glycan that is not part of a protein chain. - SO we make sure the parent residue is 0 to indicate the start point.
			if ( find_seqpos_of_saccharides_parent_residue( conf.residue(i)) == 0 ) {
				glycan_start_points[ i ] = true;
			}
		}
	}
	return glycan_start_points;
}

/// @brief Get residues further down the branch from this residue.  starting_position ->
/// @details Calls get_carbohydrate_residues_and_tips_of_branch
utility::vector1< core::Size >
get_carbohydrate_residues_of_branch(
	conformation::Conformation const & conf,
	uint const starting_position)
{
	return get_carbohydrate_residues_and_tips_of_branch(conf, starting_position).first;
}

/// @brief Get tips (end residue of linear components of branches) further down the branch from this residue.  starting_position ->
/// @details Calls get_carbohydrate_residues_and_tips_of_branch
utility::vector1< core::Size >
get_carbohydrate_tips_of_branch(
	conformation::Conformation const & conf,
	uint const starting_position)
{
	return get_carbohydrate_residues_and_tips_of_branch(conf, starting_position).second;
}

std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > >
get_carbohydrate_residues_and_tips_of_branch(
	conformation::Conformation const & conf,
	uint const starting_position,
	bool include_starting_position /*false*/)
{
	using namespace core::chemical::carbohydrates;

	utility::vector1< Size > tips;
	utility::vector1< Size > list_of_residues;
	utility::vector1< Size > children_residues;

	if ( include_starting_position ) {
		list_of_residues.push_back(starting_position);
	}


	if ( ! conf.residue( starting_position ).is_carbohydrate() && ! conf.residue( starting_position ).is_branch_point() ) {
		TR << "Delete to residue is not carbohydrate and not a branch point.  Nothing to be done." << std::endl;
		return std::make_pair( list_of_residues, tips);
	}


	Size parent_residue;
	if ( conf.residue( starting_position).is_carbohydrate() ) {
		parent_residue = find_seqpos_of_saccharides_parent_residue( conf.residue( starting_position ) );
	} else {
		parent_residue = starting_position - 1;
	}

	fill_upstream_children_res_and_tips( conf, starting_position, parent_residue, children_residues, list_of_residues, tips );

	//TR << "Children: " << utility::to_string( children_residues) << std::endl;
	get_branching_residues( conf, starting_position, children_residues, list_of_residues, tips );

	return std::make_pair( list_of_residues, tips );
}

///@brief Get the carbohydrate residue connecting the protein branch point.
core::Size
get_glycan_connecting_protein_branch_point(conformation::Conformation const & conf, core::Size const protein_branch_point_resnum){

	debug_assert(conf.residue(protein_branch_point_resnum).is_branch_point());

	core::Size parent_residue = protein_branch_point_resnum - 1;

	Size connections = conf.residue( protein_branch_point_resnum ).n_possible_residue_connections();
	for ( core::Size connection = 1; connection <= connections; ++connection ) {

		Size connecting_res = conf.residue( protein_branch_point_resnum ).connected_residue_at_resconn( connection );

		if ( connecting_res != 0 && connecting_res != parent_residue && conf.residue( connecting_res ).is_carbohydrate() ) {
			return connecting_res;
		}
	}
	return 0;

}

///@brief Get the particular resnum from a glycan position, givin the protein branch point.
/// The glycan_position is numbered 1 -> length of glycan. This is useful for easily identifying a particular glycan position.
///
core::Size
get_resnum_from_glycan_position(conformation::Conformation const & conf, core::Size const glycan_one, core::Size const glycan_position){
	using namespace utility;

	if ( ( glycan_position == 1 )| ( glycan_one == 0 ) ) {
		return glycan_one;
	} else {
		Size glycan_length = get_glycan_tree_size(conf, glycan_one);
		if ( glycan_position <= glycan_length ) {
			core::Size glycan_residue  = glycan_one + glycan_position - 1;
			return glycan_residue;
		} else {
			return 0;
		}
	}

}

core::Size
get_glycan_position_from_resnum(conformation::Conformation const & conf, core::Size const glycan_one, core::Size const glycan_residue ){

	if ( glycan_one == glycan_residue ) {
		return 1;
	} else {
		Size glycan_length = get_glycan_tree_size(conf, glycan_one);
		core::Size glycan_position = glycan_residue + 1 - glycan_one;

		if ( glycan_position > glycan_length ) {
			return 0;
		} else {
			return glycan_position;
		}
	}



}

}  // carbohydrates
}  // conformation
}  // core
