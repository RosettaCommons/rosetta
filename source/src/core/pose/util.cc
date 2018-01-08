// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/util.cc
/// @brief  Pose class utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov, Rhiju Das, Steven Lewis, Vikram K. Mulligan


// Unit header
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>

// Package headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/PositionConservedResiduesStore.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/carbohydrates/util.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/chemical/util.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/Exceptions.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/conformation/Residue.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <basic/datacache/CacheableStringMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.string.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/string_constants.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

// C/C++ headers
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <boost/functional/hash.hpp>

namespace core {
namespace pose {

static basic::Tracer TR( "core.pose.util" );

void swap_transform(Size jump_num, kinematics::RT const & xform, Pose & pose) {
	debug_assert(jump_num <= pose.num_jump());

	kinematics::FoldTree const & tree = pose.fold_tree();

	int upstream = tree.upstream_jump_residue(jump_num);
	int downstream = tree.downstream_jump_residue(jump_num);
	conformation::Residue const & upstream_res = pose.residue(upstream);
	conformation::Residue const & downstream_res = pose.residue(downstream);

	if ( upstream_res.natoms() < 3 || downstream_res.natoms() < 3 ) {
		TR.Warning << "Insufficient number of atoms for stub creation on one or more"
			<< " jump residues-- " << upstream << ", " << downstream
			<< std::endl;
		return;
	}

	id::StubID upstream_stub(
		id::AtomID(0, upstream_res.atom_index("N")),
		id::AtomID(1, upstream_res.atom_index("CA")),
		id::AtomID(2, upstream_res.atom_index("C")));

	id::StubID downstream_stub(
		id::AtomID(0, downstream_res.atom_index("N")),
		id::AtomID(1, downstream_res.atom_index("CA")),
		id::AtomID(2, downstream_res.atom_index("C")));

	debug_assert(upstream_stub.valid());
	debug_assert(downstream_stub.valid());

	pose.conformation().set_stub_transform(
		upstream_stub,
		downstream_stub,
		xform);
}

bool is_position_conserved_residue(Pose const & pose, core::Size residue) {
	using basic::datacache::BasicDataCache;
	using core::pose::datacache::PositionConservedResiduesStore;
	using core::pose::datacache::PositionConservedResiduesStoreCOP;

	debug_assert(residue > 0);
	debug_assert(residue <= pose.size());

	BasicDataCache const & cache = pose.data();
	if ( !cache.has(core::pose::datacache::CacheableDataType::POSITION_CONSERVED_RESIDUES) ) {
		return false;
	}

	PositionConservedResiduesStoreCOP store =
		utility::pointer::static_pointer_cast<PositionConservedResiduesStore const>(
		cache.get_const_ptr(core::pose::datacache::CacheableDataType::POSITION_CONSERVED_RESIDUES));

	return store->is_conserved(residue);
}


/// @brief  Create a chemical edge between two Residues
/// @param  <start_resnum>: the Residue index of the start of the chemical edge
/// @param  <end_resnum>: the Residue index of the end of the chemical edge
/// @param  <ft>: the FoldTree being modified
/// @return true if a chemical connection was found and the passed FoldTree was changed
/// @note   This function is used by set_reasonable_foldtree() for the primary purpose of rebuilding default FoldTrees
/// in cases of branching or chemical conjugation.
bool
create_chemical_edge(
	core::uint start_resnum,
	core::uint end_resnum,
	core::pose::Pose const & pose,
	core::kinematics::FoldTree & ft )
{
	using namespace core::chemical;
	using namespace conformation;

	Residue const & start_residue( pose.residue( start_resnum ) );
	Residue const & end_residue( pose.residue( end_resnum ) );
	// Don't attempt to make a connection which isn't actually in the connection map
	if ( ! end_residue.is_bonded( start_resnum ) ) {
		return false;
	}
	debug_assert( start_residue.is_bonded( end_resnum ) );

	utility::vector1< Size > const & end_connections( end_residue.connections_to_residue(start_resnum) );
	debug_assert( end_connections.size() >= 1 );

	// We arbitrarily pick the first connection - if we get double connections for a pair of residues,
	// the FoldTree only *slightly* cares about which one we use.
	core::Size end_connid( end_connections[1] );
	debug_assert( end_residue.connected_residue_at_resconn(end_connid) == start_resnum );

	// Get the partner via the connid in case there's a double connection
	core::Size start_connid( end_residue.residue_connection_conn_id(end_connid) );
	debug_assert( start_residue.connected_residue_at_resconn(start_connid) == end_resnum );

	uint start_atm_index = start_residue.residue_connection( start_connid ).atomno();
	uint end_atm_index = end_residue.residue_connection( end_connid ).atomno();
	std::string const & start_atm_name = start_residue.atom_name( start_atm_index );
	std::string const & end_atm_name = end_residue.atom_name( end_atm_index );

	if ( TR.Trace.visible() ) {
		TR.Trace << "Adding  EDGE " <<
			start_resnum << ' ' << end_resnum << ' ' << start_atm_name << ' ' << end_atm_name <<
			" to the new FoldTree." << std::endl;
	}
	ft.add_edge( start_resnum, end_resnum, start_atm_name, end_atm_name );
	return true;
}


/// @details If all ligand residues and polymer branches have been appended by a jump, this method creates a new
/// FoldTree without jumps through ligands, using CHEMICAL EDGEs instead.
void
set_reasonable_fold_tree( pose::Pose & pose )
{
	// An empty pose doesn't have jumps through ligands.
	// (Will encounter a SegFault otherwise)
	if ( pose.size() == 0 ) return;

	using namespace std;
	using namespace core::chemical;
	using namespace core::kinematics;

	if ( TR.Debug.visible() ) {
		TR.Debug << "original fold tree: " << pose.fold_tree() << endl;
	}

	FoldTree const & origft( pose.fold_tree() );
	FoldTree newft;

	// Copy the original fold tree edges to the new fold tree, possibly replacing
	// some jumps to ligand residues with ligand-ligand chemical bonds.
	// As a result, jumps must be renumbered.
	uint last_jump_id( 0 );
	bool prevent_forward_edge( false );  // Use to prevent a custom reverse edge from being overwritten.
	for ( auto i( origft.begin() ), i_end( origft.end() ); i != i_end; ++i ) {
		Edge e( *i );
		if ( TR.Trace.visible() ) {
			TR.Trace << "checking if " << e << " is reasonable for this Pose..." << endl;
		}

		if ( prevent_forward_edge ) {
			if ( TR.Trace.visible() ) {
				TR.Trace << "skipping edge, as we've already added the reverse version" << endl;
			}
			prevent_forward_edge = false;
			continue;
		}

		if ( e.is_jump() ) {
			uint const jump_stop( e.stop() );  // the residue position at the end of the jump
			uint const jump_start( e.start() );

			// Does this jump represent a direct (chemical) connection between the two residues?
			if ( pose.residue( jump_stop ).is_bonded( jump_start ) ) {
				if ( TR.Trace.visible() ) {
					TR.Trace << e << " was initially set as a Jump that matches a chemical connection." << endl;
				}
				create_chemical_edge( jump_start, jump_stop, pose, newft );
				continue;  // The Edge has already been added to newft by the function above.
			}

			// Do we have a lower connect to a residue already in the new FoldTree? If so, connect it by a chemical edge rather than a jump edge.
			if ( pose.residue( jump_stop ).has_lower_connect() && newft.residue_is_in_fold_tree( pose.residue( jump_stop ).connected_residue_at_lower() ) ) {
				if ( TR.Trace.visible() ) {
					TR.Trace << e << " was initially set as a Jump to a residue with a lower connection to a previously defined residue." << std::endl;
				}
				create_chemical_edge( pose.residue( jump_stop ).connected_residue_at_lower(), jump_stop, pose, newft );
				continue;
			}

			// Is this jump to the lower end of a polymeric edge,
			// which has it's upper end chemically connected to something that's already in the fold-tree?
			// i.e., is this a C-terminal connection?

			// RM: Assuming that it's the *next* edge which is the associated polymeric edge is slightly off,
			// But the current way we're skipping adding the polymeric edge means it's hard to do it otherwise.
			Edge const next_e( *(i + 1) );
			if ( next_e.is_polymer() && jump_stop == uint( next_e.start() ) ) {
				if ( TR.Trace.visible() ) {
					TR.Trace << e << " was initially set as a Jump to end of a polymeric connection." << endl;
				}

				uint const polymer_stop( next_e.stop() );  // the residue position at the end of the next Edge

				if ( pose.residue( polymer_stop ).has_upper_connect() && newft.residue_is_in_fold_tree( pose.residue( polymer_stop ).connected_residue_at_upper() ) ) {
					if ( TR.Trace.visible() ) {
						TR.Trace << "Adding chemical edge from " << jump_start << " to " << polymer_stop << " to the new FoldTree, along with reverse polymeric edge from " << polymer_stop << " to " << jump_stop << endl;
					}
					create_chemical_edge(  pose.residue( polymer_stop ).connected_residue_at_upper(), polymer_stop, pose, newft );
					// The CHEMICAL Edge has already been added to newft by the function above.
					// Now, we must add a reverse-direction Edge and skip making the usual forward-direction Edge.
					newft.add_edge( polymer_stop, jump_stop, Edge::PEPTIDE);
					prevent_forward_edge = true;
					continue;  // Skip the add_edge() method below; we are done here.
				}
			}

			// RM: Should we be doing something if it's a jump to a polymeric edge that connects withing the *middle* of the polymeric edge?
			// (That is, an H-shaped topology?)

			// If we made it here, it's just a normal inter-chain jump or we couldn't find a chemical bond to make. Increment the label.
			e.label() = ++last_jump_id;
		} else {
			if ( TR.Trace.visible() ) {
				TR.Trace << e << " was not a Jump." << endl;
			}
		}
		if ( TR.Trace.visible() ) {
			TR.Trace << "Adding " << e << " to the new FoldTree." << endl;
		}
		newft.add_edge( e );
	}  // next Edge

	runtime_assert( newft.size() > 0 || pose.size() == 0 );  // A valid fold tree must have >= 1 edges.

	if ( TR.Debug.visible() ) {
		TR.Debug << "new fold tree " << newft << endl;
	}

	pose.fold_tree( newft );
}

/// @brief Return the appropritate ResidueType for the virtual residue for the
/// "mode" (fullatom, centroid ...) the pose is in.
core::chemical::ResidueTypeCOP
virtual_type_for_pose(core::pose::Pose const & pose) {
	core::chemical::ResidueTypeCOP type( get_restype_for_pose( pose, "VRT" ) );
	if ( ! type ) {
		utility_exit_with_message("Cannot find VRT residue.");
	}
	return type;
}

/// @details Adds a virtual residue to a pose as the root. Jump is to the
/// residue closest to <xyz>. If the pose is already rooted on a VRT res,
/// do nothing.
void addVirtualResAsRoot(const numeric::xyzVector<core::Real>& xyz, core::pose::Pose& pose) {
	int nres = pose.size();

	// if already rooted on virtual residue, return
	if ( pose.residue( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) {
		TR.Warning << "addVirtualResAsRoot() called but pose is already rooted on a VRT residue ... continuing." << std::endl;
		return;
	}

	// return if the pose is empty (otherwise will segfault)
	if ( nres == 0 ) {
		TR.Warning << "addVirtualResAsRoot() called with empty pose!" << std::endl;
		return;
	}

	// check for terminal ligands
	int last_peptide_res = nres;
	while ( !pose.residue( last_peptide_res ).is_polymer() ) last_peptide_res--;

	// try to avoid putting the vrt too close to termini
	int i_min = 1;

#ifdef WIN32
	int r_start = static_cast< int > ( std::floor(   static_cast< double > (last_peptide_res) /3. ) );
	int r_end   = static_cast< int > ( std::ceil ( 2.* static_cast< double > (last_peptide_res)/3. ) );
#else
	auto r_start = static_cast< int > ( std::floor(   last_peptide_res/3 ) );
	auto r_end   = static_cast< int > ( std::ceil ( 2*last_peptide_res/3 ) );
#endif

	// If less than three total residues, reset starting residue
	if ( r_start == 0 ) {
		r_start = 1;
	}

	core::Length d_min = std::numeric_limits<core::Length>::infinity(), this_d;
	for ( int i=r_start; i<=r_end; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );

		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
		if ( !rsd.is_protein() ) continue;

		conformation::Atom const & atom( rsd.atom("CA") );
		this_d = (atom.xyz() - xyz).length();
		if ( this_d < d_min ) {
			d_min = this_d;
			i_min = i;
		}
	}
	core::chemical::ResidueTypeCOP rsd_type( virtual_type_for_pose(pose) );
	core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *rsd_type ) );

	// move to <xyz>
	for ( Size j=1; j<= new_res->natoms(); ++j ) {
		new_res->atom(j).xyz( new_res->atom(j).xyz()+xyz );
	}

	pose.append_residue_by_jump( *new_res , i_min );

	// update PDBinfo
	if ( pose.pdb_info() ) {
		pose.pdb_info()->chain( pose.size(), 'z' );
		pose.pdb_info()->number( pose.size(), 1 );
		pose.pdb_info()->obsolete( false );
	}

	// make the virt atom the root
	kinematics::FoldTree newF( pose.fold_tree() );
	newF.reorder( nres+1 );
	TR.Debug << "addVirtualResAsRoot() setting new fold tree to " << newF << std::endl;
	TR.Debug << "   i_min = " << i_min << "   d_min = " << d_min << std::endl;
	pose.fold_tree( newF );
}

void remove_virtual_residues(core::pose::Pose & pose) {
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue_type(i).name() == "VRT" ) {
			pose.conformation().delete_residue_slow(i);
		}
	}
}

/// @detail Get center of mass of a pose.
numeric::xyzVector< core::Real >
get_center_of_mass( core::pose::Pose const & pose ){
	numeric::xyzVector< core::Real > massSum( 0.0 );
	Size const & nres = pose.size();
	Size nAtms=0;
	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
		for ( Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );
			massSum += atom.xyz();
			nAtms++;
		}
	}
	massSum /= nAtms;
	return massSum;
}

/// @detail Find residue closest to center-of-mass
void addVirtualResAsRoot( core::pose::Pose & pose ) {
	int nres = pose.size();
	// return if the pose is empty (otherwise will segfault)
	if ( nres == 0 ) {
		TR.Warning << "addVirtualResAsRoot() called with empty pose!" << std::endl;
		return;
	}
	numeric::xyzVector< core::Real > massSum = get_center_of_mass( pose );
	return addVirtualResAsRoot(massSum, pose);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Repair pdbinfo of inserted residues that may have blank chain and zero
/// seqpos. Assumes insertions only occur _after_ a residue.
/// @details That assumption is a bad simplifying assumption, but it's all I need
/// for now.
void fix_pdbinfo_damaged_by_insertion(
	core::pose::Pose & pose
) {
	if ( !pose.pdb_info() ) return;

	PDBInfo & pdbinfo = *pose.pdb_info();

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( pdbinfo.number( ii ) != 0 || ( pdbinfo.chain( ii ) != ' ' && pdbinfo.chain( ii ) != '^' ) ) {
			// This residue is fine.
			continue;
		}

		if ( ii == 1 ) {
			utility_exit_with_message( "Error: somehow, an insertion is present at the first residue!" );
		}

		int const num_prev = pdbinfo.number( ii - 1 );
		char const chn_prev = pdbinfo.chain( ii - 1 );
		if ( ii == pose.size() ) {
			pdbinfo.number( ii, num_prev + 1 );
			pdbinfo.chain( ii, chn_prev );
			continue;
		}

		Size next_sensible_seqpos = ii + 1;
		while ( pdbinfo.number( next_sensible_seqpos ) == 0 && ( pdbinfo.chain( next_sensible_seqpos ) == ' ' || pdbinfo.chain( next_sensible_seqpos ) == '^' ) ) ++next_sensible_seqpos;

		// If this is the last residue of a chain, also easy.
		if ( pdbinfo.chain( next_sensible_seqpos ) != chn_prev ) {
			pdbinfo.number( ii, chn_prev + 1 );
			pdbinfo.chain( ii, chn_prev );
			continue;
		}

		int const num_next = pdbinfo.number( next_sensible_seqpos );
		if ( num_next - num_prev >= 2 ) {
			// Squeeze in in the middle
			pdbinfo.number( ii, num_prev + 1 );
			pdbinfo.chain( ii, chn_prev );
		} else if ( num_next - num_prev == 1 ) {
			// Insertion after num_prev, incrementing icode
			pdbinfo.number( ii, num_prev );
			pdbinfo.chain( ii, chn_prev );
			char const icode_prev = pdbinfo.icode( ii - 1 );
			if ( icode_prev == ' ' ) {
				pdbinfo.icode( ii, 'A' );
			} else {
				pdbinfo.icode( ii, char( pdbinfo.icode( ii - 1 ) + 1 ) );
			}
		} else if ( num_next == num_prev ) {
			// Inserting in the middle of an icode sequence...
			pdbinfo.number( ii, num_prev );
			pdbinfo.chain( ii, chn_prev );
			char const icode_prev = pdbinfo.icode( ii - 1 );
			char const icode_next = pdbinfo.icode( next_sensible_seqpos );
			if ( icode_prev == ' ' ) {
				pdbinfo.icode( ii, 'A' );
			} else {
				if ( int( icode_next - icode_prev ) >= 2 ) {
					// Easy, there's icode space to go around
					pdbinfo.icode( ii, char( icode_prev + 1 ) );
				} else {
					// Uh oh, gotta increment subsequent icodes.
					Size jj = ii;
					char icode = icode_prev + 1;
					while ( pdbinfo.number( jj ) == num_prev ) {
						pdbinfo.icode( jj, char( icode ) );
						jj++; icode++;
					}
				}
			}
		}
	}
}

/// @brief checks if the pose geometry is ideal
/// @param[in] pose The Pose to check.
/// @return true if all pose positions have ideal bond lengths and angles
///  up to some very small epsilon
bool is_ideal_pose(
	core::pose::Pose const & pose
) {
	bool is_ideal=true;
	for ( core::Size i=1 ; i < pose.size(); i++ ) {//leaving out last residue which always returns non-ideal for some reason
		if ( !is_ideal_position(i, pose) ) {
			is_ideal=false;
			break;
		}
	}
	return is_ideal;
}

/// @brief checks if the pose geometry is ideal in position seqpos
/// @param[in] pose The Pose to check.
/// @return true if position seqpos has ideal bond lengths and angles
///  up to some very small epsilon
bool is_ideal_position(
	core::Size seqpos,
	core::pose::Pose const & pose
) {
	return conformation::is_ideal_position( seqpos, pose.conformation() );
}

/// @brief this function removes all residues from the pose which are not protein residues.
/// @details This removal includes, but is not limited to, metals, DNA, RNA, and ligands.
/// It will NOT remove ligands which are canonical residues (for example, if a protein binds an alanine monomer,
/// the monomer will be untouched).
void remove_nonprotein_residues( core::pose::Pose & pose )
{
	core::Size i(1);
	while ( i <= pose.size() ) {
		if ( !pose.residue_type(i).is_protein() ) pose.conformation().delete_residue_slow(i);
		else ++i;
	}
}

/// @brief this function removes all residues with both UPPER and LOWER terminus types.
/// This is intended for removing ligands that are canonical residues.
void remove_ligand_canonical_residues( core::pose::Pose & pose )
{
	if ( pose.size() == 1 ) { //if we have only one residue, it cannot be removed, and this is going to crash
		throw CREATE_EXCEPTION(utility::excn::Exception, "Pose utility remove_ligand_canonical_residues: I have received a pose with only one residue but cannot delete the last residue of the pose.");
	}

	core::Size i(1);
	while ( i <= pose.size() ) {
		if ( pose.residue_type(i).is_upper_terminus() && pose.residue_type(i).is_lower_terminus() ) {
			pose.conformation().delete_residue_slow(i);
		} else ++i;
	}
}


/// @details this function compares the 3-d coordinates of two poses.
/// Along the way it is forced to check for certain other (in)equalities to prevent vector overrruns, etc.
/// These include: pose length, ResidueType, and # atoms in residue.
/// Inequalities other than 3-d coordinates result in a warning message (you shouldn't have been comparing those poses!)
/// This is NOT a complete equality operator for a pose, but I think it does a good job with the coordinates.
/// Note that it performs floating-point coordinate comparisons by floor(X*10^n_dec_places) -
/// this may cause failures if your pose is a billion angstroms from 0,0,0.
/// This comparison is preferred to an epsilon comparison std::abs( a.x - b.x ) < epsilon because it can run into
/// situations where a == b and b == c, but a != c (thanks to APL for pointing this out).
/// The last argument, n_dec_places, is the number of decimal places of precision when comparing.
bool compare_atom_coordinates(core::pose::Pose const & lhs, core::pose::Pose const & rhs, core::Size const n_dec_places){

	//number of decimal places of precision - 3 (1000) is equivalent to PDB precision.
	core::Real const n_dec(pow(static_cast< Real > (10), static_cast< int > (n_dec_places))); //this is a Real to prevent premature rounding

	//first compare pose sizes; prerequisite to iterating through length
	core::Size const lhssize(lhs.size()), rhssize(rhs.size());

	if ( lhssize != rhssize ) {
		TR.Warning << "poses of different length in compare_atom_coordinates; doomed to fail!" << std::endl;
		return false;
	}

	//now iterate through residues and make comparisons
	for ( core::Size i(1); i<=lhssize; ++i ) {
		//check equality of residue types
		core::chemical::ResidueType const & lhstype(lhs.residue_type(i)), & rhstype(rhs.residue_type(i));
		if ( lhstype.name() != rhstype.name() ) { //string matching is sufficient because ResidueType objects have unique names
			TR.Warning << "nonmatching ResidueTypes at " << i << " in compare_atom_coordinates" << std::endl;
			return false;
		}

		//get atoms vectors to compare
		core::conformation::Atoms const & lhsatoms(lhs.residue(i).atoms()), rhsatoms(rhs.residue(i).atoms());
		core::Size const lhsatmsize(lhsatoms.size()), rhsatmsize(rhsatoms.size());
		if ( lhsatmsize != rhsatmsize ) { //check vector length equality
			TR.Warning << "nonmatching numbers of atoms at residue " << i << " in compare_atom_coordinates" << std::endl;
			TR.Warning << "How did we even get here?  ResidueType comparison should have failed!" << std::endl;
			return false;
		}

		//iterate through atoms vector
		for ( core::Size atm(1); atm <= lhsatmsize; ++atm ) {
			if (  (std::floor(lhsatoms[atm].xyz().x()*n_dec) != std::floor(rhsatoms[atm].xyz().x()*n_dec))
					|| (std::floor(lhsatoms[atm].xyz().y()*n_dec) != std::floor(rhsatoms[atm].xyz().y()*n_dec))
					|| (std::floor(lhsatoms[atm].xyz().z()*n_dec) != std::floor(rhsatoms[atm].xyz().z()*n_dec)) ) {
				return false; //no warning messages, this is the "expected" failure
			}
		}//iterate through atoms vector
	}//over all residues

	return true; //whoo! we made it!
}//compare_atom_coordinates


bool
compare_binary_protein_silent_struct(
	Pose const & lhs,
	Pose const & rhs
) {
	core::io::silent::SilentFileOptions opts; // initialized from the command line
	core::io::silent::BinarySilentStruct lhs_silent_struct(opts, lhs, "" );
	std::stringstream lhs_str;
	lhs_silent_struct.print_conformation(lhs_str);

	core::io::silent::BinarySilentStruct rhs_silent_struct(opts, rhs, "" );
	std::stringstream rhs_str;
	rhs_silent_struct.print_conformation(rhs_str);

	return lhs_str.str() == rhs_str.str();
}

////////////////////////////////////////////////////////////////////////////////////////////////

id::NamedAtomID
atom_id_to_named_atom_id(
	core::id::AtomID const & atom_id,
	Pose const & pose
) {
	conformation::Residue const& rsd( pose.residue( atom_id.rsd() ) );
	return core::id::NamedAtomID( rsd.atom_name( atom_id.atomno() ), atom_id.rsd() );
}

/// @details returns an AtomID corresponding to your NamedAtomID
/// check for a valid AtomID after this.
/// following conditions return invalid ID :
/// rsd > size
/// atom not present in residue ( e.g., no CB in GLY )
id::AtomID
named_atom_id_to_atom_id(
	core::id::NamedAtomID const & named_atom_id,
	Pose const & pose,
	bool raise_exception /*default true*/
) {
	using namespace core::id;
	// work out the stubID
	if ( named_atom_id.valid() ) {
		if ( named_atom_id.rsd() <= pose.size() ) { //if in range, ... otherwise leave atomno_ on 0 --> valid() == false
			chemical::ResidueType const& rt ( pose.residue_type ( named_atom_id.rsd() ) );
			if ( rt.has( named_atom_id.atom() ) ) {
				return AtomID( rt.atom_index( named_atom_id.atom() ), named_atom_id.rsd() );
			} else {
				// tr.Error << "can't find atom " << named_atom_id.atom() << " in residue "
				//   << rt.name() << ", residue has " << rt.nheavyatoms() << " heavy atoms." << std::endl;
				//  tr.Error << "atom names are: " << std::endl;
				//rt.show_all_atom_names( tr.Error );
				if ( raise_exception ) throw CREATE_EXCEPTION(id::EXCN_AtomNotFound, named_atom_id );
				return id::GLOBAL_BOGUS_ATOM_ID;
			}
		} else {
			// tr.Error << "can't find residue " << named_atom_id.rsd()
			//  << " in pose (pose.size() = ) "
			//  << pose.size() << std::endl;
			if ( raise_exception ) throw CREATE_EXCEPTION(id::EXCN_AtomNotFound, named_atom_id );
			return id::GLOBAL_BOGUS_ATOM_ID;
		}
	} else {
		if ( raise_exception ) throw CREATE_EXCEPTION(id::EXCN_AtomNotFound, named_atom_id );
		return id::GLOBAL_BOGUS_ATOM_ID;
	}
}


id::NamedStubID
stub_id_to_named_stub_id( id::StubID const& stub_id, core::pose::Pose const& pose ) {
	using namespace core::id;
	if ( stub_id.center().valid() ) {
		return NamedStubID(
			NamedAtomID( atom_id_to_named_atom_id( stub_id.center() , pose ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 1 ), pose ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 2 ), pose ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 3 ), pose ) )
		);
	} else {
		return NamedStubID( // TODO does this make sense? if input stub is bad, what should be done?
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 1 ), pose ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 2 ), pose ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 3 ), pose ) )
		);
	}
}

id::StubID
named_stub_id_to_stub_id( id::NamedStubID const& named_stub_id, core::pose::Pose const& pose ) {
	using namespace core::id;
	if ( named_stub_id.center().valid() ) {
		return StubID(
			AtomID( named_atom_id_to_atom_id( named_stub_id.center() , pose ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 1 ), pose ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 2 ), pose ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 3 ), pose ) )
		);
	} else {
		return StubID( // TODO does this make sense? if input stub is bad, what should be done?
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 1 ), pose ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 1 ), pose ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 2 ), pose ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 3 ), pose ) )
		);
	}
}


///////////////////////////////////////////////////////////////

core::Real energy_from_pose(
	core::pose::Pose const & pose, core::scoring::ScoreType const & sc_type
) {
	return pose.energies().total_energies()[ sc_type ];
}

core::Real energy_from_pose(
	core::pose::Pose const & pose, std::string const & sc_type
) {
	return pose.energies().total_energies()[ core::scoring::score_type_from_name( sc_type ) ];
}

core::Real total_energy_from_pose( core::pose::Pose const & pose ) {
	Real total_energy = pose.energies().total_energy();
	// when initiated from a silent struct, total energy lives in "extra data".
	if ( total_energy == Real( 0.0 ) ) getPoseExtraScore( pose, "score", total_energy );
	return total_energy;
}

// criterion for sorting.
bool
sort_pose_by_score( core::pose::PoseOP const & pose1, core::pose::PoseOP const & pose2 ) {
	return ( total_energy_from_pose( *pose1 ) < total_energy_from_pose( *pose2 ) );
}

void
transfer_phi_psi( const core::pose::Pose& srcpose, core::pose::Pose& tgtpose, core::Size ir, core::Size jr )
{
	core::Size tgtlength = tgtpose.size();
	core::Size srclength = srcpose.size();
	for ( core::Size ires = ir; ires <= std::min( srclength, std::min( tgtlength, jr)) ; ires++ ) {
		if ( !srcpose.residue_type( ires ).is_protein() || !tgtpose.residue_type( ires ).is_protein() ) continue;
		tgtpose.set_phi(   ires, srcpose.phi(ires) );
		tgtpose.set_psi(   ires, srcpose.psi(ires) );
		tgtpose.set_omega( ires, srcpose.omega(ires) );
	}
}

void
transfer_phi_psi( const core::pose::Pose& srcpose, core::pose::Pose& tgtpose ){
	transfer_phi_psi( srcpose, tgtpose, 1, std::min( tgtpose.size(), srcpose.size() ) );
}

//fpd transfer the RB geometry from one pose to another
void
transfer_jumps( const core::pose::Pose& srcpose, core::pose::Pose& tgtpose )
{
	core::kinematics::FoldTree f_tgt=tgtpose.fold_tree(), f_src=srcpose.fold_tree();

	// project fold tree
	core::Size srcjmps = srcpose.num_jump();
	for ( core::Size ijmp = 1; ijmp <= srcjmps ; ijmp++ ) {
		core::kinematics::Edge srcedge_i = srcpose.fold_tree().jump_edge(ijmp);
		core::Size jjmp = tgtpose.fold_tree().jump_nr( srcedge_i.start(), srcedge_i.stop() );
		if ( jjmp != 0 ) {
			tgtpose.set_jump( jjmp, srcpose.jump(ijmp) );
		} else {
			TR.Debug << "In transfer_jumps() unable to map jump: " << srcedge_i.start() << " ,  " << srcedge_i.stop() << std::endl;
		}
	}
}


/// helper function for residue replacement/residuetype switching
/// these functions should probably move to pose/util.cc
/// @note  Will call new_rsd->fill_missing_atoms if the new residue has atoms
/// that the old one doesn't
void
replace_pose_residue_copying_existing_coordinates(
	pose::Pose & pose,
	Size const seqpos,
	core::chemical::ResidueType const & new_rsd_type
) {
	core::conformation::Residue const & old_rsd( pose.residue( seqpos ) );
	core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( new_rsd_type ) );
	conformation::copy_residue_coordinates_and_rebuild_missing_atoms( old_rsd, *new_rsd, pose.conformation() );
	pose.replace_residue( seqpos, *new_rsd, false );

}

core::chemical::ResidueTypeCOP
get_restype_for_pose(core::pose::Pose const & pose, std::string const & name) {
	core::chemical::TypeSetMode mode( pose.conformation().residue_typeset_mode() );
	if ( mode == core::chemical::MIXED_t ) {
		// Fall back to fullatom if the mode of the pose is perfectly balanced
		mode = core::chemical::FULL_ATOM_t;
	}
	return get_restype_for_pose(pose, name, mode);
}

core::chemical::ResidueTypeCOP
get_restype_for_pose(core::pose::Pose const & pose, std::string const & name, core::chemical::TypeSetMode mode) {
	core::chemical::ResidueTypeSetCOP residue_set( pose.residue_type_set_for_pose( mode ) );
	debug_assert( residue_set );
	core::chemical::ResidueTypeCOP rsd_type( residue_set->name_mapOP(name) );
	if ( ! rsd_type ) {
		TR.Error << "Can't find residue type '" << name << "' in type set of mode " << mode << std::endl;
	}
	return rsd_type;
}

/// @brief returns a Distance
core::Real
pose_max_nbr_radius( Pose const & pose )
{
	core::Real maxrad( 0.0 );
	for (  core::Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( pose.residue( ii ).nbr_radius() > maxrad ) maxrad = pose.residue_type(ii).nbr_radius();
	}
	return maxrad;
}

///////////////////////////////////////////////////////////////////////////////
void
setup_dof_to_torsion_map( pose::Pose const & pose, id::DOF_ID_Map< id::TorsionID > & dof_map )
{
	using namespace id;

	// Set DOF mask size and initialize to an invalid TorsionID
	core::pose::initialize_dof_id_map( dof_map, pose, id::TorsionID::BOGUS_TORSION_ID() );

	// Torsion angles
	Size const n_res( pose.size() );
	for ( Size i = 1; i <= n_res; ++i ) {
		// PHIL note the Residue-based helper functions you need for this
		// PHIL also note the pose.conformation() interface

		conformation::Residue const & rsd( pose.residue( i ) );

		// first the backbone torsion angles
		Size const n_bb_torsions( rsd.mainchain_atoms().size() );
		for ( uint j( 1 ); j <= n_bb_torsions; ++j ) {
			TorsionID const tor_id( i, BB, j );
			DOF_ID const & id( pose.conformation().dof_id_from_torsion_id( tor_id ) );
			if ( id.valid() ) {
				dof_map[ id ] = tor_id;
			}
		}  // j=1, n_bb_torsions

		// then the side chain torsions
		Size const n_chi_torsions( rsd.nchi() );
		for ( uint j( 1 ); j <= n_chi_torsions; ++j ) {
			TorsionID const tor_id( i, CHI, j );
			DOF_ID const & id( pose.conformation().dof_id_from_torsion_id( tor_id ) );
			if ( id.valid() ) {
				if ( dof_map[ id  ] != TorsionID::BOGUS_TORSION_ID() ) {
					// This chi angle was already assigned to the DoF map as a BB TorsionType. Do not change it.
					continue;
				}
				dof_map[ id ] = tor_id;
			}
		}  // j=1, n_chi_torsions

		// next, the internal ring torsions
		Size const n_nu_torsions( rsd.n_nus() );
		for ( uint j( 1 ); j <= n_nu_torsions; ++j ) {
			TorsionID const tor_id( i, NU, j );
			DOF_ID const & id( pose.conformation().dof_id_from_torsion_id( tor_id ) );
			if ( id.valid() ) {
				dof_map[ id ] = tor_id;
			}
		}  // j=1, n_nu_torsions

		// finally, any branch connection torsions
		Size const n_branch_torsions( rsd.n_non_polymeric_residue_connections() );
		for ( uint j( 1 ); j <= n_branch_torsions; ++j ) {
			TorsionID const tor_id( i, BRANCH, j );
			DOF_ID const & id( pose.conformation().dof_id_from_torsion_id( tor_id ) );
			if ( id.valid() ) {
				dof_map[ id ] = tor_id;
			}
		}  // j=1, n_branch_torsions

	}  // i=1, n_res

	for ( Size i = 1; i <= pose.num_jump(); ++i ) {
		for ( int j=1; j<= 6; ++j ) {
			id::TorsionID const tor_id(i,id::JUMP,j);
			id::DOF_ID const & id( pose.conformation().dof_id_from_torsion_id( tor_id ) );
			dof_map[ id ] =tor_id;
		}
	} // i=1.num_jump
}


///////////////////////////////////////////////////////////////////////////////
// Convert from allow-bb/allow-chi MoveMap to simple DOF_ID boolean mask needed by the minimizer
void
setup_dof_mask_from_move_map( kinematics::MoveMap const & mm, pose::Pose const & pose, id::DOF_ID_Mask & dof_mask )
{
	using namespace id;

	bool const PHI_default( mm.get( PHI ) );

	// Set DOF mask size and initialize to false.
	initialize_dof_id_map( dof_mask, pose, false );

	// DOF_Type defaults
	// could/should do this with a loop over types?
	// currently ignoring rb types, set these individually by jump number
	dof_mask.set( PHI  , mm.get( PHI   ) );
	dof_mask.set( THETA, mm.get( THETA ) );
	dof_mask.set( D    , mm.get( D     ) );

	// Torsion angles
	Size const n_res( pose.size() );
	for ( Size i = 1; i <= n_res; ++i ) {
		// PHIL note the Residue-based helper functions you need for this
		// PHIL also note the pose.conformation() interface

		conformation::Residue const & rsd( pose.residue( i ) );

		// first the main-chain torsion angles
		Size const n_mainchain_torsions( rsd.mainchain_atoms().size() );

		// Note: the BB setting in the MoveMap is a partial misnomer.  It really refers to main-chain torsions.  Other
		// backbone atoms may necessarily move (such as the carbonyl O of a peptide backbone carbonyl), but only
		// torsions that are part of the main chain in the AtomTree will be sampled.
		// In many (most?) cases, a ResidueType with a ring will have overlap between its internal ring torsions
		// and its main-chain torsions as defined by the AtomTree.  In such cases, while the ring atoms are part of the
		// backbone, only some of the atoms are part of the main chain.  Changing these main-chain torsions without
		// also adapting the other backbone ring atoms to form a new ring conformer will "rip open" the ring.  Thus,
		// for the purposes of the BB setting of the MoveMap, only those main-chain torsions that are not also internal
		// ring torsions will be moved.  For example, in the case of carbohydrates, phi, psi, and omega are to be
		// considered main-chain torsions, as the bonds corresponding to them are exocyclic, but the other main-chain
		// torsions should be ignored by the MoveMap BB setting.  If one wants to modify those other mainchain torsions,
		// they should be treated as nu angles and the setting for NU in the MoveMap should be used instead. ~Labonte

		for ( uint j( 1 ); j<= n_mainchain_torsions; ++j ) {
			TorsionID const torsion( i, BB, j );
			bool mm_setting;
			if ( is_mainchain_torsion_also_ring_torsion( rsd.type(), j ) ) {
				mm_setting = false;  // Do not move "backbone" torsions that are part of a ring.
			} else {
				mm_setting = mm.get( torsion );
			}
			if ( mm_setting == PHI_default ) { continue; }
			DOF_ID const & id( pose.conformation().dof_id_from_torsion_id( torsion ) );
			if ( id.valid() ) {  // If not valid, it's probably just a terminal/chainbreak torsion.
				dof_mask[ id ] = mm_setting;
			}
		} // j=1, n_bb_torsions

		// then the side chain torsions
		Size const n_chi_torsions( rsd.nchi() );
		for ( uint j = 1; j <= n_chi_torsions; ++j ) {
			TorsionID const torsion( i, CHI, j );
			bool const mm_setting( mm.get( torsion ) );
			if ( mm_setting == PHI_default ) { continue; }
			DOF_ID const & id( pose.conformation().dof_id_from_torsion_id( torsion ) );
			if ( id.valid() ) {
				if ( rsd.is_carbohydrate() && rsd.is_virtual( id.atomno() ) ) {
					dof_mask[ id ] = false;  // Do not move virtual atoms that are simply placeholders.
				} else {
					dof_mask[ id ] = mm_setting;
				}
			} else {
				TR.Warning << "Unable to find atom_tree atom for this " <<
					"Rosetta chi angle: residue " << i << " CHI " << j << std::endl;
			}
		} // j=1, n_chi_torsions

		// next, the internal ring torsions
		Size const n_nu_torsions( rsd.n_nus() );
		for ( uint j( 1 ); j <= n_nu_torsions; ++j ) {
			TorsionID const torsion( i, NU, j );
			bool const mm_setting = mm.get( torsion );
			if ( mm_setting == PHI_default ) { continue; }
			DOF_ID const & id( pose.conformation().dof_id_from_torsion_id( torsion ) );
			if ( id.valid() ) {
				dof_mask[ id ] = mm_setting;
			} else {
				TR.Warning << "Unable to find atom_tree atom for this " <<
					"Rosetta nu angle: residue " << i << " NU " << j << std::endl;
			}
		} // j=1, n_nu_torsions

		// finally, any branch connection torsions
		Size const n_branch_torsions( rsd.n_non_polymeric_residue_connections() );
		for ( uint j( 1 ); j <= n_branch_torsions; ++j ) {
			// Note: If one has multiple incoming connections on some crazy residue, she or he is still going to get a
			// warning triggered below, but it won't hurt anything, and I can't think of another way to silence such a
			// warning at the moment. ~Labonte
			TorsionID const torsion( i, BRANCH, j );
			bool const mm_setting = mm.get( torsion );
			if ( mm_setting == PHI_default ) { continue; }
			DOF_ID const & id( pose.conformation().dof_id_from_torsion_id( torsion ) );
			if ( id.valid() ) {
				dof_mask[ id ] = mm_setting;
			} else {
				TR.Warning << "Unable to find atom_tree atom for this " <<
					"Rosetta branch connection angle: residue " << i << " BRANCH " << j << std::endl;
			}
		} // j=1, n_branch_torsions
	} // i=1, n_res

	// Jumps.
	for ( Size i=1; i<= pose.num_jump(); ++i ) {
		if ( mm.get_jump(i) ) {
			for ( int j=1; j<= 6; ++j ) {
				DOF_ID const & id( pose.conformation().dof_id_from_torsion_id(TorsionID(i, JUMP, j)));
				dof_mask[ id ] = true;
			}
		}
	} // i=1, num_jump

	/////////////////////////
	// DOFs set individually
	for ( auto it=mm.dof_id_begin(), it_end=mm.dof_id_end();
			it != it_end; ++it ) {
		dof_mask[ it->first ] = it->second;
	}
} // setup_dof_mask_from_move_map

core::Size num_heavy_atoms(
	core::Size begin,
	core::Size const end,
	core::pose::Pose const & pose
) {
	core::Size total_heavy_atoms = 0;
	for ( ; begin <= end; ++begin ) {
		total_heavy_atoms += pose.residue(begin).nheavyatoms();
	}
	TR.Debug << "# of heavy atoms: "<< total_heavy_atoms << std::endl;
	return total_heavy_atoms;
}

core::Size num_atoms(
	core::Size begin,
	core::Size const end,
	core::pose::Pose const & pose
) {
	core::Size total_atoms = 0;
	for ( ; begin <= end; ++begin ) {
		total_atoms += pose.residue(begin).natoms();
	}
	TR.Debug << "# of heavy atoms: "<< total_atoms << std::endl;
	return total_atoms;
}

core::Size num_hbond_acceptors(
	core::Size begin,
	core::Size const end,
	core::pose::Pose const & pose
) {
	core::Size total_hbond_acceptors = 0;
	for ( ; begin <= end; ++begin ) {
		total_hbond_acceptors += pose.residue(begin).n_hbond_acceptors();
	}
	TR.Debug << "# of heavy atoms: "<< total_hbond_acceptors << std::endl;
	return total_hbond_acceptors;
}

core::Size num_hbond_donors(
	core::Size begin,
	core::Size const end,
	core::pose::Pose const & pose
) {
	core::Size total_hbond_donors = 0;
	for ( ; begin <= end; ++begin ) {
		total_hbond_donors += pose.residue(begin).n_hbond_donors();
	}
	TR.Debug << "# of heavy atoms: "<< total_hbond_donors << std::endl;
	return total_hbond_donors;
}

core::Size num_chi_angles(
	core::Size begin,
	core::Size const end,
	core::pose::Pose const & pose
) {
	core::Size total_chi_angles = 0;
	for ( ; begin <= end; ++begin ) {
		total_chi_angles += pose.residue(begin).nchi();
	}
	return total_chi_angles;
}

core::Real
mass(
	core::Size begin,
	core::Size const end,
	core::pose::Pose const & pose
){
	core::Real mass = 0;
	for ( ; begin <= end; ++begin ) {
		mass += pose.residue(begin).type().mass();
	}
	return mass;
}

void
initialize_disulfide_bonds(
	Pose & pose
) {
	// disulfides
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	// Fix disulfides if a file is given
	if ( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].user() ) {
		core::io::raw_data::DisulfideFile ds_file( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ]() );
		utility::vector1< std::pair<Size,Size> > disulfides;
		ds_file.disulfides(disulfides, pose);
		pose.conformation().fix_disulfides( disulfides );
	} else {
		if ( option[ in::detect_disulf ].user() ?
				option[ in::detect_disulf ]() : // detect_disulf true
				pose.is_fullatom() // detect_disulf default but fa pose
				) {
			pose.conformation().detect_disulfides();
		}
	}
}

// This version is only suitable for being called from build_pdb_from_pose_as_is1
// i.e. where the StructFileRep exists!
void
initialize_disulfide_bonds(
	Pose & pose,
	io::StructFileRep const & sfr
) {
	// disulfides
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	// Fix disulfides if a file is given
	if ( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].user() ) {
		core::io::raw_data::DisulfideFile ds_file( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ]() );
		utility::vector1< std::pair<Size,Size> > disulfides;
		ds_file.disulfides(disulfides, pose);
		pose.conformation().fix_disulfides( disulfides );
	} else if ( option[ in::detect_disulf ].user() ?
			option[ in::detect_disulf ]() : // detect_disulf true
			pose.is_fullatom() // detect_disulf default but fa pose
			) {
		//utility::vector1< std::pair< Size, Size > > disulfs;
		utility::vector1< Size > disulf_one;
		utility::vector1< Size > disulf_two;

		// Prepare a list of pose-numbered disulfides!
		for ( auto const & ssbond : sfr.ssbond_map() ) {

			// For now we really hope the vector1 is just a single element!
			if ( ssbond.second.size() != 1 ) {
				// We can salvage if it's double-entry: just take the first.
				// The length isn't actually used anyway, and so it doesn't
				// affect what conformation is preferred.

				// Are they all the same?
				std::string id1 = ssbond.second[1].resID2;
				bool identical = true;
				for ( Size i = 2; i <= ssbond.second.size(); ++i ) {
					if ( ssbond.second[ i ].resID2 != id1 ) {
						identical = false; break;
					}
				}
				if ( !identical ) {
					TR.Error << "SSBond records list multiple nonredundant disulfides for this residue!" << std::endl;
					utility_exit_with_message("Error with SSBond record.");
				}
			}

			Size seqpos_one = pose.pdb_info()->pdb2pose( ssbond.second[1].chainID1, ssbond.second[1].resSeq1, ssbond.second[1].iCode1 );
			Size seqpos_two = pose.pdb_info()->pdb2pose( ssbond.second[1].chainID2, ssbond.second[1].resSeq2, ssbond.second[1].iCode2 );

			if ( seqpos_one != 0 && seqpos_two != 0 ) {
				disulf_one.push_back( seqpos_one );
				disulf_two.push_back( seqpos_two );
			}
		}

		pose.conformation().detect_disulfides( disulf_one, disulf_two );
	}
}

core::id::SequenceMapping sequence_map_from_pdbinfo( Pose const & first, Pose const & second ) {
	core::id::SequenceMapping retval(first.size(), second.size());
	core::pose::PDBInfoCOP first_pdbinfo = first.pdb_info();
	core::pose::PDBInfoCOP second_pdbinfo = second.pdb_info();

	if ( first_pdbinfo && !first_pdbinfo->obsolete() && second_pdbinfo && !second_pdbinfo->obsolete() ) {
		for ( core::Size ii(1); ii<= first.size(); ++ii ) {
			// pdb2pose returns 0 for "not found" - 0 is also used for "not found" for SequenceMapping.
			retval[ii] = second_pdbinfo->pdb2pose( first_pdbinfo->chain(ii), first_pdbinfo->number(ii), first_pdbinfo->icode(ii) );
		}
	} else {
		TR << "One or both poses do not have usable PDBInfo, using sequence alignment instead." << std::endl;
		retval = core::sequence::map_seq1_seq2(
			core::sequence::SequenceOP( new core::sequence::Sequence(first ) ),
			core::sequence::SequenceOP( new core::sequence::Sequence(second) )
		);
	}

	return retval;
}

core::Size canonical_residue_count(core::pose::Pose const & pose)
{
	return std::count_if( pose.begin(), pose.end(),
		[&]( conformation::Residue const & resi ) {
			return resi.aa() <= core::chemical::num_canonical_aas;
		} );
}

core::Size noncanonical_residue_count(core::pose::Pose const & pose)
{
	return std::count_if( pose.begin(), pose.end(),
		[&]( conformation::Residue const & resi ) {
			return resi.aa() > core::chemical::num_canonical_aas;
		} );
}

core::Size canonical_atom_count(core::pose::Pose const & pose)
{
	return std::accumulate( pose.begin(), pose.end(), 0,
		[&]( Size const posum, conformation::Residue const & resi ) {
			if ( resi.aa() <= core::chemical::num_canonical_aas ) {
				return posum + resi.natoms();
			}
			return posum;
		} );
}

core::Size noncanonical_atom_count(core::pose::Pose const & pose)
{
	return std::accumulate( pose.begin(), pose.end(), 0,
		[&]( Size const posum, conformation::Residue const & resi ) {
			if ( resi.aa() > core::chemical::num_canonical_aas ) {
				return posum + resi.natoms();
			}
			return posum;
		} );
}

core::Size noncanonical_chi_count(core::pose::Pose const & pose)
{
	return std::accumulate( pose.begin(), pose.end(), 0,
		[&]( Size const posum, conformation::Residue const & resi ) {
			if ( resi.aa() > core::chemical::num_canonical_aas ) {
				return posum + resi.nchi();
			}
			return posum;
		} );
}

/// @brief Number of protein residues in the pose
/// @details No virtuals, membrane residues or embedding residues counted
Size nres_protein( pose::Pose const & pose ) {
	return std::count_if( pose.begin(), pose.end(),
		[&]( conformation::Residue const & resi ) {
			return resi.is_protein();
		} );
}// nres_protein


numeric::xyzVector< Real >
center_of_mass(
	pose::Pose const & pose,
	utility::vector1< bool > const & residues
) {
	using namespace numeric;
	using core::conformation::Residue;
	debug_assert( pose.size() == residues.size() );

	utility::vector1< xyzVector< Real > > coords;

	for ( Size i = residues.l(); i <= residues.u(); ++i ) {
		if ( residues[ i ] ) {
			Residue const & rsd( pose.residue( i ) );
			coords.push_back( rsd.is_protein() ? rsd.atom( "CA" ).xyz() : rsd.nbr_atom_xyz() );
		}
	}

	if ( ! coords.size() ) { utility_exit_with_message( "Cannot compute center of mass of zero residues!" ); }

	return center_of_mass( coords );
}


// anonymous function to assist in converting from start, stop to a vector1< bool >
utility::vector1< bool >
generate_vector_from_bounds(
	pose::Pose const & pose,
	int const start,
	int const stop
) {
	utility::vector1< bool > residues( pose.size(), false );
	debug_assert( (Size) stop <= residues.size() );
	debug_assert( stop > start && start > 0 );

	for ( int i = start; i <= stop; ++i ) { residues[ i ] = true; }
	return residues;
}

////////////////////////////////////////////////////////////////////////////////////
///
/// @brief calculates the center of mass of a pose
/// @details
/// the start and stop positions (or residues) within the pose are used to
/// find the starting and finishing locations
///
/// @author Monica Berrondo June 14 2007
///
/////////////////////////////////////////////////////////////////////////////////
numeric::xyzVector< core::Real>
center_of_mass(
	pose::Pose const & pose,
	int const start,
	int const stop
) {
	return center_of_mass( pose, generate_vector_from_bounds( pose, start, stop ) );
}

core::Vector
all_atom_center(
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & residues
) {

	utility::vector1< core::Vector > coords;

	for ( Size ii(1); ii <= residues.size(); ++ii ) {
		core::conformation::Residue const & rsd( pose.residue( residues[ii] ) );
		for ( Size jj(1); jj <= rsd.natoms(); ++ jj ) {
			coords.push_back( rsd.xyz( jj ) );
		}
	}

	if ( ! coords.size() ) { utility_exit_with_message( "Cannot compute center of mass for no atoms!" ); }

	return numeric::center_of_mass( coords );
}

int
residue_center_of_mass(
	pose::Pose const & pose,
	utility::vector1< bool > residues
) {
	debug_assert( pose.size() == residues.size() );

	if ( ! ( residues.has( true ) && pose.size() )  ) {
		utility_exit_with_message( "Cannot compute center of mass of zero residues!" );
	}

	Vector center = center_of_mass(pose, residues );
	return core::pose::return_nearest_residue( pose, residues, center );
}

////////////////////////////////////////////////////////////////////////////////////
///
/// @brief calculates the center of mass of a pose
/// @details
///    the start and stop positions (or residues) within the pose are used to
///    find the starting and finishing locations
///
/// @author Monica Berrondo June 14 2007
///
/////////////////////////////////////////////////////////////////////////////////
int
residue_center_of_mass(
	pose::Pose const & pose,
	int const start,
	int const stop
) {
	Vector center = center_of_mass(pose, start, stop );
	return core::pose::return_nearest_residue( pose, start, stop, center );
}

int
return_nearest_residue(
	pose::Pose const & pose,
	utility::vector1< bool > const & residues,
	Vector center
) {
	using core::conformation::Residue;

	debug_assert( pose.size() == residues.size() );

	if ( ! ( residues.has( true ) && pose.size() )  ) {
		utility_exit_with_message( "Cannot find nearest residue in empty selection!" );
	}

	Real min_dist = std::numeric_limits< Real >::infinity();
	int res = 0;
	for ( Size i = residues.l(); i <= residues.u(); ++i ) {
		if ( ! residues[ i ] ) { continue; }

		Residue const & rsd( pose.residue( i ) );
		Vector const ca_pos = rsd.is_protein() ? rsd.atom( "CA" ).xyz() : rsd.nbr_atom_xyz();

		Real const dist( ca_pos.distance_squared( center ) );
		if ( dist < min_dist ) {
			res = i;
			min_dist = dist;
		}
	}
	return res;
}

////////////////////////////////////////////////////////////////////////////////////
///
/// @brief finds the residue nearest some position passed in (normally a
///  center of mass)
/// @details
///     the start and stop positions (or residues) within the pose are used to
///     find the starting and finishing locations
///
/// @author Monica Berrondo June 14 2007
///
/////////////////////////////////////////////////////////////////////////////////
int
return_nearest_residue(
	pose::Pose const & pose,
	int const begin,
	int const end,
	Vector center
) {
	return return_nearest_residue( pose, generate_vector_from_bounds( pose, begin, end ), center );
}

// silly conversion from std::map< AtomID, AtomID> to rosetta's silly AtomID_Map class.
id::AtomID_Map< id::AtomID >
convert_from_std_map( std::map< id::AtomID, id::AtomID > const & atom_map,
	core::pose::Pose const & pose ){
	id::AtomID_Map< id::AtomID > atom_ID_map;
	initialize_atomid_map( atom_ID_map, pose, id::AtomID::BOGUS_ATOM_ID() );
	for ( auto const & it : atom_map ) {
		atom_ID_map.set( it.first, it.second );
	}
	return atom_ID_map;
}

/// @brief Create std::map from PDBPoseMap in pose (JKLeman)
std::map< std::string, core::Size >
get_pdb2pose_numbering_as_stdmap( core::pose::Pose const & pose ) {

	using namespace utility;

	// initialize empty map
	std::map< std::string, core::Size > pdb2pose_map;

	// go through protein, get chain, resn, insertion code
	// AMW: this only makes sense for poses where the protein residues come first
	// This makes sense for membranes, but not e.g. DNA/RNA/Protein complexes.
	for ( core::Size i = 1; i <= nres_protein( pose ); ++i ) {

		// get chain, PDB resnumber, and insertion code for that residue
		std::string const chain = to_string( pose.pdb_info()->chain( i ) );
		std::string const resn = to_string( pose.pdb_info()->number( i ) );
		std::string icode = to_string( pose.pdb_info()->icode( i ) );

		// replace empty insertion code with dot (spanfile has dots)
		if ( icode == " " ) {
			icode = ".";
		}

		// concatenate them into string
		std::string pdb_num = chain + resn + icode;

		// add it to local map
		pdb2pose_map[ pdb_num ] = i;

	}
	return pdb2pose_map;
}


/// @brief Create a chemical bond from lower to upper residue across CUTPOINT_LOWER/CUTPOINT_UPPER.
/// @details This will prevent steric repulsion.
/// @param[in] pose The pose to modify.
/// @param[in] cutpoint_res The index of the CUTPOINT_LOWER residue.
/// @param[in] next_res_in The index of the CUTPOINT_UPPER residue.  If not provided, or if set to 0, this defaults
/// to the cutpoint_res + 1 residue.  Must be specified for cyclic geometry.
void
declare_cutpoint_chemical_bond(
	core::pose::Pose & pose,
	Size const cutpoint_res,
	Size const next_res_in /* = 0 */
) {
	using namespace core::conformation;
	// AMW TODO: ensure this refers to the residue connected at UPPER instead of cutpoint_res + 1
	Size const next_res = ( next_res_in == 0 ) ? ( cutpoint_res + 1 ) : next_res_in; // user might specify a different "next_res" to cyclize.

	// Need to clear out any chemical bonds that might have been previously tied to upper/lower of these residues.
	Residue const & lower_rsd( pose.conformation().residue( cutpoint_res ) );

	for ( Size k = 1; k <= lower_rsd.connect_map_size(); k++ ) {
		if ( lower_rsd.has_upper_connect() && lower_rsd.residue_connect_atom_index( k ) != lower_rsd.upper_connect_atom() ) continue;
		// If there isn't an upper connect, one of the few valid options is ZO3'
		if ( !lower_rsd.has_upper_connect() && ( ! lower_rsd.has_variant_type( core::chemical::FIVEPRIME_CAP ) || lower_rsd.atom_name( lower_rsd.residue_connect_atom_index( k ) ) != "ZO3'" ) ) continue;
		Size upper( lower_rsd.connected_residue_at_resconn( k ) );
		if ( upper == 0 ) continue;
		Residue const & upper_rsd( pose.conformation().residue( upper ) ); // upper residue.
		Size const m = lower_rsd.residue_connection_conn_id( k );
		runtime_assert( upper_rsd.residue_connect_atom_index( m ) == upper_rsd.lower_connect_atom() );
		runtime_assert( upper_rsd.connected_residue_at_resconn( m ) == cutpoint_res );
		//upper_rsd.mark_connect_incomplete( m );
		//lower_rsd.mark_connect_incomplete( k );
		pose.conformation().sever_chemical_bond( cutpoint_res, k, upper, m );
	}

	// and code up analogous loop for lower/upper.
	Residue const & upper_rsd( pose.conformation().residue( next_res ) );
	for ( Size k = 1; k <= upper_rsd.connect_map_size(); k++ ) {
		if ( upper_rsd.has_lower_connect() && upper_rsd.residue_connect_atom_index( k ) != upper_rsd.lower_connect_atom() ) continue;
		Size lower( upper_rsd.connected_residue_at_resconn( k ) );
		if ( lower == 0 ) continue;
		Residue const & lower_rsd( pose.conformation().residue( lower ) ); // lower residue.
		Size const m = upper_rsd.residue_connection_conn_id( k );
		runtime_assert( lower_rsd.residue_connect_atom_index( m ) == lower_rsd.upper_connect_atom() );
		runtime_assert( lower_rsd.connected_residue_at_resconn( m ) == next_res );
		//lower_rsd.mark_connect_incomplete( m );
		//upper_rsd.mark_connect_incomplete( k );
		pose.conformation().sever_chemical_bond( next_res, k, lower, m );
	}

	// Two options. Might be ZO3'...
	if ( !lower_rsd.has_upper_connect() ) {
		debug_assert( pose.residue( cutpoint_res ).has( "ZO3'" ) );
		// Right now, assume ZO3'
		pose.conformation().declare_chemical_bond(
			cutpoint_res,
			"ZO3'",
			next_res,
			pose.residue( next_res ).atom_name( pose.residue( next_res ).lower_connect_atom() ) );
	} else {
		pose.conformation().declare_chemical_bond(
			cutpoint_res,
			pose.residue( cutpoint_res ).atom_name( pose.residue( cutpoint_res ).upper_connect_atom() ),
			next_res,
			pose.residue( next_res ).atom_name( pose.residue( next_res ).lower_connect_atom() ) );
	}
}

/// @brief Given a pose and a position that may or may not be CUTPOINT_UPPER or CUTPOINT_LOWER, determine whether this
/// position has either of these variant types, and if it does, determine whether it's connected to anything.  If it is,
/// update the C-OVL1-OVL2 bond lengths and bond angle (for CUTPOINT_LOWER) or OVU1-N bond length (for CUTPOINT_UPPER) to
/// match any potentially non-ideal geometry in the residue to which it's bonded.
/// @details Requires a little bit of special-casing for gamma-amino acids.  Throws an exception if the residue to which
/// a CUTPOINT_LOWER is bonded does not have an "N" and a "CA" or "C4".  Safe to call repeatedly, or if cutpoint variant
/// types are absent; in these cases, the function does nothing.
/// @note By default, this function calls itself again once on residues to which this residue is connected, to update their
/// geometry.  Set recurse=false to disable this.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
update_cutpoint_virtual_atoms_if_connected(
	core::pose::Pose & pose,
	core::Size const cutpoint_res,
	bool recurse /*= true*/
) {
	core::conformation::update_cutpoint_virtual_atoms_if_connected( pose.conformation(), cutpoint_res, recurse );
}

void
get_constraints_from_link_records( core::pose::Pose & pose, io::StructFileRep const & sfr )
{
	using namespace scoring::func;
	using namespace scoring::constraints;

	/*HarmonicFuncOP amide_harm_func    ( new HarmonicFunc( 1.34, 0.05 ) );
	HarmonicFuncOP thioester_harm_func( new HarmonicFunc( 1.83, 0.1 ) );
	*/
	CircularHarmonicFuncOP ang_func( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi_2_over_3(), 0.02 ) );
	CircularHarmonicFuncOP ang90_func( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi_over_2(), 0.02 ) );
	CircularHarmonicFuncOP dih_func( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi(), 0.02 ) );

	for ( auto const & it : sfr.link_map() ) {
		for ( auto const & record : it.second ) {
			TR << "|"<<record.chainID1 << "| |" << record.resSeq1 << "|" << std::endl;
			TR << "|"<<record.chainID2 << "| |" << record.resSeq2 << "|" << std::endl;

			Size id1 = pose.pdb_info()->pdb2pose( record.chainID1, record.resSeq1 );
			Size id2 = pose.pdb_info()->pdb2pose( record.chainID2, record.resSeq2 );
			conformation::Residue const & NUC = pose.residue( id1 );
			conformation::Residue const & ELEC = pose.residue( id2 );

			id::AtomID aidNUC( NUC.atom_index( record.name1 ), id1 );
			id::AtomID aidC( ELEC.atom_index( record.name2 ), id2 );

			scoring::func::HarmonicFuncOP harm_func( new scoring::func::HarmonicFunc( record.length, 0.05 ) );

			scoring::constraints::ConstraintCOP atompair(
				new scoring::constraints::AtomPairConstraint( aidNUC, aidC, harm_func ) );
			pose.add_constraint( atompair );
			TR << "Adding harmonic constraint between residue " << id1 << " atom " << record.name1;
			TR << "and residue " << id2 << " atom " << record.name2 << " with length " << record.length << std::endl;

			// Cover cyclization case first
			if ( record.name1 == " N  " && record.name2 == " C  " ) {

				id::AtomID aidCA( ELEC.atom_index( "CA" ), id2 );
				id::AtomID aidO( ELEC.atom_index( "O" ), id2 );
				id::AtomID aidCA2( NUC.atom_index( "CA" ), id1 );

				ConstraintCOP ang( new AngleConstraint( aidNUC, aidC, aidCA, ang_func ) );
				pose.add_constraint( ang );
				ConstraintCOP ang2( new AngleConstraint( aidCA2, aidNUC, aidC, ang_func ) );
				pose.add_constraint( ang2 );
				ConstraintCOP dih( new DihedralConstraint( aidNUC, aidC, aidO, aidCA, dih_func ) );
				pose.add_constraint( dih );

				// constrain omega unless N terminal residue is PRO or peptoid
				// I think we should have a residue property that describes PRO-type residue types
				if ( NUC.type().name3() == "PRO" || NUC.type().is_peptoid() ) {
					// let god sort it out
					// PERHAPS: constraints that mirror the MM torsions actually guiding this case!
				} else {
					id::AtomID aidH( NUC.atom_index( "H" ), id1 );

					ConstraintCOP omg( new DihedralConstraint( aidCA2, aidNUC, aidC, aidCA, dih_func ) );
					pose.add_constraint( omg );
					ConstraintCOP omgH( new DihedralConstraint( aidO, aidC, aidNUC, aidH, dih_func ) );
					pose.add_constraint( omgH );
					ConstraintCOP omgimp( new DihedralConstraint( aidC, aidNUC, aidH, aidCA2, dih_func ) );
					pose.add_constraint( omgimp );
				}

				TR << "Adding harmonic constraints to the angle formed by atoms N, C, O ( 120 ) and ";
				TR << "the improper torsion N, C, O, CA (180) and the dihedral CA, N, C, CA ( 180 ) " <<std::endl;
				continue;
			}

			if ( record.name1 == " C  " && record.name2 == " N  " ) {
				// swapped; nuc is actually elec and vice versa

				id::AtomID aidCA2( ELEC.atom_index( "CA" ), id2 );
				id::AtomID aidO( NUC.atom_index( "O" ), id1 );
				id::AtomID aidCA( NUC.atom_index( "CA" ), id1 );

				ConstraintCOP ang( new AngleConstraint( aidC, aidNUC, aidCA, ang_func ) );
				pose.add_constraint( ang );
				ConstraintCOP ang2( new AngleConstraint( aidCA2, aidC, aidNUC, ang_func ) );
				pose.add_constraint( ang2 );
				ConstraintCOP dih( new DihedralConstraint( aidC, aidNUC, aidO, aidCA, dih_func ) );
				pose.add_constraint( dih );

				// constrain omega unless N terminal residue is PRO or peptoid
				// I think we should have a residue property that describes PRO-type residue types
				if ( ELEC.type().name3() == "PRO" || ELEC.type().is_peptoid() ) {
					// let god sort it out
					// PERHAPS: constraints that mirror the MM torsions actually guiding this case!
				} else {
					id::AtomID aidH( ELEC.atom_index( "H" ), id1 );

					ConstraintCOP omg( new DihedralConstraint( aidCA2, aidC, aidNUC, aidCA, dih_func ) );
					pose.add_constraint( omg );
					ConstraintCOP omgH( new DihedralConstraint( aidO, aidNUC, aidC, aidH, dih_func ) );
					pose.add_constraint( omgH );
					ConstraintCOP omgimp( new DihedralConstraint( aidNUC, aidC, aidH, aidCA2, dih_func ) );
					pose.add_constraint( omgimp );
				}
				continue;
			}

			// Now let's add constraints based on residue identities and atom names. For example, let's cover
			if ( record.name2 == " CZ " && record.resName2 == "VDP" ) {
				// thiol-ene conjugation to acryl residue
				// (don't check name1 because we don't care SG/SD/SG1)
				// someday we will be fancy and check vs type.get_disulfide_atom_name()
				id::AtomID aidCB( NUC.atom_index( "CB" ), id2 );
				id::AtomID aidCE2( ELEC.atom_index( "CE2" ), id2 );
				id::AtomID aidO( ELEC.atom_index( "O" ), id2 );
				ConstraintCOP ang( new AngleConstraint( aidCB, aidNUC, aidC, ang90_func ) );
				pose.add_constraint( ang );
				//scoring::constraints::ConstraintCOP dih(
				//  new scoring::constraints::DihedralConstraint( aidCB, aidNUC, aidC, aidCE2, dih_func ) );
				//pose.add_constraint( dih );

				TR << "Assuming thiol-ene, adding harmonic constraints to the angle formed by CB, SG, CZ ( 90 )" << std::endl;// and ";
				//TR << "the dihedral CB, SG, CZ, CE2 ( 180 ) " << std::endl;
				continue;
			}

			if ( record.name2 == " C  " ) {
				// the C-terminal conjugation case:
				id::AtomID aidCA( ELEC.atom_index( "CA" ), id2 );
				id::AtomID aidO( ELEC.atom_index( "O" ), id2 );
				ConstraintCOP ang( new AngleConstraint( aidNUC, aidC, aidO, ang_func ) );
				pose.add_constraint( ang );
				ConstraintCOP dih( new scoring::constraints::DihedralConstraint( aidNUC, aidC, aidO, aidCA, dih_func ) );
				pose.add_constraint( dih );
				continue;
			}

			if ( record.name2 == " CD " ) {
				// The sidechain conjugation to glx case
				id::AtomID aidCA( ELEC.atom_index( "CG" ), id2 );
				id::AtomID aidO( ELEC.atom_index( "OE1" ), id2 );
				ConstraintCOP ang( new AngleConstraint( aidNUC, aidC, aidO, ang_func ) );
				pose.add_constraint( ang );
				ConstraintCOP dih( new DihedralConstraint( aidNUC, aidC, aidO, aidCA, dih_func ) );
				pose.add_constraint( dih );
				continue;
			}

			if ( record.name2 == " CG " ) {
				// The sidechain conjugation to asx case
				id::AtomID aidCA( ELEC.atom_index( "CB" ), id2 );
				id::AtomID aidO( ELEC.atom_index( "OD1" ), id2 );
				ConstraintCOP ang( new AngleConstraint( aidNUC, aidC, aidO, ang_func ) );
				pose.add_constraint( ang );
				ConstraintCOP dih( new DihedralConstraint( aidNUC, aidC, aidO, aidCA, dih_func ) );
				pose.add_constraint( dih );
				if ( record.name1 == " NE " ) {
					// ornithine
					id::AtomID aidH( NUC.atom_index( "1HE" ), id2 );
					id::AtomID aidCG( NUC.atom_index( "CG" ), id2 );

					pose.add_constraint( ConstraintCOP( new DihedralConstraint( aidCA, aidO, aidC, aidNUC, dih_func ) ) );
					pose.add_constraint( ConstraintCOP( new DihedralConstraint( aidO, aidC, aidNUC, aidH, dih_func ) ) );
					pose.add_constraint( ConstraintCOP( new DihedralConstraint( aidC, aidNUC, aidH, aidCG, dih_func ) ) );
					pose.add_constraint( ConstraintCOP( new AngleConstraint( aidH, aidNUC, aidC, ang_func ) ) );
					pose.add_constraint( ConstraintCOP( new AngleConstraint( aidNUC, aidC, aidO, ang_func ) ) );

				} else if ( record.name1 == " NZ " ) {
					// lysine
					id::AtomID aidH( NUC.atom_index( "1HZ" ), id2 );
					id::AtomID aidCD( NUC.atom_index( "CD" ), id2 );

					pose.add_constraint( ConstraintCOP( new DihedralConstraint( aidCA, aidO, aidC, aidNUC, dih_func ) ) );
					pose.add_constraint( ConstraintCOP( new DihedralConstraint( aidO, aidC, aidNUC, aidH, dih_func ) ) );
					pose.add_constraint( ConstraintCOP( new DihedralConstraint( aidC, aidNUC, aidH, aidCD, dih_func ) ) );
					pose.add_constraint( ConstraintCOP( new AngleConstraint( aidH, aidNUC, aidC, ang_func ) ) );
					pose.add_constraint( ConstraintCOP( new AngleConstraint( aidNUC, aidC, aidO, ang_func ) ) );

				}
				continue;
			}
		}
	}
}

/// @brief Convert PDB numbering to pose numbering. Must exist somewhere else, but I couldn't find it. -- rhiju
utility::vector1< Size > pdb_to_pose( pose::Pose const & pose, utility::vector1< int > const & pdb_res ){
	utility::vector1< Size > pose_res;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( pdb_res.has_value( pose.pdb_info()->number( n ) ) ) pose_res.push_back( n );
	}
	runtime_assert( pose_res.size() == pdb_res.size() );
	return pose_res;
}

/// @brief Convert PDB numbering to pose numbering. Must exist somewhere else, but I couldn't find it. -- rhiju
utility::vector1< Size > pdb_to_pose( pose::Pose const & pose, std::tuple< utility::vector1< int >, utility::vector1<char>, utility::vector1<std::string> > const & pdb_res ){
	utility::vector1< Size > pose_res;
	runtime_assert( std::get<0>( pdb_res ).size() == std::get<1>( pdb_res ).size() );
	for ( Size n = 1; n <= std::get<0>( pdb_res ).size(); n++ ) {
		pose_res.push_back( pdb_to_pose( pose, std::get<0>( pdb_res )[ n ], std::get<1>( pdb_res )[ n ] ) );
	}
	return pose_res;
}

/// @brief Convert PDB numbering to pose numbering. Must exist somewhere else, but I couldn't find it. -- rhiju
Size
pdb_to_pose( pose::Pose const & pose, int const res_num, char const chain /* = ' ' */ ) {
	bool found_match( false );
	Size res_num_in_full_numbering( 0 );
	for ( Size n = 1; n <= pose.size() ; n++ ) {
		if ( res_num != pose.pdb_info()->number( n ) ) continue;
		if ( chain != ' ' && pose.pdb_info()->chain( n ) != ' ' &&
				pose.pdb_info()->chain( n ) != chain ) continue;
		res_num_in_full_numbering = n;
		if ( found_match != 0 ) utility_exit_with_message( "Ambiguous match to pdb number " + ObjexxFCL::string_of(res_num ) + " and chain " + ObjexxFCL::string_of( chain ) );
		found_match = true;
	}
	if ( !found_match ) utility_exit_with_message( "Could not match residue number  " + ObjexxFCL::string_of( res_num ) + " and chain " + ObjexxFCL::string_of(chain) );
	return res_num_in_full_numbering;
}

/// @brief Convert pose numbering to pdb numbering. Must exist somewhere else, but I couldn't find it. -- rhiju
utility::vector1< Size > pose_to_pdb( pose::Pose const & pose, utility::vector1< Size > const & pose_res ){
	utility::vector1< Size > pdb_res;
	for ( Size i = 1; i <= pose_res.size(); i++ ) {
		runtime_assert( pose_res[ i ] <= pose.size() );
		pdb_res.push_back( pose.pdb_info()->number( pose_res[ i ] ) );
	}
	return pdb_res;
}

// Is the query atom in this pose residue axial or equatorial to the given ring or neither?
/// @details This function calculates an average plane and determines whether the coordinates of a given atom are
/// axial or equatorial to it (or neither).
/// @param   <pose>:       The Pose containing the Residue in question.
/// @param   <seqpos>:     The sequence position in the Pose of the Residue containing the atoms in question.
/// @param   <query_atom>: The index of the atom in question.
/// @param   <ring_atoms>: A list of indices for the atoms of a monocyclic ring system in sequence.
/// @return  An AxEqDesignation enum type value: AXIAL, EQUATORIAL, or NEITHER
/// @author  Labonte <JwLabonte@jhu.edu>
chemical::rings::AxEqDesignation
is_atom_axial_or_equatorial_to_ring(
	Pose const & pose,
	uint seqpos,
	uint query_atom,
	utility::vector1< uint > const & ring_atoms )
{
	return conformation::is_atom_axial_or_equatorial_to_ring( pose.residue( seqpos ), query_atom, ring_atoms );
}

// Is the query atom in this pose axial or equatorial to the given ring or neither?
/// @details This function calculates an average plane and determines whether the coordinates of a given atom are
/// axial or equatorial to it (or neither).
/// @param   <pose>:       The Pose containing the atoms in question.
/// @param   <query_atom>: The AtomID of the atom in question.
/// @param   <ring_atoms>: A list of AtomIDs for the atoms of a monocyclic ring system in sequence.
/// @return  An AxEqDesignation enum type value: AXIAL, EQUATORIAL, or NEITHER
/// @author  Labonte <JwLabonte@jhu.edu>
chemical::rings::AxEqDesignation
is_atom_axial_or_equatorial_to_ring(
	Pose const & pose,
	id::AtomID const & query_atom,
	utility::vector1< id::AtomID > const & ring_atoms )
{
	uint const seqpos( query_atom.rsd() );
	Size const ring_size( ring_atoms.size() );
	utility::vector1< uint > ring_atom_indices( ring_size );
	for ( uint i( 1 ); i <= ring_size; ++i ) {
		if ( ring_atoms[ i ].rsd() != seqpos ) {
			TR.Warning << "Queried atom must be in the same residue as the ring atoms." << std::endl;
			return chemical::rings::NEITHER;
		}
		ring_atom_indices[ i ] = ring_atoms[ i ].atomno();
	}
	return conformation::is_atom_axial_or_equatorial_to_ring(
		pose.residue( seqpos ), query_atom.atomno(), ring_atom_indices );
}

// Is the query atom in this pose residue axial or equatorial or neither?
/// @details This function calculates an average plane and determines whether the coordinates of a given atom are
/// axial or equatorial to it (or neither).  The ring is requested from the Residue.
/// @param   <pose>:       The Pose containing the Residue in question.
/// @param   <seqpos>:     The sequence position in the Pose of the Residue containing the atoms in question.
/// @param   <query_atom>: The index of the atom in question.
/// @return  An AxEqDesignation enum type value: AXIAL, EQUATORIAL, or NEITHER
/// @author  Labonte <JwLabonte@jhu.edu>
/*chemical::rings::AxEqDesignation
is_atom_axial_or_equatorial( Pose const & pose, uint seqpos, uint query_atom )
{
return conformation::is_atom_axial_or_equatorial( pose.residue( seqpos ), query_atom );
}*/

// Is the query atom in this pose axial or equatorial or neither?
/// @details This function calculates an average plane and determines whether the coordinates of a given atom are
/// axial or equatorial to it (or neither).  The ring is requested from the corresponding Residue.
/// @param   <pose>:       The Pose containing the atoms in question.
/// @param   <query_atom>: The AtomID of the atom in question.
/// @return  An AxEqDesignation enum type value: AXIAL, EQUATORIAL, or NEITHER
/// @author  Labonte <JwLabonte@jhu.edu>
/*chemical::rings::AxEqDesignation
is_atom_axial_or_equatorial( Pose const & pose, id::AtomID const & query_atom )
{
return conformation::is_atom_axial_or_equatorial( pose.residue( query_atom.rsd() ), query_atom.atomno() );
}*/

void
set_bb_torsion( uint torsion_id, Pose & pose, core::Size sequence_position, core::Angle new_angle){

	if ( pose.residue( sequence_position).is_carbohydrate() ) {
		carbohydrates::set_glycosidic_torsion( torsion_id, pose, sequence_position, new_angle);
		return;
	}

	auto torsion_type = static_cast< id::MainchainTorsionType >(torsion_id);

	switch (torsion_type) {

	case id::phi_dihedral :
		pose.set_phi( sequence_position, new_angle );
		break;

	case id::psi_dihedral :
		pose.set_psi( sequence_position, new_angle );
		break;

	case id::omega_dihedral :
		pose.set_omega( sequence_position, new_angle );
		break;

		//Omega2+3 undefined for bb of AA.
	case id::omega2_dihedral :
		break;
	case id::omega3_dihedral :
		break;
	}


}

core::Angle
get_bb_torsion( uint torsion_id, Pose const & pose, core::Size sequence_position )
{
	if ( pose.residue( sequence_position ).is_carbohydrate() ) {
		return carbohydrates::get_glycosidic_torsion( torsion_id, pose, sequence_position );
	}

	auto torsion_type = static_cast< id::MainchainTorsionType >( torsion_id );

	switch ( torsion_type ) {
	case id::phi_dihedral :
		return pose.phi( sequence_position );
	case id::psi_dihedral :
		return pose.psi( sequence_position );
	case id::omega_dihedral :
		return pose.omega( sequence_position );
	case id::omega2_dihedral :
		return 0.0;
	case id::omega3_dihedral :
		return 0.0;
	}
	return 0.0;  // Code cannot reach here.
}

/// @brief Set bfactors in a pose PDBInfo
void
set_bfactors_from_atom_id_map(Pose & pose, id::AtomID_Map< Real > const & bfactors){
	for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {
		for ( core::Size atom_index = 1; atom_index <= pose.residue( resnum ).natoms(); ++atom_index ) {
			id::AtomID atm = id::AtomID( atom_index, resnum );
			if ( ! bfactors.has( atm ) ) continue;
			pose.pdb_info()->bfactor( resnum, atom_index, bfactors[ atm ]);
		}
	}

}

/////////////////////////////////////////////////////////////////////////////
bool
just_modeling_RNA( std::string const & sequence ) {
	// AMW TODO: note can't distinguish an all-NCAA pose from an all-NCNT pose
	std::string const rna_letters( "acgunZX" );
	for ( Size k = 1; k <= sequence.size(); k++ ) {
		if ( rna_letters.find( sequence[k-1] ) == std::string::npos ) return false;
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// don't allow Mg(2+) or HOH yet -- must be an easier way to figure out ligand or not.
bool
stepwise_addable_pose_residue( Size const n /*in pose numbering*/, pose::Pose const & pose ) {
	using namespace pose::full_model_info;
	runtime_assert( full_model_info_defined( pose ) );
	utility::vector1< Size > const & res_list = const_full_model_info( pose ).res_list();
	return stepwise_addable_residue( res_list[ n ],
		const_full_model_info( pose ).full_model_parameters()->non_standard_residue_map() );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
stepwise_addable_residue( Size const n /* in full model numbering*/, std::map< Size, std::string > const & non_standard_residue_map )
{
	auto it = non_standard_residue_map.find( n );
	if ( it != non_standard_residue_map.end() && ( it->second == "HOH" || it->second == "MG" ) ) return false;
	return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////
bool
effective_lower_terminus_based_on_working_res( Size const i,
	utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & res_list,
	utility::vector1< Size > const & cutpoint_open_in_full_model ){

	if ( working_res.size() == 0 ) return false; // not defined

	// decrement to the nearest cutpoint -- anything planning to be sampled along the way?
	for ( Size n = res_list[ i ] - 1; n >= 1; n-- ) {
		if ( cutpoint_open_in_full_model.has_value( n ) ) {
			return true;
		} else { // not a cutpoint yet.
			if ( working_res.has_value( n ) ) return false;
		}
	}
	return true;

	// this was a problem when skip_bulge was activated.
	// return ( working_res.size() > 0 && !sample_res.has_value( res_list[ i ] - 1 ) &&
	//      ( i == 1 || ( res_list[i - 1] < res_list[i] - 1 ) ) );
}

////////////////////////////////////////////////////////////////////////////////////////////////
bool
effective_upper_terminus_based_on_working_res( Size const i,
	utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & res_list,
	utility::vector1< Size > const & cutpoint_open_in_full_model,
	Size const nres_full){

	if ( working_res.size() == 0 ) return false; // not defined

	// increment to the nearest cutpoint -- anything planning to be sampled along the way?
	for ( Size n = res_list[ i ]; n <= nres_full; n++ ) {
		if ( cutpoint_open_in_full_model.has_value( n ) ) {
			return true;
		} else { // not a cutpoint yet.
			if ( working_res.has_value( n ) ) return false;
		}
	}
	return true;

	// this was a problem when skip_bulge was activated.
	//    ( sample_res.size() > 0 && !sample_res.has_value( res_list[i] + 1 ) &&
	//     ( i == pose.size() || ( res_list[i + 1] > res_list[i] + 1 ) ) ) ){
}

////////////////////////////////////////////////////////////////////////////////////
bool
definite_terminal_root( utility::vector1< Size > const & cutpoint_open_in_full_model,
	utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & res_list,
	Size const nres,
	Size const i ) {
	if ( res_list[ i ] == 1 ||
			cutpoint_open_in_full_model.has_value( res_list[ i ] - 1 ) ||
			effective_lower_terminus_based_on_working_res( i, working_res, res_list, cutpoint_open_in_full_model ) ) {
		// great, nothing will ever get prepended here.
		return true;
	}
	if ( res_list[ i ] == nres ||
			cutpoint_open_in_full_model.has_value( res_list[ i ] ) ||
			effective_upper_terminus_based_on_working_res( i, working_res, res_list, cutpoint_open_in_full_model, nres ) ) {
		// great, nothing will ever get appended here.
		return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////
bool
definite_terminal_root( pose::Pose const & pose, Size const i ){
	using namespace core::pose::full_model_info;
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	return definite_terminal_root( full_model_info.cutpoint_open_in_full_model(),
		full_model_info.working_res(),
		full_model_info.res_list(),
		full_model_info.size(), i );
}
////////////////////////////////////////////////////////////////////////////////////////////////
Size
get_definite_terminal_root( pose::Pose const & pose,
	utility::vector1< Size > const & partition_res /* should not be empty */,
	utility::vector1< Size > const & res_list,
	utility::vector1< Size > const & fixed_domain_map /* 0 in free; 1,2,... for separate fixed domains */,
	utility::vector1< Size > const & cutpoint_open_in_full_model,
	utility::vector1< Size > const & working_res ){
	for ( Size n = 1; n <= partition_res.size(); n++ ) {
		Size const i = partition_res[ n ];
		if ( !pose.fold_tree().possible_root( i ) ) continue;
		if ( definite_terminal_root( cutpoint_open_in_full_model, working_res, res_list, fixed_domain_map.size(), i ) ) {
			return i;
		}
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
Size
get_definite_terminal_root( pose::Pose const & pose,
	utility::vector1< Size > const & partition_res /* should not be empty */ ) {
	core::pose::full_model_info::FullModelInfo const & full_model_info = core::pose::full_model_info::const_full_model_info( pose );
	utility::vector1< Size > const & res_list = full_model_info.res_list();
	utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	utility::vector1< Size > const & working_res = full_model_info.working_res();
	return get_definite_terminal_root( pose, partition_res,
		res_list, fixed_domain_map, cutpoint_open_in_full_model, working_res );
}

////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
reorder_root_partition_res(
	utility::vector1< Size > const & root_partition_res /* should not be empty */,
	utility::vector1< Size > const & res_list,
	utility::vector1< Size > const & fixed_domain_map /* 0 in free; 1,2,... for separate fixed domains */ ){

	// reorder these residues based on size of fixed domains -- prefer *big* domains first.
	utility::vector1< std::pair< Size, Size > > domain_sizes;
	for ( Size domain = 1; domain <= utility::max( fixed_domain_map ); domain++ ) {
		Size nres_with_domain( 0 );
		for ( Size n = 1; n <= root_partition_res.size(); n++ ) {
			if ( fixed_domain_map[ res_list[ root_partition_res[ n ] ] ] == domain ) {
				nres_with_domain++;
			}
		}
		domain_sizes.push_back( std::make_pair( nres_with_domain, domain ) );
	}

	// largest block to smallest.
	std::sort( domain_sizes.begin(), domain_sizes.end() );
	std::reverse( domain_sizes.begin(), domain_sizes.end() );

	utility::vector1< Size > root_partition_res_reorder;
	for ( Size k = 1; k <= domain_sizes.size(); k++ ) {
		if ( domain_sizes[ k ].first == 0 ) continue; // no residues found.
		Size const domain = domain_sizes[ k ].second;
		for ( Size n = 1; n <= root_partition_res.size(); n++ ) {
			if ( fixed_domain_map[ res_list[ root_partition_res[ n ] ] ] == domain ) {
				root_partition_res_reorder.push_back( root_partition_res[ n ] );
			}
		}
	}

	// the rest of the residues that weren't in fixed domains.
	Size const domain( 0 );
	for ( Size n = 1; n <= root_partition_res.size(); n++ ) {
		if ( fixed_domain_map[ res_list[ root_partition_res[ n ] ] ] == domain ) {
			root_partition_res_reorder.push_back( root_partition_res[ n ] );
		}
	}

	runtime_assert( root_partition_res_reorder.size() == root_partition_res.size() );
	return root_partition_res_reorder;

}

////////////////////////////////////////////////////////////////////////////////////////////////
void
reroot( pose::Pose & pose,
	utility::vector1< Size > const & root_partition_res /* should not be empty */,
	utility::vector1< Size > const & res_list,
	utility::vector1< Size > const & preferred_root_res /* can be empty */,
	utility::vector1< Size > const & fixed_domain_map /* 0 in free; 1,2,... for separate fixed domains */,
	utility::vector1< Size > const & cutpoint_open_in_full_model,
	utility::vector1< Size > const & working_res ){

	using namespace core::kinematics;
	Size new_root( 0 );
	if ( root_partition_res.size() == 0 ) return;

	runtime_assert( root_partition_res.size() > 0 );
	utility::vector1< Size > root_partition_res_ordered = reorder_root_partition_res( root_partition_res, res_list, fixed_domain_map );

	for ( Size n = 1; n <= root_partition_res_ordered.size(); n++ ) {
		Size const i = root_partition_res_ordered[ n ];
		if ( preferred_root_res.has_value( res_list[ root_partition_res_ordered[ n ] ] ) ) {
			if ( !pose.fold_tree().possible_root( i ) ) {
				TR.Warning << res_list[ root_partition_res_ordered[ n ] ] << " specified as root res but not at a pose terminal. Cannot be used." << std::endl;
				continue;
			}
			new_root = i; break;
		}
	}

	// next preference: roots that are definitely terminal -- nothing will be built past them.
	if ( new_root == 0 ) {
		new_root = get_definite_terminal_root( pose, root_partition_res_ordered,
			res_list, fixed_domain_map, cutpoint_open_in_full_model, working_res );
	}

	// if all else fails...
	if ( new_root == 0 ) {
		for ( Size n = 1; n <= root_partition_res_ordered.size(); n++ ) {
			Size const i = root_partition_res_ordered[ n ];
			if ( !pose.fold_tree().possible_root( i ) ) continue;
			new_root = i; break;
		}
	}
	runtime_assert( new_root > 0 );

	FoldTree f = pose.fold_tree();
	if ( new_root == f.root() ) return;
	f.reorder( new_root );
	pose.fold_tree( f );
}



} // pose
} // core
