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

// Package headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/PositionConservedResiduesStore.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Conformation.hh>
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

// C/C++ headers
#include <cmath>
#include <iostream>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <boost/functional/hash.hpp>

namespace core {
namespace pose {

static THREAD_LOCAL basic::Tracer TR( "core.pose.util" );

void
append_pose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	bool new_chain
){
	append_subpose_to_pose(pose1, pose2, 1, pose2.size(), new_chain);
}

void
append_subpose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	core::Size start_res,
	core::Size end_res,
	bool new_chain
){
	if ( pose2.size()<start_res ) {
		TR.Error << "Provided starting residue number " << start_res
			<< " less than number residues in appended pose. Nothing to do." << std::endl;
	}
	pose1.append_residue_by_jump(pose2.residue(start_res), pose1.size() , "", "", new_chain);
	for ( core::Size i=start_res+1; i<=end_res; ++i ) {
		if ( pose2.residue(i).is_lower_terminus() ) {
			if ( i > 1 && pose2.chain(i) == pose2.chain(i-1) ) {
				pose1.append_residue_by_jump(pose2.residue(i), pose1.size(), "","", false);
			} else {
				if ( pose2.residue(i).is_protein() ) {
					pose1.append_residue_by_bond(pose2.residue(i));
				} else {
					bool new_chain =  pose2.chain(i) != pose2.chain(i-1);
					pose1.append_residue_by_jump(pose2.residue(i), i-1, "","", new_chain);
				}
			}
		} else {
			pose1.append_residue_by_bond(pose2.residue(i));
		}
	}
}

void jumps_from_pose(core::pose::Pose const & pose, Jumps & jumps) {
	for ( Size i = 1; i <= pose.num_jump(); ++i ) {
		jumps.insert(i);
	}
}

void remove_virtual_residues(core::pose::Pose & pose) {
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue_type(i).name() == "VRT" ) {
			pose.conformation().delete_residue_slow(i);
		}
	}
}

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

void
create_subpose(
	Pose const & src,
	utility::vector1< Size > const & positions,
	kinematics::FoldTree const & f,
	Pose & pose
)
{
	Size const nres( f.nres() );
	debug_assert( nres == positions.size() );

	pose.clear();

	for ( Size i=1; i<= nres; ++i ) {
		Size const seqpos( positions[i] );
		conformation::Residue const & rsd( src.residue( seqpos ) );
		// If the residue and the previous residue are bonded in the source pose, they should be bonded in the new pose
		if ( i>1 && rsd.is_polymer_bonded( positions[ i-1 ] ) ) {
			pose.append_residue_by_bond( rsd );
		} else {
			pose.append_residue_by_jump( rsd, 1 );
		}
		if ( i>1 ) {
			// check if this residue should be in a new chain. not a perfect check...
			conformation::Residue const & prev_rsd( src.residue( positions[i-1] ) );
			if ( prev_rsd.is_upper_terminus() || rsd.is_lower_terminus() || prev_rsd.chain() != rsd.chain() ) {
				debug_assert( pose.size() == i );
				pose.conformation().insert_chain_ending( i-1 );
			}
		}
	}

	// now set the desired foldtree
	pose.fold_tree(f);

}


///////////////////////////////////////////////////////////////////////////////
// A bit like create_subpose() but does not require fold tree, and
// a little optimized to handle some weird terminal variants and RNA jumps.
void
pdbslice( core::pose::Pose & new_pose,
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & slice_res )
{
	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::chemical::rna;

	new_pose.clear();

	for ( Size i = 1; i <= slice_res.size(); i++ ) {

		ResidueOP residue_to_add = pose.residue( slice_res[ i ] ).clone() ;

		if ( (i > 1 &&  ( slice_res[i] != slice_res[i-1] + 1 )) /*new segment*/ ||
				residue_to_add->is_lower_terminus() ||
				residue_to_add->has_variant_type( N_ACETYLATION ) ||
				!residue_to_add->is_polymer() ||
				(i>1 && pose.fold_tree().is_cutpoint( slice_res[i-1] ) ) ) {

			if ( residue_to_add->is_RNA() && (i>1) && new_pose.residue_type(i-1).is_RNA() ) {

				new_pose.append_residue_by_jump(  *residue_to_add, i-1,
					chi1_torsion_atom( new_pose.residue_type(i-1) ),
					chi1_torsion_atom( residue_to_add->type() ), true /*new chain*/ );
			} else {

				new_pose.append_residue_by_jump(  *residue_to_add, i-1, "", "", true /*new chain*/ );
			}
		} else {

			new_pose.append_residue_by_bond(  *residue_to_add  ) ;
		}
	}

	using namespace core::pose::full_model_info;
	if ( full_model_info_defined( pose ) ) {
		FullModelInfoOP full_model_info = const_full_model_info( pose ).clone_info();
		utility::vector1< Size > const & res_list = full_model_info->res_list();
		utility::vector1< Size > new_res_list;
		for ( Size n = 1; n <= new_pose.size(); n++ ) {
			new_res_list.push_back( res_list[ slice_res[ n ] ] );
		}
		full_model_info->set_res_list( new_res_list );
		set_full_model_info( new_pose, full_model_info );
	}

	PDBInfoCOP pdb_info = pose.pdb_info();
	if ( pdb_info != nullptr ) {
		utility::vector1< Size > new_numbering;
		utility::vector1< char > new_chains;
		for ( Size n = 1; n <= slice_res.size(); n++ ) {
			new_numbering.push_back( pdb_info->number( slice_res[ n ] ) );
			new_chains.push_back( pdb_info->chain( slice_res[ n ] ) );
		}
		PDBInfoOP new_pdb_info( new PDBInfo( new_pose, true /*init*/ ) );
		new_pdb_info->set_numbering( new_numbering );
		new_pdb_info->set_chains( new_chains );
		new_pose.pdb_info( new_pdb_info );
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pdbslice( core::pose::Pose & pose,
	utility::vector1< core::Size > const & slice_res ){

	pose::Pose mini_pose;
	pdbslice( mini_pose, pose, slice_res );
	pose = mini_pose;

}


/// @brief  Convert a jump to the given Residue in a FoldTree to a chemical Edge.
/// @param  <end_resnum>: the Residue index of the end of the Jump\n
/// @param  <ft>: the FoldTree being modified
/// @return true if a chemical connection was found and the passed FoldTree was changed
/// @note   This function is used by set_reasonable_foldtree() for the primary purpose of rebuilding default FoldTrees
/// in cases of branching or chemical conjugation.
/// @author Labonte <JWLabonte@jhu.edu>
/// @author Joseph Harrison
bool
change_jump_to_this_residue_into_chemical_edge(
	core::uint end_resnum,
	core::pose::Pose const & pose,
	core::kinematics::FoldTree & ft )
{
	using namespace std;
	using namespace core::chemical;
	using namespace conformation;

	Residue const & end_residue( pose.residue( end_resnum ) );
	bool found_chemical_connection( false );

	if ( TR.Trace.visible() ) {
		TR.Trace << "Searching for ResidueConnection from residue " << end_resnum << endl;
	}
	// Loop through the ResidueConnections of the end Residue of the old Jump to find and attach the other end of the
	// new chemical Edge.
	Size const n_connections( end_residue.type().n_possible_residue_connections() );
	for ( uint conn_index( 1 ); conn_index <= n_connections; ++conn_index ) {
		if ( end_residue.connection_incomplete( conn_index ) ) {
			continue;  // Allow incomplete connections for design.
		}
		uint const start_resnum( end_residue.connect_map( conn_index ).resid() );
		if ( start_resnum < end_resnum ) {
			Residue const & start_residue( pose.residue( start_resnum ) );
			// Ensure that the connection is either a polymer branch or a ligand of the same chain.
			if ( ( start_residue.is_branch_point() && end_residue.is_branch_lower_terminus() ) ||
					( ( ! end_residue.is_polymer() ) && start_residue.chain() == end_residue.chain() ) ) {
				uint start_atm_index = start_residue.connect_atom( end_residue );
				uint end_atm_index = end_residue.connect_atom( start_residue );

				string start_atm_name = start_residue.atom_name( start_atm_index );
				string end_atm_name = end_residue.atom_name( end_atm_index );

				if ( TR.Trace.visible() ) {
					TR.Trace << "Adding  EDGE " <<
						start_resnum << ' ' << end_resnum << ' ' << start_atm_name << ' ' << end_atm_name <<
						" to the new FoldTree." << endl;
				}
				ft.add_edge( start_resnum, end_resnum, start_atm_name, end_atm_name );

				found_chemical_connection = true;
				break;
			}
		}
	}
	return found_chemical_connection;
}


/// @details If all ligand residues and polymer branches have been appended by a jump, this method creates a new
/// FoldTree without jumps through ligands, using CHEMICAL EDGEs instead.
void
set_reasonable_fold_tree( pose::Pose & pose )
{
	// An empty pose doesn't have jumps through ligands.
	// (Will encounter a SegFault otherwise)
	if ( pose.size() == 0 ) {
		return;
	}

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
			TR.Trace << "checking if if " << e << " is reasonable for this Pose..." << endl;
		}
		if ( e.is_jump() ) {
			uint const ii( e.stop() );  // the residue position at the end of the jump
			// Is this jump to a covalent ligand residue or forward-direction polymer branch?
			if ( ( pose.residue( ii ).is_branch_lower_terminus() ) /* polymer branch in forward direction */ ||
					( ! pose.residue( ii ).is_polymer() ) /* covalent ligand residue */ ) {
				if ( TR.Trace.visible() ) {
					TR.Trace << e << " was initially set as a Jump to a branch lower terminus or ligand." << endl;
				}
				bool jump_was_changed( change_jump_to_this_residue_into_chemical_edge( ii, pose, newft ) );

				if ( jump_was_changed ) {
					continue;  // The Edge has already been added to newft by the function above.
				} else /* We couldn't find a chemical connection. */ {
					if ( pose.residue_type( ii ).n_possible_residue_connections() > 0 ) {
						TR.Warning << "Can't find a chemical connection for residue " << ii << " " <<
							pose.residue_type( ii ).name() << endl;
					}
				}
			} else {
				if ( TR.Trace.visible() ) {
					TR.Trace << e << " was initially set as a Jump to a lower terminus." << endl;
				}
				// Figure out if the jump should connect to the LAST residue of this Edge/chain and not the 1st,
				// i.e., is this a C-terminal connection?

				Edge const next_e( *(i + 1) );
				if ( ii == uint( next_e.start() ) && pose.residue( ii ).is_lower_terminus() ) {
					uint const jj( next_e.stop() );  // the residue position at the end of the next Edge
					if ( pose.residue_type( jj ).is_branch_lower_terminus() ) {
						bool jump_was_changed( change_jump_to_this_residue_into_chemical_edge( jj, pose, newft ) );

						if ( jump_was_changed ) {
							// The CHEMICAL Edge has already been added to newft by the function above.
							// Now, we must add a reverse-direction Edge and skip making the usual forward-direction
							// Edge.
							// The Edge has already been added to newft by the function above.
							if ( TR.Trace.visible() ) {
								TR.Trace << "Adding  EDGE " << jj << ' ' << ii << " -1  to the new FoldTree." << endl;
							}
							newft.add_edge( jj, ii, Edge::PEPTIDE);
							prevent_forward_edge = true;
							continue;  // Skip the add_edge() method below; we are done here.
						} else /* We couldn't find a chemical connection. */ {
							if ( pose.residue_type( jj ).n_possible_residue_connections() > 0 ) {
								TR.Warning << "Can't find a chemical connection for residue " << jj << " " <<
									pose.residue_type( jj ).name() << endl;
							}
						}
					}
				}
			}
			// It's just a normal inter-chain jump or we couldn't find a chemical bond to make. Increment the label.
			e.label() = ++last_jump_id;
		} else {
			if ( TR.Trace.visible() ) {
				TR.Trace << e << " was not a Jump." << endl;
			}
		}
		if ( ! prevent_forward_edge ) {
			if ( TR.Trace.visible() ) {
				TR.Trace << "Adding " << e << " to the new FoldTree." << endl;
			}
			newft.add_edge( e );
		} else {
			prevent_forward_edge = false;  // Adding forward Edge was skipped; reset for next Edge.
		}
	}  // next Edge

	runtime_assert( newft.size() > 0 || pose.size() == 0 );  // A valid fold tree must have >= 1 edges.

	if ( TR.Debug.visible() ) {
		TR.Debug << "new fold tree " << newft << endl;
	}

	pose.fold_tree( newft );
}


////////////////////////////////////////////////////////////////////////////
void
partition_pose_by_jump(
	pose::Pose const & src,
	int const jump_number,
	pose::Pose & partner1,
	pose::Pose & partner2
)
{
	Size const nres( src.size() );

	// split src pose's foldtree
	kinematics::FoldTree f1, f2;
	src.fold_tree().partition_by_jump( jump_number, f1, f2 );

	TR << src.fold_tree() << '\n' << f1 << '\n' << f2 << std::endl;

	// identify residues in the two partners
	ObjexxFCL::FArray1D_bool partner1_pos( nres, false ); // FARRAY! DOH!!
	src.fold_tree().partition_by_jump( jump_number, partner1_pos );

	utility::vector1< Size > partner1_pos_list, partner2_pos_list;
	for ( Size i=1; i<= nres; ++i ) {
		if ( partner1_pos(i) ) partner1_pos_list.push_back( i );
		else partner2_pos_list.push_back( i );
	}

	create_subpose( src, partner1_pos_list, f1, partner1 );

	create_subpose( src, partner2_pos_list, f2, partner2 );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Crude way to guess secondary structure given a pose. This function
/// sets the sec_struct array of pose.conformation_ to the result of the
/// guesswork. This has been ported directly from rosetta++.
void
set_ss_from_phipsi(
	pose::Pose & pose
)
{
	// ss       :ss = 1 helix, ss = 2 sheet, ss = 3 other
	const int sstemp_offset=3;
	utility::vector1 < int > sstemp( sstemp_offset*2 + pose.size() );
	utility::vector1 < int > ss( pose.size() );

	sstemp[sstemp_offset-1] = 3; // assign loop to fictious residues at ends of chain
	sstemp[sstemp_offset+0] = 3;
	sstemp[sstemp_offset+pose.size()+1] = 3;
	sstemp[sstemp_offset+pose.size()+2] = 3;

	for ( Size i = 1; i <= pose.size(); ++i ) {

		// <murphp>
		if ( !pose.residue_type(i).is_protein() ) { // make sure we don't inquire about the phi/psi of a non-protein residue
			sstemp[sstemp_offset+i] = 3;
		} else {
			// </murphp>

			if ( pose.phi(i) < -20.0 && pose.psi(i) > -90.0 && pose.psi(i) < -10.0 ) {
				sstemp[sstemp_offset+i] = 1;
			} else if ( pose.phi(i) < -20.0 && (pose.psi(i) > 20.0 || pose.psi(i) < -170.0) ) {
				sstemp[sstemp_offset+i] = 2;
			} else {
				sstemp[sstemp_offset+i] = 3;
			}

		}
	}

	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( sstemp[sstemp_offset+i] == 2 ) {
			if ( sstemp[sstemp_offset+i-1] == 2 && sstemp[sstemp_offset+i+1] == 2 ) {
				ss[i] = 2;
			} else if ( sstemp[sstemp_offset+i-1] == 2 && sstemp[sstemp_offset+i-2] == 2 ) {
				ss[i] = 2;
			} else if ( sstemp[sstemp_offset+i+1] == 2 && sstemp[sstemp_offset+i+2] == 2 ) {
				ss[i] = 2;
			} else {
				ss[i] = 3;
			}
		} else if ( sstemp[sstemp_offset+i] == 1 ) {
			if ( sstemp[sstemp_offset+i-1] == 1 && sstemp[sstemp_offset+i+1] == 1 ) {
				ss[i] = 1;
			} else if ( sstemp[sstemp_offset+i-1] == 1 && sstemp[sstemp_offset+i-2] == 1 ) {
				ss[i] = 1;
			} else if ( sstemp[sstemp_offset+i+1] == 1 && sstemp[sstemp_offset+i+2] == 1 ) {
				ss[i] = 1;
			} else {
				ss[i] = 3;
			}
		} else {
			ss[i] = 3;
		}
	}

	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( ss[i] == 1 ) {
			pose.set_secstruct(i,'H');
		} else if ( ss[i] == 2 ) {
			pose.set_secstruct(i,'E');
		} else {
			pose.set_secstruct(i,'L');
		}
	}
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
		TR.Warning << "WARNING: addVirtualResAsRoot() called with empty pose!" << std::endl;
		return;
	}

	// check for terminal ligands
	int last_peptide_res = nres;
	while ( !pose.residue( last_peptide_res ).is_polymer() )
			last_peptide_res--;

	// try to avoid putting the vrt too close to termini
	int i_min = 1;

#ifdef WIN32
	int r_start = static_cast< int > ( std::floor(   static_cast< double > (last_peptide_res) /3. ) );
	int r_end   = static_cast< int > ( std::ceil ( 2.* static_cast< double > (last_peptide_res)/3. ) );
#else
	int r_start = static_cast< int > ( std::floor(   last_peptide_res/3 ) );
	int r_end   = static_cast< int > ( std::ceil ( 2*last_peptide_res/3 ) );
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
	bool fullatom = pose.is_fullatom();
	core::chemical::ResidueTypeSetCOP const &residue_set(
		core::chemical::ChemicalManager::get_instance()->residue_type_set
		( fullatom ? core::chemical::FA_STANDARD : core::chemical::CENTROID )
	);
	core::chemical::ResidueTypeCOP rsd_type( residue_set->get_representative_type_name3("VRT") );
	if ( ! rsd_type ) {
		utility_exit_with_message("Cannot find residue type VRT" );
	}
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
		TR.Warning << "WARNING: addVirtualResAsRoot() called with empty pose!" << std::endl;
		return;
	}
	numeric::xyzVector< core::Real > massSum = get_center_of_mass( pose );
	return addVirtualResAsRoot(massSum, pose);
}

/// @brief  Reads the comments from the pdb file and adds it into comments
void
read_comment_pdb(
	std::string const &file_name,
	core::pose::Pose & pose
) {
	utility::io::izstream data(file_name);
	if ( !data ) {
		TR<<"ERROR! PDB SCAFOLD NAME NOT FOUND!!"<<std::endl;
		utility_exit();
	}
	std::string line;
	while ( getline( data, line ) ) {
		if ( line != "##Begin comments##" ) {
			continue;
		}
		getline( data, line );
		while ( line != "##End comments##" ) {
			//TR<<"Testing read comments! :"<<line<<std::endl;
			utility::vector1<std::string> comment_line(utility::string_split(line,' '));
			core::pose::add_comment(pose,comment_line[1],comment_line[2]);
			getline( data, line );
		}
	}
}

void
dump_comment_pdb(
	std::string const &file_name,
	core::pose::Pose const& pose
) {
	std::ofstream out( file_name.c_str() );
	pose.dump_pdb(out);
	// verbose output
	out << "END\n";
	out << "##Begin comments##" << std::endl;
	using namespace std;
	map< string, string > const comments = core::pose::get_all_comments(pose);
	for ( auto const & comment : comments ) {
		out << comment.first<<" "<<comment.second << std::endl;
	}
	out << "##End comments##" << std::endl;
	out.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
hasPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const & name )
{
	using basic::datacache::CacheableStringFloatMap;
	using basic::datacache::CacheableStringFloatMapCOP;

	// make sure that the pose has one of these.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		return false;
	}

	CacheableStringFloatMapCOP data
		= utility::pointer::dynamic_pointer_cast< CacheableStringFloatMap const >
		( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
	assert( data.get() != nullptr );

	return (  data->map().find( name ) != data->map().end() );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool getPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const & name,
	core::Real & value
)
{
	using basic::datacache::CacheableStringFloatMap;
	using basic::datacache::CacheableStringFloatMapCOP;

	// make sure that the pose has one of these.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		return false;
	}

	CacheableStringFloatMapCOP data
		= utility::pointer::dynamic_pointer_cast< CacheableStringFloatMap const >
		( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
	debug_assert( data.get() != nullptr );

	auto it = data->map().find( name );
	if ( it == data->map().end() ) {
		return false;
	}
	value = it->second;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
Real
getPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const & name ) {
	Real value;
	runtime_assert( getPoseExtraScore( pose, name, value ) );
	return value;
}

bool getPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const & name,
	std::string & value
)
{
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapCOP;

	// make sure that the pose has one of these.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ) {
		return false;
	}

	CacheableStringMapCOP data
		= utility::pointer::dynamic_pointer_cast< CacheableStringMap const >
		( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) );
	debug_assert( data.get() != nullptr );

	auto it = data->map().find( name );
	if ( it == data->map().end() ) {
		return false;
	}
	value = it->second;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void setPoseExtraScore(
	core::pose::Pose & pose,
	std::string const & name,
	core::Real value
)
{
	using basic::datacache::CacheableStringFloatMap;
	using basic::datacache::CacheableStringFloatMapOP;
	using basic::datacache::DataCache_CacheableData;

	// make sure that the pose has one of these.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		pose.data().set(
			core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA,
			DataCache_CacheableData::DataOP( new CacheableStringFloatMap() )

		);
	}

	CacheableStringFloatMapOP data
		=  utility::pointer::dynamic_pointer_cast< CacheableStringFloatMap >
		( pose.data().get_ptr(core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA) );

	runtime_assert( data.get() != nullptr );
	data->map()[name] = value;
}

void setPoseExtraScore(
	core::pose::Pose & pose,
	std::string const & name,
	std::string const & value
)
{
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;
	using basic::datacache::DataCache_CacheableData;

	// make sure that the pose has one of these.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ) {
		pose.data().set(
			core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA,
			DataCache_CacheableData::DataOP( new CacheableStringMap() )

		);
	}

	CacheableStringMapOP data
		=  utility::pointer::dynamic_pointer_cast< CacheableStringMap >
		( pose.data().get_ptr(core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA) );

	runtime_assert( data.get() != nullptr );
	data->map()[name] = value;
}

void add_comment(
	core::pose::Pose & pose,
	std::string const & key,
	std::string const & val
)
{
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;
	using basic::datacache::DataCache_CacheableData;

	// make sure that the pose has a map of strings
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) ) {
		pose.data().set(
			core::pose::datacache::CacheableDataType::STRING_MAP,
			DataCache_CacheableData::DataOP( new CacheableStringMap() )

		);
	}

	CacheableStringMapOP data
		=  utility::pointer::dynamic_pointer_cast< CacheableStringMap >
		( pose.data().get_ptr(core::pose::datacache::CacheableDataType::STRING_MAP) );

	runtime_assert( data.get() != nullptr );
	data->map()[key] = val;
} // add_comment

void add_score_line_string(
	core::pose::Pose & pose,
	std::string const & key,
	std::string const & val
)
{
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;
	using basic::datacache::DataCache_CacheableData;

	// make sure that the pose has a map of strings
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS ) ) {
		pose.data().set(
			core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS,
			DataCache_CacheableData::DataOP( new CacheableStringMap() )
		);
	}

	CacheableStringMapOP data
		=  utility::pointer::dynamic_pointer_cast< CacheableStringMap >
		( pose.data().get_ptr(core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS) );

	runtime_assert( data.get() != nullptr );
	data->map()[key] = val;
}

void clearPoseExtraScores(
	core::pose::Pose & pose
) {

	using basic::datacache::DataCache_CacheableData;

	{
		using basic::datacache::CacheableStringFloatMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA,
			DataCache_CacheableData::DataOP( new CacheableStringFloatMap() )
		);
	}

	{
		using basic::datacache::CacheableStringMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA,
			DataCache_CacheableData::DataOP( new CacheableStringMap() )
		);
	}

	{
		using basic::datacache::CacheableStringMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::STRING_MAP,
			DataCache_CacheableData::DataOP( new CacheableStringMap() )

		);
	}

	{
		using basic::datacache::CacheableStringMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS,
			DataCache_CacheableData::DataOP( new CacheableStringMap() )

		);
	}
}

void clearPoseExtraScore(
	core::pose::Pose & pose,
	std::string const & name
)
{
	using basic::datacache::CacheableStringFloatMap;
	using basic::datacache::CacheableStringFloatMapOP;
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;
	using basic::datacache::DataCache_CacheableData;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		CacheableStringFloatMapOP data
			= utility::pointer::dynamic_pointer_cast< CacheableStringFloatMap >
			( pose.data().get_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
		debug_assert( data.get() != nullptr );

		if ( data->map().find( name ) != data->map().end() ) data->map().erase( name );
	}

	if ( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ) {
		CacheableStringMapOP data
			= utility::pointer::dynamic_pointer_cast< CacheableStringMap >
			( pose.data().get_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) );
		debug_assert( data.get() != nullptr );

		if ( data->map().find( name ) != data->map().end() ) data->map().erase( name );
	}
}

bool get_comment(
	core::pose::Pose const & pose,
	std::string const & key,
	std::string & val
) {
	using std::map;
	using std::string;
	map< string, string > comment_map = get_all_comments( pose );

	if ( comment_map.find( key ) == comment_map.end() ) {
		return false;
	}

	val = comment_map[ key ];
	return true;
}

bool get_score_line_string(
	core::pose::Pose const & pose,
	std::string const & key,
	std::string & val
) {
	using std::map;
	using std::string;
	map< string, string > score_line_strings_map = get_all_score_line_strings( pose );

	if ( score_line_strings_map.find( key ) == score_line_strings_map.end() ) {
		return false;
	}

	val = score_line_strings_map[ key ];
	return true;
}

void delete_comment(
	core::pose::Pose & pose,
	std::string const & key
) {
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) ) {
		CacheableStringMapOP data
			=  utility::pointer::dynamic_pointer_cast< CacheableStringMap >
			( pose.data().get_ptr(core::pose::datacache::CacheableDataType::STRING_MAP) );
		std::map< std::string, std::string >::iterator it;
		it = data->map().find(key);
		if ( it != data->map().end() ) {
			data->map().erase(it);
		}
	}
}

std::map< std::string, std::string >
get_all_score_line_strings(
	core::pose::Pose const & pose
)
{
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapCOP;

	std::map< std::string, std::string > score_line_strings;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS ) ) {
		CacheableStringMapCOP data
			= utility::pointer::dynamic_pointer_cast< CacheableStringMap const >
			( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS ) );
		score_line_strings = data->map();
		runtime_assert( data.get() != nullptr );
	}
	return score_line_strings;
}

std::map< std::string, std::string >
get_all_comments(
	core::pose::Pose const & pose
)
{
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapCOP;

	std::map< std::string, std::string > comments;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) ) {
		CacheableStringMapCOP data
			= utility::pointer::dynamic_pointer_cast< CacheableStringMap const >
			( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::STRING_MAP ) );
		comments = data->map();
		runtime_assert( data.get() != nullptr );
	}
	return comments;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< char > read_psipred_ss2_file( pose::Pose const & pose ) {
	using namespace basic::options;

	std::string filename( option[ OptionKeys::in::file::psipred_ss2 ]().name() );
	return read_psipred_ss2_file( pose, filename );
}

utility::vector1< char > read_psipred_ss2_file( pose::Pose const & pose, std::string const & filename ) {
	using namespace basic::options;

	utility::vector1< char > secstructs;
	utility::io::izstream data( filename );

	if ( !data ) {
		TR.Warning << "Cannot open psipred_ss2 file " << filename << std::endl;
		return secstructs;
	}

	std::string line;
	Size count(0);
	while ( getline( data, line ) ) {
		if ( line[0] == '#' || line == "" ) {
			continue;
		}
		std::istringstream line_stream( line );
		Size pos;
		char aa, sec;
		line_stream >> pos >> aa >> sec;
		count++;
		if ( sec != 'H' && sec != 'E' && sec != 'C' ) {
			TR.Warning << "unrecognized secstruct char : " << sec << " at seqpos " << count << std::endl;
		}
		if ( sec == 'C' ) {
			secstructs.push_back( 'L' );
		} else {
			secstructs.push_back( sec );
		}
	}

	// chu get number of protein residues
	Size nres=0;
	for ( Size i =1 ; i <= pose.size(); ++i ) {
		if ( pose.residue(i).is_protein() ) nres++;
	}

	debug_assert( secstructs.size() == nres);
	if ( secstructs.size() != nres ) {
		secstructs.clear();
	}

	return secstructs;
}
////////////////////////////////////////////////////////////////////////////////////////////////////


void conf2pdb_chain_default_map( core::pose::Pose const & pose, std::map<int,char> & chainmap ) {
	chainmap.clear();
	char letter = 'A';
	for ( core::Size i = 1; i <= pose.conformation().num_chains(); ++i ) {
		chainmap[i] = letter;
		if ( 'Z'==letter ) utility_exit_with_message("too many chains to map to letters!!!");
		letter = static_cast<char>(letter + 1);
	}
}

/// @brief get Conformation chain -> PDBInfo chain mapping
/// @remarks Any chains whose PDBInfo chain records are marked entirely as
///  PDBInfo::empty_record() will be mapped to that character.  Note that
///  Conformation -> PDBInfo is always unique, but the reverse may not be true.
/// @return the mapping if PDBInfo available and chains appear consistent,
///  otherwise returns an empty mapping
std::map< int, char > conf2pdb_chain( core::pose::Pose const & pose ) {
	using core::Size;
	using core::pose::PDBInfo;
	typedef std::map< int, char > Conf2PDB;

	Conf2PDB conf2pdb;

	if ( !pose.pdb_info().get() ) {
		TR.Warning << "WARNING: conf2pdb_chain(): PDBInfo does not exist, returning default map 1=A, 2=B, ..." << std::endl;
		conf2pdb_chain_default_map(pose,conf2pdb);
		return conf2pdb;
	}

	for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
		int const conf = pose.chain( i );
		char const pdb = pose.pdb_info()->chain( i );

		auto c2p = conf2pdb.find( conf );
		if ( c2p != conf2pdb.end() ) { // must check if existing record inconsistent
			if ( c2p->second != pdb ) {

				// three cases:
				//  (1) replace an existing empty record
				//  (2) it's an unneeded empty record, so continue
				//  (3) there is an actual problem
				if ( pdb != PDBInfo::empty_record() && c2p->second == PDBInfo::empty_record() ) {
					// replace the record
					c2p->second = pdb;
				} else if ( pdb == PDBInfo::empty_record() ) {
					continue; // skip empty record
				} else {
					// something is inconsistent
					TR.Warning << "WARNING: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
					TR.Warning << "WARNING: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! HERE BE DRAGONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
					TR.Warning << "WARNING: conf2pdb_chain(): chain mapping inconsistent, returning default map p 1=A, 2=B, ... ";
					TR.Warning << "WARNING: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
					TR.Warning << "existing " << c2p->first << " -> " << c2p->second << "  |  ";
					TR.Warning << "new " << conf << " -> " << pdb << std::endl;
					conf2pdb_chain_default_map(pose,conf2pdb);
					return conf2pdb;
				}

			}
		} else { // record doesn't exist yet
			conf2pdb[ conf ] = pdb;
		}

	} // foreach residue

	debug_assert( conf2pdb.size() == pose.conformation().num_chains() );
	return conf2pdb;
}

/// @brief Get all the chains from conformation
utility::vector1< int > get_chains( core::pose::Pose const & pose ) {

	// get map between conformation chain and PDB chainids
	std::map< int, char > chain_map = conf2pdb_chain( pose );

	// create empty vector
	utility::vector1< int > keys;

	// go through map and push back into vector
	for ( auto & it : chain_map ) {
		keys.push_back( it.first );
	}
	return keys;
}

/// @brief compute last residue number of a chain
core::Size chain_end_res( Pose const & pose, core::Size const chain ) {

	// check whether chain exists in pose
	if ( ! has_chain( chain, pose ) ) {
		TR << "??? Chain should be " << chain << "???" << std::endl;
		utility_exit_with_message("Cannot get chain end residue for chain that doesn't exist. Quitting");
	}

	int int_chain = static_cast< int >( chain );
	int end_res(0);
	int nres( static_cast< int > ( pose.size() ) );

	// go through residues
	for ( int i = 1; i <= nres ; ++i ) {

		if ( i == nres && pose.chain(i) == int_chain ) {
			end_res = nres;
		} else if ( pose.chain( i ) == int_chain && pose.chain( i+1 ) != int_chain ) {
			end_res = i;
		}
	}

	return static_cast< core::Size >( end_res );

} // chain end residue number

/// @brief compute last residue numbers of all chains
utility::vector1< core::Size > chain_end_res( Pose const & pose ) {

	utility::vector1< int > chains( get_chains( pose ) );
	utility::vector1< core::Size > end_resnumber;

	for ( core::Size i = 1; i <= chains.size(); ++i ) {
		core::Size end_res( chain_end_res( pose, static_cast< core::Size >( chains[ i ] ) ) );
		end_resnumber.push_back( end_res );
	}
	return end_resnumber;

} // chain end residue numbers


/// @brief Compute uniq chains in a complex
/// @detail Returns a vector of pose length with true/false of uniq chain
utility::vector1< bool > compute_unique_chains( Pose & pose ) {

	// initilize vector of pose length with false
	utility::vector1< bool > uniq( pose.size(), true );

	// get chains, initialize uniq chains and vector of chain sequences
	utility::vector1< core::Size > chains( get_chains( pose ) );
	utility::vector1< bool > uniq_chains( chains.size(), true );
	utility::vector1< std::string > sequences;

	// get sequences of all chains into a vector of strings
	for ( core::Size i = 1; i <= chains.size(); ++i ) {
		std::string seq( pose.chain_sequence( i ) );
		sequences.push_back( seq );
	}

	// compare sequences with respect to each other
	// this is the simplest and dirties way to compute whether the chains
	//  are similar, it does NOT do a proper sequence alignment
	// the assumptions are:
	// - if the sequences have different lengths, they are different
	//   (this obviously goes wrong if there are single residue insertions
	//   or deletions while the rest of the sequences are the same!)
	// - if the sequences have the same length and the sequence identity is
	//   above 95%, then they are "not unique"

	// go through vectors of chain sequences to compare them
	for ( core::Size i = 1; i <= sequences.size(); ++i ) {

		std::string seq1 = sequences[ i ];

		// no double counting
		for ( core::Size j = i+1; j <= sequences.size(); ++j ) {

			std::string seq2 = sequences[ j ];

			// make sure that the sequences are of same length
			core::Size num_ident_res( 0 );
			if ( seq1.size() == seq2.size() ) {

				TR << "sequences " << i << " and " << j << " have same length" << std::endl;

				// go through sequence, std::string indexing from 0
				for ( core::Size k = 0; k <= seq1.size(); ++k ) {

					if ( seq1[ k ] == seq2[ k ] ) {
						++num_ident_res;
					}
				}

				// compute sequence identity
				core::Real seqid = num_ident_res / seq1.size();
				TR << "seqid " << seqid << std::endl;

				// if sequence identity >= 95%, add a true to the uniq_chains
				// vector
				if ( seqid >= 0.95 ) {
					TR << "adding a false to chain " << i << std::endl;
					uniq_chains[ i ] = false;
				}

			} // same length sequences
		} // iterate over chains vector
	} // iterate over chains vector

	// go through uniq chains vector
	for ( core::Size i = 1; i <= uniq_chains.size(); ++i ) {

		// go through residues
		for ( core::Size j = 1; j <= pose.size(); ++j ) {

			// if residue belongs to uniq chain, set residue to true in uniq seq vector
			if ( uniq_chains[ i ] == false && chains[ i ] == static_cast< core::Size >( pose.chain( j ) ) ) {
				uniq[ j ] = false;
			}
		}
	}

	return uniq;

} // compute unique chains

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

/// @brief renumber PDBInfo based on Conformation chains; each chain starts from 1
/// @param[in,out] pose The Pose to modify.
/// @param[in] fix_chains If true, the procedure will attempt to fix any empty record
///  characters it finds in the PDBInfo. (default true)
/// @param[in] start_from_existing_numbering If true, will attempt to start each
///  chain from the existing numbering in the PDBInfo.  E.g. if the first residue
///  of chain 2 in the Conformation is 27, then the renumbering of the chain in
///  PDBInfo will start from 27. (default true)
/// @param[in] keep_insertion_codes If true, will maintain insertion codes and
///  will not increment the pdb residue numbering for those residues.  This means
///  new numbering with insertion codes will only reflect properly if the
///  old numbering included the base numbering of the insertion code residues,
///  i.e. 100 100A 100B and not just 100A 100B (with 100 never appearing).
///  (default false)
/// @param[in] rotate_chain_ids If true, allows support for more than 26 pdb chains
///  by rotating [A,Z] continuously.  WARNING: This will break the assumption
///  made by the PDBPoseMap that each pdb chain id is unique, so make sure you
///  are not using the PDBPoseMap feature downstream in your code path without
///  corrections! (default false)
/// @remarks If fixing chains and there is only one chain and the PDBInfo exists
///  but all records are marked as empty, will renumber and set the PDBInfo chain
///  to 'A'.
/// @return true if renumbering successful, false otherwise
bool renumber_pdbinfo_based_on_conf_chains(
	core::pose::Pose & pose,
	bool fix_chains,
	bool const start_from_existing_numbering,
	bool const keep_insertion_codes,
	bool const rotate_chain_ids
)
{
	using core::Size;
	using core::pose::PDBInfo;
	typedef std::map< int, char > Conf2PDB;

	if ( !pose.pdb_info().get() ) {
		TR.Warning << "WARNING: renumber_pdbinfo_based_on_conf_chains(): no PDBInfo, returning" << std::endl;
		return false;
	}

	Conf2PDB conf2pdb = conf2pdb_chain( pose );

	if ( fix_chains ) {
		if ( conf2pdb.empty() ) { // something is wrong with chain consistency
			TR.Warning << "WARNING: renumber_pdbinfo_based_on_conf_chains(): Request to fix PDBInfo chains, but ";
			TR.Warning << "chain mapping is inconsistent, so that step will be skipped." << std::endl;
			fix_chains = false;
		} else { // Try to fill in any empty record characters.

			// two different schemes: rotating and fixed length
			// WARNING: Rotating will break assumption of unique chain ids
			// inside PDBPoseMap, so make sure you are not using the PDBPoseMap
			// feature after calling this function without correcting!
			// First either remove or rotate any existing chains to the end of
			// the list.
			std::string letters( utility::UPPERCASE_LETTERS );
			for ( auto & i : conf2pdb ) {
				if ( i.second == PDBInfo::empty_record() )  continue;

				std::string::size_type const j = letters.find( i.second );
				if ( j == std::string::npos )  continue;

				if ( rotate_chain_ids ) { // rotating
					letters.push_back( letters.at( j ) );
				}

				letters.erase( j, 1 );
			}

			// Now fill in empty records.
			Size lidx = 0;
			for ( auto & i : conf2pdb ) {
				if ( i.second != PDBInfo::empty_record() )  continue;

				if ( rotate_chain_ids ) { // rotating
					i.second = letters.at( lidx % letters.size() );
				} else { // fixed length
					runtime_assert( lidx < letters.size() );
					i.second = letters.at( lidx );
				}
				++lidx;
			}

		} // if conf2pdb.empty()
	} // if fix_chains

	PDBInfo & pdbinfo = *pose.pdb_info();

	// grab all the chain endings
	utility::vector1< Size > chain_endings = pose.conformation().chain_endings();
	chain_endings.push_back( pose.size() ); // add the last res, which is not in the list

	Size res = 1;
	for ( utility::vector1< Size >::const_iterator i = chain_endings.begin(), ie = chain_endings.end(); i != ie; ++i ) {
		Size const chain_end = *i;
		int pdb_res = 0; // new chain, so reset pdb_res counter
		if ( start_from_existing_numbering && pdbinfo.chain( res ) != PDBInfo::empty_record() ) {
			pdb_res = pdbinfo.number( res ) - 1;
		}

		// find the corresponding pdb chain
		Conf2PDB::const_iterator c2p = conf2pdb.find( pose.chain( chain_end ) );
		debug_assert( ( fix_chains && c2p != conf2pdb.end() ) || !fix_chains ); // otherwise something's very wrong

		for ( ; res <= chain_end; ++res ) {
			// handle the pdb chain only if necessary
			char chain = pdbinfo.chain( res );
			if ( pdbinfo.chain( res ) == PDBInfo::empty_record() && fix_chains && c2p != conf2pdb.end() ) {
				chain = c2p->second;
			}

			// If keeping insertion codes, increment pdb_res counter only if
			// no insertion code or we're at position 1, in case there's an
			// insertion code at 1.
			char icode = pdbinfo.icode( res );
			if ( keep_insertion_codes && ( pdbinfo.icode( res ) == ' ' || res == 1 ) ) {
				++pdb_res;
			} else if ( !keep_insertion_codes ) { // always increment and clear insertion code
				icode = ' ';
				++pdb_res;
			}

			// The new pdb info for this residue must be setup in one shot.
			// The way we're redoing the info in this function can cause the
			// pdb2pose map to become out-of-sync if we attempt to make the
			// changes by separately calling chain(), icode(), and number().
			pdbinfo.set_resinfo( res, chain, pdb_res, icode );
		}

	} // foreach chain ending

	// no point updating pdb_info if it's just thrown away
	pose.pdb_info()->obsolete( false );

	debug_assert( res == pose.size() + 1 );

	return true;
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
)
{
	return conformation::is_ideal_position(
		seqpos, pose.conformation() );
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
		throw utility::excn::EXCN_Msg_Exception("Pose utility remove_ligand_canonical_residues: I have received a pose with only one residue but cannot delete the last residue of the pose.");
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
	core::io::silent::BinarySilentStruct lhs_silent_struct(lhs, "" );
	std::stringstream lhs_str;
	lhs_silent_struct.print_conformation(lhs_str);

	core::io::silent::BinarySilentStruct rhs_silent_struct(rhs, "" );
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
				// tr.Error << "Error: can't find atom " << named_atom_id.atom() << " in residue "
				//   << rt.name() << ", residue has " << rt.nheavyatoms() << " heavy atoms." << std::endl;
				//  tr.Error << "atom names are: " << std::endl;
				//rt.show_all_atom_names( tr.Error );
				if ( raise_exception ) throw id::EXCN_AtomNotFound( named_atom_id );
				return id::BOGUS_ATOM_ID;
			}
		} else {
			// tr.Error << "Error: can't find residue " << named_atom_id.rsd()
			//  << " in pose (pose.size() = ) "
			//  << pose.size() << std::endl;
			if ( raise_exception ) throw id::EXCN_AtomNotFound( named_atom_id );
			return id::BOGUS_ATOM_ID;
		}
	} else {
		if ( raise_exception ) throw id::EXCN_AtomNotFound( named_atom_id );
		return id::BOGUS_ATOM_ID;
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
std::string tag_from_pose( core::pose::Pose const & pose ) {
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;

	std::string tag( "empty_tag" );
	if ( pose.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
		tag =
			static_cast< basic::datacache::CacheableString const & >
			( pose.data().get( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
	}

	return tag;
}

void tag_into_pose( core::pose::Pose & pose, std::string const & tag ) {
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;
	using basic::datacache::CacheableString;
	using basic::datacache::DataCache_CacheableData;
	pose.data().set(
		core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG,
		DataCache_CacheableData::DataOP( new CacheableString(tag) )
	);
}

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
)
{

	core::conformation::Residue const & old_rsd( pose.residue( seqpos ) );
	core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( new_rsd_type ) );
	conformation::copy_residue_coordinates_and_rebuild_missing_atoms( old_rsd, *new_rsd, pose.conformation() );
	pose.replace_residue( seqpos, *new_rsd, false );

}


core::conformation::ResidueOP
remove_variant_type_from_residue(
	core::conformation::Residue const & old_rsd,
	core::chemical::VariantType const variant_type,
	pose::Pose const & pose )
{
	if ( !old_rsd.has_variant_type( variant_type ) ) return old_rsd.clone();

	// the type of the desired variant residue
	core::chemical::ResidueTypeSetCOP rsd_set( old_rsd.residue_type_set() );
	core::chemical::ResidueType const & new_rsd_type( rsd_set->get_residue_type_with_variant_removed( old_rsd.type(), variant_type ) );
	core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( new_rsd_type, old_rsd, pose.conformation() ) );
	core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( old_rsd, *new_rsd, pose.conformation() );
	if ( old_rsd.nchi() == new_rsd_type.nchi() ) {
		for ( Size chino=1; chino <= old_rsd.nchi(); ++chino ) {
			new_rsd->set_chi( chino, old_rsd.chi( chino ) );
		}
	} else {
		TR << "The chi angles will not be updated and your dunbrack score for this rotamer will be huge; "
			"this function is only meant to add a variant type to a residue of the same type" << std::endl;
	}

	return new_rsd;
}


conformation::ResidueOP
add_variant_type_to_residue(
	conformation::Residue const & old_rsd,
	chemical::VariantType const variant_type,
	pose::Pose const & pose )
{

	if ( old_rsd.has_variant_type( variant_type ) ) return old_rsd.clone();

	// the type of the desired variant residue
	chemical::ResidueTypeSetCOP rsd_set( old_rsd.residue_type_set() );
	chemical::ResidueType const & new_rsd_type( rsd_set->get_residue_type_with_variant_added( old_rsd.type(), variant_type ) );
	conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( new_rsd_type, old_rsd, pose.conformation() ) );
	conformation::copy_residue_coordinates_and_rebuild_missing_atoms( old_rsd, *new_rsd, pose.conformation() );
	if ( old_rsd.nchi() == new_rsd_type.nchi() ) {
		for ( Size chino=1; chino <= old_rsd.nchi(); ++chino ) {
			new_rsd->set_chi( chino, old_rsd.chi( chino ) );
		}
	} else {
		TR << "The chi angles will not be updated and your dunbrack score for this rotamer will be huge; "
			"this function is only meant to add a variant type to a residue of the same type" << std::endl;
	}

	return new_rsd;
}


/// @details E.g., make a terminus variant, and replace the original in pose.
/// @note This copies any atoms in common between old and new residues, rebuilding the others.
void
add_variant_type_to_pose_residue(
	pose::Pose & pose,
	chemical::VariantType const variant_type,
	Size const seqpos )
{
	runtime_assert( seqpos != 0 );
	if ( pose.residue( seqpos ).has_variant_type( variant_type ) ) return;

	conformation::Residue const & old_rsd( pose.residue( seqpos ) );

	// the type of the desired variant residue
	chemical::ResidueTypeSetCOP rsd_set( old_rsd.residue_type_set() );
	chemical::ResidueType const & new_rsd_type( rsd_set->get_residue_type_with_variant_added( old_rsd.type(), variant_type ) );

	core::pose::replace_pose_residue_copying_existing_coordinates( pose, seqpos, new_rsd_type );

	// update connections
	for ( Size i_con=1; i_con<=pose.conformation().residue_type(seqpos).n_possible_residue_connections(); ++i_con ) {
		if ( pose.conformation().residue(seqpos).connected_residue_at_resconn(i_con) != 0 ) {
			Size connected_seqpos = pose.conformation().residue(seqpos).connected_residue_at_resconn(i_con);
			Size connected_id = pose.residue(seqpos).connect_map(i_con).connid();
			pose.conformation().update_noncanonical_connection(seqpos, i_con, connected_seqpos, connected_id);
		}
	}
}


/// @details E.g., remove a terminus variant, and replace the original in pose.
/// @note This copies any atoms in common between old and new residues, rebuilding the others.
void
remove_variant_type_from_pose_residue(
	pose::Pose & pose,
	chemical::VariantType const variant_type,
	Size const seqpos )
{
	if ( !pose.residue( seqpos ).has_variant_type( variant_type ) ) return;

	conformation::Residue const & old_rsd( pose.residue( seqpos ) );

	// the type of the desired variant residue
	chemical::ResidueTypeSetCOP rsd_set( old_rsd.residue_type_set() );
	chemical::ResidueType const & new_rsd_type( rsd_set->get_residue_type_with_variant_removed( old_rsd.type(), variant_type ) );

	core::pose::replace_pose_residue_copying_existing_coordinates( pose, seqpos, new_rsd_type );

	// update connections
	for ( Size i_con=1; i_con<=pose.conformation().residue_type(seqpos).n_possible_residue_connections(); ++i_con ) {
		if ( pose.conformation().residue(seqpos).connected_residue_at_resconn(i_con) != 0 ) {
			Size connected_seqpos = pose.conformation().residue(seqpos).connected_residue_at_resconn(i_con);
			Size connected_id = pose.residue(seqpos).connect_map(i_con).connid();
			pose.conformation().update_noncanonical_connection(seqpos, i_con, connected_seqpos, connected_id);
		}
	}
}

void
add_lower_terminus_type_to_pose_residue(
	pose::Pose & pose,
	Size const seqpos
)
{
	add_variant_type_to_pose_residue( pose, chemical::LOWER_TERMINUS_VARIANT, seqpos );
}

void
add_upper_terminus_type_to_pose_residue(
	pose::Pose & pose,
	Size const seqpos
)
{
	add_variant_type_to_pose_residue( pose, chemical::UPPER_TERMINUS_VARIANT, seqpos );
}

void
remove_lower_terminus_type_from_pose_residue(
	pose::Pose & pose,
	Size const seqpos
)
{
	core::pose::remove_variant_type_from_pose_residue( pose, chemical::LOWER_TERMINUS_VARIANT, seqpos );
}

void
remove_upper_terminus_type_from_pose_residue(
	pose::Pose & pose,
	Size const seqpos
)
{
	core::pose::remove_variant_type_from_pose_residue( pose, chemical::UPPER_TERMINUS_VARIANT, seqpos );
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
	core::pose::initialize_dof_id_map( dof_map, pose, id::BOGUS_TORSION_ID );

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
				if ( dof_map[ id  ] != BOGUS_TORSION_ID ) {
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
				TR.Warning << "WARNING: Unable to find atom_tree atom for this " <<
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
				TR.Warning << "WARNING: Unable to find atom_tree atom for this " <<
					"Rosetta nu angle: residue " << i << " NU " << j << std::endl;
			}
		} // j=1, n_nu_torsions

		// finally, any branch connection torsions
		Size const n_branch_torsions( rsd.n_non_polymeric_residue_connections() );
		for ( uint j( 1 ); j <= n_branch_torsions; ++j ) {
			if ( j == 1 && rsd.is_branch_lower_terminus() ) { continue; }  // Only an "outgoing" connection matters.
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
				TR.Warning << "WARNING: Unable to find atom_tree atom for this " <<
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


bool
has_chain(std::string const & chain, core::pose::Pose const & pose){
	debug_assert(chain.size()==1);// chain is one char
	char chain_char= chain[0];
	return has_chain(chain_char, pose);
}

bool
has_chain(char const & chain, core::pose::Pose const & pose){
	for ( core::Size i=1; i <= pose.conformation().num_chains(); ++i ) {
		char this_char= get_chain_from_chain_id(i, pose);
		if ( this_char == chain ) {
			return true;
		}
	}
	return false;
}

bool
has_chain(core::Size chain_id, core::pose::Pose const & pose){
	for ( core::Size i=1; i <= pose.conformation().num_chains(); ++i ) {
		if ( i == chain_id ) {
			return true;
		}
	}
	return false;
}

std::set<core::Size>
get_jump_ids_from_chain_ids(std::set<core::Size> const & chain_ids, core::pose::Pose const & pose){
	std::set<core::Size> jump_ids;
	auto chain_id= chain_ids.begin();
	for ( ; chain_id != chain_ids.end(); ++chain_id ) {
		core::Size jump_id= get_jump_id_from_chain_id(*chain_id, pose);
		jump_ids.insert(jump_id);
	}
	return jump_ids;
}

core::Size
get_jump_id_from_chain_id(core::Size const & chain_id,const core::pose::Pose & pose){
	for ( core::Size jump_id=1; jump_id <= pose.num_jump(); jump_id++ ) {
		core::Size ligand_residue_id= (core::Size) pose.fold_tree().downstream_jump_residue(jump_id);
		core::Size ligand_chain_id= pose.chain(ligand_residue_id);
		if ( chain_id==ligand_chain_id ) {
			return jump_id;
		}
	}
	utility_exit();
	return 0;// this will never happen
}

core::Size
get_chain_id_from_chain(std::string const & chain, core::pose::Pose const & pose){
	debug_assert(chain.size()==1);// chain is one char
	if ( chain.size() > 1 ) utility_exit_with_message("Multiple chain_ids per chain! Are you using '-treat_residues_in_these_chains_as_separate_chemical_entities', and not using compatible movers?" );
	char chain_char= chain[0];
	return get_chain_id_from_chain(chain_char, pose);
}

core::Size
get_chain_id_from_chain(char const & chain, core::pose::Pose const & pose){
	utility::vector1<core::Size> chain_ids = get_chain_ids_from_chain(chain, pose);
	if ( chain_ids.size() == 0 ) {
		throw utility::excn::EXCN_RangeError(" chain_id "+utility::to_string(chain)+" does not exist");
	} else if ( chain_ids.size() > 1 ) {
		throw utility::excn::EXCN_RangeError("chain_id "+utility::to_string(chain)+" represents more than one chain!");
	}
	return chain_ids[1];
}

utility::vector1<core::Size>
get_chain_ids_from_chain(std::string const & chain, core::pose::Pose const & pose){
	debug_assert(chain.size()==1);// chain is one char
	char chain_char= chain[0];
	return get_chain_ids_from_chain(chain_char, pose);

}

utility::vector1<core::Size>
get_chain_ids_from_chain(char const & chain, core::pose::Pose const & pose){
	utility::vector1<core::Size> chain_ids;
	for ( core::Size i=1; i <= pose.conformation().num_chains(); i++ ) {
		char this_char= get_chain_from_chain_id(i, pose);
		if ( this_char == chain ) {
			chain_ids.push_back(i);
		}
	}
	return chain_ids;
}

utility::vector1<core::Size>
get_chain_ids_from_chains(utility::vector1<std::string> const & chains, core::pose::Pose const & pose){
	utility::vector1<core::Size> chain_ids;
	for ( Size i = 1; i <= chains.size(); ++i ) {
		Size chain_id = get_chain_id_from_chain(chains[ i ], pose);
		chain_ids.push_back(chain_id);
	}
	return chain_ids;
}

core::Size
get_jump_id_from_chain(std::string const & chain, core::pose::Pose const & pose){
	core::Size chain_id= get_chain_id_from_chain(chain, pose);
	return get_jump_id_from_chain_id(chain_id, pose);
}

core::Size
get_jump_id_from_chain(char const & chain, core::pose::Pose const & pose){
	core::Size chain_id= get_chain_id_from_chain(chain, pose);
	return get_jump_id_from_chain_id(chain_id, pose);
}

utility::vector1<core::Size>
get_jump_ids_from_chain(char const & chain, core::pose::Pose const & pose){
	utility::vector1<core::Size> jump_ids;
	utility::vector1<core::Size> chain_ids = get_chain_ids_from_chain(chain, pose);
	for ( core::Size const chain_id : chain_ids ) {
		jump_ids.push_back( get_jump_id_from_chain_id(chain_id, pose));
	}
	return jump_ids;
}

utility::vector1<core::Size>
get_jump_ids_from_chain(std::string const & chain, core::pose::Pose const & pose){
	debug_assert(chain.size()==1);// chain is one char
	char chain_char= chain[0];
	return get_jump_ids_from_chain(chain_char, pose);
}

core::Size get_chain_id_from_jump_id(core::Size const & jump_id, core::pose::Pose const & pose){
	core::Size ligand_residue_id= (core::Size) pose.fold_tree().downstream_jump_residue(jump_id);
	return pose.chain(ligand_residue_id);
}

char
get_chain_from_jump_id(core::Size const & jump_id, core::pose::Pose const & pose){
	core::Size chain_id= get_chain_id_from_jump_id(jump_id, pose);
	return get_chain_from_chain_id(chain_id, pose);
}

core::conformation::ResidueCOPs
get_chain_residues(core::pose::Pose const & pose, core::Size const chain_id){
	core::Size begin= pose.conformation().chain_begin(chain_id);
	core::Size const end= pose.conformation().chain_end(chain_id);
	core::conformation::ResidueCOPs residues;
	for ( ; begin <= end; ++begin ) {
		residues.push_back( core::conformation::ResidueOP( new core::conformation::Residue(pose.residue(begin)) ) );
	}
	return residues;
}

/// @brief Is residue number in this chain?
bool res_in_chain( core::pose::Pose const & pose, core::Size resnum, std::string chain ) {

	bool in_chain( false );

	Size chain_of_res = pose.chain( resnum );
	Size chainid = get_chain_id_from_chain( chain[0], pose );

	if ( chain_of_res == chainid ) {
		in_chain = true;
	}

	return in_chain;
}

char
get_chain_from_chain_id(core::Size const & chain_id, core::pose::Pose const & pose){
	core::Size first_chaisize= pose.conformation().chain_begin( chain_id );
	return pose.pdb_info()->chain(first_chaisize);
}

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


core::Size get_hash_from_chain(char const & chain, core::pose::Pose const & pose)
{
	std::size_t hash = 0;

	core::Size chain_id = get_chain_id_from_chain(chain,pose);
	core::Size chain_begin = pose.conformation().chain_begin(chain_id);
	core::Size chain_end = pose.conformation().chain_end(chain_id);
	for ( core::Size res_num = chain_begin; res_num <= chain_end; ++res_num ) {
		core::Size natoms = pose.conformation().residue(res_num).natoms();
		for ( core::Size atom_num = 1; atom_num <= natoms; ++atom_num ) {
			id::AtomID atom_id(atom_num,res_num);
			PointPosition current_xyz = pose.conformation().xyz(atom_id);
			boost::hash_combine(hash,current_xyz);
		}
	}

	return hash;
}

core::Size get_hash_excluding_chain(char const & chain, core::pose::Pose const & pose)
{
	std::size_t hash = 0;

	core::Size chain_id = get_chain_id_from_chain(chain,pose);

	for ( core::Size res_num = 1; res_num <= pose.size(); ++res_num ) {
		if ( (int)chain_id == pose.chain(res_num) ) {
			continue;
		}
		core::Size natoms = pose.conformation().residue(res_num).natoms();
		for ( core::Size atom_num = 1; atom_num <= natoms; ++atom_num ) {
			id::AtomID atom_id(atom_num,res_num);
			PointPosition current_xyz = pose.conformation().xyz(atom_id);
			boost::hash_combine(hash,current_xyz);
		}
	}

	return hash;
}

std::string get_sha1_hash_excluding_chain(char const & chain, core::pose::Pose const & pose)
{

	std::stringstream coord_stream;

	core::Size chain_id = get_chain_id_from_chain(chain,pose);

	for ( core::Size res_num = 1; res_num <= pose.size(); ++res_num ) {
		if ( (int)chain_id == pose.chain(res_num) ) {
			continue;
		}

		core::Size natoms = pose.conformation().residue(res_num).natoms();
		for ( core::Size atom_num = 1; atom_num <= natoms; ++atom_num ) {
			id::AtomID atom_id(atom_num,res_num);
			PointPosition current_xyz = pose.conformation().xyz(atom_id);
			coord_stream << numeric::truncate_and_serialize_xyz_vector(current_xyz,3);
		}
	}
	return utility::string_to_sha1(coord_stream.str());
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
	} else {
		if ( option[ in::detect_disulf ].user() ?
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
						TR.Error << "Error: SSBond records list multiple nonredundant disulfides for this residue!" << std::endl;
						utility_exit();
					}
				}

				Size seqpos_one = pose.pdb_info()->pdb2pose( ssbond.second[1].chainID1, ssbond.second[1].resSeq1, ssbond.second[1].iCode1 );
				Size seqpos_two = pose.pdb_info()->pdb2pose( ssbond.second[1].chainID2, ssbond.second[1].resSeq2, ssbond.second[1].iCode2 );

				if ( seqpos_one != 0 && seqpos_two != 0 ) {
					disulf_one.push_back( seqpos_one );
					disulf_two.push_back( seqpos_two );
				}
			}

			//pose.conformation().detect_disulfides( disulfs );
			pose.conformation().detect_disulfides( disulf_one, disulf_two );
		}
	}
}

std::string extract_tag_from_pose( core::pose::Pose &pose )
{
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;
	using basic::datacache::CacheableString;
	using basic::datacache::CacheableStringOP;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
		CacheableStringOP data =  utility::pointer::dynamic_pointer_cast< CacheableString > (  (pose.data().get_ptr( ( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG  )) ));
		if ( data.get() == nullptr ) return std::string("UnknownTag");
		else                      return data->str();
	}

	return std::string("UnknownTag");
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
	core::Size count = 0;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		core::conformation::Residue const & resi(pose.residue(i));
		if ( resi.aa() <= core::chemical::num_canonical_aas ) {
			++count;
		}
	}
	return count;
}

core::Size noncanonical_residue_count(core::pose::Pose const & pose)
{
	core::Size count = 0;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		core::conformation::Residue const & resi(pose.residue(i));
		if ( resi.aa() > core::chemical::num_canonical_aas ) {
			++count;
		}
	}
	return count;
}

core::Size canonical_atom_count(core::pose::Pose const & pose)
{
	core::Size count = 0;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		core::conformation::Residue const & resi(pose.residue(i));
		if ( resi.aa() <= core::chemical::num_canonical_aas ) {
			count += resi.natoms();
		}
	}
	return count;
}

core::Size noncanonical_atom_count(core::pose::Pose const & pose)
{
	core::Size count = 0;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		core::conformation::Residue const & resi(pose.residue(i));
		if ( resi.aa() > core::chemical::num_canonical_aas ) {
			count += resi.natoms();
		}
	}
	return count;
}

core::Size noncanonical_chi_count(core::pose::Pose const & pose)
{
	core::Size count = 0;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		core::conformation::Residue const & resi(pose.residue(i));
		if ( resi.aa() > core::chemical::num_canonical_aas ) {
			count += resi.nchi();
		}
	}
	return count;
}

/// @brief Number of protein residues in the pose
/// @details No virtuals, membrane residues or embedding residues counted
Size nres_protein( pose::Pose const & pose ) {

	Size cnt(0);

	// go over pose residues and ask whether is_protein
	for ( Size i = 1; i <= pose.size(); ++i ) {

		if ( pose.residue(i).is_protein() ) {
			cnt++;
		}
	}

	return cnt;
}// nres_protein


numeric::xyzVector< Real >
center_of_mass(
	pose::Pose const & pose,
	utility::vector1< bool > const & residues
)
{
	using namespace numeric;
	using core::conformation::Residue;
	assert( pose.size() == residues.size() );

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
)
{
	utility::vector1< bool > residues( pose.size(), false );
	assert( (Size) stop <= residues.size() );
	assert( stop > start && start > 0 );

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
)
{
	return center_of_mass( pose, generate_vector_from_bounds( pose, start, stop ) );
}

int
residue_center_of_mass(
	pose::Pose const & pose,
	utility::vector1< bool > residues
)
{
	assert( pose.size() == residues.size() );

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
)
{
	Vector center = center_of_mass(pose, start, stop );
	return core::pose::return_nearest_residue( pose, start, stop, center );
}

int
return_nearest_residue(
	pose::Pose const & pose,
	utility::vector1< bool > const & residues,
	Vector center
)
{
	using core::conformation::Residue;

	assert( pose.size() == residues.size() );

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
)
{
	return return_nearest_residue( pose, generate_vector_from_bounds( pose, begin, end ), center );
}

// silly conversion from std::map< AtomID, AtomID> to rosetta's silly AtomID_Map class.
id::AtomID_Map< id::AtomID >
convert_from_std_map( std::map< id::AtomID, id::AtomID > const & atom_map,
	core::pose::Pose const & pose ){
	id::AtomID_Map< id::AtomID > atom_ID_map;
	initialize_atomid_map( atom_ID_map, pose, id::BOGUS_ATOM_ID );
	for ( auto const & it : atom_map ) {
		atom_ID_map.set( it.first, it.second );
	}
	return atom_ID_map;
}

/// @brief Create std::map from PDBPoseMap in pose (JKLeman)
std::map< std::string, core::Size > get_pdb2pose_numbering_as_stdmap ( core::pose::Pose const & pose ) {

	using namespace utility;

	// initialize empty map
	std::map< std::string, core::Size > pdb2pose_map;

	// go through protein, get chain, resn, insertion code
	for ( core::Size i = 1; i <= nres_protein( pose ); ++i ) {

		// get chain, PDB resnumber, and insertion code for that residue
		std::string chain = to_string( pose.pdb_info()->chain( i ) );
		std::string resn = to_string( pose.pdb_info()->number( i ) );
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

// @brief chemical bond from lower to upper residue across CUTPOINT_LOWER/CUTPOINT_UPPER will prevent steric repulsion.
void
declare_cutpoint_chemical_bond( core::pose::Pose & pose, Size const cutpoint_res, Size const next_res_in /* = 0 */ )
{
	Size const next_res = ( next_res_in == 0 ) ? ( cutpoint_res + 1 ) : next_res_in; // user might specify a different "next_res" to cyclize.
	pose.conformation().declare_chemical_bond( cutpoint_res,
		pose.residue( cutpoint_res ).atom_name( pose.residue( cutpoint_res ).upper_connect_atom() ),
		next_res,
		pose.residue( next_res ).atom_name( pose.residue( next_res ).lower_connect_atom() ) );
}

/// @brief Add cutpoint variants to all residues annotated as cutpoints in the pose.
void
correctly_add_cutpoint_variants( core::pose::Pose & pose ) {
	const core::kinematics::FoldTree & tree(pose.fold_tree());
	for ( core::Size i = 1; i < pose.size(); ++i ) { // Less than because cutpoints are between i and i+1
		if ( tree.is_cutpoint(i) ) {
			correctly_add_cutpoint_variants( pose, i, false );
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// try to unify all cutpoint addition into this function.
void
correctly_add_cutpoint_variants( core::pose::Pose & pose,
	Size const cutpoint_res,
	bool const check_fold_tree /* = true*/,
	Size const next_res_in /* = 0 */ )
{
	using namespace core::chemical;

	Size const next_res = ( next_res_in == 0 ) ? ( cutpoint_res + 1 ) : next_res_in; // user might specify a different "next_res" to cyclize.

	if ( check_fold_tree ) runtime_assert( pose.fold_tree().is_cutpoint( cutpoint_res ) );

	remove_variant_type_from_pose_residue( pose, UPPER_TERMINUS_VARIANT, cutpoint_res );
	remove_variant_type_from_pose_residue( pose, THREE_PRIME_PHOSPHATE, cutpoint_res );
	remove_variant_type_from_pose_residue( pose, C_METHYLAMIDATION, cutpoint_res );

	remove_variant_type_from_pose_residue( pose, LOWER_TERMINUS_VARIANT, next_res );
	remove_variant_type_from_pose_residue( pose, VIRTUAL_PHOSPHATE, next_res );
	remove_variant_type_from_pose_residue( pose, FIVE_PRIME_PHOSPHATE, next_res );
	remove_variant_type_from_pose_residue( pose, N_ACETYLATION, next_res);

	if ( pose.residue_type( cutpoint_res ).is_RNA() )  rna::position_cutpoint_phosphate_torsions( pose, cutpoint_res, next_res );

	add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, cutpoint_res   );
	add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, next_res );

	// important -- to prevent artificial penalty from steric clash.
	declare_cutpoint_chemical_bond( pose, cutpoint_res, next_res );
}

void
get_constraints_from_link_records( core::pose::Pose & pose, io::StructFileRep const & sfr )
{
	using namespace scoring::func;

	/*HarmonicFuncOP amide_harm_func    ( new HarmonicFunc( 1.34, 0.05 ) );
	HarmonicFuncOP thioester_harm_func( new HarmonicFunc( 1.83, 0.1 ) );
	*/
	CircularHarmonicFuncOP ang_func( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi_2_over_3(), 0.02 ) );
	CircularHarmonicFuncOP ang90_func( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi_over_2(), 0.02 ) );
	CircularHarmonicFuncOP dih_func( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi(), 0.02 ) );

	for ( auto const & it : sfr.link_map() ) {
		for ( Size ii = 1; ii <= it.second.size(); ++ii ) {
			TR << "|"<<it.second[ii].chainID1 << "| |" << it.second[ii].resSeq1 << "|" << std::endl;
			TR << "|"<<it.second[ii].chainID2 << "| |" << it.second[ii].resSeq2 << "|" << std::endl;

			Size id1 = pose.pdb_info()->pdb2pose( it.second[ii].chainID1, it.second[ii].resSeq1 );
			Size id2 = pose.pdb_info()->pdb2pose( it.second[ii].chainID2, it.second[ii].resSeq2 );
			conformation::Residue const & NUC = pose.residue( id1 );
			conformation::Residue const & ELEC = pose.residue( id2 );

			id::AtomID aidNUC( NUC.atom_index( it.second[ii].name1 ), id1 );
			id::AtomID aidC( ELEC.atom_index( it.second[ii].name2 ), id2 );

			scoring::func::HarmonicFuncOP harm_func( new scoring::func::HarmonicFunc( it.second[ii].length, 0.05 ) );

			scoring::constraints::ConstraintCOP atompair(
				new scoring::constraints::AtomPairConstraint( aidNUC, aidC, harm_func ) );
			pose.add_constraint( atompair );
			TR << "Adding harmonic constraint between residue " << id1 << " atom " << it.second[ii].name1;
			TR << "and residue " << id2 << " atom " << it.second[ii].name2 << " with length " << it.second[ii].length << std::endl;

			// Cover cyclization case first
			if ( it.second[ii].name1 == " N  " && it.second[ii].name2 == " C  " ) {

				id::AtomID aidCA( ELEC.atom_index( "CA" ), id2 );
				id::AtomID aidO( ELEC.atom_index( "O" ), id2 );
				id::AtomID aidCA2( NUC.atom_index( "CA" ), id1 );

				scoring::constraints::ConstraintCOP ang( new scoring::constraints::AngleConstraint( aidNUC, aidC, aidCA, ang_func ) );
				pose.add_constraint( ang );
				scoring::constraints::ConstraintCOP ang2( new scoring::constraints::AngleConstraint( aidCA2, aidNUC, aidC, ang_func ) );
				pose.add_constraint( ang2 );
				scoring::constraints::ConstraintCOP dih( new scoring::constraints::DihedralConstraint( aidNUC, aidC, aidO, aidCA, dih_func ) );
				pose.add_constraint( dih );

				// constrain omega unless N terminal residue is PRO or peptoid
				// I think we should have a residue property that describes PRO-type residue types
				if ( NUC.type().name3() == "PRO" || NUC.type().is_peptoid() ) {
					// let god sort it out
					// PERHAPS: constraints that mirror the MM torsions actually guiding this case!
				} else {
					id::AtomID aidH( NUC.atom_index( "H" ), id1 );

					scoring::constraints::ConstraintCOP omg( new scoring::constraints::DihedralConstraint( aidCA2, aidNUC, aidC, aidCA, dih_func ) );
					pose.add_constraint( omg );
					scoring::constraints::ConstraintCOP omgH( new scoring::constraints::DihedralConstraint( aidO, aidC, aidNUC, aidH, dih_func ) );
					pose.add_constraint( omgH );
					scoring::constraints::ConstraintCOP omgimp( new scoring::constraints::DihedralConstraint( aidC, aidNUC, aidH, aidCA2, dih_func ) );
					pose.add_constraint( omgimp );
				}

				TR << "Adding harmonic constraints to the angle formed by atoms N, C, O ( 120 ) and ";
				TR << "the improper torsion N, C, O, CA (180) and the dihedral CA, N, C, CA ( 180 ) " <<std::endl;

			} else if ( it.second[ii].name1 == " C  " && it.second[ii].name2 == " N  " ) {
				// swapped; nuc is actually elec and vice versa

				id::AtomID aidCA2( ELEC.atom_index( "CA" ), id2 );
				id::AtomID aidO( NUC.atom_index( "O" ), id1 );
				id::AtomID aidCA( NUC.atom_index( "CA" ), id1 );

				scoring::constraints::ConstraintCOP ang( new scoring::constraints::AngleConstraint( aidC, aidNUC, aidCA, ang_func ) );
				pose.add_constraint( ang );
				scoring::constraints::ConstraintCOP ang2( new scoring::constraints::AngleConstraint( aidCA2, aidC, aidNUC, ang_func ) );
				pose.add_constraint( ang2 );
				scoring::constraints::ConstraintCOP dih( new scoring::constraints::DihedralConstraint( aidC, aidNUC, aidO, aidCA, dih_func ) );
				pose.add_constraint( dih );

				// constrain omega unless N terminal residue is PRO or peptoid
				// I think we should have a residue property that describes PRO-type residue types
				if ( ELEC.type().name3() == "PRO" || ELEC.type().is_peptoid() ) {
					// let god sort it out
					// PERHAPS: constraints that mirror the MM torsions actually guiding this case!
				} else {
					id::AtomID aidH( ELEC.atom_index( "H" ), id1 );

					scoring::constraints::ConstraintCOP omg( new scoring::constraints::DihedralConstraint( aidCA2, aidC, aidNUC, aidCA, dih_func ) );
					pose.add_constraint( omg );
					scoring::constraints::ConstraintCOP omgH( new scoring::constraints::DihedralConstraint( aidO, aidNUC, aidC, aidH, dih_func ) );
					pose.add_constraint( omgH );
					scoring::constraints::ConstraintCOP omgimp( new scoring::constraints::DihedralConstraint( aidNUC, aidC, aidH, aidCA2, dih_func ) );
					pose.add_constraint( omgimp );
				}
			}

			// Now let's add constraints based on residue identities and atom names. For example, let's cover
			if ( it.second[ii].name2 == " CZ " && it.second[ii].resName2 == "VDP" ) {
				// thiol-ene conjugation to acryl residue
				// (don't check name1 because we don't care SG/SD/SG1)
				// someday we will be fancy and check vs type.get_disulfide_atom_name()
				id::AtomID aidCB( NUC.atom_index( "CB" ), id2 );
				id::AtomID aidCE2( ELEC.atom_index( "CE2" ), id2 );
				id::AtomID aidO( ELEC.atom_index( "O" ), id2 );
				scoring::constraints::ConstraintCOP ang(
					new scoring::constraints::AngleConstraint( aidCB, aidNUC, aidC, ang90_func ) );
				pose.add_constraint( ang );
				//scoring::constraints::ConstraintCOP dih(
				//  new scoring::constraints::DihedralConstraint( aidCB, aidNUC, aidC, aidCE2, dih_func ) );
				//pose.add_constraint( dih );

				TR << "Assuming thiol-ene, adding harmonic constraints to the angle formed by CB, SG, CZ ( 90 )" << std::endl;// and ";
				//TR << "the dihedral CB, SG, CZ, CE2 ( 180 ) " << std::endl;

			} else if ( it.second[ii].name2 == " C  " ) {
				// the C-terminal conjugation case:
				id::AtomID aidCA( ELEC.atom_index( "CA" ), id2 );
				id::AtomID aidO( ELEC.atom_index( "O" ), id2 );
				scoring::constraints::ConstraintCOP ang(
					new scoring::constraints::AngleConstraint( aidNUC, aidC, aidO, ang_func ) );
				pose.add_constraint( ang );
				scoring::constraints::ConstraintCOP dih(
					new scoring::constraints::DihedralConstraint( aidNUC, aidC, aidO, aidCA, dih_func ) );
				pose.add_constraint( dih );
			} else if ( it.second[ii].name2 == " CD " ) {
				// The sidechain conjugation to glx case
				id::AtomID aidCA( ELEC.atom_index( "CG" ), id2 );
				id::AtomID aidO( ELEC.atom_index( "OE1" ), id2 );
				scoring::constraints::ConstraintCOP ang(
					new scoring::constraints::AngleConstraint( aidNUC, aidC, aidO, ang_func ) );
				pose.add_constraint( ang );
				scoring::constraints::ConstraintCOP dih(
					new scoring::constraints::DihedralConstraint( aidNUC, aidC, aidO, aidCA, dih_func ) );
				pose.add_constraint( dih );
			} else if ( it.second[ii].name2 == " CG " ) {
				// The sidechain conjugation to asx case
				id::AtomID aidCA( ELEC.atom_index( "CB" ), id2 );
				id::AtomID aidO( ELEC.atom_index( "OD1" ), id2 );
				scoring::constraints::ConstraintCOP ang(
					new scoring::constraints::AngleConstraint( aidNUC, aidC, aidO, ang_func ) );
				pose.add_constraint( ang );
				scoring::constraints::ConstraintCOP dih(
					new scoring::constraints::DihedralConstraint( aidNUC, aidC, aidO, aidCA, dih_func ) );
				pose.add_constraint( dih );
				if ( it.second[ii].name1 == " NE " ) {
					// ornithine
					id::AtomID aidH( NUC.atom_index( "1HE" ), id2 );
					id::AtomID aidCG( NUC.atom_index( "CG" ), id2 );

					/*

					id::AtomID aidH( NUC.atom_index( "H" ), id1 );

					scoring::constraints::ConstraintCOP omg( new scoring::constraints::DihedralConstraint( aidCA2, aidNUC, aidC, aidCA, dih_func ) );
					pose.add_constraint( omg );
					scoring::constraints::ConstraintCOP omgH( new scoring::constraints::DihedralConstraint( aidO, aidC, aidNUC, aidH, dih_func ) );
					pose.add_constraint( omgH );
					scoring::constraints::ConstraintCOP omgimp( new scoring::constraints::DihedralConstraint( aidC, aidNUC, aidH, aidCA2, dih_func ) );
					pose.add_constraint( omgimp );
					*/

					pose.add_constraint( scoring::constraints::ConstraintCOP( new scoring::constraints::DihedralConstraint( aidCA, aidO, aidC, aidNUC, dih_func ) ) );
					pose.add_constraint( scoring::constraints::ConstraintCOP( new scoring::constraints::DihedralConstraint( aidO, aidC, aidNUC, aidH, dih_func ) ) );
					pose.add_constraint( scoring::constraints::ConstraintCOP( new scoring::constraints::DihedralConstraint( aidC, aidNUC, aidH, aidCG, dih_func ) ) );
					pose.add_constraint( scoring::constraints::ConstraintCOP( new scoring::constraints::AngleConstraint( aidH, aidNUC, aidC, ang_func ) ) );
					pose.add_constraint( scoring::constraints::ConstraintCOP( new scoring::constraints::AngleConstraint( aidNUC, aidC, aidO, ang_func ) ) );

				} else if ( it.second[ii].name1 == " NZ " ) {
					// lysine
					id::AtomID aidH( NUC.atom_index( "1HZ" ), id2 );
					id::AtomID aidCD( NUC.atom_index( "CD" ), id2 );

					pose.add_constraint( scoring::constraints::ConstraintCOP( new scoring::constraints::DihedralConstraint( aidCA, aidO, aidC, aidNUC, dih_func ) ) );
					pose.add_constraint( scoring::constraints::ConstraintCOP( new scoring::constraints::DihedralConstraint( aidO, aidC, aidNUC, aidH, dih_func ) ) );
					pose.add_constraint( scoring::constraints::ConstraintCOP( new scoring::constraints::DihedralConstraint( aidC, aidNUC, aidH, aidCD, dih_func ) ) );
					pose.add_constraint( scoring::constraints::ConstraintCOP( new scoring::constraints::AngleConstraint( aidH, aidNUC, aidC, ang_func ) ) );
					pose.add_constraint( scoring::constraints::ConstraintCOP( new scoring::constraints::AngleConstraint( aidNUC, aidC, aidO, ang_func ) ) );

				}

			}

			//AtomID aidO( list[ ii ].second.id_O_, list[ ii ].second.resnum_ );
			//AtomID aidCA( list[ ii ].second.id_CA_, list[ ii ].second.resnum_ );
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
utility::vector1< Size > pdb_to_pose( pose::Pose const & pose, std::pair< utility::vector1< int >, utility::vector1<char> > const & pdb_res ){
	utility::vector1< Size > pose_res;
	runtime_assert( pdb_res.first.size() == pdb_res.second.size() );
	for ( Size n = 1; n <= pdb_res.first.size(); n++ ) {
		pose_res.push_back( pdb_to_pose( pose, pdb_res.first[ n ], pdb_res.second[ n ] ) );
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

/// @brief returns true if the given residue in the pose is a chain ending or has upper/lower terminal variants
bool
pose_residue_is_terminal( Pose const & pose, Size const resid )
{
	return ( is_lower_terminus( pose, resid ) || is_upper_terminus( pose, resid ) );
}

/// @brief checks to see if this is a lower chain ending more intelligently than just checking residue variants
bool
is_lower_terminus( pose::Pose const & pose, Size const resid )
{
	return ( ( resid == 1 ) || //loop starts at first residue
		( ! pose.residue(resid).is_polymer() ) || // this residue isn't a polymer
		( ! pose.residue( resid-1 ).is_protein() ) || //residue before start is not protein
		( pose.chain( resid-1 ) != pose.chain( resid ) ) || // residues before start are on another chain
		( pose.residue( resid ).is_lower_terminus() ) ); // start of residue is lower terminus
}

/// @brief checks to see if this is a lower chain ending more intelligently than just checking residue variants
bool
is_upper_terminus( pose::Pose const & pose, Size const resid )
{
	return ( ( resid == pose.size() ) || // loop end at last residue
		( !pose.residue( resid ).is_polymer() ) || // this residue isn't a polymer
		( !pose.residue( resid+1 ).is_protein() ) || // residue after end is not protein
		( pose.chain( resid+1 ) != pose.chain( resid ) ) || // residues before start is other chain
		( pose.residue( resid ).is_upper_terminus() ) ); // explicit terminus variant @ end of loop
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

	id::MainchainTorsionType torsion_type = static_cast< id::MainchainTorsionType >(torsion_id);

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

	id::MainchainTorsionType torsion_type = static_cast< id::MainchainTorsionType >( torsion_id );

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



} // pose
} // core
