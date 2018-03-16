// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/pose/carbohydrates/util.cc
/// @brief   Utility function definitions for carbohydrate-containing poses.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @author  Jared Adolf-Bryfogle <jadolfbr@gmail.com>

// Unit Headers
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/Pose.hh>

// Package Headers
#include <core/io/carbohydrates/pose_io.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/carbohydrates/database_io.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/select/residue_selector/OrResidueSelector.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/variant_util.hh>

// Project Headers
#include <core/id/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/types.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/carbohydrates/util.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/carbohydrates/GlycanTree.hh>
#include <core/conformation/carbohydrates/GlycanNode.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/residue_selector/ResidueSelector.hh>


// Utility Headers
#include <utility/string_constants.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

// Basic Headers
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/conversions.hh>
#include <numeric/angle.functions.hh>
#include <numeric/random/random.hh>

// External Headers
#include <boost/lexical_cast.hpp>
#include <basic/basic.hh>

// C++ Header
#include <list>
#include <utility>

// Construct tracer.
static basic::Tracer TR( "core.pose.carbohydrates.util" );


namespace core {
namespace pose {
namespace carbohydrates {
using std::map;
using namespace conformation::carbohydrates;

// Torsion Access /////////////////////////////////////////////////////////////
// Getters ////////////////////////////////////////////////////////////////////


///@brief Turn on/off IUPAC CHI for a particular residue number.
void
set_glycan_iupac_chi_torsions(
	core::pose::Pose const & pose,
	core::kinematics::MoveMap & movemap,
	core::Size resnum,
	bool action /*true*/)
{

	using namespace core::id;
	for ( core::Size chi_torsion_num = 1; chi_torsion_num <= pose.residue_type( resnum ).nchi(); ++ chi_torsion_num ) {
		TorsionID torsion_id = TorsionID(resnum, CHI, chi_torsion_num);
		core::Size child_resnum = get_downstream_residue_that_this_torsion_moves( pose, torsion_id);
		if ( child_resnum == 0 ) {
			TR.Debug << "Setting Non-BB CHI Torsion " << chi_torsion_num << std::endl;
			movemap.set( torsion_id, action);
			continue;
		}
	}

}

void
set_glycan_iupac_bb_torsions(
	core::pose::Pose const & pose,
	core::kinematics::MoveMap & movemap,
	core::Size resnum,
	bool action/*true*/)
{
	using namespace core::id;

	//Here is where it gets super tricky.

	// Setup MinMover-specific Glycan Movemap //
	// NOTE: THIS IS NOT IUPAC DEFINITION AND MUST BE DEALT WITH CAREFULLY.
	// 1) MM residue is the ONE BEFORE the IUPAC!
	// 2) CHI for ASN must be set to enable it's torsions to move.
	// 3) set_bb is for non-ring torsions.  Add an option to enable minimization of rings torsions (set_neus in movemap to true)
	// 4) set_branches must be on to allow the chemical-edge torsion to move.
	//    ->  This must be on for ASN and any residue that has a non-mainchain connection!
	core::Size i = resnum;
	core::Size parent_res = pose.glycan_tree_set()->get_parent( i );

	TR.Debug << "Setting up residue " << i << ", with parent " << parent_res << std::endl;
	//This residue is the protein-linkage. Set the BB up properly.
	bool branched_connection = false;
	if ( parent_res != 0 && !pose.residue( parent_res ).is_carbohydrate() ) {
		TR.Debug << "Setting ASN Linkage Torsion " << std::endl;
		movemap.set_chi( parent_res, action);
		movemap.set_branches( parent_res, action);

	} else {
		// Branch Linkage
		if ( parent_res != 0 && pose.residue(parent_res).connected_residue_at_upper() != i ) {
			core::Size total_connections = pose.residue( parent_res ).n_possible_residue_connections();
			core::Size mainchain_connections = pose.residue( parent_res).n_polymeric_residue_connections();
			for ( core::Size connection_id = 1; connection_id <= total_connections; ++connection_id ) {
				TR.Debug << "Connection: " << connection_id << std::endl;
				core::Size child = pose.residue( parent_res).residue_connection_partner( connection_id );
				core::SSize torsion_id = connection_id - mainchain_connections;

				//This is the branch.  We now turn on PHI for this connection.
				if ( child == i && torsion_id > 0 ) {
					TR.Debug << "Setting BRANCH torsion" << torsion_id << std::endl;
					movemap.set( TorsionID( parent_res, BRANCH, torsion_id ), action );
					branched_connection = true;
					break; //Only a single 'torsion' will be set here.
				}
			}

		} else if ( parent_res != 0 ) {
			//MainChain Linkage
			core::Size n_torsions_i = get_n_glycosidic_torsions_in_res(pose, i);
			core::Size phi_torsion_parent = pose.residue_type(parent_res).mainchain_atoms().size();
			TR.Debug << "Setting mainchain connection, Parent Mainchain torsions:  "<< phi_torsion_parent <<" Real Res Linkage Torsions " << n_torsions_i << std::endl;
			//n_torsions_i = 3;
			//Phi = phi_torsion_parent
			//psi = phi_torsion_parent -1
			//omega = phi_torsion_parent - 2
			//omega = phi_torsion_parent - n_torsions_i + 1
			for ( core::Size torsion = phi_torsion_parent; torsion >= phi_torsion_parent - n_torsions_i +1; --torsion ) {
				TR.Debug << "Setting BB Torsion: " << parent_res << " " << torsion << std::endl;
				movemap.set( TorsionID( parent_res, BB, torsion), action);
			}
		}
	}

	//  It either has no parent, or the parent is not carboydrate and we already set it up
	if ( parent_res == 0 || !pose.residue_type(parent_res).is_carbohydrate() ) {
		return;
	}

	//Real BB of a carbohydrate-carbohydrate linkage
	if ( branched_connection ) {
		for ( core::Size chi_torsion_num = 1; chi_torsion_num <= pose.residue_type( parent_res ).nchi(); ++ chi_torsion_num ) {
			TorsionID torsion_id = TorsionID(parent_res, CHI, chi_torsion_num);
			core::Size child_resnum = get_downstream_residue_that_this_torsion_moves( pose, torsion_id);
			//This is the BB Torsion we are after.  It is the same one as IUPAC.
			if ( child_resnum == i ) {
				TR.Debug << "Setting BB CHI Torsion " << chi_torsion_num << std::endl;
				movemap.set( torsion_id, action);
			}
		}
	}
}


///@brief
///
///  Get a pre-configured movemap from a residue selector.
///  Use the ReturnResidueSubsetSelector to obtain from a subset.
///
///  The Rosetta Movemap is VERY different from IUPAC-designated torsions for carbohydrates.
///  NEVER attempt to create a MoveMap for carbohydrates unless you know what you are doing.
///
///@details
///
/// This will create a Movemap from the residue selector for ALL residues within it.
/// including non-carbohydrates. This includes Chemical Edge Branch points, Mainchains, ASN->glycan linakge, etc.
///
///@params
///
/// include_iupac_chi:
///    Include the 'carbohydrate 'side-chains' (rotatable OH groups) and any selected non-carbohydrate side-chain
///
/// include_ring_torsions:
///    Include moveable ring torsions
///
/// include_bb_torsions:
///    Include torsions between both residues as defined by IUPAC. IE for linkage 1->5, the torsions of residue 5
///
kinematics::MoveMapOP
create_glycan_movemap_from_residue_selector(
	core::pose::Pose const & pose,
	select::residue_selector::ResidueSelectorCOP selector,
	bool include_iupac_chi /* true */,
	bool include_glycan_ring_torsions /* false */,
	bool include_bb_torsions /* true */)
{
	using namespace core::id;

	core::kinematics::MoveMapOP movemap = core::kinematics::MoveMapOP( new core::kinematics::MoveMap());

	// Make sure everything is OFF first.
	movemap->set_bb(false);
	movemap->set_chi(false);
	movemap->set_branches(false);
	movemap->set_nu(false);

	utility::vector1< bool > subset = selector->apply( pose );
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( ! subset[i] ) continue;

		TR.Debug << "Residue " << i << std::endl;
		if ( ! pose.residue( i ).is_carbohydrate() ) {
			if ( include_bb_torsions ) {
				movemap->set_bb(i, true);
			}
			if ( include_iupac_chi ) {
				movemap->set_chi(i, true);
			}

		} else {

			if ( include_glycan_ring_torsions ) {
				movemap->set_nu(i, true);
			}

			if ( include_iupac_chi ) {
				set_glycan_iupac_chi_torsions(pose, *movemap, i);
			}

			if ( include_bb_torsions ) {
				set_glycan_iupac_bb_torsions(pose, *movemap, i);
			}
		}
	}

	return movemap;
}


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
/// size, where "f" is furanose, "p" is pyranose, and "s" is septanose.\n
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
	std::string const & iupac_sequence,
	bool const idealize_linkages /*true*/,
	bool keep_pdbinfo /* true */ )
{
	using namespace std;
	using namespace utility;
	using namespace chemical;
	using namespace conformation;

	Residue const & residue( pose.residue( sequence_position ) );
	ResidueTypeSetCOP residue_set( pose.residue_type_set_for_pose( residue.type().mode() ) );

	// Get list of carbohydrate ResidueTypes from which to construct the Pose.
	ResidueTypeCOPs residue_types( residue_types_from_saccharide_sequence( iupac_sequence, *residue_set ) );

	Size const n_types( residue_types.size() );
	if ( ! n_types ) {
		TR.Warning << "No saccharide residues in sequence to append to pose!" << endl;
		return;
	}

	if ( TR.Trace.visible() ) {
		TR.Trace << "Ready to append the following residue types:  " << endl;
		for ( auto residue_type : residue_types ) {
			TR.Trace << residue_type->name() << ' ';
		}
		TR.Trace << endl;
	}

	Size const initial_nres( pose.size() );  // used below for setting up PDBInfo

	// Prepare the glycosylation site for branching.
	// TODO: Add code to attach to lipids.
	if ( residue.is_carbohydrate() ) {
		chemical::VariantType const variant(
			chemical::carbohydrates::CarbohydrateInfoManager::branch_variant_type_from_atom_name( atom_name ) );
		add_variant_type_to_pose_residue( pose, variant, sequence_position );
	} else {  // It's a typical branch point, in that it has a single type of branch point variant.
		add_variant_type_to_pose_residue( pose, SC_BRANCH_POINT, sequence_position );
	}


	// Now we can extend the Pose.
	append_pose_with_glycan_residues( pose, residue_types, sequence_position );

	// Let the Conformation know that it (now) contains sugars.
	pose.conformation().contains_carbohydrate_residues( true );

	//Size const final_nres = pose.size();
	//Size const total_glycan_residues = final_nres - initial_nres;
	Size const glycan_start_rosetta_num = initial_nres + 1;


	// Keep in an-tact PDBInfo and organize it nicely.
	if ( keep_pdbinfo && pose.pdb_info() != nullptr ) {

		char protein_chain = pose.pdb_info()->chain( sequence_position );
		int max_protein_chain_pdbnum = 0;

		//If you have a residue number that is negative - what the hell is that?

		//Find the maximum number in the chain of that the glycan is attached to.
		for ( core::Size i = 1; i <= initial_nres; ++i ) {
			if ( (pose.pdb_info()->chain(i) == protein_chain) && (pose.pdb_info()->number(i) > max_protein_chain_pdbnum) ) {
				max_protein_chain_pdbnum = pose.pdb_info()->number(i);
			}
		}

		//Set the chains
		int current_pdbnum = max_protein_chain_pdbnum + 1;
		for ( core::Size  i = glycan_start_rosetta_num; i <= pose.size(); ++i ) {

			//Actually set the PDBInfo to be correct.
			pose.pdb_info()->set_resinfo(i, protein_chain, current_pdbnum);
			current_pdbnum+=1;
		}

		//Give the Mainchain residue number
		//JAB- This is the more complicated way where the mainchain first gets numbers, and then each branch of children do.
		// This is not completely necessary (AFAIK), and this code is not yet finished.
		//utility::vector1< bool > mainchain_children = get_mainchain_children(pose, initial_nres+1, true /* include_initial */);
		//for (core::Size i = 1; i <= pose.size(); ++i){
		// if ( ! mainchain_children[ i ]) continue;
		// pose.pdb_info()->number( i, current_res );
		// current_res +=1;
		//
		//}

		pose.pdb_info()->resize_atom_records( pose );
		pose.pdb_info()->rebuild_pdb2pose(); //JAB - this may not be necessary - not sure.

		pose.pdb_info()->obsolete(false);
	} else {
		PDBInfoOP info( new PDBInfo( pose ) );
		info->name( pose.sequence() );  // Use the sequence as the default name.
		pose.pdb_info( info );
	}

	//TR << "InitialNRES: " << initial_sizes << " NRES: "<< pose.size() << " PDBINFO: "<< pose.pdb_info()->nres() << std::endl;
	TR << "Glycosylated pose with " << iupac_sequence << '-' << atom_name <<
		pose.residue( sequence_position ).name3() << sequence_position << endl;

	//Make sure that the new pose has a GlycanTreeSetObserver.
	if ( ! pose.glycan_tree_set() ) {
		pose.conformation().setup_glycan_trees();
	}

	if ( idealize_linkages ) {
		TR << "Idealizing glycosidic torsions." << endl;
		idealize_last_n_glycans_in_pose( pose, n_types );
	}
}

// Glycosylate the Pose at the given sequence position using an IUPAC sequence.
/// @details  This is a wrapper function for standard AA cases, i.e., glycosylation at Asn, Thr, Ser, or Trp.
void
glycosylate_pose(
	Pose & pose,
	uint const sequence_position,
	std::string const & iupac_sequence,
	bool const idealize_linkages /*true*/,
	bool keep_pdb_info /*true*/ )
{
	std::string const & glycosylation_site( pose.residue( sequence_position ).name3() );
	if ( glycosylation_site == "ASN" ) {
		glycosylate_pose( pose, sequence_position, "ND2", iupac_sequence, idealize_linkages, keep_pdb_info );
	} else if ( glycosylation_site == "SER" ) {
		glycosylate_pose( pose, sequence_position, "OG", iupac_sequence, idealize_linkages, keep_pdb_info );
	} else if ( glycosylation_site == "THR" ) {
		glycosylate_pose( pose, sequence_position, "OG1", iupac_sequence, idealize_linkages, keep_pdb_info );
	} else if ( glycosylation_site == "TRP" ) {
		glycosylate_pose( pose, sequence_position, "CD1", iupac_sequence, idealize_linkages, keep_pdb_info );
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
	std::string const & filename,
	bool const idealize_linkages /*true*/,
	bool keep_pdb_info /*true*/)
{
	using namespace std;
	using namespace io::carbohydrates;

	string const & sequence_file( core::chemical::carbohydrates::find_glycan_sequence_file( filename ) );
	if ( sequence_file.find( ".gws" ) != string::npos ) {
		utility_exit_with_message( "Rosetta cannot yet read .gws format sequence files." );  // TODO: Change this.
		//string const & gws_sequence( core::chemical::carbohydrates::read_glycan_sequence_file( sequence_file ) );
	} else {  // Assume any other file type contains an IUPAC sequence.
		string const & iupac_sequence( core::chemical::carbohydrates::read_glycan_sequence_file( sequence_file ) );
		glycosylate_pose( pose, sequence_position, atom_name, iupac_sequence, idealize_linkages, keep_pdb_info );
	}
}

// Glycosylate the Pose at the given sequence position using a .GWS or IUPAC sequence file.
/// @details  This is a wrapper function for standard AA cases, i.e., glycosylation at Asn, Thr, Ser, or Trp.
void
glycosylate_pose_by_file(
	Pose & pose,
	uint const sequence_position,
	std::string const & filename,
	bool const idealize_linkages /*true*/,
	bool keep_pdb_info /*true*/ )
{
	std::string const & glycosylation_site( pose.residue( sequence_position ).name3() );
	if ( glycosylation_site == "ASN" ) {
		glycosylate_pose_by_file( pose, sequence_position, "ND2", filename, idealize_linkages, keep_pdb_info );
	} else if ( glycosylation_site == "SER" ) {
		glycosylate_pose_by_file( pose, sequence_position, "OG", filename, idealize_linkages, keep_pdb_info );
	} else if ( glycosylation_site == "THR" ) {
		glycosylate_pose_by_file( pose, sequence_position, "OG1", filename, idealize_linkages, keep_pdb_info );
	} else if ( glycosylation_site == "TRP" ) {
		glycosylate_pose_by_file( pose, sequence_position, "CD1", filename, idealize_linkages, keep_pdb_info );
	} else {
		utility_exit_with_message( glycosylation_site + " is not a common site of glycosylation or else it is "
			"ambiguous; Rosetta cannot determine attachment atom.  Use glycosylate_pose_by_file( Pose & pose, uint "
			"const sequence_position, std::string const & atom_name, std::string const & filename ) instead." );
	}
}


// Return pointers to the two residues of the glycosidic bond.
/// @return  Pointers to the residue at <sequence_position> and its parent or else the same pointer twice if undefined.
std::pair< conformation::ResidueCOP, conformation::ResidueCOP >
get_glycosidic_bond_residues( Pose const & pose, uint const sequence_position )
{
	return core::conformation::carbohydrates::get_glycosidic_bond_residues( pose.conformation(), sequence_position);
}


/// @brief  Return the AtomIDs of the four phi torsion reference atoms.
utility::vector1< id::AtomID >
get_reference_atoms_for_phi( Pose const & pose, uint const sequence_position ){
	return core::conformation::carbohydrates::get_reference_atoms_for_phi( pose.conformation(), sequence_position );
}

/// @brief  Return the AtomIDs of the four psi torsion reference atoms.
utility::vector1< id::AtomID >
get_reference_atoms_for_psi( Pose const & pose, uint const sequence_position ){
	return core::conformation::carbohydrates::get_reference_atoms_for_psi( pose.conformation(), sequence_position );
}

/// @brief  Return the AtomIDs of the four omega torsion reference atoms.
utility::vector1< id::AtomID >
get_reference_atoms_for_1st_omega( Pose const & pose, uint const sequence_position ){
	return core::conformation::carbohydrates::get_reference_atoms_for_1st_omega( pose.conformation(), sequence_position );
}

/// @brief  Return the AtomIDs of the four omega2 torsion reference atoms.
utility::vector1< id::AtomID >
get_reference_atoms_for_2nd_omega( Pose const & pose, uint const sequence_position ){
	return core::conformation::carbohydrates::get_reference_atoms_for_2nd_omega( pose.conformation(), sequence_position );
}

/// @brief  Return the AtomIDs of the four reference atoms for the requested torsion.
utility::vector1< id::AtomID >
get_reference_atoms( uint const named_torsion, Pose const & pose, uint const sequence_position )
{
	return conformation::carbohydrates::get_reference_atoms( named_torsion, pose.conformation(), sequence_position );
}


// Set coordinates of virtual atoms (used as angle reference points) within a saccharide residue of the given pose.
void
align_virtual_atoms_in_carbohydrate_residue( Pose & pose, uint const sequence_position ) {
	conformation::carbohydrates::align_virtual_atoms_in_carbohydrate_residue( pose.conformation(), sequence_position );
}


// TorsionID Queries //////////////////////////////////////////////////////////
// Is this is the phi torsion angle of a glycosidic linkage?
/// @details  Carbohydrate linkages are defined as the torsion angles leading back to the previous residue.  Much
/// of Rosetta code relies on TorsionIDs and assumes TorsionID( n, BB, 1 ) is phi.  For a sugar, phi (of the next
/// residue) is the last torsion, and the number of main-chain torsions varies per saccharide residue.
bool
is_glycosidic_phi_torsion( Pose const & pose, id::TorsionID const & torsion_id )
{
	return core::conformation::carbohydrates::is_glycosidic_phi_torsion( pose.conformation(), torsion_id);
}

///@brief Get the torsion num that this torsionID moves.
core::Size
which_glycosidic_torsion(Pose const & pose, id::TorsionID const & torsion_id ){
	if      ( is_glycosidic_phi_torsion( pose, torsion_id) ) return 1;
	else if ( is_glycosidic_torsion( pose, torsion_id, core::id::psi_dihedral) ) return 2;
	else if ( is_glycosidic_torsion( pose, torsion_id, core::id::omega_dihedral) ) return 3;
	else if ( is_glycosidic_torsion( pose, torsion_id, core::id::omega2_dihedral) ) return 4;
	else return 0;
}


// Is this is the psi torsion angle of a glycosidic linkage?
/// @details  Carbohydrates linkages are defined as the torsion angles leading back to the previous residue.  Much
/// of Rosetta code relies on TorsionIDs and assumes TorsionID( n, BB, 2 ) is psi.  For a sugar, psi (of the next
/// residue) is the penultimate torsion, and the number of main-chain torsions varies per saccharide residue.
/// @note  This function will return false if the TorsionID is a CHI and the torsion is already covered by an equivalent
/// BB.  In other words, TorsionIDs for CHIs only return true for branch connections, where the ONLY way to access that
/// torsion is through CHI.
bool
is_glycosidic_psi_torsion( Pose const & pose, id::TorsionID const & torsion_id )
{
	return is_glycosidic_torsion(pose, torsion_id, core::id::psi_dihedral);
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
is_glycosidic_omega_torsion( Pose const & pose, id::TorsionID const & torsion_id )
{
	return is_glycosidic_torsion(pose, torsion_id, core::id::omega_dihedral);
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
is_glycosidic_omega2_torsion( Pose const & pose, id::TorsionID const & torsion_id )
{
	return is_glycosidic_torsion(pose, torsion_id, core::id::omega2_dihedral);
}

///@brief Base function to reduce code-duplication in torsion queries.
bool
is_glycosidic_torsion(Pose const & pose, id::TorsionID const & torsion_id, core::id::MainchainTorsionType const & torsion_type){

	using namespace id;
	core::Size mainchain_offset = 0;
	core::Size connect_atom_num = 0;

	switch( torsion_type){
	case phi_dihedral :
		return is_glycosidic_phi_torsion( pose, torsion_id);

	case psi_dihedral :
		mainchain_offset = 1;
		connect_atom_num = 3;
		break;

	case omega_dihedral :
		mainchain_offset = 2;
		connect_atom_num = 4;
		break;

	case omega2_dihedral :
		utility_exit_with_message("Omega2 currently unsupported!");

	case omega3_dihedral :
		utility_exit_with_message("Omega3 currently unsupported");

	default :
		utility_exit_with_message("Unknown behavior");
	}

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
					if ( core::Size(torsion_type) > 2 && ! info->has_exocyclic_linkage_to_child_mainchain() ) {
						return false;
					}
					return ( torsion_id.torsion() == residue.n_mainchain_atoms() - mainchain_offset );

				}
			}
			break;
		case CHI :
			{

			if ( ! residue.is_upper_terminus() ) {
				Size const n_mainchain_connections( residue.n_polymeric_residue_connections() );
				// An omega angle will always have the fourth atom of its definition be a connect atom.
				uint const N_atom( residue.chi_atoms( torsion_id.torsion() )[ connect_atom_num ] );
				Size const n_branches( residue.n_non_polymeric_residue_connections() );
				for ( uint branch_num( 1 ); branch_num <= n_branches; ++branch_num ) {
					next_rsd_num = residue.residue_connection_partner( n_mainchain_connections + branch_num );
					conformation::Residue const & next_rsd( pose.residue( next_rsd_num ) );
					if ( ( next_rsd.is_carbohydrate() ) && ( residue.connect_atom( next_rsd ) == N_atom ) ) {
						return true;
					}
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
get_downstream_residue_that_this_torsion_moves( Pose const & pose, id::TorsionID const & torsion_id )
{
	return core::conformation::carbohydrates::get_downstream_residue_that_this_torsion_moves( pose.conformation(), torsion_id );

}


///@brief Get the number of glycosidic torsions for this residue.  Up to 4 (omega2).
Size
get_n_glycosidic_torsions_in_res( Pose const & pose, uint const sequence_position )
{
	return core::conformation::carbohydrates::get_n_glycosidic_torsions_in_res( pose.conformation(), sequence_position );
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
get_glycosidic_torsion( uint const named_torsion, Pose const & pose, uint const sequence_position )
{
	return conformation::carbohydrates::get_glycosidic_torsion( named_torsion, pose.conformation(), sequence_position );
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
set_glycosidic_torsion( uint const named_torsion, Pose & pose, uint const sequence_position, core::Angle const setting )
{
	conformation::carbohydrates::set_glycosidic_torsion( named_torsion, pose.conformation(), sequence_position, setting );
}


// Glycosylation //////////////////////////////////////////////////////////////
// Idealize the glycosidic torsion angles for the last n glycan residues added or built.
void
idealize_last_n_glycans_in_pose( Pose & pose, Size const n_glycans_added )
{
	using namespace std;
	using namespace conformation;

	// Work backward through the list of glycans and fix their glycosidic bonds to be ideal, if such data is known.
	Size const sizes( pose.size() );
	uint const glycosylation_site( sizes - n_glycans_added );  // 0, if a free oligo-/polysaccharide was built.
	for ( uint i( sizes ); i > glycosylation_site; --i ) {
		Residue const & non_red_end_res( pose.residue( i ) );
		if ( ! non_red_end_res.is_carbohydrate() ) {
			TR.Warning << "Cannot idealize a non-saccharide residue!" << endl;
			continue;
		}
		uint parent_seqpos( pose.glycan_tree_set()->get_parent( i ) );
		if ( parent_seqpos == 0 ) { continue; }
		Residue const & red_end_res( pose.residue( parent_seqpos ) );
		string const & non_red_end_short_name( non_red_end_res.carbohydrate_info()->short_name() );
		string red_end_short_name;
		if ( red_end_res.is_carbohydrate() ) {
			red_end_short_name = red_end_res.carbohydrate_info()->short_name();  // 3-letter code not enough
			uint const link_pos( pose.glycan_tree_set()->get_linkage_position( i ) );
			red_end_short_name[ 2 ] = '0' + link_pos;  // Set the correct connectivity.
		} else {
			red_end_short_name = red_end_res.name3();
		}
		TR.Trace << non_red_end_short_name << "(?-" << red_end_short_name << " linkage: ";
		if ( chemical::carbohydrates::CarbohydrateInfoManager::pair_has_linkage_statistics(
				red_end_short_name, non_red_end_short_name ) ) {
			TR.Trace << "found" << endl;
			// We will assume that the 1st conformer has the highest population and use that.
			chemical::carbohydrates::LinkageConformerData const conformer(
				chemical::carbohydrates::CarbohydrateInfoManager::linkages_from_pair(
				red_end_short_name, non_red_end_short_name )[ 1 ] );
			carbohydrates::set_dihedrals_from_linkage_conformer_data(
				pose, i, conformer );
		} else {
			TR.Trace << "not found" << endl;
			// TODO: Idealize from the sugar_bb scoring term.
		}
	}
}


// Set the dihedral angles involved in a glycosidic linkage based on statistical data.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
void set_dihedrals_from_linkage_conformer_data( Pose & pose,
	uint const upper_residue,
	core::chemical::carbohydrates::LinkageConformerData const & conformer,
	bool idealize /* true */,
	bool use_prob_for_sd /* false */ )
{
	using namespace core::id;

	core::Size const n_torsions = conformer.n_torsions();

	//TR << "Sampling Bin " << conformer.conformer_bins << " which has population of " << conformer.population << std::endl;

	if ( idealize ) {
		for ( core::Size i = 1; i <= n_torsions; ++i ) {
			TR.Debug << "torsion "<< i << std::endl;
			Angle torsion_mean;
			if ( i <= 3 ) {  // phi, psi, or omega1
				auto torsion_type = static_cast< MainchainTorsionType >( i );
				torsion_mean = basic::periodic_range( conformer.get_torsion_mean( torsion_type ), 360);
			} else {
				torsion_mean = basic::periodic_range( conformer.get_torsion_mean( omega_dihedral, i - 2 ), 360);
			}
			TR.Debug << "Setting from " << get_glycosidic_torsion( i, pose, upper_residue ) <<
				" to " << torsion_mean << std::endl;  // TEMP
			set_glycosidic_torsion( i, pose, upper_residue, torsion_mean );
		}
	} else if ( use_prob_for_sd ) {
		utility_exit_with_message( " UseProbForSD in LinkageConformer sampling is not yet implemented." );

	} else {
		Real mean;
		Real sd;
		Real conformer_angle;
		for ( core::Size i = 1; i <= n_torsions; ++i ) {
			TR.Debug << "torsion: " << i << std::endl;


			if ( i <= 3 ) { // phi, psi, or omega1
				auto torsion_type = static_cast< MainchainTorsionType >( i );
				mean = conformer.get_torsion_mean( torsion_type );
				sd = conformer.get_torsion_sd( torsion_type );
				conformer_angle = basic::periodic_range( mean - sd + numeric::random::rg().uniform() * sd * 2, 360.0 );
			} else {
				mean = conformer.get_torsion_mean( omega_dihedral, i - 2 );
				sd = conformer.get_torsion_sd( omega_dihedral, i - 2);
				conformer_angle = basic::periodic_range( mean - sd + numeric::random::rg().uniform() * sd * 2, 360.0 );
			}

			//Debugging for now.
			TR.Debug << "Current: " << get_glycosidic_torsion(i, pose, upper_residue) << std::endl;
			TR.Debug << "Mean: " << mean << std::endl;
			TR.Debug << "New conformer: " << conformer_angle << std::endl;

			set_glycosidic_torsion( i, pose, upper_residue, conformer_angle );
		}
	}
}


/// @brief Remove branch points from a carbohydrate or aa residue.
void
remove_carbohydrate_branch_point_variants( Pose & pose, core::Size const seqpos)
{
	using namespace core::chemical;

	if ( pose.residue(seqpos).is_branch_point() ) {
		if ( pose.residue(seqpos).is_carbohydrate() ) {
			utility::vector1< VariantType > branch_types;

			branch_types.push_back( C1_BRANCH_POINT );
			branch_types.push_back( C2_BRANCH_POINT );
			branch_types.push_back( C3_BRANCH_POINT );
			branch_types.push_back( C4_BRANCH_POINT );
			branch_types.push_back( C5_BRANCH_POINT );
			branch_types.push_back( C6_BRANCH_POINT );
			branch_types.push_back( C7_BRANCH_POINT );
			branch_types.push_back( C8_BRANCH_POINT );
			branch_types.push_back( C9_BRANCH_POINT );


			for ( core::Size i = 1; i <= branch_types.size(); ++i ) {

				//Must check current residue since we are editing - something gets strange with valgrind
				//   since we replace the residue on variant type removal
				if ( pose.residue(seqpos).has_variant_type( branch_types[ i ] ) ) {
					remove_variant_type_from_pose_residue( pose, branch_types[ i ], seqpos);
				}
			}

		} else {
			remove_variant_type_from_pose_residue( pose, SC_BRANCH_POINT, seqpos);
		}
	}
}

/// @brief find connected atom in neighboring residue
// pair: <residue number,main chain atom number
void
find_neighbor( Pose &pose, core::Size res_pos, core::Size atom_nr, std::pair<core::Size, core::Size> &neighbor )
{
	// TODO: set proper bond length cut offs
	core::Real max_cutoff = 3.0;
	core::Real min_cutoff = 1.0;
	for ( core::Size i = res_pos + 1; i <= res_pos + 20; ++res_pos ) {
		conformation::ResidueOP res( new conformation::Residue( pose.residue( i ) ) );
		chemical::ResidueTypeOP RT( new chemical::ResidueType( pose.residue_type( i ) ) );
		// Loop through main chain atoms of potential neighbor residue
		for ( unsigned long mc_atom : RT->mainchain_atoms() ) {
			core::Real distance = pose.residue(res_pos).xyz( atom_nr ).distance(pose.residue( i ).xyz( mc_atom ));
			if ( (distance < max_cutoff) | (distance > min_cutoff) ) {
				neighbor = std::make_pair(i,mc_atom);
			}
		}
	}
	// Apparently no neighbor was found

}

/// @brief Go through pose from non-Rosetta PDB file and infer and properly set up branching
void
setup_existing_glycans(Pose &pose)
{

	//chain_map = chainID, residue number //
	map< core::Size/*chainID*/,utility::vector1<core::Size>  > chain_map;
	for ( core::Size res_pos = 1; res_pos <= pose.total_residue(); ++res_pos ) {
		conformation::ResidueOP res( new conformation::Residue( pose.residue( res_pos ) ) );
		if ( res->is_carbohydrate() ) {
			//core::Size chain = res->chain();
			chemical::ResidueTypeOP RT( new chemical::ResidueType( pose.residue_type(res_pos) ) );
			for ( auto atom = RT->mainchain_atoms().begin(); atom != RT->mainchain_atoms().end(); ++atom ) {
				// Find the fisrt exocyclic oxigen
				std::string atom_name = RT->atom_name(*atom);
				if ( atom_name.find('O') != std::string::npos ) {  // is oxigen
					std::pair< core::Size,core::Size  > neighbor = std::make_pair(0,0);
					find_neighbor(pose,res_pos,*atom,neighbor);
					if ( neighbor.first == 0 || neighbor.second == 0 ) {
						// No neighbor neighbor found
						// Must be a lower terminus
						// TODO: add terminus variant type here
					} else {
						// Main chain continues at residue neighbor.first
						// Add branch variant type
						chemical::VariantType const variant(
							chemical::carbohydrates::CarbohydrateInfoManager::branch_variant_type_from_atom_name( atom_name ) );
						add_variant_type_to_pose_residue( pose, variant, res_pos );
					}

				}
			}
		}
	}
	PDBInfoOP info( new PDBInfo( pose ) );
	info->name( pose.sequence() );  // Use the sequence as the default name.
	pose.pdb_info( info );
}





////////////////////////////////////////    Branch Information    /////////////////////////////////////////
///
///
///








/// @brief Get residues further down the branch from this residue.  starting_position ->
/// @details Calls get_carbohydrate_residues_and_tips_of_branch
utility::vector1< core::Size >
get_carbohydrate_residues_of_branch(
	Pose const & pose,
	uint const starting_position)
{
	return get_carbohydrate_residues_and_tips_of_branch(pose, starting_position).first;
}



/// @brief Get tips (end residue of linear components of branches) further down the branch from this residue.  starting_position ->
/// @details Calls get_carbohydrate_residues_and_tips_of_branch
utility::vector1< core::Size >
get_carbohydrate_tips_of_branch(
	Pose const & pose,
	uint const starting_position)
{
	return get_carbohydrate_residues_and_tips_of_branch(pose, starting_position).second;
}



/// @brief Get residues further down the branch from this residue.  starting_position ->
///  May require a better name.
///  Returns pair of all_upstream_residues, tips.
///  Tips are the ends of linear glycan branches.

std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > >
get_carbohydrate_residues_and_tips_of_branch(
	Pose const & pose,
	uint const starting_position)

{
	return conformation::carbohydrates::get_carbohydrate_residues_and_tips_of_branch(pose.conformation(), starting_position);
}





///@brief Get the carbohydrate residue connecting the protein branch point.
core::Size
get_glycan_connecting_protein_branch_point(Pose const & pose, core::Size const protein_branch_point_resnum){

	return conformation::carbohydrates::get_glycan_connecting_protein_branch_point( pose.conformation(), protein_branch_point_resnum );

}


///@brief Get the particular resnum from a glycan position, givin the protein branch point.
/// The glycan_position is numbered 1 -> length of glycan. This is useful for easily identifying a particular glycan position.
///
core::Size
get_resnum_from_glycan_position(Pose const & pose, core::Size const glycan_one, core::Size const glycan_position){
	return conformation::carbohydrates::get_resnum_from_glycan_position( pose.conformation(), glycan_one, glycan_position);

}

core::Size
get_glycan_position_from_resnum(Pose const & pose, core::Size const glycan_one, core::Size const glycan_residue ){

	return conformation::carbohydrates::get_glycan_position_from_resnum( pose.conformation(), glycan_one, glycan_residue );

}



select::residue_selector::ResidueSubset
get_resnums_from_glycan_positions(Pose const & pose, core::Size const glycan_one, utility::vector1< core::Size > const & glycan_positions){

	utility::vector1< bool > subset(pose.total_residue(), false);

	core::Size glycan_length = pose.glycan_tree_set()->get_tree(glycan_one)->get_size();
	for ( core::Size index = 1; index <= glycan_positions.size(); ++index ) {

		Size glycan_position = glycan_positions[ index ];
		if ( glycan_position <= glycan_length ) {

			Size glycan_residue = get_resnum_from_glycan_position(pose, glycan_one, glycan_position);

			if ( glycan_residue <= pose.total_residue() ) {
				subset[ glycan_residue ] = true;
			} else {
				utility_exit_with_message("Glycan residue not in pose: "+utility::to_string( glycan_residue) + " corresponding to glycan position "+utility::to_string( glycan_position));
			}

		}
	}
	return subset;
}

select::residue_selector::ResidueSubset
get_resnums_from_glycan_positions(Pose const & pose, utility::vector1< core::Size > const & glycan_positions){
	select::residue_selector::OrResidueSelector combine_subsets = select::residue_selector::OrResidueSelector();

	utility::vector1< bool > subset(pose.total_residue(), false);
	utility::vector1< bool > glycan_start_points =  pose.glycan_tree_set()->get_start_points();

	for ( Size i = 1; i <= glycan_start_points.size(); ++i ) {
		if ( ! glycan_start_points[ i ] ) continue;

		select::residue_selector::ResidueSubset single_glycan = get_resnums_from_glycan_positions(pose, i, glycan_positions);
		combine_subsets.apply_or_to_subset(single_glycan, subset);
	}
	return subset;
}


////////////////////////////////////////    Branch Deletion    /////////////////////////////////////////
///
///
///

utility::vector1< bool >
get_mainchain_children( Pose const & pose, core::Size starting_resnum, bool include_starting_resnum){

	utility::vector1< bool > mc_children( pose.size(), false);

	if ( include_starting_resnum ) mc_children[ starting_resnum ] = true;
	core::Size child = pose.glycan_tree_set()->get_node( starting_resnum )->get_mainchain_child();
	while ( child != 0 ) {
		mc_children[ child ] = true;
		child = pose.glycan_tree_set()->get_node( child )->get_mainchain_child();
	}

	return mc_children;

}

/// @brief Delete the glycan from this residue onward toward the end of the branch.  Like chopping off a tree trunk at position resnum (not including resnum). Also known as defoliating.
///  If resnum is the protein branch point, will change variant.
//   If no more carbohydrates exist in the pose, will change the pose status.
void
delete_carbohydrate_branch( Pose & pose, uint const delete_to )
{
	using namespace utility;
	using namespace chemical;
	using namespace conformation;

	if ( ! pose.residue( delete_to ).is_carbohydrate() && ! pose.residue( delete_to ).is_branch_point() ) {
		TR << "Delete to residue is not carbohydrate and not a branch point.  Nothing to be done." << std::endl;
		return;
	}

	TR << "Delete Carb Branch Starting Points " << pose.glycan_tree_set()->get_start_points() << std::endl;

	std::pair< vector1< Size >, vector1< Size > > resnums_and_tips;
	resnums_and_tips = get_carbohydrate_residues_and_tips_of_branch( pose, delete_to );

	vector1< Size > resnums_to_delete = resnums_and_tips.first;
	vector1< Size > tips              = resnums_and_tips.second;

	Size i = 0;
	std::string ref_pose_name = "temp_ref_pose" ;

	Size current_delete_to = delete_to;

	while ( tips.size() > 0 ) {

		pose.reference_pose_from_current( ref_pose_name, true );

		vector1< vector1< Size > > leafs;

		//Get linear carbohydrate leafs
		for ( core::Size x = 1; x <= tips.size(); ++x ) {
			Size tip = tips[ x ];
			vector1< Size > leaf = get_resnums_in_leaf_on_the_fly(pose, tip, delete_to);
			leafs.push_back(leaf);
		}

		//Delete the leafs
		for ( core::Size x = 1; x <= leafs.size(); ++x ) {
			vector1< Size > leaf = leafs[ x ];
			delete_leaf( pose, leaf, ref_pose_name);
		}

		//Get new tips to delete.
		current_delete_to = pose.corresponding_residue_in_current( current_delete_to, ref_pose_name);
		resnums_and_tips = get_carbohydrate_residues_and_tips_of_branch( pose, current_delete_to );

		resnums_to_delete = resnums_and_tips.first;
		tips              = resnums_and_tips.second;


		//TR << pose << std::endl; //Debugging for now.

		i+=1;
	}

	// Update branching of starting_position residue.  After all this deletion, it will not have any branches and should be a tip.
	//   -> Deal with branch of aa side chain or carbohydrate branch properly!

	remove_carbohydrate_branch_point_variants( pose, current_delete_to );


	// Update contents of Pose - does it have any carbohydrates left?
	bool found_carbohydrate = false;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue( i ).is_carbohydrate() ) {
			found_carbohydrate = true;
			break;
		}
	}
	pose.conformation().contains_carbohydrate_residues( found_carbohydrate );
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
get_branching_residues( Pose const & pose,
	Size parent_residue,
	utility::vector1< Size > & children_residues,
	utility::vector1< Size > & list_of_residues,
	utility::vector1< Size > & tips )
{
	return conformation::carbohydrates::get_branching_residues( pose.conformation(), parent_residue, children_residues, list_of_residues, tips );
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
fill_upstream_children_res_and_tips( Pose const & pose,
	Size res,
	Size parent_residue,
	utility::vector1< Size > & children_residues,
	utility::vector1< Size > & list_of_residues,
	utility::vector1< Size > & tips )
{
	return conformation::carbohydrates::fill_upstream_children_res_and_tips(pose.conformation(), res, parent_residue, children_residues, list_of_residues, tips);
}


/// @brief Get all residue numbers in order from the tip to (and not including) stop_at_residue or a branch point.
///  All residue numbers are the tip or a linear polymer of glycans.
///  Useful for glycan stripping.
utility::vector1< core::Size >
get_resnums_in_leaf( Pose const & pose, Size tip_residue, Size stop_at_residue )
{
	// This logic can be simplified even further, but this works.

	utility::vector1< Size > resnums;

	//If residue is not carbohydrate, we can't find the parent residue.
	/// We really shouldn't be here in the first place!
	if ( ! pose.residue( tip_residue ).is_carbohydrate() ) {
		return resnums;
	}

	resnums.push_back( tip_residue );
	Size res = tip_residue;

	while ( pose.residue( res ).n_current_residue_connections() < 3 ) {

		Size parent = pose.glycan_tree_set()->get_parent( res );

		if ( res != tip_residue ) { resnums.push_back( res ); }
		if ( parent == 0 || parent == stop_at_residue ) { break; }

		res = parent;
	}
	return resnums;
}

/// @brief Get all residue numbers in order from the tip to (and not including) stop_at_residue or a branch point.
///  All residue numbers are the tip or a linear polymer of glycans.
///  Useful for glycan stripping.
utility::vector1< core::Size >
get_resnums_in_leaf_on_the_fly( Pose const & pose, Size tip_residue, Size stop_at_residue )
{
	// This logic can be simplified even further, but this works.

	utility::vector1< Size > resnums;

	//If residue is not carbohydrate, we can't find the parent residue.
	/// We really shouldn't be here in the first place!
	if ( ! pose.residue( tip_residue ).is_carbohydrate() ) {
		return resnums;
	}

	resnums.push_back( tip_residue );
	Size res = tip_residue;

	while ( pose.residue( res ).n_current_residue_connections() < 3 ) {

		Size parent = find_seqpos_of_saccharides_parent_residue( pose.residue( res ) );

		if ( res != tip_residue ) { resnums.push_back( res ); }
		if ( parent == 0 || parent == stop_at_residue ) { break; }

		res = parent;
	}
	return resnums;
}

/// @brief Delete a leaf of glycan residues.  Use the ReferencePose to associate residue numbers.
/// @details
///  Uses delete_residue_slow as the debug_assert in delete_polymer_residue causes a crash.
///  This is a bug in the FoldTree where chemical edges are being treated as Jumps.
///  This is being addressed in a branch and a FoldTree.
void
delete_leaf( Pose & pose, utility::vector1< Size > leaf, std::string ref_pose_name )
{
	for ( core::Size i = 1; i <= leaf.size(); ++i ) {
		Size resnum = leaf[ i ];
		Size new_resnum = pose.corresponding_residue_in_current( resnum, ref_pose_name);

		// Uses delete_residue_slow as the debug_assert in delete_polymer_residue causes a crash.
		// This exit checks for foldtree jump, and it looks like it counts chemical edges as a foldtree jump.
		// UNTIL WE CAN RESOLVE THIS, use delete_residue_slow here.

		//pose.delete_polymer_residue( new_resnum );
		pose.delete_residue_slow( new_resnum );
	}
}

}  // namespace carbohydrates
}  // namespace pose
}  // namespace core
