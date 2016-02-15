// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/pose_to_sfr/PoseToStructFileRepConverter.cc
/// @brief A class to convert a pose to a StructFileRep.
/// @details This conversion is a first step in PDB or mmCIF output.  It could be useful for other
/// input/output, too.
/// @author Vikram K. Mulligan (vmullig@uw.edu), XRW 2016 Team
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com), XRW 2016 Team

// Unit headers
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/Remarks.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/io/HeaderInformation.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/scoring/dssp/Dssp.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// External headers
#include <ObjexxFCL/format.hh>

// C++ Headers
#include <sstream>

namespace core {
namespace io {
namespace pose_to_sfr {



static THREAD_LOCAL basic::Tracer TR( "core.io.pose_to_sfr.PoseToStructFileRepConverter" );

// Pose to StructFileRep methods ///////////////////////////////////////////////////

/// @brief Constructor.
/// @details Creates the StructFileRep object, which is subsequently accessible by owning pointer
/// (using the sfr() object).  The options_ object is created and initialized from the options system.
PoseToStructFileRepConverter::PoseToStructFileRepConverter() :
	utility::pointer::ReferenceCount(),
	options_() //Initialize from options system.
{
	new_sfr();
}

/// @brief Constructor with options input.
/// @details Creates the StructFileRep object, which is subsequently accessible by owning pointer
/// (using the sfr() object).  The options_ object is duplicated from the input options object.
PoseToStructFileRepConverter::PoseToStructFileRepConverter( StructFileRepOptions const &options_in ) :
	utility::pointer::ReferenceCount(),
	options_( options_in )
{
	new_sfr();
}


/// @brief Resets the PoseToStructFileRepConverter object, and reinitializes
/// it with a fresh StruftFileRep object, returning an owning pointer to the
/// new object.
core::io::StructFileRepOP
PoseToStructFileRepConverter::new_sfr() {
	sfr_ = StructFileRepOP( new core::io::StructFileRep() );
	return sfr_;
}


/// @details Non-const access to the StructFileRep object.
///
StructFileRepOP
PoseToStructFileRepConverter::sfr() {
	return sfr_;
}


/// @details Read atoms/residue information from Pose object and put it in StructFileRep object.
void
PoseToStructFileRepConverter::init_from_pose(core::pose::Pose const & pose)
{
	init_from_pose( pose, options_ );
}

/// @details Read atoms/residue information from Pose object and put it in StructFileRep object using options defined in
/// StructFileRepOptions.
void
PoseToStructFileRepConverter::init_from_pose(core::pose::Pose const & pose, StructFileRepOptions const & options)
{
	id::AtomID_Mask mask = id::AtomID_Mask( pose.total_residue() );
	for ( Size resnum = 1; resnum <= pose.total_residue(); ++ resnum ) {
		for ( Size atom_index = 1; atom_index <= pose.residue( resnum ).natoms(); ++ atom_index ) {
			id::AtomID atm = id::AtomID( atom_index, resnum );
			mask.set( atm, true );
		}
	}

	init_from_pose( pose, mask, options );
}

/// @details A lightweight, direct way of limiting pose pdb output to a subset of residues.
/// @note The alternative of constructing new subposes for output only would be unnecessary/less efficient (?)
void
PoseToStructFileRepConverter::init_from_pose(
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & residue_indices
)
{
	using namespace core;

	id::AtomID_Mask mask = id::AtomID_Mask( pose.total_residue() );
	for ( core::Size resnum = 1; resnum <= pose.total_residue(); ++resnum ) {

		for ( Size atom_index = 1; atom_index <= pose.residue( resnum ).natoms(); ++ atom_index ) {
			id::AtomID atm = id::AtomID( atom_index, resnum );
			if ( residue_indices.has_value( resnum ) ) {
				mask.set( atm, true );
			} else {
				mask.set( atm, false );
			}
		}

	}

	init_from_pose( pose, mask );

}

void
PoseToStructFileRepConverter::init_from_pose(const core::pose::Pose &pose, const id::AtomID_Mask &mask)
{

	init_from_pose( pose, mask, options_);

}


////////////////////////////////// Main Init Function //////////////
void
PoseToStructFileRepConverter::init_from_pose(const core::pose::Pose &pose, const id::AtomID_Mask &mask, StructFileRepOptions const & options){


	using namespace core;
	using core::pose::PDBInfo;

	new_sfr();
	if( &options_ != &options ) { // Yes, address-of comparison
		// In an if block to avoid self (re)assignment, say from the init_from_pose(pose) call.
		options_ = options;
	}

	// Get Title Section information.
	if ( (options.preserve_header() == true || options.preserve_crystinfo() == true ) && pose.pdb_info() ) {
		*(sfr_->remarks()) = pose.pdb_info()->remarks();  // Get OP to PDBInfo object for remarks.
		if ( pose.pdb_info()->header_information() ) {
			sfr_->header() = core::io::HeaderInformationOP( new core::io::HeaderInformation(*(pose.pdb_info()->header_information())) );
		} else {
			sfr_->header() = core::io::HeaderInformationOP( new core::io::HeaderInformation() );
		}
	} else {
		sfr_->header() = core::io::HeaderInformationOP( new core::io::HeaderInformation() );
	}

	// Get parametric information
	if ( options.write_pdb_parametric_info() ) {
		get_parametric_info(sfr_->remarks(), pose);
	}

	// Get Connectivity Annotation Section information.
	if ( options.write_pdb_link_records() ) {
		get_connectivity_annotation_info( pose );
	}

	// Get Crystallographic and Coordinate Transformation Section information.
	if ( options.preserve_crystinfo() && pose.pdb_info() ) {
		sfr_->crystinfo() = pose.pdb_info()->crystinfo();
	}

	// Setup options.

	bool renumber_chains(false);
	if ( options.per_chain_renumbering() ) {
		renumber_chains = true;
	}

	utility::vector1< bool > res_info_added( pose.total_residue(), false );

	core::Size new_atom_num(1);
	core::Size new_tercount(0);
	for ( Size resnum=1, resnum_max=pose.n_residue(); resnum<=resnum_max; ++resnum ) {
		conformation::Residue const & rsd( pose.residue( resnum ) );
		bool use_pdb_info = use_pdb_info_for_num( pose, resnum );

		if ( resnum > 1 && pose.chain(resnum) != pose.chain(resnum-1) ) ++new_tercount;

		for ( Size atom_index=1; atom_index<= rsd.natoms(); ++atom_index ) {

			id::AtomID atm = id::AtomID( atom_index,resnum );
			if ( ! mask[ atm ]  || ! mask.has( atm ) ) continue;
			ResidueInformation res_info;
			get_residue_information(pose, resnum, use_pdb_info, renumber_chains, new_tercount, res_info);

			///Add res information only if we havn't done so already.
			if ( ! res_info_added[ resnum ] ) {
				append_residue_info_to_sfr(pose, res_info, rsd);
				res_info_added[ resnum ] = true;
			}

			bool success = append_atom_info_to_sfr(pose, res_info, rsd, atom_index, use_pdb_info, new_atom_num, new_tercount);
			if ( success ) new_atom_num += 1;
		}
	}

}


// Append pdb information to StructFileRep for a single residue.
void
PoseToStructFileRepConverter::append_residue_to_sfr(
	core::pose::Pose const & pose,
	core::Size const resnum
) {

	core::Size new_atom_num_start = get_new_atom_serial_num();
	core::Size const new_tercount( pose.chain(resnum)-1 );
	append_residue_to_sfr( pose, resnum, new_atom_num_start, new_tercount );

}

/// @brief Append pdb information to StructFileRep for a single residue.
///  @brief Start atom numbering from given n
void
PoseToStructFileRepConverter::append_residue_to_sfr(
	pose::Pose const & pose,
	Size const resnum,
	Size & new_atom_num,
	Size const new_tercount )
{

	using namespace core;
	using namespace utility;

	conformation::Residue const & rsd = pose.residue( resnum );

	bool use_pdb_info = use_pdb_info_for_num( pose, resnum );
	bool renumber_chains(false);
	if ( options_.per_chain_renumbering() ) {
		renumber_chains = true;
	}

	ResidueInformation res_info;
	get_residue_information(pose, resnum, use_pdb_info, renumber_chains, new_tercount, res_info );
	append_residue_info_to_sfr( pose, res_info, rsd );

	// Loop through each atom in the residue and generate ATOM or HETATM data.
	for ( Size j = 1; j <= rsd.natoms(); ++j ) {
		bool success = append_atom_info_to_sfr( pose, res_info, rsd, j, use_pdb_info, new_atom_num, new_tercount);
		if ( success ) new_atom_num += 1;
	}

}


/// @brief Append just residue-based info to StructFileRep
void PoseToStructFileRepConverter::append_residue_info_to_sfr(
	pose::Pose const & /* pose */,
	ResidueInformation const & res_info,
	conformation::Residue const & rsd )
{

	// Determine residue identifier information.

	// Generate HETNAM data, if applicable.
	// TODO: For now, only output HETNAM records for saccharide residues, but in the future, outputting HETNAM records
	// for any HETATM residues could be done.
	if ( rsd.is_carbohydrate() ) {
		String const & hetID = rsd.name3();

		String resSeq( utility::pad_left(res_info.resSeq(), 4) ); //("%4d", res_info.resSeq);
		String const resID = String(1, res_info.chainID()) + resSeq + String(1, res_info.iCode());

		String const text = resID + " " + rsd.type().base_name();

		sfr_->heterogen_names().push_back(make_pair(hetID, text));
	}
}


bool
PoseToStructFileRepConverter::append_atom_info_to_sfr(
	core::pose::Pose const & pose,
	ResidueInformation const & res_info,
	core::conformation::Residue const & rsd,
	core::Size const atom_index,
	bool const use_pdb_info)
{
	return append_atom_info_to_sfr( pose, res_info, rsd, atom_index, use_pdb_info, get_new_atom_serial_num() /*atom index*/, pose.chain( rsd.seqpos() ) - 1 /*Number of termini before this atom*/);
}


bool
PoseToStructFileRepConverter::append_atom_info_to_sfr(
	core::pose::Pose const & pose,
	ResidueInformation const & res_info,
	core::conformation::Residue const & rsd,
	core::Size const atom_index,
	bool const use_pdb_info,
	core::Size const new_atom_num,
	core::Size const new_tercount)

{

	pose::PDBInfoCOP pdb_info = pose.pdb_info();
	conformation::Atom const & atom = rsd.atom( atom_index ) ;

	if ( ! add_atom_to_sfr( pose, rsd, atom_index, use_pdb_info ) ) return false;

	AtomInformation ai;
	AtomInformation orb;  //have to initialize this out here.

	ai.isHet = (!rsd.is_polymer() || rsd.is_ligand());
	ai.chainID = res_info.chainID();
	ai.resSeq = res_info.resSeq();
	ai.iCode = res_info.iCode();
	ai.serial = new_atom_num;
	ai.name = rsd.atom_name( atom_index );
	ai.resName = rsd.name3();
	ai.x = atom.xyz()(1);
	ai.y = atom.xyz()(2);
	ai.z = atom.xyz()(3);
	ai.occupancy = 1.0; // dummy occupancy, can be overridden by PDBInfo
	ai.segmentID = res_info.segmentID();
	ai.terCount = new_tercount;

	// Output with pdb-specific info if possible.
	if ( use_pdb_info ) {
		if ( pdb_info->is_het( rsd.seqpos(), atom_index ) ) { // override standard het only if .is_het() is true
			ai.isHet = true;
		}
		ai.altLoc = pdb_info->alt_loc( rsd.seqpos(), atom_index );
		ai.occupancy = pdb_info->occupancy( rsd.seqpos(), atom_index );
		ai.temperature = pdb_info->temperature( rsd.seqpos(), atom_index );
		ai.segmentID = pdb_info->segmentID( rsd.seqpos() );
	}

	if ( options_.output_virtual_zero_occ() ) {
		if ( rsd.atom_type( atom_index ).is_virtual() ) {
			ai.occupancy = 0.0;
		}
	}

	grab_conect_records_for_atom( pose, rsd.seqpos(), atom_index, ai );

	// Element
	// (written by fpd; moved here by Labonte)
	core::chemical::AtomTypeSet const &ats = rsd.type().atom_type_set();
	ai.element = ats[atom.type()].element();
	if ( ai.element.length() == 1 ) ai.element = " "+ai.element;

	// 'chains' is member data
	if ( sfr_->chains().size() < static_cast <core::Size> (rsd.chain() + 1) ) sfr_->chains().resize( rsd.chain() + 1 );
	sfr_->chains()[rsd.chain()].push_back(ai);

	return true;
}

core::Size
PoseToStructFileRepConverter::get_new_atom_serial_num() const {
	core::Size new_atom_num;
	if ( total_sfr_atoms( *sfr_ ) == 0 ) {
		new_atom_num = 1;
	} else {
		core::Size total_chains = sfr_->chains().size();
		core::Size total_entries_last_chain = total_sfr_atoms( *sfr_, total_chains - 1);
		new_atom_num = sfr_->chains()[ total_chains - 1][ total_entries_last_chain - 1].serial;
	}

	return new_atom_num;
}

bool
PoseToStructFileRepConverter::add_atom_to_sfr(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	core::Size atom_index,
	bool use_pdb_info ) const

{

	//skip outputting virtual atom unless specified
	if ( !options_.output_virtual() &&
			rsd.atom_type( atom_index ).is_virtual() ) return false;

	//fpd optionally don't output centroids
	if ( options_.no_output_cen() &&
			rsd.atom_name( atom_index ) == " CEN" ) return false;

	// skip outputting zero occupancy atoms if specified
	if ( use_pdb_info && options_.suppress_zero_occ_pdb_output() &&
			( rsd.seqpos() <= pose.pdb_info()->nres() ) ) {
		if ( pose.pdb_info()->occupancy( rsd.seqpos(), atom_index ) < 0.0001 ) return false;
	}

	return true;
}

bool
PoseToStructFileRepConverter::use_pdb_info_for_num(
	pose::Pose const & pose,
	Size resnum)
{
	// Setup options.
	pose::PDBInfoCOP pdb_info = pose.pdb_info();
	if (
			pdb_info
			&& !(pdb_info->obsolete())
			&& resnum <= pdb_info->nres()
			&& !( options_.renumber_pdb() ) ) {

		return true;
	} else {

		return false;
	}

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


LinkInformation
PoseToStructFileRepConverter::get_link_record( core::pose::Pose const & pose, core::Size ii, core::Size conn ) {

	using namespace id;
	LinkInformation link;
	Size jj = pose.residue( ii ).connected_residue_at_resconn( conn );
	Size jj_conn = pose.residue( ii ).residue_connection_conn_id( conn );

	link.name1 = pose.residue( ii ).atom_name( pose.residue( ii ).residue_connect_atom_index( conn ) );
	link.resName1 = pose.residue( ii ).name3();
	link.chainID1 = pose.pdb_info()->chain( ii );
	link.resSeq1 = pose.pdb_info()->number( ii );
	link.iCode1 = pose.pdb_info()->icode( ii );
	std::stringstream ss;
	ss.width(6);
	ss << std::right << pose.pdb_info()->number( ii );
	link.resID1 = ss.str() + link.iCode1 + link.chainID1;

	link.name2 =  pose.residue( jj ).atom_name( pose.residue( jj ).residue_connect_atom_index( jj_conn ) );
	link.resName2 = pose.residue( jj ).name3();
	link.chainID2 = pose.pdb_info()->chain( jj );
	link.resSeq2 = pose.pdb_info()->number( jj );
	link.iCode2 = pose.pdb_info()->icode( jj );
	std::stringstream ss2;
	ss2.width(6);
	ss2 << std::right << pose.pdb_info()->number( jj );
	link.resID2 = ss2.str() + link.iCode2 + link.chainID2;

	// Calculate bond distance.
	uint start_atom_index = pose.residue( ii ).atom_index( link.name1 );
	uint stop_atom_index = pose.residue( jj ).atom_index( link.name2 );
	link.length = pose.conformation().bond_length(
		AtomID( start_atom_index, ii ),
		AtomID( stop_atom_index, jj ) );
	return link;
}

SSBondInformation
PoseToStructFileRepConverter::get_ssbond_record( core::pose::Pose const & pose, core::Size ii, core::Size conn ) {

	using namespace id;
	SSBondInformation ssbond;

	Size jj = pose.residue( ii ).connected_residue_at_resconn( conn );

	ssbond.resName1 = pose.residue( ii ).name3();
	ssbond.chainID1 = pose.pdb_info()->chain( ii );
	ssbond.resSeq1 = pose.pdb_info()->number( ii );
	ssbond.iCode1 = pose.pdb_info()->icode( ii );
	std::stringstream ss;
	ss.width(6);
	ss << std::right << pose.pdb_info()->number( ii );
	ssbond.resID1 = ss.str() + ssbond.iCode1 + ssbond.chainID1;

	ssbond.resName2 = pose.residue( jj ).name3();
	ssbond.chainID2 = pose.pdb_info()->chain( jj );
	ssbond.resSeq2 = pose.pdb_info()->number( jj );
	ssbond.iCode2 = pose.pdb_info()->icode( jj );
	std::stringstream ss2;
	ss2.width(6);
	ss2 << std::right << pose.pdb_info()->number( jj );
	ssbond.resID2 = ss2.str() + ssbond.iCode2 + ssbond.chainID2;

	// Calculate bond distance.
	uint start_atom_index = pose.residue( ii ).atom_index( pose.residue( ii ).type().get_disulfide_atom_name() );
	uint stop_atom_index = pose.residue( jj ).atom_index( pose.residue( ii ).type().get_disulfide_atom_name() );
	ssbond.length = pose.conformation().bond_length(
		AtomID( start_atom_index, ii ),
		AtomID( stop_atom_index, jj ) );

	return ssbond;
}

/// @brief Get connectivity annotation information from the Pose object and create LinkInformation and
/// SSBondInformation data as appropriate.
/// @author Watkins
void PoseToStructFileRepConverter::get_connectivity_annotation_info( core::pose::Pose const & pose ) {

	using namespace utility;
	using namespace id;
	using namespace kinematics;
	using namespace conformation;

	// In the past, we walked through the fold_tree and found the termini of all
	// the chemical edges.
	// This technique is not very general because many complicated branched
	// topologies would end up creating cycles and thus aren't included.
	// Similarly, it doesn't permit the use of SSBOND records, as the
	// FoldTree doesn't go through those (except in exotic and generally
	// temporary circumstances.

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {

		if ( pose.residue( ii ).has_lower_connect() ) {
			Size lower = pose.residue( ii ).lower_connect().index();

			// if bonded to not ii - 1 or bonded to not ii - 1's upper
			if ( pose.residue( ii ).connected_residue_at_resconn( lower ) != ii - 1 ||
					pose.residue( ii ).residue_connection_conn_id( lower ) != static_cast<Size>(pose.residue( ii - 1 ).upper_connect().index()) ) {

				vector1<LinkInformation> links;
				LinkInformation link = get_link_record( pose, ii, lower );

				// If key is found in the links map, add this new linkage information to the links already keyed to
				// this residue.
				if ( sfr_->link_map().count(link.resID1) ) {
					links = sfr_->link_map()[link.resID1];
				}
				links.push_back(link);

				sfr_->link_map()[link.resID1] = links;

			}
		}

		if ( pose.residue( ii ).has_upper_connect() ) {
			Size upper = pose.residue( ii ).upper_connect().index();

			// if bonded to not ii + 1 or bonded to not ii + 1's lower
			if ( pose.residue( ii ).connected_residue_at_resconn( upper ) != ii + 1 ||
					pose.residue( ii ).residue_connection_conn_id( upper ) != static_cast<Size>( pose.residue( ii + 1 ).lower_connect().index()) ) {

				// Escape if it's bonded to residue 1's lower--we don't want to double-count cyclization here.
				// If jj < ii, we already counted it
				if ( pose.residue( ii ).connected_residue_at_resconn( upper ) < ii ) continue;

				LinkInformation link = get_link_record( pose, ii, upper );
				vector1<LinkInformation> links;

				// If key is found in the links map, add this new linkage information to the links already keyed to
				// this residue.
				if ( sfr_->link_map().count(link.resID1) ) {
					links = sfr_->link_map()[link.resID1];
				}
				links.push_back(link);

				sfr_->link_map()[link.resID1] = links;

			}
		}

		if ( pose.residue( ii ).n_non_polymeric_residue_connections() == 0 ) continue;

		for ( Size conn = pose.residue( ii ).n_polymeric_residue_connections()+1; conn <= pose.residue( ii ).n_residue_connections(); ++conn ) {

			Size jj = pose.residue( ii ).connected_residue_at_resconn( conn );
			Size jj_conn = pose.residue( ii ).residue_connection_conn_id( conn );

			// Either LINK or SSBOND
			// Note that it's not an SSBOND if you are a cysteine bonded
			// to the UPPER of another cysteine!
			// Are the two atoms both the get_disulfide_atom_name() of the
			// connected residue types?
			if ( ( pose.residue( ii ).type().has_property( chemical::DISULFIDE_BONDED ) || pose.residue( ii ).type().has_property( chemical::SIDECHAIN_THIOL ) ) &&
					( pose.residue( jj ).type().has_property( chemical::DISULFIDE_BONDED ) || pose.residue( jj ).type().has_property( chemical::SIDECHAIN_THIOL ) ) &&
					pose.residue( ii ).residue_connect_atom_index( conn ) ==
					pose.residue( ii ).atom_index( pose.residue( ii ).type().get_disulfide_atom_name() ) &&
					pose.residue( jj ).residue_connect_atom_index( jj_conn ) ==
					pose.residue( ii ).atom_index( pose.residue( jj ).type().get_disulfide_atom_name() ) ) {
				// Disulfide.

				// If jj < ii, we already counted it
				if ( jj < ii ) continue;

				SSBondInformation ssbond = get_ssbond_record( pose, ii, conn );
				vector1<SSBondInformation> ssbonds;

				// If key is found in the links map, add this new linkage information to the links already keyed to
				// this residue.
				if ( sfr_->ssbond_map().count(ssbond.resID1) ) {
					ssbonds = sfr_->ssbond_map()[ssbond.resID1];
				}
				ssbonds.push_back(ssbond);

				sfr_->ssbond_map()[ssbond.resID1] = ssbonds;

			} else {

				LinkInformation link = get_link_record( pose, ii, conn );
				vector1<LinkInformation> links;

				// If key is found in the links map, add this new linkage information to the links already keyed to
				// this residue.
				if ( sfr_->link_map().count(link.resID1) ) {
					links = sfr_->link_map()[link.resID1];
				}
				// If this link is found under the OTHER record...
				bool skip = false;
				if ( sfr_->link_map().count( link.resID2 ) ) {
					for ( Size i = 1; i <= sfr_->link_map()[link.resID2].size(); ++i ) {
						if ( link.resID1 == sfr_->link_map()[link.resID2][i].resID2 ) {
							skip = true;
							break;
						}
					}
				}
				if ( skip )  continue;
				// Make sure it isn't a dupe--for example, make sure that
				// we didn't already push this back as a lower to upper thing.
				// Right now, we assume you can't have two noncanonical connections to the same residue.
				// This is not strictly necessarily true, but until we implement
				// LinkInformation ==, it's good enough.
				bool push_it = true;
				for ( Size i = 1; i <= links.size(); ++i ) {
					if ( ( link.resID1 == links[i].resID1 && link.resID2 == links[i].resID2 )
							|| ( link.resID2 == links[i].resID1 && link.resID1 == links[i].resID2 ) ) {
						push_it = false;
					}
				}
				if ( push_it ) {
					links.push_back(link);
				}

				sfr_->link_map()[link.resID1] = links;
			}
		}
	}
}

/// @brief Get parametric information from the Pose object and add it to the PDB remarks.
///
void PoseToStructFileRepConverter::get_parametric_info(
	core::io::RemarksOP remarks,
	core::pose::Pose const & pose
) {

	using namespace core::conformation::parametric;

	core::Size const nsets(pose.conformation().n_parameters_sets()); //How many ParametersSet objects are there in the pose?
	if ( nsets==0 ) return; //No need to proceed if this isn't a parametric conformation.

	for ( core::Size iset=1; iset<=nsets; ++iset ) { //Loop through all of the ParametersSet objects.
		ParametersSetCOP curset( pose.conformation().parameters_set(iset) );
		std::stringstream curset_summary;
		curset->get_pdb_remark( curset_summary );

		int cur_remark_number(1); //int instead of core::Size, to match the Remarks class.
		if ( remarks->size() > 0 ) cur_remark_number = (*sfr_->remarks())[remarks->size()-1].num + 1;

		std::string cur_remark_str;
		while ( std::getline(curset_summary,cur_remark_str) ) {
			core::io::RemarkInfo cur_remark;
			cur_remark.value = cur_remark_str; // Ugh.  The RemarkInfo class provides no getters or setters -- only public access to its members.
			cur_remark.num=cur_remark_number;
			remarks->push_back( cur_remark );
			cur_remark_str.clear();
		}
	}

	return;
}

/// @brief Get additional pose data.
/// @details This is rewritten from the old core/io/pdb/file_data.cc:write_additional_pdb_data() function.
/// This was a catch-all for dumping out all sorts of protocol-specific stuff (e.g. membrane info).
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
PoseToStructFileRepConverter::grab_additional_pose_data( core::pose::Pose const &pose )
{
	grab_membrane_info( pose, options_.normalize_to_thk() );
	grab_foldtree( pose, options_.fold_tree_io() );
	grab_pdb_parents( pose, options_.pdb_parents() );
	grab_pdb_comments( pose, options_.pdb_comments() );
	grab_torsion_records( pose, options_.output_torsions() );
	grab_pdbinfo_labels( pose );
}





/// @brief Set whether to write the fold tree, in the
/// StructFileRepOptions object (options_).
void
PoseToStructFileRepConverter::set_fold_tree_io(
	bool const setting
) {
	options_.set_fold_tree_io( setting );
}

/********* PRIVATE FUNCTIONS ************/

/// @brief Get the membrane information from the pose and store it in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
/// @param [in] pose The pose.
/// @param [in] normalize_to_thk Normalized MEM lines, useful for visualizing the boundaries
/// of the membrane by coupling the NORM and THK coordinates.  Added by Rebecca.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
PoseToStructFileRepConverter::grab_membrane_info(
	core::pose::Pose const &pose,
	bool const normalize_to_thk
) {

	if ( pose.conformation().is_membrane() && normalize_to_thk == true ) {

		// Grab membrane residue & current data
		core::Real const thkn( pose.conformation().membrane_info()->membrane_thickness() );
		core::Vector const cntr( pose.conformation().membrane_info()->membrane_center(pose.conformation()) );
		core::Vector norm( pose.conformation().membrane_info()->membrane_normal(pose.conformation()) );

		// Actually normalize the membrane residue to thk
		norm.normalize();
		norm.x( norm.x() * thkn);
		norm.y( norm.y() * thkn);
		norm.z( norm.z() * thkn);

		// Get rescount, chain count
		core::Size const chaincount( sfr_->chains().size() - 1 );
		core::Size const rescount( sfr_->chains()[chaincount][ sfr_->chains()[chaincount].size() - 1 ].serial );
		char const newchain( sfr_->chains()[chaincount][ sfr_->chains()[chaincount].size() - 1 ].chainID + 1 );

		AtomInformation ai1, ai2, ai3;
		ai1.isHet = true;
		ai2.isHet = true;
		ai3.isHet = true;
		ai1.serial = rescount + 1;
		ai2.serial = rescount + 2;
		ai3.serial = rescount + 3;
		ai1.name = "THKN";
		ai2.name = "CNTR";
		ai3.name = "NORM";
		ai1.resName = "MEM";
		ai2.resName = "MEM";
		ai3.resName = "MEM";
		ai1.chainID = newchain;
		ai2.chainID = newchain;
		ai3.chainID = newchain;
		ai1.resSeq = rescount + 1;
		ai2.resSeq = rescount + 2;
		ai3.resSeq = rescount + 3;
		ai1.x = thkn; ai1.y = 0; ai1.z = 0;
		ai2.x = cntr.x(); ai2.y = cntr.y(); ai2.z = cntr.z();
		ai3.x = norm.x(); ai3.y = norm.y(); ai3.z = norm.z();
		ai1.occupancy=1.0;
		ai2.occupancy=1.0;
		ai3.occupancy=1.0;
		ai1.element=" H";
		ai2.element=" H";
		ai3.element=" H";

		sfr_->chains().push_back( utility::vector0< AtomInformation >() );
		sfr_->chains()[ chaincount + 1 ].push_back( ai1 );
		sfr_->chains()[ chaincount + 1 ].push_back( ai2 );
		sfr_->chains()[ chaincount + 1 ].push_back( ai3 );

		//out << "HETATM XXXX THKN MEM " << new_chain << I(4,resid+1) << "    " << F(8, 3, thkn) << "   0.000   0.000 \n";
		//out << "HETATM XXXX CNTR MEM " << new_chain << I(4,resid+1) << "    " << F(8, 3, cntr.x()) << F(8, 3, cntr.y()) << F(8, 3, cntr.z()) << "\n";
		//out << "HETATM XXXX NORM MEM " << new_chain << I(4,resid+1) << "    " << F(8, 3, norm.x()) << F(8, 3, norm.y()) << F(8, 3, norm.z()) << "\n";

		//              HETATM XXXX THKN MEM X  81      15.000   0.000   0.000
		//    HETATM XXXX CNTR MEM X  81       0.000   0.000   0.000
		//              HETATM XXXX NORM MEM X  81       0.000   0.000  15.000

	}

}

/// @brief Get the conect record information from the pose and store it in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
/// @param [in] pose The pose, for determining what's bonded to what.
/// @param [in] atom_index_in_pose The index (in the whole pose) of the atom that we're setting up.
/// @param [in] res_index The current residue index.
/// @param [in] atom_index_rsd The index (in the residue) of the atom that we're setting up.
/// @param [in] dist_cutoff The atom separation, above which a CONECT record is written.
/// @param [in] write_virtuals Are virtual atoms being written out?
/// @param [in/out] ai The AtomInformation object that's being set up.  Modified by this function
void
PoseToStructFileRepConverter::grab_conect_records_for_atom(
	core::pose::Pose const &pose,
	core::Size const res_index,
	core::Size const atom_index_in_rsd,
	core::io::AtomInformation &ai
) const {

	if ( options_.skip_connect_info() ) return;
	if ( pose.n_residue() == 0 ) return; //Probably unnecesary, but why not?

	debug_assert( res_index <= pose.n_residue() );
	debug_assert( atom_index_in_rsd <= pose.residue( res_index ).natoms() );

	bool const write_virtuals( options_.output_virtual() );

	//Return if this is virtual and we're not writing virtuals:
	if ( ( pose.residue(res_index).is_virtual_residue() || pose.residue(res_index).atom_type(atom_index_in_rsd).is_virtual() ) && !write_virtuals ) return;

	bool const writeall( options_.write_all_connect_info() );
	bool const this_res_is_canonical_or_solvent( pose.residue(res_index).type().is_canonical() || pose.residue(res_index).type().is_solvent() );
	core::Real const dist_cutoff_sq( options_.connect_info_cutoff()*options_.connect_info_cutoff() );
	core::Size count(0);
	core::id::AtomID const this_atom_id( atom_index_in_rsd, res_index ); //The AtomID of this atom.
	utility::vector1<core::id::AtomID> const & bonded_ids(  pose.conformation().bonded_neighbor_all_res( this_atom_id, write_virtuals ) ); //List of AtomIDs of atoms bound to this atom.

	for ( core::Size i=1; i<=res_index; ++i ) { //Loop through all residues
		for ( core::Size j=1, jmax=pose.residue(i).natoms(); j<=jmax; ++j ) { //Loop through all atoms in this residue
			if ( i == res_index && j == atom_index_in_rsd ) {
				break;
			}

			if ( !write_virtuals && pose.residue(i).atom_type(j).is_virtual() ) continue; //Skip writing virtuals if we should do so.
			++count;
			if ( !writeall ) {
				if ( i == res_index ) { //If this is an intra-residue bond:
					if ( pose.residue(i).type().is_canonical() || pose.residue(i).type().is_solvent() ) continue; //Skip canonical and solvent intra-res bonds.
				} else { //If this is an inter-residue bond:
					if ( this_res_is_canonical_or_solvent && ( pose.residue(i).type().is_canonical() || pose.residue(i).type().is_solvent() ) ) continue; //Skip canonical or solvent inter-res bonds.
				}
			}

			core::id::AtomID const other_atom_id( j, i ); //Candidate other atom to which this one might be bonded.
			for ( core::Size n=1, nmax=bonded_ids.size(); n<=nmax; ++n ) {
				if ( bonded_ids[n] == other_atom_id ) { //If the candidate atom is in the list of bonded atoms for this atom, add CONECT data.
					if ( pose.xyz( this_atom_id ).distance_squared( pose.xyz( other_atom_id ) ) >= dist_cutoff_sq ) { //Are we over the distance cutoff?
						ai.connected_indices.push_back( count );
					}
					break;
				}
			} //Loop through bonded atoms

		} //Loop through all atoms
	} //Loop through all residues

}

/// @brief Get the foldtree from the pose and store it in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
void
PoseToStructFileRepConverter::grab_foldtree(
	core::pose::Pose const &pose,
	bool const output_foldtree
) {
	if ( !output_foldtree ) return;
	std::stringstream sstr;
	sstr << pose.fold_tree();
	sfr_->foldtree_string() = sstr.str();
}

/// @brief Get the pdb parents from the pose and store it in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
void
PoseToStructFileRepConverter::grab_pdb_parents(
	core::pose::Pose const &pose,
	bool const output_parents
) {
	if ( !output_parents ) return;
	std::string value;
	bool has_parents = core::pose::get_comment( pose, "parents", value );
	if ( has_parents ) {
		sfr_->additional_string_output() = sfr_->additional_string_output() + "REMARK PARENT    " + value.substr(0,5) + "\n";
	}
}

/// @brief Get the pdb comments from the pose and store it in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
void
PoseToStructFileRepConverter::grab_pdb_comments(
	core::pose::Pose const &pose,
	bool const output_comments
) {
	if ( !output_comments ) return;
	sfr_->pdb_comments() = core::pose::get_all_comments( pose );
}

/// @brief Get the torsion information from the pose and store it in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
/// @details This is an incredibly ugly way to handle the output of these
/// records.  This abuses the REMARK lines in the PDB format, and should be
/// refactored.  Keeping as-is for now to preserve legacy behaviour (VKM
/// 29 January 2016, Chemical XRW).
void
PoseToStructFileRepConverter::grab_torsion_records(
	core::pose::Pose const &pose,
	bool const output_torsions
) {
	using namespace ObjexxFCL::format;

	if ( !output_torsions ) return;

	std::stringstream out("");
	core::pose::PDBInfoCOP pdb_info(pose.pdb_info());

	if ( !core::pose::is_ideal_pose(pose) ) {
		TR << "Ignoring out::file::output_torsions option because pose is non-ideal!" << std::endl;
	} else {
		ObjexxFCL::FArray1D_char dssp_reduced_secstruct(pose.n_residue());
		scoring::dssp::Dssp(pose).dssp_reduced(dssp_reduced_secstruct);
		if ( pdb_info ) {
			out << "REMARK torsions: res pdbres pdbchain seq dssp phi psi omega" << std::endl;
		} else {
			out << "REMARK torsions: res    res    chain seq dssp phi psi omega" << std::endl;
		}
		for ( core::Size i=1; i<=pose.n_residue(); ++i ) {
			if ( pdb_info ) {
				out << "REMARK " << I( 4, i ) << " " << I( 4, pose.pdb_info()->number(i)) << " " << pose.pdb_info()->chain(i) << " " << pose.residue( i ).name1() << " " <<
					dssp_reduced_secstruct(i) << " " << F( 9, 3, pose.phi(i)) << " " << F( 9, 3, pose.psi(i)) << " " << F( 9, 3, pose.omega(i)) << std::endl;
			} else {
				out << "REMARK " << I( 4, i ) << " " << I( 4, i) << " " << pose.chain(i) << " " << pose.residue( i ).name1() << " " <<
					dssp_reduced_secstruct(i) << " " << F( 9, 3, pose.phi(i)) << " " << F( 9, 3, pose.psi(i)) << " " << F( 9, 3, pose.omega(i)) << std::endl;
			}
		}
	}

	//Append the results of all of this to the additional_string_output in the StructFileRep:
	sfr_->additional_string_output() = sfr_->additional_string_output() + out.str();
}


/// @brief Get the pdbinfo labels from the pose and store them in the
/// StructFileRep for output to pdbs/mmCIF/whatnot.
/// @details This also abuses REMARK lines, and should be rewritten.
void
PoseToStructFileRepConverter::grab_pdbinfo_labels(
	core::pose::Pose const &pose
) {
	using namespace ObjexxFCL::format;

	// Added by Daniel-Adriano Silva, used to write the PDBInfoLabels to the REMARK
	// First test that the pdb_info() is not empty
	if ( pose.pdb_info() ) {
		// Then output the labels
		std::stringstream out;
		for ( core::Size i=1; i<=pose.n_residue(); ++i ) {
			utility::vector1 < std::string > const tmp_v_reslabels( pose.pdb_info()->get_reslabels(i) ); //Ugh.  Passing a vector of strings by return value.
			core::Size const numLables( tmp_v_reslabels.size() );
			//Only write if the residue has any label (keep the file as small as possible)
			if ( numLables > 0 ) {
				out << "REMARK PDBinfo-LABEL: " << I( 4, i );
				for ( core::Size lndx=1; lndx <= numLables; ++lndx ) {
					out << " " << tmp_v_reslabels[lndx];
				}
				out << std::endl;
			}
		}
		sfr_->additional_string_output() = sfr_->additional_string_output() + out.str();
	}
}

/// @brief Get the total number of atoms in the SFR.
///
core::Size
PoseToStructFileRepConverter::total_sfr_atoms(
	StructFileRep const & sfr
) const {
	core::Size total = 0;
	for ( core::Size chain = 0, chain_max=sfr.chains().size(); chain < chain_max ; ++chain ) {
		total += total_sfr_atoms( sfr, chain );
	}
	return total;
}

/// @brief Get the total number of atoms in a chain in the SFR.
///
core::Size
PoseToStructFileRepConverter::total_sfr_atoms(
	StructFileRep const & sfr,
	core::Size const chain_num
) const {
	runtime_assert_string_msg( chain_num < sfr.chains().size(), "Error in core::io::pose_to_sfr::PoseToStructFileRepConverter::total_sfr_atoms(): The chain index is out of range.");
	return sfr.chains()[ chain_num ].size();
}

/// @brief Return the PDB resName, chainID, resSeq, and iCode for the given Rosetta sequence position.
/// @details Output is res_info.
void
PoseToStructFileRepConverter::get_residue_information(
	core::pose::Pose const & pose,
	core::uint const seqpos,
	bool const use_PDB,
	bool const renumber_chains,
	core::Size const new_tercount,
	ResidueInformation &res_info
) const {
	using namespace std;
	using namespace utility;
	using namespace core;
	using namespace core::pose;
	using core::chemical::chr_chains;

	res_info.resName( pose.residue(seqpos).name3() );

	// Use PDB-specific information?
	if ( use_PDB ) {
		PDBInfoCOP pdb_info = pose.pdb_info();

		res_info.chainID( pdb_info->chain(seqpos) );
		if ( res_info.chainID() == PDBInfo::empty_record() ) {  // safety
			TR.Warning << "PDBInfo chain ID was left as character '" << PDBInfo::empty_record()
				<< "', denoting an empty record; for convenience, replacing with space." << endl;
			res_info.chainID( ' ' );
		}
		res_info.resSeq(    pdb_info->number(seqpos));
		res_info.iCode(     pdb_info->icode(seqpos));
		res_info.segmentID( pdb_info->segmentID(seqpos));
		res_info.terCount( new_tercount );
		// ...or not?
	} else {
		uint chain_num = pose.chain(seqpos);
		runtime_assert(chain_num > 0);

		res_info.chainID( chr_chains[(chain_num - 1) % chr_chains.size()] );
		res_info.resSeq(  seqpos );
		res_info.iCode( ' ' );
		res_info.terCount( new_tercount );

		// If option is specified, renumber per-chain.
		if ( renumber_chains ) {
			vector1<uint> const & chn_ends = pose.conformation().chain_endings();
			for ( uint i = 1; i <= chn_ends.size(); ++i ) {
				if ( chn_ends[i] < seqpos ) {
					res_info.resSeq( seqpos - chn_ends[i] );
				}
			}
		}

		// Fix for >10k residues.
		res_info.resSeq( res_info.resSeq() % 10000 );
	}
}

} // namespace pose_to_sfr
} // namespace io
} // namespace core
