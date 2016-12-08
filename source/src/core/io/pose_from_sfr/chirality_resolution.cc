// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/pose_from_sfr/chirality_resolution.cc
/// @brief
/// @author Rhiju Das
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Unit header
#include <core/io/pose_from_sfr/chirality_resolution.hh>

// Package headers
#include <core/io/ResidueInformation.hh>
#include <core/io/AtomInformation.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>
#include <numeric/numeric.functions.hh>

// External headers
#include <ObjexxFCL/format.hh>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/mcgregor_common_subgraphs.hpp>

// C++ headers
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <utility>

using namespace basic::options;
using namespace basic::options::OptionKeys;

// Tracer instance for this file
static THREAD_LOCAL basic::Tracer TR( "core.io.pose_from_sfr.chirality_resolution" );

////////////////////////////////////////////////////////////////////////////////////////////
// useful utility scripts for fixing up residue and atom names.
//
// Originally implemented for nucleic acids, to permit backward-compatibility with
// old Rosetta atom & residue names.
//
// Some  useful functions below for fixing up H5'<-->H5'' ambiguities based on geometric
//  comparisons to ideal coordinates.
//
//
namespace core {
namespace io {
namespace pose_from_sfr {


//////////////////////////////////////////////////////////////////////////////////////////
// @brief due to differences in different crystallography/NMR/modeling packages, labeling of sister atoms
//  (like OP1 <--> OP2, or H41 <--> H42) in PDBs is totally wacky. This is an attempt to regularize...
//  and it can actually make a difference since sometimes partial charges on sister hydrogens
//  can be different. Right now only set up for nucleic acids, but could probably generalize.
void
check_and_correct_sister_atoms( core::conformation::ResidueOP & rsd ){

	// in all nucleic acids
	check_and_correct_sister_atom_based_on_chirality( rsd, " OP1", " OP2", " P  ", " O5'" );
	check_and_correct_sister_atom_based_on_chirality( rsd, " H5'", "H5''", " C5'", " C4'" );
	// in DNA
	// AMW:
	// REMOVED because we no longer blanket apply the same
	// 2H2' => H2'' mapping. Do we still have to do this?
	//check_and_correct_sister_atom_based_on_chirality( rsd, " H2'", "H2''", " C2'", " C3'" );

	// in adenosine
	check_and_correct_sister_atom_based_on_outgroup( rsd, " H61", " H62", " N1 " );
	// in guanosine
	check_and_correct_sister_atom_based_on_outgroup( rsd, " H21", " H22", " N1 " );
	// in cytidine
	check_and_correct_sister_atom_based_on_outgroup( rsd, " H41", " H42", " N3 " );
}


//////////////////////////////////////////////////////////////////////////////////////////
// sisters sprout off the same parent, and outer_ref is something else bonded to the parent. a cousin, i guess.
void
check_and_correct_sister_atom_based_on_chirality( core::conformation::ResidueOP & rsd,
	std::string const & sister1_name,
	std::string const & sister2_name,
	std::string const & parent_name,
	std::string const & cousin_name ){

	if ( !rsd->has( sister1_name ) ) return;
	if ( !rsd->has( sister2_name ) ) return;
	if ( !rsd->has( parent_name ) ) return;
	if ( !rsd->has( cousin_name ) ) return;

	Vector const current_xyz_sister1        = rsd->xyz( sister1_name );
	Vector const current_xyz_sister2        = rsd->xyz( sister2_name );
	Vector const current_xyz_parent         = rsd->xyz( parent_name );
	Vector const current_xyz_cousin      = rsd->xyz( cousin_name );
	int current_sign = get_chirality_sign( current_xyz_sister1, current_xyz_sister2, current_xyz_parent, current_xyz_cousin );

	core::chemical::ResidueType const & rsd_type = rsd->type();
	Vector const ideal_xyz_sister1        = rsd_type.atom( sister1_name ).ideal_xyz();
	Vector const ideal_xyz_sister2        = rsd_type.atom( sister2_name ).ideal_xyz();
	Vector const ideal_xyz_parent         = rsd_type.atom( parent_name ).ideal_xyz();
	Vector const ideal_xyz_cousin      = rsd_type.atom( cousin_name ).ideal_xyz();
	int ideal_sign = get_chirality_sign( ideal_xyz_sister1, ideal_xyz_sister2, ideal_xyz_parent, ideal_xyz_cousin );

	if ( current_sign != ideal_sign ) flip_atom_xyz( rsd, sister1_name, sister2_name );

}

//////////////////////////////////////////////////////////////////////////////////
void
check_and_correct_sister_atom_based_on_outgroup( core::conformation::ResidueOP & rsd,
	std::string const & sister1_name,
	std::string const & sister2_name,
	std::string const & outgroup_name ){

	if ( !rsd->has( sister1_name ) ) return;
	if ( !rsd->has( sister2_name ) ) return;
	if ( !rsd->has( outgroup_name ) ) return;

	Vector const current_xyz_sister1        = rsd->xyz( sister1_name );
	Vector const current_xyz_sister2        = rsd->xyz( sister2_name );
	Vector const current_xyz_outgroup       = rsd->xyz( outgroup_name );

	int current_closest_sister = get_closest_sister( current_xyz_sister1, current_xyz_sister2, current_xyz_outgroup );

	core::chemical::ResidueType const & rsd_type = rsd->type();
	Vector const ideal_xyz_sister1        = rsd_type.atom( sister1_name ).ideal_xyz();
	Vector const ideal_xyz_sister2        = rsd_type.atom( sister2_name ).ideal_xyz();
	Vector const ideal_xyz_outgroup       = rsd_type.atom( outgroup_name ).ideal_xyz();
	int ideal_closest_sister = get_closest_sister( ideal_xyz_sister1, ideal_xyz_sister2, ideal_xyz_outgroup );

	if ( current_closest_sister != ideal_closest_sister ) flip_atom_xyz( rsd, sister1_name, sister2_name );

}

//////////////////////////////////////////////////////////////////////////////////
void
flip_atom_xyz( core::conformation::ResidueOP & rsd,
	std::string const & sister1_name,
	std::string const & sister2_name ) {
	// following is to show warnings or cap number.
	static Size nfix( 0 );
	static Size const max_fix( 2 );
	static bool const show_all_fixup( option[ in::show_all_fixes ]() );
	static bool showed_warning( false );

	if ( ++nfix <= max_fix || show_all_fixup ) {
		TR << "Flipping atom xyz for " << sister1_name << " and " << sister2_name << " for residue " << rsd->name3() <<  std::endl;
	}
	Vector const temp_xyz = rsd->xyz( sister1_name );
	rsd->set_xyz( sister1_name, rsd->xyz( sister2_name ) );
	rsd->set_xyz( sister2_name, temp_xyz );

	if ( nfix > max_fix && !show_all_fixup && !showed_warning ) {
		TR << "Number of flip-atom fixups exceeds output limit. Rerun with -show_all_fixes to show everything." << std::endl;
		showed_warning = true;
	}
}

//////////////////////////////////////////////////////////////////////////////////
int
sgn( Real const & x ){
	return ( x > 0 ) - ( x < 0 );
}

//////////////////////////////////////////////////////////////////////////////////
int
get_chirality_sign(  Vector const & xyz_sister1,
	Vector const & xyz_sister2,
	Vector const & xyz_parent,
	Vector const & xyz_cousin ) {
	int const sign = sgn( dot( xyz_cousin - xyz_parent, cross( xyz_sister1 - xyz_parent, xyz_sister2 - xyz_parent ) ) );
	if ( sign == 0 ) utility_exit_with_message( "unexpected sign error when checking chirality" );
	return sign;
}

//////////////////////////////////////////////////////////////////////////////////
// returns 1 or 2 based on which sister is closest to outgroup.
int
get_closest_sister(  Vector const & xyz_sister1,
	Vector const & xyz_sister2,
	Vector const & xyz_outgroup ) {
	return ( xyz_sister1.distance( xyz_outgroup ) < xyz_sister2.distance( xyz_outgroup) ) ? 1 : 2;
}


//////////////////////////////////////////////////////////////////////////////////
// Utility functions for name remapping.


/// @brief Get theshold distance below which two atoms are considered bonded. (1.2*covalent)
/// @details The closest distance of a non-bonded contact is likely to be something
/// like the opposite atoms cyclobutane. This would be sqrt(2)*covalent bond distance.
/// We thus set the contact distance threshold to 1.2*covalent bond distance to allow
/// bond length flexibility. In initial tests this looks to give a clean decision.
/// @details Pass-by-value is deliberate, as we want to strip the elements of whitespace
core::Real
bonding_distance_threshold( std::string element1, std::string element2 ) {
	core::chemical::ElementSetCOP element_types = core::chemical::ChemicalManager::get_instance()->element_set("default");

	utility::strip_whitespace( element1 );
	utility::strip_whitespace( element2 );
	if ( ! element_types->contains_element_type( element1 ) ||
			! element_types->contains_element_type( element2 ) ) {
		utility_exit_with_message("Cannot find bonding distance threshold for elements '"+element1+"' and '"+element2+"'");
	}
	core::chemical::ElementCOP elem1( element_types->element(element1) );
	core::chemical::ElementCOP elem2( element_types->element(element2) );
	core::Real rad1( elem1->get_property(  core::chemical::Element::CovalentRadius ) );
	core::Real rad2( elem2->get_property(  core::chemical::Element::CovalentRadius ) );
	return 1.2 * (rad1 + rad2) ;
}

core::Real score_mapping( NameBimap const & mapping, ResidueInformation const & rinfo, chemical::ResidueType const & rsd_type ) {
	// Scoring:
	// -5 point for each chirality mis-match
	// -2 points for each bond that shouldn't be made
	// 1/natoms point for each atom present
	// 1/natoms^2 point for each name match
	//
	// This scoring scheme is set up such that we sort primarily on certain features,
	// and tie-break on the next feature.

	core::Size const natoms = rsd_type.natoms();
	core::Size const natoms2 = natoms*natoms;
	//Inverse mapping - from Pose names to rinfo names
	NameBimap::right_map const & inverse( mapping.right );

	// Sum the scores for each of the atoms individually
	core::Real score(0);
	for ( core::Size ii(1); ii <= natoms; ++ii ) {
		std::string const & name( rsd_type.atom_name(ii) );
		if ( inverse.count( name ) == 0 ) {
			continue;
		}
		score += 1.0/natoms; // 1/natoms point for each atom present.
		if ( utility::stripped_whitespace(name) == utility::stripped_whitespace( inverse.find(name)->second ) ) {
			score += 1.0 / natoms2; // 1/natoms^2 point for each name identity match.
		}
		// Check the chirality of the atoms
		core::chemical::AtomIndices const & nbrs( rsd_type.bonded_neighbor(ii) );
		// Make a name mapping for the bonded atoms, if present.
		utility::vector1< std::string > type_nbrs;
		for ( core::Size jj(1); jj <= nbrs.size(); ++jj ) {
			std::string const & posename( rsd_type.atom_name(nbrs[jj]) );
			if ( inverse.count( posename ) &&
					rinfo.xyz().count( inverse.find(posename)->second ) ) {
				// Only look at neighbors with correspondences and coordinates
				type_nbrs.push_back( posename );
			}
		}
		// Make sure we have enough information to compute the chirality
		if ( type_nbrs.size() == 3 ) {
			// We can compute chirality with three atoms if we use the current atom as a reference
			type_nbrs.push_back( name );
		}
		if ( type_nbrs.size() >= 4 ) {
			// We hope that the first four are sufficient for a chirality signature if there's more than 4 bonded neighbors
			core::Real const rsd_dhd = numeric::dihedral_degrees(
				rsd_type.atom( type_nbrs[1] ).ideal_xyz(),
				rsd_type.atom( type_nbrs[2] ).ideal_xyz(),
				rsd_type.atom( type_nbrs[3] ).ideal_xyz(),
				rsd_type.atom( type_nbrs[4] ).ideal_xyz() );
			//These will be found, because we only pushed back names which were present and had coordinates.
			core::Real const rinfo_dhd = numeric::dihedral_degrees(
				rinfo.xyz().find( inverse.find(type_nbrs[1])->second )->second,
				rinfo.xyz().find( inverse.find(type_nbrs[2])->second )->second,
				rinfo.xyz().find( inverse.find(type_nbrs[3])->second )->second,
				rinfo.xyz().find( inverse.find(type_nbrs[4])->second )->second );
			//TR << "Dihedral difference: " << rsd_dhd - rinfo_dhd << std::endl;
			if ( numeric::sign( rsd_dhd ) != numeric::sign( rinfo_dhd ) && // The chiral orientation is different,
					numeric::abs_difference( rsd_dhd, rinfo_dhd ) > 10.0 ) { // and they aren't insignificantly close to each other (e.g. planar).
				score -= 5.0; // -5.0 point for each chirality mis-match.
			}
		} // else too few bonded neighbors to do chirality checks.
	}

	// For each bond in the residue type, make sure that the mapped atoms are "close enough" to also be bonded.
	// This check will keep the heurisitic from smashing together distal parts of the molecule
	core::chemical::EIter bonditr, bonditr_end;
	for ( boost::tie( bonditr, bonditr_end ) = rsd_type.bond_iterators(); bonditr != bonditr_end; ++bonditr ) {
		core::chemical::VD source( boost::source(*bonditr,rsd_type.graph()) ), target( boost::target(*bonditr,rsd_type.graph()) );
		std::string const & name1( rsd_type.atom_name(source) );
		std::string const & name2( rsd_type.atom_name(target) );
		if ( inverse.count(name1) && inverse.count(name2) &&
				rinfo.xyz().count( inverse.find(name1)->second ) &&
				rinfo.xyz().count( inverse.find(name2)->second ) ) {
			Vector const & pos1(rinfo.xyz().find( inverse.find(name1)->second )->second);
			Vector const & pos2(rinfo.xyz().find( inverse.find(name2)->second )->second);
			std::string const & elem1( rsd_type.atom(source).element_type()->get_chemical_symbol() ) ;
			std::string const & elem2( rsd_type.atom(target).element_type()->get_chemical_symbol() ) ;
			core::Real bond_thresh( bonding_distance_threshold(elem1,elem2) );
			if ( pos1.distance_squared(pos2) > (bond_thresh*bond_thresh) ) {
				score -= 2.0;
			}
		}
	}

	return score;
}

typedef boost::undirected_graph<AtomInformation /*Node information only*/ > AtomInfoGraph;
typedef AtomInfoGraph::vertex_descriptor AIVD;

class GeometricRenameIsomorphismCallback {
public:
	GeometricRenameIsomorphismCallback(AtomInfoGraph const & aigraph,
		ResidueInformation const & rinfo,
		core::chemical::ResidueType const & rsdtype,
		NameBimap & mapping,
		core::Real &  mapscore ):
		aigraph_( aigraph ),
		rinfo_(rinfo),
		rsdtype_( rsdtype ),
		mapping_( mapping ),
		best_score_( mapscore ),
		n_mappings_( new core::Size(0) )
	{
	}

	template< class VD2VDmap1, class VD2VDmap2 >
	bool operator()(VD2VDmap1 map_1_to_2, VD2VDmap2, core::Size = 0 /*needed for McGregor*/ ) {
		NameBimap newmap;
		AtomInfoGraph::vertex_iterator iter, iter_end;
		for ( boost::tie( iter, iter_end) = boost::vertices(aigraph_); iter != iter_end; ++iter ) {
			if ( map_1_to_2[ *iter ] != core::chemical::ResidueGraph::null_vertex() ) {
				newmap.insert( NameBimap::value_type( aigraph_[*iter].name , rsdtype_.atom_name( map_1_to_2[ *iter ] ) ) );
			}
		}
		core::Real newscore( score_mapping( newmap, rinfo_, rsdtype_ ) );
		++(*n_mappings_);
		if ( *n_mappings_ % 1000 == 0 ) {
			TR.Debug << "Mapping " << *n_mappings_ << " has a score of " << newscore << std::endl;
		}
		if ( newscore > best_score_ ) {
			if ( TR.Trace.visible() ) {
				TR.Trace << "Found new best mapping, with score " << newscore << std::endl;
			}
			best_score_ = newscore;
			mapping_ = newmap;
		}
		// We don't want to go forever - if there are too many mappings, just go with a "good enough" one.
		// RM: The "too many" criteria here are arbitrary, based on wait times on my machine.
		core::Size natoms( numeric::min(rsdtype_.natoms(),boost::num_vertices(aigraph_)) );
		if ( (*n_mappings_ >=  100000 && mapping_.left.size() > (3*natoms)/4 ) ||
				(*n_mappings_ >=  300000 && mapping_.left.size() >   natoms/2 ) ||
				(*n_mappings_ >= 1000000 && mapping_.left.size() >= 3 ) ) {
			TR << "Too many mappings to consider in geometric name remapping - truncating search after " << *n_mappings_
				<< " mappings, with " << mapping_.left.size() << " of " << natoms << " possible atoms matched." << std::endl;
			return false;
		}
		return true;
	}

private:
	AtomInfoGraph const & aigraph_;
	ResidueInformation const & rinfo_;
	core::chemical::ResidueType const & rsdtype_;
	NameBimap & mapping_;
	core::Real & best_score_;
	// This needs to be a shared pointer as the callback is passed around by value
	boost::shared_ptr< core::Size > n_mappings_;
};

/// @brief Will consider two verticies equivalent if they have the same element.
class GeometricRenameVerticiesEquivalent {
public:

	GeometricRenameVerticiesEquivalent(AtomInfoGraph const & aigraph,
		core::chemical::ResidueGraph const & rsdtype ):
		aigraph_( aigraph ),
		rsdtype_( rsdtype )
	{}

	bool operator() ( AIVD vd1, core::chemical::VD vd2 ) {
		std::string pdb_elem( aigraph_[vd1].element );
		utility::strip_whitespace( pdb_elem );
		utility::uppercase( pdb_elem );
		std::string rsdtype_elem( rsdtype_[vd2].element_type()->get_chemical_symbol() );
		utility::strip_whitespace( rsdtype_elem );
		utility::uppercase( rsdtype_elem );
		//TR << "Element match '" << pdb_elem << "' '" << rsdtype_elem << "'" << std::endl;
		return pdb_elem == rsdtype_elem;
	}

private:
	AtomInfoGraph const & aigraph_;
	core::chemical::ResidueGraph const & rsdtype_;
};

/// @brief Attempt to use element identity and connectivity to map atom names from the rinfo object onto the rsd_type object names.
void
remap_names_on_geometry( NameBimap & mapping,
	ResidueInformation const & rinfo,
	chemical::ResidueType const & rsd_type) {
	// Set up the graph on the PDB side
	AtomInfoGraph aigraph;
	std::map< std::string, AIVD > name_aivd_map;
	utility::vector1< AIVD > small_order;

	for ( utility::vector1< AtomInformation >::const_iterator iter=rinfo.atoms().begin(), iter_end=rinfo.atoms().end();
			iter != iter_end; ++iter ) {
		if ( ! rinfo.xyz().count( iter->name ) ) {
			continue; // Only look at atoms with coordinates
		}
		AIVD aivd = boost::add_vertex(*iter,aigraph);
		//Fix up the element name, if we have an old-style short ATOM line.
		//Atom info in the graph is a copy - so we can modify it if we need.
		std::string & elem( aigraph[ aivd ].element );
		if ( elem == "  " ||  elem == "" ) {
			aigraph[ aivd ].element = aigraph[ aivd ].name.substr(0,2); // First two letters of the atom name - typically this is the element code.
			// AMW: WHAT?! It could also be a number! (e.g. 1H in NtermProteinFull)
			if ( isdigit( int( aigraph[ aivd ].name.at(0) ) ) ) {
				aigraph[ aivd ].element = aigraph[ aivd ].name.substr(1,1);
			}

		}
		name_aivd_map[ iter->name ] = aivd;
		small_order.push_back( aivd );
	}

	// Add bonds to all atoms which are within the contact distance of each other.
	for ( core::Size ii(1); ii <= rinfo.atoms().size(); ++ii ) {
		std::string const & atomname( rinfo.atoms()[ii].name );
		Vector const & atomxyz( rinfo.xyz().find( atomname )->second );
		// Need to pull the element out of the graph, as we may have adjusted it above
		std::string const & elem1( aigraph[ name_aivd_map[ atomname ] ].element );
		for ( core::Size jj(ii+1); jj <= rinfo.atoms().size(); ++jj ) {
			std::string const & atom2name( rinfo.atoms()[jj].name );
			Vector const &atom2xyz( rinfo.xyz().find( atom2name )->second );
			std::string const & elem2( aigraph[ name_aivd_map[ atom2name ] ].element );
			core::Real bond_thresh( bonding_distance_threshold(elem1,elem2) );
			if ( atomxyz.distance_squared( atom2xyz ) < (bond_thresh*bond_thresh) ) {
				boost::add_edge( name_aivd_map[atomname], name_aivd_map[atom2name], aigraph);
			}
		}
	}

	TR.Debug << "Graph sizes: Ainfo " << boost::num_vertices(aigraph)
		<< " ResidueType " <<  boost::num_vertices(rsd_type.graph()) << std::endl;
	TR.Debug << "      sizes: order " << small_order.size() << std::endl;
	TR.Debug << " Number of edges: Ainfo " << boost::num_edges(aigraph)
		<< " ResidueType " <<  boost::num_edges(rsd_type.graph()) << std::endl;

	core::Real const no_match_found_score( -999999 );
	core::Real best_score( no_match_found_score );
	GeometricRenameIsomorphismCallback callback( aigraph, rinfo, rsd_type, mapping, best_score );
	GeometricRenameVerticiesEquivalent vertices_equivalent( aigraph, rsd_type.graph() );

	boost::vf2_subgraph_mono( aigraph, rsd_type.graph(), callback, small_order,
		boost::vertices_equivalent( vertices_equivalent ) );

	if ( best_score == no_match_found_score ) {
		// The McGregor approach takes longer, but picks up additional mapping correspondences.
		// (In particular, it allows for ignoring extraneous atoms on the PDB side.)
		// It also tends to output a *lot* more mappings. In fact, we're only considering connected
		// subgraphs to speed things up. This means that we're going to match on just the largest
		// connected component, rather than multiple disconnected subgraphs.
		TR.Debug << "Using the McGregor fall-back approach." << std::endl;
		boost::mcgregor_common_subgraphs( aigraph, rsd_type.graph(), true /*only_connected_subgraphs*/,
			callback, boost::vertices_equivalent( vertices_equivalent ) );
	}

	if ( best_score == no_match_found_score ) {
		TR.Error << "ERROR: Difficulties mapping atom names from geometry for " << rinfo.resName() << " " << rinfo.chainID() << rinfo.resSeq()
			<< rinfo.iCode() << " onto " << rsd_type.name() << std::endl;
		utility_exit_with_message("Can't find good mapping between input residue and residue type.");
	}

	TR.Debug << "After geometric name remapping, missing " << (rsd_type.natoms()-mapping.left.size()) << " of " << rsd_type.natoms() << " atoms." << std::endl;
}


} // namespace pose_from_sfr
} // namespace io
} // namespace core
