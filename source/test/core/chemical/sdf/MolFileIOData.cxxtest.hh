// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/sdf/MolfileIOData.cxxtest.hh
/// @brief unit tests for the MolfileIOData file
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/chemical/sdf/MolFileIOData.hh>
#include <core/chemical/sdf/MolFileIOReader.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/Atom.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/MMAtomType.hh>


#include <core/chemical/icoor_support.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/util.hh>

#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/sdf/mol_writer.hh>

#include <numeric/model_quality/rms.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/FArray2.io.hh>

#include <boost/graph/mcgregor_common_subgraphs.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

// C++ Headers
#include <string>
#include <set>
#include <cmath>

static basic::Tracer TR("core.chemical.sdf.MolfileIOData.cxxtest");

using namespace core::chemical;
using namespace core::chemical::sdf;

typedef std::map< core::chemical::VD, core::chemical::VD > VDVDmap;

template< typename Graph1, typename Graph2 >
class IsomorphismCallback {
	typedef std::map< typename Graph1::vertex_descriptor, typename Graph2::vertex_descriptor > Mapping;
	IsomorphismCallback();

public:

	IsomorphismCallback(Graph1 const & rsd1, Graph2 const & rsd2, utility::vector1< Mapping > & vertex_maps ):
		rsd1_( rsd1 ),
		rsd2_( rsd2 ),
		mappings_( vertex_maps )
	{}

	template< class VD2VDmap1, class VD2VDmap2 >
	bool operator()(VD2VDmap1 map_1_to_2, VD2VDmap2 ) {
		Mapping newmap;
		//TR << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Found isomorphism" << std::endl;
		typename Graph1::vertex_iterator iter, iter_end;
		for ( boost::tie( iter, iter_end) = boost::vertices(rsd1_); iter != iter_end; ++iter ) {
			if ( map_1_to_2[ *iter ] != Graph2::null_vertex() ) {
				newmap[ *iter ] = map_1_to_2[ *iter ]; // Because of random type conflicts
				//TR << rsd1_[ *iter].name() << "  " << rsd2_[ map_1_to_2[ *iter ] ].name()  << std::endl;
			}
		}
		mappings_.push_back( newmap );
		//TR << "--------------------------" << std::endl;
		return true;
	}

public:
	Graph1 const & rsd1_;
	Graph2 const & rsd2_;
	utility::vector1< Mapping > & mappings_;
};

bool
atom_equivalent( Atom const & one, Atom const & other, bool exact, bool verbose ) {
	// Only compare names if exact
	if ( exact && one.name() != other.name() ) {
		if ( verbose ) TR << "Bad name. " << one.name() << " " << other.name() << std::endl;
		return false;
	}
	// Elements must be set and match to be considered equivalent.
	if ( ! one.element_type() || ! other.element_type() ) {
		if ( verbose ) TR << "No element. " << one.name() << " " << other.name() << std::endl;
		return false;
	}
	if ( one.element_type()->get_chemical_name() != other.element_type()->get_chemical_name() ) {
		if ( verbose ) TR << "Bad element. " << one.name() << " " << other.name() << std::endl;
		return false;
	}
	// If the type indicies are set, they must match. (For exact, they must always match)
	if ( ( exact || (one.atom_type_index() && other.atom_type_index()) )
			&& one.atom_type_index() != other.atom_type_index() ) {
		if ( verbose ) {
			TR << "Bad atom type. " << one.name() << " (" << one.atom_type_index() << ") "
				<< other.name() << " (" << other.atom_type_index() << ") " << std::endl;
		}
		return false;
	}
	// VIRT/X issue with testing runs -- only care about mm type if exact
	if ( ( exact ) && one.mm_atom_type_index() != other.mm_atom_type_index() )  {
		if ( verbose ) TR << "Bad mm atom type. " << one.name() << " " << other.name()  << std::endl;
		return false;
	}
	if ( one.gasteiger_atom_type() && other.gasteiger_atom_type() &&
			one.gasteiger_atom_type()->get_name() != other.gasteiger_atom_type()->get_name() ) {
		if ( verbose ) TR << "Bad gast atom type. " << one.name() << " " << other.name()  << std::endl;
		return false;
	} else if ( exact ) {
		if ( ( one.gasteiger_atom_type() && ! other.gasteiger_atom_type() ) ||
				( ! one.gasteiger_atom_type() && other.gasteiger_atom_type() ) ) {
			if ( verbose ) TR << "No gast atom type. " << one.name() << " " << other.name()  << std::endl;
			return false;
		}
	}
	if ( one.formal_charge() != other.formal_charge() ) {
		if ( verbose ) TR << "Bad formal charge. " << one.name() << " " << other.name() << ": " << one.formal_charge() << " vs  " << other.formal_charge() << std::endl;
		return false;
	}
	if ( (one.charge() < other.charge() - 0.01) || (one.charge() > other.charge() + 0.01) ) { // 0.01 being the typical resolution in the params files.
		if ( verbose ) TR << "Bad partial charge. " << one.name() << " " << other.name() << ": " << one.charge() << " vs " << other.charge() << std::endl;
		return false;
	}
	if ( exact && one.ideal_xyz() != other.ideal_xyz() ) {
		if ( verbose ) TR << "Bad coord. " << one.name() << " " << other.name()  << std::endl;
		return false;
	}

	if ( one.is_hydrogen() != other.is_hydrogen() ) {
		if ( verbose ) TR << "Bad is_hydrogen " << one.name() << " " << other.name()  << std::endl;
		return false;
	}
	if ( one.is_polar_hydrogen() != other.is_polar_hydrogen() ) {
		if ( verbose ) TR << "Bad is_polar_hydrogen " << one.name() << " " << other.name()  << std::endl;
		return false;
	}
	if ( one.is_haro() != other.is_haro() ) {
		if ( verbose ) TR << "Bad is_haro " << one.name() << " " << other.name()  << std::endl;
		return false;
	}
	if ( one.is_acceptor() != other.is_acceptor() ) {
		if ( verbose ) TR << "Bad is_acceptor " << one.name() << " " << other.name()  << std::endl;
		return false;
	}
	if ( one.is_virtual() != other.is_virtual() ) {
		if ( verbose ) TR << "Bad is_virtual " << one.name() << " " << other.name()  << std::endl;
		return false;
	}
	if ( one.heavyatom_has_polar_hydrogens() != other.heavyatom_has_polar_hydrogens() ) {
		if ( verbose ) TR << "Bad heavyatom_has_polar_hydrogens " << one.name() << " " << other.name()  << std::endl;
		return false;
	}
	if ( one.has_orbitals() != other.has_orbitals() ) {
		if ( verbose ) TR << "Bad has_orbitals " << one.name() << " " << other.name()  << std::endl;
		return false;
	}

	return true;
}

template< typename Graph1, typename Graph2 >
class VerticiesEquivalent {
public:

	VerticiesEquivalent( Graph1 const & rsd1, Graph2 const & rsd2 ) :
		rsd1_( rsd1 ),
		rsd2_( rsd2 )
	{}

	bool operator() ( typename Graph1::vertex_descriptor vd1, typename Graph2::vertex_descriptor vd2 ) {
		return rsd1_[vd1].element_type() == rsd2_[vd2].element_type(); // part of equivalency
		//return atom_equivalent( rsd1_[vd1], rsd2_[vd2], false, false );
	}

public:
	Graph1 const & rsd1_;
	Graph2 const & rsd2_;
};

template< typename Graph1, typename Graph2 >
class EdgesEquivalent {
public:

	EdgesEquivalent( Graph1 const & rsd1, Graph2 const & rsd2 ) :
		rsd1_( rsd1 ),
		rsd2_( rsd2 )
	{}

	bool operator() ( typename Graph1::edge_descriptor /*ed1*/, typename Graph2::edge_descriptor /*ed2*/ ) {
		return true;
		//core::chemical::Bond const & bond1( rsd1_[ed1] );
		//core::chemical::Bond const & bond2( rsd2_[ed2] );
		//return bond1.bond_name() != bond2.bond_name() && bond1.bond_name() != 4 && bond2.bond_name() != 4;
	}

public:
	Graph1 const & rsd1_;
	Graph2 const & rsd2_;
};
/// @details map is a mapping of VDs of rsd1 to rsd2
bool compare_residues_mapping( ResidueType const & rsd1, ResidueType const & rsd2, VDVDmap & map) {
	// First, compare the non-mapping dependant options.
	if ( rsd1.atom_type_set_ptr() != rsd2.atom_type_set_ptr() ) {
		TR << "Bad atom type set match." << std::endl;
		return false;
	}
	if ( rsd1.name() != rsd2.name() ) {
		TR << "Name mismatch" << std::endl;
		return false;
	}
	if ( rsd1.name3() != rsd2.name3() ) {
		TR << "Name3 mismatch" << std::endl;
		return false;
	}
	if ( rsd1.name1() != rsd2.name1() ) {
		TR << "Name1 mismatch" << std::endl;
		return false;
	}
	if ( rsd1.aa() != rsd2.aa() ) {
		TR << "AA mismatch" << std::endl;
		return false;
	}
	if ( rsd1.interchangeability_group() != rsd2.interchangeability_group() ) {
		TR << "Interchangeability group mismatch" << std::endl;
		return false;
	}
	if ( rsd1.is_polymer() != rsd2.is_polymer()  ||
			rsd1.is_protein() != rsd2.is_protein()  ||
			rsd1.is_alpha_aa() != rsd2.is_alpha_aa()  ||
			rsd1.is_beta_aa() != rsd2.is_beta_aa()  ||
			rsd1.is_d_aa() != rsd2.is_d_aa()  ||
			rsd1.is_l_aa() != rsd2.is_l_aa()  ||
			rsd1.is_DNA() != rsd2.is_DNA()  ||
			rsd1.is_RNA() != rsd2.is_RNA()  ||
			rsd1.is_NA() != rsd2.is_NA()  ||
			rsd1.is_coarse() != rsd2.is_coarse()  ||
			rsd1.is_ligand() != rsd2.is_ligand()  ||
			rsd1.is_metal() != rsd2.is_metal()  ||
			rsd1.is_metalbinding() != rsd2.is_metalbinding()  ||
			rsd1.is_surface() != rsd2.is_surface()  ||
			rsd1.is_polar() != rsd2.is_polar()  ||
			rsd1.is_charged() != rsd2.is_charged()  ||
			rsd1.is_aromatic() != rsd2.is_aromatic()  ||
			rsd1.is_cyclic() != rsd2.is_cyclic()  ||
			rsd1.is_terminus() != rsd2.is_terminus()  ||
			rsd1.is_lower_terminus() != rsd2.is_lower_terminus()  ||
			rsd1.is_upper_terminus() != rsd2.is_upper_terminus()  ||
			rsd1.is_branch_lower_terminus() != rsd2.is_branch_lower_terminus()  ||
			rsd1.is_acetylated_nterminus() != rsd2.is_acetylated_nterminus()  ||
			rsd1.is_virtual_residue() != rsd2.is_virtual_residue()  ||
			rsd1.is_adduct() != rsd2.is_adduct() ) {
		TR << "Properties mismatch!" << std::endl;
		TR << rsd1.name() << ": ";
		utility::vector1< std::string > properties( rsd1.properties().get_list_of_properties() );
		for ( core::Size p1( 1 ); p1 <= properties.size(); ++p1 ) {
			TR << properties[p1] << " ";
		}
		TR << std::endl;
		TR << rsd2.name() << ": ";
		properties = rsd2.properties().get_list_of_properties();
		for ( core::Size p2( 1 ); p2 <= properties.size(); ++p2 ) {
			TR << properties[p2] << " ";
		}
		TR << std::endl;
		return false;
	}
	if ( ! variants_match( rsd1, rsd2 ) ) {
		TR << "Variants Mismatch" << std::endl;
		return false;
	}

	// For testing purposes we want to have a strict comparison.
	if ( rsd1.natoms() != rsd2.natoms() || map.size() != rsd1.natoms() ) {
		TR << "Wrong size!" << std::endl;
		return false;
	}

	if ( rsd1.nheavyatoms() != rsd2.nheavyatoms() ) {
		TR << "Wrong nheavyatoms!" << rsd1.nheavyatoms() << " " << rsd2.nheavyatoms() << std::endl;
		return false;
	}
	if ( rsd1.n_hbond_acceptors() != rsd2.n_hbond_acceptors() ) {
		TR << "Wrong n_hbond_acceptors!" << rsd1.n_hbond_acceptors() << " " << rsd2.n_hbond_acceptors() << std::endl;
		return false;
	}
	if ( rsd1.n_hbond_donors() != rsd2.n_hbond_donors() ) {
		TR << "Wrong n_hbond_donors!" << rsd1.n_hbond_donors() << " " << rsd2.n_hbond_donors() << std::endl;
		return false;
	}
	if ( rsd1.nbonds() != rsd2.nbonds() ) {
		TR << "Wrong nbonds!" << rsd1.nbonds() << " " << rsd2.nbonds() << std::endl;
		return false;
	}
	if ( rsd1.last_backbone_atom() != rsd2.last_backbone_atom() ) {
		TR << "Wrong last_backbone_atom!" << rsd1.last_backbone_atom() << " " << rsd2.last_backbone_atom() << std::endl;
		return false;
	}
	if ( rsd1.first_sidechain_atom() != rsd2.first_sidechain_atom() ) {
		TR << "Wrong first_sidechain_atom!" << rsd1.first_sidechain_atom() << " " << rsd2.first_sidechain_atom() << std::endl;
		return false;
	}
	if ( rsd1.first_sidechain_hydrogen() != rsd2.first_sidechain_hydrogen() ) {
		TR << "Wrong first_sidechain_hydrogen!" << rsd1.first_sidechain_hydrogen() << " " << rsd2.first_sidechain_hydrogen() << std::endl;
		return false;
	}
	if ( rsd1.n_orbitals() != rsd2.n_orbitals() ) {
		TR << "Wrong n_orbitals!" << rsd1.n_orbitals() << " " << rsd2.n_orbitals() << std::endl;
		return false;
	}

	//compare the Atoms
	for ( VDVDmap::const_iterator iter(map.begin()); iter != map.end(); ++iter ) {
		if ( ! atom_equivalent( rsd1.atom( iter->first ), rsd2.atom( iter->second ), false, true ) ) {
			TR << "Bad Atoms" << std::endl;
			return false;
		}
	}

	//compare the Bonds
	core::chemical::EIter eiter, eiter_end;
	for ( boost::tie( eiter, eiter_end ) = boost::edges( rsd1.graph() ); eiter != eiter_end; ++eiter ) {
		core::chemical::VD source1( boost::source( *eiter, rsd1.graph()) ), target1( boost::target( *eiter, rsd1.graph()) );
		if ( map.count( source1 ) == 0 || map.count( target1 ) == 0 ) {
			TR << "Missing in map " << rsd1.atom( source1 ).name() << " or " << rsd1.atom( target1 ).name() << std::endl;
			return false;
		}
		core::chemical::ED edge2;
		bool found(true);
		boost::tie( edge2, found ) = boost::edge( map[ source1 ], map[ target1 ], rsd2.graph() );
		if ( ! found ) {
			TR << "Edge not present in both !!!! " << rsd1.atom( source1 ).name() << "  " << rsd1.atom( target1 ).name() << std::endl;
			return false;
		}
		core::chemical::Bond const & bond1( rsd1.graph()[ *eiter ] );
		core::chemical::Bond const & bond2( rsd2.graph()[ edge2 ] );
		if ( bond1.bond_name() != bond2.bond_name() && bond1.bond_name() != 4 && bond2.bond_name() != 4 ) { // Ignore aro-single/double mismatch
			TR << "Not same bond type!!!! " << rsd1.atom( source1 ).name() << "  " << rsd1.atom( target1 ).name() << bond1.bond_name() << " " << bond2.bond_name() << std::endl;
			return false;
		}
		// The ring can be split in different locations
		//  if( bond1.cut_bond() != bond2.cut_bond() ) {
		//   TR << "Not same cut bond type!!!! " << rsd1.atom( source1 ).name() << "  " << rsd1.atom( target1 ).name() << std::endl;
		//   return false;
		//  }
	}

	//Nbr atom:
	core::chemical::VD nbr1( rsd1.atom_vertex( rsd1.nbr_atom() ) );
	core::chemical::VD nbr2( rsd2.atom_vertex( rsd2.nbr_atom() ) );
	if ( map[ nbr1 ] != nbr2 ) {
		TR << "Nbr atom mismatch:" << rsd1.atom_name( nbr1 ) << " -> should be " << rsd2.atom_name( map[ nbr1 ] ) << " is " << rsd2.atom_name( nbr2 ) << std::endl;
		VD nbr2_rev( ResidueGraph::null_vertex() );
		for ( VDVDmap::const_iterator itr( map.begin() ); itr != map.end(); ++itr ) {
			if ( itr->second == nbr2 ) {
				nbr2_rev = itr->first;
			}
		}
		std::string maxatom, maxatom1, maxatom2;
		core::Real maxdist( 0 ), maxdist1( 0 ), maxdist2( 0 );
		core::Vector nbrxyz( rsd1.atom( nbr1 ).ideal_xyz() );
		core::Vector nbr2xyz( rsd1.atom( nbr2_rev ).ideal_xyz() );
		for ( core::Size ii(1); ii <= rsd1.natoms(); ++ii ) {
			if ( rsd1.atom_is_hydrogen( ii ) ) continue;
			core::Real dist( nbrxyz.distance( rsd1.atom(ii).ideal_xyz() ) );
			if ( dist > maxdist ) {
				maxdist = dist;
				maxatom = rsd1.atom_name(ii);
			} else if ( dist == maxdist ) {
				maxatom = maxatom + " " + rsd1.atom_name(ii);
			}
			core::Real dist2( nbr2xyz.distance( rsd1.atom(ii).ideal_xyz() ) );
			if ( dist2 > maxdist2 ) {
				maxdist2 = dist2;
				maxatom2 = rsd1.atom_name(ii);
			} else if ( dist == maxdist2 ) {
				maxatom2 = maxatom2 + ", " + rsd1.atom_name(ii);
			}
		}
		TR << "Residue 1: nbr atom " << rsd1.atom_name( nbr1 ) << " distance to " << maxatom << " is " << maxdist << std::endl;
		TR << "Residue 1: residue 2 eq neighbor " << rsd1.atom_name( nbr2_rev ) << " distance to " << maxatom2 << " is " << maxdist2 << std::endl;
		maxatom = "";
		maxdist = 0;
		nbrxyz = rsd2.atom( rsd2.nbr_atom() ).ideal_xyz();
		core::Vector nbr1xyz( rsd2.atom( map[ nbr1 ] ).ideal_xyz() );
		for ( core::Size ii(1); ii <= rsd2.natoms(); ++ii ) {
			if ( rsd2.atom_is_hydrogen( ii ) ) continue;
			core::Real dist( nbrxyz.distance( rsd2.atom(ii).ideal_xyz() ) );
			if ( dist > maxdist ) {
				maxdist = dist;
				maxatom = rsd2.atom_name(ii);
			} else if ( dist == maxdist ) {
				maxatom = maxatom + ", " + rsd2.atom_name(ii);
			}
			core::Real dist1( nbr1xyz.distance( rsd2.atom(ii).ideal_xyz() ) );
			if ( dist1 > maxdist1 ) {
				maxdist1 = dist1;
				maxatom1 = rsd2.atom_name(ii);
			} else if ( dist == maxdist1 ) {
				maxatom1 = maxatom1 + ", " + rsd2.atom_name(ii);
			}
		}
		TR << "Residue 2: residue 1 eq neighbor " << rsd2.atom_name( map[ nbr1 ] ) << " distance to " << maxatom1 << " is " << maxdist1 << std::endl;
		TR << "Residue 2: nbr atom " << rsd2.atom_name( rsd2.nbr_atom() ) << " distance to " << maxatom << " is " << maxdist << std::endl;
		return false;

		//Nbr radius:
		if ( rsd1.nbr_radius() != rsd2.nbr_radius() ) {
			TR << "Nbr radius mismatch." << std::endl;
			return false;
		}
	} // if map[ nbr1 ] != nbr2

	if ( rsd1.nchi() != rsd2.nchi() ) {
		TR << "Wrong nchi! " << rsd1.nchi() << " " << rsd2.nchi() << std::endl;
		TR << "Residue 1:" << std::endl;
		print_chis( TR, rsd1);
		TR << "Residue 2:" << std::endl;
		print_chis( TR, rsd2);
		return false;
	}
	if ( rsd1.n_proton_chi() != rsd2.n_proton_chi() ) {
		TR << "Wrong n_proton_chi!" << rsd1.n_proton_chi() << " " << rsd2.n_proton_chi() << std::endl;
		return false;
	}
	// Chis -- a chi representing each rotatable bond must be present, and the proton states must match.
	for ( core::Size chi1(1); chi1 <= rsd1.nchi(); ++chi1 ) {
		AtomIndices const & chiatoms1( rsd1.chi_atoms( chi1 ) );
		// The center two atoms in the chi -- different outside references may exist.
		VD mapped2a( map[ rsd1.atom_vertex( chiatoms1[2] ) ] ), mapped2b( map[ rsd1.atom_vertex( chiatoms1[3] ) ] );
		bool found(false), isproton(false);
		for ( core::Size chi2(1); chi2 <= rsd2.nchi(); ++chi2 ) {
			AtomIndices const & chiatoms2( rsd2.chi_atoms( chi2 ) );
			VD vd2a( rsd2.atom_vertex( chiatoms2[2] ) ), vd2b( rsd2.atom_vertex(chiatoms2[3]) );
			// The orientation of the bond is allowed to flip
			if ( ( mapped2a == vd2a && mapped2b == vd2b ) || ( mapped2a == vd2b && mapped2b == vd2a ) ) {
				found=true;
				isproton= rsd2.is_proton_chi( chi2 );
				break;
			}
		}
		if ( ! found ) {
			TR << "Couldn't find chi matching: " << rsd1.atom_name( chiatoms1[2] ) << " -- " << rsd1.atom_name( chiatoms1[3] ) << std::endl;
			for ( core::Size chi2(1); chi2 <= rsd2.nchi(); ++chi2 ) {
				AtomIndices const & chiatoms2( rsd2.chi_atoms( chi2 ) );
				TR << rsd2.atom_name( chiatoms2[2] ) << " -- " << rsd2.atom_name( chiatoms2[3] ) << std::endl;
			}
			return false;

		}
		if ( isproton != rsd1.is_proton_chi( chi1 ) ) {
			TR << "Proton chi states don't match: " << rsd1.atom_name( chiatoms1[2] ) << " -- " << rsd1.atom_name( chiatoms1[3] ) << std::endl;
			return false;
		}
	}

	// Icoords -- here we just check that the internal coordinates match up, even after bouncing through the icoor->xyz calculator.

	// First, straight.
	{ // scoping
		ObjexxFCL::FArray2D< core::Real > p1a( 3, map.size() );
		ObjexxFCL::FArray2D< core::Real > p2a( 3, map.size() );
		core::Size natoms(0);
		for ( VDVDmap::const_iterator iter( map.begin() ); iter != map.end(); ++iter ) {
			++natoms;
			core::Vector const & coord1( rsd1.graph()[ iter->first ].ideal_xyz() );
			core::Vector const & coord2( rsd2.graph()[ iter->second ].ideal_xyz() );
			for ( core::Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
				p1a(k,natoms) = coord1(k); // Call is 1-based indexing of xyzVector
				p2a(k,natoms) = coord2(k);
			}
		}
		core::Real rmsd( numeric::model_quality::rms_wrapper( map.size(), p1a, p2a ) );

		if ( rmsd > 0.01 ) { // Somewhat arbitrary choice
			TR << "Ideal coordinates don't match. Rmsd: " << rmsd << std::endl;
			//TR << "Coords 1:" << std::endl;
			//TR << p1a << std::endl;
			//TR << "Coords 2:" << std::endl;
			//TR << p2a << std::endl;
			//core::conformation::ResidueCOP residue1( new core::conformation::Residue( rsd1, true ) );
			//core::conformation::ResidueCOP residue2( new core::conformation::Residue( rsd2, true ) );
			//TR << "Residue 1 icoor: " << formatted_icoord_tree( rsd1 ) << std::endl;
			//TR << "Residue 2 icoor: " << formatted_icoord_tree( rsd2 ) << std::endl;
			//print_chis( TR, rsd1);
			//core::chemical::sdf::MolWriter write;
			//write.output_residue( rsd1.name() + "_rsd1.sdf", residue1 );
			//write.output_residue( rsd2.name() + "_rsd2.sdf", residue2 );
			//core::pose::Pose test_pose;
			//test_pose.append_residue_by_jump( *residue1, 1 );
			//test_pose.dump_pdb( rsd1.name() + "_rsd1.pdb" );
			//test_pose.replace_residue( 1, *residue2, false );
			//test_pose.dump_pdb( rsd2.name() + "_rsd2.pdb" );
			//test_pose.replace_residue( 1, *residue1, true );
			//test_pose.dump_pdb( rsd1.name() + "_rsd1_orient.pdb" );
			return false;
		}
	}

	// Now we call icoor->xyz (on both) and see if they still match
	ResidueType copy1(rsd1), copy2(rsd2); // Make copy for alteration.

	//TR << "------------------- rsd1 ------------------------" <<std::endl;
	//rsd1.dump_vd_info();
	//TR << " ~~~~ " << std::endl;
	//pretty_print_atomicoor(TR, rsd1.icoor( 1 /*rsd1.nbr_atom()*/ ), rsd1);
	//TR << "------------------- rsd2 ------------------------" <<std::endl;
	//rsd2.dump_vd_info();
	//TR << " ~~~~ " << std::endl;
	//pretty_print_atomicoor(TR, rsd2.icoor( 1 /*rsd2.nbr_atom()*/ ), rsd2);
	//TR << "------------------- copy1 ------------------------" <<std::endl;
	//copy1.dump_vd_info();
	//TR << " ~~~~ " << std::endl;
	//pretty_print_atomicoor(TR, copy1.icoor( 1 /*copy1.nbr_atom()*/ ), copy1);
	//TR << "------------------- copy2 ------------------------" <<std::endl;
	//copy2.dump_vd_info();
	//TR << " ~~~~ " << std::endl;
	//pretty_print_atomicoor(TR, copy2.icoor( 1 /*copy2.nbr_atom()*/ ), copy2);
	//TR << "--------------------------------------------------" <<std::endl;

	fill_ideal_xyz_from_icoor( copy1, copy1.graph() );
	fill_ideal_xyz_from_icoor( copy2, copy2.graph() );

	{ // scoping
		ObjexxFCL::FArray2D< core::Real > p1a( 3, map.size() );
		ObjexxFCL::FArray2D< core::Real > p2a( 3, map.size() );
		core::Size natoms(0);
		for ( VDVDmap::const_iterator iter( map.begin() ); iter != map.end(); ++iter ) {
			++natoms;
			// Need to go through names as the vds for the copies have changed.
			std::string const & name1( rsd1.atom_name( iter->first ) );
			std::string const & name2( rsd2.atom_name( iter->second ) );
			VD cpvd1( copy1.atom_vertex( copy1.atom_index( name1 ) ) ), cpvd2( copy2.atom_vertex( copy2.atom_index( name2 ) ) );
			core::Vector const & coord1( copy1.graph()[ cpvd1 ].ideal_xyz() );
			core::Vector const & coord2( copy2.graph()[ cpvd2 ].ideal_xyz() );
			for ( core::Size k = 1; k <= 3; ++k ) { // k = X, Y and Z
				p1a(k,natoms) = coord1(k); // Call is 1-based indexing of xyzVector
				p2a(k,natoms) = coord2(k);
			}
		}
		core::Real rmsd( numeric::model_quality::rms_wrapper( map.size(), p1a, p2a ) );

		if ( rmsd > 0.01 ) { // Somewhat arbitrary choice
			TR << "Ideal coordinates don't match after icoor->xyz. rmsd: " << rmsd << std::endl;
			//core::chemical::sdf::MolWriter write;
			//core::conformation::ResidueCOP residue1( new core::conformation::Residue( copy1, true ) );
			//core::conformation::ResidueCOP residue2( new core::conformation::Residue( copy2, true ) );
			//write.output_residue( copy1.name() + "_copy1.sdf", residue1 );
			//write.output_residue( copy2.name() + "_copy2.sdf", residue2 );
			return false;
		}
	}

	// If we've gotten here, the residues are pretty similar
	return true;
}

/// @brief utility function for seeing if two residue types are "equivalent"
bool match_restype( ResidueType const & rsd1, ResidueType const & rsd2 ) {
	if ( rsd1.natoms() != rsd2.natoms() ) {
		return false;
	} // Otherwise we need to assure that rsd1 has fewer atoms than rsd2 for the isomorphism functions.

	// If we were doing heavy atoms only ...
	//core::chemical::HeavyAtomGraph const ha1( rsd1.heavy_atoms() );
	//core::chemical::ResidueGraph const & ha_only( rsd1.graph() );
	//core::chemical::CopyVertex<core::chemical::HeavyAtomGraph,core::chemical::ResidueGraph> copy_atom( ha1, ha_only );
	//core::chemical::CopyEdge<core::chemical::HeavyAtomGraph,core::chemical::ResidueGraph> copy_edge( ha1, ha_only );
	//boost::copy_graph( ha1, ha_only, boost::vertex_copy( copy_atom ).edge_copy( copy_edge ) );

	utility::vector1< core::chemical::VD > small_order;
	core::chemical::VIter iter, iter_end;
	for ( boost::tie( iter, iter_end) = boost::vertices(rsd1.graph()); iter != iter_end; ++iter ) {
		small_order.push_back( *iter );
	}

	utility::vector1< VDVDmap > mappings;
	IsomorphismCallback<core::chemical::ResidueGraph, core::chemical::ResidueGraph> callback( rsd1.graph(), rsd2.graph(), mappings );
	VerticiesEquivalent<core::chemical::ResidueGraph, core::chemical::ResidueGraph> vertices_equivalent( rsd1.graph(), rsd2.graph() );
	EdgesEquivalent<core::chemical::ResidueGraph, core::chemical::ResidueGraph> edges_equivalent( rsd1.graph(), rsd2.graph() );
	//TR << "$$$ FINDING SUBGRAPHS" << std::endl;
	boost::vf2_subgraph_mono( rsd1.graph(), rsd2.graph(), callback, small_order,  boost::vertices_equivalent( vertices_equivalent ).edges_equivalent( edges_equivalent) );
	//TR << "$$$ DONE FINDING SUBGRAPHS" << std::endl;
	TR << "Found " << mappings.size() << " subgraph mappings between " << rsd1.name() << " and " << rsd2.name() << std::endl;
	bool one_matches( false );
	for ( core::Size n(1); n <= mappings.size(); ++n ) {
		VDVDmap & map( mappings[n] );
		for ( VDVDmap::const_iterator itr( map.begin() ); itr != map.end(); ++itr ) {
			TR.Debug << rsd1.atom_name( itr->first ) << ":" << rsd2.atom_name( itr->second ) << "  ";
		}
		TR.Debug << std::endl;
		if ( compare_residues_mapping( rsd1, rsd2, map ) ) {
			one_matches=true;
			TR << "Mapping " << n << " Passes" << std::endl;
			return true;
		} else {
			TR << "Mapping " << n << " Failed" << std::endl;
		}
	}
	return one_matches;
}


class MolFileIODataTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_sdfreader() {
		using namespace core::chemical;
		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		ElementSetCOP element_types = cm->element_set("default");
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
		orbitals::OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);
		ResidueTypeSetOP rsd_types( new ResidueTypeSet );

		sdf::MolFileIOReader molfile_reader;

		//This file should contain names for paired sdf/params files.
		utility::io::izstream paramslist("core/chemical/sdf/paired_list.txt");
		std::string molfile, paramsfile;
		paramslist >> molfile >> paramsfile;
		while ( paramslist ) {
			if ( molfile[0] != '#' ) {
				TR << "------- Comparing  " << molfile << " and " << paramsfile << std::endl;

				// Read reference
				core::chemical::ResidueTypeOP rsd_ref = read_topology_file(paramsfile,
					atom_types, element_types, mm_atom_types, orbital_types, ResidueTypeSetCAP(rsd_types));

				// Read molfile (reader has sensible defaults for typesets in use)
				utility::vector1< sdf::MolFileIOMoleculeOP > data( molfile_reader.parse_file( molfile ) );
				utility::vector1< ResidueTypeOP > rtvec( sdf::convert_to_ResidueTypes( data, false ) );
				TS_ASSERT( rtvec.size() == 1 ); // These should all have a single entry.
				if ( rtvec.size() > 0 ) {
					bool restypes_match( match_restype( *rtvec[1], *rsd_ref ) );
					if ( ! restypes_match ) {
						core::chemical::write_topology_file( *rtvec[1] );
						core::chemical::ResidueType const & rsd1( *rtvec[1] );
						core::chemical::ResidueType const & rsd2( *rsd_ref );
						core::chemical::sdf::MolWriter write;
						core::conformation::ResidueCOP residue1( core::conformation::ResidueOP( new core::conformation::Residue( rsd1, true ) ) );
						core::conformation::ResidueCOP residue2( core::conformation::ResidueOP( new core::conformation::Residue( rsd2, true ) ) );
						write.output_residue( rsd1.name() + "_rsd1.sdf", residue1 );
						write.output_residue( rsd2.name() + "_rsd2.sdf", residue2 );
						core::pose::Pose test_pose;
						test_pose.append_residue_by_jump( *residue1, 1 );
						test_pose.dump_pdb( rsd1.name() + "_rsd1.pdb" );
						test_pose.replace_residue( 1, *residue2, false );
						test_pose.dump_pdb( rsd2.name() + "_rsd2.pdb" );
						test_pose.replace_residue( 1, *residue1, true );
						test_pose.dump_pdb( rsd1.name() + "_rsd1_orient.pdb" );
					} else {
						TR << ">>>>>>" << rtvec[1]->name() << " PASSES <<<<<<<<" << std::endl;
					}
					TS_ASSERT( restypes_match );
				}
			}
			paramslist >> molfile >> paramsfile;
		} // while( paramslist )
		paramslist.close();
	} //  test_sdf_reader

	void test_extra_data() {
		sdf::MolFileIOReader molfile_reader;
		utility::vector1< sdf::MolFileIOMoleculeOP > data( molfile_reader.parse_file( "core/chemical/sdf/example.sdf" ) );
		TS_ASSERT_EQUALS( data.size(), 1 );
		if ( data.size() < 1 ) return;
		sdf::MolFileIOMoleculeOP molecule( data[1] );

		sdf::StrStrMap const & extras( molecule->get_str_str_data() );
		TS_ASSERT_EQUALS( extras.size(), 8 );
		TS_ASSERT( extras.count( "Atom_EffectivePolarizability" ) );
		TS_ASSERT( extras.count( "Atom_LonePairEN" ) );
		TS_ASSERT_EQUALS( extras.find("Atom_LonePairEN")->second, "0 0 0.860678 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 " );
		TS_ASSERT( extras.count( "Atom PiCharge" ) );
		TS_ASSERT( extras.count( "Atom PiEN" ) );
		TS_ASSERT( extras.count( "Atom Polarizability" ) );
		TS_ASSERT( extras.count( "DummyActivity" ) );
		TS_ASSERT_EQUALS( extras.find("DummyActivity")->second, "0" );
		TS_ASSERT( extras.count( "EMOL_PARENT_ID" ) );
		TS_ASSERT_EQUALS( extras.find("EMOL_PARENT_ID")->second, "569781" );
		TS_ASSERT( extras.count( "MatchedFragmentIndices" ) );
	}
};
