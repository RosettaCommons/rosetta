// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/rdkit/ConformationSelectionToRDMol.cc
/// @brief  This class takes a part of a conformation from Rosetta and converts it into a RDKit RWMol object
/// @author Andy Watkins (watkina6@gene.com)

#include <protocols/drug_design/ConformationSelectionToRDMol.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/rdkit/util.hh>
// #include <core/chemical/MutableResidueType.hh>
// #include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/Bond.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueConnection.hh>

#include <core/id/AtomID.hh>

#include <utility/numbers.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>

#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SanitException.h>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h> // For MolToSmiles

#include <core/chemical/Element.hh> // AUTO IWYU For Element



namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.ConformationSelectionToRDMol");

ConformationSelectionToRDMol::ConformationSelectionToRDMol(
	core::conformation::Conformation const & conf,
	utility::vector1< core::Size > const & res,
	bool neutralize /* = true*/,
	bool keep_hydro /*= false*/
):
	conf_( conf ),
	res_( res ),
	neutralize_( neutralize ),
	keep_hydro_( keep_hydro )
{
	if ( neutralize_ && keep_hydro_ ) {
		utility_exit_with_message("In ConformationSelectionToRDMol, enabling neutralization while keeping the hydrogens doesn't make sense.");
	}
}


::RDKit::RWMOL_SPTR
ConformationSelectionToRDMol::Mol() {

	/// @brief Mapping of AtomID to indices of the rdkit object
	// Uses utility::get_undefined_size() for non-represented RDKit indices
	std::map< core::id::AtomID, core::Size > aid_to_index;

	::RDKit::RWMOL_SPTR rdmol( new ::RDKit::RWMol );
	rdmol->setProp("_Name", "PEPTIDE");
	rdmol->setProp("_3DConf",true);
	//rdmol->setProp("_MolFileChiralFlag",true);

	// Atoms
	std::map< core::Size, ::RDGeom::Point3D > pos_map;

	// core::chemical::VIter vertex_start, vertex_end; //iterators to start traversing the graph
	// boost::tie(vertex_start, vertex_end) = boost::vertices(res_.graph());
	for ( auto const & resi : res_ ) {
		for ( core::Size iiatom = 1; iiatom <= conf_.residue_type( resi ).natoms(); ++iiatom ) {
			//RDKit handling of virtual atoms is a little iffy at this point.
			if ( conf_.residue_type( resi ).is_virtual( iiatom ) ) {
				continue;
			}

			core::Size atomic_number = conf_.residue_type( resi ).element_type( iiatom )->get_atomic_number();
			debug_assert( atomic_number > 0 );

			// We keep hydrogens on the graph for now, for Kekulization reasons
			// -- with the explicit hydrogen settings, we can probably avoid that, but for now ...

			numeric::xyzVector<core::Real> coords = conf_.residue( resi ).xyz( iiatom );
			int formal_charge = conf_.residue_type( resi ).formal_charge( iiatom );

			// Raw pointer here is intentional - the RWMol will take exclusive ownership of the Atom
			::RDKit::Atom *rd_atom( new ::RDKit::Atom( (unsigned int) atomic_number ) );
			if ( formal_charge != 0 ) { // We take care of neutralizing later on
				rd_atom->setFormalCharge(formal_charge);
			}

			// Set atom name and index information (how well this is preserved depends on protocol).
			// Raw pointer here is intentional - the Atom will take exclusive ownership of the AtomPDBResidueInfo
			::RDKit::AtomPDBResidueInfo *info( new ::RDKit::AtomPDBResidueInfo(
				conf_.residue_type( resi ).atom_name( iiatom ), iiatom ) );
			rd_atom->setMonomerInfo(info); // takes ownership

			unsigned int aid = rdmol->addAtom(rd_atom,true,true); //Will take ownership of atom

			::RDGeom::Point3D pos( coords[0], coords[1], coords[2] );
			pos_map[ aid ] = pos;

			aid_to_index[core::id::AtomID( iiatom, resi ) ] = core::Size(aid);
		}
	}

	// Raw pointer here is intentional - the RWMol will take exclusive ownership of the conformer
	::RDKit::Conformer *conf = new ::RDKit::Conformer( rdmol->getNumAtoms() );
	conf->setId(0);
	conf->set3D(true);

	for ( auto const & elem : pos_map ) {
		conf->setAtomPos(elem.first, elem.second);
	}
	rdmol->addConformer(conf, true); // Takes ownership

	// Bonds, intra residue
	std::set< std::pair< unsigned int, unsigned int > > added;
	for ( auto const & resi : res_ ) {
		for ( core::Size iiatom = 1; iiatom <= conf_.residue_type( resi ).natoms(); ++iiatom ) {
			auto const & nbrs = conf_.residue_type( resi ).bonded_neighbor( iiatom );
			auto const & bond_types = conf_.residue_type( resi ).bonded_neighbor_types( iiatom );
			auto const & bond_ringnesses = conf_.residue_type( resi ).bonded_neighbor_ringnesses( iiatom );

			for ( core::Size jj = 1; jj <= nbrs.size(); ++jj ) {
				auto const nbr = nbrs[ jj ];
				auto const & bt = bond_types[ jj ];
				auto const & br = bond_ringnesses[ jj ];
				//Do not look at bonds between virtual atoms or virtual bonds
				if ( conf_.residue_type( resi ).is_virtual( iiatom )
						|| conf_.residue_type( resi ).is_virtual( nbr ) ) {
					// || bond.order() == core::chemical::PseudoBondOrder ) {
					continue;
				}

				auto idx1_it = aid_to_index.find( core::id::AtomID( iiatom, resi ) ); //core::Size index1 = aid_to_index_[];
				auto idx2_it = aid_to_index.find( core::id::AtomID( nbr, resi) );

				// Do not add bonds if the atoms aren't present in the fragment
				if ( idx1_it == aid_to_index.end() || idx2_it == aid_to_index.end() ) {
					TR.Trace << "Skipping bond between residue " << resi << " atoms " << iiatom << " and " << nbr << " due to skipped atoms." << std::endl;
					continue;
				}

				// Do not double count
				if ( added.find( std::make_pair( idx1_it->second, idx2_it->second ) ) != added.end() ) {
					continue;
				}
				if ( added.find( std::make_pair( idx2_it->second, idx1_it->second ) ) != added.end() ) {
					continue;
				}

				added.insert( std::make_pair( idx1_it->second, idx2_it->second ) );

				// Aromaticity only possible in aromatic side chains. iiatom and nbr must
				// both be side chain atoms (C-terminus of an aromatic side chain, also not
				// aromatic.) Note that this will fail for RARE NCAA residue types,
				// like 4-carboxyl-Phe.
				// AMW: supplemental check via iteration through bonded_neighbor_ringnesses
				// and looking for BondInRing

				::RDKit::Bond::BondType type( core::chemical::rdkit::convert_to_rdkit_bondtype( bt ) );

				if ( type == ::RDKit::Bond::AROMATIC ) {
					if ( iiatom >= conf_.residue_type( resi ).first_sidechain_atom()
							&& nbr >= conf_.residue_type( resi ).first_sidechain_atom()
							&& br == core::chemical::BondInRing // only endocyclic bonds!
							&& conf_.residue_type( resi ).has_property( core::chemical::AROMATIC ) ) {

						rdmol->getAtomWithIdx(idx1_it->second)->setIsAromatic(true);
						rdmol->getAtomWithIdx(idx2_it->second)->setIsAromatic(true);
						// TR << "Adding bond: " << (unsigned int) idx1_it->second << " " << (unsigned int) idx2_it->second << std::endl;
						rdmol->addBond( (unsigned int) idx1_it->second, (unsigned int) idx2_it->second, type);
					} else {
						// definitely not aromatic.
						rdmol->addBond( (unsigned int) idx1_it->second, (unsigned int) idx2_it->second, ::RDKit::Bond::ONEANDAHALF  );
					}
				} else {
					rdmol->addBond( (unsigned int) idx1_it->second, (unsigned int) idx2_it->second, type);
				}
			}
		}
		// bonded_neighbor_
	}

	// Now INTER-residue connections.
	for ( auto const & resi : res_ ) {
		core::conformation::Residue const & ii_res( conf_.residue( resi ) );
		core::Size const ii_nresconn = ii_res.type().n_possible_residue_connections();
		for ( core::Size jj = 1; jj <= ii_nresconn; ++jj ) {
			core::Size const jj_atid = ii_res.residue_connection( jj ).atomno();
			core::Size const jj_conn_res( ii_res.residue_connection_partner( jj ));
			core::Size const jj_conn_id(  ii_res.residue_connection_conn_id( jj ));
			if ( ! ( jj_atid && jj_conn_res && jj_conn_id ) ) continue;
			core::Size const jj_conn_atom( conf_.residue( jj_conn_res ).residue_connection( jj_conn_id ).atomno() );

			// OK, we have a single bond from an atom in our selection
			if ( !res_.has_value( jj_conn_res ) ) continue;

			if ( conf_.residue_type( resi ).is_virtual(jj_atid)
					|| conf_.residue_type( jj_conn_res ).is_virtual(jj_conn_atom) ) {
				continue;
			}

			auto idx1_it = aid_to_index.find( core::id::AtomID( jj_atid, resi ) ); //core::Size index1 = aid_to_index_[];
			auto idx2_it = aid_to_index.find( core::id::AtomID( jj_conn_atom, jj_conn_res) );

			// Do not add bonds if the atoms aren't present in the fragment
			if ( idx1_it == aid_to_index.end() || idx2_it == aid_to_index.end() ) {
				TR.Trace << "Skipping bond between residue " << resi << " atom " << jj_atid << " and residue " << jj_conn_res << " atom " << jj_conn_atom << " due to skipped atoms." << std::endl;
				continue;
			}

			// Do not double count
			if ( added.find( std::make_pair( idx1_it->second, idx2_it->second ) ) != added.end() ) {
				continue;
			}
			// Do not double count
			if ( added.find( std::make_pair( idx2_it->second, idx1_it->second ) ) != added.end() ) {
				continue;
			}

			added.insert( std::make_pair( idx1_it->second, idx2_it->second ) );
			// For now -- since bond types do not live on in Conformation -- we treat
			// all interresidue bonds as single. They certainly are for nucleic acids
			// and proteins, and that's good enough for now.
			::RDKit::Bond::BondType type( core::chemical::rdkit::convert_to_rdkit_bondtype( core::chemical::SingleBond ) );

			if ( type == ::RDKit::Bond::AROMATIC ) {
				rdmol->getAtomWithIdx(idx1_it->second)->setIsAromatic(true);
				rdmol->getAtomWithIdx(idx2_it->second)->setIsAromatic(true);
			}

			// TR << "Adding bond: " << (unsigned int) idx1_it->second << " " << (unsigned int) idx2_it->second << std::endl;
			rdmol->addBond( (unsigned int) idx1_it->second, (unsigned int) idx2_it->second, type);
		}
	}



	// Now attempt to fill out the other parameters for this molecule

	// (Nothing Here Yet)

	// Sanitize the molecule -- the one from Rosetta should be a decent molecule - probably.
	try {
		if ( neutralize_ ) {
			core::chemical::rdkit::neutralize_rdmol( *rdmol, /*addHs*/ false ); // Attempt to clean up silly charges.
			// Neutralization has an internal sanitization.
		} else {
			::RDKit::MolOps::sanitizeMol(*rdmol);
		}
	} catch (::RDKit::MolSanitizeException &se){
		TR.Error << "Cannot Sanitize molecule with RDKit: " << se.message() << std::endl;
		TR.Error << "    molecule: " << ::RDKit::MolToSmiles( *rdmol ) << std::endl;
		utility_exit_with_message("Encountered molecule which cannot properly be represented in RDKit.");
	}

	if ( ! keep_hydro_ ) {
		// Convert the existing physical hydrogens into explicit annotations.
		::RDKit::MolOps::removeHs(*rdmol, /*implicitOnly,*/ false, /*updateExplicitCount*/ true);
	}

	// Redo protonation on charged items.
	if ( neutralize_ ) {
		core::chemical::rdkit::final_neutralize( rdmol );
	}

	::RDKit::MolOps::assignChiralTypesFrom3D(*rdmol);

	return rdmol;
}

}
}

