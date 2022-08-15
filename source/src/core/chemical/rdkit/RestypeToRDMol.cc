// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/rdkit/RestypeToRDMol.cc
/// @brief  This class takes a residuetype from Rosetta and converts it into a RDKit RWMol object
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/rdkit/util.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/Bond.hh>

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



namespace core {
namespace chemical {
namespace rdkit {

static basic::Tracer TR("core.chemical.rdkit.RestypeToRDMol");

RestypeToRDMol::RestypeToRDMol(MutableResidueType const & res, bool neutralize /* = true*/, bool keep_hydro /*= false*/ ):
	res_( res ),
	neutralize_( neutralize ),
	keep_hydro_( keep_hydro ),
	vd_to_index_( ResidueGraph::null_vertex(), utility::get_undefined_size() )
{
	if ( neutralize_ && keep_hydro_ ) {
		utility_exit_with_message("In RestypeToRDMol, enabling neutralization while keeping the hydrogens doesn't make sense.");
	}
}


::RDKit::RWMOL_SPTR
RestypeToRDMol::Mol() {
	vd_to_index_.clear();

	::RDKit::RWMOL_SPTR rdmol( new ::RDKit::RWMol );
	rdmol->setProp("_Name", res_.name());
	rdmol->setProp("_3DConf",true);
	//rdmol->setProp("_MolFileChiralFlag",true);

	// Atoms
	std::map< core::Size, ::RDGeom::Point3D > pos_map;

	core::chemical::VIter vertex_start, vertex_end; //iterators to start traversing the graph
	boost::tie(vertex_start, vertex_end) = boost::vertices(res_.graph());
	for ( core::chemical::VIter iter_vd = vertex_start; iter_vd != vertex_end; ++iter_vd ) {
		core::chemical::Atom const & atom = res_.atom(*iter_vd);
		//RDKit handling of virtual atoms is a little iffy at this point.
		if ( atom.is_virtual() ) {
			continue;
		}

		core::Size atomic_number = atom.element_type()->get_atomic_number();

		// We keep hydrogens on the graph for now, for Kekulization reasons
		// -- with the explicit hydrogen settings, we can probably avoid that, but for now ...

		numeric::xyzVector<core::Real> coords = atom.ideal_xyz();
		int formal_charge = atom.formal_charge();

		debug_assert( atomic_number > 0 );
		// Raw pointer here is intentional - the RWMol will take exclusive ownership of the Atom
		::RDKit::Atom *rd_atom( new ::RDKit::Atom( (unsigned int) atomic_number ) );
		if ( formal_charge != 0 ) { // We take care of neutralizing later on
			rd_atom->setFormalCharge(formal_charge);
		}

		// Set atom name and index information (how well this is preserved depends on protocol).
		// Raw pointer here is intentional - the Atom will take exclusive ownership of the AtomPDBResidueInfo
		::RDKit::AtomPDBResidueInfo *info( new ::RDKit::AtomPDBResidueInfo(atom.name(),res_.atom_index(*iter_vd) ) );
		rd_atom->setMonomerInfo(info); // takes ownership

		unsigned int aid = rdmol->addAtom(rd_atom,true,true); //Will take ownership of atom

		::RDGeom::Point3D pos( coords[0], coords[1], coords[2] );
		pos_map[ aid ] = pos;

		vd_to_index_[*iter_vd] = core::Size(aid);
	}

	// Raw pointer here is intentional - the RWMol will take exclusive ownership of the conformer
	::RDKit::Conformer *conf = new ::RDKit::Conformer( rdmol->getNumAtoms() );
	conf->setId(0);
	conf->set3D(true);

	for ( std::map< core::Size, ::RDGeom::Point3D >::const_iterator itr( pos_map.begin() ), itr_end(pos_map.end()); itr != itr_end; ++itr ) {
		conf->setAtomPos(itr->first, itr->second);
	}
	rdmol->addConformer(conf, true); // Takes ownership

	// Bonds
	core::chemical::EIter edge_start, edge_end;
	boost::tie(edge_start, edge_end) = boost::edges(res_.graph());
	for ( core::chemical::EIter itr_edge = edge_start; itr_edge != edge_end; ++itr_edge ) {
		core::chemical::Bond const & bond = res_.bond(*itr_edge);

		core::chemical::VD source = boost::source( *itr_edge, res_.graph() );
		core::chemical::VD target = boost::target(*itr_edge, res_.graph());
		//Do not look at bonds between virtual atoms or virtual bonds
		if ( res_.atom(source).is_virtual() || res_.atom(target).is_virtual() || bond.order() == core::chemical::PseudoBondOrder ) {
			continue;
		}
		core::Size index1 = vd_to_index_[source];
		core::Size index2 = vd_to_index_[target];

		// Do not add bonds if the atoms aren't present in the fragment
		if ( index1 == vd_to_index_.invalid_entry() || index2 == vd_to_index_.invalid_entry() ) {
			TR.Trace << "Skipping bond between " << res_.graph()[source].name() << " and " << res_.graph()[target].name() << " due to skipped atoms." << std::endl;
			continue;
		}

		::RDKit::Bond::BondType type( convert_to_rdkit_bondtype( bond.bond_name() ) );

		if ( type == ::RDKit::Bond::AROMATIC ) {
			rdmol->getAtomWithIdx(index1)->setIsAromatic(true);
			rdmol->getAtomWithIdx(index2)->setIsAromatic(true);
		}

		rdmol->addBond( (unsigned int) index1, (unsigned int) index2, type);
	}
	// Now attempt to fill out the other parameters for this molecule

	// (Nothing Here Yet)

	// Sanitize the molecule -- the one from Rosetta should be a decent molecule - probably.
	try {
		if ( neutralize_ ) {
			neutralize_rdmol( *rdmol, /*addHs*/ false ); // Attempt to clean up silly charges.
			// Neutralization has an internal sanitization.
		} else {
			::RDKit::MolOps::sanitizeMol(*rdmol);
		}
	} catch (::RDKit::MolSanitizeException &se){
		TR.Error << "Cannot Sanitize molecule with RDKit: " << se.what() << std::endl;
		TR.Error << "    molecule: " << ::RDKit::MolToSmiles( *rdmol ) << std::endl;
		utility_exit_with_message("Encountered molecule which cannot properly be represented in RDKit.");
	}

	if ( ! keep_hydro_ ) {
		// Convert the existing physical hydrogens into explicit annotations.
		::RDKit::MolOps::removeHs(*rdmol, /*implicitOnly,*/ false, /*updateExplicitCount*/ true);

		vd_to_index_.clear();
		for ( core::Size ii(0); ii < rdmol->getNumAtoms(); ++ii ) {
			::RDKit::AtomPDBResidueInfo* ami( dynamic_cast< ::RDKit::AtomPDBResidueInfo* >(rdmol->getAtomWithIdx(ii)->getMonomerInfo()) ); // Raw pointer reference - don't delete.
			if ( ami != nullptr && ami->getSerialNumber() > 0 && core::Size(ami->getSerialNumber()) <= res_.all_atoms().size() ) {
				vd_to_index_[ res_.all_atoms()[ ami->getSerialNumber() ] ] = ii;
			} else {
				TR.Warning << "Issue with mapping atoms before/after Restype->RDKit! Atom " << ii << std::endl;
			}
		}
	}

	// Redo protonation on charged items.
	if ( neutralize_ ) {
		final_neutralize( rdmol );
	}

	::RDKit::MolOps::assignChiralTypesFrom3D(*rdmol);

	return rdmol;
}

}
}
}

