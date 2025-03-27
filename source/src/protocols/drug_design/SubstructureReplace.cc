// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/SubstructureReplace.hh
/// @brief use RDKit to replace a substructure in a ResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/drug_design/SubstructureReplace.hh>
#include <protocols/drug_design/SubstructureReplaceCreator.hh>
#include <protocols/drug_design/ReactionChemistry.hh>
#include <protocols/drug_design/substitution_support.hh>

#include <protocols/chemistries/Chemistry.hh>

#include <core/chemical/rdkit/RDKit.fwd.hh>
#include <core/chemical/rdkit/RDMolToRestype.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/rdkit/util.hh>
#include <core/chemical/MutableResidueType.hh>
#include <protocols/chemistries/util.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h> // For MolToSmiles
#include <rdkit/GraphMol/MolAlign/AlignMolecules.h>
#include <rdkit/GraphMol/SanitException.h>
#include <rdkit/ForceField/ForceField.h>

namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.SubstructureReplace");

//------------------------- Creator -----------------------------

protocols::chemistries::ChemistryOP
SubstructureReplaceCreator::create_chemistry() const {
	return protocols::chemistries::ChemistryOP( new SubstructureReplace );
}

std::string
SubstructureReplaceCreator::keyname() const {
	return SubstructureReplace::class_name();
}

void
SubstructureReplaceCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SubstructureReplace::provide_xml_schema( xsd );
}

//------------------------- Chemistry -----------------------------

SubstructureReplace::SubstructureReplace():
	protocols::chemistries::Chemistry( class_name() ),
	H_as_dummy_( false ),
	V_as_dummy_( false ),
	dist_threshold_(1.0)
{}

/// @brief The file which contains the fragments to add to input residue type.
void
SubstructureReplace::substructure_database( std::string filename, bool append /*=false*/ ) {
	if ( ! append ) {
		substructures_.clear();
	}

	core::Size new_item_start( substructures_.size() + 1 );

	core::chemical::rdkit::load_sdf( filename, substructures_, /*removeHs=*/ false );

	if ( substructures_.size() < new_item_start ) {
		TR.Warning << "No molecule fragments found in file " << filename << std::endl;
	}

	for ( core::Size ii(new_item_start); ii <= substructures_.size(); ++ii ) {
		// Do H/V mangling, if needed.
		if ( H_as_dummy_ || V_as_dummy_ ) {
			::RDKit::RWMolOP mod_mol( new ::RDKit::RWMol( *substructures_[ii] ) );
			::RDKit::AtomOP qatom( new ::RDKit::QueryAtom(0) );
			qatom->setQuery( ::RDKit::makeAtomNullQuery() );
			for ( unsigned int jj(0); jj < mod_mol->getNumAtoms(); ++jj ) {
				if ( (H_as_dummy_ && mod_mol->getAtomWithIdx(jj)->getAtomicNum() == 1 ) ||
						(V_as_dummy_ && mod_mol->getAtomWithIdx(jj)->getAtomicNum() == 23 ) // Vanadium is 23
						) {
					mod_mol->replaceAtom( jj, qatom.get() ); // Atom will be copied.
				}
			}
			substructures_[ii] = mod_mol;
		}
		// Check that dummy atoms are only bonded to one other atom.
		utility::vector1< unsigned int > substruct_dummies;
		find_dummies( *substructures_[ii], substruct_dummies );
		for ( core::Size jj(1); jj <= substruct_dummies.size(); ++jj ) {
			::RDKit::ROMol::OBOND_ITER_PAIR bond_itrs( substructures_[ii]->getAtomBonds( substructures_[ii]->getAtomWithIdx( substruct_dummies[jj] ) ) );
			if ( (bond_itrs.second - bond_itrs.first) != 1 ) {
				TR.Error << "In file " << filename << " molecule '" << core::chemical::rdkit::get_name( *substructures_[ii] ) << ": ";
				TR.Error << ::RDKit::MolToSmiles( *substructures_[ii] ) << std::endl;
				TR.Error << "Dummy atom " << substruct_dummies[jj] << " has " << (bond_itrs.second - bond_itrs.first) << " bonds, should only be one." << std::endl;
				utility_exit_with_message("Dummy atom in template structure has too many bonds.");
			}
		}
	}

}

void
SubstructureReplace::apply( core::chemical::MutableResidueType & rsdtype )
{

	if ( substructures_.size() <= 1 ) {
		utility_exit_with_message("Not enough substructures found for SubstructureReplace!");
	}

	// We need hydrogens as physical entities, as we may be swapping them.
	core::chemical::rdkit::RestypeToRDMol to_converter(rsdtype, /*neutralize=*/ false, /*keep_hydro=*/ true);
	::RDKit::RWMolOP rdmol( to_converter.Mol() ); // Convert

	TR.Debug << "Number of atoms in ResType: " << rsdtype.natoms() << " in converted: " << rdmol->getNumAtoms() << std::endl;

	// Select location to do the replacement on
	MoleculeSubstitutionOP template_molsub( pick_template( rdmol, substructures_ ) );

	if ( ! template_molsub ) {
		TR.Warning << "No suitable template structure found. Doing nothing." << std::endl;
		mapping_.clear();
		mapping_.identity( true );
		set_last_status( core::chemical::modifications::FAIL_DO_NOT_RETRY );
		return;
	}

	MoleculeSubstitutionOP replacement_molsub( pick_replacement( template_molsub, substructures_, dist_threshold_, property_name_ ) );

	if ( ! replacement_molsub ) {
		TR.Warning << "No suitable replacement structure found. Doing nothing." << std::endl;
		mapping_.clear();
		mapping_.identity( true );
		set_last_status( core::chemical::modifications::FAIL_RETRY ); // RETRY, as a different template match might work.
		return;
	}

	TR.Debug << "Number of atoms in template: " << replacement_molsub->templt()->getNumAtoms() << " in replacement " << replacement_molsub->replace()->getNumAtoms() << std::endl;

	// Compute the transform from the template/replacement frame to the original molecule frame
	::RDKit::MatchVectType atom_mapping( replacement_molsub->make_match_vector() );
	RDGeom::Transform3D transform_replace_to_mol;
	::RDKit::MolAlign::getAlignmentTransform( *replacement_molsub->templt(), *replacement_molsub->mol(),
		transform_replace_to_mol, -1, -1, &atom_mapping); // Needs raw pointer

	// Now build a new molecule
	::RDKit::RWMolOP new_mol( new ::RDKit::RWMol );
	new_mol->addConformer( new ::RDKit::Conformer(0) ); // Make new zero-atom conformer - RDKit takes ownership!

	// If there is any global-level information that would need to be copied, here's the place to do it.

	// Copy over atoms from the replacement
	::RDKit::ROMolOP replace( replacement_molsub->replace() );

	for ( core::Size rdx(0); rdx < replace->getNumAtoms(); ++rdx ) { // 0 start!
		// All replacement atoms should have atom sub entries
		AtomSubstitution & atom_sub( replacement_molsub->substitution_for_rdx( rdx ) );
		if ( !replacement_molsub->rdx_is_dummy( rdx ) ) {
			// Regular atom - transfer it
			unsigned int ndx( new_mol->addAtom( replace->getAtomWithIdx( rdx ) ) );
			::RDGeom::Point3D old_pos( replace->getConformer().getAtomPos(rdx) );
			new_mol->getConformer().setAtomPos( ndx, transform_replace_to_mol * old_pos ); // transformed position
			atom_sub.set_ndx( ndx );
		} else if ( atom_sub.mdx() == AtomSubstitution::invalid_index ) {
			// Stub we don't have matched with the original molecule - turn it into a hydrogen
			unsigned int ndx( new_mol->addAtom( replace->getAtomWithIdx( rdx ) ) );
			new_mol->getAtomWithIdx(ndx)->setAtomicNum(1);
			::RDGeom::Point3D old_pos( replace->getConformer().getAtomPos(rdx) );
			new_mol->getConformer().setAtomPos( ndx, transform_replace_to_mol * old_pos ); // transformed position
			atom_sub.set_ndx( ndx );
		} // else - a dummy atom that we're replacing with a sidechain
	}

	// Copy over bonds from the replacement
	for ( core::Size bnd(0); bnd < replace->getNumBonds(); ++bnd ) { // 0 start!
		::RDKit::Bond* r_bond( replace->getBondWithIdx(bnd) );
		unsigned int bgn_rdx( r_bond->getBeginAtomIdx() );
		unsigned int end_rdx( r_bond->getEndAtomIdx() );
		unsigned int bgn_ndx( replacement_molsub->substitution_for_rdx(bgn_rdx).ndx() );
		unsigned int end_ndx( replacement_molsub->substitution_for_rdx(end_rdx).ndx() );
		if ( bgn_ndx != AtomSubstitution::invalid_index && end_ndx != AtomSubstitution::invalid_index ) {
			new_mol->addBond( bgn_ndx, end_ndx, r_bond->getBondType() );
		}
	}

	// Now we need to copy over all the sidechains from the original molecule
	// Sidechains are atoms that match a dummy atom in the replacement, and everything connected to it that doesn't match a template.
	utility::vector1< unsigned int > const & template_dummies( replacement_molsub->template_dummies());
	utility::vector1< unsigned int > mdx_atoms_to_copy;
	for ( core::Size ii(1); ii <= template_dummies.size(); ++ii ) {
		AtomSubstitution & atom_sub( replacement_molsub->substitution_for_tdx( template_dummies[ii] ) );
		if ( atom_sub.mdx() != AtomSubstitution::invalid_index && atom_sub.rdx() != AtomSubstitution::invalid_index ) {
			// Matches a replacement stub.
			mdx_atoms_to_copy.push_back( atom_sub.mdx() );
		}
	}

	for ( core::Size anum(1); anum <= mdx_atoms_to_copy.size(); ++anum ) { // Add to list while iterating - don't cache size!
		// Copy atom.
		unsigned int mdx( mdx_atoms_to_copy[ anum ] );
		AtomSubstitution & atom_sub( replacement_molsub->substitution_for_mdx( mdx ) );
		if ( atom_sub.ndx() != AtomSubstitution::invalid_index ) { // Already made a copy - skip it.
			continue;
		}
		unsigned int ndx( new_mol->addAtom( rdmol->getAtomWithIdx( mdx ) ) );
		::RDGeom::Point3D old_pos( rdmol->getConformer().getAtomPos(mdx) );
		new_mol->getConformer().setAtomPos( ndx, old_pos ); // not transformed!
		atom_sub.set_ndx( ndx );

		// Now add any bonded atoms.
		for ( ::RDKit::ROMol::OBOND_ITER_PAIR bonds( rdmol->getAtomBonds( rdmol->getAtomWithIdx(mdx) ) );
				bonds.first != bonds.second; ++bonds.first ) {
			unsigned int mdx_bnd( (*rdmol)[ *bonds.first ]->getOtherAtomIdx(mdx) );
			if ( replacement_molsub->substitution_for_mdx( mdx_bnd ).ndx() != AtomSubstitution::invalid_index // Already added
					|| replacement_molsub->substitution_for_mdx( mdx_bnd ).tdx() != AtomSubstitution::invalid_index ) { // matches the template
				continue;
			}
			mdx_atoms_to_copy.push_back( mdx_bnd ); // Copy the bonded atom
		}
	}

	// for( core::Size mdx(0); mdx < rdmol->getNumAtoms(); ++mdx ) { // 0 start!
	//  // All mol atoms should have infos
	//  AtomSubstitution & atom_sub( replacement_molsub->substitution_for_mdx( mdx ) );
	//  if( atom_sub.tdx() == AtomSubstitution::invalid_index || // Not matched to template - it's a sidechain
	//    (atom_sub.rdx() != AtomSubstitution::invalid_index &&
	//     replacement_molsub->rdx_is_dummy( atom_sub.rdx() ) )// We match a dummy in the replacement
	//   ) {
	//   // Quick spot check for hydrogens which are bonded to the core, but not matched to dummies in the template - ignore these.
	//   if( rdmol->getAtomWithIdx( mdx )->getAtomicNum() == 1 ) {
	//    unsigned int mdx_bnd( get_bonded_atom( *rdmol, mdx ) ); // Should only be a single bond to hydrogen.
	//    unsigned int tdx_bnd( replacement_molsub->substitution_for_mdx( mdx_bnd ).tdx() );
	//    if( tdx_bnd != AtomSubstitution::invalid_index && replacement_molsub->templt()->getAtomWithIdx(tdx_bnd)->getAtomicNum() != 0 ) {
	//     continue; // Ignoring hydrogen bonded to core, but not matched to template
	//    }
	//   }
	//
	//   unsigned int ndx( new_mol->addAtom( rdmol->getAtomWithIdx( mdx ) ) );
	//   ::RDGeom::Point3D old_pos( rdmol->getConformer().getAtomPos(mdx) );
	//   debug_assert( atom_sub.ndx() == AtomSubstitution::invalid_index ); // We shouldn't already have an associated atom
	//   new_mol->getConformer().setAtomPos( ndx, old_pos ); // not transformed!
	//   atom_sub.set_ndx( ndx );
	//  } // else // We're matched to the template without a correspondence to the replacement - ignore
	// }

	// Copy over bonds from original molecule
	for ( core::Size bnd(0); bnd < rdmol->getNumBonds(); ++bnd ) { // 0 start!
		::RDKit::Bond* m_bond( rdmol->getBondWithIdx(bnd) );
		unsigned int bgn_mdx( m_bond->getBeginAtomIdx() );
		unsigned int end_mdx( m_bond->getEndAtomIdx() );
		unsigned int bgn_ndx( replacement_molsub->substitution_for_mdx(bgn_mdx).ndx() );
		unsigned int end_ndx( replacement_molsub->substitution_for_mdx(end_mdx).ndx() );
		if ( bgn_ndx != AtomSubstitution::invalid_index && end_ndx != AtomSubstitution::invalid_index ) {
			debug_assert( new_mol->getBondBetweenAtoms(bgn_ndx, end_ndx) == nullptr ); // Bond should not already exist
			new_mol->addBond( bgn_ndx, end_ndx, m_bond->getBondType() );
		}
	}

	// Okay, now need to bond the sidechains from the original molecule with the core of the new one.
	// We're looking for t_dummies which are matched with r_dummies
	//utility::vector1< unsigned int > const & template_dummies( replacement_molsub->template_dummies());
	for ( core::Size ii(1); ii <= template_dummies.size(); ++ii ) {
		AtomSubstitution & atom_sub( replacement_molsub->substitution_for_tdx( template_dummies[ii] ) );
		if ( atom_sub.rdx() != AtomSubstitution::invalid_index && replacement_molsub->rdx_is_dummy( atom_sub.rdx() ) ) {
			unsigned int bnd_rdx( get_bonded_atom( *replacement_molsub->replace(), atom_sub.rdx() ) );
			// We attach the core atom on the new replacement to the original equivalent of the template dummy
			unsigned int ndx_orig( atom_sub.ndx() );
			unsigned int ndx_bnd( replacement_molsub->substitution_for_rdx(bnd_rdx).ndx() );
			::RDKit::Bond::BondType bt( get_first_bondtype(*replacement_molsub->replace(), atom_sub.rdx()) );
			debug_assert( new_mol->getBondBetweenAtoms(ndx_orig, ndx_bnd) == nullptr ); // Bond should not already exist
			new_mol->addBond( ndx_orig, ndx_bnd, bt );
		}
	}

	TR.Debug << "Pre sanitation new molecule size " << new_mol->getNumAtoms() << std::endl;
	// Now clean up the molecule.
	try {
		::RDKit::MolOps::sanitizeMol(*new_mol);
	} catch (::RDKit::MolSanitizeException &se){
		TR.Warning << "Cannot clean up molecule made by replacement - skipping." << std::endl;
		mapping_.clear();
		mapping_.identity( true );
		set_last_status( core::chemical::modifications::FAIL_RETRY ); // RETRY as a different template/replacement might work
		return;
	}

	TR.Debug << "Post sanitation new molecule size " << new_mol->getNumAtoms() << std::endl;

	// Check to make sure everything got matched up correctly. (Rosetta can't handle disconnected Restypes.)
	std::vector<int> frag_mapping;
	unsigned int num_frags( ::RDKit::MolOps::getMolFrags(*new_mol, frag_mapping) );
	if ( num_frags != 1 ) {
		TR.Error << "When substituting " << ::RDKit::MolToSmiles( *rdmol ) << std::endl;
		TR.Error << "with template " << ::RDKit::MolToSmiles( *replacement_molsub->templt() ) << std::endl;
		TR.Error << "and replacement " << ::RDKit::MolToSmiles( *replacement_molsub->replace() ) << std::endl;
		TR.Error << "disconnected molecule obtained: " << ::RDKit::MolToSmiles( *new_mol ) << std::endl;
		utility_exit_with_message("SubstructureReplace made a molecule that was disconnected!");
	}

	// Minimize with somewhat loose tolerances to fix up bond geometries
	::RDKit::ForceFieldOP ff( core::chemical::rdkit::get_forcefield( *new_mol ) );
	if ( ff ) {
		ff->minimize(200, 1e-2, 1e-4);
	} else {
		TR.Warning << "Cannot find appropriate forcefield for minimization - skipping min." << std::endl;
	}

	// We need to find a mapping from the molecule to the post molecule
	core::chemical::IndexIndexMapping sub_map( replacement_molsub->find_mdx_to_ndx_mapping() );

	// Should be good. Now convert the residue into a Rosetta residue type.

	core::chemical::VDIndexMapping restype_prod_map( combine( to_converter.vd_to_index(), sub_map ) );

	core::chemical::rdkit::RDMolToRestype from_converter(*new_mol);
	from_converter.set_nbr( restype_prod_map[rsdtype.nbr_vertex()] );

	core::chemical::MutableResidueTypeOP new_resop( from_converter.generate_restype(rsdtype,restype_prod_map) );

	TR << "Replaced the core of " << rsdtype.name() << std::endl;
	TR << "using template " << ::RDKit::MolToSmiles( *replacement_molsub->templt() ) << std::endl;
	TR << "and replacement " << ::RDKit::MolToSmiles( *replacement_molsub->replace() ) << std::endl;
	TR << "to convert " << ::RDKit::MolToSmiles( *replacement_molsub->mol() ) << std::endl;
	TR << "to " << ::RDKit::MolToSmiles( *new_mol ) << std::endl;

	mapping_ = combine( restype_prod_map, from_converter.index_to_vd() );
	rsdtype = *new_resop;
	mapping_ = combine( mapping_, combine( core::chemical::VDStringMapping(*new_resop), core::chemical::StringVDMapping(rsdtype)) );

	// That's it, we're successful
	set_last_status( core::chemical::modifications::SUCCESS );
	return;

}

core::chemical::VDVDMapping
SubstructureReplace::get_mapping() const {
	return mapping_;
}

std::string
SubstructureReplace::class_name() {
	return "SubstructureReplace";
}

void
SubstructureReplace::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("substructures", xs_string,
		"The name of the SDF file which contains the *aligned* substructures to swap.")
		+ XMLSchemaAttribute::attribute_w_default("weight_by_property", xs_string,
		"When randomly picking the substructure from the file, weight by the given property from the SDF", "")
		+ XMLSchemaAttribute::attribute_w_default("distance_threshold", xsct_real,
		"When replacing, how far apart can two attachment point atoms be and still be considered the same point", "")
		+ XMLSchemaAttribute::attribute_w_default("H_as_dummy", xsct_rosetta_bool,
		"If true, use hydrogens in the input file as attachment points", "0")
		+ XMLSchemaAttribute::attribute_w_default("V_as_dummy", xsct_rosetta_bool,
		"If true, use vanadium atoms in the input file as attachment points", "0");

	protocols::chemistries::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Replace a given substructure of a ResidueType with another, grafting substituents.",
		attlist );
}

/// @brief Initialize any data members of this instance from an input tag
/// and a DataMap object
void
SubstructureReplace::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &)
{
	// These need to be set prior to loading the database
	H_as_dummy( tag->getOption<bool>("H_as_dummy", false) );
	V_as_dummy( tag->getOption<bool>("V_as_dummy", false) );

	substructure_database( tag->getOption<std::string>("substructures") );
	weight_by_property( tag->getOption<std::string>("weight_by_property", "") );
	distance_threshold( tag->getOption<core::Real>("distance_threshold", 1.0) );
}


}
}
