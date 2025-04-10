// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/SubstituentReplace.hh
/// @brief use RDKit to replace a subsituent on a substructure in a ResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/drug_design/SubstituentReplace.hh>
#include <protocols/drug_design/SubstituentReplaceCreator.hh>
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
#include <numeric/random/random.hh>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h> // For MolToSmiles
#include <rdkit/GraphMol/MolAlign/AlignMolecules.h>
#include <rdkit/GraphMol/SanitException.h>
#include <rdkit/ForceField/ForceField.h>

namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.SubstituentReplace");

//------------------------- Creator -----------------------------

protocols::chemistries::ChemistryOP
SubstituentReplaceCreator::create_chemistry() const {
	return protocols::chemistries::ChemistryOP( new SubstituentReplace );
}

std::string
SubstituentReplaceCreator::keyname() const {
	return SubstituentReplace::class_name();
}

void
SubstituentReplaceCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SubstituentReplace::provide_xml_schema( xsd );
}

//------------------------- Chemistry -----------------------------

SubstituentReplace::SubstituentReplace():
	protocols::chemistries::Chemistry( class_name() ),
	H_as_dummy_( false ),
	V_as_dummy_( false )
{}

void
SubstituentReplace::template_database( std::string filename, bool append /*=false*/ ) {
	if ( ! append ) {
		templates_.clear();
	}

	core::Size new_item_start( templates_.size() + 1 );

	core::chemical::rdkit::load_sdf( filename, templates_, /*removeHs=*/ false );

	if ( templates_.size() < new_item_start ) {
		TR.Warning << "No molecule fragments found in file " << filename << std::endl;
	}

	for ( core::Size ii(new_item_start); ii <= templates_.size(); ++ii ) {
		// Do H/V mangling, if needed.
		if ( H_as_dummy_ || V_as_dummy_ ) {
			::RDKit::RWMolOP mod_mol( new ::RDKit::RWMol( *templates_[ii] ) );
			::RDKit::AtomOP qatom( new ::RDKit::QueryAtom(0) );
			qatom->setQuery( ::RDKit::makeAtomNullQuery() );
			for ( unsigned int jj(0); jj < mod_mol->getNumAtoms(); ++jj ) {
				if ( (H_as_dummy_ && mod_mol->getAtomWithIdx(jj)->getAtomicNum() == 1 ) ||
						(V_as_dummy_ && mod_mol->getAtomWithIdx(jj)->getAtomicNum() == 23 ) // Vanadium is 23
						) {
					mod_mol->replaceAtom( jj, qatom.get() ); // Atom will be copied.
				}
			}
			templates_[ii] = mod_mol;
		}
		// Check that dummy atoms are only bonded to one other atom.
		utility::vector1< unsigned int > substruct_dummies;
		find_dummies( *templates_[ii], substruct_dummies );
		for ( core::Size jj(1); jj <= substruct_dummies.size(); ++jj ) {
			::RDKit::ROMol::OBOND_ITER_PAIR bond_itrs( templates_[ii]->getAtomBonds( templates_[ii]->getAtomWithIdx( substruct_dummies[jj] ) ) );
			if ( (bond_itrs.second - bond_itrs.first) != 1 ) {
				TR.Error << "In file " << filename << " molecule '" << core::chemical::rdkit::get_name( *templates_[ii] ) << ": ";
				TR.Error << ::RDKit::MolToSmiles( *templates_[ii] ) << std::endl;
				TR.Error << "Dummy atom " << substruct_dummies[jj] << " has " << (bond_itrs.second - bond_itrs.first) << " bonds, should only be one." << std::endl;
				utility_exit_with_message("Dummy atom in template structure has too many bonds.");
			}
		}
	}

}

/// @brief The file which contains the fragments to add to input residue type.
void
SubstituentReplace::substituents_database( std::string filename, bool append /*=false*/ ) {
	if ( ! append ) {
		substituents_.clear();
	}

	core::Size new_item_start( substituents_.size() + 1 );

	core::chemical::rdkit::load_sdf( filename, substituents_, /*removeHs=*/ false );

	if ( substituents_.size() < new_item_start ) {
		TR.Warning << "No molecule fragments found in file " << filename << std::endl;
	}

	// There isn't necessarily any dummies in the substituents database.
}

void
SubstituentReplace::apply( core::chemical::MutableResidueType & rsdtype )
{

	if ( templates_.size() < 1 ) {
		utility_exit_with_message("Not enough substructures found for SubstituentReplace!");
	}

	if ( substituents_.size() < 1 ) {
		utility_exit_with_message("Not enough substructures found for SubstituentReplace!");
	}

	// We need hydrogens as physical entities, as we may be swapping them.
	core::chemical::rdkit::RestypeToRDMol to_converter(rsdtype, /*neutralize=*/ false, /*keep_hydro=*/ true);
	::RDKit::RWMolOP rdmol( to_converter.Mol() ); // Convert

	// Select location to do the replacement on
	MoleculeSubstitutionOP template_molsub( pick_template( rdmol, templates_, false ) );

	if ( ! template_molsub ) {
		TR.Warning << "No suitable template structure found. Doing nothing." << std::endl;
		mapping_.clear();
		mapping_.identity( true );
		set_last_status( core::chemical::modifications::FAIL_DO_NOT_RETRY );
		return;
	}

	unsigned int tdx( AtomSubstitution::invalid_index );
	MoleculeSubstitutionOP replace_molsub( pick_replacement( template_molsub, tdx ) );

	if ( ! replace_molsub ) {
		TR.Warning << "No suitable substituent sources found for template " << ::RDKit::MolToSmiles( *template_molsub->templt() ) <<std::endl;
		TR.Warning << "  Doing nothing." << std::endl;
		mapping_.clear();
		mapping_.identity( true );
		set_last_status( core::chemical::modifications::FAIL_RETRY );
		return;
	}

	unsigned int mdx( replace_molsub->substitution_for_tdx( tdx ).mdx() );
	unsigned int rdx( replace_molsub->substitution_for_tdx( tdx ).rdx() );

	unsigned int tdx_bnd( get_bonded_atom( *replace_molsub->templt(), tdx ) );
	unsigned int mdx_bnd( replace_molsub->substitution_for_tdx( tdx_bnd ).mdx() );
	unsigned int rdx_bnd( replace_molsub->substitution_for_tdx( tdx_bnd ).rdx() );

	// Now build a new molecule
	::RDKit::RWMolOP new_mol( new ::RDKit::RWMol );
	new_mol->addConformer( new ::RDKit::Conformer(0) ); // Make new zero-atom conformer - RDKit takes ownership!

	// If there is any global-level information that would need to be copied, here's the place to do it.

	replace_molsub->add_newmol( new_mol );

	// We want all the atoms from the original molecule which are connected to the non-dummy side of the bond,
	// and everything connected to it, skipping the dummy side of the bond

	copy_attached_atoms( *replace_molsub, OriginalMol, ::RDGeom::Transform3D(), mdx_bnd, mdx );

	// Compute the transform from the replacement frame to the original molecule frame
	// The easiest is probably just to get the transform across all mapped atoms.
	::RDKit::MatchVectType atom_mapping;
	for ( unsigned int mm(0); mm < replace_molsub->mol()->getNumAtoms(); ++mm ) { // 0 Based!!
		unsigned int rr( replace_molsub->substitution_for_mdx( mm ).rdx() );
		if ( rr != AtomSubstitution::invalid_index ) {
			atom_mapping.push_back( std::pair< unsigned int, unsigned int >( rr, mm ) );
		}
	}
	::RDGeom::Transform3D transform_replace_to_mol;
	::RDKit::MolAlign::getAlignmentTransform( *replace_molsub->replace(), *replace_molsub->mol(),
		transform_replace_to_mol, -1, -1, &atom_mapping); // Needs raw pointer

	// Now we want to copy the replace atoms from the dummy side of the bond.

	copy_attached_atoms( *replace_molsub, ReplaceMol, transform_replace_to_mol, rdx, rdx_bnd );

	// Add the bond between the two fragments
	::RDKit::Bond::BondType bt( replace_molsub->mol()->getBondBetweenAtoms(mdx, mdx_bnd)->getBondType() );
	unsigned int bgn_ndx( replace_molsub->substitution_for_mdx(mdx_bnd).ndx() );
	unsigned int end_ndx( replace_molsub->substitution_for_rdx(rdx).ndx() );
	TR << "Bridging fragments with bond between " << bgn_ndx << " and " << end_ndx << std::endl;
	debug_assert( new_mol->getBondBetweenAtoms(bgn_ndx, end_ndx) == nullptr ); // Bond should not already exist
	new_mol->addBond( bgn_ndx, end_ndx, bt );

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

	// Check to make sure everything got matched up correctly. (Rosetta can't handle disconnected Restypes.)
	std::vector<int> frag_mapping;
	unsigned int num_frags( ::RDKit::MolOps::getMolFrags(*new_mol, frag_mapping) );
	if ( num_frags != 1 ) {
		TR.Error << "When substituting " << ::RDKit::MolToSmiles( *rdmol ) << std::endl;
		TR.Error << "with template " << ::RDKit::MolToSmiles( *replace_molsub->templt() ) << std::endl;
		TR.Error << "and replacement " << ::RDKit::MolToSmiles( *replace_molsub->replace() ) << std::endl;
		TR.Error << "disconnected molecule obtained: " << ::RDKit::MolToSmiles( *new_mol ) << std::endl;
		utility_exit_with_message("SubstituentReplace made a molecule that was disconnected!");
	}

	// Minimize with somewhat loose tolerances to fix up bond geometries
	::RDKit::ForceFieldOP ff( core::chemical::rdkit::get_forcefield( *new_mol ) );
	if ( ff ) {
		ff->minimize(200, 1e-2, 1e-4);
	} else {
		TR.Warning << "Cannot find appropriate forcefield for minimization - skipping min." << std::endl;
	}

	// We need to find a mapping from the molecule to the post molecule
	core::chemical::IndexIndexMapping sub_map( replace_molsub->find_mdx_to_ndx_mapping() );

	// Should be good. Now convert the residue into a Rosetta residue type.

	core::chemical::VDIndexMapping restype_prod_map( combine( to_converter.vd_to_index(), sub_map ) );

	core::chemical::rdkit::RDMolToRestype from_converter(*new_mol);
	from_converter.set_nbr( restype_prod_map[rsdtype.nbr_vertex()] );

	core::chemical::MutableResidueTypeOP new_resop( from_converter.generate_restype(rsdtype,restype_prod_map) );

	TR << "Replaced a substituent of '" << rsdtype.name() << "'" << std::endl;
	TR << "using template " << ::RDKit::MolToSmiles( *replace_molsub->templt() ) << std::endl;
	TR << "and replacement " << ::RDKit::MolToSmiles( *replace_molsub->replace() ) << std::endl;
	TR << "to convert " << ::RDKit::MolToSmiles( *replace_molsub->mol() ) << std::endl;
	TR << "to " << ::RDKit::MolToSmiles( *new_mol ) << std::endl;

	mapping_ = combine( restype_prod_map, from_converter.index_to_vd() );
	rsdtype = *new_resop;
	mapping_ = combine( mapping_, combine( core::chemical::VDStringMapping(*new_resop), core::chemical::StringVDMapping(rsdtype)) );

	// That's it, we're successful
	set_last_status( core::chemical::modifications::SUCCESS );
	return;

}

MoleculeSubstitutionOP
SubstituentReplace::pick_replacement(MoleculeSubstitutionOP template_molsub, unsigned int & tdx_out) const {

	// Subset the possible substitutuent sources to those which are applicable.
	utility::vector1< ::RDKit::ROMolOP > poss_subst;
	numeric::random::WeightedSampler subst_sampler;
	for ( core::Size jj(1); jj <= substituents_.size(); ++jj ) {
		::RDKit::MatchVectType match_vect;
		// Only need pass/fail at this point.
		if ( ::RDKit::SubstructMatch(*substituents_[jj], *template_molsub->templt(), match_vect) ) {
			core::Real weight( 1.0 );
			if ( property_name_.size() && substituents_[jj]->hasProp(property_name_) ) {
				weight = substituents_[jj]->getProp<core::Real>(property_name_);
			}
			subst_sampler.add_weight( weight );
			poss_subst.push_back( substituents_[jj] );
		}
	}

	if ( poss_subst.size() == 0 ) {
		TR.Warning << "Found no matching cores for substituent replacment. Are you missing hydrogens?" << std::endl;
		return nullptr;
	}

	while ( subst_sampler.update_cumulative_distribution() ) { // While there are valid samples.
		// Pick possible swaps from the possibilities
		core::Size selected_subst( subst_sampler.random_sample() );
		::RDKit::ROMolOP subst_source( poss_subst[ selected_subst ] );

		// Pick a template match from the subst_source (this deals with orientation issues.)
		std::vector< ::RDKit::MatchVectType > matches;
		::RDKit::SubstructMatch(*subst_source, *template_molsub->templt(), matches, /*uniquify=*/ false); // Want redundancy in orientation
		while ( matches.size() > 0 ) {
			core::Size picked_match_index( numeric::random::random_range(0,matches.size()-1) ); // All are equivalent weight
			::RDKit::MatchVectType picked_match( matches[ picked_match_index ] );
			std::map< unsigned int, unsigned int > r_to_t_map( convert_matchvect_to_map(picked_match) );
			MoleculeSubstitutionOP replace_molsub( template_molsub->add_replacement( subst_source, r_to_t_map ) );

			utility::vector1< unsigned int > possible_dummies( replace_molsub->template_dummies() ); // Make copy - we'll modify this
			while ( possible_dummies.size() > 0 ) {
				// Now pick a random dummy atom from the template
				core::Size dummy_index( numeric::random::random_range(1,possible_dummies.size()) );

				unsigned int tdx( possible_dummies[ dummy_index ] );

				unsigned int mdx( replace_molsub->substitution_for_tdx( tdx ).mdx() );
				debug_assert( mdx != AtomSubstitution::invalid_index );
				unsigned int rdx( replace_molsub->substitution_for_tdx( tdx ).rdx() );
				debug_assert( rdx != AtomSubstitution::invalid_index );

				unsigned int tdx_bnd( get_bonded_atom( *replace_molsub->templt(), tdx ) );
				debug_assert( tdx_bnd != AtomSubstitution::invalid_index );
				unsigned int mdx_bnd( replace_molsub->substitution_for_tdx( tdx_bnd ).mdx() );
				debug_assert( mdx_bnd != AtomSubstitution::invalid_index );
				unsigned int rdx_bnd( replace_molsub->substitution_for_tdx( tdx_bnd ).rdx() );
				debug_assert( rdx_bnd != AtomSubstitution::invalid_index );

				TR.Debug << "Bond on original is "    << mdx << "  -- (" << mdx_bnd <<")* " << std::endl;
				TR.Debug << "Bond on replacement is " << rdx << "* -- (" << rdx_bnd <<") " << std::endl;

				if ( ! bond_is_in_ring( *replace_molsub->mol(), mdx, mdx_bnd ) &&
						! bond_is_in_ring( *replace_molsub->replace(), rdx, rdx_bnd) &&
						( replace_molsub->mol()->getBondBetweenAtoms(mdx, mdx_bnd)->getBondType() ==
						replace_molsub->replace()->getBondBetweenAtoms(rdx, rdx_bnd)->getBondType() ) ) {
					tdx_out = tdx;
					return replace_molsub;
				}
				// This stub won't work - remove it from the list
				// (the erase-remove idiom - oh, for proper delete methods in C++)
				possible_dummies.erase( std::remove( possible_dummies.begin(), possible_dummies.end(), tdx), possible_dummies.end() );
			}
			// This template match won't work - remove it from the list
			matches.erase( matches.begin() + picked_match_index ); //For a decent index delete ...
		}
		// This molecule won't work - turn it off.
		subst_sampler.set_weight(selected_subst , 0);
	}

	TR.Warning << "Cannot find possible replacement in SubstituentReplace." << std::endl;
	return nullptr; // Can't return anything - it's not possible.
}

core::chemical::VDVDMapping
SubstituentReplace::get_mapping() const {
	return mapping_;
}

std::string
SubstituentReplace::class_name() {
	return "SubstituentReplace";
}

void
SubstituentReplace::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("templates", xs_string,
		"The name of the SDF file which contains the templates which define the substituents.")
		+ XMLSchemaAttribute::required_attribute("substituent", xs_string,
		"The name of the SDF file which contains the database of substituents to use in the replacement.")
		+ XMLSchemaAttribute::attribute_w_default("weight_by_property", xs_string,
		"When randomly picking the substructure from the file, weight by the given property from the SDF", "")
		+ XMLSchemaAttribute::attribute_w_default("H_as_dummy", xsct_rosetta_bool,
		"If true, use hydrogens in the input file as attachment points", "0")
		+ XMLSchemaAttribute::attribute_w_default("V_as_dummy", xsct_rosetta_bool,
		"If true, use vanadium atoms in the input file as attachment points", "0");

	protocols::chemistries::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Replace a substituent on a given substructure of a ResidueType with another.",
		attlist );
}

/// @brief Initialize any data members of this instance from an input tag
/// and a DataMap object
void
SubstituentReplace::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &)
{
	// These need to be set prior to loading the database
	H_as_dummy( tag->getOption<bool>("H_as_dummy", false) );
	V_as_dummy( tag->getOption<bool>("V_as_dummy", false) );

	template_database( tag->getOption<std::string>("templates") );
	substituents_database( tag->getOption<std::string>("substituent") );
	weight_by_property( tag->getOption<std::string>("weight_by_property", "") );

	TR << "Defined SubstituentReplace Chemistry with " << templates_.size() << " templates and " << substituents_.size() << " substituents, weighting by property '" << weight_by_property() << "'" << std::endl;
}


}
}
