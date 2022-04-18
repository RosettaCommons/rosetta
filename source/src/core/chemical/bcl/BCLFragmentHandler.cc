// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//////////////////////////////////////////////////////////////////////
/// @file
///
/// @brief
/// A class for converting between BCl and Rosetta small molecule objects
///
/// @details
/// This class converts between BCL and Rosetta small molecule objects and sets
/// small molecule BCL conformers as rotamers in a Rosetta pose residue.
/// This class is meant to deprecate FragmentToRestype and RestypeToFragment
/// and make one umbrella class that can be extended and/or generalized.
///
/// @author
/// Sandeep Kothiwale
/// Rocco Moretti
/// Steven Combs
/// Benjamin P. Brown (benjamin.p.brown17@gmail.com)
///
////////////////////////////////////////////////////////////////////////

// unit include
#include <core/chemical/bcl/BCLFragmentHandler.hh>

// core includes
#include <core/chemical/Atom.hh>
#include <core/chemical/atomtype_support.hh>
#include <core/chemical/bcl/util.hh>
#include <core/chemical/Bond.hh>
#include <core/chemical/bond_support.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/Elements.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>
#include <core/chemical/icoor_support.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/MutableResidueConnection.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/residue_support.hh>

// util includes
#include <utility/exit.hh>
#include <utility/numbers.hh>
#include <utility/string_util.hh>

// basic includes
#include <basic/Tracer.hh>

// numeric includes
#include <numeric/xyzVector.hh>

// BCL includes
#ifdef USEBCL
#include <bcl/include/linal/bcl_linal_vector_3d.h>
#include <bcl/include/storage/bcl_storage_vector.h>
#include <bcl/include/storage/bcl_storage_list.h>
#endif

namespace core
{
namespace chemical
{
namespace bcl
{

static basic::Tracer TR("core.chemical.bcl.BCLFragmentHandler");

//////////////////////////////////
// construction and destruction //
//////////////////////////////////

BCLFragmentHandler::BCLFragmentHandler() :
	utility::VirtualBase(),
	index_to_vd_(utility::get_undefined_size(), ResidueGraph::null_vertex()),
	vd_to_index_( ResidueGraph::null_vertex(), utility::get_undefined_size()),
	nbr_(utility::get_undefined_size())
{
#ifdef USEBCL
	bcl_fragment_ = ::bcl::chemistry::FragmentComplete();
#endif
}

#ifdef USEBCL
/// @brief construct with BCL fragment
BCLFragmentHandler::BCLFragmentHandler
(
	::bcl::chemistry::FragmentComplete const &fragment
) :
	utility::VirtualBase(),
	bcl_fragment_( fragment),
	index_to_vd_( utility::get_undefined_size(), ResidueGraph::null_vertex()),
	vd_to_index_( ResidueGraph::null_vertex(), utility::get_undefined_size()),
	nbr_( utility::get_undefined_size())
{
	// extras=bcl required for construction of this class object
	core::chemical::bcl::require_bcl();
}
#endif

#ifdef USEBCL
/// @brief construct with Rosetta residue
BCLFragmentHandler::BCLFragmentHandler
(
	MutableResidueTypeCOP restype
) :
	utility::VirtualBase(),
	bcl_fragment_( ::bcl::chemistry::FragmentComplete()),
	rosetta_restype_( restype),
	index_to_vd_( utility::get_undefined_size(), ResidueGraph::null_vertex()),
	vd_to_index_( ResidueGraph::null_vertex(), utility::get_undefined_size()),
	nbr_( utility::get_undefined_size())
{
	// extras=bcl required for construction of this class object
	core::chemical::bcl::require_bcl();
}
#endif


/////////////////
// data access //
/////////////////

#ifdef USEBCL
//! Return the BCL molecule as a standard FragmentComplete
::bcl::chemistry::FragmentComplete const &BCLFragmentHandler::get_bcl_fragment() const
{
	return bcl_fragment_;
}
#endif

//! @brief Return a Rosetta molecule as a MutableResidueType
MutableResidueTypeCOP BCLFragmentHandler::get_rosetta_restype() const
{
	return rosetta_restype_;
}

/// @brief Get how the most recently created ResidueType corresponds to the underlying fragment.
IndexVDMapping const &BCLFragmentHandler::get_index_to_vd() const
{
	return index_to_vd_;
}

//! @brief Get mapping of restype vertex descriptors to indices of the bcl fragments
VDIndexMapping const &BCLFragmentHandler::get_vd_to_index() const
{
	return vd_to_index_;
}

//! @brief Get the atom used as the neighbor atom when restype is generated
core::Size BCLFragmentHandler::get_nbr() const
{
	return nbr_;
}

////////////////
// operations //
////////////////

#ifdef USEBCL
//! @brief Set the BCL molecule
void BCLFragmentHandler::set_bcl_fragment( ::bcl::chemistry::FragmentComplete const &fragment)
{
	bcl_fragment_ = fragment;
}
#endif

//! @brief Set the Rosetta molecule
void BCLFragmentHandler::set_rosetta_restype( MutableResidueTypeCOP restype)
{
	rosetta_restype_ = restype;
}

/// @brief Which atom in the fragment to use as the neighbor atom when the a restype is generated.
void BCLFragmentHandler::set_nbr( core::Size nbr )
{
	nbr_ = nbr;
}

#ifdef USEBCL
//! @brief Convert BCL member Fragment to Rosetta MutableResidueType
//! @return mutable residue type object for Rosetta
MutableResidueTypeOP BCLFragmentHandler::fragment_to_restype()
{
	return fragment_to_restype( bcl_fragment_);
}
#endif

#ifdef USEBCL
//! @brief Convert BCL Fragment to Rosetta MutableResidueType using default restype information
//! @param FRAGMENT the BCL FragmentComplete to be converted
//! @param CONFS fragment conformers to save as rotamers
//! @return mutable residue type object for Rosetta
MutableResidueTypeOP BCLFragmentHandler::fragment_to_restype
(
	::bcl::chemistry::FragmentComplete const &fragment,
	::bcl::chemistry::FragmentEnsemble const &confs
)
{
	// obtain default standard atom typing and related Rosetta atom info
	ChemicalManager * cm(core::chemical::ChemicalManager::get_instance());
	std::string const tag("fa_standard");
	AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
	ElementSetCOP element_types = cm->element_set("default");
	MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);
	orbitals::OrbitalTypeSetCOP orbital_types = cm->orbital_type_set(tag);

	// setup mutable residue type object to which BCL fragment will be converted
	MutableResidueType restype(atom_types, element_types, mm_atom_types, orbital_types);
	restype.name( "LIGAND" );
	restype.name3( "LIG" );
	restype.name1( 'X' );
	restype.interchangeability_group( "LIG" );

	// initialize an empty mapping
	VDIndexMapping mapping;

	// convert
	return fragment_to_restype(fragment, restype, mapping, confs);
}
#endif

#ifdef USEBCL
//! @brief Convert BCL Fragment to Rosetta MutableResidueType with additional information from provided restype
//! @param FRAGMENT the BCL FragmentComplete to be converted
//! @param RESTYPE the the restype from which to obtain additional restype info
//! @param MAPPING the VD index mapping object
//! @param CONFS fragment conformers to save as rotamers
//! @return mutable residue type object for Rosetta
MutableResidueTypeOP BCLFragmentHandler::fragment_to_restype
(
	::bcl::chemistry::FragmentComplete const &fragment,
	MutableResidueType const &restype,
	VDIndexMapping const &mapping,
	::bcl::chemistry::FragmentEnsemble const &confs
)
{
	// obtain restype information from provided ResidueType
	MutableResidueTypeOP residue( utility::pointer::make_shared< core::chemical::MutableResidueType>(
		restype.atom_type_set_ptr(),
		restype.element_set_ptr(),
		restype.mm_atom_types_ptr(),
		restype.orbital_types_ptr()) );

	// set to gasteiger atom types
	residue->set_gasteiger_atom_typeset( restype.gasteiger_atom_typeset());
	if ( ! residue->gasteiger_atom_typeset() ) {
		ChemicalManager * cm(core::chemical::ChemicalManager::get_instance());
		residue->set_gasteiger_atom_typeset( cm->gasteiger_atom_type_set("default") );
	}

	// set residue name
	residue->name(  restype.name() );
	residue->name3( restype.name3() );
	residue->name1( restype.name1() );
	residue->interchangeability_group( restype.interchangeability_group() );

	// clear index matching
	index_to_vd_.clear();

	// map BCL fragment atominfo to each restype atom
	size_t atom_index( 0); // size_t used for compatibility with BCL
	::bcl::storage::Vector< ::bcl::sdf::AtomInfo> atom_info( fragment.GetAtomInfo());
	for
		(
				::bcl::storage::Vector< ::bcl::sdf::AtomInfo>::const_iterator
				itr_atom( atom_info.Begin()), itr_atom_end( atom_info.End());
				itr_atom != itr_atom_end;
				++itr_atom, ++atom_index
				) {
		const std::string element_name = itr_atom->GetAtomType()->GetElementType()->GetChemicalSymbol();
		const short charge = itr_atom->GetAtomType()->GetFormalCharge(); // return type on GetFormalCharge is a short
		std::string gasteiger_name = itr_atom->GetAtomType().GetName();

		// Try copying the atom name
		std::string name = "";
		VD orig_vd( mapping.reverse_lookup(atom_index) );
		if ( orig_vd != mapping.invalid_key() ) {
			name = restype.atom_name( orig_vd );
		}

		// obtain atomic coordinates
		VD atom_vd = residue->add_atom(name);
		Atom & restype_atom( residue->atom( atom_vd ) );
		numeric::xyzVector<core::Real> xyz_coords;
		xyz_coords.x() = itr_atom->GetCoordinates().X();
		xyz_coords.y() = itr_atom->GetCoordinates().Y();
		xyz_coords.z() = itr_atom->GetCoordinates().Z();

		// assign element charge info to restype
		restype_atom.element_type( residue->element_set().element(element_name) );
		restype_atom.charge( charge );
		restype_atom.formal_charge( charge );
		restype_atom.ideal_xyz( xyz_coords );
		restype_atom.mm_name( "VIRT" ); // We need to do better on this typing.

		// allow unknown atom types
		if ( ! residue->gasteiger_atom_typeset()->contains_atom_type(gasteiger_name ) ) {
			gasteiger_name = "";
		}
		residue->set_gasteiger_atom_type(atom_vd, gasteiger_name);
		index_to_vd_[atom_index] = atom_vd;
	} // end setting atom info


	// map BCL fragment atominfo to between each restype atom pair
	::bcl::storage::Vector< ::bcl::sdf::BondInfo> bond_info( fragment.GetBondInfo());
	for (
			::bcl::storage::Vector< ::bcl::sdf::BondInfo>::const_iterator
			itr_bond( bond_info.Begin()), itr_bond_end( bond_info.End());
			itr_bond != itr_bond_end;
			++itr_bond, ++atom_index
			) {
		// create the bond type object
		// Previously we were only using raw bond order to go between Rosetta and BCL bond types, but in practice this results in some
		// distorted amide bonds and sub-optimal ring geometries. Rosetta does not currently have an AmideBond type, so until we figure out
		// an alternative plan to mapping bonds, I am going to just have this little logic thing here that manually assigns
		// Rosetta bond types based on the BCL enum name; in other words, DO NOT change the BCL bond type enum
		// I repeat, DO NOT change the BCL bond type enum used in this class
		// NOTE - this may end up being score function-dependent. This works reasonably well for now from what I can tell
		// (and much better than before for amide bonds), but ultimately this is a scientific question we should discuss so that
		// we can improve it further.
		std::string bcl_bond_enum_name( ::bcl::util::Format()( itr_bond->GetConfigurationalBondType().GetName()));
		core::chemical::BondName bond_name;
		// single bonds
		if
				(
						bcl_bond_enum_name == "AromaticSingleBond" ||
						bcl_bond_enum_name == "NonConjugatedSingleBond" || bcl_bond_enum_name == "ConjugatedSingleBond" ||
						bcl_bond_enum_name == "NonConjugatedSingleBondInRing" || bcl_bond_enum_name == "ConjugatedSingleBondInRing" // ||
						//   bcl_bond_enum_name == "AmideSingleBond"
						) {
			// create the bond
			residue->add_bond
				(
				index_to_vd_[ itr_bond->GetAtomIndexLow()],
				index_to_vd_[ itr_bond->GetAtomIndexHigh()],
				bond_name = core::chemical::BondName( SingleBond)
			);
		} else if
		// double bonds
				(
						bcl_bond_enum_name == "AromaticDoubleBond" ||
						bcl_bond_enum_name == "ConjugatedDoubleBond" || bcl_bond_enum_name == "ConjugatedDoubleBondInRing" ||
						bcl_bond_enum_name == "ConjugatedDoubleBond_X" || bcl_bond_enum_name == "ConjugatedDoubleBond_E" ||
						bcl_bond_enum_name == "ConjugatedDoubleBond_Z" || bcl_bond_enum_name == "ConjugatedDoubleBondInRing_X" ||
						bcl_bond_enum_name == "ConjugatedDoubleBondInRing_E" || bcl_bond_enum_name == "ConjugatedDoubleBondInRing_Z"
						) {
			// create the bond
			residue->add_bond
				(
				index_to_vd_[ itr_bond->GetAtomIndexLow()],
				index_to_vd_[ itr_bond->GetAtomIndexHigh()],
				bond_name = core::chemical::BondName( DoubleBond)
			);
		} else if
		// triple bonds
				(
						bcl_bond_enum_name == "AromaticTripleBond" ||
						bcl_bond_enum_name == "ConjugatedTripleBond" || bcl_bond_enum_name == "ConjugatedTripleBondInRing"
						) {
			// create the bond
			residue->add_bond
				(
				index_to_vd_[ itr_bond->GetAtomIndexLow()],
				index_to_vd_[ itr_bond->GetAtomIndexHigh()],
				bond_name = core::chemical::BondName( TripleBond)
			);
		} else if
		// kekulizing everything; treat amide bonds as aromatic to enforce planar geometry
				(
						//   bcl_bond_enum_name == "AromaticSingleBond" || bcl_bond_enum_name == "AromaticDoubleBond" ||
						//   bcl_bond_enum_name == "AromaticTripleBond" || bcl_bond_enum_name == "AromaticBond" ||
						bcl_bond_enum_name == "AmideSingleBond"
						) {
			// create the bond
			residue->add_bond
				(
				index_to_vd_[ itr_bond->GetAtomIndexLow()],
				index_to_vd_[ itr_bond->GetAtomIndexHigh()],
				bond_name = core::chemical::BondName( AromaticBond)
			);
		} else {
			// if it is something else then assign UnknownBond and throw a warning;
			// the remaining bond types in BCL are generic e_ConjugatedBond and
			// e_ConjugatedBondInRing with no specified order; this will be treated as Unknown in Rosetta
			// but then recalculated in the BCL
			// create the bond
			TR << "Warning! Bond type conversion from BCL to Rosetta may be inaccurate!" << std::endl;
			residue->add_bond
				(
				index_to_vd_[ itr_bond->GetAtomIndexLow()],
				index_to_vd_[ itr_bond->GetAtomIndexHigh()],
				bond_name = core::chemical::BondName( UnknownBond)
			);
		}
	}

	// finalize
	core::chemical::rename_atoms(*residue, /*preserve=*/true);
	rosetta_retype_fullatom(*residue, false);
	rosetta_recharge_fullatom(*residue);
	core::chemical::find_bonds_in_rings( *residue );
	VD nbr_atom = MutableResidueType::null_vertex;
	if ( ! utility::is_undefined( nbr_ ) ) {
		nbr_atom = index_to_vd_[ nbr_ ];
	}
	residue->nbr_radius( find_nbr_dist( *residue, nbr_atom ) );
	residue->nbr_atom( nbr_atom );
	residue->assign_internal_coordinates(); // Also sets atom base. Needs nbr atom assignment
	residue->autodetermine_chi_bonds(); // Must be after internal coordinate setup

	// if there are no conformers to save as rotamers, just return the new restype
	if ( !confs.GetSize() ) {
		return residue;
	}

	// save conformers as rotamers on the new restype
	core::chemical::rotamers::StoredRotamerLibrarySpecificationOP
		rotamers_spec( utility::pointer::make_shared< core::chemical::rotamers::StoredRotamerLibrarySpecification >() );

	// iterate over all conformers
	for
		(
				auto itr_mols( confs.Begin()), itr_mols_end( confs.End());
				itr_mols != itr_mols_end;
				++itr_mols
				) {
		std::map< std::string, core::Vector > single_rotamer_spec;

		::bcl::util::SiPtrVector< const ::bcl::linal::Vector3D> atom_coords( itr_mols->GetAtomCoordinates());

		//keep track of what the atom number we are on
		core::Size i=0;
		for
			(
					::bcl::util::SiPtrVector< const ::bcl::linal::Vector3D>::const_iterator
					itr_coords( atom_coords.Begin()), itr_coords_end( atom_coords.End());
					itr_coords != itr_coords_end;
					++itr_coords, ++i
					) {
			//convert BCL fragment index base to Rosetta's VD based atoms
			//This relies on the fact that the BCL manipulation doesn't change indices
			VD const vd = index_to_vd_[i];

			//get the xyz coordinates
			core::Real x, y, z;
			x = ( *itr_coords)->X();
			y = ( *itr_coords)->Y();
			z = ( *itr_coords)->Z();

			//set the xyz coordinates in Rosetta
			single_rotamer_spec[ residue->atom( vd).name() ] = core::Vector( x, y, z );
		}
		// add the current conformer as a rotamer
		rotamers_spec->add_rotamer( single_rotamer_spec );
	}
	// add the whole Rosetta rotamer library to the residue type
	residue->rotamer_library_specification( rotamers_spec );
	return residue;
}
#endif

#ifdef USEBCL
//! @brief Convert ResidueType to a BCL Fragment.
//! @param RESTYPE the restype to be converted into a FragmentComplete
//! @param NORMALIZE if true, remove charges and reprotonate based on BCL heuristics.
//! Note, charges due to heavy atoms, like quaternary amines, should still be present.
::bcl::chemistry::FragmentComplete BCLFragmentHandler::restype_to_fragment
(
	MutableResidueType const &restype,
	bool normalize
)
{
	// initialize atominfo vector that we will fill with atom data from restype
	std::vector< ::bcl::sdf::AtomInfo > atominfo_vector;

	//iterators to start traversing the graph
	core::chemical::VIter vertex_start, vertex_end;
	boost::tie(vertex_start, vertex_end) = boost::vertices( restype.graph());

	// 0 based indexing in BCL
	core::Size index( 0);

	// iterate over all atoms
	for ( core::chemical::VIter iter_vd( vertex_start); iter_vd != vertex_end; ++iter_vd ) {
		// current atom
		core::chemical::Atom const &atom = restype.atom(*iter_vd);

		// BCL does not use virtual atoms, so skip them
		if ( atom.is_virtual() ) {
			continue;
		}

		// no need to normalize hydrogen atoms
		std::string element = atom.element_type()->get_chemical_symbol();
		if ( normalize && element == "H" ) {
			continue;
		}

		// set atom coordinates and charge
		vd_to_index_[*iter_vd] = index;
		numeric::xyzVector<core::Real> coords = atom.ideal_xyz();
		core::Real charge = atom.charge();
		if ( normalize ) {
			charge = 0;
		}

		// fill the atom info vector using the data generated above
		const ::bcl::chemistry::ElementType element_type( ::bcl::chemistry::GetElementTypes().ElementTypeLookup(element) );
		::bcl::chemistry::AtomType type( ::bcl::chemistry::AtomTypes::GetAtomType( element_type, charge));
		::bcl::linal::Vector3D coordinates( coords[0], coords[1], coords[2] );
		::bcl::sdf::AtomInfo atom_info(type, ::bcl::chemistry::e_UnknownChirality, coordinates);
		atominfo_vector.push_back( atom_info);
		++index;
	}

	// note the use of std::vector for BCL objects; 0-indexing is standard for BCL, kept here
	std::vector< ::bcl::sdf::BondInfo> bondinfo_vector;
	utility::vector1<core::chemical::Bond> bonds;
	core::chemical::EIter edge_start, edge_end;
	boost::tie(edge_start, edge_end) = boost::edges(restype.graph());
	for ( core::chemical::EIter itr_edge( edge_start); itr_edge != edge_end; ++itr_edge ) {
		core::chemical::Bond const & bond = restype.bond( *itr_edge);
		bonds.push_back(bond);
		core::chemical::VD source = boost::source( *itr_edge, restype.graph() );
		core::chemical::VD target = boost::target(*itr_edge, restype.graph());
		// Do not look at bonds between virtual atoms or virtual bonds
		if
				(
						restype.atom(source).is_virtual() ||
						restype.atom(target).is_virtual() || bond.order() == core::chemical::PseudoBondOrder
						) {
			continue;
		}
		core::Size const index1 = vd_to_index_[source];
		core::Size const index2 = vd_to_index_[target];

		// Do not add bonds if the atoms aren't present in the fragment
		if
				(
						index1 == vd_to_index_.invalid_entry() ||
						index2 == vd_to_index_.invalid_entry() ) {
			TR.Debug << "Skipping bond between " << restype.graph()[source].name() << " and "
				<< restype.graph()[target].name() << " due to skipped atoms." << std::endl;
			continue;
		}

		// convert bond order to size_t for the bcl
		int order = (int)bond.bond_name();
		if ( order < 6 ) {
			// in the bcl all of the bond order/conjugation/ringness/aromaticity stuff is contained within
			// a single enum, while in Rosetta they are split out. while not an explicit item to-do,
			// this function could afford more logic to allow more detailed bond type parsing
			// NOTE - these will match the BondOrderAmideOrAromatic enum, but are only initial types;
			// atom standardization will adjust atom types and the corresponding bond types
			::bcl::chemistry::ConfigurationalBondType name_to_bcl_type[6] = {
				::bcl::chemistry::ConfigurationalBondType(),
				::bcl::chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
				::bcl::chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond_X,
				::bcl::chemistry::GetConfigurationalBondTypes().e_ConjugatedTripleBond,
				::bcl::chemistry::GetConfigurationalBondTypes().e_AmideSingleBond, // Rosetta does not have Amide bond
				::bcl::chemistry::GetConfigurationalBondTypes().e_AromaticBond
				};
			::bcl::sdf::BondInfo bond_info(index1, index2, name_to_bcl_type[ order]);
			bondinfo_vector.push_back(bond_info);
		}
	}

	// build our complete atomvector from the atom and bond info vectors
	::bcl::chemistry::AtomVector< ::bcl::chemistry::AtomComplete> atoms( atominfo_vector, bondinfo_vector);

	// determine atom and bond types
	const std::string mol_id( "Rosetta ResidueType from BCL " );
	::bcl::chemistry::AtomsCompleteStandardizer standardizer( atoms, mol_id, true); // force recalculation

	// add bond stereocenter/isometry info then build full molecule
	::bcl::chemistry::BondIsometryHandler::AddIsometryInformation( atoms, true);
	::bcl::chemistry::StereocentersHandler::AddChiralityFromConformation( atoms);
	::bcl::chemistry::FragmentComplete new_molecule(atoms, "");

	if ( normalize ) {
		new_molecule.SaturateWithH();
	}

	// done, return the BCL fragment
	return new_molecule;
}
#endif

} // bcl
} // chemical
} // core

