// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
// This file is based off of code from the Biochemistry Library (BCL).
// The BCL is copyright Vanderbilt University (Meiler Lab), a Rosetta Commons Member Institution


/// @file   core/chemical/gasteiger/GasteigerAtomTyper.cc
/// @brief  The type assigner for gasteiger type data.
/// @author Rosetta conversion: Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/gasteiger/GasteigerAtomTyper.hh>

#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/ElectronConfiguration.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/Bond.hh>
#include <core/chemical/ResidueConnection.hh>

#include <core/chemical/sdf/mol_writer.hh>

#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <utility/graph/ring_detection.hh>

#include <utility/numbers.hh>

#include <basic/Tracer.hh>

#include <map>
#include <string>

namespace core {
namespace chemical {
namespace gasteiger {

static THREAD_LOCAL basic::Tracer TR( "core.chemical.gasteiger.GasteigerAtomTypeSet" );

void
PossibleAtomTypesForAtom::gasteiger_atom_type_set( GasteigerAtomTypeSetCOP GASTEIGER_ATOM_TYPE_SET ) {
	debug_assert( GASTEIGER_ATOM_TYPE_SET );
	gasteiger_atom_type_set_ = GASTEIGER_ATOM_TYPE_SET;
}

/// @brief store the atomic environment for each atom type in a map for fast lookup when
/// determining atom types
std::map< std::string, PossibleAtomTypesForAtom >
PossibleAtomTypesForAtom::CreateAtomicEnvironmentToTypesMap ( const bool IN_AROMATIC_RING, GasteigerAtomTypeSetCOP GASTEIGER_ATOM_TYPE_SET ) {
	std::map< std::string, PossibleAtomTypesForAtom > atomic_environment_to_possible_types;
	for ( core::Size ii(1); ii <= GASTEIGER_ATOM_TYPE_SET->n_types(); ++ii ) {
		GasteigerAtomTypeDataCOP const & type( (*GASTEIGER_ATOM_TYPE_SET)[ii] );
		if ( ! type->is_gasteiger_atom_type() ) { // skip non-gasteiger type
			TR << "Skipping Non-Gasteiger type." << std::endl;
			continue;
		}

		if ( !IN_AROMATIC_RING ) {
			const std::string search_str
				(
				type->get_element_type()->get_chemical_symbol()
				+ utility::to_string( type->get_number_electrons_in_bonds())
				+ utility::to_string( type->get_number_bonds())
			);

			// add the type to the map, ignoring charge, since it is redundant with the # bonds and e- in bonds
			atomic_environment_to_possible_types[ search_str].AddAtomType( type);
		} else if ( type->is_conjugated() && type->get_hybrid_orbital_type() != GasteigerAtomTypeData::Unhybridized ) {
			//  else if( type->get_number_bonds() < 2 )
			//  {
			//   // must be at least two bonds to be inside a ring
			// Not true - open babel can type aromatic even if it's not in a ring.
			//   continue;
			//  }
			// depending on what the atom in the aromatic ring is connected to, it may or may not have an extra
			// e- based considering aromatic bond types to have 3 e- in bonds
			const std::string search_str
				(
				type->get_element_type()->get_chemical_symbol() + utility::to_string( type->get_number_bonds() )
			);

			// add search strings for type in an aromatic ring with charges -1 - 1  to the map
			for
				(
						core::Size min_double_bonds( 0), max_double_bonds( type->get_number_electrons_in_bonds() - type->get_number_bonds());
						min_double_bonds <= max_double_bonds;
						++min_double_bonds
						) {
				// the type may be chosen when an anion was requested if its charge is <= 0
				// an atom type with < 0 charge will be preferred, however
				if ( type->get_formal_charge() <= 0 ) {
					atomic_environment_to_possible_types[ search_str + utility::to_string( min_double_bonds) + "N"].AddAromaticAtomType( type, short( -1));
				}

				// since charges are often omitted, all atom types are considered if a neutral atom type was requested
				// in a ring with declared aromatic bonds
				atomic_environment_to_possible_types[ search_str + utility::to_string( min_double_bonds) + "O"].AddAromaticAtomType( type, short( 0));

				// the type may be chosen when a cation was requestion if its charge is > 0
				// an atom type with > 0 charge will be preferred, however
				if ( type->get_formal_charge() >= 0 ) {
					atomic_environment_to_possible_types[ search_str + utility::to_string( min_double_bonds) + "P"].AddAromaticAtomType( type, short( 1));
				}
			}
		}
	}

	TR << "Removing unhybridized atoms for " << (IN_AROMATIC_RING?"aromatic":"non-aromatic")
		<< " map with " << atomic_environment_to_possible_types.size() << " entries for GasteigerAtomTypeSet with " << GASTEIGER_ATOM_TYPE_SET->n_types() << " entries. " << std::endl;

	// Remove unhybridized types wherever possible; the unhybridized types are only used when there is no
	// hybridized alternative
	for
		(
				std::map< std::string, PossibleAtomTypesForAtom >::iterator
				itr( atomic_environment_to_possible_types.begin()),
				itr_end( atomic_environment_to_possible_types.end());
				itr != itr_end;
				++itr
				) {
		itr->second.gasteiger_atom_type_set( GASTEIGER_ATOM_TYPE_SET );
		if ( IN_AROMATIC_RING ) {
			// get the last character, which indicates the desired charge
			const char last_char( itr->first[ itr->first.size() - 1]);
			const int charge( last_char - 'O');
			itr->second.FinalizeAromatic( charge);
		} else {
			// finalize; only keep the best type
			itr->second.Finalize();
		}
	}

	TR << "Initializing Gasteiger atom type map for " << (IN_AROMATIC_RING?"aromatic":"non-aromatic")
		<< " with " << atomic_environment_to_possible_types.size() << " entries for GasteigerAtomTypeSet with " << GASTEIGER_ATOM_TYPE_SET->n_types() << " entries. " << std::endl;
	return atomic_environment_to_possible_types;
}

////////////////////////////////////
//// construction and destruction //
////////////////////////////////////

PossibleAtomTypesForAtom::PossibleAtomTypesForAtom() :
	m_NumberAtomTypesWithHybridization( GasteigerAtomTypeData::NumberHybridOrbitalType, 0 ),
	m_AtomTypesByDecreasingStability(),
	m_NumberConjugatedTypes( 0),
	m_FinalizeFunction( 0 )
{
}

//! @brief constructor from the known information about the atom
//! @param ELEMENT element type,
//! @param NUMBER_ELECTRONS_IN_BONDS number of electrons in bonds for the atom type; for atom in declared aromatic environment, number of exocyclic bonds
//! @param NUMBER_BONDS number of bonds for the atom
//! @param SUSPECTED_CHARGE; expected charge, ignored if no atom type matching the other criteria if found
//! @param IN_AROMATIC_RING true iff the atom has bonds of the aromatic unspecified type
PossibleAtomTypesForAtom
PossibleAtomTypesForAtom::FindPossibleAtomTypesForAtom
(
	GasteigerAtomTypeSetCOP gasteiger_atom_type_set,
	const Element &ELEMENT,
	const core::Size NUMBER_ELECTRONS_IN_BONDS,
	const core::Size NUMBER_BONDS,
	const int SUSPECTED_CHARGE,
	const bool IN_AROMATIC_RING
) {
	// create maps from atomic environment to possible types
	typedef std::map< std::string, PossibleAtomTypesForAtom> EnvTypesMap;
	typedef std::map< GasteigerAtomTypeSetCOP, EnvTypesMap > SetToEnvTypesMap;
	static SetToEnvTypesMap non_aro_typing, aro_typing;

	// The less than and equivalent semantics for owning pointers are for the underlying raw pointers,
	// so this should key on object identity
	if ( non_aro_typing.find(gasteiger_atom_type_set) == non_aro_typing.end() ) {
		non_aro_typing.insert( SetToEnvTypesMap::value_type( gasteiger_atom_type_set, CreateAtomicEnvironmentToTypesMap( false, gasteiger_atom_type_set  )));
		aro_typing.insert( SetToEnvTypesMap::value_type( gasteiger_atom_type_set, CreateAtomicEnvironmentToTypesMap( true, gasteiger_atom_type_set )));
	}

	EnvTypesMap const & s_atomic_env_outside_arom_ring_to_types_map = non_aro_typing.find(gasteiger_atom_type_set)->second;
	EnvTypesMap const & s_element_bonds_in_arom_ring_to_types_map = aro_typing.find(gasteiger_atom_type_set)->second;

	//static const std::map< std::string, PossibleAtomTypesForAtom>
	//s_atomic_env_outside_arom_ring_to_types_map( CreateAtomicEnvironmentToTypesMap( false ) );
	//static const std::map< std::string, PossibleAtomTypesForAtom>
	//s_element_bonds_in_arom_ring_to_types_map( CreateAtomicEnvironmentToTypesMap( true ) );

	std::map< std::string, PossibleAtomTypesForAtom>::const_iterator itr;
	std::string primary_search_string;

	if ( IN_AROMATIC_RING ) {
		// when aromatic bonds (sdf id = 4) are given, 1 of 2 different numbers of electrons in bonds are possible,
		// depending on the charge.  Thus, first search for the type with the specified charge
		primary_search_string =
			(
			ELEMENT.get_chemical_symbol()
			+ utility::to_string( NUMBER_BONDS)
			+ utility::to_string( NUMBER_ELECTRONS_IN_BONDS - NUMBER_BONDS - 1)
			+ std::string( 1, char( 'O' + SUSPECTED_CHARGE))
		);

		itr = s_element_bonds_in_arom_ring_to_types_map.find( primary_search_string);
		if ( itr != s_element_bonds_in_arom_ring_to_types_map.end() ) {
			// use the type with known charge
			return itr->second;
		}
	} else {

		primary_search_string =
			(
			ELEMENT.get_chemical_symbol()
			+ utility::to_string( NUMBER_ELECTRONS_IN_BONDS)
			+ utility::to_string( NUMBER_BONDS)
		);
		// look for the type in the type-outside-ring map
		itr = s_atomic_env_outside_arom_ring_to_types_map.find( primary_search_string);
		if ( itr != s_atomic_env_outside_arom_ring_to_types_map.end() ) {
			return itr->second;
		}
	}
	// TODO: We probably need to have a better fall-back catching system here.
	TR.Error << "Error when getting all gasteiger atom types for " << ELEMENT.get_chemical_name() << " with " << NUMBER_BONDS << " bonds." << std::endl;
	TR.Error << "Search string '" << primary_search_string << "' not found in "
		<< (IN_AROMATIC_RING?"aromatic":"non-aromatic") << " database having "
		<< (IN_AROMATIC_RING?s_element_bonds_in_arom_ring_to_types_map.size():s_atomic_env_outside_arom_ring_to_types_map.size()) << " entries." << std::endl;
	//utility_exit_with_message("No gasgteiger atom types exist for atom.");
	return PossibleAtomTypesForAtom(); // To silence warnings
}

///////////////////
//// data access //
///////////////////

bool PossibleAtomTypesForAtom::CouldHaveHybridization( const GasteigerAtomTypeData::HybridOrbitalType HYBRID) const
{
	return m_NumberAtomTypesWithHybridization[ HYBRID ];
}

core::Size PossibleAtomTypesForAtom::GetNumberPossibleTypes() const
{
	return m_AtomTypesByDecreasingStability.size();
}

GasteigerAtomTypeDataCOP PossibleAtomTypesForAtom::GetMostStableType() const
{
	return !m_AtomTypesByDecreasingStability.empty()
		? m_AtomTypesByDecreasingStability.front() : GasteigerAtomTypeDataCOP(0);
}

////! @brief get the alternate atom types
////! @return the alternative atom types
//storage::Vector< AtomType> PossibleAtomTypesForAtom::GetAlternateTypes() const
//{
//  return
//    storage::Vector< AtomType>( ++m_AtomTypesByDecreasingStability.Begin(), m_AtomTypesByDecreasingStability.End());
//}
//
////! @brief get the alternate atom type with the given charge
////! @param CHARGE the charge desired
////! @return an alternative atom type
//AtomType PossibleAtomTypesForAtom::GetAlternateTypeWithCharge( const short &CHARGE) const
//{
//  for
//  (
//    storage::List< AtomType>::const_iterator
//      itr( m_AtomTypesByDecreasingStability.Begin()), itr_end( m_AtomTypesByDecreasingStability.End());
//    itr != itr_end;
//    ++itr
//  )
//  {
//    if( ( *itr)->get_formal_charge() == CHARGE)
//    {
//      return *itr;
//    }
//  }
//  return GetAtomTypes().Undefined;
//}
//
//bool PossibleAtomTypesForAtom::CouldBeConjugated() const
//{
//  // return true if any conjugated types are available
//  return m_NumberConjugatedTypes;
//}
//
//bool PossibleAtomTypesForAtom::MustBeConjugated() const
//{
//  // return true if there is at least one conjugated type and the number of conjugated types
//  // is the same as the number of types
//  return m_NumberConjugatedTypes && m_NumberConjugatedTypes == m_AtomTypesByDecreasingStability.GetSize();
//}
//
////! @brief determine the maximal # of pi-electrons in the pi-electron system
////! @return the maximal # of pi-electrons in the pi-electron system
//size_t PossibleAtomTypesForAtom::GetMaxElectronsParticipatingInPiSystem() const
//{
//  storage::List< AtomType>::const_iterator itr( m_AtomTypesByDecreasingStability.Begin()),
//                                           itr_end( m_AtomTypesByDecreasingStability.End());
//
//  while( itr != itr_end && !( *itr)->is_conjugated())
//  {
//    ++itr;
//  }
//
//  // if we found a non-conjugated type before the end, then return false, otherwise, return true
//  return itr == itr_end ? 0 : ( *itr)->get_maxE_contribution_to_pi_system();
//}
//
//////////////////
//// operations //
//////////////////

//! @brief add an atom type to be considered
//! @param ATOM_TYPE the type of atom to consider
void PossibleAtomTypesForAtom::AddAtomType( GasteigerAtomTypeDataCOP ATOM_TYPE)
{
	const GasteigerAtomTypeData::HybridOrbitalType hybrid( ATOM_TYPE->get_hybrid_orbital_type());

	// handle the common case where there are currently no atom types
	// also always add atom types for unhybridized types and atom types with only 1 bond
	// which are handled separately when Finalize is called
	if
			(
					m_AtomTypesByDecreasingStability.empty()
					|| hybrid == GasteigerAtomTypeData::Unhybridized
					|| ATOM_TYPE->get_number_bonds() == 1
					) {
		++m_NumberAtomTypesWithHybridization[ hybrid ];
		if ( ATOM_TYPE->is_conjugated() ) {
			++m_NumberConjugatedTypes;
		}
		m_AtomTypesByDecreasingStability.push_back( ATOM_TYPE);
		return;
	}

	// handle the complicated case
	const core::Real stability( ATOM_TYPE->get_stability_metric());

	if ( CouldHaveHybridization( hybrid) ) { // find the type with that hybridization, test its stability
		// test stability criterion to decide whether we should insert this hybrid
		std::list< GasteigerAtomTypeDataCOP >::iterator itr( m_AtomTypesByDecreasingStability.begin());

		// find the hybrid orbital type in the most stable type vector
		// we know this type is there, otherwise CouldHaveHybridization would have returned false
		while ( ( *itr)->get_hybrid_orbital_type() != hybrid )
				{
			++itr;
		}

		// test the stability
		if ( stability > ( *itr)->get_stability_metric() ) {
			// ATOM_TYPE is more stable than the existing type with that hybridization, so replace it
			if ( ( *itr)->is_conjugated() ) {
				--m_NumberConjugatedTypes;
			}
			m_AtomTypesByDecreasingStability.erase( itr);
			--m_NumberAtomTypesWithHybridization[ hybrid ];
		}
	}

	// only insert if there are no types remaining with that hybridization
	if ( !CouldHaveHybridization( hybrid) ) {
		++m_NumberAtomTypesWithHybridization[ hybrid ];
		std::list< GasteigerAtomTypeDataCOP >::iterator itr( m_AtomTypesByDecreasingStability.begin()),
			itr_end( m_AtomTypesByDecreasingStability.end());

		while ( itr != itr_end && ( *itr)->get_stability_metric() > stability )
				{
			++itr;
		}

		m_AtomTypesByDecreasingStability.insert( itr, ATOM_TYPE);
		if ( ATOM_TYPE->is_conjugated() ) {
			++m_NumberConjugatedTypes;
		}
	}
}

//! @brief set this object to only consider the given atom type
//! @param ATOM_TYPE the atom type desired
void PossibleAtomTypesForAtom::SetToType( GasteigerAtomTypeDataCOP ATOM_TYPE)
{
	if ( std::find (
			m_AtomTypesByDecreasingStability.begin(),
			m_AtomTypesByDecreasingStability.end(),
			ATOM_TYPE) == m_AtomTypesByDecreasingStability.end() ) {
		TR.Warning << "WARNING: Tried to limit atom type to " + ATOM_TYPE->get_name()
			+ " but could not find it " << std::endl;
		TR.Warning << "Possible types are: ";
		for ( std::list< GasteigerAtomTypeDataCOP >::iterator itr(m_AtomTypesByDecreasingStability.begin()),
				end(m_AtomTypesByDecreasingStability.end()); itr != end; ++itr ) {
			TR.Warning << " " << (*itr)->get_name();
		}
		if ( m_AtomTypesByDecreasingStability.size() == 0 ) {
			TR.Warning << " (None) ";
		}
		TR.Warning << std::endl;
	}

	// reset this object
	m_NumberAtomTypesWithHybridization.clear(); //
	m_NumberAtomTypesWithHybridization.resize(GasteigerAtomTypeData::NumberHybridOrbitalType, 0);
	m_AtomTypesByDecreasingStability.clear();
	m_NumberConjugatedTypes = 0;

	// add the atom type
	++m_NumberAtomTypesWithHybridization[ ATOM_TYPE->get_hybrid_orbital_type()];
	if ( ATOM_TYPE->is_conjugated() ) {
		++m_NumberConjugatedTypes;
	}
	m_AtomTypesByDecreasingStability.push_back( ATOM_TYPE);
	m_FinalizeFunction = NULL;
}

//! @brief set the final type based on the given atom and smallest ring size
//! @param ATOM the atom of interest
//! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
void PossibleAtomTypesForAtom::Finalize(const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD)
{
	// check whether there is any need for finalization
	if ( m_FinalizeFunction ) {
		// use the defined finalize function
		( this->*m_FinalizeFunction)(graph, atomVD);
	}
}

/////////////////
//// operators //
/////////////////
//
////////////////////////
//// input and output //
////////////////////////
//
//std::istream &PossibleAtomTypesForAtom::Read( std::istream &ISTREAM)
//{
//  io::Serialize::Read( m_AtomTypesByDecreasingStability, ISTREAM);
//  m_NumberAtomTypesWithHybridization.SetAllElements( 0);
//  m_NumberConjugatedTypes = 0;
//  for
//  (
//    storage::List< AtomType>::const_iterator
//      itr( m_AtomTypesByDecreasingStability.Begin()), itr_end( m_AtomTypesByDecreasingStability.End());
//    itr != itr_end;
//    ++itr
//  )
//  {
//    ++m_NumberAtomTypesWithHybridization( ( *itr)->get_hybrid_orbital_type().GetIndex());
//    if( ( *itr)->is_conjugated())
//    {
//      ++m_NumberConjugatedTypes;
//    }
//  }
//
//  return ISTREAM;
//}
//
//std::ostream &PossibleAtomTypesForAtom::Write( std::ostream &OSTREAM, const size_t INDENT) const
//{
//  io::Serialize::Write( m_AtomTypesByDecreasingStability, OSTREAM, INDENT);
//
//  return OSTREAM;
//}
//
////////////////////////
//// helper functions //
////////////////////////

core::Size PossibleAtomTypesForAtom::hybridization_rank( GasteigerAtomTypeData::HybridOrbitalType const hybrid ) {
	switch ( hybrid ) {
	case GasteigerAtomTypeData::Unhybridized :
		return 3;
	case GasteigerAtomTypeData::SP :
		return 2;
	case GasteigerAtomTypeData::SP2 :
		return 0;
	case GasteigerAtomTypeData::SP3 :
		return 1;
	default :
		TR << "Attempting to get rank for insufficiently specified bond." << std::endl;
		return 4;
	}
}

//! @brief remove a particular hybrid orbital type from the possible types, unless that would remove all possibilities
//! @param HYBRID the type of hybrid orbital to remove
void PossibleAtomTypesForAtom::RemoveHybridization( const GasteigerAtomTypeData::HybridOrbitalType HYBRID)
{
	// if there are no alternatives to HYBRID for this type, just return
	if ( m_NumberAtomTypesWithHybridization[ HYBRID ] == GetNumberPossibleTypes() ) {
		return;
	}

	for
		(
				std::list< GasteigerAtomTypeDataCOP >::iterator
				itr( m_AtomTypesByDecreasingStability.begin()), itr_end( m_AtomTypesByDecreasingStability.end());
				itr != itr_end;
				// iteration in loop
				) {
		if ( ( *itr)->get_hybrid_orbital_type() == HYBRID ) {
			if ( ( *itr)->is_conjugated() ) {
				--m_NumberConjugatedTypes;
			}
			itr = m_AtomTypesByDecreasingStability.erase( itr);
			--m_NumberAtomTypesWithHybridization[ HYBRID ];
		} else { // incorrect hybridization,
			++itr;
		}
	}
}

void PossibleAtomTypesForAtom::AddAromaticAtomType( GasteigerAtomTypeDataCOP ATOM_TYPE, const int DESIRED_CHARGE)
{
	const GasteigerAtomTypeData::HybridOrbitalType hybrid( ATOM_TYPE->get_hybrid_orbital_type());
	const int charge_diff( std::abs( ATOM_TYPE->get_formal_charge() - DESIRED_CHARGE));

	// only add types with smaller than 2 difference in charge; changing charge by 2 would mean that the atom
	// was given a -1 charge and that we would replace it with a type that would have a +1 charge, or vice versa
	if ( charge_diff > 1 ) {
		return;
	}

	++m_NumberConjugatedTypes;
	++m_NumberAtomTypesWithHybridization[ hybrid ];

	// handle the complicated case
	const core::Real stability( ATOM_TYPE->get_stability_metric());

	// find out where this atom type should be placed
	std::list< GasteigerAtomTypeDataCOP >::iterator itr( m_AtomTypesByDecreasingStability.begin());

	const core::Size hybrid_rank( hybridization_rank( hybrid ) );

	// place itr at the first atom type that ATOM_TYPE should be considered before in search order
	for
		(
				std::list< GasteigerAtomTypeDataCOP>::const_iterator itr_end( m_AtomTypesByDecreasingStability.end());
				itr != itr_end;
				++itr
				) {
		const int itr_charge_diff( std::abs( ( *itr)->get_formal_charge() - DESIRED_CHARGE));

		if ( charge_diff < itr_charge_diff ) {
			break;
		}
		if ( charge_diff > itr_charge_diff ) {
			continue;
		}

		const size_t itr_hybrid_rank( hybridization_rank( ( *itr)->get_hybrid_orbital_type() ));
		if ( hybrid_rank < itr_hybrid_rank ) {
			break;
		}
		if ( hybrid_rank > itr_hybrid_rank ) {
			continue;
		}

		// equal hybridization ranks and charge difference; go for stability
		if ( stability > ( *itr)->get_stability_metric() ) {
			break;
		}
	}
	m_AtomTypesByDecreasingStability.insert( itr, ATOM_TYPE);
}

//! @brief Select the best choice for the atom type wherever possible
//! @see @link https://structbio.vanderbilt.edu:8443/display/MeilerLab/RethinkingAtomTypeDetection @endlink
//! @details the link above contains the statistics and models used to select the current set of rules
void PossibleAtomTypesForAtom::Finalize()
{
	using namespace core::chemical;

	if ( m_AtomTypesByDecreasingStability.size() <= 1 ) {
		return;
	}

	// extract # bonds, element type, main group
	const core::Size n_bonds( GetMostStableType()->get_number_bonds());
	const ElementCOP element_type( GetMostStableType()->get_element_type());
	const core::Size main_group( element_type->get_main_group());

	// handle group I and II elements; prefer unhybridized types with bonding S orbitals
	// handle group VII elements; prefer unhybridized types with bonding P orbitals
	if ( n_bonds <= 2 && ( main_group == 1 || main_group == 2 || main_group == 7) ) {
		// Hybridization generally does not benefit these elements, which generally only form ionic bonds
		FinalizeUnhybridized();
		return;
	}

	// handle remaining elements.  These will all prefer to be hybridized, so remove unhybridized choices first
	RemoveHybridization( GasteigerAtomTypeData::Unhybridized );
	if ( m_AtomTypesByDecreasingStability.size() == 1 ) {
		return;
	}

	// get # of electrons in bonds
	const core::Size n_e_in_bonds( GetMostStableType()->get_number_electrons_in_bonds());

	// determine preferred hybridization
	GasteigerAtomTypeData::HybridOrbitalType preferred_hybridization = GasteigerAtomTypeData::Unhybridized;

	if ( n_bonds == 1 ) { // handle remaining cases with a single bond
		// Formally charged atoms with a single bond
		// choose hybrid orbital type according to VSEPR number (bonds + lone pairs)
		preferred_hybridization = GasteigerAtomTypeData::HybridOrbitalType( 4 - n_e_in_bonds);
	} else if ( n_bonds >= 4 || ( element_type->get_period() >= 3 && n_e_in_bonds == n_bonds) ) {
		// with 4+ bonds, only Te types are available under the standard gasteiger types
		// period 3+ atom types (except group 1, 2, and 7) engage in d-orbital bonding, which is not represented
		// in gasteiger atom types.  Choose the most similar type by bond angles (SP3)
		preferred_hybridization = GasteigerAtomTypeData::SP3;
	} else if ( n_e_in_bonds - n_bonds >= 2 ) {
		// at least one triple bond or two double bonds
		// prefer digonal types if there are two bonds, trigonal if there are three
		if ( n_bonds == 2 ) {
			preferred_hybridization = GasteigerAtomTypeData::SP;
		} else { // if( n_bonds == 3)
			preferred_hybridization = GasteigerAtomTypeData::SP2;
		}
	} else if ( element_type->get_period() >= size_t( 3) ) {
		// period 3+, group 3-6 elements with 2-3 bonds (with the current atom types, this can only be P, S, and Si),
		// exactly one double bond
		// period 3+ atom types (except group 1, 2, and 7) engage in d-orbital bonding, which is not represented
		// in gasteiger atom types.

		// under the assumption that similar bond angles reflect the most similar types, all types in this category
		// are closest to tetrahedral except Si with 3 bonds (1 double), which is closest to trigonal
		// On the other hand, none of the candidate atom types are tetrahedral since Te cannot form unsaturated bonds,
		// it is only through d-orbitals that this is achieved for period 3+ elements.
		// Likewise, prefer to eliminate all types except Tr and Te;  will assert if any new types are added
		// that violate this rule.
		preferred_hybridization = GasteigerAtomTypeData::SP3;
	} else if ( n_bonds == 3 && n_e_in_bonds == 4 ) {
		// 1 unsaturated bond; two single bonds.  The pi-orbital assures planarity
		preferred_hybridization = GasteigerAtomTypeData::SP2;
	} else if ( element_type->element() == element::B || element_type->element() == element::C ) {
		// down to B, C, N, O, with either 2 single bonds, 1 single & 1 double bond, or 3 single bonds
		// O, and N atom types depend heavily on their environment
		// boron and carbon with 2 - 3 bonds are always trigonal
		preferred_hybridization = GasteigerAtomTypeData::SP2;
	} else if ( element_type->element() == core::chemical::element::N ) {
		// down to N and O, with either 2 single bonds, 1 single & 1 double bond, or 3 single bonds
		if ( n_bonds == 3 && n_e_in_bonds == 3 ) {
			m_FinalizeFunction = &PossibleAtomTypesForAtom::FinalizeNitrogenThreeSingle;
		} else if ( n_bonds == 2 && n_e_in_bonds == 3 ) {
			m_FinalizeFunction = &PossibleAtomTypesForAtom::FinalizeNitrogenSingleDouble;
		} else if ( n_bonds == 2 && n_e_in_bonds == 2 ) {
			m_FinalizeFunction = &PossibleAtomTypesForAtom::FinalizeNitrogenTwoSingle;
		}
	} else if ( element_type->element() == element::O ) {
		if ( n_bonds == 3 && n_e_in_bonds == 3 ) {
			m_FinalizeFunction = &PossibleAtomTypesForAtom::FinalizeOxygenThreeSingle;
		} else if ( n_bonds == 2 && n_e_in_bonds == 3 ) {
			m_FinalizeFunction = &PossibleAtomTypesForAtom::FinalizeOxygenSingleDouble;
		} else if ( n_bonds == 2 && n_e_in_bonds == 2 ) {
			m_FinalizeFunction = &PossibleAtomTypesForAtom::FinalizeOxygenTwoSingle;
		}
	}

	if ( preferred_hybridization == GasteigerAtomTypeData::SP ) {
		RemoveHybridization( GasteigerAtomTypeData::SP3);
		RemoveHybridization( GasteigerAtomTypeData::SP2);
	} else if ( preferred_hybridization == GasteigerAtomTypeData::SP2 ) {
		RemoveHybridization( GasteigerAtomTypeData::SP);
		RemoveHybridization( GasteigerAtomTypeData::SP3);
	} else if ( preferred_hybridization == GasteigerAtomTypeData::SP3 ) {
		RemoveHybridization( GasteigerAtomTypeData::SP);
		RemoveHybridization( GasteigerAtomTypeData::SP2);
	}

	if ( !m_FinalizeFunction && m_AtomTypesByDecreasingStability.size() > 1 ) {
		TR.Error << "Cannot assign Gasteiger atom types: residual unhandled degeneracy: ";
		for
			(
					std::list< GasteigerAtomTypeDataCOP >::const_iterator
					itr( m_AtomTypesByDecreasingStability.begin()), itr_end( m_AtomTypesByDecreasingStability.end());
					itr != itr_end;
					++itr
					) {
			TR.Error << ' ' << (*itr)->get_name();
		}
		TR.Error << std::endl;
		utility_exit_with_message( "Cannot assign gasteiger atom types: residual unhandled degeneracy" );
	}
}

//! @brief choose the preferred atom type (using VSEPR theory) assuming that the orbitals do not hybridize
//! @details This is used for elements in group 1, 2, & 7, which do hybridize in the gasteiger scheme
void PossibleAtomTypesForAtom::FinalizeUnhybridized()
{
	RemoveHybridization( GasteigerAtomTypeData::SP3);
	RemoveHybridization( GasteigerAtomTypeData::SP2);
	RemoveHybridization( GasteigerAtomTypeData::SP);
	if ( m_AtomTypesByDecreasingStability.size() <= 1 ) {
		return;
	}

	// ensure that only two types remain that differ only based on whether they bond with a sigma or just p orbitals
	runtime_assert_string_msg( m_AtomTypesByDecreasingStability.size() == 2,
		"Should have been exactly two choices remaining for unhybridized types");

	// choose the type with bonding sigma orbitals, based on normal valence electron / Lewis shell model
	GasteigerAtomTypeDataCOP & atom_type_0_bonding_s_orbitals( m_AtomTypesByDecreasingStability.front());
	GasteigerAtomTypeDataCOP & atom_type_1_bonding_s_orbital( m_AtomTypesByDecreasingStability.back());
	if ( atom_type_0_bonding_s_orbitals->get_number_unhybridized_sigma_orbitals() ) {
		runtime_assert_string_msg(
			!atom_type_1_bonding_s_orbital->get_number_unhybridized_sigma_orbitals(),
			"One of the two unhybridized types should have not had a bonding sigma orbital"
		);
		std::swap( atom_type_0_bonding_s_orbitals, atom_type_1_bonding_s_orbital);
	} else {
		runtime_assert_string_msg(
			atom_type_1_bonding_s_orbital->get_number_unhybridized_sigma_orbitals(),
			"One of the two unhybridized types should have had a bonding sigma orbital"
		);
	}
	if ( atom_type_0_bonding_s_orbitals->get_element_type()->get_main_group() == 7 ) {
		// the two hybridizations available are SP2P2P2 and S2P2P2P; prefer the latter,
		// based on normal valence electron / Lewis shell model
		SetToType( atom_type_0_bonding_s_orbitals);
	} else {
		// the hybridizations available are S/P or SP/PP.  Prefer the type with that bonds with s-orbitals (S/SP), which
		// are lower energy
		SetToType( atom_type_1_bonding_s_orbital);
	}
}

//! @brief only keep the most stable types for the atom that span the set of desired pi orbital electrons (0-2)
//! @param DESIRED_CHARGE the preferred charge
//! used during construction of the maps when there is no part of standardization that
//! should edit this class
void PossibleAtomTypesForAtom::FinalizeAromatic( const int DESIRED_CHARGE)
{
	// to work in any aromatic system, one of several different types may be necessary to satisfy hueckels rule
	// thus, always choose the first type in the list with the given number of pi electrons that has the best
	// hybridization order of any of the atom types with that hybridization
	std::map< core::Size, GasteigerAtomTypeDataCOP > rank_to_best_type;

	// foreach x electrons in pi system
	//   - Prefer getting something closer to desired charge
	//   - If distance from desired charge is the same, prefer getting something with better hybridization rank
	//   - All else being equal, prefer the more stable compound
	for ( core::Size pi_electrons( 0), limit_pi_electrons( 3); pi_electrons < limit_pi_electrons; ++pi_electrons ) {
		core::Size best_hybrid_orbital_rank( 10);
		int best_charge_diff( 10);
		GasteigerAtomTypeDataCOP best_type( 0 );
		core::Size index( 0), best_index( utility::get_undefined_size() );
		std::list< GasteigerAtomTypeDataCOP > equivalent_choices;
		for ( std::list< GasteigerAtomTypeDataCOP >::const_iterator
				itr( m_AtomTypesByDecreasingStability.begin()), itr_end( m_AtomTypesByDecreasingStability.end());
				itr != itr_end;
				++itr, ++index ) {
			if ( ( *itr)->get_maxE_contribution_to_pi_system() != pi_electrons &&
					( pi_electrons != 0 || ( *itr)->get_pi_electron_contribution_type() != GasteigerAtomTypeData::ZeroOrTwo) ) {
				continue;
			}

			const int itr_charge_diff( std::abs( ( *itr)->get_formal_charge() - DESIRED_CHARGE));
			const core::Size hybrid_orbital_rank( hybridization_rank( ( *itr)->get_hybrid_orbital_type() ));
			if ( hybrid_orbital_rank < best_hybrid_orbital_rank && itr_charge_diff <= best_charge_diff ) {
				best_hybrid_orbital_rank = hybrid_orbital_rank;
				best_charge_diff = itr_charge_diff;
				best_type = *itr;
				best_index = index;
				equivalent_choices.clear();
			} else if ( hybrid_orbital_rank == best_hybrid_orbital_rank && itr_charge_diff == best_charge_diff ) {
				equivalent_choices.push_back( *itr);
			}
		}

		if ( best_type ) {
			rank_to_best_type[ best_index] = best_type;
		}
	}

	m_AtomTypesByDecreasingStability.clear();
	for
		(
				std::map< core::Size, GasteigerAtomTypeDataCOP>::const_iterator
				itr( rank_to_best_type.begin()), itr_end( rank_to_best_type.end());
				itr != itr_end;
				++itr
				) {
		m_AtomTypesByDecreasingStability.push_back( itr->second);
	}

	m_NumberConjugatedTypes = m_AtomTypesByDecreasingStability.size();
}

//! @brief choose the final atom type for Nitrogen with two single bonds
void PossibleAtomTypesForAtom::FinalizeNitrogenTwoSingle( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD ) {
	using namespace core::chemical;

	// in ring of size 3 - 5 -> Tetrahedral (51 cases)
	// Bound to N or S -> Tetrahedral (90 cases)
	// Else -> Trigonal (53 cases)
	core::Size SMALLEST_RING_SIZE = utility::graph::smallest_ring_size(atomVD, graph, 5);
	if ( SMALLEST_RING_SIZE >= 3 && SMALLEST_RING_SIZE <= 5 ) {
		SetToType( gasteiger_atom_type_set_->atom_type("N_Te2Te2TeTe") );
	} else {
		// check for bonds to N or S
		RealResidueAdjacentIter iter, iter_end;
		for ( boost::tie(iter,iter_end) = boost::adjacent_vertices( atomVD, graph ); iter != iter_end; ++iter ) {
			const ElementCOP & element_data = graph[ *iter ].element_type();
			debug_assert( element_data );
			if ( element_data->element() == element::S || element_data->element() == element::N ) {
				SetToType( gasteiger_atom_type_set_->atom_type("N_Te2Te2TeTe"));
				return;
			}
		}
		SetToType( gasteiger_atom_type_set_->atom_type("N_Tr2TrTrPi2"));
	}
}

//! @brief choose the final atom type for Nitrogen with three single bonds
void PossibleAtomTypesForAtom::FinalizeNitrogenThreeSingle( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD ) {
	// Node 1 In a 3 membered ring Tetrahedral (Accuracy: 100%, 293 cases)
	// Bond orbital overlap likely precludes trigonal
	core::Size SMALLEST_RING_SIZE = utility::graph::smallest_ring_size(atomVD, graph, 5);
	if ( SMALLEST_RING_SIZE == 3 ) {
		SetToType( gasteiger_atom_type_set_->atom_type("N_Te2TeTeTe"));
		return;
	}
	const core::Size unsaturated_neighbor_count( CountUnsaturatedNeighbors( graph, atomVD ));
	if ( unsaturated_neighbor_count >= 2 ) {
		// Node 2 2+ unsaturated neighbors Trigonal (Accuracy: 99.5%, 11070 cases)
		// 2+ unsaturated neighbors guarantees that maximal orbital overlap will be achieved with trigonal geometry
		SetToType( gasteiger_atom_type_set_->atom_type("N_TrTrTrPi2"));
		return;
	} else if ( IsBondedToAHalogen( graph, atomVD ) ) {
		// Node 3 Connected to any halogen Tetrahedral (Accuracy: 88.9%, 45 cases)
		// Lone pair repulsion and strong electronegativity induces extra P character on N
		SetToType( gasteiger_atom_type_set_->atom_type("N_Te2TeTeTe"));
		return;
	} else if ( unsaturated_neighbor_count == 1 ) {
		// Node 4 1 unsaturated neighbor Trigonal (Accuracy: 94.1%, 14440 cases)
		// A single unsaturated neighbor is otherwise sufficient to grant trigonality in most cases
		// Exceptions are generally small (4-5 membered), fused or bridged rings, in which case bond angle strain can be severe, but presumably does not change the hybridization of N
		SetToType( gasteiger_atom_type_set_->atom_type("N_TrTrTrPi2"));
		return;
	}

	bool all_C = true;
	RealResidueAdjacentIter iter, iter_end;
	for ( boost::tie(iter,iter_end) = boost::adjacent_vertices( atomVD, graph ); iter != iter_end; ++iter ) {
		const ElementCOP & element_data = graph[ *iter ].element_type();
		if ( element_data->element() != element::C ) {
			all_C = false;
		}
		// Node 5 Connected to any element period > 2 except those in group 6 Trigonal (Accuracy: 98.1%, 327 cases)
		// This is likely due to the relative similarity of d & p orbitals, since elements in period > 2 appear to engage
		// heavily in d-orbital bonding
		if ( element_data->get_period() > 2 && ( element_data->get_main_group() != 6 ) ) {
			// as a rule, H-O-H is tetrahedral; also holds for any halogens, N, and O (Node 2)
			SetToType( gasteiger_atom_type_set_->atom_type("N_TrTrTrPi2"));
			return;
		}
	}
	// Node 6a All C neighbors Tetrahedral (Accuracy: 96.4%, 5981 cases) No impetus to become planar
	if ( all_C ) {
		SetToType( gasteiger_atom_type_set_->atom_type("N_Te2TeTeTe"));
		return;
	}

	// check for being in multiple, non-aromatic rings (aromatic rings were already checked
	// for by looking for unsaturated neighbors)
	if ( false ) { // TODO: We actually do need a ring enumeration check: ATOM.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)) == size_t( 3))
		// Node 6b in multiple rings (Accuracy: 96.4%, 5981 cases) No impetus to become planar
		SetToType( gasteiger_atom_type_set_->atom_type("N_Te2TeTeTe"));
	} else if ( SMALLEST_RING_SIZE == 6 ) {
		RealResidueAdjacentIter iter, iter_end;
		for ( boost::tie(iter,iter_end) = boost::adjacent_vertices( atomVD, graph ); iter != iter_end; ++iter ) {
			const ElementCOP & element_data = graph[ *iter ].element_type();
			// Node 7: In a 6 membered ring (Accuracy: 70.4%, 252 cases)
			if ( element_data->element() == element::O ) {
				// Node 7a at least 1 O Trigonal (181 cases)
				SetToType( gasteiger_atom_type_set_->atom_type("N_TrTrTrPi2"));
				return;
			}
		}
		// Node 7b no O Tetrahedral (71 cases)
		SetToType( gasteiger_atom_type_set_->atom_type("N_Te2TeTeTe"));
	} else {
		// Node 8 Everything else (Accuracy: 78.6%, 564 cases)
		size_t heteroatom_count( 0);
		RealResidueAdjacentIter iter, iter_end;
		for ( boost::tie(iter,iter_end) = boost::adjacent_vertices( atomVD, graph ); iter != iter_end; ++iter ) {
			const ElementCOP & element_data = graph[ *iter ].element_type();
			if ( element_data->element() != element::C ) { // TODO: Do we need a Hydrogen check here, too?
				++heteroatom_count;
			}
		}
		if ( heteroatom_count >= 2 ) {
			// Node 8a: Connected to 2+ heteroatoms Trigonal (71 cases)
			// Heteroatoms are generally quite electronegative and can form conjugated systems
			SetToType( gasteiger_atom_type_set_->atom_type("N_TrTrTrPi2"));
		} else {
			// Node 8b: Connected to 1 heteroatom (493 cases)
			// Very little impetus to become planar
			SetToType( gasteiger_atom_type_set_->atom_type("N_Te2TeTeTe"));
		}
	}
}

//! @brief choose the final atom type for a single and a double bond
//! @param ATOM the atom of interest
//! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
void PossibleAtomTypesForAtom::FinalizeNitrogenSingleDouble( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD ) {
	// use the digonal type only if bound to an atom w/ period 3+ and
	// whose main group (5 for N, 6 for O) is <= the bonded atom, and is at least 4
	// For example, Al should not induce digonality
	RealResidueAdjacentIter iter, iter_end;
	for ( boost::tie(iter,iter_end) = boost::adjacent_vertices( atomVD, graph ); iter != iter_end; ++iter ) {
		const ElementCOP & element_data = graph[ *iter ].element_type();
		if ( element_data->get_period() > 2 && ( element_data->get_main_group() == 4 || element_data->get_main_group() == 5) ) {
			SetToType( gasteiger_atom_type_set_->atom_type("N_DiDiPi2Pi"));
			return;
		}
	}
	SetToType( gasteiger_atom_type_set_->atom_type("N_Tr2TrTrPi"));
}

//! @brief choose the final atom type for Oxygen with two single bonds
//! @param ATOM the atom of interest
//! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
void PossibleAtomTypesForAtom::FinalizeOxygenTwoSingle( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD ) {
	// any empty valences imply hydrogens; hydrogens imply tetrahedral
	if ( boost::out_degree( atomVD, graph ) < 2 ) {
		SetToType( gasteiger_atom_type_set_->atom_type("O_Te2Te2TeTe"));
		return;
	}
	// check for aromatic ring membership
	// only the first bond needs to be checked
	if ( graph[ *(boost::out_edges( atomVD, graph).first) ].aromaticity() == IsAromaticBond ) {
		SetToType( gasteiger_atom_type_set_->atom_type("O_Tr2TrTrPi2"));
		return;
	}
	//RM: Because aromatic typing may not be complete, check for unsaturated neighbors.
	// If we have two unsaturated neighbors, we're likely to be Trigonal.
	if ( CountUnsaturatedNeighbors( graph, atomVD ) >= 2 ) {
		SetToType( gasteiger_atom_type_set_->atom_type("O_Tr2TrTrPi2"));
		return;
	}
	//
	// Node 1 In a 3-5 membered ring Tetrahedral (Accuracy: 99.5%, 10210 cases)
	// Bond orbital overlap likely precludes trigonal
	core::Size SMALLEST_RING_SIZE = utility::graph::smallest_ring_size(atomVD, graph, 5);
	if ( SMALLEST_RING_SIZE >= 3 && SMALLEST_RING_SIZE <= 5 ) {
		SetToType( gasteiger_atom_type_set_->atom_type("O_Te2Te2TeTe"));
	} else {
		bool has_saturated_sulfur = false;
		bool only_carbon = true;
		bool only_silicon = true;
		RealResidueAdjacentIter iter, iter_end;
		for ( boost::tie(iter,iter_end) = boost::adjacent_vertices( atomVD, graph ); iter != iter_end; ++iter ) {
			const ElementCOP & element_data = graph[ *iter ].element_type();
			element::Elements elem( element_data->element() );
			if ( elem == element::H || elem == element::N ||
					elem == element::O || element_data->get_main_group() == 7 /*halogen */ ) {
				SetToType( gasteiger_atom_type_set_->atom_type("O_Te2Te2TeTe"));
				return;
			}
			// special case for sulfur, which needs to be checked for saturation
			if ( elem == element::S && !IsUnsaturated( graph, *iter ) ) {
				has_saturated_sulfur = true;
			}
			if ( elem != element::C ) {
				only_carbon = false;
			}
			if ( elem != element::Si ) {
				only_silicon = false;
			}
		}
		if ( has_saturated_sulfur ) {
			// part of node 2
			SetToType( gasteiger_atom_type_set_->atom_type("O_Te2Te2TeTe"));
		} else if ( only_carbon && ( CountUnsaturatedNeighbors( graph, atomVD ) == 0 ) ) {
			// only bonded to carbon, have to count unsaturated neighbors
			SetToType( gasteiger_atom_type_set_->atom_type("O_Te2Te2TeTe"));
		} else if ( only_silicon ) {
			// only connected to Si; Si - O - Si has strong digonal character (typically > 130 degrees)
			SetToType( gasteiger_atom_type_set_->atom_type("O_DiDiPi2Pi2"));
		} else {
			// other hetero atoms, choose the trigonal type (Node 3)
			// Node 3
			SetToType( gasteiger_atom_type_set_->atom_type("O_Tr2TrTrPi2"));
		}
	}
}

//! @brief choose the final atom type for Oxygen with a single and a double bond
//! @param ATOM the atom of interest
//! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
void PossibleAtomTypesForAtom::FinalizeOxygenSingleDouble( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD ) {
	// use the digonal type only if bound to an atom w/ period 3+ and
	// whose main group (5 for N, 6 for O) is <= the bonded atom, and is at least 4
	// For example, Al should not induce digonality
	RealResidueAdjacentIter iter, iter_end;
	for ( boost::tie(iter,iter_end) = boost::adjacent_vertices( atomVD, graph ); iter != iter_end; ++iter ) {
		const ElementCOP & element_data = graph[ *iter ].element_type();
		if ( element_data->get_period() > 2 && ( element_data->get_main_group() >= 4 && element_data->get_main_group() <= 6) ) {
			SetToType( gasteiger_atom_type_set_->atom_type("O_DiDiPi2Pi"));
			return;
		}
	}
	SetToType( gasteiger_atom_type_set_->atom_type("O_Tr2TrTrPi"));
}

//! @brief choose the final atom type for Oxygen with three single bonds
void PossibleAtomTypesForAtom::FinalizeOxygenThreeSingle( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD ) {
	// Trigonal if bound to anything group 1 - 3, else tetrahedral
	RealResidueAdjacentIter iter, iter_end;
	for ( boost::tie(iter,iter_end) = boost::adjacent_vertices( atomVD, graph ); iter != iter_end; ++iter ) {
		const ElementCOP & element_data = graph[ *iter ].element_type();
		if ( element_data->get_main_group() >= 1 && element_data->get_main_group() <= 3 ) {
			SetToType( gasteiger_atom_type_set_->atom_type("O_TrTrTrPi2"));
			return;
		}
	}
	SetToType( gasteiger_atom_type_set_->atom_type("O_Te2TeTeTe"));
}

////! @brief get connected element types
////! @param ATOM the atom of interest
////! @return a set of the connected element types
//storage::Set< ElementType> PossibleAtomTypesForAtom::GetConnectedElementTypes
//(
//  const AtomConformationalInterface &ATOM
//)
//{
//  storage::Set< ElementType> element_types;
//  for
//  (
//    storage::Vector< BondConformational>::const_iterator
//      itr( ATOM.GetBonds().Begin()), itr_end( ATOM.GetBonds().End());
//    itr != itr_end;
//    ++itr
//  )
//  {
//    element_types.Insert( itr->GetTargetAtom().get_element_type());
//  }
//  return element_types;
//}

//! @brief test whether a particular atom is unsaturated
//! @return true if atom has no A. unsaturated bonds or B. is part of an aromatic ring or C. has empty orbitals
//! @details Virtual atoms are treated like ordinary atoms.
bool PossibleAtomTypesForAtom::IsUnsaturated( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD ) const
{
	RealResidueOutEdgeIter iter, iter_end;
	for ( boost::tie(iter,iter_end) = boost::out_edges(atomVD, graph); iter != iter_end; ++iter ) {
		BondOrder const & order = graph[ *iter ].order();
		if ( order == DoubleBondOrder || order == TripleBondOrder || graph[ *iter ].aromaticity() == IsAromaticBond ) {
			return true;
		}
	}
	// no aromatic or unsaturated bonds; last resort, check for empty orbitals
	const ElementCOP & element_data = graph[ atomVD ].element_type();
	// TODO: Need to normalize this with the formal charge resetting
	if ( element_data->get_main_group() - graph[ atomVD ].formal_charge() != 3 ) {
		return false;
	}
	return true;
}

//! @brief count unsaturated neighbors
//! @param ATOM the atom of interest
//! @return the number of unsaturated neighbors around ATOM
//! @details Virtual atoms are treated like ordinary atoms.
core::Size PossibleAtomTypesForAtom::CountUnsaturatedNeighbors( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD ) const
{
	core::Size unsaturated_counts( 0);
	RealResidueAdjacentIter iter, iter_end;
	for ( boost::tie(iter,iter_end) = boost::adjacent_vertices( atomVD, graph ); iter != iter_end; ++iter ) {
		if ( IsUnsaturated( graph, *iter ) ) {
			++unsaturated_counts;
		}
	}
	return unsaturated_counts;
}

//! @brief test whether atom is bonded to any halogens
//! @param ATOM the atom of interest
//! @return true if the atom is bonded to any halogens
bool PossibleAtomTypesForAtom::IsBondedToAHalogen( const core::chemical::RealResidueGraph & graph, const core::chemical::RealResidueVD & atomVD ) const
{
	RealResidueAdjacentIter iter, iter_end;
	for ( boost::tie(iter,iter_end) = boost::adjacent_vertices( atomVD, graph ); iter != iter_end; ++iter ) {
		const ElementCOP & element_data = graph[ *iter ].element_type();
		if ( element_data->get_main_group() == 7 ) {
			return true;
		}
	}
	return false;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////

void assign_gasteiger_atom_types( core::chemical::ResidueType & restype, bool keep_existing /*= true*/, bool allow_unknown /*= false*/) {
	GasteigerAtomTypeSetCOP typeset( restype.gasteiger_atom_typeset() );
	if ( ! typeset ) {
		typeset = core::chemical::ChemicalManager::get_instance()->gasteiger_atom_type_set();
		assert( typeset );
		restype.set_gasteiger_typeset( typeset );
	}
	assign_gasteiger_atom_types( restype, typeset, keep_existing, allow_unknown );
}

/// @brief Assign gasteiger atom types to the atoms in restype.
/// If "keep_existing" is true keep any already-assigned types
/// Otherwise assign all new types.
///
/// @details  Remaining issues:
///
/// * Input with explicitly mentioned aromatic bonds is not necessarily handled correctly.
/// ** (I haven't extensively tested these yet.)
/// * Oxygen typing in aromatic rings is still an issue.
/// ** While most work, N-O-N bridges like in PDB ligands G98 and 3C3 are missed.
/// ** Conversely, methylated/oxidized dibenzofurans like in KRC and FCM are mistyped as trigonal when they're not.
/// * Rosetta needs to have formal charges notated. (Cannot calculate them from the number of hydrogens present.
///
/// Differences from gasteiger types which aren't necessarily errors:
/// * Phosphooxime oxygens like SEH & CGT are typed trig versus gasteiger tet.
/// * Carbamoyl oxime oxygens like GHR and VPU are typed trig versus gasteiger tet.
/// * 5 member ring oxygens where there's an exocyclic double bond like PHX and FUR are typed trig versus gasteiger tet.
/// * 5 member carbamoyl oxygens like CLW and A91 are typed trig versus gasteiger tet.
/// * The nitroperoxide oxygens of 16X are typed tet versus gasteiger trig
/// * The bridge nitrogens of fused ring structures like PO5, C7M and AZQ are typed trig versus gasteiger tet.

void
assign_gasteiger_atom_types( core::chemical::ResidueType & restype, GasteigerAtomTypeSetCOP gasteiger_atom_type_set, bool keep_existing, bool allow_unknown /* = false */) {
	debug_assert( gasteiger_atom_type_set );
	// This functionality was taken from AtomsCompleteStandardizer

	//# RemoveObviousIonicBonds( ATOMS )
	//# graph::EdgeCoverRingPerception( m_Graph)
	//# AddRingInformationToBondTypes()

	// Initialize();
	std::map< VD, PossibleAtomTypesForAtom > PossibleAtomTypes;

	// A subgraph of the residue type with all the fake bonds and edges removed.
	// The Gasteiger typing assumes what you're typing is all real bonds and atoms.
	RealFilter filter( restype.graph() );
	RealResidueGraph real_graph( restype.graph(), filter );

	//Make sure that the internal graph representation is up to date.
	core::chemical::regenerate_graph_vertex_index( real_graph);

	VIter iter, iter_end;
	for ( boost::tie(iter,iter_end) = restype.atom_iterators(); iter != iter_end; ++iter ) {
		if ( keep_existing && restype.atom(*iter).gasteiger_atom_type() &&
				restype.atom(*iter).gasteiger_atom_type()->is_gasteiger_atom_type() ) {
			; // Ignore vertex
		} else if ( restype.atom(*iter).is_fake() ) {
			restype.atom(*iter).gasteiger_atom_type( gasteiger_atom_type_set->type_for_fake_atoms() );
		} else {
			core::Size connections( restype.n_residue_connections_for_atom( restype.atom_index(*iter) ) );
			// This makes the (probably justifilable) assumption that VDs for the main graph can be converted simply to VDs of the filtered graph.
			PossibleAtomTypes[ *iter ] = GetPossibleTypesForAtom( real_graph, *iter, gasteiger_atom_type_set, connections);
			if ( ! allow_unknown && ! PossibleAtomTypes[ *iter ].GetMostStableType() ) {
				TR.Error << "Cannot find appropriate gasteiger type for atom" << restype.atom_name(*iter) << std::endl;
				sdf::output_residue( std::cerr, restype );
				utility_exit_with_message("Cannot find appropriate gasteiger type for atom.");
			}
		}
	}

	//# SplitRingsByAromaticity();

	// SetAtomTypes();
	std::map< VD, PossibleAtomTypesForAtom >::iterator miter, miter_end;
	for ( miter = PossibleAtomTypes.begin(), miter_end = PossibleAtomTypes.end(); miter != miter_end; ++miter ) {
		if ( miter->second.GetMostStableType() ) {
			// determine final atom type, taking into account neighboring atom types
			miter->second.Finalize( real_graph, miter->first );
			restype.atom(miter->first).gasteiger_atom_type( miter->second.GetMostStableType() );
			// TODO: reset the formal charge if there isn't one already.
		} else if ( allow_unknown ) {
			restype.atom(miter->first).gasteiger_atom_type( 0 ); //
		} else {
			TR.Error << "Cannot find appropriate gasteiger type for atom " << restype.atom_name(miter->first) << std::endl;
			sdf::output_residue( std::cerr, restype );
			utility_exit_with_message("Cannot find appropriate gasteiger type for atom.");
		}
	}

	//# SetConjugationOfBondTypes();
	//# DetermineUnknownBondOrders();

	// Special case for N-terminal lower connection.
	// We assume that connections are going to be primarily used by connection to a C-terminal carbonyl
	// This will turn what is typed separately as an amine (N_Te2TeTeTe) atom into an amide (N_TrTrTrPi2).
	if ( restype.is_polymer() && restype.lower_connect_id() != 0 ) {
		VD lower_connect_vd = restype.lower_connect().vertex();
		assert( lower_connect_vd != ResidueType::null_vertex );
		if ( restype.atom( lower_connect_vd ).gasteiger_atom_type()->get_name() == "N_Te2TeTeTe" ) {
			restype.atom( lower_connect_vd ).gasteiger_atom_type( gasteiger_atom_type_set->atom_type("N_TrTrTrPi2") );
		}
	}
}

//! @brief get all atom types matching a given atom considering its bonds
//! @param ATOM atom of interest
//! @return PossibleAtomTypesForAtom
PossibleAtomTypesForAtom GetPossibleTypesForAtom(
	const core::chemical::RealResidueGraph & graph,
	VD const & atomVD,
	GasteigerAtomTypeSetCOP gasteiger_atom_type_set,
	core::Size connections )
{
	// determine number of sigma orbitals (assume one per bond)
	core::Size nr_bonds( 0 );

	// count conjugated bonds, that is, that had a partial bond order
	core::Size number_aromatic_bonds( 0);
	core::Size number_bonds_unknown_order( 0);

	// count electrons in bonds
	core::Size nr_e_in_bonds( 0);
	// core::Size min_e_in_bonds( 0);
	// core::Size max_e_in_bonds( 0);

	RealResidueOutEdgeIter itr_bonds, iter_end;
	for ( boost::tie(itr_bonds,iter_end) = boost::out_edges(atomVD, graph); itr_bonds != iter_end; ++itr_bonds ) {
		++nr_bonds;
		nr_e_in_bonds += graph[ *itr_bonds ].GetNumberOfElectrons();
		//  min_e_in_bonds += graph[ *itr_bonds ].GetMinimumElectrons();
		//  max_e_in_bonds += graph[ *itr_bonds ].GetMaximumElectrons();
		number_bonds_unknown_order += ! graph[ *itr_bonds ].IsBondOrderKnown();
		number_aromatic_bonds += (graph[ *itr_bonds ].aromaticity() == IsAromaticBond);
	}

	if ( number_bonds_unknown_order == 1 ) {
		// discontinuous aromatic ring, nr_e_in_bonds by 1 to compensate
		++nr_e_in_bonds;
	}

	// Handle connections - each is a (single) bond
	nr_bonds += connections;
	nr_e_in_bonds += 2*connections;

	// divide by 2 to count electrons in orbitals
	nr_e_in_bonds /= 2;
	// min_e_in_bonds /= 2;
	// max_e_in_bonds /= 2;

	// handle implicit hydrogens and charge
	const ElementCOP & element_data = graph[ atomVD ].element_type();
	const ElectronConfiguration &e_config( element_data->get_electron_configuration());
	int formal_charge = graph[ atomVD ].formal_charge();

	int implicit_hydrogens = 0 ;

	//Things like transition metals have undefined sp valence electrons.
	// We'll just assume that they don't have implict hydrogens.
	// (For Rosetta any hydrogens should be explicit anyway.)
	if ( ! utility::is_undefined( e_config.valence_electrons_sp() ) ) {
		const int valence_e( e_config.valence_electrons_sp() - formal_charge);

		// Find the total number of unpaired valence electrons of *itr_atom
		const int unpaired_valence_e(
			nr_e_in_bonds > 4
			? std::max( valence_e, int( e_config.max_valence_electrons_sp()) - valence_e)
			: std::min( valence_e, int( e_config.max_valence_electrons_sp()) - valence_e)
		);
		implicit_hydrogens = std::min( int( 4 - nr_bonds), unpaired_valence_e - int(nr_e_in_bonds));
	}

	// Note: Most implicit hydrogen issues I've seen have been due to misplaced/missing formal charges on atoms,
	// Or atoms simply lacking substituents (e.g. a conjugated residue treated as a stand-alone residue).
	if ( implicit_hydrogens > 0 ) {
		TR.Warning << "Gasteiger atom typing thinks you have implicit Hydrogens!" << std::endl;
		TR.Warning << "   " << core::Size(implicit_hydrogens) << " missing hydrogens on atom " << graph[ atomVD ].name() << std::endl;
		nr_e_in_bonds += implicit_hydrogens;
		nr_bonds += implicit_hydrogens;
	}

	return PossibleAtomTypesForAtom::FindPossibleAtomTypesForAtom(
		gasteiger_atom_type_set,
		*element_data,
		nr_e_in_bonds,
		nr_bonds,
		formal_charge,
		number_aromatic_bonds
	);
} // GetPossibleTypesForAtom


} // gasteiger
} // chemical
} // core
