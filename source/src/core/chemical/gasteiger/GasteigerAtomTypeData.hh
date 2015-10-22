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
// The BCL is copyright Vanderbilt University (Meiler Lab), a RosettaCommons member

/// @file   core/chemical/gasteiger/GasteigerAtomTypeData.hh
/// @brief  The data for the gasteiger atom types
/// @author Rosetta conversion: Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_gasteiger_GasteigerAtomTypeData_hh
#define INCLUDED_core_chemical_gasteiger_GasteigerAtomTypeData_hh

#include <core/chemical/gasteiger/GasteigerAtomTypeData.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/Element.fwd.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <set>
#include <vector>

// Only the Windows PyRosetta build complains about this. But somehow,
// WIN32 is not a useful ifdef here. What.
//#ifdef WIN32
#include <string>
//#else
//#include <iosfwd>
//#endif

namespace core {
namespace chemical {
namespace gasteiger {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!
//! @class AtomTypeData
//! @brief contains hybridization and bond geometry data, which is used in Atom
//!
//! @see @link example_chemistry_atom_type_data.cpp @endlink
//! @author mueller, woetzen, mendenjl
//! @date 08/23/2009
//!
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GasteigerAtomTypeData : public utility::pointer::ReferenceCount
{
	//#      friend class AtomTypes;
	//#      friend class BondLengths;

public:

	///////////
	// Enums //
	///////////

	//! enum properties for atom types
	// see Hinze, Jaffe: Electronegativity. I. Orbital Electronegativity of Neutral Atoms
	// each atom type has up to two different parameter sets for sigma- and pi-orbitals
	enum Properties
	{
		SigmaValenceStateIonizationPotential,    //!< sigma ionization potential
		SigmaValenceStateElectronAffinity,       //!< sigma electron affinity
		SigmaOrbitalElectronegativityMulliken,   //!< sigma orbital electronegativity, Mulliken scale
		SigmaOrbitalElectronegativityPauling,    //!< sigma orbital electronegativity, Pauli scale
		PiValenceStateIonizationPotential,       //!< pi ionization potential
		PiValenceStateElectronAffinity,          //!< pi electron affinity
		PiOrbitalElectronegativityMulliken,      //!< pi orbital electronegativity, Mulliken scale
		PiOrbitalElectronegativityPauling,       //!< pi orbital electronegativity, Pauli scale
		LonePairIonizationPotential,             //!< lone pair ionization potential
		LonePairElectronAffinity,                //!< lone pair electron affinity
		LonePairElectronegativityMulliken,       //!< lone pair electronegativity
		AdditiveAtomicPolarizability,            //!< additive atomic polarizability, see J. Am. Chem. Soc. 1990, 112, 8533-8542
		VdWaalsRadiusCSD,                        //!< Vdw radius for the atom type using data from the CSD; see notes below
		CovalentRadiusSingleBond,                //!< Covalent radius for single bond
		CovalentRadiusDoubleBond,                //!< Covalent radius for double bond
		CovalentRadiusTripleBond,                //!< Covalent radius for triple bond
		CovalentRadiusAromaticBond,              //!< Covalent radius for an aromatic bond
		NumberOfProperties                       //!< Number of properties
	};

	//! @note VdWaalsRadiusCSD is only valid for atoms without hydrogens, and can only be used to
	//! @note detect bad geometries if both atoms have no H (otherwise H bonding may lead to close contacts, also
	//! @note could imply missing bonds), and which do not have opposite charges
	//! @note VdWaalsRadiusCSD was calculated using data from the cambridge structural database
	//! @note Vdw radii may be violated in certain bridged ring systems

	//! @brief element type property as string
	//! @param PROPERTY the property desired
	//! @return the property as string
	static const std::string &GetPropertyName( const Properties PROPERTY);

	//#      //! PropertyEnum simplifies the usage of the Properties enum of this class
	//#      typedef util::WrapperEnum< Properties, &GetPropertyName, NumberOfProperties> PropertyEnum;

	//! how the atom type can contribute to a pi system
	enum PiContributionType
	{
		Zero,       //!< no pi electrons will be contributed; e.g. atom type with no non-singular bonds, no lone pairs, like B_TrTrTr
		One,        //!< exactly one pi electron will be contributed by this atom type (atom type with a double bond)
		Two,        //!< exactly two pi electrons will be contributed by this atom type (atom type with a triple bond or two 2x bonds)
		ZeroOrTwo   //!< Either zero or two pi electrons will be contributed by this atom type (atom type with lone pairs but only singular bonds)
	};

	enum AtomicOrbitalTypes
	{
		S,
		Px,
		Py,
		Pz,
		Dxy,
		Dxz,
		Dyz,
		Dz2,
		Dx2y2,
		NumberOfAtomicOrbitalTypes
	};

	static std::vector<std::string> const & AtomicOrbitalTypes_strings();

	enum HybridOrbitalType  // Set of AtomicOrbitalTypes
	{
		Unhybridized, // empty
		SP,           // S, Px
		SP2,          // S, Px, Py
		SP3,           // S, Px, Py, Pz
		NumberHybridOrbitalType
	}; // Was a more complex data structure: GetNumberOfPossibleSigmaBondingPartners() method is basically the conversion to Size (0,1,2,3)


	//#      //! Base atom type; only contains element type and charge
	//#      util::SiPtr< const AtomType> m_BaseType;

public:

	//////////
	// data //
	//////////

	//#      //! single instance of that class
	//#      static const util::SiPtr< const util::ObjectInterface> s_Instance;

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	//! @brief default constructor
	GasteigerAtomTypeData();

	//! @brief destructor
	~GasteigerAtomTypeData();

	//! @brief the usual constructor
	GasteigerAtomTypeData
	(
		const std::string &NAME,
		const ElementOP &ELEMENT_TYPE,
		const HybridOrbitalType &HYBRIDIZATION,
		const core::Size &HYBRID_ORBITALS_IN_SIGMA_BONDS,
		const core::Size &HYBRID_ORBITALS_NONBINDING,
		const std::set< AtomicOrbitalTypes > &PI_ORBITALS_IN_BONDS,
		const std::set< AtomicOrbitalTypes > &ATOMIC_ORBITALS_NONBINDING,
		const core::Real &SIGMA_VALENCE_STATE_IONIZATION_POTENTIAL,
		const core::Real &SIGMA_VALENCE_STATE_ELECTRON_AFFINITY,
		const core::Real &PI_VALENCE_STATE_IONIZATION_POTENTIAL,
		const core::Real &PI_VALENCE_STATE_ELECTRON_AFFINITY,
		const core::Real &LONE_PAIR_IONIZATION_POTENTIAL,
		const core::Real &LONE_PAIR_ELECTRON_AFFINITY,
		const core::Real &ATOMIC_POLARIZABILITY
	);

	//! @brief the usual constructor, with element set.
	GasteigerAtomTypeData
	(
		const std::string &NAME,
		const std::string &ELEMENT,
		const HybridOrbitalType &HYBRIDIZATION,
		const core::Size &HYBRID_ORBITALS_IN_SIGMA_BONDS,
		const core::Size &HYBRID_ORBITALS_NONBINDING,
		const std::set< AtomicOrbitalTypes > &PI_ORBITALS_IN_BONDS,
		const std::set< AtomicOrbitalTypes > &ATOMIC_ORBITALS_NONBINDING,
		const core::Real &SIGMA_VALENCE_STATE_IONIZATION_POTENTIAL,
		const core::Real &SIGMA_VALENCE_STATE_ELECTRON_AFFINITY,
		const core::Real &PI_VALENCE_STATE_IONIZATION_POTENTIAL,
		const core::Real &PI_VALENCE_STATE_ELECTRON_AFFINITY,
		const core::Real &LONE_PAIR_IONIZATION_POTENTIAL,
		const core::Real &LONE_PAIR_ELECTRON_AFFINITY,
		const core::Real &ATOMIC_POLARIZABILITY,
		ElementSetCOP ele_set
	);

	//! @brief constructor from just an element type and charge
	GasteigerAtomTypeData
	(
		const std::string &NAME,
		ElementOP ELEMENT_TYPE,
		const short &CHARGE
	);


private:

	/// @brief Initialize derived values.
	void initialize();

public:
	//! @brief Clone function
	//! @return pointer to new AtomTypeData
	GasteigerAtomTypeDataOP Clone() const;

	/////////////////
	// data access //
	/////////////////

	//#      //! @brief returns class name of the object behind a pointer or the current object
	//#      //! @return the class name
	//#      const std::string &GetClassIdentifier() const;

	//! return Name
	std::string const & get_name() const;

	//! return ElementType
	ElementCOP get_element_type() const;

	//! @brief returns the hybridization of the atom type
	//! @return the type of hybrid orbital
	HybridOrbitalType get_hybrid_orbital_type() const;

	//! @brief returns the number of hybridized orbitals
	//! @return the number of hybridized orbitals
	core::Size get_number_hybrid_orbitals() const;

	//! @brief returns the number of lone pairs in hybrid orbitals
	//! @return the number of lone pairs in hybrid orbitals
	core::Size get_nNumber_hybrid_lone_pairs() const;

	//! @return Number of hybridized bonds
	core::Size get_number_hybrid_bonds() const;

	//! @return Number of bonds
	core::Size get_number_bonds() const;

	//! @return Number of Sigma orbitals that are not hybridized
	core::Size get_number_unhybridized_sigma_orbitals() const;

	//! @return Number of Sigma orbitals
	core::Size get_number_sigma_orbitals() const;

	//! @return Number of electrons in p orbitals (whether hybridized or not)
	core::Size get_number_electrons_in_p_orbitals() const;

	//! @return Number of pi-orbitals
	core::Size get_number_pi_orbitals() const;

	//! @return Charge
	short get_formal_charge() const;

	//! @return valence electrons in sp orbitals
	core::Size get_valence_electrons_sp() const;

	//! @brief get_number_electrons_in_bonds calculates the total number of electrons in pi-orbital and sigma bonds
	core::Size get_number_electrons_in_bonds() const;

	//! @brief return Number of unhybridized lone pairs
	core::Size get_number_unhybridized_lone_pairs() const;

	//! @brief atom type property as core::Real
	//! @param PROPERTY the property desired
	//! @return the property as core::Real
	core::Real get_atom_type_property( const GasteigerAtomTypeData::Properties PROPERTY) const;

	//! @return the orbital electronegativity associated with the charged state
	core::Real get_orbital_E_neg_pos() const;

	//! @brief determine if this atom type can participate in pi-bond conjugation
	//! @return true iff this atom type has any non-single bonds or lone pairs
	bool is_conjugated() const;

	//! @brief is this a well characterized gasteiger atom type
	//! @return true iff this atom type is this a well characterized gasteiger atom type
	bool is_gasteiger_atom_type() const;

	//#      //! @brief get the base atom type
	//#      //! @return the atom type with only element type and charge information
	//#      const AtomType &GetBaseAtomType() const;

	//! @brief Get the max number of electrons available for contribution to an aromatic ring
	//! @return the max electrons contributed by this atom type to a pi system
	core::Size get_maxE_contribution_to_pi_system() const;

	//! @brief Get the type of contribution this atom type can make to a pi system
	//! @return the type of contribution this atom type can make to a pi system
	PiContributionType get_pi_electron_contribution_type() const;

	//! @brief Get the stability metric.  Electronic stability is indicated by a larger number
	//! This is used to decide between atom types when no other means are possible
	core::Real get_stability_metric() const;

	////////////////
	// operations //
	////////////////

	///////////////
	// operators //
	///////////////

	//////////////////////
	// input and output //
	//////////////////////

public:

	//! @brief read from std::istream
	//! @param ISTREAM input stream
	//! @return istream which was read from
	std::istream &read( std::istream &ISTREAM, ElementSetCAP ele_set);

	//! @brief write to std::ostream
	//! @param OSTREAM output stream to write to
	//! @param INDENT number of indentations
	//! @return output stream which was written to
	std::ostream &write( std::ostream &OSTREAM ) const;

	//////////////////////
	// helper functions //
	//////////////////////

	//! Type difference specifies the difference between two atom types
	enum TypeDifference
	{
		None,                    //!< atom types are identical
		NumberBondingSOrbitals,  //!< different # of electrons in bonding s-orbitals
		NumberBondingPOrbitals,  //!< different # of electrons in bonding p-orbitals
		NumberLonePairOrbitals,  //!< different # of electrons in lone pairs
		Other,                   //!< multiple differences or other differences (e.g. element type),
		NumberTypeDifferences
	};

	//! @brief calculate the stability metric.  Electronic stability is indicated by a larger number
	//! This is used to decide between atom types when no other means are possible
	core::Real calculate_stability_metric() const;

	//#      //! @brief Estimate electronegativities of the anionic, cationic, or ground state given values from at least one other state
	//#      //! @param ELECTRONEGATIVITIES electronegativities for the charge = -1,0, and +1 state
	//#      //! @param TYPE_DIFFERENCE what orbital the electronegativity relates to
	//#      void EstimateUndefinedElectronegativities
	//#      (
	//#        linal::Vector< core::Real> &ELECTRONEGATIVITIES,
	//#        const TypeDifference &TYPE_DIFFERENCE
	//#      );
	//#
	//#      //! @brief set orbital electronegativity for the charged species
	//#      //! @param TYPE_DIFFERENCE whether the electronegativity is S/P Bonding/LonePair orbitals
	//#      //! @param ELECTRONEGATIVITY electronegativity of the charged or neutral species, if available
	//#      //! Used only in AtomTypes constructor
	//#      void SetSimilarOrbitalElectronegativity
	//#      (
	//#        const TypeDifference &TYPE_DIFFERENCE,
	//#        const core::Real &ELECTRONEGATIVITY = numeric::get_undefined_real()
	//#      );
	//#
	//#      //! @brief set orbital electronegativity from the neutral type
	//#      //! @param TYPE_DIFFERENCE whether the electronegativity is S/P Bonding/LonePair orbitals
	//#      //! @param NEUTRAL_TYPE neutral type of the atom, if available
	//#      //! Used only in AtomTypes constructor
	//#      void SetSimilarOrbitalElectronegativity
	//#      (
	//#        const TypeDifference &TYPE_DIFFERENCE,
	//#        const GasteigerAtomTypeData &NEUTRAL_TYPE
	//#      );

	//! @brief get the average ionization potential ratio between cation and neutral atom type that differ by TYPE_DIFFERENCE
	//! @param TYPE_DIFFERENCE the type difference to get the corresponding ratio for
	//! @return the ratio
	core::Real get_average_ip_change_cation_to_neutral( const TypeDifference TYPE_DIFFERENCE) const;

	//! @brief get the average ionization potential ratio between neutral and cation atom type that differ by TYPE_DIFFERENCE
	//! @param TYPE_DIFFERENCE the type difference to get the corresponding ratio for
	//! @return the ratio
	core::Real get_average_ip_change_neutral_to_anion( const TypeDifference TYPE_DIFFERENCE) const;

	//! @brief type difference as string
	//! @param TYPE_DIFFERENCE the type difference for which a string is desired
	//! @return the type difference as a string
	static const std::string &get_type_difference_name( const TypeDifference TYPE_DIFFERENCE);

	//! @brief determine the difference betweent his atom type data and another
	//! @param OTHER the atom type data to compare this atom type data to
	//! @return the corresponding TypeDifference
	TypeDifference difference_from( const GasteigerAtomTypeData &OTHER);

	//! @brief get the electronegativity type corresponding to a TypeDifference
	//! @param TYPE_DIFFERENCE the type difference to get the corresponding electronegativity for
	//! @return the electronegativity type corresponding to TypeDifference
	core::Real get_electronegativity( const TypeDifference TYPE_DIFFERENCE) const;

	//! @brief get the ionization potential type corresponding to a TypeDifference
	//! @param TYPE_DIFFERENCE the type difference to get the corresponding ionization potential for
	//! @return the ionization potential type corresponding to TypeDifference
	core::Real get_ionization_potential( const TypeDifference TYPE_DIFFERENCE) const;

	//! @brief get the electron affinity type corresponding to a TypeDifference
	//! @param TYPE_DIFFERENCE the type difference to get the corresponding electron affinity for
	//! @return the electron affinity type corresponding to TypeDifference
	core::Real get_electron_affinity( const TypeDifference TYPE_DIFFERENCE) const;

	//#      //! @brief get the electronegativity from charge corresponding to a TypeDifference
	//#      //! @param TYPE_DIFFERENCE the type difference to get the corresponding function for
	//#      //! @return the electronegativity as a function of charge polynomial corresponding to TypeDifference
	//#      math::Polynomial &get_electronegativityFromChargeFunction( const TypeDifference &TYPE_DIFFERENCE);

	//! @brief set a particular data
	//! @param DATA the property to set
	//! @param VALUE the value to set the property to
	void set_property( const Properties DATA, const core::Real VALUE);

	//! @brief GetAverageNeutralSigmaIVToEARatio helper function for AtomTypes::CalculateElectronegativityValues
	//! @return reference to a core::Real, which returns the ratio of Average(SigmaValenceStateIonizationPotential) for neutral atoms vs. anions
	static core::Real get_average_neutral_sigma_ip_to_anion_ip_ratio();

	//! @brief get_average_neutral_pi_ip_to_anion_ip_ratio helper function for AtomTypes::CalculateElectronegativityValues
	//! @return reference to a core::Real, which returns the ratio of Average(PiValenceStateIonizationPotential) for neutral atoms vs. anions
	static core::Real get_average_neutral_pi_ip_to_anion_ip_ratio();

	//! @brief get_average_cation_sigma_ip_to_neutral_ip_ratio helper function for AtomTypes::CalculateElectronegativityValues
	//! @return reference to a core::Real, which returns the ratio of Average(SigmaValenceStateIonizationPotential) for cations vs. neutral atoms
	static core::Real get_average_cation_sigma_ip_to_neutral_ip_ratio();

	//! @brief get_average_cation_pi_ip_to_neutral_ip_ratio helper function for AtomTypes::CalculateElectronegativityValues
	//! @return reference to a core::Real, which returns the ratio of Average(PiValenceStateIonizationPotential) for cations vs. neutral atoms
	static core::Real get_average_cation_pi_ip_to_neutral_ip_ratio();

private:

	//////////
	// data //
	//////////

	std::string name_; //! name of the atom type

	ElementCOP element_type_;                 //!< element type

	HybridOrbitalType hybridization_;         //!< Unhybridized, SP, SP2, or SP3

	core::Size number_hybrid_orbitals_sigma_binding_; //!< hybrid orbitals that are binding
	core::Size number_hybrid_orbitals_nonbinding_;   //!< hybrid orbitals that are non binding with their electrons
	core::Size number_electrons_in_bonds_;           //!< Number of e- bonds (calculated)
	core::Size number_bonds_;                      //!< Number of bonds (calculated)
	core::Size max_e_contribution_to_pi_system_;       //!< Number of e- in pi system (calculated)

	std::set< AtomicOrbitalTypes > pi_orbitals_binding_; //!< pi orbitals involved in binding

	//! non hybridized non binding orbitals with their electrons
	std::set< AtomicOrbitalTypes > atomic_orbitals_nonbinding_;

	core::Real properties_[ int( NumberOfProperties)]; //!< real-valued properties
	core::Real orbital_e_neg_pos_; //!< orbital electronegativity associated with the charged state

	//! estimated stability of the atom type; only used to resolve clashes in atom types
	core::Real stability_metric_;

	//! charge of atom type = difference between electrons in bonds and non-bonds minus valence electrons
	short charge_;

	//! whether this atom is conjugated
	bool conjugated_;


	// inflection points are needed for the above charge equations for the rare cases that the charge can end up on
	// the other side of the quadratic polynomial curve, in which case charges run off towards infinity
	// To avoid this problem, if the desired charge falls below the inflection point of the above quadratic or linear
	// functions, we use the value at the inflection point
	core::Real sigma_charge_to_en_inflection_;         //!< Inflection point for sigma charge to inflection
	core::Real pi_charge_to_en_inflection_;            //!< Inflection point for pi charge to inflection
	core::Real sigma_charge_to_lone_pair_en_inflection_; //!< Inflection point for sigma charge to inflection

	bool is_gasteiger_type_; //!< Is the type a proper gasteiger type
	// before SmallMoleculeStandardizer is used on a molecule all atoms in it will not have a proper gasteiger type
	// after SmallMoleculeStandardizer is used, all atoms in it will have have a proper gasteiger atom type, provided
	// one exists

	PiContributionType pi_contribution_; //!< Contribution of this atom type to a pi-conjugated system


}; // class AtomTypeData

inline std::ostream &
operator<< (std::ostream & out, GasteigerAtomTypeData const & obj ){
	return obj.write( out );
}

// No operator>> because we also need to pass a ElementSet in order to parse the element types.


} // gasteiger
} // chemical
} // core

#endif
