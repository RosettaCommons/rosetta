// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
// This file is based off of code from the Biochemistry Library (BCL).
// The BCL is copyright Vanderbilt University (Meiler Lab), a RosettaCommons member

/// @file   core/chemical/gasteiger/GasteigerAtomTypeData.cc
/// @brief  The data for the gasteiger atom types
/// @author Rosetta conversion: Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>

#include <core/chemical/ElementSet.hh>
#include <core/chemical/Element.hh>

#include <core/chemical/gasteiger/util.hh>

#include <utility>
#include <utility/numbers.hh>
#include <utility/exit.hh>
#include <utility/tools/make_vector.hh>
#include <utility/thread/backwards_thread_local.hh>

namespace core {
namespace chemical {
namespace gasteiger {

using utility::tools::make_vector;

std::vector<std::string> const & GasteigerAtomTypeData::AtomicOrbitalTypes_strings()
{
	static const std::vector<std::string> AtomicOrbitalTypes_strings_ =
		make_vector<std::string> (
		"S",
		"Px",
		"Py",
		"Pz",
		"Dxy",
		"Dxz",
		"Dyz",
		"Dz2",
		"Dx2y2",
		"AtomicOrbitalTypes");

	return AtomicOrbitalTypes_strings_;
}


//#    // instantiate s_Instance
//#    const util::SiPtr< const util::ObjectInterface> GasteigerAtomTypeData::s_Instance
//#    (
//#      GetObjectInstances().AddInstance( new GasteigerAtomTypeData())
//#    );

///////////
// Enums //
///////////

//! @brief element type property as string
//! @param PROPERTY the property desired
//! @return the property as string
const std::string &GasteigerAtomTypeData::GetPropertyName( const GasteigerAtomTypeData::Properties PROPERTY)
{
	static const std::string s_properties[] =
		{
		"SigmaValenceStateIonizationPotential",
		"SigmaValenceStateElectronAffinity",
		"SigmaOrbitalElectronegativityMulliken",
		"SigmaOrbitalElectronegativityPauling",
		"PiValenceStateIonizationPotential",
		"PiValenceStateElectronAffinity",
		"PiOrbitalElectronegativityMulliken",
		"PiOrbitalElectronegativityPauling",
		"LonePairIonizationPotential",
		"LonePairElectronAffinity",
		"LonePairElectronegativity",
		"AdditiveAtomicPolarizability",
		"VdWaalsRadiusCSD",
		"CovalentRadiusSingleBond",
		"CovalentRadiusDoubleBond",
		"CovalentRadiusTripleBond",
		"CovalentRadiusAromaticBond",
		"GasteigerAtomTypeDataProperties" //GetStaticClassName< Properties>()
		};

	return s_properties[ PROPERTY];
}

//////////////////////////////////
// construction and destruction //
//////////////////////////////////

//! @brief default constructor
GasteigerAtomTypeData::GasteigerAtomTypeData() :
	name_("UNKNOW_GASTEIGER_ATOM_TYPE"),
	hybridization_(Unhybridized), // amw cppcheck initialize this!
	number_hybrid_orbitals_sigma_binding_( 0),
	number_hybrid_orbitals_nonbinding_( 0),
	number_electrons_in_bonds_( 0),
	number_bonds_( 0),
	max_e_contribution_to_pi_system_( 0),
	orbital_e_neg_pos_( 0),
	stability_metric_( 0),
	charge_( 0),
	conjugated_( false),
	sigma_charge_to_en_inflection_( -std::numeric_limits< core::Real>::infinity()),
	pi_charge_to_en_inflection_( -std::numeric_limits< core::Real>::infinity()),
	sigma_charge_to_lone_pair_en_inflection_( -std::numeric_limits< core::Real>::infinity()),
	is_gasteiger_type_( false),
	pi_contribution_( Zero)
{
	// make all the properties undefined
	for ( double & propertie : properties_ ) {
		propertie = utility::get_undefined_real();
	}
}

//! @brief destructor
GasteigerAtomTypeData::~GasteigerAtomTypeData() = default;

//! @brief the usual constructor
GasteigerAtomTypeData::GasteigerAtomTypeData
(
	std::string const & NAME,
	const ElementOP &ELEMENT_TYPE,
	const HybridOrbitalType &HYBRIDIZATION,
	const core::Size &HYBRID_ORBITALS_IN_SIGMA_BONDS,
	const core::Size &HYBRID_ORBITALS_NONBINDING,
	std::set< AtomicOrbitalTypes > const & PI_ORBITALS_IN_BONDS,
	std::set< AtomicOrbitalTypes > const & ATOMIC_ORBITALS_NONBINDING,
	const core::Real &SIGMA_VALENCE_STATE_IONIZATION_POTENTIAL,
	const core::Real &SIGMA_VALENCE_STATE_ELECTRON_AFFINITY,
	const core::Real &PI_VALENCE_STATE_IONIZATION_POTENTIAL,
	const core::Real &PI_VALENCE_STATE_ELECTRON_AFFINITY,
	const core::Real &LONE_PAIR_IONIZATION_POTENTIAL,
	const core::Real &LONE_PAIR_ELECTRON_AFFINITY,
	const core::Real &ATOMIC_POLARIZABILITY
) :
	name_( NAME ),
	element_type_( ELEMENT_TYPE),
	hybridization_( HYBRIDIZATION),
	number_hybrid_orbitals_sigma_binding_( HYBRID_ORBITALS_IN_SIGMA_BONDS),
	number_hybrid_orbitals_nonbinding_( HYBRID_ORBITALS_NONBINDING),
	pi_orbitals_binding_( PI_ORBITALS_IN_BONDS),
	atomic_orbitals_nonbinding_( ATOMIC_ORBITALS_NONBINDING)
{

	properties_[ SigmaValenceStateIonizationPotential] = SIGMA_VALENCE_STATE_IONIZATION_POTENTIAL;
	properties_[ SigmaValenceStateElectronAffinity] = SIGMA_VALENCE_STATE_ELECTRON_AFFINITY;
	properties_[ PiValenceStateIonizationPotential] = PI_VALENCE_STATE_IONIZATION_POTENTIAL;
	properties_[ PiValenceStateElectronAffinity] = PI_VALENCE_STATE_ELECTRON_AFFINITY;
	properties_[ LonePairIonizationPotential] = LONE_PAIR_IONIZATION_POTENTIAL;
	properties_[ LonePairElectronAffinity] = LONE_PAIR_ELECTRON_AFFINITY;
	properties_[ AdditiveAtomicPolarizability] = ATOMIC_POLARIZABILITY;

	initialize();
}

//! @brief the usual constructor
GasteigerAtomTypeData::GasteigerAtomTypeData
(
	std::string const & NAME,
	const std::string &ELEMENT,
	const HybridOrbitalType &HYBRIDIZATION,
	const core::Size &HYBRID_ORBITALS_IN_SIGMA_BONDS,
	const core::Size &HYBRID_ORBITALS_NONBINDING,
	std::set< AtomicOrbitalTypes > const & PI_ORBITALS_IN_BONDS,
	std::set< AtomicOrbitalTypes > const & ATOMIC_ORBITALS_NONBINDING,
	const core::Real &SIGMA_VALENCE_STATE_IONIZATION_POTENTIAL,
	const core::Real &SIGMA_VALENCE_STATE_ELECTRON_AFFINITY,
	const core::Real &PI_VALENCE_STATE_IONIZATION_POTENTIAL,
	const core::Real &PI_VALENCE_STATE_ELECTRON_AFFINITY,
	const core::Real &LONE_PAIR_IONIZATION_POTENTIAL,
	const core::Real &LONE_PAIR_ELECTRON_AFFINITY,
	const core::Real &ATOMIC_POLARIZABILITY,
	ElementSetCOP ele_set
) :
	name_( NAME ),
	hybridization_( HYBRIDIZATION),
	number_hybrid_orbitals_sigma_binding_( HYBRID_ORBITALS_IN_SIGMA_BONDS),
	number_hybrid_orbitals_nonbinding_( HYBRID_ORBITALS_NONBINDING),
	pi_orbitals_binding_( PI_ORBITALS_IN_BONDS),
	atomic_orbitals_nonbinding_( ATOMIC_ORBITALS_NONBINDING)
{
	element_type_ = ele_set->element( ELEMENT );
	properties_[ SigmaValenceStateIonizationPotential] = SIGMA_VALENCE_STATE_IONIZATION_POTENTIAL;
	properties_[ SigmaValenceStateElectronAffinity] = SIGMA_VALENCE_STATE_ELECTRON_AFFINITY;
	properties_[ PiValenceStateIonizationPotential] = PI_VALENCE_STATE_IONIZATION_POTENTIAL;
	properties_[ PiValenceStateElectronAffinity] = PI_VALENCE_STATE_ELECTRON_AFFINITY;
	properties_[ LonePairIonizationPotential] = LONE_PAIR_IONIZATION_POTENTIAL;
	properties_[ LonePairElectronAffinity] = LONE_PAIR_ELECTRON_AFFINITY;
	properties_[ AdditiveAtomicPolarizability] = ATOMIC_POLARIZABILITY;

	initialize();
}


//! @brief constructor from just an element type and charge
GasteigerAtomTypeData::GasteigerAtomTypeData
(
	std::string const & NAME,
	const ElementOP ELEMENT_TYPE,
	const short CHARGE
) :
	name_( NAME ),
	element_type_( ELEMENT_TYPE),
	number_hybrid_orbitals_sigma_binding_( utility::get_undefined_size()),
	number_hybrid_orbitals_nonbinding_( utility::get_undefined_size()),
	number_electrons_in_bonds_( 0),
	number_bonds_( 0),
	max_e_contribution_to_pi_system_( 0),
	orbital_e_neg_pos_( 0),
	stability_metric_( 0),
	charge_( CHARGE),
	conjugated_( false),
	sigma_charge_to_en_inflection_( -std::numeric_limits< core::Real >::infinity()),
	pi_charge_to_en_inflection_( -std::numeric_limits< core::Real >::infinity()),
	sigma_charge_to_lone_pair_en_inflection_( -std::numeric_limits< core::Real >::infinity()),
	is_gasteiger_type_( false),
	pi_contribution_( Zero)
{
	// make all the properties undefined
	for ( double & propertie : properties_ ) {
		propertie = utility::get_undefined_real();
	}
}

void GasteigerAtomTypeData::initialize() {
	orbital_e_neg_pos_ = 0;
	stability_metric_ = 0;
	charge_ = element_type_->get_electron_configuration().valence_electrons_sp()
		- number_hybrid_orbitals_sigma_binding_ - pi_orbitals_binding_.size()
		- 2 * ( number_hybrid_orbitals_nonbinding_ + atomic_orbitals_nonbinding_.size());
	sigma_charge_to_en_inflection_ = -std::numeric_limits< core::Real>::infinity();
	pi_charge_to_en_inflection_ = -std::numeric_limits< core::Real>::infinity();
	sigma_charge_to_lone_pair_en_inflection_ = -std::numeric_limits< core::Real>::infinity();
	is_gasteiger_type_ = true ;

	// make all remaining properties undefined
	for ( int property_number( VdWaalsRadiusCSD); property_number < NumberOfProperties; ++property_number ) {
		properties_[ property_number] = utility::get_undefined_real();
	}

	// calculate and store the mulliken electronegativities, (ElectronAffinity+IonizationPotential)/2
	properties_[ PiOrbitalElectronegativityMulliken]
		= 0.5 * ( properties_[ PiValenceStateElectronAffinity] + properties_[ PiValenceStateIonizationPotential]);
	properties_[ SigmaOrbitalElectronegativityMulliken]
		= 0.5 * ( properties_[ SigmaValenceStateElectronAffinity] + properties_[ SigmaValenceStateIonizationPotential] );
	properties_[ LonePairElectronegativityMulliken]
		= 0.5 * ( 3.0 * properties_[ LonePairElectronAffinity] - properties_[ LonePairIonizationPotential] );

	// convert from the mulliken electronegativity scale to pauling, e.g. see Hinze, Jaffe, 1963
	properties_[ SigmaOrbitalElectronegativityPauling]
		= 0.336 * ( properties_[ SigmaOrbitalElectronegativityMulliken] - 0.615);
	properties_[ PiOrbitalElectronegativityPauling]
		= 0.336 * ( properties_[ PiOrbitalElectronegativityMulliken] - 0.615);

	number_electrons_in_bonds_ = number_hybrid_orbitals_sigma_binding_ + pi_orbitals_binding_.size();
	number_bonds_ = get_number_hybrid_bonds() == 0 ? number_electrons_in_bonds_ : get_number_hybrid_bonds();
	if ( get_number_pi_orbitals() > get_number_unhybridized_sigma_orbitals() ) {
		// case where there is at least one core::Real or triple bond
		// then the # contributed is exactly = this expression since the electrons are directly on the bonds
		// in the ring
		max_e_contribution_to_pi_system_ =
			std::min( core::Size( 2), get_number_pi_orbitals() - get_number_unhybridized_sigma_orbitals());
	} else {
		// all singular bonds; contribution to pi system is 2 if there are any lone pairs, 0 otherwise
		max_e_contribution_to_pi_system_ =
			2 * std::min( get_number_unhybridized_lone_pairs() + get_nNumber_hybrid_lone_pairs(), core::Size( 1));
	}
	stability_metric_ = calculate_stability_metric();

	conjugated_ =
		get_number_pi_orbitals() > get_number_unhybridized_sigma_orbitals()
		|| !atomic_orbitals_nonbinding_.empty()
		|| number_hybrid_orbitals_nonbinding_
		|| hybridization_ == SP2;
	pi_contribution_ = PiContributionType( max_e_contribution_to_pi_system_);
	if ( max_e_contribution_to_pi_system_ == core::Size( 2) && number_electrons_in_bonds_ == number_bonds_ ) {
		// lone pair, no pi electrons in bonds, so the lone pair electrons are either all in or all out of the aromatic
		// system
		pi_contribution_ = ZeroOrTwo;
	}
}

//! @brief Clone function
//! @return pointer to new GasteigerAtomTypeData
GasteigerAtomTypeDataOP GasteigerAtomTypeData::Clone() const
{
	return GasteigerAtomTypeDataOP( new GasteigerAtomTypeData( *this) );
}

/////////////////
// data access //
/////////////////

//#    //! @brief returns class name of the object behind a pointer or the current object
//#    //! @return the class name
//#    const std::string &GasteigerAtomTypeData::GetClassIdentifier() const
//#    {
//#      return GetStaticClassName( *this);
//#    }

//! return Name
std::string const & GasteigerAtomTypeData::get_name() const {
	return name_;
}

//! return ElementType
ElementCOP GasteigerAtomTypeData::get_element_type() const
{
	return element_type_;
}

//! @brief returns the hybridization of the atom type
//! @return the type of hybrid orbital
GasteigerAtomTypeData::HybridOrbitalType GasteigerAtomTypeData::get_hybrid_orbital_type() const
{
	return hybridization_;
}

//! @brief returns the number of hybridized orbitals
//! @return the number of hybridized orbitals
core::Size GasteigerAtomTypeData::get_number_hybrid_orbitals() const
{
	return number_hybrid_orbitals_sigma_binding_ + number_hybrid_orbitals_nonbinding_;
}

//! @brief returns the number of lone pairs in hybrid orbitals
//! @return the number of lone pairs in hybrid orbitals
core::Size GasteigerAtomTypeData::get_nNumber_hybrid_lone_pairs() const
{
	return number_hybrid_orbitals_nonbinding_;
}

//! return Number of bonds
core::Size GasteigerAtomTypeData::get_number_hybrid_bonds() const
{
	return number_hybrid_orbitals_sigma_binding_ == 0
		? pi_orbitals_binding_.size()
		: number_hybrid_orbitals_sigma_binding_;
}

//! return Number of Sigma orbitals that are not hybridized
core::Size GasteigerAtomTypeData::get_number_unhybridized_sigma_orbitals() const
{
	return core::Size( pi_orbitals_binding_.find(S) != pi_orbitals_binding_.end() );
}

//! return Number of Sigma orbitals
core::Size GasteigerAtomTypeData::get_number_sigma_orbitals() const
{
	// start with the number of pi orbitals binding
	return get_number_unhybridized_sigma_orbitals() + number_hybrid_orbitals_sigma_binding_;
}

//! @brief return Number of unhybridized lone pairs
core::Size GasteigerAtomTypeData::get_number_unhybridized_lone_pairs() const
{
	return atomic_orbitals_nonbinding_.size();
}

//! return Number of occupied p orbitals
core::Size GasteigerAtomTypeData::get_number_electrons_in_p_orbitals() const
{
	return number_electrons_in_bonds_ - get_number_unhybridized_sigma_orbitals();
}

//! return Number of pi-orbitals
core::Size GasteigerAtomTypeData::get_number_pi_orbitals() const
{
	return pi_orbitals_binding_.size();
}

//! return Charge
short GasteigerAtomTypeData::get_formal_charge() const
{
	return charge_;
}

//! @brief get_number_electrons_in_bonds calculates the total number of electrons in pi-orbital and sigma bonds
core::Size GasteigerAtomTypeData::get_number_electrons_in_bonds() const
{
	return number_electrons_in_bonds_;
}

//! @brief atom type property as core::Real
//! @param PROPERTY the property desired
//! @return the property as core::Real
core::Real GasteigerAtomTypeData::get_atom_type_property( const GasteigerAtomTypeData::Properties PROPERTY) const
{
	return properties_[ PROPERTY];
}


//! @return the orbital electronegativity associated with the charged state
core::Real GasteigerAtomTypeData::get_orbital_E_neg_pos() const
{
	return orbital_e_neg_pos_;
}

//! @return valence electrons in sp orbitals
core::Size GasteigerAtomTypeData::get_valence_electrons_sp() const
{
	return element_type_->get_electron_configuration().valence_electrons_sp() - get_formal_charge();
}

//! @return Number of bonds
core::Size GasteigerAtomTypeData::get_number_bonds() const
{
	return number_bonds_;
}

//! @brief determine if this atom type can participate in pi-bond conjugation
//! @return true iff this atom type has any e- in p-orbitals or lone pairs
bool GasteigerAtomTypeData::is_conjugated() const
{
	return conjugated_;
}

//! @brief is this a well characterized gasteiger atom type
//! @return true iff this atom type is this a well characterized gasteiger atom type
bool GasteigerAtomTypeData::is_gasteiger_atom_type() const
{
	return is_gasteiger_type_;
}

//#    //! @brief get the base atom type
//#    //! @return the atom type with only element type and charge information
//#    const AtomType &GasteigerAtomTypeData::GetBaseAtomType() const
//#    {
//#      return *m_BaseType;
//#    }

//! @brief Get the max number of electrons available for contribution to an aromatic ring
//! @return the max electrons contributed by this atom type to a pi system
core::Size GasteigerAtomTypeData::get_maxE_contribution_to_pi_system() const
{
	return max_e_contribution_to_pi_system_;
}

//! @brief Get the type of contribution this atom type can make to a pi system
//! @return the type of contribution this atom type can make to a pi system
GasteigerAtomTypeData::PiContributionType GasteigerAtomTypeData::get_pi_electron_contribution_type() const
{
	return pi_contribution_;
}

//! @brief Get the stability metric.  Electronic stability is indicated by a larger number
//! This is used to decide between atom types when no other means are possible
core::Real GasteigerAtomTypeData::get_stability_metric() const
{
	return stability_metric_;
}

////////////////
// operations //
////////////////

///////////////
// operators //
///////////////

//////////////////////
// input and output //
//////////////////////

//! @brief read from std::istream
//! @param ISTREAM input stream
//! @return istream which was read from
std::istream &GasteigerAtomTypeData::read( std::istream &ISTREAM, ElementSetCAP ele_set_ap)
{
	ElementSetCOP ele_set( ele_set_ap );
	// runtime_assert( ele_set ); // ^ above will throw an exception is ele_set_ap is invalid or null
	// write member
	ISTREAM >> name_;
	std::string symbol;
	ISTREAM >> symbol;
	if ( ! ele_set->contains_element_type( symbol ) ) {
		utility_exit_with_message("Unrecognized element symbol '" + symbol + "' reading Gasteiger atom type data.");
	}
	element_type_ = ele_set->element( symbol );
	int hyb;
	ISTREAM >> hyb;
	hybridization_ = HybridOrbitalType(hyb);
	safe_read( ISTREAM, number_hybrid_orbitals_sigma_binding_ );
	safe_read( ISTREAM, number_hybrid_orbitals_nonbinding_ );
	std::string enumset;
	ISTREAM >> enumset;
	if ( enumset[0] != 'p' || enumset[1] != 'i' ) utility_exit_with_message( "Unrecognized compact enum set representation in gasteiger atom type file" );
	pi_orbitals_binding_ = parse_enum_set< AtomicOrbitalTypes >( enumset );
	ISTREAM >> enumset;
	if ( enumset[0] != 'n' || enumset[1] != 'b' ) utility_exit_with_message( "Unrecognized compact enum set representation in gasteiger atom type file" );
	atomic_orbitals_nonbinding_ = parse_enum_set< AtomicOrbitalTypes >( enumset );

	safe_read(ISTREAM, properties_[SigmaValenceStateIonizationPotential]);
	safe_read(ISTREAM, properties_[SigmaValenceStateElectronAffinity]);
	safe_read(ISTREAM, properties_[PiValenceStateIonizationPotential]);
	safe_read(ISTREAM, properties_[PiValenceStateElectronAffinity]);
	safe_read(ISTREAM, properties_[LonePairIonizationPotential]);
	safe_read(ISTREAM, properties_[LonePairElectronAffinity]);
	safe_read(ISTREAM, properties_[AdditiveAtomicPolarizability]);

	initialize();

	return ISTREAM;
}

//! @brief write to std::ostream
//! @param OSTREAM output stream to write to
//! @param INDENT number of indentations
//! @return output stream which was written to
std::ostream &GasteigerAtomTypeData::write( std::ostream &OSTREAM) const
{
	// write member
	OSTREAM << name_ << ' ';
	OSTREAM << element_type_->get_chemical_symbol() << ' ';
	OSTREAM << hybridization_ << ' ';
	safe_write( OSTREAM, number_hybrid_orbitals_sigma_binding_ );
	safe_write( OSTREAM, number_hybrid_orbitals_nonbinding_ );

	OSTREAM << compact_enum_set( pi_orbitals_binding_, "pi") << ' ';
	OSTREAM << compact_enum_set( atomic_orbitals_nonbinding_, "nb") << ' ';

	safe_write(OSTREAM, properties_[SigmaValenceStateIonizationPotential]);
	safe_write(OSTREAM, properties_[SigmaValenceStateElectronAffinity]);
	safe_write(OSTREAM, properties_[PiValenceStateIonizationPotential]);
	safe_write(OSTREAM, properties_[PiValenceStateElectronAffinity]);
	safe_write(OSTREAM, properties_[LonePairIonizationPotential]);
	safe_write(OSTREAM, properties_[LonePairElectronAffinity]);
	safe_write(OSTREAM, properties_[AdditiveAtomicPolarizability]);

	OSTREAM << std::endl;
	return OSTREAM;
}

//////////////////////
// helper functions //
//////////////////////

//! @brief calculate the stability metric.  Electronic stability is indicated by a larger number
//! This is used to decide between atom types when no other means are possible
core::Real GasteigerAtomTypeData::calculate_stability_metric() const
{
	// Consider the markov chain model with this atom type
	// Allowed transitions are
	// 1. Gaining an electron in a sigma orbital
	//    This releases SigmaValenceStateElectronAffinity eV
	// 2. Removing an electron from a sigma orbital
	//    This requires SigmaValenceStateIonizationPotential eV
	// 3. Gaining an electron in a pi orbital
	//    This releases PiValenceStateElectronAffinity eV
	// 4. Removing an electron from a pi orbital
	//    This releases PiValenceStateIonizationPotential eV
	// Given this markov chain model, the effective time constant is
	// 1 / TimeConstant
	// = 1 / ( 1/SigmaValenceStateElectronAffinity + 1/SigmaValenceStateIonizationPotential
	//         + 1/PiValenceStateElectronAffinity + 1/PiValenceStateIonizationPotential)
	// The atom type with the longest time constant should be the most stable, provided that they both have the same
	// charge

	core::Real sigma_ea( properties_[ GasteigerAtomTypeData::SigmaValenceStateElectronAffinity]);
	core::Real sigma_ip( properties_[ GasteigerAtomTypeData::SigmaValenceStateIonizationPotential]);
	core::Real pi_ea( properties_[ GasteigerAtomTypeData::PiValenceStateElectronAffinity]);
	core::Real pi_ip( properties_[ GasteigerAtomTypeData::PiValenceStateIonizationPotential]);

	core::Real reciprocal_time_constant( 0.0);
	if ( sigma_ea > 0.0 && !utility::is_undefined(sigma_ea) ) { // add data on the energy released by adding 1 sigma-electron
		reciprocal_time_constant += 1.0 / sigma_ea;
	}
	if ( sigma_ip > 0.0 && !utility::is_undefined(sigma_ip) ) { // add data on the energy required to remove a sigma-electron
		reciprocal_time_constant += 1.0 / sigma_ip;
	}
	if ( pi_ea > 0.0 && !utility::is_undefined(pi_ea) ) { // add data on the energy released by adding 1 pi-electron
		reciprocal_time_constant += 1.0 / pi_ea;
	}
	if ( pi_ip > 0.0 && !utility::is_undefined(pi_ip) ) { // add data on the energy required to remove a pi-electron
		reciprocal_time_constant += 1.0 / pi_ip;
	}

	// if the time constant is 0.0 or less, then the type is very unstable, so we should return 0 for the stability
	return ( reciprocal_time_constant <= 0.0 ? 0.0 : 1.0 / reciprocal_time_constant);
}

//! @brief GetAverageNeutralSigmaIVToEARatio helper function for AtomTypes::CalculateElectronegativityValues
//! @return reference to a core::Real, which returns the ratio of Average(SigmaValenceStateIonizationPotential) for neutral atoms vs. anions
core::Real GasteigerAtomTypeData::get_average_neutral_sigma_ip_to_anion_ip_ratio()
{
	static const core::Real s_Data( utility::get_undefined_real() );
	return s_Data;
}

//! @brief get_average_neutral_pi_ip_to_anion_ip_ratio helper function for AtomTypes::CalculateElectronegativityValues
//! @return reference to a core::Real, which returns the ratio of Average(PiValenceStateIonizationPotential) for neutral atoms vs. anions
core::Real GasteigerAtomTypeData::get_average_neutral_pi_ip_to_anion_ip_ratio()
{
	static const core::Real s_Data( utility::get_undefined_real() );
	return s_Data;
}

//! @brief get_average_cation_sigma_ip_to_neutral_ip_ratio helper function for AtomTypes::CalculateElectronegativityValues
//! @return reference to a core::Real, which returns the ratio of Average(SigmaValenceStateIonizationPotential) for cations vs. neutral atoms
core::Real GasteigerAtomTypeData::get_average_cation_sigma_ip_to_neutral_ip_ratio()
{
	static const core::Real s_Data( utility::get_undefined_real() );
	return s_Data;
}

//! @brief get_average_cation_pi_ip_to_neutral_ip_ratio helper function for AtomTypes::CalculateElectronegativityValues
//! @return reference to a core::Real, which returns the ratio of Average(PiValenceStateIonizationPotential) for cations vs. neutral atoms
core::Real GasteigerAtomTypeData::get_average_cation_pi_ip_to_neutral_ip_ratio()
{
	static const core::Real s_Data( utility::get_undefined_real() );
	return s_Data;
}
//! @brief type difference as string
//! @param TYPE_DIFFERENCE the type difference for which a string is desired
//! @return the type difference as a string
const std::string &GasteigerAtomTypeData::get_type_difference_name( const GasteigerAtomTypeData::TypeDifference TYPE_DIFFERENCE)
{
	static const std::string s_Properties[] =
		{
		"None",
		"NumberBondingSOrbitals",
		"NumberBondingPOrbitals",
		"NumberLonePairOrbitals",
		"Other",
		"GasteigerAtomTypeData_TypeDifference" //GetStaticClassName< TypeDifference>()
		};

	return s_Properties[ TYPE_DIFFERENCE];
}

//! @brief get the electronegativity type corresponding to TypeDifference
//! @param TYPE_DIFFERENCE the type difference to get the corresponding electronegativity for
//! @return the electronegativity type corresponding to TypeDifference
core::Real GasteigerAtomTypeData::get_electronegativity( const GasteigerAtomTypeData::TypeDifference TYPE_DIFFERENCE) const
{
	switch( TYPE_DIFFERENCE)
			{
			case NumberBondingPOrbitals : return properties_[ PiOrbitalElectronegativityMulliken];
			case NumberBondingSOrbitals : return properties_[ SigmaOrbitalElectronegativityMulliken];
			case NumberLonePairOrbitals : return properties_[ LonePairElectronegativityMulliken];
			default :
				return utility::get_undefined_real();
			}
	return utility::get_undefined_real();
}

//! @brief set a particular data
//! @param DATA the property to set
//! @param VALUE the value to set the property to
void GasteigerAtomTypeData::set_property( const Properties DATA, const core::Real VALUE)
{
	properties_[ DATA] = VALUE;
}

//! @brief get the average ionization potential ratio between cation and neutral atom type that differ by TYPE_DIFFERENCE
//! @param TYPE_DIFFERENCE the type difference to get the corresponding ratio for
//! @return the ratio
core::Real GasteigerAtomTypeData::get_average_ip_change_cation_to_neutral( const TypeDifference TYPE_DIFFERENCE) const
{
	switch( TYPE_DIFFERENCE)
			{
			case NumberBondingPOrbitals : return get_average_cation_pi_ip_to_neutral_ip_ratio();
			case NumberBondingSOrbitals : return get_average_cation_sigma_ip_to_neutral_ip_ratio();
			case NumberLonePairOrbitals :
				return 0.5 * ( get_average_cation_pi_ip_to_neutral_ip_ratio() + get_average_cation_sigma_ip_to_neutral_ip_ratio());
			default :
				return utility::get_undefined_real();
			}
	return utility::get_undefined_real();
}

//! @brief get the average ionization potential ratio between neutral and cation atom type that differ by TYPE_DIFFERENCE
//! @param TYPE_DIFFERENCE the type difference to get the corresponding ratio for
//! @return the ratio
core::Real GasteigerAtomTypeData::get_average_ip_change_neutral_to_anion( const TypeDifference TYPE_DIFFERENCE) const
{
	switch( TYPE_DIFFERENCE)
			{
			case NumberBondingPOrbitals : return get_average_neutral_pi_ip_to_anion_ip_ratio();
			case NumberBondingSOrbitals : return get_average_neutral_sigma_ip_to_anion_ip_ratio();
			case NumberLonePairOrbitals :
				return 0.5 * ( get_average_neutral_pi_ip_to_anion_ip_ratio() + get_average_neutral_sigma_ip_to_anion_ip_ratio());
			default :
				return utility::get_undefined_real();
			}
	return utility::get_undefined_real();
}

//! @brief get the ionization potential type corresponding to TypeDifference
//! @param TYPE_DIFFERENCE the type difference to get the corresponding ionization potential for
//! @return the ionization potential type corresponding to TypeDifference
core::Real GasteigerAtomTypeData::get_ionization_potential( const TypeDifference TYPE_DIFFERENCE) const
{
	switch( TYPE_DIFFERENCE)
			{
			case NumberBondingPOrbitals : return properties_[ PiValenceStateIonizationPotential];
			case NumberBondingSOrbitals : return properties_[ SigmaValenceStateIonizationPotential];
			case NumberLonePairOrbitals : return properties_[ LonePairIonizationPotential];
			default :
				return utility::get_undefined_real();
			}
	return utility::get_undefined_real();
}

//! @brief get the electron affinity type corresponding to TypeDifference
//! @param TYPE_DIFFERENCE the type difference to get the corresponding electron affinity for
//! @return the electron affinity type corresponding to TypeDifference
core::Real GasteigerAtomTypeData::get_electron_affinity( const TypeDifference TYPE_DIFFERENCE) const
{
	switch( TYPE_DIFFERENCE)
			{
			case NumberBondingPOrbitals : return properties_[ PiValenceStateElectronAffinity];
			case NumberBondingSOrbitals : return properties_[ SigmaValenceStateElectronAffinity];
			case NumberLonePairOrbitals : return properties_[ LonePairElectronAffinity];
			default :
				return utility::get_undefined_real();
			}
	return utility::get_undefined_real();
}

//#    //! @brief get the electronegativity from charge corresponding to a TypeDifference
//#    //! @param TYPE_DIFFERENCE the type difference to get the corresponding function for
//#    //! @return the electronegativity as a function of charge polynomial corresponding to TypeDifference
//#    math::Polynomial &GasteigerAtomTypeData::get_electronegativityFromChargeFunction( const TypeDifference &TYPE_DIFFERENCE)
//#    {
//#      static math::Polynomial s_undefined;
//#      switch( TYPE_DIFFERENCE)
//#      {
//#        case NumberBondingPOrbitals: return m_PiChargeToEN;
//#        case NumberBondingSOrbitals: return m_SigmaChargeToEN;
//#        case NumberLonePairOrbitals: return m_SigmaChargeToLonePairEN;
//#        default:
//#          return s_undefined;
//#      }
//#      return s_undefined;
//#    }

//! @brief determine the difference betweent his atom type data and another
//! @param OTHER the atom type data to compare this atom type data to
//! @return the corresponding TypeDifference
GasteigerAtomTypeData::TypeDifference GasteigerAtomTypeData::difference_from( const GasteigerAtomTypeData &OTHER)
{
	if ( &OTHER == this ) {
		return None;
	}

	if ( charge_ == OTHER.charge_ || element_type_ != OTHER.element_type_ || hybridization_ != OTHER.hybridization_ ) {
		return Other;
	}

	const core::Size number_sigma_orbitals( get_number_sigma_orbitals());
	const core::Size number_sigma_orbitals_b( OTHER.get_number_sigma_orbitals());

	if
			(
					get_nNumber_hybrid_lone_pairs() + get_number_unhybridized_lone_pairs()
					!= OTHER.get_nNumber_hybrid_lone_pairs() + OTHER.get_number_unhybridized_lone_pairs()
					) {
		const core::Size number_unhybridized_lone_pairs( get_number_unhybridized_lone_pairs());
		const core::Size number_unhybridized_lone_pairs_b( OTHER.get_number_unhybridized_lone_pairs());
		const core::Size number_hybridized_lone_pairs( get_nNumber_hybrid_lone_pairs());
		const core::Size number_hybridized_lone_pairs_b( OTHER.get_nNumber_hybrid_lone_pairs());

		const core::Size lone_pairs( number_unhybridized_lone_pairs + number_hybridized_lone_pairs);
		const core::Size lone_pairs_b( number_unhybridized_lone_pairs_b + number_hybridized_lone_pairs_b);

		// is the difference between the types == the difference in the # of lone pairs?
		if
				(
						lone_pairs + charge_ != lone_pairs_b + OTHER.charge_
						|| get_number_hybrid_orbitals() != OTHER.get_number_hybrid_orbitals()
						) {
			// the type is different in multiple ways
			return Other;
		}

		if
				(
						get_number_bonds() - charge_ != OTHER.get_number_bonds() - OTHER.charge_
						&& ( number_electrons_in_bonds_ - get_number_bonds() - charge_) != ( OTHER.number_electrons_in_bonds_ - OTHER.get_number_bonds() - OTHER.charge_)
						) {
			return Other;
		}

		// yes
		return NumberLonePairOrbitals;
	}

	if ( number_sigma_orbitals + charge_ == number_sigma_orbitals_b + OTHER.charge_ ) {
		return NumberBondingSOrbitals;
	}

	// if the number of sigma orbitals remained the same, then the electrons must have exclusively went into the p-orbitals
	if ( number_sigma_orbitals == number_sigma_orbitals_b ) {
		return NumberBondingPOrbitals;
	}

	// the type is different in multiple ways
	return Other;
}

} // gasteiger
} // chemical
} // core
