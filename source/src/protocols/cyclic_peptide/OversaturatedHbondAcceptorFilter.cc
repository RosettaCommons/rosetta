// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/OversaturatedHbondAcceptorFilter.cc
/// @brief This filter flags poses containing more than two hydrogen bonds to an oxygen atom, a common pathology that results
/// from Rosetta's pairwise-decomposible scorefunction, which can't penalize excessive hydrogen bonds.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#include <protocols/cyclic_peptide/OversaturatedHbondAcceptorFilter.hh>
#include <protocols/cyclic_peptide/OversaturatedHbondAcceptorFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Atom.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyMap.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

static basic::Tracer TR( "protocols.cyclic_peptide.OversaturatedHbondAcceptorFilter" );

namespace protocols {
namespace cyclic_peptide {

/// @brief Default constructor.
///
OversaturatedHbondAcceptorFilter::OversaturatedHbondAcceptorFilter():
	protocols::filters::Filter( "OversaturatedHbondAcceptorFilter" ),
	max_allowed_oversaturated_(0),
	hbond_energy_cutoff_(-0.1),
	consider_mainchain_only_(true),
	donor_selector_(),
	acceptor_selector_(),
	scorefxn_()
	//TODO Initialize all private member variables here.
{}

/// @brief Copy constructor.
///
OversaturatedHbondAcceptorFilter::OversaturatedHbondAcceptorFilter(
	OversaturatedHbondAcceptorFilter const &src
):
	protocols::filters::Filter( "OversaturatedHbondAcceptorFilter" ),
	max_allowed_oversaturated_( src.max_allowed_oversaturated_ ),
	hbond_energy_cutoff_( src.hbond_energy_cutoff_ ),
	consider_mainchain_only_( src.consider_mainchain_only_ ),
	donor_selector_(), //Cloned below
	acceptor_selector_(), //Cloned below
	scorefxn_() //Cloned below
	//TODO Copy all private member variables here.
{
	if ( src.donor_selector_ ) donor_selector_ = src.donor_selector_->clone();
	if ( src.acceptor_selector_ ) acceptor_selector_ = src.acceptor_selector_->clone();
	if ( src.scorefxn_ ) scorefxn_ = src.scorefxn_->clone();
}

/// @brief Destructor (important for properly forward-declaring smart-pointer members)
///
OversaturatedHbondAcceptorFilter::~OversaturatedHbondAcceptorFilter() = default;


/// @brief Required in the context of the parser/scripting scheme.
/// @details Make a copy of this object and return an owning pointer to the copy.
protocols::filters::FilterOP
OversaturatedHbondAcceptorFilter::clone() const
{
	return protocols::filters::FilterOP( new OversaturatedHbondAcceptorFilter( *this ) );
}

/// @brief Required in the context of the parser/scripting scheme.
///
protocols::filters::FilterOP
OversaturatedHbondAcceptorFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new OversaturatedHbondAcceptorFilter );
}

// ---------- PUBLIC FUNCTIONS ------------------

std::string
OversaturatedHbondAcceptorFilter::get_name() const
{
	return OversaturatedHbondAcceptorFilter::class_name();
}

/// @brief Returns true if the structure passes the filter, false otherwise.
///
bool
OversaturatedHbondAcceptorFilter::apply(
	core::pose::Pose const &pose
) const {
	core::Size const oversaturated_atoms( compute(pose) );
	return ( oversaturated_atoms <= max_allowed_oversaturated() );
}

/// @brief Required for reporting score values.
///
core::Real
OversaturatedHbondAcceptorFilter::report_sm(
	core::pose::Pose const &pose
) const {
	return static_cast<core::Real>( compute(pose) );
}

/// @brief Allows printing a summary of this filter to a stream.
///
void
OversaturatedHbondAcceptorFilter::report(
	std::ostream &os,
	core::pose::Pose const &pose
) const {
	core::Size const oversaturated_atoms( compute(pose) );
	os << "The OversaturatedHbondAcceptorFilter reports that " << oversaturated_atoms
		<< " hydrogen bond acceptors in the pose are receiving more than the allowed number of hydrogen bonds.  Since the cutoff is "
		<< max_allowed_oversaturated() << ", the filter "
		<< (oversaturated_atoms > max_allowed_oversaturated() ? "fails." : "passes.");
}

/// @brief Parse XML tag (to use this Mover in RosettaScripts).
///
void
OversaturatedHbondAcceptorFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &data,
	protocols::filters::Filters_map const & /*filters_map*/,
	protocols::moves::Movers_map const & /*movers_map*/,
	core::pose::Pose const & /*pose*/)
{
	using namespace core::select::residue_selector;

	set_max_allowed_oversaturated( tag->getOption<core::Size>( "max_allowed_oversaturated", max_allowed_oversaturated() ) );
	set_hbond_energy_cutoff( tag->getOption<core::Real>( "hbond_energy_cutoff", hbond_energy_cutoff() ) );
	set_consider_mainchain_only( tag->getOption<bool>( "consider_mainchain_only", consider_mainchain_only() ) );

	//Parse the residue selectors:
	ResidueSelectorCOP donor( protocols::rosetta_scripts::parse_residue_selector( tag, data, "donor_selector" ) );
	ResidueSelectorCOP acceptor( protocols::rosetta_scripts::parse_residue_selector( tag, data, "acceptor_selector" ) );
	if ( donor ) set_donor_selector(donor);
	if ( acceptor ) set_acceptor_selector(acceptor);

	//Parse the scorefunction:
	if ( tag->hasOption("scorefxn") ) { //This "if" is unfortunately necessary.  The parse function defaults to giving talaris back, rather than returning a null pointer, if nothing is specified in the tag.
		set_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, "scorefxn", data ) );
	}

}

// ---------- SETTERS AND GETTERS ---------------

/// @brief Set the maximum allowed number of instances of an oversaturated
/// hydrogen bond acceptor.
void
OversaturatedHbondAcceptorFilter::set_max_allowed_oversaturated(
	core::Size const max_allowed
) {
	max_allowed_oversaturated_ = max_allowed;
}

/// @brief Get the maximum allowed number of instances of an oversaturated
/// hydrogen bond acceptor.
core::Size
OversaturatedHbondAcceptorFilter::max_allowed_oversaturated() const {
	return max_allowed_oversaturated_;
}

/// @brief Set the threshold for considering something to be a
/// hydrogen bond.
void
OversaturatedHbondAcceptorFilter::set_hbond_energy_cutoff(
	core::Real const &input_value
) {
	hbond_energy_cutoff_ = input_value;
}

/// @brief Get the threshold for considering something to be a
/// hydrogen bond.
core::Real const &
OversaturatedHbondAcceptorFilter::hbond_energy_cutoff() const {
	return hbond_energy_cutoff_;
}

/// @brief Set whether we only consider mainchain hydrogen bond
/// donors and acceptors.
void
OversaturatedHbondAcceptorFilter::set_consider_mainchain_only (
	bool const input_setting
) {
	consider_mainchain_only_ = input_setting;
}

/// @brief Get whether we only consider mainchain hydrogen bond
/// donors and acceptors.
bool
OversaturatedHbondAcceptorFilter::consider_mainchain_only() const {
	return consider_mainchain_only_;
}

/// @brief Set the residue selector for donor residues.
/// @details Clones the input.
void
OversaturatedHbondAcceptorFilter::set_donor_selector(
	core::select::residue_selector::ResidueSelectorCOP donor_selector_in
) {
	runtime_assert_string_msg( donor_selector_in, "Error in protocols::cyclic_peptide::OversaturatedHbondAcceptorFilter::set_donor_selector(): A null pointer was passed to this function." ); //Can't be null.
	donor_selector_ = donor_selector_in->clone();
}

/// @brief Get the residue selector for donor residues.
///
core::select::residue_selector::ResidueSelectorCOP
OversaturatedHbondAcceptorFilter::donor_selector() const {
	return donor_selector_;
}

/// @brief Set the residue selector for acceptor residues.
/// @details Clones the input.
void
OversaturatedHbondAcceptorFilter::set_acceptor_selector(
	core::select::residue_selector::ResidueSelectorCOP acceptor_selector_in
) {
	runtime_assert_string_msg( acceptor_selector_in, "Error in protocols::cyclic_peptide::OversaturatedHbondAcceptorFilter::set_acceptor_selector(): A null pointer was passed to this function." ); //Can't be null.
	acceptor_selector_ = acceptor_selector_in->clone();
}

/// @brief Get the residue selector for acceptor residues.
///
core::select::residue_selector::ResidueSelectorCOP
OversaturatedHbondAcceptorFilter::acceptor_selector() const {
	return acceptor_selector_;
}

/// @brief Set the scorefunction.
/// @details Clones the input.
void
OversaturatedHbondAcceptorFilter::set_scorefxn(
	core::scoring::ScoreFunctionCOP sfxn_in
) {
	runtime_assert_string_msg( sfxn_in, "Error in protocols::cyclic_peptide::OversaturatedHbondAcceptorFilter::set_scorefxn(): A null pointer was passed to this function." ); //Can't be null.
	scorefxn_ = sfxn_in->clone();
}

/// @brief Get the scorefunction.
///
core::scoring::ScoreFunctionCOP OversaturatedHbondAcceptorFilter::scorefxn() const { return scorefxn_; }

// ---------- PRIVATE FUNCTIONS -----------------

/// @brief Given an atom, return the maximum number of hydrogen bonds that an atom of this type
/// is allowed to receive.
/// @details For now, this function returns 0 if this is not an acceptor, 2 if it is.  This is crude,
/// of course.  Ultimately, it should figure out how many hbonds an atom can actually receive.
core::Size
OversaturatedHbondAcceptorFilter::max_allowed_hbonds(
	core::pose::Pose const &pose,
	core::Size const res_index,
	core::Size const atom_index
) const {
	if ( !pose.residue(res_index).type().atom(atom_index).is_acceptor() ) return 0; //If it's not an acceptor, it can have only 0.
	return 2; //By default, assume that 2 is the max.
}

/// @brief The function that actually calculates the value that this filter returns, called by the apply(),
/// report(), and report_sm() functions.
/// @details Returns the number of atoms receiving more than the allowed number of hydrogen bonds.
core::Size
OversaturatedHbondAcceptorFilter::compute(
	core::pose::Pose const &pose
) const {
	using namespace core::select::residue_selector;
	using namespace core::scoring::hbonds;

	//The scorefunction to use:
	core::scoring::ScoreFunctionOP scorefxn( (scorefxn_ ? scorefxn_ : core::scoring::get_score_function()->clone() ) );
	core::scoring::methods::EnergyMethodOptionsOP energy_options( new core::scoring::methods::EnergyMethodOptions(scorefxn->energy_method_options()) );
	energy_options->hbond_options().decompose_bb_hb_into_pair_energies(true);
	scorefxn->set_energy_method_options(*energy_options);

	//Ugh.  I have to copy the pose to score it:
	core::pose::Pose pose_copy( pose );
	(*scorefxn)(pose_copy);

	//Number of residues in the pose:
	core::Size const nres(pose_copy.size());

	//The hydrogen bond set:
	HBondSet hbond_set( energy_options->hbond_options() );
	hbond_set.setup_for_residue_pair_energies( pose_copy, false, consider_mainchain_only() );
	HBondDatabaseCOP hb_database( HBondDatabase::get_database( hbond_set.hbond_options().params_database_tag()));
	core::scoring::EnergyMap hbond_emap;

	// Select donor and acceptor subsets:
	ResidueSubset donor_selection, acceptor_selection;
	if ( donor_selector_ ) {
		donor_selection = donor_selector_->apply( pose_copy );
	} else {
		donor_selection.resize( nres, true);
	}
	if ( acceptor_selector_ ) {
		acceptor_selection = acceptor_selector_->apply( pose_copy );
	} else {
		acceptor_selection.resize( nres, true);
	}


	utility::vector1< utility::vector1 < core::Size > > hbond_counts_by_atom;
	hbond_counts_by_atom.resize( nres );

	for ( core::Size ir=1; ir<=nres; ++ir ) { //Loop through all potential acceptor residues
		if ( !acceptor_selection[ir] ) continue; //Skip unselected residues
		hbond_counts_by_atom[ir].resize( pose_copy.residue(ir).natoms(), 0 );

		for ( core::Size jr=1; jr<=nres; ++jr ) { //Loop through all potential donor residues
			if ( !donor_selection[jr] ) continue; //Skip unselected residues

			HBondSet pair_hbond_set( energy_options->hbond_options(), ir == jr ? 1 : 2);
			if ( ir == jr ) {
				if ( core::scoring::hbonds::calculate_intra_res_hbonds( pose_copy.residue(ir), energy_options->hbond_options() ) ) {
					identify_intra_res_hbonds( *hb_database, pose_copy.residue(ir), hbond_set.nbrs(ir), false, pair_hbond_set,
						false /*never exclude bb-bb*/, consider_mainchain_only_, consider_mainchain_only_, consider_mainchain_only_);
				} else {
					continue;
				}
			} else {
				identify_hbonds_1way( *hb_database, pose_copy.residue(jr), pose_copy.residue(ir), hbond_set.nbrs(jr), hbond_set.nbrs(ir), false,
					false /*never exclude bb-bb*/, consider_mainchain_only_, consider_mainchain_only_, consider_mainchain_only_, pair_hbond_set);
			}
			if ( pair_hbond_set.nhbonds() == 0 ) continue; //Continue if this pair of residues has no hydrogen bonds.

			for ( core::Size ii=1, iimax=pair_hbond_set.nhbonds(); ii<=iimax; ++ii ) {
				if ( !pair_hbond_set.allow_hbond(ii) ) continue;
				HBond cur_hbond( pair_hbond_set.hbond(ii) );
				core::Real const cur_energy( cur_hbond.energy() * cur_hbond.weight() );
				if ( cur_energy > hbond_energy_cutoff() ) continue; //Ignore hydrogen bonds too weak to count.
				++hbond_counts_by_atom[ir][cur_hbond.acc_atm()]; //Add to the acceptor count.
			}
		}
	}

	//Loop through hbond_counts_by_atom and count the number of atoms with counts greater than should be allowed for the type:
	core::Size oversaturated_atoms(0); //Accumulator
	if ( TR.Debug.visible() ) {
		TR.Debug << "Atoms receiving too many hydrogen bonds:\nRES\tATOM\tHBONDS\tALLOWED\n";
	}
	for ( core::Size i=1, imax=hbond_counts_by_atom.size(); i<=imax; ++i ) {
		for ( core::Size j=1, jmax=hbond_counts_by_atom[i].size(); j<=jmax; ++j ) {
			if ( hbond_counts_by_atom[i][j] > 0 ) { //Don't bother looking up max allowed counts if not needed.
				if ( hbond_counts_by_atom[i][j] > max_allowed_hbonds( pose_copy, i, j ) ) {
					++oversaturated_atoms;
					if ( TR.Debug.visible() ) {
						TR.Debug << i << "\t" << pose_copy.residue(i).atom_name(j) << "\t" << hbond_counts_by_atom[i][j] << "\t" << max_allowed_hbonds( pose_copy, i, j ) << "\n";
					}
				}
			}
		}
	}
	if ( TR.Debug.visible() ) {
		TR.Debug << "TOTAL NUMBER OF ATOMS WITH TOO MANY HBONDS: " << oversaturated_atoms << std::endl;
		TR.Debug.flush();
	}

	return oversaturated_atoms;
}


/////////////// Creator ///////////////

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP OversaturatedHbondAcceptorFilterCreator::create_filter() const
// XRW TEMP {
// XRW TEMP  return protocols::filters::FilterOP( new OversaturatedHbondAcceptorFilter );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP OversaturatedHbondAcceptorFilterCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return OversaturatedHbondAcceptorFilter::class_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP OversaturatedHbondAcceptorFilter::class_name()
// XRW TEMP {
// XRW TEMP  return "OversaturatedHbondAcceptorFilter";
// XRW TEMP }

std::string OversaturatedHbondAcceptorFilter::name() const {
	return class_name();
}

std::string OversaturatedHbondAcceptorFilter::class_name() {
	return "OversaturatedHbondAcceptorFilter";
}

void OversaturatedHbondAcceptorFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute(
		"max_allowed_oversaturated", xsct_non_negative_integer,
		"How many oversaturated acceptors are allowed before the filter fails? "
		"Default 0 (filter fails if any oversaturated acceptors are found)." )
		+ XMLSchemaAttribute(
		"hbond_energy_cutoff", xsct_real,
		"A hydrogen bond must have energy less than or equal to this threshold "
		"in order to be counted. Default -0.1 Rosetta energy units." )
		+ XMLSchemaAttribute(
		"consider_mainchain_only", xsct_rosetta_bool,
		"If true (the default), only mainchain-mainchain hydrogen bonds are "
		"considered. If false, all hydrogen bonds are considered." );
	core::select::residue_selector::attributes_for_parse_residue_selector(
		attlist, "acceptor_selector",
		"Selector that defines the hydrogen bond acceptor" );
	core::select::residue_selector::attributes_for_parse_residue_selector(
		attlist, "donor_selector",
		"Selector that defines the hydrogen bond donor" );
	rosetta_scripts::attributes_for_parse_score_function( attlist, "scorefxn" );

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"This filter counts the number of hydrogen bond acceptors that are "
		"receiving hydrogen bonds from more than the allowed number of donors. "
		"If the count is greater than a threshold value (default 0), the filter fails. "
		"This filter is intended to address a limitation of Rosetta's pairwise-"
		"decomposible hydrogen bond score terms: it is not possible for a score term "
		"that is examining the interaction between only two residues to know that a "
		"third is also interacting with one of the residues, creating artifacts in "
		"which too many hydrogen bond donors are all forming hydrogen bonds to the same acceptor.",
		attlist );
}

std::string OversaturatedHbondAcceptorFilterCreator::keyname() const {
	return OversaturatedHbondAcceptorFilter::class_name();
}

protocols::filters::FilterOP
OversaturatedHbondAcceptorFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new OversaturatedHbondAcceptorFilter );
}

void OversaturatedHbondAcceptorFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	OversaturatedHbondAcceptorFilter::provide_xml_schema( xsd );
}


} //protocols
} //cyclic_peptide
