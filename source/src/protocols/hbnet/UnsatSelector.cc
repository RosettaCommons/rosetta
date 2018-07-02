// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/UnsatSelector.hh
/// @brief  A ResidueSelector that selects hydrogen bond acceptors or selectors that are not satisfied with h-bond
/// @author Parisa Hosseinzadeh (parisah@uw.edu)

// Unit headers
#include <protocols/hbnet/UnsatSelectorCreator.hh>
#include <protocols/hbnet/UnsatSelector.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

//Tracer headers
#include <basic/Tracer.hh>

// Package headers
#include <core/chemical/AtomType.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Atom.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/pose/selection.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <protocols/hbnet/HBNet_util.hh>
#include <protocols/hbnet/HBNet.hh>
#include <core/id/AtomID.hh>

// Utility Headers
#include <utility>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>
#include <iostream>
#include <string>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.hbnet.UnsatSelector" );


namespace protocols {
namespace hbnet {

/// @brief selector_specific.
///
core::select::residue_selector::ResidueSelectorOP
UnsatSelectorCreator::create_residue_selector() const
{
	return core::select::residue_selector::ResidueSelectorOP( new UnsatSelector );
}


std::string
UnsatSelectorCreator::keyname() const
{
	return UnsatSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
UnsatSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	UnsatSelector::provide_xml_schema( xsd );
}

/// @brief Constructor.
///
UnsatSelector::UnsatSelector() :
	hbond_energy_cutoff_(-0.5), //a relatviely mild cutoff as a default
	consider_mainchain_only_(true),
	acceptors_(true),
	scorefxn_(),
	legacy_(false)
	//TODO -- initialize all vars here.
{}

/// @brief Destructor.
///
UnsatSelector::~UnsatSelector() = default;

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the new object.
core::select::residue_selector::ResidueSelectorOP
UnsatSelector::clone() const
{
	return core::select::residue_selector::ResidueSelectorOP( new UnsatSelector(*this) );
}

// ---------- APPLY FUNCTION --------------------

core::select::residue_selector::ResidueSubset UnsatSelector::apply (core::pose::Pose const & pose ) const
{
	core::Size const nres( pose.total_residue() ); //The number of residues in the pose.

	utility::vector1< utility::vector1 < core::Size > > hbond_counts_by_atom (compute (pose));//calculates the hbond numbers

	core::select::residue_selector::ResidueSubset selected( nres, false ); //Output, initialized to a vector of "false".

	for ( core::Size i=1; i<=nres; ++i ) {//goes through all atoms and residues in the pose
		core::conformation::Residue const & rsd = pose.residue( i );
		for ( core::Size j=1; j <= rsd.natoms(); ++j ) {
			if ( acceptors_ ) {
				if ( rsd.atom_type(j).is_acceptor() ) { //skips the unrelevant atoms
					if ( hbond_counts_by_atom[i][j] == 0 ) {
						if ( pose.residue(i).atom_is_backbone(j) ) {
							selected[i]=true;
						}
					}
				}
			} else {
				if ( rsd.atom_is_polar_hydrogen(j) ) { //skips the unrelevant atoms
					if ( hbond_counts_by_atom[i][j] == 0 ) {
						if ( pose.residue(i).atom_is_backbone(j) ) {
							selected[i]=true;
						}
					}
				}
			}
		}
	}

	return selected;

}
// ---------- PUBLIC FUNCTIONS ------------------


/// @brief Parse XML tag (to use this Mover in RosettaScripts).
///
void UnsatSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &data)
{
	using namespace core::select::residue_selector;

	set_hbond_energy_cutoff( tag->getOption<core::Real>( "hbond_energy_cutoff", hbond_energy_cutoff() ) );
	set_consider_mainchain_only( tag->getOption<bool>( "consider_mainchain_only", consider_mainchain_only() ) );
	set_mode( tag->getOption<bool>( "check_acceptors", mode() ) );
	set_legacy( tag->getOption<bool>( "legacy", legacy() ) );

	//Parse the scorefunction:
	if ( tag->hasOption("scorefxn") ) { //This "if" is unfortunately necessary.  The parse function defaults to giving talaris back, rather than returning a null pointer, if nothing is specified in the tag.
		set_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, "scorefxn", data ) );
	}

}

// ---------- SETTERS AND GETTERS ---------------


/// @brief Set the threshold for considering something to be a
/// hydrogen bond.
void UnsatSelector::set_hbond_energy_cutoff(core::Real const &input_value)
{
	hbond_energy_cutoff_ = input_value;
}

/// @brief Get the threshold for considering something to be a
/// hydrogen bond.
core::Real const & UnsatSelector::hbond_energy_cutoff() const
{
	return hbond_energy_cutoff_;
}

/// @brief Set whether we only consider mainchain hydrogen bond
/// donors and acceptors.
void UnsatSelector::set_consider_mainchain_only ( bool const input_setting )
{
	consider_mainchain_only_ = input_setting;
}

/// @brief Get whether we only consider mainchain hydrogen bond
/// donors and acceptors.
bool UnsatSelector::consider_mainchain_only() const
{
	return consider_mainchain_only_;
}

/// @brief Set whether we are checking acceptors or selectors
void UnsatSelector::set_mode ( bool const input_setting )
{
	acceptors_ = input_setting;
}

/// @brief Get whether we are checking acceptors or selectors
bool UnsatSelector::mode() const
{
	return acceptors_;
}

/// @brief Set whether we are checking acceptors or selectors
void UnsatSelector::set_legacy ( bool const input_setting )
{
	legacy_ = input_setting;
}

/// @brief Get the legacy mode
/// default=false
bool UnsatSelector::legacy() const
{
	return legacy_;
}

/// @brief Set the scorefunction.
/// @details Clones the input.
void UnsatSelector::set_scorefxn(core::scoring::ScoreFunctionCOP sfxn_in)
{
	runtime_assert_string_msg( sfxn_in, "Error in core::select::residue_selector::UnsatSelector::set_scorefxn(): A null pointer was passed to this function." ); //Can't be null.
	scorefxn_ = sfxn_in->clone();
}

/// @brief Get the scorefunction.
///
core::scoring::ScoreFunctionCOP UnsatSelector::scorefxn() const { return scorefxn_; }

// ---------- PRIVATE FUNCTIONS -----------------
/// @brief The function that actually calculates the value that this filter returns, called by the apply(),
/// report(), and report_sm() functions.
/// @details Returns the number of atoms receiving more than the allowed number of hydrogen bonds.
utility::vector1< utility::vector1 < core::Size > > UnsatSelector::compute(core::pose::Pose const &pose) const
{
	using namespace core::select::residue_selector;
	using namespace core::scoring::hbonds;


	//The scorefunction to use:
	core::scoring::ScoreFunctionOP scorefxn( (scorefxn_ ? scorefxn_ : core::scoring::get_score_function()->clone() ) );
	core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose( pose, scorefxn );

	core::scoring::methods::EnergyMethodOptionsOP energy_options( new core::scoring::methods::EnergyMethodOptions(scorefxn->energy_method_options()) );

	core::pose::Pose pose_copy( pose );
	//Number of residues in the pose:
	core::Size const nres(pose_copy.total_residue());

	utility::vector1< utility::vector1 < core::Size > > hbond_counts_by_atom;
	hbond_counts_by_atom.resize( nres );

	core::select::residue_selector::ResidueSubset donor_selection, acceptor_selection;
	donor_selection.resize( nres, true);
	acceptor_selection.resize( nres, true);

	//Ugh.  I have to copy the pose to score it:
	(*scorefxn)(pose_copy);

	if ( legacy_ ) {
		energy_options->hbond_options().decompose_bb_hb_into_pair_energies(true);
		scorefxn->set_energy_method_options(*energy_options);

		//The hydrogen bond set:
		HBondSet hbond_set( energy_options->hbond_options() );
		hbond_set.setup_for_residue_pair_energies( pose_copy, false, consider_mainchain_only() );
		HBondDatabaseCOP hb_database( HBondDatabase::get_database( hbond_set.hbond_options().params_database_tag()));
		core::scoring::EnergyMap hbond_emap;

		// Select donor and acceptor subsets:note that these are right now defaulted to be the main pose but can be later changed.
		/*if ( donor_selector_ ) {
		donor_selection = donor_selector_->apply( pose_copy );
		} else {
		donor_selection.resize( nres, true);
		}
		if ( acceptor_selector_ ) {
		acceptor_selection = acceptor_selector_->apply( pose_copy );
		} else {
		acceptor_selection.resize( nres, true);
		}*/

		for ( core::Size ir=1; ir<=nres; ++ir ) { //Loop through all potential acceptor residues
			if ( !acceptor_selection[ir] ) continue; //Skip unselected residues
			hbond_counts_by_atom[ir].resize( pose_copy.residue(ir).natoms(), 0 );

			for ( core::Size jr=1; jr<=nres; ++jr ) { //Loop through all potential donor residues
				hbond_counts_by_atom[jr].resize( pose_copy.residue(jr).natoms(), 0 );
				if ( !donor_selection[jr] ) continue; //Skip unselected residues

				HBondSet pair_hbond_set( energy_options->hbond_options(), ir == jr ? 1 : 2);
				if ( ir == jr ) {
					if ( core::scoring::hbonds::calculate_intra_res_hbonds( pose_copy.residue(ir), energy_options->hbond_options() ) ) {
						identify_intra_res_hbonds( *hb_database, pose_copy.residue(ir), hbond_set.nbrs(ir), false, pair_hbond_set,false /*never exclude bb-bb*/, consider_mainchain_only_, consider_mainchain_only_, consider_mainchain_only_);
					} else {
						continue;
					}
				} else {
					identify_hbonds_1way( *hb_database, pose_copy.residue(jr), pose_copy.residue(ir), hbond_set.nbrs(jr), hbond_set.nbrs(ir), false,false /*never exclude bb-bb*/, consider_mainchain_only_, consider_mainchain_only_, consider_mainchain_only_, pair_hbond_set);
				}
				if ( pair_hbond_set.nhbonds() == 0 ) continue; //Continue if this pair of residues has no hydrogen bonds.

				for ( core::Size ii=1, iimax=pair_hbond_set.nhbonds(); ii<=iimax; ++ii ) {
					if ( !pair_hbond_set.allow_hbond(ii) ) continue;
					HBond cur_hbond( pair_hbond_set.hbond(ii) );
					core::Real const cur_energy( cur_hbond.energy() * cur_hbond.weight() );
					if ( cur_energy > hbond_energy_cutoff() ) continue; //Ignore hydrogen bonds too weak to count.
					if ( acceptors_ ) {
						if ( cur_hbond.acc_atm_is_protein_backbone() ) { //because this selector only cares about bb satisfaction
							++hbond_counts_by_atom[ir][cur_hbond.acc_atm()]; //Add to the acceptor count.
						}
					} else {
						if ( cur_hbond.don_hatm_is_protein_backbone() ) { //because this selector only cares about bb satisfaction
							++hbond_counts_by_atom[jr][cur_hbond.don_hatm()]; // Add to donor H atom
						}
					}
				}
			}
		}
	} else {

		HBondSet temp_hbond_set;

		core::scoring::hbonds::HBondOptions new_options( temp_hbond_set.hbond_options() );
		new_options.use_hb_env_dep(false);
		new_options.decompose_bb_hb_into_pair_energies(false);
		//new_options.bb_donor_acceptor_check(bb_exclusion); // don't use bb exclusion logic when penalizing unsatisfied -- ideally would only eclude N-H donors and not exclude C=O with only 1 h-bond
		new_options.exclude_intra_res_protein( false ); // we want to count this for unsat calc, by default they are excluded

		HBondSetOP full_hbond_set( new HBondSet(new_options) );

		core::scoring::hbonds::fill_hbond_set( pose_copy, false /* deriv */, *full_hbond_set, false /* exclude bb-bb */, consider_mainchain_only_ /* exclude bb-sc */, consider_mainchain_only_ /* exclude sc-bb */, true /* exclude sc-sc */ ); //currently we just consider bb

		for ( core::Size ir=1; ir<=nres; ++ir ) { //Loop through all potential acceptor residues
			//if ( !acceptor_selection[ir] ) continue; //Skip unselected residues
			hbond_counts_by_atom[ir].resize( pose_copy.residue(ir).natoms(), 0 );
			for ( core::Size jr=1; jr <= pose_copy.residue(ir).natoms(); ++jr ) {//looping through atoms
				if ( pose_copy.residue(ir).atom_is_backbone(jr) ) { //care only about backbones
					if ( acceptors_ ) { //checking acceptors
						if ( pose_copy.residue(ir).atom_type(jr).is_acceptor() ) {
							utility::vector1< core::scoring::hbonds::HBondCOP > hbonds_for_atom_jr(full_hbond_set->atom_hbonds(core::id::AtomID( jr, ir ),true));//getting hbonds to this particular atom
							utility::vector1< core::scoring::hbonds::HBondCOP > hbond_vec(0);
							for ( auto & j : hbonds_for_atom_jr ) {
								if ( !( protocols::hbnet::hbond_exists_in_vector( hbond_vec,j )) ) { //check that it's not counted twice
									hbond_vec.push_back( j );
								}
							}
							for ( core::Size i=0; i< hbond_vec.size(); ++i ) {//looping through all hbonds
								core::Real hb_energy((hbond_vec[i+1])->energy());//unweighted h-bond energy
								if ( hb_energy < hbond_energy_cutoff() ) { //only count that if the energy is fine
									++hbond_counts_by_atom[ir][jr];
								}
							}
						}
					} else { //checking donors now
						if ( pose_copy.residue(ir).atom_type(jr).is_polar_hydrogen() ) {
							utility::vector1< core::scoring::hbonds::HBondCOP > hbonds_for_atom_jr(full_hbond_set->atom_hbonds(core::id::AtomID( jr, ir ),true));//getting hbonds to this particular atom
							utility::vector1< core::scoring::hbonds::HBondCOP > hbond_vec(0);
							for ( auto & j : hbonds_for_atom_jr  ) {
								if ( !( protocols::hbnet::hbond_exists_in_vector( hbond_vec,j )) ) { //check that it's not counted twice
									hbond_vec.push_back(j);
								}
							}
							for ( core::Size i=0; i< hbond_vec.size(); ++i ) {//looping through all hbonds
								core::Real hb_energy((hbond_vec[i+1])->energy());//unweighted h-bond energy
								if ( hb_energy < hbond_energy_cutoff() ) { //only count that if the energy is fine
									++hbond_counts_by_atom[ir][jr];
								}
							}
						}
					}
				}
			}
		}

	}

	return hbond_counts_by_atom;

}
/// @brief Get the mover class name.
///
std::string
UnsatSelector::get_name() const {
	return UnsatSelector::class_name();
}

/// @brief Get the mover class name.
///
std::string
UnsatSelector::class_name() {
	return "Unsat";
}
/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
///
void
UnsatSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute(  "check_acceptors", xsct_rosetta_bool, "whether you want to count acceptors or donors" )
		+ XMLSchemaAttribute(  "legacy", xsct_rosetta_bool, "do you want the hbnet style hbond detection or legacy style" )
		+ XMLSchemaAttribute(  "hbond_energy_cutoff", xsct_real, "what is the hbond energy cutoff" )
		+ XMLSchemaAttribute(  "consider_mainchain_only", xsct_rosetta_bool, "should we consider just mainchains?" );

	protocols::rosetta_scripts::attributes_for_parse_score_function (attributes);
	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(),"selects hbond acceptors or donors that are not satisfied", attributes );
}

} //hbnet
} //protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::hbnet::UnsatSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( hbond_energy_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( consider_mainchain_only_ ) ); // bool
	arc( CEREAL_NVP( acceptors_ ) ); // bool
	arc( CEREAL_NVP( scorefxn_ ) ); // core::scoring::ScoreFunctionOP
	arc( CEREAL_NVP( legacy_ ) ); // bool

}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::hbnet::UnsatSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( hbond_energy_cutoff_ ); // core::Real
	arc( consider_mainchain_only_ ); // bool
	arc( acceptors_ ); // bool
	arc( scorefxn_ ); // core::scoring::ScoreFunctionOP
	arc( legacy_ ); // bool
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::hbnet::UnsatSelector );
CEREAL_REGISTER_TYPE( protocols::hbnet::UnsatSelector )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_hbnet_UnsatSelector )
#endif // SERIALIZATION
