// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/PeptideInternalHbondsMetric.cc
/// @brief A simple metric for measuring the hydrogen bonds within a peptide backbone, ignoring hydrogen
/// bonds to other parts of a pose.  Note that this metric does not count intra-residue hydrogen bonds.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Unit headers
#include <protocols/cyclic_peptide/PeptideInternalHbondsMetric.hh>
#include <protocols/cyclic_peptide/PeptideInternalHbondsMetricCreator.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/xml_util.hh>
#include <core/id/AtomID.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>
#include <utility/graph/Graph.hh>
#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

// Ugh.  ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.cyclic_peptide.PeptideInternalHbondsMetric" );


namespace protocols {
namespace cyclic_peptide {

/// @brief Duplicate this object, and return a smart pointer to the copy.
core::simple_metrics::SimpleMetricOP
PeptideInternalHbondsMetric::clone() const {
	return utility::pointer::make_shared< PeptideInternalHbondsMetric >( *this );
}

///@brief Calculate the metric.
core::Real
PeptideInternalHbondsMetric::calculate(
	core::pose::Pose const & pose
) const {
	core::select::residue_selector::ResidueSubset const selection(
		residue_selector_ == nullptr ?
		utility::vector1< bool >( pose.total_residue(), true ) :
		residue_selector_->apply(pose)
	);

	core::Size const count( count_hbonds(pose, selection) );
	TR << "Counted " << count << " hbonds in " << (residue_selector_ == nullptr ? "pose" : "selection") << "." << std::endl;
	return static_cast< core::Real >( count );
}

std::string
PeptideInternalHbondsMetric::name() const {
	return name_static();
}

std::string
PeptideInternalHbondsMetric::name_static() {
	return "PeptideInternalHbondsMetric";

}

/// @brief Return text indicating what we're measuring.
std::string
PeptideInternalHbondsMetric::metric() const {
	return "internal_hbonds";
}

void
PeptideInternalHbondsMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &datamap  )
{
	configure_from_filter_tag( tag, datamap );
	// Add any options that are unique to this metric but not shared by the filter here.
}

/// @brief Configure this metric from the tag passed to the companion filter (the
/// PeptideInternalHbondsFilter).
void
PeptideInternalHbondsMetric::configure_from_filter_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) {
	SimpleMetric::parse_base_tag( tag );

	set_hbond_types(
		tag->getOption< bool >( "backbone_backbone", backbone_backbone() ),
		tag->getOption< bool >( "backbone_sidechain", backbone_sidechain() ),
		tag->getOption< bool >( "sidechain_sidechain", sidechain_sidechain() )
	);
	set_exclusion_distance(
		tag->getOption< core::Size >( "exclusion_distance", exclusion_distance() )
	);
	set_hbond_energy_cutoff(
		tag->getOption< core::Real >( "hbond_energy_cutoff", hbond_energy_cutoff() )
	);
	if ( tag->hasOption("scorefxn") ) {
		set_scorefxn( core::scoring::parse_score_function( tag, "scorefxn", datamap ) );
	}
	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector( core::select::residue_selector::parse_residue_selector( tag, datamap, "residue_selector" ) );
	}
}

void
PeptideInternalHbondsMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	AttributeList attlist;
	provide_shared_xml_schema_elements(attlist);

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for counting the number of hydrogen bonds that are internal to a cyclic peptide or other residue selection.  Note that this metric does not count intra-residue hydrogen bonds.", attlist);
}

/// @brief Provide XML interface elements shared between this metric and the companion filter.
void
PeptideInternalHbondsMetric::provide_shared_xml_schema_elements(
	utility::tag::AttributeList & attlist
) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	attlist + XMLSchemaAttribute::attribute_w_default( "backbone_backbone", xsct_rosetta_bool, "If true, backbone-backbone hydrogen bonds are counted.  True by default.", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "backbone_sidechain", xsct_rosetta_bool, "If true, backbone-sidechain hydrogen bonds are counted.  False by default.", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "sidechain_sidechain", xsct_rosetta_bool, "If true, sidechain-sidechain hydrogen bonds are counted.  False by default.", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "exclusion_distance", xsct_non_negative_integer, "Hydrogen bonds between residues that are N positions apart in terms of covalent connectivity will not be counted if this is set higher than 0.  Set to 1 by default to exclude hydrogen bonds between neighbouring residues.", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "hbond_energy_cutoff", xsct_real, "The maximum energy of a hydrogen bond, if that hydrogen bond is to be counted.  Defaults to -0.25", "-0.25" );

	core::scoring::attributes_for_parse_score_function_w_description( attlist,  "scorefxn", "The scorefunction to use for computing hydrogen bonds.  If not provided, the default scoring function is used." );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "An optional residue selector that selects the peptide.  If provided, hydrogen bonds within this selection will be counted.  If not, the entire pose will be considered." );
}

/// @brief allows printing data to a stream.
void
PeptideInternalHbondsMetric::report(
	std::ostream &outstream,
	core::pose::Pose const &pose
) const {
	core::select::residue_selector::ResidueSubset const selection(
		residue_selector_ == nullptr ?
		utility::vector1< bool >( pose.total_residue(), true ) :
		residue_selector_->apply(pose)
	);
	outstream << "PeptideInternalHbondsMetric summary:\n";
	outstream << "\tSelected residues:" << (residue_selector_ == nullptr ? " ALL" : "" );
	if ( residue_selector_ != nullptr ) {
		for ( core::Size i(1), imax(selection.size()); i<=imax; ++i ) {
			if ( selection[i] ) outstream << " " << i;
		}
	}
	outstream << "\n\tInternal hydrogen bond count:" << count_hbonds( pose, selection ) << std::endl;
}

////////// Setters: //////////

/// @brief Set whether we're counting backbone-backbone, backbone-sidechain, and sidechain-sidechain hydrogen bonds.
/// @details Throws if all are set to false.  Defaults to only backbone-backbone being true.
void
PeptideInternalHbondsMetric::set_hbond_types(
	bool const backbone_backbone_setting,
	bool const backbone_sidechain_setting,
	bool const sidechain_sidechain_setting
) {
	runtime_assert_string_msg(
		backbone_backbone_setting || backbone_sidechain_setting || sidechain_sidechain_setting,
		"Error in PeptideInternalHbondsMetric::set_hbond_types(): At least one of backbone_backbone, backbone_sidechain, or sidechain_sidechain must be set to true!"
	);
	backbone_backbone_ = backbone_backbone_setting;
	backbone_sidechain_ = backbone_sidechain_setting;
	sidechain_sidechain_ = sidechain_sidechain_setting;
}

/// @brief Set the number of residues apart in terms of covalent connectivity that two residues have to be in order for their
/// hydrogen bonds to be counted.
/// @details Defaults to 1.
void
PeptideInternalHbondsMetric::set_exclusion_distance(
	core::Size const setting
) {
	exclusion_distance_ = setting;
}

/// @brief Set the energy cutoff for counting a hydrogen bond.
/// @details Defaults to -0.25.
void
PeptideInternalHbondsMetric::set_hbond_energy_cutoff(
	core::Real const setting
) {
	hbond_energy_cutoff_ = setting;
}

/// @brief Set the scorefunction for identifying hydrogen bonds.  If this is not called, the default scorefunction is used.
void
PeptideInternalHbondsMetric::set_scorefxn(
	core::scoring::ScoreFunctionCOP const & sfxn_in
) {
	scorefxn_ = sfxn_in;
}

/// @brief Set the scorefunction for selecting the peptide.  If this is not called, whole pose is counted.
void
PeptideInternalHbondsMetric::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP const & res_selector_in
) {
	residue_selector_ = res_selector_in;
}

////////// CITATION MANAGER FUNCTIONS //////////

/// @brief Does this simple metric provide information about how to cite it?
/// @details Returns false.
bool
PeptideInternalHbondsMetric::simple_metric_provides_citation_info() const {
	return false;
}

/// @brief Provide the citation.
/// @returns An empty vector for this simple metric, since it's unpublished.  Provides citations for the scorefunction
/// and residue selector, if available, however.
utility::vector1< basic::citation_manager::CitationCollectionCOP >
PeptideInternalHbondsMetric::provide_citation_info() const {
	utility::vector1< basic::citation_manager::CitationCollectionCOP > returnvec;

	if ( residue_selector_ != nullptr ) {
		basic::citation_manager::merge_into_citation_collection_vector( residue_selector_->provide_citation_info(), returnvec );
	}
	if ( scorefxn_ != nullptr ) {
		basic::citation_manager::merge_into_citation_collection_vector( scorefxn_->provide_citation_info(), returnvec );
	}
	return returnvec;
}

/// @brief Does this simple metric indicate that it is unpublished (and, by extension, that the author should be
/// included in publications resulting from it)?
/// @returns True, since this is unpublished.
bool
PeptideInternalHbondsMetric::simple_metric_is_unpublished() const {
	return true;
}

/// @brief Provide a list of authors and their e-mail addresses, as strings.
/// @details Also provides unpublished author information for residue selectors and scorefunctions used.
utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >
PeptideInternalHbondsMetric::provide_authorship_info_for_unpublished() const {
	utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP > returnvec{
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		name_static(),
		basic::citation_manager::CitedModuleType::SimpleMetric,
		"Vikram K. Mulligan",
		"Systems Biology, Center for Computational Biology, Flatiron Institute",
		"vmulligan@flatironinstitute.org"
		)
		};

	if ( residue_selector_ != nullptr ) {
		basic::citation_manager::merge_into_unpublished_collection_vector( residue_selector_->provide_authorship_info_for_unpublished(), returnvec );
	}
	if ( scorefxn_ != nullptr ) {
		basic::citation_manager::merge_into_unpublished_collection_vector( scorefxn_->provide_authorship_info_for_unpublished(), returnvec );
	}
	return returnvec;
}

//////////////// PRIVATE METHODS: ////////////////

/// @brief Given a pose and a selection, count the internal hydrogen bonds.
core::Size
PeptideInternalHbondsMetric::count_hbonds(
	core::pose::Pose const & original_pose,
	core::select::residue_selector::ResidueSubset const & selection
) const {
	std::string const errmsg( "Error in PeptideInternalHbondsMetric::count_hbonds(): " );

	// Make sure that we have a scorefunction.
	if ( scorefxn_ == nullptr ) {
		TR << "No scorefunction was provided to the PeptideInternalHbondsMetric.  Fetching default scorefunction." << std::endl;
		core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );
		core::scoring::methods::EnergyMethodOptions energy_options( sfxn->energy_method_options() );
		energy_options.hbond_options().decompose_bb_hb_into_pair_energies(true);
		energy_options.hbond_options().use_hb_env_dep(false);
		energy_options.exclude_intra_res_protein(true);
		energy_options.hbond_options().exclude_intra_res_protein(true);
		sfxn->set_energy_method_options(energy_options);
		scorefxn_ = sfxn;
	}

	// Check the scorefunction to ensure that appropriate hbonding terms are on.
	runtime_assert_string_msg(
		scorefxn_->energy_method_options().hbond_options().decompose_bb_hb_into_pair_energies(),
		errmsg + "The scoring function provided to the PeptideInternalHbondsMetric was not set to decompose backbone hydrogen bonds into pair energies!  This must be set."
	);
	if ( backbone_backbone_ ) {
		runtime_assert_string_msg(
			scorefxn_->has_nonzero_weight( core::scoring::hbond_lr_bb ) &&
			scorefxn_->has_nonzero_weight( core::scoring::hbond_sr_bb ),
			errmsg + "The scoring function used for counting hydrogen bonds must have the hbond_lr_bb and hbond_sr_bb terms turned on if we are to count backbone-backbone hydrogen bonds!"
		);
	}
	if ( backbone_sidechain_ ) {
		runtime_assert_string_msg(
			scorefxn_->has_nonzero_weight( core::scoring::hbond_bb_sc ),
			errmsg + "The scoring function used for counting hydrogen bonds must have the hbond_bb_sc term turned on if we are to count backbone-sidechain hydrogen bonds!"
		);
	}
	if ( sidechain_sidechain_ ) {
		runtime_assert_string_msg(
			scorefxn_->has_nonzero_weight( core::scoring::hbond_sc ),
			errmsg + "The scoring function used for counting hydrogen bonds must have the hbond_sc term turned on if we are to count sidechain-sidechain hydrogen bonds!"
		);
	}

	// Ugh.  Need to copy the pose to score it.
	core::pose::Pose pose( original_pose );

	// Score the pose:
	(*scorefxn_)(pose);

	// Check that the selection is sensible.
	runtime_assert( selection.size() == pose.total_residue() );

	// Generate a list of residue indices from the selection.
	utility::vector1 < core::Size > const resindices( core::select::residue_selector::selection_positions( selection ) );

	// Generate list of allowed partners for each residue in my list.
	// The outer index matches resindices, while the inner list is a set of possible partners.
	utility::vector1< std::list < core::Size > > const allowed_partners_for_each_residue( generate_allowed_partners( pose, resindices ) );

	// Loop through and count hydrogen bonds.
	core::Size hbond_counter(0); //Note: each hydrogen bond will be counted twice (once when examining the donor residue and once when examining the receiver).
	core::scoring::hbonds::HBondSet hbond_set( scorefxn_->energy_method_options().hbond_options() );
	core::scoring::hbonds::fill_hbond_set( pose, false, hbond_set, !backbone_backbone(), !backbone_sidechain(), !backbone_sidechain(), !sidechain_sidechain() );

	for ( core::Size i(1), imax(resindices.size()); i<=imax; ++i ) { //Iterating over potential donor and acceptor residues.
		core::Size const resindex( resindices[i] );
		std::list< core::Size > const & allowed_partners_for_this_residue( allowed_partners_for_each_residue[i] );

		core::chemical::ResidueType const & restype( pose.residue_type(resindex) );
		for ( core::Size ia(1), iamax(restype.natoms()); ia<=iamax; ++ia ) {
			bool const atom_is_bb( restype.atom_is_backbone(ia) );
			bool const atom_is_sc( restype.atom_is_sidechain(ia) );
			if ( (!backbone_backbone()) && (!backbone_sidechain()) && atom_is_bb ) continue;
			if ( (!backbone_sidechain()) && !(sidechain_sidechain()) && atom_is_sc ) continue;
			if ( restype.heavyatom_is_an_acceptor(ia) ) {
				//ACCEPTOR ATOMS:
				utility::vector1< core::scoring::hbonds::HBondCOP > const hbonds_for_atom(hbond_set.atom_hbonds(core::id::AtomID( ia, resindex ), false));
				for ( auto const & hbond : hbonds_for_atom ) {
					//Looping through all the hbonds that this atom makes.
					if ( std::find( allowed_partners_for_this_residue.begin(), allowed_partners_for_this_residue.end(), hbond->don_res() ) != allowed_partners_for_this_residue.end() ) {
						//We've got an hbond from a partner residue.  Now to check if it's kosher:
						if ( atom_is_bb ) {
							if ( (!backbone_sidechain()) && (!hbond->don_hatm_is_backbone()) ) continue;
							if ( (!backbone_backbone()) && hbond->don_hatm_is_backbone() ) continue;
						}
						if ( atom_is_sc ) {
							if ( (!sidechain_sidechain()) && (!hbond->don_hatm_is_backbone()) ) continue;
							if ( (!backbone_sidechain()) && hbond->don_hatm_is_backbone() ) continue;
						}
						//If we reach here, we've got a good hbond from another residue in the set to this residue.  Count it.
						++hbond_counter;
					}
				}
			}
			if ( restype.atom_is_polar_hydrogen(ia) ) {
				//DONOR ATOMS:
				utility::vector1< core::scoring::hbonds::HBondCOP > const hbonds_for_atom(hbond_set.atom_hbonds(core::id::AtomID( ia, resindex ), false));
				for ( auto const & hbond : hbonds_for_atom ) {
					//Looping through all the hbonds that this atom makes.
					if ( std::find( allowed_partners_for_this_residue.begin(), allowed_partners_for_this_residue.end(), hbond->acc_res() ) != allowed_partners_for_this_residue.end() ) {
						//We've got an hbond to a partner residue.  Now to check if it's kosher:
						if ( atom_is_bb ) {
							if ( (!backbone_sidechain()) && (!hbond->acc_atm_is_backbone()) ) continue;
							if ( (!backbone_backbone()) && hbond->acc_atm_is_backbone() ) continue;
						}
						if ( atom_is_sc ) {
							if ( (!sidechain_sidechain()) && (!hbond->acc_atm_is_backbone()) ) continue;
							if ( (!backbone_sidechain()) && hbond->acc_atm_is_backbone() ) continue;
						}
						//If we reach here, we've got a good hbond to another residue in the set from this residue.  Count it.
						++hbond_counter;
					}
				}
			}
		}

	}

	runtime_assert( hbond_counter % 2 == 0 );
	return hbond_counter / 2;
}

/// @brief Given a vector of selected resiudes and a pose, construct a vector of lists of allowed residue partners.
/// @details This is an O(N^2) operation, unfortunately.  However, the selected_residues list should be small, and N
/// is its length.
utility::vector1< std::list < core::Size > >
PeptideInternalHbondsMetric::generate_allowed_partners(
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & selected_residues
) const {
	utility::vector1< std::list <core::Size> > returnvec( selected_residues.size() );

	/// @brief Get the shortest paths between all residues.
	utility::graph::Graph connectivity_graph;
	connectivity_graph.set_num_nodes( pose.total_residue() );
	for ( core::Size ir(1), irmax( pose.total_residue()); ir<irmax; ++ir ) {
		core::conformation::Residue const & curres( pose.residue(ir));
		for ( core::Size iconn(1), iconnmax( curres.n_possible_residue_connections() ); iconn<=iconnmax; ++iconn ) {
			core::Size const other_res( curres.connected_residue_at_resconn(iconn) );
			if ( other_res > ir ) {
				connectivity_graph.add_edge( ir, other_res );
			}
		}
	}
	ObjexxFCL::FArray2D_int const path_lengths( connectivity_graph.all_pairs_shortest_paths() );  //Ugh.  O(N^3 operation).

	for ( core::Size i(1), imax(selected_residues.size()); i<=imax; ++i ) {
		std::list< core::Size > & curlist( returnvec[i] );
		core::Size const curres_i( selected_residues[i] );
		for ( core::Size j(1), jmax(selected_residues.size()); j<=jmax; ++j ) {
			core::Size const curres_j( selected_residues[j] );
			if ( curres_i == curres_j ) continue;
			if ( static_cast< core::Size >( path_lengths( curres_i, curres_j ) ) <= exclusion_distance() ) continue;
			curlist.push_back( curres_j );
		}
	}
	return returnvec;
}

//////////////////// CREATOR FUNCTIONS: ////////////////////

void
PeptideInternalHbondsMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PeptideInternalHbondsMetric::provide_xml_schema( xsd );
}

std::string
PeptideInternalHbondsMetricCreator::keyname() const {
	return PeptideInternalHbondsMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
PeptideInternalHbondsMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< PeptideInternalHbondsMetric >();
}

} //cyclic_peptide
} //protocols


#ifdef    SERIALIZATION

template< class Archive >
void
protocols::cyclic_peptide::PeptideInternalHbondsMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::RealMetric>( this ) );
	arc( CEREAL_NVP( backbone_backbone_ ) );
	arc( CEREAL_NVP( backbone_sidechain_ ) );
	arc( CEREAL_NVP( sidechain_sidechain_ ) );
	arc( CEREAL_NVP( exclusion_distance_ ) );
	arc( CEREAL_NVP( hbond_energy_cutoff_ ) );
	arc( CEREAL_NVP( scorefxn_ ) );
	arc( CEREAL_NVP( residue_selector_ ) );
}

template< class Archive >
void
protocols::cyclic_peptide::PeptideInternalHbondsMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::RealMetric >( this ) );
	arc( backbone_backbone_ );
	arc( backbone_sidechain_ );
	arc( sidechain_sidechain_ );
	arc( exclusion_distance_ );
	arc( hbond_energy_cutoff_ );

	//Have to do a little bit of a tap-dance around const OPs:
	core::scoring::ScoreFunctionOP sfxn_temp;
	arc( sfxn_temp );
	scorefxn_ = sfxn_temp; //Copy OP to COP.
	core::select::residue_selector::ResidueSelectorOP ressel_temp;
	arc( ressel_temp );
	residue_selector_ = ressel_temp; //Copy OP to COP.
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::cyclic_peptide::PeptideInternalHbondsMetric );
CEREAL_REGISTER_TYPE( protocols::cyclic_peptide::PeptideInternalHbondsMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_cyclic_peptide_PeptideInternalHbondsMetric )
#endif // SERIALIZATION




