// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/PerResidueProbabilitiesMetric.cc
///
/// @brief Main class for a per residue probability simple metric.
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)

// Unit Headers
#include <core/simple_metrics/PerResidueProbabilitiesMetric.hh>
#include <core/simple_metrics/util.hh>

// Protocol Headers

// Core headers
#include <core/simple_metrics/SimpleMetricData.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/conformation/Residue.hh>

// // utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/ozstream.hh>

// STL headers
#include <iostream>
#include <boost/format.hpp>
#include <boost/range/algorithm/copy.hpp>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {
using namespace core::select::residue_selector;

PerResidueProbabilitiesMetric::PerResidueProbabilitiesMetric():
	SimpleMetric("PerResidueProbabilitiesMetric")
{
	selector_ = utility::pointer::make_shared< TrueResidueSelector >();
}

PerResidueProbabilitiesMetric::~PerResidueProbabilitiesMetric() = default;

PerResidueProbabilitiesMetric::PerResidueProbabilitiesMetric( PerResidueProbabilitiesMetric const & src ):
	SimpleMetric( src ),
	output_as_pdb_nums_( src.output_as_pdb_nums_ )
{
	if ( src.selector_ ) selector_ = src.selector_->clone();
}


PerResidueProbabilitiesMetric &
PerResidueProbabilitiesMetric::operator=( PerResidueProbabilitiesMetric const & ot ) {
	SimpleMetric::operator=( ot );
	if ( ot.selector_ ) selector_ = ot.selector_->clone();
	output_as_pdb_nums_ = ot.output_as_pdb_nums_;

	return *this;
}

utility::vector1< std::string >
PerResidueProbabilitiesMetric::get_metric_names() const {
	utility::vector1< std::string > names;
	names.push_back( metric() );
	return names;
}

void
PerResidueProbabilitiesMetric::set_residue_selector(select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

void
PerResidueProbabilitiesMetric::set_output_as_pdb_nums(bool output_as_pdb_nums){
	output_as_pdb_nums_ = output_as_pdb_nums;
}

select::residue_selector::ResidueSelectorCOP
PerResidueProbabilitiesMetric::get_selector() const {
	return selector_;
}

void
PerResidueProbabilitiesMetric::parse_per_residue_tag(utility::tag::TagCOP tag, basic::datacache::DataMap & datamap){
	set_output_as_pdb_nums(tag->getOption< bool >("output_as_pdb_nums", false));
	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(select::residue_selector::parse_residue_selector( tag, datamap ));
	}
}

///@brief Add options to the schema from this base class.
void
PerResidueProbabilitiesMetric::add_schema( utility::tag::XMLSchemaComplexTypeGeneratorOP ct_gen) {
	add_per_residue_simple_metric_schema( ct_gen );
}

std::map< core::Size, std::map< core::chemical::AA, core::Real>>
PerResidueProbabilitiesMetric::cached_calculate(
	pose::Pose const & pose,
	bool use_cache,
	std::string const & prefix,
	std::string const & suffix,
	bool fail_on_missing_cache,
	bool use_ref_pose_for_cache) const {

	std::string name = prefix + get_final_sm_type() + suffix;

	if ( use_cache && has_sm_data( pose ) ) {
		std::map< core::Size, std::map< core::chemical::AA, core::Real >> value;
		bool data_found = get_sm_data(pose)->get_value(name, value, pose, use_ref_pose_for_cache);
		if ( data_found ) {
			return value;
		} else if ( fail_on_missing_cache ) {
			utility_exit_with_message("Could not find PerResidueProbabilitiesMetric: "+name+" in pose");
		} else {
			return calculate(pose);
		}
	} else {
		return calculate(pose);
	}
}

void
PerResidueProbabilitiesMetric::apply( std::string const & out_tag, pose::Pose & pose, bool override_existing ) const {
	std::map< core::Size, std::map< core::chemical::AA, core::Real >> const values = calculate( pose ); //Index to value map

	MetricKey mk;

	std::map< core::Size, std::map< core::chemical::AA, core::Real >> stored_value;

	if ( ( ! override_existing ) && get_sm_data(pose)->get_value(out_tag, stored_value) ) {
		throw_sm_override_error(out_tag, name());
	}

	get_sm_data(pose)->set_value(mk, pose, out_tag, values, output_as_pdb_nums_);

}

/// @brief Return the probabilities in psi-blast position-specific-scoring-matrix (PSSM) format
/// @param[in] sequence The sequence of the pose
/// @param[in] logit_map A map containing the predicted logits for each position
/// @param[in] output_filename A string defining the name of the output file
void
PerResidueProbabilitiesMetric::output_sequence_profile( std::string const & sequence, std::map< core::Size, std::map< core::chemical::AA, core::Real >> const & logit_map,
	std::string const &output_filename) {
	// define the alphabet for the sequence profile
	utility::vector1<std::string> seqprof_alphabet_vec = {
		"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};
	core::sequence::SequenceProfile sequence_profile;
	sequence_profile.alphabet( seqprof_alphabet_vec ); // set alphabet
	// Same as above but as string now
	std::string const psiblast_alphabet = "ARNDCQEGHILKMFPSTWYV";
	// vector to fill with all probability rows for each position
	utility::vector1<utility::vector1<core::Real> > all_pssm_rows_vec(logit_map.size(),
		utility::vector1<core::Real>(20, 0));
	// go through positions and fill sequence profile
	for ( auto const & position_map_pair : logit_map ) {
		utility::vector1< core::Real > pssm_row_vec( 20, 0);
		for ( auto const & aa_logits_pair : position_map_pair.second ) {
			char cur_aa = core::chemical::oneletter_code_from_aa( aa_logits_pair.first );
			core::Size index = psiblast_alphabet.find( cur_aa ) + 1;
			pssm_row_vec[index] = aa_logits_pair.second;
		}
		all_pssm_rows_vec[position_map_pair.first] = pssm_row_vec;
	}
	sequence_profile.profile( all_pssm_rows_vec );
	sequence_profile.sequence( sequence );
	write_profile( sequence_profile, output_filename );
}

/// @brief Output the sequence_profile
/// @param[in] profile A SequenceProfile filled with the logits
/// @param[in] output_filename A string defining the name of the output file
void
PerResidueProbabilitiesMetric::write_profile(
	core::sequence::SequenceProfile & profile,
	std::string const & output_filename ) {

	utility::io::ozstream out( output_filename, std::ios::out | std::ios::binary );
	out << std::endl;
	out << "Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts" << std::endl;
	out << "            ";
	boost::copy( profile.alphabet(), std::ostream_iterator< std::string >( out, "  ") );
	out << std::endl;

	boost::format line_format(
		"%5i %1s   %3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f%3.0f");

	for ( core::Size i = 0; i < profile.length(); i++ ) {
		line_format.clear_binds();
		line_format.bind_arg( 1, i+1 );
		line_format.bind_arg( 2, profile.sequence()[i]);
		for ( core::Size j = 0; j < profile.width(); ++j ) {
			line_format.bind_arg(
				line_format.cur_arg(),
				profile.prof_row( i+1 )[ j+1 ]
			);
		}
		out << line_format << std::endl;
	}
	// rest of psi-blast pssm formatting, but its just placeholders
	out << std::endl;
	out << "                      K         Lambda" << std::endl;
	out << "Standard Ungapped    0.0000     0.0000" << std::endl;
	out << "Standard Gapped      0.0000     0.0000" << std::endl;
	out << "PSI Ungapped         0.0000     0.0000" << std::endl;
	out << "PSI Gapped           0.0000     0.0000" << std::endl;

	out.close();
}

} // simple_metrics
} // core

#ifdef    SERIALIZATION
template< class Archive >
void
core::simple_metrics::PerResidueProbabilitiesMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::SimpleMetric>( this ) );
	arc( CEREAL_NVP( selector_) );
	arc( CEREAL_NVP( output_as_pdb_nums_ ) );
}

template< class Archive >
void
core::simple_metrics::PerResidueProbabilitiesMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::SimpleMetric >( this ) );
	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector ); // ResidueSelectorCOP
	selector_ = local_selector; // copy the non-const pointer(s) into the const pointer(s)
	arc( output_as_pdb_nums_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::PerResidueProbabilitiesMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::PerResidueProbabilitiesMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_PerResidueProbabilitiesMetric )
#endif // SERIALIZATION


