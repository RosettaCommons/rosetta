// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/simple_metrics/SequenceRecoveryMetric.cc
/// @brief Calculate sequence recovery statistics on a protein, relative to a reference.
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <protocols/analysis/simple_metrics/SequenceRecoveryMetric.hh>
#include <protocols/analysis/simple_metrics/SequenceRecoveryMetricCreator.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/pose/ref_pose.hh>
#include <core/pack/task/xml_util.hh>

#include <protocols/residue_selectors/TaskSelector.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.analysis.simple_metrics.SequenceRecoveryMetric" );


namespace protocols {
namespace analysis {
namespace simple_metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SequenceRecoveryMetric::SequenceRecoveryMetric():
	core::simple_metrics::RealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SequenceRecoveryMetric::~SequenceRecoveryMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
SequenceRecoveryMetric::SequenceRecoveryMetric( SequenceRecoveryMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
SequenceRecoveryMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new SequenceRecoveryMetric( *this ) );

}

std::string
SequenceRecoveryMetric::name() const {
	return name_static();
}

std::string
SequenceRecoveryMetric::name_static() {
	return "SequenceRecoveryMetric";

}

std::string
SequenceRecoveryMetric::metric() const {
	return "seqrec";
}

void
SequenceRecoveryMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	SimpleMetric::parse_base_tag( tag );

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(core::select::residue_selector::parse_residue_selector( tag, datamap ));
	}
	if ( tag->hasOption("task_operations") ) {
		if ( tag->hasOption("residue_selector") ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Cannot set both residue_selector and task_operations in  SequenceRecoveryMetric.");
		}

		core::pack::task::TaskFactoryOP tf( core::pack::task::parse_task_operations( tag, datamap ) );
		core::select::residue_selector::ResidueSelectorOP sele( new protocols::residue_selectors::TaskSelector( tf, true, false, false ) );
		set_residue_selector( sele );
	}

	if ( tag->hasOption("residue_selector_ref") ) {
		set_residue_selector_ref(core::select::residue_selector::parse_residue_selector( tag, datamap, "residue_selector_ref" ) );
	}

	// Comparison pose.
	if ( tag->hasOption("reference_name") ) {
		if ( tag->getOption<bool>("use_native", false) || tag->hasOption("pssm") ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Cannot use use_native with reference_name in SequenceRecoveryMetric.");
		}
		ref_pose_ = core::pose::saved_reference_pose(tag, datamap, "reference_name");
		TR<<"Loaded reference pose: "<<tag->getOption< std::string >( "reference_name" )<<" with "<<ref_pose_->size()<<" residues"<<std::endl;
	} else if ( tag->getOption<bool>("use_native", false) ) {
		if ( datamap.has_resource("native_pose") ) {
			ref_pose_ = core::pose::saved_native_pose(datamap);
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Use native specified with SequenceRecoveryMetric, but no native structure has been specified.");
		}
	}

	if ( tag->hasOption("pssm") ) {
		load_pssm( tag->getOption< std::string >("pssm") );
	} else {
		if ( ref_pose_ == nullptr ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "SequenceRecoveryMetric must have a reference set. Use one of reference_name, use_native, or pssm.");
		}
	}

	set_use_ave_pssm( tag->getOption<bool>("use_ave_pssm", false ) );
}

void
SequenceRecoveryMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "use_native", xsct_rosetta_bool, "Use the native if present on the cmd-line.", "false" )
		+ XMLSchemaAttribute( "pssm", xs_string, "a filename of a blast formatted pssm file containing the sequence profile to use" )
		+ XMLSchemaAttribute::attribute_w_default( "use_ave_pssm", xsct_rosetta_bool, "Use an average value metric for PSSM, rather than a pass/fail", "false" );

	core::pose::attributes_for_saved_reference_pose( attlist );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector",
		"Calculate the sequence recovery for these residues" );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector_ref",
		"Selector for the reference pose. If not specified, assume there's a 1-to-1 correspondence between the active pose and the reference pose.");

	core::pack::task::attributes_for_parse_task_operations( attlist,
		"As a convenience, instead of a residue_selector for selecting which residues to count, you can pass task operations "
		"and the design residues will be used. The use of the residue_selector attribute is prefered." );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for measuring sequence recovery and adding it to the resulting score file.\n"
		"\n"
		"There's several options for how the sequence recovery is calculated, depending on what parameters are set.\n"
		"Each metric is only calculated over the set of residues specified by the residue selector.\n"
		"\n"
		"* Standard - Used when PSSM isn't set but reference pose is\n"
		"     This is a strict match/no-match fraction.\n"
		"* Pass/Fail PSSM - Used when PSSM is set and use_ave_pssm is not (does not use reference pose)\n"
		"     This is the PSSM recovery metric from DeLuca, Dorr and Meiler 2011 Biochem 50(40):8521\n"
		"     Residue identities with positive (or zero) values in the PSSM count as a match, those with negative vales as no-match.\n"
		"* Ave PSSM - Used when PSSM is set, use_ave_pssm is true, and no reference PDB is given\n"
		"     This value is the average of the values in the PSSM matrix for the residue identities.\n"
		"* Delta PSSM - Used when PSSM is set, use_ave_pssm is true, and a reference PDB is provided.\n"
		"     This value is the average of the change in value of the PSSM matrix (mut - ref)\n"
		"\n"
		"For PSSM metrics, it's assumed that the Pose numbering of both the main and reference structure matches the numbering of the PSSM.\n"
		, attlist);

}

core::Real
SequenceRecoveryMetric::calculate(core::pose::Pose const & pose) const {

	if ( res_select_ == nullptr ) {
		utility_exit_with_message("Must specify a residue subset over which to calculate the SequenceRecoveryMetric.");
	}

	utility::vector1< core::Size > reslist = core::select::get_residues_from_subset( res_select_->apply(pose) );

	utility::vector1< core::Size > reslist_ref = reslist;
	if ( (ref_pose_ != nullptr) && (res_select_ref_ != nullptr) ) {
		reslist_ref = core::select::get_residues_from_subset( res_select_ref_->apply(*ref_pose_) );
		if ( reslist.size() != reslist_ref.size() ) {
			utility_exit_with_message("Unequal number of residues specified in SequenceRecoveryMetric. " +
				std::to_string( reslist.size() ) + " versus " + std::to_string( reslist_ref.size() ) + " in the reference");
		}
	}

	core::Real numerator = 0.0;
	core::Real denominator = 0.0;

	std::string metric_type = "";

	for ( core::Size ii(1); ii <= reslist.size(); ++ii ) {
		core::Size pos( reslist[ii] );
		core::Size ref_pos( reslist_ref[ii] );

		if ( ref_profile_ != nullptr ) {
			// Don't attempt to do UNK AAs - there's no consistent mapping on PSSM inputs.
			if ( pose.aa(pos) == core::chemical::aa_unk ) { continue; }

			utility::vector1< std::string > const & alphabet( ref_profile_->alphabet() );
			if ( pos > ref_profile_->length() ) {
				utility_exit_with_message("Attempted to access position " + std::to_string(pos) + " which doesn't exist in profile.");
			}

			std::string const & pos_aa{ pose.residue(pos).name1() };

			int index = alphabet.index( pos_aa );
			if ( index == 0 ) {
				TR.Warning << "Couldn't find residue `" << pos_aa << "` for position " << pos << " in PSSM profile." << std::endl;
				continue;
			}

			core::Real prof_value = ref_profile_->prof_row(pos)[index];
			TR.Debug << "POS " << pos << " ident " << pos_aa << " value " << prof_value << std::endl;
			debug_assert( ref_profile_->negative_better() == false );

			if ( use_ave_pssm_ && ref_pose_ != nullptr ) {
				debug_assert( metric_type == "" || metric_type == "Delta PSSM" );
				metric_type = "Delta PSSM";

				std::string const & ref_aa{ ref_pose_->residue(ref_pos).name1() };
				int ref_index = alphabet.index( ref_aa );
				if ( ref_index == 0 ) {
					TR.Warning << "Couldn't find residue `" << pos_aa << "` for reference position " << pos << " in PSSM profile." << std::endl;
					continue;
				}
				core::Real ref_value = ref_profile_->prof_row(ref_pos)[ref_index];
				TR.Debug << "REF " << ref_pos << " ident " << ref_aa << " value " << ref_value << std::endl;

				numerator += prof_value - ref_value;

			} else if ( use_ave_pssm_ ) {
				debug_assert( metric_type == "" || metric_type == "Ave PSSM" );
				metric_type = "Ave PSSM";

				numerator += prof_value;

			} else {
				// Use DeLuca-style match/no-match
				debug_assert( metric_type == "" || metric_type == "Pass/Fail PSSM" );
				metric_type = "Pass/Fail PSSM";

				if ( prof_value >= 0 ) {
					numerator += 1.0;
				}
			}
			denominator += 1.0;

		} else if ( ref_pose_ != nullptr ) {
			debug_assert( metric_type == "" || metric_type == "Standard" );
			metric_type = "Standard";

			core::chemical::AA aa( pose.aa( pos ) );
			core::chemical::AA aa_ref( ref_pose_->aa( ref_pos ) );
			if ( aa == aa_ref ) {
				numerator += 1.0;
			}
			denominator += 1.0;

		} else {
			utility_exit_with_message("Must specify either a reference pose or pssm for SequenceRecoveryMetric.");
		}
	}

	TR << "Calculating sequence metric using the " << metric_type << " settings." << std::endl;
	return numerator/denominator;

}

void
SequenceRecoveryMetric::set_comparison_pose( core::pose::PoseCOP pose ) {
	ref_pose_ = pose;
}

void
SequenceRecoveryMetric::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP residue_selector ) {
	res_select_ = residue_selector;
}

void
SequenceRecoveryMetric::set_residue_selector_ref( core::select::residue_selector::ResidueSelectorCOP residue_selector ) {
	res_select_ref_ = residue_selector;
}

void
SequenceRecoveryMetric::load_pssm( std::string const & pssm_filename ) {
	ref_profile_ = core::sequence::SequenceProfileOP( new core::sequence::SequenceProfile );
	ref_profile_->read_from_file( utility::file::FileName( pssm_filename ) );
}

void
SequenceRecoveryMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SequenceRecoveryMetric::provide_xml_schema( xsd );

}

std::string
SequenceRecoveryMetricCreator::keyname() const {
	return SequenceRecoveryMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
SequenceRecoveryMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new SequenceRecoveryMetric );

}

} //protocols
} //analysis
} //simple_metrics


#ifdef    SERIALIZATION



template< class Archive >
void
protocols::analysis::simple_metrics::SequenceRecoveryMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::RealMetric>( this ) );
	arc( CEREAL_NVP( res_select_ ) );
	arc( CEREAL_NVP( ref_pose_ ) );
	arc( CEREAL_NVP( res_select_ref_ ) );
	arc( CEREAL_NVP( ref_profile_ ) );
	arc( CEREAL_NVP( use_ave_pssm_ ) );
}

template< class Archive >
void
protocols::analysis::simple_metrics::SequenceRecoveryMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::RealMetric >( this ) );
	arc( res_select_ );
	core::pose::PoseOP refpose;
	arc( refpose );
	ref_pose_ = refpose;
	arc( res_select_ref_ );
	arc( ref_profile_ );
	arc( use_ave_pssm_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::analysis::simple_metrics::SequenceRecoveryMetric );
CEREAL_REGISTER_TYPE( protocols::analysis::simple_metrics::SequenceRecoveryMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_analysis_simple_metrics_SequenceRecoveryMetric )
#endif // SERIALIZATION




