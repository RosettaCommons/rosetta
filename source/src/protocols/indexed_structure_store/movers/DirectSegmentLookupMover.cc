// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/indexed_structure_store/movers/DirectSegmentLookupMover.cc

#include <boost/format.hpp>

#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <core/select/residue_selector/CachedResidueSubset.hh>

#include <protocols/indexed_structure_store/StructureStore.hh>
#include <protocols/indexed_structure_store/StructureStoreManager.hh>
#include <protocols/indexed_structure_store/search/QueryDatabase.hh>
#include <protocols/indexed_structure_store/vector_tools.hh>
#include <protocols/indexed_structure_store/pose_utility.hh>

#include <protocols/indexed_structure_store/movers/DirectSegmentLookupMover.hh>
#include <protocols/indexed_structure_store/movers/DirectSegmentLookupMoverCreator.hh>

#include <protocols/indexed_structure_store/utility.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols::indexed_structure_store::movers::DirectSegmentLookupMover" );

namespace protocols { namespace indexed_structure_store { namespace  movers {

moves::MoverOP
DirectSegmentLookupMover::clone() const {
	return utility::pointer::make_shared< DirectSegmentLookupMover >( *this );
}

void
DirectSegmentLookupMover::apply(core::pose::Pose & pose) {
	source_pose_ = utility::pointer::make_shared< core::pose::Pose >();
	source_pose_->detached_copy( pose );

	runtime_assert(from_chain_ <= source_pose_->conformation().num_chains() && from_chain_ >= 1);
	runtime_assert(to_chain_ <= source_pose_->conformation().num_chains() && to_chain_ >= 1);
	runtime_assert(from_chain_ != to_chain_);

	int n_term_pos = source_pose_->conformation().chain_end(from_chain_);
	int c_term_pos = source_pose_->conformation().chain_begin(to_chain_);

	DirectSegmentLookup lookup(lookup_config_);

	lookup_results_ = lookup.segment_lookup(
		structure_store()->residue_entries,
		*structure_database(),
		*source_pose_,
		n_term_pos - (lookup_context_ - 1), n_term_pos + 1,
		c_term_pos, c_term_pos + lookup_context_
	);
	lookup_result_index_ = 0;

	// Sort results to return "most dense" centers first.
	std::sort(
		lookup_results_.begin(), lookup_results_.end(),
		[](DirectSegmentLookupResult const & a, DirectSegmentLookupResult const & b) {
			return a.query_results.size() > b.query_results.size();
		});

	if ( lookup_results_.size() == 0 ) {
		set_last_move_status( moves::MS_FAIL_DO_NOT_RETRY );
		return;
	}

	pose = *generate_result_pose(*source_pose_, lookup_results_[lookup_result_index_]);

	runtime_assert(
		!stored_subset_name().empty() ?
		core::select::residue_selector::CachedResidueSubset::from_pose_datacache(pose).has_subset(stored_subset_name()) : true);
}

core::pose::PoseOP
DirectSegmentLookupMover::get_additional_output() {
	if ( !(lookup_result_index_ < lookup_results_.size() - 1) ) {
		return NULL;
	}

	if ( !(lookup_result_index_ < max_num_results_) && (max_num_results_ > 0) ) {
		return NULL;
	}

	lookup_result_index_ += 1;

	return generate_result_pose(*source_pose_, lookup_results_[lookup_result_index_]);
}

core::pose::PoseOP
DirectSegmentLookupMover::generate_result_pose(core::pose::Pose & source_pose, DirectSegmentLookupResult lookup_result) {
	auto pose_chains = source_pose.split_by_chain();

	core::pose::PoseOP segment_pose = residue_entries_to_pose(v_to_a(lookup_result.result_residues));
	core::pose::PoseOP upstream_pose = pose_chains.at(from_chain_);
	core::pose::PoseOP downstream_pose = pose_chains.at(to_chain_);
	core::pose::PoseOP joined_pose = append_pose_with_overlap(
		*append_pose_with_overlap(*upstream_pose, *segment_pose, lookup_context_, DELETE_DOWNSTREAM),
		*downstream_pose, lookup_context_, DELETE_UPSTREAM
	);

	pose_chains.at(from_chain_).reset();
	pose_chains.at(to_chain_).reset();

	core::Size insertion_point = from_chain_ ? from_chain_ <= to_chain_ : to_chain_;
	pose_chains.at(insertion_point) = joined_pose;

	core::Size insertion_start_res = 1;
	for ( core::Size c = 1; c < insertion_point; ++c ) {
		if ( pose_chains[c] ) {
			insertion_start_res += pose_chains[c]->total_residue();
		}
	}
	insertion_start_res += upstream_pose->total_residue() - lookup_context_;
	core::Size insertion_end_res = insertion_start_res + segment_pose->total_residue();

	core::pose::PoseOP result = pose_chains.front();
	for ( core::pose::PoseOP chain_pose : boost::make_iterator_range(pose_chains.begin() + 1, pose_chains.end()) ) {
		if ( chain_pose ) {
			result->append_pose_by_jump(*chain_pose, result->total_residue());
		}
	}

	if ( !stored_subset_name_.empty() ) {
		using namespace core::select::residue_selector;
		CachedResidueSubset & stored_subsets = CachedResidueSubset::from_pose_datacache(*result);
		if ( !overwrite_stored_subset_ && stored_subsets.has_subset( stored_subset_name_ ) ) {
			utility_exit_with_message( "A stored residue subset with the name " + stored_subset_name_ + " already exists; you must set overwrite_stored_subset flag to true to overwrite." );
		}

		ResidueSubset inserted_subset(result->total_residue(), false);
		for ( core::Size r = insertion_start_res; r < insertion_end_res; ++r ) {
			inserted_subset.at(r) = true;
		}

		TR << boost::format("Storing residue subset %s as range [%i, %i)") % stored_subset_name_ % insertion_start_res % insertion_end_res << std::endl;

		stored_subsets.set_subset( utility::pointer::make_shared< ResidueSubset >(inserted_subset), stored_subset_name_);
	}

	if ( !label_insertion_.empty() ) {
		if ( !result->pdb_info() ) {
			TR << "Result pose does not have PDBInfo, skipping label insertion." << std::endl;
		}
		TR << boost::format("Seting residue label %s on range [%i, %i)") % label_insertion_ % insertion_start_res % insertion_end_res << std::endl;

		for ( core::Size r = insertion_start_res; r < insertion_end_res; ++r ) {
			result->pdb_info()->add_reslabel(r, label_insertion_);
		}
	}

	if ( !output_lookup_count_.empty() ) {
		TR.Info << boost::format("Storing segment lookup hit count %s : %i") % output_lookup_count_ % lookup_result.query_results.size() << std::endl;
		core::pose::setPoseExtraScore(*result, output_lookup_count_, lookup_result.query_results.size());
	}

	if ( !output_lookup_length_.empty() ) {
		TR.Info << boost::format("Storing segment lookup length %s : %i") % output_lookup_length_ % lookup_result.result_residues.size() << std::endl;
		core::pose::setPoseExtraScore(*result, output_lookup_length_, lookup_result.result_residues.size());
	}

	return result;
}

void DirectSegmentLookupMover::structure_store_path(std::string path) {
	if ( path != structure_store_path_ ) {
		structure_store_.reset();
		structure_database_.reset();
		structure_store_path_ = path;
	}
}

std::string DirectSegmentLookupMover::structure_store_path() {
	return structure_store_path_;
}

void DirectSegmentLookupMover::lookup_config(DirectSegmentLookupConfig config) {
	runtime_assert(config.rmsd_tolerance > 0);
	runtime_assert(config.segment_cluster_tolerance > 0);
	runtime_assert(config.max_insertion_length >= 0);
	lookup_config_ = config;
}

DirectSegmentLookupConfig
DirectSegmentLookupMover::lookup_config() {
	return lookup_config_;
}

StructureStoreOP DirectSegmentLookupMover::structure_store() {
	init_structure_store();
	return structure_store_;
}

search::StructureDatabaseOP DirectSegmentLookupMover::structure_database() {
	init_structure_store();
	return structure_database_;
}

void
DirectSegmentLookupMover::init_structure_store() {
	if ( !structure_store_ ) {
		runtime_assert(!structure_database_);
		structure_store_ = StructureStoreManager::get_instance()->load_structure_store(structure_store_path());
		structure_database_ = utility::pointer::make_shared< search::StructureDatabase >();
		structure_database_->initialize(structure_store_->residue_entries);
	}
}

void
DirectSegmentLookupMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
	structure_store_path( tag->getOption< std::string >( "structure_store" ) );

	DirectSegmentLookupConfig config;
	config.rmsd_tolerance = tag->getOption<numeric::Real>( "rmsd_tolerance", .5 );
	config.segment_cluster_tolerance = tag->getOption<numeric::Real>( "segment_cluster_tolerance", 1.5 );
	config.max_insertion_length = tag->getOption<search::Index>( "max_insertion_length", 5 );
	config.max_insertion_length = tag->getOption<search::Index>( "max_insertion_length", 5 );
	lookup_config(config);

	max_num_results_ = tag->getOption<search::Index>( "max_num_results", 0 );

	from_chain(tag->getOption<Size>("from_chain", 1));
	to_chain(tag->getOption<Size>("to_chain", 2));

	stored_subset_name_ = tag->getOption< std::string >( "stored_subset_name", "");
	overwrite_stored_subset_ = tag->getOption< bool >( "overwrite_stored_subset", overwrite_stored_subset_ );

	label_insertion_ = tag->getOption< std::string >( "label_insertion", "");
	output_lookup_count_ = tag->getOption< std::string >( "output_lookup_count", "");
	output_lookup_length_ = tag->getOption< std::string >( "output_lookup_length", "");
}

void
DirectSegmentLookupMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute("structure_store" , xs_string , "Path of target struture store.");
	attlist + XMLSchemaAttribute::attribute_w_default("rmsd_tolerance" , xsct_real, "Lookup endpoint rmsd tolerance.", ".5");
	attlist + XMLSchemaAttribute::attribute_w_default("segment_cluster_tolerance" , xsct_real, "Segment clustering rmsd tolerance.", "1.5");
	attlist + XMLSchemaAttribute::attribute_w_default("max_insertion_length" , xsct_non_negative_integer, "Maximum segment insertion length.", "5");
	attlist + XMLSchemaAttribute::attribute_w_default("max_num_results" , xsct_non_negative_integer, "Maximum number of result segments returned, no limit if 0.", "0");

	attlist + XMLSchemaAttribute::attribute_w_default("from_chain" , xsct_positive_integer, "Upstream chain for connection.", "1");
	attlist + XMLSchemaAttribute::attribute_w_default("to_chain" , xsct_positive_integer, "Downstream chain for connection.", "2");

	attlist + XMLSchemaAttribute::attribute_w_default( "stored_subset_name", xs_string, "Store inserted residues as residue subset, recall with StoredResidueSubset.", "");
	attlist + XMLSchemaAttribute::attribute_w_default( "overwrite_stored_subset", xsct_rosetta_bool, "Overwrite stored subset if already present.", "false");

	attlist + XMLSchemaAttribute::attribute_w_default( "label_insertion", xs_string, "Label inserted residues in pdbinfo.", "");
	attlist + XMLSchemaAttribute::attribute_w_default( "output_lookup_count", xs_string, "Report lookup result cluster size as output pose score.", "");
	attlist + XMLSchemaAttribute::attribute_w_default( "output_lookup_length", xs_string, "Report lookup result length as output pose score.", "");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Perform alignment-based segment lookup in target protein structure db. TODO fordas", attlist );
}

protocols::moves::MoverOP
DirectSegmentLookupMoverCreator::create_mover() const {
	return utility::pointer::make_shared< DirectSegmentLookupMover >();
}

std::string
DirectSegmentLookupMover::get_name() const {
	return mover_name();
}

std::string
DirectSegmentLookupMover::mover_name() {
	return "DirectSegmentLookupMover";
}

std::string
DirectSegmentLookupMoverCreator::keyname() const {
	return DirectSegmentLookupMover::mover_name();
}

void
DirectSegmentLookupMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	DirectSegmentLookupMover::provide_xml_schema(xsd);
}

} } }
