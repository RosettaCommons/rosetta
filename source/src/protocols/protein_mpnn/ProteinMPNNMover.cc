// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/protein_mpnn/ProteinMPNNMover.cc
/// @brief  Mover for the ProteinMPNN PyTorch model developed by Dauparas et al.
/// @author Frederick Chan (fredchan@uw.edu)

// Unit Headers
#include <protocols/protein_mpnn/ProteinMPNNMover.hh>
#include <protocols/protein_mpnn/ProteinMPNNMoverCreator.hh>

// Package Headers
#include <protocols/protein_mpnn/ProteinMPNN.hh>
#include <protocols/protein_mpnn/util.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResidueLevelTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/SimpleThreadingMover.hh>

// Utility Headers
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/Tracer.hh>
#include <core/pose/symmetry/util.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace protein_mpnn {

static basic::Tracer TR("protocols.protein_mpnn.ProteinMPNNMover");

/// @brief Destructor.
ProteinMPNNMover::~ProteinMPNNMover() = default;

void
ProteinMPNNMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {
	deterministic_flag_ = tag->getOption< bool >( "deterministic_flag", deterministic_flag_ );

	temperature_ = tag->getOption< core::Real >( "temperature", temperature_ );

	std::string omit_AAs_string = tag->getOption< std::string >( "omit_AAs", "X" );
	omit_AAs_ = utility::vector1< char >(omit_AAs_string.begin(), omit_AAs_string.end());

	// need to get bool mask vectors from these residue selectors after we have the pose,
	// so we need to do that in apply. However, we need to get the residue selectors here while
	// we have the DataMap
	std::string design_selector_rs_name = tag->getOption< std::string >( "design_selector", "" );
	if ( design_selector_rs_name == "" ) {
		design_selector_rs_ = nullptr;
	} else {
		design_selector_rs_ = core::select::residue_selector::get_residue_selector( design_selector_rs_name, data );
	}

	std::string coord_selector_rs_name = tag->getOption< std::string >( "coord_selector", "" );
	if ( coord_selector_rs_name == "" ) {
		coord_selector_rs_ = nullptr;
	} else {
		coord_selector_rs_ = core::select::residue_selector::get_residue_selector( coord_selector_rs_name, data );
	}

	tied_pos_rs_.clear();
	for ( utility::tag::TagCOP tied_def : tag->getTags( "TiedPositions" ) ) {
		if ( tied_def->hasOption( "residue_selectors" ) ) {
			utility::vector1< std::string > const res_selector_names =
				utility::string_split( tied_def->getOption< std::string >( "residue_selectors" ), ',' );

			utility::vector1< core::select::residue_selector::ResidueSelectorCOP > res_seles;
			for ( auto rs : res_selector_names ) {
				res_seles.emplace_back( core::select::residue_selector::get_residue_selector( rs, data ) );
			}

			tied_pos_rs_.emplace_back( res_seles );
		}
	}

	/// Parse TaskOps
	core::pack::task::TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	if ( new_task_factory != nullptr ) {
		task_factory_ = new_task_factory;
	}
}

void
ProteinMPNNMover::apply( core::pose::Pose & pose ) {
#ifndef USE_TORCH
	(void)pose; // avert error: unused parameter
	utility_exit_with_message( "ProteinMPNNMover is disabled in this build of Rosetta. To use ProteinMPNN you must compile Rosetta with extras=pytorch. After compiling with PyTorch, make sure you run a version of Rosetta with 'pytorch' in the name (rosetta_scripts.pytorch.linuxgccrelease, for example)." );
#else
	using namespace protocols::simple_moves;
	using namespace core::select;

	// Pass pose into ProteinMPNN
	ProteinMPNNOptions pose_options( pose );

	pose_options.deterministic_flag = deterministic_flag_;
	pose_options.temperature = temperature_;
	pose_options.omit_AAs = omit_AAs_;

	if ( design_selector_rs_ != nullptr ) {
		utility::vector1< bool > pos_mask = design_selector_rs_->apply( pose );
		pose_options.pos_mask = pos_mask;
	}

	if ( coord_selector_rs_ != nullptr ) {
		utility::vector1< bool > coord_mask = coord_selector_rs_->apply( pose );
		pose_options.coord_mask = coord_mask;
	}

	utility::vector1< utility::vector1< core::Size > > tied_positions;

	// Convert tied positions as specified in TiedPositions tags to ProteinMPNN tied positions
	for ( auto tied_res_selectors : tied_pos_rs_ ) {
		// getting indices from res selectors can occasionally be slow.
		// we store them here so we need to only compute them once
		utility::vector1 < utility::vector1< core::Size > > selection_positions = residue_selectors_to_indices(pose, tied_res_selectors);

		core::Size rs_size = selection_positions[1].size();
		validate_residue_selectors( selection_positions, rs_size );

		// for every residue of the selectors
		for ( core::Size residue_idx = 1; residue_idx <= rs_size; ++residue_idx ) {
			utility::vector1< core::Size > tied_resnums;

			// look at all the residue selectors in the tied_res_selectors and tie them together in order
			for ( auto pos_list : selection_positions ) {
				// put the residue number at index residue_idx of the residue selector into the tied_resnums
				tied_resnums.emplace_back( pos_list[ residue_idx ] );
			}

			tied_positions.emplace_back( tied_resnums );
		}
	}

	//// BEGIN TASK OPERATIONS SUPPORT
	if ( task_factory_ != nullptr ) {
		core::pack::task::PackerTaskCOP packer_task = task_factory_->create_task_and_apply_taskoperations( pose );

		for( core::Size ii = 1; ii <= pose.total_residue(); ++ii ){
			core::pack::task::ResidueLevelTask const & rlt( packer_task->residue_task( ii ) );

			// if no types are allowed, simply mask position
			if( rlt.allowed_residue_types().empty() ){
				pose_options.pos_mask[ ii ] = false;
				continue;
			}

			// create vector with name1 of allowed types
			utility::vector1< char > allowed_types;
			for( auto & restype : rlt.allowed_residue_types() ){
				allowed_types.push_back( restype->name1() );
			}

			// add disallowed AA types to omit_AAs_pos
			for( auto it = AA_ALPHABET.begin(); it != AA_ALPHABET.end(); ++it ){
				if( find( allowed_types.begin(), allowed_types.end(), *it ) == allowed_types.end() ){
					pose_options.omit_AAs_pos[ ii ].push_back( *it );
				}
			}
		}
	}

	//// END TASK OPERATIONS SUPPORT

	//// BEGIN SYMMETRY SUPPORT
	// Convert tied positions as specified in pose symmetry to ProteinMPNN tied positions
	// Also apply taskops to symmetric positions
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::conformation::symmetry::SymmetryInfoOP sym_info = core::pose::symmetry::symmetry_info( pose )->clone();

		// iterate over residues in asymmetric unit
		for ( core::Size ii = 1; ii <= sym_info->num_independent_residues(); ++ii ) {
			utility::vector1< core::Size > tied_set;
			tied_set.push_back( ii );

			// iterate over symmetric copies
			for ( core::Size subunit = 2; subunit <= sym_info->subunits(); ++subunit ) {
				core::Size clone = sym_info->equivalent_residue_on_subunit( subunit, ii );

				// tie symmetric positions
				tied_set.push_back( clone );

				// replicate mask and taskops at symmetric positions
				pose_options.pos_mask[ clone ] = pose_options.pos_mask[ ii ];
				pose_options.omit_AAs_pos[ clone ] = pose_options.omit_AAs_pos[ ii ];
			}

			tied_positions.push_back( tied_set );
		}
	}
	pose_options.tied_positions = tied_positions;

	//// END SYMMETRY SUPPORT

	std::string sample_seq = ProteinMPNN::get_instance()->sample( pose, pose_options );

	// Replace virtual residues with '-' to tell SimpleThreadingMover to skip the residue
	for ( core::Size ii = 0; ii < sample_seq.length(); ++ii ) {
		if ( sample_seq[ii] == 'X' ||
		     sample_seq[ii] == pose.residue( ii+1 ).name1() ) {
			sample_seq[ ii ] = '-';
		}
	}

	// Call SimpleThreadingMover with predicted sequence and input pose
	SimpleThreadingMoverOP threader( utility::pointer::make_shared< SimpleThreadingMover >(sample_seq, 1) );
	threader->set_pack_neighbors(false); // want to completely skip packing
	threader->set_pack_rounds(0);

	threader->apply( pose );

	TR << "Complete" << std::endl;
#endif //USE_TORCH
}

protocols::moves::MoverOP
ProteinMPNNMover::clone() const {
	return utility::pointer::make_shared< ProteinMPNNMover >( *this );
}

protocols::moves::MoverOP
ProteinMPNNMover::fresh_instance() const {
	return utility::pointer::make_shared< ProteinMPNNMover >();
}

void
ProteinMPNNMover::validate_residue_selectors(
	utility::vector1< utility::vector1< core::Size > > const & selection_positions,
	core::Size const expected_size
) {
	if ( selection_positions.size() < 2 ) {
		utility_exit_with_message(
			"Invalid TiedPositions with less than 2 residue selectors provided. "
			"A single TiedPositions tag should specify at least 2 residue selectors to tie together. ");
	}

	for ( auto pos_list : selection_positions ) {
		if ( pos_list.size() != expected_size ) {
			utility_exit_with_message(
				"Invalid TiedPositions with resiude selectors of different size provided. "
				"Tied residue selectors must have the same number of residues. ");
		}
	}
}

utility::vector1< utility::vector1< core::Size > >
ProteinMPNNMover::residue_selectors_to_indices(
	core::pose::Pose const & pose,
	utility::vector1< core::select::residue_selector::ResidueSelectorCOP > const & residue_selectors
) {
	utility::vector1< utility::vector1< core::Size > > selection_positions;
	for ( auto rs : residue_selectors ) {
		selection_positions.emplace_back(rs->selection_positions(pose));
	}
	return selection_positions;
}

void
ProteinMPNNMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
#ifndef USE_TORCH
	(void)citations; // avert error: unused parameter
#else
	using namespace basic::citation_manager;

	citations.add( ProteinMPNN::get_ProteinMPNN_neural_net_citation() );

	protocols::simple_moves::SimpleThreadingMover threader;
	threader.provide_citation_info(citations);
#endif //USE_TORCH
}

////////////// Creator /////////

std::string ProteinMPNNMover::get_name() const {
	return mover_name();
}

std::string ProteinMPNNMover::mover_name() {
	return "ProteinMPNNMover";
}

void
ProteinMPNNMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	// Attributes for ProteinMPNNMover
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"deterministic_flag", xsct_rosetta_bool,
		"Whether to enable deterministic mode. Setting to true overrides temperature. Defaults to false.", "false");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"temperature", xsct_real,
		"Sampling temperature. Defaults to 0.1.", "0.1");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"omit_AAs", xs_string,
		"List of globally disallowed AA types (plus unknown 'X'). Defaults to all AAs allowed except X.", "X");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"design_selector", xs_string,
		"Name of residue selector that selects designable positions. If none provided, all residues are designable.", "");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"coord_selector", xs_string,
		"Name of residue selector that selects per-residue coordinates to pass to ProteinMPNN. If none selected, all coordinates are passed.", "");

	// Attributes for ProteinMPNNMover>TiedPositions
	AttributeList tied_pos_attrs;
	tied_pos_attrs + XMLSchemaAttribute( "residue_selectors", xs_string,
		"Comma separated list of residue selectors to tie together. "
		"The first residues of each selector will be tied together, then the second, etc. "
		"Each residue selector must have the same number of residues.");

	XMLSchemaSimpleSubelementList attrs_subelements;
	attrs_subelements.add_simple_subelement( "TiedPositions", tied_pos_attrs,
		"Used to define tied sets of residues." );

	rosetta_scripts::attributes_for_parse_task_operations( attlist );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & protocols::moves::complex_type_name_for_mover )
		.element_name( mover_name() )
		.description( "Author: Frederick Chan (fredchan@uw.edu)\n"
		"Predict (design) AA identities using the ProteinMPNN deep learning model, and thread the predicted sequence onto the pose. Attempting to use this will result in a runtime error unless Rosetta was successfully compiled with `extras=pytorch`." )
		.add_attributes( attlist )
		.add_optional_name_attribute()
		.set_subelements_repeatable( attrs_subelements, 0 )
		.write_complex_type_to_schema( xsd );
}

std::string
ProteinMPNNMoverCreator::keyname() const {
	return ProteinMPNNMover::mover_name();
}

protocols::moves::MoverOP
ProteinMPNNMoverCreator::create_mover() const {
	return utility::pointer::make_shared< ProteinMPNNMover >();
}

void
ProteinMPNNMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ProteinMPNNMover::provide_xml_schema( xsd );
}

////////////// Option setters /////////

void
ProteinMPNNMover::set_deterministic_flag( bool deterministic_flag ) {
	deterministic_flag_ = deterministic_flag;
}

void
ProteinMPNNMover::set_temperature( core::Real temperature ) {
	temperature_ = temperature;
}

void
ProteinMPNNMover::set_omit_AAs( utility::vector1< char > omit_AAs ) {
	omit_AAs_ = omit_AAs;
}

void
ProteinMPNNMover::set_design_selector_rs( core::select::residue_selector::ResidueSelectorCOP design_selector_rs ) {
	design_selector_rs_ = design_selector_rs;
}

void
ProteinMPNNMover::set_coord_selector_rs( core::select::residue_selector::ResidueSelectorCOP coord_selector_rs ) {
	coord_selector_rs_ = coord_selector_rs;
}

void
ProteinMPNNMover::set_tied_pos_rs( utility::vector1< utility::vector1< core::select::residue_selector::ResidueSelectorCOP > > tied_pos_rs ) {
	tied_pos_rs_ = tied_pos_rs;
}

} // protein_mpnn
} // protocols
