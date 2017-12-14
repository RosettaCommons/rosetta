// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/protein_interface_design/movers/DockAndRetrieveSidechains.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/DockAndRetrieveSidechains.hh>
#include <protocols/protein_interface_design/movers/DockAndRetrieveSidechainsCreator.hh>

#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <utility/tag/Tag.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/symmetric_docking/SymDockProtocol.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <core/pack/task/TaskFactory.hh>
#include <utility/string_util.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.DockAndRetrieveSidechains" );

// XRW TEMP std::string
// XRW TEMP DockAndRetrieveSidechainsCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return DockAndRetrieveSidechains::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP DockAndRetrieveSidechainsCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new DockAndRetrieveSidechains );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP DockAndRetrieveSidechains::mover_name()
// XRW TEMP {
// XRW TEMP  return "Docking";
// XRW TEMP }

DockAndRetrieveSidechains::DockAndRetrieveSidechains() :
	protocols::moves::Mover( DockAndRetrieveSidechains::mover_name() )
{}

DockAndRetrieveSidechains::~DockAndRetrieveSidechains() = default;

protocols::moves::MoverOP
DockAndRetrieveSidechains::clone() const{
	return( protocols::moves::MoverOP( new DockAndRetrieveSidechains( *this ) ) );
}

void
DockAndRetrieveSidechains::apply( core::pose::Pose & pose )
{
	// If the pose is not symmetric, then make it so
	if ( symmetry_ ) {
		protocols::simple_moves::symmetry::SetupForSymmetryMoverOP setup_mover( new protocols::simple_moves::symmetry::SetupForSymmetryMover );
		setup_mover->apply( pose );
	}

	core::pose::PoseCOP saved_pose( core::pose::PoseOP( new core::pose::Pose( pose ) ) );
	core::kinematics::FoldTree saved_ft( pose.fold_tree() );

	if ( symmetry_ ) {
		sym_docking_mover_->set_native_pose( saved_pose );
		sym_docking_mover_->set_input_pose( saved_pose );
		sym_docking_mover_->apply( pose );
		//set_last_move_status( sym_docking_mover_->get_last_move_status() ); might be needed
	} else {
		docking_mover_->set_native_pose( saved_pose );
		docking_mover_->set_input_pose( saved_pose );
		docking_mover_->apply( pose );
		//Allow this mover to see the status of DockingProtocol because it can fail and that failure needs to be known by this mover.
		set_last_move_status( docking_mover_->get_last_move_status() );
	}

	if ( low_res_protocol_only_ ) {
		protocols::simple_moves::SwitchResidueTypeSetMover to_all_atom( core::chemical::FA_STANDARD );
		protocols::simple_moves::ReturnSidechainMover recover_sidechains( *saved_pose );
		to_all_atom.apply( pose );
		recover_sidechains.apply( pose );
		pose.update_residue_neighbors();
	}
	if ( conserve_foldtree_ ) pose.fold_tree( saved_ft );

}


// XRW TEMP std::string
// XRW TEMP DockAndRetrieveSidechains::get_name() const {
// XRW TEMP  return DockAndRetrieveSidechains::mover_name();
// XRW TEMP }

void
DockAndRetrieveSidechains::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	std::string const score_low( tag->getOption<string>( "score_low", "score_docking_low" ) );
	std::string const score_high( protocols::rosetta_scripts::get_score_function_name( tag, "score_high" ) );
	low_res_protocol_only_ = !tag->getOption< bool >( "fullatom", false );
	conserve_foldtree_ = tag->getOption< bool >( "conserve_foldtree", false );
	bool const local_refine( tag->getOption<bool>( "local_refine", false ));
	bool const view( tag->getOption<bool>( "view", false ) );
	bool const design( tag->getOption<bool>( "design", false ) );
	symmetry_ = tag->getOption<bool>( "symmetry", false );
	//allow replacement of the default DockingTask with whatever you give it. Does not work with symmetry yet.
	bool const ignore_default_docking_task( tag->getOption<bool>( "ignore_default_docking_task", false ) );

	if ( symmetry_ ) {
		using namespace core::scoring::symmetry;
		ScoreFunctionOP scorelo = core::scoring::symmetry::symmetrize_scorefunction( *data.get_ptr< ScoreFunction >( "scorefxns", score_low ) );
		ScoreFunctionOP scorehi = core::scoring::symmetry::symmetrize_scorefunction( *data.get_ptr< ScoreFunction >( "scorefxns", score_high ));

		sym_docking_mover_ = protocols::symmetric_docking::SymDockProtocolOP( new protocols::symmetric_docking::SymDockProtocol( !low_res_protocol_only_, local_refine, view, scorelo, scorehi ) );

		sym_docking_mover_->task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
		sym_docking_mover_->design( design );
		TR << "symmetric docking mover with parameters low_res_protocol_only_ " << low_res_protocol_only_ << " local_refine " << local_refine
			<< " view "<< view << "  lowres_scorefxn= " << score_low
			<< "  highres_scorefxn= " << score_high << std::endl;
		return;
	}

	using namespace core::scoring;
	ScoreFunctionOP scorelo = data.get< ScoreFunction * >( "scorefxns", score_low )->clone();
	ScoreFunctionOP scorehi = data.get< ScoreFunction * >( "scorefxns", score_high )->clone();

	utility::vector1<std::string> jumps_str = utility::string_split( tag->getOption<string>( "jumps", "1" ), ',' );
	utility::vector1<core::Size> movable_jumps;
	for ( utility::vector1<std::string>::const_iterator it = jumps_str.begin(); it != jumps_str.end(); ++it ) {
		//movable_jumps.push_back( std::strtoul( *it, NULL, 0 ));
		movable_jumps.push_back( std::atoi( it->c_str() ));
	}
	//core::Size const rb_jump( tag->getOption< core::Size >( "rb_jump", 1 ) );
	bool const optimize_foldtree = tag->getOption<bool>( "optimize_fold_tree", true );
	docking_mover_ = protocols::docking::DockingProtocolOP( new protocols::docking::DockingProtocol( movable_jumps, low_res_protocol_only_, local_refine, optimize_foldtree, scorelo, scorehi ) );

	docking_mover_->set_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	docking_mover_->set_ignore_default_docking_task( ignore_default_docking_task );
	// //debugging
	// core::pack::task::TaskFactoryOP tf_debug(protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	// std::cout << "DockAndRetrieveSidechains task : \n" << *(tf_debug->create_task_and_apply_taskoperations( pose ) ) << std::endl;

	docking_mover_->set_design( design );
	TR<<"docking mover with parameters low_res_protocol_only_ "<< low_res_protocol_only_ <<" local_refine "<<local_refine<<" view "<<view<< "  lowres_scorefxn= " << score_low <<
		"  highres_scorefxn= " << score_high<<" over rb_jumps ";
	for ( utility::vector1<core::Size>::const_iterator it = movable_jumps.begin(); it != movable_jumps.end(); ++it ) {
		TR << *it << ",";
	}
	TR << " optimize fold tree="<<optimize_foldtree<<std::endl;
}

std::string DockAndRetrieveSidechains::get_name() const {
	return mover_name();
}

std::string DockAndRetrieveSidechains::mover_name() {
	return "Docking";
}

void DockAndRetrieveSidechains::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "score_low", xs_string, "Low-resolution scorefunction", "score_docking_low" )
		+ XMLSchemaAttribute( "score_high", xs_string, "High-resolution scorefunction" )
		+ XMLSchemaAttribute::attribute_w_default( "fullatom", xsct_rosetta_bool, "Run the high-resolution phase of the protocol", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "conserve_foldtree", xsct_rosetta_bool, "Keep the foldtree the same through the whole protocol", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "local_refine", xsct_rosetta_bool, "Do local refinement", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "view", xsct_rosetta_bool, "XRW TODO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "design", xsct_rosetta_bool, "XRW TODO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "symmetry", xsct_rosetta_bool, "Account for symmetry and associated information", "0" )

		+ XMLSchemaAttribute::attribute_w_default( "ignore_default_docking_task", xsct_rosetta_bool, "Ignore default docking task, instead using whatever is provided", "0" )
		;

	rosetta_scripts::attributes_for_parse_task_operations( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default( "jumps", xsct_int_cslist, "list of jumps for docking", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "optimize_fold_tree", xsct_rosetta_bool, "Obtain an optimal foldtree given the desired docking jumps", "1" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string DockAndRetrieveSidechainsCreator::keyname() const {
	return DockAndRetrieveSidechains::mover_name();
}

protocols::moves::MoverOP
DockAndRetrieveSidechainsCreator::create_mover() const {
	return protocols::moves::MoverOP( new DockAndRetrieveSidechains );
}

void DockAndRetrieveSidechainsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DockAndRetrieveSidechains::provide_xml_schema( xsd );
}



} //movers
} //protein_interface_design
} //protocols
