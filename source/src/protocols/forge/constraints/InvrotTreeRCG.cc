// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/constraints/InvrotTreeRCG.cc
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, may 2012
/// @modified Tom Linsky, tlinsky@uw.edu, nov 2012

// unit headers
#include <protocols/forge/constraints/InvrotTreeRCG.hh>
#include <protocols/forge/constraints/InvrotTreeCstGeneratorCreator.hh>

// project headers
#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/toolbox/match_enzdes_util/AlignPoseToInvrotTreeMover.hh>
#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTree.hh>

// utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

// boost headers
#include <boost/assign.hpp>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer tr( "protocols.forge.constraints.InvrotTreeRCG" );

namespace protocols {
namespace forge {
namespace constraints {

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP InvrotTreeCstGeneratorCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new InvrotTreeRCG() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP InvrotTreeCstGeneratorCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return InvrotTreeRCG::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP InvrotTreeRCG::mover_name()
// XRW TEMP {
// XRW TEMP  return "InvrotTreeCstGenerator";
// XRW TEMP }

InvrotTreeRCG::InvrotTreeRCG()
: RemodelConstraintGenerator(),
	add_ligand_to_pose_( false ),
	invrot_tree_( /* NULL */ ),
	enzcst_io_( /* NULL */ ),
	geomcst_seqpos_( /* NULL */ ),
	setup_align_pose_( /* NULL */ )
{}

InvrotTreeRCG::InvrotTreeRCG( InvrotTreeRCG const & rval )
: RemodelConstraintGenerator( rval ),
	add_ligand_to_pose_( rval.add_ligand_to_pose_ ),
	invrot_tree_( rval.invrot_tree_ ),
	enzcst_io_( rval.enzcst_io_ ),
	geomcst_seqpos_( rval.geomcst_seqpos_ ),
	setup_align_pose_( rval.setup_align_pose_ )
{}

InvrotTreeRCG::InvrotTreeRCG(
	toolbox::match_enzdes_util::InvrotTreeOP invrot_tree,
	toolbox::match_enzdes_util::AllowedSeqposForGeomCstOP geomcst_seqpos
)
: invrot_tree_(invrot_tree),
	geomcst_seqpos_(geomcst_seqpos)
{}

InvrotTreeRCG::~InvrotTreeRCG(){}

void
InvrotTreeRCG::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	RemodelConstraintGenerator::parse_my_tag( tag, data, filters, movers, pose );
	//in case we'ref folding up around a ligand
	std::string cstfilename = tag->getOption<std::string>( "cstfile", "" );
	if ( cstfilename == "" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "XML tag for invrot tree does not contain a cstfile option... this will probably fail when the InvrotTreeCstGenerator mover runs." );
	}
	set_cstfile( cstfilename );

	// should we add in the ligand, or is it already there
	set_add_ligand_to_pose( tag->getOption<bool>( "add_ligand_to_pose", add_ligand_to_pose_ ) );
}

// XRW TEMP std::string
// XRW TEMP InvrotTreeRCG::get_name() const
// XRW TEMP {
// XRW TEMP  return InvrotTreeRCG::mover_name();
// XRW TEMP }

protocols::moves::MoverOP
InvrotTreeRCG::fresh_instance() const
{
	return protocols::moves::MoverOP( new InvrotTreeRCG() );
}

protocols::moves::MoverOP
InvrotTreeRCG::clone() const
{
	return protocols::moves::MoverOP( new InvrotTreeRCG( *this ) );
}

void
InvrotTreeRCG::apply( core::pose::Pose & pose )
{
	//tr << "IN apply, id=" << id() << std::endl;
	init( pose );
	setup_align_pose_->apply( pose );

	runtime_assert( invrot_tree_ != 0 );
	runtime_assert( geomcst_seqpos_ != 0 );
	runtime_assert( setup_align_pose_ != 0 );

	//tr << "now id=" << id() << std::endl;
	// generate and add constraints
	this->add_remodel_constraints_to_pose( pose );
}

void
InvrotTreeRCG::generate_remodel_constraints(
	core::pose::Pose const & pose )
{
	tr << "in invrottreercg" << std::endl;

	//for now
	runtime_assert( invrot_tree_->num_target_states() == 1 );
	Size target_state = 1;

	tr << "generating inverse_rotamer_constraints" << std::endl;
	invrot_tree_->generate_inverse_rotamer_constraints( pose, geomcst_seqpos_ );
	tr << "done generating inverse_rotamer_constraints" << std::endl;
	core::scoring::constraints::ConstraintCOP cst_to_add( invrot_tree_->get_constraint_for_target_state( target_state ) );

	if ( !cst_to_add ) utility_exit_with_message("InvrotTree failed to generate anything but NULL pointer csts. Something is wrong somewhere. Check your starting structure for whether it already has every interaction in the cstfile." );
	add_constraint( cst_to_add );
}

/// @brief tells the mover whether it should add the ligand to the pose
void
InvrotTreeRCG::set_add_ligand_to_pose( bool const add_lig )
{
	add_ligand_to_pose_ = add_lig;
}

/// @brief sets up the invrot_tree_ and enzcst_io_ from an enzdes constraint filename.
void
InvrotTreeRCG::set_cstfile( std::string const & cstfilename )
{
	toolbox::match_enzdes_util::EnzConstraintIOOP enzcst_io( new toolbox::match_enzdes_util::EnzConstraintIO( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ) );
	enzcst_io->read_enzyme_cstfile( cstfilename );
	invrot_tree_ = toolbox::match_enzdes_util::InvrotTreeOP( new protocols::toolbox::match_enzdes_util::TheozymeInvrotTree( enzcst_io ) );
	invrot_tree_->generate_targets_and_inverse_rotamers();
	enzcst_io_ = enzcst_io;
} //set_cstfile()

/// @brief initializes the rcg object
void
InvrotTreeRCG::init( core::pose::Pose const & pose )
{
	// this function is invoked from BluePrintBDR after the VLB object is attached, but before we enter centroid mode
	// it is also invoked when the mover is applied, if it is called from somewhere else
	geomcst_seqpos_ = toolbox::match_enzdes_util::AllowedSeqposForGeomCstOP( new protocols::toolbox::match_enzdes_util::AllowedSeqposForGeomCst() );
	//stupid: apparently we have to make a copy of the pose on the heap for
	//the initialization of allowed_seqpos to work
	core::pose::PoseOP posecopy( new core::pose::Pose( pose ) );
	geomcst_seqpos_->initialize_from_command_line( posecopy );

	// setup the align_pose movers
	// this one is called once before adding constraints
	setup_align_pose_ = toolbox::match_enzdes_util::AlignPoseToInvrotTreeMoverOP( new toolbox::match_enzdes_util::AlignPoseToInvrotTreeMover( invrot_tree_, geomcst_seqpos_ ) );
	setup_align_pose_->set_add_target_to_pose( add_ligand_to_pose_ );
	setup_align_pose_->set_geomcst_for_superposition_from_enz_io( enzcst_io_ );

	// if a VLB is going on, this one can be called as a user-provided mover
	run_align_pose_ = toolbox::match_enzdes_util::AlignPoseToInvrotTreeMoverOP( new toolbox::match_enzdes_util::AlignPoseToInvrotTreeMover( invrot_tree_, geomcst_seqpos_ ) );
	run_align_pose_->set_geomcst_for_superposition_from_enz_io( enzcst_io_ );

	protocols::forge::components::VarLengthBuildOP vlbop = vlb().lock();
	if ( vlbop ) {
		vlbop->loop_mover_fold_tree_constant( true ); //we're taking care of the fold tree through the above align movers
		vlbop->add_setup_mover( setup_align_pose_ );
		vlbop->add_user_provided_mover( run_align_pose_ );
	}
}

std::string InvrotTreeRCG::get_name() const {
	return mover_name();
}

std::string InvrotTreeRCG::mover_name() {
	return "InvrotTreeCstGenerator";
}

void InvrotTreeRCG::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	RemodelConstraintGenerator::attributes_for_remodel_constraint_generator( attlist );
	attlist
		+ XMLSchemaAttribute::required_attribute( "cstfile", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "add_ligand_to_pose", xsct_rosetta_bool, "XRW TO DO" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string InvrotTreeCstGeneratorCreator::keyname() const {
	return InvrotTreeRCG::mover_name();
}

protocols::moves::MoverOP
InvrotTreeCstGeneratorCreator::create_mover() const {
	return protocols::moves::MoverOP( new InvrotTreeRCG );
}

void InvrotTreeCstGeneratorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InvrotTreeRCG::provide_xml_schema( xsd );
}


} //namespace remodel
} //namespace forge
} //namespace protocols
