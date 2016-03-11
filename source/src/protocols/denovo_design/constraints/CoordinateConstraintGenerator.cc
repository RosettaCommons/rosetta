// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/constraints/CoordinateConstraintGenerator.cc
/// @brief Coordinate constraint remodel constraint generator
/// @detailed
/// @author Tom Linsky

//Unit Headers
#include <protocols/denovo_design/constraints/CoordinateConstraintGenerator.hh>
#include <protocols/denovo_design/constraints/CoordinateConstraintGeneratorCreator.hh>

//Project Headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/loops/Loops.hh>

//Protocol Headers
#include <protocols/rosetta_scripts/util.hh>

//Core Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>

//Basic Headers
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers

//C++ Headers

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.constraints.CoordinateConstraintGenerator" );

namespace protocols {
namespace denovo_design {
namespace constraints {

std::string
CoordinateConstraintGeneratorCreator::keyname() const
{
	return CoordinateConstraintGeneratorCreator::mover_name();
}

protocols::moves::MoverOP
CoordinateConstraintGeneratorCreator::create_mover() const
{
	return protocols::moves::MoverOP( new CoordinateConstraintGenerator() );
}

std::string
CoordinateConstraintGeneratorCreator::mover_name()
{
	return "CoordinateConstraintGenerator";
}

CoordinateConstraintGenerator::CoordinateConstraintGenerator() :
	protocols::moves::ConstraintGenerator(),
	reference_pose_(),
	selector_( core::select::residue_selector::ResidueSelectorCOP( new core::select::residue_selector::TrueResidueSelector() ) ),
	weight_( 1.0 )
{
}

CoordinateConstraintGenerator::~CoordinateConstraintGenerator() {}

std::string
CoordinateConstraintGenerator::get_name() const
{
	return "CoordinateConstraintGenerator";
}

protocols::moves::MoverOP
CoordinateConstraintGenerator::clone() const
{
	return protocols::moves::MoverOP( new CoordinateConstraintGenerator( *this ) );
}

void
CoordinateConstraintGenerator::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	protocols::moves::ConstraintGenerator::parse_my_tag( tag, data, filters, movers, pose );
	if ( tag->hasOption( "reference_name" ) ) {
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose( tag, data );
	}
	if ( tag->hasOption( "weight" ) ) {
		weight_ = tag->getOption< core::Real >( "weight" );
	}
	core::select::residue_selector::ResidueSelectorCOP selector = protocols::rosetta_scripts::parse_residue_selector( tag, data );
	if ( selector ) selector_ = selector;
}

core::scoring::constraints::ConstraintCOPs
CoordinateConstraintGenerator::generate_constraints( core::pose::Pose const & pose )
{
	core::scoring::constraints::ConstraintCOPs csts;
	if ( reference_pose_ ) {
		if ( reference_pose_->total_residue() != pose.total_residue() ) {
			std::stringstream msg;
			msg << "CoordinateConstraintGenerator: reference pose size (" << reference_pose_->total_residue()
				<< ") does not match input pose size (" << pose.total_residue() << ")" << std::endl;
			throw utility::excn::EXCN_BadInput( msg.str() );
		}
		TR.Debug << "Adding coordinate constraints based on reference pose" << std::endl;
		add_constraints( csts, *reference_pose_ );
	} else {
		TR.Debug << "Adding coordinate constraints based on input pose" << std::endl;
		add_constraints( csts, pose );
	}
	return csts;
}

void
CoordinateConstraintGenerator::add_constraints( core::scoring::constraints::ConstraintCOPs & csts, core::pose::Pose const & pose ) const
{
	TR.Debug << "Selecting residues" << std::endl;
	debug_assert( selector_ );
	core::select::residue_selector::ResidueSubset subset = selector_->apply( pose );
	for ( core::Size resid=1; resid<=subset.size(); ++resid ) {
		if ( subset[ resid ] ) {
			TR.Debug << "Adding constraint for " << resid << std::endl;
			csts.push_back( create_coordinate_cst( pose, resid ) );
		}
	}
}

void
CoordinateConstraintGenerator::set_selector( core::select::residue_selector::ResidueSelectorCOP selector )
{
	selector_ = selector;
}

/// @brief creates a Ca coordinate constraint for residue resi
core::scoring::constraints::ConstraintOP
CoordinateConstraintGenerator::create_coordinate_cst(
	core::pose::Pose const & pose,
	core::Size const resi ) const
{
	core::Size atom( pose.residue_type(resi).nbr_atom() );
	if ( pose.residue_type(resi).has("CA") ) {
		atom = pose.residue_type(resi).atom_index("CA");
	}

	return core::scoring::constraints::ConstraintOP(
		new core::scoring::constraints::CoordinateConstraint(
		core::id::AtomID(atom,resi),
		core::id::AtomID(pose.residue(1).nbr_atom(),1),
		pose.residue(resi).xyz(atom),
		core::scoring::func::FuncOP( new core::scoring::func::HarmonicFunc(0.0, 0.5*weight_) ) ) );
}

}
}
}

