// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Restrict design to residues with electron density correlation above threshold value
/// @author Patrick Conway

// Unit Headers
#include <protocols/toolbox/task_operations/SelectByDensityFitOperation.hh>
#include <protocols/toolbox/task_operations/SelectByDensityFitOperationCreator.hh>

// Project Headers
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>

#include <basic/options/util.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <utility/exit.hh>

// C++ Headers

static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.task_operations.SelectByDensityFitOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

using namespace core;
using namespace basic;
using namespace utility;
using namespace protocols;
using namespace toolbox;
//using namespace task_operations;

core::pack::task::operation::TaskOperationOP
SelectByDensityFitOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new SelectByDensityFitOperation );
}

void SelectByDensityFitOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SelectByDensityFitOperation::provide_xml_schema( xsd );
}

std::string SelectByDensityFitOperationCreator::keyname() const
{
	return SelectByDensityFitOperation::keyname();
}

SelectByDensityFitOperation::SelectByDensityFitOperation( core::Real threshold, bool invert ):
	threshold_(threshold),
	invert_(invert)
{}

SelectByDensityFitOperation::~SelectByDensityFitOperation() {}

core::pack::task::operation::TaskOperationOP SelectByDensityFitOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new SelectByDensityFitOperation( *this ) );
}

void
SelectByDensityFitOperation::apply( core::pose::Pose const & const_pose, core::pack::task::PackerTask & task ) const
{
	using namespace protocols::moves;
	using namespace scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::pose::Pose pose = const_pose;
	core::Size nres = pose.size();

	protocols::electron_density::SetupForDensityScoringMoverOP dockindens( new protocols::electron_density::SetupForDensityScoringMover );
	dockindens->apply( pose );

	core::scoring::electron_density::getDensityMap().set_nres( nres );
	core::scoring::electron_density::getDensityMap().setScoreWindowContext( true );
	if ( option[ edensity::sliding_window ].user() ) {
		core::scoring::electron_density::getDensityMap().setWindow( option[ edensity::sliding_window ] );
	}

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	scorefxn->set_weight( core::scoring::elec_dens_window, 1.0 );
	(*scorefxn)(pose);

	for ( Size r=1; r<=nres; ++r ) {
		if ( task.nonconst_residue_task(r).being_designed() || task.nonconst_residue_task(r).being_packed() ) {
			Real score = core::scoring::electron_density::getDensityMap().matchRes( r , pose.residue(r), pose, NULL , false);
			TR.Debug << pose.pdb_info()->name() << " residue: " << r << " density_score: " << score << std::endl;
			if ( (score < threshold_) != invert_ ) {         // != invert_ flips the Boolean when true
				task.nonconst_residue_task(r).prevent_repacking();
			}
		}
	}
}

void
SelectByDensityFitOperation::parse_tag( TagCOP tag, DataMap & )
{
	threshold_ = tag->getOption<core::Real>("threshold", 0.72);
	invert_ = tag->getOption<bool>("invert", 0);
}

void SelectByDensityFitOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "threshold", xsct_real, "threshold value, electron density correlation must be above to pass",  "0.72"  )
		+ XMLSchemaAttribute::attribute_w_default(  "invert", xsct_rosetta_bool, "invert_ flips the Boolean when true",  "false"  );

	task_op_schema_w_attributes( xsd, keyname(), attributes, "Restrict design to residues with electron density correlation above threshold value." );
}


} //namespace task_operations
} //namespace toolbox
} //namespace protocols
