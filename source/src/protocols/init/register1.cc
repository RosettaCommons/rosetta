// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/init/register1.cc
/// @brief  Declare WidgetRegistrators as static (global) variables in this .cc file
///         so that at load time, they will be initialized, and the Creator classes
///         they register will be handed to the appropriate WidgetFactory.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/init/register1.hh>

#include <protocols/init/init.LegacyRequirementCreators.ihh>
#include <protocols/init/init.LegacyRequirementRegistrators.ihh>

#include <protocols/init/init.AssemblyRequirementCreators.ihh>
#include <protocols/init/init.AssemblyRequirementRegistrators.ihh>

#include <protocols/init/init.AssemblyScorerCreators.ihh>
#include <protocols/init/init.AssemblyScorerRegistrators.ihh>

#include <protocols/init/init.ConstraintCreators.ihh>
#include <protocols/init/init.ConstraintRegistrators.ihh>

#include <protocols/init/init.EnergyMethodCreators.ihh>
#include <protocols/init/init.EnergyMethodRegistrators.ihh>

#include <protocols/init/init.RDFFunctionCreators.ihh>
#include <protocols/init/init.RDFFunctionRegistrators.ihh>

#include <protocols/init/init.ResourceLoaderCreators.ihh>
#include <protocols/init/init.ResourceLoaderRegistrators.ihh>

#include <protocols/init/init.RotamerRecoveryCreators.ihh>
#include <protocols/init/init.RotamerRecoveryRegistrators.ihh>

#include <protocols/init/init.DataLoaderCreators.ihh>
#include <protocols/init/init.DataLoaderRegistrators.ihh>

#include <protocols/init/init.GridCreators.ihh>
#include <protocols/init/init.GridRegistrators.ihh>

#include <protocols/init/init.JobInputterCreators.ihh>
#include <protocols/init/init.JobInputterRegistrators.ihh>

#include <protocols/init/init.JobOutputterCreators.ihh>
#include <protocols/init/init.JobOutputterRegistrators.ihh>

#include <protocols/init/init.FeaturesReporterCreators.ihh>
#include <protocols/init/init.FeaturesReporterRegistrators.ihh>

namespace protocols {
namespace init {

void register1()
{
	// A no-op, but needed to make sure that the Creators/Registrators above get statically allocated.
	// (As they might not before a function in this compilation unit is called.
}

} //namespace init
} //namespace protocols


