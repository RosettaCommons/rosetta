// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/init/init.cc
/// @brief  Declare WidgetRegistrators as static (global) variables in this .cc file
///         so that at load time, they will be initialized, and the Creator classes
///         they register will be handed to the appropriate WidgetFactory.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <protocols/init/init.hh>
#include <core/init/init.hh>

//Note protocols/init/init has been split into included-headers for organizational purposes.
//
#include <protocols/init/init.RequirementCreators.ihh>
#include <protocols/init/init.RequirementRegistrators.ihh>

#include <protocols/init/init.ConstraintCreators.ihh>
#include <protocols/init/init.ConstraintRegistrators.ihh>

#include <protocols/init/init.EnergyMethodCreators.ihh>
#include <protocols/init/init.EnergyMethodRegistrators.ihh>

#include <protocols/init/init.FallbackConfigurationCreators.ihh>
#include <protocols/init/init.FallbackConfigurationRegistrators.ihh>

#include <protocols/init/init.RDFFunctionCreators.ihh>
#include <protocols/init/init.RDFFunctionRegistrators.ihh>

#include <protocols/init/init.ResourceLoaderCreators.ihh>
#include <protocols/init/init.ResourceLoaderRegistrators.ihh>

#include <protocols/init/init.ResourceLocatorCreators.ihh>
#include <protocols/init/init.ResourceLocatorRegistrators.ihh>

#include <protocols/init/init.ResourceOptionsCreators.ihh>
#include <protocols/init/init.ResourceOptionsRegistrators.ihh>

#include <protocols/init/init.ResourceManagerCreators.ihh>
#include <protocols/init/init.ResourceManagerRegistrators.ihh>

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

#include <protocols/init/init.EvaluatorCreators.ihh>
#include <protocols/init/init.EvaluatorRegistrators.ihh>

#include <protocols/init/init.LoopsDefinerCreators.ihh>
#include <protocols/init/init.LoopsDefinerRegistrators.ihh>

#include <protocols/init/init.TaskOperationCreators.ihh>
#include <protocols/init/init.TaskOperationRegistrators.ihh>

#include <protocols/init/init.FilterCreators.ihh>
#include <protocols/init/init.FilterRegistrators.ihh>

#include <protocols/init/init.MoverCreators.ihh>
#include <protocols/init/init.MoverRegistrators.ihh>

#include <protocols/init/init.PosePropertyReporterCreators.ihh>
#include <protocols/init/init.PosePropertyReporterRegistrators.ihh>

#include <protocols/init/init.PoseSelectorCreators.ihh>
#include <protocols/init/init.PoseSelectorRegistrators.ihh>

#include <protocols/init/init.WriteableCacheableDataCreators.ihh>
#include <protocols/init/init.WriteableCacheableDataRegistrators.ihh>

namespace protocols {
namespace init {
void init( int argc, char * argv [] )
{
	core::init::init( argc, argv );
}

void init( utility::vector1< std::string > const & args )
{
	core::init::init( args );
}

} //namespace init
} //namespace protocols


