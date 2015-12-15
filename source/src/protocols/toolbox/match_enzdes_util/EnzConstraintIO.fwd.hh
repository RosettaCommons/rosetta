// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file IO-functionality for enzyme Constraints
/// @brief
/// @author Florian Richter, floric@u.washington.edu

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_EnzConstraintIO_fwd_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_EnzConstraintIO_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {


class EnzConstraintParameters;
class EnzCstTemplateRes;
class EnzCstTemplateResAtoms;
class EnzConstraintIO;
class CovalentConnectionReplaceInfo;

typedef utility::pointer::shared_ptr< EnzConstraintParameters > EnzConstraintParametersOP;
typedef utility::pointer::shared_ptr< EnzConstraintParameters const > EnzConstraintParametersCOP;
typedef utility::pointer::weak_ptr< EnzConstraintParameters > EnzConstraintParametersAP;
typedef utility::pointer::weak_ptr< EnzConstraintParameters const > EnzConstraintParametersCAP;

typedef utility::pointer::shared_ptr< EnzCstTemplateRes > EnzCstTemplateResOP;
typedef utility::pointer::shared_ptr< EnzCstTemplateRes const > EnzCstTemplateResCOP;

typedef utility::pointer::shared_ptr< EnzCstTemplateResAtoms > EnzCstTemplateResAtomsOP;
typedef utility::pointer::shared_ptr< EnzCstTemplateResAtoms const > EnzCstTemplateResAtomsCOP;

typedef utility::pointer::shared_ptr< EnzConstraintIO > EnzConstraintIOOP;
typedef utility::pointer::shared_ptr< EnzConstraintIO const > EnzConstraintIOCOP;
typedef utility::pointer::weak_ptr< EnzConstraintIO const > EnzConstraintIOCAP;
typedef utility::pointer::weak_ptr< EnzConstraintIO > EnzConstraintIOAP;

typedef utility::pointer::shared_ptr< CovalentConnectionReplaceInfo >CovalentConnectionReplaceInfoOP;
typedef utility::pointer::shared_ptr< CovalentConnectionReplaceInfo const > CovalentConnectionReplaceInfoCOP;

}
} // enzdes
} //protocols


#endif
