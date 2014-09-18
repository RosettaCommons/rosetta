// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SAXSScoreFilter.cc
/// @brief runs reject or accept filters on pose
/// @detailed
///	  Contains currently: SAXSScoreFilter
///
///
/// @author Dominik Gront

// Unit Headers
#include <protocols/simple_filters/SAXSScoreFilter.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/options/keys/casp.OptionKeys.gen.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <basic/options/option.hh>
#include <basic/options/keys/filters.OptionKeys.gen.hh>

#include <core/scoring/saxs/SAXSEnergyFA.hh>

// Utility headers
#include <basic/Tracer.hh>

#include <core/chemical/ResidueType.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/scoring/saxs/FormFactorManager.hh>



//// C++ headers
static thread_local basic::Tracer tr( "protocols.simple_filters.SAXSScoreFilter" );

namespace protocols {
namespace simple_filters {

SAXSScoreFilter::SAXSScoreFilter() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

//    score_ = new protocols::scoring::methods::saxs::SAXSEnergyFA();
    cutoff_ = basic::options::option[ basic::options::OptionKeys::filters::set_saxs_filter ]();
    score_value_ = cutoff_ + 1;
}


bool SAXSScoreFilter::apply( core::pose::Pose const & pose ) const {

    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    core::pose::Pose fa_pose ( pose );
    core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
    if ( !pose.is_fullatom() )
		core::util::switch_to_residue_type_set( fa_pose, core::chemical::FA_STANDARD);
    if ( option[ casp::repack ].user() ) {
	core::pack::task::PackerTaskOP task(
                            core::pack::task::TaskFactory::create_packer_task( fa_pose ));
	task->initialize_from_command_line();
        task->restrict_to_repacking();
        core::pack::pack_rotamers(fa_pose, (*scorefxn), task);
    }
    if ( option[ casp::sc_min ] ) {
//	std::string const min_type( option[ run::min_type ]() );
	std::string min_type("dfpmin");
        core::kinematics::MoveMap final_mm;
        final_mm.set_chi( true );
        final_mm.set_bb( false );
        core::optimization::AtomTreeMinimizer().run( fa_pose, final_mm, *scorefxn,
    		core::optimization::MinimizerOptions( min_type, 0.001, true ) );
    }
    score_value_ = score.total_energy( fa_pose );

    if ( score_value_ < cutoff_ ) {
	tr.Info << " Passed with score " << score_value_ << std::endl;
	return true;
    }

    tr.Info << " Failed with score " << score_value_ << std::endl;
    return false;
}

} // filters
} // protocols
