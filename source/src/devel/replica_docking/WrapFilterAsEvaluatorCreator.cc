// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Zhe Zhang

// Unit Headers
#include <devel/replica_docking/WrapFilterAsEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>
#include <protocols/filters/Filter.hh>
#include <devel/replica_docking/InteractionScoreFilter.hh>
#include <devel/replica_docking/IrmsdFilter.hh>
#include <devel/replica_docking/FnatFilter.hh>
#include <devel/replica_docking/LrmsdFilter.hh>
#include <devel/replica_docking/FnonnatFilter.hh>
#include <devel/replica_docking/CaIrmsdFilter.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <devel/replica_docking/WrapFilterAsEvaluator.hh>
#include <core/io/silent/silent.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <utility/file/FileName.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// due to template function
#include <core/io/silent/SilentStruct.hh>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>

#include <utility/vector0.hh>

//Auto Headers


static THREAD_LOCAL basic::Tracer tr( "devel.replica_docking.WrapFilterAsEvaluatorCreator" );

namespace devel {
namespace replica_docking {

WrapFilterAsEvaluatorCreator::~WrapFilterAsEvaluatorCreator() {}

void WrapFilterAsEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::I_sc );
	OPT( evaluation::Irms );
	OPT( evaluation::Fnat );
	OPT( evaluation::Lrmsd );
	OPT( evaluation::Ca_Irms );
	OPT( evaluation::DockMetrics );

}

void WrapFilterAsEvaluatorCreator::add_evaluators( protocols::evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;
	using protocols::filters::FilterOP;


	if ( option[ OptionKeys::evaluation::I_sc ].user() || option[ OptionKeys::evaluation::DockMetrics ].user() ) {
		tr << "Add I_sc evaluator " << std::endl;
		protocols::filters::FilterOP f( new devel::replica_docking::InteractionScoreFilter() );
		eval.add_evaluation( PoseEvaluatorOP( new devel::replica_docking::WrapFilterAsEvaluator( f , "I_sc" ) ) );
	}
	if ( option[ OptionKeys::evaluation::Irms ].user() || option[ OptionKeys::evaluation::DockMetrics ].user() ) {
		tr << "Add Irms evaluator " << std::endl;
		eval.add_evaluation( PoseEvaluatorOP( new WrapFilterAsEvaluator( FilterOP( new devel::replica_docking::IrmsdFilter() ), "Irms" ) ) );
	}

	if ( option[ OptionKeys::evaluation::Ca_Irms ].user() || option[ OptionKeys::evaluation::DockMetrics ].user() ) {
		tr << "Add Irms evaluator " << std::endl;
		eval.add_evaluation( PoseEvaluatorOP( new WrapFilterAsEvaluator( FilterOP( new devel::replica_docking::CaIrmsdFilter() ), "Ca_Irms" ) ) );
	}
	if ( option[ OptionKeys::evaluation::Fnat ].user() || option[ OptionKeys::evaluation::DockMetrics ].user() ) {
		tr << "Add evaluator of Fnat " << std::endl;
		eval.add_evaluation( PoseEvaluatorOP( new WrapFilterAsEvaluator( FilterOP( new FnatFilter() ), "Fnat_n" ) ) );
	}
	if ( option[ OptionKeys::evaluation::Lrmsd ].user() || option[ OptionKeys::evaluation::DockMetrics ].user() ) {
		tr << "Add Lrmsd evaluator " << std::endl;
		eval.add_evaluation( PoseEvaluatorOP( new WrapFilterAsEvaluator( FilterOP( new devel::replica_docking::LrmsdFilter() ), "Lrmsd" ) ) );
	}
	if ( option[ OptionKeys::evaluation::Fnonnat ].user() || option[ OptionKeys::evaluation::DockMetrics ].user() ) {
		tr << "Add evaluator of Fnonnat " << std::endl;
		eval.add_evaluation( PoseEvaluatorOP( new WrapFilterAsEvaluator( FilterOP( new FnonnatFilter() ), "Fnonnat" ) ) );
	}
}

std::string WrapFilterAsEvaluatorCreator::type_name() const {
	return "WrapFilterAsEvaluator";
}

} //namespace
} //namespace
