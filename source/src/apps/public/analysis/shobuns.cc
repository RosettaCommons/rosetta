// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief computes, using the SHO solvation model, the buried unsatisfied atoms
///  among a target set of polar atoms
///
/// @author Andrea Bazzoli (ndrbzz@gmail.com)
///
/// The program accepts only the following combinations of program-specific flags:
///
/// 1. -tgt_amino <AMINO> -tgt_atom <ATOM> ,
///    to select all polar atoms named <ATOM> in all residues of type <AMINO>,
///    where <AMINO> is a one-letter amino acid code; if <AMINO> is equal to
///    "any", then the program selects atom <ATOM> from all residues that
///    have it;
///
/// 2. -tgt_res <TGTFIL> ,
///    to select all polar atoms from the residues specified in file <TGTFIL>.
///    The format of <TGTFIL> is specified in the comments to function
///    protocols::toolbox::pose_metric_calculators::residue_subset() in file
///    "protocols/toolbox/pose_metric_calculators/SHOBuriedUnsatisfiedPolarsCalculator.cc";
///
/// 3. [NO FLAGS],
///    to select all polar atoms in the pose.
///
/// NOTE: flags -tgt_amino, -tgt_atom, -tgt_res are all in namespace
///       pose_metrics:shobuns.
///
/// NOTE: flags that activate SHO within a scoring function, like
///       '-score:patch occ_Hbond_sol_exact_talaris', are irrelevant to
///       determining whether a polar atom is buried unsatisfied or not.


#include <protocols/toolbox/pose_metric_calculators/SHOBuriedUnsatisfiedPolarsCalculator.fwd.hh>
#include <protocols/toolbox/pose_metric_calculators/SHOBuriedUnsatisfiedPolarsCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/after_opts.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

static basic::Tracer TR( "apps.public.analysis.shobuns.main" );

using basic::options::option;

/// MAIN

int main( int argc, char * argv [] )
{
	try {

		//// build pose and score function
		devel::init(argc, argv);

		core::pose::Pose ps;
		std::string const input_pdb_name( basic::options::start_file() );
		core::import_pose::pose_from_file( ps, input_pdb_name , core::import_pose::PDB_file);

		core::scoring::ScoreFunctionOP scorefxn(core::scoring::get_score_function());
		(*scorefxn)(ps);

		//// initialize SHO buried unsatisfied calculator from command line
		using protocols::toolbox::pose_metric_calculators::SHOBuriedUnsatisfiedPolarsCalculator;
		using protocols::toolbox::pose_metric_calculators::SHOBuriedUnsatisfiedPolarsCalculatorOP;
		SHOBuriedUnsatisfiedPolarsCalculatorOP shobuns_calculator;

		using basic::options::OptionKeys::pose_metrics::shobuns::sho_cutoff;
		core::Real cutoff = option[sho_cutoff];

		using basic::options::OptionKeys::pose_metrics::shobuns::tgt_amino;
		if ( option[tgt_amino].user() ) {

			// select by amino acid type and atom name

			using basic::options::OptionKeys::pose_metrics::shobuns::tgt_atom;

			std::string amino = option[tgt_amino];
			std::string atom = option[tgt_atom];
			shobuns_calculator = SHOBuriedUnsatisfiedPolarsCalculatorOP(
				new SHOBuriedUnsatisfiedPolarsCalculator(
				cutoff, amino, atom, scorefxn));
		} else {

			// select by residue indexes

			using protocols::toolbox::pose_metric_calculators::residue_subset;
			using basic::options::OptionKeys::pose_metrics::shobuns::tgt_res;

			utility::vector1<core::Size> res;
			if ( option[tgt_res].user() ) {
				std::string resfil = option[tgt_res];
				residue_subset(resfil, res, ps);
			}
			shobuns_calculator = SHOBuriedUnsatisfiedPolarsCalculatorOP(
				new SHOBuriedUnsatisfiedPolarsCalculator(cutoff, res, scorefxn));
		}

		core::pose::metrics::CalculatorFactory::Instance().register_calculator(
			"shobuns", shobuns_calculator);

		//// compute and print SHO buried unsatisfied polar atoms
		shobuns_calculator->recompute_and_print(ps);

	} // try
catch (utility::excn::Exception const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
}

}
