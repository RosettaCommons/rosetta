// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/smlewis/ChemXRW/PDB_diagnostic.cc
/// @brief This app is meant to report if Rosetta can successfully read a PDB, and if not, help with diagnostics on why it failed. Spiritual child of old AnchorFinder.

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>

#include <core/chemical/ResidueType.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <basic/options/option.hh>
//#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>

#include <devel/init.hh>

#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

//#include <utility/vector0.hh>
//#include <utility/vector1.hh>

OPT_1GRP_KEY( Boolean, PDB_diagnostic, reading_only ) //basic::options::OptionKeys::PDB_diagnostic::reading_only
OPT_1GRP_KEY( Boolean, PDB_diagnostic, skip_pack_and_min ) //basic::options::OptionKeys::PDB_diagnostic::skip_pack_and_min

static THREAD_LOCAL basic::Tracer TR( "PDB_diagnostic" );

///local mover for testing purposes
class PDBDiagnosticMover : public protocols::moves::Mover {
public:
	PDBDiagnosticMover()
	{
		TR << "PDBDiagnosticMover ctor" << std::endl;

	}

	virtual ~PDBDiagnosticMover(){};

	///@details this function sums statistics: for each is_?? function in ResidueType, it accumulates how many residues of the Pose have that quality, and dumps that data to the scorefile via the Job's interface.
	void residue_type_statistics( core::pose::Pose const & pose, protocols::jd2::JobOP job_me, core::Size const nres ){
		//I wonder if it's faster to pass nres in the function signature or just recalculate it?  I wonder if I've already spent more time thinking about it than the computer would spend, over all of human history, on the slower option?  I wonder if the wasted disk space for this comment outweighs having thought about it?

		//counter variable for each thing to track
		//core::Size is_??(0);
		core::Size is_polymer(0);
		core::Size is_sidechain_thiol(0);
		core::Size is_disulfide_bonded(0);
		core::Size is_sidechain_amine(0);
		core::Size is_protein(0);
		core::Size is_alpha_aa(0);
		core::Size is_beta_aa(0);
		core::Size is_gamma_aa(0);
		core::Size is_sri(0);
		core::Size is_triazolemer(0);
		core::Size is_d_aa(0);
		core::Size is_l_aa(0);
		core::Size is_achiral_backbone(0);
		core::Size is_DNA(0);
		core::Size is_RNA(0);
		core::Size is_coarse(0);
		core::Size is_NA(0);
		core::Size is_peptoid(0);
		core::Size is_carbohydrate(0);
		core::Size is_ligand(0);
		core::Size is_lipid(0);
		core::Size is_metal(0);
		core::Size is_metalbinding(0);
		core::Size is_membrane(0);
		core::Size is_surface(0);
		core::Size has_sc_orbitals(0);
		core::Size is_polar(0);
		core::Size is_charged(0);
		core::Size is_aromatic(0);
		core::Size is_cyclic(0);
		core::Size is_terminus(0);
		core::Size is_lower_terminus(0);
		core::Size is_upper_terminus(0);
		core::Size is_branch_point(0);
		core::Size is_branch_lower_terminus(0);
		core::Size is_acetylated_nterminus(0);
		core::Size is_methylated_cterminus(0);
		core::Size is_virtual_residue(0);
		core::Size is_adduct(0);

		//loop over pose residues, accumulating properties
		for ( core::Size i(1); i<=nres; ++i ) {
			//if(pose.residue_type(i).??) ++??;
			if ( pose.residue_type(i).is_polymer() ) ++is_polymer;
			if ( pose.residue_type(i).is_sidechain_thiol() ) ++is_sidechain_thiol;
			if ( pose.residue_type(i).is_disulfide_bonded() ) ++is_disulfide_bonded;
			if ( pose.residue_type(i).is_sidechain_amine() ) ++is_sidechain_amine;
			if ( pose.residue_type(i).is_protein() ) ++is_protein;
			if ( pose.residue_type(i).is_alpha_aa() ) ++is_alpha_aa;
			if ( pose.residue_type(i).is_beta_aa() ) ++is_beta_aa;
			if ( pose.residue_type(i).is_gamma_aa() ) ++is_gamma_aa;
			if ( pose.residue_type(i).is_sri() ) ++is_sri;
			if ( pose.residue_type(i).is_triazolemer() ) ++is_triazolemer;
			if ( pose.residue_type(i).is_d_aa() ) ++is_d_aa;
			if ( pose.residue_type(i).is_l_aa() ) ++is_l_aa;
			if ( pose.residue_type(i).is_achiral_backbone() ) ++is_achiral_backbone;
			if ( pose.residue_type(i).is_DNA() ) ++is_DNA;
			if ( pose.residue_type(i).is_RNA() ) ++is_RNA;
			if ( pose.residue_type(i).is_coarse() ) ++is_coarse;
			if ( pose.residue_type(i).is_NA() ) ++is_NA;
			if ( pose.residue_type(i).is_peptoid() ) ++is_peptoid;
			if ( pose.residue_type(i).is_carbohydrate() ) ++is_carbohydrate;
			if ( pose.residue_type(i).is_ligand() ) ++is_ligand;
			if ( pose.residue_type(i).is_lipid() ) ++is_lipid;
			if ( pose.residue_type(i).is_metal() ) ++is_metal;
			if ( pose.residue_type(i).is_metalbinding() ) ++is_metalbinding;
			if ( pose.residue_type(i).is_membrane() ) ++is_membrane;
			if ( pose.residue_type(i).is_surface() ) ++is_surface;
			if ( pose.residue_type(i).has_sc_orbitals() ) ++has_sc_orbitals;
			if ( pose.residue_type(i).is_polar() ) ++is_polar;
			if ( pose.residue_type(i).is_charged() ) ++is_charged;
			if ( pose.residue_type(i).is_aromatic() ) ++is_aromatic;
			if ( pose.residue_type(i).is_cyclic() ) ++is_cyclic;
			if ( pose.residue_type(i).is_terminus() ) ++is_terminus;
			if ( pose.residue_type(i).is_lower_terminus() ) ++is_lower_terminus;
			if ( pose.residue_type(i).is_upper_terminus() ) ++is_upper_terminus;
			if ( pose.residue_type(i).is_branch_point() ) ++is_branch_point;
			if ( pose.residue_type(i).is_branch_lower_terminus() ) ++is_branch_lower_terminus;
			if ( pose.residue_type(i).is_acetylated_nterminus() ) ++is_acetylated_nterminus;
			if ( pose.residue_type(i).is_methylated_cterminus() ) ++is_methylated_cterminus;
			if ( pose.residue_type(i).is_virtual_residue() ) ++is_virtual_residue;
			if ( pose.residue_type(i).is_adduct() ) ++is_adduct;

		}

		//dump all these counts into scorefile
		//job_me->add_string_real_pair("??", ??);

		job_me->add_string_real_pair("nres", nres);

		job_me->add_string_real_pair("is_polymer", is_polymer);
		job_me->add_string_real_pair("is_sidechain_thiol", is_sidechain_thiol);
		job_me->add_string_real_pair("is_disulfide_bonded", is_disulfide_bonded);
		job_me->add_string_real_pair("is_sidechain_amine", is_sidechain_amine);
		job_me->add_string_real_pair("is_protein", is_protein);
		job_me->add_string_real_pair("is_alpha_aa", is_alpha_aa);
		job_me->add_string_real_pair("is_beta_aa", is_beta_aa);
		job_me->add_string_real_pair("is_gamma_aa", is_gamma_aa);
		job_me->add_string_real_pair("is_sri", is_sri);
		job_me->add_string_real_pair("is_triazolemer", is_triazolemer);
		job_me->add_string_real_pair("is_d_aa", is_d_aa);
		job_me->add_string_real_pair("is_l_aa", is_l_aa);
		job_me->add_string_real_pair("is_achiral_backbone", is_achiral_backbone);
		job_me->add_string_real_pair("is_DNA", is_DNA);
		job_me->add_string_real_pair("is_RNA", is_RNA);
		job_me->add_string_real_pair("is_coarse", is_coarse);
		job_me->add_string_real_pair("is_NA", is_NA);
		job_me->add_string_real_pair("is_peptoid", is_peptoid);
		job_me->add_string_real_pair("is_carbohydrate", is_carbohydrate);
		job_me->add_string_real_pair("is_ligand", is_ligand);
		job_me->add_string_real_pair("is_lipid", is_lipid);
		job_me->add_string_real_pair("is_metal", is_metal);
		job_me->add_string_real_pair("is_metalbinding", is_metalbinding);
		job_me->add_string_real_pair("is_membrane", is_membrane);
		job_me->add_string_real_pair("is_surface", is_surface);
		job_me->add_string_real_pair("has_sc_orbitals", has_sc_orbitals);
		job_me->add_string_real_pair("is_polar", is_polar);
		job_me->add_string_real_pair("is_charged", is_charged);
		job_me->add_string_real_pair("is_aromatic", is_aromatic);
		job_me->add_string_real_pair("is_cyclic", is_cyclic);
		job_me->add_string_real_pair("is_terminus", is_terminus);
		job_me->add_string_real_pair("is_lower_terminus", is_lower_terminus);
		job_me->add_string_real_pair("is_upper_terminus", is_upper_terminus);
		job_me->add_string_real_pair("is_branch_point", is_branch_point);
		job_me->add_string_real_pair("is_branch_lower_terminus", is_branch_lower_terminus);
		job_me->add_string_real_pair("is_acetylated_nterminus", is_acetylated_nterminus);
		job_me->add_string_real_pair("is_methylated_cterminus", is_methylated_cterminus);
		job_me->add_string_real_pair("is_virtual_residue", is_virtual_residue);
		job_me->add_string_real_pair("is_adduct", is_adduct);

		return;
	}


	virtual
	void
	apply( core::pose::Pose & pose ){

		//Diagnostic plan:
		//1) try reading PDB in (done by JD, already complete at this step)
		//2) calculate statistics based on the constellation of ResidueTypes in the Pose
		//3) try scoring
		//4) try packing
		//5) try minimizing

		using protocols::jd2::JobDistributor;
		protocols::jd2::JobOP job_me( JobDistributor::get_instance()->current_job() );
		std::string this_pdb_name(JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );
		core::Size const nres(pose.total_residue());
		TR << this_pdb_name << " nres " << nres << std::endl;

		//step 2: accumulate residue type stats
		residue_type_statistics(pose, job_me, nres);

		//if option, cut and run! Useful for not-mega-clusters where packing a virus/ribosome will crush memory.
		if ( basic::options::option[ basic::options::OptionKeys::PDB_diagnostic::reading_only ].value() ) {
			return;
		}

		//create a ScoreFunction from commandline options
		core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();

		//step 3: scoring
		using protocols::simple_moves::ScoreMoverOP;
		using protocols::simple_moves::ScoreMover;
		ScoreMoverOP score_mover(new ScoreMover(score_fxn));
		score_mover->set_verbose(false);
		score_mover->apply(pose);

		if ( !basic::options::option[ basic::options::OptionKeys::PDB_diagnostic::skip_pack_and_min ].value() ) {

			//packing crashes on zero nres
			if ( nres > 0 ) {

				//step 4: packing
				//note the main function forces -repack_only
				//set up SF and TF defaults (manually for some reason)
				//create a task factory: this will create a new PackerTask for each input pose
				using core::pack::task::operation::TaskOperationCOP;
				core::pack::task::TaskFactoryOP main_task_factory( new core::pack::task::TaskFactory );
				main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
				using protocols::simple_moves::PackRotamersMoverOP;
				using protocols::simple_moves::PackRotamersMover;
				PackRotamersMoverOP pack_rotamers(new protocols::simple_moves::PackRotamersMover());
				pack_rotamers->task_factory( main_task_factory );
				pack_rotamers->score_function( score_fxn );
				pack_rotamers->apply(pose);
			} //skip all this if no residues!

			//step 5: minimizing
			using protocols::simple_moves::MinMoverOP;
			using protocols::simple_moves::MinMover;
			MinMoverOP min_mover(new protocols::simple_moves::MinMover());
			min_mover->score_function( score_fxn );
			min_mover->apply(pose);
		}//skip_pack_and_min

		return;
	}

	virtual
	std::string
	get_name() const {
		return "PDBDiagnosticMover";
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return protocols::moves::MoverOP( new PDBDiagnosticMover );
	}

	virtual
	bool
	reinitialize_for_each_job() const { return true; }

	virtual
	bool
	reinitialize_for_new_input() const { return true; }

private:

};

typedef utility::pointer::shared_ptr< PDBDiagnosticMover > PDBDiagnosticMoverOP;


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		NEW_OPT( PDB_diagnostic::reading_only, "Skip scoring/packing/minimization: only test read-in functions", false );
		NEW_OPT( PDB_diagnostic::skip_pack_and_min, "Skip packing/minimization: only test read-in and scoring", false );

		devel::init(argc, argv);

		//force these options: don't leak memory, don't pack on read-in, and never ever design
		basic::options::option[ basic::options::OptionKeys::jd2::delete_old_poses ].value(true);
		basic::options::option[ basic::options::OptionKeys::packing::pack_missing_sidechains ].value(false);
		basic::options::option[ basic::options::OptionKeys::packing::repack_only ].value(true);

		//force no-output, otherwise we are wasting time
		bool const score_only(basic::options::option[ basic::options::OptionKeys::out::file::score_only ].user());
		if ( !score_only ) { //if not, do nothing
			//do nothing
			TR << "PDB_diagnostic requires -score_only (or you'll write god knows how many PDBs out!)" << std::endl;
		} else { //do work

			PDBDiagnosticMoverOP test_mover(new PDBDiagnosticMover);

			protocols::jd2::JobDistributor::get_instance()->go(test_mover);
		} //if score_only
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}


	TR << "*****************COMPLETED******************" << std::endl;
	return 0;
}
