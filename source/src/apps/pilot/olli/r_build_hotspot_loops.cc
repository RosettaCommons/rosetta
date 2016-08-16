// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief demo program for implementing loop relax + FA relax
/// @author Srivatsan Raman
/// @author James Thompson
/// @author Mike Tyka
/// @author Daniel J. Mandell

//unit headers
#include <protocols/loop_build/LoopBuild.hh>

// include these first for building on Visual Studio

#include <protocols/moves/MoverStatus.hh>
#include <core/types.hh>

#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/util.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/Loops.hh>
#include <utility/excn/Exceptions.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

#include <core/kinematics/FoldTree.hh>


#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <protocols/protein_interface_design/movers/SetupHotspotConstraintsMover.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/simple_filters/ConstraintScoreCutoffFilter.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/BackboneStubConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <core/fragment/FragSet.hh>
#include <utility/exit.hh>

#include <core/scoring/electron_density/util.hh>

#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/MoverContainer.hh>
// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

//silly using/typedef
#include <core/id/AtomID_Map.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/chemical/AA.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/pack_rotamers.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <protocols/comparative_modeling/LoopRelaxMover.hh>

#include <core/import_pose/import_pose.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
// Utility headers
#include <basic/options/option_macros.hh>

#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <boost/foreach.hpp>
#include <devel/init.hh>


#include <protocols/filters/Filter.hh>

using basic::T;
using basic::Error;
using basic::Warning;


static THREAD_LOCAL basic::Tracer tr( "main" );

OPT_KEY( File, hotspots )
OPT_KEY( File, filter_hotspots1 )
OPT_KEY( File, filter_hotspots2 )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( hotspots, "load hotspots as constraints","" );
	NEW_OPT( filter_hotspots1, "filter1", "" );
	NEW_OPT( filter_hotspots2, "filter2", "" );
}

using namespace protocols;
using namespace hotspot_hashing;
using namespace protocols::protein_interface_design;
using namespace protocols::simple_filters;
/// @brief utility function to add hotspot bbcst's to a pose
class ApplyFilterMover : public  protocols::moves::Mover {
private:
	utility::vector1<protocols::filters::FilterOP> filters_;
};

class SetupHotspotMover : public protocols::moves::Mover {

public:
	SetupHotspotMover( HotspotStubSetOP hotspot_stub_set )
	: Mover ("SetupHotspotMover" ),
		hotspot_stub_set_( hotspot_stub_set ),
		resfile_( "NONE" )
	{}
	void apply( core::pose::Pose & pose );
	core::Size generate_csts( core::pose::Pose const& pose,  core::scoring::constraints::ConstraintCOPs& constraints );
	std::string get_name() const { return "SetupHotspotMover"; };
	void set_resfile( std::string const& setting ) {
		resfile_ = setting;
	}

private:
	HotspotStubSetOP hotspot_stub_set_;
	std::string resfile_;
};


core::Size SetupHotspotMover::generate_csts( core::pose::Pose const& pose,  core::scoring::constraints::ConstraintCOPs& constraints ) {
core::id::AtomID fixed_atom(1, 1);
		//	core::pack::task::PackerTaskCOP const packer_task,
	core::Real CB_force_constant = 0.5;
	core::Real worst_allowed_stub_bonus = 0;
	bool apply_self_energies = false;
	bool apply_ambiguous_constraints = true;

	// Take in a scaffold pose (with PackerTask, for its DesignMap), and a set of stubs.
	// Each repacked residue will get one "AmbiguousConstraint".
	// This AmbiguousConstraint will contain a series of BackboneStubConstraints (one for each valid stub)
	runtime_assert( CB_force_constant > -1E-6 ); // these can't be negative
	runtime_assert( worst_allowed_stub_bonus < 1E-6 ); // these can't be positive

	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
	if ( resfile_ =="NONE" && basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
		core::pack::task::parse_resfile(pose, *task);
	} else if ( resfile_ != "NONE" ) {
		core::pack::task::parse_resfile(pose, *task, resfile_ );
	}
	// *****associate the stub set with the unbound pose
	// *****hotspot_stub_set_->pair_with_scaffold( pose, partner );

	protocols::filters::FilterCOP true_filter( new protocols::filters::TrueFilter );
	for ( core::Size resnum=55; resnum <= 71; ++resnum ) {
		//		if ( task->pack_residue(resnum) )	{
		hotspot_stub_set_->pair_with_scaffold( pose, pose.chain( resnum ), true_filter );
		break;
		//}
	}

	tr.Info << "Making hotspot constraints..." << std::endl;
	Size scaffold_seqpos(0);
	Size ct_cst( 0 );
	for ( core::Size resnum=55; resnum <= 71; ++resnum ) {

		// Check that this position is allowed to be used for stub constraints
		//		if ( ! task->pack_residue(resnum) ) continue;

		// sets the index used by the hotspot for its associated scaffold
		scaffold_seqpos = resnum - pose.conformation().chain_begin( pose.chain( resnum ) );

		// Start the vector which will become a single AmbiguousConstraint, if apply_ambiguous_constraints is true
		utility::vector1< core::scoring::constraints::ConstraintCOP > ambig_csts;
		// Loop over all allowed AAs at this position
		std::list< core::chemical::ResidueTypeCOP > allowed_aas = task->residue_task( resnum ).allowed_residue_types();
		for (std::list< core::chemical::ResidueTypeCOP >::const_iterator restype = allowed_aas.begin();
				 restype != allowed_aas.end(); ++restype) {

			// Loop over all stubs with this restype
			HotspotStubSet::Hotspots res_stub_set( hotspot_stub_set_->retrieve( (*restype )->name3() ) );
			for (std::multimap<core::Real,HotspotStubOP >::iterator hs_stub = res_stub_set.begin();
					 hs_stub != res_stub_set.end(); ++hs_stub) {

				// prevent Gly/Pro constraints
				if ( (hs_stub->second->residue()->aa() == core::chemical::aa_gly) || (hs_stub->second->residue()->aa() == core::chemical::aa_pro && !basic::options::option[basic::options::OptionKeys::hotspot::allow_proline] ) ) {
					tr.Info << "ERROR - Gly/Pro stubs cannot be used for constraints." << std::endl;
					continue;
				}

				// prevent Gly/Pro constraints
				if ( (pose.residue(resnum).aa() == core::chemical::aa_gly) || (pose.residue(resnum).aa() == core::chemical::aa_pro && !basic::options::option[basic::options::OptionKeys::hotspot::allow_proline]) ) {
									tr.Debug << "ERROR - Position " << resnum << " is currently Gly/Pro and cannot be used for stub constraints." << std::endl;
									continue;
				}

				core::Real stub_bonus_value = hs_stub->second->bonus_value();
				if ( stub_bonus_value < worst_allowed_stub_bonus ) {
					hs_stub->second->set_scaffold_status( resnum, protocols::hotspot_hashing::accept );
					//tr.Info << " SuccSelfEnergy=" << stub_bonus_value << std::endl;
					// ****** accept the pairing -- do we really want this? better to just reject, since bb fit doesn't necessarily mean good pair
					// ****** hs_stub->scaffold_status( resnum, accept );

					// Build a BackboneStubConstraint from this stub
					if ( apply_ambiguous_constraints ) {
						// Push it onto ambig_csts for this residue
						ct_cst++;
						ambig_csts.push_back( new core::scoring::constraints::BackboneStubConstraint( pose, resnum, fixed_atom, *(hs_stub->second->residue()), stub_bonus_value, CB_force_constant ) );
					} else {
						// Apply it directly
						constraints.push_back( new core::scoring::constraints::BackboneStubConstraint( pose, resnum, fixed_atom, *(hs_stub->second->residue()), stub_bonus_value, CB_force_constant ) );
					}
				}	else hs_stub->second->set_scaffold_status( resnum, protocols::hotspot_hashing::reject );
				//else tr.Info << " FailSelfEnergy=" << stub_bonus_value << std::endl;
				// ****** reject the pairing
				// ******else hs_stub->scaffold_status( resnum, reject );
			}
		}

		// Finally, add the constraint corresponding to this resnum to the main set
		if ( ( apply_ambiguous_constraints ) && ( ambig_csts.size() > 0 ) )
			constraints.push_back( new core::scoring::constraints::AmbiguousConstraint(ambig_csts) );
	}
	return ct_cst;
}

void SetupHotspotMover::apply( core::pose::Pose & pose ) {

	core::scoring::constraints::ConstraintCOPs constraints_;
	Size ct_cst = generate_csts( pose, constraints_ );
	constraints_ = pose.add_constraints( constraints_ );
	tr.Info << "Applied " << ct_cst << " hotspots in " << constraints_.size() << " constraints to the pose." << std::endl;
	return;
}

ConstraintScoreCutoffFilterOP setup_hotspot_filter( core::pose::Pose const& pose, std::string filename ) {
	tr.Info << "setup filter for " << filename << " read hotspot set..." << std::endl;
	HotspotStubSetOP hotspot_set = new HotspotStubSet;
	hotspot_set->read_data( filename );
	SetupHotspotMover hspmover( hotspot_set );
	core::scoring::constraints::ConstraintCOPs constraints;
	core::Size ncst = hspmover.generate_csts( pose, constraints );
	ConstraintScoreCutoffFilterOP cst_filter = new ConstraintScoreCutoffFilter();
	cst_filter->set_constraints( constraints );
	tr.Info << "created " << ncst << " hotspot constraints from file " << filename << " for filtering... "<< std::endl;
	return cst_filter;
}

int
LoopBuild_main() {

	basic::Tracer tr( "protocols.loop_build.LoopBuild" );

	using namespace basic::options;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::chemical;
	using namespace core::id;
	using namespace protocols;
	using namespace protocols::loops;

	std::string remodel            ( option[ OptionKeys::loops::remodel ]() );
	std::string const intermedrelax( option[ OptionKeys::loops::intermedrelax ]() );
	std::string const refine       ( option[ OptionKeys::loops::refine ]() );
	std::string const relax        ( option[ OptionKeys::loops::relax ]() );
	bool const keep_time      ( option[ OptionKeys::loops::timer ]() );

	tr << "==== Loop protocol: ================================================="
		<< std::endl;
	tr << " remodel        " << remodel        << std::endl;
	tr << " intermedrelax  " << intermedrelax  << std::endl;
	tr << " refine         " << refine         << std::endl;
	tr << " relax          " << relax          << std::endl;

	core::pose::Pose start_pose, pose, init_pose_obj;
	core::chemical::ResidueTypeSetCAP rsd_set;

	// DJM: must first load the structure in fullatom mode. Otherwise
	// disulfide-bonded cysteines will not have their residue types set
	// correctly, which causes crashes upon switching to full-atom and repacking.
	if ( option[ OptionKeys::in::file::fullatom ]() ) {
		// if full-atom load starting structure as full-atom to recover
		// sidechains later
		rsd_set
			= ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		core::import_pose::pose_from_file(
			start_pose, *rsd_set, option[ OptionKeys::loops::input_pdb ]().name()
		);
	} else { // no full-atom, create a centroid PDB
		rsd_set
			= ChemicalManager::get_instance()->residue_type_set( "centroid" );
		core::import_pose::pose_from_file(
			start_pose, *rsd_set, option[ OptionKeys::loops::input_pdb ]().name()
		);
	}

	bool psipred_ss2_ok = loops::set_secstruct_from_psipred_ss2( start_pose );
	if ( !psipred_ss2_ok ) {
		std::string dssp_name( option[ OptionKeys::in::file::dssp ]().name() );
		bool dssp_ok = loops::set_secstruct_from_dssp(start_pose, dssp_name);
		if ( !dssp_ok ) {
			core::pose::set_ss_from_phipsi( start_pose );
		}
	}

	// symmetrize start pose & loopfile
	if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::simple_moves::symmetry::SetupForSymmetryMover pre_mover;
			pre_mover.apply( start_pose );
	}

	// bit of a hack for looprelax-into-density
	// set pose for density scoring if a map was input
	// (potentially) dock map into density -- do this here so we only need to
	// dock the pose once
	if ( option[ OptionKeys::edensity::mapfile ].user() ) {
		protocols::electron_density::SetupForDensityScoringMover pre_mover;
		pre_mover.apply( start_pose );
	}

	core::pose::Pose native_pose;
	if ( option[ OptionKeys::in::file::native ].user() ) {
		core::import_pose::pose_from_file(
			native_pose, option[ OptionKeys::in::file::native ]()
		);
		core::pose::set_ss_from_phipsi( native_pose );
	} else	{
		native_pose = start_pose;
	}

	// fragment initialization
	// is there any way to clean this up? This logic is very convoluted.
	utility::vector1< core::fragment::FragSetOP > frag_libs;
	if ( remodel == "perturb_ccd" || remodel == "quick_ccd" ||
			remodel == "quick_ccd_moves" || remodel == "old_loop_relax" ||
			remodel == "sdwindow" ||
			option[ OptionKeys::loops::build_initial ].value() ||
	   	( option[ OptionKeys::loops::frag_files ].user()
	      	&& (refine == "refine_ccd" || intermedrelax != "no" || relax != "no")
			)
	) {
		// these protocols optionally take a fragment set .. only load if
		// specified
		loops::read_loop_fragments( frag_libs );
	}

	evaluation::MetaPoseEvaluatorOP evaluator = new evaluation::MetaPoseEvaluator;
	evaluation::EvaluatorFactory::get_instance()->add_all_evaluators(*evaluator);
	evaluator->add_evaluation(
		new simple_filters::SelectRmsdEvaluator( native_pose, "_native" )
	);

	// load loopfile
	protocols::loops::LoopsOP my_loops = new protocols::loops::Loops( true );

	comparative_modeling::LoopRelaxMoverOP mover = new comparative_modeling::LoopRelaxMover;
	mover->frag_libs( frag_libs );
	mover->loops( my_loops );
	mover->relax( relax );
	mover->refine( refine );
	mover->remodel( remodel );
	mover->intermedrelax( intermedrelax );

	// add density wts from cmd line to looprelax scorefunctions
	if ( option[ OptionKeys::edensity::mapfile ].user() ) {
		ScoreFunctionOP lr_cen_scorefxn = loops::get_cen_scorefxn();
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *lr_cen_scorefxn );
		ScoreFunctionOP lr_fa_scorefxn = loops::get_fa_scorefxn();
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *lr_fa_scorefxn );
		mover->scorefxns( lr_cen_scorefxn, lr_fa_scorefxn );
	}
	moves::MoverOP final_mover = mover;
	if ( option[ OptionKeys::hotspots ].user() ) {

		tr.Info << "read hotspot set..." << std::endl;
		HotspotStubSetOP hotspot_set = new HotspotStubSet;
		hotspot_set->read_data( option[ OptionKeys::hotspots ]() );
		moves::MoverOP hotspot_mover =	new SetupHotspotMover( hotspot_set );
		moves::SequenceMoverOP seq_mover = new moves::SequenceMover;
		seq_mover->add_mover( hotspot_mover );
		seq_mover->add_mover( mover );
		final_mover = seq_mover;
		if ( option[ OptionKeys::filter_hotspots1 ].user() ) {
			ConstraintScoreCutoffFilterOP cst_filter( setup_hotspot_filter( start_pose, option[ OptionKeys::filter_hotspots1 ]() ) );
			cst_filter->set_user_defined_name( "hotspot1" );
			seq_mover->add_mover( cst_filter );
		}
		if ( option[ OptionKeys::filter_hotspots2 ].user() ) {
			ConstraintScoreCutoffFilterOP cst_filter( setup_hotspot_filter( start_pose, option[ OptionKeys::filter_hotspots2 ]() ) );
			cst_filter->set_user_defined_name( "hotspot2" );
			seq_mover->add_mover( cst_filter );
		}
	}
	tr.Info << "start job-distribution..." << std::endl;
	jd2::JobDistributor::get_instance()->go( final_mover );

	return 0;
} // Looprelax_main

int
main( int argc, char * argv [] )
{
	try{
	register_options();
	devel::init( argc, argv );
	LoopBuild_main();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
