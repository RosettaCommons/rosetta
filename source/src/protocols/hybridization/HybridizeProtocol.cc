// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Yifan Song

#include <core/init/score_function_corrections.hh>

#include <protocols/hybridization/HybridizeProtocolCreator.hh>
#include <protocols/hybridization/HybridizeProtocol.hh>
#include <protocols/hybridization/FoldTreeHybridize.hh>
#include <protocols/hybridization/CartesianHybridize.hh>
#include <protocols/hybridization/TemplateHistory.hh>
#include <protocols/hybridization/util.hh>
#include <protocols/hybridization/DomainAssembly.hh>
#include <protocols/hybridization/DDomainParse.hh>
#include <protocols/hybridization/TMalign.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/util.hh>

#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>

#include <protocols/electron_density/SetupForDensityScoringMover.hh>

// dynamic fragpick
#include <protocols/moves/DsspMover.hh>
#include <core/fragment/picking_old/vall/util.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/util.hh>

#include <core/sequence/util.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// task operation
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

#include <basic/datacache/DataMap.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// utility
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <numeric/random/WeightedSampler.hh>
#include <ObjexxFCL/format.hh>
#include <boost/foreach.hpp>

// evaluation
#include <core/scoring/rms_util.hh>
#include <protocols/comparative_modeling/coord_util.hh>

// option
#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh> // strand pairings
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh> // check pre talaris

//docking
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>

#include <string>

static basic::Tracer TR( "protocols.hybridization.HybridizeProtocol" );
static numeric::random::RandomGenerator RG(541938);

namespace protocols {
namespace hybridization {

//
using namespace core;
using namespace kinematics;
using namespace sequence;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;
using namespace constraints;


/////////////
// creator
std::string
HybridizeProtocolCreator::keyname() const {
	return HybridizeProtocolCreator::mover_name();
}

protocols::moves::MoverOP
HybridizeProtocolCreator::create_mover() const {
	return new HybridizeProtocol;
}

std::string
HybridizeProtocolCreator::mover_name() {
	return "Hybridize";
}



/////////////
// mover
HybridizeProtocol::HybridizeProtocol() :
	template_weights_sum_(0)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	init();

	// initialization may come from command line or from RS
	if (option[cm::hybridize::template_list].user()) {
		read_template_structures( option[cm::hybridize::template_list]() );
	}
}

// sets default options
void
HybridizeProtocol::init() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	stage1_probability_ = option[cm::hybridize::stage1_probability]();
	stage1_increase_cycles_ = option[cm::hybridize::stage1_increase_cycles]();
	stage1_1_cycles_ = option[cm::hybridize::stage1_1_cycles]();
	stage1_2_cycles_ = option[cm::hybridize::stage1_2_cycles]();
	stage1_3_cycles_ = option[cm::hybridize::stage1_3_cycles]();
	stage1_4_cycles_ = option[cm::hybridize::stage1_4_cycles]();
	domain_assembly_ = false;
	add_hetatm_ = false;
	realign_domains_ = option[cm::hybridize::realign_domains]();
	realign_domains_stage2_ = option[cm::hybridize::realign_domains_stage2]();
	add_non_init_chunks_ = option[cm::hybridize::add_non_init_chunks]();
	frag_weight_aligned_ = option[cm::hybridize::frag_weight_aligned]();
	auto_frag_insertion_weight_ = option[cm::hybridize::auto_frag_insertion_weight]();
	max_registry_shift_ = option[cm::hybridize::max_registry_shift]();
	frag_1mer_insertion_weight_ = option[cm::hybridize::frag_1mer_insertion_weight]();
	small_frag_insertion_weight_ = option[cm::hybridize::small_frag_insertion_weight]();
	big_frag_insertion_weight_ = option[cm::hybridize::big_frag_insertion_weight]();
	chunk_insertion_weight_ = 1.;
	hetatm_self_cst_weight_ = 10.;
	hetatm_prot_cst_weight_ = 0.;
	cartfrag_overlap_ = 2;
	seqfrags_only_ = false;
	nofragbias_ = false;
	skip_long_min_ = true;   //fpd  this is no longer necessary and seems to hurt model accuracy
	keep_pose_constraint_ = false;   //fpd PLEASE INITIALIZE NEW VARIABLES

	jump_move_ = false;
	jump_move_repeat_ = 1;

	if (option[cm::hybridize::starting_template].user()) {
		starting_templates_ = option[cm::hybridize::starting_template]();
	}

	stage2_increase_cycles_ = option[cm::hybridize::stage2_increase_cycles]();
	stage25_increase_cycles_ = option[cm::hybridize::stage2min_increase_cycles]();
	no_global_frame_ = option[cm::hybridize::no_global_frame]();
	linmin_only_ = option[cm::hybridize::linmin_only]();

	// default scorefunction
	stage1_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(
		option[cm::hybridize::stage1_weights](), option[cm::hybridize::stage1_patch]() );

	if (!option[mistakes::restore_pre_talaris_2013_behavior] && !option[corrections::score::cenrot]) {
	  ///////////////////////////////////////////////////////////////////////////////////////
	  // restore hb term to score12
	  option[ corrections::score::hb_sp2_chipen ].value( false );
	  option[ corrections::score::hb_fade_energy ].value( false );
	  option[ corrections::score::hbond_measure_sp3acc_BAH_from_hvy ].value( false );
	  option[ corrections::score::hb_sp2_outer_width ].value( 0.33333 );
	  option[ score::hbond_params ].value( "score12_params" );
	}

	stage2_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(
		option[cm::hybridize::stage2_weights](), option[cm::hybridize::stage2_patch]() );

	if (!option[mistakes::restore_pre_talaris_2013_behavior] && !option[corrections::score::cenrot]) {
	  core::scoring::methods::EnergyMethodOptions options2(stage2_scorefxn_->energy_method_options());
	  core::scoring::hbonds::HBondOptions hbopt;
	  hbopt.params_database_tag("score12_params");
	  options2.hbond_options(hbopt);
	  stage2_scorefxn_->set_energy_method_options(options2);

    ///////////////////////////////////////////////////////////////////////////////////////
	  // take talaris2013 back
	  option[ corrections::score::hb_sp2_chipen ].value( true );
	  option[ corrections::score::hb_fade_energy ].value( true );
	  option[ corrections::score::hbond_measure_sp3acc_BAH_from_hvy ].value( true );
	  option[ corrections::score::hb_sp2_outer_width ].value( 0.357 );
	  option[ score::hbond_params ].value( "sp2_elec_params" );
	}

	fa_scorefxn_ = core::scoring::getScoreFunction();

	if (!option[mistakes::restore_pre_talaris_2013_behavior] && !option[corrections::score::cenrot]) {
		core::scoring::methods::EnergyMethodOptions optionsfa(fa_scorefxn_->energy_method_options());
		core::scoring::hbonds::HBondOptions hboptfa;
		hboptfa.params_database_tag("sp2_elec_params");
		optionsfa.hbond_options(hboptfa);
		fa_scorefxn_->set_energy_method_options(optionsfa);
	}

	core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *fa_scorefxn_ );


	if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
		utility::vector1< std::string > cst_files = option[ OptionKeys::constraints::cst_fa_file ]();
		core::Size choice = core::Size( RG.random_range( 1, cst_files.size() ) );
		fa_cst_fn_ = cst_files[choice];
		TR.Info << "Fullatom Constraint choice: " << fa_cst_fn_ << std::endl;
	}

	batch_relax_ = option[ cm::hybridize::relax ]();
	relax_repeats_ = option[ basic::options::OptionKeys::relax::default_repeats ]();

	// disulfide file
	if (option[ in::fix_disulf ].user()) {
		disulf_file_ = option[ in::fix_disulf ]();
	}

	// read fragments
	if ( option[ in::file::frag9 ].user() ) {
		using namespace core::fragment;
		FragSetOP frags = FragmentIO().read_data( option[ in::file::frag9 ]() );
		fragments_big_.push_back(frags);
	}
	if ( option[ in::file::frag3 ].user() ) {
		using namespace core::fragment;
		FragSetOP frags = FragmentIO().read_data( option[ in::file::frag3 ]() );
		fragments_small_.push_back(frags);
	}

	// native
	if ( option[ in::file::native ].user() ) {
		native_ = new core::pose::Pose;
		if (option[in::file::fullatom]()){
			core::import_pose::pose_from_pdb( *native_, option[ in::file::native ]() );
		}
		else {
			core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
			core::import_pose::pose_from_pdb( *native_, *residue_set, option[ in::file::native ]()  );
		}
	} else if ( option[ evaluation::align_rmsd_target ].user() ) {
		native_ = new core::pose::Pose;
		utility::vector1< std::string > const & align_rmsd_target( option[ evaluation::align_rmsd_target ]() );
		core::import_pose::pose_from_pdb( *native_, align_rmsd_target[1] ); // just use the first one for now
	}

	// strand pairings
 	pairings_file_ = option[jumps::pairing_file]();
	if ( option[jumps::sheets].user() ) {
		sheets_ = option[jumps::sheets]();
	} else {
		random_sheets_ = option[jumps::random_sheets]();
	}
	filter_templates_ = option[jumps::filter_templates]();
}

void
HybridizeProtocol::check_and_create_fragments( core::pose::Pose & pose ) {
	if (fragments_big_.size() > 0 && fragments_small_.size() > 0) return;

	core::Size fragbiglen=9;

	if (!fragments_big_.size()) {
		// number of residues
		core::Size nres_tgt = pose.total_residue();
		core::conformation::symmetry::SymmetryInfoCOP symm_info;
		if ( core::pose::symmetry::is_symmetric(pose) ) {
			core::conformation::symmetry::SymmetricConformation & SymmConf (
				dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
			symm_info = SymmConf.Symmetry_Info();
			nres_tgt = symm_info->num_independent_residues();
		}
		if (pose.residue(nres_tgt).aa() == core::chemical::aa_vrt) nres_tgt--;
		while (!pose.residue(nres_tgt).is_protein()) nres_tgt--;

		fragbiglen = std::min( fragbiglen, nres_tgt );

		core::fragment::FragSetOP frags = new core::fragment::ConstantLengthFragSet(  );

		// sequence
		std::string tgt_seq = pose.sequence();
		std::string tgt_ss(nres_tgt, '0');

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		if (option[ OptionKeys::in::file::psipred_ss2 ].user()) {
			utility::vector1< char > psipred = read_psipred_ss2_file( pose );
			for ( core::Size j=1; j<=nres_tgt; ++j ) {
				tgt_ss[j-1] = psipred[j];
			}
		} else {
			// templates vote on secstruct
			for (core::Size i=1; i<=templates_.size(); ++i) {
				for (core::Size j=1; j<=templates_[i]->total_residue(); ++j ) {
					if (!templates_[i]->residue(j).is_protein()) continue;
					core::Size tgt_pos = templates_[i]->pdb_info()->number(j);

					runtime_assert( tgt_pos<=nres_tgt );
					char tgt_ss_j = templates_[i]->secstruct(j);

					if (tgt_ss[tgt_pos-1] == '0') {
						tgt_ss[tgt_pos-1] = tgt_ss_j;
					} else if (tgt_ss[tgt_pos-1] != tgt_ss_j) {
						tgt_ss[tgt_pos-1] = 'D'; // templates disagree
					}
				}
			}
			for ( core::Size j=1; j<=nres_tgt; ++j ) {
				if (tgt_ss[j-1] == '0') tgt_ss[j-1] = 'D';
			}
		}

		// pick from vall based on template SS + target sequence
		for ( core::Size j=1; j<=nres_tgt-fragbiglen+1; ++j ) {
			core::fragment::FrameOP frame = new core::fragment::Frame( j, fragbiglen );
			frame->add_fragment(
				core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa( tgt_ss.substr( j-1, fragbiglen ), tgt_seq.substr( j-1, fragbiglen ), 25, true, core::fragment::IndependentBBTorsionSRFD() ) );
			frags->add( frame );
		}
		fragments_big_.push_back( frags );
	}
	if (!fragments_small_.size()) {
		core::fragment::FragSetOP frags = new core::fragment::ConstantLengthFragSet( 3 );

		// make them from big fragments
		core::fragment::chop_fragments( *fragments_big_[1], *frags );
		fragments_small_.push_back( frags );
	}
}


void
HybridizeProtocol::initialize_and_sample_loops(
		core::pose::Pose &pose,
		core::pose::PoseOP chosen_templ,
		protocols::loops::Loops template_contigs_icluster,
		core::scoring::ScoreFunctionOP scorefxn)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// xyz copy starting model
	for (Size i=1; i<=chosen_templ->total_residue(); ++i)
		for (Size j=1; j<=chosen_templ->residue(i).natoms(); ++j) {
			core::id::AtomID src(j,i), tgt(j, chosen_templ->pdb_info()->number(i));
			pose.set_xyz( tgt, chosen_templ->xyz( src ) );
		}

	// make loops as inverse of template_contigs_icluster
	TR << "CONTIGS" << std::endl << template_contigs_icluster << std::endl;
	//core::Size ncontigs = template_contigs_icluster.size();
	core::Size nres_tgt = pose.total_residue();

	//symmetry
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres_tgt = symm_info->num_independent_residues();
	}

	if (pose.residue(nres_tgt).aa() == core::chemical::aa_vrt) nres_tgt--;

	protocols::loops::LoopsOP loops = new protocols::loops::Loops;
	utility::vector1< bool > templ_coverage(nres_tgt, false);

	for (Size i=1; i<=chosen_templ->total_residue(); ++i) {
		core::Size cres = chosen_templ->pdb_info()->number(i);
		templ_coverage[cres] = true;
	}

	// remove 1-3 residue "segments"
	for (Size i=1; i<=nres_tgt-2; ++i) {
		if (!templ_coverage[i] && templ_coverage[i+1] && !templ_coverage[i+2]) {
			templ_coverage[i+1]=false;
		} else if(i<=nres_tgt-3 && !templ_coverage[i] && templ_coverage[i+1] && templ_coverage[i+2] && !templ_coverage[i+3]) {
			templ_coverage[i+1]=false;
			templ_coverage[i+2]=false;
		} else if(i<=nres_tgt-4 &&
		          !templ_coverage[i] && templ_coverage[i+1] && templ_coverage[i+2] && templ_coverage[i+3] && !templ_coverage[i+4]) {
			templ_coverage[i+1]=false;
			templ_coverage[i+2]=false;
			templ_coverage[i+3]=false;
		}
	}
	// make loopfile
	bool inloop=!templ_coverage[1];
	core::Size loopstart=1, loopstop;
	for (Size i=2; i<=nres_tgt; ++i) {
		if (templ_coverage[i] && inloop && allowed_to_move_[i]==true) {
		//if (templ_coverage[i] && inloop) {
			// end loop
			inloop = false;
			loopstop = i;
			if (loopstop < loopstart + 2) {
				if (loopstart>1) loopstart--;
				if (loopstop<nres_tgt) loopstop++;
			}
			loops->add_loop( loopstart,loopstop );
		} else if (!templ_coverage[i] && !inloop) {
			// start loop
			inloop = true;
			loopstart = i-1;
		}
	}
	if (inloop) {
		// end loop
		loopstop = nres_tgt;
		while (loopstop < loopstart + 2) { loopstart--; }
		loops->add_loop( loopstart,loopstop );
	}
	TR << "LOOPS" << std::endl << *loops << std::endl;

	// now do insertions
	if (loops->size() != 0) {
		// set foldtree + variants
		loops->auto_choose_cutpoints(pose);
		core::kinematics::FoldTree f_in=pose.fold_tree(), f_new;
		protocols::loops::fold_tree_from_loops( pose, *loops, f_new);
		pose.fold_tree( f_new );
		protocols::loops::add_cutpoint_variants( pose );

		// set movemap
		core::kinematics::MoveMapOP mm_loop = new core::kinematics::MoveMap();
		for( protocols::loops::Loops::const_iterator it=loops->begin(), it_end=loops->end(); it!=it_end; ++it ) {
			for( Size i=it->start(); i<=it->stop(); ++i ) {
				mm_loop->set_bb(i, true);
				mm_loop->set_chi(i, true); // chi of loop residues
			}
		}

		// setup fragment movers
		protocols::simple_moves::ClassicFragmentMoverOP frag3mover, frag9mover;
		core::fragment::FragSetOP frags_small = fragments_small_[RG.random_range(1,fragments_small_.size())];
		core::fragment::FragSetOP frags_big = fragments_big_[RG.random_range(1,fragments_big_.size())];
		TR.Info << "FRAGMENTS small max length: " << frags_small->max_frag_length() << std::endl;
		TR.Info << "FRAGMENTS big max length: " << frags_big->max_frag_length() << std::endl;
		frag3mover = new protocols::simple_moves::ClassicFragmentMover( frags_small, mm_loop );
		frag3mover->set_check_ss( false ); frag3mover->enable_end_bias_check( false );
		frag9mover = new protocols::simple_moves::ClassicFragmentMover( frags_big, mm_loop );
		frag9mover->set_check_ss( false ); frag9mover->enable_end_bias_check( false );

		// extend + idealize loops
		for( protocols::loops::Loops::const_iterator it=loops->begin(), it_end=loops->end(); it!=it_end; ++it ) {
			protocols::loops::Loop to_idealize( *it );
			protocols::loops::set_extended_torsions( pose, *it );
		}

		// setup MC
		scorefxn->set_weight( core::scoring::linear_chainbreak, 0.5 );
		(*scorefxn)(pose);
		protocols::moves::MonteCarloOP mc1 = new protocols::moves::MonteCarlo( pose, *scorefxn, 2.0 );

		core::Size neffcycles = (core::Size)(1000*option[cm::hybridize::stage1_increase_cycles]());
		for (Size n=1; n<=neffcycles; ++n) {
			frag9mover->apply( pose ); (*scorefxn)(pose); mc1->boltzmann( pose , "frag9" );
			frag3mover->apply( pose ); (*scorefxn)(pose); mc1->boltzmann( pose , "frag3" );

			if (n%100 == 0) {
				mc1->show_scores();
				mc1->show_counters();
			}
		}
		mc1->recover_low( pose );

		scorefxn->set_weight( core::scoring::linear_chainbreak, 2.0 );
		(*scorefxn)(pose);
		protocols::moves::MonteCarloOP mc2 = new protocols::moves::MonteCarlo( pose, *scorefxn, 2.0 );
		for (Size n=1; n<=neffcycles; ++n) {
			frag9mover->apply( pose ); (*scorefxn)(pose); mc2->boltzmann( pose , "frag9" );
			frag3mover->apply( pose ); (*scorefxn)(pose); mc2->boltzmann( pose , "frag3" );

			if (n%100 == 0) {
				mc2->show_scores();
				mc2->show_counters();
			}
		}

		// restore input ft
		mc2->recover_low( pose );
		protocols::loops::remove_cutpoint_variants( pose );
		pose.fold_tree( f_in );
	}
}



void HybridizeProtocol::add_template(
	std::string template_fn,
	std::string cst_fn,
	std::string symm_file,
	core::Real weight,
	core::Real domain_assembly_weight,
	core::Size cluster_id,
	utility::vector1<core::Size> cst_reses)
{
	core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	core::pose::PoseOP template_pose = new core::pose::Pose();
	core::import_pose::pose_from_pdb( *template_pose, *residue_set, template_fn );

	add_template( template_pose, cst_fn, symm_file, weight, domain_assembly_weight, cluster_id, cst_reses, template_fn );
}


void HybridizeProtocol::add_template(
	core::pose::PoseOP template_pose,
	std::string cst_fn,
	std::string symm_file,
	core::Real weight,
	core::Real domain_assembly_weight,
	core::Size cluster_id,
	utility::vector1<core::Size> cst_reses,
	std::string filename)
{
	// add secondary structure information to the template pose
	core::scoring::dssp::Dssp dssp_obj( *template_pose );
	dssp_obj.insert_ss_into_pose( *template_pose );

	// find ss chunks in template
	protocols::loops::Loops contigs = protocols::loops::extract_continuous_chunks(*template_pose);
	protocols::loops::Loops chunks = protocols::loops::extract_secondary_structure_chunks(*template_pose, "HE", 3, 6, 3, 4);

	if (chunks.num_loop() == 0)
		chunks = contigs;

	TR.Debug << "Chunks from template\n" << chunks << std::endl;
	TR.Debug << "Contigs from template\n" << contigs << std::endl;

	template_fn_.push_back(filename);
	templates_.push_back(template_pose);
	template_cst_fn_.push_back(cst_fn);
	symmdef_files_.push_back(symm_file);
	template_weights_.push_back(weight);
	domain_assembly_weights_.push_back(domain_assembly_weight);
	template_clusterID_.push_back(cluster_id);
	template_chunks_.push_back(chunks);
	template_contigs_.push_back(contigs);
	template_cst_reses_.push_back(cst_reses);
}


///@brief  Old way of parsing hybrid config files; RosettaScripts is now preferred.
void HybridizeProtocol::read_template_structures(utility::file::FileName template_list)
{
	utility::io::izstream f_stream( template_list );
	std::string line;
	while (!f_stream.eof()) {
		getline(f_stream, line);

		if (line.size() == 0) continue;

		std::istringstream str_stream(line);
		std::string template_fn;
		std::string cst_fn;
		core::Size cluster_id(0);
		core::Real weight(1.);
		core::Real domain_assembly_weight(0.);
		if (!str_stream.eof()) {
			str_stream >> template_fn;
			if (template_fn.empty()) continue;
			if (template_fn[0] == '#') continue;

			if (!str_stream.eof()) str_stream >> cst_fn;
			if (!str_stream.eof()) str_stream >> cluster_id;
			if (!str_stream.eof()) str_stream >> weight;

			TR << template_fn << " " << cst_fn << " " << cluster_id << " " << weight << std::endl;

			std::string cst_reses_str;
			utility::vector1<core::Size> cst_reses;
			if ( str_stream >> cst_reses_str ) {
				utility::vector1<std::string> cst_reses_parsed = utility::string_split( cst_reses_str , ',' ) ;
				for (Size i=1; i<= cst_reses_parsed.size(); ++i ) {
					cst_reses.push_back( (core::Size) std::atoi( cst_reses_parsed[i].c_str() ) );
				}
			}
			add_template(template_fn, cst_fn, "", weight,domain_assembly_weight, cluster_id, cst_reses);
		}
	}
	f_stream.close();
}

void HybridizeProtocol::read_template_structures(utility::vector1 < utility::file::FileName > const & template_filenames)
{
	templates_.clear();
	templates_.resize(template_filenames.size());

	core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

	for (core::Size i_ref=1; i_ref<= template_filenames.size(); ++i_ref) {
		templates_[i_ref] = new core::pose::Pose();
		core::import_pose::pose_from_pdb( *(templates_[i_ref]), *residue_set, template_filenames[i_ref].name() );
		core::scoring::dssp::Dssp dssp_obj( *templates_[i_ref] );
		dssp_obj.insert_ss_into_pose( *templates_[i_ref] );
	}
}

void HybridizeProtocol::pick_starting_template(core::Size & initial_template_index,
							core::Size & initial_template_index_icluster,
							utility::vector1 < core::Size > & template_index_icluster,
							utility::vector1 < core::pose::PoseOP > & templates_icluster,
							utility::vector1 < core::Real > & weights_icluster,
							utility::vector1 < protocols::loops::Loops > & template_chunks_icluster,
							utility::vector1 < protocols::loops::Loops > & template_contigs_icluster )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	template_index_icluster.clear();
	templates_icluster.clear();
	weights_icluster.clear();
	template_chunks_icluster.clear();
	template_contigs_icluster.clear();

	if (starting_templates_.size() > 0) {
		numeric::random::WeightedSampler weighted_sampler;
		for (Size i=1; i<=starting_templates_.size(); ++i) {
			weighted_sampler.add_weight(template_weights_[starting_templates_[i]]);
		}
		Size k = weighted_sampler.random_sample(RG);
		initial_template_index = starting_templates_[k];
	}
	else {
		numeric::random::WeightedSampler weighted_sampler;
		weighted_sampler.weights(template_weights_);
		initial_template_index = weighted_sampler.random_sample(RG);
	}
	Size cluster_id = template_clusterID_[initial_template_index];

	for (Size i_template = 1; i_template <= template_clusterID_.size(); ++i_template) {
		if ( cluster_id == template_clusterID_[i_template] ) {
			// add templates from the same cluster
			template_index_icluster.push_back(i_template);
			templates_icluster.push_back(templates_[i_template]);
			weights_icluster.push_back(template_weights_[i_template]);
			template_chunks_icluster.push_back(template_chunks_[i_template]);
			template_contigs_icluster.push_back(template_contigs_[i_template]);

			if (i_template == initial_template_index) {
				initial_template_index_icluster = template_index_icluster.size();
			}
		}
	}
}

void HybridizeProtocol::apply( core::pose::Pose & pose )
{
	using namespace protocols::moves;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose::datacache;
	using namespace core::io::silent;
	using namespace ObjexxFCL::format;

	//save necessary constraint in pose
  core::scoring::constraints::ConstraintSetOP save_pose_constraint_set = new core::scoring::constraints::ConstraintSet() ;
  if ( keep_pose_constraint_ ) {
		save_pose_constraint_set = pose.constraint_set()->clone();
		core::scoring::constraints::remove_nonbb_constraints(pose);
	}

	if (pose.is_fullatom()) {
		protocols::moves::MoverOP tocen = new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );
		tocen->apply(pose);
	}

	// make fragments if we don't have them at this point
	check_and_create_fragments( pose );

	// select random fragments if given variable lengths
	core::fragment::FragSetOP frags_small = fragments_small_[RG.random_range(1,fragments_small_.size())];
	core::fragment::FragSetOP frags_big = fragments_big_[RG.random_range(1,fragments_big_.size())];
	TR.Info << "FRAGMENTS small max length: " << frags_small->max_frag_length() << std::endl;
	TR.Info << "FRAGMENTS big max length: " << frags_big->max_frag_length() << std::endl;

	// starting structure pool
	std::vector < SilentStructOP > post_centroid_structs;
	bool need_more_samples = true;

	// number of residues in asu without VRTs
	core::Size nres_tgt = pose.total_residue();
	core::Size nres_protein_tgt = pose.total_residue();
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres_tgt = symm_info->num_independent_residues();
		nres_protein_tgt = symm_info->num_independent_residues();
	}
	if (pose.residue(nres_tgt).aa() == core::chemical::aa_vrt) nres_tgt--;
	while (!pose.residue(nres_protein_tgt).is_protein()) nres_protein_tgt--;

	core::Real gdtmm = 0.0;
	core::sequence::SequenceAlignmentOP native_aln;
	while(need_more_samples) {
		need_more_samples = false;

		// pick starting template
		core::Size initial_template_index;
		core::Size initial_template_index_icluster; // index in the cluster
		utility::vector1 < core::Size > template_index_icluster; // index look back
		utility::vector1 < core::pose::PoseOP > templates_icluster;
		utility::vector1 < core::Real > weights_icluster;
		utility::vector1 < protocols::loops::Loops > template_chunks_icluster;
		utility::vector1 < protocols::loops::Loops > template_contigs_icluster;
		pick_starting_template(
			initial_template_index, initial_template_index_icluster,
			template_index_icluster, templates_icluster, weights_icluster,
			template_chunks_icluster, template_contigs_icluster );
		TR << "Using initial template: " << I(4,initial_template_index) << " " << template_fn_[initial_template_index] << std::endl;

		// Three opts for template preprocessing:
		//  -> adding hetero residues
		//  -> domain assembly
		//  -> local realignment
		// Currently they are mutually exclusive
		utility::vector1< std::pair< core::Size,core::Size > > hetatms;
		if (add_hetatm_) {
			for ( Size ires=1; ires <= templates_[initial_template_index]->total_residue(); ++ires ) {
				if (templates_[initial_template_index]->pdb_info()->number(ires) > (int)nres_tgt) {
					TR.Debug << "Insert hetero residue: " << templates_[initial_template_index]->residue(ires).name3() << std::endl;
					if ( templates_[initial_template_index]->residue(ires).is_polymer() && !templates_[initial_template_index]->residue(ires).is_lower_terminus() ) {
						pose.append_residue_by_bond(templates_[initial_template_index]->residue(ires));
					} else {
						pose.append_residue_by_jump(templates_[initial_template_index]->residue(ires), 1);
					}
					hetatms.push_back( std::make_pair( ires, pose.total_residue() ) );
				}
			}
		}
		else if (domain_assembly_) {
			DomainAssembly domain_assembly(templates_, domain_assembly_weights_);
			domain_assembly.run();
		} else if (realign_domains_ || realign_domains_stage2_) {
			// realign each template to the starting template by domain
			// does not to domain realignment if in domain assembly mode
			// domain parsing
			DDomainParse ddom(pcut_,hcut_,length_);
			utility::vector1< utility::vector1< loops::Loops > > domains_all_templ;
			domains_all_templ.resize( templates_.size() );
			for (Size i_template=1; i_template<=templates_.size(); ++i_template) {
				if (template_clusterID_[i_template] != template_clusterID_[initial_template_index]) continue;
				domains_all_templ[i_template] = ddom.split( *templates_[i_template], nres_protein_tgt );

				//protocols::loops::Loops my_chunks(template_chunks_[initial_template_index_]);
				// convert domain numbering to target pose numbering
				for (Size iloops=1; iloops<=domains_all_templ[i_template].size(); ++iloops) {
					for (Size iloop=1; iloop<=domains_all_templ[i_template][iloops].num_loop(); ++iloop) {
						Size seqpos_start_pose = templates_[i_template]->pdb_info()->number(domains_all_templ[i_template][iloops][iloop].start());
						domains_all_templ[i_template][iloops][iloop].set_start( seqpos_start_pose );

						Size seqpos_stop_pose = templates_[i_template]->pdb_info()->number(domains_all_templ[i_template][iloops][iloop].stop());
						domains_all_templ[i_template][iloops][iloop].set_stop( seqpos_stop_pose );
					}
				}

				TR << "Found " << domains_all_templ[i_template].size() << " domains for template " << template_fn_[i_template] << std::endl;
				for (Size i=1; i<=domains_all_templ[i_template].size(); ++i) {
					TR << "domain " << i << ": " << domains_all_templ[i_template][i] << std::endl;
				}
			}

			// combine domains that are not in the initial template
			domains_ = expand_domains_to_full_length(domains_all_templ, initial_template_index, nres_tgt);
			TR << "Final decision: " << domains_.size() << " domains" << std::endl;
			for (Size i=1; i<= domains_.size(); ++i) {
				TR << "domain " << i << ": " << domains_[i] << std::endl;
			}

			// local align
			align_by_domain(templates_, domains_, templates_[initial_template_index]);

			// update chunk, contig informations
			for (Size i_template=1; i_template<=templates_.size(); ++i_template) {
				// default minimum length is 3 and CA distance is 4
				template_contigs_[i_template] = protocols::loops::extract_continuous_chunks(*templates_[i_template]); // for chunk insertions
				template_chunks_[i_template] = protocols::loops::extract_secondary_structure_chunks(*templates_[i_template], "HE", 3, 6, 3, 4); // for fold tree setup
				if (template_chunks_[i_template].num_loop() == 0)
					template_chunks_[i_template] = template_contigs_[i_template];
			}
		}

		// symmetrize
		std::string symmdef_file = symmdef_files_[initial_template_index];
		if (!symmdef_file.empty() && symmdef_file != "NULL") {
			protocols::simple_moves::symmetry::SetupForSymmetryMover makeSymm( symmdef_file );
			makeSymm.apply(pose);

			//fpd   to get the right rotamer set we need to do this
			basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].value( "dummy" );

			// xyz copy hetatms -- this makes certain symmetries a little easier to setup
			if (add_hetatm_) {
				for ( Size ihet=1; ihet <= hetatms.size(); ++ihet ) {
					core::conformation::Residue const &res_in = templates_[initial_template_index]->residue(hetatms[ihet].first);

					for (Size iatm=1; iatm<=res_in.natoms(); ++iatm) {
						core::id::AtomID tgt(iatm,hetatms[ihet].second);
						pose.set_xyz( tgt, res_in.xyz( iatm ) );
					}
				}
			}
		}

		// set pose for density scoring if a map was input
		// >> keep this after symmetry
		if ( option[ OptionKeys::edensity::mapfile ].user() || user_csts_.size() > 0) {
			MoverOP dens( new protocols::electron_density::SetupForDensityScoringMover );
			dens->apply( pose );
		}

		// initialize template history
		// >> keep this after symmetry
		TemplateHistoryOP history = new TemplateHistory(pose);
		history->setall( initial_template_index_icluster );
		pose.data().set( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY, history );

    allowed_to_move_.clear();
    allowed_to_move_.resize(pose.total_residue(),true);

		// if a task factory is given, use it to set allowable residues
		if( task_factory_ ){
			task_ = task_factory_->create_task_and_apply_taskoperations( pose );
			for( core::Size resi = 1; resi <= get_num_residues_nonvirt(pose); ++resi )
				allowed_to_move_[resi] = ( task_->residue_task( resi ).being_designed() || task_->residue_task( resi ).being_packed()) ;
		}

		// STAGE 1
		//fpd constraints are handled a little bit weird
		//  * foldtree hybridize sets chainbreak variants then applies constraints (so c-beta csts are treated properly)
		//  * after chainbreak variants are removed, csts are removed and reapplied in this function
		//  * finally after switch to fullatom CSTs are reapplied again (optionally using a different CST file)
		// If "AUTO" is specified for cen_csts automatic centroid csts are generated
		// If "AUTO" is specified for fa_csts automatic fa csts are generated
		// If "NONE" is specified for fa_csts, cen csts are used (AUTO or otherwise)
		// An empty string is treated as equivalent to "NONE"

		// fold tree hybridize
		if (RG.uniform() < stage1_probability_) {
			for ( core::Size repeatstage1=0; repeatstage1 < jump_move_repeat_; ++repeatstage1 ) {
				std::string cst_fn = template_cst_fn_[initial_template_index];

				FoldTreeHybridizeOP ft_hybridize(
					new FoldTreeHybridize(
						initial_template_index_icluster, templates_icluster, weights_icluster,
						template_chunks_icluster, template_contigs_icluster, frags_small, frags_big) ) ;
				ft_hybridize->set_constraint_file( cst_fn );
				ft_hybridize->set_scorefunction( stage1_scorefxn_ );
				ft_hybridize->set_increase_cycles( stage1_increase_cycles_ );
				ft_hybridize->set_stage1_1_cycles( stage1_1_cycles_ );
				ft_hybridize->set_stage1_2_cycles( stage1_2_cycles_ );
				ft_hybridize->set_stage1_3_cycles( stage1_3_cycles_ );
				ft_hybridize->set_stage1_4_cycles( stage1_4_cycles_ );
				ft_hybridize->set_add_non_init_chunks( add_non_init_chunks_ );
				ft_hybridize->set_domain_assembly( domain_assembly_ );
				ft_hybridize->set_add_hetatm( add_hetatm_, hetatm_self_cst_weight_, hetatm_prot_cst_weight_ );
				ft_hybridize->set_frag_1mer_insertion_weight( frag_1mer_insertion_weight_ );
				ft_hybridize->set_small_frag_insertion_weight( small_frag_insertion_weight_ );
				ft_hybridize->set_big_frag_insertion_weight( big_frag_insertion_weight_ );
				ft_hybridize->set_chunk_insertion_weight( chunk_insertion_weight_ );
				ft_hybridize->set_frag_weight_aligned( frag_weight_aligned_ );
				ft_hybridize->set_auto_frag_insertion_weight( auto_frag_insertion_weight_ );
				ft_hybridize->set_max_registry_shift( max_registry_shift_ );

				// strand pairings
				ft_hybridize->set_pairings_file( pairings_file_ );
				ft_hybridize->set_sheets( sheets_ );
				ft_hybridize->set_random_sheets( random_sheets_ );
				ft_hybridize->set_filter_templates( filter_templates_ );

				// other cst stuff
				//ft_hybridize->set_movable_region( allowed_to_move_ );
				ft_hybridize->set_task_factory( task_factory_ );
				ft_hybridize->set_user_csts( user_csts_ );

				// finally run stage 1
				ft_hybridize->apply(pose);

				// get strand pairings if they exist for constraints in later stages
				strand_pairs_ = ft_hybridize->get_strand_pairs();

				//jump perturbation and minimization
				if ( jump_move_ ) {
					// call docking or symm docking
					if (core::pose::symmetry::is_symmetric(pose) ) {
						protocols::symmetric_docking::SymDockingLowResOP docking_lowres_mover = new protocols::symmetric_docking::SymDockingLowRes(stage1_scorefxn_);
						docking_lowres_mover->apply(pose);

						core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
						mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( true );
						core::pose::symmetry::make_symmetric_movemap( pose, *mm );

						protocols::simple_moves::symmetry::SymMinMoverOP min_mover =
							new protocols::simple_moves::symmetry::SymMinMover( mm, stage1_scorefxn_, "lbfgs_armijo_nonmonotone", 0.01, true );
						min_mover->apply(pose);

						core::pose::PoseOP stage1pose = new core::pose::Pose();
						core::pose::symmetry::extract_asymmetric_unit(pose, *stage1pose);
						align_by_domain(templates_, domains_, stage1pose);
					} else {
						core::Size const rb_move_jump = 1; // use the first jump as the one between partners <<<< fpd: MAKE THIS A PARSIBLE OPTION
						protocols::docking::DockingLowResOP docking_lowres_mover = new protocols::docking::DockingLowRes( stage1_scorefxn_, rb_move_jump );
						docking_lowres_mover->apply(pose);

						core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
						mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( rb_move_jump, true );
						protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm, stage1_scorefxn_, "lbfgs_armijo_nonmonotone", 0.01, true );
						min_mover->apply(pose);

						core::pose::PoseOP stage1pose = new core::pose::Pose( pose );
						align_by_domain(templates_, domains_, stage1pose);
					}
				}
			} //end of repeatstage1
		} else {
			// just do frag insertion in unaligned regions
			core::pose::PoseOP chosen_templ = templates_icluster[initial_template_index_icluster];
			protocols::loops::Loops chosen_contigs = template_contigs_icluster[initial_template_index_icluster];
			initialize_and_sample_loops(pose, chosen_templ, chosen_contigs, stage1_scorefxn_);
		}

		if (realign_domains_stage2_) {
			// realign domains to the output of stage 1
			TR << "Realigning template domains to stage1 pose." << std::endl;
			core::pose::PoseOP stage1pose = new core::pose::Pose( pose );
			align_by_domain(templates_, domains_, stage1pose);
		}


		//write gdtmm to output
		if (native_ && native_->total_residue()) {
			gdtmm = get_gdtmm(*native_, pose, native_aln);
			core::pose::setPoseExtraScores( pose, "GDTMM_after_stage1", gdtmm);
			TR << "GDTMM_after_stage1" << F(8,3,gdtmm) << std::endl;
		}

		// STAGE 2
		// apply constraints
		if ( stage2_scorefxn_->get_weight( core::scoring::atom_pair_constraint ) != 0 ) {
			std::string cst_fn = template_cst_fn_[initial_template_index];
			if (!keep_pose_constraint_ )
					setup_centroid_constraints( pose, templates_, template_weights_, cst_fn );
			if (add_hetatm_)
				add_non_protein_cst(pose, hetatm_self_cst_weight_, hetatm_prot_cst_weight_);
			if (strand_pairs_.size())
				add_strand_pairs_cst(pose, strand_pairs_);
			if ( task_factory_ )
				setup_interface_atompair_constraints(pose,allowed_to_move_);
		}

		if ( stage2_scorefxn_->get_weight( core::scoring::coordinate_constraint ) != 0) {
			if ( user_csts_.size() > 0 )
				setup_user_coordinate_constraints(pose,user_csts_);
			if ( task_factory_ )
				setup_interface_coordinate_constraints(pose,allowed_to_move_);
		}

		if (!option[cm::hybridize::skip_stage2]()) {
			core::scoring::ScoreFunctionOP stage2_scorefxn_clone = stage2_scorefxn_->clone();

			core::scoring::methods::EnergyMethodOptions lowres_options(stage2_scorefxn_clone->energy_method_options());
			lowres_options.set_cartesian_bonded_linear(true);
			stage2_scorefxn_clone->set_energy_method_options(lowres_options);

			CartesianHybridizeOP cart_hybridize (
				new CartesianHybridize(
					templates_icluster, weights_icluster,
					template_chunks_icluster,template_contigs_icluster, frags_big ) );
			cart_hybridize->set_scorefunction( stage2_scorefxn_);
			cart_hybridize->set_increase_cycles( stage2_increase_cycles_ );
			cart_hybridize->set_no_global_frame( no_global_frame_ );
			cart_hybridize->set_linmin_only( linmin_only_ );
			cart_hybridize->set_nofragbias( nofragbias_ );
			cart_hybridize->set_seqfrags_only( seqfrags_only_ );
			cart_hybridize->set_cartfrag_overlap( cartfrag_overlap_ );
			cart_hybridize->set_skip_long_min( skip_long_min_ );

			// finally run stage 2
			cart_hybridize->apply(pose);
		}

		//write gdtmm to output
		if (native_ && native_->total_residue()) {
			gdtmm = get_gdtmm(*native_, pose, native_aln);
			core::pose::setPoseExtraScores( pose, "GDTMM_after_stage2", gdtmm);
			TR << "GDTMM_after_stage2" << ObjexxFCL::format::F(8,3,gdtmm) << std::endl;
		}
		// get fragment history
		runtime_assert( pose.data().has( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY ) );
		history = *( static_cast< TemplateHistory* >( pose.data().get_ptr( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY )() ));

		TR << "History :";
		for (Size i=1; i<= history->size(); ++i ) { TR << I(4,i); }
		TR << std::endl;
		TR << "History :";
		for (Size i=1; i<= history->size(); ++i ) { TR << I(4, history->get(i)); }
		TR << std::endl;

		core::kinematics::MoveMapOP mm=new core::kinematics::MoveMap;
		mm->set_bb  ( true );
		mm->set_chi ( true );
		mm->set_jump( true );

		// stage "2.5" .. minimize with centroid energy + full-strength cart bonded
		if (!option[cm::hybridize::skip_stage2]()) {
			core::optimization::MinimizerOptions options_lbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
			core::optimization::CartesianMinimizer minimizer;
			core::kinematics::MoveMapOP mm=new core::kinematics::MoveMap;
			if (core::pose::symmetry::is_symmetric(pose) )
				core::pose::symmetry::make_symmetric_movemap( pose, *mm );
			Size n_min_cycles =(Size) (200.*stage25_increase_cycles_);
			options_lbfgs.max_iter(n_min_cycles);
			(*stage2_scorefxn_)(pose); minimizer.run( pose, *mm, *stage2_scorefxn_, options_lbfgs );
		}

		// STAGE 3: RELAX
		if (batch_relax_ > 0) {
			// set disulfides before going to FA
			if (disulf_file_.length() > 0) {
				// manual disulfide
				TR << " add disulfide: " << std::endl;
				basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].value(disulf_file_);
				core::pose::initialize_disulfide_bonds(pose);
				// must reset this since initialize_disulfide_bonds is used in pose io
				basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].deactivate();
				basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].to_default(); // reset to the default value
			} else {
				pose.conformation().detect_disulfides();
			}

			protocols::moves::MoverOP tofa = new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD );
			tofa->apply(pose);

			// apply fa constraints
			std::string cst_fn = template_cst_fn_[initial_template_index];
			if ( stage2_scorefxn_->get_weight( core::scoring::atom_pair_constraint ) != 0 ) {
					if (!keep_pose_constraint_ ) {
							setup_fullatom_constraints( pose, templates_, template_weights_, cst_fn, fa_cst_fn_ );
					} else {
							pose.constraint_set(save_pose_constraint_set);
					}
				if (add_hetatm_) {
					add_non_protein_cst(pose, hetatm_self_cst_weight_, hetatm_prot_cst_weight_);
				}
				if (strand_pairs_.size()) {
					add_strand_pairs_cst(pose, strand_pairs_);
				}
				if ( task_factory_ ) {
			  	setup_interface_atompair_constraints(pose,allowed_to_move_);
				}
			}

			if ( stage2_scorefxn_->get_weight( core::scoring::coordinate_constraint ) != 0) {
				// note that CSTs are updated based on stage 1 movement
				//   this may or may not be desirable
				if ( user_csts_.size() > 0 )
					setup_user_coordinate_constraints(pose,user_csts_);
				if ( task_factory_ )
					setup_interface_coordinate_constraints(pose,allowed_to_move_);
			}

			if (batch_relax_ == 1) {
				// standard relax
				// do relax _without_ ramping down coordinate constraints
				TR << " batch_relax 1 : " << std::endl;
				protocols::relax::FastRelax relax_prot( fa_scorefxn_, relax_repeats_ ,"NO CST RAMPING" );
				relax_prot.min_type("lbfgs_armijo_nonmonotone");
				relax_prot.apply(pose);
			} else {
				// batch relax
				// add stucture to queue
				SilentStructOP new_struct = SilentStructFactory::get_instance()->get_silent_struct("binary");
				new_struct->fill_struct( pose );
				new_struct->energies_from_pose( pose );
				post_centroid_structs.push_back( new_struct );

				if (post_centroid_structs.size() == batch_relax_) {
					protocols::relax::FastRelax relax_prot( fa_scorefxn_ );
					relax_prot.min_type("lbfgs_armijo_nonmonotone");
					relax_prot.set_force_nonideal(true);
					relax_prot.set_script_to_batchrelax_default( relax_repeats_ );

					// need to use a packer task factory to handle poses with different disulfide patterning
					core::pack::task::TaskFactoryOP tf = new core::pack::task::TaskFactory;
					tf->push_back( new core::pack::task::operation::InitializeFromCommandline );
					tf->push_back( new core::pack::task::operation::IncludeCurrent );
					tf->push_back( new core::pack::task::operation::RestrictToRepacking );
					relax_prot.set_task_factory( tf );

					// notice! this assumes all poses in a set have the same constraints!
					relax_prot.batch_apply(post_centroid_structs, pose.constraint_set()->clone());

					// reinflate pose
					post_centroid_structs[0]->fill_pose( pose );
				} else {
					protocols::moves::MoverOP tocen = new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );
					tocen->apply(pose);
					need_more_samples = true;
				}
			}
		} else {
			// no fullatom sampling
			(*stage2_scorefxn_)(pose);
		}
	}
	if (native_ && native_->total_residue()) {
		gdtmm = get_gdtmm(*native_, pose, native_aln);
		core::pose::setPoseExtraScores( pose, "GDTMM_final", gdtmm);
		TR << "GDTMM_final" << ObjexxFCL::format::F(8,3,gdtmm) << std::endl;
	}
}

utility::vector1 <Loops>
HybridizeProtocol::expand_domains_to_full_length(utility::vector1 < utility::vector1 < Loops > > all_domains, Size ref_domains_index, Size n_residues)
{
	utility::vector1 < Loops > domains = all_domains[ref_domains_index];
	if (all_domains[ref_domains_index].size() == 0) return domains;

	utility::vector1 < bool > residue_mask(n_residues, false);
	residue_mask.resize(n_residues);

	//mask residues from all domains not to be cut
	for (Size i_pose=1; i_pose <= all_domains.size(); ++i_pose) {
		for (Size idomain=1; idomain <= all_domains[i_pose].size(); ++idomain) {
			for (core::Size iloop=1; iloop<=all_domains[i_pose][idomain].num_loop(); ++iloop) {
				for (core::Size ires=all_domains[i_pose][idomain][iloop].start()+1; ires<=all_domains[i_pose][idomain][iloop].stop(); ++ires) {
					residue_mask[ires] = true;
				}
			}
		}
	}

	// assuming loops in domains are sequential
	for (Size idomain=1; idomain <= all_domains[ref_domains_index].size(); ++idomain) {
		for (core::Size iloop=1; iloop<=all_domains[ref_domains_index][idomain].num_loop(); ++iloop) {
			if (idomain == 1 && iloop == 1) {
				domains[idomain][iloop].set_start(1);
			}
			if (idomain == all_domains[ref_domains_index].size() && iloop == all_domains[ref_domains_index][idomain].num_loop()) {
				domains[idomain][iloop].set_stop(n_residues);
			}

			// the next loop
			Size jdomain = idomain;
			Size jloop = iloop+1;
			if (jloop > all_domains[ref_domains_index][idomain].num_loop()) {
				++jdomain;
				jloop = 1;
				if (jdomain > all_domains[ref_domains_index].size()) continue;

				Size gap_start = all_domains[ref_domains_index][idomain][iloop].stop()+1;
				Size gap_stop  = all_domains[ref_domains_index][jdomain][jloop].start();
				utility::vector1<Size> cut_options;
				for (Size ires=gap_start; ires<=gap_stop; ++ires) {
					if (residue_mask[ires]) continue;
					cut_options.push_back(ires);
				}
				if (cut_options.size() == 0) {
					for (Size ires=gap_start; ires<=gap_stop; ++ires) {
						cut_options.push_back(ires);
					}
				}
				Size cut = cut_options[RG.random_range(1,cut_options.size())];

				domains[idomain][iloop].set_stop(cut-1);
				domains[jdomain][jloop].set_start(cut);
			}
		}
	}
	return domains;
}

void
HybridizeProtocol::align_by_domain(utility::vector1<core::pose::PoseOP> & poses, utility::vector1 < Loops > domains, core::pose::PoseOP & ref_pose)
{
	for (Size i_pose=1; i_pose <= poses.size(); ++i_pose) {
		if (poses[i_pose] == ref_pose) continue;
		align_by_domain(*poses[i_pose], *ref_pose, domains);

		//std::string out_fn = template_fn_[i_pose] + "_realigned.pdb";
		//poses[i_pose]->dump_pdb(out_fn);
	}
}


void
HybridizeProtocol::align_by_domain(core::pose::Pose & pose, core::pose::Pose const & ref_pose, utility::vector1 <Loops> domains)
{
	for (Size i_domain = 1; i_domain <= domains.size() ; ++i_domain) {
		core::id::AtomID_Map< core::id::AtomID > atom_map;
		core::pose::initialize_atomid_map( atom_map, pose, core::id::BOGUS_ATOM_ID );
		std::list <Size> residue_list;
		core::Size n_mapped_residues=0;
		for (core::Size ires=1; ires<=pose.total_residue(); ++ires) {
			if (!pose.residue_type(ires).is_protein()) continue;
			int pose_res = (pose.pdb_info()) ? pose.pdb_info()->number(ires) : ires;
			for (core::Size iloop=1; iloop<=domains[i_domain].num_loop(); ++iloop) {
				if ( pose_res < (int)domains[i_domain][iloop].start() || pose_res > (int)domains[i_domain][iloop].stop() ) continue;
				residue_list.push_back(ires);
			}
		}
		std::list <Size> ref_residue_list;
		for (core::Size jres=1; jres<=ref_pose.total_residue(); ++jres) {
			if (!ref_pose.residue_type(jres).is_protein()) continue;
			int ref_pose_res = (ref_pose.pdb_info()) ? ref_pose.pdb_info()->number(jres) : jres;
			for (core::Size iloop=1; iloop<=domains[i_domain].num_loop(); ++iloop) {
				if ( ref_pose_res < (int)domains[i_domain][iloop].start() || ref_pose_res > (int)domains[i_domain][iloop].stop() ) continue;
				ref_residue_list.push_back(jres);
			}
		}

		TMalign tm_align;
		std::string seq_pose, seq_ref, aligned;
		int reval = tm_align.apply(pose, ref_pose, residue_list, ref_residue_list);
		if (reval == 0) {
			tm_align.alignment2AtomMap(pose, ref_pose, residue_list, ref_residue_list, n_mapped_residues, atom_map);
			tm_align.alignment2strings(seq_pose, seq_ref, aligned);

			using namespace ObjexxFCL::format;
			Size norm_length = residue_list.size() < ref_residue_list.size() ? residue_list.size():ref_residue_list.size();
			TR << "Align domain with TMscore of " << F(8,3,tm_align.TMscore(norm_length)) << std::endl;
			TR << seq_pose << std::endl;
			TR << aligned << std::endl;
			TR << seq_ref << std::endl;

			if (n_mapped_residues >= 6) {
				utility::vector1< core::Real > aln_cutoffs;
				aln_cutoffs.push_back(6);
				aln_cutoffs.push_back(4);
				aln_cutoffs.push_back(3);
				aln_cutoffs.push_back(2);
				aln_cutoffs.push_back(1.5);
				aln_cutoffs.push_back(1);
				core::Real min_coverage = 0.2;
				partial_align(pose, ref_pose, atom_map, residue_list, true, aln_cutoffs, min_coverage); // iterate_convergence = true
			}
		}
		//else {
		//	TR << "This domain cannot be aligned: " << n_mapped_residues<< std::endl;
		//	TR << domains[i_domain];
		//}
	}
}

protocols::moves::MoverOP HybridizeProtocol::clone() const { return new HybridizeProtocol( *this ); }
protocols::moves::MoverOP HybridizeProtocol::fresh_instance() const { return new HybridizeProtocol; }

std::string
HybridizeProtocol::get_name() const {
	return "HybridizeProtocol";
}

void
HybridizeProtocol::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const & pose )
{
	// config file
	if( tag->hasOption( "config_file" ) )
		read_template_structures( tag->getOption< std::string >( "config_file" ) );

	// basic options
	stage1_increase_cycles_ = tag->getOption< core::Real >( "stage1_increase_cycles", 1. );
	stage2_increase_cycles_ = tag->getOption< core::Real >( "stage2_increase_cycles", 1. );
	stage25_increase_cycles_ = tag->getOption< core::Real >( "stage2.5_increase_cycles", 1. );
	fa_cst_fn_ = tag->getOption< std::string >( "fa_cst_file", "" );
	batch_relax_ = tag->getOption< core::Size >( "batch" , 1 );
	jump_move_= tag->getOption< bool >( "jump_move" , false );
	jump_move_repeat_= tag->getOption< core::Size >( "jump_move_repeat" , 1 );
	keep_pose_constraint_= tag->getOption< bool >( "keep_pose_constraint" , false );

	if( tag->hasOption( "task_operations" ) ){
		task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data );
		task_ = task_factory_->create_task_and_apply_taskoperations( pose );
	}	else {
		task_factory_ = NULL;
	}

	// force starting template
	if( tag->hasOption( "starting_template" ) ) {
		utility::vector1<std::string> buff = utility::string_split( tag->getOption<std::string>( "starting_template" ), ',' );
		BOOST_FOREACH(std::string field, buff){
			Size const value = std::atoi( field.c_str() );
			starting_templates_.push_back(value);
		}
	}

	// tons of ab initio options
	if( tag->hasOption( "stage1_1_cycles" ) )
		stage1_1_cycles_ = tag->getOption< core::Size >( "stage1_1_cycles" );
	if( tag->hasOption( "stage1_2_cycles" ) )
		stage1_2_cycles_ = tag->getOption< core::Size >( "stage1_2_cycles" );
	if( tag->hasOption( "stage1_3_cycles" ) )
		stage1_3_cycles_ = tag->getOption< core::Size >( "stage1_3_cycles" );
	if( tag->hasOption( "stage1_4_cycles" ) )
		stage1_4_cycles_ = tag->getOption< core::Size >( "stage1_4_cycles" );
	if( tag->hasOption( "stage1_probability" ) )
		stage1_probability_ = tag->getOption< core::Real >( "stage1_probability" );
	if( tag->hasOption( "add_hetatm" ) )
		add_hetatm_ = tag->getOption< bool >( "add_hetatm" );
	if( tag->hasOption( "hetatm_cst_weight" ) )
		hetatm_self_cst_weight_ = tag->getOption< core::Real >( "hetatm_cst_weight" );
	if( tag->hasOption( "hetatm_to_protein_cst_weight" ) )
		hetatm_prot_cst_weight_ = tag->getOption< core::Real >( "hetatm_to_protein_cst_weight" );
	if( tag->hasOption( "domain_assembly" ) )
		domain_assembly_ = tag->getOption< bool >( "domain_assembly" );
	if( tag->hasOption( "realign_domains" ) )
		realign_domains_ = tag->getOption< bool >( "realign_domains" );
	if( tag->hasOption( "realign_domains_stage2" ) )
		realign_domains_stage2_ = tag->getOption< bool >( "realign_domains_stage2" );
	if( tag->hasOption( "add_non_init_chunks" ) )
		add_non_init_chunks_ = tag->getOption< bool >( "add_non_init_chunks" );
	if( tag->hasOption( "frag_1mer_insertion_weight" ) )
		frag_1mer_insertion_weight_ = tag->getOption< core::Real >( "frag_1mer_insertion_weight" );
	if( tag->hasOption( "small_frag_insertion_weight" ) )
		small_frag_insertion_weight_ = tag->getOption< core::Real >( "small_frag_insertion_weight" );
	if( tag->hasOption( "big_frag_insertion_weight" ) )
		big_frag_insertion_weight_ = tag->getOption< core::Real >( "big_frag_insertion_weight" );
	if( tag->hasOption( "chunk_insertion_weight" ) )
		chunk_insertion_weight_ = tag->getOption< core::Real >( "chunk_insertion_weight" );
	if( tag->hasOption( "frag_weight_aligned" ) )
		frag_weight_aligned_ = tag->getOption< core::Real >( "frag_weight_aligned" );
	if( tag->hasOption( "auto_frag_insertion_weight" ) )
		auto_frag_insertion_weight_ = tag->getOption< bool >( "auto_frag_insertion_weight" );
	if( tag->hasOption( "max_registry_shift" ) )
		max_registry_shift_ = tag->getOption< core::Size >( "max_registry_shift" );
	if( tag->hasOption( "repeats" ) )
		relax_repeats_ = tag->getOption< core::Size >( "repeats" );
	if( tag->hasOption( "disulf_file" ) )
		disulf_file_ = tag->getOption< std::string >( "disulf_file" );


	// stage 2-specific options
	if( tag->hasOption( "no_global_frame" ) )
		no_global_frame_ = tag->getOption< bool >( "no_global_frame" );
	if( tag->hasOption( "linmin_only" ) )
		linmin_only_ = tag->getOption< bool >( "linmin_only" );
	if( tag->hasOption( "cartfrag_overlap" ) )
		cartfrag_overlap_ = tag->getOption< core::Size >( "cartfrag_overlap" );
	if( tag->hasOption( "seqfrags_only" ) )
		seqfrags_only_ = tag->getOption< core::Size >( "seqfrags_only" );
	if( tag->hasOption( "nofragbias" ) )
		nofragbias_ = tag->getOption< core::Size >( "nofragbias" );
	if( tag->hasOption( "skip_long_min" ) )
		skip_long_min_ = tag->getOption< core::Size >( "skip_long_min" );


	// scorefxns
	if( tag->hasOption( "stage1_scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "stage1_scorefxn" ) );
		stage1_scorefxn_ = (data.get< ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	}
	if( tag->hasOption( "stage2_scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "stage2_scorefxn" ) );
		stage2_scorefxn_ = (data.get< ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	}
	if( tag->hasOption( "fa_scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "fa_scorefxn" ) );
		fa_scorefxn_ = (data.get< ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	}

	// ddomain options
	hcut_ = tag->getOption< core::Real >( "domain_hcut" , 0.18);
	pcut_ = tag->getOption< core::Real >( "domain_pcut" , 0.81);
	length_ = tag->getOption< core::Size >( "domain_length" , 38);

	// user constraints
	if( tag->hasOption( "coord_cst_res" ) ) {
		user_csts_ = core::pose::get_resnum_list_ordered( tag->getOption<std::string>( "coord_cst_res" ), pose );
	}

	// if user constraints are defined, make sure coord_csts are defined in at least one stage
	if (user_csts_.size() > 0) {
		runtime_assert(
			stage1_scorefxn_->get_weight( core::scoring::coordinate_constraint ) > 0 ||
			stage2_scorefxn_->get_weight( core::scoring::coordinate_constraint ) > 0 ||
			fa_scorefxn_->get_weight( core::scoring::coordinate_constraint ) > 0 );
	}

	// fragments
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for (tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it) {
		if ( (*tag_it)->getName() == "Fragments" ) {
			using namespace core::fragment;
			if ( (*tag_it)->hasOption( "3mers" ) ) {
				core::fragment::FragSetOP frags = core::fragment::FragmentIO().read_data( (*tag_it)->getOption<std::string>( "3mers" )  );
				fragments_small_.push_back(frags);
			} else if ( (*tag_it)->hasOption( "small" ) ) {
				utility::vector1<std::string> frag_files = utility::string_split((*tag_it)->getOption<std::string>( "small" ), ',');
				for (core::Size i=1; i<= frag_files.size(); ++i )
					fragments_small_.push_back(core::fragment::FragmentIO().read_data( frag_files[i] ));
			}
			if ( (*tag_it)->hasOption( "9mers" ) ) {
				core::fragment::FragSetOP frags = core::fragment::FragmentIO().read_data( (*tag_it)->getOption<std::string>( "9mers" ) );
				fragments_big_.push_back(frags);
			} else if ( (*tag_it)->hasOption( "big" ) ) {
				utility::vector1<std::string> frag_files = utility::string_split((*tag_it)->getOption<std::string>( "big" ), ',');
				for (core::Size i=1; i<= frag_files.size(); ++i )
					fragments_big_.push_back(core::fragment::FragmentIO().read_data( frag_files[i] ));
			}
		}

		if ( (*tag_it)->getName() == "Template" ) {
			std::string template_fn = (*tag_it)->getOption<std::string>( "pdb" );
			std::string cst_fn = (*tag_it)->getOption<std::string>( "cst_file", "AUTO" );
			core::Real weight = (*tag_it)->getOption<core::Real>( "weight", 1 );
			core::Real domain_assembly_weight = (*tag_it)->getOption<core::Real>( "domain_assembly_weight", 0. );
			core::Size cluster_id = (*tag_it)->getOption<core::Size>( "cluster_id", 1 );
			std::string symm_file = (*tag_it)->getOption<std::string>( "symmdef", "" );
			utility::vector1<core::Size> cst_reses;
			add_template(template_fn, cst_fn, symm_file, weight, domain_assembly_weight, cluster_id, cst_reses);
		}

		// strand pairings
		if ( (*tag_it)->getName() == "Pairings" ) {
			pairings_file_ = (*tag_it)->getOption< std::string >( "file", "" );
			if ( (*tag_it)->hasOption("sheets") ) {
				core::Size sheets = (*tag_it)->getOption< core::Size >( "sheets" );
				sheets_.clear();
				random_sheets_.clear();
				sheets_.push_back(sheets);
			} else if ( (*tag_it)->hasOption("random_sheets") ) {
				core::Size random_sheets = (*tag_it)->getOption< core::Size >( "random_sheets" );
				sheets_.clear();
				random_sheets_.clear();
				random_sheets_.push_back(random_sheets);
			}
			filter_templates_ = (*tag_it)->getOption<bool>( "filter_templates" , 0 );
		}

	}
}

} // hybridization
} // protocols
