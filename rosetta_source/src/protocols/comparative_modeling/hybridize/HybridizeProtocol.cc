// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#include <protocols/comparative_modeling/hybridize/HybridizeProtocol.hh>
#include <protocols/comparative_modeling/hybridize/FoldTreeHybridize.hh>
#include <protocols/comparative_modeling/hybridize/CartesianHybridize.hh>
#include <protocols/comparative_modeling/hybridize/TemplateHistory.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/FragmentMover.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/util.hh>

#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.hh>

#include <core/sequence/util.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>

#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// task operation
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/rosetta_scripts/util.hh>

// utility
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <numeric/random/WeightedSampler.hh>
#include <ObjexxFCL/format.hh>

// evaluation
#include <core/scoring/rms_util.hh>
#include <protocols/comparative_modeling/coord_util.hh>

// option
#include <basic/options/option.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <string>

static basic::Tracer TR( "protocols.comparative_modeling.hybridize.HybridizeProtocol" );
static numeric::random::RandomGenerator RG(541938);

namespace protocols {
namespace comparative_modeling {
namespace hybridize {

using namespace core;
using namespace kinematics;
using namespace sequence;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;
using namespace constraints;

HybridizeProtocol::HybridizeProtocol() :
	template_weights_sum_(0)
{
	check_options();
	
	//read templates
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
	read_template_structures( option[cm::hybridize::template_list]() );
	
	//read fragments
	if ( option[ in::file::frag9 ].user() ) {
		using namespace core::fragment;
		fragments9_ = new ConstantLengthFragSet( 9 );
		fragments9_ = FragmentIO().read_data( option[ in::file::frag9 ]() );
		for (core::fragment::FrameIterator i = fragments9_->begin(); i != fragments9_->end(); ++i) {
			core::Size position = (*i)->start();
			//library_[position] = **i;
		}
	}
	if ( option[ in::file::frag3 ].user() ) {
		using namespace core::fragment;
		fragments3_ = new ConstantLengthFragSet( 3 );
		fragments3_ = FragmentIO().read_data( option[ in::file::frag3 ]() );
		for (core::fragment::FrameIterator i = fragments9_->begin(); i != fragments9_->end(); ++i) {
			core::Size position = (*i)->start();
			//library_[position] = **i;
		}
	}
	
	// native
	if ( option[ in::file::native ].user() ) {
		native_ = new core::pose::Pose;
		core::import_pose::pose_from_pdb( *native_, option[ in::file::native ]() );
	}

}
	
core::Real HybridizeProtocol::get_gdtmm( core::pose::Pose &pose ) {
	core::Real gdtmm = 0;
	if ( native_ && native_->total_residue() > 0) {
		if ( !aln_ ) {
			core::sequence::SequenceOP model_seq ( new core::sequence::Sequence( pose.sequence(),  "model",  1 ) );
			core::sequence::SequenceOP native_seq( new core::sequence::Sequence( native_->sequence(), "native", 1 ) );
			aln_ = new core::sequence::SequenceAlignment;
			*aln_ = align_naive(model_seq,native_seq);
		}
		
		int n_atoms;
		ObjexxFCL::FArray2D< core::Real > p1a, p2a;
		protocols::comparative_modeling::gather_coords( pose, *native_, *aln_, n_atoms, p1a, p2a );
		
		core::Real m_1_1, m_2_2, m_3_3, m_4_3, m_7_4;
		gdtmm = core::scoring::xyz_gdtmm( p1a, p2a, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
	}
	return gdtmm;
}


void
HybridizeProtocol::initialize_and_sample_loops(
		core::pose::Pose &pose,
		core::pose::PoseOP chosen_templ,
		protocols::loops::Loops template_contigs_icluster,
		core::scoring::ScoreFunctionOP scorefxn)
{
	// xyz copy starting model
	for (int i=1; i<=chosen_templ->total_residue(); ++i)
		for (int j=1; j<=chosen_templ->residue(i).natoms(); ++j) {
			core::id::AtomID src(j,i), tgt(j, chosen_templ->pdb_info()->number(i));
			pose.set_xyz( tgt, chosen_templ->xyz( src ) );
		}

	// make loops as inverse of template_contigs_icluster
	TR << "CONTIGS" << std::endl << template_contigs_icluster << std::endl;
	core::Size ncontigs = template_contigs_icluster.size();
	core::Size nres_tgt = pose.total_residue();
	protocols::loops::LoopsOP loops = new protocols::loops::Loops;
	utility::vector1< bool > templ_coverage(nres_tgt, false);
//	for (int i=1; i<=ncontigs; ++i) {
//		core::Size cstart = chosen_templ->pdb_info()->number(template_contigs_icluster[i].start());
//		core::Size cstop = chosen_templ->pdb_info()->number(template_contigs_icluster[i].stop());
//		for (int j=cstart; j<=cstop; ++j) templ_coverage[j] = true;
//	}
	for (int i=1; i<=chosen_templ->total_residue(); ++i) {
		core::Size cres = chosen_templ->pdb_info()->number(i);
		templ_coverage[cres] = true;
	}

	// remove 1-3 residue "segments"
	for (int i=1; i<=nres_tgt-2; ++i) {
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
	for (int i=2; i<=nres_tgt; ++i) {
		if (templ_coverage[i] && inloop) {
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
		frag3mover = new protocols::simple_moves::ClassicFragmentMover( fragments3_, mm_loop );
		frag3mover->set_check_ss( false ); frag3mover->enable_end_bias_check( false );
		frag9mover = new protocols::simple_moves::ClassicFragmentMover( fragments9_, mm_loop );
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

		core::Size neffcycles = 1000;  //fpd  to do: base on # of loop residues
		for (int n=1; n<=neffcycles; ++n) {
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
		for (int n=1; n<=neffcycles; ++n) {
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



HybridizeProtocol::~HybridizeProtocol(){}

void HybridizeProtocol::add_template(
	std::string template_fn,
	std::string cst_fn,
	core::Real weight,
	core::Size cluster_id,
	utility::vector1<core::Size> cst_reses)
{
	core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	core::pose::PoseOP template_pose = new core::pose::Pose();
	core::import_pose::pose_from_pdb( *template_pose, *residue_set, template_fn );

	// add secondary structure information to the template pose
	core::scoring::dssp::Dssp dssp_obj( *template_pose );
	dssp_obj.insert_ss_into_pose( *template_pose );
	
	// find ss chunks in template
	//protocols::loops::Loops chunks = protocols::loops::extract_secondary_structure_chunks(*template_pose); 
	protocols::loops::Loops chunks = protocols::loops::extract_secondary_structure_chunks(*template_pose, "HE", 3, 6, 3, 4);
	TR.Debug << "Chunks from template\n" << chunks << std::endl;

	// break templates into contigs
	protocols::loops::Loops contigs = protocols::loops::extract_continuous_chunks(*template_pose); 
	TR.Debug << "Contigs from template\n" << contigs << std::endl;

	template_fn_.push_back(template_fn);
	templates_.push_back(template_pose);
	template_cst_fn_.push_back(cst_fn);
	template_weights_.push_back(weight);
	template_clusterID_.push_back(cluster_id);
	template_chunks_.push_back(chunks);
	template_contigs_.push_back(contigs);
	template_cst_reses_.push_back(cst_reses);
}

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
				std::vector<std::string> cst_reses_parsed = utility::string_split( cst_reses_str , ',' ) ;
				for (int i=0; i< cst_reses_parsed.size(); ++i ) {
					cst_reses.push_back( (core::Size) std::atoi( cst_reses_parsed[i].c_str() ) );
				}
			}
			add_template(template_fn, cst_fn, weight, cluster_id, cst_reses);
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

void HybridizeProtocol::check_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !basic::options::option[ cm::hybridize::template_list ].user() ) {
		utility_exit_with_message("Error! Need the -cm::hybridize::template_list for input template structures.");
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
	template_index_icluster.clear();
	templates_icluster.clear();
	weights_icluster.clear();
	template_chunks_icluster.clear();
	template_contigs_icluster.clear();
	
	numeric::random::WeightedSampler weighted_sampler;
	weighted_sampler.weights(template_weights_);

	initial_template_index = weighted_sampler.random_sample(RG);
	Size cluster_id = template_clusterID_[initial_template_index];
	
	for (Size i_template = 1; i_template <= template_clusterID_.size(); ++i_template) {
		if ( cluster_id == template_clusterID_[i_template] ) {
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
	
	// pick starting template
	core::Size initial_template_index;
	core::Size initial_template_index_icluster; // index in the cluster
	utility::vector1 < core::Size > template_index_icluster; // index look back
	utility::vector1 < core::pose::PoseOP > templates_icluster;
	utility::vector1 < core::Real > weights_icluster;
	utility::vector1 < protocols::loops::Loops > template_chunks_icluster;
	utility::vector1 < protocols::loops::Loops > template_contigs_icluster;
	pick_starting_template(initial_template_index, initial_template_index_icluster,
	                      template_index_icluster, templates_icluster, weights_icluster, template_chunks_icluster, template_contigs_icluster);

	using namespace ObjexxFCL::fmt;
	TR << "Using initial template: " << I(4,initial_template_index) << " " << template_fn_[initial_template_index] << std::endl;

	// initialize template history
	TemplateHistoryOP history = new TemplateHistory(pose);
	history->setall( initial_template_index_icluster );
	pose.data().set( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY, history );
	
	// apply constraints
	core::scoring::constraints::ConstraintSetOP constraint_set;
	std::string cst_fn = template_cst_fn_[initial_template_index];
	if (!cst_fn.empty() && cst_fn != "NONE") {
		constraint_set = ConstraintIO::get_instance()->read_constraints_new( cst_fn, new ConstraintSet, pose );
	}
	pose.constraint_set( constraint_set );
	
	// fold tree hybridize
	core::scoring::ScoreFunctionOP scorefxn_stage1 = 
		core::scoring::ScoreFunctionFactory::create_score_function(option[cm::hybridize::stage1_weights](), option[cm::hybridize::stage1_patch]());

	if (RG.uniform() < option[cm::hybridize::stage1_probability]()) {
		core::pose::PoseOP chosen_templ = templates_icluster[initial_template_index_icluster];
		protocols::loops::Loops chosen_contigs = template_contigs_icluster[initial_template_index_icluster];
		initialize_and_sample_loops(pose, chosen_templ, chosen_contigs, scorefxn_stage1);
	} else {
		// package up fragments
		utility::vector1 < core::fragment::FragSetOP > frag_libs;
		frag_libs.push_back(fragments3_);
		frag_libs.push_back(fragments9_);

		FoldTreeHybridizeOP ft_hybridize(
			new FoldTreeHybridize(initial_template_index_icluster, templates_icluster, weights_icluster, template_chunks_icluster, template_contigs_icluster,  frag_libs) ) ;
		ft_hybridize->set_scorefunction(scorefxn_stage1);
		ft_hybridize->apply(pose);
	}

	//write gdtmm to output
	using namespace ObjexxFCL::fmt;
	core::Real gdtmm = get_gdtmm(pose);
	core::pose::setPoseExtraScores( pose, "GDTMM_after_stage1", gdtmm);
	TR << "GDTMM_after_stage1" << F(8,3,gdtmm) << std::endl;

	//pose.dump_pdb( "after_stage1.pdb" );
	
	// cartesian fragment hybridize
	pose.constraint_set( constraint_set );  //reset constraints
	CartesianHybridizeOP cart_hybridize ( new CartesianHybridize( templates_icluster, weights_icluster,template_chunks_icluster,template_contigs_icluster, fragments9_ ) );
	core::scoring::ScoreFunctionOP scorefxn_stage2 =
	core::scoring::ScoreFunctionFactory::create_score_function(option[cm::hybridize::stage2_weights](), option[cm::hybridize::stage2_patch]());
	cart_hybridize->set_scorefunction(scorefxn_stage2);

	if (!option[cm::hybridize::skip_stage2]()) {
		cart_hybridize->apply(pose);
	}

	// get fragment history
	runtime_assert( pose.data().has( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY ) );
	history = *( static_cast< TemplateHistory* >( pose.data().get_ptr( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY )() ));

	//write gdtmm to output
	gdtmm = get_gdtmm(pose);
	core::pose::setPoseExtraScores( pose, "GDTMM_after_stage2", gdtmm);
	TR << "GDTMM_after_stage2" << F(8,3,gdtmm) << std::endl;
	
	//look back and apply constraints
	TR << "History :";
	for (int i=1; i<= history->size(); ++i ) {
		TR << I(4,i);
	}
	TR << std::endl;
	TR << "History :";
	for (int i=1; i<= history->size(); ++i ) {
		TR << I(4, history->get(i));
	}
	TR << std::endl;

	// optional relax
	if (option[ cm::hybridize::relax ]()) {
		protocols::moves::MoverOP tofa = new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD );
		tofa->apply(pose);

		// add constraints
		core::pose::addVirtualResAsRoot(pose);
		core::Size rootres = pose.fold_tree().root();
		Real const coord_sdev( 1 ); // this should be an option perhaps
		for ( core::Size i=1; i<pose.total_residue(); ++i ) {
			if ( !pose.residue(i).is_polymer() ) continue;
			utility::vector1<core::Size> const &source_list = template_cst_reses_[ template_index_icluster[ history->get(i) ] ];
			if ( std::find( source_list.begin(), source_list.end(), i ) != source_list.end() ) {
				TR << "Constrain residue " << i << std::endl;
				core::Size ii = pose.residue(i).atom_index(" CA ");
				pose.add_constraint(
					new CoordinateConstraint(
						core::id::AtomID(ii,i), core::id::AtomID(1,rootres), pose.xyz( core::id::AtomID(ii,i) ), new HarmonicFunc( 0.0, coord_sdev ) ) );
			}
		}

		// get scorefunction, add cst weight from commandline
		core::scoring::ScoreFunctionOP fa_scorefxn = core::scoring::getScoreFunction();
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *fa_scorefxn );

		// do relax _without_ ramping down coordinate constraints
		protocols::relax::FastRelax relax_prot( fa_scorefxn, option[ basic::options::OptionKeys::relax::default_repeats ]() ,"NO CST RAMPING" );
		relax_prot.set_min_type("lbfgs_armijo_nonmonotone");

		core::pose::Pose pose_pre_relax = pose;

relaxing:
		// fpd more nan problems with torsion derivs
		try {
			relax_prot.apply(pose);
		} catch( utility::excn::EXCN_Base& excn ) {
			//fpd hbond fail? start over
			pose = pose_pre_relax;
			goto relaxing;
		}

		gdtmm = get_gdtmm(pose);
		core::pose::setPoseExtraScores( pose, "GDTMM_final", gdtmm);
		TR << "GDTMM_final" << F(8,3,gdtmm) << std::endl;
	} else {
		(*scorefxn_stage2)(pose);
	}
}


protocols::moves::MoverOP HybridizeProtocol::clone() const { return new HybridizeProtocol( *this ); }
protocols::moves::MoverOP HybridizeProtocol::fresh_instance() const { return new HybridizeProtocol; }

std::string
HybridizeProtocol::get_name() const {
	return "HybridizeProtocol";
}

} // hybridize 
} // comparative_modeling 
} // protocols

