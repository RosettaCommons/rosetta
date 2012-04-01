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

#include <protocols/comparative_modeling/hybridize/HybridizeProtocolCreator.hh>
#include <protocols/comparative_modeling/hybridize/HybridizeProtocol.hh>
#include <protocols/comparative_modeling/hybridize/FoldTreeHybridize.hh>
#include <protocols/comparative_modeling/hybridize/CartesianHybridize.hh>
#include <protocols/comparative_modeling/hybridize/TemplateHistory.hh>
#include <protocols/comparative_modeling/hybridize/util.hh>
#include <protocols/comparative_modeling/hybridize/DDomainParse.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/util.hh>

#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>

#include <protocols/electron_density/SetupForDensityScoringMover.hh>

// dynamic fragpick
#include <protocols/moves/DsspMover.hh>
#include <core/fragment/picking_old/vall/util.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>

#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>

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
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
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

#include <protocols/moves/DataMap.hh>

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

#include <string>

#define foreach BOOST_FOREACH

static basic::Tracer TR( "protocols.comparative_modeling.hybridize.HybridizeProtocol" );
static numeric::random::RandomGenerator RG(541938);

namespace protocols {
namespace comparative_modeling {
namespace hybridize {

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
	template_weights_sum_(0),
	fragments3_(NULL),
	fragments9_(NULL)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	init();

	// initialization may come from command line or from RS
	if (option[cm::hybridize::template_list].user()) {
		read_template_structures( option[cm::hybridize::template_list]() );
	}
}

HybridizeProtocol::HybridizeProtocol(std::string template_list_file) :
	template_weights_sum_(0),
	fragments3_(NULL),
	fragments9_(NULL)
{
	init();
	read_template_structures(template_list_file);
}

// sets default options
void
HybridizeProtocol::init() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	stage1_probability_ = option[cm::hybridize::stage1_probability]();
	stage1_increase_cycles_ = option[cm::hybridize::stage1_increase_cycles]();
	realign_domains_ = option[cm::hybridize::realign_domains]();
	add_non_init_chunks_ = option[cm::hybridize::add_non_init_chunks]();
	frag_weight_aligned_ = option[cm::hybridize::frag_weight_aligned]();
	max_registry_shift_ = option[cm::hybridize::max_registry_shift]();
	cartfrag_overlap_ = 1;

	if (option[cm::hybridize::starting_template].user()) {
		starting_templates_ = option[cm::hybridize::starting_template]();
	}
	
	stage2_increase_cycles_ = option[cm::hybridize::stage2_increase_cycles]();
	no_global_frame_ = option[cm::hybridize::no_global_frame]();
	linmin_only_ = option[cm::hybridize::linmin_only]();

	// default scorefunction
	stage1_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( 
		option[cm::hybridize::stage1_weights](), option[cm::hybridize::stage1_patch]() );
	stage2_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( 
		option[cm::hybridize::stage2_weights](), option[cm::hybridize::stage2_patch]() );

	fa_scorefxn_ = core::scoring::getScoreFunction();
	core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *fa_scorefxn_ );

	if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
		utility::vector1< std::string > cst_files = option[ OptionKeys::constraints::cst_fa_file ]();
		core::Size choice = core::Size( RG.random_range( 1, cst_files.size() ) );
		fa_cst_fn_ = cst_files[choice];
		TR.Info << "Fullatom Constraint choice: " << fa_cst_fn_ << std::endl;
	}
	
	batch_relax_ = option[ cm::hybridize::relax ]();
	relax_repeats_ = option[ basic::options::OptionKeys::relax::default_repeats ]();

	// read fragments
	if ( option[ in::file::frag9 ].user() ) {
		using namespace core::fragment;
		fragments9_ = new ConstantLengthFragSet( 9 );
		fragments9_ = FragmentIO().read_data( option[ in::file::frag9 ]() );
	}
	if ( option[ in::file::frag3 ].user() ) {
		using namespace core::fragment;
		fragments3_ = new ConstantLengthFragSet( 3 );
		fragments3_ = FragmentIO().read_data( option[ in::file::frag3 ]() );
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
HybridizeProtocol::check_and_create_fragments( core::pose::Pose & pose ) {
	if (fragments9_ && fragments3_) return;

	if (!fragments9_) {
		fragments9_ = new core::fragment::ConstantLengthFragSet( 9 );

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
	
		// sequence
		std::string tgt_seq = pose.sequence();
		std::string tgt_ss(nres_tgt, '0');

		// templates vote on secstruct
		for (core::Size i=1; i<=templates_.size(); ++i) {
			for (core::Size j=1; j<=templates_[i]->total_residue(); ++j ) {
				core::Size tgt_pos = templates_[i]->pdb_info()->number(j);
				char tgt_ss_j = templates_[i]->secstruct(j);

				if (tgt_ss[tgt_pos] == '0') {
					tgt_ss[tgt_pos] = tgt_ss_j;
				} else if (tgt_ss[tgt_pos] != tgt_ss_j) {
					tgt_ss[tgt_pos] = 'D'; // templates disagree
				}
			}
		}
		for ( core::Size j=1; j<=nres_tgt; ++j ) {
			if (tgt_ss[j] == '0') tgt_ss[j] = 'D';
		}

		// pick from vall based on template SS + target sequence
		for ( core::Size j=1; j<=nres_tgt-8; ++j ) {
			std::string ss_sub = tgt_ss.substr( j-1, 9 );
			std::string aa_sub = tgt_seq.substr( j-1, 9 );

			core::fragment::FrameOP frame = new core::fragment::Frame( j, 9 );
			frame->add_fragment( 
				core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa( ss_sub, aa_sub, 25, true, core::fragment::IndependentBBTorsionSRFD() ) );
			fragments9_->add( frame );
		}
	}
	if (!fragments3_) {
		fragments3_ = new core::fragment::ConstantLengthFragSet( 3 );

		// make them from 9mers
		core::fragment::chop_fragments( *fragments9_, *fragments3_ );
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
	for (int i=1; i<=chosen_templ->total_residue(); ++i)
		for (int j=1; j<=chosen_templ->residue(i).natoms(); ++j) {
			core::id::AtomID src(j,i), tgt(j, chosen_templ->pdb_info()->number(i));
			pose.set_xyz( tgt, chosen_templ->xyz( src ) );
		}

	// make loops as inverse of template_contigs_icluster
	TR << "CONTIGS" << std::endl << template_contigs_icluster << std::endl;
	core::Size ncontigs = template_contigs_icluster.size();
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

		core::Size neffcycles = (core::Size)(1000*option[cm::hybridize::stage1_increase_cycles]());
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



void HybridizeProtocol::add_template(
	std::string template_fn,
	std::string cst_fn,
	std::string symm_file,
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
	protocols::loops::Loops contigs = protocols::loops::extract_continuous_chunks(*template_pose); 
	protocols::loops::Loops chunks = protocols::loops::extract_secondary_structure_chunks(*template_pose, "HE", 3, 6, 3, 4);

	if (chunks.num_loop() == 0)
		chunks = contigs;

	TR.Debug << "Chunks from template\n" << chunks << std::endl;
	TR.Debug << "Contigs from template\n" << contigs << std::endl;

	template_fn_.push_back(template_fn);
	templates_.push_back(template_pose);
	template_cst_fn_.push_back(cst_fn);
	symmdef_files_.push_back(symm_file);
	template_weights_.push_back(weight);
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
			add_template(template_fn, cst_fn, "", weight, cluster_id, cst_reses);
		}
	}
	f_stream.close();
}

void HybridizeProtocol::read_template_structures(utility::vector1 < utility::file::FileName > const & template_filenames)
{
	templates_.clear();
	templates_.resize(template_filenames.size());

	//core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

	//fpd remember sidechains from input structure
	core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

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
	using namespace ObjexxFCL::fmt;

	// make fragments if we don't have them at this point
	check_and_create_fragments( pose );

	// set pose for density scoring if a map was input
	//fpd eventually this should be moved to after pose initialization
	if ( option[ OptionKeys::edensity::mapfile ].user() ) {
		MoverOP dens( new protocols::electron_density::SetupForDensityScoringMover );
		dens->apply( pose );
	}

	// starting structure pool
	std::vector < SilentStructOP > post_centroid_structs;
	bool need_more_samples = true;

	// number of residues in asu without VRTs
	core::Size nres_tgt = pose.total_residue();
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres_tgt = symm_info->num_independent_residues();
	}
	if (pose.residue(nres_tgt).aa() == core::chemical::aa_vrt) nres_tgt--;
	
	core::Real gdtmm = 0.0;
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

		// realign each template to the starting template by domain
		if (realign_domains_) {
			// domain parsing
			DDomainParse ddom;
			utility::vector1< utility::vector1< loops::Loops > > domains_all_templ;
			domains_all_templ.resize( templates_.size() );
			for (Size i_template=1; i_template<=templates_.size(); ++i_template) {
				if (template_clusterID_[i_template] != template_clusterID_[initial_template_index]) continue;
				domains_all_templ[i_template] = ddom.split( *templates_[i_template], nres_tgt );
				
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
				for (int i=1; i<=domains_all_templ[i_template].size(); ++i) {
					TR << "domain " << i << ": " << domains_all_templ[i_template][i] << std::endl;
				}
			}
			
			// combine domains that are not in the initial template
			utility::vector1< loops::Loops > domains = expand_domains_to_full_length(domains_all_templ, initial_template_index, nres_tgt);
			TR << "Final decision: " << domains.size() << " domains" << std::endl;
			for (int i=1; i<= domains.size(); ++i) {
				TR << "domain " << i << ": " << domains[i] << std::endl;
			}
			
			// local align
			align_by_domain(templates_, domains, initial_template_index);
			
			// update chunk, contig informations
			for (Size i_template=1; i_template<=templates_.size(); ++i_template) {
				template_contigs_[i_template] = protocols::loops::extract_continuous_chunks(*templates_[i_template]); 
				template_chunks_[i_template] = protocols::loops::extract_secondary_structure_chunks(*templates_[i_template], "HE", 3, 6, 3, 4);
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
		}

		// initialize template history
		//fpd this must be done after symmetrization!
		TemplateHistoryOP history = new TemplateHistory(pose);
		history->setall( initial_template_index_icluster );
		pose.data().set( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY, history );

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
   		std::string cst_fn = template_cst_fn_[initial_template_index];
	
			FoldTreeHybridizeOP ft_hybridize(
				new FoldTreeHybridize(
					initial_template_index_icluster, templates_icluster, weights_icluster,
					template_chunks_icluster, template_contigs_icluster, fragments3_, fragments9_) ) ;
			ft_hybridize->set_constraint_file( cst_fn );
			ft_hybridize->set_scorefunction( stage1_scorefxn_ );
			ft_hybridize->set_increase_cycles( stage1_increase_cycles_ );
			ft_hybridize->set_add_non_init_chunks( add_non_init_chunks_ );
			ft_hybridize->set_frag_weight_aligned( frag_weight_aligned_ );
			ft_hybridize->set_max_registry_shift( max_registry_shift_ );
			ft_hybridize->apply(pose);
		} else {
			// just do frag insertion in unaligned regions
			core::pose::PoseOP chosen_templ = templates_icluster[initial_template_index_icluster];
			protocols::loops::Loops chosen_contigs = template_contigs_icluster[initial_template_index_icluster];
			initialize_and_sample_loops(pose, chosen_templ, chosen_contigs, stage1_scorefxn_);
		}

		//write gdtmm to output
		gdtmm = get_gdtmm(pose);
		core::pose::setPoseExtraScores( pose, "GDTMM_after_stage1", gdtmm);
		TR << "GDTMM_after_stage1" << F(8,3,gdtmm) << std::endl;

		// apply constraints
		std::string cst_fn = template_cst_fn_[initial_template_index];
		setup_centroid_constraints( pose, templates_, template_weights_, cst_fn );

		if (!option[cm::hybridize::skip_stage2]()) {
			CartesianHybridizeOP cart_hybridize (
				new CartesianHybridize(
					templates_icluster, weights_icluster, 
					template_chunks_icluster,template_contigs_icluster, fragments9_ ) );
			cart_hybridize->set_scorefunction( stage2_scorefxn_);
			cart_hybridize->set_increase_cycles( stage2_increase_cycles_ );
			cart_hybridize->set_no_global_frame( no_global_frame_ );
			cart_hybridize->set_linmin_only( linmin_only_ );
			cart_hybridize->set_cartfrag_overlap( cartfrag_overlap_ );
			bool linbonded_old = option[ score::linear_bonded_potential ]();
			option[ score::linear_bonded_potential ].value( true );  //fpd hack
			cart_hybridize->apply(pose);
			option[ score::linear_bonded_potential ].value( linbonded_old ); //fpd hack
		}

		//write gdtmm to output
		gdtmm = get_gdtmm(pose);
		core::pose::setPoseExtraScores( pose, "GDTMM_after_stage2", gdtmm);
		TR << "GDTMM_after_stage2" << F(8,3,gdtmm) << std::endl;

		// get fragment history
		runtime_assert( pose.data().has( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY ) );
		history = *( static_cast< TemplateHistory* >( pose.data().get_ptr( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY )() ));

		TR << "History :";
		for (int i=1; i<= history->size(); ++i ) { TR << I(4,i); }
		TR << std::endl;
		TR << "History :";
		for (int i=1; i<= history->size(); ++i ) { TR << I(4, history->get(i)); }
		TR << std::endl;

		// stage "2.5" .. minimize with centroid energy + full-strength cart bonded
		if (!option[cm::hybridize::skip_stage2]()) {
		 	core::optimization::MinimizerOptions options_lbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
		 	core::optimization::CartesianMinimizer minimizer;
		 	core::scoring::ScoreFunctionOP stage2_scorefxn_copy = stage2_scorefxn_->clone();
		 	core::Real fa_cart_bonded_wt = fa_scorefxn_->get_weight(core::scoring::cart_bonded);
		 	if (fa_cart_bonded_wt == 0) fa_cart_bonded_wt = 0.5;
		 	stage2_scorefxn_copy->set_weight( core::scoring::cart_bonded, fa_cart_bonded_wt );
		 	core::kinematics::MoveMap mm;
		 	mm.set_bb( true ); mm.set_chi( true ); mm.set_jump( true );
			if (core::pose::symmetry::is_symmetric(pose) )
				core::pose::symmetry::make_symmetric_movemap( pose, mm );
		 	options_lbfgs.max_iter(200);
		 	(*stage2_scorefxn_copy)(pose); minimizer.run( pose, mm, *stage2_scorefxn_copy, options_lbfgs );
		}

		// optional relax
		if (batch_relax_ > 0) {
			protocols::moves::MoverOP tofa = new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD );
			tofa->apply(pose);

			// apply fa constraints
			std::string cst_fn = template_cst_fn_[initial_template_index];
			setup_fullatom_constraints( pose, templates_, template_weights_, fa_cst_fn_, cst_fn );

			if (batch_relax_ == 1) {
				// add additional _CALPHA_ constraints
				//fpd  generally this is unused
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

				// standard relax
				// do relax _without_ ramping down coordinate constraints
				protocols::relax::FastRelax relax_prot( fa_scorefxn_, relax_repeats_ ,"NO CST RAMPING" );
				relax_prot.set_min_type("lbfgs_armijo_nonmonotone");
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
					relax_prot.set_min_type("lbfgs_armijo_nonmonotone");
					relax_prot.set_force_nonideal(true);
					relax_prot.set_script_to_batchrelax_default( relax_repeats_ );

					// need to use a packer task factory to handle poses with different disulfide patterning
					core::pack::task::TaskFactoryOP tf = new core::pack::task::TaskFactory;
					tf->push_back( new core::pack::task::operation::InitializeFromCommandline );
					tf->push_back( new core::pack::task::operation::IncludeCurrent );
					tf->push_back( new core::pack::task::operation::RestrictToRepacking );
					//tf->push_back( new core::pack::task::operation::NoRepackDisulfides );
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

	gdtmm = get_gdtmm(pose);
	core::pose::setPoseExtraScores( pose, "GDTMM_final", gdtmm);
	TR << "GDTMM_final" << F(8,3,gdtmm) << std::endl;
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
HybridizeProtocol::align_by_domain(utility::vector1<core::pose::PoseOP> & poses, utility::vector1 < Loops > domains, core::Size const ref_index)
{
	for (Size i_pose=1; i_pose <= poses.size(); ++i_pose) {
		if (i_pose == ref_index) continue;
		align_by_domain(*poses[i_pose], *poses[ref_index], domains);
		
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
			for (core::Size iloop=1; iloop<=domains[i_domain].num_loop(); ++iloop) {
				if ( pose.pdb_info()->number(ires) < domains[i_domain][iloop].start() || pose.pdb_info()->number(ires) > domains[i_domain][iloop].stop() ) continue;
				
				residue_list.push_back(ires);
			
				for (core::Size jres=1; jres<=ref_pose.total_residue(); ++jres) {
					if ( pose.pdb_info()->number(ires) != ref_pose.pdb_info()->number(jres) ) continue;
					if ( !pose.residue_type(ires).is_protein() ) continue;
					core::id::AtomID const id1( pose.residue_type(ires).atom_index("CA"), ires );
					core::id::AtomID const id2( ref_pose.residue_type(jres).atom_index("CA"), jres );

					atom_map[ id1 ] = id2;
					++n_mapped_residues;
				}
			}
		}
		if (n_mapped_residues >= 6) {
			partial_align(pose, ref_pose, atom_map, residue_list, true); // iterate_convergence = true
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
	utility::tag::TagPtr const tag, moves::DataMap & data, filters::Filters_map const & , moves::Movers_map const & , core::pose::Pose const & pose )
{
	// config file
	if( tag->hasOption( "config_file" ) )
		read_template_structures( tag->getOption< std::string >( "config_file" ) );

	// basic options
	stage1_increase_cycles_ = tag->getOption< core::Real >( "stage1_increase_cycles", 1. );
	stage2_increase_cycles_ = tag->getOption< core::Real >( "stage2_increase_cycles", 1. );
	fa_cst_fn_ = tag->getOption< std::string >( "fa_cst_file", "" );
	batch_relax_ = tag->getOption< core::Size >( "batch" , 1 );

	if( tag->hasOption( "starting_template" ) ) {
		std::vector<std::string> buff = utility::string_split( tag->getOption<std::string>( "starting_template" ), ',' );
		foreach(std::string field, buff){
			Size const value = std::atoi( field.c_str() ); // convert to C string, then convert to integer, then set a Size (phew!)
			starting_templates_.push_back(value);
		}
	}
	
	if( tag->hasOption( "stage1_probability" ) )
		stage1_probability_ = tag->getOption< core::Real >( "stage1_probability" );
	if( tag->hasOption( "realign_domains" ) )
		realign_domains_ = tag->getOption< bool >( "realign_domains" );
	if( tag->hasOption( "add_non_init_chunks" ) )
		add_non_init_chunks_ = tag->getOption< bool >( "add_non_init_chunks" );
	if( tag->hasOption( "frag_weight_aligned" ) )
		frag_weight_aligned_ = tag->getOption< core::Real >( "frag_weight_aligned" );
	if( tag->hasOption( "max_registry_shift" ) )
		max_registry_shift_ = tag->getOption< core::Size >( "max_registry_shift" );
	if( tag->hasOption( "no_global_frame" ) )
		no_global_frame_ = tag->getOption< bool >( "no_global_frame" );
	if( tag->hasOption( "linmin_only" ) )
		linmin_only_ = tag->getOption< bool >( "linmin_only" );
	if( tag->hasOption( "repeats" ) )
		relax_repeats_ = tag->getOption< core::Size >( "repeats" );
	if( tag->hasOption( "cartfrag_overlap" ) )
		cartfrag_overlap_ = tag->getOption< core::Size >( "cartfrag_overlap" );

	// scorfxns
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

	// fragments
	utility::vector1< utility::tag::TagPtr > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagPtr >::const_iterator tag_it;
	for (tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it) {
		if ( (*tag_it)->getName() == "Fragments" ) {
			using namespace core::fragment;

			fragments3_ = new ConstantLengthFragSet( 3 );
			fragments3_ = FragmentIO().read_data( (*tag_it)->getOption<std::string>( "3mers" )  );
			fragments9_ = new ConstantLengthFragSet( 9 );
			fragments9_ = FragmentIO().read_data( (*tag_it)->getOption<std::string>( "9mers" ) );
		}

		if ( (*tag_it)->getName() == "Template" ) {
			std::string template_fn = (*tag_it)->getOption<std::string>( "pdb" );
			std::string cst_fn = (*tag_it)->getOption<std::string>( "cst_file", "AUTO" );
			core::Real weight = (*tag_it)->getOption<core::Real>( "weight", 1 );
			core::Size cluster_id = (*tag_it)->getOption<core::Size>( "cluster_id", 1 );
			utility::vector1<core::Size> cst_reses;
			if ((*tag_it)->hasOption( "constrain_res" ))
				 cst_reses = protocols::rosetta_scripts::get_resnum_list_ordered( (*tag_it)->getOption<std::string>("constrain_res"), pose );

			std::string symm_file = (*tag_it)->getOption<std::string>( "symmdef", "" );

			add_template(template_fn, cst_fn, symm_file, weight, cluster_id, cst_reses);
		}

	}
}

} // hybridize 
} // comparative_modeling 
} // protocols

