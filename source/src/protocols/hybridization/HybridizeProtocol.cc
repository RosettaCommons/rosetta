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
#include <protocols/loops/Loop.hh>

#include <protocols/rigid/RB_geometry.hh>

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
#include <core/io/CrystInfo.hh>
#include <core/import_pose/import_pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/util.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/func/CircularSplineFunc.hh>
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

#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>

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

static THREAD_LOCAL basic::Tracer TR( "protocols.hybridization.HybridizeProtocol" );

namespace protocols {
namespace hybridization {


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
	return protocols::moves::MoverOP( new HybridizeProtocol );
}

std::string
HybridizeProtocolCreator::mover_name() {
	return "Hybridize";
}


/////////////
// mover
HybridizeProtocol::HybridizeProtocol()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	init();
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
	stage2_temperature_ = option[cm::hybridize::stage2_temperature]();
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
	skip_long_min_ = true;   //fpd  this is no longer necessary and seems to hurt model accuracy
	keep_pose_constraint_ = false;   //fpd PLEASE INITIALIZE NEW VARIABLES

	cenrot_ = option[corrections::score::cenrot]();

	csts_from_frags_ = false;  // generate dihedral constraints from fragments
	max_contig_insertion_ = -1;  // don't insert contigs larger than this size (-1 ==> don't limit)
	min_after_stage1_ = false;   // tors min after stage1
	fragprob_stage2_ = 0.3;  // ratio of fragment vs. template moves
	randfragprob_stage2_ = 0.5; // given a fragmove, how often is it applied to a random position (as opposed to a chainbreak position)


	// domain parsing options
	hcut_ = 0.18;
	pcut_ = 0.81;
	length_ = 38;


	jump_move_ = false;
	jump_move_repeat_ = 1;

	stage2_increase_cycles_ = option[cm::hybridize::stage2_increase_cycles]();
	stage25_increase_cycles_ = option[cm::hybridize::stage2min_increase_cycles]();
	no_global_frame_ = option[cm::hybridize::no_global_frame]();
	linmin_only_ = option[cm::hybridize::linmin_only]();

	// default scorefunctions
	//    stage2 scorefunctions are initialized in CartesianHybridize
	stage1_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	fa_scorefxn_ = core::scoring::get_score_function();

	core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *fa_scorefxn_ );


	if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
		utility::vector1< std::string > cst_files = option[ OptionKeys::constraints::cst_fa_file ]();
		core::Size choice = core::Size( numeric::random::rg().random_range( 1, cst_files.size() ) );
		fa_cst_fn_ = cst_files[choice];
		TR.Info << "Fullatom Constraint choice: " << fa_cst_fn_ << std::endl;
	}

	batch_relax_ = option[ cm::hybridize::relax ]();
	relax_repeats_ = option[ basic::options::OptionKeys::relax::default_repeats ]();

	// disulfide file
	if ( option[ in::fix_disulf ].user() ) {
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
		native_ = core::pose::PoseOP( new core::pose::Pose );
		if ( option[in::file::fullatom]() ) {
			core::import_pose::pose_from_file( *native_, option[ in::file::native ]() , core::import_pose::PDB_file);
		} else {
			core::chemical::ResidueTypeSetCOP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
			core::import_pose::pose_from_file( *native_, *residue_set, option[ in::file::native ]()  , core::import_pose::PDB_file);
		}
	} else if ( option[ evaluation::align_rmsd_target ].user() ) {
		native_ = core::pose::PoseOP( new core::pose::Pose );
		utility::vector1< std::string > const & align_rmsd_target( option[ evaluation::align_rmsd_target ]() );
		core::import_pose::pose_from_file( *native_, align_rmsd_target[1] , core::import_pose::PDB_file); // just use the first one for now
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
	if ( fragments_big_.size() > 0 && fragments_small_.size() > 0 ) return;

	core::Size fragbiglen=9;

	if ( !fragments_big_.size() ) {
		// number of residues
		core::Size nres_tgt = pose.total_residue();
		core::conformation::symmetry::SymmetryInfoCOP symm_info;
		if ( core::pose::symmetry::is_symmetric(pose) ) {
			core::conformation::symmetry::SymmetricConformation & SymmConf (
				dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
			symm_info = SymmConf.Symmetry_Info();
			nres_tgt = symm_info->num_independent_residues();
		}
		if ( pose.residue(nres_tgt).aa() == core::chemical::aa_vrt ) nres_tgt--;
		while ( !pose.residue(nres_tgt).is_protein() ) nres_tgt--;

		fragbiglen = std::min( fragbiglen, nres_tgt );

		core::fragment::FragSetOP frags( new core::fragment::ConstantLengthFragSet(  ) );

		// sequence
		std::string tgt_seq = pose.sequence();
		std::string tgt_ss(nres_tgt, '0');

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		if ( option[ OptionKeys::in::file::psipred_ss2 ].user() ) {
			utility::vector1< char > psipred = read_psipred_ss2_file( pose );
			for ( core::Size j=1; j<=nres_tgt; ++j ) {
				tgt_ss[j-1] = psipred[j];
			}
		} else {
			// templates vote on secstruct
			for ( core::Size i=1; i<=templates_.size(); ++i ) {
				for ( core::Size j=1; j<=templates_[i]->total_residue(); ++j ) {
					if ( !templates_[i]->residue(j).is_protein() ) continue;
					core::Size tgt_pos = templates_[i]->pdb_info()->number(j);

					runtime_assert( tgt_pos<=nres_tgt );
					char tgt_ss_j = templates_[i]->secstruct(j);

					if ( tgt_ss[tgt_pos-1] == '0' ) {
						tgt_ss[tgt_pos-1] = tgt_ss_j;
					} else if ( tgt_ss[tgt_pos-1] != tgt_ss_j ) {
						tgt_ss[tgt_pos-1] = 'D'; // templates disagree
					}
				}
			}
			for ( core::Size j=1; j<=nres_tgt; ++j ) {
				if ( tgt_ss[j-1] == '0' ) tgt_ss[j-1] = 'D';
			}
		}

		// pick from vall based on template SS + target sequence
		for ( core::Size j=1; j<=nres_tgt-fragbiglen+1; ++j ) {
			core::fragment::FrameOP frame( new core::fragment::Frame( j, fragbiglen ) );

			if ( j > residue_sample_abinitio_.size() || residue_sample_abinitio_[j] ) {
				frame->add_fragment(
					core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa(
					tgt_ss.substr( j-1, fragbiglen ), tgt_seq.substr( j-1, fragbiglen ), 25, true, core::fragment::IndependentBBTorsionSRFD() ) );
				frags->add( frame );
			}
		}
		fragments_big_.push_back( frags );
	}
	if ( !fragments_small_.size() ) {
		core::fragment::FragSetOP frags( new core::fragment::ConstantLengthFragSet( 3 ) );

		// make them from big fragments
		core::fragment::chop_fragments( *fragments_big_[1], *frags );
		fragments_small_.push_back( frags );
	}
}

//fpd add fragment-derived constraints
void
HybridizeProtocol::add_fragment_csts( core::pose::Pose &pose ) {
	TR << "Adding fragment constraints" << std::endl;

	core::fragment::FragSetOP frags = fragments_small_[1];

	// stolen from SecondaryStructure.cc
	if ( frags->global_offset() != 0 ) {
		TR.Error << "[ERROR] SecondaryStructure computations must be carried out with local coordinates (global offset of fragments must be 0)." << std::endl;
		runtime_assert( false );
	}

	// 1 - collect fragment statistics
	Size frag_nres = frags->max_pos();
	utility::vector1< utility::vector1< core::Real > > phi_distr( frag_nres, utility::vector1< core::Real >(36,0.0) );
	utility::vector1< utility::vector1< core::Real > > psi_distr( frag_nres, utility::vector1< core::Real >(36,0.0) );
	utility::vector1< int > N( frag_nres, 0 );

	for ( core::fragment::FragID_Iterator it=frags->begin(), eit=frags->end(); it!=eit; ++it ) { //carefully checked that I don't change FrameData
		core::Size loop_start = 1;
		core::Size loop_end = it->frame().length();
		for ( core::Size fpos = loop_start; fpos <= loop_end; ++fpos ) {
			core::fragment::BBTorsionSRFDCOP res_i =
				utility::pointer::dynamic_pointer_cast<const core::fragment::BBTorsionSRFD> (it->fragment().get_residue( fpos ) );

			Size pos = it->frame().seqpos( fpos );
			core::Real phi = std::fmod( res_i->torsion(1), 360.0 ); if ( phi<0 ) phi +=360.0;
			core::Real psi = std::fmod( res_i->torsion(2), 360.0 ); if ( psi<0 ) psi +=360.0;

			core::Size phibin = (core::Size) std::floor( phi/10.0 ); if ( phibin == 36 ) phibin=0;
			core::Size psibin = (core::Size) std::floor( psi/10.0 ); if ( psibin == 36 ) psibin=0;

			phi_distr[pos][phibin+1]+=1.0;
			psi_distr[pos][psibin+1]+=1.0;
			N[pos]++;
		}
	}

	// smoothing
	//core::Real smoothing[7] = {0.09, 0.14, 0.175, 0.19, 0.175, 0.14, 0.09};  // very soft
	core::Real smoothing[7] = {0.01, 0.05, 0.24, 0.40, 0.24, 0.05, 0.01}; // less soft
	for ( int i=1; i<=(int)frag_nres; ++i ) {
		utility::vector1< core::Real > phi_i = phi_distr[i];
		utility::vector1< core::Real > psi_i = psi_distr[i];
		for ( int j=1; j<=36; ++j ) {
			phi_distr[i][j] = smoothing[0]*phi_i[1+((j+32)%36)]
				+ smoothing[1]*phi_i[1+((j+33)%36)]
				+ smoothing[2]*phi_i[1+((j+34)%36)]
				+ smoothing[3]*phi_i[1+((j+35)%36)]
				+ smoothing[4]*phi_i[1+((j+0)%36)]
				+ smoothing[5]*phi_i[1+((j+1)%36)]
				+ smoothing[6]*phi_i[1+((j+2)%36)];

			psi_distr[i][j] = smoothing[0]*psi_i[1+((j+32)%36)]
				+ smoothing[1]*psi_i[1+((j+33)%36)]
				+ smoothing[2]*psi_i[1+((j+34)%36)]
				+ smoothing[3]*psi_i[1+((j+35)%36)]
				+ smoothing[4]*psi_i[1+((j+0)%36)]
				+ smoothing[5]*psi_i[1+((j+1)%36)]
				+ smoothing[6]*psi_i[1+((j+2)%36)];

			phi_distr[i][j] /= N[i];
			psi_distr[i][j] /= N[i];
		}
	}

	// convert to energy
	core::Real mest=0.1;
	for ( int i=1; i<=(int)frag_nres; ++i ) {
		for ( int j=1; j<=36; ++j ) {
			// shift so max is 0
			phi_distr[i][j] = -std::log( mest/36.0 + (1-mest)*phi_distr[i][j] ) + std::log( mest/36.0 ) ;
			psi_distr[i][j] = -std::log( mest/36.0 + (1-mest)*psi_distr[i][j] ) + std::log( mest/36.0 ) ;
		}

		//TR << "phi " << i;
		//for (int j=1; j<=36; ++j) TR << " " << phi_distr[i][j];
		//TR << std::endl;

		//TR << "psi " << i;
		//for (int j=1; j<=36; ++j) TR << " " << psi_distr[i][j];
		//TR << std::endl;
	}

	// finally .. add constraints
	for ( int i=1; i<=(int)frag_nres; ++i ) {
		using namespace core::scoring::func;
		using namespace core::scoring::constraints;

		if ( !pose.residue(i).is_protein() ) continue;

		if ( i>1 && pose.residue(i-1).is_protein() ) {
			FuncOP phi_func( new CircularSplineFunc( 1.0, phi_distr[i] ) );
			ConstraintOP phi_cst( new DihedralConstraint(
				core::id::AtomID(3,i-1),core::id::AtomID(1,i),core::id::AtomID(2,i),core::id::AtomID(3,i), phi_func  ) );
			pose.add_constraint( scoring::constraints::ConstraintCOP( phi_cst ) );
		}

		if ( i<(int)frag_nres && pose.residue(i+1).is_protein() ) {
			FuncOP psi_func( new CircularSplineFunc( 1.0, psi_distr[i] ) );
			ConstraintOP psi_cst( new DihedralConstraint(
				core::id::AtomID(1,i),core::id::AtomID(2,i),core::id::AtomID(3,i),core::id::AtomID(1,i+1), psi_func  ) );
			pose.add_constraint( scoring::constraints::ConstraintCOP( psi_cst ) );
		}
	}
}

void
HybridizeProtocol::initialize_and_sample_loops(
	core::pose::Pose &pose,
	core::pose::PoseOP chosen_templ,
	protocols::loops::Loops template_contigs,
	core::scoring::ScoreFunctionOP scorefxn)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// xyz copy starting model
	for ( Size i=1; i<=chosen_templ->total_residue(); ++i ) {
		for ( Size j=1; j<=chosen_templ->residue(i).natoms(); ++j ) {
			core::id::AtomID src(j,i), tgt(j, chosen_templ->pdb_info()->number(i));
			pose.set_xyz( tgt, chosen_templ->xyz( src ) );
		}
	}

	// make loops as inverse of template_contigs
	TR << "CONTIGS" << std::endl << template_contigs << std::endl;
	//core::Size ncontigs = template_contigs.size();
	core::Size nres_tgt = pose.total_residue();

	//symmetry
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres_tgt = symm_info->num_independent_residues();
	}

	if ( pose.residue(nres_tgt).aa() == core::chemical::aa_vrt ) nres_tgt--;

	protocols::loops::LoopsOP loops( new protocols::loops::Loops );
	utility::vector1< bool > templ_coverage(nres_tgt, false);

	for ( Size i=1; i<=chosen_templ->total_residue(); ++i ) {
		core::Size cres = chosen_templ->pdb_info()->number(i);
		templ_coverage[cres] = true;
	}

	// remove 1-3 residue "segments"
	for ( Size i=1; i<=nres_tgt-2; ++i ) {
		if ( !templ_coverage[i] && templ_coverage[i+1] && !templ_coverage[i+2] ) {
			templ_coverage[i+1]=false;
		} else if ( i<=nres_tgt-3 && !templ_coverage[i] && templ_coverage[i+1] && templ_coverage[i+2] && !templ_coverage[i+3] ) {
			templ_coverage[i+1]=false;
			templ_coverage[i+2]=false;
		} else if ( i<=nres_tgt-4 &&
				!templ_coverage[i] && templ_coverage[i+1] && templ_coverage[i+2] && templ_coverage[i+3] && !templ_coverage[i+4] ) {
			templ_coverage[i+1]=false;
			templ_coverage[i+2]=false;
			templ_coverage[i+3]=false;
		}
	}
	// make loopfile
	bool inloop=!templ_coverage[1];
	core::Size loopstart=1, loopstop;
	for ( Size i=2; i<=nres_tgt; ++i ) {
		if ( templ_coverage[i] && inloop ) {
			inloop = false;
			loopstop = i;
			if ( loopstop < loopstart + 2 ) {
				if ( loopstart>1 ) loopstart--;
				if ( loopstop<nres_tgt ) loopstop++;
			}
			loops->add_loop( loopstart,loopstop );
		} else if ( !templ_coverage[i] && !inloop ) {
			// start loop
			inloop = true;
			loopstart = i-1;
		}
	}
	if ( inloop ) {
		// end loop
		loopstop = nres_tgt;
		while ( loopstop < loopstart + 2 ) { loopstart--; }
		loops->add_loop( loopstart,loopstop );
	}
	TR << "LOOPS" << std::endl << *loops << std::endl;

	// now do insertions
	if ( loops->size() != 0 ) {
		// set foldtree + variants
		loops->auto_choose_cutpoints(pose);
		core::kinematics::FoldTree f_in=pose.fold_tree(), f_new;
		protocols::loops::fold_tree_from_loops( pose, *loops, f_new);
		pose.fold_tree( f_new );
		protocols::loops::add_cutpoint_variants( pose );

		// set movemap
		core::kinematics::MoveMapOP mm_loop( new core::kinematics::MoveMap() );
		for ( protocols::loops::Loops::const_iterator it=loops->begin(), it_end=loops->end(); it!=it_end; ++it ) {
			for ( Size i=it->start(); i<=it->stop(); ++i ) {
				mm_loop->set_bb(i, true);
				mm_loop->set_chi(i, true); // chi of loop residues
			}
		}

		// setup fragment movers
		protocols::simple_moves::ClassicFragmentMoverOP frag3mover, frag9mover;
		core::fragment::FragSetOP frags_small = fragments_small_[numeric::random::rg().random_range(1,fragments_small_.size())];
		core::fragment::FragSetOP frags_big = fragments_big_[numeric::random::rg().random_range(1,fragments_big_.size())];
		TR.Info << "FRAGMENTS small max length: " << frags_small->max_frag_length() << std::endl;
		TR.Info << "FRAGMENTS big max length: " << frags_big->max_frag_length() << std::endl;
		frag3mover = protocols::simple_moves::ClassicFragmentMoverOP( new protocols::simple_moves::ClassicFragmentMover( frags_small, mm_loop ) );
		frag3mover->set_check_ss( false ); frag3mover->enable_end_bias_check( false );
		frag9mover = protocols::simple_moves::ClassicFragmentMoverOP( new protocols::simple_moves::ClassicFragmentMover( frags_big, mm_loop ) );
		frag9mover->set_check_ss( false ); frag9mover->enable_end_bias_check( false );

		// extend + idealize loops
		for ( protocols::loops::Loops::const_iterator it=loops->begin(), it_end=loops->end(); it!=it_end; ++it ) {
			protocols::loops::Loop to_idealize( *it );
			protocols::loops::set_extended_torsions( pose, *it );
		}

		// setup MC
		scorefxn->set_weight( core::scoring::linear_chainbreak, 0.5 );
		(*scorefxn)(pose);
		protocols::moves::MonteCarloOP mc1( new protocols::moves::MonteCarlo( pose, *scorefxn, 2.0 ) );

		core::Size neffcycles = (core::Size)(1000*option[cm::hybridize::stage1_increase_cycles]());
		for ( Size n=1; n<=neffcycles; ++n ) {
			frag9mover->apply( pose ); (*scorefxn)(pose); mc1->boltzmann( pose , "frag9" );
			frag3mover->apply( pose ); (*scorefxn)(pose); mc1->boltzmann( pose , "frag3" );

			if ( n%100 == 0 ) {
				mc1->show_scores();
				mc1->show_counters();
			}
		}
		mc1->recover_low( pose );

		scorefxn->set_weight( core::scoring::linear_chainbreak, 2.0 );
		(*scorefxn)(pose);
		protocols::moves::MonteCarloOP mc2( new protocols::moves::MonteCarlo( pose, *scorefxn, 2.0 ) );
		for ( Size n=1; n<=neffcycles; ++n ) {
			frag9mover->apply( pose ); (*scorefxn)(pose); mc2->boltzmann( pose , "frag9" );
			frag3mover->apply( pose ); (*scorefxn)(pose); mc2->boltzmann( pose , "frag3" );

			if ( n%100 == 0 ) {
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
	utility::vector1<char> randchains)
{
	core::chemical::ResidueTypeSetCOP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

	if ( template_fn == "extended" ) {
		//
		core::pose::PoseOP template_pose( new core::pose::Pose() );

		// auto constraints make no sense
		if ( cst_fn == "AUTO" ) {
			TR.Error << "Warning!  Turning off auto constraints for extended pose" << std::endl;
			cst_fn = "NONE";
		}

		add_null_template( template_pose, cst_fn, symm_file, weight );
	} else {
		core::pose::PoseOP template_pose( new core::pose::Pose() );
		core::import_pose::pose_from_file( *template_pose, *residue_set, template_fn , core::import_pose::PDB_file);

		add_template( template_pose, cst_fn, symm_file, weight, randchains, template_fn );
	}
}


void HybridizeProtocol::add_null_template(
	core::pose::PoseOP template_pose,
	std::string cst_fn,
	std::string symm_file,
	core::Real weight)
{
	template_fn_.push_back("null");
	templates_.push_back(template_pose);
	templates_aln_.push_back(template_pose); // shallow copy
	template_cst_fn_.push_back(cst_fn);
	symmdef_files_.push_back(symm_file);
	template_weights_.push_back(weight);
	template_chunks_.push_back(protocols::loops::Loops());
	template_contigs_.push_back(protocols::loops::Loops());
	randomize_chains_.push_back(utility::vector1<core::Size>(0));
}

void HybridizeProtocol::add_template(
	core::pose::PoseOP template_pose,
	std::string cst_fn,
	std::string symm_file,
	core::Real weight,
	utility::vector1<char> rand_chains,
	std::string filename)
{
	// add secondary structure information to the template pose
	core::scoring::dssp::Dssp dssp_obj( *template_pose );
	dssp_obj.insert_ss_into_pose( *template_pose );

	// find ss chunks in template
	protocols::loops::Loops contigs = protocols::loops::extract_continuous_chunks(*template_pose);
	protocols::loops::Loops chunks = protocols::loops::extract_secondary_structure_chunks(*template_pose, "HE", 3, 6, 3, 4);

	// if there are no SS elts in the pose, use contigs
	if ( chunks.num_loop() == 0 ) chunks = contigs;

	TR.Debug << "Chunks from template\n" << chunks << std::endl;
	TR.Debug << "Contigs from template\n" << contigs << std::endl;

	template_fn_.push_back(filename);
	templates_.push_back(template_pose);
	templates_aln_.push_back(template_pose); // shallow copy
	template_cst_fn_.push_back(cst_fn);
	symmdef_files_.push_back(symm_file);
	template_weights_.push_back(weight);
	template_chunks_.push_back(chunks);
	template_contigs_.push_back(contigs);
	randomize_chains_.push_back(rand_chains);

	non_null_template_indices_.push_back( templates_.size() );
}

void HybridizeProtocol::update_last_template()
{
	core::pose::PoseOP template_pose = templates_[templates_.size()];

	// add secondary structure information to the template pose
	core::scoring::dssp::Dssp dssp_obj( *template_pose );
	dssp_obj.insert_ss_into_pose( *template_pose );

	// find ss chunks in template
	protocols::loops::Loops contigs = protocols::loops::extract_continuous_chunks(*template_pose);
	protocols::loops::Loops chunks = protocols::loops::extract_secondary_structure_chunks(*template_pose, "HE", 3, 6, 3, 4);

	if ( chunks.num_loop() == 0 ) chunks = contigs;

	template_chunks_[templates_.size()] = chunks;
	template_contigs_[templates_.size()] = contigs;
}


// validate input templates match input sequence
// TO DO: if only sequences mismatch try realigning
void HybridizeProtocol::validate_template(
	std::string filename,
	std::string fasta,
	core::pose::PoseOP template_pose,
	bool & align_pdb_info)
{
	core::Size nres_fasta = fasta.length(), nres_templ = template_pose->total_residue();
	core::Size next_ligand = nres_fasta+1; // ligands must be sequentially numbered

	bool requires_alignment = false;
	for ( core::Size i=1; i<=nres_templ; ++i ) {
		char templ_aa = template_pose->residue(i).name1();
		core::Size i_fasta = template_pose->pdb_info()->number(i);

		// ligands
		if ( i_fasta > nres_fasta ) {
			bool number_correct = (i_fasta==next_ligand);
			bool is_ligand = !template_pose->residue(i).is_protein();
			if ( !(number_correct && is_ligand) ) {
				TR.Error << "Error in input numbering in template " << filename
					<< " while reading residue " << i << std::endl;
				if ( !is_ligand ) {
					TR.Error << "   Residue number " << i_fasta << " is outside of input fasta range." << std::endl;
					utility_exit();
				} else {
					TR.Error << "   Ligands must be numbered sequentially after protein residues." << std::endl;
					TR.Error << "   Expected: " << next_ligand << "   Saw: " << i_fasta << std::endl;
					utility_exit();
				}
			}
			next_ligand++;
		} else {
			char fasta_aa = fasta[i_fasta-1];
			if ( fasta_aa != templ_aa ) {
				TR.Error << "Sequence mismatch between input fasta and template " << filename
					<< " at residue " << i_fasta << std::endl;
				TR.Error << "   Expected: " << fasta_aa << "   Saw: " << templ_aa << std::endl;
				if ( !align_pdb_info ) utility_exit();
				requires_alignment = true;
				break;
			}
		}
	}
	if ( requires_alignment ) {
		TR.Error << "THE PDB INFO IS NOT IN SYNC WITH THE FASTA. ATTEMPTING TO RESOLVE THE ISSUE AUTOMATICALLY" << std::endl;
		//this block code shouldn't be needed since hybrid needs asymmetric poses anyway but it doesn't hurt
		core::pose::Pose pose_for_seq;
		core::pose::symmetry::extract_asymmetric_unit(*template_pose, pose_for_seq, false);

		core::sequence::SequenceOP full_length_seq( new core::sequence::Sequence( fasta, "target" ));
		core::sequence::SequenceOP t_pdb_seq( new core::sequence::Sequence( pose_for_seq.sequence(), "pose_seq" ));
		core::sequence::SWAligner sw_align;
		core::sequence::ScoringSchemeOP ss(  new core::sequence::SimpleScoringScheme(120, 0, -100, 0));
		core::sequence::SequenceAlignment fasta2template;

		fasta2template = sw_align.align(full_length_seq, t_pdb_seq, ss);

		TR << fasta2template << std::endl;


		core::id::SequenceMapping sequencemap = fasta2template.sequence_mapping(1,2);
		sequencemap.reverse();
		core::Size ndel = 0;
		core::Size currres = 1;
		core::Size nres = template_pose->total_residue();
		for ( Size i=1; i<=nres; i++ ) {
			Size pdbnumber = template_pose->pdb_info()->number(i);
			if ( sequencemap[i] == 0 ) { // extra residues in template
				template_pose->delete_residue_range_slow( currres,currres );
				ndel++;
				continue;
			}
			if ( pdbnumber != sequencemap[i] ) {
				Size fastanumber = sequencemap[i];
				template_pose->pdb_info()->number(currres,fastanumber);
			}
			currres++;
		}
		if ( ndel > 0 ) {
			TR.Error << "WARNING! Realignment removed " << ndel << " residues!" << std::endl;
		}
	} else {
		align_pdb_info = false;
	}
}


core::Size HybridizeProtocol::pick_starting_template( ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	numeric::random::WeightedSampler weighted_sampler;
	weighted_sampler.weights(template_weights_);
	return weighted_sampler.random_sample(numeric::random::rg());
}

void HybridizeProtocol::domain_parse_templates(core::Size nres) {
	if ( domains_all_templ_.size() == templates_.size() ) {
		return;
	}

	DDomainParse ddom(pcut_,hcut_,length_);

	domains_all_templ_.resize( templates_.size() );
	for ( Size i_template=1; i_template<=templates_.size(); ++i_template ) {
		if ( templates_[i_template]->total_residue() < 3 ) {
			continue;  // ???
		}

		domains_all_templ_[i_template] = ddom.split( *templates_[i_template] );

		// convert domain numbering to target pose numbering
		protocols::loops::Loops all_domains;
		for ( Size iloops=1; iloops<=domains_all_templ_[i_template].size(); ++iloops ) {
			for ( Size iloop=1; iloop<=domains_all_templ_[i_template][iloops].num_loop(); ++iloop ) {
				Size seqpos_start_pose = templates_[i_template]->pdb_info()->number(domains_all_templ_[i_template][iloops][iloop].start());

				if ( seqpos_start_pose > nres ) continue; // sometimes dna can do this

				domains_all_templ_[i_template][iloops][iloop].set_start( seqpos_start_pose );

				Size seqpos_stop_pose = templates_[i_template]->pdb_info()->number(domains_all_templ_[i_template][iloops][iloop].stop());
				domains_all_templ_[i_template][iloops][iloop].set_stop( seqpos_stop_pose );

				all_domains.add_loop( seqpos_start_pose, seqpos_stop_pose );
			}
		}

		// extend to cover full-length seq
		// assumes no overlap in domain defs!
		all_domains.sequential_order();
		if ( all_domains.num_loop() == 0 ) continue;

		std::map<core::Size,core::Size> remap_endpoints;
		remap_endpoints [ all_domains[1].start() ] = 1;
		remap_endpoints [ all_domains[all_domains.num_loop()].stop() ] = nres;
		for ( Size iloops=2; iloops<=all_domains.size(); ++iloops ) {
			core::Size start_i = (all_domains[iloops-1].stop() + all_domains[iloops].start() + 1)/2;
			remap_endpoints [ all_domains[iloops-1].stop() ] = start_i-1;
			remap_endpoints [ all_domains[iloops].start() ] = start_i;
		}

		for ( Size iloops=1; iloops<=domains_all_templ_[i_template].size(); ++iloops ) {
			for ( Size iloop=1; iloop<=domains_all_templ_[i_template][iloops].num_loop(); ++iloop ) {
				Size seqpos_start_pose = remap_endpoints[ domains_all_templ_[i_template][iloops][iloop].start() ];
				domains_all_templ_[i_template][iloops][iloop].set_start( seqpos_start_pose );

				Size seqpos_stop_pose = remap_endpoints[ domains_all_templ_[i_template][iloops][iloop].stop() ];
				domains_all_templ_[i_template][iloops][iloop].set_stop( seqpos_stop_pose );
			}
		}

		TR << "Found " << domains_all_templ_[i_template].size() << " domains using template " << template_fn_[i_template] << std::endl;
		for ( Size i=1; i<=domains_all_templ_[i_template].size(); ++i ) {
			TR << "domain " << i << ": " << domains_all_templ_[i_template][i] << std::endl;
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
	if ( pose.residue(nres_tgt).aa() == core::chemical::aa_vrt ) nres_tgt--;
	while ( !pose.residue(nres_protein_tgt).is_protein() ) nres_protein_tgt--;

	//save necessary constraint in pose
	core::scoring::constraints::ConstraintSetOP save_pose_constraint_set( new core::scoring::constraints::ConstraintSet() ) ;
	if ( keep_pose_constraint_ ) {
		save_pose_constraint_set = pose.constraint_set()->clone();
		core::scoring::constraints::remove_nonbb_constraints(pose);
	}

	if ( pose.is_fullatom() ) {
		protocols::moves::MoverOP tocen( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );
		tocen->apply(pose);
	}

	// make fragments if we don't have them at this point
	check_and_create_fragments( pose );

	// domain parse templates if we do not have domain definitions
	domain_parse_templates(nres_tgt);

	// select random fragment sets if >2 are provided
	core::fragment::FragSetOP frags_small = fragments_small_[numeric::random::rg().random_range(1,fragments_small_.size())];
	core::fragment::FragSetOP frags_big = fragments_big_[numeric::random::rg().random_range(1,fragments_big_.size())];
	TR.Info << "FRAGMENTS small max length: " << frags_small->max_frag_length() << std::endl;
	TR.Info << "FRAGMENTS big max length: " << frags_big->max_frag_length() << std::endl;

	// starting structure pool (for batch relax)
	std::vector < SilentStructOP > post_centroid_structs;
	bool need_more_samples = true;

	// main centroid generation loop
	core::Real gdtmm = 0.0;
	while ( need_more_samples ) {
		need_more_samples = false;

		// pick starting template
		core::Size initial_template_index = pick_starting_template();
		TR << "Using initial template: " << I(4,initial_template_index) << " " << template_fn_[initial_template_index] << std::endl;

		// ensure
		//    1)that no CONTIGS cross multiple CHAINS
		//    2)that every input CHAIN has a chunk in it
		// if not, add the corresponding contig as a chunk
		utility::vector1< int > cuts = pose.fold_tree().cutpoints();
		protocols::loops::Loops const &contigs_old = template_contigs_[initial_template_index];
		protocols::loops::Loops contigs_new;
		for ( core::Size j=1; j<=contigs_old.num_loop(); ++j ) {
			protocols::loops::Loop contig_j = contigs_old[j];

			for ( core::Size i=1; i<=cuts.size(); ++i ) {
				int nextcut = cuts[i];
				if ( nextcut > (int)nres_protein_tgt ) break;

				core::Size start = templates_[initial_template_index]->pdb_info()->number(contig_j.start());
				core::Size stop  = templates_[initial_template_index]->pdb_info()->number(contig_j.stop());

				if ( (int)start <= nextcut && (int)stop > nextcut ) {
					// find template res corresponding to cut
					for ( core::Size k=contig_j.start(); k<=contig_j.stop(); ++k ) {
						if ( templates_[initial_template_index]->pdb_info()->number(k) == nextcut ) {
							// ensure contig >= 3 res
							if ( k > contig_j.start()+1 ) {
								contigs_new.push_back( protocols::loops::Loop(contig_j.start(), k) );
							}
							TR << "Split contig ("<<contig_j.start()<<","<<contig_j.stop()<< ") at " << k << " [pdbnum: "<<start<<"/"<<stop<<"/" << nextcut << "]" << std::endl;
							contig_j.set_start( k+1 );
						}
					}

				}
			}
			// ensure contig >= 3 res
			if ( contig_j.stop() > contig_j.start()+1 ) {
				contigs_new.push_back( contig_j );
			}
		}
		template_contigs_[initial_template_index] = contigs_new;

		protocols::loops::Loops &chunks = template_chunks_[initial_template_index];
		for ( core::Size i=0; i<=cuts.size(); ++i ) {
			int prevcut = (i==0) ? 1 : cuts[i];
			int nextcut = (i==cuts.size()) ? nres_protein_tgt : cuts[i+1];
			if ( nextcut > (int)nres_protein_tgt ) break;
			bool haschunk = false;
			for ( core::Size j=1; j<=chunks.num_loop() && !haschunk; ++j ) {
				protocols::loops::Loop const &chunk_j = chunks[j];
				core::Size start = templates_[initial_template_index]->pdb_info()->number(chunk_j.start());
				core::Size stop  = templates_[initial_template_index]->pdb_info()->number(chunk_j.stop());
				if ( (int)start > prevcut && (int)stop <= nextcut ) haschunk=true;
			}

			if ( ! haschunk ) {
				TR << "Segment ("<<prevcut<<","<<nextcut<<") has no chunks!  Adding contigs!" << std::endl;
				bool hascontig=false;
				for ( core::Size j=1; j<=contigs_new.num_loop(); ++j ) {
					protocols::loops::Loop const &contig_j = contigs_new[j];
					core::Size start = templates_[initial_template_index]->pdb_info()->number(contig_j.start());
					core::Size stop  = templates_[initial_template_index]->pdb_info()->number(contig_j.stop());
					if ( (int)start > prevcut && (int)stop <= nextcut ) {
						hascontig = true;
						chunks.push_back( contig_j );
					}
				}
				if ( !hascontig ) {
					TR << "Warning!  No contigs found for segment!" << std::endl;
					TR << "If you are not using the 'randomize=X' option there is likely something wrong with your input and models will not be reasonable!" << std::endl;
				}
			}
		}

		// (0) randomize chains
		if ( randomize_chains_[initial_template_index].size() > 0 ) {
			runtime_assert ( templates_[initial_template_index]->pdb_info() );
			numeric::xyzVector< core::Real > comFixed(0,0,0), comMoving(0,0,0);
			core::Real nFixed=0, nMoving=0;

			for ( core::Size i=1; i<=templates_[initial_template_index]->total_residue(); ++i ) {
				char chain = templates_[initial_template_index]->pdb_info()->chain(i);
				if ( std::find(
						randomize_chains_[initial_template_index].begin(),
						randomize_chains_[initial_template_index].end(),
						chain) != randomize_chains_[initial_template_index].end()
						) {
					comMoving += templates_[initial_template_index]->residue(i).xyz(2);
					nMoving++;
				} else {
					comFixed += templates_[initial_template_index]->residue(i).xyz(2);
					nFixed++;
				}
			}
			comMoving /= nMoving;
			comFixed /= nFixed;

			// randomize orientation
			numeric::xyzMatrix< core::Real > R = protocols::geometry::random_reorientation_matrix( 360, 360 );

			// randomize offset
			numeric::xyzVector< core::Real > T = numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() );

			// perturb & slide into contact
			bool done=false, compacting=false;
			core::Real DIST = 100.0;
			core::scoring::ScoreFunction sf_vdw;
			core::Real vdw_score = -1.0, step_size = -2.0, min_step_size = 0.25;
			sf_vdw.set_weight( core::scoring::vdw, 1.0 );
			core::pose::Pose template_orig = *(templates_[initial_template_index]);
			while ( !done && DIST > 0 ) {
				utility::vector1< id::AtomID > atm_ids;
				utility::vector1< numeric::xyzVector< core::Real> > atm_xyzs;

				for ( core::Size i=1; i<=templates_[initial_template_index]->total_residue(); ++i ) {
					core::conformation::Residue const & rsd_i = template_orig.residue(i);
					char chain = templates_[initial_template_index]->pdb_info()->chain(i);
					if ( std::find(
							randomize_chains_[initial_template_index].begin(),
							randomize_chains_[initial_template_index].end(),
							chain) != randomize_chains_[initial_template_index].end()
							) {
						for ( core::Size j = 1; j <= rsd_i.natoms(); ++j ) {
							id::AtomID id( j,i );
							atm_ids.push_back( id );
							numeric::xyzVector< core::Real > new_ij = R*(rsd_i.xyz(j) - comMoving) + comFixed + DIST * T;
							atm_xyzs.push_back( new_ij );
						}
					}
				}
				templates_[initial_template_index]->batch_set_xyz( atm_ids, atm_xyzs );
				core::Real score_d = sf_vdw( *(templates_[initial_template_index]) );
				if ( vdw_score < 0 ) vdw_score = score_d;

				if ( !compacting && score_d > vdw_score + 2.0 ) {
					compacting = true;
					step_size = min_step_size;
				} else if ( compacting  && vdw_score <= score_d + 2.0 ) {
					done = true;
				}

				TR << "D=" << DIST << " score=" << score_d << std::endl;

				DIST += step_size;
			}
		}

		// (1) steal hetatms from template
		utility::vector1< std::pair< core::Size,core::Size > > hetatms;
		if ( add_hetatm_ ) {
			for ( Size ires=1; ires <= templates_[initial_template_index]->total_residue(); ++ires ) {
				if ( templates_[initial_template_index]->pdb_info()->number(ires) > (int)nres_tgt ) {
					TR.Debug << "Insert hetero residue: " << templates_[initial_template_index]->residue(ires).name3() << std::endl;
					if ( templates_[initial_template_index]->residue(ires).is_polymer()
							&& !templates_[initial_template_index]->residue(ires).is_lower_terminus()
							&& !pose.residue(pose.total_residue()).is_upper_terminus() ) {
						pose.append_residue_by_bond(templates_[initial_template_index]->residue(ires));
					} else {
						pose.append_residue_by_jump(templates_[initial_template_index]->residue(ires), 1);
					}
					hetatms.push_back( std::make_pair( ires, pose.total_residue() ) );
				}
			}
		}

		// (2) realign structures per-domain
		if ( realign_domains_ ) {
			//
			if ( std::find( non_null_template_indices_.begin(), non_null_template_indices_.end(), initial_template_index )
					!= non_null_template_indices_.end() ) {
				align_templates_by_domain(templates_[initial_template_index], domains_all_templ_[initial_template_index]);
			} else {
				if ( non_null_template_indices_.size() == 0 ) {
					TR << "No non-extended templates.  Skipping alignment." << std::endl;
				} else {
					TR << "Cannot align to extended template! Choosing a random template instead." << std::endl;
					core::Size aln_target = numeric::random::random_range(1, non_null_template_indices_.size() );
					align_templates_by_domain(templates_[aln_target], domains_all_templ_[aln_target]);
				}
			}
		}

		// (3) apply symmetry
		std::string symmdef_file = symmdef_files_[initial_template_index];
		if ( !symmdef_file.empty() && symmdef_file != "NULL" ) {
			protocols::simple_moves::symmetry::SetupForSymmetryMover makeSymm( symmdef_file );
			makeSymm.apply(pose);

			//fpd   to get the right rotamer set we need to do this
			basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].value( "dummy" );

			// xyz copy hetatms (properly handle cases where scoring subunit is not the first)
			if ( add_hetatm_ ) {
				for ( Size ihet=1; ihet <= hetatms.size(); ++ihet ) {
					core::conformation::Residue const &res_in = templates_[initial_template_index]->residue(hetatms[ihet].first);
					for ( Size iatm=1; iatm<=res_in.natoms(); ++iatm ) {
						core::id::AtomID tgt(iatm,hetatms[ihet].second);
						pose.set_xyz( tgt, res_in.xyz( iatm ) );
					}
				}
			}
		}

		// (4) add a virtual if we have a map or coord csts
		//     keep after symm
		//     should we check if a map is loaded (or density scoring is on) instead??
		if ( option[ OptionKeys::edensity::mapfile ].user() || user_csts_.size() > 0 ) {
			MoverOP dens( new protocols::electron_density::SetupForDensityScoringMover );
			dens->apply( pose );
		}

		// (5) initialize template history
		//     keep after symm
		TemplateHistoryOP history( new TemplateHistory(pose) );
		history->setall( initial_template_index );
		pose.data().set( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY, history );

		// (6) allowed movement
		utility::vector1<bool> allowed_to_move;
		allowed_to_move.resize(pose.total_residue(),true);
		for ( int i=1; i<=(int)residue_sample_template_.size(); ++i ) allowed_to_move[i] = allowed_to_move[i] && residue_sample_template_[i];
		for ( int i=1; i<=(int)residue_sample_abinitio_.size(); ++i ) allowed_to_move[i] = allowed_to_move[i] && residue_sample_abinitio_[i];

		// (7) steal crystal parameters (if specified)
		if ( templates_[initial_template_index]->pdb_info() ) {
			core::io::CrystInfo ci = templates_[initial_template_index]->pdb_info()->crystinfo();
			if ( ci.A()*ci.B()*ci.C() != 0 ) {
				pose.pdb_info()->set_crystinfo( ci );
			}
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
		if ( numeric::random::rg().uniform() < stage1_probability_ ) {

			for ( core::Size repeatstage1=0; repeatstage1 < jump_move_repeat_; ++repeatstage1 ) {

				std::string cst_fn = template_cst_fn_[initial_template_index];

				FoldTreeHybridizeOP ft_hybridize(
					new FoldTreeHybridize(
					initial_template_index, templates_aln_, template_weights_,
					template_chunks_, frags_small, frags_big) ) ;

				ft_hybridize->set_constraint_file( cst_fn );
				ft_hybridize->set_scorefunction( stage1_scorefxn_ );

				// options
				ft_hybridize->set_increase_cycles( stage1_increase_cycles_ );
				ft_hybridize->set_stage1_1_cycles( stage1_1_cycles_ );
				ft_hybridize->set_stage1_2_cycles( stage1_2_cycles_ );
				ft_hybridize->set_stage1_3_cycles( stage1_3_cycles_ );
				ft_hybridize->set_stage1_4_cycles( stage1_4_cycles_ );
				ft_hybridize->set_add_non_init_chunks( add_non_init_chunks_ );
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

				// allowed movement
				ft_hybridize->set_per_residue_controls( residue_sample_template_, residue_sample_abinitio_ );
				ft_hybridize->set_minimize_at_end( min_after_stage1_ );
				ft_hybridize->set_minimize_sf( stage2_scorefxn_ );

				// max insertion
				ft_hybridize->set_max_insertion( max_contig_insertion_ );

				// other cst stuff
				ft_hybridize->set_user_csts( user_csts_ );

				// finally run stage 1
				ft_hybridize->apply(pose);

				// get strand pairings if they exist for constraints in later stages
				strand_pairs_ = ft_hybridize->get_strand_pairs();

				//jump perturbation and minimization
				//fpd!  these docking moves should be incorporated as part of stage 1!
				if ( jump_move_ ) {
					do_intrastage_docking( pose );
				}
			} //end of repeatstage1

		} else {
			// just do frag insertion in unaligned regions
			core::pose::PoseOP chosen_templ = templates_aln_[initial_template_index];
			protocols::loops::Loops chosen_contigs = template_contigs_[initial_template_index];
			initialize_and_sample_loops(pose, chosen_templ, chosen_contigs, stage1_scorefxn_);
		}

		// realign domains to the output of stage 1
		if ( jump_move_ || realign_domains_stage2_ ) {
			TR << "Realigning template domains to stage1 pose." << std::endl;
			core::pose::PoseOP stage1pose( new core::pose::Pose( pose ) );
			align_templates_by_domain(stage1pose); // don't use stored domains; recompute parse from model
		}

		//write gdtmm to output
		if ( native_ && native_->total_residue() ) {
			gdtmm = get_gdtmm(*native_, pose, aln_);
			core::pose::setPoseExtraScore( pose, "GDTMM_after_stage1", gdtmm);
			TR << "GDTMM_after_stage1" << F(8,3,gdtmm) << std::endl;
		}

		// STAGE 2
		// apply constraints
		if ( stage2_scorefxn_->get_weight( core::scoring::atom_pair_constraint ) != 0 ) {
			std::string cst_fn = template_cst_fn_[initial_template_index];
			if ( !keep_pose_constraint_ ) {
				setup_centroid_constraints( pose, templates_aln_, template_weights_, cst_fn );
			}
			if ( add_hetatm_ ) {
				add_non_protein_cst(pose, *templates_aln_[initial_template_index], hetatm_self_cst_weight_, hetatm_prot_cst_weight_);
			}
			if ( strand_pairs_.size() ) {
				add_strand_pairs_cst(pose, strand_pairs_);
			}
		}

		// add torsion constraints derived from fragments
		if ( csts_from_frags_ ) {
			if ( stage2_scorefxn_->get_weight( core::scoring::dihedral_constraint ) != 0 ) {
				add_fragment_csts( pose );
			} else {
				TR << "Warning! csts_from_frags is on but dihedral_constraint weight=0.  Ignoring!" << std::endl;
			}
		}

		if ( stage2_scorefxn_->get_weight( core::scoring::coordinate_constraint ) != 0 ) {
			if ( user_csts_.size() > 0 ) {
				setup_user_coordinate_constraints(pose,user_csts_);
			}
		}

		if ( !option[cm::hybridize::skip_stage2]() ) {
			core::scoring::ScoreFunctionOP stage2_scorefxn_clone = stage2_scorefxn_->clone();

			core::scoring::methods::EnergyMethodOptions lowres_options(stage2_scorefxn_clone->energy_method_options());
			lowres_options.set_cartesian_bonded_linear(true);
			stage2_scorefxn_clone->set_energy_method_options(lowres_options);

			CartesianHybridizeOP cart_hybridize(
				new CartesianHybridize( templates_aln_, template_weights_, template_chunks_,template_contigs_, frags_big ) );

			// default scorefunctions (cenrot-compatable)
			if ( stage2_scorefxn_!=NULL ) cart_hybridize->set_scorefunction( stage2_scorefxn_ );
			if ( stage2pack_scorefxn_!=NULL ) cart_hybridize->set_pack_scorefunction( stage2pack_scorefxn_ );
			if ( stage2min_scorefxn_!=NULL ) cart_hybridize->set_min_scorefunction( stage2min_scorefxn_ );

			// set options
			cart_hybridize->set_increase_cycles( stage2_increase_cycles_ );
			cart_hybridize->set_no_global_frame( no_global_frame_ );
			cart_hybridize->set_linmin_only( linmin_only_ );
			cart_hybridize->set_seqfrags_only( seqfrags_only_ );
			cart_hybridize->set_cartfrag_overlap( cartfrag_overlap_ );
			cart_hybridize->set_skip_long_min( skip_long_min_ );
			cart_hybridize->set_cenrot( cenrot_ );
			cart_hybridize->set_fragment_probs(fragprob_stage2_, randfragprob_stage2_);

			// max insertion
			cart_hybridize->set_max_insertion( max_contig_insertion_ );

			// per-residue controls
			cart_hybridize->set_per_residue_controls( residue_sample_template_, residue_sample_abinitio_ );

			// finally run stage 2
			cart_hybridize->apply(pose);
		}

		//write gdtmm to output
		if ( native_ && native_->total_residue() ) {
			gdtmm = get_gdtmm(*native_, pose, aln_);
			core::pose::setPoseExtraScore( pose, "GDTMM_after_stage2", gdtmm);
			TR << "GDTMM_after_stage2" << ObjexxFCL::format::F(8,3,gdtmm) << std::endl;
		}
		// get fragment history
		runtime_assert( pose.data().has( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY ) );
		history = utility::pointer::static_pointer_cast< TemplateHistory >( pose.data().get_ptr( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY ) );

		TR << "History :";
		for ( Size i=1; i<= history->size(); ++i ) { TR << I(4,i); }
		TR << std::endl;
		TR << "History :";
		for ( Size i=1; i<= history->size(); ++i ) { TR << I(4, history->get(i)); }
		TR << std::endl;

		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb  ( true );
		mm->set_chi ( true );
		mm->set_jump( true );

		//fpd -- in positions where no fragment insertions are allowed, also allow no minimization
		if ( residue_sample_template_.size() > 0 && residue_sample_abinitio_.size() > 0 ) {
			for ( int i=1; i<=(int)nres_tgt; ++i ) {
				if ( !residue_sample_template_[i] && !residue_sample_abinitio_[i] ) {
					TR.Trace << "locking residue " << i << std::endl;
					mm->set_bb  ( i, false );
					mm->set_chi ( i, false );
				}
			}
		}

		if ( core::pose::symmetry::is_symmetric(pose) ) {
			core::pose::symmetry::make_symmetric_movemap( pose, *mm );
		}

		// stage "2.5" .. minimize with centroid energy + full-strength cart bonded
		if ( !option[cm::hybridize::skip_stage2]() ) {
			core::optimization::MinimizerOptions options_lbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
			core::optimization::CartesianMinimizer minimizer;
			Size n_min_cycles =(Size) (200.*stage25_increase_cycles_);
			options_lbfgs.max_iter(n_min_cycles);
			(*stage2_scorefxn_)(pose); minimizer.run( pose, *mm, *stage2_scorefxn_, options_lbfgs );
		}

		// STAGE 3: RELAX
		if ( batch_relax_ > 0 ) {
			// set disulfides before going to FA
			if ( disulf_file_.length() > 0 ) {
				// delete current disulfides
				for ( int i=1; i<=(int)pose.total_residue(); ++i ) {
					if ( !pose.residue(i).is_protein() ) continue;
					if ( pose.residue(i).type().is_disulfide_bonded() ) {
						core::conformation::change_cys_state( i, "CYS", pose.conformation() );
					}
				}

				// manually add new ones
				TR << " add disulfide: " << std::endl;
				basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].value(disulf_file_);
				core::pose::initialize_disulfide_bonds(pose);
				// must reset this since initialize_disulfide_bonds is used in pose io
				basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].deactivate();
				basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].to_default(); // reset to the default value
			} else {
				pose.conformation().detect_disulfides();
			}

			protocols::moves::MoverOP tofa( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD ) );
			tofa->apply(pose);

			// apply fa constraints
			std::string cst_fn = template_cst_fn_[initial_template_index];
			if ( fa_scorefxn_->get_weight( core::scoring::atom_pair_constraint ) != 0 ) {
				if ( !keep_pose_constraint_ ) {
					setup_fullatom_constraints( pose, templates_aln_, template_weights_, cst_fn, fa_cst_fn_ );
				} else {
					pose.constraint_set(save_pose_constraint_set);
				}
				if ( add_hetatm_ ) {
					add_non_protein_cst(pose, *templates_aln_[initial_template_index], hetatm_self_cst_weight_, hetatm_prot_cst_weight_);
				}
				if ( strand_pairs_.size() ) {
					add_strand_pairs_cst(pose, strand_pairs_);
				}
			}

			// add torsion constraints derived from fragments
			if ( csts_from_frags_ ) {
				if ( fa_scorefxn_->get_weight( core::scoring::dihedral_constraint ) != 0 ) {
					add_fragment_csts( pose );
				} else {
					TR << "Warning! csts_from_frags is on but dihedral_constraint weight=0.  Ignoring!" << std::endl;
				}
			}

			if ( stage2_scorefxn_->get_weight( core::scoring::coordinate_constraint ) != 0 ) {
				// note that CSTs are updated based on stage 1 movement
				//   this may or may not be desirable
				if ( user_csts_.size() > 0 ) {
					setup_user_coordinate_constraints(pose,user_csts_);
				}
			}

			if ( batch_relax_ == 1 ) {
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

				if ( post_centroid_structs.size() == batch_relax_ ) {
					protocols::relax::FastRelax relax_prot( fa_scorefxn_ );
					relax_prot.min_type("lbfgs_armijo_nonmonotone");
					relax_prot.set_force_nonideal(true);
					relax_prot.set_script_to_batchrelax_default( relax_repeats_ );

					// need to use a packer task factory to handle poses with different disulfide patterning
					core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory );
					tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
					tf->push_back( TaskOperationCOP( new core::pack::task::operation::IncludeCurrent ) );
					tf->push_back( TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
					relax_prot.set_task_factory( tf );

					// notice! this assumes all poses in a set have the same constraints!
					relax_prot.batch_apply(post_centroid_structs, pose.constraint_set()->clone());

					// reinflate pose
					post_centroid_structs[0]->fill_pose( pose );
				} else {
					protocols::moves::MoverOP tocen( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );
					tocen->apply(pose);
					need_more_samples = true;
				}
			}
		} else {
			// no fullatom sampling
			(*stage2_scorefxn_)(pose);
		}
	}
	if ( native_ && native_->total_residue() ) {
		gdtmm = get_gdtmm(*native_, pose, aln_);
		core::pose::setPoseExtraScore( pose, "GDTMM_final", gdtmm);
		TR << "GDTMM_final" << ObjexxFCL::format::F(8,3,gdtmm) << std::endl;
	}
}


void
HybridizeProtocol::align_templates_by_domain(core::pose::PoseOP & ref_pose, utility::vector1 <Loops> domains)
{
	// get domain parse if not specified
	core::pose::PoseOP working_pose;
	if ( core::pose::symmetry::is_symmetric(*ref_pose) ) {
		working_pose = core::pose::PoseOP(new core::pose::Pose);
		core::pose::symmetry::extract_asymmetric_unit(*ref_pose, *working_pose);
	} else {
		working_pose = ref_pose;
	}

	if ( domains.size() == 0 ) {
		DDomainParse ddom(pcut_,hcut_,length_);
		utility::vector1 <Loops> domains = ddom.split( *working_pose );
	}

	// clone original poses -> copies for alignment
	templates_aln_.clear();
	for ( Size i_pose=1; i_pose <= templates_.size(); ++i_pose ) {
		templates_aln_.push_back(
			core::pose::PoseOP( new core::pose::Pose( *templates_[i_pose] ) ) );
	}

	for ( Size i_pose=1; i_pose <= templates_aln_.size(); ++i_pose ) {
		if ( templates_aln_[i_pose] == ref_pose ) continue; // compare pointers
		align_by_domain(*templates_aln_[i_pose], *working_pose, domains);
		//std::string out_fn = template_fn_[i_pose] + "_realigned.pdb";
		//poses[i_pose]->dump_pdb(out_fn);
	}
}


void
HybridizeProtocol::align_by_domain(core::pose::Pose & pose, core::pose::Pose const & ref_pose, utility::vector1 <Loops> domains) {
	// ASSUME:
	//    domains cover the full pose (no gaps)
	utility::vector1< core::Real > aln_cutoffs;
	aln_cutoffs.push_back(6);
	aln_cutoffs.push_back(4);
	aln_cutoffs.push_back(3);
	aln_cutoffs.push_back(2);
	aln_cutoffs.push_back(1.5);
	aln_cutoffs.push_back(1);
	core::Real min_coverage = 0.2;

	// find corresepondance of ligands to nearest residues (in moving pose)
	core::Size last_protein_residue=pose.total_residue();
	while ( last_protein_residue>0 && !pose.residue_type(last_protein_residue).is_protein() ) last_protein_residue--;

	if ( last_protein_residue == 0 ) return;

	std::map< core::Size, core::Size > ligres_map;
	for ( core::Size i=last_protein_residue+1; i<=pose.total_residue(); ++i ) {
		core::Real mindist = 999.0;
		core::Size minj = 1;
		for ( core::Size j=1; j<=last_protein_residue; ++j ) {
			core::Real dist_ij = (pose.residue(i).xyz(1) - pose.residue(j).xyz(1)).length(); // hacky (using atom 1 only)
			if ( dist_ij < mindist ) {
				mindist = dist_ij;
				minj = (pose.pdb_info()) ? pose.pdb_info()->number(j) : j;
			}
		}
		ligres_map[i] = minj;
	}


	for ( Size i_domain = 1; i_domain <= domains.size() ; ++i_domain ) {
		core::id::AtomID_Map< core::id::AtomID > atom_map;
		core::pose::initialize_atomid_map( atom_map, pose, core::id::BOGUS_ATOM_ID );

		// find all residues in moving pose corresponding to this domain
		std::list <Size> residue_list;
		core::Size n_mapped_residues=0;
		for ( core::Size ires=1; ires<=pose.total_residue(); ++ires ) {
			if ( !pose.residue_type(ires).is_protein() ) continue;
			int pose_res = (pose.pdb_info()) ? pose.pdb_info()->number(ires) : ires;
			for ( core::Size iloop=1; iloop<=domains[i_domain].num_loop(); ++iloop ) {
				if ( pose_res < (int)domains[i_domain][iloop].start() || pose_res > (int)domains[i_domain][iloop].stop() ) continue;
				residue_list.push_back(ires);
			}
		}

		// find all residues in reference pose corresponding to this domain
		std::list <Size> ref_residue_list;
		for ( core::Size jres=1; jres<=ref_pose.total_residue(); ++jres ) {
			if ( !ref_pose.residue_type(jres).is_protein() ) continue;
			int ref_pose_res = (ref_pose.pdb_info()) ? ref_pose.pdb_info()->number(jres) : jres;
			for ( core::Size iloop=1; iloop<=domains[i_domain].num_loop(); ++iloop ) {
				if ( ref_pose_res < (int)domains[i_domain][iloop].start() || ref_pose_res > (int)domains[i_domain][iloop].stop() ) continue;
				ref_residue_list.push_back(jres);
			}
		}

		// get tmalign sequence mapping
		TMalign tm_align;
		std::string seq_pose, seq_ref, aligned;
		int reval = tm_align.apply(pose, ref_pose, residue_list, ref_residue_list);
		if ( reval != 0 ) continue;  // TO DO: remove this domain

		tm_align.alignment2AtomMap(pose, ref_pose, residue_list, ref_residue_list, n_mapped_residues, atom_map);
		tm_align.alignment2strings(seq_pose, seq_ref, aligned);

		using namespace ObjexxFCL::format;
		Size norm_length = residue_list.size() < ref_residue_list.size() ? residue_list.size():ref_residue_list.size();
		TR << "Align domain with TMscore of " << F(8,3,tm_align.TMscore(norm_length)) << std::endl;
		TR << seq_pose << std::endl;
		TR << aligned << std::endl;
		TR << seq_ref << std::endl;

		if ( n_mapped_residues < 6 ) continue;  // TO DO: remove this domain

		// add in ligand residues
		for ( core::Size i=last_protein_residue+1; i<=pose.total_residue(); ++i ) {
			core::Size res_controlling_i = ligres_map[i];
			for ( core::Size iloop=1; iloop<=domains[i_domain].num_loop(); ++iloop ) {
				if ( res_controlling_i < domains[i_domain][iloop].start() || res_controlling_i > domains[i_domain][iloop].stop() ) continue;
				residue_list.push_back(i);
			}
		}

		partial_align(pose, ref_pose, atom_map, residue_list, true, aln_cutoffs, min_coverage); // iterate_convergence = true
	}
}


void
HybridizeProtocol::do_intrastage_docking(core::pose::Pose & pose) {
	// call docking or symm docking
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		/////
		/////  SYMM LOGIC
		protocols::symmetric_docking::SymDockingLowResOP docking_lowres_mover( new protocols::symmetric_docking::SymDockingLowRes(stage1_scorefxn_) );
		docking_lowres_mover->apply(pose);

		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( true );
		core::pose::symmetry::make_symmetric_movemap( pose, *mm );

		protocols::simple_moves::symmetry::SymMinMoverOP min_mover(
			new protocols::simple_moves::symmetry::SymMinMover( mm, stage1_scorefxn_, "lbfgs_armijo_nonmonotone", 0.01, true ) );
		min_mover->apply(pose);
	} else {
		/////
		/////  ASYMM LOGIC
		core::Size const rb_move_jump = 1; // use the first jump as the one between partners <<<< fpd: MAKE THIS A PARSIBLE OPTION
		protocols::docking::DockingLowResOP docking_lowres_mover( new protocols::docking::DockingLowRes( stage1_scorefxn_, rb_move_jump ) );
		docking_lowres_mover->apply(pose);

		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb( false ); mm->set_chi( false ); mm->set_jump( rb_move_jump, true );
		protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( mm, stage1_scorefxn_, "lbfgs_armijo_nonmonotone", 0.01, true ) );
		min_mover->apply(pose);
	}
}

protocols::moves::MoverOP HybridizeProtocol::clone() const { return protocols::moves::MoverOP( new HybridizeProtocol( *this ) ); }
protocols::moves::MoverOP HybridizeProtocol::fresh_instance() const { return protocols::moves::MoverOP( new HybridizeProtocol ); }

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
	// basic options
	stage1_increase_cycles_ = tag->getOption< core::Real >( "stage1_increase_cycles", 1. );
	stage2_increase_cycles_ = tag->getOption< core::Real >( "stage2_increase_cycles", 1. );
	stage2_temperature_ = tag->getOption< core::Real >( "stage2_temperature", 2. );
	stage25_increase_cycles_ = tag->getOption< core::Real >( "stage2.5_increase_cycles", 1. );
	fa_cst_fn_ = tag->getOption< std::string >( "fa_cst_file", "" );
	batch_relax_ = tag->getOption< core::Size >( "batch" , 1 );
	jump_move_= tag->getOption< bool >( "jump_move" , false );
	jump_move_repeat_= tag->getOption< core::Size >( "jump_move_repeat" , 1 );
	keep_pose_constraint_= tag->getOption< bool >( "keep_pose_constraint" , false );
	cenrot_= tag->getOption< bool >( "cenrot" , false );

	//if( tag->hasOption( "task_operations" ) ){
	//FPD   remove this option, it was used as a residue selector and not as options controlling packing
	//      the 'DetailedControls' tag should replace one part of this, and a separate mover should generate interface constraints

	// force starting template
	//FPD   removed this option; it is redundant with weights! (set weight to 0 for templates we do not start with)


	if ( tag->hasOption( "csts_from_frags" ) ) {
		csts_from_frags_ = tag->getOption< bool >( "csts_from_frags" );
	}
	if ( tag->hasOption( "max_contig_insertion" ) ) {
		max_contig_insertion_ = tag->getOption< int >( "max_contig_insertion" );
	}

	if ( tag->hasOption( "min_after_stage1" ) ) {
		min_after_stage1_ = tag->getOption< bool >( "min_after_stage1" );
	}
	if ( tag->hasOption( "fragprob_stage2" ) ) {
		fragprob_stage2_ = tag->getOption< core::Real >( "fragprob_stage2" );
	}
	if ( tag->hasOption( "randfragprob_stage2" ) ) {
		randfragprob_stage2_ = tag->getOption< core::Real >( "randfragprob_stage2" );
	}


	// tons of ab initio options
	if ( tag->hasOption( "stage1_1_cycles" ) ) {
		stage1_1_cycles_ = tag->getOption< core::Size >( "stage1_1_cycles" );
	}
	if ( tag->hasOption( "stage1_2_cycles" ) ) {
		stage1_2_cycles_ = tag->getOption< core::Size >( "stage1_2_cycles" );
	}
	if ( tag->hasOption( "stage1_3_cycles" ) ) {
		stage1_3_cycles_ = tag->getOption< core::Size >( "stage1_3_cycles" );
	}
	if ( tag->hasOption( "stage1_4_cycles" ) ) {
		stage1_4_cycles_ = tag->getOption< core::Size >( "stage1_4_cycles" );
	}
	if ( tag->hasOption( "stage1_probability" ) ) {
		stage1_probability_ = tag->getOption< core::Real >( "stage1_probability" );
	}
	if ( tag->hasOption( "add_hetatm" ) ) {
		add_hetatm_ = tag->getOption< bool >( "add_hetatm" );
	}
	if ( tag->hasOption( "hetatm_cst_weight" ) ) {
		hetatm_self_cst_weight_ = tag->getOption< core::Real >( "hetatm_cst_weight" );
	}
	if ( tag->hasOption( "hetatm_to_protein_cst_weight" ) ) {
		hetatm_prot_cst_weight_ = tag->getOption< core::Real >( "hetatm_to_protein_cst_weight" );
	}
	if ( tag->hasOption( "realign_domains" ) ) {
		realign_domains_ = tag->getOption< bool >( "realign_domains" );
	}
	if ( tag->hasOption( "realign_domains_stage2" ) ) {
		realign_domains_stage2_ = tag->getOption< bool >( "realign_domains_stage2" );
	}
	if ( tag->hasOption( "add_non_init_chunks" ) ) {
		add_non_init_chunks_ = tag->getOption< core::Real >( "add_non_init_chunks" );
	}
	if ( tag->hasOption( "frag_1mer_insertion_weight" ) ) {
		frag_1mer_insertion_weight_ = tag->getOption< core::Real >( "frag_1mer_insertion_weight" );
	}
	if ( tag->hasOption( "small_frag_insertion_weight" ) ) {
		small_frag_insertion_weight_ = tag->getOption< core::Real >( "small_frag_insertion_weight" );
	}
	if ( tag->hasOption( "big_frag_insertion_weight" ) ) {
		big_frag_insertion_weight_ = tag->getOption< core::Real >( "big_frag_insertion_weight" );
	}
	if ( tag->hasOption( "chunk_insertion_weight" ) ) {
		chunk_insertion_weight_ = tag->getOption< core::Real >( "chunk_insertion_weight" );
	}
	if ( tag->hasOption( "frag_weight_aligned" ) ) {
		frag_weight_aligned_ = tag->getOption< core::Real >( "frag_weight_aligned" );
	}
	if ( tag->hasOption( "auto_frag_insertion_weight" ) ) {
		auto_frag_insertion_weight_ = tag->getOption< bool >( "auto_frag_insertion_weight" );
	}
	if ( tag->hasOption( "max_registry_shift" ) ) {
		max_registry_shift_ = tag->getOption< core::Size >( "max_registry_shift" );
	}
	if ( tag->hasOption( "repeats" ) ) {
		relax_repeats_ = tag->getOption< core::Size >( "repeats" );
	}
	if ( tag->hasOption( "disulf_file" ) ) {
		disulf_file_ = tag->getOption< std::string >( "disulf_file" );
	}

	// stage 2-specific options
	if ( tag->hasOption( "no_global_frame" ) ) {
		no_global_frame_ = tag->getOption< bool >( "no_global_frame" );
	}
	if ( tag->hasOption( "linmin_only" ) ) {
		linmin_only_ = tag->getOption< bool >( "linmin_only" );
	}
	if ( tag->hasOption( "cartfrag_overlap" ) ) {
		cartfrag_overlap_ = tag->getOption< core::Size >( "cartfrag_overlap" );
	}
	if ( tag->hasOption( "seqfrags_only" ) ) {
		seqfrags_only_ = tag->getOption< core::Size >( "seqfrags_only" );
	}
	if ( tag->hasOption( "skip_long_min" ) ) {
		skip_long_min_ = tag->getOption< core::Size >( "skip_long_min" );
	}


	// scorefxns
	if ( tag->hasOption( "stage1_scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "stage1_scorefxn" ) );
		stage1_scorefxn_ = (data.get< ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	}
	if ( tag->hasOption( "stage2_scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "stage2_scorefxn" ) );
		stage2_scorefxn_ = (data.get< ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	}
	if ( tag->hasOption( "stage2_min_scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "stage2_min_scorefxn" ) );
		stage2min_scorefxn_ = (data.get< ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	}
	if ( tag->hasOption( "stage2_pack_scorefxn" ) ) {
		if ( !cenrot_ ) {
			TR << "Warning! Ignoring stage2_pack_scorefxn declaration since cenrot is not set." << std::endl;
		} else {
			std::string const scorefxn_name( tag->getOption<std::string>( "stage2_pack_scorefxn" ) );
			stage2pack_scorefxn_ = (data.get< ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
		}
	}
	if ( tag->hasOption( "fa_scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "fa_scorefxn" ) );
		fa_scorefxn_ = (data.get< ScoreFunction * >( "scorefxns", scorefxn_name ))->clone();
	}

	// ddomain options
	hcut_ = tag->getOption< core::Real >( "domain_hcut" , 0.18);
	pcut_ = tag->getOption< core::Real >( "domain_pcut" , 0.81);
	length_ = tag->getOption< core::Size >( "domain_length" , 38);

	// user constraints
	if ( tag->hasOption( "coord_cst_res" ) ) {
		user_csts_ = core::pose::get_resnum_list_ordered( tag->getOption<std::string>( "coord_cst_res" ), pose );
	}

	// if user constraints are defined, make sure coord_csts are defined in at least one stage
	if ( user_csts_.size() > 0 ) {
		runtime_assert(
			stage1_scorefxn_->get_weight( core::scoring::coordinate_constraint ) > 0 ||
			stage2_scorefxn_->get_weight( core::scoring::coordinate_constraint ) > 0 ||
			fa_scorefxn_->get_weight( core::scoring::coordinate_constraint ) > 0 );
	}

	// fragments
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for ( tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() == "Fragments" ) {
			using namespace core::fragment;
			if ( (*tag_it)->hasOption( "3mers" ) ) {
				core::fragment::FragSetOP frags = core::fragment::FragmentIO().read_data( (*tag_it)->getOption<std::string>( "3mers" )  );
				fragments_small_.push_back(frags);
			} else if ( (*tag_it)->hasOption( "small" ) ) {
				utility::vector1<std::string> frag_files = utility::string_split((*tag_it)->getOption<std::string>( "small" ), ',');
				for ( core::Size i=1; i<= frag_files.size(); ++i ) {
					fragments_small_.push_back(core::fragment::FragmentIO().read_data( frag_files[i] ));
				}
			}
			if ( (*tag_it)->hasOption( "9mers" ) ) {
				core::fragment::FragSetOP frags = core::fragment::FragmentIO().read_data( (*tag_it)->getOption<std::string>( "9mers" ) );
				fragments_big_.push_back(frags);
			} else if ( (*tag_it)->hasOption( "big" ) ) {
				utility::vector1<std::string> frag_files = utility::string_split((*tag_it)->getOption<std::string>( "big" ), ',');
				for ( core::Size i=1; i<= frag_files.size(); ++i ) {
					fragments_big_.push_back(core::fragment::FragmentIO().read_data( frag_files[i] ));
				}
			}
		}

		std::string fasta = pose.sequence();
		if ( (*tag_it)->getName() == "Template" ) {
			std::string template_fn = (*tag_it)->getOption<std::string>( "pdb" );
			std::string cst_fn = (*tag_it)->getOption<std::string>( "cst_file", "AUTO" );  // should this be NONE?
			core::Real weight = (*tag_it)->getOption<core::Real>( "weight", 1 );
			std::string symm_file = (*tag_it)->getOption<std::string>( "symmdef", "" );

			// randomize some chains
			std::string rand_chain_str = (*tag_it)->getOption<std::string>( "randomize", "" );
			utility::vector1< std::string > rand_chain_strs = utility::string_split( rand_chain_str, ',');
			utility::vector1< char > rand_chains;
			for ( core::Size j = 1; j<=rand_chain_strs.size(); ++j ) {
				rand_chains.push_back( rand_chain_strs[j][0] );
			}

			bool align_pdb_info = (*tag_it)->getOption<bool>( "auto_align", true );

			add_template(template_fn, cst_fn, symm_file, weight, rand_chains);

			// validate here since we have the pose (add_template does not) ... could do this in apply as well (?)
			validate_template( template_fn, fasta, templates_[templates_.size()], align_pdb_info );

			if ( align_pdb_info ) {
				align_pdb_info = false;
				validate_template( template_fn, fasta, templates_[templates_.size()], align_pdb_info );
				update_last_template();
			}
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

		// per-residue control
		if ( (*tag_it)->getName() == "DetailedControls" ) {
			core::Size nres_nonvirt = get_num_residues_nonvirt(pose);

			residue_sample_template_.resize(nres_nonvirt, true);
			residue_sample_abinitio_.resize(nres_nonvirt, true);

			if ( (*tag_it)->hasOption( "task_operations" ) ) {
				core::pack::task::TaskFactoryOP task_factory = protocols::rosetta_scripts::parse_task_operations( *tag_it, data );
				core::pack::task::PackerTaskOP task = task_factory->create_task_and_apply_taskoperations( pose );

				for ( core::Size ires = 1; ires <= nres_nonvirt; ++ires ) {
					if ( !task->residue_task( ires ).being_packed() ) {
						residue_sample_template_[ires] = false;
						residue_sample_abinitio_[ires] = false;
					}
				}
			} else {
				core::Size start_res = (*tag_it)->getOption<core::Size>( "start_res", 1 );
				core::Size stop_res = (*tag_it)->getOption<core::Size>( "stop_res", nres_nonvirt );
				if ( (*tag_it)->hasOption( "sample_template" ) ) {
					bool sample_template = (*tag_it)->getOption<bool>( "sample_template", true );
					for ( core::Size ires=start_res; ires<=stop_res; ++ires ) {
						residue_sample_template_[ires] = sample_template;
					}
				}
				if ( (*tag_it)->hasOption( "sample_abinitio" ) ) {
					bool sample_abinitio = (*tag_it)->getOption<bool>( "sample_abinitio", true );
					for ( core::Size ires=start_res; ires<=stop_res; ++ires ) {
						residue_sample_abinitio_[ires] = sample_abinitio;
					}
				}
			} // if tag == DetailedControls
		} //forach tag
	}
}

} // hybridization
} // protocols
