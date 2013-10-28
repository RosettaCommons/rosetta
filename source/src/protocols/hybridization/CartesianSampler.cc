// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Yifan Song
/// @author Frank DiMaio

#include <protocols/hybridization/HybridizeProtocolCreator.hh>
#include <protocols/hybridization/CartesianSampler.hh>
#include <protocols/hybridization/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/relax/util.hh>

#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/symmetry/SetupNCSMover.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/comparative_modeling/coord_util.hh>
#include <basic/datacache/DataMap.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>

#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/cryst/PhenixInterface.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/Energies.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Remarks.hh>
#include <core/pose/selection.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>

#include <core/pose/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ScalarWeightedFunc.hh>
#include <core/scoring/constraints/SOGFunc.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/USOGFunc.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/conformation/util.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>


#include <utility/excn/Exceptions.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/WeightedSampler.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>

#include <boost/unordered/unordered_map.hpp>

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

static basic::Tracer TR("protocols.hybridization.CartesianSampler");
static numeric::random::RandomGenerator RG(8403155);

/////////////
// creator
std::string
CartesianSamplerCreator::keyname() const {
	return CartesianSamplerCreator::mover_name();
}

protocols::moves::MoverOP
CartesianSamplerCreator::create_mover() const {
	return new CartesianSampler;
}

std::string
CartesianSamplerCreator::mover_name() {
	return "CartesianSampler";
}

////////////

CartesianSampler::CartesianSampler( ) {
	init();
}

CartesianSampler::CartesianSampler(
		utility::vector1 < core::fragment::FragSetOP > fragments_in ) {
	init();

	fragments_ = fragments_in;
	update_fragment_library_pointers();
}

void
CartesianSampler::init() {
	temp_ = 2.0;
	rms_cutoff_ = 1.5;
	overlap_ = 2;
	ncycles_ = 250;
	nminsteps_ = 10;
	ref_cst_weight_ = 1.0;
	input_as_ref_ = false;
	fullatom_ = false;
	bbmove_ = false;

	fragment_bias_strategy_ = "uniform";
	selection_bias_ = "none";

	// default scorefunction
	if (!fullatom_)
		set_scorefunction ( core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" ) );
	else
		set_scorefunction ( core::scoring::ScoreFunctionFactory::create_score_function( "soft_rep" ) );

	set_fa_scorefunction ( core::scoring::getScoreFunction() );
}

void
CartesianSampler::update_fragment_library_pointers() {
	core::Size nfragsets = fragments_.size();

	// map positions to fragments
	library_.resize( nfragsets );
	for (int i=1; i<=(int)nfragsets; ++i) {
		for (core::fragment::ConstFrameIterator j = fragments_[i]->begin(); j != fragments_[i]->end(); ++j) {
			core::Size position = (*j)->start();
			library_[i][position] = **j;
		}
	}
}

protocols::moves::MoverOP CartesianSampler::clone() const { return new CartesianSampler( *this ); }
protocols::moves::MoverOP CartesianSampler::fresh_instance() const { return new CartesianSampler; }

// get frag->pose transform, return RMS
core::Real
CartesianSampler::get_transform(
					core::pose::Pose const &pose, core::pose::Pose const &frag, core::Size startpos,
					core::Vector &preT, core::Vector &postT, numeric::xyzMatrix< core::Real > &R) {
	int aln_len = overlap_;
	int len = frag.total_residue();
	if (frag.residue(len).aa() == core::chemical::aa_vrt) len--;
	ObjexxFCL::FArray1D< core::Real > ww( 2*4*aln_len, 1.0 );
	ObjexxFCL::FArray2D< core::Real > uu( 3, 3, 0.0 ), init_coords( 3, 2*4*aln_len ), final_coords( 3, 2*4*aln_len );
	preT = numeric::xyzVector< core::Real >(0,0,0);
	postT = numeric::xyzVector< core::Real >(0,0,0);

	bool nterm = ( (startpos == 1) || pose.fold_tree().is_cutpoint(startpos-1) );
	bool cterm = ( pose.fold_tree( ).is_cutpoint(startpos+len-1) );

	// grab coords
	for (int ii=-aln_len; ii<aln_len; ++ii) {
		int i = (ii>=0) ? (nterm?len-ii-1:ii) : (cterm?-ii-1:len+ii);
		numeric::xyzVector< core::Real > x_1 = pose.residue(startpos+i).atom(1).xyz();
		numeric::xyzVector< core::Real > x_2 = pose.residue(startpos+i).atom(2).xyz();
		numeric::xyzVector< core::Real > x_3 = pose.residue(startpos+i).atom(3).xyz();
		numeric::xyzVector< core::Real > x_4 = pose.residue(startpos+i).atom(4).xyz();
		postT += x_1+x_2+x_3+x_4;
		for (int j=0; j<3; ++j) {
			init_coords(j+1,4*(ii+aln_len)+1) = x_1[j];
			init_coords(j+1,4*(ii+aln_len)+2) = x_2[j];
			init_coords(j+1,4*(ii+aln_len)+3) = x_3[j];
			init_coords(j+1,4*(ii+aln_len)+4) = x_4[j];
		}
	}
	postT /= 2.0*4.0*aln_len;
	for (int ii=0; ii<2*4*aln_len; ++ii) {
		for ( int j=0; j<3; ++j ) init_coords(j+1,ii+1) -= postT[j];
	}

	// grab new coords
	for (int ii=-aln_len; ii<aln_len; ++ii) {
		int i = (ii>=0) ? (nterm?len-ii-1:ii) : (cterm?-ii-1:len+ii);
		numeric::xyzVector< core::Real > x_1 = frag.residue(i+1).atom(1).xyz();
		numeric::xyzVector< core::Real > x_2 = frag.residue(i+1).atom(2).xyz();
		numeric::xyzVector< core::Real > x_3 = frag.residue(i+1).atom(3).xyz();
		numeric::xyzVector< core::Real > x_4 = frag.residue(i+1).atom(4).xyz();
		preT += x_1+x_2+x_3+x_4;
		for (int j=0; j<3; ++j) {
			final_coords(j+1,4*(ii+aln_len)+1) = x_1[j];
			final_coords(j+1,4*(ii+aln_len)+2) = x_2[j];
			final_coords(j+1,4*(ii+aln_len)+3) = x_3[j];
			final_coords(j+1,4*(ii+aln_len)+4) = x_4[j];
		}
	}
	preT /= 2.0*4.0*aln_len;
	for (int ii=0; ii<2*4*aln_len; ++ii) {
		for ( int j=0; j<3; ++j ) final_coords(j+1,ii+1) -= preT[j];
	}

	numeric::Real ctx;
	float rms;
	numeric::model_quality::findUU( final_coords, init_coords, ww, 2*4*aln_len, uu, ctx );
	numeric::model_quality::calc_rms_fast( rms, final_coords, init_coords, ww, 2*4*aln_len, ctx );

	R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
	R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
	R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );

	return ((core::Real)rms);
}

// apply csts to frag
void
CartesianSampler::apply_fragcsts( core::pose::Pose &working_frag,	core::pose::Pose const &pose, core::Size startpos ) {
	using namespace core::scoring::constraints;

	working_frag.remove_constraints();
	int len = working_frag.total_residue();
	if (working_frag.residue(len).aa() == core::chemical::aa_vrt) len--;

	bool nterm = ( (startpos == 1) || pose.fold_tree().is_cutpoint(startpos-1) );
	bool cterm = ( pose.fold_tree( ).is_cutpoint(startpos+len-1) );

	if (!nterm) {
		for (int j=0; j<overlap_; ++j) {
			for (int i=1; i<=3; ++i) {
				working_frag.add_constraint(
					new CoordinateConstraint(
						core::id::AtomID(i,j+1),
						core::id::AtomID(2,working_frag.total_residue()),
						pose.residue(startpos+j).atom(i).xyz(),
						new HarmonicFunc( 0.0, 1.0 )
					)
				);
			}
		}
	}
	if (!cterm) {
		for (int j=len-overlap_; j<len; ++j) {
			for (int i=1; i<=3; ++i) {
				working_frag.add_constraint(
					new CoordinateConstraint(
						core::id::AtomID(i,j+1),
						core::id::AtomID(2,working_frag.total_residue()),
						pose.residue(startpos+j).atom(i).xyz(),
						new HarmonicFunc( 0.0, 1.0 )
					)
				);
			}
		}
	}
}


// transform fragment
void
CartesianSampler::apply_transform( core::pose::Pose &frag, core::Vector const &preT, core::Vector const &postT, numeric::xyzMatrix< core::Real > const &R) {
	// apply rotation to ALL atoms
	// x_i' <- = R*x_i + com1;
	for ( Size i = 1; i <= frag.total_residue(); ++i ) {
		for ( Size j = 1; j <= frag.residue_type(i).natoms(); ++j ) {
			core::id::AtomID id( j, i );
			frag.set_xyz( id, R * ( frag.xyz(id) - preT) + postT );
		}
	}
}



bool
CartesianSampler::apply_frame( core::pose::Pose & pose, core::fragment::Frame &frame ) {
	core::Size start = frame.start(),len = frame.length();
	core::Size end = start + len - 1;
	runtime_assert( overlap_>=1 && overlap_<=len/2);

	// set the frame's insert point to 1
	frame.shift_to(1);

	// see if the pose has NCS
	simple_moves::symmetry::NCSResMappingOP ncs;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ) {
		ncs = ( static_cast< simple_moves::symmetry::NCSResMapping* >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING )() ));
	}

	// make subpose at this position
	// we assume the frame does not cross a jump
	core::pose::Pose frag;
	for (int i=0; i<(int)len; ++i) frag.append_residue_by_bond( pose.residue( start+i ) );
	for (int i=0; i<(int)len; ++i) core::conformation::idealize_position(i+1, frag.conformation());

	core::Vector preT, postT;
	numeric::xyzMatrix< core::Real > R;
	core::Size maxtries,frag_toget;
	core::Real rms;

	if (selection_bias_ == "none") {
		maxtries = frame.nr_frags(); // bias-dependent
		int tries;
		bool frag_chosen=false;
		for (tries = 0; tries<maxtries && !frag_chosen; ++tries) {
			frag_toget = numeric::random::random_range( 1, frame.nr_frags() ); // bias-dependent
			frame.apply( frag_toget, frag );
			rms = get_transform( pose,  frag,  start,	preT, postT, R);
			frag_chosen = (rms < 0.5) || (tries >= maxtries/4 && rms < 1) || (tries >= maxtries/2 && rms < 2); // bias-dependent
		}
		TR << "after " << tries+1 << " tries, fragment " << frag_toget << "chosen with RMS = " << rms << std::endl;

		apply_transform( frag,	preT, postT, R);
		for ( Size i = 0; i < len; ++i )
			pose.replace_residue( start+i, frag.residue(i+1), false );
	} else if (selection_bias_ == "density") {

		// set up minimizer
		core::optimization::AtomTreeMinimizer rbminimizer;
		core::optimization::MinimizerOptions options_rb( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
		options_rb.max_iter(20);
		core::kinematics::MoveMap mm_rb;
		mm_rb.set_bb ( bbmove_ );
		mm_rb.set_chi ( false );  //?
		mm_rb.set_jump ( true );

		// scorefunctions
		core::scoring::ScoreFunctionOP nonsymm_fa_scorefxn = new core::scoring::ScoreFunction(*fa_scorefxn_);
		if (nonsymm_fa_scorefxn->get_weight( core::scoring::elec_dens_fast ) == 0)
			nonsymm_fa_scorefxn->set_weight( core::scoring::elec_dens_fast , 20 );
		core::scoring::ScoreFunctionOP densonly = new core::scoring::ScoreFunction();
		if (bbmove_ && nonsymm_fa_scorefxn->get_weight( core::scoring::coordinate_constraint ) == 0)
			nonsymm_fa_scorefxn->set_weight( core::scoring::coordinate_constraint , 1 );

		densonly->set_weight( core::scoring::elec_dens_fast, 5.0 );

		// set up packer
		core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
		main_task_factory->push_back( new core::pack::task::operation::RestrictToRepacking );
		protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover;
		pack_mover->task_factory( main_task_factory );
		pack_mover->score_function( nonsymm_fa_scorefxn );

		// prepare fragment
		core::pose::addVirtualResAsRoot(frag);
		protocols::simple_moves::SwitchResidueTypeSetMover to_fa("fa_standard");
		protocols::simple_moves::SwitchResidueTypeSetMover to_cen("centroid");
		if (!fullatom_) to_fa.apply( frag );

		core::Size nattempts = 0;
		core::Real best_dens_score = 1e30, best_rms;
		for (int i=1; i<=frame.nr_frags(); ++i) {
			core::pose::Pose working_frag = frag;
			frame.apply( i, working_frag );
			rms = get_transform( pose,  working_frag,  start,	preT, postT, R);
			if (rms<=rms_cutoff_) {
				nattempts++;

				// orient
				apply_transform( working_frag,	preT, postT, R);

				// cenmin+repack+min
				to_cen.apply( working_frag );
				(*densonly)(working_frag);
				rbminimizer.run( working_frag, mm_rb, *densonly, options_rb );
				to_fa.apply( working_frag );
				pack_mover->apply( working_frag );

				core::Real dens_score;
				if (bbmove_) {
					apply_fragcsts( working_frag,	pose, start );
					(*nonsymm_fa_scorefxn)(working_frag);
					rbminimizer.run( working_frag, mm_rb, *nonsymm_fa_scorefxn, options_rb );
					dens_score = (*nonsymm_fa_scorefxn) (working_frag);
				} else {
					(*densonly)(working_frag);
					rbminimizer.run( working_frag, mm_rb, *densonly, options_rb );
					dens_score = (*densonly) (working_frag);
				}

				if (dens_score<best_dens_score) {
					best_rms = rms;
					best_dens_score = dens_score;
					frag = working_frag;
				}
				if (nattempts >=25) break;  //fpd within helices often all fragments will match
			}
		}

		if (nattempts > 0) {
			TR << "Best frag ( out of "<< nattempts << ") with score="  << best_dens_score << "  rms=" << best_rms << std::endl;

			// best scoring fragment saved
			if (!fullatom_) to_cen.apply( frag );
			for ( Size i = 0; i < len; ++i )
				pose.replace_residue( start+i, frag.residue(i+1), false );
		} else {
			TR << "No acceptable fragments" << std::endl;
			frame.shift_to(start);
			return false;
		}
	} else {
		utility_exit_with_message("Unrecognized fragbias!");
	}

	// apply to NCS-symmetric copies
	if (ncs) {
		for (int j=1; j<=ncs->ngroups(); ++j ) {
			bool all_are_mapped = true;
			for ( Size k= 0; k< len && all_are_mapped; ++k )
				all_are_mapped &= (ncs->get_equiv( j,start+k )!=0);
			if (!all_are_mapped) continue;

			core::Size remap_start = ncs->get_equiv( j, start );
			get_transform( pose,  frag,  remap_start,	preT, postT, R);
			apply_transform( frag,	preT, postT, R);
			for ( Size i = 0; i < len; ++i )
				pose.replace_residue( remap_start+i, frag.residue(i+1), false );
		}
	}

	// restore frame
	frame.shift_to(start);
	return true;
}


// apply constraints from ref_pose into current pose
void
CartesianSampler::apply_constraints( core::pose::Pose &pose )
{
	using namespace core::scoring::constraints;

	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	core::Size MINSEQSEP = 8;
	core::Real MAXDIST = 12.0;
	core::Size GAPBUFFER = 3;
	core::Real COORDDEV = 0.39894;

	bool user_rebuild = (fragment_bias_strategy_ == "user");

	core::Size nres_tgt = get_num_residues_nonvirt( pose );

	utility::vector1< utility::vector1< core::Real > > tgt_dists(nres_tgt);
	utility::vector1< utility::vector1< core::Real > > tgt_weights(nres_tgt);

	for (core::Size j=1; j < ref_model_.total_residue(); ++j ) {
		if (!ref_model_.residue_type(j).is_protein()) continue;

		if (user_rebuild && user_pos_.find(j) != user_pos_.end() ) continue;
		if (user_rebuild && loops_ && loops_->is_loop_residue(j) ) continue;

		for (core::Size k=j+1; k < ref_model_.total_residue(); ++k ) {
			if (!ref_model_.residue_type(k).is_protein()) continue;
			if (ref_model_.pdb_info()->number(k) - ref_model_.pdb_info()->number(j) < (int)MINSEQSEP) continue;

			if (user_rebuild && user_pos_.find(k) != user_pos_.end() ) continue;
			if (user_rebuild && loops_ && loops_->is_loop_residue(k) ) continue;

			core::Real dist = ref_model_.residue(j).xyz(2).distance( ref_model_.residue(k).xyz(2) );

			if ( dist <= MAXDIST ) {
				core::Size resid_j = ref_model_.pdb_info()->number(j);
				core::Size resid_k = ref_model_.pdb_info()->number(k);
				char chnid_j = ref_model_.pdb_info()->chain(j);
				char chnid_k = ref_model_.pdb_info()->chain(k);

				core::Size tgt_resid_j = pose.pdb_info()->pdb2pose(chnid_j, resid_j);
				core::Size tgt_resid_k = pose.pdb_info()->pdb2pose(chnid_k, resid_k);

				runtime_assert ( tgt_resid_j > 0 );
				runtime_assert ( tgt_resid_k > 0 );

				// ???
				if (symm_info && !symm_info->bb_is_independent( tgt_resid_j ) ) continue;
				if (symm_info && !symm_info->bb_is_independent( tgt_resid_k ) )	continue;

				pose.add_constraint(
						new AtomPairConstraint( core::id::AtomID(2,tgt_resid_j), core::id::AtomID(2,tgt_resid_k),
							new ScalarWeightedFunc( ref_cst_weight_, new USOGFunc( dist, COORDDEV ) )
						)
					);
			}
		}
	}
}


void
CartesianSampler::compute_fragment_bias(Pose & pose) {
	using namespace core::scoring;

	utility::vector1< core::Real > fragmentProbs;
	frag_bias_.resize( fragments_.size() );

	// get nres, accounting for symmetry/vrt/ligands
	core::Size nres = pose.total_residue();
	core::conformation::symmetry::SymmetryInfoCOP symminfo(0);
	if (core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symminfo = SymmConf.Symmetry_Info();
		nres = symminfo->num_independent_residues();
	}
	while (!pose.residue(nres).is_protein()) nres--;

	fragmentProbs.resize(nres);

	if (fragment_bias_strategy_ == "density") {
		// find segments with worst agreement to the density
		core::scoring::electron_density::ElectronDensity &edm = core::scoring::electron_density::getDensityMap();
		edm.setScoreWindowContext( true );
		edm.setWindow( 3 );  // smoother to use 3-res window

		// score the pose
		core::scoring::ScoreFunctionOP myscore = new core::scoring::ScoreFunction();
		myscore->set_weight( core::scoring::elec_dens_window, 1.0 );

		// make sure interaction graph gets computed
		if (pose.is_fullatom())
			myscore->set_weight( core::scoring::fa_rep, 1.0 );
		else
			myscore->set_weight( core::scoring::vdw, 1.0 );

		if (core::pose::symmetry::is_symmetric(pose) )
			myscore = new core::scoring::symmetry::SymmetricScoreFunction(*myscore);

		utility::vector1<core::Real> per_resCC;
		per_resCC.resize(nres);
		core::Real CCsum=0, CCsum2=0;
		(*myscore)(pose);

		for (int r=1; r<=(int)nres; ++r) {
			per_resCC[r] = core::scoring::electron_density::getDensityMap().matchRes( r , pose.residue(r), pose, symminfo , false);
			CCsum += per_resCC[r];
			CCsum2 += per_resCC[r]*per_resCC[r];
		}
		CCsum /= nres;
		CCsum2 = sqrt( CCsum2/nres-CCsum*CCsum );

		for (int r=1; r<=(int)nres; ++r) {
			if (per_resCC[r]<0.6)
				fragmentProbs[r] = 1.0;
			else if (per_resCC[r]<0.8)
				fragmentProbs[r] = 0.25;
			else
				fragmentProbs[r] = 0.01;
			TR << "residue " << r << ": " << " rscc=" << per_resCC[r] << " weight=" <<fragmentProbs[r] << std::endl;
		}
	} else if (fragment_bias_strategy_ == "bfactors") {
		// find segments with highest bfactors
		core::Real Btemp=25;  // no idea what value makes sense here
		                      // with Btemp = 25, a B=100 is ~54 times more likely to be sampled than B=0

		runtime_assert( pose.pdb_info() );
		for (int r=1; r<=(int)nres; ++r) {
			core::Real Bsum=0;
			core::Size nbb = pose.residue_type(r).last_backbone_atom();
			for (core::Size atm=1; atm<=nbb; ++atm) {
				Bsum += pose.pdb_info()->temperature( r, atm );
			}
			Bsum /= nbb;
			fragmentProbs[r] = exp( Bsum/Btemp );
			TR << "Prob_dens_bfact( " << r << " ) = " << fragmentProbs[r] << " ; B=" << Bsum << std::endl;
		}
	} else if (fragment_bias_strategy_ == "rama") {
		// find segments with worst rama score
		core::Real Rtemp=1;  // again, this is a guess

		// score the pose
		core::scoring::ScoreFunctionOP myscore = new core::scoring::ScoreFunction();
		myscore->set_weight( core::scoring::rama, 1.0 );
		if (core::pose::symmetry::is_symmetric(pose) ) {
			myscore = new core::scoring::symmetry::SymmetricScoreFunction(*myscore);
		}

		Energies & energies( pose.energies() );
		(*myscore)(pose);

		for (int r=1; r<=(int)nres; ++r) {
			EnergyMap & emap( energies.onebody_energies( r ) );
					// i dont think this will work for symmetric systems where 1st subunit is not the scoring one
			core::Real ramaScore = emap[ rama ];
			fragmentProbs[r] = exp( ramaScore / Rtemp );
			TR << "Prob_dens_rama( " << r << " ) = " << fragmentProbs[r] << " ; rama=" << ramaScore << std::endl;
		}
	} else if (fragment_bias_strategy_ == "chainbreak") {
		for (int r=1; r<(int)nres; ++r) {
			if (!pose.residue_type(r).is_protein()) continue;
			if (pose.fold_tree().is_cutpoint(r+1)) continue;

			numeric::xyzVector< core::Real > c0 , n1;
			c0 = pose.residue(r).atom(" C  ").xyz();
			n1 = pose.residue(r+1).atom(" N  ").xyz();
			core::Real d2 = c0.distance( n1 );
			if ( (d2-1.328685)*(d2-1.328685) > 0.1*0.1 ) {
				fragmentProbs[r] = 1.0;
			} else {
				fragmentProbs[r] = 0.001;
			}
			TR << "Prob_dens_cb( " << r << " ) = " << fragmentProbs[r] << std::endl;
		}
	} else if (fragment_bias_strategy_ == "user") {
		// user defined segments to rebuild
		runtime_assert( user_pos_.size()>0 || (loops_ && !loops_->empty()) );

		for (int r=1; r<=(int)nres; ++r) {
			fragmentProbs[r] = 0.0;
			if ( user_pos_.find(r) != user_pos_.end() ) fragmentProbs[r] = 1.0;
			if ( loops_ && loops_->is_loop_residue(r) ) fragmentProbs[r] = 1.0;
			TR << "Prob_dens_user( " << r << " ) = " << fragmentProbs[r] << std::endl;
		}
	} else {
		// default to uniform
		for (int r=1; r<=(int)nres; ++r) {
			fragmentProbs[r] = 1.0;
			TR << "Prob_dens_uniform( " << r << " ) = " << 1.0 << std::endl;
		}
	}


	// for each fragment size, smooth over the fragment window
	//   - handle mapping from frame->seqpos
	//   - don't allow any insertions that cross cuts
	for (Size i_frag_set = 1; i_frag_set<=fragments_.size(); ++i_frag_set) {
		utility::vector1< core::Real > frame_weights(fragments_[i_frag_set]->nr_frames(), 0.0);
		for (Size i_frame = 1; i_frame <= fragments_[i_frag_set]->nr_frames(); ++i_frame) {
			core::fragment::ConstFrameIterator frame_it = fragments_[i_frag_set]->begin(); // first frame of the fragment library
			advance(frame_it, i_frame-1);  // point frame_it to the i_frame of the library
			core::Size seqpos_start = (*frame_it)->start();  // find starting and ending residue seqpos of the inserted fragment
			core::Size seqpos_end   = (*frame_it)->end();

			frame_weights[i_frame] = 0;
			for(int i_pos = (int)seqpos_start; i_pos<=(int)seqpos_end; ++i_pos)
				frame_weights[i_frame] += fragmentProbs[i_pos];
			frame_weights[i_frame] /= (seqpos_end-seqpos_start+1);

			for(int i_pos = (int)seqpos_start; i_pos<=(int)seqpos_end-1; ++i_pos)
				if (pose.fold_tree().is_cutpoint(i_pos))
					frame_weights[i_frame] = 0.0;
		}
		frag_bias_[i_frag_set].weights(frame_weights);
	}
}

void
CartesianSampler::apply( Pose & pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose::datacache;

	// autogenerate fragments if they are not loaded yet
	if (fragments_.size() == 0) {
		fragments_.push_back( create_fragment_set_no_ssbias(pose, 9, 25) );
		update_fragment_library_pointers();
	}

	if (!mc_scorefxn_) mc_scorefxn_ = scorefxn_;

	// derived scorefunctions
	scorefxn_dens_ = new core::scoring::ScoreFunction();
	scorefxn_dens_->set_weight( core::scoring::elec_dens_fast, 1.0 );
	scorefxn_xray_ = new core::scoring::ScoreFunction();
	scorefxn_xray_->set_weight( core::scoring::xtal_ml, 100.0 );
	if (core::pose::symmetry::is_symmetric(pose) ) {
		scorefxn_dens_ = new core::scoring::symmetry::SymmetricScoreFunction(*scorefxn_dens_);
		scorefxn_xray_ = new core::scoring::symmetry::SymmetricScoreFunction(*scorefxn_xray_);
	}

	// see if the pose has NCS
	simple_moves::symmetry::NCSResMappingOP ncs;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ) {
		ncs = ( static_cast< simple_moves::symmetry::NCSResMapping* >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING )() ));
	}

	// using the current fragment_bias_strategy_, compute bias
	// do this from fullatom (if the input pose is fullatom)
	compute_fragment_bias(pose);

	// number of protein residues in ASU
	core::Size nres = pose.total_residue();
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
	 	core::conformation::symmetry::SymmetricConformation & SymmConf (
	 		dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
	 	symm_info = SymmConf.Symmetry_Info();
	 	nres = symm_info->num_independent_residues();
	}
	while (!pose.residue(nres).is_protein()) nres--;

	//
	bool fullatom_input = pose.is_fullatom();
	protocols::moves::MoverOP restore_sc;
	if (fullatom_input && !fullatom_) {
	  restore_sc = new protocols::simple_moves::ReturnSidechainMover( pose );
		protocols::moves::MoverOP tocen = new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );
		tocen->apply( pose );
	} else if (!fullatom_input && fullatom_) {
		utility_exit_with_message("ERROR! Expected fullatom input.");
	}

	// constraints to reference model
	// save current constraint set
	core::scoring::constraints::ConstraintSetOP saved_csts = pose.constraint_set()->clone();
	if (input_as_ref_) ref_model_ = pose;
	apply_constraints( pose );

	// stepwise minimizer
	core::optimization::MinimizerOptions options_minilbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
	options_minilbfgs.max_iter(nminsteps_);

	// to do ... make this parsable
	core::optimization::CartesianMinimizer minimizer;
	core::kinematics::MoveMap mm;
	mm.set_bb  ( true ); mm.set_chi ( true ); mm.set_jump( true );

	if (core::pose::symmetry::is_symmetric(pose) ) {
		core::pose::symmetry::make_symmetric_movemap( pose, mm );
	}

	Pose pose_in = pose;
	(*scorefxn_)(pose);
	protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( pose, *mc_scorefxn_, temp_ );

	if (TR.Debug.visible()) {
		scorefxn_->show_line_headers(TR);
		scorefxn_->show_line(TR,pose);
	}

	for (int n=1; n<=(int)ncycles_; ++n) {
		bool success=false;
		int try_count=1000;
		int insert_pos;
		core::Size i_frag_set;

		while (!success && --try_count>0) {
			// pick fragment set
			i_frag_set = numeric::random::random_range(1, fragments_.size());

			// pick insertion position
			insert_pos = (int)frag_bias_[i_frag_set].random_sample(RG);
			int ntries=50;
			while (library_[i_frag_set].find(insert_pos) == library_[i_frag_set].end() && --ntries>0)
				insert_pos = (int)frag_bias_[i_frag_set].random_sample(RG);

			if (library_[i_frag_set].find(insert_pos) == library_[i_frag_set].end()) {
				TR << "ERROR! unable to find allowed fragment inserts after " << ntries << " attempts.  Continuing.";
				continue;
			}

			success = apply_frame (pose, library_[i_frag_set][insert_pos]);
		}

		// restricted movemap
		// TO DO we should check of chainbreaks
		// TO DO min window extension (curr 6) should be a parameter
		core::kinematics::MoveMap mm_local;
		int start_move = std::max(1,insert_pos-6);
		int stop_move = std::min(nres,insert_pos+library_[i_frag_set][insert_pos].length()+5);
		for (int i=start_move; i<=stop_move; ++i) {
			mm_local.set_bb(i,true);
			mm_local.set_chi(i,true);
			// ncs
			if (ncs) {
				for (int j=1; j<=ncs->ngroups(); ++j ) {
					core::Size remap_i = ncs->get_equiv( j, i );
					if (remap_i!=0) {
						mm_local.set_bb(remap_i,true);
						mm_local.set_chi(remap_i,true);
					}
				}
			}
		}

		// min + MC
		core::Real scoreinit = (*scorefxn_)(pose);
		if (TR.Debug.visible()) {
			scorefxn_->show_line(TR,pose);
			TR << std::endl;
		}
		minimizer.run( pose, mm_local, *scorefxn_, options_minilbfgs );
		core::Real scorefinal = (*scorefxn_)(pose);
		if (TR.Debug.visible()) {
			scorefxn_->show_line(TR,pose);
			TR << std::endl;
		}

		(*mc_scorefxn_)(pose);
		bool accept = mc->boltzmann( pose );

		mc->show_scores();
		if (accept) {
			TR << "Insert at " << insert_pos << " accepted!" << std::endl;
		} else {
			TR << "Insert at " << insert_pos << " rejected!" << std::endl;
		}
	}
	mc->show_scores();
	mc->show_counters();
	mc->recover_low(pose);

	// final minimization (make optional?)
	//(*scorefxn_)(pose);
	//minimizer.run( pose, mm, *scorefxn_, options_lbfgs );
	//(*scorefxn_)(pose);

	if (fullatom_input && !fullatom_) {
		protocols::moves::MoverOP tofa = new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD );
		tofa->apply( pose );
		restore_sc->apply( pose );
	}

	// reapply old constraints
	pose.remove_constraints();
	pose.constraint_set( saved_csts );
}


// parse_my_tag
void
CartesianSampler::parse_my_tag(
	utility::tag::TagCOP const tag, basic::datacache::DataMap & data, filters::Filters_map const & , moves::Movers_map const & , core::pose::Pose const & pose )
{
	using namespace core::scoring;

	// scorefunction
	if( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn" ) );
		set_scorefunction ( data.get< ScoreFunction * >( "scorefxns", scorefxn_name ) );
	}

	// fullatom scorefunction
	if( tag->hasOption( "fascorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "fascorefxn" ) );
		set_fa_scorefunction ( data.get< ScoreFunction * >( "scorefxns", scorefxn_name ) );
	}

	// mc scorefunction
	if( tag->hasOption( "mcscorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "mcscorefxn" ) );
		mc_scorefxn_ = data.get< ScoreFunction * >( "scorefxns", scorefxn_name );
	}

	// options
	if( tag->hasOption( "fullatom" ) ) {
		fullatom_ = tag->getOption<bool >( "fullatom" );
	}
	if( tag->hasOption( "bbmove" ) ) {
		bbmove_ = tag->getOption<bool >( "bbmove" );
	}
	if( tag->hasOption( "temp" ) ) {
		temp_ = tag->getOption<core::Real  >( "temp" );
	}
	if( tag->hasOption( "rms" ) ) {
		rms_cutoff_ = tag->getOption<core::Real  >( "rms" );
	}
	if( tag->hasOption( "nminsteps" ) ) {
		nminsteps_ = tag->getOption<core::Size  >( "nminsteps" );
	}
	if( tag->hasOption( "overlap" ) ) {
		overlap_ = tag->getOption<core::Size>( "overlap" );
	}
	if( tag->hasOption( "ncycles" ) ) {
		ncycles_ = tag->getOption<core::Size>( "ncycles" );
	}
	if( tag->hasOption( "strategy" ) ) {
		fragment_bias_strategy_ = tag->getOption<std::string>( "strategy" );
	}
	if( tag->hasOption( "fragbias" ) ) {
		selection_bias_ = tag->getOption<std::string>( "fragbias" );
	}

	if( tag->hasOption( "reference_model" ) ) {
		std::string ref_model_pdb = tag->getOption<std::string>( "reference_model" );
		if (ref_model_pdb != "input")
			core::import_pose::pose_from_pdb( ref_model_, ref_model_pdb );
		else
			input_as_ref_ = true;
	}
	if( tag->hasOption( "reference_cst_wt" ) ) {
		ref_cst_weight_ = tag->getOption<core::Real  >( "reference_cst_wt" );
	}

	// user-specified residues
	runtime_assert( !( tag->hasOption( "residues" ) && tag->hasOption( "loops_in" ) ));  // one or the other
	if( tag->hasOption( "residues" ) ) {
		user_pos_ = core::pose::get_resnum_list(	tag->getOption<std::string>( "residues" ), pose );
	}
	if (tag->hasOption( "loops_in" ) ) {
		std::string looptag = tag->getOption<std::string>( "loops_in" );
		runtime_assert( data.has( "loops", looptag ) );
		loops_ = data.get< protocols::loops::Loops * >( "loops", looptag );
	}

	// fragments
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for (tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it) {
		if ( (*tag_it)->getName() == "Fragments" ) {
			using namespace core::fragment;
			fragments_.push_back( FragmentIO().read_data( (*tag_it)->getOption<std::string>( "fragfile" )  ) );
		}
	}

	core::Size nfrags = tag->getOption< core::Size >( "nfrags", 25 );

	// autofragments
	if( tag->hasOption( "fraglens" ) ) {
		utility::vector1<std::string> fraglens = utility::string_split( tag->getOption< std::string >( "fraglens" ), ',' );
		for (core::Size i=1; i<=fraglens.size(); ++i) {
			fragments_.push_back( create_fragment_set_no_ssbias(pose, atoi(fraglens[i].c_str()), nfrags) );
		}
	}
	update_fragment_library_pointers();
}

}
//}
}
