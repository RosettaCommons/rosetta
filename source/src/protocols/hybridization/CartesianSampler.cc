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

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/comparative_modeling/coord_util.hh>
#include <protocols/moves/DataMap.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>

#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/rms_util.hh>
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
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	overlap_ = 2;
	ncycles_ = 250;
	input_as_ref_ = false;

	fragment_bias_strategy_ = "uniform";

	// default scorefunction
	set_scorefunction ( core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" ) );
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

void
CartesianSampler::apply_frame( core::pose::Pose & pose, core::fragment::Frame &frame ) {
	core::Size start = frame.start(),len = frame.length();
	core::Size end = start + len - 1;

	int aln_len = overlap_;
	runtime_assert( overlap_>=1 && overlap_<=len/2);

	// number of protein residues in the asymmetric unit
	core::Size nres = pose.total_residue();
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres = symm_info->num_independent_residues();
	}
	while (!pose.residue(nres).is_protein()) nres--;

	bool nterm = ( (start == 1) || pose.fold_tree( ).is_cutpoint(start-1) );
	bool cterm = ( (end == nres) || pose.fold_tree( ).is_cutpoint(end) );

	// insert frag
	core::pose::Pose pose_copy = pose;

	ObjexxFCL::FArray1D< numeric::Real > ww( 2*4*aln_len, 1.0 );
	ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
	numeric::xyzVector< core::Real > com1(0,0,0), com2(0,0,0);

	for (int i=0; i<(int)len; ++i) {
		core::conformation::idealize_position(start+i, pose_copy.conformation());
	}

	int maxtries = frame.nr_frags() / 2;
	for (int tries = 0; tries<maxtries; ++tries) {
		ww = 1.0;
		uu = 0.0;
		com1 = numeric::xyzVector< core::Real >(0,0,0);
		com2 = numeric::xyzVector< core::Real >(0,0,0);

		// grab coords
		ObjexxFCL::FArray2D< core::Real > init_coords( 3, 2*4*aln_len );
		for (int ii=-aln_len; ii<aln_len; ++ii) {
			int i = (ii>=0) ? (nterm?len-ii-1:ii) : (cterm?-ii-1:len+ii);
			numeric::xyzVector< core::Real > x_1 = pose.residue(start+i).atom(" C  ").xyz();
			numeric::xyzVector< core::Real > x_2 = pose.residue(start+i).atom(" O  ").xyz();
			numeric::xyzVector< core::Real > x_3 = pose.residue(start+i).atom(" CA ").xyz();
			numeric::xyzVector< core::Real > x_4 = pose.residue(start+i).atom(" N  ").xyz();
			com1 += x_1+x_2+x_3+x_4;
			for (int j=0; j<3; ++j) {
				init_coords(j+1,4*(ii+aln_len)+1) = x_1[j];
				init_coords(j+1,4*(ii+aln_len)+2) = x_2[j];
				init_coords(j+1,4*(ii+aln_len)+3) = x_3[j];
				init_coords(j+1,4*(ii+aln_len)+4) = x_4[j];
			}
		}
		com1 /= 2.0*4.0*aln_len;
		for (int ii=0; ii<2*4*aln_len; ++ii) {
			for ( int j=0; j<3; ++j ) init_coords(j+1,ii+1) -= com1[j];
		}

		core::Size toget = numeric::random::random_range( 1, frame.nr_frags() );
		frame.apply( toget, pose_copy );

		// grab new coords
		ObjexxFCL::FArray2D< core::Real > final_coords( 3, 2*4*aln_len );
		for (int ii=-aln_len; ii<aln_len; ++ii) {
			int i = (ii>=0) ? (nterm?len-ii-1:ii) : (cterm?-ii-1:len+ii);
			numeric::xyzVector< core::Real > x_1 = pose_copy.residue(start+i).atom(" C  ").xyz();
			numeric::xyzVector< core::Real > x_2 = pose_copy.residue(start+i).atom(" O  ").xyz();
			numeric::xyzVector< core::Real > x_3 = pose_copy.residue(start+i).atom(" CA ").xyz();
			numeric::xyzVector< core::Real > x_4 = pose_copy.residue(start+i).atom(" N  ").xyz();
			com2 += x_1+x_2+x_3+x_4;
			for (int j=0; j<3; ++j) {
			final_coords(j+1,4*(ii+aln_len)+1) = x_1[j];
				final_coords(j+1,4*(ii+aln_len)+2) = x_2[j];
				final_coords(j+1,4*(ii+aln_len)+3) = x_3[j];
				final_coords(j+1,4*(ii+aln_len)+4) = x_4[j];
			}
		}
		com2 /= 2.0*4.0*aln_len;
		for (int ii=0; ii<2*4*aln_len; ++ii) {
			for ( int j=0; j<3; ++j ) final_coords(j+1,ii+1) -= com2[j];
		}

		numeric::Real ctx; float rms;
		numeric::model_quality::findUU( final_coords, init_coords, ww, 2*4*aln_len, uu, ctx );
		numeric::model_quality::calc_rms_fast( rms, final_coords, init_coords, ww, 2*4*aln_len, ctx );

		if (rms < 0.5) break;
		if (tries >= maxtries/4 && rms < 1) break;
		if (tries >= maxtries/2 && rms < 1.5) break;
	}

	numeric::xyzMatrix< core::Real > R;
	R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
	R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
	R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );

	// apply rotation to ALL atoms
	// x_i' <- = R*x_i + com1;
	for ( Size i = 0; i < len; ++i ) {
		for ( Size j = 1; j <= pose.residue_type(start+i).natoms(); ++j ) {
			core::id::AtomID id( j, start+i );
			pose.set_xyz( id, R * ( pose_copy.xyz(id) - com2) + com1 );
		}
	}
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
	core::Real COORDDEV = 1.0;

	pose.remove_constraints();

	core::Size nres_tgt = get_num_residues_nonvirt( pose );

	utility::vector1< utility::vector1< core::Real > > tgt_dists(nres_tgt);
	utility::vector1< utility::vector1< core::Real > > tgt_weights(nres_tgt);

	for (core::Size j=1; j < ref_model_.total_residue(); ++j ) {
		if (!ref_model_.residue_type(j).is_protein()) continue;

		for (core::Size k=j+1; k < ref_model_.total_residue(); ++k ) {
			if (!ref_model_.residue_type(k).is_protein()) continue;
			if (ref_model_.pdb_info()->number(k) - ref_model_.pdb_info()->number(j) < (int)MINSEQSEP) continue;

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
							new USOGFunc( dist, COORDDEV )
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

	if (fragment_bias_strategy_ == "density") {
		// find segments with worst agreement to the density
		core::scoring::electron_density::ElectronDensity &edm = core::scoring::electron_density::getDensityMap();
		edm.setScoreWindowContext( true );
		edm.setWindow( 3 );

		// score the pose
		core::scoring::ScoreFunctionOP myscore = new core::scoring::ScoreFunction();
		myscore->set_weight( core::scoring::elec_dens_window, 1.0 );
		if (pose.is_fullatom())
			myscore->set_weight( core::scoring::fa_rep, 1.0 );
		else
			myscore->set_weight( core::scoring::vdw, 1.0 );

		if (core::pose::symmetry::is_symmetric(pose) ) {
			myscore = new core::scoring::symmetry::SymmetricScoreFunction(*myscore);
		}

		utility::vector1<core::Real> per_resCC;
		per_resCC.resize(nres);
		fragmentProbs.resize(nres);
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
			fragmentProbs[r] = exp( (CCsum-per_resCC[r])/CCsum2 );
			TR << "Prob_dens_density( " << r << " ) = " << fragmentProbs[r] << " ; CC=" << per_resCC[r] << " (Z=" << (CCsum-per_resCC[r])/CCsum2 << ")" << std::endl;
		}
	} else if (fragment_bias_strategy_ == "bfactors") {
		// find segments with highest bfactors
		fragmentProbs.resize(nres);
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
	} else if (fragment_bias_strategy_ == "user") {
		// user defined segments to rebuild
		runtime_assert( user_pos_.size()>0 );

		for (int r=1; r<=(int)nres; ++r) {
			fragmentProbs[r] = 0.0;
			if ( user_pos_.find(r) != user_pos_.end() ) fragmentProbs[r] = 1.0;
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
		fragments_.push_back( create_fragment_set(pose, 9, 25) );
		update_fragment_library_pointers();
	}

	// using the current fragment_bias_strategy_, compute bias
	// do this from fullatom (if the input pose is fullatom)
	compute_fragment_bias(pose);

	//
	bool fullatom_input = pose.is_fullatom();
	protocols::moves::MoverOP restore_sc;
	if (fullatom_input) {
	  restore_sc = new protocols::simple_moves::ReturnSidechainMover( pose );
		protocols::moves::MoverOP tocen = new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );
		tocen->apply( pose );
	}

	// constraints to reference model
	// save current constraint set
	core::scoring::constraints::ConstraintSetOP saved_csts = pose.constraint_set()->clone();
	if (input_as_ref_) ref_model_ = pose;
	apply_constraints( pose );

	// minimizer
	core::optimization::MinimizerOptions options_minilbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
	options_minilbfgs.max_iter(5);
	core::optimization::MinimizerOptions options_lbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
	options_lbfgs.max_iter(200);

	// to do ... make this parsable
	core::optimization::CartesianMinimizer minimizer;
	core::kinematics::MoveMap mm;
	mm.set_bb  ( true );
	mm.set_chi ( true );
	mm.set_jump( true );

	if (core::pose::symmetry::is_symmetric(pose) ) {
		core::pose::symmetry::make_symmetric_movemap( pose, mm );
	}

	Pose pose_in = pose;
	core::Size nres = pose.total_residue();
	if (pose.residue(nres).aa() == core::chemical::aa_vrt) nres--;
	core::Size n_prot_res = pose.total_residue();
	while (!pose.residue(n_prot_res).is_protein()) n_prot_res--;

	(*scorefxn_)(pose);
	protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( pose, *scorefxn_, 2.0 );

	for (int n=1; n<=(int)ncycles_; ++n) {
		// pick fragment set
		core::Size i_frag_set = numeric::random::random_range(1, fragments_.size());

		// pick insertion position
		core::Size insert_pos = frag_bias_[i_frag_set].random_sample(RG);

		if (library_[i_frag_set].find(insert_pos) != library_[i_frag_set].end())
			apply_frame (pose, library_[i_frag_set][insert_pos]);

		// MC
		(*scorefxn_)(pose);
		minimizer.run( pose, mm, *scorefxn_, options_minilbfgs );
		mc->boltzmann( pose );
	}
	mc->recover_low(pose);

	// final minimization
	(*scorefxn_)(pose); minimizer.run( pose, mm, *scorefxn_, options_lbfgs );
	(*scorefxn_)(pose);

	if (fullatom_input) {
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
	utility::tag::TagPtr const tag, moves::DataMap & data, filters::Filters_map const & , moves::Movers_map const & , core::pose::Pose const & pose )
{
	using namespace core::scoring;

	// scorefunction
	if( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn" ) );
		set_scorefunction ( data.get< ScoreFunction * >( "scorefxns", scorefxn_name ) );
	}

	// options
	if( tag->hasOption( "overlap" ) ) {
		overlap_ = tag->getOption<core::Size>( "overlap" );
	}
	if( tag->hasOption( "ncycles" ) ) {
		ncycles_ = tag->getOption<core::Size>( "ncycles" );
	}
	if( tag->hasOption( "strategy" ) ) {
		fragment_bias_strategy_ = tag->getOption<std::string>( "strategy" );
	}

	if( tag->hasOption( "reference_model" ) ) {
		std::string ref_model_pdb = tag->getOption<std::string>( "reference_model" );
		if (ref_model_pdb != "none")
			core::import_pose::pose_from_pdb( ref_model_, ref_model_pdb );
	} else {
		input_as_ref_ = true;
	}

	if( tag->hasOption( "residues" ) ) {
		user_pos_ = core::pose::get_resnum_list(	tag->getOption<std::string>( "residues" ), pose );
	}

	// fragments
	utility::vector1< utility::tag::TagPtr > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagPtr >::const_iterator tag_it;
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
			fragments_.push_back( create_fragment_set(pose, atoi(fraglens[i].c_str()), nfrags) );
		}
	}
	update_fragment_library_pointers();
}

}
//}
}
