// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/seeded_abinitio/movers/
/// @brief repurposing logic and some functions from CartesianHybridze protocols for segment insertions and chain closure
/// @author Eva-Maria Strauch (evas01@uw.edu)

// Unit headers
#include <protocols/seeded_abinitio/SegmentHybridizer.hh>
#include <protocols/seeded_abinitio/SegmentHybridizerCreator.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <core/scoring/ScoreFunction.hh>
#include <boost/foreach.hpp>
#include <core/scoring/dssp/Dssp.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <numeric/xyzVector.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <numeric/random/random.hh>
#include <core/pose/selection.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/simple_moves/FragmentMover.hh>

#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/comparative_modeling/coord_util.hh>
#include <core/scoring/rms_util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/util.hh>

#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/conformation/util.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/picking_old/vall/util.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/util.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>

#include <basic/Tracer.hh>
#include <boost/unordered/unordered_map.hpp>

namespace protocols {
namespace seeded_abinitio {

static THREAD_LOCAL basic::Tracer TR( "protocols.seeded_abinitio.movers.SegmentHybridizer" );
static THREAD_LOCAL basic::Tracer TR_ccd( "protocols.seeded_abinitio.movers.SegmentHybridizer_ccd" );

std::string
SegmentHybridizerCreator::keyname() const
{
	return SegmentHybridizerCreator::mover_name();
}

protocols::moves::MoverOP
SegmentHybridizerCreator::create_mover() const {
	return protocols::moves::MoverOP( new SegmentHybridizer );
}

std::string
SegmentHybridizerCreator::mover_name()
{
	return "SegmentHybridizer";
}


SegmentHybridizer::SegmentHybridizer() :
	Mover( SegmentHybridizerCreator::mover_name() ),
	cartfrag_overlap_(2)
	/*big_( 9 ),
	small_( 3 ),
	nfrags_( 50 ),
	cartfrag_overlap_( 2 ),
	use_seq_( 0 ),
	rms_( 0.5 ),
	auto_mm_( 1 ),
	tries_( 100 ),
	mc_cycles_( 4 ),
	temp_(2.0),
	use_frags_(1),
	min_cycles_(2)
	*/
{
	//torsion_database_.clear();
	//delta_lengths_.clear();
}


SegmentHybridizer::~SegmentHybridizer() = default;

protocols::moves::MoverOP SegmentHybridizer::clone() const { return protocols::moves::MoverOP( new SegmentHybridizer( *this ) ); }
protocols::moves::MoverOP SegmentHybridizer::fresh_instance() const { return protocols::moves::MoverOP( new SegmentHybridizer ); }


void
SegmentHybridizer::init() {
	fragments_big_ = core::fragment::FragSetOP( new core::fragment::ConstantLengthFragSet( big_ ) );
	fragments_small_ = core::fragment::FragSetOP( new core::fragment::ConstantLengthFragSet( small_ ) );

	// default scorefunction
	set_scorefunction ( core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" ) );
}

void
SegmentHybridizer::set_scorefunction(core::scoring::ScoreFunctionOP scorefxn_in) {
	lowres_scorefxn_ = scorefxn_in->clone();

	min_scorefxn_ = scorefxn_in->clone();

	//bonds_scorefxn_ = new core::scoring::symmetry::SymmetricScoreFunction();
	bonds_scorefxn_ = scorefxn_in->clone();
	bonds_scorefxn_->reset();
	bonds_scorefxn_->set_weight( core::scoring::vdw, lowres_scorefxn_->get_weight( core::scoring::vdw ) );
	bonds_scorefxn_->set_weight( core::scoring::cart_bonded, lowres_scorefxn_->get_weight( core::scoring::cart_bonded ) );
	bonds_scorefxn_->set_weight( core::scoring::cart_bonded_angle, lowres_scorefxn_->get_weight( core::scoring::cart_bonded_angle ) );
	bonds_scorefxn_->set_weight( core::scoring::cart_bonded_length, lowres_scorefxn_->get_weight( core::scoring::cart_bonded_length ) );
	bonds_scorefxn_->set_weight( core::scoring::cart_bonded_torsion, lowres_scorefxn_->get_weight( core::scoring::cart_bonded_torsion ) );

	nocst_scorefxn_ = lowres_scorefxn_->clone();
	nocst_scorefxn_->set_weight( core::scoring::atom_pair_constraint, 0.0 );
}

void
SegmentHybridizer::apply_frame( core::pose::Pose & pose, core::fragment::Frame &frame ) {

	core::Size start = frame.start();
	core::Size len = frame.length();

	//std::cout << "len: " << len << std::endl;

	int aln_len = cartfrag_overlap_;
	runtime_assert( cartfrag_overlap_>=1 &&  cartfrag_overlap_<= len/2 + 1);//cartfrag_overlap_<=len/2);
	core::Size nres = pose.total_residue();

	if ( pose.residue(nres).aa() == core::chemical::aa_vrt ) nres--;
	bool nterm = (start == 1);
	bool cterm = (start == nres- len -1 );

	/// insert frag:

	core::pose::Pose pose_copy = pose;

	ObjexxFCL::FArray1D< numeric::Real > ww( 2*4*aln_len, 1.0 );
	ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
	numeric::xyzVector< core::Real > com1(0,0,0), com2(0,0,0);

	for ( int i=0; i<(int)len; ++i ) {
		core::conformation::idealize_position(start+i, pose_copy.conformation());
	}
	for ( int tries = 0; tries< tries_; ++tries ) {
		ww = 1.0;
		uu = 0.0;
		com1 = numeric::xyzVector< core::Real >(0,0,0);
		com2 = numeric::xyzVector< core::Real >(0,0,0);

		// grab coords
		ObjexxFCL::FArray2D< core::Real > init_coords( 3, 2*4*aln_len );
		for ( int ii=-aln_len; ii<aln_len; ++ii ) {
			int i = (ii>=0) ? (nterm?len-ii-1:ii) : (cterm?-ii-1:len+ii);
			numeric::xyzVector< core::Real > x_1 = pose.residue(start+i).atom(" C  ").xyz();
			numeric::xyzVector< core::Real > x_2 = pose.residue(start+i).atom(" O  ").xyz();
			numeric::xyzVector< core::Real > x_3 = pose.residue(start+i).atom(" CA ").xyz();
			numeric::xyzVector< core::Real > x_4 = pose.residue(start+i).atom(" N  ").xyz();
			com1 += x_1+x_2+x_3+x_4;
			for ( int j=0; j<3; ++j ) {
				init_coords(j+1,4*(ii+aln_len)+1) = x_1[j];
				init_coords(j+1,4*(ii+aln_len)+2) = x_2[j];
				init_coords(j+1,4*(ii+aln_len)+3) = x_3[j];
				init_coords(j+1,4*(ii+aln_len)+4) = x_4[j];
			}
		}
		com1 /= 2.0*4.0*aln_len;
		for ( int ii=0; ii<2*4*aln_len; ++ii ) {
			for ( int j=0; j<3; ++j ) init_coords(j+1,ii+1) -= com1[j];
		}

		core::Size toget = numeric::random::random_range( 1, frame.nr_frags() );
		frame.apply( toget, pose_copy );

		// grab new coords
		ObjexxFCL::FArray2D< core::Real > final_coords( 3, 2*4*aln_len );
		for ( int ii=-aln_len; ii<aln_len; ++ii ) {
			int i = (ii>=0) ? (nterm?len-ii-1:ii) : (cterm?-ii-1:len+ii);
			numeric::xyzVector< core::Real > x_1 = pose_copy.residue(start+i).atom(" C  ").xyz();
			numeric::xyzVector< core::Real > x_2 = pose_copy.residue(start+i).atom(" O  ").xyz();
			numeric::xyzVector< core::Real > x_3 = pose_copy.residue(start+i).atom(" CA ").xyz();
			numeric::xyzVector< core::Real > x_4 = pose_copy.residue(start+i).atom(" N  ").xyz();
			com2 += x_1+x_2+x_3+x_4;
			for ( int j=0; j<3; ++j ) {
				final_coords(j+1,4*(ii+aln_len)+1) = x_1[j];
				final_coords(j+1,4*(ii+aln_len)+2) = x_2[j];
				final_coords(j+1,4*(ii+aln_len)+3) = x_3[j];
				final_coords(j+1,4*(ii+aln_len)+4) = x_4[j];
			}
		}
		com2 /= 2.0*4.0*aln_len;
		for ( int ii=0; ii<2*4*aln_len; ++ii ) {
			for ( int j=0; j<3; ++j ) final_coords(j+1,ii+1) -= com2[j];
		}

		numeric::Real ctx;
		float rms;

		numeric::model_quality::findUU( final_coords, init_coords, ww, 2*4*aln_len, uu, ctx );
		numeric::model_quality::calc_rms_fast( rms, final_coords, init_coords, ww, 2*4*aln_len, ctx );

		TR.Debug  << "fragments rms:\t" << rms << std::endl;

		if ( rms < rms_ ) break;
		if ( tries >= 50 && rms < 1 ) break; //20
		if ( tries >= 70 && rms < 2 ) break; //40
		if ( tries >= 100 && rms < 4 ) break; //60
	}

	numeric::xyzMatrix< core::Real > R;
	R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
	R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
	R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );

	for ( Size i = 0; i < len; ++i ) {
		for ( Size j = 1; j <= pose.residue_type(start+i).natoms(); ++j ) {
			core::id::AtomID id( j, start+i );
			pose.set_xyz( id, R * ( pose_copy.xyz(id) - com2) + com1 );
		}
	}
}


void
SegmentHybridizer::check_and_create_fragments( core::pose::Pose & pose, core::Size insert_start, core::Size insert_stop ) {
	//if (fragments_big_ && fragments_small_) return;
	// how do I reset the fragments best?

	//if (!fragments_big_present) {
	/// want to use the native ss here!

	std::string tgt_seq = pose.sequence();
	core::scoring::dssp::Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );
	std::string tgt_ss = pose.secstruct();
	//std::cout << "fragment picking start position " << insert_start << " and stop: " << insert_stop <<std::endl;

	// pick from vall based on template SS or target sequence
	for ( core::Size j=insert_start; j<= insert_stop - big_-1 ; ++j ) {
		std::string ss_sub = tgt_ss.substr( j-1, big_ );
		std::string aa_sub = tgt_seq.substr( j-1, big_ );

		core::fragment::FrameOP frame( new core::fragment::Frame( j, big_ ) );
		if ( use_seq_ ) {
			frame->add_fragment( core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa( ss_sub, aa_sub,
				nfrags_, true, core::fragment::IndependentBBTorsionSRFD() ) );
			//std::cout << "picking " << nfrags_ << " fragments based on ss and sequence" <<  std::cout;
		} else {
			frame->add_fragment( core::fragment::picking_old::vall::pick_fragments_by_ss( ss_sub, nfrags_, true,
				core::fragment::IndependentBBTorsionSRFD() ) );
			//std::cout << "picking " << nfrags_ << " fragments based on ss only " <<  std::endl;
		}

		fragments_big_->add( frame );
	}
}

void
SegmentHybridizer::hybridize( core::pose::Pose & pose , core::Size insert_pos_start, core::Size insert_pos_stop){

	Pose pose_in = pose;
	core::Size nres = pose.total_residue();
	if ( pose.residue(nres).aa() == core::chemical::aa_vrt ) nres--;

	// pick an insert position based on gap
	utility::vector1<core::Real> residuals( nres , 0.0 );
	utility::vector1<core::Real> max_residuals(3,0);
	utility::vector1<int> max_poses(4,-1);


	for ( Size /*int*/ i=1; i<nres; ++i ) {
		if ( pose.fold_tree().is_cutpoint(i+1) ) {
			residuals[i] = -1;
		} else {
			numeric::xyzVector< core::Real > c0 , n1;
			c0 = pose.residue(i).atom(" C  ").xyz();
			n1 = pose.residue(i+1).atom(" N  ").xyz();
			core::Real d2 = c0.distance( n1 );
			residuals[i] = (d2-1.328685)*(d2-1.328685);
			if ( residuals[i] > max_residuals[1] ) {
				max_residuals[3] = max_residuals[2]; max_residuals[2] = max_residuals[1]; max_residuals[1] = residuals[i];
				max_poses[3] = max_poses[2]; max_poses[2] = max_poses[1]; max_poses[1] = i;
			} else if ( residuals[i] > max_residuals[2] ) {
				max_residuals[3] = max_residuals[2]; max_residuals[2] = residuals[i];
				max_poses[3] = max_poses[2]; max_poses[2] = i;
			} else if ( residuals[i] > max_residuals[3] ) {
				max_residuals[3] = residuals[i];
				max_poses[3] = i;
			}
		}
	}

	int select_position = numeric::random::random_range(1,3); //4);
	core::Size max_pos = max_poses[ select_position ];

	//std::cout << "selection position: "<< select_position <<
	//"\nmaxpos: " << max_pos << std::endl;

	// select random pos in the middle depending on fragment size
	core::Size insert_pos = max_pos - numeric::random::random_range(big_/2-1, big_/2); //+1);
	//std::cout << "insert position before apply frame : " << insert_pos << std::endl;

	//insert_pos = std::min( insert_pos, nres - big_-1);
	insert_pos = std::min( insert_pos, insert_pos_stop - big_-1);
	insert_pos = std::max( (int)insert_pos, (int)insert_pos_start);

	// for debugging of frames
	//for (boost::unordered_map<core::Size, core::fragment::Frame>::iterator iter = library_.begin() ; iter != library_.end() ; ++iter){
	//std::cout << " frame  " << (*iter).first << " len: " << library_[insert_pos].length() << std::endl;
	//}

	if ( library_.find(insert_pos) != library_.end() ) {
		apply_frame (pose, library_[insert_pos]);
		//TR<< "applying fragments on position: " << insert_pos << std::endl;
	}

}

void
SegmentHybridizer::apply( core::pose::Pose & pose ){
	//using namespace basic::options;
	//using namespace basic::options::OptionKeys;
	//using namespace core::pose::datacache;
	init();

	/// 1. setting weights -- needs some cleaning

	core::Real max_cart = lowres_scorefxn_->get_weight( core::scoring::cart_bonded );
	core::Real max_cart_angle = lowres_scorefxn_->get_weight( core::scoring::cart_bonded_angle );
	core::Real max_cart_length = lowres_scorefxn_->get_weight( core::scoring::cart_bonded_length );
	core::Real max_cart_torsion = lowres_scorefxn_->get_weight( core::scoring::cart_bonded_torsion );
	core::Real max_cst  = lowres_scorefxn_->get_weight( core::scoring::atom_pair_constraint );
	core::Real max_vdw  = lowres_scorefxn_->get_weight( core::scoring::vdw );

	core::Real bonded_weight = 0.1*max_cart;
	core::Real bonded_weight_angle = 0.1*max_cart_angle;
	core::Real bonded_weight_length = 0.1*max_cart_length;
	core::Real bonded_weight_torsion = 0.1*max_cart_torsion;
	core::Real cst_weight = 2*max_cst;
	core::Real vdw_weight = 0.1*max_vdw;

	lowres_scorefxn_->set_weight( core::scoring::cart_bonded, bonded_weight );
	lowres_scorefxn_->set_weight( core::scoring::cart_bonded_angle, bonded_weight_angle );
	lowres_scorefxn_->set_weight( core::scoring::cart_bonded_length, bonded_weight_length );
	lowres_scorefxn_->set_weight( core::scoring::cart_bonded_torsion, bonded_weight_torsion );
	lowres_scorefxn_->set_weight( core::scoring::atom_pair_constraint, cst_weight );
	lowres_scorefxn_->set_weight( core::scoring::vdw, vdw_weight );

	using namespace core::pose::datacache;

	protocols::moves::MoverOP tocen( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );
	tocen->apply( pose );

	/// 2. set up minimizer

	//core::optimization::MinimizerOptions options( "linmin", 0.01, true, false, false );
	core::optimization::MinimizerOptions options_minilbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
	options_minilbfgs.max_iter(5);
	core::optimization::MinimizerOptions options_lbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
	options_lbfgs.max_iter(200);
	core::optimization::CartesianMinimizer minimizer;


	//// 3. iterate through segments of interest

	for ( Size iter = 1 ; iter <= seg_vector_.size() ; ++ iter ) {

		///  runtime parsing of input information
		//core::Size insert_pos_start = protocols::rosetta_scripts::parse_resnum( seg_vector_[iter].first, pose ) - extend_outside_ ;
		//core::Size insert_pos_stop  = extend_outside_ + protocols::rosetta_scripts::parse_resnum( seg_vector_[iter].second, pose );
		core::Size insert_pos_start( core::pose::parse_resnum( seg_vector_[iter].first, pose ) - extend_outside_ );
		core::Size insert_pos_stop( extend_outside_ + core::pose::parse_resnum( seg_vector_[iter].second, pose ) );

		core::Size nterm_mm = insert_pos_start;
		core::Size cterm_mm = insert_pos_stop;
		//std::cout << "nterm before adjustment : "<< nterm_mm << " cterm: " << cterm_mm << std::endl;

		// assertionsssss ...

		/// 3.a. define and set movemap

		extended_mm_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
		extended_mm_->set_bb  ( false );
		extended_mm_->set_chi ( false );
		extended_mm_->set_jump( true );


		//if there is more than 1 chain, use only the last one
		if ( all_movable_ ) {
			TR<< "allowing all elements to move, if there is more than 1 chain, it will be the last chain that is allowed to move  " <<std::endl;
			core::Size chain_start  = 1 ;
			core::Size chain_num = pose.conformation().num_chains();
			if ( chain_num > 1 ) {
				chain_start = pose.conformation().chain_begin( chain_num );
			}
			for ( core::Size i = chain_start; i <= pose.total_residue(); i++ ) {
				extended_mm_->set_bb(i, true);
				extended_mm_->set_chi(i,true);
			}
		}

		if ( !all_movable_ ) {
			if ( auto_mm_ ) {
				core::scoring::dssp::Dssp dssp( pose );
				dssp.insert_ss_into_pose( pose );
				std::string tgt_ss = pose.secstruct();
				//extend the movemap on both sides until it hits a loop
				for ( core::Size it = insert_pos_start - 2 ; it > 0 ; it-- ) { //2 instead of 1 because it is typically misassigned as loop
					//std::cout << "ss: " << tgt_ss[it-1] << std::endl;
					if ( tgt_ss[it-1] == 'L' ) break;
					else {
						nterm_mm = it;
						//std::cout << "ss: " << tgt_ss[it-1] << std::endl;
					}
				}
				for ( Size it = insert_pos_stop + 1  ; it <= pose.total_residue() ; it++ ) { // 1 off because ss is typically misassigned
					//std::cout << "ss: " << tgt_ss[it-1] << std::endl;
					if ( tgt_ss[it-1] == 'L' ) break;
					else {
						cterm_mm = it;
						//std::cout << "ss: " << tgt_ss[it-1] << std::endl;
					}
				}
			}

			for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
				if ( i >= nterm_mm  && i <= cterm_mm ) {
					extended_mm_->set_bb( i, true );
					extended_mm_->set_chi( i, true );
					if ( i >= insert_pos_start + extend_outside_ + extend_inside_ && i <= insert_pos_stop - extend_outside_ ) {
						//reset to false
						extended_mm_->set_bb( i, false );
						extended_mm_->set_chi( i, false );
					}
				}
			}
		} //specific movemap

		/// 3.b fragment hybridize

		if ( use_frags_ ) {
			///  get fragments -- should sample more than 1
			// TODO check that the insert position does not go into the first chain!
			check_and_create_fragments( pose, insert_pos_start, insert_pos_stop );


			///  map resids to frames, keeping track of positions
			core::Size insert_frags_pos = fragments_big_->min_pos();
			//std::cout << "start frags = fragments->min_pos:"<< insert_frags_pos <<
			//"\nmax_pos " << fragments_big_->max_pos() <<
			//"\nfragset for positions = nr_frames " << fragments_big_->nr_frames() << std::endl;

			for ( core::fragment::ConstFrameIterator i = fragments_big_->begin(); i != fragments_big_->end(); ++i ) {
				core::Size position = (*i)->start();
				//std::cout << "position after iterator: " << position << std::endl;
				library_[position] =  **i;
				//library_[insert_frags_pos]= **i;
				//adjust position counter (for debugging stuff)
				//std::cout<< "position counter: " << insert_frags_pos << std::endl;
				insert_frags_pos++;

			}

			(*lowres_scorefxn_)(pose);
			protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *lowres_scorefxn_, temp_ ) );

			for ( core::Size i=1; i <= mc_cycles_; ++i ) {
				hybridize(pose, insert_pos_start, insert_pos_stop);
				(*lowres_scorefxn_)(pose);
				mc->boltzmann(pose);
			}
			mc->show_scores();
			mc->show_counters();
			mc->recover_low(pose);
		}

		TR.Debug << "final rms for fragment align: " << rms_ << std::endl;

		/// 3.c. finish with minimizing
		minimizer.run( pose, *extended_mm_, *min_scorefxn_, options_minilbfgs );

	}//finish interating through pieces that need to be smoothened

	/// 4. final minimization (optional)

	if ( extra_min_ ) {
		TR<< "final minimization" << std::endl;
		(*min_scorefxn_)(pose);   minimizer.run( pose, *extended_mm_, *min_scorefxn_, options_lbfgs );
	}
}

std::string
SegmentHybridizer::get_name() const {
	return SegmentHybridizerCreator::mover_name();
}


void
SegmentHybridizer::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap &data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*pose*/ ){

	TR << "initialized SegmentHybridizer mover " << std::endl;

	//jump_num_ = tag->getOption<core::Size>("jump_number", 1 );
	highres_scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, "fa_scorefxn", data )->clone();

	TR << "initialized SegmentHybridizer mover " << std::endl;

	big_   = tag->getOption <core::Size>( "frag_big" , 9 );
	//currently uses only 1 size
	small_ = tag->getOption <core::Size>( "frag_smal", 3 );
	nfrags_= tag->getOption<core::Size>( "nfrags", 50 );
	cartfrag_overlap_ = tag->getOption<core::Size>("cartfrag_overlap" , 2 );
	use_seq_ = tag->getOption<bool> ("use_seq", 0 );
	rms_ = tag->getOption<core::Real> ("rms_frags" , 0.5 );
	auto_mm_= tag->getOption<bool>("auto_mm", 1 );
	tries_=tag->getOption<int>("frag_tries" , 100 );
	mc_cycles_=tag->getOption<core::Size>("mc_cycles", 200 );
	temp_=tag->getOption<core::Real>("temp", 2.0);
	use_frags_=tag->getOption<bool>("use_frags", 1 );
	min_cycles_=tag->getOption<core::Size>("min_cycles", 5 );
	all_movable_=tag->getOption<bool>("all_movable", 0 );
	extra_min_=tag->getOption<bool>("extra_min", 0 );

	/// read areas that are supposed to be remodeled
	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );
	BOOST_FOREACH ( TagCOP const btag, branch_tags ) {

		if ( btag->getName() == "Span" ) { //need an assertion for the presence of these or at least for the option file
			std::string const beginS( btag->getOption<std::string>( "begin" ) );
			std::string const endS( btag->getOption<std::string>( "end" ) );

			extend_outside_ = btag->getOption<core::Size>( "extend_outside", 5);
			extend_inside_  =  btag->getOption<core::Size>( "extend_inside" , 1 );
			//   N_mm_ = btag->getOption<core::Size>( "N_mm", 2 );
			//   C_mm_ = btag->getOption<core::Size>( "C_mm", 2 );

			std::pair <std::string,std::string> segpair;
			segpair.first  = beginS;
			//std::cout  <<"parsing spans: " << beginS << " " <<endS << " and remdoel "<< std::endl;
			segpair.second = endS;
			seg_vector_.push_back( segpair ); // parse at runtime for possible length changes

		}//end seeds
	}//end b-tags
}//end parse tag

} //seeded_abinitio
} //protocols
