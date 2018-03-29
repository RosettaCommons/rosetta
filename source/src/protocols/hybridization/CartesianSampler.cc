// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Yifan Song
/// @author Frank DiMaio
/// @author Ray Wang

//keep first
#include <core/scoring/cryst/PhenixInterface.hh>

#include <protocols/hybridization/HybridizeProtocolCreator.hh>
#include <protocols/hybridization/CartesianSampler.hh>
#include <protocols/hybridization/FragmentBiasAssigner.hh>
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
#include <protocols/symmetry/SetupNCSMover.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/util.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/comparative_modeling/coord_util.hh>
#include <basic/datacache/DataMap.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/jd2/util.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>

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
#include <core/io/Remarks.hh>
#include <core/pose/selection.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>

// dump intermediate pdb
#include <core/io/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <utility/io/ozstream.hh>

#include <core/pose/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/func/SOGFunc.hh>

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
#include <core/scoring/func/USOGFunc.hh>

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
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>

#include <boost/unordered/unordered_map.hpp>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace hybridization {

using namespace core;

static basic::Tracer TR( "protocols.hybridization.CartesianSampler" );

/////////////
// creator
// XRW TEMP std::string
// XRW TEMP CartesianSamplerCreator::keyname() const {
// XRW TEMP  return CartesianSampler::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP CartesianSamplerCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new CartesianSampler );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CartesianSampler::mover_name() {
// XRW TEMP  return "CartesianSampler";
// XRW TEMP }

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
	nfrags_=25;
	ref_cst_weight_ = 1.0;
	ref_cst_maxdist_ = 12.0;
	debug_ = false;
	input_as_ref_ = false;
	fullatom_ = false;
	bbmove_ = false;
	recover_low_ = true;
	restore_csts_ = true;
	force_ss_='D';
	// ray added
	exclude_residues_ = false;
	include_residues_ = false;
	dump_pdb_ = false;
	dump_pdb_tag_ = "1";
	automode_scorecut_ = -0.5;
	rsd_wdw_to_refine_ = 0; //
	wdw_to_freeze_ = 0;
	freeze_endpoints_ = true;
	score_threshold_ = 9999;

	selection_bias_ = "none";
	cumulate_prob_ = false;

	// default scorefunction
	set_scorefunction ( core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" ) );
	set_fa_scorefunction ( core::scoring::ScoreFunctionFactory::create_score_function( "soft_rep" ) );
}

void
CartesianSampler::update_fragment_library_pointers() {
	core::Size nfragsets = fragments_.size();

	// map positions to fragments
	library_.resize( nfragsets );
	for ( int i=1; i<=(int)nfragsets; ++i ) {
		for ( core::fragment::ConstFrameIterator j = fragments_[i]->begin(); j != fragments_[i]->end(); ++j ) {
			core::Size position = (*j)->start();
			library_[i][position] = **j;
		}
	}
}

protocols::moves::MoverOP CartesianSampler::clone() const { return protocols::moves::MoverOP( new CartesianSampler( *this ) ); }
protocols::moves::MoverOP CartesianSampler::fresh_instance() const { return protocols::moves::MoverOP( new CartesianSampler ); }

// get frag->pose transform, return RMS
// default overlap_ = 2
core::Real
CartesianSampler::get_transform(
	core::pose::Pose const &pose,
	core::pose::Pose const &frag,
	core::Size startpos,
	core::Vector &preT,
	core::Vector &postT,
	numeric::xyzMatrix< core::Real > &R
){
	int aln_len = overlap_;
	int len = frag.size();
	if ( frag.residue(len).aa() == core::chemical::aa_vrt ) len--;
	ObjexxFCL::FArray1D< core::Real > ww( 2*4*aln_len, 1.0 );
	ObjexxFCL::FArray2D< core::Real > uu( 3, 3, 0.0 ), init_coords( 3, 2*4*aln_len ), final_coords( 3, 2*4*aln_len );
	preT = numeric::xyzVector< core::Real >(0,0,0);
	postT = numeric::xyzVector< core::Real >(0,0,0);

	bool nterm = ! freeze_endpoints_ && ( (startpos == 1) || pose.fold_tree().is_cutpoint(startpos-1) );
	bool cterm = ! freeze_endpoints_ && ( pose.fold_tree( ).is_cutpoint(startpos+len-1) );

	// grab coords from input pose
	for ( int ii=-aln_len; ii<aln_len; ++ii ) {
		int i = (ii>=0) ? (nterm?len-ii-1:ii) : (cterm?-ii-1:len+ii);
		numeric::xyzVector< core::Real > x_1 = pose.residue(startpos+i).atom(1).xyz();
		numeric::xyzVector< core::Real > x_2 = pose.residue(startpos+i).atom(2).xyz();
		numeric::xyzVector< core::Real > x_3 = pose.residue(startpos+i).atom(3).xyz();
		numeric::xyzVector< core::Real > x_4 = pose.residue(startpos+i).atom(4).xyz();
		postT += x_1+x_2+x_3+x_4;
		for ( int j=0; j<3; ++j ) {
			init_coords(j+1,4*(ii+aln_len)+1) = x_1[j];
			init_coords(j+1,4*(ii+aln_len)+2) = x_2[j];
			init_coords(j+1,4*(ii+aln_len)+3) = x_3[j];
			init_coords(j+1,4*(ii+aln_len)+4) = x_4[j];
		}
	}
	postT /= 2.0*4.0*aln_len;
	for ( int ii=0; ii<2*4*aln_len; ++ii ) {
		for ( int j=0; j<3; ++j ) init_coords(j+1,ii+1) -= postT[j];
	}

	// grab new coords from frag
	for ( int ii=-aln_len; ii<aln_len; ++ii ) {
		int i = (ii>=0) ? (nterm?len-ii-1:ii) : (cterm?-ii-1:len+ii);
		numeric::xyzVector< core::Real > x_1 = frag.residue(i+1).atom(1).xyz();
		numeric::xyzVector< core::Real > x_2 = frag.residue(i+1).atom(2).xyz();
		numeric::xyzVector< core::Real > x_3 = frag.residue(i+1).atom(3).xyz();
		numeric::xyzVector< core::Real > x_4 = frag.residue(i+1).atom(4).xyz();
		preT += x_1+x_2+x_3+x_4;
		for ( int j=0; j<3; ++j ) {
			final_coords(j+1,4*(ii+aln_len)+1) = x_1[j];
			final_coords(j+1,4*(ii+aln_len)+2) = x_2[j];
			final_coords(j+1,4*(ii+aln_len)+3) = x_3[j];
			final_coords(j+1,4*(ii+aln_len)+4) = x_4[j];
		}
	}
	preT /= 2.0*4.0*aln_len;
	for ( int ii=0; ii<2*4*aln_len; ++ii ) {
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

// apply endpoint csts from pose to frag
void
CartesianSampler::apply_fragcsts(
	core::pose::Pose &working_frag,
	core::pose::Pose const &pose,
	core::Size startpos
){
	using namespace core::scoring::constraints;

	working_frag.remove_constraints();
	int len = working_frag.size();
	if ( working_frag.residue(len).aa() == core::chemical::aa_vrt ) len--;

	bool nterm = ! freeze_endpoints_ && ( (startpos == 1) || pose.fold_tree().is_cutpoint(startpos-1) );
	bool cterm = ! freeze_endpoints_ && ( pose.fold_tree( ).is_cutpoint(startpos+len-1) );

	if ( !nterm ) {
		for ( core::uint j = 0; j < overlap_; ++j ) {
			for ( core::uint i = 1; i <= 3; ++i ) {
				core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
				working_frag.add_constraint(
					scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new CoordinateConstraint(
					core::id::AtomID(i,j+1),
					core::id::AtomID(2,working_frag.size()),
					pose.residue(startpos+j).atom(i).xyz(),
					fx
					) ) )
				);
			}
		}
	}
	if ( !cterm ) {
		for ( int j=len-overlap_; j<len; ++j ) {
			for ( int i=1; i<=3; ++i ) {
				core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
				working_frag.add_constraint(
					scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new CoordinateConstraint(
					core::id::AtomID(i,j+1),
					core::id::AtomID(2,working_frag.size()),
					pose.residue(startpos+j).atom(i).xyz(),
					fx
					) ) )
				);
			}
		}
	}
}


// transform fragment
void
CartesianSampler::apply_transform(
	core::pose::Pose &frag,
	core::Vector const &preT,
	core::Vector const &postT,
	numeric::xyzMatrix< core::Real > const &R
){
	// apply rotation to ALL atoms
	// x_i' <- = R*x_i + com1;
	for ( Size i = 1; i <= frag.size(); ++i ) {
		for ( Size j = 1; j <= frag.residue_type(i).natoms(); ++j ) {
			core::id::AtomID id( j, i );
			frag.set_xyz( id, R * ( frag.xyz(id) - preT) + postT );
		}
	}
}


bool
CartesianSampler::apply_frame(
	core::pose::Pose & pose,
	core::fragment::Frame &frame
){
	using core::pack::task::operation::TaskOperationCOP;

	core::Size start = frame.start(),len = frame.length();
	runtime_assert( overlap_>=1 && overlap_<=len/2); // at least two-residue overlap, 4 mer at least

	// set the frame's insert point to 1
	frame.shift_to(1);

	// see if the pose has NCS
	symmetry::NCSResMappingOP ncs;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ) {
		ncs = ( utility::pointer::static_pointer_cast< symmetry::NCSResMapping > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ));
	}

	protocols::simple_moves::SwitchResidueTypeSetMover to_fa("fa_standard");
	protocols::simple_moves::SwitchResidueTypeSetMover to_cen("centroid");

	// make subpose at this position
	// we assume the frame does not cross a jump
	core::pose::Pose frag;
	for ( int i=0; i<(int)len; ++i ) frag.append_residue_by_bond( pose.residue( start+i ) );
	for ( int i=0; i<(int)len; ++i ) {
		if ( frag.residue(i+1).aa() == chemical::aa_cys && frag.residue(i+1).has_variant_type( chemical::DISULFIDE ) ) {
			TR << "temp remove dslf at " << i+1 << std::endl;
			conformation::change_cys_state( i+1, "CYS", frag.conformation() );
		}
	}

	for ( int i=0; i<(int)len; ++i ) core::conformation::idealize_position(i+1, frag.conformation());

	core::Vector preT, postT;
	numeric::xyzMatrix< core::Real > R;
	core::Size maxtries,frag_toget=0;
	core::Real rms=rms_cutoff_;

	// set up minimizer
	core::optimization::AtomTreeMinimizer rbminimizer;
	core::optimization::MinimizerOptions options_rb( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
	options_rb.max_iter(20);
	core::kinematics::MoveMap mm_rb;
	mm_rb.set_bb ( bbmove_ );
	mm_rb.set_chi ( true );
	mm_rb.set_jump ( true );

	// set up packer
	core::scoring::ScoreFunctionOP nonsymm_fa_scorefxn = core::scoring::symmetry::asymmetrize_scorefunction(*fa_scorefxn_);
	core::scoring::ScoreFunctionOP nonsymm_cen_scorefxn = core::scoring::symmetry::asymmetrize_scorefunction(*scorefxn_);

	if ( bbmove_ && nonsymm_fa_scorefxn->get_weight( core::scoring::coordinate_constraint ) == 0 ) {
		nonsymm_fa_scorefxn->set_weight( core::scoring::coordinate_constraint , 1 );
	}
	core::pack::task::TaskFactoryOP main_task_factory( new core::pack::task::TaskFactory );
	main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
	protocols::minimization_packing::PackRotamersMoverOP pack_mover( new protocols::minimization_packing::PackRotamersMover );
	pack_mover->task_factory( main_task_factory );
	pack_mover->score_function( nonsymm_fa_scorefxn );

	if ( selection_bias_ == "density" ) {
		// scorefunctions
		if ( nonsymm_fa_scorefxn->get_weight( core::scoring::elec_dens_fast ) == 0 ) {
			nonsymm_fa_scorefxn->set_weight( core::scoring::elec_dens_fast , 20 );
		}
		core::scoring::ScoreFunctionOP densonly( new core::scoring::ScoreFunction() );
		densonly->set_weight( core::scoring::elec_dens_fast, 5.0 );

		// prepare fragment
		core::pose::addVirtualResAsRoot(frag);
		if ( !fullatom_ ) to_fa.apply( frag );

		core::Size nattempts = 0;
		core::Real best_dens_score = 1e30, best_rms = 1e30, best_rms_aftermin = 1e30;
		for ( core::uint i = 1; i <= frame.nr_frags(); ++i ) {
			core::pose::Pose working_frag = frag;
			frame.apply( i, working_frag );
			rms = get_transform( pose,  working_frag,  start, preT, postT, R); // rms straight from comparing frag to pose at region with nres overlap_
			if ( rms <= rms_cutoff_ ) {
				nattempts++;

				// orient
				apply_transform( working_frag, preT, postT, R);

				// cenmin+repack+min
				to_cen.apply( working_frag );
				(*densonly)(working_frag);
				rbminimizer.run( working_frag, mm_rb, *densonly, options_rb );
				to_fa.apply( working_frag );
				pack_mover->apply( working_frag );

				core::Real dens_score;
				if ( bbmove_ ) {
					apply_fragcsts( working_frag, pose, start ); // apply csts from pose to frag at endpoints
					(*nonsymm_fa_scorefxn)(working_frag);
					rbminimizer.run( working_frag, mm_rb, *nonsymm_fa_scorefxn, options_rb );
					dens_score = (*nonsymm_fa_scorefxn) (working_frag);
				} else {
					(*densonly)(working_frag);
					rbminimizer.run( working_frag, mm_rb, *densonly, options_rb );
					dens_score = (*densonly) (working_frag);
				}

				if ( dens_score<best_dens_score ) {
					best_rms_aftermin = get_transform( pose,  working_frag,  start, preT, postT, R); // rms after minimization
					best_rms = rms; // this rms is before minimization, should return after minimization
					best_dens_score = dens_score;
					frag = working_frag;
				}
				if ( nattempts >=25 ) break;  //fpd within helices often all fragments will match

			} // if rms <= rms_cutoff_
		} // try all the frags in the position until some have rms < cutoff

		if ( nattempts > 0 ) {
			TR << "Best frag ( out of "<< nattempts << ") with score="  << best_dens_score << "  rms=" << best_rms << "  rms_aftermin=" << best_rms_aftermin << std::endl;
			if ( !fullatom_ ) to_cen.apply( frag );
			for ( Size i = 0; i < len; ++i ) {
				bool dslf_i = ( pose.residue(start+i).aa() == chemical::aa_cys && pose.residue(start+i).has_variant_type( chemical::DISULFIDE ) );

				// remember dslf connectivity
				if ( dslf_i ) {
					for ( int j=1; j<=(int)pose.residue(start+i).natoms(); ++j ) {  // copy everything but HG
						TR << "restore dslf at " << start+i << std::endl;
						pose.set_xyz( core::id::AtomID( j,start+i ), frag.residue(i+1).xyz( j ) );
					}
				} else {
					pose.replace_residue( start+i, frag.residue(i+1), false );
				}
			}
		} else { // best scoring fragment saved
			TR << "No acceptable fragments" << std::endl;
			frame.shift_to(start);
			return false;
		}
	} else if ( selection_bias_ == "none" ) {
		maxtries = frame.nr_frags(); // bias-dependent
		core::uint tries;
		bool frag_chosen=false;
		for ( tries = 0; tries < maxtries && !frag_chosen; ++tries ) {
			frag_toget = numeric::random::random_range( 1, frame.nr_frags() );
			frame.apply( frag_toget, frag );
			rms = get_transform( pose,  frag,  start, preT, postT, R);
			frag_chosen = (rms < rms_cutoff_);

			// minimize backbone fragment
			if ( frag_chosen && bbmove_ ) {
				apply_fragcsts( frag, pose, start );
				(*nonsymm_cen_scorefxn)(frag);
				rbminimizer.run( frag, mm_rb, *nonsymm_cen_scorefxn, options_rb );
				frag_chosen = (rms < 0.5);
			}
		}
		TR << "after " << tries+1 << " tries, fragment " << frag_toget << " chosen with RMS = " << rms << std::endl;

		apply_transform( frag, preT, postT, R);

		// set up packer
		core::pack::task::TaskFactoryOP main_task_factory( new core::pack::task::TaskFactory );
		main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
		protocols::minimization_packing::PackRotamersMoverOP pack_mover( new protocols::minimization_packing::PackRotamersMover );
		pack_mover->task_factory( main_task_factory );
		pack_mover->score_function( nonsymm_fa_scorefxn );

		for ( Size i = 0; i < len; ++i ) {
			bool dslf_i = ( pose.residue(start+i).aa() == chemical::aa_cys && pose.residue(start+i).has_variant_type( chemical::DISULFIDE ) );

			// remember dslf connectivity
			if ( dslf_i ) {
				for ( int j=1; j<=(int)pose.residue(start+i).natoms(); ++j ) {  // copy everything but HG
					TR << "restore dslf at " << start+i << std::endl;
					pose.set_xyz( core::id::AtomID( j,start+i ), frag.residue(i+1).xyz( j ) );
				}
			} else {
				pose.replace_residue( start+i, frag.residue(i+1), false );
			}
		}

		// if fullatom do a repack here
		if ( fullatom_ ) {
			if ( core::pose::symmetry::is_symmetric(pose) ) {
				protocols::minimization_packing::PackRotamersMoverOP sympack_mover( new protocols::minimization_packing::symmetry::SymPackRotamersMover );
				sympack_mover->task_factory( main_task_factory );
				sympack_mover->score_function( fa_scorefxn_ );
				sympack_mover->apply( pose );
			} else {
				pack_mover->apply( pose );
			}
		}
	} else {
		utility_exit_with_message("Unrecognized fragbias!");
	}

	// apply to NCS-symmetric copies
	if ( ncs ) {
		for ( core::uint j = 1; j <= ncs->ngroups(); ++j ) {
			bool all_are_mapped = true;
			for ( Size k= 0; k< len && all_are_mapped; ++k ) {
				all_are_mapped &= (ncs->get_equiv( j,start+k )!=0);
			}
			if ( !all_are_mapped ) continue;

			core::Size remap_start = ncs->get_equiv( j, start );
			get_transform( pose,  frag,  remap_start, preT, postT, R);
			apply_transform( frag, preT, postT, R);
			for ( Size i = 0; i < len; ++i ) {
				bool dslf_i = ( pose.residue(remap_start+i).aa() == chemical::aa_cys && pose.residue(remap_start+i).has_variant_type( chemical::DISULFIDE ) );
				if ( dslf_i ) {
					for ( int k=1; k<=(int)pose.residue(remap_start+i).natoms(); ++k ) {  // copy everything but HG
						TR << "(ncs) restore dslf at " << remap_start+i << std::endl;
						pose.set_xyz( core::id::AtomID( k,remap_start+i ), frag.residue(i+1).xyz( k ) );
					}
				} else {
					pose.replace_residue( remap_start+i, frag.residue(i+1), false );
				}
			}
		}
	}

	// restore frame
	frame.shift_to(start);
	return true;
}


// when "reference_model" specified in the options, apply constraints from ref_pose into current pose
// for residues in user selection, ignore <- does it make sense though (ray), shouldn't we ignore it by using all the selected residues in FragmentBiasAssigner?
void
CartesianSampler::apply_constraints(
	core::pose::Pose &pose
){
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		auto & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	core::Size MINSEQSEP = 8;
	core::Real MAXDIST = ref_cst_maxdist_;
	core::Real COORDDEV = 0.39894;

	// user_pos_ rewrite rw
	bool user_rebuild = ( std::find( fragment_bias_strategies_.begin(), fragment_bias_strategies_.end(), "user" ) !=fragment_bias_strategies_.end() );

	core::Size nres_tgt = get_num_residues_nonvirt( pose );

	utility::vector1< utility::vector1< core::Real > > tgt_dists(nres_tgt);
	utility::vector1< utility::vector1< core::Real > > tgt_weights(nres_tgt);

	std::set< core::Size > user_pos;
	if ( user_pos_ ) {
		user_pos = core::select::get_residue_set_from_subset( user_pos_->apply( pose ) );
	}
	for ( core::Size j=1; j < ref_model_.size(); ++j ) {
		if ( !ref_model_.residue_type(j).is_protein() ) continue;

		if ( user_rebuild && user_pos.find(j) != user_pos.end() ) continue;
		if ( user_rebuild && loops_ && loops_->is_loop_residue(j) ) continue;

		for ( core::Size k=j+1; k < ref_model_.size(); ++k ) {
			if ( !ref_model_.residue_type(k).is_protein() ) continue;
			if ( ref_model_.pdb_info()->number(k) - ref_model_.pdb_info()->number(j) < (int)MINSEQSEP ) continue;

			if ( user_rebuild && user_pos.find(k) != user_pos.end() ) continue;
			if ( user_rebuild && loops_ && loops_->is_loop_residue(k) ) continue;

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
				if ( symm_info && !symm_info->bb_is_independent( tgt_resid_j ) ) continue;
				if ( symm_info && !symm_info->bb_is_independent( tgt_resid_k ) ) continue;

				FuncOP fx( new ScalarWeightedFunc( ref_cst_weight_, FuncOP( new USOGFunc( dist, COORDDEV ) ) ) );
				pose.add_constraint(
					scoring::constraints::ConstraintCOP( new AtomPairConstraint( core::id::AtomID(2,tgt_resid_j), core::id::AtomID(2,tgt_resid_k), fx ) )
				);
			}
		}
	}
}


void
CartesianSampler::
compute_fragment_bias(
	Pose & pose
){
	TR << "compute fragment bias" << std::endl;
	FragmentBiasAssigner frag_bias_assigner( pose );

	if ( fragment_bias_strategies_.empty() ) {
		fragment_bias_strategies_.push_back("uniform");
	}

	frag_bias_assigner.set_rsd_wdw_to_assign_prob( rsd_wdw_to_refine_ );
	frag_bias_assigner.set_score_threshold( score_threshold_ );

	// select residues
	if ( std::find( fragment_bias_strategies_.begin(), fragment_bias_strategies_.end(), "uniform" ) !=fragment_bias_strategies_.end() ) {
		frag_bias_assigner.uniform();
	} else if ( std::find( fragment_bias_strategies_.begin(), fragment_bias_strategies_.end(), "auto" ) !=fragment_bias_strategies_.end() ) {
		frag_bias_assigner.automode( pose, automode_scorecut_ );
	} else {
		if ( cumulate_prob_ ) {
			frag_bias_assigner.cumulate_probability();
		}

		if ( std::find( fragment_bias_strategies_.begin(), fragment_bias_strategies_.end(), "user" ) !=fragment_bias_strategies_.end() ) {
			std::set< core::Size > user_pos;
			if ( user_pos_ ) {
				user_pos = core::select::get_residue_set_from_subset( user_pos_->apply(pose) );
			}
			frag_bias_assigner.user( user_pos, loops_ );
		}

		if ( std::find( fragment_bias_strategies_.begin(), fragment_bias_strategies_.end(), "density" ) !=fragment_bias_strategies_.end() ) {
			frag_bias_assigner.density( pose );
		}

		if ( std::find( fragment_bias_strategies_.begin(), fragment_bias_strategies_.end(), "density_nbr" ) !=fragment_bias_strategies_.end() ) {
			frag_bias_assigner.density_nbr( pose );
		}

		if ( std::find( fragment_bias_strategies_.begin(), fragment_bias_strategies_.end(), "geometry" ) !=fragment_bias_strategies_.end() ) {
			frag_bias_assigner.geometry( pose );
		}

		if ( std::find( fragment_bias_strategies_.begin(), fragment_bias_strategies_.end(), "rama" ) !=fragment_bias_strategies_.end() ) {
			frag_bias_assigner.rama( pose );
		}

		if ( std::find( fragment_bias_strategies_.begin(), fragment_bias_strategies_.end(), "bfactors" ) !=fragment_bias_strategies_.end() ) {
			frag_bias_assigner.bfactors( pose );
		}

		if ( std::find( fragment_bias_strategies_.begin(), fragment_bias_strategies_.end(), "chainbreak" ) !=fragment_bias_strategies_.end() ) {
			frag_bias_assigner.chainbreak( pose );
		}
	}

	// remove residues form the frag_bias container - should probably consider window as well
	if ( include_residues_ ) {
		std::set< core::Size > residues_to_include;
		if ( user_pos_ ) {
			residues_to_include = core::select::get_residue_set_from_subset( residues_to_include_->apply(pose) );
		}
		frag_bias_assigner.include_residues( residues_to_include );
	}

	// remove residues form the frag_bias container - should probably consider window as well
	if ( exclude_residues_ ) {
		std::set< core::Size > residues_to_exclude;
		if ( user_pos_ ) {
			residues_to_exclude = core::select::get_residue_set_from_subset( residues_to_exclude_->apply(pose) );
		}
		frag_bias_assigner.exclude_residues( residues_to_exclude );
	}

	// set residues to freeze (residues to chew back from the missing density jump)
	frag_bias_assigner.set_wdw_to_freeze( wdw_to_freeze_ );

	// compute frag bias
	frag_bias_assigner.compute_frag_bias( frag_bias_, pose, fragments_ );

}

void
CartesianSampler::
apply(
	core::pose::Pose & pose
){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose::datacache;

	// dump pdb right before doing anything
	if ( dump_pdb_ ) {
		std::string outfile = ("intermediate_" + dump_pdb_tag_ + ".pdb");
		utility::io::ozstream out( outfile );
		core::io::pdb::dump_pdb( pose, out );
	}

	// autogenerate fragments if they are not loaded yet
	if ( fragments_.size() == 0 ) {
		if ( frag_sizes_.size() == 0 ) frag_sizes_.push_back(9); // default is 9-mers only

		for ( core::uint i = 1; i <= frag_sizes_.size(); ++i ) {
			/*if ( (std::find( fragment_bias_strategies_.begin(), fragment_bias_strategies_.end(), "user" ) !=fragment_bias_strategies_.end())
			&& user_pos_.size() > 0) { // pick only on user selected positions
			fragments_.push_back( create_fragment_set_no_ssbias(pose, user_pos_, frag_sizes_[i], nfrags_, force_ss_) );
			}
			else {
			fragments_.push_back( create_fragment_set_no_ssbias(pose, frag_sizes_[i], nfrags_, force_ss_) );
			}
			fragments_.push_back( create_fragment_set_no_ssbias(pose, frag_sizes_[i], nfrags_, force_ss_) );
			}*/

			TR << "frag_sizes_: " << frag_sizes_[i] << " nfrags: " << nfrags_ <<  std::endl;
			fragments_.push_back( create_fragment_set_no_ssbias(pose, frag_sizes_[i], nfrags_, force_ss_) );
			update_fragment_library_pointers();
		}
	}

	if ( fullatom_ ) scorefxn_ = fa_scorefxn_;
	if ( !mc_scorefxn_ ) mc_scorefxn_ = scorefxn_;

	// see if the pose has NCS
	symmetry::NCSResMappingOP ncs;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ) {
		ncs = ( utility::pointer::static_pointer_cast< symmetry::NCSResMapping > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ));
	}

	// using the current fragment_bias_strategy_, compute bias
	// do this from fullatom (if the input pose is fullatom)
	compute_fragment_bias(pose);

	// add heteroatom constraint into it
	//if (add_hetatm_)
	//add_non_protein_cst(pose, hetatm_self_cst_weight_, hetatm_prot_cst_weight_);

	// number of protein residues in ASU
	core::Size nres = pose.size();
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		auto & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres = symm_info->num_independent_residues();
	}
	while ( !pose.residue(nres).is_protein() ) nres--;


	bool fullatom_input = pose.is_fullatom();
	protocols::moves::MoverOP restore_sc;
	if ( fullatom_input && !fullatom_ ) {
		restore_sc = protocols::moves::MoverOP( new protocols::simple_moves::ReturnSidechainMover( pose ) );
		protocols::moves::MoverOP tocen( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );
		tocen->apply( pose );
	} else if ( !fullatom_input && fullatom_ ) {
		utility_exit_with_message("ERROR! Expected fullatom input.");
	}

	// constraints to reference model
	// save current constraint set
	core::scoring::constraints::ConstraintSetOP saved_csts = pose.constraint_set()->clone();
	if ( input_as_ref_ ) ref_model_ = pose;
	apply_constraints( pose ); // rw: this is not in used unless you specified reference_model

	// stepwise minimizer
	core::optimization::MinimizerOptions options_minilbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
	options_minilbfgs.max_iter(nminsteps_);

	// ray, why do we need this mm?
	// to do ... make this parsable
	core::optimization::CartesianMinimizer minimizer;
	core::kinematics::MoveMap mm;
	mm.set_bb( true ); mm.set_chi( true ); mm.set_jump( true );

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::pose::symmetry::make_symmetric_movemap( pose, mm );
	}

	Pose pose_in = pose;
	(*scorefxn_)(pose);
	protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *mc_scorefxn_, temp_ ) );

	if ( TR.Debug.visible() ) {
		scorefxn_->show_line_headers(TR);
		scorefxn_->show_line(TR,pose);
	}

	std::string base_name = protocols::jd2::current_input_tag();
	utility::vector1< std::string > temp_out_names= utility::split( base_name );
	utility::file::FileName out_name = utility::file::combine_names( temp_out_names );
	base_name = out_name.base();
	std::string outname = option[ out::prefix ]()+base_name+option[ out::suffix ]();

	for ( int n=1; n<=(int)ncycles_; ++n ) {
		// apply_frame to a random selected insert pos
		bool success=false;
		int try_count=1000;
		int insert_pos=-999; //Adding initialization value for quick fix of compiler errors - this might not be the most sensible way of handling this.
		core::Size i_frag_set=999; //Adding initialization value for quick fix of compiler errors - this might not be the most sensible way of handling this.

		while ( !success && --try_count>0 ) {
			// pick fragment set
			i_frag_set = numeric::random::random_range(1, fragments_.size());

			// pick insertion position
			insert_pos = (int)frag_bias_[i_frag_set].random_sample(numeric::random::rg());
			int ntries=50;
			while ( library_[i_frag_set].find(insert_pos) == library_[i_frag_set].end() && --ntries>0 )
					insert_pos = (int)frag_bias_[i_frag_set].random_sample(numeric::random::rg());

			if ( library_[i_frag_set].find(insert_pos) == library_[i_frag_set].end() ) {
				// ray: why does this happen though?
				// frag_bias_ gives positions where there is no fragments selected?
				TR << "ERROR! Unable to find allowed fragment insert at " << insert_pos << " after 50 attempts.  Continuing." << std::endl;
				continue;
			}

			success = apply_frame(pose, library_[i_frag_set][insert_pos]) ;
		}// pick and insert_pos and apply_frame (with the best density_score fragment saved)

		if ( insert_pos == -999 ) {
			TR.Error << "Cannot find insert position!" << std::endl;
		}

		// restricted movemap
		// TODO we should check of chainbreaks
		// TODO min window extension (curr 6) should be a parameter
		core::kinematics::MoveMap mm_local;
		int start_move = std::max(1,insert_pos-6);
		int stop_move = std::min(nres,insert_pos+library_[i_frag_set][insert_pos].length()+5);
		for ( int i=start_move; i<=stop_move; ++i ) {
			mm_local.set_bb(i,true);
			mm_local.set_chi(i,true);
			// ncs
			if ( ncs ) {
				for ( core::uint j = 1; j <= ncs->ngroups(); ++j ) {
					core::Size remap_i = ncs->get_equiv( j, i );
					if ( remap_i!=0 ) {
						mm_local.set_bb(remap_i,true);
						mm_local.set_chi(remap_i,true);
					}
				}
			}
		} // assigned movemap

		// min + MC
		(*scorefxn_)(pose);
		if ( TR.Debug.visible() ) {
			scorefxn_->show_line(TR,pose);
			TR << std::endl;
		}
		minimizer.run( pose, mm_local, *scorefxn_, options_minilbfgs );
		(*scorefxn_)(pose);
		if ( TR.Debug.visible() ) {
			scorefxn_->show_line(TR,pose);
			TR << std::endl;
		}

		(*mc_scorefxn_)(pose);
		bool accept = mc->boltzmann( pose );

		mc->show_scores();
		if ( accept ) {
			TR << "Insert at " << insert_pos << " accepted!" << std::endl;
		}
		if ( debug_ ) {
			static int ctr=1;
			pose.dump_pdb( outname+"_acc_"+ObjexxFCL::right_string_of( ctr, 6, '0' )+".pdb" );
			ctr++;
		} else {
			TR << "Insert at " << insert_pos << " rejected!" << std::endl;
		}

	}// ncycles_ to apply_frame

	mc->show_scores();
	mc->show_counters();
	if ( recover_low_ ) mc->recover_low(pose); // default: true

	if ( fullatom_input && !fullatom_ ) {
		protocols::moves::MoverOP tofa( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD ) );
		tofa->apply( pose );
		restore_sc->apply( pose );
	}

	// reapply old constraints
	if ( restore_csts_ ) {
		pose.remove_constraints();
		pose.constraint_set( saved_csts );
	}
}


// parse_my_tag
void
CartesianSampler::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &
){
	using namespace core::scoring;

	// scorefunction
	if ( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn" ) );
		set_scorefunction ( data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_name ) );
	}

	// fullatom scorefunction
	if ( tag->hasOption( "fascorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "fascorefxn" ) );
		set_fa_scorefunction ( data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_name ) );
	}

	// mc scorefunction
	if ( tag->hasOption( "mcscorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "mcscorefxn" ) );
		mc_scorefxn_ = data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_name );
	}

	// options
	if ( tag->hasOption( "debug" ) ) {
		debug_ = tag->getOption< bool >( "debug" );
	}
	if ( tag->hasOption( "fullatom" ) ) {
		fullatom_ = tag->getOption< bool >( "fullatom" );
	}
	if ( tag->hasOption( "bbmove" ) ) {
		bbmove_ = tag->getOption< bool >( "bbmove" );
	}
	if ( tag->hasOption( "temp" ) ) {
		temp_ = tag->getOption< core::Real >( "temp" );
	}
	if ( tag->hasOption( "rms" ) ) {
		rms_cutoff_ = tag->getOption< core::Real >( "rms" );
	}
	if ( tag->hasOption( "nminsteps" ) ) {
		nminsteps_ = tag->getOption< core::Size >( "nminsteps" );
	}
	if ( tag->hasOption( "overlap" ) ) {
		overlap_ = tag->getOption<core::Size>( "overlap" );
	}
	if ( tag->hasOption( "ncycles" ) ) {
		ncycles_ = tag->getOption<core::Size>( "ncycles" );
	}
	if ( tag->hasOption( "fragbias" ) ) {
		selection_bias_ = tag->getOption< std::string >( "fragbias" );
	}

	recover_low_ = tag->getOption< bool >( "recover_low" , true );
	force_ss_    = tag->getOption< char >( "force_ss" , 'D' );

	////////////////////////////////////////////////////////////////////////////////
	// residue selection strategy
	if ( tag->hasOption( "strategy" ) ) {
		fragment_bias_strategies_ = utility::string_split( tag->getOption<std::string>("strategy"), ',' );
	}
	// to control window of residues to refine for strategy: user, rama, geometry, auto and "residues_to_include"
	if ( tag->hasOption( "rsd_wdw_to_refine" ) ) {
		rsd_wdw_to_refine_ = tag->getOption<core::Size>( "rsd_wdw_to_refine" );
	}
	if ( tag->hasOption( "score_threshold" ) ) {
		score_threshold_ = tag->getOption<core::Real>( "score_threshold" );
	}
	if ( tag->hasOption( "cumulate_prob" ) ) {
		cumulate_prob_ = tag->getOption<bool>( "cumulate_prob" );
	}
	if ( tag->hasOption( "automode_scorecut" ) ) {
		automode_scorecut_ = tag->getOption<core::Real>( "automode_scorecut" );
	}
	if ( tag->hasOption( "dump_pdb" ) ) {
		dump_pdb_ = tag->getOption<bool>( "dump_pdb" );
	}
	if ( tag->hasOption( "dump_pdb_tag" ) ) {
		dump_pdb_tag_ = tag->getOption<std::string>( "dump_pdb_tag" );
	}

	if ( tag->hasOption( "wdw_to_freeze" ) ) {
		wdw_to_freeze_ = tag->getOption<int>( "wdw_to_freeze" );
	}
	if ( tag->hasOption( "freeze_endpoints" ) ) {
		freeze_endpoints_ = tag->getOption<bool>( "freeze_endpoints" );
	}

	// user-specified residues
	runtime_assert( !( tag->hasOption( "residues" ) && tag->hasOption( "loops_in" ) ));  // one or the other
	if ( tag->hasOption( "residues" ) ) {
		user_pos_ = core::pose::get_resnum_selector( tag, "residues" );
	}
	if ( tag->hasOption( "loops_in" ) ) {
		std::string looptag = tag->getOption<std::string>( "loops_in" );
		runtime_assert( data.has( "loops", looptag ) );
		loops_ = data.get_ptr<protocols::loops::Loops>( "loops", looptag );
	}
	// residues_to_exclude duing modeling, prob assigned to be zero, should check frame as well
	if ( tag->hasOption( "residues_to_exclude" ) ) {
		residues_to_exclude_ = core::pose::get_resnum_selector( tag, "residues_to_exclude" );
		exclude_residues_=true;
	}
	// residues_to_include duing modeling, prob assigned to be one
	if ( tag->hasOption( "residues_to_include" ) ) {
		residues_to_include_ = core::pose::get_resnum_selector( tag, "residues_to_include" );
		include_residues_=true;
	}

	////////////////////////////////////////////////////////////////////////////////
	// csts from reference model
	if ( tag->hasOption( "reference_model" ) ) {
		std::string ref_model_pdb = tag->getOption<std::string>( "reference_model" );
		if ( ref_model_pdb != "input" ) {
			core::import_pose::pose_from_file( ref_model_, ref_model_pdb , core::import_pose::PDB_file);
		} else {
			input_as_ref_ = true;
		}
	}
	if ( tag->hasOption( "reference_cst_wt" ) ) {
		ref_cst_weight_ = tag->getOption< core::Real >( "reference_cst_wt" );
	}
	if ( tag->hasOption( "reference_cst_maxdist" ) ) {
		ref_cst_maxdist_ = tag->getOption< core::Real >( "reference_cst_maxdist" );
	}

	////////////////////////////////////////////////////////////////////////////////
	// fragments
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for ( tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() == "Fragments" ) {
			using namespace core::fragment;
			if ( tag->hasOption( "nfrags" ) ) {
				nfrags_ = tag->getOption< core::Size >( "nfrags" );
				fragments_.push_back( FragmentIO( nfrags_, 1, true ).read_data( (*tag_it)->getOption<std::string>( "fragfile" )  ) );
				TR << "top " << nfrags_ << " has been used while reading fragment file: "
					<< (*tag_it)->getOption<std::string>( "fragfile" ) << std::endl;
			} else {
				fragments_.push_back( FragmentIO().read_data( (*tag_it)->getOption<std::string>( "fragfile" )  ) );
				TR << "all fragments has been used while reading fragment file: "
					<< (*tag_it)->getOption<std::string>( "fragfile" ) << std::endl;
			}
		}
	}

	// if no Fragments tag is given, nfrags is defining how many fragments per pos for picking autofrag
	nfrags_ = tag->getOption< core::Size >( "nfrags", 25 );

	// autofragments
	if ( tag->hasOption( "fraglens" ) ) {
		utility::vector1<std::string> fraglens = utility::string_split( tag->getOption< std::string >( "fraglens" ), ',' );
		for ( core::Size i=1; i<=fraglens.size(); ++i ) {
			frag_sizes_.push_back( atoi(fraglens[i].c_str()) );
		}
	}
	update_fragment_library_pointers();


	////////////////////////////////////////////////////////////////////////////////
	// hetatm
	/*
	if( tag->hasOption( "add_hetatm" ) )
	add_hetatm_ = tag->getOption< bool >( "add_hetatm" );
	if( tag->hasOption( "hetatm_cst_weight" ) )
	hetatm_self_cst_weight_ = tag->getOption< core::Real >( "hetatm_cst_weight" );
	if( tag->hasOption( "hetatm_to_protein_cst_weight" ) )
	hetatm_prot_cst_weight_ = tag->getOption< core::Real >( "hetatm_to_protein_cst_weight" );
	*/
}

std::string CartesianSampler::get_name() const {
	return mover_name();
}

std::string CartesianSampler::mover_name() {
	return "CartesianSampler";
}

std::string CartesianSampler_subelement_ct_name( std::string const & name ) {
	return "CartesianSampler_subelement_" + name + "Type";
}

void CartesianSampler::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "scorefxn", xs_string, "(centroid) scorefunction")
		+ XMLSchemaAttribute( "fascorefxn", xs_string, "(fullatom) scorefunction")
		+ XMLSchemaAttribute( "mcscorefxn", xs_string, "Monte Carlo (?) scorefunction")
		+ XMLSchemaAttribute( "debug", xsct_rosetta_bool, "XRW TODO")
		+ XMLSchemaAttribute( "fullatom", xsct_rosetta_bool, "XRW TODO")
		+ XMLSchemaAttribute( "bbmove", xsct_rosetta_bool, "XRW TODO")
		+ XMLSchemaAttribute( "temp", xsct_real, "XRW TODO")
		+ XMLSchemaAttribute( "rms", xsct_real, "XRW TODO")
		+ XMLSchemaAttribute( "nminsteps", xsct_non_negative_integer, "XRW TODO")
		+ XMLSchemaAttribute( "overlap", xsct_non_negative_integer, "XRW TODO")
		+ XMLSchemaAttribute( "ncycles", xsct_non_negative_integer, "XRW TODO")
		+ XMLSchemaAttribute( "fragbias", xs_string, "XRW TODO")
		+ XMLSchemaAttribute::attribute_w_default( "recover_low", xsct_rosetta_bool, "XRW TODO", "true")
		+ XMLSchemaAttribute::attribute_w_default( "force_ss", xsct_char, "XRW TODO", "D");
	//take a break from this fun list to restrict the next attribute.

	//argument legality checker for "strategy":
	// AMW: could sanely do this as an enumeration instead
	std::string const strategy_possibles("user|uniform|auto|density|density_nbr|geometry|rama|bfactors|chainbreak");
	XMLSchemaRestriction strategy_type;
	strategy_type.name("strategy_type");
	strategy_type.base_type( xs_string );
	strategy_type.add_restriction( xsr_pattern, strategy_possibles + "(," + strategy_possibles + ")+" );
	xsd.add_top_level_element( strategy_type );

	attlist
		+ XMLSchemaAttribute( "strategy", "strategy_type", "fragment bias strategies.  This string is a comma-separated list; allowed values include 'user', 'uniform', 'auto', 'density', 'density_nbr', 'geometry', 'rama', 'bfactors', 'chainbreak' ")
		//now we return to the FUN LIST OF FUN
		+ XMLSchemaAttribute( "rsd_wdw_to_refine", xsct_non_negative_integer, "residue window to refine")
		+ XMLSchemaAttribute( "score_threshold", xsct_real, "XRW TODO")
		+ XMLSchemaAttribute( "cumulate_prob", xsct_rosetta_bool, "XRW TODO")
		+ XMLSchemaAttribute( "automode_scorecut", xsct_real, "XRW TODO")
		+ XMLSchemaAttribute( "dump_pdb", xsct_rosetta_bool, "XRW TODO")
		+ XMLSchemaAttribute( "dump_pdb_tag", xs_string, "XRW TODO")
		+ XMLSchemaAttribute( "wdw_to_freeze", xs_integer, "XRW TODO") //is stored as int
		+ XMLSchemaAttribute( "freeze_endpoints", xsct_rosetta_bool, "XRW TODO");

	std::string const residues_loops_warning(" 'residues' and 'loops_in' are mutually exclusive.");

	core::pose::attributes_for_get_resnum_selector( attlist, xsd, "residues"); //XRW TODO, fix function to take a docstring; here is partial docstring:
	//+ XMLSchemaAttribute( "residues", xsct_refpose_enabled_residue_number_cslist, "XRW TODO" + residues_loops_warning) // XRW TODO this is the wrong type, use get_attributes_for_get_resnum_selector instead

	attlist + XMLSchemaAttribute( "loops_in", xs_string, "XRW TODO.  Looks for a Loops object in the DataMap with this name." + residues_loops_warning);

	core::pose::attributes_for_get_resnum_selector( attlist, xsd, "residues_to_include");
	core::pose::attributes_for_get_resnum_selector( attlist, xsd, "residues_to_exclude");

	attlist
		+ XMLSchemaAttribute( "reference_model", xs_string, "XRW TODO; if 'input' uses the Pose from apply(); otherwise loads this pdb file")
		+ XMLSchemaAttribute( "reference_cst_wt", xsct_real, "XRW TODO")
		+ XMLSchemaAttribute( "reference_cst_maxdist", xsct_real, "XRW TODO");

	// attributes for Fragments subelement
	AttributeList fragment_subelement_attributes;
	fragment_subelement_attributes
		+ XMLSchemaAttribute::required_attribute( "fragfile", xs_string, "file to read fragments from")
		+ XMLSchemaAttribute( "nfrags", xsct_non_negative_integer, "read this many fragments from the top of 'fragfile'");

	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( & CartesianSampler_subelement_ct_name );
	subelements.add_simple_subelement( "Fragments", fragment_subelement_attributes, "Instructions for fragments files.");

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "nfrags", xsct_non_negative_integer, "XRW TODO; unclear if this should be used with the Fragments subelement", "25")
		+ XMLSchemaAttribute( "fraglens", xsct_positive_integer_cslist, "XRW TODO; comma separated list of positive integers, unclear if this should be used with the Fragments subelement");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TODO", attlist, subelements );
}

std::string CartesianSamplerCreator::keyname() const {
	return CartesianSampler::mover_name();
}

protocols::moves::MoverOP
CartesianSamplerCreator::create_mover() const {
	return protocols::moves::MoverOP( new CartesianSampler );
}

void CartesianSamplerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CartesianSampler::provide_xml_schema( xsd );
}


} // hybridization
}// protocols
