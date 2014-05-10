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

#include <protocols/hybridization/CartesianHybridize.hh>
#include <protocols/hybridization/CartesianHybridizeCreator.hh>
#include <protocols/hybridization/HybridizeSetup.hh>

#include <protocols/hybridization/TemplateHistory.hh>
#include <protocols/hybridization/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

//#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/util.hh>

#include <protocols/simple_moves/FragmentMover.hh>

#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/comparative_modeling/coord_util.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>

#include <core/scoring/rms_util.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Remarks.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/pose/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
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
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>


#include <utility/excn/Exceptions.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <basic/Tracer.hh>
#include <boost/unordered/unordered_map.hpp>

// parser
#include <protocols/rosetta_scripts/util.hh>
//#include <protocols/moves/DataMap.hh>
#include <utility/tag/Tag.hh>

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

const core::Size DEFAULT_NCYCLES=400;

static basic::Tracer TR("protocols.hybridization.CartesianHybridize");

CartesianHybridize::CartesianHybridize( ) :
	hybridize_setup_(NULL),
	ncycles_(DEFAULT_NCYCLES)
{
	init();
}

CartesianHybridize::CartesianHybridize(
		utility::vector1 < core::pose::PoseOP > const & templates_in,
		utility::vector1 < core::Real > const & template_wts_in,
		utility::vector1 < protocols::loops::Loops > const & template_chunks_in,
		utility::vector1 < protocols::loops::Loops > const & template_contigs_in,
		core::fragment::FragSetOP fragments9_in ) : ncycles_(DEFAULT_NCYCLES) {
	init();

	templates_ = templates_in;
	template_wts_ = template_wts_in;
	template_contigs_ = template_contigs_in;
	fragments9_ = fragments9_in;

	// make sure all data is there
	runtime_assert( templates_.size() == template_wts_.size() );
	runtime_assert( templates_.size() == template_contigs_.size() );

	// normalize weights
	core::Real weight_sum = 0.0;
	for (int i=1; i<=(int)templates_.size(); ++i) weight_sum += template_wts_[i];
	for (int i=1; i<=(int)templates_.size(); ++i) template_wts_[i] /= weight_sum;

	// map resids to frames
	for (core::fragment::ConstFrameIterator i = fragments9_->begin(); i != fragments9_->end(); ++i) {
		core::Size position = (*i)->start();
		library_[position] = **i;
	}

	// use chunks to subdivide contigs
	core::Size ntempls = templates_.size();
	for( core::Size tmpl = 1; tmpl <= ntempls; ++tmpl) {
		core::Size ncontigs = template_contigs_[tmpl].size();  // contigs to start
		for (int i=1; i<=(int)ncontigs; ++i) {
			core::Size cstart = template_contigs_[tmpl][i].start(), cstop = template_contigs_[tmpl][i].stop();
			bool spilt_chunk=false;

			// assumes sorted
			for (int j=2; j<=(int)template_chunks_in[tmpl].size(); ++j) {
				core::Size j0start = template_chunks_in[tmpl][j-1].start(), j0stop = template_chunks_in[tmpl][j-1].stop();
				core::Size j1start = template_chunks_in[tmpl][j].start(), j1stop = template_chunks_in[tmpl][j].stop();

				bool j0incontig = ((j0start>=cstart) && (j0stop<=cstop));
				bool j1incontig = ((j1start>=cstart) && (j1stop<=cstop));

				if (j0incontig && j1incontig) {
					spilt_chunk=true;
					core::Size cutpoint = (j0stop+j1start)/2;
					template_contigs_[tmpl].add_loop( cstart, cutpoint-1 );
					TR.Debug << "Make subfrag " << cstart << " , " << cutpoint-1 << std::endl;
					cstart=cutpoint;
				} else if(spilt_chunk && j0incontig && !j1incontig) {
					spilt_chunk=false;
					template_contigs_[tmpl].add_loop( cstart, cstop );
					TR.Debug << "Make subfrag " << cstart << " , " << cstop << std::endl;
				}
			}
		}
		template_contigs_[tmpl].sequential_order();
	}
	TR.Debug << "template_contigs:" << std::endl;
	for (int i=1; i<= (int)template_contigs_.size(); ++i) {
		TR.Debug << "templ. " << i << std::endl << template_contigs_[i] << std::endl;
	}
}

void CartesianHybridize::setup_for_parser()
{
	// RosettaScripts uses default constructor. Need to set up everything using datamap at the beginning of apply()
	templates_ = hybridize_setup_->template_poses();
	template_wts_ = hybridize_setup_->template_wts();
	template_contigs_ = hybridize_setup_->template_contigs();
	fragments9_ = hybridize_setup_->fragments_big()[numeric::random::random_range(1,hybridize_setup_->fragments_big().size())];
	
	// make sure all data is there
	runtime_assert( templates_.size() == template_wts_.size() );
	runtime_assert( templates_.size() == template_contigs_.size() );
	
	// normalize weights
	core::Real weight_sum = 0.0;
	for (int i=1; i<=(int)templates_.size(); ++i) weight_sum += template_wts_[i];
	for (int i=1; i<=(int)templates_.size(); ++i) template_wts_[i] /= weight_sum;
	
	// map resids to frames
	for (core::fragment::ConstFrameIterator i = fragments9_->begin(); i != fragments9_->end(); ++i) {
		core::Size position = (*i)->start();
		library_[position] = **i;
	}
	
	// use chunks to subdivide contigs
	core::Size ntempls = templates_.size();
	utility::vector1 < protocols::loops::Loops > template_chunks_in (hybridize_setup_->template_chunks());
	for( core::Size tmpl = 1; tmpl <= ntempls; ++tmpl) {
		core::Size ncontigs = template_contigs_[tmpl].size();  // contigs to start
		for (int i=1; i<=(int)ncontigs; ++i) {
			core::Size cstart = template_contigs_[tmpl][i].start(), cstop = template_contigs_[tmpl][i].stop();
			bool spilt_chunk=false;
			
			// assumes sorted
			for (int j=2; j<=(int)template_chunks_in[tmpl].size(); ++j) {
				core::Size j0start = template_chunks_in[tmpl][j-1].start(), j0stop = template_chunks_in[tmpl][j-1].stop();
				core::Size j1start = template_chunks_in[tmpl][j].start(), j1stop = template_chunks_in[tmpl][j].stop();
				
				bool j0incontig = ((j0start>=cstart) && (j0stop<=cstop));
				bool j1incontig = ((j1start>=cstart) && (j1stop<=cstop));
				
				if (j0incontig && j1incontig) {
					spilt_chunk=true;
					core::Size cutpoint = (j0stop+j1start)/2;
					template_contigs_[tmpl].add_loop( cstart, cutpoint-1 );
					TR.Debug << "Make subfrag " << cstart << " , " << cutpoint-1 << std::endl;
					cstart=cutpoint;
				} else if(spilt_chunk && j0incontig && !j1incontig) {
					spilt_chunk=false;
					template_contigs_[tmpl].add_loop( cstart, cstop );
					TR.Debug << "Make subfrag " << cstart << " , " << cstop << std::endl;
				}
			}
		}
		template_contigs_[tmpl].sequential_order();
	}
	TR.Debug << "template_contigs:" << std::endl;
	for (int i=1; i<= (int)template_contigs_.size(); ++i) {
		TR.Debug << "templ. " << i << std::endl << template_contigs_[i] << std::endl;
	}
}
	
void
CartesianHybridize::init() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	increase_cycles_ = option[cm::hybridize::stage2_increase_cycles]();
	no_global_frame_ = option[cm::hybridize::no_global_frame]();
	linmin_only_ = option[cm::hybridize::linmin_only]();

	// only adjustable via methods (for now)
	cartfrag_overlap_ = 2;
	align_templates_to_pose_ = false;
	
	seqfrags_only_ = false;
	nofragbias_ = false;
	skip_long_min_ = false;

	// default scorefunction
	set_scorefunction ( core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" ) );
}

void
CartesianHybridize::set_scorefunction(core::scoring::ScoreFunctionOP scorefxn_in) {
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
CartesianHybridize::apply_frag( core::pose::Pose &pose, core::pose::Pose &templ, protocols::loops::Loop &frag, bool superpose) {
	// superimpose frag
	numeric::xyzMatrix< core::Real > R;
	numeric::xyzVector< core::Real > preT(0,0,0), postT(0,0,0);
	R.xx() = R.yy() = R.zz() = 1;
	R.xy() = R.yx() = R.zx() = R.zy() = R.yz() = R.xz() = 0;

	if (superpose) {
		core::Size len = frag.size();
		core::Size aln_len = std::min( (core::Size)9999, len );   //fpd  can change 9999 to some max alignment sublength
		core::Size aln_start = numeric::random::random_range(frag.start(), len-aln_len+frag.start() );

		// don't try to align really short frags
		if (len > 2) {
			ObjexxFCL::FArray2D< core::Real > final_coords( 3, 4*aln_len );
			ObjexxFCL::FArray2D< core::Real > init_coords( 3, 4*aln_len );

			for (int ii=0; ii<(int)aln_len; ++ii) {
				int i=aln_start+ii;
				numeric::xyzVector< core::Real > x_1 = templ.residue(i).atom(" C  ").xyz();
				numeric::xyzVector< core::Real > x_2 = templ.residue(i).atom(" O  ").xyz();
				numeric::xyzVector< core::Real > x_3 = templ.residue(i).atom(" CA ").xyz();
				numeric::xyzVector< core::Real > x_4 = templ.residue(i).atom(" N  ").xyz();
				preT += x_1+x_2+x_3+x_4;

				numeric::xyzVector< core::Real > y_1 = pose.residue(templ.pdb_info()->number(i)).atom(" C  ").xyz();
				numeric::xyzVector< core::Real > y_2 = pose.residue(templ.pdb_info()->number(i)).atom(" O  ").xyz();
				numeric::xyzVector< core::Real > y_3 = pose.residue(templ.pdb_info()->number(i)).atom(" CA ").xyz();
				numeric::xyzVector< core::Real > y_4 = pose.residue(templ.pdb_info()->number(i)).atom(" N  ").xyz();
				postT += y_1+y_2+y_3+y_4;

				for (int j=0; j<3; ++j) {
					init_coords(j+1,4*ii+1) = x_1[j];
					init_coords(j+1,4*ii+2) = x_2[j];
					init_coords(j+1,4*ii+3) = x_3[j];
					init_coords(j+1,4*ii+4) = x_4[j];
					final_coords(j+1,4*ii+1) = y_1[j];
					final_coords(j+1,4*ii+2) = y_2[j];
					final_coords(j+1,4*ii+3) = y_3[j];
					final_coords(j+1,4*ii+4) = y_4[j];
				}
			}
			preT /= 4*len;
			postT /= 4*len;
			for (int i=1; i<=4*(int)len; ++i) {
				for ( int j=0; j<3; ++j ) {
					init_coords(j+1,i) -= preT[j];
					final_coords(j+1,i) -= postT[j];
				}
			}

			// get optimal superposition
			// rotate >init< to >final<
			ObjexxFCL::FArray1D< numeric::Real > ww( 4*len, 1.0 );
			ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
			numeric::Real ctx;

			numeric::model_quality::findUU( init_coords, final_coords, ww, 4*len, uu, ctx );
			R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
			R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
			R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );
		}
	}

	// xyz copy fragment to pose
	for (int i=(int)frag.start(); i<=(int)frag.stop(); ++i) {
		for (int j=1; j<=(int)templ.residue(i).natoms(); ++j) {
			core::id::AtomID src(j,i), tgt(j, templ.pdb_info()->number(i));
			pose.set_xyz( tgt, postT + (R*(templ.xyz( src )-preT)) );
		}
	}
}


void
CartesianHybridize::apply_frame( core::pose::Pose & pose, core::fragment::Frame &frame ) {
	core::Size start = frame.start(),len = frame.length();

	// we might want to tune this
	// it might make sense to change this based on gap width
	// for really large gaps make it one sided to emulate fold-tree fragment insertion
	int aln_len = cartfrag_overlap_;
	runtime_assert( cartfrag_overlap_>=1 && cartfrag_overlap_<=len/2);

	core::Size nres = pose.total_residue();

	//symmetry
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres = symm_info->num_independent_residues();
	}

	while (!pose.residue(nres).is_protein()) nres--;
	bool nterm = (start == 1);
	bool cterm = (start == nres-8);

	// insert frag
	core::pose::Pose pose_copy = pose;

	ObjexxFCL::FArray1D< numeric::Real > ww( 2*4*aln_len, 1.0 );
	ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
	numeric::xyzVector< core::Real > com1(0,0,0), com2(0,0,0);

	for (int i=0; i<(int)len; ++i) {
		core::conformation::idealize_position(start+i, pose_copy.conformation());
	}
	for (int tries = 0; tries<80; ++tries) {
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

		numeric::Real ctx;
		float rms;

		numeric::model_quality::findUU( final_coords, init_coords, ww, 2*4*aln_len, uu, ctx );
		numeric::model_quality::calc_rms_fast( rms, final_coords, init_coords, ww, 2*4*aln_len, ctx );

		//fpd  another place where we might want to tune parameters
		if (rms < 0.5) break;
		if (tries >= 20 && rms < 1) break;
		if (tries >= 40 && rms < 2) break;
		if (tries >= 60 && rms < 4) break;
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


void
CartesianHybridize::apply( Pose & pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose::datacache;

	if (hybridize_setup_) {
		if (align_templates_to_pose_) {
			TR << "Realigning template domains to the current pose." << std::endl;
			core::pose::PoseOP pose_copy = new core::pose::Pose( pose );
			hybridize_setup_->realign_templates(pose_copy);
		}
		
		setup_for_parser();
	}
    else {
        core::Size nres_nonvirt = get_num_residues_nonvirt(pose);
        residue_sample_template_.resize(nres_nonvirt, true);
        residue_sample_abinitio_.resize(nres_nonvirt, true);
    }

	//protocols::viewer::add_conformation_viewer(  pose.conformation(), "hybridize" );
	
	///////////////////////////////
	// added by yuan 06-28-2013
	// packer
	using namespace core::pack::task;
	simple_moves::PackRotamersMoverOP pack_rotamers;
	pack_rotamers = new protocols::simple_moves::PackRotamersMover();
	TaskFactoryOP main_task_factory = new TaskFactory;
	main_task_factory->push_back( new operation::RestrictToRepacking );
	pack_rotamers->task_factory(main_task_factory);
	pack_rotamers->score_function(lowres_scorefxn_);
	
	if (option[corrections::score::cenrot]) {
		protocols::moves::MoverOP tocenrot =
			new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID_ROT );
		tocenrot->apply( pose );
		pack_rotamers->apply(pose);
	}
	else {
		protocols::moves::MoverOP tocen =
			new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );
		tocen->apply( pose );
	}

	// minimizer
	core::optimization::MinimizerOptions options( "linmin", 0.01, true, false, false );
	core::optimization::MinimizerOptions options_minilbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
	options_minilbfgs.max_iter(5);
	core::optimization::MinimizerOptions options_lbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
	if (increase_cycles_ < 1.) {
		Size niter = (Size) (200*increase_cycles_);
		options_lbfgs.max_iter(niter);
	}	else {
		options_lbfgs.max_iter(200);
	}
	core::optimization::CartesianMinimizer minimizer;
	core::kinematics::MoveMap mm;

	mm.set_bb  ( true );
	mm.set_chi ( true );
	mm.set_jump( true );

	//fpd  --  this should really be automated somewhere
	if (core::pose::symmetry::is_symmetric(pose) ) {
		core::pose::symmetry::make_symmetric_movemap( pose, mm );
	}

	core::Real max_cart = lowres_scorefxn_->get_weight( core::scoring::cart_bonded );
	core::Real max_cart_angle = lowres_scorefxn_->get_weight( core::scoring::cart_bonded_angle );
	core::Real max_cart_length = lowres_scorefxn_->get_weight( core::scoring::cart_bonded_length );
	core::Real max_cart_torsion = lowres_scorefxn_->get_weight( core::scoring::cart_bonded_torsion );
	core::Real max_cst  = lowres_scorefxn_->get_weight( core::scoring::atom_pair_constraint );
	core::Real max_vdw  = lowres_scorefxn_->get_weight( core::scoring::vdw );

	// for i = 1 to n cycles
	core::Size NMACROCYCLES = 4;
	TR << "RUNNING FOR " << NMACROCYCLES << " MACROCYCLES" << std::endl;

	Pose pose_in = pose;

	core::Size nres = pose.total_residue();
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres = symm_info->num_independent_residues();
	}
	while (pose.residue(nres).aa() == core::chemical::aa_vrt) nres--;

	core::Size n_prot_res = pose.total_residue();
	while (!pose.residue(n_prot_res).is_protein()) n_prot_res--;

	// 10% of the time, skip moves in the global frame
	bool no_ns_moves = no_global_frame_; // (numeric::random::uniform() <= 0.1);

sampler:
	for (core::Size m=1; m<=NMACROCYCLES; ++m) {
		core::Real bonded_weight = max_cart;
		if (m==1) bonded_weight = 0.0*max_cart;
		if (m==2) bonded_weight = 0.01*max_cart;
		if (m==3) bonded_weight = 0.1*max_cart;

		core::Real bonded_weight_angle = max_cart_angle;
		if (m==1) bonded_weight_angle = 0.0*max_cart_angle;
		if (m==2) bonded_weight_angle = 0.01*max_cart_angle;
		if (m==3) bonded_weight_angle = 0.1*max_cart_angle;

		core::Real bonded_weight_length = max_cart_length;
		if (m==1) bonded_weight_length = 0.0*max_cart_length;
		if (m==2) bonded_weight_length = 0.01*max_cart_length;
		if (m==3) bonded_weight_length = 0.1*max_cart_length;

		core::Real bonded_weight_torsion = max_cart_torsion;
		if (m==1) bonded_weight_torsion = 0.0*max_cart_torsion;
		if (m==2) bonded_weight_torsion = 0.01*max_cart_torsion;
		if (m==3) bonded_weight_torsion = 0.1*max_cart_torsion;

		core::Real cst_weight = max_cst;
		if (m==1)  cst_weight = 2*max_cst;
		if (m==2)  cst_weight = 2*max_cst;
		if (m==3)  cst_weight = 2*max_cst;

		core::Real vdw_weight = max_vdw;
		if (m==1) vdw_weight = 0.1*max_vdw;
		if (m==2) vdw_weight = 0.1*max_vdw;
		if (m==3) vdw_weight = 0.1*max_vdw;

		TR << "CYCLE " << m << std::endl;
		TR << "  setting bonded weight = " << bonded_weight << std::endl;
		TR << "  setting bonded angle weight = " << bonded_weight_angle << std::endl;
		TR << "  setting bonded length weight = " << bonded_weight_length << std::endl;
		TR << "  setting bonded torsion weight = " << bonded_weight_torsion << std::endl;
		TR << "  setting cst    weight = " << cst_weight << std::endl;
		TR << "  setting vdw    weight = " << vdw_weight << std::endl;

		lowres_scorefxn_->set_weight( core::scoring::cart_bonded, bonded_weight );
		lowres_scorefxn_->set_weight( core::scoring::cart_bonded_angle, bonded_weight_angle );
		lowres_scorefxn_->set_weight( core::scoring::cart_bonded_length, bonded_weight_length );
		lowres_scorefxn_->set_weight( core::scoring::cart_bonded_torsion, bonded_weight_torsion );
		lowres_scorefxn_->set_weight( core::scoring::atom_pair_constraint, cst_weight );
		lowres_scorefxn_->set_weight( core::scoring::vdw, vdw_weight );

		(*lowres_scorefxn_)(pose);
		protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( pose, *lowres_scorefxn_, 
			option[cm::hybridize::stage2_temperature]() ); //cenrot may use higher temp

		core::Size neffcycles = (core::Size)(ncycles_*increase_cycles_);
		if (m==4 || seqfrags_only_) neffcycles /=2;
		for (int n=1; n<=(int)neffcycles; ++n) {
			// possible actions:
			//  1 - insert homologue frag, global frame
			//  2 - insert homologue frag, local frame
			//  3 - insert sequence frag
			core::Real action_picker = numeric::random::uniform();
			core::Size action = 0;
			if (m==1) {
				action = no_ns_moves?2:1;
				if (action_picker < 0.2) action = 3;
			} else if (m==2) {
				action = 2;
				if (action_picker < 0.2) action = 3;
			} else if (m==3) {
				action = 2;
				if (action_picker < 0.2) action = 3;
			} else if (m==4) {
				action = 3;
			}

			if ( seqfrags_only_ ) action = 3;

			std::string action_string;
			if (action == 1) action_string = "fragNS";
			if (action == 2) action_string = "frag";
			if (action == 3) action_string = "picker";

			if (action == 1 || action == 2) {
				core::Size max_templates_trial=templates_.size()*5;
				bool movable_loop=false;
				protocols::loops::LoopOP frag;
                
				core::Size templ_id=1;
				for (core::Size i=1; i<=max_templates_trial; ++i) {
                    templ_id = numeric::random::random_range( 1, templates_.size() );
                    core::Size nfrags = template_contigs_[templ_id].num_loop();
                    // remove non-protein frags
                    while (templates_[templ_id]->pdb_info()->number(template_contigs_[templ_id][nfrags].start()) >
                           (int)n_prot_res) {
                        --nfrags;
                    }
                    
                    core::Size frag_id;
                    core::Size max_frag_trial=nfrags*3;
                    for (core::Size ii=1; ii<=max_frag_trial; ++ii) {
                        frag_id = numeric::random::random_range( 1, nfrags );
                        frag =  new protocols::loops::Loop ( template_contigs_[templ_id][frag_id] );
                        for (core::Size iii=frag->start(); iii<=frag->stop(); ++iii) {
                            if ( residue_sample_template_[templates_[templ_id]->pdb_info()->number(iii)]==true ) {
                                movable_loop=true;
                                break;
                            }
                        }
                        if ( movable_loop==true)
                            break;
                    } //trial different frags in a tempalte
				    if ( movable_loop==true)
                        break;
                } //end of trial different templates
                
				if (frag->size() > 14)
					action_string = action_string+"_15+";
				else if (frag->size() <= 4)
					action_string = action_string+"_0-4";
				else
					action_string = action_string+"_5-14";
                
				if ( frag->size() > 0 )
					apply_frag( pose, *templates_[templ_id], *frag, (action==2) );
                
				if (action == 1) {
					//fpd assume this was initialized elsewhere
					runtime_assert( pose.data().has( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY ) );
					TemplateHistory &history =
                    *( static_cast< TemplateHistory* >( pose.data().get_ptr( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY )() ));
					history.set( frag->start(), frag->stop(), templ_id );
				}
			}

			if (action == 3) {
				// pick an insert position based on gap
				utility::vector1<core::Real> residuals( n_prot_res , 0.0 );
				utility::vector1<core::Real> max_residuals(3,0);
				utility::vector1<int> max_poses(4,-1);
				if (!nofragbias_) {
					for (int i=1; i<(int)n_prot_res; ++i) {
						if (!pose.residue_type(i).is_protein()){
							residuals[i] = -1;
						} else if (pose.fold_tree().is_cutpoint(i+1)) {
							residuals[i] = -1;
						} else if ( i > static_cast <int> (residue_sample_abinitio_.size()) ) {
							residuals[i] = -1;
                        } else if ( residue_sample_abinitio_[i]==false ) {
							residuals[i] = -1;
						} else {
							numeric::xyzVector< core::Real > c0 , n1;
							c0 = pose.residue(i).atom(" C  ").xyz();
							n1 = pose.residue(i+1).atom(" N  ").xyz();
							core::Real d2 = c0.distance( n1 );
							residuals[i] = (d2-1.328685)*(d2-1.328685);
							if ( residuals[i] > max_residuals[1]) {
								max_residuals[3] = max_residuals[2]; max_residuals[2] = max_residuals[1]; max_residuals[1] = residuals[i];
								max_poses[3] = max_poses[2]; max_poses[2] = max_poses[1]; max_poses[1] = i;
							} else if ( residuals[i] > max_residuals[2]) {
								max_residuals[3] = max_residuals[2]; max_residuals[2] = residuals[i];
								max_poses[3] = max_poses[2]; max_poses[2] = i;
							} else if ( residuals[i] > max_residuals[3]) {
								max_residuals[3] = residuals[i];
								max_poses[3] = i;
							}
						}
					}
				} else {
					max_poses[3] = max_poses[2] = max_poses[1] = 0;
				}

				// 25% chance of random position
				int random_residue_move=numeric::random::random_range(1, residue_sample_abinitio_.size());
                int ntrials=500;
				while (!residue_sample_abinitio_[random_residue_move] && --ntrials>0) {
                    random_residue_move=numeric::random::random_range(1, residue_sample_abinitio_.size());
				}
                if (ntrials<=0) {
                    TR << "Warning! Fail to find a free residue for sampling." << std::endl;
                    continue;
                }
				max_poses[ 4 ] = random_residue_move;
				int select_position = numeric::random::random_range(1,4);

				//fpd  hack
				if (nofragbias_) select_position=4;

				if (select_position == 4)
					action_string = action_string+"_rand";
				core::Size max_pos = max_poses[ select_position ];

				// select random pos in [i-8,i]
				core::Size insert_pos = max_pos - numeric::random::random_range(3,5);
				insert_pos = std::min( insert_pos, n_prot_res-8);
				insert_pos = std::max( (int)insert_pos, 1);

				if (library_.find(insert_pos) != library_.end())
					apply_frame (pose, library_[insert_pos]);
			}

			//////////////////////////////
			// added by yuan 06-28-2013
			// repack after fragment insertion
			//////////////////////////////
			if (option[corrections::score::cenrot]) {
				pack_rotamers->apply(pose);
			}

			// MC
			try {
				(*min_scorefxn_)(pose);
				if ( m<4 || linmin_only_ )
					minimizer.run( pose, mm, *min_scorefxn_, options );
				else
					minimizer.run( pose, mm, *min_scorefxn_, options_minilbfgs );

				mc->boltzmann( pose , action_string );
			} catch( utility::excn::EXCN_Base& excn ) {
				//fpd hbond fail? start over
				pose = pose_in;
				goto sampler;
			}

			if (n%100 == 0) {
				mc->show_scores();
				mc->show_counters();
			}

		}
		mc->recover_low(pose);
	}

	//for (int i=1; i<=templates_.size(); ++i) {
	//	std::ostringstream oss;
	//	oss << "templ"<<i<<".pdb";
	//	templates_[i]->dump_pdb( oss.str() );
	//}

	// final minimization
	try {
		if (!skip_long_min_) {
				(*min_scorefxn_)(pose); minimizer.run( pose, mm, *min_scorefxn_, options_lbfgs );
				(*nocst_scorefxn_)(pose); minimizer.run( pose, mm, *nocst_scorefxn_, options_lbfgs );
				(*bonds_scorefxn_)(pose); minimizer.run( pose, mm, *bonds_scorefxn_, options_lbfgs );
				(*nocst_scorefxn_)(pose); minimizer.run( pose, mm, *nocst_scorefxn_, options_lbfgs );
		}
	} catch( utility::excn::EXCN_Base& excn ) {
		//fpd hbond fail? start over
		pose = pose_in;
		goto sampler;
	}

	lowres_scorefxn_->set_weight( core::scoring::cart_bonded, max_cart );
	lowres_scorefxn_->set_weight( core::scoring::cart_bonded_angle, max_cart_angle );
	lowres_scorefxn_->set_weight( core::scoring::cart_bonded_length, max_cart_length );
	lowres_scorefxn_->set_weight( core::scoring::cart_bonded_torsion, max_cart_torsion );
	lowres_scorefxn_->set_weight( core::scoring::atom_pair_constraint, max_cst );
	lowres_scorefxn_->set_weight( core::scoring::vdw, max_vdw );
	(*lowres_scorefxn_)(pose);
}

protocols::moves::MoverOP CartesianHybridize::clone() const { return new CartesianHybridize( *this ); }
protocols::moves::MoverOP CartesianHybridize::fresh_instance() const { return new CartesianHybridize; }

void
CartesianHybridize::parse_my_tag(
								utility::tag::TagCOP tag,
								basic::datacache::DataMap & data,
								filters::Filters_map const &,
								moves::Movers_map const &,
								core::pose::Pose const & pose)
{
	if( tag->hasOption( "hybridize_setup" ) ) {
		std::string const setup_data_name( tag->getOption<std::string>( "hybridize_setup" ) );
		hybridize_setup_ = (data.get< protocols::hybridization::HybridizeSetup * >( "HybridizeSetup", setup_data_name ))->clone();
	}
	else {
		utility_exit_with_message("Fatal error: hybridize_setup needs to be defined!");
	}

	if( tag->hasOption( "increase_cycles" ) ) {
		set_increase_cycles( tag->getOption< core::Real >( "increase_cycles" ) );
	}
	
	if( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn" ) );
		set_scorefunction( (data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ))->clone() );
	}

	//task operations
    /*
    allowed_to_move_.clear();
	allowed_to_move_.resize(pose.total_residue(),true);
	if( tag->hasOption( "task_operations" ) ){
		core::pack::task::TaskFactoryOP task_factory = protocols::rosetta_scripts::parse_task_operations( tag, data );
		core::pack::task::PackerTaskOP task = task_factory->create_task_and_apply_taskoperations( pose );
		for( core::Size resi = 1; resi <= get_num_residues_prot(pose); ++resi ){
			if( task->residue_task( resi ).being_designed() || task->residue_task( resi ).being_packed())
				allowed_to_move_[resi]=true;
    		else
				allowed_to_move_[resi]=false;
		}
	}
	*/
    
	if( tag->hasOption( "no_global_frame" ) )
		set_no_global_frame( tag->getOption< bool >( "no_global_frame" ) );
	if( tag->hasOption( "linmin_only" ) )
		set_linmin_only( tag->getOption< bool >( "linmin_only" ) );
	if( tag->hasOption( "cartfrag_overlap" ) )
		set_cartfrag_overlap( tag->getOption< core::Size >( "cartfrag_overlap" ) );
	if( tag->hasOption( "align_templates_to_pose" ) ) {
		align_templates_to_pose_ = tag->getOption< bool >( "align_templates_to_pose" );
	}
    
    residue_sample_template_.resize(hybridize_setup_->nres_tgt_asu(), true);
    residue_sample_abinitio_.resize(hybridize_setup_->nres_tgt_asu(), true);
    
    utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for (tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it) {
        // per-residue control
        if ( (*tag_it)->getName() == "DetailedControls" ) {
            if( (*tag_it)->hasOption( "task_operations" ) ){
                core::pack::task::TaskFactoryOP task_factory = protocols::rosetta_scripts::parse_task_operations( *tag_it, data );
                core::pack::task::PackerTaskOP task = task_factory->create_task_and_apply_taskoperations( pose );
                
                for( core::Size ires = 1; ires <= hybridize_setup_->nres_tgt_asu(); ++ires ){
                    
                    if( task->residue_task( ires ).being_designed() && task->residue_task( ires ).being_packed() ) {
                        // residue_cst_cross_chain_[ires] = true;
                    }
                    else {
                        residue_sample_template_[ires] = false;
                        residue_sample_abinitio_[ires] = false;
                    }
                }
            }
            else {
                core::Size start_res = (*tag_it)->getOption<core::Size>( "start_res", 1 );
                core::Size stop_res = (*tag_it)->getOption<core::Size>( "stop_res", hybridize_setup_->nres_tgt_asu() );
                if( (*tag_it)->hasOption( "sample_template" ) ){
                    bool sample_template = (*tag_it)->getOption<bool>( "sample_template", true );
                    for (core::Size ires=start_res; ires<=stop_res; ++ires) {
                        residue_sample_template_[ires] = sample_template;
                    }
                }
                if( (*tag_it)->hasOption( "sample_abinitio" ) ){
                    bool sample_abinitio = (*tag_it)->getOption<bool>( "sample_abinitio", true );
                    for (core::Size ires=start_res; ires<=stop_res; ++ires) {
                        residue_sample_abinitio_[ires] = sample_abinitio;
                    }
                }
            }
        }
        
    }
}
	
/////////////
// creator
std::string
CartesianHybridizeCreator::keyname() const {
	return CartesianHybridizeCreator::mover_name();
}

protocols::moves::MoverOP
CartesianHybridizeCreator::create_mover() const {
	return new CartesianHybridize;
}

std::string
CartesianHybridizeCreator::mover_name() {
	return "CartesianHybridize";
}

}
}
