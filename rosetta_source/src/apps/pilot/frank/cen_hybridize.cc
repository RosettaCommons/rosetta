/// @file
/// @brief

#include <devel/init.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <protocols/moves/PackRotamersMover.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/SwitchResidueTypeSetMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/rbsegment_moves/util.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/ConstraintSetMover.hh>
#include <protocols/jumping/Dssp.hh>
#include <protocols/comparative_modeling/coord_util.hh>

#include <protocols/evaluation/Align_RmsdEvaluator.hh>
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

#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

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

#include <boost/unordered/unordered_map.hpp>

OPT_1GRP_KEY(FileVector, fpd, templates)
OPT_1GRP_KEY(FileVector, fpd, fragments)
OPT_1GRP_KEY(File, fpd, frag9)
OPT_1GRP_KEY(Integer, fpd, ncycles)
OPT_1GRP_KEY(Integer, fpd, nmacrocycles)
OPT_1GRP_KEY(Integer, fpd, subfraglen)
OPT_1GRP_KEY(Integer, fpd, minfraglen)
OPT_1GRP_KEY(Boolean, fpd, movie)


///////////////////////////////////////////////////////////////////////////////

class CustomMover : public protocols::moves::Mover {
public:
	CustomMover( ) {
		using namespace protocols::moves;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// load all templates for insertion
		utility::vector1< utility::file::FileName > files = option[ fpd::templates ]();
		core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
		for (int i=1; i<= files.size(); ++i) {
			core::pose::PoseOP pose_i = new core::pose::Pose();
			core::import_pose::pose_from_pdb( *pose_i, *residue_set, files[i].name() );
			pose1s_.push_back( pose_i );
		}

		// load all the fragments for insertion
		files = option[ fpd::fragments ]();
		int MINFRAGLEN = option[ fpd::minfraglen ]();
		for (int i=1; i<= files.size(); ++i) {
			core::pose::PoseOP pose_i = new core::pose::Pose();
			core::import_pose::pose_from_pdb( *pose_i, *residue_set, files[i].name() );
			if (pose_i->total_residue() >= MINFRAGLEN) {
				pose2s_.push_back( pose_i );
				std::cout << "PDB " << files[i].name() << " insert point " << pose_i->pdb_info()->number(1) << std::endl;
			}
		}

		// default chunking on remaining long fragments
		// int FRAGLEN = option[ fpd::subfraglen ]();
		// core::Size nfrags = pose1s_.size();
		// for (int i=1; i<=nfrags; ++i) {
		// 	core::pose::PoseOP pose_i = pose1s_[i];
		// 	if (pose_i->total_residue() > FRAGLEN+4) {
		// 		int nres=pose_i->total_residue();
		// 		for (int j=1; j<=pose_i->total_residue()-FRAGLEN; j+=1) {
		// 			core::pose::PoseOP pose_j = new core::pose::Pose( *pose_i, j, j+FRAGLEN );
		// 			pose2s_.push_back( pose_j );
		// 			pose2_insert_points_.push_back( pose1_insert_points_[i]+j-1 );
		// 			std::cout << "Created " << pose_j->total_residue() << "-residue subfragment"
		// 					  << " from PDB " << files[i].name() 
		// 					  << " insert point " << pose_i->pdb_info()->number(j) << std::endl;
		// 		}
		// 	}
		// }

		// screfunction init
		lowres_scorefxn_ = core::scoring::getScoreFunction();

		min1_scorefxn_ = new core::scoring::ScoreFunction();
		min1_scorefxn_->set_weight( core::scoring::vdw, lowres_scorefxn_->get_weight( core::scoring::vdw ) );
		min1_scorefxn_->set_weight( core::scoring::cart_bonded, lowres_scorefxn_->get_weight( core::scoring::cart_bonded ) );

		min2_scorefxn_ = new core::scoring::ScoreFunction();
		min2_scorefxn_->set_weight( core::scoring::vdw, lowres_scorefxn_->get_weight( core::scoring::vdw ) );
		min2_scorefxn_->set_weight( core::scoring::hbond_lr_bb, lowres_scorefxn_->get_weight( core::scoring::hbond_lr_bb ) );
		min2_scorefxn_->set_weight( core::scoring::hbond_sr_bb, lowres_scorefxn_->get_weight( core::scoring::hbond_sr_bb ) );
		min2_scorefxn_->set_weight( core::scoring::cart_bonded, lowres_scorefxn_->get_weight( core::scoring::cart_bonded ) );
		min2_scorefxn_->set_weight( core::scoring::rama, lowres_scorefxn_->get_weight( core::scoring::rama ) );
		min2_scorefxn_->set_weight( core::scoring::omega, lowres_scorefxn_->get_weight( core::scoring::omega ) );

		min3_scorefxn_ = new core::scoring::ScoreFunction();
		min3_scorefxn_->set_weight( core::scoring::vdw, lowres_scorefxn_->get_weight( core::scoring::vdw ) );
		min3_scorefxn_->set_weight( core::scoring::hbond_lr_bb, lowres_scorefxn_->get_weight( core::scoring::hbond_lr_bb ) );
		min3_scorefxn_->set_weight( core::scoring::hbond_sr_bb, lowres_scorefxn_->get_weight( core::scoring::hbond_sr_bb ) );
		min3_scorefxn_->set_weight( core::scoring::cart_bonded, lowres_scorefxn_->get_weight( core::scoring::cart_bonded ) );
		min3_scorefxn_->set_weight( core::scoring::rama, lowres_scorefxn_->get_weight( core::scoring::rama ) );
		min3_scorefxn_->set_weight( core::scoring::omega, lowres_scorefxn_->get_weight( core::scoring::omega ) );
		min3_scorefxn_->set_weight( core::scoring::atom_pair_constraint, lowres_scorefxn_->get_weight( core::scoring::atom_pair_constraint ) );

		// change VDW set
		core::scoring::methods::EnergyMethodOptions lowres_options(lowres_scorefxn_->energy_method_options());
		lowres_options.atom_vdw_atom_type_set_name("centroid_min");
		min1_scorefxn_->set_energy_method_options(lowres_options);
		min2_scorefxn_->set_energy_method_options(lowres_options);
		min3_scorefxn_->set_energy_method_options(lowres_options);

		if ( option[ OptionKeys::fpd::frag9 ].user() ) {
			using namespace core::fragment;
			fragments_ = new ConstantLengthFragSet( 9 );
			fragments_ = FragmentIO().read_data( option[ OptionKeys::fpd::frag9 ]().name() );

			// code shamelessly stolen from nonlocal/SingleFragmentMover.cc
			for (core::fragment::FrameIterator i = fragments_->begin(); i != fragments_->end(); ++i) {
				core::Size position = (*i)->start();
				library_[position] = **i;
			}
		}

		// native
		if ( option[ in::file::native ].user() ) {
			native_ = new core::pose::Pose;
			core::import_pose::pose_from_pdb( *native_, option[ in::file::native ]() );
		}
	}

	void superpose( core::pose::Pose &frag, core::pose::Pose &pose, 
	                numeric::xyzMatrix< core::Real > &R, numeric::xyzVector< core::Real > &preT, numeric::xyzVector< core::Real > &postT) {
		// com of both
		core::Size len = frag.total_residue();
		core::Size aln_len = std::min( (core::Size)9, len );
		core::Size aln_start = numeric::random::random_range(1, len-aln_len+1 );

		preT = postT = numeric::xyzVector< core::Real >(0,0,0);
		if (len <= 2) {
			R.xx() = R.yy() = R.zz() = 1;
			R.xy() = R.yx() = R.zx() = R.zy() = R.yz() = R.xz() = 0;
			return;
		}

		ObjexxFCL::FArray2D< core::Real > final_coords( 3, 4*len );
		ObjexxFCL::FArray2D< core::Real > init_coords( 3, 4*len );
		preT = postT = numeric::xyzVector< core::Real >(0,0,0);
		for (int ii=1; ii<=(int)aln_len; ++ii) {
			int i=aln_start+ii-1;
			numeric::xyzVector< core::Real > x_1 = frag.residue(i).atom(" C  ").xyz();
			numeric::xyzVector< core::Real > x_2 = frag.residue(i).atom(" O  ").xyz();
			numeric::xyzVector< core::Real > x_3 = frag.residue(i).atom(" CA ").xyz();
			numeric::xyzVector< core::Real > x_4 = frag.residue(i).atom(" N  ").xyz();
			preT += x_1+x_2+x_3+x_4;
			numeric::xyzVector< core::Real > y_1 = pose.residue(frag.pdb_info()->number(i)).atom(" C  ").xyz();
			numeric::xyzVector< core::Real > y_2 = pose.residue(frag.pdb_info()->number(i)).atom(" O  ").xyz();
			numeric::xyzVector< core::Real > y_3 = pose.residue(frag.pdb_info()->number(i)).atom(" CA ").xyz();
			numeric::xyzVector< core::Real > y_4 = pose.residue(frag.pdb_info()->number(i)).atom(" N  ").xyz();
			postT += x_1+x_2+x_3+x_4;
			for (int j=0; j<3; ++j) { 
				init_coords(j+1,4*(i-1)+1) = x_1[j];
				init_coords(j+1,4*(i-1)+2) = x_2[j];
				init_coords(j+1,4*(i-1)+3) = x_3[j];
				init_coords(j+1,4*(i-1)+4) = x_4[j];
				final_coords(j+1,4*(i-1)+1) = y_1[j];
				final_coords(j+1,4*(i-1)+2) = y_2[j];
				final_coords(j+1,4*(i-1)+3) = y_3[j];
				final_coords(j+1,4*(i-1)+4) = y_4[j];
			}
		}
		preT /= 4*len;
		postT /= 4*len;
		for (int i=1; i<=(int)4*len; ++i) {
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

	///////////////////////////
	///////////////////////////
	void apply_frame( core::pose::Pose & pose, core::fragment::Frame &frame ) {
		core::Size start = frame.start(),len = frame.length();

		int aln_len = 4;
	
		bool nterm = (start == 1);
		bool cterm = (start == pose.total_residue()-8);

		// insert frag
		core::pose::Pose pose_copy = pose;

		ObjexxFCL::FArray1D< numeric::Real > ww( 2*4*aln_len, 1.0 );
		ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
		numeric::xyzVector< core::Real > com1(0,0,0), com2(0,0,0);

		for (int i=0; i<(int)len; ++i) {
			core::conformation::idealize_position(start+i, pose_copy.conformation());
		}
		for (int tries = 0; tries<100; ++tries) {
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
	
			// get optimal superposition
			// rotate >final< to >init<
			numeric::Real ctx;
			float rms;

			numeric::model_quality::findUU( final_coords, init_coords, ww, 2*4*aln_len, uu, ctx );
			numeric::model_quality::calc_rms_fast( rms, final_coords, init_coords, ww, 2*4*aln_len, ctx );

			//std::cout << "try " << tries << " rms " << rms << std::endl;

			if (rms < 0.5) break;
			if (tries >= 20 && rms < 1) break;
			if (tries >= 40 && rms < 2) break;
			if (tries >= 60 && rms < 3) break;
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

	///////////////////////////
	///////////////////////////
	void apply( Pose & pose ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// minimizer
		core::optimization::MinimizerOptions options( "linmin", 0.01, true, false, false );
		core::optimization::MinimizerOptions options_minilbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
		options_minilbfgs.max_iter(10);
		core::optimization::MinimizerOptions options_lbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
		options_lbfgs.max_iter(200);
		core::optimization::CartesianMinimizer minimizer;
		core::kinematics::MoveMap mm;
		mm.set_bb  ( true ); mm.set_chi ( true ); mm.set_jump( true );

		protocols::moves::MoverOP tocen = new protocols::moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );
		tocen->apply( pose );

		core::Real max_cart = lowres_scorefxn_->get_weight( core::scoring::cart_bonded );
		core::Real max_cst  = lowres_scorefxn_->get_weight( core::scoring::atom_pair_constraint );
		core::Real max_vdw  = lowres_scorefxn_->get_weight( core::scoring::vdw );

		// for i = 1 to n cycles
		core::Size ncycles = option[ fpd::ncycles ]();
		core::Size nmacrocycles = option[ fpd::nmacrocycles ]();
		std::cout << "RUNNING FOR " << nmacrocycles << " MACROCYCLES" << std::endl;
		
		if (nmacrocycles > 0) {
			for (int m=1; m<=std::min((int)nmacrocycles,5); m+=1) {
				core::Real bonded_weight = 0;
				if (m==1) bonded_weight = 0;
				if (m==2) bonded_weight = 0;
				if (m==3) bonded_weight = 0.01*max_cart;
				if (m==4) bonded_weight = 0.1*max_cart;
				if (m==5) bonded_weight = max_cart;
	
				core::Real cst_weight = max_cst;
				if (m<5)  cst_weight = 2*max_cst;

				core::Real vdw_weight = max_vdw;
				if (m==1)  vdw_weight = max_vdw;
				if (m==2)  vdw_weight = 0.1*max_vdw;
				if (m==3)  vdw_weight = 0.1*max_vdw;
				if (m==4)  vdw_weight = 0.1*max_vdw;
	
				std::cout << "CYCLE " << m << std::endl;
				std::cout << "  setting bonded weight = " << bonded_weight << std::endl;
				std::cout << "  setting cst    weight = " << cst_weight << std::endl;
				std::cout << "  setting vdw    weight = " << vdw_weight << std::endl;
				lowres_scorefxn_->set_weight( core::scoring::cart_bonded, bonded_weight );
				lowres_scorefxn_->set_weight( core::scoring::atom_pair_constraint, cst_weight );
				lowres_scorefxn_->set_weight( core::scoring::vdw, vdw_weight );
	
				(*lowres_scorefxn_)(pose);
				protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( pose, *lowres_scorefxn_, 2.0 );
				protocols::moves::MonteCarloOP mc_inner
					= new protocols::moves::MonteCarlo( pose, *lowres_scorefxn_, 0.0 );
	
				core::Size neffcycles = option[ fpd::ncycles ]();

				if (m==1) neffcycles = 5;

				for (int n=1; n<=neffcycles; ++n) {
					// superimpose frag
					numeric::xyzMatrix< core::Real > R;
					numeric::xyzVector< core::Real > preT(0,0,0), postT(0,0,0);
					R.xx() = R.yy() = R.zz() = 1;
					R.xy() = R.yx() = R.zx() = R.zy() = R.yz() = R.xz() = 0;
	
					// 1 - insert homologue frag, don't superpose
					// 2 - insert homologue subfrag, don't superpose
					// 3 - insert homologue frag, superpose
					// 4 - insert homologue subfrag, superpose
					// 5 - insert sequence frag
					core::Real action_picker = numeric::random::uniform();
					core::Size action = 0;
	
					if (m==1) {
						//if (n%STAGE1STEP == 1)
							action = 1;
						//else
						//	if (action_picker < 0.5) action = 2;
						//	else action = 4;
					} else if (m==2) {
						action = 4;
						if (action_picker < 0.2) action = 5;
					} else if (m==3) {
						action = 4;
						if (action_picker < 0.2) action = 5;
					} else if (m==4) {
						action = 4;
						if (action_picker < 0.2) action = 5;
					} else if (m==5) {
						action = 5;
					}

					std::string action_string;
					if (action == 1) action_string = "fragNS";
					if (action == 2) action_string = "subfragNS";
					if (action == 3) action_string = "frag";
					if (action == 4) action_string = "subfrag";
					if (action == 5) action_string = "picker";
	
					if (action == 5) {
						// pick an insert position
						utility::vector1<core::Real> residuals( pose.total_residue() , 0.0 );
						utility::vector1<core::Real> max_residuals(3,0);
						utility::vector1<int> max_poses(4,-1);
						for (int i=1; i<pose.total_residue(); ++i) {
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
	
						// 25% chance of random position
						max_poses[ 4 ] = numeric::random::random_range(1,pose.total_residue());
						int select_position = numeric::random::random_range(1,4);
						if (select_position == 4)
							action_string = action_string+"_rand";
						core::Size max_pos = max_poses[ select_position ];
	
						// select random pos in [i-8,i]
						core::Size insert_pos = max_pos - numeric::random::random_range(3,5);
						insert_pos = std::min( insert_pos, pose.total_residue()-8);
						insert_pos = std::max( (int)insert_pos, 1);
	
						//core::Size insert_pos = numeric::random::random_range(1, pose.total_residue()-8);
						//std::cerr << "type 3 insert at pos " << insert_pos << std::endl;
						apply_frame (pose, library_[insert_pos]);
	
					} else {
						// pick a fragment at random
						core::Size frag_set = (action == 2 || action == 4) ? 2 : 1;
						core::Size frag_id = numeric::random::random_range(1, frag_set==1 ? pose1s_.size():pose2s_.size() );
						core::pose::PoseOP frag = (frag_set==1) ? pose1s_[frag_id] : pose2s_[frag_id];
						//int insert_pt = (frag_set==1) ? pose1_insert_points_[frag_id] : pose2_insert_points_[frag_id];
	
						if (frag_set == 2)
							if (frag->total_residue() > 14)
								action_string = action_string+"_15+";
							else if (frag->total_residue() <= 4)
								action_string = action_string+"_0-4";
							else
								action_string = action_string+"_5-14";
						else
							action_string = action_string+"_MEGA";

						if (action == 3 || action == 4)
							superpose( *frag, pose, R, preT, postT);
	
						// xyz copy fragment to pose
						for (int i=1; i<=frag->total_residue(); ++i)
						for (int j=1; j<=frag->residue(i).natoms(); ++j) {
							core::id::AtomID src(j,i), tgt(j, frag->pdb_info()->number(i));
							pose.set_xyz( tgt, postT + (R*(frag->xyz( src )-preT)) );
						}
					}


					// //////   //
					// MC stuff //
					// //////   //
					try {
						(*min3_scorefxn_)(pose);
						minimizer.run( pose, mm, *min3_scorefxn_, options );
	
						if (m>1)
							mc->boltzmann( pose , action_string );
						else { // m==1
							minimizer.run( pose, mm, *min3_scorefxn_, options_minilbfgs );
							mc->boltzmann( pose , action_string );
						}
					} catch( utility::excn::EXCN_Base& excn ) {  // bad hbond shit
						mc->recover_low(pose);
					}
					if (n%100 == 0) {
						mc->show_scores();
						mc->show_counters();
					}

				}
				mc->recover_low(pose);
	
	
				// evaluator
				if ( option[ in::file::native ].user() && native_) {
					// align first time
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
					core::Real gdtmm = core::scoring::xyz_gdtmm( p1a, p2a, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
					std::cout << "CYCLE " << m << "  GDTMM = " << gdtmm
							  << " (" << m_1_1 << "," << m_2_2 << "," << m_3_3 << "," << m_4_3 << "," << m_7_4 << ")" << std::endl;
				}
	
				if (m==5) {
					(*min3_scorefxn_)(pose); minimizer.run( pose, mm, *min3_scorefxn_, options_lbfgs );
					(*min2_scorefxn_)(pose); minimizer.run( pose, mm, *min2_scorefxn_, options_lbfgs );
					(*min1_scorefxn_)(pose); minimizer.run( pose, mm, *min1_scorefxn_, options_lbfgs );
					(*min2_scorefxn_)(pose); minimizer.run( pose, mm, *min2_scorefxn_, options_lbfgs );
					(*min3_scorefxn_)(pose); minimizer.run( pose, mm, *min3_scorefxn_, options_lbfgs );
	
					// evaluator
					if ( option[ in::file::native ].user() && native_) {
						// align first time
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
						core::Real gdtmm = core::scoring::xyz_gdtmm( p1a, p2a, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
						std::cout << "CYCLE " << m << "+min  GDTMM = " << gdtmm
								  << " (" << m_1_1 << "," << m_2_2 << "," << m_3_3 << "," << m_4_3 << "," << m_7_4 << ")" << std::endl;
					}
				}
			}
		} else { // nmacrocycles==0
			if ( option[ in::file::native ].user() && native_) {
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
				core::Real gdtmm = core::scoring::xyz_gdtmm( p1a, p2a, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
				std::cout << "FINAL GDTMM = " << gdtmm
						  << " (" << m_1_1 << "," << m_2_2 << "," << m_3_3 << "," << m_4_3 << "," << m_7_4 << ")" << std::endl;
			}

			(*min3_scorefxn_)(pose); minimizer.run( pose, mm, *min3_scorefxn_, options_lbfgs );
			(*min2_scorefxn_)(pose); minimizer.run( pose, mm, *min2_scorefxn_, options_lbfgs );
			(*min1_scorefxn_)(pose); minimizer.run( pose, mm, *min1_scorefxn_, options_lbfgs );
			(*min2_scorefxn_)(pose); minimizer.run( pose, mm, *min2_scorefxn_, options_lbfgs );
			(*min3_scorefxn_)(pose); minimizer.run( pose, mm, *min3_scorefxn_, options_lbfgs );

			if ( option[ in::file::native ].user() && native_) {
				int n_atoms;
				ObjexxFCL::FArray2D< core::Real > p1a, p2a;
				protocols::comparative_modeling::gather_coords( pose, *native_, *aln_, n_atoms, p1a, p2a );

				core::Real m_1_1, m_2_2, m_3_3, m_4_3, m_7_4;
				core::Real gdtmm = core::scoring::xyz_gdtmm( p1a, p2a, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
				std::cout << "FINAL+min  GDTMM = " << gdtmm
						  << " (" << m_1_1 << "," << m_2_2 << "," << m_3_3 << "," << m_4_3 << "," << m_7_4 << ")" << std::endl;
			}
		}

		lowres_scorefxn_->set_weight( core::scoring::cart_bonded, max_cart );
		(*lowres_scorefxn_)(pose);
	}

	virtual std::string get_name() const {
		return "CustomMover";
	}

private:
	utility::vector1< core::pose::PoseOP > pose1s_,pose2s_;
	//utility::vector1< int > pose1_insert_points_, pose2_insert_points_;
	core::scoring::ScoreFunctionOP lowres_scorefxn_, min1_scorefxn_, min2_scorefxn_, min3_scorefxn_;

	core::fragment::FragSetOP fragments_;
	boost::unordered_map<core::Size, core::fragment::Frame> library_;

	// native pose, aln
	core::pose::PoseOP native_;
	core::sequence::SequenceAlignmentOP aln_;
};

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	SequenceMoverOP seq( new SequenceMover() );
	//seq->add_mover( new protocols::moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );
	seq->add_mover( new protocols::moves::ConstraintSetMover() );
	seq->add_mover( new CustomMover() );

	try{
		protocols::jd2::JobDistributor::get_instance()->go( seq );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
	NEW_OPT(fpd::templates, "templates", utility::vector1<utility::file::FileName >(0));
	NEW_OPT(fpd::fragments, "fragments", utility::vector1<utility::file::FileName >(0));
	NEW_OPT(fpd::frag9, "frag9", "");
	NEW_OPT(fpd::ncycles, "ncycles", 500);
	NEW_OPT(fpd::nmacrocycles, "nmacrocycles", 5);
	NEW_OPT(fpd::subfraglen, "subfraglen", 9);
	NEW_OPT(fpd::minfraglen, "minfraglen", 4);
	NEW_OPT(fpd::movie, "movie", false);

	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
}


