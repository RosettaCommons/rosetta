// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/pack/min_pack.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/CrystInfo.hh>
#include <core/pose/motif/reference_frames.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Remarks.hh>
#include <core/pose/util.hh>
#include <core/pose/selection.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/vdwaals/VDW_Energy.hh>
#include <core/types.hh>
#include <devel/init.hh>


#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/random/random.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/fourier/FFT.hh>
#include <numeric/constants.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>

#include <boost/foreach.hpp>
#define foreach_ BOOST_FOREACH

// option includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>


using namespace core;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using utility::vector1;
using utility::file::FileName;
using basic::T;
using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "calc_ssm_energies" );

OPT_1GRP_KEY(Boolean, ssm, interface)
OPT_1GRP_KEY(Boolean, ssm, packing)
OPT_1GRP_KEY(Boolean, ssm, verbose)
OPT_1GRP_KEY(Boolean, ssm, with_ss)
OPT_1GRP_KEY(IntegerVector, ssm, parallel)


//
// helper function, get interface residues (directional)
void
get_interface_residues( core::pose::Pose & pose, utility::vector1< bool > &interface, core::Real K) {
	using namespace core;
	using namespace core::scoring;
	using namespace core::conformation::symmetry;

	core::Real b=0.28;

	core::Size nres = pose.total_residue();

	// make pose polyA
	core::pose::Pose pose_working = pose;
	utility::vector1< Size > protein_residues;
	for ( Size i=1, i_end = nres; i<= i_end; ++i )
		if ( pose.residue( i ).is_protein() )
			protein_residues.push_back( i );
	protocols::toolbox::pose_manipulation::construct_poly_XXX_pose( "ALA", pose_working, protein_residues, false, false, false );

	interface.clear();
	interface.resize( pose.total_residue(), false );

	// initialize energy graph
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	(*sfxn)(pose_working);
	Energies & energies( pose_working.energies() );
	EnergyGraph & energy_graph( energies.energy_graph() );

	Size ninterface=0;
	for ( Size i=1, i_end = nres; i<= i_end; ++i ) {
		conformation::Residue const & rsd1( pose_working.residue( i ) );

		if ( pose.pdb_info()->chain(i) != pose.pdb_info()->chain(1) ) continue;

		for ( graph::Graph::EdgeListIter
				iru  = energy_graph.get_node(i)->edge_list_begin(),
				irue = energy_graph.get_node(i)->edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge & edge( static_cast< EnergyEdge & > (**iru) );

			Size const j = edge.get_other_ind( i );
			conformation::Residue const & rsd2( pose_working.residue( j ) );

			if ( i==j ) continue;  // don't think this is necessary
			if ( !rsd1.is_protein() || !rsd2.is_protein() ) continue;
			if ( pose.pdb_info()->chain(i) == pose.pdb_info()->chain(j) ) continue;

			// if CB-CB distance < 8A and CA-CB-CB angles are >75 deg then design
			core::Real dist = (rsd1.atom("CB").xyz() - rsd2.atom("CB").xyz()).length();
			core::Real angle1 = numeric::angle_degrees(rsd1.atom("CA").xyz(), rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz() ) ;
			core::Real angle2 = numeric::angle_degrees(rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz(), rsd2.atom("CA").xyz() ) ;
			core::Real angle_tgt = K*exp(b*dist);
			if (angle_tgt < 180 && angle1 > angle_tgt && angle2 > angle_tgt) {
				interface[i] = interface[j] = true;
				ninterface++;
			}
		}
	}

	if (ninterface==0) {
		utility_exit_with_message( "No interface residues found!" );
	}
}


//
// helper function, get neighbor residues (directional)
void
get_neighbor_residues(
		core::pose::Pose & pose,
		core::Size res_i,
		utility::vector1< bool > &neighbor,
		core::Real K) {
	using namespace core;
	using namespace core::scoring;

	core::Real b=0.28;

	core::Size nres = pose.total_residue();

	// make pose polyA
	core::pose::Pose pose_working = pose;
	utility::vector1< Size > protein_residues;
	for ( Size i=1, i_end = nres; i<= i_end; ++i )
		if ( pose.residue( i ).is_protein() )
			protein_residues.push_back( i );
	protocols::toolbox::pose_manipulation::construct_poly_XXX_pose( "ALA", pose_working, protein_residues, false, false, true );

	neighbor.clear();
	neighbor.resize( nres, false );
	neighbor[res_i] = true;

	for ( Size res_j=1, j_end = nres; res_j<= j_end; ++res_j ) {
		conformation::Residue const & rsd1( pose_working.residue( res_i ) );
		conformation::Residue const & rsd2( pose_working.residue( res_j ) );

		if ( res_i==res_j ) continue;
		if ( !rsd1.is_protein() || !rsd2.is_protein() ) continue;

		core::Real dist = (rsd1.atom("CB").xyz() - rsd2.atom("CB").xyz()).length();
		core::Real angle1 = numeric::angle_degrees(rsd1.atom("CA").xyz(), rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz() ) ;
		core::Real angle2 = numeric::angle_degrees(rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz(), rsd2.atom("CA").xyz() ) ;

		core::Real angle_tgt = K*exp(b*dist);

		if (angle_tgt < 180 && angle1 > angle_tgt && angle2 > angle_tgt) {
			neighbor[res_j] = true;
		}
	}
}


class Packing_energies : public protocols::moves::Mover {
public:
	Packing_energies() {
		// get scorefunction, set ref to 0
		sfxn_ = core::scoring::get_score_function();
		//sfxn_->set_weight( core::scoring::ref, 0.0 );

		K_=12.0;
		NCYC_=1;
		NEXP_=0;

		// quick and dirty parallelization
		if (basic::options::option[basic::options::OptionKeys::ssm::parallel].user()){
			utility::vector1< core::Size > parallel_args = basic::options::option[basic::options::OptionKeys::ssm::parallel]();
			runtime_assert( parallel_args.size() == 2);
			runtime_assert( parallel_args[1]<=parallel_args[2] );
			i_ = parallel_args[1];
			j_ = parallel_args[2];
		} else {
			i_ = j_ = 1;
		}
	}


	///
	void
	apply(core::pose::Pose &pose) {
		// load packer task from command line
		core::Size nres = pose.total_residue();

		// read resfile
		core::pack::task::TaskFactoryOP task ( new core::pack::task::TaskFactory );
		if (basic::options::option[basic::options::OptionKeys::packing::resfile].user()){
			task->push_back( core::pack::task::operation::TaskOperationOP( new core::pack::task::operation::ReadResfile ) );
		}
		core::pack::task::PackerTaskOP ptask_resfile = task->create_task_and_apply_taskoperations( pose );

		// restrict to interface if requested
		if (basic::options::option[basic::options::OptionKeys::ssm::interface].user()){
			utility::vector1< bool > interface;
			get_interface_residues( pose, interface, K_);
			ptask_resfile->restrict_to_residues(interface);
		}

		for ( Size i_res=1; i_res <= nres; ++i_res ) {
			if (j_>1 && ( (i_res%j_) != (i_%j_) )) continue;

			bool design_i = ptask_resfile->design_residue( i_res );
			if (!design_i) continue;

			// find neighbor residues
			utility::vector1<bool> neighbor;
			get_neighbor_residues( pose, i_res, neighbor, K_);

			// find neighbors' neighbors
			utility::vector1<bool> neighbor_neighbor = neighbor;
			for ( Size j_res=1; j_res <= nres; ++j_res ) {
				if (neighbor[j_res]) {
					utility::vector1<bool> neighbor_j;
					get_neighbor_residues( pose, j_res, neighbor_j, K_);

					for ( Size k_res=1; k_res <= nres; ++k_res ) {
						neighbor_neighbor[k_res] = (neighbor_neighbor[k_res] || neighbor_j[k_res]);
					}
				}
			}

			TR << i_res << " ";
			if (  basic::options::option[basic::options::OptionKeys::ssm::verbose] ) {
				TR << pose.residue(i_res).aa() << " ";
			}

			for ( Size i_aa=1; i_aa <= (Size)core::chemical::num_canonical_aas; ++i_aa ) {
				utility::vector1<bool> allowed_aas( core::chemical::num_canonical_aas, false );
				allowed_aas[i_aa] = true;

				// make a one-res packer task
				core::pack::task::PackerTaskOP ptask_working (core::pack::task::TaskFactory::create_packer_task( pose ));
				ptask_working->or_include_current(false);
				ptask_working->restrict_to_residues(neighbor_neighbor);
				for ( Size j_res=1; j_res <= nres; ++j_res ) {
					if ( i_res == j_res ) {
						ptask_working->nonconst_residue_task( j_res ).restrict_absent_canonical_aas( allowed_aas );
					} else if (neighbor[j_res]) {
						; //do nothing
					} else if (neighbor_neighbor[j_res]) {
						ptask_working->nonconst_residue_task( j_res ).restrict_to_repacking();
					}
				}

				core::pose::Pose pose_copy = pose;
				core::pack::pack_rotamers( pose_copy, *sfxn_, ptask_working );

				core::Real score_ij = (*sfxn_)(pose_copy);
				score_ij -= pose_copy.energies().onebody_energies( i_res )[core::scoring::ref]; // subtract reference weight of center residue
				TR << score_ij << " ";
			} // foreach aa
			TR << std::endl;
		} // foreach res
	}

	virtual std::string get_name() const {
		return "Packing_energies";
	}

private:
	core::scoring::ScoreFunctionOP sfxn_;
	core::pack::task::PackerTaskOP ptask_;
	core::Real K_,NCYC_;
	core::Size NEXP_;

	// parallelize
	core::Size i_,j_;
};



////////////
////////////
////////////

class SSM_energies : public protocols::moves::Mover {
public:
	SSM_energies() {
		// get scorefunction, set ref to 0
		sfxn_ = core::scoring::get_score_function();
		sfxn_->set_weight( core::scoring::ref, 0.0 );

		if (  basic::options::option[basic::options::OptionKeys::ssm::with_ss] ) {
			sfxn_->set_weight( core::scoring::rama_prepro, 0.0 );
			sfxn_->set_weight( core::scoring::rama, 0.0 );
			sfxn_->set_weight( core::scoring::p_aa_pp, 0.0 );
		}

		K_=12.0;
		NCYC_=1;
	}

	void
	optimization_loop(
			core::pose::Pose & pose,
			core::scoring::ScoreFunctionOP sf,
			core::pack::task::PackerTaskOP ptask,
			core::kinematics::MoveMapOP mm) {

		core::pack::pack_rotamers( pose, *sf, ptask );

		core::optimization::MinimizerOptions options( "lbfgs_armijo", 1e-2, true, false, false );

		//options.nblist_auto_update( true ); //?
		core::optimization::AtomTreeMinimizer minimizer;
		(*sf)(pose);  // this needs to be here
		minimizer.run( pose, *mm, *sf, options );

	}


	///
	void
	apply(core::pose::Pose &pose) {
    core::scoring::dssp::Dssp dssp( pose );
    dssp.insert_ss_into_pose( pose );

		// load packer task from command line
		core::Size nres = pose.total_residue();

		// read resfile
		core::pack::task::TaskFactoryOP task ( new core::pack::task::TaskFactory );
		if (basic::options::option[basic::options::OptionKeys::packing::resfile].user()){
			task->push_back( core::pack::task::operation::TaskOperationOP( new core::pack::task::operation::ReadResfile ) );
		}
		core::pack::task::PackerTaskOP ptask_resfile = task->create_task_and_apply_taskoperations( pose );

		// restrict to interface if requested
		if (basic::options::option[basic::options::OptionKeys::ssm::interface].user()){
			utility::vector1< bool > interface;
			get_interface_residues( pose, interface, K_);
			ptask_resfile->restrict_to_residues(interface);
		}

		core::scoring::ScoreFunction Srama, Spaapp;
		Srama.set_weight(core::scoring::rama_prepro,1.0);
		Spaapp.set_weight(core::scoring::p_aa_pp,1.0);

		for ( Size i_res=1; i_res <= nres; ++i_res ) {
			bool design_i = ptask_resfile->design_residue( i_res );
			if (!design_i) continue;

			// find neighbor residues
			utility::vector1<bool> neighbor;
			get_neighbor_residues( pose, i_res, neighbor, K_);

			// set up movemap
			core::kinematics::MoveMapOP mm(new core::kinematics::MoveMap);
			mm->set_jump(true); mm->set_chi(false); mm->set_bb(false);
			for ( Size j_res=1; j_res <= nres; ++j_res ) {
				if (i_res == j_res || neighbor[j_res] ) mm->set_chi(j_res, true);
			}

			TR << i_res << " ";

			if (  basic::options::option[basic::options::OptionKeys::ssm::verbose] ) {
				TR << pose.residue(i_res).aa() << " ";
			}
			if (  basic::options::option[basic::options::OptionKeys::ssm::with_ss] ) {
				TR << pose.secstruct(i_res) << " ";
			}

			utility::vector1<core::Real> Es, ramas,paapps;
			for ( Size i_aa=1; i_aa <= (Size)core::chemical::num_canonical_aas; ++i_aa ) {
				utility::vector1<bool> allowed_aas( core::chemical::num_canonical_aas, false );
				allowed_aas[i_aa] = true;

				// make a one-res packer task
				core::pack::task::PackerTaskOP ptask_working (core::pack::task::TaskFactory::create_packer_task( pose ));
				ptask_working->or_include_current(false);
				ptask_working->restrict_to_residues(neighbor);
				for ( Size j_res=1; j_res <= nres; ++j_res ) {
					if ( i_res == j_res ) {
						ptask_working->nonconst_residue_task( j_res ).restrict_absent_canonical_aas( allowed_aas );
					} else if (neighbor[j_res]) {
						ptask_working->nonconst_residue_task( j_res ).restrict_to_repacking();
					}
				}

				// optimize
				core::pose::Pose pose_copy = pose;
				optimization_loop( pose_copy, sfxn_, ptask_working, mm );
				Es.push_back( (*sfxn_)(pose_copy) );
				ramas.push_back( Srama(pose_copy) );
				paapps.push_back( Spaapp(pose_copy) );
			} // foreach aa

			for ( Size i_aa=1; i_aa <= (Size)core::chemical::num_canonical_aas; ++i_aa ) TR << Es[i_aa] << " ";
			if (  basic::options::option[basic::options::OptionKeys::ssm::with_ss] ) {
				for ( Size i_aa=1; i_aa <= (Size)core::chemical::num_canonical_aas; ++i_aa ) TR << ramas[i_aa] << " ";
				for ( Size i_aa=1; i_aa <= (Size)core::chemical::num_canonical_aas; ++i_aa ) TR << paapps[i_aa] << " ";
			}
			TR << std::endl;
		} // foreach res
	}

	virtual std::string get_name() const {
		return "SSM_energies";
	}

private:
	core::scoring::ScoreFunctionOP sfxn_;
	core::pack::task::PackerTaskOP ptask_;
	core::Real K_,NCYC_;
};



///
int main( int argc, char * argv [] )
{
	using namespace protocols::moves;
	using namespace protocols;
	using namespace protocols::jd2;

	try {
		NEW_OPT(ssm::interface, "interface", false);
		NEW_OPT(ssm::packing, "packing", false);
		NEW_OPT(ssm::with_ss, "with_ss", false);
		NEW_OPT(ssm::verbose, "verbose", false);
		NEW_OPT(ssm::with_ss, "with_ss", false);
		NEW_OPT(ssm::parallel, "parallel", utility::vector1<core::Size>());

		devel::init(argc, argv);

		SequenceMoverOP seq( new SequenceMover() );

		if (basic::options::option[basic::options::OptionKeys::ssm::packing]()) {
			seq->add_mover( MoverOP(new Packing_energies()) );
		} else {
			seq->add_mover( MoverOP(new SSM_energies()) );
		}

		// main loop
		protocols::jd2::JobDistributor::get_instance()->go( seq );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

