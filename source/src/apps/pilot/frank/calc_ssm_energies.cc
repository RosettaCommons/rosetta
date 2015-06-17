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
static thread_local basic::Tracer TR( "calc_ssm_energies" );

OPT_1GRP_KEY(Boolean, ssm, interface)


class SSM_energies : public protocols::moves::Mover {
public:
	SSM_energies() {
		// get scorefunction, set ref to 0
		sfxn_ = core::scoring::get_score_function();
		sfxn_->set_weight( core::scoring::ref, 0.0 );

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

		// pack+relax
		//protocols::relax::FastRelax relax_prot( sf, NCYC_ );
		//relax_prot.min_type("lbfgs_armijo_nonmonotone");
		//relax_prot.cartesian( true );
		//relax_prot.set_movemap( mm );
		//relax_prot.apply(pose);
	}


	//
	void
	get_neighbor_residues(
			core::pose::Pose & pose,
			core::Size res_i,
			utility::vector1< bool > &neighbor) {
		using namespace core;
		using namespace core::scoring;

		core::Real K=K_;
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

	//
	void
	get_interface_residues( Pose & pose, utility::vector1< bool > &interface) {
		using namespace core;
		using namespace core::scoring;
		using namespace core::conformation::symmetry;

		core::Real K=K_;
		core::Real b=0.28;

		core::Size nres = pose.total_residue();

		// make pose polyA
		Pose pose_working = pose;
		utility::vector1< Size > protein_residues;
		for ( Size i=1, i_end = nres; i<= i_end; ++i )
			if ( pose.residue( i ).is_protein() )
				protein_residues.push_back( i );
		protocols::toolbox::pose_manipulation::construct_poly_XXX_pose( "ALA", pose_working, protein_residues, false, false, false );

		interface.clear();
		interface.resize( pose.total_residue(), false );

		// initialize energy graph
		(*sfxn_)(pose_working);
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
			get_interface_residues( pose, interface);
			ptask_resfile->restrict_to_residues(interface);
		}

		for ( Size i_res=1; i_res <= nres; ++i_res ) {
			bool design_i = ptask_resfile->design_residue( i_res );
			if (!design_i) continue;

			// find neighbor residues
			utility::vector1<bool> neighbor;
			get_neighbor_residues( pose, i_res, neighbor);

			// set up movemap
			core::kinematics::MoveMapOP mm(new core::kinematics::MoveMap);
			mm->set_jump(true); mm->set_chi(false); mm->set_bb(false);
			for ( Size j_res=1; j_res <= nres; ++j_res ) {
				if (i_res == j_res || neighbor[j_res] ) mm->set_chi(j_res, true);
			}

			TR << i_res << " ";
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
				TR << (*sfxn_)(pose_copy) << " ";
			} // foreach aa
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

		devel::init(argc, argv);

		SequenceMoverOP seq( new SequenceMover() );
		seq->add_mover( MoverOP(new SSM_energies()) );

		// main loop
		protocols::jd2::JobDistributor::get_instance()->go( seq );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

