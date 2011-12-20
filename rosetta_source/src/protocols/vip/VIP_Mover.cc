//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file VIP_Mover.cc
/// @brief Void Identification and Packing Mover. Identifies buried voids in an input structure
/// and attempts to fill these using fixed backbone packing and the GOE energy
/// @author ben bborgo@genetics.wustl.edu

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/relax/MiniRelax.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <core/scoring/Energies.hh>

#include <core/init.hh>
#include <core/types.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/packstat/types.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/packstat/packing_score_params.hh>
#include <core/scoring/packstat/AtomRadiusMap.hh>
#include <core/scoring/packstat/SimplePDB.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <core/conformation/Residue.hh>
#include <numeric/all.hh>
#include <utility/all.hh>

#include <basic/Tracer.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/MoverContainer.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
//#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/pdb/file_data.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cp.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <utility/options/keys/OptionKey.hh>
#include <protocols/simple_moves/AddCavitiesMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <utility/io/izstream.hh>
#include <ObjexxFCL/format.hh>
#include <utility/string_util.hh>

#include <protocols/vip/VIP_Mover.hh>
#include <protocols/vip/VIP_Report.hh>

#include <cmath>




namespace protocols {
namespace vip {


static basic::Tracer TR("VIP");

                VIP_Mover::VIP_Mover() {}
		VIP_Mover::VIP_Mover(
        	        core::pose::Pose p,
	                core::pose::Pose cp,
                	core::pose::Pose fp,
                	utility::vector1<core::pose::Pose> tp,
                	utility::vector1<core::pose::Pose> vp,
                	core::Size nc,
                	utility::vector1<core::Size> cb,
                	utility::vector1<std::string> fm,
                	utility::vector1<core::Size> vn,
                	utility::vector1<core::Size> vm,
			core::Real fe){
				initial_pose = p;
		                cavity_pose = cp;
                		final_pose = fp;
                		temp_poses = tp;
				favorable_poses = vp;
                		number_cavities = nc;
				cavity_balls = cb;
                		favorable_mutations = fm;
				void_neighbors = vn;
				void_mutatables = vm;
				final_energy = fe;}
		VIP_Mover::~VIP_Mover(){}

		void VIP_Mover::set_initial_pose( core::pose::Pose pose ){
		        using namespace basic::options;
			using namespace basic::options::OptionKeys;
//			core::pose::Pose pose;
//			core::io::pdb::build_pose_from_pdb_as_is( pose, option[ OptionKeys::in::file::s ]().vector().front() );
			initial_pose = pose;}


		void VIP_Mover::minimize_conformation(){
			using namespace basic::options;
			using namespace basic::options::OptionKeys;

			core::pose::Pose pose = initial_pose;
			core::scoring::ScoreFunctionOP sf2 = core::scoring::ScoreFunctionFactory::create_score_function( option[cp::minimizer_score_fxn] );
			core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
        		movemap->set_jump(false);
        		movemap->set_chi(true);
        		movemap->set_bb(true);
        		protocols::moves::MoverOP min_native = new protocols::simple_moves::MinMover( movemap, sf2, "dfpmin_armijo_nonmonotone", 1e-2, true );
        		min_native->apply( pose );
			initial_pose = pose;}

                void VIP_Mover::compute_number_cavities(){
			core::Size numsuck = 0;
			for( core::Size a = 1 ; a < (cavity_pose.total_residue()+1) ; a++ ){
                		if( cavity_pose.residue(a).name() == "SUCK" ){
                        	numsuck = numsuck + 1;}};
       			number_cavities = numsuck;}

                void VIP_Mover::apply_holes(){
			core::pose::Pose pose = initial_pose;
        		protocols::simple_moves::AddCavitiesMover cavget;
        		cavget.apply( pose );
			cavity_pose = pose;}

                void VIP_Mover::get_cavity_positions(){
			utility::vector1<core::Size> cav_positions;
			for( core::Size i = 1; i <= cavity_pose.total_residue(); i++ ){
				if( cavity_pose.residue(i).name() == "SUCK" ){
					cav_positions.push_back(i);}}
			cavity_balls = cav_positions;}

                void VIP_Mover::dump_pdb_to_file( core::pose::Pose posey, std::string filename ){
		        posey.dump_pdb(filename);}

                void VIP_Mover::get_neighbors(){
                        using namespace basic::options;
                        using namespace basic::options::OptionKeys;
		    utility::vector1<core::Size> neighbors;

			if ( cavity_balls.size() < 1 ) {
				initial_pose.dump_pdb("intermediate.pdb");
                        	TR.Error << "No more cavities! Dumping structure to intermediate.pdb" << std::endl;
                        return;}

			for( core::Size i = 1; i <= cavity_balls.size(); i++ ){
			  numeric::xyzVector<core::Real> cav_center = cavity_pose.residue( cavity_balls[i]).xyz(1);
			    for( core::Size a = 1; a <= initial_pose.total_residue(); a++ ){
				for( core::Size b = 1; b <= initial_pose.residue(a).nheavyatoms(); b++ ){
			numeric::xyzVector<core::Real> test_position = initial_pose.residue(a).xyz(b);
			if( test_position.distance(cav_center) <= option[cp::cutoff] ){
				neighbors.push_back( a );
				neighbors.push_back( b );}}}}
			void_neighbors = neighbors;}


		void VIP_Mover::cull_mutatable_residues(){
			utility::vector1<core::Size> mutatable_residues;
			for( core::Size i = 1; i <= void_neighbors.size(); i+=2 ){
		if( (initial_pose.residue(void_neighbors[i]).is_surface() == false) && (initial_pose.residue(void_neighbors[i]).is_polar() == false) && (initial_pose.residue(void_neighbors[i]).atom_is_backbone(void_neighbors[i+1]) == false ) ){
					mutatable_residues.push_back(void_neighbors[i]);}}
			utility::vector1<core::Size> mm_res;
		for( core::Size ii = 1; ii < mutatable_residues.size(); ii++ ){
			if( mutatable_residues[ii] != mutatable_residues[ii+1] ){
				mm_res.push_back( mutatable_residues[ii] );}}

			void_mutatables = mm_res;}


                void VIP_Mover::try_point_mutants(){
                        using namespace basic::options;
                        using namespace basic::options::OptionKeys;

			for( core::Size i = 1; i <= void_mutatables.size(); i++ ){
        			temp_poses.push_back(initial_pose);}

        		core::pack::task::ResfileCommandOP command = new core::pack::task::APOLAR;

		for( core::Size aa = 1; aa <= void_mutatables.size(); aa++ ){
			core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( temp_poses[aa]));
	        	core::scoring::ScoreFunctionOP score_fxn = core::scoring::ScoreFunctionFactory::create_score_function( option[cp::pack_sfxn] );
			core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
		for( core::Size j = 1; j <= temp_poses[aa].total_residue(); j++ ){
			if( j != void_mutatables[aa] ){
        	task->nonconst_residue_task(j).prevent_repacking();}
			else{
        	task->nonconst_residue_task(void_mutatables[aa]).or_ex1(true);
        	task->nonconst_residue_task(void_mutatables[aa]).or_ex2(true);
        	task->nonconst_residue_task(void_mutatables[aa]).or_ex3(true);
        	command->residue_action(*task,void_mutatables[aa]);}}
        protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover(score_fxn, task);
        pack_mover->apply(temp_poses[aa]);}
}

		void VIP_Mover::print_pack_report(){
			VIP_Report();
			VIP_Report vip_report;

			vip_report.set_goe_native( initial_pose );
			vip_report.set_goe_repack( favorable_poses );
			vip_report.get_GOE_repack_report();}

                void VIP_Mover::sort_fill_energies(){
                        using namespace basic::options;
                        using namespace basic::options::OptionKeys;

		         core::scoring::ScoreFunctionOP score_fxn = core::scoring::ScoreFunctionFactory::create_score_function( option[cp::pack_sfxn] );
			protocols::simple_moves::ScoreMoverOP score_em = new protocols::simple_moves::ScoreMover(score_fxn);
			score_em->apply( initial_pose );

			core::Real baseE = initial_pose.energies().total_energy();
			core::pose::Pose basePose = initial_pose;
			for( core::Size i = 1; i <= temp_poses.size(); i++ ){
				core::Real tempE = temp_poses[i].energies().total_energy();
				if( tempE < baseE ){
					favorable_poses.push_back( temp_poses[i] );}}

			if( option[ cp::print_reports ] ){
				print_pack_report();}}

		void VIP_Mover::print_relax_report(){
			VIP_Report();
			VIP_Report vip_report;

			vip_report.set_goe_native( initial_pose );
			vip_report.set_goe_relax( favorable_poses );
			vip_report.get_GOE_relaxed_report();
                	vip_report.get_GOE_packstat_report();
}

		void VIP_Mover::relax_favorable_poses(){
                        using namespace basic::options;
                        using namespace basic::options::OptionKeys;
			core::scoring::ScoreFunctionOP relax_score_fxn = core::scoring::ScoreFunctionFactory::create_score_function( option[cp::relax_sfxn] );

			std::string rmover = option[ cp::relax_mover ];

if( rmover == "relax" ){
			for( core::Size i = 1; i <= favorable_poses.size(); i++ ){
				protocols::moves::MoverOP relaxmover = new protocols::relax::FastRelax( relax_score_fxn );
				relaxmover->apply(favorable_poses[i]);}}
if( rmover == "cst_relax" ){
			for( core::Size i = 1; i <= favorable_poses.size(); i++ ){
				protocols::moves::MoverOP cstrelaxmover = new protocols::relax::MiniRelax( relax_score_fxn );
				cstrelaxmover->apply(favorable_poses[i]);}
}}

		void VIP_Mover::sort_relaxed_poses(){
                        using namespace basic::options;
                        using namespace basic::options::OptionKeys;

			core::Real bestE = 99999;
			core::Size bestP = 0;
   			for( core::Size i = 1; i <= favorable_poses.size(); i++ ){
				if( favorable_poses[i].energies().total_energy() < bestE ){
					bestE = favorable_poses[i].energies().total_energy();
					bestP = i;}}
			final_pose = favorable_poses[bestP];
			dump_pdb_to_file( final_pose, "final.pdb" );
			final_energy = bestE;

			if( option[ cp::print_reports ] ){
				print_relax_report();}}

		void VIP_Mover::nook_finder(){
			minimize_conformation();
			apply_holes();
			get_cavity_positions();
			get_neighbors();
			cull_mutatable_residues();}

		void VIP_Mover::cranny_packer(){
			try_point_mutants();
			sort_fill_energies();
			relax_favorable_poses();
			sort_relaxed_poses();}

		void VIP_Mover::apply(){
			nook_finder();
			cranny_packer();}
}}
