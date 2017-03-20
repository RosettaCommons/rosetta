// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/ddg.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/database/open.hh>

#include <devel/init.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <ObjexxFCL/format.hh>

// C++ headers
#include <protocols/ddg/ddGMover.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <basic/Tracer.hh>

OPT_1GRP_KEY(Integer, ddg, bbnbr)
OPT_1GRP_KEY(Boolean, ddg, cartesian)
OPT_1GRP_KEY(Boolean, ddg, fd_mode)

//Auto Headers
using basic::T;
using basic::Error;
using basic::Warning;
static basic::Tracer TR("apps.pilot.wendao.ddg");

using namespace core;
using namespace pose;
using namespace scoring;
using namespace conformation;
using namespace kinematics;
using namespace ObjexxFCL::format;

typedef utility::vector1< core::chemical::AA > mutations;
typedef utility::vector1< bool > bools;

void
find_neighbors(
	bools const & is_mutated,
	Pose const & pose,
	bools & is_flexible,
	Real const heavyatom_distance_threshold = 6.0 )
{
	Size const nres( pose.size() );
	is_flexible = is_mutated;

	for ( Size i=1; i<= nres; ++i ) {
		Residue const & rsd1( pose.residue(i) );
		if ( rsd1.is_virtual_residue() ) continue;
		for ( Size j=1; j<= nres; ++j ) {
			if ( !is_mutated[j] ) continue;
			if ( is_flexible[i] ) break;
			Residue const & rsd2( pose.residue(j) );
			if ( rsd2.is_virtual_residue() ) continue;
			if ( rsd1.nbr_atom_xyz().distance_squared( rsd2.nbr_atom_xyz() ) <=
					numeric::square( rsd1.nbr_radius() + rsd2.nbr_radius() + heavyatom_distance_threshold ) ) {
				is_flexible[i] = true;
			}
		}
	}
}

void
find_neighbors_directional(
	bools const & neighborin,
	Pose const & pose,
	bools & neighbor,
	Real const K = 8.0
) {
	using namespace core;
	using namespace core::scoring;

	// interface parameters
	//   if (angle_degrees > K*exp(b*distance_angstrom) we have an interface interaction
	//   reduce K to increase # detected interface residues
	//      max dist = (1/b)*ln(180/k)
	//   at K= 10, maxdist = 10.32
	//         12             9.67
	//         14             9.12
	//         16             8.64
	core::Real b=0.28;

	core::Size nres = pose.total_residue();

	// make pose polyA
	Pose pose_working = pose;
	utility::vector1< Size > protein_residues;
	for ( Size i=1, i_end = nres; i<= i_end; ++i )
		if ( pose.residue( i ).is_protein() )
			protein_residues.push_back( i );
	protocols::toolbox::pose_manipulation::construct_poly_XXX_pose( "ALA", pose_working, protein_residues, false, false, false );

	neighbor = neighborin;

	for ( Size i=1, i_end = nres; i<= i_end; ++i ) {
		if ( ! neighborin[i] ) continue;
		neighbor[i] = true;

		for ( Size j=1, j_end = nres; j<= j_end; ++j ) {
			conformation::Residue const & rsd1( pose_working.residue( i ) );
			conformation::Residue const & rsd2( pose_working.residue( j ) );

			if ( i==j ) continue;
			if ( !rsd1.is_protein() || !rsd2.is_protein() ) continue;

			core::Real dist = (rsd1.atom("CB").xyz() - rsd2.atom("CB").xyz()).length();
			core::Real angle1 = numeric::angle_degrees(rsd1.atom("CA").xyz(), rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz() ) ;
			core::Real angle2 = numeric::angle_degrees(rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz(), rsd2.atom("CA").xyz() ) ;

			core::Real angle_tgt = K*exp(b*dist);

			if (angle_tgt < 180 && angle1 > angle_tgt && angle2 > angle_tgt) {
				neighbor[j] = true;
			}
		}
	}
}


/// @brief The input file is a list of mutation blocks.  Usually, this will be a set of point mutations.
/// where each "block" lists a single mutation.  However, it is possible to specify multiple mutations
/// together in a single block.
///
/// The file format is:
/// "total N"
/// followed by N blocks, where each block is
/// "M"
/// specifying followed by M lines of wt/resid/mutaa triples
/// "wtaa resid mutaa"
/// N, M and resid are all supposed to be integers.
/// wtaa, and mutaa are supposed to be 1-letter amino acid codes.
void
read_in_mutations(
	utility::vector1< mutations > & res_to_mut,
	bools & is_mutated, std::string filename,
	Pose & pose )
{
	std::ifstream inputstream;
	inputstream.open(filename.c_str());
	if ( inputstream.is_open() ) {
		int total;
		std::string total_keyword;
		inputstream >> total_keyword;
		assert(total_keyword.compare("total") == 0);

		inputstream >> total; //keep for cross-checking
		while ( !inputstream.eof() && total>0 ) {
			mutations current_mutation( pose.size(), core::chemical::aa_unk );
			int num_mutations;
			inputstream >> num_mutations;
			runtime_assert(num_mutations>0);
			while ( num_mutations > 0 ) {
				char wt; int resnum; char mut;
				inputstream >> wt >> resnum >> mut;
				TR.Debug << "wt is " << wt << " resnum is " << resnum << " and mut is " << mut << std::endl;
				runtime_assert( pose.residue(resnum).name1() == wt ); /// APL -- never use regular asserts when it comes to user input
				runtime_assert( core::chemical::oneletter_code_specifies_aa( mut ) ); /// APL -- input should specify the 1-letter code for an amino acid.
				core::chemical::AA mutation = core::chemical::aa_from_oneletter_code( mut );
				current_mutation[resnum] = mutation;
				is_mutated[resnum] = true;
				num_mutations--; total--;
			}
			TR.Debug << "end reading mutations for this" << std::endl;
			if ( num_mutations < 0 ) {
				TR.Error << "number of mutations mismatch! num_mutations < 0" << std::endl;
				runtime_assert(num_mutations==0);
			} else {
				res_to_mut.push_back(current_mutation);
			}
		}

		if ( total < 0 ) {
			TR.Error << "total number of mutations mismatch! total < 0" << std::endl;
			return;
		}
	}
}

void
compute_folding_energies(
	ScoreFunctionOP fa_scorefxn,
	Pose & pose,
	bools const & is_flexible,
	bools const & is_mutpos,
	Size bbnbrs=0 )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );
	fastrelax.cartesian( option[ ddg::cartesian ]()  );

	MoveMapOP movemap(new MoveMap);
	movemap->set_bb( false );
	movemap->set_chi( false );
	movemap->set_jump( false );

	bool flexbb = true;
	bool cartmin = true;

	if ( flexbb ) {
		runtime_assert( cartmin );
		for ( Size i=1; i<= pose.size(); ++i ) {
			if ( is_flexible[i] ) {
				movemap->set_chi( i, true );

				for ( Size j=0; j<=bbnbrs; j++ ) {
					if ( is_mutpos[i] ||
							( i+j <= pose.size() && is_mutpos[i+j] ) ||
							( i-j >= 1 && is_mutpos[i-j] ) ) {
						movemap->set_bb ( i, true );
					}
				}
			}
		}
	}

	fastrelax.set_movemap( movemap );
	fastrelax.apply(pose);
}


int
main( int argc, char * argv [] )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::task;
	using namespace protocols::moves;

	try {
		NEW_OPT(ddg::bbnbr, "bb neighbor", 0);
		NEW_OPT(ddg::cartesian, "cartesian", true);
		NEW_OPT(ddg::fd_mode, "mode", false);

		//init
		devel::init(argc, argv);
		Size num_iterations = option[ OptionKeys::ddg::iterations ]();
		bool fd_mode = option[ OptionKeys::ddg::fd_mode ]();
		double cutoff = fd_mode? 8.0 : 6.0;
		if ( basic::options::option[ OptionKeys::ddg::opt_radius].user() ) {
			cutoff = basic::options::option[ OptionKeys::ddg::opt_radius ]();
		}
		Size bbnbr = basic::options::option[ OptionKeys::ddg::bbnbr ]();

		// read the pose
		Pose pose;
		core::import_pose::pose_from_file( pose, basic::options::start_file(),
			core::import_pose::PDB_file );

		core::scoring::ScoreFunctionOP score_fxn;
		score_fxn = core::scoring::get_score_function();
		//TR << (*score_fxn)(pose) << std::endl;

		//interface mode? setting = jump number to use for interface
		Size interface_ddg = option[ OptionKeys::ddg::interface_ddg ]();

		//fd try to be smart .. look for interchain jump
		if( interface_ddg > pose.num_jump() ) {
			for (int i=1; i<=(int)pose.fold_tree().num_jump(); ++i) {
				kinematics::Edge jump_i = pose.fold_tree().jump_edge( i ) ;
				if (pose.pdb_info()->chain(jump_i.start()) != pose.pdb_info()->chain(jump_i.stop())) {
					runtime_assert( interface_ddg > pose.num_jump() );
					interface_ddg = i;
				}
			}
			TR << "Autosetting interface jump to " << interface_ddg << std::endl;
		}
		runtime_assert ( interface_ddg <= pose.num_jump() );

		bools is_mutated( pose.size(), false );   //mutated position
		bools is_flexible( pose.size(), false );  //repackable

		if ( !option[ OptionKeys::ddg::mut_file ].user() ) {
			//quit!
			runtime_assert(false);
		}
		TR.Debug << "reading in mutfile" << std::endl;
		std::string filename = option[OptionKeys::ddg::mut_file]();
		utility::vector1< mutations > res_to_mut;
		read_in_mutations( res_to_mut, is_mutated, filename, pose);
		TR.Debug << "mutations:" << res_to_mut.size() << std::endl;
		//the logic here is:
		//to be comparable, all the mutations should have the save dof!
		//so take the union set
		if (fd_mode) {
			find_neighbors_directional( is_mutated, pose, is_flexible, cutoff );
		} else {
			find_neighbors( is_mutated, pose, is_flexible, cutoff );
		}

		//output file, mark if rerun
		std::map< std::string, Size > before_jump_done_list;
		std::map< std::string, Size > after_jump_done_list;
		utility::io::ozstream ofp;
		utility::file::FileName parsefn(filename);
		std::string ofn = parsefn.base()+".ddg";
		if ( utility::file::file_exists(ofn) ) {
			TR << "This is a restart run!" << std::endl;
			utility::io::izstream ifp;
			ifp.open(ofn);
			if ( !ifp.good() ) utility_exit_with_message( "can't open cluster file " + ofn );
			while ( !ifp.eof() ) {
				std::string line;
				utility::io::getline( ifp, line );
				//parse the line
				std::istringstream istr(line);
				std::string cat, rd, mut;
				istr >> cat >> rd >> mut;
				mut = mut.substr(0, mut.size()-1);

				//skip if no keyword
				const size_t pos1 = rd.find("Round");
				const size_t pos2 = rd.find(':');
				if ( std::string::npos == pos1 || std::string::npos == pos2 ) continue;
				std::string nrdstr = rd.substr(pos1+5, pos2-5);
				Size nrd = std::atoi(nrdstr.c_str());

				if ( cat=="COMPLEX:" ) {
					if ( before_jump_done_list.count(mut)>0 ) {
						if ( before_jump_done_list[mut]<nrd ) before_jump_done_list[mut]=nrd;
					} else {
						before_jump_done_list[mut]=nrd;
					}
				} else if ( cat=="OPT_APART:" ) {
					if ( after_jump_done_list.count(mut) > 0 ) {
						if ( after_jump_done_list[mut] < nrd ) after_jump_done_list[mut]=nrd;
					} else {
						after_jump_done_list[mut]=nrd;
					}
				}
			}
			ifp.close();

			//append in previous file
			ofp.open_append(ofn);
		} else {
			//create new file
			ofp.open(ofn);
		}

		//start from 0(WT), unless "mut_only" is specifyied
		Size mstart=0;
		if ( option[OptionKeys::ddg::mut_only].user() ) {
			if ( option[OptionKeys::ddg::mut_only]() ) mstart=1;
		}
		Size mend=res_to_mut.size();
		if ( option[OptionKeys::ddg::wt_only].user() ) {
			if ( option[OptionKeys::ddg::wt_only]() ) mend=0;
		}
		for ( Size mi(mstart); mi<=mend; mi++ ) {
			Pose work_pose(pose);

			//generate tag
			std::ostringstream tag;
			if ( mi==0 ) {
				//WT, keep the pose as is
				tag << "WT";
			} else {
				tag << "MUT";
				for ( Size j=1; j <= is_mutated.size(); j++ ) {
					if ( res_to_mut[mi][j] != core::chemical::aa_unk ) {
						tag << "_" << j << res_to_mut[mi][j];
					}
				}
			}

			//figure out where to start, if restart
			Size bfj_rd_start(1);
			Size afj_rd_start(1);
			if ( before_jump_done_list.count(tag.str())>0 ) bfj_rd_start = before_jump_done_list[tag.str()]+1;
			if ( after_jump_done_list.count(tag.str())>0 ) afj_rd_start = after_jump_done_list[tag.str()]+1;

			//debug
			//std::cout << tag.str() << " bfj: " << bfj_rd_start << std::endl;
			//std::cout << tag.str() << " afj: " << afj_rd_start << std::endl;
			if ( bfj_rd_start>num_iterations ) {
				if ( interface_ddg==0 || afj_rd_start>num_iterations ) {
					TR << "Job has been done for " << tag.str() << std::endl;
					continue;
				}
			}

			TR << "Job is running for: " << tag.str() << std::endl;

			//MUT, point(s) mutation, mute_file
			if ( mi>0 ) {
				pack::task::PackerTaskOP packer_task(pack::task::TaskFactory::create_packer_task(work_pose));
				for ( Size j=1; j <= is_mutated.size(); j++ ) {
					if ( res_to_mut[mi][j] != core::chemical::aa_unk ) {
						bools restrict_to_aa( 20, false );
						restrict_to_aa[res_to_mut[mi][j]] = true;
						packer_task->nonconst_residue_task(j).restrict_absent_canonical_aas(restrict_to_aa);
					} else {
						packer_task->nonconst_residue_task(j).prevent_repacking();
					}
				}
				protocols::simple_moves::PackRotamersMoverOP packer(new protocols::simple_moves::PackRotamersMover( score_fxn, packer_task ));
				packer->apply(work_pose);
			}

			//monomer folding energy
			for ( Size r=bfj_rd_start; r<=num_iterations; r++ ) {
				Pose local_pose(work_pose);
				compute_folding_energies( score_fxn, local_pose, is_flexible, is_mutated, bbnbr );

				//output
				Real final_score( (*score_fxn)( local_pose ) );
				if ( option[OptionKeys::ddg::dump_pdbs]() ) {
					std::ostringstream dump_fn;
					dump_fn << tag.str() << "_bj" << r << ".pdb";
					local_pose.dump_pdb(dump_fn.str());
				}
				ofp << "COMPLEX:   Round" << r << ": " << tag.str() << ": " << F(9,3,final_score) << " "
					<< local_pose.energies().total_energies().weighted_string_of( score_fxn->weights() ) << std::endl;

				Size rb_jump(interface_ddg);
				protocols::rigid::RigidBodyTransMoverOP separate_partners( new protocols::rigid::RigidBodyTransMover( local_pose, rb_jump ) );
				separate_partners->step_size(1000.0);
				separate_partners->apply(local_pose);
				final_score = (*score_fxn)( local_pose );
				ofp << "APART:     Round" << r << ": " << tag.str() << ": " << F(9,3,final_score) << " "
					<< local_pose.energies().total_energies().weighted_string_of( score_fxn->weights() ) << std::endl;
			}

			//interface mode, seperate and score
			if ( interface_ddg > 0 ) {
				Size rb_jump(interface_ddg);
				protocols::rigid::RigidBodyTransMoverOP separate_partners( new protocols::rigid::RigidBodyTransMover( work_pose, rb_jump ) );
				separate_partners->step_size(1000.0);
				separate_partners->apply(work_pose);

				//repack or not?
				//seperate partners energy
				for ( Size r=afj_rd_start; r<=num_iterations; r++ ) {
					Pose local_pose(work_pose);
					compute_folding_energies( score_fxn, local_pose, is_flexible, is_mutated, bbnbr );

					//output
					Real const final_score( (*score_fxn)( local_pose ) );
					if ( option[OptionKeys::ddg::dump_pdbs]() ) {
						std::ostringstream dump_fn;
						dump_fn << tag.str() << "_aj" << r << ".pdb";
						local_pose.dump_pdb(dump_fn.str());
					}
					ofp << "OPT_APART: Round" << r << ": " << tag.str() << ": " << F(9,3,final_score) << " "
						<< local_pose.energies().total_energies().weighted_string_of( score_fxn->weights() ) << std::endl;
				}
			}
		}

		ofp.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

