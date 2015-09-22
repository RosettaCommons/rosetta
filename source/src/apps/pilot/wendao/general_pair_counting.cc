// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file general_pair_counting.cc
/// @brief generate octave "load_table" scripts for calculating env and pair term
/// @author Yuan Liu (wendao@uw.edu)

// Core Headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>

#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
//#include <core/scoring/TenANeighborGraph.hh>

// 
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.hh>

#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>

#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

// namespace
using namespace core;
using namespace core::import_pose;
using namespace core::io::pdb;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace core::pose;
using namespace core::chemical;
using namespace core::scoring;
using namespace core::kinematics;
using namespace protocols::moves;

// create a TaskFactory with the resfile
using namespace core::pack::task;
using namespace core::pack::dunbrack;
using namespace core::pack::dunbrack::cenrot;

// declare
static THREAD_LOCAL basic::Tracer TR( "pilot.wendao.cenrot" );
utility::vector1<std::string> rot_type_list;
utility::vector1<Size> rot_list_num;

void my_main();

// define
OPT_KEY(Boolean, restype_rot_dependent)
OPT_KEY(IntegerVector, env_dependent_pair)
OPT_KEY(Integer, pair_bin_number)
OPT_KEY(Real, pair_bin_width)
OPT_KEY(Integer, env_bin_number)
OPT_KEY(Integer, sequence_seperation)

OPT_KEY(Boolean, anisovdw)

int main( int argc, char * argv [] )
{
	NEW_OPT(restype_rot_dependent, 
		"make the restype rotamer dependent", false);
	NEW_OPT(env_dependent_pair, 
		"make the restype environment dependent by given list", 
		utility::vector1<core::Size>());
	NEW_OPT(pair_bin_number,
		"bin number of pair term", 30);
	NEW_OPT(pair_bin_width,
		"bin width of pair term", 0.5);
	NEW_OPT(env_bin_number,
		"bin number of environment", 40);
	NEW_OPT(sequence_seperation,
		"sequence distance exlusion cutoff", 9);
	NEW_OPT(anisovdw, 
		"output anisovdw parameter", false);

	devel::init(argc, argv);
	my_main();
	return 0;
}

Size get_env_bin_number(Size nbr_num)
{
	if (option[env_dependent_pair].user()) {
		for (Size i=1; i<=option[env_dependent_pair].size();i++) {
			if ((int)nbr_num<=option[env_dependent_pair][i]) {
				return i;
			}
		}
		//out of range
		return option[env_dependent_pair].size();
	}
	else {
		return 1;
	}
}

Size get_best_rotamer_index(Pose &p, Size ndx)
{
	if (!option[restype_rot_dependent]) return 0;

	RotamerLibrary const & rotamerlibCAP = * RotamerLibrary::get_instance();

	SingleResidueRotamerLibraryCAP residue_rotamer_library(
		rotamerlibCAP.get_rsd_library(p.residue(ndx).type()));

	if (residue_rotamer_library==0) return 1;

	SingleResidueCenrotLibraryCAP residue_cenrot_library(
		dynamic_cast< SingleResidueCenrotLibrary const * >(residue_rotamer_library.get()));

	Size closest_rot;
	Real closest_dis;
	//CentroidRotamerSampleData const &sampledata ( residue_cenrot_library->get_closest_rotamer( p.residue(ndx), closest_rot, closest_dis));
	residue_cenrot_library->get_closest_rotamer( p.residue(ndx), closest_rot, closest_dis);
	return closest_rot;
}

Size get_rotamer_type_index(std::string const & aa, Size nrot, Size nenv=0)
{
		std::stringstream rot_type;
		std::string aastr = aa.substr(0,3);
		if (aa.compare(0,3,"CYD")==0) aastr="CYS";

		//TR.Debug << "|" << aa.substr(0,3) << "| -->" << aastr << std::endl;

		if (option[restype_rot_dependent]) {
				rot_type << aastr << nrot;
		}
		else {
				rot_type << aastr;
		}

		if (option[env_dependent_pair].user()) {
				rot_type << " " << nenv;
		}

		Size i;
		for(i=1; i<=rot_type_list.size(); i++) {
				if (rot_type_list[i]==rot_type.str()) break;
				//TR.Debug << "rot_list " << i << ": " << rot_type_list[i] << std::endl;
		}

		if (i>rot_type_list.size()) {
				//TR.Debug << "Adding new type: " << rot_type.str() << std::endl;
				rot_type_list.push_back(rot_type.str());
				rot_list_num.push_back(0);
		}

		return i;
}

void fill_1D_vector(utility::vector1<Size> &t1, Size nx)
{
	for (Size i=1; i<=nx; i++) {
		t1.push_back(0);
	}
}

void fill_2D_vector(
	utility::vector1< utility::vector1<Size> > &t2,
	Size nx, Size ny
) {
	utility::vector1<Size> t1;
	fill_1D_vector(t1, ny);
	for (Size i=1; i<=nx; i++)	{
		t2.push_back(t1);
	}
}

void fill_3D_vector(
	utility::vector1< 
		utility::vector1< 
			utility::vector1<
				Size
			> 
		> 
	> &t3,
	Size nx, Size ny, Size nz
) {
	utility::vector1<utility::vector1<Size> > t2;
	fill_2D_vector(t2, ny, nz);
	for (Size i=1; i<=nx; i++)	{
		t3.push_back(t2);
	}
}

void fill_4D_vector(
	utility::vector1<
		utility::vector1<
			utility::vector1<
				utility::vector1<
					Size
				>
			>
		>
	> &t4,
	Size nx, Size ny, Size nz, Size no
) {
	utility::vector1 <utility::vector1<utility::vector1<Size> > > t3;
	fill_3D_vector(t3, ny, nz, no);
	for (Size i=1; i<=nx; i++)	{
		t4.push_back(t3);
	}
}

void my_main()
{
	ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set( "centroid_rot" );

	Size npdbs = option[ in::file::l ]().size();
	if (npdbs<=0) {
		TR.Error << "No pdb file found, use -in::file::l to specify a list of pdb!" << std::endl;
	}

	// build the count table
	utility::vector1 <utility::vector1<utility::vector1<Size> > > pair_table;
	utility::vector1 <utility::vector1<utility::vector1<Size> > > contact_table;
	utility::vector1 <utility::vector1<Size> > env_table;
	//theta 0~180, 18bin
	utility::vector1 <utility::vector1 <utility::vector1<utility::vector1<Size> > > > th1_table;
	utility::vector1 <utility::vector1 <utility::vector1<utility::vector1<Size> > > > th2_table;
	//dih -180~180, 36bin
	utility::vector1 <utility::vector1 <utility::vector1<utility::vector1<Size> > > > dih_table;

	Real wbin = option[pair_bin_width];
	Size nbin_pair = option[pair_bin_number];
	//Real Lbin = wbin * nbin_pair;
	Size nbin_env = option[env_bin_number]; //max neighbour number
	Size nrottype = option[restype_rot_dependent]?86:20; //magic number
	Size nenv_pair = get_env_bin_number(nbin_env);
	if (option[env_dependent_pair].user()) {
		nrottype *= nenv_pair;
	}
	Size nbin_contact = 150;
	Real wbin_c = 0.1;

	Size nbin_ang = 16;
	Size nbin_dih = 18;
	Real cos_ang0 = -1;
	Real d_ang = 2.0/16;
	Real dih0 = -180.0;
	Real d_dih = 20.0;

	//fill pair table
	fill_3D_vector(pair_table, nrottype, nrottype, nbin_pair);
	fill_3D_vector(contact_table, nrottype, nrottype, nbin_contact);
	//fill env table
	fill_2D_vector(env_table, nrottype, nbin_env);
	//fill angle table
	fill_4D_vector(th1_table, nrottype, nrottype, nbin_pair, nbin_ang);
	fill_4D_vector(th2_table, nrottype, nrottype, nbin_pair, nbin_ang);
	fill_4D_vector(dih_table, nrottype, nrottype, nbin_pair, nbin_dih);

	core::scoring::ScoreFunctionOP score_fxn = new core::scoring::ScoreFunction();
	score_fxn->set_weight(core::scoring::cen_rot_pair, 1.0);

	// go through each pdb file
	for (Size npdb=1; npdb<=npdbs; npdb++)
	{
		PoseOP pose = new Pose();
		Pose &p(*pose);
		pose_from_pdb( p, *rsd_set, option[ in::file::l ]()[npdb] );
		std::cerr << option[ in::file::l ]()[npdb] << std::endl;

		utility::vector1<Size> rotliblst; //residue type index for each residue

		(*score_fxn)(p);
		EnergyGraph const & energy_graph( p.energies().energy_graph() );

		for (Size i=1; i<=p.n_residue(); i++) {
			Size nbr_within_ten=0;

			conformation::Residue const & rsd1 ( p.residue(i) );
			for ( graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node(i)->const_edge_list_begin(),
					irue = energy_graph.get_node(i)->const_edge_list_end();
					iru != irue; ++iru ) {
				EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
				Size const j( edge->get_other_ind(i) );
				conformation::Residue const & rsd2 ( p.residue(j) );
				Real const cendist = rsd1.atom("CEN").xyz().distance(
								rsd2.atom(rsd2.nbr_atom()).xyz());
				if ( cendist <= 10.0 ) {
					nbr_within_ten++;
				}
			}

			//debug
			//if (nbr_within_ten<5 || nbr_within_ten>35) {
			//	std::cerr<< p.residue(i).name() << " " << i << ": nbr=" << nbr_within_ten << std::endl;
			//}

			//save its type index
			Size res_index = get_rotamer_type_index(
							p.residue(i).name(),
							get_best_rotamer_index(p,i),
							get_env_bin_number(nbr_within_ten));
			rotliblst.push_back( res_index );
			rot_list_num[res_index]++; //save number of each type

			//TR.Debug << rotliblst[i] << ": nbr= " <<  nbr_within_ten << std::endl;

			//add to env_table
			if (nbr_within_ten>=1 && (int)nbr_within_ten<=option[env_bin_number])
			{
					env_table[rotliblst[i]][nbr_within_ten]++;
			}
		}

		/// for each res-pair, cal dist
		for (Size i=1; i<=p.n_residue(); i++) {
			// for(Size j=i+option[sequence_seperation]; j<=p.n_residue(); j++) {
			// 	// pair
			// 	Real d = p.residue(i).atom("CEN").xyz().distance(p.residue(j).atom("CEN").xyz());
			// 	Size bin_p = int(d/wbin)+1;
			// 	if ( bin_p <= nbin_pair ) {
			// 		pair_table[rotliblst[i]][rotliblst[j]][bin_p]++;
			// 		pair_table[rotliblst[j]][rotliblst[i]][bin_p]++;
			// 	}
			// }

			// use energy graph may lost some of counting where r > 12 (very unlikely)
			// but scince our cutoff is 12, it's fine, i guess

			conformation::Residue const & rsd1 ( p.residue(i) );
		  if (rsd1.name3()=="CYS" &&  rsd1.type().is_disulfide_bonded() ) continue;

			for ( graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
					irue = energy_graph.get_node(i)->const_upper_edge_list_end();
					iru != irue; ++iru ) {
				//only for i<j
				EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
				Size const j( edge->get_second_node_ind() );

				//check the sequence seperation, skip local interaction
				if ((int)(j-i)<option[sequence_seperation]) continue;

				conformation::Residue const & rsd2 ( p.residue(j) );
		    if (rsd2.name3()=="CYS" &&  rsd2.type().is_disulfide_bonded()) continue;

				//save r
				Real d = rsd1.atom("CEN").xyz().distance(rsd2.atom("CEN").xyz());
				Size bin_r = int(d/wbin)+1;
				Size bin_c_r = int(d/wbin_c)+1;
				if ( bin_r <= nbin_pair ) {
			 		pair_table[rotliblst[i]][rotliblst[j]][bin_r]++;
			 		pair_table[rotliblst[j]][rotliblst[i]][bin_r]++;
					if ( bin_c_r <= nbin_contact ) {
						// save 0:0.1:15.0 bin counting
			 		  contact_table[rotliblst[i]][rotliblst[j]][bin_c_r]++;
			 		  contact_table[rotliblst[j]][rotliblst[i]][bin_c_r]++;
					}

				 	if (rsd1.name3()=="GLY" || rsd1.name3()=="ALA"
				 		|| rsd2.name3()=="GLY" || rsd2.name3()=="ALA") continue;

				 	//save angle
				 	core::kinematics::Stub::Vector ra, rb, rc, rd;
					ra = rsd1.atom("CB").xyz();
					rb = rsd1.atom("CEN").xyz();
					rc = rsd2.atom("CEN").xyz();
					rd = rsd2.atom("CB").xyz();

					Real ang1 = numeric::angle_radians(ra,rb,rc);
					Real cos_ang1 = cos(ang1);
					Real ang2 = numeric::angle_radians(rb,rc,rd);
					Real cos_ang2 = cos(ang2);
					Real dih = numeric::dihedral_degrees(ra,rb,rc,rd);
					Size bin_ang1 = int((cos_ang1-cos_ang0)/d_ang)+1;
					Size bin_ang2 = int((cos_ang2-cos_ang0)/d_ang)+1;
					Size bin_dih = int((dih-dih0)/d_dih)+1;

					th1_table[rotliblst[i]][rotliblst[j]][bin_r][bin_ang1]++;
					th2_table[rotliblst[i]][rotliblst[j]][bin_r][bin_ang2]++;
					dih_table[rotliblst[i]][rotliblst[j]][bin_r][bin_dih]++;

					//test aniso-vdw parameter
					if (option[anisovdw] && bin_r<=nbin_pair/2) {
							//close enough
							core::kinematics::Stub::Vector rCA1, rCA2;
							rCA1 = rsd1.atom("CA").xyz();
							rCA2 = rsd2.atom("CA").xyz();
							Real dih1 = numeric::dihedral_degrees(rc,rb,ra,rCA1);
							Real dih2 = numeric::dihedral_degrees(rb,rc,rd,rCA2);
							std::cout << "#: " << rsd1.name3() << " " << rsd2.name3() << " " << d << " " 
									<< ang1*numeric::constants::r::radians_to_degrees << " " << dih1 << std::endl;
							std::cout << "#: " << rsd2.name3() << " " << rsd1.name3() << " " << d << " " 
									<< ang2*numeric::constants::r::radians_to_degrees << " " << dih2 << std::endl;
					}
				}
			}
		}
	}


	//output contact table
	std::ofstream out_contact_dat("contact_table.txt");
	for (Size i=1; i<=rot_type_list.size(); i++) {
		for (Size j=i; j<=rot_type_list.size(); j++) {
			out_contact_dat << rot_type_list[i] << " " << rot_type_list[j] << " ";
			for (Size k=1; k<=nbin_contact; k++) {
				out_contact_dat << contact_table[i][j][k] << " ";
			}
			out_contact_dat << std::endl;
		}
	}

	///////////////////////////////////////////////////////////////////////////////////
	// output octave scritps for pair and env data
	///////////////////////////////////////////////////////////////////////////////////
	std::ofstream octave_scr("load_pair_env_table.m");

	//title
	octave_scr << "function [list, pair, env, numlist] = load_pair_env_table()" << std::endl;

	// output restype list
	octave_scr << "list = [ " << std::endl;
	for (Size i=1; i<rot_type_list.size(); i++) {
		octave_scr << " \"" << rot_type_list[i] << "\"; %" << i << std::endl;
	}
	octave_scr << " \"" << rot_type_list[rot_type_list.size()] << "\" ];" << std::endl;

	octave_scr << "numlist = [ " << std::endl;
	for (Size i=1; i<rot_list_num.size(); i++) {
		octave_scr << rot_list_num[i] << "; %" << i << std::endl;
	}
	octave_scr << rot_list_num[rot_list_num.size()] << " ];" << std::endl;

	//output pair table
	octave_scr << std::endl;
	octave_scr << "pair = zeros(" << nrottype << "," << nrottype 
		<< "," << nbin_pair << ");" << std::endl;
	for (Size i=1; i<=rot_type_list.size(); ++i)
	{
		for (Size j=1; j<=rot_type_list.size(); ++j)
		{
			octave_scr << "pair(" << i << "," << j << ",:) = [ ";
			for (Size k=1; k<=nbin_pair; ++k)
			{
				octave_scr << pair_table[i][j][k] << " ";
			}
			octave_scr << "];" << std::endl;
		}
	}

	//output env table
	octave_scr << std::endl;
	octave_scr << "env = zeros(" << rot_type_list.size() 
		<< "," << nbin_env << ");" << std::endl; 
	for (Size i=1; i<=rot_type_list.size(); i++) {
		octave_scr << "env(" << i << ",:)= [ ";
		for (Size j=1; j<=nbin_env; j++) {
			octave_scr << env_table[i][j] << " ";
		}
		octave_scr << "];" << std::endl;
	}

	///////////////////////////////////////////////////////////////////////////////////
	// output octave scritps for orientation term
	///////////////////////////////////////////////////////////////////////////////////
	std::ofstream octave_scr2("load_angle_r_table.m");
	octave_scr2 << "function [list, th1_table, th2_table, dih_table, numlist] = load_angle_r_table()" << std::endl;

	// output restype list
	octave_scr2 << "list = [ " << std::endl;
	for (Size i=1; i<rot_type_list.size(); i++) {
		octave_scr2 << " \"" << rot_type_list[i] << "\"; %" << i << std::endl;
	}
	octave_scr2 << " \"" << rot_type_list[rot_type_list.size()] << "\" ];" << std::endl;

	octave_scr2 << "numlist = [ " << std::endl;
	for (Size i=1; i<rot_list_num.size(); i++) {
		octave_scr2 << rot_list_num[i] << "; %" << i << std::endl;
	}
	octave_scr2 << rot_list_num[rot_list_num.size()] << " ];" << std::endl;

	octave_scr2 << "th1_table = zeros(" << nrottype << "," 
		<< nrottype << "," << nbin_pair << ","
		<< nbin_ang << ");" << std::endl;
	for (Size i=1; i<=rot_type_list.size(); ++i) {
		for (Size j=1; j<=rot_type_list.size(); ++j) {
			for (Size k=1; k<=nbin_pair; k++) {
				octave_scr2 << "th1_table(" << i << "," << j << "," << k << ",:) = [ ";
				for (Size a=1; a<=nbin_ang; a++) {
					octave_scr2 << th1_table[i][j][k][a] << " ";
				}
				octave_scr2 << "];" << std::endl;
			}
		}
	}

	octave_scr2 << "th2_table = zeros(" << nrottype << "," 
		<< nrottype << "," << nbin_pair << ","
		<< nbin_ang << ");" << std::endl;
	for (Size i=1; i<=rot_type_list.size(); ++i) {
		for (Size j=1; j<=rot_type_list.size(); ++j) {
			for (Size k=1; k<=nbin_pair; k++) {
				octave_scr2 << "th2_table(" << i << "," << j << "," << k << ",:) = [ ";
				for (Size a=1; a<=nbin_ang; a++) {
					octave_scr2 << th2_table[i][j][k][a] << " ";
				}
				octave_scr2 << "];" << std::endl;
			}
		}
	}

	octave_scr2 << "dih_table = zeros(" << nrottype << ","
		<< nrottype << "," << nbin_pair << ","
		<< nbin_dih << ");" << std::endl;
	for (Size i=1; i<=rot_type_list.size(); ++i) {
		for (Size j=1; j<=rot_type_list.size(); ++j) {
			for (Size k=1; k<=nbin_pair; k++) {
				octave_scr2 << "dih_table(" << i << "," << j << "," << k << ",:) = [ ";
				for (Size a=1; a<=nbin_dih; a++) {
					octave_scr2 << dih_table[i][j][k][a] << " ";
				}
				octave_scr2 << "];" << std::endl;
			}
		}
	}
}

