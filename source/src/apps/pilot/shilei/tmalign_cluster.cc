
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <protocols/toolbox/KCluster.hh>

#include <core/init/init.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/Residue.hh>
#include <protocols/docking/metrics.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <iostream>

#include <protocols/cluster/cluster.hh>
#include <core/scoring/rms_util.hh>

#include <basic/options/option_macros.hh>
#include <map>
#include <sstream>
#include <iomanip>

//use TM align for sequence independent
#include <protocols/hybridization/TMalign.hh>
#include <protocols/hybridization/util.hh>
#include <boost/lexical_cast.hpp>

//using boost threads
#ifdef USE_BOOST_THREAD
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#endif

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::toolbox;
using namespace core::io::silent;
//using namespace protocols::moves::mc_convergence_checks;

OPT_1GRP_KEY(Real,tmalign_cluster,cluster_radius)

static basic::Tracer TR("tmalign_cluster");

void do_tmscore( core::pose::Pose const &pose1, core::pose::Pose const &pose2, core::Real &tmscore) {
	protocols::hybridization::TMalign tm_align;
	if (pose1.total_residue()==0 && pose2.total_residue()==0)
		tmscore=0;
	else {
		int reval = tm_align.apply(pose1,pose2);
                if (reval == 0) {
                	tmscore=tm_align.TMscore(std::min(pose1.total_residue(),pose2.total_residue()));
		} else {
			tmscore=0;
		}
	}
}

void do_tmalign(core::pose::Pose & aligned_pose,core::pose::Pose const & ref_pose) {
	protocols::hybridization::TMalign tm_align;

        std::string seq_pose, seq_ref, aligned;
        core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::Size n_mapped_residues=0;

        core::pose::initialize_atomid_map( atom_map, aligned_pose, core::id::BOGUS_ATOM_ID );

	std::list <core::Size> residue_list;
    	std::list <core::Size> ref_residue_list;
    	for ( Size ires=1; ires<= aligned_pose.total_residue(); ++ires ) {
      		if ( !aligned_pose.residue(ires).is_protein() ) continue;
      			residue_list.push_back(ires);
    	}
    	for ( Size ires=1; ires<= ref_pose.total_residue(); ++ires ) {
      		if ( !ref_pose.residue(ires).is_protein() ) continue;
      			ref_residue_list.push_back(ires);
    	}

        tm_align.apply(aligned_pose, ref_pose, residue_list, ref_residue_list);
        tm_align.alignment2AtomMap(aligned_pose, ref_pose, residue_list, ref_residue_list, n_mapped_residues, atom_map);
        tm_align.alignment2strings(seq_pose, seq_ref, aligned);

        std::list <Size> full_residue_list;
        for (Size ires=1; ires<=aligned_pose.total_residue(); ++ires) {
          full_residue_list.push_back(ires);
        }

  	if (n_mapped_residues >= 6) {
                utility::vector1< core::Real > aln_cutoffs;
                aln_cutoffs.push_back(6);
                aln_cutoffs.push_back(4);
                aln_cutoffs.push_back(3);
                aln_cutoffs.push_back(2);
                aln_cutoffs.push_back(1.5);
                aln_cutoffs.push_back(1);
                core::Real min_coverage = 0.2;
 	        protocols::hybridization::partial_align(aligned_pose, ref_pose, atom_map, full_residue_list, true, aln_cutoffs, min_coverage);
	}
}

int main(int argc, char *argv[])
{
  try {
	NEW_OPT(tmalign_cluster::cluster_radius, "radius for cluster",0.8);

	core::init::init(argc, argv);

	//read a list of pdbs
	Size input_pdb_size=basic::options::option[ basic::options::OptionKeys::in::file::s ]().size();
	assert(input_pdb_size>1);
	utility::vector1< core::pose::PoseOP > PoseVec;
	core::pose::PoseOP pose;
        for ( Size i = 1; i <= input_pdb_size; ++i ) {
		std::string pdbnamei=basic::options::option[ basic::options::OptionKeys::in::file::s ]()[i];
//		TR << pdbnamei << endl;
		pose=new core::pose::Pose;
		core::import_pose::pose_from_pdb( *pose, pdbnamei.c_str() );
		PoseVec.push_back(pose);
	}

	//std::vector< core::Real > ddglist;
	//std::vector< std::string > taglist;
	//std::vector< std::string > outtaglist;
	//utility::vector1< std::string > outtaglist;
        //std::string scorecolumn=basic::options::option[ basic::options::OptionKeys::tmalign_cluster::column];
        //core::scoring::ScoreFunctionOP sfxn = core::scoring::getScoreFunction();
	//Real score;

	TR.Info << "PoseVec.size(): "<< PoseVec.size() << std::endl;

	//compute similarity matrix
	//use boost threads frank/rr_opt.cc
	//protocols/frag_picker/FragmentPicker.cc
        ObjexxFCL::FArray2D< core::Real > sc_matrix( PoseVec.size(),PoseVec.size(), 1.0 );
	#ifdef USE_BOOST_THREAD
		TR << " Running multi-thread version " << std::endl;
	boost::thread_group threads;
	for ( Size i = 1; i <= PoseVec.size(); ++i ) {
    		for ( Size j = i ; j <= PoseVec.size(); ++j ) {
			TR << "i " << i << " j " << j << " ";
			//boost::thread*  threadij = new boost::thread( &do_tmscore, boost::cref(*PoseVec[i]), boost::cref(*PoseVec[j]), boost::ref(sc_matrix(i,j)));
			//threads.add_thread(threadij);
			threads.create_thread( boost::bind( &do_tmscore, boost::cref(*PoseVec[i]), boost::cref(*PoseVec[j]), boost::ref(sc_matrix(i,j)) )) ;
			TR << "thread created"<< std::endl;
		}
	}
        threads.join_all();
	#else
	//TR << " use extras=boost_thread to compile boost version, seg fault, why? " << endl;
        for ( Size i = 1; i <= PoseVec.size(); ++i ) {
                for ( Size j = i+1 ; j <= PoseVec.size(); ++j ) {
                        do_tmscore(*PoseVec[i], *PoseVec[j], sc_matrix(i,j));
						TR << "TMscore("<<i<<","<<j<<"): " << sc_matrix(i,j) << std::endl;
                }
        }
	//std::cerr << "compile with extras=boost_thread-mt!" << std::endl;
	#endif

	/*
	//print similarity matrix
	for ( Size i = 1; i <= PoseVec.size(); ++i ) {
	    for ( Size j = 1 ; j <= PoseVec.size(); ++j ) {
			TR.Info << sc_matrix(i,j) << " ";
		}
		TR.Info << std::endl;
	}
	*/

	//assign the other half of the matrix
	for ( Size i = 1; i <= PoseVec.size(); ++i ) {
	    for ( Size j = 1 ; j < i; ++j ) {
			sc_matrix(i,j)=sc_matrix(j,i);
		}
	}

	//print similarity matrix
	for ( Size i = 1; i <= PoseVec.size(); ++i ) {
	    for ( Size j = 1 ; j <= PoseVec.size(); ++j ) {
			TR.Info << sc_matrix(i,j) << " ";
		}
		TR.Info << std::endl;
	}

	//cluster similarity matrix
	//what algorithms will you use?
	//the one in Rosetta
	core::Real cluster_radius_=basic::options::option[ basic::options::OptionKeys::tmalign_cluster::cluster_radius];
	core::Size listsize=PoseVec.size();
	std::vector < int > neighbors ( PoseVec.size(), 0 );
	std::vector < int > clusternr ( PoseVec.size(), -1 );
	std::vector < int > clustercenter;
	core::Size mostneighbors;
	core::Size nclusters=0;
	core::Size i,j;
	std::vector <int> clustercentre;

        TR.Info << "Clustering of " << listsize << " structures with radius (TMscore) " <<  cluster_radius_ <<  std::endl;

  	// now assign groupings
  	while(true) {
    	// count each's neighbors
        for (i=0;i<listsize;++i ) {
          neighbors[i] = 0;
          if (clusternr[i]>=0) continue; // ignore ones already taken
          for (j=0;j<listsize;++j ) {
            if (clusternr[j]>=0) continue; // ignore ones already taken
	    //TR.Info << "sc_matrix("<<i+1<<","<<j+1<<") = " << sc_matrix( i+1, j+1 ) << std::endl;
            if ( sc_matrix( i+1, j+1 ) > cluster_radius_ ) neighbors[i]++;
          }
	  TR.Info << "i: " << i << " " << neighbors[i] << std::endl;
        }

        mostneighbors = 0;
        for (i=0;i<listsize;++i ) {
          if (neighbors[i]>neighbors[mostneighbors]) {
		mostneighbors=i;
		}
        }
	TR.Info << "mostneighbors =" << mostneighbors <<std::endl;
        if (neighbors[mostneighbors] <= 0) break;  // finished!

	for (i=0;i<listsize;++i ) {
      		if (clusternr[i]>=0) continue; // ignore ones already taken
      		if ( sc_matrix( i+1, mostneighbors+1 ) > cluster_radius_) {
        		clusternr[i] = mostneighbors;
      		}
    	}

    	clustercentre.push_back(mostneighbors);
    	nclusters++;

  	}

	TR.Info << "ncluster: " << nclusters << std::endl;
	for (i=0;i<clustercentre.size();++i ) {
		TR.Info << "CLUSTER " << i << ", " << clustercentre[i] <<" : ";
    		for (j=0;j<listsize;++j ) {
      			if (clusternr[j] == clustercentre[i]) {
				TR.Info << j << " " ;
				//std::string tag( tag_from_pose( *PoseVec[j+1] ) );
				std::string tag=basic::options::option[ basic::options::OptionKeys::in::file::s ]()[j+1];
				std::string fn( "c_"+boost::lexical_cast<std::string>(i)+"_" + boost::lexical_cast<std::string>(j) +"_" + tag );
				do_tmalign(*PoseVec[j+1],*PoseVec[clustercentre[i]+1]);
				PoseVec[j+1]->dump_pdb(fn);
			}
		}
		TR.Info << std::endl;
	}

	TR.Info << "done tmalign_cluster" << std::endl;

  } catch ( utility::excn::EXCN_Base const & e ) {
    std::cout << "caught exception " << e.msg() << std::endl;
    return -1;
  }
}
