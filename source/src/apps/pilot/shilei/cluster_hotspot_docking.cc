// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <protocols/toolbox/KCluster.hh>

#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/id/AtomID.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
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
#include <boost/lexical_cast.hpp>
#include <sstream>
#include <iomanip>

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::toolbox;
using namespace core::io::silent;
//using namespace protocols::moves::mc_convergence_checks;

OPT_1GRP_KEY(Real,cluster_hotspot_docking,cluster_radius)
OPT_1GRP_KEY(Integer,cluster_hotspot_docking,max_clusters)
OPT_1GRP_KEY(Integer,cluster_hotspot_docking,output_ddg_clusters)
OPT_1GRP_KEY(Boolean,cluster_hotspot_docking,heavyatom)
OPT_1GRP_KEY(Boolean,cluster_hotspot_docking,ca)
OPT_1GRP_KEY(String,cluster_hotspot_docking,prefix)
OPT_1GRP_KEY(String,cluster_hotspot_docking,column)

static basic::Tracer TR( "cluster_hotspot_docking" );

core::Real aa2sim_sc(core::pose::Pose pose1, core::pose::Pose pose2) {
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	core::Real rms(0.0);
	runtime_assert( pose1.conformation().size() == 1 );
	runtime_assert( pose2.conformation().size() == 1 );
	rms=core::scoring::all_scatom_rmsd_nosuper(pose1,pose2);
	return rms;
}


core::Real aa2sim_all(core::pose::Pose pose1, core::pose::Pose pose2) {
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	core::Real rms(0.0);
	runtime_assert( pose1.conformation().size() == 1 );
	runtime_assert( pose2.conformation().size() == 1 );
	rms=core::scoring::all_atom_rmsd_nosuper(pose1,pose2);
	return rms;
}

core::Real aa2sim_ca(core::pose::Pose pose1, core::pose::Pose pose2) {
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	core::Real rms(0.0);
	std::map< core::id::AtomID, core::id::AtomID > CA_map;
	core::scoring::setup_matching_CA_atoms(pose1, pose2, CA_map );
	rms=core::scoring::rms_at_corresponding_atoms_no_super( pose1, pose2, CA_map );
	return rms;
}

int main(int argc, char *argv[])
{
	try {

		NEW_OPT(cluster_hotspot_docking::cluster_radius, "radius for cluster",1.0);
		NEW_OPT(cluster_hotspot_docking::max_clusters, "radius for cluster",10);
		NEW_OPT(cluster_hotspot_docking::output_ddg_clusters, "output the lowest ddg N of the clusters",5);
		NEW_OPT(cluster_hotspot_docking::heavyatom, "use all heavy atom for clustering",false);
		NEW_OPT(cluster_hotspot_docking::ca, "use CA atom for clustering (proteins)",false);
		NEW_OPT(cluster_hotspot_docking::prefix, "Prefix for output","");
		NEW_OPT(cluster_hotspot_docking::column, "Score column for clustering","ddg");
		devel::init(argc, argv);

		//read in silent file
		std::string silentin = option[ in::file::silent ]()[1];
		SilentFileOptions opts; // initialized from the command line
		SilentFileData sfd(opts), sfd2(opts);
		sfd.read_file( silentin );

		//Get access to individual pose in the silent file
		utility::vector1< std::string > tags = sfd.tags();

		pose::Pose pose;
		std::vector< core::pose::Pose > poselist;
		std::vector< core::Real > ddglist;
		std::vector< std::string > taglist;
		//std::vector< std::string > outtaglist;
		utility::vector1< std::string > outtaglist;
		std::string scorecolumn=basic::options::option[ basic::options::OptionKeys::cluster_hotspot_docking::column];

		core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
		Real score;

		for ( core::Size ii = 1; ii <= tags.size(); ii++ ) {
			SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_in( opts );
			ss = sfd[ tags[ ii ] ];
			//save score column from the silent files

			ss->fill_pose(pose);

			//control over the columns for ranking
			//getPoseExtraScore( pose, "ddg", score);
			if ( getPoseExtraScore( pose, scorecolumn, score) ) {
				TR.Info << "Pose: "<< scorecolumn << ": " << score << std::endl;
			} else if ( scorecolumn == "Isc" ) {
				score=protocols::docking::calc_interaction_energy(pose,sfxn,utility::tools::make_vector1<core::Size>(1));
				TR.Info << "Compute " << scorecolumn << ": " << score << std::endl;
			} else {
				utility_exit_with_message("scorecolumn not in silent header, only those in header or Isc is supported");
			}

			//score = (*sfxn)(pose) ;
			// remove the chain one from the pose

			//get chain 1 sequence from pose and remove chain 1
			Size inputchain=1;
			Size chainstart(pose.conformation().chain_begin( inputchain ));
			Size chainend(pose.conformation().chain_end( inputchain ) );
			pose.conformation().delete_residue_range_slow( chainstart, chainend);
			TR.Info << "PoseSize: "<< pose.conformation().size() << " " << scorecolumn<< "Score: " << score << " description "<< tags[ ii ] << std::endl;

			poselist.push_back( pose );
			ddglist.push_back( score);
			taglist.push_back( tags[ ii ] );
		}
		TR.Info << "poselist.size(): "<< poselist.size() << std::endl;

		//compute similarity matrix of chain 2 and use sidechain allatom or all heavy atom (including bb)
		ObjexxFCL::FArray2D< core::Real > sc_matrix( poselist.size(),poselist.size(), 0.0 );
		if ( basic::options::option[ basic::options::OptionKeys::cluster_hotspot_docking::heavyatom] ) {
			for ( Size i = 0; i < poselist.size(); i++ ) {
				for ( Size j = i+1 ; j < poselist.size(); j++ ) {
					sc_matrix(i+1,j+1)=aa2sim_all(poselist[i],poselist[j]);
				}
			}
		} else if ( basic::options::option[ basic::options::OptionKeys::cluster_hotspot_docking::ca] ) {
			for ( Size i = 0; i < poselist.size(); i++ ) {
				for ( Size j = i+1 ; j < poselist.size(); j++ ) {
					sc_matrix(i+1,j+1)=aa2sim_ca(poselist[i],poselist[j]);
				}
			}
		} else {
			for ( Size i = 0; i < poselist.size(); i++ ) {
				for ( Size j = i+1 ; j < poselist.size(); j++ ) {
					sc_matrix(i+1,j+1)=aa2sim_sc(poselist[i],poselist[j]);
				}
			}
		}


		//assign the other half of the matrix
		for ( Size i = 0; i < poselist.size(); i++ ) {
			for ( Size j = 0 ; j < i; j++ ) {
				sc_matrix(i+1,j+1)=sc_matrix(j+1,i+1);
			}
		}

		//print similarity matrix
		/*
		for ( Size i = 0; i < poselist.size(); i++ ) {
		for ( Size j = 0 ; j < poselist.size(); j++ ) {
		TR.Info << sc_matrix(i+1,j+1) << " ";
		}
		TR.Info << std::endl;
		}
		*/

		//cluster similarity matrix
		//what algorithms will you use?
		//the one in Rosetta
		core::Real cluster_radius_=basic::options::option[ basic::options::OptionKeys::cluster_hotspot_docking::cluster_radius];
		core::Size max_total_cluster=basic::options::option[ basic::options::OptionKeys::cluster_hotspot_docking::max_clusters];
		core::Size output_ddg_clusters=basic::options::option[ basic::options::OptionKeys::cluster_hotspot_docking::output_ddg_clusters];
		std::string prefix=basic::options::option[ basic::options::OptionKeys::cluster_hotspot_docking::prefix];
		core::Size listsize=poselist.size();
		std::vector < int > neighbors ( poselist.size(), 0 );
		std::vector < int > clusternr ( poselist.size(), -1 );
		std::vector < int > clustercenter;
		core::Size mostneighbors;
		core::Size nclusters=0;
		core::Size i,j;
		std::vector <int> clustercentre;

		TR.Info << "Clustering of " << listsize << " structures with radius " <<  cluster_radius_ << " with max of  " << max_total_cluster << " cluster "<< " and output lowest " <<  output_ddg_clusters << " " << scorecolumn << " cluster centers"<<  std::endl;

		// now assign groupings
		while ( true ) {
			// count each's neighbors
			for ( i=0; i<listsize; i++ ) {
				neighbors[i] = 0;
				if ( clusternr[i]>=0 ) continue; // ignore ones already taken
				for ( j=0; j<listsize; j++ ) {
					if ( clusternr[j]>=0 ) continue; // ignore ones already taken
					if ( sc_matrix( i+1, j+1 ) < cluster_radius_ ) neighbors[i]++;
				}
				//TR.Info << "i: " << i << " " << neighbors[i] << std::endl;
			}

			mostneighbors = 0;
			for ( i=0; i<listsize; i++ ) {
				if ( neighbors[i]>neighbors[mostneighbors] ) mostneighbors=i;
			}
			if ( neighbors[mostneighbors] <= 0 ) break;  // finished!

			for ( i=0; i<listsize; i++ ) {
				if ( clusternr[i]>=0 ) continue; // ignore ones already taken
				if ( sc_matrix( i+1, mostneighbors+1 ) < cluster_radius_ ) {
					clusternr[i] = mostneighbors;
				}
			}

			clustercentre.push_back(mostneighbors);
			nclusters++;

			if ( nclusters > max_total_cluster ) break;  // currently fixed but ought to be a paraemter

			//if ((nclusters%10)==0) TR.Info << ".";

			//TR.Info.flush();
		}

		TR.Info << "ncluster: " << nclusters << std::endl;
		core::Real bestddg;
		core::Size bestj;
		core::Size sizei=0;

		std::multimap< core::Real, core::Size> ddgclustermap;

		for ( i=0; i<clustercentre.size(); i++ ) {
			TR.Info << "CLUSTER " << i << " : ";
			bestddg=999;
			bestj=999;
			for ( j=0; j<listsize; j++ ) {
				// if that struture belongs to a given cluster
				if ( clusternr[j] == clustercentre[i] ) {
					if ( ddglist[j] < bestddg ) {
						bestddg=ddglist[j];
						bestj=j;
					}
					sizei++;
					//new_cluster.add_member(j);         // add structure
					//TR.Info << " "<< taglist[j]<<".pdb" << " " << ddglist[j] << " " ;
					TR.Info << " "<< taglist[j] << " ";
				}
			}
			TR.Info << std::endl;

			TR.Info << "best"<<scorecolumn<<": " << bestddg << " from: " << taglist[bestj] << " of Cluster " << i << " Size: " << sizei << std::endl;
			ddgclustermap.insert( std::make_pair(bestddg, bestj) );

			sizei=0;
		}

		core::Size outsize=0;
		for ( std::multimap< core::Real, core::Size>::iterator ddg_index_pair = ddgclustermap.begin(); ddg_index_pair != ddgclustermap.end(); ++ddg_index_pair ) {
			if ( outsize<output_ddg_clusters ) {
				outtaglist.push_back(taglist[ddg_index_pair->second]);
				TR.Info << ddg_index_pair->second<< " : " << taglist[ddg_index_pair->second] << " : " << ddg_index_pair->first << std::endl;
			}
			outsize++;
			if ( outsize>=output_ddg_clusters ) break;
		}

		core::import_pose::pose_stream::PoseInputStreamOP input( new core::import_pose::pose_stream::SilentFilePoseInputStream( option[ in::file::silent ](), outtaglist) );

		while ( input->has_another_pose() ) {
			input->fill_pose( pose );
			if ( getPoseExtraScore( pose, scorecolumn, score) ) {
				TR.Info << "Pose: "<< scorecolumn << ": " << score << std::endl;
			} else if ( scorecolumn == "Isc" ) {
				score=protocols::docking::calc_interaction_energy(pose,sfxn,utility::tools::make_vector1<core::Size>(1));
			} else {
				utility_exit_with_message("scorecolumn not in silent header, only those in header or Isc is supported");
			}
			//getPoseExtraScore( pose, "ddg", score);
			std::string tag( tag_from_pose( pose ) );

			std::ostringstream ss;
			ss << std::fixed << std::setprecision(3);
			ss << score;
			std::string out_prefix = prefix+"_"+scorecolumn+"_"+ss.str()+"_";
			//std::string out_prefix = prefix+"_ddg_"+ss.str()+"_";
			//std::string out_prefix = prefix+"_ddg_"+boost::lexical_cast<std::string>(score)+"_";
			std::string fn( out_prefix + tag + ".pdb" );
			pose.dump_pdb(fn);
		}

		TR.Info << "done cluster_hotspot_docking" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
