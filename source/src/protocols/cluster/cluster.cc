// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
///
/// @author Mike Tyka


#include <core/conformation/Residue.hh>
#include <core/pose/symmetry/util.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/cluster/cluster.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/moves/Mover.hh>
#include <utility/file/FileName.hh>
#include <utility/string_util.hh>
#include <core/scoring/constraints/util.hh>
#include <protocols/loops/loops_main.hh>
#include <utility/io/izstream.hh>
#include <numeric/random/random.hh>

// C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <deque>
#include <map>
#include <cmath>

// option key includes

#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>


//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end

namespace protocols {
namespace cluster {

using namespace core;
using namespace ObjexxFCL;
using namespace pose;
using namespace basic::options;
using namespace evaluation;

static THREAD_LOCAL basic::Tracer tr( "protocols.cluster" );


std::map< std::string, core::Real > read_template_scores( const std::string & filename ){

	utility::io::izstream data( filename.c_str() );
	if ( !data.good() ) {
		utility_exit_with_message(
			"ERROR: Unable to open extra_scores file: '" + filename + "'"
		);
	}

	std::map< std::string, core::Real > score_data;
	std::string line;//,tag;
	while ( getline(data,line) ) {
		std::istringstream l( line );
		std::string template_name;
		core::Real  template_score;
		l >> template_name;
		l >> template_score;

		score_data[template_name] = template_score;
	}

	return score_data;
}


bool compareIndexEnergyPair(
	const std::pair< int, Real > & p1, const std::pair< int, Real > & p2 )
{
	return p1.second  < p2.second;
}

/*
class PoseComparator {
public:
PoseComparator() { }
virtual Real measure( Pose &pose1, Pose &pose2 ) = 0;
private:
};


class PoseComparator_RMSD : public PoseComparator {
public:
PoseComparator() { }
ual Real measure( Pose &pose1, Pose &pose2 ) = 0;
private:
};

*/

GatherPosesMover::GatherPosesMover() : Mover()
{
	filter_ = 1000000.0;
	cluster_by_protein_backbone_ = false;
	cluster_by_all_atom_ = false;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ OptionKeys::cluster::template_scores].user() ) {
		template_scores = read_template_scores( option[ OptionKeys::cluster::template_scores]() );

		tr.Info << "Read template scores: " << std::endl;
		for ( std::map<std::string, core::Real >::iterator ii=template_scores.begin(); ii!=template_scores.end(); ++ii ) {
			tr.Info << (*ii).first << ": " << (*ii).second << std::endl;
		}
	}

}

// This should probably be replaced with an object -- "PosePoseRMSD" or something.
Real
GatherPosesMover::get_distance_measure(
	const Pose & pose1,
	const Pose & pose2
) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ basic::options::OptionKeys::symmetry::symmetric_rmsd ]() &&
			!core::pose::symmetry::is_symmetric( pose1 ) ) {
		tr.Info << "Warning!!! For symmetric rmsd selected but pose is not symmetric. Ignoring symmetry" << std::endl;
	}

	if ( option[ OptionKeys::cluster::hotspot_hash ]() ) {
		Size const resnum1 = pose1.total_residue();
		Size const resnum2 = pose2.total_residue();
		conformation::ResidueCOP res1( pose1.residue(resnum1).get_self_ptr() );
		conformation::ResidueCOP res2( pose2.residue(resnum2).get_self_ptr() );
		//return protocols::hotspot_hashing::residue_sc_rmsd_no_super( res1, res2 );
		return core::scoring::residue_sc_rmsd_no_super( res1, res2 );
	}

	if ( option[ OptionKeys::cluster::gdtmm ]() ) {
		if ( option[ basic::options::OptionKeys::symmetry::symmetric_rmsd ]() ) utility_exit_with_message( "No symmetric gdtmm available!!!!\n" ) ;

		core::scoring::ResidueSelection residues;
		core::scoring::invert_exclude_residues( pose1.total_residue(), option[ OptionKeys::cluster::exclude_res ](), residues );
		// Make a new measure that is like an inverse gdtmm running from 10 (=gdtm of 0) to
		// 0 (= gdtmm of 1)
		return 10 -10.0 * scoring::CA_gdtmm( pose1, pose2, residues );
	}

	if ( option[ OptionKeys::cluster::exclude_res ].user() ) {
		core::scoring::ResidueSelection residues;
		core::scoring::invert_exclude_residues( pose1.total_residue(), option[ OptionKeys::cluster::exclude_res ](), residues );
		if ( pose1.residue(1).is_RNA() ) utility_exit_with_message( "Hey put in all atom rmsd code for residue subset!\n" ) ;
		if ( option[ basic::options::OptionKeys::symmetry::symmetric_rmsd ]() &&
				core::pose::symmetry::is_symmetric( pose1 ) ) {
			tr.Info << "Warning!!! For symmetric clustering currently only all CA clustering is available. Calculating all CA rmsd instead..." << std::endl;
			return scoring::CA_rmsd_symmetric( pose1, pose2 );
		}
		return scoring::CA_rmsd( pose1, pose2, residues );
	} else {
		// no residues excluded from the native.
		if ( option[ basic::options::OptionKeys::symmetry::symmetric_rmsd ]() &&
				core::pose::symmetry::is_symmetric( pose1 ) ) return scoring::CA_rmsd_symmetric( pose1, pose2 );
		if ( option[ basic::options::OptionKeys::cluster::skip_align ].user() && pose1.residue(1).is_RNA() ) return scoring::all_atom_rmsd_nosuper( pose1, pose2 );
		if ( option[ basic::options::OptionKeys::cluster::skip_align ].user() ) return scoring::rmsd_no_super( pose1, pose2, scoring::is_protein_backbone );
		if ( pose1.residue(1).is_RNA() ) return scoring::all_atom_rmsd( pose1, pose2 );
		if ( cluster_by_all_atom_ ) return scoring::all_atom_rmsd( pose1, pose2 );
		if ( cluster_by_protein_backbone_ ) return scoring::rmsd_with_super( pose1, pose2, scoring::is_protein_backbone );
		if ( option[ basic::options::OptionKeys::cluster::skip_align ].user() ) return scoring::rmsd_no_super( pose1, pose2, scoring::is_protein_backbone );
		return scoring::CA_rmsd( pose1, pose2 );
	}

} // native_CA_rmsd


void GatherPosesMover::set_score_function(  scoring::ScoreFunctionOP sfxn ) {
	sfxn_ = sfxn;
}

void GatherPosesMover::set_filter( Real filter ) {
	filter_ = filter;
}

core::Real GatherPosesMover::get_filter() const {
	return filter_;
}

void GatherPosesMover::apply( Pose & pose ) {

	core::Real score(0.0);
	// does pose already have score set ?
	if ( !getPoseExtraScore( pose, "silent_score", score ) ) {
		// if not do we have scoring function ?
		// (Should we be falling back to score or total_score as well?)

		if ( !sfxn_ ) {
			//oh well - can't do anything ... just leave it at 0.
		} else {
			// if so - try and score the structure and set that energy column
			////////////////// add constraints if specified by user.
			scoring::constraints::add_constraints_from_cmdline( pose, *sfxn_ );
			score = (*sfxn_)(pose);
			tr.Info << "RESCORING: " << core::pose::extract_tag_from_pose( pose )<< std::endl;
			setPoseExtraScore( pose, "silent_score", score );

			if ( get_native_pose() ) {
				setPoseExtraScore(
					pose, "rms", get_distance_measure( pose, *get_native_pose() )
				);
			}
			pose.energies().clear();
		}

		// remember tag!
		// RM: Should this be outside the getPoseExtraScore block?
		tag_list.push_back( core::pose::extract_tag_from_pose( pose ) );
	}

	if ( template_scores.size() > 0 ) {

		std::map< std::string, std::string > score_line_strings( core::pose::get_all_score_line_strings( pose ) );

		std::string template_name = score_line_strings["aln_id"];
		template_name = template_name.substr(0,5);

		std::cout << "Found template name: " << template_name << std::endl;
		core::Real template_score = template_scores[ template_name ];

		score += template_score;
		setPoseExtraScore( pose, "silent_score", score );
	}

	// filter structures if the filter is active
	if ( score > filter_ ) {
		tr.Info << "Ignoring structure " << core::pose::extract_tag_from_pose( pose ) << " because score " << score << " is greater than threshold " << filter_ << std::endl;
		return; // ignore pose if filter value is too high
	}

	tr.Info << "Adding struc: " << score << std::endl;

	// now save the structure.
	poselist.push_back( pose );
	//tr.Info << "Read " << poselist.size() << " structures" << std::endl;
}

std::string
GatherPosesMover::get_name() const {
	return "GatherPosesMover";
}

bool GatherPosesMover::check_tag( const std::string &query_tag ) {
	using std::find;
	if ( find( tag_list.begin(), tag_list.end(), query_tag ) == tag_list.end() ) {
		return false;
	}
	return true;
}

ClusterBase::ClusterBase()
: GatherPosesMover(),
	export_only_low_(false),
	median_rms_( 0.0 ),
	population_weight_( 0.09 ),
	cluster_radius_( 2.0 )
{}

std::string
ClusterBase::get_name() const {
	return "ClusterBase";
}

void ClusterBase::set_cluster_radius( Real cluster_radius ) {
	cluster_radius_ = cluster_radius;
}

void ClusterBase::set_population_weight( Real population_weight ) {
	population_weight_ = population_weight ;
}

Real ClusterBase::get_cluster_radius() {
	return cluster_radius_;
}

void ClusterBase::calculate_distance_matrix() {
	FArray2D< Real > p1a, p2a;
	distance_matrix = FArray2D< Real > ( poselist.size(), poselist.size(), 0.0 );
	int count = 0;
	tr.Info << "Calculating RMS matrix: " << std::endl;

	Real hist_resolution = 0.25;
	int   hist_size = int(20.0 / hist_resolution);

	std::vector<int> histcount(hist_size,0);

	for ( Size i = 0; i < poselist.size(); i++ ) {
		for ( Size j = i; j < poselist.size(); j++ ) {
			// only do comparisons once
			if ( i < j ) {
				// get the similarity between structure with index i and with index j
				Real dist  = get_distance_measure( poselist[i],poselist[j]);
				distance_matrix( i+1, j+1 ) = dist;
				distance_matrix( j+1, i+1 ) = dist;
				//    tr.Info << "( " << i+1 << ";" << j+1 << " ) " << dist << std::endl;
				int histbin = int(dist/hist_resolution);
				if ( histbin < hist_size ) histcount[histbin]+=1;

				// print some stats of progress
				count ++;
				if ( count % 5000 == 0 ) {
					Real const percent_done ( 200.0 * static_cast< Real > ( count ) / ( (poselist.size() - 1) * poselist.size()  ) );
					tr.Info << count
						<< "/" << ( poselist.size() - 1 )*( poselist.size() )/2
						<< " ( " << F(8,1,percent_done) << "% )"
						<< std::endl;
				}
			}
		} // for it2
	} // for it1


	tr.Info << "Histogram of pairwise similarity values for the initial clustering set" << std::endl;
	int maxcount_count = 0;
	int maxcount_i = 0;
	for ( Size i = 0; i <(Size)hist_size; i++ ) {
		tr.Info << "hist " <<  Real(i)*hist_resolution << "   " << histcount[i] << std::endl;
		if ( histcount[i] > maxcount_count ) {
			maxcount_count = histcount[i];
			maxcount_i = i;
		}
	}

	median_rms_ = maxcount_i * hist_resolution;
	tr.Info << "Median RMS " << get_median_rms() << std::endl;

} // calculate_distance_matrix

// PostProcessing ---------------------------------------------------------
void ClusterBase::add_structure( Pose & pose ) {

	poselist.push_back( pose );
	int nexindex = poselist.size() - 1;
	int lowrmsi=-1;
	Real lowrms=1000000.0;
	for ( int m=0; m<(int)clusterlist.size(); m++ ) {
		Real rms;
		rms =   get_distance_measure(
			poselist[ nexindex ],
			poselist[ clusterlist[m].get_cluster_center() ]  ); // clustercentre m
		if ( rms < lowrms ) {
			lowrms  = rms;
			lowrmsi = m;
		}
	}

	if ( lowrms <= 0.001 ) {
		tr.Info << "Structure identical to existing structure - ignoring" << std::endl;
		return;
	}

	if ( lowrmsi >= 0 ) {
		if ( lowrms < get_cluster_radius() ) {  // if within our radius - then add to cluster
			clusterlist[lowrmsi].add_member(nexindex);
			tr.Info << "Adding to cluster " << lowrmsi << " Cluster_rad: " << get_cluster_radius() << std::endl;
		} else {                                 // else make a new cluster !!
			Cluster new_cluster( nexindex );
			clusterlist.push_back( new_cluster );
			tr.Info << "Adding as new cluster " << " Cluster_rad: " << get_cluster_radius() << std::endl;
		}
	}
}

void Cluster::shuffle(){
	// fisher-yates
	for ( int i=member.size()-1; i>-1; i-- ) {
		int j = numeric::random::rg().random_range(0, 100000000) % (i + 1);
		if ( i != j ) {
			std::swap(member[j], member[i]);
		}
	}
}
// PostProcessing ---------------------------------------------------------

void ClusterBase::sort_each_group_by_energy( ) {
	tr.Info << "Sorting each cluster's structures by energy: " << std::endl;
	int i,j;
	for ( i=0; i<(int)clusterlist.size(); i++ ) {

		std::vector< std::pair< int, Real > > cluster_energies;
		for ( j=0; j<(int)clusterlist[i].size(); j++ ) {
			Real score=0.0;
			if ( !getPoseExtraScore( poselist[ clusterlist[i][j] ], "silent_score", score ) ) {
				tr.Error << "Warning: no score available for " << std::endl;
			}
			cluster_energies.push_back( std::pair< int, Real > ( clusterlist[i][j], score ) );
		}

		std::sort( cluster_energies.begin(), cluster_energies.end(), compareIndexEnergyPair );
		clusterlist[i].clear();  // This retains the cluster center!
		for ( j=0; j<(int)cluster_energies.size(); j++ ) clusterlist[i].push_back( cluster_energies[j].first );

	}

}

// PostProcessing ---------------------------------------------------------

void ClusterBase::remove_highest_energy_member_of_each_group() {
	int i,j;
	for ( i=0; i<(int)clusterlist.size(); i++ ) {
		if ( clusterlist[i].size() <= 1 ) continue;

		// sort group by energy
		std::vector< std::pair< int, Real > > cluster_energies;
		for ( j=0; j<(int)clusterlist[i].size(); j++ ) {
			Real score=0.0;
			if ( !getPoseExtraScore( poselist[ clusterlist[i][j] ], "silent_score", score ) ) {
				tr.Error << "Warning: no score available for " << std::endl;
			}

			cluster_energies.push_back( std::pair< int, Real > ( clusterlist[i][j], score ) );
		}
		std::sort( cluster_energies.begin(), cluster_energies.end(), compareIndexEnergyPair );
		clusterlist[i].clear();
		for ( j=0; j<(int)(cluster_energies.size()-1); j++ ) clusterlist[i].push_back( cluster_energies[j].first );
	}

}

void ClusterBase::sort_groups_by_energy() {
	tr.Info << "Sorting clsuters by energy: " << std::endl;

	int i;
	std::vector < Cluster >   temp = clusterlist;

	std::vector< std::pair< int, Real > > cluster_energies;
	for ( i=0; i<(int)clusterlist.size(); i++ ) {
		Real score;

		// This assumes that the first member is already the lowest energy one!
		int energy_member = 0;
		if ( !getPoseExtraScore( poselist[ clusterlist[i][energy_member] ], "silent_score", score ) ) {
			tr.Error << "Warning: no score available for " << std::endl;
		}

		Real combo_score = (1.0 - population_weight_ ) * score -  population_weight_  * clusterlist[i].group_size() ;  // group size keeps trakc of the size of the clsuter, even with pruning
		cluster_energies.push_back( std::pair< int, Real > ( i, combo_score ) );
	}

	std::sort( cluster_energies.begin(), cluster_energies.end(), compareIndexEnergyPair );

	// clear the cluster list
	clusterlist.clear();

	// now put back every cluster in the order of energies
	for ( i=0; i<(int)cluster_energies.size(); i++ ) clusterlist.push_back( temp[cluster_energies[i].first ] );

}

void ClusterBase::remove_singletons() {
	int i;
	std::vector < Cluster >   clusterlist_copy = clusterlist;
	clusterlist.clear();

	//std::vector< std::pair< int, Real > > cluster_energies;
	for ( i=0; i<(int)clusterlist_copy.size(); i++ ) {
		if ( clusterlist_copy[i].group_size() > 1 ) {
			clusterlist.push_back( clusterlist_copy[i] );
		}
	}

	//In case nothing is left, better at least return one cluster?
	if ( clusterlist.size() < 1 ) {
		clusterlist.push_back( clusterlist_copy[ 0 ] );
	}


}

// take first 'limit' from each cluster
void ClusterBase::limit_groupsize( int limit ) {
	tr.Info << "Limiting each cluster to a total size of : " << limit << std::endl;
	int i,j;
	if ( limit < 1 ) return;
	for ( i=0; i<(int)clusterlist.size(); i++ ) {
		Cluster temp = clusterlist[i];
		clusterlist[i].clear();
		for ( j=0; (j<(int)temp.size()) && (j<limit); j++ ) {
			clusterlist[i].push_back( temp[j] );
		}
	}
}

// take first 'limit'%  from each cluster
void ClusterBase::limit_groupsize( core::Real percent_limit ) {
	int i,j;//, limit;
	if ( percent_limit > 1.0 ) return;
	for ( i=0; i<(int)clusterlist.size(); i++ ) {
		Cluster temp = clusterlist[i];
		clusterlist[i].clear();
		int limit = static_cast<core::Size>( std::floor(percent_limit*temp.size()) );
		tr << "truncating from " << temp.size() << " to " << limit << std::endl;
		for ( j=0; j<limit; j++ ) {
			clusterlist[i].push_back( temp[j] );
		}
	}
}

// take random 'limit'% from each cluster
void ClusterBase::random_limit_groupsize( core::Real percent_limit ) {
	int i,j;//,limit;
	if ( percent_limit >= 1.0 ) return;
	for ( i=0; i<(int)clusterlist.size(); i++ ) {
		Cluster temp = clusterlist[i];
		clusterlist[i].clear();
		int limit = static_cast<core::Size> (std::floor(percent_limit*temp.size()) );
		tr << "truncating from " << temp.size() << " to " << limit << std::endl;
		temp.shuffle();
		for ( j=0; j<limit; j++ ) {
			clusterlist[i].push_back( temp[j] );
		}
	}
}

// take first 'limit' clusters
void ClusterBase::limit_groups( int limit ) {
	tr.Info << "Limiting total number of clusters to : " << limit << std::endl;
	int i;
	if ( limit < 0 ) return;
	std::vector < Cluster >   temp = clusterlist;
	clusterlist.clear();
	for ( i=0; i<(int)limit; i++ ) {
		if ( (int)i >= (int)temp.size() ) break;
		clusterlist.push_back( temp[i] );
	}

	// Remove the Poses from the removed clusters!
	for ( i=limit; i<(int)temp.size() ; i++ ) {
		for ( int j=0; j<(int)temp[i].size(); j++ ) {
			poselist[ temp[i].member[j] ] = Pose();
		}
	}
}

// take first 'limit' total_structures
void ClusterBase::limit_total_structures( int limit) {
	tr.Info << "Limiting total structure count to : " << limit << std::endl;
	int i,j;
	if ( limit < 0 ) return;
	std::vector < Cluster >   temp = clusterlist;
	clusterlist.clear();
	int count=0;
	for ( i=0; i<(int)temp.size(); i++ ) {
		Cluster newcluster( temp[i].get_cluster_center() );
		newcluster.clear();
		for ( j=0; j<(int)temp[i].size(); j++ ) {
			if ( count < limit ) {
				newcluster.push_back( temp[i][j] );
			}
			count ++;
		}
		if ( newcluster.size() > 0 ) {
			clusterlist.push_back( newcluster );
		}
	}
}


void ClusterBase::clean_pose_store() {
	tr.Info << "Cleaning Pose store to save memory ... " << std::endl;
	// for each structure in the pose list
	for ( Size index=0; index < poselist.size(); index ++ ) {
		bool ispresent = false;
		// try and find it in the cluster assigments
		Size i,j;
		for ( i=0; i<clusterlist.size(); i++ ) {
			for ( j=0; j<clusterlist[i].size(); j++ ) {
				if ( (int)index == clusterlist[i][j] ) { ispresent = true; break; }
			}
			if ( (int)index == clusterlist[i].get_cluster_center() )  ispresent = true;
			if ( ispresent ) break;
		}
		if ( !ispresent ) {
			// zap the pose (but retain a DUMMY POSE so that the indexing stays the smae !!! )
			// i know this is retarded, PoseOP would be better at least.
			poselist[index] = Pose() ;
		}
	}
}


//  ------ OUTPUT -----------------------------------------------------------------------

void ClusterBase::print_summary() {
	tr.Info << "---------- Summary ---------------------------------" << std::endl;
	int i,j;

	int count=0;

	for ( i=0; i<(int)clusterlist.size(); i++ ) {
		tr.Info << "Cluster:  " << i << "  N: " << (int)clusterlist[i].size() <<  "GN: " <<  (int)clusterlist[i].group_size() << "   c." << i <<".*.pdb " ;
		tr.Info << std::endl;
		count += clusterlist[i].size();
		for ( j=0; j<(int)clusterlist[i].size(); j++ ) {
			tr.Info << "    ";
			tr.Info << core::pose::extract_tag_from_pose( poselist[ clusterlist[i][j] ]  ) << "  " ;
			Real score = 0.0;
			if ( !getPoseExtraScore( poselist[ clusterlist[i][j] ], "silent_score", score ) ) {
				tr.Info << "----" ;
			} else {
				tr.Info << score ;
			}

			tr.Info << std::endl;
		}
	}
	tr.Info << "----------------------------------------------------" << std::endl;
	tr.Info << "  Clusters: " << clusterlist.size() << std::endl;
	tr.Info << "  Structures: " << count << std::endl;
	tr.Info << "----------------------------------------------------" << std::endl;
}

void ClusterBase::print_raw_numbers() {
	int i,j;
	for ( i=0; i<(int)clusterlist.size(); i++ ) {
		tr.Info << i << " : ";
		for ( j=0; j<(int)clusterlist[i].size(); j++ ) {
			tr.Info << clusterlist[i][j] << "  ";
		}
		tr.Info << std::endl;
	}
}

void ClusterBase::print_cluster_assignment() {
	int i,j;
	for ( i=0; i<(int)clusterlist.size(); i++ ) {
		for ( j=0; j<(int)clusterlist[i].size(); j++ ) {
			tr.Info << clusterlist[i][j]
				<< " " << core::pose::extract_tag_from_pose( poselist[ clusterlist[i][j] ]  )
				<< "  " << i << "  " << j << std::endl;
		}
	}
}

void ClusterBase::print_cluster_PDBs( std::string prefix ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	bool idealize_final = option[ OptionKeys::cluster::idealize_final_structures ]();
	int i,j;
	for ( i=0; i<(int)clusterlist.size(); i++ ) {
		//std::vector< int > new_cluster;
		for ( j=0; j<(int)clusterlist[i].size(); j++ ) {
			if ( export_only_low_ ) {
				if ( clusterlist[i].size() >= 1 && j!=0 ) continue; // print only clustercenter
			}
			//std::string output_name = "c." + string_of( i ) + "." + string_of( j ) + "." +
			//                         core::pose::extract_tag_from_pose( poselist[ clusterlist[i][j] ] ) + ".pdb";
			std::string output_name = "c." + string_of( i ) + "." + string_of( j ) + "." + "pdb";
			utility::replace_in( output_name, '/', "_" );

			// Idealize ?
			Pose pose;
			pose = poselist[ clusterlist[i][j] ];

			// make sure the parent/template name gets included in the REMARK line of the PDB
			using std::map;
			using std::string;
			map< string, string > score_line_strings( core::pose::get_all_score_line_strings( pose ) );
			for ( map< string, string >::const_iterator it = score_line_strings.begin(),
					end = score_line_strings.end();
					it != end; ++it ) {
				if ( it->first != "aln_id" ) continue;
				core::pose::add_comment( pose, "parents", it->second );
			}

			map< string, string > comments = core::pose::get_all_comments( pose );
			for ( map< string, string >::const_iterator it = comments.begin(),
					end = comments.end(); it != end; ++it
					) {
				if ( it->first != "aln_id" ) continue;
				core::pose::add_comment( pose, "parents", it->second );
			}


			if ( idealize_final ) {
				protocols::idealize::IdealizeMover idealizer;
				idealizer.fast( false );
				idealizer.apply( pose );
			}

			pose.dump_pdb( prefix + output_name );
		}
	}
}

std::vector< PoseOP >  ClusterBase::return_lowest_poses_in_clusters() {
	int i;
	std::vector< PoseOP > templist;
	for ( i=0; i<(int)clusterlist.size(); i++ ) {
		PoseOP tempPointer;
		if ( clusterlist[i].size() > 0 ) templist.push_back( core::pose::PoseOP( new Pose(poselist[ clusterlist[i][0] ]) ) );
	}
	return templist;
}

std::vector< core::pose::PoseOP >  ClusterBase::return_top_poses_in_clusters( core::Size count) {
	int i;
	std::vector< core::pose::PoseOP > templist;
	for ( i=0; i<(int)clusterlist.size(); i++ ) {

		core::pose::PoseOP tempPointer;
		if ( clusterlist[i].size() < 2 ) {
			//if singleton return the cluster center
			tempPointer = core::pose::PoseOP( new core::pose::Pose(poselist[ clusterlist[i][0]]) );
			templist.push_back(tempPointer);
		} else if ( clusterlist[i].size() >= 2 &&
				clusterlist[i].size() <= count+1 ) { //+1 because of cluster center
			for ( Size j = 1; j < clusterlist[i].size(); j++ ) {
				tempPointer = core::pose::PoseOP( new core::pose::Pose(poselist[ clusterlist[i][j]]) );
				templist.push_back(tempPointer);
			}
		} else if ( clusterlist[i].size() > count+1 ) {
			for ( Size j = 1; j <= count; j++ ) {
				tempPointer = core::pose::PoseOP( new core::pose::Pose(poselist[ clusterlist[i][j]]) );
				templist.push_back(tempPointer);
			}
		} else { }
	}
	return templist;
}

void ClusterBase::print_clusters_silentfile( std::string prefix ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	bool idealize_final = option[ OptionKeys::cluster::idealize_final_structures ]();

	using io::silent::SilentStructFactory;
	using io::silent::SilentStructOP;

	int i,j;
	io::silent::SilentStructOP ss, ss_center;
	io::silent::SilentFileData sfd, sfd_center;

	ss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
	std::string silent_file_ = option[ OptionKeys::out::file::silent ]();

	if ( option[ OptionKeys::cluster::write_centers ].user() ) {
		ss_center = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
		std::string silent_file_center_ = "centers_" + silent_file_;
		for ( i=0; i<(int)clusterlist.size(); ++i ) {
			std::string tag = "center_" + string_of( i );
			utility::replace_in( tag, '/', "_" );
			Pose pose;
			pose = poselist[ clusterlist[i].get_cluster_center() ];
			ss_center->fill_struct( pose, tag );
			sfd_center.write_silent_struct( *ss_center, silent_file_center_ );
		}
	}

	for ( i=0; i<(int)clusterlist.size(); i++ ) {
		//std::vector< int > new_cluster;
		for ( j=0; j<(int)clusterlist[i].size(); j++ ) {

			if ( export_only_low_ ) {
				if ( clusterlist[i].size() >= 1 && j!=0 ) continue; // print only clustercenter
			}

			std::string tag = prefix + "c." + string_of( i ) + "." + string_of( j ) + "." + "pdb";
			utility::replace_in( tag, '/', "_" );

			// Idealize ?
			Pose pose;
			pose = poselist[ clusterlist[i][j] ];
			if ( idealize_final ) {
				protocols::idealize::IdealizeMover idealizer;
				idealizer.fast( false );
				idealizer.apply( pose );
			}

			ss->fill_struct( pose, tag );

			if ( option[ OptionKeys::cluster::write_centers ].user() && i<10 ) { // this could be its own option
				std::string separate_file_name = "c" + string_of( i ) + "_" + silent_file_;
				sfd.write_silent_struct( *ss, separate_file_name );
			}

			sfd.write_silent_struct( *ss, silent_file_ );

		}
	}
}

void ClusterBase::create_constraints(
	std::string prefix,
	EnsembleConstraints &constraint_maker)
{
	int i,j;
	tr.Info << "Making constraints .. " << std::endl;
	for ( i=0; i<(int)clusterlist.size(); i++ ) {
		constraint_maker.clear();
		for ( j=0; j<(int)clusterlist[i].size(); j++ ) {
			constraint_maker.push_back(  poselist[ clusterlist[i][j] ]  );
		}

		std::string output_name = "c." + string_of( i ) + "." + "cst";
		utility::replace_in( output_name, '/', "_" );
		std::ofstream outf(  std::string(prefix + output_name).c_str() );
		constraint_maker.createConstraints( outf );
		outf.close();
	}

}


ClusterPhilStyle::ClusterPhilStyle() : ClusterBase()
{

}

std::string
ClusterPhilStyle::get_name() const {
	return "ClusterPhilStyle";
}

void ClusterPhilStyle::do_clustering( Size max_total_cluster ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	int listsize = poselist.size();
	if ( listsize <= 0 ) return;

	tr.Info << "Clustering an initial set of " << listsize << " structures " << std::endl;
	calculate_distance_matrix();

	int i,j;
	std::vector < int > neighbors ( poselist.size(), 0 );
	std::vector < int > clusternr ( poselist.size(), -1 );
	//std::vector < int > clustercenter;
	int    mostneighbors;
	Size    nclusters=0;

	if ( listsize == 0 ) {
		utility_exit_with_message( "Error: no Poses to cluster! Try -in:file:s or -in:file:silent!" );
	}

	if ( get_cluster_radius() < 0 ) {

		cluster_radius_ = get_median_rms()*1.1;

		// slightly different rule for gdtmm clsutering
		if ( option[ OptionKeys::cluster::gdtmm ]() ) {
			cluster_radius_ = get_median_rms() * 0.7;
			if ( cluster_radius_ > 5.0 ) cluster_radius_ = 5.0; // Cluster radius of GDTMM < 0.5 is just too coarse.
			if ( cluster_radius_ < 1.0 ) cluster_radius_ = 1.0; // Cluster radius of GDTMM > 0.9 is too fine for clustering
		}

		tr.Info << "Clustering of " << listsize << "structures with radius " <<  get_cluster_radius() << " (auto) " << std::endl;
	} else {
		tr.Info << "Clustering of " << listsize << "structures with radius " <<  get_cluster_radius() << std::endl;
	}

	std::vector <int> clustercentre;

	tr.Info << "Assigning initial cluster centres " << std::endl;
	// now assign groupings
	while ( true ) {
		// count each's neighbors
		for ( i=0; i<listsize; i++ ) {
			neighbors[i] = 0;
			if ( clusternr[i]>=0 ) continue; // ignore ones already taken
			for ( j=0; j<listsize; j++ ) {
				if ( clusternr[j]>=0 ) continue; // ignore ones already taken
				if ( distance_matrix( i+1, j+1 ) < get_cluster_radius() ) neighbors[i]++;
			}
		}

		mostneighbors = 0;
		for ( i=0; i<listsize; i++ ) {
			if ( neighbors[i]>neighbors[mostneighbors] ) mostneighbors=i;
		}
		if ( neighbors[mostneighbors] <= 0 ) break;  // finished!


		for ( i=0; i<listsize; i++ ) {
			if ( clusternr[i]>=0 ) continue; // ignore ones already taken
			if ( distance_matrix( i+1, mostneighbors+1 ) < get_cluster_radius() ) {
				clusternr[i] = mostneighbors;
			}
		}

		clustercentre.push_back(mostneighbors);
		nclusters++;

		if ( nclusters > max_total_cluster ) break;  // currently fixed but ought to be a paraemter
		if ( (nclusters%10)==0 ) tr.Info << ".";
		tr.Info.flush();
	}
	tr.Info << std::endl;

	for ( i=0; i<(int)clustercentre.size(); i++ ) {
		Cluster new_cluster( clustercentre[i] );
		new_cluster.clear();
		for ( j=0; j<listsize; j++ ) {
			// if that struture belongs to a given cluster
			if ( clusternr[j] == clustercentre[i] ) {
				new_cluster.add_member(j);         // add structure
			}
		}

		// add the new list to the clusterlist
		clusterlist.push_back(new_cluster);
	}


}

// this ensures every structure is in the cluster to who's cluster center it is most similar too
void ClusterPhilStyle::do_redistribution() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	int i,j;

	tr.Info << "Redistributing groups ..." << clusterlist.size() << " cluster centers";

	// redistribute groups - i.e. take each structure, calculate the rms to each cluster centre.
	for ( i=0; i<(int)clusterlist.size(); i++ ) {
		for ( j=1; j<(int)clusterlist[i].size(); j++ ) {
			int lowrmsi=i;
			Real lowrms=10000.0;
			for ( int m=0; m<(int)clusterlist.size(); m++ ) {
				Real rms;
				rms = distance_matrix( clusterlist[i][j]+1,                    // current structure vs
					clusterlist[m].get_cluster_center()+1); // clustercentre of cluster m

				if ( rms < lowrms ) {
					lowrms  = rms;
					lowrmsi = m;
				}
			}
			if ( lowrmsi != i ) { // is a different cluster centre closer to us than our current centre ?
				tr.Info << "Switched " << lowrmsi << "<--" << i << std::endl;
				// switch over;
				clusterlist[lowrmsi].add_member(clusterlist[i][j]);
				clusterlist[i].remove_member( j );
			}
		}

	}


}


std::string
ClusterPhilStyle_Loop::get_name() const {
	return "ClusterPhilStyle_Loop";
}

Real
ClusterPhilStyle_Loop::
get_distance_measure(
	const Pose & pose1,
	const Pose & pose2
) const
{
	return protocols::loops::loop_rmsd_with_superimpose(pose2, pose1, loop_def_, false );
}


//////////////////////////////////////////////////////////////////////////////
// Hook so ClusterPhilStyle can be used with PoseSelectors

std::string
ClusterPhilStyle_PoseReporter::get_name() const {
	return "ClusterPhilStyle_PoseReporter";
}

Real
ClusterPhilStyle_PoseReporter::
get_distance_measure(
	const Pose & pose1,
	const Pose & pose2
) const
{
	// NOTE: reporter_ takes non-const poses because one of the implemented reporters is the
	// EnergyReporter, which scores the pose. The score function alters the Energies object on the
	// pose so it takes non-const poses.
	// It's kind of rediculous that a pose copy is required here -- or making get_distance_measure
	// thing (and get_native_pose()) non-cost. For now, it's a copy unless it proves to be too slow.
	Pose pose1_copy(pose1);
	Pose pose2_copy(pose2);

	return reporter_->report_property(pose1_copy, pose2_copy);
}

//////////////////////////////////////////////////////////////////////////////


AssignToClustersMover::AssignToClustersMover( ClusterBaseOP cluster_base): GatherPosesMover()
{
	cluster_base_ = cluster_base;
}

void AssignToClustersMover::apply( Pose & pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	poselist.clear();
	GatherPosesMover::apply( pose );

	// check we hvnt added this structure before !
	if ( cluster_base_->check_tag( core::pose::extract_tag_from_pose( pose ) ) ) {
		tr.Info << "Already added: " << core::pose::extract_tag_from_pose( pose )  << std::endl;
		return;
	}

	//silent_score shoud be assigned in the GatherPosesMover::apply() step.
	Real score;
	getPoseExtraScore( pose , "silent_score", score );

	if ( score > cluster_base_->get_filter() ) {
		tr.Info << "Ignoring structure " << core::pose::extract_tag_from_pose( pose ) << " because score " << score << " is greater than threshold " << cluster_base_->get_filter() << std::endl;
		return; // ignore pose if filter value is too high
	}

	// find the best group by looping round all the existing groups and comparing
	// protected:
	cluster_base_->add_structure( pose );

	// to prevent memory crazyness, periodically limit back to the number of
	// desired structures!  giving preference to low energy clusters
	static int count=0;

	tr.Info << "Adding a "
		<< core::pose::extract_tag_from_pose( pose )
		<< "  " << score
		<< std::endl;

	if ( ( (count+1) %  150 ) == 0 ) {
		cluster_base_->print_summary();
		cluster_base_->sort_each_group_by_energy();
		cluster_base_->print_summary();
		cluster_base_->sort_groups_by_energy();
		cluster_base_->print_summary();
		cluster_base_->limit_groupsize( option[ OptionKeys::cluster::limit_cluster_size] );
		cluster_base_->limit_groups( option[OptionKeys::cluster::limit_clusters] );
		cluster_base_->limit_total_structures( option[ OptionKeys::cluster::limit_total_structures] );
		cluster_base_->clean_pose_store();
		cluster_base_->print_summary();
	}

	count++;
}

std::string
AssignToClustersMover::get_name() const {
	return "AssignToClustersMover";
}

std::string
EnsembleConstraints::get_name() const {
	return "EnsembleConstraints";
}

std::string
EnsembleConstraints_Simple::get_name() const {
	return "EnsembleConstraints_Simple";
}

void EnsembleConstraints_Simple::createConstraints( std::ostream &out) {
	int nres = poselist[0].total_residue();
	int residuesep = 5;
	Real strength = 1.0;

	out << "[ atompairs ]" << std::endl;

	for ( int ir = 1; ir <= nres; ir ++ ) {
		for ( int jr = 1; jr <= nres; jr ++ ) {
			if ( ir >= (jr - residuesep) ) continue;

			Real lowdist=1000000.0;
			Real highdist=-100000.0;
			for ( Size i=0; i<poselist.size(); i++ ) {
				Vector ir_CA = poselist[i].residue(ir).xyz( "CA" );
				Vector jr_CA = poselist[i].residue(jr).xyz( "CA" );
				Real dist = ir_CA.distance( jr_CA );
				if ( dist < lowdist ) lowdist = dist;
				if ( dist > highdist ) highdist = dist;
			}

			if ( lowdist > 11.0 ) continue;
			if ( highdist > 11.0 ) continue;

			if ( ( highdist -  lowdist ) < minimum_width_ ) {
				highdist += 0.5 * minimum_width_;
				lowdist  -= 0.5 * minimum_width_;
			}

			out << "     CA" << right_string_of(ir,7,' ')
				<< "     CA" << right_string_of(jr,7,' ')
				<< " BOUNDED"
				<< F( 12, 3, lowdist)
				<< F( 12, 3, highdist)
				<< F(4,1,strength )
				<< "  na"
				<< "; " <<  F( 12, 3, highdist - lowdist)
				<< std::endl;
		}
	}
}


} // cluster
} // protocols
