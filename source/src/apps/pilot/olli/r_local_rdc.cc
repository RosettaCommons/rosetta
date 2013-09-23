// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/r_frag_quality.cc
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange


//#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/FragID.hh>
//#include <core/fragment/BBTorsionSRFD.hh>

#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>

//#include <core/fragment/OrderedFragSet.hh>

//#include <core/kinematics/MoveMap.hh>
//#include <core/kinematics/Stub.hh>
//#include <core/kinematics/RT.hh>


#include <core/fragment/util.hh>

//#include <protocols/jumping/JumpSample.hh>
//#include <protocols/jumping/JumpSetup.hh>
//#include <protocols/jumping/SecondaryStructure.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/rms_util.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>
#include <devel/init.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/chemical/ChemicalManager.hh>


#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>

#include <basic/options/option.hh>
#include <core/scoring/ResidualDipolarCoupling.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/options/option_macros.hh>


#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <protocols/toolbox/Cluster.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

// option key includes

#include <basic/options/keys/cluster.OptionKeys.gen.hh>
//#include <basic/options/keys/OptionKeys.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>



static basic::Tracer tr("main");

class ThisApplication  {
public:
  ThisApplication();
  static void register_options();
private:
};

ThisApplication::ThisApplication()
{}

OPT_1GRP_KEY( File, in, frags )
OPT_1GRP_KEY( File, out, qual )
OPT_1GRP_KEY( File, cluster, in )

using namespace core;
using namespace fragment;
using namespace pose;
using namespace kinematics;
using namespace protocols::toolbox;

//using namespace protocols::jumping;


using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace scoring;
using namespace ObjexxFCL::format;

void ThisApplication::register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  OPT( in::file::native );
	OPT( in::file::rdc );
  NEW_OPT( in::frags, "fragment file", "" );
  NEW_OPT( out::qual, "quality_data", "" );
	NEW_OPT( cluster::in, "cluster file which determines which frags RDC scores are averaged", "");
}




//typedef utility::vector1< FragID > FragID_List;
typedef utility::vector1< FragID_List > FragID_ClusterList;

void recover_clustered_frags( FragSetOP predicted_frags, FragID_ClusterList& clustered_frags ) {
	typedef ClusterBase::ClusterList ClusterList;
	clustered_frags.clear();
	//if all clusters together have at least this value, we don't allow any std-fragset stuff
	FrameList fra;
	predicted_frags->frames( 1, fra);
	assert( fra.size() );

	//Size const nres( fill_frags->max_pos() );
	Size max_frag_id( 0 );

	utility::io::izstream cluster_in;
	cluster_in.open( option[ cluster::in ]() );
	if ( !cluster_in.good() ) utility_exit_with_message( "can't open cluster file " + std::string( option[ cluster::in ]()) );

	//as many clusters as residues
	typedef utility::vector1< ClusterList > Clusters;
	Clusters clusters( predicted_frags->max_pos() );

	typedef utility::vector1< Size > SizeVector;
	SizeVector clust_sum( predicted_frags->max_pos(), 0 );

	//iterate over frames to read all clusters
	for ( FrameIterator itframe = predicted_frags->begin(), eitframe=predicted_frags->end();
				itframe != eitframe; ++itframe ) {
		Frame const& frame ( **itframe );

		if ( frame.nr_frags() == 0 ) continue;
		if ( frame.nr_frags() > max_frag_id ) max_frag_id = frame.nr_frags();

		//retrieve corresponding cluster from file

		//first check position
		std::string line;
		getline( cluster_in, line );
		std::istringstream str( line );
		Size pos;
		str >> pos;
		if ( pos != frame.start() ) utility_exit_with_message(" didn't find correct cluster in file " );

		//now read clusters for this position
		ClusterList cluster;
		str >> cluster;
		clusters[ pos ] = cluster;

	} // read clusters

	for ( FrameIterator itframe = predicted_frags->begin(), eitframe=predicted_frags->end();
				itframe != eitframe; ++itframe ) {
		Size pos = (*itframe)->start();
		ClusterList all_cluster_at_pos = clusters[ pos ];
		tr.Debug << (*itframe)->start() << " " << std::endl;
		for ( ClusterList::const_iterator cluster_it = all_cluster_at_pos.begin(); cluster_it != all_cluster_at_pos.end(); ++cluster_it ) {
			FragID_List frags_in_cluster;
			for ( ClusterBase::IntraClusterIterator frag_in_cluster_it = cluster_it->begin(); frag_in_cluster_it != cluster_it->end(); ++frag_in_cluster_it ) {
				if ( *frag_in_cluster_it > (*itframe)->nr_frags() ) {
					tr.Warning << "IGNORE: more elements in cluster than in frame... frame:" << (*itframe)->nr_frags()
										 << " in cluster: " << *frag_in_cluster_it << std::endl;
					tr.Warning << " position: " << pos << std::endl;
				} else {
					tr.Debug << *frag_in_cluster_it << " ";
					FragID frag( *itframe, *frag_in_cluster_it );
					frags_in_cluster.push_back( frag );
				}
			}
			tr.Debug << std::endl;
			clustered_frags.push_back( frags_in_cluster );
		}
	}
}


ResidualDipolarCouplingOP filter_rdcs_for_frame( Frame const& frame, ResidualDipolarCoupling const& orig_rdcs ) {
  residual_dipolar_coupling::RDC_lines const& rdcs = orig_rdcs.get_RDC_data();
  residual_dipolar_coupling::RDC_lines filtered;

  for ( residual_dipolar_coupling::RDC_lines::const_iterator it = rdcs.begin(); it != rdcs.end(); ++it ) {
    if ( frame.moves_residue( it->res1() ) ) {
			if ( rdc_it->res1() < it
	filtered.push_back( *it );
    }
  }
  return new ResidualDipolarCoupling( filtered );
}

Real compare_cartesian_rmsd( Pose const &orig_frag, Pose const &pred_frag, Frame const& frame ) {
  if ( !frame.is_continuous() ) utility_exit_with_message("can only determine rmsd for cont. frames");
	Size const nres( orig_frag.total_residue() );
	Size const nres2( pred_frag.total_residue() );
	assert( nres == nres2 );

	return scoring::CA_rmsd( orig_frag, pred_frag, frame.start(), frame.end() );
}

void score_clustered_frags( FragSetOP frags, Pose& test_pose, Pose& native_pose ) {

  ResidualDipolarCoupling original_rdcs; //reads automatically from file
  Size max_pdb( 0 );

  scoring::ScoreFunction scorefxn;
  scorefxn.set_weight( scoring::rdc, 1 );

  utility::io::ozstream out( option[ out::qual ] );

	FragID_ClusterList clustered_frags;
	recover_clustered_frags( frags, clustered_frags );

	for ( FragID_ClusterList::const_iterator cluster_it = clustered_frags.begin();
				cluster_it != clustered_frags.end(); ++cluster_it ) {
		Real mean_rmsd( 0.0 );
		Real min_rmsd( 100.0 );
		Real max_rmsd( 0.0 );
		Real mean_rdc_score( 0.0 );
		Size len;
		Size pos;
		Size ct;
		for ( FragID_List::const_iterator frag_it = cluster_it->begin(); frag_it != cluster_it->end(); ++frag_it ) {
			Frame const& frame( frag_it->frame() );
			scoring::store_RDC_in_pose( filter_rdcs_for_frame( frame, original_rdcs ), test_pose );
			len = frame.length();
			if ( frame.start() != pos ) ct=0;
			pos = frame.start();
      frag_it->apply( test_pose );
      Real const score( scorefxn( test_pose ) );
			Real rmsd;
			rmsd = compare_cartesian_rmsd( native_pose, test_pose, frame );
			if ( min_rmsd > rmsd ) min_rmsd = rmsd;
			if ( max_rmsd < rmsd ) max_rmsd = rmsd;
			mean_rmsd += rmsd/cluster_it->size();
			mean_rdc_score += score/cluster_it->size();
		}
		++ct;
		out << RJ(6, len) << RJ(6,pos) << " " << RJ(4,ct)<< " " << F(10,4, mean_rdc_score ) << " "
				<< F(10,4, mean_rmsd ) <<" " << F(10,4,min_rmsd)	<< " " << RJ(4,cluster_it->size()) << std::endl;
	}
}

int main( int argc, char** argv ) {
	try{
  ThisApplication::register_options();
  devel::init( argc, argv );

  Pose native;
  //read it
  std::string const native_pdb ( option[ in::file::native ]() );
  core::import_pose::pose_from_pdb( native, native_pdb );
  core::util::switch_to_residue_type_set( native, chemical::CENTROID );

  Pose test_pose;
	Pose clean_pose;
  core::pose::make_pose_from_sequence(
				    test_pose,
				    native.sequence(),
				    *( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ))
  );
	clean_pose = test_pose;
  FragSetOP orig_frags;
  orig_frags = FragmentIO().read( option[ in::frags ]() );

	if ( option[ cluster::in ].user() ) {
		score_clustered_frags( orig_frags, test_pose, native );
		return 0;
	}


  ResidualDipolarCoupling original_rdcs; //reads automatically from file
  Size max_pdb( 0 );

  scoring::ScoreFunction scorefxn;
  scorefxn.set_weight( scoring::rdc, 1 );

  utility::io::ozstream out( option[ out::qual ] );

	Size frame_ct=1;
  for ( FrameIterator frame = orig_frags->begin(), eframe=orig_frags->end();
				frame != eframe; ++frame, ++frame_ct ) {
    scoring::store_RDC_in_pose( filter_rdcs_for_frame( **frame, original_rdcs ), test_pose );
		std::list< Size > sel_residues;

		if ( !frame->is_continuous() ) {
			kinematics::FoldTree fold_tree;
			test_pose = clean_pose;
			make_simple_fold_tree_from_jump_frame( **frame, test_pose.total_residue(), fold_tree );
			fold_tree.put_jump_stubs_intra_residue();
			test_pose.fold_tree( fold_tree );
			for ( Size i = 1; i<=frame->length(); ++i ) {
				sel_residues.push_back( frame->seqpos( i ) );
			}
			tr.Trace << " frame: "; frame->show_header( tr.Trace );
			tr.Trace << fold_tree << std::endl;
		}
		test_pose.dump_pdb("original.pdb");

		for ( Size i=1; i<=frame->nr_frags(); i++ ) {
      Size len = frame->length();
      frame->apply( i, test_pose );
      Real const score( scorefxn( test_pose ) );
			if ( tr.Trace.visible() ) test_pose.dump_pdb("test_pose"+string_of( i )+".pdb");
			Real rmsd;
			if ( frame->is_continuous() ) rmsd = compare_cartesian_rmsd( native, test_pose, **frame );
			else rmsd = scoring::CA_rmsd( test_pose, native, sel_residues );

      out << RJ(6, len) << RJ(6, frame_ct ) << RJ(6,frame->start()) << RJ(6,i) << " " << F(10,4, score ) << " "
					<< F(10,4, rmsd )	<< std::endl;
    }
  }
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl; 
	} 
}

