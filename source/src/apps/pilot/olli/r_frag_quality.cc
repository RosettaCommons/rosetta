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

#include <core/conformation/Conformation.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FragID_Iterator.hh>
// AUTO-REMOVED #include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/OrderedFragSet.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>


#include <core/fragment/util.hh>
#include <basic/Tracer.hh>

#include <protocols/jumping/JumpSample.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <core/fragment/SecondaryStructure.hh>

#include <core/scoring/rms_util.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>
#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>


#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>

#include <basic/options/option.hh>

#include <ObjexxFCL/string.functions.hh>
#include <basic/options/option_macros.hh>
#include <utility/excn/Exceptions.hh>


#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>


#include <protocols/toolbox/Cluster.hh>

// option key includes

#include <basic/options/keys/cluster.OptionKeys.gen.hh>
//#include <basic/options/keys/OptionKeys.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/fragment/FragData.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>

//Auto Headers
#include <core/kinematics/AtomTree.hh>


static basic::Tracer tr("main");


class ThisApplication  {
public:
	ThisApplication();
	static void register_options();
private:
};

ThisApplication::ThisApplication()
{}

OPT_KEY( File, f )
OPT_KEY( Boolean, torsion )
OPT_KEY( Integer, chop )
OPT_KEY( File, jumpss )
OPT_KEY( File, fold_tree )
OPT_KEY( Integer, ref_size )
OPT_KEY( File, exclude )
OPT_KEY( File, filter )
OPT_KEY( File, write )
OPT_KEY( Boolean, intrinsic )
OPT_KEY( File, write_big_cluster)
OPT_2GRP_KEY( FileVector, out, file, torsions )
//OPT_KEY( File, custer_check )
OPT_KEY( Integer, min_chop_in_quality_check )
OPT_KEY( File, fill_frags)
OPT_KEY( Integer, cluster_size)
OPT_1GRP_KEY( File, cluster, out )
OPT_1GRP_KEY( File, cluster, in )
OPT_1GRP_KEY( Real, cluster, acc_size )
OPT_1GRP_KEY( Real, cluster, acc_sum )
OPT_1GRP_KEY( Real, cluster, nmax )
OPT_1GRP_KEY( IntegerVector, cluster, range )
OPT_1GRP_KEY( File, out, qual )
OPT_KEY( File, ss_content )

using namespace core;
using namespace fragment;
using namespace pose;
using namespace kinematics;


using namespace protocols::jumping;
using namespace protocols::toolbox;


using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace ObjexxFCL::format;

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::native );
	OPT( in::file::s );
	NEW_OPT( f, "fragment file", "" );
	NEW_OPT( torsion, "test torsions of fragments", false );
	NEW_OPT( chop, "chop into 1mer frags first", 1 );
	NEW_OPT( jumpss, "jump-definition  -- makes us to check also the SS-library generated fragments","jump.def");
	NEW_OPT( fold_tree, "read jump-definitio as FOLD_TREE ","fold_tree.dat");
	NEW_OPT( ref_size,"make all comparisons with Xmers",0);
	NEW_OPT( exclude, "exclude the CA positions in this file","");
	NEW_OPT( filter, "filter nfrag file, according to entries in filter-filee","");
	NEW_OPT( intrinsic, "compute the std-deviation of the torsion dofs",false);
	NEW_OPT( write,"write fragments to file","fragout.dat");
	NEW_OPT( ss_content, "write ss_content to file","");
	NEW_OPT( min_chop_in_quality_check, "just ignore less than 9mers for the final quality check",5 );
	NEW_OPT( write_big_cluster,"write big cluster fragments to file","");
	//	NEW_OPT( cluster_check, "write analysis of big cluster frags to this file", "
	NEW_OPT( cluster_size,"complement frags from largest cluster with frags from fill_frags until this number is reached",30);
	NEW_OPT( fill_frags,"fill with these fragmenTs where no big cluster is found","");
	NEW_OPT( cluster::out, "file to write cluster definitions to", "" );
	NEW_OPT( cluster::in, "file to read cluster definitions from","");
	NEW_OPT( cluster::nmax, "how many clusters are maximally created ?",10);

	NEW_OPT( cluster::acc_sum, "fraction of fragments that has to be clustered (sum of all)",0.8);
	NEW_OPT( cluster::acc_size,"minimum fraction of fragments in winning cluster",0.7);
	NEW_OPT( cluster::range, "work only on these residues", 0);
	NEW_OPT( out::qual, "quality_data", "frag_qual.dat" );
	NEW_OPT( out::file::torsions, "torsion angles first file name takes the phi, the second takes the psi angles", "" );
	OPT( cluster::radius );
}

Real check_jump( pose::Pose const& pose, pose::Pose const& native, JumpSample const& , int downstream_res_nr );


Real compare_cartesian_rmsd( Pose const &orig_frag, Pose const &pred_frag ) {
	Size const nres( orig_frag.total_residue() );
	assert( nres == pred_frag.total_residue() );
	Size const cmp ( 9 ); //compare 9mers
	Size const ncmp ( nres - cmp + 1 );
	Real total ( 0 );
	for ( Size i = 1; i <= ncmp; ++i ) {
		tr.Debug << "get rmsd " << i << " " << i+cmp-1 << std::endl;
		total += scoring::CA_rmsd( orig_frag, pred_frag, i, i+cmp-1);
	}
	return total/ncmp;
	//	return scoring::rmsd_with_super( orig_frag, pred_frag, scoring::is_protein_backbone );
}

inline Real sqr ( Real x ) {
	return x*x;
}

Real compare_torsion_rmsd( Pose const &orig_frag, Pose const &pred_frag ) {
	Real err ( 0.0 );
	for ( Size pos = 1; pos <= orig_frag.total_residue(); ++pos ) {
		for ( Size dof = 1; dof <= 2; ++dof ) { //check phi and psi
			Real orig = orig_frag.torsion( id::TorsionID( pos, id::BB, dof ) );
			Real pred = pred_frag.torsion( id::TorsionID( pos, id::BB, dof ) );
			//	std::cout << orig << ' ' << pred << ' ' << orig-pred << " "
			//					<< numeric::nearest_angle_degrees(orig-pred,0.0)  << std::endl;
			err += sqr( numeric::nearest_angle_degrees(orig-pred,0.0) );
		}
	}
	return sqrt( err / 2*orig_frag.total_residue() );
}

Real compare_frags_pose( Pose const &native, Pose const &test_pose, Frame const& frame, utility::vector1< Size > const& excl ) {
	if ( !frame.is_continuous() ) utility_exit_with_message("can only determine rmsd for cont. frames");

	Size const start( frame.start() );

	Size cmp( frame.length() );
	if ( option[ ref_size ].user() ) {
		cmp = option[ ref_size ];
	}
	Size const ncmp ( frame.length() - cmp + 1 );
	Real total ( 0 );
	for ( Size i = start; i <= start+ncmp-1; ++i ) {
		tr.Debug << "get rmsd " << i << " " << i+cmp-1 << std::endl;
		total += scoring::CA_rmsd( native, test_pose, start, start+cmp-1, excl);
	}
	return total/ncmp;
	//	return scoring::rmsd_with_super( orig_frag, pred_frag, scoring::is_protein_backbone );
}


bool compute_min_mean_rmsd_frag( FrameList& frames, Pose const& native_pose, Real& min_rmsd, Real& mean_rmsd, Real &max_rmsd ) {
	Pose test_pose( native_pose );
	mean_rmsd = 0.0;
	max_rmsd = 0;
	min_rmsd = 1000;
	bool chops( false );
	Size ct( 0 );
	for ( FragID_Iterator it = frames.begin(); it != frames.end(); ++it ) {
		it->apply( test_pose );
		utility::vector1< Size > excl;
		Real r =  compare_frags_pose( native_pose, test_pose, it->frame(), excl );
		if ( it->frame().length() < 9 ) chops = true;
		if ( (Size) option[ min_chop_in_quality_check ]() <= it->frame().length() ) {
			r *= ( 9 / it->frame().length() ); //total heuristic!
			if ( min_rmsd > r ) min_rmsd = r;
			if ( max_rmsd < r ) max_rmsd = r;
			mean_rmsd += r;
			++ct;
		}
	}
	mean_rmsd /= ct;
	return chops;
}


void check_quality_of_cluster_frags( Pose const& native_pose, FragSetOP decoy_frags, FragSetOP fill_frags, FragSetOP new_frags ) {

	utility::io::ozstream os( std::string( static_cast<std::string>(option[ write_big_cluster]()))+".quality" ); //I had to do a static_cast here to stop gcc 4.4 from complaining about option being an array of pointers
	//	std::ostream& os( std::cout );
	for ( Size pos = 1; pos <= new_frags->max_pos(); pos++ ) {
		os << RJ(6,pos) << " ";
		FrameList frames;
		Real min_rmsd, mean_rmsd, max_rmsd;

		new_frags->frames( pos, frames );
		bool chops=compute_min_mean_rmsd_frag( frames, native_pose, min_rmsd, mean_rmsd, max_rmsd );
		os << RJ( 10, mean_rmsd ) << " " << RJ( 10, min_rmsd ) << " " << RJ( 10, max_rmsd ) << " ";

		frames.clear();
		decoy_frags->frames( pos, frames );
		compute_min_mean_rmsd_frag( frames, native_pose, min_rmsd, mean_rmsd, max_rmsd );
		os << RJ( 10, mean_rmsd ) << " " << RJ( 10, min_rmsd ) << " " << RJ( 10, max_rmsd ) << " ";

		frames.clear();
		fill_frags->frames( pos, frames );
		compute_min_mean_rmsd_frag( frames, native_pose, min_rmsd, mean_rmsd, max_rmsd );
		os << RJ( 10, mean_rmsd ) << " " << RJ( 10, min_rmsd ) << " " << RJ( 10, max_rmsd ) << " ";
		os << ( chops ? 1 : 0 );
		os << std::endl;
	}
}

void write_cluster_frags( FragSetOP predicted_frags, FragSetOP fill_frags, FragSetOP new_frags ) {
	typedef ClusterBase::ClusterList ClusterList;

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

	//movemap to keep track where we have GOOD clusters  -- will be set to FALSE where GOOD fragments are found
	// we use the movemap to track this, because the core/fragment supports easy fragment chopping at the fringes ( 9mer .....xxxx --> 5mer )
	kinematics::MoveMap mm;
	mm.set_bb( true );

	//iterate over frames to read all clusters
	for ( ConstFrameIterator itframe = predicted_frags->begin(), eitframe=predicted_frags->end();
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

		//compute cluster sum for this position
		Size sum( 0 );
		for ( ClusterBase::ClusterIterator it = cluster.begin(), eit =cluster.end();
					it!=eit; ++it ) {
			sum += it->size();
		}
		clust_sum[ pos ] = sum;

		Size const N( frame.nr_frags() );
		Size const clust_sum_value( static_cast< Size >(option[ cluster::acc_sum ]()*N ) );
		Size const min_cluster_size( static_cast< Size >( option[ cluster::acc_size ]()*N) );
		tr.Info << "pos: " << pos << " min_cluster_size: " << min_cluster_size << " cluster_sum "<< clust_sum_value << std::endl;

		//decide whether to use these fragments
		if ( sum >= clust_sum_value && cluster.begin()->size() >= min_cluster_size) mm.set_bb( pos, false );
	} // read clusters

	//create insert_map -- > fast and painfree access to the fragments we like
	InsertMap insert_map;
	InsertSize insert_size;
	pose::Pose pose; //empty pose --> jump-frags won't be applicable: pose(s) no problem here.
	fill_frags->generate_insert_map( mm, insert_map, insert_size );

// 	// we need complete numbering:
// 	InsertMap insert_map( nres, 0 );
// 	InsertSize insert_size( nres, 0 );
// 	for ( Size i = 1; i <= _insert_map.size(); i++ ) {
// 		insert_map[ _insert_map[ pos ] ] = i;
// 		insert_size[ _insert_map[ pos ] ] = _insert_size[ i ];
// 	}
	//DEBUG output
	Size const total_insert = insert_map.size();
	tr.Trace << "size of insertmap: " << total_insert << " -- ";
	for ( Size i = 1; i<=total_insert; i++ ) tr.Trace << " " << insert_map[ i ];
	tr.Trace << "insert_size: \n";
	for ( Size i = 1; i<=insert_map[ total_insert ]; i++ ) tr.Trace << " " << RJ(3,i);
	tr.Trace <<"\n";
	for ( Size i = 1; i<=insert_map[ total_insert ]; i++ ) tr.Trace << " " << RJ(3,insert_size[ i ]);
	tr.Trace << std::endl;

	//if insert_size is 9 take all fragments ( decoys and fill_frags )
	//	in insert_size smaller 9 take fragments that are taken in positions to the right ( clusters ) and add chopped std_frags
	//	if insert_size 0 take cluster frags

	for ( ConstFrameIterator itframe = predicted_frags->begin(), eitframe=predicted_frags->end();
				itframe != eitframe; ++itframe ) {
		Size pos = (*itframe)->start();
		ClusterList cluster = clusters[ pos ];

		Size const N( (*itframe)->nr_frags() );
		Size const clust_sum_value( static_cast< Size >(option[ cluster::acc_sum ]()*N ) );
		Size const min_cluster_size( static_cast< Size >( option[ cluster::acc_size ]()*N) );


		if ( !mm.get_bb( pos ) ) { //we have confidence in our decoys: no fill_frags

			//get all the fragments from the first cluster
			ClusterList::const_iterator it = cluster.begin();

			for ( ClusterBase::IntraClusterIterator fit = it->begin(); fit != it->end(); ++fit ) {
				if ( *fit > (*itframe)->nr_frags() ) {
					tr.Warning << "IGNORE: more elements in cluster than in frame... frame:" << (*itframe)->nr_frags()
										 << " in cluster: " << *fit << std::endl;
					tr.Warning << " position: " << pos << std::endl;
				} else {
					FragID frag( *itframe, *fit );
					new_frags->add( frag );
				}
			}
			continue; // don't add anyting else at this position
		}

		if ( insert_size[ pos ] == 9 ) { //no confident cluster anywhere to the right
			new_frags->add( *itframe );//add all fragments
			//and fill up with standard fragments and fragments from the decoys
			FrameList decoy_frames;
			predicted_frags->frames( pos, decoy_frames );
			new_frags->add( decoy_frames );
			if ( fill_frags ) {
				FrameList fill_frames;
				fill_frags->frames( pos, fill_frames); //put roughly the same amount of fill_frags in there!
				Size const count ( std::max( 1, static_cast< int > ( decoy_frames.size() / fill_frames.size() ) ) );
				for ( Size i=1; i<= count; i++ ) new_frags->add( fill_frames );
			}
			continue; //don't add anything else... we already have both: decoy_fragments and original_fragments
		}

		// somewhere less than 9 residues from here we have convergence:
		// 1.) add full-length frags for those frag_id's those frags are in a downstream winning cluster
		// 2.) fill up with chopped frags that do not overlap with the converged region

		//vector of boolean, to track which frag_ids we want to keep
		typedef utility::vector1< bool> InsertSet;
		InsertSet insert_set( max_frag_id, false );

		//at insert_size[pos] the first good cluster starts,
		// off - from here (pos ) to the cluster ( cpos ).. go up to 8 since a 9mer can't overlap with a cluster that is 9 down from here
		//---- with old core/fragment: off = insert_size[ pos ]
		for ( Size off = 1; off <= 8; off++ ) {
			Size cpos = pos + off;

			// set insert_set[]=true for fragments from decoys whose fragments are in a good downstream cluster.
			ClusterList cluster = clusters[ cpos ];
			if ( clust_sum[ cpos ] > clust_sum_value ) {
				for ( ClusterList::const_iterator it = cluster.begin(), eit = cluster.end(); it != eit; ++it ) {
					if ( it->size() > min_cluster_size ) {
						//keep only fragments from cluster
						for ( ClusterBase::IntraClusterIterator fit = it->begin(), efit = it->end();
									fit != efit; ++fit ) {
							insert_set[ *fit ] = true;
						}
					}
				}
			}
		} //for ( Size off... ) loop to find insertable frag_ids
		if ( tr.Trace.visible() ) {
			Size ct( 0 );
			for ( Size i = 1; i<=max_frag_id; i++ ) ct += insert_set[ i ];
			tr.Trace << "mixed bag at pos " << pos << " found " << ct << " fragments from good downstream cluster " << std::endl;
		}
		// now insert frames
		FrameList fill_frames; // if not an accepted decoy frame, it can be still used as fill_frame

		// get decoys fragments
		FrameList decoy_frames;
		predicted_frags->frames( pos, decoy_frames);
		runtime_assert( decoy_frames.size() == 1);
		FrameCOP frame( decoy_frames[ 1 ] );

		// collect fragments from bad decoys here
		FrameOP decoy_fill_frame = new Frame( pos, decoy_frames[ 1 ]->length() );

		for ( Size nr = 1; nr <= frame->nr_frags(); ++nr ) {
			//from good decoy ?
			if ( insert_set[ nr ] ) { //from good decoy->add straight to new_frags
				FragID frag( frame, nr );
				new_frags->add( frag );
			} else { //not from good decoy--> add to decoy_fill for chopping
				decoy_fill_frame->add_fragment( frame->fragment_ptr( nr ) );
			}
		}
		//fill_frames: combine frags from bad decoys...
		if ( decoy_fill_frame->nr_frags() ) {
			fill_frames.push_back( decoy_fill_frame );
		}
		//... and from the standard library fill_frags
		fill_frags->frames( pos, fill_frames );
		if ( tr.Trace.visible() ) {
			tr.Trace << " fill_frames have " << fill_frames.size() << " frames with: ";
			for ( FrameList::iterator it = fill_frames.begin(), eit = fill_frames.end(); it != eit; ++it ) 	tr.Trace << (*it)->nr_frags() << " ";
			tr.Trace << "number of fragments" << std::endl;
		}

		//chop everything in fill_frames and add to new_frags
		Size const chop( insert_size[ pos ] );
		for ( FrameList::iterator it = fill_frames.begin(), eit = fill_frames.end(); it != eit; ++it ) {
			FrameOP chop_frame = new Frame( pos, chop );
			for ( Size nr = 1; nr <= (*it)->nr_frags(); nr ++ ) {
				chop_frame->add_fragment( (*it)->fragment( nr ).generate_sub_fragment( 1, chop ) );
			}
			new_frags->add( chop_frame );
		}

	} //big loop over all positions


	//at this point we should have a complete fragment set in new_frags...
	runtime_assert( new_frags );
	tr.Info << " fill_frags: " << fill_frags->max_frag_length() << " new_frag: " << fill_frags->max_frag_length() << std::endl;


	//now decide if we want to chop all fragments to match up with the size of the standard fragments
	if ( fill_frags->max_frag_length() < new_frags->max_frag_length() ) {
		//ups fill_frags are e.g, 3mers but our new_frags where collected as 9mers...
		//chop them to 3mers
		Size const chops( fill_frags->max_frag_length() );
		tr.Info << " fill_frags are only " << fill_frags->max_frag_length() <<"mers. Going to chop all longer frags "  << std::endl;

		//chopping...
		FragSetOP chopped_frags = new OrderedFragSet;
		for ( ConstFrameIterator it = new_frags->begin(), eit = new_frags->end(); it != eit; ++it ) {
			Frame const& fr( **it );

			//already small enough --- don't chop
			if ( fr.length() <= chops ) {
				chopped_frags->add( *it );
				continue;
			}

			//chop...
			for ( Size pos = fr.start(); pos <= fr.end() - chops + 1; pos ++ ) {
				FrameOP cf = new Frame( pos, chops );
				for ( Size nr = 1; nr <= fr.nr_frags(); ++nr ) {
					cf->add_fragment( fr.fragment( nr ).generate_sub_fragment( pos-fr.start() + 1, pos-fr.start()+ chops ) );
				}
				chopped_frags->add( cf );
			}

		} //finished chopping .. only tail left to sort out

		//fill up the tail with fill_frags. the starting postions will never reach that region
		for ( Size pos = chopped_frags->max_pos() + 1; pos <= fill_frags->max_pos(); pos++ ) {
			FrameList fr;
			fill_frags->frames( pos, fr );
			chopped_frags->add( fr );
		}

		//done --- replace new_frags with the chopped fragments
		new_frags = chopped_frags;
	}// if chopping

	//write new fragment library
	FragmentIO().write_data( option[ write_big_cluster](), *new_frags );
}

void compute_intrinsic_deviation( Pose& test_pose, FragSetOP predicted_frags, Pose const& native_pose ) {
	/* ---            how to read output CLUSTER  ----

		 pos iRMSD(all_frames)  CLUSTER_1 CLUSTER_2 ... CLUSTER_20

		 each CLUSTER_i has 4 numbers
		 		 nr_elem meanRMSD(to native) minRMSD( to native) intrinsicRMSD(in_cluster)

 	 */
	utility::io::ozstream cluster_out;
	if ( option[ cluster::out ].user() ) {
		cluster_out.open( option[ cluster::out ]() );
	}

	utility::io::izstream cluster_in;
	if ( option[ cluster::in ].user() &&  !option[ cluster::out ].user() ) {
		cluster_in.open( option[ cluster::in ]() );
		if ( !cluster_in.good() ) utility_exit_with_message( "can't open cluster file " + std::string( option[ cluster::in ]()) );
	}

	for ( ConstFrameIterator itframe = predicted_frags->begin(), eitframe=predicted_frags->end();
				itframe != eitframe; ++itframe ) {
		Frame const& frame ( **itframe );

		//testing
		//		if ( frame.start() != 100 ) continue;

		if ( frame.nr_frags() == 0 ) continue;
		if ( option[ cluster::range ].user() ) {
			utility::vector1< int > range = option[ cluster::range ]();
			utility::vector1< int >::const_iterator iter = find( range.begin(), range.end(), (int) frame.start() );
			if ( iter == range.end() ) {
				tr.Info << "skip residue " << frame.start() << std::endl;
				continue;
			}
		}

		if ( option[ torsion ]() ) {
			runtime_assert( frame.length() == 1 );
			Real phi_av( 0.0 ), psi_av( 0.0 );
			Real phi_std( 0.0 ), psi_std( 0.0 );

			Size pos = frame.start();
			id::TorsionID phiID( pos, id::BB, 1 );
			id::TorsionID psiID( pos, id::BB, 2 );

			for ( Size i = 1; i<= frame.nr_frags(); i++ ) {
				frame.apply( i, test_pose );
				Real phi = test_pose.torsion( phiID );
				Real psi = test_pose.torsion( psiID );
				phi_av+=phi;
				psi_av+=psi;
				phi_std+=phi*phi;
				psi_std+=psi*psi;
			}
			Real invn = 1.0/frame.nr_frags();
			Real const dev ( std::sqrt( invn*(phi_std-phi_av*phi_av*invn + psi_std-psi_av*psi_av*invn) ) );
			std::cout << RJ(6,frame.start()) << RJ(10,dev) << std::endl;
		} else { //rmsd


			ClusterPhilStyle cluster( frame.nr_frags() );
			if ( option[ cluster::in ].user() &&  !option[ cluster::out ].user() ) {
				std::string line;
				getline( cluster_in, line );
				std::istringstream str( line );
				Size pos;
				str >> pos;
				if ( pos != frame.start() ) utility_exit_with_message(" didn't find correct cluster in file " );
				str >> cluster;
				std::cout << RJ(6,frame.start()) << RJ(10,"nan");
			} else {
				//compute dist matrix
				//compute rmsd for each pair of fragments at this position
				Pose pose2( test_pose );
				Real total( 0.0 );
				for ( Size i = 1; i<=frame.nr_frags(); i++ ) {
					frame.apply( i, test_pose );
					for ( Size j = 1; j<=frame.nr_frags(); j++ ) {
						frame.apply( j, pose2 );
						Real dist = scoring::CA_rmsd( pose2, test_pose, frame.start(), frame.end() );
						cluster.dist( i, j ) = dist;
						total += dist;
					}
				}
				Real invn = 1.0/frame.nr_frags();
				Real const dev ( total*invn*invn );
				std::cout << RJ(6,frame.start()) << RJ(10,dev);
				cluster.set_radius( option[ cluster::radius ]() ); //flex_rad was of no help for T374... //std::min( option[ cluster::max_rad ](), dev ) );
				cluster.set_n_max_cluster( (core::Size)option[ cluster::nmax ]() );
				cluster.compute();
			};

			if ( cluster_out.good() ) cluster_out << frame.start() << " " << cluster << std::endl;
			//		cluster.print_cluster_assignment( tr.Info );
			Size const nout( option[ cluster::nmax ] );
			for ( Size ncl = 1; ncl <= (Size) std::min( (int) nout, (int) cluster.size() ); ncl ++ ) {
				Real total( 0 );
				Real min_rms( 1000 );
				for ( ClusterPhilStyle::IntraClusterIterator it = cluster.cluster( ncl ).begin(), eit = cluster.cluster( ncl ).end();
							it != eit; ++it ) {
					Size frag_nr = *it;
					frame.apply( frag_nr, test_pose );
					utility::vector1< Size > excl;
					Real rms =  compare_frags_pose( native_pose, test_pose, frame, excl );
					if ( min_rms > rms ) min_rms = rms;
					total += rms;
				}
				Size const nelem( cluster.cluster( ncl ).size() );
				std::cout << RJ( 4, nelem ) << " " << F( 6, 3, total /nelem ) << " " << F( 6, 3, min_rms )<< " ";
				// get intrinsic cluster rms
				if ( !option[ cluster::in ].user() || option[ cluster::out ].user() ) { //if cluster_in... no distance matrix
					//compute average distance within cluster -- intrinsic RMSD
					Real itotal ( 0 );
					for ( Size k1 = 0; k1 < nelem; k1++ ) { //clusters are std::container
						for ( Size k2 = 0; k2 < nelem; k2++ ) {
							itotal+=cluster.dist( cluster.cluster( ncl )[ k1 ], cluster.cluster( ncl )[ k2 ] );
						}
					}
					std::cout << F( 6, 3, itotal*(1.0/nelem)*(1.0/nelem));
				} else {
					std::cout << RJ( 6, "nan");
				}
			}
			for ( Size ncl = cluster.size()+1; ncl <= nout; ncl ++ ) {
				std::cout << RJ( 4, 0 ) << " " << RJ( 5, "nan" ) << " " << RJ( 5, "nan") << " " << RJ( 5, "nan");
			}
			std::cout << std::endl;
		} //torsion or rmsd
	} //for frame
}

Real compare_frags( Pose const &orig_frag, Pose const &pred_frag ) {
	if ( option[ torsion ]() ) {
		return compare_torsion_rmsd( orig_frag, pred_frag );
	}
	return compare_cartesian_rmsd( orig_frag, pred_frag );
}

FragSetOP filter_frags( FragSet const& frags_in, std::string const& filter_file ) {
	utility::io::izstream filter( filter_file );
	std::string line;
	FragSetOP new_frags = frags_in.empty_clone();
	while ( getline( filter, line ) ) {
		std::istringstream line_str( line );
		Size size, pos, frag_nr;
		line_str >> size >> pos >> frag_nr;
		if ( line_str.fail() ) {
			utility_exit_with_message("read error in file: " + filter_file + " at line " + line );
		}
		tr.Info << "use " << size << "mer fragment: " << pos << " "  << frag_nr << std::endl;
		FrameList frames;
		frags_in.frames( pos, frames );

		// find right size
		FrameList::iterator it = frames.begin(), eit = frames.end();
		while( it != eit ) {
			if ( (*it)->length() == size ) break;
			++it;
		}
		if ( it == eit ) {
			tr.Error << "ERROR: no frame with size " << size << " was found at position " <<  pos << std::endl;
		} else {
			tr.Info << "frame with " << (*it)->nr_frags() << " found" << std::endl;
			if ( (*it)->nr_frags() >= frag_nr ) {
				FrameOP new_frame = (*it)->clone();
				new_frame->add_fragment( (*it)->fragment_ptr( frag_nr ));
				new_frags->add( new_frame );
			} else {
				tr.Error << "ERROR: not enough fragments in frame " << (**it) << std::endl;
			}
		}
	}
	return new_frags;
}




int main( int argc, char** argv ) {
	try{
	ThisApplication::register_options();
	devel::init( argc, argv );

	//NEW_OPT( cluster::radius, "radius for clustering", 1.0);
	if ( !basic::options::option[ basic::options::OptionKeys::cluster::radius ].user() ) {
		basic::options::option[ basic::options::OptionKeys::cluster::radius ].def( 1.0 );
	}

	if ( option[ out::file::torsions ].user() && option[ in::file::s ].user() ) {
		utility::io::ozstream out_phi;
		utility::io::ozstream out_psi;
		out_phi.open( option[out::file::torsions ]()[1] );
		out_psi.open( option[out::file::torsions ]()[2] );
		for ( Size ct=1; ct <= option[ in::file::s ]().size(); ct++ ) {
			Pose pose;
			//read it
			core::import_pose::pose_from_pdb( pose, option[ in::file::s ]()[ ct ] );
			for ( Size pos=1; pos <= pose.total_residue(); pos++ ) {
				out_phi << RJ( 6, pos ) <<  RJ(6, ct) << F(10,4, pose.phi( pos ) ) << std::endl;
				out_psi << RJ( 6, pos ) <<  RJ(6, ct) << F(10,4, pose.psi( pos ) ) << std::endl;
			}
		}
	}

	kinematics::MoveMap move_all;

	std::string const native_pdb ( option[ in::file::native ]() );

	utility::vector1< Size > excl;
	if ( option[ exclude ].user() ) {
		utility::io::izstream file( option[ exclude ] );
		Size pos;
		while ( file >> pos ) {
			excl.push_back( pos );
		}
	}

	Pose native;
	//read it
	core::import_pose::pose_from_pdb( native, native_pdb );
	core::util::switch_to_residue_type_set( native, chemical::CENTROID );

	Pose test_pose;
	core::pose::make_pose_from_sequence(
		test_pose,
		native.sequence(),
		*( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ))
	);

	bool const bJumps ( option[ jumpss ].user() || option[ fold_tree ].user() );

	if ( !bJumps ) {

		FragSetOP orig_frags;
		orig_frags = FragmentIO().read_data( option[ f ]() );

		if ( option[ filter ].user() ) {
			orig_frags = filter_frags( *orig_frags, option[ filter ] );
			FragmentIO().write_data( option[ OptionKeys::write ](), *orig_frags );
		}

		if ( option[ ss_content ].user() ) {
			core::fragment::SecondaryStructure ss_def( *orig_frags, true /*no JustUseCentralResidue */ );
			utility::io::ozstream out( option[ ss_content ]() );
			ss_def.write_psipred_ss2( out, native.sequence() ); // << ss_def << std::endl;
		}



		FragSetOP predicted_frags = NULL;
		if ( ( option[ intrinsic ] && option[ torsion ] ) || option[ chop ].user() ) {
			ConstantLengthFragSetOP short_frags = new ConstantLengthFragSet( option[ chop ] );
			chop_fragments( *orig_frags, *short_frags );
			predicted_frags = short_frags;
			FragmentIO().write_data( "dump_chop.dat", *short_frags);
		} else {
			predicted_frags = orig_frags;
		}

			Size max_pdb( 0 );
		typedef std::map< std::string, Size > PDB_IDS;
		PDB_IDS pdb_code;

		if ( option [ intrinsic ] ) {
			compute_intrinsic_deviation( test_pose, predicted_frags, native );
		}
		if ( option [ write_big_cluster ].user() ) {

			FragSetOP fill_frags( FragmentIO( 25 ).read_data( option[ OptionKeys::fill_frags ]() ) );
			FragSetOP new_frags( new OrderedFragSet );

			write_cluster_frags( predicted_frags, fill_frags, new_frags );

			check_quality_of_cluster_frags( native, predicted_frags, fill_frags, new_frags );

			return 0;

		}
		utility::io::ozstream out_phi;
		utility::io::ozstream out_psi;
		if ( option[ out::file::torsions ].user() ) {
			out_phi.open( option[out::file::torsions ]()[1] );
			out_psi.open( option[out::file::torsions ]()[2] );
			if ( option[ in::file::native ].user() ) {
				for ( Size pos=1; pos <= native.total_residue(); pos++ ) {
					out_phi << RJ( 6, pos ) <<  RJ(6, 0) << F(10,4, native.phi( pos ) ) << std::endl;
					out_psi << RJ( 6, pos ) <<  RJ(6, 0) << F(10,4, native.psi( pos ) ) << std::endl;
				}
			}
		}
		if ( !option [ intrinsic ] ) {
			utility::io::ozstream out( option[ out::qual ] );
			for ( ConstFrameIterator frame = predicted_frags->begin(), eframe=predicted_frags->end();
						frame != eframe; ++frame ) {
				for ( Size i=1; i<=frame->nr_frags(); i++ ) {
					Size len = frame->length();
					//					if ( len < 9 ) continue;
					std::string const pdb (frame->fragment( i ).pdbid() );
					Size pdb_id;
					if ( pdb.size() > 4 ) {
						if ( !pdb_code[ pdb ] ) pdb_code[ pdb ] = ++max_pdb;
						pdb_id = pdb_code[ pdb ] ;
					} else pdb_id = 0;
					frame->apply( i, test_pose );
					if ( option[ out::file::torsions ].user() ) {
						for ( Size pos = frame->start(); pos <= frame->end(); pos ++ ) {
							out_phi << RJ( 6, pos ) <<  RJ(6,i) << F(10,4, test_pose.phi( pos ) ) << std::endl;
							out_psi << RJ( 6, pos ) <<  RJ(6,i) << F(10,4, test_pose.psi( pos ) ) << std::endl;
						}
					} else {
						out << RJ(6, len) << RJ(6,frame->start()) << RJ(6,i) << F(10,4,compare_frags_pose( native, test_pose, **frame, excl))
								<< RJ(6, pdb_id ) << RJ( 6, frame->fragment( i ).pdbpos() )
								<< std::endl;
					}
				}
			}
		}
		for ( PDB_IDS::const_iterator it= pdb_code.begin(), eit = pdb_code.end(); it!=eit; ++it ) {
			std::cout << "MAPPING: " << it->second << " " << it->first << std::endl;
		}
//

// 	Size max_pdb( 0 );
// 		typedef std::map< std::string, Size > PDB_IDS;
// 		PDB_IDS pdb_code;
// 		for ( FrameIterator frame = predicted_frags->begin(), eframe=predicted_frags->end();
// 					frame != eframe; ++frame ) {
// 			Pose orig_frag;
// 			Pose pred_frag;
// 			//			frame->fragment_as_pose( 1, orig_frag);
// 			frame->fragment_as_pose( 1, pred_frag);
// 			FrameOP native_frame = frame->clone_with_template();
// 			native_frame->steal( native );
// 			native_frame->fragment_as_pose( 1, orig_frag );

// 			for ( Size i=1; i<=frame->nr_frags(); i++ ) {
// 				Size len = frame->fragment( i ).apply( pred_frag, 1, frame->length() );
// 				if ( len < 9 ) continue;
// 				std::string const pdb (frame->fragment( i ).pdbid() );
// 				Size pdb_id;
// 				if ( pdb.size() > 4 ) {
// 					if ( !pdb_code[ pdb ] ) pdb_code[ pdb ] = ++max_pdb;
// 					pdb_id = pdb_code[ pdb ] ;
// 				} else pdb_id = 0;
// 				std::cout << RJ(6, len) << RJ(6,frame->start()) << RJ(6,i) << RJ(10,compare_frags( orig_frag, pred_frag))
// 									<< RJ(6, pdb_id ) << RJ( 6, frame->fragment( i ).pdbpos() )
// 									<< std::endl;
// 			}
// 		}
// 		for ( PDB_IDS::const_iterator it= pdb_code.begin(), eit = pdb_code.end(); it!=eit; ++it ) {
// 			std::cout << "MAPPING: " << it->second << " " << it->first << std::endl;
// 		}
// //


// 		ConstantLengthFragSet native_frags( predicted_frags->max_frag_length() );
// 		steal_constant_length_frag_set_from_pose( native , native_frags );

// 		Pose orig_frag; //an original fragment
// 		Pose pred_frag; //a predicted fragment

// 		FrameIterator fr_nat = native_frags.begin();

// 		//initialize poses
// 		fr_nat->fragment_as_pose( 1, orig_frag);
// 		fr_nat->fragment_as_pose( 1, pred_frag);

// 		for ( FrameIterator
// 						efr_nat=native_frags.end(),
// 						fr_pred=predicted_frags->begin(),
// 						efr_pred=predicted_frags->end();
// 					fr_nat!=efr_nat && fr_pred!=efr_pred;
// 					++fr_pred, ++fr_nat
// 		) {
// 			fr_nat->fragment( 1 ).apply( orig_frag, 1, fr_nat->length() );
// 			for ( Size i=1; i<=fr_pred->nr_frags(); i++ ) {
// 				Size len = fr_pred->fragment( i ).apply( pred_frag, 1, fr_nat->length() );
// 				std::cout << RJ(6, len) << RJ(6,fr_pred->start()) << RJ(6,i) << RJ(10,compare_frags( orig_frag, pred_frag)) << std::endl;
// 			}
// 		}
	} else { //bJumps
		JumpSetup jump_def( native.total_residue() );
		JumpSample jump_setup;
		if ( option[ jumpss ].user() ) {
			jump_def.read_file( option[ jumpss ] );
			//if ( jump_def.size() != 1 ) { //don't deal with other cases
			//			utility_exit_with_message("give me a single jump each time you run this: found "+string_of( jump_def.size() ) + " jumps");
			//		}
			jump_setup = JumpSample ( jump_def );
		}	else if ( option[ fold_tree ].user() ) {
			utility::io::izstream file( option[ fold_tree ]() );
			core::kinematics::FoldTree f;
			file >> f;
			jump_setup = JumpSample( f );
		}

		jump_setup.set_fold_tree_in_pose( native );
		assert( native.fold_tree().num_jump() == jump_setup.size() );
		FrameList jumps;
		jump_setup.steal_orientation_and_pleating( native );
		kinematics::MoveMap mm;
		mm.set_bb( true );
		jump_setup.generate_jump_frags( *StandardPairingLibrary::get_instance(), mm, false /* bWithTorsion */, jumps );
		tr.Info << "JUMPS: " << jump_setup << std::endl;
		int nrj = 1;
		for ( FrameList::const_iterator it = jumps.begin(), eit = jumps.end(); it!=eit; ++it, ++nrj ) {
			pose::Pose pose( native );
			for ( Size i=1; i<=(*it)->nr_frags(); ++i ) { // (*it)->nr_frags(); i++) {
				(*it )->fragment( i ).apply( pose, **it );
				kinematics::Edge jump_edge = pose.fold_tree().get_residue_edge( (*it)->end() );
				tr.Info << "jump: " << jump_edge.label() << " " << (*it)->start() << " " << (*it)->end() << " jump_frag " << RJ(3,i) << " ";
				check_jump( pose, native, jump_setup, (*it)->end() ); //);
			}
		}
		{
			// pose::Pose pose( native );
// 			for ( FrameList::const_iterator it = jumps.begin(), eit = jumps.end(); it!=eit; ++it ) {
// 				(*it )->fragment( 123 ).apply( pose, **it );
// 				//				int ii = nr_jumps;
// 			}
			//torsion frag setup
			//			ConstantLengthFragSetOP orig_frags = new ConstantLengthFragSet;
			//			orig_frags->read_fragment_file( frag_file );
			//			kinematics::MoveMapOP mm = new kinematics::MoveMap;
			//			mm ->set_bb( true );

			//			ClassicFragmentMover mover( orig_frags, mm );

			//			for ( Size pos = 1; pos<=15; pos++) {
			//					mm->set_bb( pos, false );
			//			}

			//			pose.set_psi( 1, -45 );

			// native = pose;
// 			for ( Size ii = 1; ii<=jump_setup.size();  ii++ ) {
// 				tr.Info << "START " <<  " jump: " << ii <<" ";
// 				check_jump( pose, native, jump_setup, ii ); );
// 			}


// 			for ( Size pos = 1; pos <= pose.total_residue(); pos++ ) {
// 								if ( pos == 21 || pos == 20 || pos == 22 ) continue;
// 				pose.set_phi( 128, -45 );
// 				pose.set_psi( pos, -45 );
// 				pose.set_omega( pos, 180 );
// 			}


// 			for ( int cycl = 1; cycl <= 5000; cycl++ ) {
// 				for ( Size ii = 1; ii<=jump_setup.size();  ii++ ) {
// 					tr.Info << "MOVED " << cycl << " jump: " << ii <<" ";
// 					check_jump( pose, native, jump_setup, ii ); );
// 				}
// 				torsion-frag apply
// 				mover.apply( pose );
// 			}
		}
		// now try some frag insertion and look if jump-qual stays invariant
	}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}








/* =================================================*/
Real check_jump( pose::Pose const& pose, pose::Pose const& native, JumpSample const& , int downstream_res_nr ) {

	/*
	{
		kinematics::Jump nat_jump = native.jump( 1 );
		kinematics::Jump pred_jump = pose.jump( 1 );
		Real dist, theta;
		jump_distance( nat_jump, pred_jump, dist, theta );
		tr.Info << F(10,4,dist) << F( 10,4,numeric::conversions::degrees( theta )) ;
		//	return -1;
	}
	*/
	//		kinematics::Edge jump_edge = pose.fold_tree().jump_edge( jump_nr );
	kinematics::Edge jump_edge = pose.fold_tree().get_residue_edge( downstream_res_nr );
		//work out upstream and downstream residue
		Size res1=jump_edge.start();
		Size res2=jump_edge.stop();
		tr.Info << RJ(5,res1) << RJ(5,res2)<<" ";
		// work out the stubID
		chemical::ResidueType const& rt1 ( pose.residue_type ( res1 ) );
		chemical::ResidueType const& rt2 ( pose.residue_type ( res2 ) );

		id::AtomID a1( rt1.atom_index ("N") , res1 );
		id::AtomID a2( rt1.atom_index ("CA") , res1 );
		id::AtomID a3( rt1.atom_index ("C") , res1 );
		id::StubID down_stub( a1, a2, a3 );

		id::AtomID b1( rt2.atom_index ("N") , res2 );
		id::AtomID b2( rt2.atom_index ("CA") , res2 );
		id::AtomID b3( rt2.atom_index ("C") , res2 );
		id::StubID up_stub( b1, b2, b3 );

	// 	Real dist = distance( pose.xyz( a1 ), pose.xyz( b3 ));
// 		Real nat_dist = distance( native.xyz( a1 ), native.xyz( b3 ));
// 		tr.Info <<  " dist: " << F( 10, 4, dist-nat_dist );

// 		dist = distance( pose.xyz( a3 ), pose.xyz( b3 ));
// 		nat_dist = distance( native.xyz( a3 ), native.xyz( b3 ));
// 		tr.Info <<  " dist: " << F( 10, 4, dist-nat_dist );

// 		dist = distance( pose.xyz( b1 ), pose.xyz( a3 ));
// 		nat_dist = distance( native.xyz( b1 ), native.xyz( a3 ));
// 		tr.Info <<  " dist: " << F( 10, 4, dist-nat_dist );

// 		dist = distance( pose.xyz( b1 ), pose.xyz( b3 ));
// 		nat_dist = distance( native.xyz( b1 ), native.xyz( b3 ));
// 		tr.Info <<  " dist: " << F( 10, 4, dist-nat_dist );

		Stub up = pose.conformation().atom_tree().stub_from_id( up_stub );
		Stub down = pose.conformation().atom_tree().stub_from_id( down_stub );
		RT rt(up, down);

		Stub native_up = native.conformation().atom_tree().stub_from_id( up_stub );
		Stub native_down = native.conformation().atom_tree().stub_from_id( down_stub );
		RT rt_native(  native_up,native_down );
	// 	tr.Info << " rtdist " << F(10, 4, distance(rt, rt_native) );
// 		kinematics::Jump nat_jump ( rt_native );
// 		kinematics::Jump pred_jump ( rt );
		//		Real dist, theta;
		//		jump_distance( nat_jump, pred_jump, dist, theta );
		//		tr.Info << F(10,4,dist) << F( 10,4,numeric::conversions::degrees( theta )) ;
		//		return -1;


		//apply rt to native_up
		Stub zero_test_down;
		rt_native.make_jump( native_up, zero_test_down );
		//	tr.Info << distance(native_down, zero_test_down );
		Stub test_down;
		rt.make_jump( native_up, test_down );

		Real zero_rms ( 0.0 );
		Real rms ( 0.0 );

		for ( Size i=1; i<=3; i++ ) {
			Vector tv = test_down.build_fake_xyz( i );
			Vector ztv = zero_test_down.build_fake_xyz( i );
			Vector nv = native_down.build_fake_xyz( i );
			Vector d = nv-tv;
			Vector zd = nv-ztv;
			rms += d.length();
			zero_rms += zd.length();
		}
		tr.Info << " : " << F(10,3,rms) <<  std::endl;
		return rms;
}
