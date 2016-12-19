// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>  //make_pose_from_sequence
#include <core/pose/util.hh>

#include <core/import_pose/import_pose.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/util/SwitchResidueTypeSet.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <core/conformation/ResidueFactory.hh>

// minimize pose into density
#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/DensitySymmInfo.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/electron_density/DockIntoDensityMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/idealize/IdealizeMover.fwd.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OptH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>


#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>

#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <basic/Tracer.hh>

#include <iostream>
#include <string>
#include <list>
#include <algorithm>

#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

OPT_KEY( String, mode )
OPT_KEY( String, ca_positions )
OPT_KEY( String, startmodel )
OPT_KEY( String, scorefile )
OPT_KEY( String, symm_type )
OPT_KEY( File, fragfile )
OPT_KEY( Integer, num_frags )
OPT_KEY( IntegerVector, pos )
OPT_KEY( IntegerVector, designated_rank )
OPT_KEY( Integer, movestep )
OPT_KEY( Integer, ncyc )
OPT_KEY( Boolean, min_pack_min )
OPT_KEY( Boolean, min_bb )
OPT_KEY( Boolean, norm_scores )
OPT_KEY( Integer, bw )
OPT_KEY( Integer, n_to_search )
OPT_KEY( Integer, n_filtered )
OPT_KEY( Integer, n_output )
OPT_KEY( Boolean, verbose )
OPT_KEY( Boolean, native_placements )
OPT_KEY( Real, delR )
OPT_KEY( Real, clust_radius )
OPT_KEY( Real, scale_cycles )
OPT_KEY( Real, point_radius )
OPT_KEY( Integer, clust_oversample )
OPT_KEY( Integer, n_matches )
OPT_KEY( Real, frag_dens )
OPT_KEY( IntegerVector, frag_len )
OPT_KEY( RealVector, assembly_weights )
OPT_KEY( Real, null_weight )
OPT_KEY( Real, consensus_frac )
OPT_KEY( Real, consensus_stdev )
OPT_KEY( Real, energy_cut )
OPT_KEY( Boolean, cheat )

using namespace ObjexxFCL::format;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::pose;
using namespace core;

static basic::Tracer TR("denovo_density");

// helper function: read CAs from a PDB
void
ReadCAsFromPDB( std::string pdbfile, utility::vector1< numeric::xyzVector<core::Real> > &cas ) {
	cas.clear();

	std::ifstream inpdb(pdbfile.c_str());
	std::string buf;

	while ( std::getline(inpdb, buf ) ) {
		if ( buf.substr(0,4)!="ATOM" && buf.substr(0,6)!="HETATM" ) continue;
		if ( buf.substr(12,4) != " CA " ) continue;

		numeric::xyzVector<core::Real> ca_i(
			atof(buf.substr(30,8).c_str()),
			atof(buf.substr(38,8).c_str()),
			atof(buf.substr(46,8).c_str())
		);

		cas.push_back( ca_i );
	}
}

// store a CA trace
struct CAtrace {
	CAtrace() {
		dens_score_ = 0;
		rms_ = 0;
		tag_ = "";
		idx_ = 0;
	}

	CAtrace(core::pose::Pose pose, std::string tag, core::Real score, core::Real rms) {
		cas_.reserve( pose.size() );
		for ( int i=1; i<=(int)pose.size(); ++i ) {
			if ( pose.residue(i).is_protein() ) cas_.push_back( pose.residue(i).xyz(2) );
		}
		dens_score_ = score;
		rms_ = rms;
		tag_ = tag;
		idx_ = 0;
	}

	utility::vector1< numeric::xyzVector< core::Real> > cas_;
	core::Real dens_score_, rms_;
	std::string tag_;
	core::Size idx_;
};

// per fragment info
struct FragID {
	FragID( int pos, std::string tag) {
		pos_=pos;
		tag_=tag;
	}
	int pos_;
	std::string tag_;
};

// ignores idx (pos/tag pair is unique)
bool operator<(const FragID& x, const FragID& y)
{
	return boost::make_tuple(x.pos_,x.tag_) < boost::make_tuple(y.pos_,y.tag_);
}

/// fpd
///   -- none of these are movers, since they do not follow the 1 struct in / 1 struct out model
///   -- they may be at some point ...
class DockFragmentsMover  {
public:
	//
	void run( );

	// remove density around pose from map
	void cut_from_map( core::pose::Pose const &pose );

	// read and preprocess fragment file
	void process_fragfile();

	// cluster fragfile
	core::fragment::FragSetOP
	cluster_frags(core::fragment::FragSetOP fragments);

	virtual std::string get_name() const {
		return "DockFragments";
	}

private:
	std::map<core::Size, core::fragment::Frame> library_, libsmall_;
	core::Size nmer_, nmer_small_;
};

///
class ScoreFragmentSetMover {
public:
	//
	ScoreFragmentSetMover();

	void run( );

	core::Real overlap_score(
		CAtrace &pose1,
		CAtrace &pose2,
		core::Size offset,
		protocols::electron_density::DensitySymmInfo const &symminfo);

	core::Real clash_score(
		CAtrace &pose1,
		CAtrace &pose2,
		core::Size offset,
		protocols::electron_density::DensitySymmInfo const &symminfo);


	core::Real closability_score(
		CAtrace &pose1,
		CAtrace &pose2,
		core::Size offset,
		protocols::electron_density::DensitySymmInfo const &symminfo);

	virtual std::string get_name() const {
		return "ScoreFragmentSet";
	}

private:
	// parameters
	core::Real steepness_, clash_dist_, overlap_width_, unclosable_penalty_;
	utility::vector1< core::Real > gap_lengths_, gap_weights_;
};

///
class FragmentAssemblyMover {
public:
	//
	FragmentAssemblyMover();

	void run( );

	// compute the score for the current assignment
	core::Real score(bool);

	virtual std::string get_name() const {
		return "FragmentAssembly";
	}

private:
	// parameters
	core::Real overlap_wt_, clash_wt_, close_wt_, dens_wt_;
	core::Real null_frag_;

	// data structures
	utility::vector1< FragID > allfrags;  // all fragment indices
	std::map< FragID, int > frag2idx;  // all fragment indices
	utility::vector1< utility::vector1<int> > pos2frags; // frag_candidates as fragidx
	utility::vector1<int> assigned_frags; // current MC assignment
	utility::vector1< utility::vector1< float > > scores_2b;  // keep as float
	utility::vector1<core::Real> scores_1b;
	utility::vector1<core::Real> rmses;

	// debug
	utility::vector1< utility::vector1< float > > clash2b, overlap2b, close2b;
};

///
class ConsensusFragmentMover {
public:
	void run( );

	virtual std::string get_name() const {
		return "ConsensusFragment";
	}
};

///
class SolutionRescoreMover {
public:
	void run( );

	virtual std::string get_name() const {
		return "SolutionRescore";
	}
};


///
core::fragment::FragSetOP
DockFragmentsMover::cluster_frags(core::fragment::FragSetOP fragments) {
	core::Size nfrags = basic::options::option[ num_frags ];
	core::Real fragfilter = 10;

	core::fragment::FragSetOP fragments_clust = fragments->empty_clone();

	for ( core::fragment::ConstFrameIterator it = fragments->begin(); it != fragments->end(); ++it ) {
		core::fragment::Frame frame = **it;
		core::fragment::FrameOP filteredframe = frame.clone();

		core::Size nfrags_i = std::min( nfrags, frame.nr_frags() );
		for ( core::Size i_frag=1; i_frag<=nfrags_i; i_frag++ ) {
			core::fragment::FragDataCOP oldfrag = frame.fragment_ptr(i_frag);

			bool addfrag=true;
			core::Size storedfragcount = filteredframe->nr_frags();
			for ( core::Size j_frag=1; j_frag<=storedfragcount && addfrag; j_frag++ ) {
				core::fragment::FragDataCOP storedfrag = filteredframe->fragment_ptr(j_frag);

				core::Real error = 0;
				core::Size ntors = 0;
				for ( core::Size i_res=1; i_res<=oldfrag->size(); ++i_res ) {
					core::fragment::BBTorsionSRFD const & frag_i =
						dynamic_cast< core::fragment::BBTorsionSRFD const &> ( *(oldfrag->get_residue(i_res)) );
					core::fragment::BBTorsionSRFD const & frag_j =
						dynamic_cast< core::fragment::BBTorsionSRFD const &> ( *(storedfrag->get_residue(i_res)) );

					for ( core::Size i_tors=1; i_tors<=frag_i.nbb(); ++i_tors ) {
						core::Real err_i = (frag_i.torsion(i_tors) - frag_j.torsion(i_tors));
						error += err_i*err_i;
						ntors++;
					}
				}
				core::Real RMS = std::sqrt(error/ntors);
				if ( RMS <= fragfilter ) addfrag=false;
			}

			if ( addfrag ) filteredframe->add_fragment(oldfrag);
		}

		TR.Debug << "pos " << (*it)->start() << " mer " << (*it)->length()
			<< "  nfrags " << frame.nr_frags() << " -> " << filteredframe->nr_frags() << std::endl;

		fragments_clust->add(filteredframe);
	}

	return fragments_clust;
}

///
void DockFragmentsMover::process_fragfile() {
	/////
	// read in fragments file
	// process

	core::fragment::FragSetOP fragments, fragments_small;
	fragments = core::fragment::FragmentIO().read_data( basic::options::option[ fragfile ]() );
	nmer_ = fragments->max_frag_length(); // assumes constant length fragments!

	utility::vector1< int > nmer_target = basic::options::option[ frag_len ]();
	if ( nmer_target.size() == 0 ) {
		nmer_target.push_back(nmer_);
		nmer_target.push_back(nmer_);
	}
	runtime_assert( nmer_target.size() == 2);

	core::Size nmer_target_big = (core::Size) std::max( nmer_target[1],nmer_target[2] );
	core::Size nmer_target_small = (core::Size) std::min( nmer_target[1],nmer_target[2] );

	if ( nmer_target_big<nmer_ ) {
		TR << "Chopping to " << nmer_target_big << " mers" << std::endl;
		fragments_small = core::fragment::FragSetOP( new core::fragment::ConstantLengthFragSet( nmer_target_big ) );
		core::fragment::chop_fragments( *fragments, *fragments_small );
		fragments = fragments_small->clone();
	}

	if ( nmer_target_small<nmer_target_big ) {
		TR << "Chopping to " << nmer_target_small << " mers" << std::endl;
		fragments_small = core::fragment::FragSetOP( new core::fragment::ConstantLengthFragSet( nmer_target_small ) );
		core::fragment::chop_fragments( *fragments, *fragments_small );
	} else {
		fragments_small = fragments->clone();
	}

	nmer_ = nmer_target_big;
	nmer_small_ = nmer_target_small;

	core::fragment::FragSetOP fragclust = cluster_frags( fragments );
	core::fragment::FragSetOP fragsmallclust = cluster_frags( fragments_small );

	// map resids to frames
	for ( core::fragment::ConstFrameIterator i = fragclust->begin(); i != fragclust->end(); ++i ) {
		core::Size position = (*i)->start();
		library_[position] = **i;
	}
	for ( core::fragment::ConstFrameIterator i = fragsmallclust->begin(); i != fragsmallclust->end(); ++i ) {
		core::Size position = (*i)->start();
		libsmall_[position] = **i;
	}
}



///
void DockFragmentsMover::run() {
	// get sequence from fasta or native pdb if provided
	std::string sequence;
	core::pose::Pose native_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_file( native_pose, option[ in::file::native ]().name() , core::import_pose::PDB_file);
		sequence = native_pose.sequence();
	} else {
		sequence = core::sequence::read_fasta_file( option[ in::file::fasta ]()[1])[1]->sequence();
	}

	// read fragments, preprocess
	process_fragfile();

	// setup symmetry
	protocols::electron_density::DensitySymmInfo symminfo( basic::options::option[ symm_type ]() );
	symminfo.detect_axes( core::scoring::electron_density::getDensityMap() );

	// figure out: which fragments to use, what positions to steal
	utility::vector1<bool> use_big( sequence.length() , true );
	utility::vector1<bool> search_positions( sequence.length() , true );
	utility::vector1<bool> steal_positions( sequence.length() , false );

	if ( basic::options::option[ pos ].user() ) {
		utility::vector1< core::Size > inpos = basic::options::option[ pos ]();
		for ( core::Size i=1; i<=search_positions.size(); ++i ) {
			search_positions[i] = (std::find( inpos.begin(), inpos.end(), i ) != inpos.end());
		}
	}

	core::Real STRAND_LONGFRAG_CUT = 0.5; // if at least this percent of fragdata is extended, use shortfrags

	for ( std::map<core::Size, core::fragment::Frame>::iterator it=library_.begin(); it!=library_.end(); ++it ) {
		core::Size idx = it->first;
		core::fragment::Frame &f = it->second;

		// foreach fragment
		core::Size nres=0,nstrand=0;
		for ( core::Size i=1; i<=f.nr_frags(); ++i ) {
			core::fragment::FragData const &f_i = f.fragment( i );
			for ( core::Size j=1; j<=f_i.size(); ++j ) {
				nres++;
				if ( f_i.secstruct(j) == 'E' ) { nstrand++; }
			}
		}

		core::Real strand_frac = ((core::Real)nstrand)/((core::Real)nres);
		if ( strand_frac > STRAND_LONGFRAG_CUT ) {
			use_big[idx] = false;
		}
	}

	bool have_initial_pose = option[ startmodel ].user();
	core::pose::Pose initial_pose;
	std::map<int,int> initial_pose_seqmap;
	core::Size nsteal=0;
	if ( have_initial_pose ) {
		core::import_pose::pose_from_file( initial_pose, option[ startmodel ]() , core::import_pose::PDB_file);

		// trim map
		cut_from_map( initial_pose );

		// adjust search positions
		utility::vector1<bool> init_pose_covers(search_positions.size());
		for ( int i=1; i<=(int)initial_pose.size(); ++i ) {
			init_pose_covers[initial_pose.pdb_info()->number(i)] = true;
			initial_pose_seqmap[initial_pose.pdb_info()->number(i)] = i;
		}
		for ( int i=1; i<=(int)search_positions.size(); ++i ) {
			if ( !search_positions[i] ) continue;
			core::Size nmer = use_big[i]? nmer_ : nmer_small_;

			bool frag_i_covered = init_pose_covers[i];
			for ( int j=i+1; j<i+(int)nmer && frag_i_covered; ++j ) {
				frag_i_covered = frag_i_covered && init_pose_covers[j];
			}

			search_positions[i] = !frag_i_covered;
			steal_positions[i] = frag_i_covered;
			nsteal++;
		}
	}

	//
	TR << "Searching positions:";
	for ( int i=1; i<=(int)search_positions.size(); ++i ) {
		if ( search_positions[i] ) TR << " " << i;
	}
	TR << std::endl;
	if ( nsteal>0 ) {
		TR << "Stealing positions:";
		for ( int i=1; i<=(int)search_positions.size(); ++i ) {
			if ( steal_positions[i] ) TR << " " << i;
		}
		TR << std::endl;
	}

	// set up docking
	protocols::electron_density::DockIntoDensityMover dock;
	dock.setDelR(option[ delR ]); // option?
	dock.setB( option[ bw ]() );
	dock.setTopN( option[ n_to_search ]() , option[ n_filtered ]() , option[ n_output ]() ); //y
	dock.setGridStep(option[ movestep ]()); //y
	dock.setMinBackbone(option[ min_bb ]()); //y
	dock.setNCyc(option[ ncyc ]()); //y
	dock.setClusterRadius(option[ clust_radius ]()); //y
	dock.setPointRadius(option[ point_radius ]()); //y
	dock.setFragDens(option[ frag_dens ]()); //y
	dock.setNormScores(option[ norm_scores ]());
	dock.setClusterOversamp(option[ clust_oversample ]()); //y

	int maxRotPerTrans = (int)std::ceil( (core::Real)option[ n_filtered ]() / (core::Real)option[ n_to_search ]() );
	dock.setMaxRotPerTrans( maxRotPerTrans );

	if ( option[ out::file::silent ].user() ) {
		std::string silent_fn = option[ out::file::silent ]();
		dock.setOutputSilent( silent_fn );
	}

	// read CA positions (if specified)
	if ( option[ ca_positions ].user() ) {
		utility::vector1< numeric::xyzVector<core::Real> > cas;
		ReadCAsFromPDB( option[ ca_positions ](), cas );
		dock.predefine_search(cas);
		dock.setCenterOnMiddleCA(true);
	}

	// use symmetry (if specified)
	if ( symminfo.enabled() ) {
		dock.setSymminfo(symminfo);
	}

	// for each position
	for ( std::map<core::Size, core::fragment::Frame>::iterator it=library_.begin(); it!=library_.end(); ++it ) {
		core::Size idx = it->first;
		core::fragment::Frame frame;
		if ( use_big[idx] ) {
			frame = library_[idx];
		} else {
			frame = libsmall_[idx];
		}

		core::Size seq_pos = frame.start();
		core::Size nmer_len = frame.length();
		if ( !search_positions[seq_pos] && !steal_positions[seq_pos] ) continue;

		// create native_frag_pose from native_pose if native are provided
		core::pose::PoseOP native_frag_pose ( new core::pose::Pose() ) ;
		if ( native_pose.size() > 0 ) {
			utility::vector1< core::Size > positions;
			core::kinematics::FoldTree fold_tree( nmer_len );
			for ( core::Size irsd=seq_pos; irsd<seq_pos+nmer_len; ++irsd ) { positions.push_back( irsd ); }
			core::pose::create_subpose( native_pose, positions, fold_tree, *native_frag_pose ); // make native_frag_pose

			// pass native to docking
			dock.setNative( native_frag_pose );
		}

		std::string frag_seq = sequence.substr( seq_pos-1, nmer_len );
		TR << "search for fragment: [" << seq_pos << "] " << frag_seq << std::endl;

		utility::vector1< core::pose::PoseOP > frags;

		if ( search_positions[seq_pos] ) {
			dock.setPassThrough( false );
			dock.setDoRefine(option[ min_pack_min ]());

			// each frag candidate at a designated position
			core::pose::Pose frag_ref;
			core::pose::make_pose_from_sequence( frag_ref, frag_seq,
				*( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID )) );
			for ( core::Size j=1; j<=frame.nr_frags(); ++j ) {
				core::pose::PoseOP frag_i(new core::pose::Pose(frag_ref));
				frame.fragment( j ).apply( *frag_i, 1, nmer_len );
				//std::string const pdbid ( frame->fragment( j ).pdbid() );
				frags.push_back( frag_i );
			}
		} else if ( steal_positions[seq_pos] ) {
			dock.setPassThrough( true );
			dock.setDoRefine( false );

			core::pose::PoseOP steal_pose ( new core::pose::Pose() );

			utility::vector1< core::Size > positions;
			core::kinematics::FoldTree fold_tree( nmer_len );
			for ( core::Size irsd=seq_pos; irsd<seq_pos+nmer_len; ++irsd ) {
				positions.push_back( initial_pose_seqmap[irsd] );
			}
			core::pose::create_subpose( initial_pose, positions, fold_tree, *steal_pose ); // make native_frag_pose

			dock.setNative( steal_pose );
			frags.push_back( steal_pose );
		}

		dock.setTag( "S_"+utility::to_string(seq_pos) );
		dock.apply_multi(frags);
	}
}


void
DockFragmentsMover::cut_from_map( core::pose::Pose const &pose ) {
	using core::scoring::electron_density::poseCoords;
	using core::scoring::electron_density::poseCoord;

	core::Size edge_trim = 5;
	core::Real mask_radius = 2; //?
	poseCoords litePose;

	for ( int i = 1; i <= (int)pose.size(); ++i ) {
		core::conformation::Residue const & rsd_i ( pose.residue(i) );
		bool skipres = ( rsd_i.aa() == core::chemical::aa_vrt );
		for ( int j = i-(int)edge_trim; j < (i+(int)edge_trim) && !skipres; ++j ) {
			if ( j>=1 && j<=(int)pose.size() && pose.fold_tree().is_cutpoint( j ) ) skipres=true;
		}

		if ( skipres ) continue;

		core::Size natoms = rsd_i.nheavyatoms();
		for ( core::Size j = 1; j <= natoms; ++j ) {
			core::chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );
			poseCoord coord_j;
			coord_j.x_ = rsd_i.xyz( j );
			coord_j.B_ = core::scoring::electron_density::getDensityMap().getEffectiveBfactor();
			coord_j.elt_ = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			litePose.push_back( coord_j );
		}
	}
	ObjexxFCL::FArray3D< double > rhoC, rhoMask;
	core::scoring::electron_density::getDensityMap().calcRhoC( litePose, 0, rhoC, rhoMask, -1, 600, mask_radius );

	// apply mask to map
	ObjexxFCL::FArray3D< float > densnew = core::scoring::electron_density::getDensityMap().get_data();
	for ( int z=1; z<=(int)densnew.u3(); z++ ) {
		for ( int y=1; y<=(int)densnew.u2(); y++ ) {
			for ( int x=1; x<=(int)densnew.u1(); x++ ) {
				densnew(x,y,z) *= (1-rhoMask(x,y,z));
			}
		}
	}
	core::scoring::electron_density::getDensityMap().set_data( densnew );

	if ( basic::options::option[ basic::options::OptionKeys::edensity::debug ]() ) {
		core::scoring::electron_density::getDensityMap().writeMRC( "trimmed.mrc" );
	}
}


///
ScoreFragmentSetMover::ScoreFragmentSetMover() {
	steepness_ = 8.0;
	overlap_width_ = 3.0;
	clash_dist_ = 2.0;
	unclosable_penalty_ = 10.0;

	gap_lengths_.push_back(6.0); // 1
	gap_lengths_.push_back(9.6); // 2
	gap_lengths_.push_back(13.0); // 3
	gap_lengths_.push_back(16.0); // 4
	gap_lengths_.push_back(19.5); // 5
	gap_lengths_.push_back(22.0); // 6
	gap_lengths_.push_back(26.0); // 7
	gap_lengths_.push_back(29.0); // 8
	gap_lengths_.push_back(33.0); // 9
	gap_lengths_.push_back(35.0); // 10
	gap_lengths_.push_back(38.0); // 11
	gap_lengths_.push_back(40.0); // 12
	gap_lengths_.push_back(43.0); // 13
	gap_lengths_.push_back(45.0); // 14
	gap_lengths_.push_back(47.0); // 15

	gap_weights_.push_back(1.0000); // 1
	gap_weights_.push_back(0.5563); // 2
	gap_weights_.push_back(0.3529); // 3
	gap_weights_.push_back(0.2907); // 4
	gap_weights_.push_back(0.1735); // 5
	gap_weights_.push_back(0.1149); // 6
	gap_weights_.push_back(0.1032); // 7
	gap_weights_.push_back(0.0580); // 8
	gap_weights_.push_back(0.0469); // 9
	gap_weights_.push_back(0.0545); // 10
	gap_weights_.push_back(0.0457); // 11
	gap_weights_.push_back(0.0504); // 12
	gap_weights_.push_back(0.0510); // 13
	gap_weights_.push_back(0.0416); // 14
	gap_weights_.push_back(0.0422); // 15
}

///
void ScoreFragmentSetMover::run() {
	// read silent file
	utility::vector1<utility::file::FileName> insilent = option[ in::file::silent ]();
	utility::vector1< utility::vector1< CAtrace > > all_frags;

	// setup symmetry
	protocols::electron_density::DensitySymmInfo symminfo( basic::options::option[ symm_type ]() );
	symminfo.detect_axes( core::scoring::electron_density::getDensityMap() );

	core::io::silent::SilentFileData sfd;
	for ( core::Size n=1; n<=insilent.size(); ++n ) {
		sfd.read_file( insilent[n] );
	}

	for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
		std::string tag = iter->decoy_tag();
		utility::vector1< std::string > fields = utility::string_split( tag, '_' );
		core::Real score_i = iter->get_energy( "dens_score" );
		core::Real rms_i = iter->get_energy( "rms" );
		core::Real rank_i = iter->get_energy( "dens_rank" );

		if ( option[ n_matches ].user() && (option[ n_matches ]+0.5) < rank_i ) continue;
		if ( option[ cheat ].user() && rms_i > 2 ) continue;

		if ( score_i<-999 ) score_i=1000.0; // hack

		runtime_assert( score_i != 0 ); //??

		TR.Debug << "Processing " << tag << std::endl;

		// NOTE: order must match what is written!
		std::string fragtag = fields[1];
		core::Size resid = atoi( fields[2].c_str() );

		// get pose
		core::pose::Pose frag;
		iter->fill_pose( frag );

		CAtrace ca_i(frag,tag, -score_i, rms_i);

		//if (mer_size==0) mer_size=ca_i.cas_.size();
		//runtime_assert( mer_size == ca_i.cas_.size() );

		if ( all_frags.size() < resid ) all_frags.resize( resid );
		if ( all_frags[resid].size() < rank_i ) all_frags[resid].resize( rank_i );
		all_frags[resid][rank_i] = ca_i;
	}

	// score: 1b first (keep them first!)
	//   make idx to tag mapping that 2b energies will use
	core::Size nres = all_frags.size();
	core::Size idx=1, type=1;
	std::ofstream outscore( option[ scorefile ]().c_str(), std::ios::out | std::ios::binary );
	for ( core::Size i_res=1; i_res<=nres; ++i_res ) {
		TR << "1b: [" << i_res << "/" << nres << "]" << "\r" << std::flush;
		core::Size nfrag_i = all_frags[i_res].size();
		for ( core::Size i_frag=1; i_frag<=nfrag_i; ++i_frag ) {
			core::Real scalefactor = 9.0 / all_frags[i_res][i_frag].cas_.size();

			core::Real dens_sc = scalefactor * all_frags[i_res][i_frag].dens_score_;
			core::Real rms = all_frags[i_res][i_frag].rms_;
			std::string tag = all_frags[i_res][i_frag].tag_;

			char tagstr[32];
			std::strncpy( tagstr, tag.c_str(), 31 );
			tagstr[31]='\0';

			outscore.write( (char*)&type , sizeof(type) );
			outscore.write( (char*)tagstr , 32*sizeof(char) ); // truncates tag at 32 chars (should be ok ....)
			outscore.write( (char*)&idx , sizeof(idx) );
			outscore.write( (char*)&i_res , sizeof(i_res) );
			outscore.write( (char*)&dens_sc , sizeof(dens_sc) );
			outscore.write( (char*)&rms , sizeof(rms) );

			all_frags[i_res][i_frag].idx_ = idx++; // make idx->tag mapping now
		}
	}
	TR << "1b: [" << nres << "/" << nres << "]" << std::endl;

	// then 2b
	type = 2;
	for ( core::Size i_res=1; i_res<=nres; ++i_res ) {
		TR << "2b: [" << i_res << "/" << nres << "]" << "\r" << std::flush;
		core::Size nfrag_i = all_frags[i_res].size();

		for ( core::Size j_res=i_res+1; j_res<=nres; ++j_res ) {
			core::Size nfrag_j = all_frags[j_res].size();
			core::Size offset = j_res-i_res;
			for ( core::Size i_frag=1; i_frag<=nfrag_i; ++i_frag ) {
				for ( core::Size j_frag=1; j_frag<=nfrag_j; ++j_frag ) {
					core::Real clash_sc = clash_score( all_frags[i_res][i_frag], all_frags[j_res][j_frag], offset, symminfo );
					core::Real overlap_sc = overlap_score( all_frags[i_res][i_frag], all_frags[j_res][j_frag], offset, symminfo );
					core::Real closability_sc = closability_score( all_frags[i_res][i_frag], all_frags[j_res][j_frag], offset, symminfo );

					if ( std::abs(clash_sc)+std::abs(overlap_sc)+std::abs(closability_sc) > 1e-4 ) {
						core::Size idxi = all_frags[i_res][i_frag].idx_;
						core::Size idxj = all_frags[j_res][j_frag].idx_;

						outscore.write( (char*)&type , sizeof(type) );
						outscore.write( (char*)&i_res , sizeof(i_res) );
						outscore.write( (char*)&idxi , sizeof(idxi) );
						outscore.write( (char*)&j_res , sizeof(j_res) );
						outscore.write( (char*)&idxj , sizeof(idxj) );
						outscore.write( (char*)&clash_sc , sizeof(clash_sc) );
						outscore.write( (char*)&overlap_sc , sizeof(overlap_sc) );
						outscore.write( (char*)&closability_sc , sizeof(closability_sc) );
					}
				}
			}
		}
	}
	TR << "2b: [" << nres << "/" << nres << "]" << std::endl;
}

core::Real
ScoreFragmentSetMover::overlap_score(
	CAtrace &pose1,
	CAtrace &pose2,
	core::Size offset,
	protocols::electron_density::DensitySymmInfo const &symminfo)
{
	core::Size mersize = pose1.cas_.size();
	core::Real overlap_ij = 0.0;

	if ( offset<mersize ) {
		for ( int i=(int)(offset+1); i<=(int)mersize; ++i ) {
			int j = i-offset;
			//core::Real dist = (pose1.cas_[i] - pose2.cas_[j]).length();
			core::Real dist = std::sqrt( symminfo.min_symm_dist2( pose1.cas_[i] , pose2.cas_[j] ) );

			// for bad violations, give full penalty
			if ( dist > 5.0 ) return 8.0; // should be (mer size)-1?

			// check for clashes
			//  if there are any clashes give a large penalty
			for ( int ii=1; ii<=(int)mersize; ++ii ) {
				for ( int jj=1; jj<=(int)mersize; ++jj ) {
					int seqsep = offset-ii+jj;
					if ( seqsep >= 5 || seqsep <= -5 ) {
						//core::Real distC = (pose1.cas_[ii] - pose2.cas_[jj]).length_squared();
						core::Real distC = symminfo.min_symm_dist2( pose1.cas_[ii] , pose2.cas_[jj] );
						if ( distC < clash_dist_*clash_dist_ ) {
							return 8.0; // should be nmer-1?
						}
					}
				}
			}

			// sigmoid
			core::Real overlap_ij_xyz = 2.0 / ( 1 + std::exp( -steepness_*(dist-overlap_width_)) ) - 1;
			//core::Real overlap_ij_xyz = std::exp( -1*dist*steepness_);
			overlap_ij += overlap_ij_xyz;
		}
	}
	return overlap_ij;
}

core::Real
ScoreFragmentSetMover::clash_score(
	CAtrace &pose1,
	CAtrace &pose2,
	core::Size offset,
	protocols::electron_density::DensitySymmInfo const &symminfo)
{
	core::Size mersize = pose1.cas_.size();
	core::Real clash_ij = 0.0;
	if ( offset >= mersize ) {
		// fragments do not overlap -- use standard penalty
		for ( int i=1; i<=(int)mersize; ++i ) {
			for ( int j=1; j<=(int)mersize; ++j ) {
				//if (i==(int)mersize && j==1 && offset==mersize) continue;
				//core::Real dist = (pose1.cas_[i] - pose2.cas_[j]).length_squared();

				core::Real dist2 = symminfo.min_symm_dist2( pose1.cas_[i] , pose2.cas_[j] );
				if ( dist2 < clash_dist_*clash_dist_ ) clash_ij += 1.0;
			}
		}
	}
	return clash_ij;
}

core::Real
ScoreFragmentSetMover::closability_score(
	CAtrace &pose1,
	CAtrace &pose2,
	core::Size offset,
	protocols::electron_density::DensitySymmInfo const &symminfo)
{
	runtime_assert( gap_lengths_.size() == gap_weights_.size() );

	core::Size mersize = pose1.cas_.size();
	core::Real close_ij = 0.0;

	if ( offset < mersize ) return 0.0;

	core::Size gap_size = offset - mersize + 1;
	if ( gap_size < gap_lengths_.size() ) {
		//core::Real dist = (pose1.cas_[mersize] - pose2.cas_[1]).length();
		core::Real dist = std::sqrt( symminfo.min_symm_dist2( pose1.cas_[mersize] , pose2.cas_[1] ) );

		if ( dist < gap_lengths_[gap_size] ) {
			close_ij -= gap_weights_[gap_size];
		} else {
			close_ij += unclosable_penalty_;
		}
	}
	return close_ij;
}



FragmentAssemblyMover::FragmentAssemblyMover() {
	overlap_wt_ = 4.0;
	clash_wt_ = 20.0;
	close_wt_ = 6.0;
	dens_wt_ = 0.85;  //fd

	if ( option[ assembly_weights ].user() ) {
		overlap_wt_ = option[ assembly_weights ]()[1];
		clash_wt_ = option[ assembly_weights ]()[2];
		close_wt_ = option[ assembly_weights ]()[3];
	}

	null_frag_ = option[ null_weight ]();
}

core::Real
FragmentAssemblyMover::score( bool verbose=false ) {
	core::Size nres = assigned_frags.size();
	core::Real score_total = 0.0;

	for ( core::Size i=1; i<=nres; ++i ) {
		core::Real clash_i=0, overlap_i=0, close_i=0, dens_i=0, score_i=0, rms_i=0;

		if ( assigned_frags[i] == 0 ) {
			score_total+=null_frag_; //null frag
			dens_i = null_frag_;
		} else {
			score_total += scores_1b[ assigned_frags[i] ];
			dens_i = scores_1b[ assigned_frags[i] ];
			score_i = scores_1b[ assigned_frags[i] ];
			rms_i = rmses[ assigned_frags[i] ];

			for ( core::Size j=1; j<=nres; ++j ) {
				if ( assigned_frags[j] != 0 ) {
					score_i += 0.5*scores_2b[ assigned_frags[i] ][ assigned_frags[j] ];
					score_total += 0.5*scores_2b[ assigned_frags[i] ][ assigned_frags[j] ];

					//clash_i += 0.5*clash2b[ assigned_frags[i] ][ assigned_frags[j] ];
					//overlap_i += 0.5*overlap2b[ assigned_frags[i] ][ assigned_frags[j] ];
					//close_i += 0.5*close2b[ assigned_frags[i] ][ assigned_frags[j] ];
				}
			}
		}
		if ( verbose ) TR << "res " << i << ": " << rms_i << " " << score_i << " = " << clash_i << "/" << overlap_i << "/" << close_i << "/" << dens_i << std::endl;
	}
	return score_total;
}


void
FragmentAssemblyMover::run( ) {
	if ( !option[ out::file::silent ].user() ) {
		TR << "Must specify an output file with -out::file::silent!" << std::endl;
		return;
	}

	// read scorefile
	std::ifstream inscore( option[ scorefile ]().c_str(), std::ios::in | std::ios::binary );
	TR << "reading scorefile" << std::endl;
	while ( !inscore.eof() ) {
		std::string tag;
		core::Real dens_sc,clash_sc,overlap_sc,closability_sc, rms_i;
		core::Size ires,jres,idxi,idxj;
		std::string itag;

		// ???  prevent too many reallocations
		pos2frags.reserve(1000);
		scores_1b.reserve(50*1000);
		rmses.reserve(50*1000);

		core::Size type;
		inscore.read( (char*)&type , sizeof(type) );

		if ( inscore.eof() ) break;

		if ( type==1 ) {
			char tagstr[32];
			inscore.read( (char*)tagstr , 32*sizeof(char) ); // assumes max tag == 32 bytes
			inscore.read( (char*)&idxi , sizeof(idxi) );
			inscore.read( (char*)&ires , sizeof(ires) );
			inscore.read( (char*)&dens_sc , sizeof(dens_sc) );
			inscore.read( (char*)&rms_i , sizeof(rms_i) );

			itag = tagstr;

			FragID frag_i( ires, itag );
			allfrags.push_back( frag_i );

			frag2idx[frag_i] = idxi;

			// have to assume fragments are not ordered
			if ( ires > pos2frags.size() ) pos2frags.resize(ires);
			pos2frags[ires].push_back( idxi );

			if ( idxi > scores_1b.size() ) scores_1b.resize( idxi );
			scores_1b[idxi] = dens_wt_ * dens_sc;

			if ( idxi > rmses.size() ) rmses.resize( idxi );
			rmses[idxi] = rms_i;

		} else if ( type==2 ) {
			// assumption: all 1body energies have been seen
			if ( scores_2b.size() == 0 ) {
				TR << "read " << allfrags.size() << " fragments" << std::endl;
				TR << "reading 2b scores" << std::endl;
				scores_2b.resize( allfrags.size(), utility::vector1<core::Real>(allfrags.size(),0) );

				// debug
				//clash2b.resize( allfrags.size(), utility::vector1<core::Real>(allfrags.size(),0) );
				//overlap2b.resize( allfrags.size(), utility::vector1<core::Real>(allfrags.size(),0) );
				//close2b.resize( allfrags.size(), utility::vector1<core::Real>(allfrags.size(),0) );
			}

			inscore.read( (char*)&ires , sizeof(ires) );
			inscore.read( (char*)&idxi , sizeof(idxi) );
			inscore.read( (char*)&jres , sizeof(jres) );
			inscore.read( (char*)&idxj , sizeof(idxj) );
			inscore.read( (char*)&clash_sc , sizeof(clash_sc) );
			inscore.read( (char*)&overlap_sc , sizeof(overlap_sc) );
			inscore.read( (char*)&closability_sc , sizeof(closability_sc) );

			scores_2b[idxi][idxj] =
				overlap_wt_ * overlap_sc + clash_wt_ * clash_sc + close_wt_ * closability_sc;

			// debug
			//clash2b[idxi][idxj] = clash_wt_ * clash_sc;
			//overlap2b[idxi][idxj] = overlap_wt_ * overlap_sc;
			//close2b[idxi][idxj] = close_wt_ * closability_sc;

			// mirror
			scores_2b[idxj][idxi] = scores_2b[idxi][idxj];
			//clash2b[idxj][idxi] = clash2b[idxi][idxj];
			//overlap2b[idxj][idxi] = overlap2b[idxi][idxj];
			//close2b[idxj][idxi] = close2b[idxi][idxj];
		}
	}

	core::Size nres = pos2frags.size();
	core::Size nstruct = option [out::nstruct]();

	// read silent files
	utility::vector1< core::pose::PoseOP > frags_sel(nres);
	core::io::silent::SilentFileData sfd;
	utility::vector1<utility::file::FileName> insilent = option[ in::file::silent ]();
	for ( core::Size n=1; n<=insilent.size(); ++n ) {
		sfd.read_file( insilent[n] );
	}


	// parameters
	core::Real sa_start_temp = 500.0, sa_end_temp=1.0;
	core::Size sa_nsteps = 200, mc_nsteps = 25*nres*option[ scale_cycles ];

	for ( int rd=1; rd<=(int)nstruct; ++rd ) {
		// MC
		// initialize
		core::Size nres = pos2frags.size();

		// add null frags, initialize to it
		for ( int i=1; i<=(int)nres; ++i ) {
			pos2frags[i].push_back( 0 );
		}
		assigned_frags.resize( nres, 0 );


		// run
		core::Real temp = sa_start_temp;
		core::Real temp_scale = std::pow( sa_end_temp/sa_start_temp, 1.0/((core::Real)sa_nsteps) );

		for ( core::Size temp_ctr=1; temp_ctr<=sa_nsteps; ++temp_ctr ) {
			for ( core::Size step=1; step<=mc_nsteps; ++step ) {
				int select_pos = numeric::random::random_range( 1, nres );

				// calculate compatibility scores for all the cadidate placements at the given pos
				utility::vector1<int> const &frag_cands = pos2frags[select_pos];
				core::Size n_frag_cands = frag_cands.size();

				utility::vector1<core::Real> frag_scores( n_frag_cands );
				core::Real best_frag_score = 1e30, prob_sum = 0.0;
				for ( int i=1; i<=(int)n_frag_cands; ++i ) {
					frag_scores[i] = 0;
					if ( frag_cands[i] == 0 ) {
						frag_scores[i] = null_frag_;
					} else {
						frag_scores[i] = scores_1b[ frag_cands[i] ];
						for ( core::Size pos=1; pos<=nres; ++pos ) {
							int assigned_fragidx = assigned_frags[pos];
							if ( assigned_fragidx != 0 && (int)pos != select_pos ) {
								frag_scores[i] += scores_2b[ frag_cands[i] ][ assigned_fragidx ];
							}
						}
					}

					best_frag_score = std::min( frag_scores[i], best_frag_score );
				}

				for ( core::Size i=1; i<=n_frag_cands; ++i ) {
					frag_scores[i] -= best_frag_score;
					frag_scores[i] = std::exp( -std::min( frag_scores[i]/temp, 100.0 ) );

					//if (option[ cheat_assem ]() && frag_cands[i] != 0 && rmses[ frag_cands[i] ]<4) {
					// frag_scores[i] *= 1e6;
					//}

					prob_sum += frag_scores[i];
				}

				// select fragment, normalize prob
				core::Real fragpicker = prob_sum*numeric::random::uniform();
				int picked = 0;
				while ( fragpicker>=0 && picked<int(n_frag_cands) ) {
					picked++;
					fragpicker -= frag_scores[picked];
				}
				assert( fragpicker < 0 );
				assigned_frags[select_pos] = frag_cands[picked];
			}
			core::Real score_total = score( false );
			core::Real nassigned = 0;
			for ( core::Size i=1; i<=nres; ++i ) { if ( assigned_frags[i] != 0 ) nassigned++; }

			TR << "Finished temp = " << temp << " : score = " << score_total << " [" << nassigned << " assigned]" << std::endl;

			temp *= temp_scale;
		}

		// report
		core::Real score_total = score( true );

		utility::vector1< std::string > tags_to_fetch;
		utility::vector1< core::Size > tag_indices;
		TR << "Total score: " << score_total << std::endl;
		for ( core::Size i=1; i<=nres; ++i ) {
			//TR << "   res " << i << "  frag " << assigned_frags[i];
			if ( assigned_frags[i] != 0 ) {
				std::ostringstream oss;
				oss << allfrags[assigned_frags[i]].tag_; // << "_" << i;
				//TR << " (" << oss.str() << ")";

				tags_to_fetch.push_back( allfrags[assigned_frags[i]].tag_ );
				tag_indices.push_back( i );
			}
			//TR << std::endl;
		}

		for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
			std::string tag = iter->decoy_tag();
			utility::vector1< std::string >::iterator tag_it = std::find( tags_to_fetch.begin(), tags_to_fetch.end(), tag );
			if ( tag_it == tags_to_fetch.end() ) continue;
			int index = std::distance (tags_to_fetch.begin(), tag_it)+1; // +1 for 1-indexing

			// do something ...
			core::pose::PoseOP frag_i ( new core::pose::Pose );
			iter->fill_pose( *frag_i );
			frags_sel[ tag_indices[index] ] = frag_i;
		}


		// write silent file of saved fragments + averaged model
		std::string outfile = option[ out::file::silent ]()+"_"+ObjexxFCL::right_string_of( rd, 4, '0' )+".silent";
		core::io::silent::SilentFileData sfd_out( outfile, false, false, "binary" ); //true to store argv in silent file
		for ( core::Size i=1; i<=tags_to_fetch.size(); ++i ) {
			core::Size idx = tag_indices[i];

			core::io::silent::BinarySilentStruct silent_stream( *(frags_sel[idx]), tags_to_fetch[i] );
			silent_stream.add_energy( "pos", idx );
			silent_stream.add_energy( "total_score", score_total );
			sfd.write_silent_struct( silent_stream, outfile );
		}

		// averaged pose ... to do
	}
}


void
SolutionRescoreMover::run() {
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	// setup symmetry
	protocols::electron_density::DensitySymmInfo symminfo( basic::options::option[ symm_type ]() );
	symminfo.detect_axes( core::scoring::electron_density::getDensityMap() );

	// read silent files
	utility::vector1<utility::file::FileName> insilent = option[ in::file::silent ]();

	// energy cut
	utility::vector1<core::Real> allscores;
	for ( core::Size n=1; n<=insilent.size(); ++n ) {
		core::io::silent::SilentFileData sfd;
		sfd.read_file( insilent[n] );

		// map resids to fragments
		std::map< core::Size, CAtrace > all_frags;

		for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
			std::string tag = iter->decoy_tag();
			utility::vector1< std::string > fields = utility::string_split( tag, '_' );

			core::Real score_i = iter->get_energy( "dens_score" );
			core::Real rms_i = iter->get_energy( "rms" );

			if ( fields.size() < 2 || fields[1] == "empty" ) { // weird bug
				//TR << "Skipping tag >" << tag << "<" << std::endl;
				continue;
			}

			core::Size resid = atoi( fields[2].c_str() );

			TR.Debug << "Processing " << tag << std::endl;
			core::pose::PoseOP frag_i (new core::pose::Pose);
			iter->fill_pose( *frag_i );

			CAtrace ca_i(*frag_i, tag, -score_i, rms_i);

			//all_frags[resid] = ca_i;
			all_frags.insert ( std::make_pair(resid,ca_i) );
		}

		// we have our structure, now score it
		std::map< core::Size, CAtrace >::iterator iter1, iter2;

		// Directly call scorefunctions
		ScoreFragmentSetMover scoring;

		core::Real overlap_wt = 4.0, clash_wt = 20.0, close_wt = 6.0, dens_wt = 0.85;
		if ( option[ assembly_weights ].user() ) {
			overlap_wt = option[ assembly_weights ]()[1];
			clash_wt = option[ assembly_weights ]()[2];
			close_wt = option[ assembly_weights ]()[3];
		}

		core::Real score=0.0, clash_score=0.0, close_score=0.0, overlap_score=0.0, dens_score=0.0;
		core::Size ncorr = 0;
		for ( iter1 = all_frags.begin(); iter1 != all_frags.end(); iter1++ ) {
			core::Real overlap_i=0.0;
			core::Real clash_i=0.0;
			core::Real close_i=0.0;
			core::Real dens_i = dens_wt * iter1->second.dens_score_;
			core::Real rms = iter1->second.rms_;

			if ( rms < 2.0 ) ncorr++;

			dens_score += dens_i;

			//TR << iter1->second.tag_ << "  " << dens_i << std::endl;
			for ( iter2 = all_frags.begin(); iter2 != all_frags.end(); iter2++ ) {
				if ( iter1->first == iter2->first ) continue;

				core::Real overlap_ij, clash_ij, close_ij;
				if ( iter2->first > iter1->first ) {
					core::Size offset = iter2->first-iter1->first;
					overlap_ij = scoring.overlap_score(iter1->second, iter2->second, offset, symminfo);
					clash_ij = scoring.clash_score(iter1->second, iter2->second, offset, symminfo);
					close_ij = scoring.closability_score(iter1->second, iter2->second, offset, symminfo);
				} else {
					core::Size offset = iter1->first-iter2->first;
					overlap_ij = scoring.overlap_score(iter2->second, iter1->second, offset, symminfo);
					clash_ij = scoring.clash_score(iter2->second, iter1->second, offset, symminfo);
					close_ij = scoring.closability_score(iter2->second, iter1->second, offset, symminfo);
				}

				overlap_i += 0.5*overlap_wt * overlap_ij;
				clash_i += 0.5*clash_wt * clash_ij;
				close_i += 0.5*close_wt * close_ij;

				//if (overlap_ij != 0.0 || clash_ij != 0.0 || close_ij != 0.0) {
				// TR << "   " << iter1->second.tag_ << " : " << iter2->second.tag_ << "  " << overlap_ij << " / " << clash_ij << " / " << close_ij << std::endl;
				//}
			}
			overlap_score += overlap_i;
			clash_score += clash_i;
			close_score += close_i;
			TR << iter1->second.tag_ << ": " << dens_i + overlap_i + clash_i + close_i << "  " << dens_i << " / " << overlap_i << " / " << clash_i << " / " << close_i << std::endl;
		}
		score = clash_score + close_score + overlap_score + dens_score;
		TR << "TOTAL: " << score << " [" << dens_score << " / " << overlap_score << " / " << clash_score << " / " << close_score << "] "
			<< ncorr << " " << all_frags.size() << std::endl;

	}
}

void
ConsensusFragmentMover::run() {
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	// read silent files
	utility::vector1<utility::file::FileName> insilent = option[ in::file::silent ]();

	// map resids to fragments
	std::map< core::Size, utility::vector1< core::pose::PoseOP > > all_frags;

	// energy cut
	utility::vector1<core::Real> allscores;
	core::io::silent::SilentFileData sfd;
	for ( core::Size n=1; n<=insilent.size(); ++n ) {
		sfd.read_file( insilent[n] );
	}

	for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
		core::Real score_i = iter->get_energy( "total_score" );
		allscores.push_back(score_i); // stored for each fragment, this is wasteful
	}
	std::sort( allscores.begin(), allscores.end() );
	core::Real score_cut = allscores[ (int)std::floor( option[energy_cut]*allscores.size() )+1 ];
	TR << "Score cut = " << score_cut << " (" << (int)std::floor( option[energy_cut]*allscores.size() ) << "/" << allscores.size() << ")" << std::endl;

	for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
		core::Real score_i = iter->get_energy( "total_score" );
		if ( score_i > score_cut ) continue;

		std::string tag = iter->decoy_tag();
		utility::vector1< std::string > fields = utility::string_split( tag, '_' );

		if ( fields.size() < 2 || fields[1] == "empty" ) { // weird bug
			TR << "Skipping tag >" << tag << "<" << std::endl;
			continue;
		}

		core::Size resid = atoi( fields[2].c_str() );

		TR.Debug << "Processing " << tag << std::endl;
		core::pose::PoseOP frag_i (new core::pose::Pose);
		iter->fill_pose( *frag_i );

		all_frags[resid].push_back( frag_i );
	}

	// find largest fragment occupancy
	core::Size maxcount = 0;
	std::map< core::Size, utility::vector1< core::pose::PoseOP > >::iterator iter;
	for ( iter = all_frags.begin(); iter != all_frags.end(); iter++ ) {
		maxcount = std::max( maxcount, iter->second.size() );
	}
	TR << "Highest frequency = " << maxcount << std::endl;

	// cutoff
	maxcount = (core::Size) std::floor( maxcount*option[consensus_frac] );

	// now average all positions with low stdev (<2A?)
	std::map< core::Size , numeric::xyzVector< core::Real > > ca_sum, c_sum, n_sum, o_sum;
	std::map< core::Size , numeric::xyzVector< core::Real > > ca_sum2;
	std::map< core::Size , core::Size > rescounts;

	for ( iter = all_frags.begin(); iter != all_frags.end(); iter++ ) {
		core::Size resid=iter->first;
		if ( iter->second.size() < maxcount ) continue;

		for ( int i=1; i<=(int)iter->second.size(); ++i ) {
			core::pose::PoseOP frag_i = iter->second[i];
			debug_assert( frag_i );
			for ( int j=1; j<=(int)frag_i->size(); ++j ) {
				if ( !frag_i->residue(j).is_protein() ) continue;
				core::Size resid_j = resid+j-1;
				if ( rescounts.find( resid_j ) == rescounts.end() ) {
					rescounts[resid_j] = 0;
					ca_sum[resid_j] = c_sum[resid_j] = n_sum[resid_j] = o_sum[resid_j]
						= numeric::xyzVector<core::Real>(0,0,0);
					ca_sum2[resid_j] = numeric::xyzVector<core::Real>(0,0,0);
				}
				numeric::xyzVector< core::Real > ca_j, c_j, n_j, o_j;
				n_j = frag_i->residue(j).atom(1).xyz();
				ca_j = frag_i->residue(j).atom(2).xyz();
				c_j = frag_i->residue(j).atom(3).xyz();
				o_j = frag_i->residue(j).atom(4).xyz();

				rescounts[resid_j]++;
				n_sum[resid_j] += n_j;
				ca_sum[resid_j] += ca_j;
				c_sum[resid_j] += c_j;
				o_sum[resid_j] += o_j;

				ca_sum2[resid_j] += numeric::xyzVector<core::Real>(
					ca_j[0]*ca_j[0], ca_j[1]*ca_j[1], ca_j[2]*ca_j[2]);
			}
		}
	}

	core::pose::Pose averaged_pose;
	std::map< core::Size , core::Size >::iterator res_iter;
	core::Size lastres = 0;

	utility::vector1< int > pdb_numbering;
	utility::vector1< char > pdb_chains;
	for ( res_iter = rescounts.begin(); res_iter != rescounts.end(); res_iter++ ) {
		core::Size resid = res_iter->first;

		n_sum[resid] /= rescounts[resid];
		ca_sum[resid] /= rescounts[resid];
		c_sum[resid] /= rescounts[resid];
		o_sum[resid] /= rescounts[resid];

		for ( int j=0; j<3; ++j ) {
			ca_sum2[resid][j] = ( ca_sum2[resid][j] / rescounts[resid] - ca_sum[resid][j]*ca_sum[resid][j] );
		}
		core::Real ca_var = ca_sum2[resid].length();
		TR << "Residue " << resid << " stdev = " << ca_var << " pose = " << ca_sum[resid] << " count = " << rescounts[resid] << std::endl;
		if ( ca_var > option[ consensus_stdev ] ) continue;  // TODO: make this a parameter

		// find a fragment containing this residue
		// this could be more efficient if we assume constant size
		core::conformation::ResidueOP res_to_add;
		bool done=false;
		for ( iter = all_frags.begin(); iter != all_frags.end() && !done; iter++ ) {
			core::Size resid_tgt=iter->first;
			for ( int i=1; i<=(int)iter->second.size() && !done; ++i ) {
				core::pose::PoseOP frag_i = iter->second[i];
				debug_assert( frag_i );
				for ( int j=1; j<=(int)frag_i->size() && !done; ++j ) {
					if ( !frag_i->residue(j).is_protein() ) continue;
					int resid_j = resid_tgt+j-1;
					if ( resid_j == (int)resid ) {
						done = true;

						core::conformation::Residue old_rsd = frag_i->residue(j);
						chemical::ResidueTypeSetCOP rsd_set( frag_i->residue_type_set_for_pose( old_rsd.type().mode() ) );
						chemical::ResidueType const & new_rsd_type(
							rsd_set->get_residue_type_with_variant_removed( old_rsd.type(),
							chemical::LOWERTERM_TRUNC_VARIANT ) );
						chemical::ResidueType const & new2_rsd_type(
							rsd_set->get_residue_type_with_variant_removed( new_rsd_type,
							chemical::UPPERTERM_TRUNC_VARIANT ) );
						chemical::ResidueType const & new3_rsd_type(
							rsd_set->get_residue_type_with_variant_removed( new2_rsd_type,
							chemical::UPPER_TERMINUS_VARIANT ) );
						chemical::ResidueType const & new4_rsd_type(
							rsd_set->get_residue_type_with_variant_removed( new3_rsd_type,
							chemical::LOWER_TERMINUS_VARIANT ) );
						chemical::ResidueType const & new5_rsd_type(
							rsd_set->get_residue_type_with_variant_removed( new4_rsd_type,
							chemical::DISULFIDE ) );

						res_to_add = conformation::ResidueFactory::create_residue( new5_rsd_type );
						res_to_add->atom("N").xyz(n_sum[resid_j]);
						res_to_add->atom("CA").xyz(ca_sum[resid_j]);
						res_to_add->atom("C").xyz(c_sum[resid_j]);
						res_to_add->atom("O").xyz(o_sum[resid_j]);
					}
				}
			}
		}

		if ( averaged_pose.size() == 0 || resid != lastres-1 ) {
			averaged_pose.append_residue_by_bond( *res_to_add );
		} else {
			core::conformation::add_variant_type_to_conformation_residue(
				averaged_pose.conformation(), chemical::UPPER_TERMINUS_VARIANT, averaged_pose.size() );
			averaged_pose.append_residue_by_jump( *res_to_add, averaged_pose.size() );
			core::conformation::add_variant_type_to_conformation_residue(
				averaged_pose.conformation(), chemical::LOWER_TERMINUS_VARIANT, averaged_pose.size() );
		}
		lastres = resid;
		pdb_numbering.push_back(resid);
		pdb_chains.push_back('A');
	}

	if ( averaged_pose.size() == 0 ) {
		TR << "NO RESIDUES IN CONSENSUS ASSIGNMENT." << std::endl;
	} else {
		// make pdbinfo, set residues
		core::pose::PDBInfoOP new_pdb_info( new core::pose::PDBInfo(averaged_pose,true) );
		new_pdb_info->set_numbering( pdb_numbering );
		new_pdb_info->set_chains( pdb_chains );
		averaged_pose.pdb_info( new_pdb_info );
		averaged_pose.pdb_info()->obsolete( false );

		// terminal variants
		if ( !averaged_pose.residue_type(1).has_variant_type( chemical::LOWER_TERMINUS_VARIANT ) ) {
			core::pose::add_variant_type_to_pose_residue( averaged_pose, chemical::LOWER_TERMINUS_VARIANT, 1 );
		}
		core::Size nres = averaged_pose.size();
		while ( !averaged_pose.residue(nres).is_protein() ) nres--;
		if ( !averaged_pose.residue_type(nres).has_variant_type( chemical::UPPER_TERMINUS_VARIANT ) ) {
			core::pose::add_variant_type_to_pose_residue( averaged_pose, chemical::UPPER_TERMINUS_VARIANT, nres );
		}

		// sidechain + Hydrogens
		id::AtomID_Mask missing( false );
		core::pose::initialize_atomid_map( missing, averaged_pose ); // dimension the missing-atom mask
		for ( Size i=1; i<= averaged_pose.size(); ++i ) {
			core::conformation::Residue const & rsd( averaged_pose.residue(i) );
			for ( Size j=5; j<= rsd.natoms(); ++j ) {
				core::id::AtomID atom_id( j, i );
				missing[ atom_id ] = true;
			}
		}
		averaged_pose.conformation().fill_missing_atoms( missing );

		// repack
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		TaskFactoryOP tf( new TaskFactory() );
		tf->push_back(TaskOperationCOP( new RestrictToRepacking() ));
		protocols::simple_moves::PackRotamersMoverOP pack_full_repack( new protocols::simple_moves::PackRotamersMover( scorefxn ) );
		pack_full_repack->task_factory(tf);
		pack_full_repack->apply( averaged_pose );

		// cartmin
		kinematics::MoveMap mm;
		mm.set_bb  ( true );
		mm.set_chi ( true );
		mm.set_jump( true );
		if ( scorefxn->get_weight( core::scoring::cart_bonded ) == 0 ) {
			scorefxn->set_weight( core::scoring::cart_bonded, 0.5 );
		}
		if ( scorefxn->get_weight( core::scoring::pro_close ) != 0 ) {
			scorefxn->set_weight( core::scoring::pro_close, 0.0 );
		}
		core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.00001, true, false, false );
		core::optimization::CartesianMinimizer minimizer;
		//minimizer.run( averaged_pose, mm, *scorefxn, options );

		averaged_pose.dump_pdb( "S_0001.pdb" );
	}
}


//////////////////////////
//////////////////////////
/// main
int main(int argc, char* argv[]) {
	try {

		NEW_OPT( mode , "What mode to run (place, score, assemble, or consensus)", "place" );
		NEW_OPT( verbose , "Be verbose?", false );

		// general options
		NEW_OPT( symm_type, "Symmetry type of system", "C1" );

		// search options
		NEW_OPT( ca_positions, "CA positions to limit search", "" );
		NEW_OPT( startmodel, "A starting model for matching", "" );
		NEW_OPT( fragfile, "Fragfile name", "" );
		NEW_OPT( num_frags, "Number of fragments to use for each residue position", 25 );
		NEW_OPT( pos, "Only calculate at these position", utility::vector1<int>() );
		NEW_OPT( designated_rank , "Only calculate at these ranks", utility::vector1<int>() );
		NEW_OPT( min_pack_min, "Min fragment hits into density?", true );
		NEW_OPT( movestep,  "The grid stepsize for translational search", 2 );
		NEW_OPT( ncyc, "Min cycles", 1 );
		NEW_OPT( min_bb, "Allow backbone minimization?", false );
		NEW_OPT( norm_scores, "Normalize scores over the map", false );
		NEW_OPT( bw, "spharm bandwidth", 16 );
		NEW_OPT( n_to_search, "how many translations to search", 4000 );
		NEW_OPT( n_filtered,  "how many solutions to take to refinement", 2000 );
		NEW_OPT( n_output, "how many solutions to output", 50 );
		NEW_OPT( clust_radius, "Cluster radius", 3.0 );
		NEW_OPT( clust_oversample, "Cluster oversampling", 2 );
		NEW_OPT( frag_dens, "Radius to use for docking fragments", 0.7 );
		NEW_OPT( frag_len, "Trim fragments to this length (0=use input length)", utility::vector1<int>() );
		NEW_OPT( native_placements , "Generate native placements only", false );
		NEW_OPT( point_radius , "Filters the grid points to be searched to be at least this distance appart", 0 );

		// scoring and assembly options
		NEW_OPT( scorefile, "Scorefile name", "fragscores.sc" );
		NEW_OPT( scale_cycles , "scale the number of cycles", 1.0 );
		NEW_OPT( assembly_weights , "Weights in assembly (dens=1): <overlap> <clash> <close>", utility::vector1<core::Real>() );
		NEW_OPT( null_weight , "Weight on null fragment", -150 );
		NEW_OPT( n_matches , "number of fragment matches", 0 );
		NEW_OPT( cheat , "Cheat in scoring/assembly", false );

		// consensus options
		NEW_OPT( delR , "del R for rotational search", 2.0 );
		NEW_OPT( consensus_frac , "Fraction of assigned models needed for consensus", 1.0 );
		NEW_OPT( consensus_stdev , "Max stdev of CAs for consensus", 2.5 );
		NEW_OPT( energy_cut , "Energy cut before consensus assembly", 0.05 );

		devel::init(argc, argv);

		// force some options
		option[ in::missing_density_to_jump ].value(true);
		option[ out::nooutput ].value(true);

		if ( option[mode]() == "place" ) {
			DockFragmentsMover().run();
		} else if ( option[mode]() == "score" ) {
			ScoreFragmentSetMover().run();
		} else if ( option[mode]() == "assemble" ) {
			FragmentAssemblyMover().run();
		} else if ( option[mode]() == "rescore" ) {
			SolutionRescoreMover().run();
		} else if ( option[mode]() == "consensus" ) {
			ConsensusFragmentMover().run();
		} else {
			TR << "Unknown mode " << option[mode]() << std::endl;
			return -1;
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

