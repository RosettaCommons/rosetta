// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SwapAlignedFragments might be a better name
/// @author Eva-Maria Strauch (evas01@u.washington.edu)
/// @brief	mover for pulling fragments into a pose so that they can be aligned onto stubs or other segements

// Unit headers
#include <protocols/seeded_abinitio/PlaceFragments.hh>
#include <protocols/seeded_abinitio/PlaceFragmentsCreator.hh>

//#include <protocols/protein_interface_design/util.hh>
#include <utility/string_util.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/vector1.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/picking_old/vall/util.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/util.hh>

#include <protocols/forge/build/Interval.hh>

// unit headers
#include <protocols/jd2/parser/BluePrint.hh>
//#include <protocols/fldsgn/BluePrintBDR.hh>
//#include <protocols/fldsgn/BluePrintBDRCreator.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/methods/pose_mod.hh>

#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/picking_old/FragmentLibraryManager.hh>
/////#include <protocols/jd2/parser/BluePrint.fwd.hh>


#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

//#include <protocols/evaluation/Align_RmsdEvaluator.hh>
#include <core/scoring/rms_util.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>


// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/tag/Tag.hh>


namespace protocols {
namespace seeded_abinitio {

using namespace core;
using namespace scoring::constraints;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.seeded_abinitio.PlaceFragments" );

	std::string
	PlaceFragmentsCreator::keyname() const
	{
		return PlaceFragmentsCreator::mover_name();
	}

	protocols::moves::MoverOP
	PlaceFragmentsCreator::create_mover() const {
		return new PlaceFragments;
	}

	std::string
	PlaceFragmentsCreator::mover_name()
	{
		return "PlaceFragments";
	}

	PlaceFragments::~PlaceFragments() {}

	PlaceFragments::PlaceFragments() :
	protocols::moves::Mover( PlaceFragmentsCreator::mover_name() )
	{
		frags_onflight_ = true;
	}




void
PlaceFragments::initialize_fragments( const core::fragment::FragSetOP & fragments ) {
	using namespace core::fragment;
	assert(fragments);
	fragments_ = fragments;

	core::fragment::FrameIterator i;
	for (i = fragments_->begin(); i != fragments_->end(); ++i) {
		Size position = (*i)->start();
		library_[position] = **i;
	}
}


void
PlaceFragments::create_fragments( core::pose::Pose & pose, core::Size insert_start, core::Size insert_stop ) {

	// pick from vall based on template SS, and target sequence if specified
	// randomly pick a sequence from the input that fits to the specified fragment size, this of course only makes sense
    // if the fragment size is shorter than the specified secondary structure or amino acid stretch

    Size seq_start = numeric::random::random_range( 0, ss_.size() - fsize_ );
    std::string ss_sub = "" ;
    std::string aa_sub = "" ;

    for( Size i = seq_start ; i < ss_.size(); ++ i )
    	 ss_sub+= ss_[i];

    std::cout << "picking fragments for secondary structure "<< ss_sub << std::endl;


	for ( core::Size j=insert_start; j<= insert_stop - fsize_-1 ; ++j ) {
		//std::string ss_sub = tgt_ss.substr( j-1, fsize_ );
		//std::string aa_sub = tgt_seq.substr( j-1, fsize_ );
		core::fragment::FrameOP frame = new core::fragment::Frame( j, fsize_ );

		if( use_seq_ )
			frame->add_fragment( core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa( ss_sub, aa_sub,
															nfrags_, true, core::fragment::IndependentBBTorsionSRFD() ) );

		else
			frame->add_fragment( core::fragment::picking_old::vall::pick_fragments_by_ss( ss_sub, nfrags_, true,
															core::fragment::IndependentBBTorsionSRFD() ) );

		fragments_->add( frame );
	}
}


void
//PlaceFragments::apply_frag( core::pose::Pose &pose, core::pose::Pose &templ, protocols::loops::Loop &frag, bool superpose) {
PlaceFragments::apply_frag( core::pose::Pose &pose, core::fragment::Frame &frame ){//core::pose::Pose &templ, protocols::loops::Loop &frag, bool superpose) {

	/// superimpose fragment

	numeric::xyzMatrix< core::Real > R;
	numeric::xyzVector< core::Real > preT(0,0,0), postT(0,0,0);
	R.xx() = R.yy() = R.zz() = 1;
	R.xy() = R.yx() = R.zx() = R.zy() = R.yz() = R.xz() = 0;

  core::Size fstart = frame.start();
  core::Size len = frame.length();
  core::pose::Pose pose_copy = pose;

	// in case that the overlap is ever more than 1 residues
	//int cartfrag_overlap = aln_len_;
	runtime_assert( cartfrag_overlap_>=1 &&  cartfrag_overlap_<= len/2 + 1);//cartfrag_overlap_<=len/2);
	core::Size nres = pose.total_residue();


	// pick a random fragment
	core::Size toget = numeric::random::random_range( 1, frame.nr_frags() );
	frame.apply( toget, pose_copy );

	core::Size aln_len = std::min( (core::Size)9999, len );   //can change 9999 to some max alignment sublength
	core::Size aln_start = numeric::random::random_range(frag.start(), len-aln_len+frag.start() );

	// don't try to align really short frags
	if (len > 2) {
		ObjexxFCL::FArray2D< core::Real > final_coords( 3, 4*aln_len );
		ObjexxFCL::FArray2D< core::Real > init_coords( 3, 4*aln_len );

		for (int ii=0; ii<(int)aln_len; ++ii) {
			int i=aln_start+ii;
			numeric::xyzVector< core::Real > x_1 = pose_copy.residue(i).atom(" C  ").xyz();
			numeric::xyzVector< core::Real > x_2 = pose_copy.residue(i).atom(" O  ").xyz();
			numeric::xyzVector< core::Real > x_3 = pose_copy.residue(i).atom(" CA ").xyz();
			numeric::xyzVector< core::Real > x_4 = pose_copy.residue(i).atom(" N  ").xyz();
			preT += x_1+x_2+x_3+x_4;

			numeric::xyzVector< core::Real > y_1 = pose.residue(templ.pdb_info()->number(i)).atom(" C  ").xyz();
			numeric::xyzVector< core::Real > y_2 = pose.residue(templ.pdb_info()->number(i)).atom(" O  ").xyz();
			numeric::xyzVector< core::Real > y_3 = pose.residue(templ.pdb_info()->number(i)).atom(" CA ").xyz();
			numeric::xyzVector< core::Real > y_4 = pose.residue(templ.pdb_info()->number(i)).atom(" N  ").xyz();
			postT += y_1+y_2+y_3+y_4;

			for (int j=0; j<3; ++j) {
				init_coords(j+1,4*ii+1) = x_1[j];
				init_coords(j+1,4*ii+2) = x_2[j];
				init_coords(j+1,4*ii+3) = x_3[j];
				init_coords(j+1,4*ii+4) = x_4[j];
				final_coords(j+1,4*ii+1) = y_1[j];
				final_coords(j+1,4*ii+2) = y_2[j];
				final_coords(j+1,4*ii+3) = y_3[j];
				final_coords(j+1,4*ii+4) = y_4[j];
			}
		}
		preT /= 4*len;
		postT /= 4*len;
		for (int i=1; i<=(int)4*len; ++i) {
			for ( int j=0; j<3; ++j ) {
				init_coords(j+1,i) -= preT[j];
				final_coords(j+1,i) -= postT[j];
			}
		}

		// get optimal superposition
		// rotate >init< to >final<
		ObjexxFCL::FArray1D< numeric::Real > ww( 4*len, 1.0 );
		ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
		numeric::Real ctx;

		numeric::model_quality::findUU( init_coords, final_coords, ww, 4*len, uu, ctx );
		R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
		R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
		R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );
	}


/// xyz copy fragment to pose

for (int i=frag.start(); i<=frag.stop(); ++i) {
	for (int j=1; j<=pose_copy.residue(i).natoms(); ++j) {
		core::id::AtomID src(j,i), tgt(j, pose_copy.pdb_info()->number(i));
		pose.set_xyz( tgt, postT + (R*(pose_copy.xyz( src )-preT)) );
	}
}
}


void
PlaceFragments::apply_frame( 	core::pose::Pose & pose,
								core::fragment::Frame &frame,
								int aln_len,
								core::Size seq_start,
								core::Size max_frag_len  ) {



	core::Size start = frame.start() + seq_start,len = frame.length();
	bool nterm = (start == seq_start ); //seq_position );
	bool cterm = (start == pose.total_residue()-max_frag_len );

	// insert frag
	core::pose::Pose pose_copy = pose;

	//compare two, use 4 atoms from alnlen defined residues, 1 is the initialized number
	ObjexxFCL::FArray1D< numeric::Real > ww( 2*4*aln_len, 1.0 );
	ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
	numeric::xyzVector< core::Real > com1(0,0,0), com2(0,0,0);

	for (int tries = 0; tries<100; ++tries) {
		ww = 1.0;
		uu = 0.0;
		com1 = numeric::xyzVector< core::Real >(0,0,0);
		com2 = numeric::xyzVector< core::Real >(0,0,0);

		// grab coords
		ObjexxFCL::FArray2D< core::Real > init_coords( 3, 2*4*aln_len );
		for (int ii=-aln_len; ii<aln_len; ++ii) {
			int i = (ii>=0) ? (nterm?len-ii-1:ii) : (cterm?-ii-1:len+ii);

			numeric::xyzVector< core::Real > x_1 = pose.residue(start+i).atom(" C  ").xyz();
			numeric::xyzVector< core::Real > x_2 = pose.residue(start+i).atom(" O  ").xyz();
			numeric::xyzVector< core::Real > x_3 = pose.residue(start+i).atom(" CA ").xyz();
			numeric::xyzVector< core::Real > x_4 = pose.residue(start+i).atom(" N  ").xyz();
			com1 += x_1+x_2+x_3+x_4;

			for (int j=0; j<3; ++j) {
				init_coords(j+1,4*(ii+aln_len)+1) = x_1[j];
				init_coords(j+1,4*(ii+aln_len)+2) = x_2[j];
				init_coords(j+1,4*(ii+aln_len)+3) = x_3[j];
				init_coords(j+1,4*(ii+aln_len)+4) = x_4[j];
			}
		}
		com1 /= 2.0*4.0*aln_len;
		for (int ii=0; ii<2*4*aln_len; ++ii) {
			for ( int j=0; j<3; ++j ) init_coords(j+1,ii+1) -= com1[j];
		}

		core::Size toget = numeric::random::random_range( 1, frame.nr_frags() );
		frame.apply( toget, pose_copy );

		// grab new coords
		ObjexxFCL::FArray2D< core::Real > final_coords( 3, 2*4*aln_len );
		for (int ii=-aln_len; ii<aln_len; ++ii) {
			int i = (ii>=0) ? (nterm?len-ii-1:ii) : (cterm?-ii-1:len+ii);
			numeric::xyzVector< core::Real > x_1 = pose_copy.residue(start+i).atom(" C  ").xyz();
			numeric::xyzVector< core::Real > x_2 = pose_copy.residue(start+i).atom(" O  ").xyz();
			numeric::xyzVector< core::Real > x_3 = pose_copy.residue(start+i).atom(" CA ").xyz();
			numeric::xyzVector< core::Real > x_4 = pose_copy.residue(start+i).atom(" N  ").xyz();
			com2 += x_1+x_2+x_3+x_4;
			for (int j=0; j<3; ++j) {
				final_coords(j+1,4*(ii+aln_len)+1) = x_1[j];
				final_coords(j+1,4*(ii+aln_len)+2) = x_2[j];
				final_coords(j+1,4*(ii+aln_len)+3) = x_3[j];
				final_coords(j+1,4*(ii+aln_len)+4) = x_4[j];
			}
		}
		com2 /= 2.0*4.0*aln_len;
		for (int ii=0; ii<2*4*aln_len; ++ii) {
			for ( int j=0; j<3; ++j ) final_coords(j+1,ii+1) -= com2[j];
		}

		// get optimal superposition
		// rotate >final< to >init<
		numeric::Real ctx;
		float rms;

		numeric::model_quality::findUU( final_coords, init_coords, ww, 2*4*aln_len, uu, ctx );
		numeric::model_quality::calc_rms_fast( rms, final_coords, init_coords, ww, 2*4*aln_len, ctx );

		std::cout << "try " << tries << " rms " << rms << std::endl;

		if (rms < 0.5) break;
		if (tries >= 20 && rms < 1) break;
		if (tries >= 40 && rms < 2) break;
		if (tries >= 60 && rms < 3) break;
	}
	numeric::xyzMatrix< core::Real > R;
	R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
	R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
	R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );

	// apply rotation to ALL atoms
	// x_i' <- = R*x_i + com1;
	for ( Size i = 0; i < len; ++i ) {
		for ( Size j = 1; j <= pose.residue_type(start+i).natoms(); ++j ) {
			core::id::AtomID id( j, start+i );
			pose.set_xyz( id, R * ( pose_copy.xyz(id) - com2) + com1 );
		}
	}
}

///adjustment since parse time specified residues are different numbered than run time residues
utility::vector1< core::Size >
parse_residues( pose::Pose & pose,
					std::string resid ){

	utility::vector1< std::string > const design_keys( utility::string_split( resid, ',' ) );
	utility::vector1< core::Size > res;

	foreach( std::string const key, design_keys ){
		core::Size const resnum( core::pose::parse_resnum( key, pose ));
		res.push_back( resnum);
		TR<<"parsed: "<<key<<std::endl;
	}
	return res;
}

void
PlaceFragments::apply( pose::Pose & pose ){

	/// overall scheme:
	/// grow around stub (can be done outside )
	/// decide on fragment length (or a window of it)
	/// make fragments
	/// align and copy fragment coordinates
	/// modify with mover
	/// filter

	// need to keep track of starting position or loop for alignement, start with 1 res for simplicity

	//protocols::loops::Loops parse_seeds( );

	utility::vector1<core::Size> parsed_residues( parse_residues( input_stubs ));
	Size stub = parsed_residues[1];
	std::cout << "stub: " << stub <<std::endl;

	/// 1. define fragment insert/starting position and get and assign fragments

	// get length of the last chain
	Size const num_chains( pose.conformation().num_chains() );
	Size chainB_len = pose.split_by_chain( num_chains ).total_residue();

	Size insert_start = pose.chain_begin( num_chains );
	Size insert_stop = chainB_len - fsize_ + 1;
	runtime_assert( insert_start <= insert_stop );

	create_fragments( pose, insert_start, insert_stop);

  	// map resids to frames
  	core::Size insert_frags_pos = fragments_->min_pos();
  	std::cout << "start frags = fragments->min_pos:"<< insert_frags_pos <<
  	"\nmax_pos " << fragments_->max_pos() <<
  	"\nfragset for positions = nr_frames " << fragments_->nr_frames() << std::endl;

  	for (core::fragment::ConstFrameIterator i = fragments_->begin(); i != fragments_->end(); ++i){
	  core::Size position = (*i)->start();
	  std::cout << "position after iterator: " << position << std::endl;
	  library_[position] =  **i;
  	}

    /// 2. alignment

	  //int select_position = numeric::random::random_range(1,3); //4);
	  //core::Size max_pos = max_poses[ select_position ];
	  //int select_position = numeric::random::random_range(insert_start,insert_stop);

	  // insert more complex alignment other than a single stub (eg. loop for segement) here...
	  // for now just simple single
	  int in_position = stub_pos;

	  // select random pos around the middle depending on fragment size
	  core::Size insert_pos = max_pos - numeric::random::random_range(fsize_/2 -1, fsize_/2); //+1); /////////?
	  std::cout << "insert position before apply frags : " << insert_pos << std::endl;

	  //insert_pos = std::min( insert_pos, nres - big_-1);
	  insert_pos = std::min( insert_pos, insert_stop - fsize_-1);
	  insert_pos = std::max( (int)insert_pos, (int)insert_start);

	  // for debugging of frames
	  for (boost::unordered_map<core::Size, core::fragment::Frame>::iterator iter = library_.begin() ; iter != library_.end() ; ++iter){
	  	std::cout << " frame 	" << (*iter).first << " len: " << library_[insert_pos].length() << std::endl;
	  }

	  if (library_.find(insert_pos) != library_.end()){
		  apply_frag (pose, library_[insert_pos]);
		  std::cout << "applying fragments on position: " << insert_pos << std::endl;
	  }

}

std::string
PlaceFragments::get_name() const {
		return PlaceFragmentsCreator::mover_name();
}

void
PlaceFragments::parse_my_tag( 	TagCOP const tag,
								basic::datacache::DataMap & /*data*/,
								protocols::filters::Filters_map const & filters,
								Movers_map const & movers,
								Pose const & pose){

	using core::fragment::FragmentIO;
	using core::fragment::FragSetOP;
	using std::string;
	using namespace filters;

	if( tag->hasOption ( "fragments" )){
		string fragments_file = tag->getOption<string>("fragments");
		FragSetOP fragments = FragmentIO().read_data(fragments_file);
		initialize_fragments(fragments);
		frags_onflight_ = false;
	}

	fsize_ = tag->getOption < core::Size >("frag_length", 6 );
	ss_ = tag->getOption<std::string>( "secstr", "" );
	if( tag->!hasOption("secstr") || !hasOption("fragments") )
	   throw utility::excn::EXCN_RosettaScriptsOption("either need to specify secondary structure or supply fragements!!");

    // option for amino acid sequence not yet specified
    use_seq_ = false

	nfrags_= tag->getOption<core::Size>( "nfrags", 50 );
	cartfrag_overlap_ = tag->getOption < int >("aln_len", 1);

	/// to simplify for now, just take in one residue
    if( tag->hasOption("stubs") )
       input_stubs_ = tag->getOption< std::string >( "stubs" );

    else{
		/// read input seeds
    	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );
    	foreach( TagCOP const btag, branch_tags ){

	   		if( btag->getName() == "Seeds" ) { //need an assertion for the presence of these or at least for the option file

		   		std::string const beginS( btag->getOption<std::string>( "begin" ) );
		   		std::string const endS( btag->getOption<std::string>( "end" ) );
		   		std::pair <std::string,std::string> seedpair;
		   		seedpair.first 	= beginS;
		   		TR.Debug <<"parsing seeds: " << beginS << " " <<endS <<std::endl;
		   		seedpair.second = endS;
		   		seed_vector_.push_back( seedpair );
 	   		}//end seeds
            else {
                throw utility::excn::EXCN_RosettaScriptsOption("need to either specify a stub residue or a seed/segment");
            }
    	}//end b-tags
	}

    /// get movers
	Movers_map::const_iterator find_mover( movers.find( mover_name ) );
    std::string const mover_name( tag->getOption< std::string >( "mover", "null" ) );
    protocols::moves::Movers_map::const_iterator mover_it( movers.find( mover_name ) );
    if( mover_it == movers.end() )
       throw utility::excn::EXCN_RosettaScriptsOption( "mover "+ mover_name+" not found" );
	mover_ = mover_it->second ;

    /// get filters
    Filters_map::const_iterator find_filter( filters.find( filter_name ));
    if( find_filter == filters.end() ) {
	   TR<<"WARNING WARNING!!! filter not found in map. skipping: \n"<<tag<<"defaulting to truefilter "<<std::endl;
	   //runtime_assert( find_filter == filters.end() );
    }
	else
	   find_filter_ = new protocols::filters::TrueFilter;

    filter_ = find_filter->second->clone();

    TR << "with mover \"" << mover_name << "\" and filter \"" << filter_name << std::endl ;
    TR.flush();
}



}//seeded_abinitio
}//protocol
