// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/r_frag_quality.cc
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange


#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/BBTorsionSRFD.hh>

#include <core/id/SequenceMapping.hh>
#include <core/sequence/util.hh>


#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>


#include <core/fragment/util.hh>
#include <protocols/jumping/JumpSample.hh>
#include <protocols/jumping/JumpSetup.hh>

#include <core/scoring/rms_util.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>


#include <protocols/idealize.hh>

#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>

#include <basic/options/option.hh>


#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/options/option_macros.hh>


#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/excn/Exceptions.hh>


// option key includes

//#include <basic/options/keys/pick.OptionKeys.gen.hh>


static thread_local basic::Tracer tr( "main" );


OPT_1GRP_KEY( File, pick, f )
OPT_1GRP_KEY( File, pick, a )
OPT_1GRP_KEY( File, pick, o )
OPT_1GRP_KEY( Boolean, pick, no_idealize )
OPT_1GRP_KEY( Integer, pick, size )

using namespace core;
using namespace fragment;
using namespace pose;
using namespace kinematics;
//using namespace protocols::abinitio; // for util
//using namespace protocols::jumping;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace ObjexxFCL::format;

class ThisApplication  {
public:
  ThisApplication();
  static void register_options();
  void setup();
  void run();
private:
  fragment::FragSetOP fragset_;
  Pose aligned_;
  std::string aligned_seq_;
  std::string target_seq_;
  Size frag_length_;
  core::id::SequenceMapping mapping_;
};

ThisApplication::ThisApplication()
{}

void ThisApplication::register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  NEW_OPT( pick::f ,"the structure of the aligned sequence", "aligned.pdb" );
	NEW_OPT( pick::a ,"the alignment file ", "alignment.hhpred");
  NEW_OPT( pick::o ,"the fragments for the target sequence", "fragout.dat" );
  NEW_OPT( pick::no_idealize, "idealize structure", false );
  NEW_OPT( pick::size, "size of fragments to pick", 9 );
}


void make_your_own_alignment( id::SequenceMapping &mapping ) {
  Size starts1[ 5 ] = {1, 78, 81, 136, 138 };
  Size stops1[ 5 ] = {75, 80, 133, 136, 143 };
  Size starts2[ 5 ] = {4, 79, 82, 135,136 };
  Size stops2[ 5 ] = {78, 81, 134, 135, 141 };
  for ( Size n = 0; n< 5; n++ ) {
    for ( Size s1 = starts1[ n ], s2 = starts2[ n ]; s1<=stops1[ n ]; s1++,s2++ ) {
      assert( s2 <= stops2[ n ] );
      mapping.insert_aligned_residue_safe( s1, s2 );
    }
  }
}

void ThisApplication::setup() {
  //read it
	std::string pdb_file( option[ pick::f ]() );
  core::import_pose::pose_from_pdb( aligned_, pdb_file );
  //core::util::switch_to_residue_type_set( aligned_, chemical::CENTROID );

  //idealize
  if ( !option[ pick::no_idealize ] ) {
    std::cerr << "idealization not implemented yet" << std::endl;
		protocols::IdealizeMover idealizer;
		idealizer.fast( true );
		idealizer.apply( aligned_ );
		aligned_.dump_pdb(  +"_ideal.pdb");
  } else {
    tr.Info << "stealing from structure without idealization ... hope it is already " << std::endl;
  }

  std::string const align_file( option[ pick::a ]() );
  read_alignment_file(align_file,target_seq_, aligned_seq_, mapping_);
	//  make_your_own_alignment( mapping_ );
  tr.Info << "no alignment reading implemented ,,, use constant alignment " << std::endl;
  frag_length_ = option[ pick::size ];
  fragset_ = new ConstantLengthFragSet( frag_length_ );
}

void
dump_movemap( kinematics::MoveMap const& mm, Size nres, std::ostream& out ) {
  for ( Size i = 1; i<=nres; i++ ) {
    if ( (i-1)%10 == 0 ) { out << i; continue; }
    //large numbers take several characters... skip appropriate
    if ( (i>=10) && (i-2)%10 == 0 ) { continue; }
    if ( (i>=100) && (i-3)%10 == 0 ) { continue; }
    if ( (i>=1000) && (i-4)%10 == 0 ) { continue; }
    out << ".";
  }
  out << std::endl;
  for ( Size i = 1; i<=nres; i++ ) {
    if ( mm.get_bb( i ) ) out << 'F'; //cuttable
    else out << '.';
  }
  out << std::endl;
}

void dump_sequences( utility::vector1< std::string > seq, std::ostream& out ) {
  Size nres( 0 );
  for ( Size i = 1; i<=seq.size(); i++ ){
    nres = std::max( nres, seq[ i ].size() );
  }
  for ( Size i = 1; i<=nres; i++ ) {
    if ( (i-1)%10 == 0 ) { out << i; continue; }
    //large numbers take several characters... skip appropriate
    if ( (i>=10) && (i-2)%10 == 0 ) { continue; }
    if ( (i>=100) && (i-3)%10 == 0 ) { continue; }
    if ( (i>=1000) && (i-4)%10 == 0 ) { continue; }
    out << ".";
  }
  out << std::endl;
  for ( Size s = 1; s<=seq.size(); s++ ) {
    for ( Size i = 1; i<=nres; i++ ) {
      out << seq[ s][ i ];
    }
  }
  out << std::endl;
}


void ThisApplication::run() {
  std::string pose_seq = aligned_.sequence();
  if ( pose_seq.substr(0,10) != aligned_seq_ ) {
		// tr.Error << "aligned sequence and sequence of structure are not compatible! " << std::endl;
    tr.Error << aligned_seq_ << std::endl << pose_seq << std::endl;
    //    return;
  }

  //create movemap that is true for every properly aligned residue
  //but not for phi right after gap and not for psi right before gap --> this will mean that no fragment is picked for these residues
  // this will already be enough to take care of gaps and insertions
  // after mm is setup we go throug target sequence and if a aligned pos is present we attempt to pick fragment at this position
  // if it hits a gap or insertion this picking will fail.
  kinematics::MoveMap mm;
  mm.set_bb( false );
  bool bConnected( true );
  for ( Size pos = 1; pos<=mapping_.size1(); pos++ ) {
    Size tpos = mapping_[ pos ];
    tr.Trace << "pos: " << pos << "tpos: " << tpos << std::endl;
    if ( tpos > 0 ) {
      if ( bConnected ) {
				mm.set_bb( tpos, true );
      }
      bConnected = true;
    } else {
      if ( bConnected ) {
				mm.set_bb( tpos-1, false );
      }
      bConnected = false;
    }
  }
  tr.Trace << "insertion mask " << std::endl;
  dump_movemap( mm, mapping_.size1(), tr.Trace );
	pose::set_ss_from_phipsi( aligned_);
	for ( Size pos = 1; pos <= 20; pos ++ ) {
		tr.Trace << aligned_.torsion( id::TorsionID( pos, id::BB, 1 )) << " ";
	}
	tr.Trace << std::endl;
  FragData frag( new BBTorsionSRFD, frag_length_ );
  for ( Size pos = 1; pos <= mapping_.size1(); pos++ ) {
		Size tpos = mapping_[ pos ];
    if ( tpos > 0 && frag.is_applicable( mm, tpos, tpos+frag_length_ - 1/*, aligned_*/) == frag_length_ ) {
			frag.steal( aligned_, tpos, tpos + frag_length_ - 1 );
      FrameOP frame = new Frame( pos, frag_length_ );
      bool success = frame->add_fragment( frag.clone() );
      assert( success );
      fragset_->add( frame );
    } else {
			tr.Info << "no frag at pos " << pos << std::endl;
		}
  }

  utility::io::ozstream dump_frag( option[ pick::o ]() );
  Size ct( 1 );
  for ( FrameIterator it=fragset_->begin(), eit=fragset_->end(); it!=eit; ++it, ++ct ) {
    tr.Trace << ct << std::endl;
    (*it)->show( dump_frag );
  }
}

int main( int argc, char** argv ) {
	try{
  ThisApplication::register_options();
  devel::init( argc, argv );

  ThisApplication app;
  app.setup();

  app.run();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
