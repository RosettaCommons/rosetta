// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/init.hh>
#include <core/types.hh>

#include <protocols/moves/Mover.hh>

#include <basic/Tracer.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>
//#include <protocols/loops/Loops.hh>
//#include <protocols/loops/util.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>
#include <protocols/docking/types.hh>
#include <protocols/docking/metrics.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/io/pdb/file_data.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/util.hh>
#include <core/sequence/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/rms_util.hh>

#include <numeric/model_quality/rms.hh>
#include <numeric/xyz.functions.hh>



#include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/format.hh>

#include <utility/io/util.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

#include <string>
#include <utility>
#include <iterator>
#include <list>

#include <math.h>

//Auto Headers
#include <protocols/simple_filters/DdgFilter.hh>


// option key includes

static basic::Tracer TR("frag_contacts");

bool
score_filter(
  const float score
)
{

  static std::vector< float > scores;
  if (scores.size() == 0) {
     scores.push_back( 9999 );
  }
  int max_rank = 0;
  if( int(scores.size()) < 500 ) {
     scores.push_back( score );
     std::sort(scores.begin(), scores.end());
     max_rank =  int( scores.size() * 0.01 );
  }
  return ( score <= scores[max_rank]  );

}


/// @brief Finds the fold tree boundaries to the left and right of <pos>.
void FindBoundaries(const core::kinematics::FoldTree& tree,
                    core::Size pos,
                    core::Size* left,
                    core::Size* right) {
  using core::Size;
  assert(left);
  assert(right);

  Size lower_cut = 0;
  Size upper_cut = tree.nres();
  for (Size i = 1; i <= tree.num_cutpoint(); ++i) {
    // find the upper boundary (inclusive)
    if (tree.cutpoint(i) >= pos && tree.cutpoint(i) < upper_cut)
      upper_cut = tree.cutpoint(i);

    // find the lower boundary (exclusive)
    if (tree.cutpoint(i) < pos && tree.cutpoint(i) > lower_cut)
      lower_cut = tree.cutpoint(i);
  }

  // set output parameters
  *left = lower_cut + 1;
  *right = upper_cut;
}

core::kinematics::Stub getxform(numeric::xyzVector<core::Real> m1,
                                numeric::xyzVector<core::Real> m2,
                                numeric::xyzVector<core::Real> m3,
                                numeric::xyzVector<core::Real> f1,
                                numeric::xyzVector<core::Real> f2,
                                numeric::xyzVector<core::Real> f3 ) {
  core::kinematics::Stub s;
  s.M = alignVectorSets(m1-m2, m3-m2, f1-f2, f3-f2);
  s.v = f2-s.M*m2;
  return s;
}

void xform_pose(core::pose::Pose& pose, const core::kinematics::Stub& s, core::Size sres=1, core::Size eres=0 ) {
  using core::Size;
  using core::id::AtomID;

  if(eres==0) eres = pose.n_residue();
  for (Size ir = sres; ir <= eres; ++ir) {
    for (Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      const AtomID aid(AtomID(ia,ir));
      pose.set_xyz(aid, s.local2global(pose.xyz(aid)));
    }
  }
}


///local mover for testing purposes
class JDtestmover : public protocols::moves::Mover {


struct fragpairdata {
  core::Size rank_a;
  core::Size rank_b;
  core::fragment::FragDataOP frag_a;
  core::fragment::FragDataOP frag_b;
};


public:
	JDtestmover()
	{}

	virtual ~JDtestmover(){};

	virtual
	void
	apply( core::pose::Pose & pose ){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::fragment;
		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;
		using namespace core;
		using namespace core::fragment;
		using namespace core::kinematics;
		using namespace std;
		using core::id::AtomID;
		using core::id::StubID;
		using core::conformation::Residue;

		int repeats = option[ OptionKeys::relax::default_repeats ]();
    std::string relaxmethod = "fast";
    if ( option[ OptionKeys::relax::quick ]() ){
      relaxmethod = "quick";
    } else if ( option[ OptionKeys::relax::thorough ]() ){
      relaxmethod = "thorough";
			repeats = 15;
    } else if ( option[ OptionKeys::relax::sequence ]() ) {
      relaxmethod = "sequence";
    } else if ( option[ OptionKeys::relax::fast ]() ) {
      relaxmethod = "fast";
    } else if ( option[ OptionKeys::relax::classic ]() ) {
      relaxmethod = "classic";
    } else if ( option[ OptionKeys::relax::mini ]() ) {
      relaxmethod = "mini";
    }

		std::stringstream relaxiter;
		relaxiter << repeats;
		relaxmethod = relaxmethod + relaxiter.str();

		std::string frag_file;
		if (option[ in::file::fragA ].user()) {
			frag_file  = option[ in::file::fragA ]();
		} else {
			utility_exit_with_message(
				"Error: in:file:fragA option required!"
			);
		}

		int MIN_SEQ_SEP = 12;
		int MIN_DIST_SQUARED = 81;
		int MIN_CONTACTS_PER_RES = 3;

		if (option[ frags::nonlocal::min_seq_sep ].user()) {
			MIN_SEQ_SEP = option[ frags::nonlocal::min_seq_sep ]();
		}

    if (option[ frags::nonlocal::ca_dist ].user()) {
			int min_dist = option[ frags::nonlocal::ca_dist ]();
      MIN_DIST_SQUARED = min_dist*min_dist;
    }

		if (option[ frags::nonlocal::min_contacts_per_res ].user()) {
			MIN_CONTACTS_PER_RES = option[ frags::nonlocal::min_contacts_per_res ]();
		}

		Size N_FRAGS = option[ frags::n_frags ]();

		std::stringstream out;
		out << frag_file << ".contact_pair_rmsds." << relaxmethod << "." << N_FRAGS << "." << MIN_CONTACTS_PER_RES << "." << MIN_SEQ_SEP << "." << MIN_DIST_SQUARED << ".filt";
		std::string outfile = out.str();

		typedef std::pair<Size,Size> pospair;

// read in a contact file - this is just the svmcon contact prediction file
/*
		std::string contacts_file;
		std::map<pospair,bool> contact_map;
		if (option[ in::file::contacts ].user()) {
			contacts_file = option[ in::file::contacts ]();
		} else {
			contacts_file = "contacts.txt";
		}
		utility::io::izstream is( contacts_file );
		int a,b;
		std::string line;
		if ( is.good() ) {
			while( getline( is, line ) ) {
				std::istringstream iss(line);
				iss >> a >> b;
				pospair contact = ( a < b ) ? std::make_pair(a,b) : std::make_pair(b,a);
				contact_map[contact]=true;
			}
		} else {
			utility_exit_with_message(
				"Error: in:file:contacts option required!"
			);
		}
*/

// get all fragment pairs that are from the same pdb that contain contacting residues
		std::multimap<pospair,fragpairdata> pairmap;
		std::multimap<pospair,fragpairdata>::iterator it;

		FragSetOP fragments = FragmentIO().read_data( frag_file ); // read the fragment file
		for ( FrameIterator ita=fragments->begin(), eita=fragments->end(); ita!=eita; ++ita ) { // iterate through frames (position X) for a
			const Frame* framea = *ita;
			Size fraga_len = framea->length();
			for ( FrameIterator itb=fragments->begin(), eitb=fragments->end(); itb!=eitb; ++itb ) { // iterate through frames (position X) for b
				const Frame* frameb = *itb;
				if (framea->start() >= frameb->start()) continue; // skip same pairs
				if (std::abs(int(framea->start()-frameb->start())) < fraga_len) continue; // skip overlapping query
				Size fragb_len = frameb->length();
				for ( Size ranka = 1; ranka <= framea->nr_frags(); ++ranka ) { // iterate through fragments (neighbors, ie. rank) in frame (position X) for a
					if (ranka > N_FRAGS) break;
					const FragData& fraga = framea->fragment( ranka );
					const std::string pdba = fraga.pdbid() + fraga.chain();






//if (fraga.pdbid() != "1wp9") continue;






					for ( Size rankb = 1; rankb <= frameb->nr_frags(); ++rankb ) { // iterate through fragments (neighbors, ie. rank) in frame (position X) for b
						bool next_bfrag = false;
						if (rankb > N_FRAGS) break;
						const FragData& fragb = frameb->fragment( rankb );
						const std::string pdbb = fragb.pdbid() + fragb.chain();

						if (pdba != pdbb) continue; // skip if not from the same pdb and chain
						if (std::abs(int(fraga.pdbpos()-fragb.pdbpos())) < fraga_len) continue; // skip overlapping fragments
						for ( Size i=1; i<=fraga_len; ++i ) { // iterate through each residue position in fragment for a
							if (next_bfrag) break;
							BBTorsionSRFD* resa = (BBTorsionSRFD*)(fraga.get_residue(i))();
							core::Vector xyza(resa->x(),resa->y(),resa->z());
							for ( Size j=1; j<=fragb_len; ++j ) { // iterate through each residue position in fragment for b
								if (std::abs(int(framea->start()+i-frameb->start()+j)) < MIN_SEQ_SEP) continue;  // skip local contacts relative to query
								if (std::abs(int(fraga.pdbpos()+i-fragb.pdbpos()+j)) < MIN_SEQ_SEP) continue; // skip local contacts relative to fragments
								BBTorsionSRFD* resb = (BBTorsionSRFD*)(fragb.get_residue(j))();
								core::Vector xyzb(resb->x(),resb->y(),resb->z());
								core::Real dist_squared = xyza.distance_squared(xyzb);
								if (dist_squared > MIN_DIST_SQUARED || dist_squared <= 0) continue;  // skip if not within min dist
								// found a contact! save the frag pair data
								pospair querypair = std::make_pair(framea->start(),frameb->start());
								fragpairdata frag_pair;
								frag_pair.rank_a = ranka;
								frag_pair.rank_b = rankb;
								frag_pair.frag_a = fraga.clone();
								frag_pair.frag_a->set_valid();
								frag_pair.frag_b = fragb.clone();
								frag_pair.frag_b->set_valid();
								pairmap.insert(std::pair<pospair,fragpairdata>(querypair,frag_pair));
								next_bfrag = true; // move on to the next b fragment
								break;
							}
						}
					}
				}
			}
		}

		// at this point we have fragment pairs that have at least one non-local contact

		// now iterate the pairmap
		for ( it=pairmap.begin() ; it != pairmap.end(); it++ ) {
			pospair querypair = (*it).first; // frame/frag start query position
			fragpairdata fragment_pair = (*it).second; // frag data
			//std::cout << "Query pair: " << (int)querypair.first << " " << (int)querypair.second << std::endl;
			//std::cout << "frag rank: " << (int)fragment_pair.rank_a << " " << (int)fragment_pair.rank_b << std::endl;
			//std::cout << *(fragment_pair.frag_a) << std::endl;
			//std::cout << *(fragment_pair.frag_b) << std::endl;

			// get pose CA coords
			std::vector< core::Vector > input_pose_coords;
			std::vector< core::Vector > relaxed_pose_coords;
			std::vector< core::Vector > fragment_coords;


			// frag a coords
			Size contact_file_contacts = 0;
			Size frag_contacts = 0;
			Size helix_rescnt = 0;
			for (Size i=1; i<=fragment_pair.frag_a->size(); i++) {
				// get pose residue CA
				Size resposa = querypair.first+i-1;
				input_pose_coords.push_back( pose.residue(resposa).xyz("CA") );
				// get fragment CA
				BBTorsionSRFD* res = (BBTorsionSRFD*)(fragment_pair.frag_a->get_residue(i))();
				core::Vector xyza(res->x(),res->y(),res->z());
				fragment_coords.push_back( xyza );
				SecstructSRFD* resss = (SecstructSRFD*)(fragment_pair.frag_a->get_residue(i))();
				if (resss->secstruct() == 'H') helix_rescnt++;
				// contact file contacts?
				for (Size j=1; j<=fragment_pair.frag_b->size(); j++) {
					Size resposb = querypair.second+j-1;
					pospair contactpair = std::make_pair(resposa,resposb);
				  if (std::abs(int(resposa-resposb)) < MIN_SEQ_SEP) continue;  // skip local contacts relative to query
					if (std::abs(int(fragment_pair.frag_a->pdbpos()+i-fragment_pair.frag_b->pdbpos()+j)) < MIN_SEQ_SEP) continue;  // skip local contacts relative to fragments
					// get distance
					BBTorsionSRFD* resb = (BBTorsionSRFD*)(fragment_pair.frag_b->get_residue(j))();
				  SecstructSRFD* resbss = (SecstructSRFD*)(fragment_pair.frag_b->get_residue(j))();
					if (i == 1 && resbss->secstruct() == 'H') helix_rescnt++;
					core::Vector xyzb(resb->x(),resb->y(),resb->z());
					core::Real dist_squared = xyza.distance_squared(xyzb);
					if (dist_squared > MIN_DIST_SQUARED || dist_squared <= 0) continue;  // skip if not within min dist or if distance is 0
//					if (contact_map.count(std::make_pair(resposa,resposb)) == 1) contact_file_contacts++;
					frag_contacts++;
				}
			}



			// are there enough contacts for the length of the fragment?
			if (frag_contacts/fragment_pair.frag_a->size() < MIN_CONTACTS_PER_RES) {
				continue;
			}
	//		if (fragment_pair.frag_a->size() > 3 && (core::Real)helix_rescnt/(core::Real)(fragment_pair.frag_a->size()+fragment_pair.frag_b->size()) > 0.4) {
	//		  if (frag_contacts/fragment_pair.frag_a->size() < 2.0 ) continue; // if there's at least 50% helix residues, reduce the min contacts per res
	//		} else {
	//		  if (fragment_pair.frag_a->size() <= 3) {
	//			  if (frag_contacts < 8) continue; // allow less contacts for 3mers
	//		  } else if (frag_contacts/fragment_pair.frag_a->size() < MIN_CONTACTS_PER_RES) {
	//				continue;
	//			}
	//		}

//		std::cout << "size: " << fragment_pair.frag_a->size() << " H: " << helix_rescnt << " contacts: " << frag_contacts << " helixfrac: " << (core::Real)helix_rescnt/(core::Real)(fragment_pair.frag_a->size()+fragment_pair.frag_b->size()) << std::endl;

// make a pose of fragment pair
      std::string sequence;
      for (Size i=1; i<=fragment_pair.frag_a->size(); i++) {
        Size resposa = querypair.first+i-1;
        sequence += pose.residue(resposa).name1();
        std::cout << "adding aa: " << pose.residue(resposa).name1() << " at pos " << resposa << std::endl;
      }
      for (Size i=1; i<=fragment_pair.frag_b->size(); i++) {
        Size resposb = querypair.second+i-1;
        sequence += pose.residue(resposb).name1();
        std::cout << "adding aa: " << pose.residue(resposb).name1() << " at pos " << resposb << std::endl;
      }
      utility::vector1<FragDataCOP> fragdatapair;
      fragdatapair.push_back(fragment_pair.frag_a);
      fragdatapair.push_back(fragment_pair.frag_b);

			core::pose::Pose fpose;
			make_pose_from_frags( fpose, sequence, fragdatapair, true );

      std::cout << "Pose reorientation complete!" << std::endl;
      std::stringstream outputpdb;
      outputpdb << "frags_pose_dump_" << fragment_pair.frag_a->size() << "_" << querypair.first << "_" << querypair.second << "_" << fragment_pair.frag_a->pdbpos() << "_" << fragment_pair.frag_b->pdbpos() << "_" <<
      fragment_pair.rank_a << "_" << fragment_pair.rank_b << ".pdb";

//			static int totaldumpcntfrag = 0;
//			if (totaldumpcntfrag < 1500) {
//				fpose.dump_pdb(outputpdb.str());
//				totaldumpcntfrag++;
//			}

			// save pose after reorientation
			core::pose::Pose startpose = fpose;

			core::scoring::ScoreFunctionOP scorefxn( NULL );
			scorefxn = core::scoring::getScoreFunction();

			// RELAX!
			protocols::relax::RelaxProtocolBaseOP protocol = protocols::relax::generate_relax_from_cmd();
      std::cout << "run mover... " << std::endl;
			core::kinematics::MoveMapOP mm = protocol->get_movemap();
			mm->set_jump(true); // set jumps movable
			protocol->set_current_tag( outputpdb.str() );
      protocol->apply( fpose );
			std::string newoutputpdb = outputpdb.str() + "." + relaxmethod + ".relaxed.pdb";

			// DDG!
			// might want to try  calc_interaction_energy in docking/metrics.cc
			protocols::simple_filters::DdgFilter ddg = protocols::simple_filters::DdgFilter( 1000, scorefxn, 1, 5);
			core::Real ddgval = ddg.compute( fpose );
			//core::Real ddgval_startpose = ddg.compute( startpose );

			// INTERACTION ENERGY!
			protocols::docking::DockJumps dockjumps = utility::tools::make_vector1<core::Size>(1);
			core::Real interaction_energy_val = protocols::docking::calc_interaction_energy( fpose, scorefxn, dockjumps);
			//core::Real interaction_energy_val_startpose = protocols::docking::calc_interaction_energy( startpose, scorefxn, dockjumps);

//			fpose.dump_pdb(newoutputpdb);

			// SCORE!
      float current_score = (*scorefxn)( fpose );

			// SCORE WITHOUT TERM RESIDUES
			utility::vector1<core::Size> excluderes;
			excluderes.push_back(1);
			excluderes.push_back(fragment_pair.frag_a->size());
			excluderes.push_back(fragment_pair.frag_a->size()+1);
			excluderes.push_back(fragment_pair.frag_a->size()+fragment_pair.frag_b->size());
			float truncated_score = scorefxn()->get_sub_score_exclude_res( fpose, excluderes );

      // calculate rms of fragment pair with relaxed fragment pair
			core::Real rms_fragment_vs_relaxed = core::scoring::CA_rmsd( fpose, startpose );

      for (Size i=1; i<=fragment_pair.frag_a->size(); i++) {
        relaxed_pose_coords.push_back( fpose.residue(i).xyz("CA") );
			}

			// frag b coords
			for (Size i=1; i<=fragment_pair.frag_b->size(); i++) {
				// get pose residue CA
				input_pose_coords.push_back( pose.residue(querypair.second+i-1).xyz("CA") );
				relaxed_pose_coords.push_back( fpose.residue(fragment_pair.frag_a->size()+i).xyz("CA") );
				// get fragment CA
				BBTorsionSRFD* res = (BBTorsionSRFD*)(fragment_pair.frag_b->get_residue(i))();
				core::Vector xyzb(res->x(),res->y(),res->z());
				fragment_coords.push_back( xyzb );
			}
			assert( input_pose_coords.size() == fragment_coords.size() );

			int const natoms = input_pose_coords.size();
			FArray2D< core::Real > p1a( 3, natoms );  // input pose
			FArray2D< core::Real > p2a( 3, natoms );  // fragment original
			FArray2D< core::Real > p3a( 3, natoms );  // relaxed fragment
			for ( int i = 0; i < natoms; ++i ) {
				for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
					p1a(k+1,i+1) = input_pose_coords[i][k];
					p2a(k+1,i+1) = fragment_coords[i][k];
					p3a(k+1,i+1) = relaxed_pose_coords[i][k];
				}
			}
			std::string pdb = fragment_pair.frag_a->pdbid() + fragment_pair.frag_a->chain();
      // calculate rms of native to original fragment pair
			core::Real rms = numeric::model_quality::rms_wrapper( natoms, p1a, p2a );
      // calculate rms of native to relaxed fragment pair
			core::Real rms_relaxed_vs_native = numeric::model_quality::rms_wrapper( natoms, p1a, p3a );

//			static int totaldumpcnt = 0;
			if (rms < 3.5) { //score_filter( current_score ) || rms_fragment_vs_relaxed < 0.5 || rms < 1.5 || rms_relaxed_vs_native < 1.5 ) {
				// only dump interesting pairs
//			  if (totaldumpcnt < 1000) {
			  fpose.dump_pdb(newoutputpdb);
			  startpose.dump_pdb(outputpdb.str());
//					totaldumpcnt++;
//				}
			}

		  utility::io::ozstream outz;
			outz.open_append(outfile.c_str());
			outz << pdb << " " << I(4,(int)querypair.first) << " " << I(4,(int)querypair.second) << " " <<
			  I(4,(int)fragment_pair.frag_a->pdbpos()) << " " << I(4,(int)fragment_pair.frag_b->pdbpos()) << " " <<
			  I(4,(int)fragment_pair.rank_a) << " " << I(4,(int)fragment_pair.rank_b) << " " <<
				I(4, frag_contacts) << " " << I(4,contact_file_contacts) << " " << I(4,contact_file_contacts+frag_contacts) << " " << F(6,3,rms) << " " << F(6,3,rms_fragment_vs_relaxed) << " " <<
				F(6,3,rms_relaxed_vs_native) << " " << F(8,3,current_score) << " " << F(8,3,truncated_score) << " " <<
				F(8,3,ddgval) << " " << F(8,3,interaction_energy_val) << //" " << F(8,3,ddgval_startpose) << " " << F(8,3,interaction_energy_val_startpose) <<
				std::endl;
		  outz.close();
		}

		return;
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return new JDtestmover;
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return false; }


	virtual std::string get_name() const { return "frag_contacts"; }

private:

};

typedef utility::pointer::owning_ptr< JDtestmover > JDtestmoverOP;




///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	core::init(argc, argv);

	JDtestmoverOP test_mover(new JDtestmover);

	protocols::jd2::JobDistributor::get_instance()->go(test_mover);

	TR << "*********************successful completion**************************" << std::endl;
}

