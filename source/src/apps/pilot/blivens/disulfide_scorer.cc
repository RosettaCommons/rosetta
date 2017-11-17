// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file disulfide_scorer.cc
/// @brief Looks through a PDB file and outputs disulfide formation scores
/// @details Outputs scores for all detected disulfide bonds.
/// Also outputs non-bonded pairs with probability 'nds_prob'. Thus if nds_prob=0 only
/// disulfide bonds will be output, while if prob=1 the output will span
/// all n*(n-1)/2 pairs of residues in the protein.
///
/// Similarly, non-disulfide bonded cysteines are output with probability
/// 'cys_prob'. Note that cysteines are excluded from the nds_prob output
/// because of this second option.
///
/// The expected number of nonbonded pairs output is prob*n*(n-1)/2
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date Created November 2008
/// @details
/// @section cli Command Line
/// @code disulfide_scorer -s input.pdb -o output.txt -nds_prob 0.0 -cys_prob 0.0 -database db @endcode


//Utility
#include <devel/init.hh>
#include <utility/vector1.hh>
#include <fstream>

//Options
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/blivens.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//Poses
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/disulfides/CentroidDisulfidePotential.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <apps/pilot/blivens/disulfides.hh>

// Numeric headers
#include <numeric/random/random.hh>

using namespace core;
using namespace std;
using utility::vector1;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using core::conformation::Residue;

#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/jumping/StrandPairing.hh>

static basic::Tracer TR( "pilot_apps.blivens.disulfide_scorer" );

int
usage(char* msg)
{
	TR.Error << "usage: disulfide_scorer -s in.pdb -o out.txt -nds_prob prob -database db" << endl
		<< msg << endl;
	exit(1);
}

class PairScore
{
public:
	PairScore( Size positionA, Size positionB):
		filename_(""),
		posA_(positionA),
		resA_name_(""),
		pdb_posA_(0),
		chainA_(' '),
		posB_(positionB),
		resB_name_(""),
		pdb_posB_(0),
		chainB_(' '),

		//default unset values
		isDS_(' '),
		sameSecondary_(' '),
		cb_score_factor_(0.0),
		cbcb_distance_sq_(0.0),
		centroid_distance_sq_(0.0),
		cacbcb_angle_1_(0.0),
		cacbcb_angle_2_(0.0),
		cacbcbca_dihedral_(0.0),
		backbone_dihedral_(0.0),
		cbcb_distance_score_(0.0),
		centroid_distance_score_(0.0),
		cacbcb_angle_1_score_(0.0),
		cacbcb_angle_2_score_(0.0),
		cacbcbca_dihedral_score_(0.0),
		backbone_dihedral_score_(0.0)
	{}
	~PairScore() {}

	void score(pose::Pose const& pose, scoring::disulfides::CentroidDisulfidePotential const& potential) {
		Residue const& resA(pose.residue(posA_));
		Residue const& resB(pose.residue(posB_));
		resA_name_ = resA.name1();
		resB_name_ = resB.name1();

		pose::PDBInfoCOP pdbinfo = pose.pdb_info();
		pdb_posA_ = pdbinfo->number(posA_);
		chainA_   = pdbinfo->chain(posA_);
		pdb_posB_ = pdbinfo->number(posB_);
		chainB_   = pdbinfo->chain(posB_);

		potential.score_disulfide(resA,resB,
			cbcb_distance_sq_,
			centroid_distance_sq_,
			cacbcb_angle_1_,
			cacbcb_angle_2_,
			cacbcbca_dihedral_,
			backbone_dihedral_,
			cbcb_distance_score_,
			centroid_distance_score_,
			cacbcb_angle_1_score_,
			cacbcb_angle_2_score_,
			cacbcbca_dihedral_score_,
			backbone_dihedral_score_,
			cb_score_factor_
		);

		if ( isDS_ == ' ' ) {
			isDS_ = actual_disulfide(pose,posA_,posB_)?'T':'F';
		}
		if ( sameSecondary_ == ' ' ) {
			sameSecondary_ = same_secondary_structure(pose,posA_,posB_)?'T':'F';
		}

	}

	//Be sure to score() before printing
	void print(ofstream &out) {
		out << filename_ << sep;
		out << posA_ << sep;
		out << resA_name_ << sep;
		out << pdb_posA_ << sep;
		out << chainA_ << sep;
		out << posB_ << sep;
		out << resB_name_ << sep;
		out << pdb_posB_ << sep;
		out << chainB_ << sep;
		out << isDS_ << sep;
		out << std::abs(static_cast<int>(posA_-posB_+1)) << sep;
		out << sameSecondary_ << sep;
		out << cb_score_factor_ << sep;

		out << std::sqrt(cbcb_distance_sq_) << sep;
		out << std::sqrt(centroid_distance_sq_) << sep;
		out << cacbcb_angle_1_ << sep;
		out << cacbcb_angle_2_ << sep;
		out << cacbcbca_dihedral_ << sep;
		out << backbone_dihedral_ << sep;

		out << cbcb_distance_score_ << sep;
		out << centroid_distance_score_ << sep;
		out << cacbcb_angle_1_score_ << sep;
		out << cacbcb_angle_2_score_ << sep;
		out << cacbcbca_dihedral_score_ << sep;
		out << backbone_dihedral_score_ << sep;
		out <<  cbcb_distance_score_ +
			centroid_distance_score_ +
			cacbcb_angle_1_score_*.5 +
			cacbcb_angle_2_score_*.5 +
			cacbcbca_dihedral_score_ +
			backbone_dihedral_score_ << endl;
	}

	static void printHeader(ofstream &out) {
		out << "filename" << sep;
		out << "posL" << sep;
		out << "resL" << sep;
		out << "pdb_posL" << sep;
		out << "pdb_chnL" << sep;
		out << "posU" << sep;
		out << "resU" << sep;
		out << "pdb_posU" << sep;
		out << "pdb_chnU" << sep;
		out << "isDS" << sep;
		out << "Seq_dist" << sep;
		out << "sameSecondary" << sep;

		out << "Score_factor" << sep;
		out << "Cb_dist" << sep;
		out << "Cen_dist" << sep;
		out << "Cb_ang1" << sep;
		out << "Cb_ang2" << sep;
		out << "Cb_dihed" << sep;
		out << "bb_dihed" << sep;

		out << "Cb_dist_scr" << sep;
		out << "Cen_dist_scr" << sep;
		out << "Cb_ang1_scr" << sep;
		out << "Cb_ang2_scr" << sep;
		out << "Cb_dihed_scr" << sep;
		out << "bb_dihed_scr" << sep;
		out << "total" << endl;
	}


	void setFilename(string fn) { filename_ = fn; }
	void setDS(bool isDS) { isDS_ = isDS?'T':'F'; }

private:
	static const string sep;

	//Metadata
	string filename_;
	Size posA_;
	string resA_name_;
	int pdb_posA_;
	char chainA_;
	Size posB_;
	string resB_name_;
	int pdb_posB_;
	char chainB_;
	char isDS_;// rather than bool store 'T', 'F', or ' ' (unset)
	char sameSecondary_;// rather than bool store 'T', 'F', or ' ' (unset)

	//Scoring Info
	Real cb_score_factor_;

	Real cbcb_distance_sq_;
	Real centroid_distance_sq_;
	Real cacbcb_angle_1_;
	Real cacbcb_angle_2_;
	Real cacbcbca_dihedral_;
	Real backbone_dihedral_;

	Energy cbcb_distance_score_;
	Energy centroid_distance_score_;
	Energy cacbcb_angle_1_score_;
	Energy cacbcb_angle_2_score_;
	Energy cacbcbca_dihedral_score_;
	Energy backbone_dihedral_score_;
};
const string PairScore::sep = "\t";
//end PairScore

int main( int argc, char * argv [] )
{
	try {

		//init options system
		option.add(blivens::disulfide_scorer::nds_prob,"Probability of outputing a bond");
		devel::init(argc, argv);

		vector1<string> pdbs;
		if ( option[ in::file::s ].user() || option[in::file::l].user() ) {
			pdbs = basic::options::start_files();
		} else return usage("No in file given: Use -s or -l to designate pdb files to search for disulfides");

		string outfile;
		ofstream out;
		if ( option[ out::file::o ].user() ) {
			outfile = option[ out::file::o ]();
			out.open(outfile.c_str());
		} else return usage("No out file given: Use -o to designate an output file");


		Real prob = option[blivens::disulfide_scorer::nds_prob]();
		Real cys_prob = option[blivens::disulfide_scorer::cys_prob]();

		if ( prob<0.0 || 1.0<prob ) {
			prob = 0.0;
			TR.Warning << "nds_prob out of range; using "<<prob <<std::endl;
		}
		if ( cys_prob<0.0 || 1.0<cys_prob ) {
			cys_prob = prob;
			TR.Warning << "cys_prob out of range; using "<<cys_prob <<std::endl;
		}


		PairScore::printHeader(out);

		core::scoring::disulfides::CentroidDisulfidePotential potential;

		for ( vector1<string>::const_iterator infile_it = pdbs.begin();   infile_it != pdbs.end(); ++infile_it ) {

			pose::Pose pose;
			core::import_pose::centroid_pose_from_pdb( pose, *infile_it , core::import_pose::PDB_file);

			//Assign secodary structure
			core::scoring::dssp::Dssp dssp(pose);
			dssp.insert_ss_into_pose(pose);

			TR << pose.secstruct() << endl;


			//for each residue pair
			for ( Size i(1); i<= pose.size()-1; ++i ) {
				char i_name = pose.residue_type(i).name1();
				for ( Size j(i+1); j<= pose.size(); ++j ) {
					bool isDS = actual_disulfide(pose,i,j);
					char j_name = pose.residue_type(j).name1();

					double coin = numeric::random::uniform();

					//Decide whether to score this pair
					//Score all disulfides
					//Score everything (except Gly) that passes cointoss
					//Score Cys pairs that pass cointoss
					if ( isDS ||
							coin < prob && i_name != 'G' && j_name != 'G' ||
							coin < cys_prob && i_name == 'C' && j_name == 'C' ) {
						//score the pair
						PairScore scores(i,j);
						scores.setFilename(*infile_it);
						scores.setDS(isDS);
						scores.score(pose, potential);
						scores.print(out);
					}
				}
			}//end residue pair
		}

		out.close();

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // end main

