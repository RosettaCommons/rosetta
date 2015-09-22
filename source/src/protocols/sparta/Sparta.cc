// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
/*                                                                            */
/*                           ----  SPARTA   ----                              */
/*     Shifts Prediction from Analogue of Residue type and Torsion Angle      */
/*                           Yang Shen and Ad Bax                             */
/*                    J. Biomol. NMR, xx, xxx-xxx (2010)                      */
/*                 NIH, NIDDK, Laboratory of Chemical Physics                 */
/*                     version, 1.00 (build 2010.0607.00)                     */
/*                                                                            */
/*                      for any problem, please contact                       */
/*                          shenyang@niddk.nih.gov                            */
/*                                                                            */
/******************************************************************************/
///  modified for use inside CS-Rosetta  by Oliver Lange
///
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


///Problems found during porting:
// ANN is loaded with residues r1+1 .. rN-1
// but PRED_SUM is loaded from r1 .. rN --> first residue can be uninitialized (came only up in MSPARTA_PI runs for some reason... )
// -- HA3 -> means aN.size -> 9 need to have extra memory alloacated
// string functions were weird ... potential memory problems... replaced in util.cc where fishy...


/// @author Oliver Lange

// Unit Headers
#include <protocols/sparta/Sparta.hh>
#include <protocols/sparta/SpartaUtil.hh>

#include <core/pose/Pose.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>


//// C++ headers
#include <cstdlib>
#include <string>
#include <algorithm>    // std::min

#ifdef WIN32
#include <direct.h>
#include <ctime>
#else
#endif


#include <protocols/sparta/constants.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/database/open.hh>

#include <numeric/NumericTraits.hh>
#include <boost/algorithm/string/erase.hpp>


namespace protocols {
namespace sparta {

static THREAD_LOCAL basic::Tracer tr( "protocols.sparta" );

bool protocols::sparta::Sparta::options_registered_( false );

void protocols::sparta::Sparta::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( options_registered_ ) return;
	options_registered_ = true;
}


Sparta::~Sparta() {
	//deallocate_arrays();
}

Sparta::SpartaLib* Sparta::lib_instance_;

using namespace core;
using namespace std;
Sparta::SpartaLib::SpartaLib() {
	//string libvar;
	if ( getenv( "SPARTA_DIR" ) == NULL ) {
		SPARTA_DIR = ".";
	} else {
		SPARTA_DIR = getenv( "SPARTA_DIR" );
	}

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	slash_char = "/"; //default Unix
	SPARTA_DIR=basic::database::full_name( "external/SPARTA+" );
	if ( option[ OptionKeys::evaluation::sparta_dir ].user() ) SPARTA_DIR=option[ OptionKeys::evaluation::sparta_dir ]();

	if ( SPARTA_DIR.find("/") != string::npos ) slash_char = "/"; // unix
	else if ( SPARTA_DIR.find("\\") != string::npos ) slash_char = "\\"; // Windows
	else SPARTA_DIR = ".";

	string temp;
	if ( getenv( "PATH" ) != NULL ) {
		temp = getenv( "PATH" );
		if ( temp.find("/") != string::npos ) slash_char = "/"; // unix
		else if ( temp.find("\\") != string::npos ) slash_char = "\\"; // Windows
	};

	aN[1]="N"; aN[2]="HA"; aN[3]="C"; aN[4]="CA"; aN[5]="CB"; aN[6]="HN"; //aN[7]="H3";
	/*
	aN_ALL[1]="N"; aN_ALL[2]="HA"; aN_ALL[3]="C"; aN_ALL[4]="CA"; aN_ALL[5]="CB"; aN_ALL[6]="HN"; //aN_ALL[7]="H3";
	*/

	//if ( option[ OptionKeys::sparta::dir ].user() ) SPARTA_DIR = option[ OptionKeys::sparta::dir ]();

	init();
}


Sparta::Sparta( std::string const & cs_file ) :
	REF_CS_Tab( cs_file ),
	bCreateOutput_( false )
{
	if ( !lib_instance_ ) lib_instance_ = new SpartaLib;
	refCSFileName=cs_file;
}

void Sparta::SpartaLib::setup_for_scoring( core::pose::Pose const & pose ) {
	inName = "INTERNAL";
	inPDB.loadPDB( pose );
	residList = inPDB.residListOne;

	r1 = inPDB.r1;
	rN = inPDB.rN;

	if ( firstRes < r1 ) firstRes = r1;
	if ( lastRes  < r1 ) lastRes  = r1;

	if ( firstRes > rN ) firstRes = rN;
	if ( lastRes  > rN ) lastRes  = rN;

	if ( firstRes > lastRes ) {
		int itemp = firstRes;
		firstRes  = lastRes;
		lastRes   = itemp;
	}
	tr.Info << "run ANN Sparta for pose with " << rN-r1+1 << " residues " << std::endl;
}


void Sparta::check_pose( core::pose::Pose const & pose ) {
	std::string pose_sequence( pose.sequence() );
	std::string tab_sequence( REF_CS_Tab.getData("SEQUENCE") );
	//std::replace(tab_sequence.begin(),tab_sequence.end(),' ','');
	boost::erase_all(tab_sequence," ");
	core::Size compare_len( min( pose_sequence.length(),tab_sequence.length() ) );
	if ( pose_sequence.compare(0,compare_len, tab_sequence,0,compare_len ) != 0 ) {
		bool match(true);
		for ( std::string::size_type i=0; i<=compare_len; ++i ) {
			if ( pose_sequence[i] != tab_sequence[i] && tab_sequence[i] != '_' ) {
				match=false;
				break;
			}
		}
		if ( match==false ) {
			tr.Debug << "The sequence of pose is:  " << pose_sequence << std::endl;
			tr.Debug << "The sequence of tab is:   " << tab_sequence << std::endl;
			utility_exit_with_message("the pose and tab don't start with the same residues");
		}
	}
}


Real Sparta::score_pose( core::pose::Pose const & pose ) {

	check_pose( pose);
	lib().setup_for_scoring( pose );
	return run_A_ANN_Prediction();
}

utility::vector1< float > Sparta::score_pose_per_residue(
	core::pose::Pose const & pose
) {
	lib().setup_for_scoring(pose);

	GDB PRED_SUM = lib().get_ANN_data( bCreateOutput_ );
	lib().deallocate_arrays();

	GDB COMP_TAB;
	utility::vector1< float > scores( pose.total_residue(), 0.0 );
	calc_per_residue_scores( lib().aN, PRED_SUM, REF_CS_Tab, COMP_TAB, scores );
	return scores;
}

//preset the args form command line SHIFT_DIR
void Sparta::SpartaLib::setup_defaults() {

	TAB_DIR = SPARTA_DIR + slash_char+ "tab";
	SHIFT_DIR = SPARTA_DIR + slash_char+ "shifts";
	PDB_DIR = SPARTA_DIR + slash_char+ "pdb";

	//later use Evaluator to determine scratch dir as in ExternalEvaluator...
	PRED_DIR = "pred";
	inName = "INTERNAL";
	// if( args["in"].length() > 0 ) inName = args["in"];
	//  if( args["ins"].length() > 0 ) inNames = args["ins"];

	tripFileName = TAB_DIR + slash_char+ "sparta.tab";
	weightFileName = TAB_DIR + slash_char + "weight.tab";
	homoFileName = TAB_DIR + slash_char + "homology.tab";
	fitFileName = TAB_DIR + slash_char + "fitting.tab";
	sumName = PRED_DIR + slash_char +"pred.tab";

	// if( args["ref"].length() > 0 ) refCSFileName = args["ref"];

	rcFileName = TAB_DIR + slash_char + "randcoil.tab";
	adjFileName = TAB_DIR + slash_char + "rcadj.tab";
	prevFileName = TAB_DIR + slash_char + "rcprev.tab";
	nextFileName = TAB_DIR + slash_char + "rcnext.tab";

	//Other Options
	EXCLUDED="";

	// if(args["atom"].length() > 0 ) {
	//   aN.clear();
	//   utility::vector0< string > temp = GDB::split(" ", args["atom"]);
	//   int cnt = 1;
	//   for(int i = 0; i < temp.size(); i++)
	//    {
	//     if( temp[i]!="N" && temp[i]!="HA"&& temp[i]!="C"&& temp[i]!="CA"&& temp[i]!="CB"&& temp[i]!="HN")
	//      {
	//       cerr << "\tInvalid atom -" << temp[i] << endl;
	//       exit(0);
	//      }
	//     aN[cnt++] = temp[i];
	//    }
	//  }

	matchCount = 20;
	tVal = 500.0; // not used
	firstRes = -9999;
	lastRes = 9999;

}

void Sparta::SpartaLib::init() {
	setup_defaults();
	tr.Info << "Reading Random Coil Shifts from " << rcFileName << endl;
	RC_Tab.loadGDB( rcFileName );

	tr.Info << "Reading RC Adjustments from " << adjFileName << endl;
	ADJ_Tab.loadGDB( adjFileName );

	//load BLOSUM62 table
	AAlist = "ACDEFGHIKLMNPQRSTVWY";
	GDB B62;
	string B62_fname = TAB_DIR + slash_char+ "BLOSUM62.tab";
	tr.Info << "Reading BLOSUM62 Table from " << B62_fname << endl;
	B62.loadGDB( B62_fname );
	for ( GDB::EntryList::iterator it = B62.Entries.begin(), end = B62.Entries.end(); it != end; ++it ) {
		string aa = (it->second)["RESNAME"];
		BLOSUM_62[aa]=BlosumMatrix::mapped_type( AAlist.size(), 0 );
	}
	for ( GDB::EntryList::iterator it = B62.Entries.begin(), end = B62.Entries.end(); it != end; ++it ) {
		//int index=it->first;
		string aa = (it->second)["RESNAME"];
		for ( GDB::GDB_Entry::iterator itS = (it->second).begin(), endS = (it->second).end(); itS != endS; ++itS ) {
			if ( itS->first == "RESNAME" ) continue;
			size_t index( AAlist.find( itS->first ) );
			runtime_assert( index  != string::npos);
			BLOSUM_62[aa][ index ] = atof( ( itS->second ).c_str() )/10.0;
		}
	} // end of assigning sequence homology vector (using blosum62 matrix)

	tr.Info << "Load ANN parameters ... ";
	for ( AtomNameList::iterator itN = aN.begin(), end = aN.end(); itN != end; ++itN ) {
		string atomName = itN->second;
		if ( atomName == "H" ) atomName="HN";
		SPARTA_ANN[atomName].init(113,30,1,9,6,3,TAB_DIR,atomName);
	}
	init_PredErrorSurface();
	tr.Info << "done " << std::endl;
}

//Get the list of angles\ring shifts\h-bond information from coordinates for all possible residues
//**************** NOT ABLE TO HANDLE PROTEIN WITH MULTIPLE CHAINS ****************
void Sparta::SpartaLib::getResInfo( bool create_output )
{
	const Real SPARTA_PI = numeric::NumericTraits<Real>::pi();
	const Real SPARTA_RADS_PER_DEG = SPARTA_PI / 180.0;
	const Real SIN_PI = sin(SPARTA_PI);
	const Real COS_PI = cos(SPARTA_PI);

	inTab.Entries.clear();

	// allocation
	int n = rN-r1+1, m = 10;
	U_ANGLES = new float* [n];
	U_ANGLES[0] = new float [n*m];
	for ( int i = 1; i < n; ++i ) {
		U_ANGLES[i] = U_ANGLES[i-1] + m;
	}

	U_RING_SHIFTS = new float* [n];
	U_RING_SHIFTS[0] = new float [n*(aN.size()+1)];
	for ( int i = 1; i < n; ++i ) {
		U_RING_SHIFTS[i] = U_RING_SHIFTS[i-1] + aN.size()+1;
	}

	n = rN-r1+1; m = 4;
	U_NAME = new string* [n];
	U_NAME[0] = new string [n*m];
	for ( int i = 1; i < n; ++i ) {
		U_NAME[i] = U_NAME[i-1] + m;
	}

	U_HN_HB = new float [n];
	U_HA_HB = new float [n];
	U_CO_HB = new float [n];

	int pos0 = inName.find_last_of(slash_char)+1;
	int pos1 = inName.find_last_of(".");

	sourceName=inName.substr(pos0,pos1-pos0);

	int cnt = 0;
	// format the sequence read from PDB coordinates
	sequence="";
	for ( ResidList::iterator itN = residList.begin(), end = residList.end(); itN != end; ++itN ) {
		sequence += itN->second;
		cnt++;
		if ( cnt%10 == 0 ) sequence += " "; //separator for each 10 residues
		++itN;
		if ( itN != residList.end() ) { //add "?" if sequence numbers are not consecutive
			int j = itN->first;
			--itN;
			for ( int i = 1; i< j - itN->first; i++ ) {
				sequence += "?"; cnt++;
				if ( cnt%10 == 0 ) sequence += " ";
			}
		} else --itN;
	}

	//   clock_t start, finish;
	//   start = clock();

	inPDB.initOrbitalShift();
	//finish = clock();
	//tr.Info << "\n\t initOrbitalShift running time: " << (float)(finish - start)/ CLOCKS_PER_SEC << " seconds" << endl;
	inPDB.initHBond();
	//finish = clock();
	//tr.Info << "\n\t initHBond running time: " << (float)(finish - start)/ CLOCKS_PER_SEC << " seconds" << endl;
	inPDB.collect_HN_S2_and_EF();
	//inPDB.calc_HN_S2();
	//finish = clock();
	//tr.Info << "\n\t calc_HN_S2 running time: " << (float)(finish - start)/ CLOCKS_PER_SEC << " seconds" << endl;

	inTab.setData("SEQUENCE", sequence);
	inTab.VARS_str_parser("  RESID_R1 RESNAME_R1 PHI_R1 PSI_R1 CHI1_R1 RESID_R2 RESNAME_R2 PHI_R2 PSI_R2 CHI1_R2 RESID_R3 RESNAME_R3 PHI_R3 PSI_R3 CHI1_R3 N_HM HA_HM C_HM CA_HM CB_HM H_HM H_HB HA_HB CO_HB SOURCE");
	inTab.FORMAT_str_parser("%4d %s %8.3f %8.3f %8.3f %4d %s %8.3f %8.3f %8.3f %4d %s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %s");

	//calculate H-bond for the first residue
	float dist = inPDB.HBDistList[r1]["HN"];
	U_HN_HB[0] = dist;
	dist = inPDB.HBDistList[r1]["HA"];
	U_HA_HB[0] = dist;
	dist = inPDB.HBDistList[r1]["O"];
	U_CO_HB[0] = dist;


	for ( int i = r1; i <= rN; i++ ) {
		CHI2_ANGLES[i] = inPDB.getChi2(1,i); // chi2 angle for residue r1, index by i-r1 (confusing..., but consistent with the loop)
		//OMEGA_ANGLES[i] = inPDB.getOmega(1,i); // chi2 angle for residue r1, index by i-r1 (confusing..., but consistent with the loop)
	}

	float shift;
	//loop for the polypeptide chain
	for ( int i = r1+1; i < rN; i++ ) {
		int index = i-r1;

		if ( residList.find(i-1) == residList.end() ||
				residList.find(i) == residList.end() ||
				residList.find(i+1) == residList.end() ) continue;


		shift = inPDB.getPhi(1,i+1);
		U_ANGLES[index][7] = shift;
		shift = inPDB.getPsi(1,i+1);
		U_ANGLES[index][8] = shift;
		shift = inPDB.getChi1(1,i+1);
		U_ANGLES[index][9] = shift;

		U_NAME[index][3] = residList[i+1];

		inTab.Entries[index]["PHI_R3"] = ftoa(U_ANGLES[index][7], buf);
		inTab.Entries[index]["PSI_R3"] = ftoa(U_ANGLES[index][8], buf);
		inTab.Entries[index]["CHI1_R3"] = ftoa(U_ANGLES[index][9], buf);
		inTab.Entries[index]["RESID_R3"] = itoa(i+1, buf);
		inTab.Entries[index]["RESNAME_R3"] = residList[i+1];
		inTab.Entries[index]["SOURCE"] = inName.substr(pos0,pos1-pos0);

		//Ring current shifts
		for ( AtomNameList::iterator itN_unordered = aN.begin(), end = aN.end(); itN_unordered != end; ++itN_unordered ) {
			string name = itN_unordered->second;
			if ( name == "H" ) {
				name = "HN";
				if ( residList[i] == "P" ) continue;
			} else if ( name == "HA" && residList[i] == "G" ) {
				U_RING_SHIFTS[index][7-1] = inPDB.getOrbitalShift(1,i,"HA3"); // change to use standard HA2/3 names
				name = "HA2";
			} else if ( name == "CB" && residList[i] == "G" ) continue;

			U_RING_SHIFTS[index][itN_unordered->first-1] = inPDB.getOrbitalShift(1,i,name) ;
			inTab.Entries[index][itN_unordered->second+"_HM"] = ftoa(U_RING_SHIFTS[index][itN_unordered->first-1], buf);
		}

		//H-Honds
		dist = inPDB.HBDistList[i]["HN"];
		inTab.Entries[index]["H_HB"] = ftoa(dist, buf);
		U_HN_HB[index] = dist;

		dist = inPDB.HBDistList[i]["HA"];
		inTab.Entries[index]["HA_HB"] = ftoa(dist, buf);
		U_HA_HB[index] = dist;

		dist = inPDB.HBDistList[i]["O"];
		inTab.Entries[index]["CO_HB"] = ftoa(dist, buf);
		U_CO_HB[index] = dist;

		if ( inTab.Entries.find(index-1) != inTab.Entries.end() ) {
			//if tripet i-1 exist
			//assign the values of positions 1 and 2 of triplet i using the postions 2 and 3 of triplet i-1
			U_ANGLES[index][1] = U_ANGLES[index-1][4];
			U_ANGLES[index][2] = U_ANGLES[index-1][5];
			U_ANGLES[index][3] = U_ANGLES[index-1][6];
			U_ANGLES[index][4] = U_ANGLES[index-1][7];
			U_ANGLES[index][5] = U_ANGLES[index-1][8];
			U_ANGLES[index][6] = U_ANGLES[index-1][9];
			U_NAME[index][1] = U_NAME[index-1][2];
			U_NAME[index][2] = U_NAME[index-1][3];

			inTab.Entries[index]["PHI_R1"] = inTab.Entries[index-1]["PHI_R2"];
			inTab.Entries[index]["PSI_R1"] = inTab.Entries[index-1]["PSI_R2"];
			inTab.Entries[index]["CHI1_R1"] = inTab.Entries[index-1]["CHI1_R2"];
			inTab.Entries[index]["RESID_R1"] = inTab.Entries[index-1]["RESID_R2"];
			inTab.Entries[index]["RESNAME_R1"] = inTab.Entries[index-1]["RESNAME_R2"];

			inTab.Entries[index]["PHI_R2"] = inTab.Entries[index-1]["PHI_R3"];
			inTab.Entries[index]["PSI_R2"] = inTab.Entries[index-1]["PSI_R3"];
			inTab.Entries[index]["CHI1_R2"] = inTab.Entries[index-1]["CHI1_R3"];
			inTab.Entries[index]["RESID_R2"] = inTab.Entries[index-1]["RESID_R3"];
			inTab.Entries[index]["RESNAME_R2"] = inTab.Entries[index-1]["RESNAME_R3"];
		} else { //else, calculate the values from coordinates

			shift = inPDB.getPhi(1,i-1);
			U_ANGLES[index][1] = shift;
			shift = inPDB.getPsi(1,i-1);
			U_ANGLES[index][2] = shift;
			shift = inPDB.getChi1(1,i-1);
			U_ANGLES[index][3] = shift;
			U_NAME[index][1] = residList[i-1];

			inTab.setEntry(index, "PHI_R1", ftoa(U_ANGLES[index][1], buf) );
			inTab.setEntry(index, "PSI_R1", ftoa(U_ANGLES[index][2], buf) );
			inTab.setEntry(index, "CHI1_R1", ftoa(U_ANGLES[index][3], buf) );
			inTab.setEntry(index, "RESID_R1", itoa(i-1, buf) );
			inTab.setEntry(index, "RESNAME_R1", residList[i-1] );

			shift = inPDB.getPhi(1,i);
			U_ANGLES[index][4] = shift;
			shift = inPDB.getPsi(1,i);
			U_ANGLES[index][5] = shift;
			shift = inPDB.getChi1(1,i);
			U_ANGLES[index][6] = shift;
			U_NAME[index][2] = residList[i];

			inTab.Entries[index]["PHI_R2"] = ftoa(U_ANGLES[index][4], buf);
			inTab.Entries[index]["PSI_R2"] = ftoa(U_ANGLES[index][5], buf);
			inTab.Entries[index]["CHI1_R2"] = ftoa(U_ANGLES[index][6], buf);
			inTab.Entries[index]["RESID_R2"] = itoa(i, buf);
			inTab.Entries[index]["RESNAME_R2"] = residList[i];

			if ( tr.Trace.visible() ) {
				tr.Trace << std::endl;
			}
		}

		//ANN input preparation
		// (20 BLOSSUM + 2 PHI + 2 PSI + 2 CHI1 + 2 CHI2)*3 + (4 ASA)*3 -chi2_c_asa
		// (20 BLOSSUM + 2 PHI + 2 PSI + 2 CHI1 + 2 CHI2 + 2 Oemga)*3 + (4 H-bond)*5 [O(i-1),HN,HA,O,HN(i+1)]
		utility::vector0<float> temp;
		//add ANN input for residue i-1
		string resName=residList[i-1]; if ( resName=="c" ) resName="C";
		//1-20
		temp.insert(temp.end(), BLOSUM_62[resName].begin(), BLOSUM_62[resName].end());
		//21-22
		float phi = U_ANGLES[index][1], psi = U_ANGLES[index][2], chi1 = U_ANGLES[index][3], chi2 = CHI2_ANGLES[i-1];//, omega=OMEGA_ANGLES[i-1];
		if ( phi<999 ) { temp.push_back(sin(phi*SPARTA_RADS_PER_DEG)); temp.push_back(cos(phi*SPARTA_RADS_PER_DEG));}//phi
		else { temp.push_back(SIN_PI); temp.push_back(COS_PI); }
		//23-24
		if ( psi<999 ) { temp.push_back(sin(psi*SPARTA_RADS_PER_DEG)); temp.push_back(cos(psi*SPARTA_RADS_PER_DEG));}//psi
		else { temp.push_back(SIN_PI); temp.push_back(COS_PI); }
		//25-26
		if ( chi1<999 ) { temp.push_back(sin(chi1*SPARTA_RADS_PER_DEG)); temp.push_back(cos(chi1*SPARTA_RADS_PER_DEG));}//chi1
		else { temp.push_back(0); temp.push_back(0); }
		//27
		temp.push_back(chi1<999);
		//28-29
		if ( chi2<999 ) { temp.push_back(sin(chi2*SPARTA_RADS_PER_DEG)); temp.push_back(cos(chi2*SPARTA_RADS_PER_DEG));}//chi2
		else { temp.push_back(0); temp.push_back(0); }
		//30
		temp.push_back(chi2<999);

		//add ANN input for residue i
		//31-50
		resName=residList[i]; if ( resName=="c" ) resName="C";
		temp.insert(temp.end(), BLOSUM_62[resName].begin(), BLOSUM_62[resName].end());
		phi = U_ANGLES[index][4]; psi = U_ANGLES[index][5]; chi1 = U_ANGLES[index][6]; chi2 = CHI2_ANGLES[i];//, omega=OMEGA_ANGLES[i];
		//51-52
		if ( phi<999 ) { temp.push_back(sin(phi*SPARTA_RADS_PER_DEG)); temp.push_back(cos(phi*SPARTA_RADS_PER_DEG));}//phi
		else { temp.push_back(SIN_PI); temp.push_back(COS_PI); }
		//53-54
		if ( psi<999 ) { temp.push_back(sin(psi*SPARTA_RADS_PER_DEG)); temp.push_back(cos(psi*SPARTA_RADS_PER_DEG));}//psi
		else { temp.push_back(SIN_PI); temp.push_back(COS_PI); }
		//55-56
		if ( chi1<999 ) { temp.push_back(sin(chi1*SPARTA_RADS_PER_DEG)); temp.push_back(cos(chi1*SPARTA_RADS_PER_DEG));}//chi1
		else { temp.push_back(0); temp.push_back(0); }
		//57
		temp.push_back(chi1<999);
		//58-59
		if ( chi2<999 ) { temp.push_back(sin(chi2*SPARTA_RADS_PER_DEG)); temp.push_back(cos(chi2*SPARTA_RADS_PER_DEG));}//chi2
		else { temp.push_back(0); temp.push_back(0); }
		//60
		temp.push_back(chi2<999);

		//add ANN input for residue i+1
		//61-80
		resName=residList[i+1]; if ( resName=="c" ) resName="C";
		temp.insert(temp.end(), BLOSUM_62[resName].begin(), BLOSUM_62[resName].end());
		//81-82
		phi = U_ANGLES[index][7]; psi = U_ANGLES[index][8]; chi1 = U_ANGLES[index][9]; chi2 = CHI2_ANGLES[i+1];//, omega=OMEGA_ANGLES[i+1];
		if ( phi<999 ) { temp.push_back(sin(phi*SPARTA_RADS_PER_DEG)); temp.push_back(cos(phi*SPARTA_RADS_PER_DEG));}//phi
		else { temp.push_back(SIN_PI); temp.push_back(COS_PI); }
		//83-84
		if ( psi<999 ) { temp.push_back(sin(psi*SPARTA_RADS_PER_DEG)); temp.push_back(cos(psi*SPARTA_RADS_PER_DEG));}//psi
		else { temp.push_back(SIN_PI); temp.push_back(COS_PI); }
		//85-86
		if ( chi1<999 ) { temp.push_back(sin(chi1*SPARTA_RADS_PER_DEG)); temp.push_back(cos(chi1*SPARTA_RADS_PER_DEG));}//chi1
		else { temp.push_back(0); temp.push_back(0); }
		//87
		temp.push_back(chi1<999);
		//88-89
		if ( chi2<999 ) { temp.push_back(sin(chi2*SPARTA_RADS_PER_DEG)); temp.push_back(cos(chi2*SPARTA_RADS_PER_DEG));}//chi2
		else { temp.push_back(0); temp.push_back(0); }
		//90
		temp.push_back(chi2<999);

		//91-94
		float hb = inPDB.HBDistList[i-1]["O"];
		if ( hb>0 ) {
			temp.push_back(1.0); temp.push_back(hb); temp.push_back( cos(inPDB.HB_DHO_AngleList[i-1]["O"]*SPARTA_RADS_PER_DEG) ); temp.push_back( cos(inPDB.HB_HOA_AngleList[i-1]["O"]*SPARTA_RADS_PER_DEG) );
		} else { temp.push_back(0.0); temp.push_back(0.0); temp.push_back(0.0); temp.push_back(0.0); }
		//95-98
		hb = inPDB.HBDistList[i]["HN"];
		if ( hb>0 ) {
			temp.push_back(1.0); temp.push_back(hb); temp.push_back( cos(inPDB.HB_DHO_AngleList[i]["HN"]*SPARTA_RADS_PER_DEG) ); temp.push_back( cos(inPDB.HB_HOA_AngleList[i]["HN"]*SPARTA_RADS_PER_DEG) );
		} else { temp.push_back(0.0); temp.push_back(0.0); temp.push_back(0.0); temp.push_back(0.0); }
		//99-102
		hb = inPDB.HBDistList[i]["HA"];
		if ( hb>0 ) {
			temp.push_back(1.0); temp.push_back(hb); temp.push_back( cos(inPDB.HB_DHO_AngleList[i]["HA"]*SPARTA_RADS_PER_DEG) ); temp.push_back( cos(inPDB.HB_HOA_AngleList[i]["HA"]*SPARTA_RADS_PER_DEG) );
		} else { temp.push_back(0.0); temp.push_back(0.0); temp.push_back(0.0); temp.push_back(0.0); }
		//103-106
		hb = inPDB.HBDistList[i]["O"];
		if ( hb>0 ) {
			temp.push_back(1.0); temp.push_back(hb); temp.push_back( cos(inPDB.HB_DHO_AngleList[i]["O"]*SPARTA_RADS_PER_DEG) ); temp.push_back( cos(inPDB.HB_HOA_AngleList[i]["O"]*SPARTA_RADS_PER_DEG) );
		} else { temp.push_back(0.0); temp.push_back(0.0); temp.push_back(0.0); temp.push_back(0.0); }
		//107-110
		hb = inPDB.HBDistList[i+1]["HN"];
		if ( hb>0 ) {
			temp.push_back(1.0); temp.push_back(hb); temp.push_back( cos(inPDB.HB_DHO_AngleList[i+1]["HN"]*SPARTA_RADS_PER_DEG) ); temp.push_back( cos(inPDB.HB_HOA_AngleList[i+1]["HN"]*SPARTA_RADS_PER_DEG) );
		} else { temp.push_back(0.0); temp.push_back(0.0); temp.push_back(0.0); temp.push_back(0.0); }

		//111-113
		temp.push_back(inPDB.HN_S2[i-1]);
		temp.push_back(inPDB.HN_S2[i]);
		temp.push_back(inPDB.HN_S2[i+1]);

		//done store
		runtime_assert( temp.size() == 113);
		ANN_IN_MTX[i]=temp;
	}


	//calculate H-bond for the last residue
	dist = inPDB.HBDistList[rN]["HN"];
	U_HN_HB[rN-r1] = dist;
	dist = inPDB.HBDistList[rN]["HA"];
	U_HA_HB[rN-r1] = dist;
	dist = inPDB.HBDistList[rN]["O"];
	U_CO_HB[rN-r1] = dist;

	if ( create_output ) inTab.saveGDB(PRED_DIR+slash_char+inName.substr(pos0,pos1-pos0) + "_in.tab");

	if ( tr.Trace.visible() ) {
		ANN::ANN_Matrix::iterator itX, end;
		for ( itX = ANN_IN_MTX.begin(), end = ANN_IN_MTX.end(); itX != end; ++itX ) {
			for ( int i=0; i< (int)(itX->second).size(); i++ ) {
				tr.Trace << (itX->first) << " " << (itX->second)[i] << std::endl;
			}
		}
	}
}


// run ANN prediction for a single protein chain
//void Sparta::runANN_Prediction() {
//  clock_t start/*, finish*/;
//  start = clock();

//  //  init(); now in constructor
//  if ( bCreateOutput_ ) {
//   // mkdir for prediction
//   if (PRED_DIR.find_last_of(slash_char) == PRED_DIR.length()-1 ) {
//    PRED_DIR = PRED_DIR.substr(0,PRED_DIR.length()-1);
//   }
//   mkdir_pred(PRED_DIR);
//  }
//   //for( itN = aN.begin(); itN != aN.end(); itN++ )
//   // mkdir_pred(PRED_DIR+slash_char+itN->second);

//  tr.Info << "Reading PDB Coordinates from " << inName << endl;
//  inPDB.loadPDB(inName);

//  residList = inPDB.residListOne;

//  r1 = inPDB.r1;
//  rN = inPDB.rN;

//  if (firstRes < r1) firstRes = r1;
//  if (lastRes  < r1) lastRes = r1;

//  if (firstRes > rN) firstRes = rN;
//  if (lastRes  > rN) lastRes = rN;

//  if (firstRes > lastRes) {
//   int itemp = firstRes;
//   firstRes  = lastRes;
//   lastRes   = itemp;
//  }

//   run_A_ANN_Prediction();
// }

// run ANN prediction for a single protein chain
Real Sparta::run_A_ANN_Prediction() {

	GDB COMP_Tab;
	GDB PRED_SUM = lib().get_ANN_data( bCreateOutput_ );
	Real score( compareRef_fxn( lib().aN, PRED_SUM, REF_CS_Tab, COMP_Tab ) );

	if ( bCreateOutput_ ) {
		COMP_Tab.addRemark( "Observed chemical shift from: " + refCSFileName );
		COMP_Tab.saveGDB( lib().sumName );
		REF_CS_Tab.saveGDB( lib().PRED_DIR + lib().slash_char + "ref.tab" );
	}

	lib().deallocate_arrays();
	return score;
}

void Sparta::SpartaLib::deallocate_arrays() {
	// deallocation - structures created in getResInfo
	delete [] U_ANGLES[0];
	delete [] U_ANGLES;
	delete [] U_NAME[0];
	delete [] U_NAME;
	delete [] U_RING_SHIFTS[0];
	delete [] U_RING_SHIFTS;
	delete [] U_HN_HB;
	delete [] U_HA_HB;
	delete [] U_CO_HB;
}

GDB Sparta::SpartaLib::get_ANN_data( bool create_output ) {
	/*clock_t start, finish*/;
	//start = clock();

	GDB PRED_SUM;

	tr.Info << "Analyzing " << inName << " " ;
	getResInfo( create_output ); // loads info from PDB
	tr.Info << residList.size() << " residues read " << endl;

	//finish = clock();
	//tr.Info << "\t getResInfo() running time: " << (float)(finish - start)/ CLOCKS_PER_SEC << " seconds" << endl;

	//start = clock();
	tr.Info << "ANN prediction ..." << endl;
	for ( AtomNameList::iterator itN = aN.begin(), end = aN.end(); itN != end; ++itN ) {
		string atomName = itN->second;
		if ( atomName == "H" ) atomName="HN";

		SPARTA_ANN[atomName].ANN_OUT_MTX_LEVEL1.clear();
		SPARTA_ANN[atomName].runSpartaANN(ANN_IN_MTX);

		ANN_CS_OUTPUT_FULL[atomName] = SPARTA_ANN[atomName].ANN_OUT_MTX_LEVEL1;
	}

	//finish = clock();
	//tr.Info << "\t ANNPredict() running time: " << (float)(finish - start)/ CLOCKS_PER_SEC << " seconds" << endl;

	//GDB PRED_SUM;
	PRED_SUM.VARS_str_parser("  RESID RESNAME ATOMNAME SS_SHIFT SHIFT RC_SHIFT HM_SHIFT EF_SHIFT SIGMA SOURCE");
	PRED_SUM.FORMAT_str_parser(" %4d %4s %4s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %s");
	string str = itoa(r1, buf);
	PRED_SUM.setData("FIRST_RESID", str+"\n");
	PRED_SUM.setData("SEQUENCE", sequence);

	float RC, RCadj, pred_2nd_shift, pred_shift/*, HB*/;
	for ( int i = r1+1; i <= rN-1; i++ ) { //olange: we have not loaded the ANN with stuff for residue 1 or rN as it would be the 0,1,2 triplett.. ignore here TOO!
		for ( AtomNameList::iterator itN = aN.begin(), end = aN.end(); itN != end; ++itN ) {
			string atomName = itN->second;
			if ( atomName == "H" ) atomName="HN";
			int index = PRED_SUM.Entries.size()+1;

			if ( residList[i].empty() ) continue;
			if ( residList[i] == "P" && (atomName == "HN" || atomName == "N") ) continue;
			if ( residList[i] == "G" && atomName == "CB" ) continue;
			if ( i==r1 && (atomName == "HN"|| atomName == "N") ) continue; //added from email Yang Shen/ Aug 6th.
			if ( i==rN && atomName == "C" ) continue; //added from email Yang Shen/ Aug 6th

			PRED_SUM.setEntry(index, "RESID",  itoa(i,buf));
			PRED_SUM.setEntry(index, "RESNAME",  residList[i]);
			if ( atomName=="HA" && residList[i] == "G" ) PRED_SUM.setEntry(index, "ATOMNAME",  "HA2"); //added from email YangShen Aug 6th.
			else PRED_SUM.setEntry(index, "ATOMNAME",  atomName);


			RC = getRC(residList[i],atomName);
			RCadj = getRCadj(residList[i],atomName);
			if ( i==r1 || i==rN ) {
				pred_2nd_shift = 0.0; //may not good for the last residue, for which the neighoring residue effect is not considered.
			} else {
				pred_2nd_shift = 0.0;
				if ( ANN_CS_OUTPUT_FULL[atomName].size() >= static_cast< Size > (i) ) {
					pred_2nd_shift = ANN_CS_OUTPUT_FULL[atomName][i][0];
				}

				if      ( atomName == "HA" ) pred_2nd_shift /= 4.0;
				else if ( atomName == "HN" ) pred_2nd_shift /= 2.0;
				else if ( atomName == "N" ) pred_2nd_shift  *= 2.5;

				if ( pred_2nd_shift > 20.0 || pred_2nd_shift < -20.0 ) pred_shift = 0.0;

				/*
				if(atomName == "HA") pred_2nd_shift = ANN_CS_OUTPUT_FULL[atomName][i][0]/4.0;
				else if(atomName == "HN") pred_2nd_shift = ANN_CS_OUTPUT_FULL[atomName][i][0]/2.0;
				else if(atomName == "N") pred_2nd_shift = ANN_CS_OUTPUT_FULL[atomName][i][0]*2.5;
				else pred_2nd_shift = ANN_CS_OUTPUT_FULL[atomName][i][0];
				if( pred_2nd_shift > 20.0 || pred_2nd_shift < -20.0 ) pred_shift = 0.0;
				*/
			}

			pred_shift = pred_2nd_shift + RC + RCadj; // + PrevRCadj + NextRCadj;
			if ( pred_shift > 999.0 ) pred_shift = SPARTA_MAX_NUM;


			PRED_SUM.setEntry(index, "SS_SHIFT",  ftoa(pred_2nd_shift,buf));

			pred_shift += 0.6*atof(inTab.Entries[i-r1][atomName+"_HM"].c_str());
			if ( atomName == "HN" || atomName == "HA" ) pred_shift-= inPDB.ElectricField[i][atomName]; //marked off to exclude shifts from "global" contacts and to test MFR
			PRED_SUM.setEntry(index, "SHIFT",  ftoa(pred_shift,buf));

			PRED_SUM.setEntry(index, "RC_SHIFT", ftoa(RC+RCadj ,buf) );
			PRED_SUM.setEntry(index, "SOURCE", sourceName );
			PRED_SUM.setEntry(index, "SIGMA", ftoa(getANN_PredError(U_ANGLES[i-r1][4],U_ANGLES[i-r1][5],residList[i],atomName),buf) );
			//tr.Info << U_ANGLES[i-r1][4] << "\t" << U_ANGLES[i-r1][5] << "\t" << getANN_PredError(U_ANGLES[i-r1][4],U_ANGLES[i-r1][5],residList[i],atomName) << endl;
			PRED_SUM.setEntry(index, "HM_SHIFT", inTab.Entries[i-r1][atomName+"_HM"] );
			PRED_SUM.setEntry(index, "EF_SHIFT", ftoa(inPDB.ElectricField[i][atomName],buf) );

			if ( atomName=="HA" && residList[i] == "G" ) { // for GLY HA3
				index++;
				PRED_SUM.setEntry(index, "RESID",  itoa(i,buf));
				PRED_SUM.setEntry(index, "RESNAME",  residList[i]);
				PRED_SUM.setEntry(index, "ATOMNAME",  "HA3");
				PRED_SUM.setEntry(index, "SS_SHIFT",  ftoa(pred_2nd_shift,buf));
				pred_shift = pred_2nd_shift + RC + RCadj + 0.6*U_RING_SHIFTS[i-r1][0]; // atof(inTab.Entries[i-r1][atomName+"_HM"].c_str());;
				pred_shift-= inPDB.ElectricField[i]["HA"];
				PRED_SUM.setEntry(index, "SHIFT",  ftoa(pred_shift,buf));
				PRED_SUM.setEntry(index, "RC_SHIFT", ftoa(RC+RCadj ,buf) );
				PRED_SUM.setEntry(index, "SOURCE", sourceName );
				PRED_SUM.setEntry(index, "SIGMA", ftoa(getANN_PredError(U_ANGLES[i-r1][4],U_ANGLES[i-r1][5],residList[i],atomName),buf) );
				PRED_SUM.setEntry(index, "HM_SHIFT", ftoa(U_RING_SHIFTS[i-r1][0],buf) );
				PRED_SUM.setEntry(index, "EF_SHIFT", ftoa(inPDB.ElectricField[i]["HA"],buf) );
			}
		}
	}

	if ( tr.Debug.visible() ) {
		tr.Debug << " ============== PRED_SUM ==================== " << std::endl;
		PRED_SUM.showGDB( tr.Debug );
		tr.Debug << " ============== END_ PRED_SUM ==================== " << std::endl;
	}
	//finish = clock();
	//tr.Info << "\t ANNPredict() running time: " << (float)(finish - start)/ CLOCKS_PER_SEC << " seconds" << endl;
	if ( create_output ) {
		PRED_SUM.saveGDB(sumName);
		// snippet moved from compareRef. Not sure why this is done twice
		// in slightly different ways ...
		int pos = sumName.find_last_of(".");
		PRED_SUM.saveGDB( sumName.substr(0, pos) + "_full.tab" ); //save the original prediction summary file to a new name
	}

	return PRED_SUM;
} // get_ANN_data

// // run ANN prediction for multiple protein chains
// void Sparta::runANN_Predictions() {
//   //init(); //now in constructor
//   // mkdir for prediction
//  if ( bCreateOutput_ ) {

//   if (PRED_DIR.find_last_of(slash_char) == PRED_DIR.length()-1 ) {
//    PRED_DIR = PRED_DIR.substr(0,PRED_DIR.length()-1);
//   }
//   mkdir_pred(PRED_DIR);
//   // for( itN = aN.begin(); itN != aN.end(); itN++ )
//   //  mkdir_pred(PRED_DIR+slash_char+itN->second);
//  }

//   utility::vector0< string > temp = split(" ", inNames);
//   string  outName = sumName;

//   tr.Info << inNames << endl;

//  for ( Size i = 0; i < temp.size(); i++) {
//   inName = temp[i];
//   //tr.Info << "Reading PDB Coordinates from " << inName << endl;
//   inPDB.loadPDB(inName);

//   residList = inPDB.residListOne;

//   r1 = inPDB.r1;
//   rN = inPDB.rN;

//   if (firstRes < r1) firstRes = r1;
//   if (lastRes  < r1) lastRes = r1;

//   if (firstRes > rN) firstRes = rN;
//   if (lastRes  > rN) lastRes = rN;

//   if (firstRes > lastRes) {
//    int itemp = firstRes;
//    firstRes  = lastRes;
//    lastRes   = itemp;
//   }

//   ANN_IN_MTX.clear();
//   ANN_CS_OUTPUT_FULL.clear();

//   int pos0 = inName.find_last_of(slash_char)+1;
//   int pos1 = inName.find_last_of(".");

//   sourceName=inName.substr(pos0,pos1-pos0);
//   sumName = PRED_DIR + slash_char + sourceName + "_pred.tab";

//   run_A_ANN_Prediction();
//   tr.Info << "\tPrediction file " << sumName << " is ready for protein " << inName << endl;

//  }

// }


// Initiate an ANN prediction for a single protein using its file name
//void Sparta::runANN_Prediction(const string& pName)
//{
//init(); now in constructor
//  if ( bCreateOutput_ ) {

//   // mkdir for prediction
//   if (PRED_DIR.find_last_of(slash_char) == PRED_DIR.length()-1 )
//    PRED_DIR = PRED_DIR.substr(0,PRED_DIR.length()-1);
//   for( itN = aN.begin(); itN != aN.end(); itN++ )
//    mkdir_pred(PRED_DIR+slash_char+itN->second);
//  }
//  inName = pName;
//  inPDB.loadPDB(inName);

// run_A_ANN_Prediction();
//}


void Sparta::SpartaLib::init_PredErrorSurface()
{
	int step = 5;

	for ( AtomNameList::iterator itN = aN.begin(), end = aN.end(); itN != end; ++itN ) {
		string atomName = itN->second;
		if ( atomName == "H" ) atomName="HN";

		for ( Size i=0; i<AAlist.length(); i++ ) {
			string AA = AAlist.substr(i,1);
			if ( AA == " " ) continue;

			if ( AA == "G" && atomName == "CB" ) continue;
			if ( AA == "P" && atomName == "HN" ) continue;

			string surfName = TAB_DIR + slash_char + "errorSurface" + slash_char + atomName + slash_char + AA + "..A450.S5.RMS.tab";
			GDB surf(surfName);

			for ( GDB::EntryList::iterator it = surf.Entries.begin(); it != surf.Entries.end(); ++it ) {
				int phi = atoi( it->second["PHI"].c_str() );
				for ( int y=-180; y<180; y+=step ) {
					string psi=itoa(y,buf);
					SPARTA_ERR_SURF[AA][atomName][phi][y] = atof( it->second[psi].c_str() );
				}
			}
		}
	}
}


float Sparta::SpartaLib::getANN_PredError(float phi, float psi, string aa, string aName)
{
	return SPARTA_ERR_SURF[aa][aName][5*int(phi/5)][5*int(psi/5)];
}


// get random coil chemical shift for atom 'aName' of residue 'resName'
float Sparta::SpartaLib::getRC(const string& resName, const string& aName)
{
	GDB::GDB_Entry temp = RC_Tab.getEntry("RESNAME",resName,1);

	if ( temp.size() != 0 ) return atof( temp[aName].c_str() );

	return 9999.0;
}


float Sparta::SpartaLib::getRCadj(const string& resName, const string& aName)
{
	GDB::GDB_Entry temp = ADJ_Tab.getEntry("RESNAME",resName,1);

	if ( temp.size() != 0 ) return atof( temp[aName].c_str() );

	return 0.0;
}


float Sparta::SpartaLib::getPrevRCadj(const string& prev_rName, const string& aName)
{
	GDB::GDB_Entry temp = PREV_Tab.getEntry("RESNAME",prev_rName,1);

	if ( temp.size() != 0 ) return atof( temp[aName].c_str() );

	return 0.0;
}


float Sparta::SpartaLib::getNextRCadj(const string& next_rName, const string& aName)
{
	GDB::GDB_Entry temp = NEXT_Tab.getEntry("RESNAME",next_rName,1);

	if ( temp.size() != 0 ) return atof( temp[aName].c_str() );

	return 0.0;
}

float Sparta::SpartaLib::getWeight(const string& Name, const string& aName)
{
	GDB::GDB_Entry temp = WEIGHT_Tab.getEntry("RESNAME",Name,1);

	if ( temp.size() != 0 ) return atof( temp[aName].c_str() );

	return 9999.0;

}

void Sparta::SpartaLib::mkdir_pred(const string& d)// create a directory for prediction results
{
	if ( ! isDirExists( d ) ) {
		int i = MKDIR(d.c_str())+1;
		if ( !i ) { // if not success
			string parentD = d.substr(0,d.find_last_of(slash_char)+1);
			mkdir_pred(parentD.c_str());
			if ( !(MKDIR(d.c_str())+1) ) {
				string msg("\tCan't create prediction directory " + d);
				utility_exit_with_message(msg);
			}
		}
	}
}

} // namespace sparta
} // namespace protocols
