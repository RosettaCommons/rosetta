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
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author James Thompson

// Unit Headers
#include <boost/unordered_map.hpp>
#include <protocols/sparta/SpartaUtil.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>

#include <cstdio>
#include <string>


#ifdef WIN32
#include <direct.h>
#else
#include <dirent.h>
#include <sys/stat.h>

//Auto Headers
#endif
static THREAD_LOCAL basic::Tracer tr( "protocols.sparta" );

namespace protocols {
namespace sparta {

using namespace core;
using namespace std;

void
calc_per_residue_scores(
	Sparta::SpartaLib::AtomNameList& atom_names,
	GDB& Pred_Sum,
	GDB& REF_CS_Tab,
	GDB& COMP_Tab,
	utility::vector1< float >& per_residue_scores
) {
	using boost::unordered_map;

	COMP_Tab.REMARKS = Pred_Sum.REMARKS;
	COMP_Tab.DATA = Pred_Sum.DATA;

	// to compute chemshift score take CS_DIFF^2 / SIGMA ^2

	COMP_Tab.VARS_str_parser("  RESID RESNAME ATOMNAME CS_OBS SHIFT RC_SHIFT CS_DIFF SIGMA W");
	COMP_Tab.FORMAT_str_parser("  %4d %4s %4s %9.3f %9.3f %9.3f %9.3f %9.3f %.2f");

	//boost::unordered_map< int, string >::iterator itN;
	// boost::unordered_map< int, boost::unordered_map< string, string > >::iterator it;

	//std::cout << "begin_scoring: scores.size() = " << per_residue_scores.size() << std::endl;

	utility::vector0< float > OBS_V, PRED_V, DIFF_V, OBS_V_CORRECTED;
	for (auto & atom_name : atom_names) {
		string aName = atom_name.second;
		if ( aName == "H" ) aName="HN";
		bool floating_sign( REF_CS_Tab.isVarFloat("SHIFT2") );
		if ( floating_sign ) tr.Info << " use floating sign1 " << std::endl;
		else tr.Info << " no floating sign " << std::endl;
		for ( auto it = REF_CS_Tab.Entries.begin(), end2 = REF_CS_Tab.Entries.end(); it != end2; ++it ) {
			float obs_shift, pred_shift, obs_shift2( 0.0 );
			string aName_ref = it->second["ATOMNAME"];
			if ( aName_ref == "H" ) aName_ref = "HN";

			if ( aName == "HA" && aName_ref == "G" ) {
				// AMW: cppcheck claims that compare would be faster than find here
				if ( aName_ref.find("HA") != 0 ) continue;
			} else if ( aName_ref != aName ) continue;

			GDB::GDB_Entry temp = Pred_Sum.getEntry("RESID",it->second["RESID"],"ATOMNAME",aName,1);

			if ( temp["SHIFT"].length() <= 0 ) continue;

			obs_shift = atof( (it->second["SHIFT"]).c_str() );
			if ( floating_sign ) obs_shift2 = atof( (it->second["SHIFT2"]).c_str() );

			pred_shift = atof( (temp["SHIFT"]).c_str() );

			if ( aName == "HA" && it->second["RESNAME"] == "G" ) { //for HA of Gly
				if ( aName_ref == "HA2" ) { // assign HA2 to the one with smaller shift
					float shift_HA3 = 9999.000;
					float shift_HA3_2 = 9999;
					GDB::GDB_Entry HA3 = REF_CS_Tab.getEntry("RESID",it->second["RESID"],"ATOMNAME","HA3",1);
					if ( HA3["ATOMNAME"] == "HA3" ) shift_HA3 = atof( (HA3["SHIFT"]).c_str() );
					if ( floating_sign && HA3["ATOMNAME"] == "HA3" ) shift_HA3_2 = atof( (HA3["SHIFT2"]).c_str() );

					if ( shift_HA3 < obs_shift ) obs_shift = shift_HA3;
					if ( floating_sign && shift_HA3_2 < obs_shift2 ) obs_shift2 = shift_HA3_2;

					shift_HA3 = 9999.000;
					HA3 = Pred_Sum.getEntry("RESID",it->second["RESID"],"ATOMNAME","HA3",1);
					if ( HA3["ATOMNAME"] == "HA3" ) shift_HA3 = atof( (HA3["SHIFT"]).c_str() );

					if ( shift_HA3 < pred_shift ) pred_shift = shift_HA3;
				} else if ( aName_ref == "HA3" ) { // assign HA3 to the one with larger shift

					float shift_HA2 = -9999.000;
					float shift_HA2_2 = -9999;
					GDB::GDB_Entry HA2 = REF_CS_Tab.getEntry("RESID",it->second["RESID"],"ATOMNAME","HA2",1);
					if ( HA2["ATOMNAME"] == "HA2" ) shift_HA2 = atof( (HA2["SHIFT"]).c_str() );
					if ( floating_sign && HA2["ATOMNAME"] == "HA2" ) shift_HA2_2 = atof( (HA2["SHIFT2"]).c_str() );

					if ( shift_HA2 > obs_shift ) obs_shift = shift_HA2;
					if ( floating_sign && shift_HA2_2 > obs_shift2 ) obs_shift2 = shift_HA2_2;

					shift_HA2 = -9999.000;
					HA2 = Pred_Sum.getEntry("RESID",it->second["RESID"],"ATOMNAME","HA2",1);
					if ( HA2["ATOMNAME"] == "HA2" ) shift_HA2 = atof( (HA2["SHIFT"]).c_str() );

					if ( shift_HA2 > pred_shift ) pred_shift = shift_HA2;
				}
			}

			if ( pred_shift > 999.0 ) continue;

			if ( std::fabs( pred_shift-obs_shift ) > std::fabs( pred_shift-obs_shift2 ) ) {
				obs_shift=obs_shift2;
			}
			OBS_V.push_back( obs_shift );
			PRED_V.push_back( pred_shift );
			DIFF_V.push_back( pred_shift-obs_shift );

			Real sigma( atof( (temp["SIGMA"]).c_str() ) );
			if ( sigma > 0.1 ) {
				Size const res_id( ObjexxFCL::ulong_of( it->second["RESID"] ) );
				Real const score(
					(obs_shift-pred_shift)*(obs_shift-pred_shift)/(sigma*sigma)
				);
				//std::cout << "score(" << res_id << "," << aName << ") = " << score << std::endl;
				//if ( per_residue_scores.size() < res_id ) {
				// per_residue_scores.resize( res_id, 0.0 );
				//}

				// don't run off the end of the per_residue_scores, as the per_residue_scores tell
				// you how many residues actually exist in your Pose.
				if ( per_residue_scores.size() > res_id ) {
					//std::cout << per_residue_scores[res_id] << " += " << score << std::endl;
					per_residue_scores[res_id] += score;
				}
			}
		} // REF_CS_Tab

		if ( tr.Info.visible() ) {
			tr.Info << "RMS(OBS, PRED) for " << aName << " (n=" << OBS_V.size() << "): "
				<<  getRMS(OBS_V, PRED_V) <<  " ppm " << std::endl;
		}

		float avg_diff = getAVG(DIFF_V);
		for ( Size i = 0; i< OBS_V.size(); i++ ) {
			OBS_V_CORRECTED.push_back( OBS_V[i] + avg_diff );
		}

		OBS_V.clear(); PRED_V.clear(); DIFF_V.clear(); OBS_V_CORRECTED.clear();
	} // for ( itN = atom_names.begin(); itN != atom_names.end(); itN++ )

	if ( tr.Debug.visible() ) {
		tr.Debug << " ============== COMP_Tab ==================== " << std::endl;
		COMP_Tab.showGDB( tr.Debug );
		tr.Debug << " ============== END COMP_Tab ==================== " << std::endl;
		tr.Debug << " ============== REF_CS_TAB ==================== " << std::endl;
		REF_CS_Tab.showGDB( tr.Debug );
		tr.Debug << " ============== END_CS_TAB ==================== " << std::endl;
	}
}

Real compareRef_fxn(
	Sparta::SpartaLib::AtomNameList& names,
	GDB & Pred_Sum,
	GDB & REF_CS_Tab,
	GDB & COMP_Tab
) {
	using utility::vector1;
	vector1< float > scores( Pred_Sum.Entries.size(), 0.0 );
	//vector1< float > scores;
	calc_per_residue_scores( names, Pred_Sum, REF_CS_Tab, COMP_Tab, scores );

	typedef vector1< float >::const_iterator iter;
	Real SCORE_SUM(0.0);
	for ( iter it = scores.begin(), end = scores.end(); it != end; ++it ) {
		SCORE_SUM += *it;
	}
	return SCORE_SUM/4;
} // compareRef_fxn

float getDiff( float ang1, float ang2 ) {
	float a = fabs(ang1-ang2);

	if ( a < 180 ) return a;
	else if ( a > 180 && a < 360.1 && ang1*ang2 < 0 ) {
		return 360.0-a;
	}

	return 0;
}


float getAVG(utility::vector0<float> &v1) {
	int n = v1.size();

	float sum = 0.0;

	for ( int i = 0; i < n; i++ ) {
		sum += v1[i];
	}

	return sum/(float)n;
}

float getSTD(utility::vector0<float> &v1) {
	int n = v1.size();
	float avg = getAVG(v1);

	float sum = 0.0;

	for ( int i = 0; i < n; i++ ) {
		sum += (v1[i]-avg)*(v1[i]-avg);
	}

	return sqrt(sum/(float)(n-1));
}

float getRMS(utility::vector0<float> &v1, utility::vector0<float> &v2) {
	if ( v1.size() != v2.size() ) return -1.0;

	int n = v1.size();

	float sum = 0.0;
	// using iterators should speed this up ...
	for ( int i = 0; i < n; i++ ) {
		sum += ( (v1[i]-v2[i])*(v1[i]-v2[i]) );
	}

	return sqrt(sum/(float)(n-1));
}


char * itoa( int n, char *buff, int /*base*/ )
{
	sprintf(buff, "%d", n);
	return buff;
}


char * ftoa( float n, char *buff, char f, int prec )
{
	if ( !(f=='f' || f=='F' || f=='e' || f=='E' || f=='g' || f=='G') ) {
		f = 'f';
	}
	char format[20];
	char *fs = format;    // generate format string
	*fs++ = '%';     // "%.<prec>l<f>"
	if ( prec >= 0 ) {
		if ( prec > 99 ) {   // buf big enough for precision?
			prec = 99;
		}
		*fs++ = '.';
		if ( prec >= 10 ) {
			*fs++ = prec / 10 + '0';
			*fs++ = prec % 10 + '0';
		} else {
			*fs++ = prec + '0';
		}
	}
	*fs++ = 'l';
	*fs++ = f;
	*fs = '\0';
	sprintf( buff, format, n );

	return buff;
}

int MKDIR(const char *dirName)
{
#ifndef  __native_client__
#ifdef WIN32
	return mkdir(dirName);
#else
	// 777? dangerous ...
	return mkdir(dirName, 0777);
#endif
#endif
}

bool isDirExists(const string &Dir) {
#ifdef WIN32
	char oldDir[100];
	if ( NULL == getcwd(oldDir, _MAX_PATH) ) return false;

	if ( chdir( Dir.c_str() ) == -1 ) return false;

	if ( chdir( oldDir ) == -1 ) return false;

	return true;
#else
	DIR *dp;
	dp = opendir(Dir.c_str());
	if ( dp != nullptr ) {
		(void) closedir (dp);
		return true;
	} else {
		return false;
	}
#endif
}

} // namespace sparta
} // namespace protocols
