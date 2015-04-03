// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packstat/types.hh
///
/// @brief
/// @author will sheffler


#ifndef INCLUDED_core_scoring_packstat_sasa_dot_locations_hh
#define INCLUDED_core_scoring_packstat_sasa_dot_locations_hh

#include <core/scoring/packstat/types.hh>

namespace core {
namespace scoring {
namespace packstat {

inline utility::vector1<numeric::xyzVector<PackstatReal> > const get_sasa_dot_locations() {

	utility::vector1<numeric::xyzVector<PackstatReal> > sasa_dot_locations(162);

  sasa_dot_locations[1].x(0.000000);
  sasa_dot_locations[1].y(0.000000);
  sasa_dot_locations[1].z(1.000000);
  sasa_dot_locations[2].x(0.276393);
  sasa_dot_locations[2].y(0.850651);
  sasa_dot_locations[2].z(0.447214);
  sasa_dot_locations[3].x(0.894427);
  sasa_dot_locations[3].y(0.000000);
  sasa_dot_locations[3].z(0.447214);
  sasa_dot_locations[4].x(0.162460);
  sasa_dot_locations[4].y(0.500000);
  sasa_dot_locations[4].z(0.850651);
  sasa_dot_locations[5].x(0.525731);
  sasa_dot_locations[5].y(0.000000);
  sasa_dot_locations[5].z(0.850651);
  sasa_dot_locations[6].x(0.077609);
  sasa_dot_locations[6].y(0.238856);
  sasa_dot_locations[6].z(0.967949);
  sasa_dot_locations[7].x(0.251148);
  sasa_dot_locations[7].y(0.000000);
  sasa_dot_locations[7].z(0.967949);
  sasa_dot_locations[8].x(0.361803);
  sasa_dot_locations[8].y(0.262866);
  sasa_dot_locations[8].z(0.894427);
  sasa_dot_locations[9].x(0.688191);
  sasa_dot_locations[9].y(0.500000);
  sasa_dot_locations[9].z(0.525731);
  sasa_dot_locations[10].x(0.232827);
  sasa_dot_locations[10].y(0.716567);
  sasa_dot_locations[10].z(0.657513);
  sasa_dot_locations[11].x(0.447214);
  sasa_dot_locations[11].y(0.525731);
  sasa_dot_locations[11].z(0.723607);
  sasa_dot_locations[12].x(0.483974);
  sasa_dot_locations[12].y(0.716567);
  sasa_dot_locations[12].z(0.502295);
  sasa_dot_locations[13].x(0.638197);
  sasa_dot_locations[13].y(0.262866);
  sasa_dot_locations[13].z(0.723607);
  sasa_dot_locations[14].x(0.753443);
  sasa_dot_locations[14].y(0.000000);
  sasa_dot_locations[14].z(0.657513);
  sasa_dot_locations[15].x(0.831052);
  sasa_dot_locations[15].y(0.238856);
  sasa_dot_locations[15].z(0.502295);
  sasa_dot_locations[16].x(-0.723607);
  sasa_dot_locations[16].y(0.525731);
  sasa_dot_locations[16].z(0.447214);
  sasa_dot_locations[17].x(-0.425325);
  sasa_dot_locations[17].y(0.309017);
  sasa_dot_locations[17].z(0.850651);
  sasa_dot_locations[18].x(-0.203183);
  sasa_dot_locations[18].y(0.147621);
  sasa_dot_locations[18].z(0.967949);
  sasa_dot_locations[19].x(-0.138197);
  sasa_dot_locations[19].y(0.425325);
  sasa_dot_locations[19].z(0.894427);
  sasa_dot_locations[20].x(-0.262866);
  sasa_dot_locations[20].y(0.809017);
  sasa_dot_locations[20].z(0.525731);
  sasa_dot_locations[21].x(-0.609548);
  sasa_dot_locations[21].y(0.442863);
  sasa_dot_locations[21].z(0.657513);
  sasa_dot_locations[22].x(-0.361803);
  sasa_dot_locations[22].y(0.587785);
  sasa_dot_locations[22].z(0.723607);
  sasa_dot_locations[23].x(-0.531939);
  sasa_dot_locations[23].y(0.681718);
  sasa_dot_locations[23].z(0.502295);
  sasa_dot_locations[24].x(-0.052786);
  sasa_dot_locations[24].y(0.688191);
  sasa_dot_locations[24].z(0.723607);
  sasa_dot_locations[25].x(0.029644);
  sasa_dot_locations[25].y(0.864188);
  sasa_dot_locations[25].z(0.502295);
  sasa_dot_locations[26].x(-0.723607);
  sasa_dot_locations[26].y(-0.525731);
  sasa_dot_locations[26].z(0.447214);
  sasa_dot_locations[27].x(-0.425325);
  sasa_dot_locations[27].y(-0.309017);
  sasa_dot_locations[27].z(0.850651);
  sasa_dot_locations[28].x(-0.203183);
  sasa_dot_locations[28].y(-0.147621);
  sasa_dot_locations[28].z(0.967949);
  sasa_dot_locations[29].x(-0.447214);
  sasa_dot_locations[29].y(0.000000);
  sasa_dot_locations[29].z(0.894427);
  sasa_dot_locations[30].x(-0.850651);
  sasa_dot_locations[30].y(0.000000);
  sasa_dot_locations[30].z(0.525731);
  sasa_dot_locations[31].x(-0.609548);
  sasa_dot_locations[31].y(-0.442863);
  sasa_dot_locations[31].z(0.657513);
  sasa_dot_locations[32].x(-0.670820);
  sasa_dot_locations[32].y(-0.162460);
  sasa_dot_locations[32].z(0.723607);
  sasa_dot_locations[33].x(-0.812731);
  sasa_dot_locations[33].y(-0.295242);
  sasa_dot_locations[33].z(0.502295);
  sasa_dot_locations[34].x(-0.670820);
  sasa_dot_locations[34].y(0.162460);
  sasa_dot_locations[34].z(0.723607);
  sasa_dot_locations[35].x(-0.812731);
  sasa_dot_locations[35].y(0.295242);
  sasa_dot_locations[35].z(0.502295);
  sasa_dot_locations[36].x(0.276393);
  sasa_dot_locations[36].y(-0.850651);
  sasa_dot_locations[36].z(0.447214);
  sasa_dot_locations[37].x(0.162460);
  sasa_dot_locations[37].y(-0.500000);
  sasa_dot_locations[37].z(0.850651);
  sasa_dot_locations[38].x(0.077609);
  sasa_dot_locations[38].y(-0.238856);
  sasa_dot_locations[38].z(0.967949);
  sasa_dot_locations[39].x(-0.138197);
  sasa_dot_locations[39].y(-0.425325);
  sasa_dot_locations[39].z(0.894427);
  sasa_dot_locations[40].x(-0.262866);
  sasa_dot_locations[40].y(-0.809017);
  sasa_dot_locations[40].z(0.525731);
  sasa_dot_locations[41].x(0.232827);
  sasa_dot_locations[41].y(-0.716567);
  sasa_dot_locations[41].z(0.657513);
  sasa_dot_locations[42].x(-0.052786);
  sasa_dot_locations[42].y(-0.688191);
  sasa_dot_locations[42].z(0.723607);
  sasa_dot_locations[43].x(0.029644);
  sasa_dot_locations[43].y(-0.864188);
  sasa_dot_locations[43].z(0.502295);
  sasa_dot_locations[44].x(-0.361803);
  sasa_dot_locations[44].y(-0.587785);
  sasa_dot_locations[44].z(0.723607);
  sasa_dot_locations[45].x(-0.531939);
  sasa_dot_locations[45].y(-0.681718);
  sasa_dot_locations[45].z(0.502295);
  sasa_dot_locations[46].x(0.361803);
  sasa_dot_locations[46].y(-0.262866);
  sasa_dot_locations[46].z(0.894427);
  sasa_dot_locations[47].x(0.688191);
  sasa_dot_locations[47].y(-0.500000);
  sasa_dot_locations[47].z(0.525731);
  sasa_dot_locations[48].x(0.638197);
  sasa_dot_locations[48].y(-0.262866);
  sasa_dot_locations[48].z(0.723607);
  sasa_dot_locations[49].x(0.831052);
  sasa_dot_locations[49].y(-0.238856);
  sasa_dot_locations[49].z(0.502295);
  sasa_dot_locations[50].x(0.447214);
  sasa_dot_locations[50].y(-0.525731);
  sasa_dot_locations[50].z(0.723607);
  sasa_dot_locations[51].x(0.483974);
  sasa_dot_locations[51].y(-0.716567);
  sasa_dot_locations[51].z(0.502295);
  sasa_dot_locations[52].x(0.723607);
  sasa_dot_locations[52].y(0.525731);
  sasa_dot_locations[52].z(-0.447214);
  sasa_dot_locations[53].x(0.951057);
  sasa_dot_locations[53].y(0.309017);
  sasa_dot_locations[53].z(0.000000);
  sasa_dot_locations[54].x(0.956626);
  sasa_dot_locations[54].y(0.147621);
  sasa_dot_locations[54].z(0.251148);
  sasa_dot_locations[55].x(0.861803);
  sasa_dot_locations[55].y(0.425325);
  sasa_dot_locations[55].z(0.276393);
  sasa_dot_locations[56].x(0.587785);
  sasa_dot_locations[56].y(0.809017);
  sasa_dot_locations[56].z(0.000000);
  sasa_dot_locations[57].x(0.670820);
  sasa_dot_locations[57].y(0.688191);
  sasa_dot_locations[57].z(0.276393);
  sasa_dot_locations[58].x(0.436009);
  sasa_dot_locations[58].y(0.864188);
  sasa_dot_locations[58].z(0.251148);
  sasa_dot_locations[59].x(0.809017);
  sasa_dot_locations[59].y(0.587785);
  sasa_dot_locations[59].z(0.000000);
  sasa_dot_locations[60].x(0.860696);
  sasa_dot_locations[60].y(0.442863);
  sasa_dot_locations[60].z(-0.251148);
  sasa_dot_locations[61].x(0.687157);
  sasa_dot_locations[61].y(0.681718);
  sasa_dot_locations[61].z(-0.251148);
  sasa_dot_locations[62].x(-0.276393);
  sasa_dot_locations[62].y(0.850651);
  sasa_dot_locations[62].z(-0.447214);
  sasa_dot_locations[63].x(0.000000);
  sasa_dot_locations[63].y(1.000000);
  sasa_dot_locations[63].z(0.000000);
  sasa_dot_locations[64].x(0.155218);
  sasa_dot_locations[64].y(0.955423);
  sasa_dot_locations[64].z(0.251148);
  sasa_dot_locations[65].x(-0.138197);
  sasa_dot_locations[65].y(0.951057);
  sasa_dot_locations[65].z(0.276393);
  sasa_dot_locations[66].x(-0.587785);
  sasa_dot_locations[66].y(0.809017);
  sasa_dot_locations[66].z(0.000000);
  sasa_dot_locations[67].x(-0.447214);
  sasa_dot_locations[67].y(0.850651);
  sasa_dot_locations[67].z(0.276393);
  sasa_dot_locations[68].x(-0.687157);
  sasa_dot_locations[68].y(0.681718);
  sasa_dot_locations[68].z(0.251148);
  sasa_dot_locations[69].x(-0.309017);
  sasa_dot_locations[69].y(0.951057);
  sasa_dot_locations[69].z(0.000000);
  sasa_dot_locations[70].x(-0.155218);
  sasa_dot_locations[70].y(0.955423);
  sasa_dot_locations[70].z(-0.251148);
  sasa_dot_locations[71].x(-0.436009);
  sasa_dot_locations[71].y(0.864188);
  sasa_dot_locations[71].z(-0.251148);
  sasa_dot_locations[72].x(-0.894427);
  sasa_dot_locations[72].y(0.000000);
  sasa_dot_locations[72].z(-0.447214);
  sasa_dot_locations[73].x(-0.951057);
  sasa_dot_locations[73].y(0.309017);
  sasa_dot_locations[73].z(0.000000);
  sasa_dot_locations[74].x(-0.860696);
  sasa_dot_locations[74].y(0.442863);
  sasa_dot_locations[74].z(0.251148);
  sasa_dot_locations[75].x(-0.947214);
  sasa_dot_locations[75].y(0.162460);
  sasa_dot_locations[75].z(0.276393);
  sasa_dot_locations[76].x(-0.951057);
  sasa_dot_locations[76].y(-0.309017);
  sasa_dot_locations[76].z(0.000000);
  sasa_dot_locations[77].x(-0.947214);
  sasa_dot_locations[77].y(-0.162460);
  sasa_dot_locations[77].z(0.276393);
  sasa_dot_locations[78].x(-0.860696);
  sasa_dot_locations[78].y(-0.442863);
  sasa_dot_locations[78].z(0.251148);
  sasa_dot_locations[79].x(-1.000000);
  sasa_dot_locations[79].y(0.000000);
  sasa_dot_locations[79].z(0.000000);
  sasa_dot_locations[80].x(-0.956626);
  sasa_dot_locations[80].y(0.147621);
  sasa_dot_locations[80].z(-0.251148);
  sasa_dot_locations[81].x(-0.956626);
  sasa_dot_locations[81].y(-0.147621);
  sasa_dot_locations[81].z(-0.251148);
  sasa_dot_locations[82].x(-0.276393);
  sasa_dot_locations[82].y(-0.850651);
  sasa_dot_locations[82].z(-0.447214);
  sasa_dot_locations[83].x(-0.587785);
  sasa_dot_locations[83].y(-0.809017);
  sasa_dot_locations[83].z(0.000000);
  sasa_dot_locations[84].x(-0.687157);
  sasa_dot_locations[84].y(-0.681718);
  sasa_dot_locations[84].z(0.251148);
  sasa_dot_locations[85].x(-0.447214);
  sasa_dot_locations[85].y(-0.850651);
  sasa_dot_locations[85].z(0.276393);
  sasa_dot_locations[86].x(0.000000);
  sasa_dot_locations[86].y(-1.000000);
  sasa_dot_locations[86].z(0.000000);
  sasa_dot_locations[87].x(-0.138197);
  sasa_dot_locations[87].y(-0.951057);
  sasa_dot_locations[87].z(0.276393);
  sasa_dot_locations[88].x(0.155218);
  sasa_dot_locations[88].y(-0.955423);
  sasa_dot_locations[88].z(0.251148);
  sasa_dot_locations[89].x(-0.309017);
  sasa_dot_locations[89].y(-0.951057);
  sasa_dot_locations[89].z(0.000000);
  sasa_dot_locations[90].x(-0.436009);
  sasa_dot_locations[90].y(-0.864188);
  sasa_dot_locations[90].z(-0.251148);
  sasa_dot_locations[91].x(-0.155218);
  sasa_dot_locations[91].y(-0.955423);
  sasa_dot_locations[91].z(-0.251148);
  sasa_dot_locations[92].x(0.723607);
  sasa_dot_locations[92].y(-0.525731);
  sasa_dot_locations[92].z(-0.447214);
  sasa_dot_locations[93].x(0.587785);
  sasa_dot_locations[93].y(-0.809017);
  sasa_dot_locations[93].z(0.000000);
  sasa_dot_locations[94].x(0.436009);
  sasa_dot_locations[94].y(-0.864188);
  sasa_dot_locations[94].z(0.251148);
  sasa_dot_locations[95].x(0.670820);
  sasa_dot_locations[95].y(-0.688191);
  sasa_dot_locations[95].z(0.276393);
  sasa_dot_locations[96].x(0.951057);
  sasa_dot_locations[96].y(-0.309017);
  sasa_dot_locations[96].z(0.000000);
  sasa_dot_locations[97].x(0.861803);
  sasa_dot_locations[97].y(-0.425325);
  sasa_dot_locations[97].z(0.276393);
  sasa_dot_locations[98].x(0.956626);
  sasa_dot_locations[98].y(-0.147621);
  sasa_dot_locations[98].z(0.251148);
  sasa_dot_locations[99].x(0.809017);
  sasa_dot_locations[99].y(-0.587785);
  sasa_dot_locations[99].z(0.000000);
  sasa_dot_locations[100].x(0.687157);
  sasa_dot_locations[100].y(-0.681718);
  sasa_dot_locations[100].z(-0.251148);
  sasa_dot_locations[101].x(0.860696);
  sasa_dot_locations[101].y(-0.442863);
  sasa_dot_locations[101].z(-0.251148);
  sasa_dot_locations[102].x(0.262866);
  sasa_dot_locations[102].y(0.809017);
  sasa_dot_locations[102].z(-0.525731);
  sasa_dot_locations[103].x(0.531939);
  sasa_dot_locations[103].y(0.681718);
  sasa_dot_locations[103].z(-0.502295);
  sasa_dot_locations[104].x(0.447214);
  sasa_dot_locations[104].y(0.850651);
  sasa_dot_locations[104].z(-0.276393);
  sasa_dot_locations[105].x(0.309017);
  sasa_dot_locations[105].y(0.951057);
  sasa_dot_locations[105].z(0.000000);
  sasa_dot_locations[106].x(0.138197);
  sasa_dot_locations[106].y(0.951057);
  sasa_dot_locations[106].z(-0.276393);
  sasa_dot_locations[107].x(-0.029644);
  sasa_dot_locations[107].y(0.864188);
  sasa_dot_locations[107].z(-0.502295);
  sasa_dot_locations[108].x(-0.688191);
  sasa_dot_locations[108].y(0.500000);
  sasa_dot_locations[108].z(-0.525731);
  sasa_dot_locations[109].x(-0.483974);
  sasa_dot_locations[109].y(0.716567);
  sasa_dot_locations[109].z(-0.502295);
  sasa_dot_locations[110].x(-0.670820);
  sasa_dot_locations[110].y(0.688191);
  sasa_dot_locations[110].z(-0.276393);
  sasa_dot_locations[111].x(-0.809017);
  sasa_dot_locations[111].y(0.587785);
  sasa_dot_locations[111].z(0.000000);
  sasa_dot_locations[112].x(-0.861803);
  sasa_dot_locations[112].y(0.425325);
  sasa_dot_locations[112].z(-0.276393);
  sasa_dot_locations[113].x(-0.831052);
  sasa_dot_locations[113].y(0.238856);
  sasa_dot_locations[113].z(-0.502295);
  sasa_dot_locations[114].x(-0.688191);
  sasa_dot_locations[114].y(-0.500000);
  sasa_dot_locations[114].z(-0.525731);
  sasa_dot_locations[115].x(-0.831052);
  sasa_dot_locations[115].y(-0.238856);
  sasa_dot_locations[115].z(-0.502295);
  sasa_dot_locations[116].x(-0.861803);
  sasa_dot_locations[116].y(-0.425325);
  sasa_dot_locations[116].z(-0.276393);
  sasa_dot_locations[117].x(-0.809017);
  sasa_dot_locations[117].y(-0.587785);
  sasa_dot_locations[117].z(0.000000);
  sasa_dot_locations[118].x(-0.670820);
  sasa_dot_locations[118].y(-0.688191);
  sasa_dot_locations[118].z(-0.276393);
  sasa_dot_locations[119].x(-0.483974);
  sasa_dot_locations[119].y(-0.716567);
  sasa_dot_locations[119].z(-0.502295);
  sasa_dot_locations[120].x(0.262866);
  sasa_dot_locations[120].y(-0.809017);
  sasa_dot_locations[120].z(-0.525731);
  sasa_dot_locations[121].x(-0.029644);
  sasa_dot_locations[121].y(-0.864188);
  sasa_dot_locations[121].z(-0.502295);
  sasa_dot_locations[122].x(0.138197);
  sasa_dot_locations[122].y(-0.951057);
  sasa_dot_locations[122].z(-0.276393);
  sasa_dot_locations[123].x(0.309017);
  sasa_dot_locations[123].y(-0.951057);
  sasa_dot_locations[123].z(0.000000);
  sasa_dot_locations[124].x(0.447214);
  sasa_dot_locations[124].y(-0.850651);
  sasa_dot_locations[124].z(-0.276393);
  sasa_dot_locations[125].x(0.531939);
  sasa_dot_locations[125].y(-0.681718);
  sasa_dot_locations[125].z(-0.502295);
  sasa_dot_locations[126].x(0.850651);
  sasa_dot_locations[126].y(0.000000);
  sasa_dot_locations[126].z(-0.525731);
  sasa_dot_locations[127].x(0.812731);
  sasa_dot_locations[127].y(-0.295242);
  sasa_dot_locations[127].z(-0.502295);
  sasa_dot_locations[128].x(0.947214);
  sasa_dot_locations[128].y(-0.162460);
  sasa_dot_locations[128].z(-0.276393);
  sasa_dot_locations[129].x(1.000000);
  sasa_dot_locations[129].y(0.000000);
  sasa_dot_locations[129].z(0.000000);
  sasa_dot_locations[130].x(0.947214);
  sasa_dot_locations[130].y(0.162460);
  sasa_dot_locations[130].z(-0.276393);
  sasa_dot_locations[131].x(0.812731);
  sasa_dot_locations[131].y(0.295242);
  sasa_dot_locations[131].z(-0.502295);
  sasa_dot_locations[132].x(0.000000);
  sasa_dot_locations[132].y(0.000000);
  sasa_dot_locations[132].z(-1.000000);
  sasa_dot_locations[133].x(0.425325);
  sasa_dot_locations[133].y(0.309017);
  sasa_dot_locations[133].z(-0.850651);
  sasa_dot_locations[134].x(0.609548);
  sasa_dot_locations[134].y(0.442863);
  sasa_dot_locations[134].z(-0.657513);
  sasa_dot_locations[135].x(0.361803);
  sasa_dot_locations[135].y(0.587785);
  sasa_dot_locations[135].z(-0.723607);
  sasa_dot_locations[136].x(-0.162460);
  sasa_dot_locations[136].y(0.500000);
  sasa_dot_locations[136].z(-0.850651);
  sasa_dot_locations[137].x(0.052786);
  sasa_dot_locations[137].y(0.688191);
  sasa_dot_locations[137].z(-0.723607);
  sasa_dot_locations[138].x(-0.232827);
  sasa_dot_locations[138].y(0.716567);
  sasa_dot_locations[138].z(-0.657513);
  sasa_dot_locations[139].x(0.138197);
  sasa_dot_locations[139].y(0.425325);
  sasa_dot_locations[139].z(-0.894427);
  sasa_dot_locations[140].x(0.203183);
  sasa_dot_locations[140].y(0.147621);
  sasa_dot_locations[140].z(-0.967949);
  sasa_dot_locations[141].x(-0.077609);
  sasa_dot_locations[141].y(0.238856);
  sasa_dot_locations[141].z(-0.967949);
  sasa_dot_locations[142].x(-0.447214);
  sasa_dot_locations[142].y(0.525731);
  sasa_dot_locations[142].z(-0.723607);
  sasa_dot_locations[143].x(-0.525731);
  sasa_dot_locations[143].y(0.000000);
  sasa_dot_locations[143].z(-0.850651);
  sasa_dot_locations[144].x(-0.638197);
  sasa_dot_locations[144].y(0.262866);
  sasa_dot_locations[144].z(-0.723607);
  sasa_dot_locations[145].x(-0.753443);
  sasa_dot_locations[145].y(0.000000);
  sasa_dot_locations[145].z(-0.657513);
  sasa_dot_locations[146].x(-0.361803);
  sasa_dot_locations[146].y(0.262866);
  sasa_dot_locations[146].z(-0.894427);
  sasa_dot_locations[147].x(-0.251148);
  sasa_dot_locations[147].y(0.000000);
  sasa_dot_locations[147].z(-0.967949);
  sasa_dot_locations[148].x(-0.638197);
  sasa_dot_locations[148].y(-0.262866);
  sasa_dot_locations[148].z(-0.723607);
  sasa_dot_locations[149].x(-0.162460);
  sasa_dot_locations[149].y(-0.500000);
  sasa_dot_locations[149].z(-0.850651);
  sasa_dot_locations[150].x(-0.447214);
  sasa_dot_locations[150].y(-0.525731);
  sasa_dot_locations[150].z(-0.723607);
  sasa_dot_locations[151].x(-0.232827);
  sasa_dot_locations[151].y(-0.716567);
  sasa_dot_locations[151].z(-0.657513);
  sasa_dot_locations[152].x(-0.361803);
  sasa_dot_locations[152].y(-0.262866);
  sasa_dot_locations[152].z(-0.894427);
  sasa_dot_locations[153].x(-0.077609);
  sasa_dot_locations[153].y(-0.238856);
  sasa_dot_locations[153].z(-0.967949);
  sasa_dot_locations[154].x(0.052786);
  sasa_dot_locations[154].y(-0.688191);
  sasa_dot_locations[154].z(-0.723607);
  sasa_dot_locations[155].x(0.425325);
  sasa_dot_locations[155].y(-0.309017);
  sasa_dot_locations[155].z(-0.850651);
  sasa_dot_locations[156].x(0.361803);
  sasa_dot_locations[156].y(-0.587785);
  sasa_dot_locations[156].z(-0.723607);
  sasa_dot_locations[157].x(0.609548);
  sasa_dot_locations[157].y(-0.442863);
  sasa_dot_locations[157].z(-0.657513);
  sasa_dot_locations[158].x(0.138197);
  sasa_dot_locations[158].y(-0.425325);
  sasa_dot_locations[158].z(-0.894427);
  sasa_dot_locations[159].x(0.203183);
  sasa_dot_locations[159].y(-0.147621);
  sasa_dot_locations[159].z(-0.967949);
  sasa_dot_locations[160].x(0.670820);
  sasa_dot_locations[160].y(-0.162460);
  sasa_dot_locations[160].z(-0.723607);
  sasa_dot_locations[161].x(0.670820);
  sasa_dot_locations[161].y(0.162460);
  sasa_dot_locations[161].z(-0.723607);
  sasa_dot_locations[162].x(0.447214);
  sasa_dot_locations[162].y(0.000000);
  sasa_dot_locations[162].z(-0.894427);

  // locations are negative so reverse them!!!
  for(int ii = 1; ii <= 162; ii++) {
		sasa_dot_locations[ii] *= -1;
	}

	return sasa_dot_locations;

	// PackstatReal mn = 0.26;
	//
	// if ( false ) {
	// 	for (int ii=1; ii<=1; ii++) {
	// 		for( int jj = 1; jj <= 162; jj += 1 ) {
	// 			PackstatReal const xx = sasa_dot_locations(1,ii) + sasa_dot_locations(1,jj);
	// 			PackstatReal const yy = sasa_dot_locations(2,ii) + sasa_dot_locations(2,jj);
	// 			PackstatReal const zz = sasa_dot_locations(3,ii) + sasa_dot_locations(3,jj);
	// 			PackstatReal const dis = sqrt(xx*xx + yy*yy + zz*zz);
	// 			if ( 0 < dis && dis < mn ) {
	// 				cout << ii << ' ' << jj << ' ' << dis << endl;
	// 			}
	// 		}
	// 	}
	// }

	// bool perturb_dots = false;
	// if ( perturb_dots ) {
	//
	// 	FArray2D_float newsdl(3,162,-1234);
	//
	// 	PackstatReal rx,ry,rz;
	// 	PackstatReal r2 = 9999;
	// 	while ( r2 > 1.0f )  {
	// 		rx = ran3();
	// 		ry = ran3();
	// 		rz = ran3();
	// 		r2 = rx*rx + ry*ry + rz*rz;
	// 	}
	// 	PackstatReal r = sqrt(r2);
	// 	rx = rx/r;
	// 	ry = ry/r;
	// 	rz = rz/r;
	// 	PackstatReal const phi = acos( rz );
	// 	numeric::xyzVector_float const axis( rx,ry,0 );
	// 	numeric::xyzMatrix_float rot( rotation_matrix(axis,phi) );
	// 	for (int ii=1; ii <= 162; ++ii ) {
	// 		xyzVector_float x(sasa_dot_locations(1,ii));
	// 		x = rot * x;
	// 		newsdl(1,ii) = x.x();
	// 		newsdl(2,ii) = x.y();
	// 		newsdl(3,ii) = x.z();
	// 	}
	//
	// 	if ( false ) {
	// 		for (int ii=1; ii<=1; ii++) {
	// 			PackstatReal const xx = sasa_dot_locations(1,ii);
	// 			PackstatReal const yy = sasa_dot_locations(2,ii);
	// 			PackstatReal const zz = sasa_dot_locations(3,ii);
	// 			cout << F( 7,2, xx) <<  F( 7,2, yy) <<  F( 7,2, zz)
	// 				<< F( 7,2, newsdl(1,ii)) <<  F( 7,2, newsdl(2,ii)) <<  F( 7,2, newsdl(3,ii) )
	// 				<< endl;
	// 		}
	// 	}


	}


} // namespace packstat
} // namespace scoring
} // namespace core

#endif
