// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/clusters/CDRClusterEnumManager.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/clusters/CDRClusterEnumManager.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>

#include <map>
#include <string>

#include <utility/vector1.hh>
#include <core/types.hh>
#include <utility/py/PyAssert.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace antibody {
namespace clusters {


CDRClusterEnumManager::CDRClusterEnumManager() {
	setup();
}

CDRClusterEnumManager::~CDRClusterEnumManager() = default;

void
CDRClusterEnumManager::setup() {

	enum_to_string_.resize(CDRClusterEnum_total);
	enum_to_string_[H1_10_1] = "H1-10-1";
	enum_to_string_[H1_12_1] = "H1-12-1";
	enum_to_string_[H1_13_1] = "H1-13-1";
	enum_to_string_[H1_13_2] = "H1-13-2";
	enum_to_string_[H1_13_3] = "H1-13-3 ";
	enum_to_string_[H1_13_4] = "H1-13-4";
	enum_to_string_[H1_13_5] = "H1-13-5 ";
	enum_to_string_[H1_13_6] = "H1-13-6";
	enum_to_string_[H1_13_7] = "H1-13-7";
	enum_to_string_[H1_13_8] = "H1-13-8";
	enum_to_string_[H1_13_9] = "H1-13-9";
	enum_to_string_[H1_13_10]="H1-13-10";
	enum_to_string_[H1_13_11]="H1-13-11";
	enum_to_string_[H1_13_cis9_1] = "H1-13-cis9-1";
	enum_to_string_[H1_14_1] = "H1-14-1";
	enum_to_string_[H1_15_1] = "H1-15-1";
	enum_to_string_[H1_16_1] = "H1-16-1";
	enum_to_string_[H2_8_1] = "H2-8-1";
	enum_to_string_[H2_10_1] = "H2-10-1";
	enum_to_string_[H2_10_2] = "H2-10-2";
	enum_to_string_[H2_10_3] = "H2-10-3";
	enum_to_string_[H2_10_4] = "H2-10-4";
	enum_to_string_[H2_10_5] = "H2-10-5";
	enum_to_string_[H2_10_6] = "H2-10-6";
	enum_to_string_[H2_10_7] = "H2-10-7";
	enum_to_string_[H2_10_8] = "H2-10-8";
	enum_to_string_[H2_10_9] = "H2-10-9";
	enum_to_string_[H2_12_1] = "H2-12-1";
	enum_to_string_[H2_15_1] = "H2-15-1";
	enum_to_string_[H2_9_1] = "H2-9-1";
	enum_to_string_[H2_9_2] = "H2-9-2";
	enum_to_string_[H2_9_3] = "H2-9-3";
	enum_to_string_[H3_13_1] = "H3-13-1";
	enum_to_string_[H3_10_1] = "H3-10-1";
	enum_to_string_[H3_11_1] = "H3-11-1";
	enum_to_string_[H3_12_1] = "H3-12-1";
	enum_to_string_[H3_13_2] = "H3-13-2";
	enum_to_string_[H3_14_1] = "H3-14-1";
	enum_to_string_[H3_16_2] = "H3-16-2";
	enum_to_string_[H3_12_2] = "H3-12-2";
	enum_to_string_[H3_15_2] = "H3-15-2";
	enum_to_string_[H3_14_3] = "H3-14-3";
	enum_to_string_[H3_5_1] = "H3-5-1";
	enum_to_string_[H3_11_2] = "H3-11-2";
	enum_to_string_[H3_12_cis9_1] = "H3-12-cis9-1";
	enum_to_string_[H3_9_1] = "H3-9-1";
	enum_to_string_[H3_7_1] = "H3-7-1";
	enum_to_string_[H3_13_3] = "H3-13-2";
	enum_to_string_[H3_19_1] = "H3-19-1";
	enum_to_string_[H3_16_1] = "H3-16-1";
	enum_to_string_[H3_7_cis4_1] = "H3-7-cis4-1";
	enum_to_string_[H3_6_1] = "H3-6-1";
	enum_to_string_[H3_18_2] = "H3-18-2";
	enum_to_string_[H3_21_1] = "H3-21-1";
	enum_to_string_[H3_8_1] = "H3-8-1";
	enum_to_string_[H3_17_1] = "H3-17-1";
	enum_to_string_[H3_9_2] = "H3-9-2";
	enum_to_string_[H3_7_3] = "H3-7-3";
	enum_to_string_[H3_5_2] = "H3-5-2";
	enum_to_string_[H3_10_2] = "H3-10-2";
	enum_to_string_[H3_13_cis7_1] = "H3-13-cis7-1";
	enum_to_string_[H3_10_cis5_1] = "H3-10-cis5-1";
	enum_to_string_[H3_14_2] = "H3-14-2";
	enum_to_string_[H3_7_2] = "H3-7-2";
	enum_to_string_[H3_15_1] = "H3-15-1";
	enum_to_string_[H3_10_3] = "H3-10-3";
	enum_to_string_[H3_20_1] = "H3-20-1";
	enum_to_string_[H3_14_cis7_1] = "H3-14-cis7-1";
	enum_to_string_[H3_18_1] = "H3-18-1";
	enum_to_string_[H3_24_1] = "H3-24-1";
	enum_to_string_[H3_26_1] = "H3-26-1";
	enum_to_string_[H3_8_2] = "H3-8-2";
	enum_to_string_[H3_14_4] = "H3-14-4";
	enum_to_string_[H3_13_cis8_1] = "H3-13-cis8-1";
	enum_to_string_[H3_22_1] = "H3-22-1";
	enum_to_string_[L1_10_1] = "L1-10-1";
	enum_to_string_[L1_10_2] = "L1-10-2";
	enum_to_string_[L1_11_1] = "L1-11-1";
	enum_to_string_[L1_11_2] = "L1-11-2";
	enum_to_string_[L1_11_3] = "L1-11-3";
	enum_to_string_[L1_12_1] = "L1-12-1";
	enum_to_string_[L1_12_2] = "L1-12-2";
	enum_to_string_[L1_12_3] = "L1-12-3";
	enum_to_string_[L1_13_1] = "L1-13-1";
	enum_to_string_[L1_13_2] = "L1-13-2";
	enum_to_string_[L1_14_1] = "L1-14-1";
	enum_to_string_[L1_14_2] = "L1-14-2";
	enum_to_string_[L1_15_1] = "L1-15-1";
	enum_to_string_[L1_15_2] = "L1-15-2";
	enum_to_string_[L1_16_1] = "L1-16-1";
	enum_to_string_[L1_17_1] = "L1-17-1";
	enum_to_string_[L2_8_1] = "L2-8-1";
	enum_to_string_[L2_8_2] = "L2-8-2";
	enum_to_string_[L2_8_3] = "L2-8-3";
	enum_to_string_[L2_8_4] = "L2-8-4";
	enum_to_string_[L2_8_5] = "L2-8-5";
	enum_to_string_[L2_12_1] = "L2-12-1";
	enum_to_string_[L2_12_2] = "L2-12-2";
	enum_to_string_[L3_10_1] = "L3-10-1";
	enum_to_string_[L3_11_1] = "L3-11-1";
	enum_to_string_[L3_12_1] = "L3-12-1";
	enum_to_string_[L3_13_1] = "L3-13-1";
	enum_to_string_[L3_7_1] = "L3-7-1";
	enum_to_string_[L3_8_1] = "L3-8-1";
	enum_to_string_[L3_8_2] = "L3-8-2";
	enum_to_string_[L3_8_cis6_1] = "L3-8-cis6-1";
	enum_to_string_[L3_9_1] = "L3-9-1";
	enum_to_string_[L3_9_2] = "L3-9-2";
	enum_to_string_[L3_9_cis6_1] = "L3-9-cis6-1";
	enum_to_string_[L3_9_cis7_1] = "L3-9-cis7-1";
	enum_to_string_[L3_9_cis7_2] = "L3-9-cis7-2";
	enum_to_string_[L3_9_cis7_3] = "L3-9-cis7-3";
	enum_to_string_[L3_10_cis8_1] = "L3-10-cis8-1";
	enum_to_string_[L3_10_cis7_8_1] = "L3-10-cis7,8-1";
	enum_to_string_[L3_11_cis7_1] = "L3-11-cis7-1";
	enum_to_string_[NA] = "NA";
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	string_to_enum_["H1-10-1"] = H1_10_1;
	string_to_enum_["H1-12-1"] = H1_12_1;
	string_to_enum_["H1-13-1"] = H1_13_1;
	string_to_enum_["H1-13-2"] = H1_13_2;
	string_to_enum_["H1-13-3"] = H1_13_3 ;
	string_to_enum_["H1-13-4"] = H1_13_4;
	string_to_enum_["H1-13-5"] = H1_13_5 ;
	string_to_enum_["H1-13-6"] = H1_13_6;
	string_to_enum_["H1-13-7"] = H1_13_7;
	string_to_enum_["H1-13-8"] = H1_13_8;
	string_to_enum_["H1-13-9"] = H1_13_9;
	string_to_enum_["H1-13-10"] = H1_13_10;
	string_to_enum_["H1-13-11"] = H1_13_11;
	string_to_enum_["H1-13-cis9-1"] = H1_13_cis9_1;
	string_to_enum_["H1-13-CIS9-1"] = H1_13_cis9_1;
	string_to_enum_["H1-14-1"] = H1_14_1;
	string_to_enum_["H1-15-1"] = H1_15_1;
	string_to_enum_["H1-16-1"] = H1_16_1;
	string_to_enum_["H2-8-1"] = H2_8_1;
	string_to_enum_["H2-10-1"] = H2_10_1;
	string_to_enum_["H2-10-2"] = H2_10_2;
	string_to_enum_["H2-10-3"] = H2_10_3;
	string_to_enum_["H2-10-4"] = H2_10_4;
	string_to_enum_["H2-10-5"] = H2_10_5;
	string_to_enum_["H2-10-6"] = H2_10_6;
	string_to_enum_["H2-10-7"] = H2_10_7;
	string_to_enum_["H2-10-8"] = H2_10_8;
	string_to_enum_["H2-10-9"] = H2_10_9;
	string_to_enum_["H2-12-1"] = H2_12_1;
	string_to_enum_["H2-15-1"] = H2_15_1;
	string_to_enum_["H2-9-1"] = H2_9_1;
	string_to_enum_["H2-9-2"] = H2_9_2;
	string_to_enum_["H2-9-3"] = H2_9_3;
	string_to_enum_["H3-13-1"] = H3_13_1;
	string_to_enum_["H3-10-1"] = H3_10_1;
	string_to_enum_["H3-11-1"] = H3_11_1;
	string_to_enum_["H3-12-1"] = H3_12_1;
	string_to_enum_["H3-13-2"] = H3_13_2;
	string_to_enum_["H3-14-1"] = H3_14_1;
	string_to_enum_["H3-16-2"] = H3_16_2;
	string_to_enum_["H3-12-2"] = H3_12_2;
	string_to_enum_["H3-15-2"] = H3_15_2;
	string_to_enum_["H3-14-3"] = H3_14_3;
	string_to_enum_["H3-5-1"] = H3_5_1;
	string_to_enum_["H3-11-2"] = H3_11_2;
	string_to_enum_["H3-12-cis9-1"] = H3_12_cis9_1;
	string_to_enum_["H3-12-CIS9-1"] = H3_12_cis9_1;
	string_to_enum_["H3-9-1"] = H3_9_1;
	string_to_enum_["H3-7-1"] = H3_7_1;
	string_to_enum_["H3-13-3"] = H3_13_2;
	string_to_enum_["H3-19-1"] = H3_19_1;
	string_to_enum_["H3-16-1"] = H3_16_1;
	string_to_enum_["H3-7-cis4-1"] = H3_7_cis4_1;
	string_to_enum_["H3-7-CIS4-1"] = H3_7_cis4_1;
	string_to_enum_["H3-6-1"] = H3_6_1;
	string_to_enum_["H3-18-2"] = H3_18_2;
	string_to_enum_["H3-21-1"] = H3_21_1;
	string_to_enum_["H3-8-1"] = H3_8_1;
	string_to_enum_["H3-17-1"] = H3_17_1;
	string_to_enum_["H3-9-2"] = H3_9_2;
	string_to_enum_["H3-7-3"] = H3_7_3;
	string_to_enum_["H3-5-2"] = H3_5_2;
	string_to_enum_["H3-10-2"] = H3_10_2;
	string_to_enum_["H3-13-cis7-1"] = H3_13_cis7_1;
	string_to_enum_["H3-13-CIS7-1"] = H3_13_cis7_1;
	string_to_enum_["H3-10-cis5-1"] = H3_10_cis5_1;
	string_to_enum_["H3-10-CIS5-1"] = H3_10_cis5_1;
	string_to_enum_["H3-14-2"] = H3_14_2;
	string_to_enum_["H3-7-2"] = H3_7_2;
	string_to_enum_["H3-15-1"] = H3_15_1;
	string_to_enum_["H3-10-3"] = H3_10_3;
	string_to_enum_["H3-20-1"] = H3_20_1;
	string_to_enum_["H3-14-cis7-1"] = H3_14_cis7_1;
	string_to_enum_["H3-14-CIS7-1"] = H3_14_cis7_1;
	string_to_enum_["H3-18-1"] = H3_18_1;
	string_to_enum_["H3-24-1"] = H3_24_1;
	string_to_enum_["H3-26-1"] = H3_26_1;
	string_to_enum_["H3-8-2"] = H3_8_2;
	string_to_enum_["H3-14-4"] = H3_14_4;
	string_to_enum_["H3-13-cis8-1"] = H3_13_cis8_1;
	string_to_enum_["H3-13-CIS8-1"] = H3_13_cis8_1;
	string_to_enum_["H3-22-1"] = H3_22_1;
	string_to_enum_["L1-10-1"] = L1_10_1;
	string_to_enum_["L1-10-2"] = L1_10_2;
	string_to_enum_["L1-11-1"] = L1_11_1;
	string_to_enum_["L1-11-2"] = L1_11_2;
	string_to_enum_["L1-11-3"] = L1_11_3;
	string_to_enum_["L1-12-1"] = L1_12_1;
	string_to_enum_["L1-12-2"] = L1_12_2;
	string_to_enum_["L1-12-3"] = L1_12_3;
	string_to_enum_["L1-13-1"] = L1_13_1;
	string_to_enum_["L1-13-2"] = L1_13_2;
	string_to_enum_["L1-14-1"] = L1_14_1;
	string_to_enum_["L1-14-2"] = L1_14_2;
	string_to_enum_["L1-15-1"] = L1_15_1;
	string_to_enum_["L1-15-2"] = L1_15_2;
	string_to_enum_["L1-16-1"] = L1_16_1;
	string_to_enum_["L1-17-1"] = L1_17_1;
	string_to_enum_["L2-8-1"] = L2_8_1;
	string_to_enum_["L2-8-2"] = L2_8_2;
	string_to_enum_["L2-8-3"] = L2_8_3;
	string_to_enum_["L2-8-4"] = L2_8_4;
	string_to_enum_["L2-8-5"] = L2_8_5;
	string_to_enum_["L2-12-1"] = L2_12_1;
	string_to_enum_["L2-12-2"] = L2_12_2;
	string_to_enum_["L3-10-1"] = L3_10_1;
	string_to_enum_["L3-11-1"] = L3_11_1;
	string_to_enum_["L3-12-1"] = L3_12_1;
	string_to_enum_["L3-13-1"] = L3_13_1;
	string_to_enum_["L3-7-1"] = L3_7_1;
	string_to_enum_["L3-8-1"] = L3_8_1;
	string_to_enum_["L3-8-2"] = L3_8_2;
	string_to_enum_["L3-8-cis6-1"] = L3_8_cis6_1;
	string_to_enum_["L3-8-CIS6-1"] = L3_8_cis6_1;
	string_to_enum_["L3-9-1"] = L3_9_1;
	string_to_enum_["L3-9-2"] = L3_9_2;
	string_to_enum_["L3-9-cis6-1"] = L3_9_cis6_1;
	string_to_enum_["L3-9-cis7-1"] = L3_9_cis7_1;
	string_to_enum_["L3-9-cis7-2"] = L3_9_cis7_2;
	string_to_enum_["L3-9-cis7-3"] = L3_9_cis7_3;
	string_to_enum_["L3-10-cis8-1"] = L3_10_cis8_1;
	string_to_enum_["L3-10-cis7,8-1"] = L3_10_cis7_8_1;
	string_to_enum_["L3-11-cis7-1"] = L3_11_cis7_1;
	string_to_enum_["L3-9-CIS6-1"] = L3_9_cis6_1;
	string_to_enum_["L3-9-CIS7-1"] = L3_9_cis7_1;
	string_to_enum_["L3-9-CIS7-2"] = L3_9_cis7_2;
	string_to_enum_["L3-9-CIS7-3"] = L3_9_cis7_3;
	string_to_enum_["L3-10-CIS8-1"] = L3_10_cis8_1;
	string_to_enum_["L3-10-CIS7,8-1"] = L3_10_cis7_8_1;
	string_to_enum_["L3-11-CIS7-1"] = L3_11_cis7_1;
	string_to_enum_["NA"] = NA;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	string_to_enum_["H1_10_1"] = H1_10_1;
	string_to_enum_["H1_12_1"] = H1_12_1;
	string_to_enum_["H1_13_1"] = H1_13_1;
	string_to_enum_["H1_13_2"] = H1_13_2;
	string_to_enum_["H1_13_3"] = H1_13_3 ;
	string_to_enum_["H1_13_4"] = H1_13_4;
	string_to_enum_["H1_13_5"] = H1_13_5 ;
	string_to_enum_["H1_13_6"] = H1_13_6;
	string_to_enum_["H1_13_7"] = H1_13_7;
	string_to_enum_["H1_13_8"] = H1_13_8;
	string_to_enum_["H1_13_9"] = H1_13_9;
	string_to_enum_["H1_13_10"] = H1_13_10;
	string_to_enum_["H1_13_11"] = H1_13_11;
	string_to_enum_["H1_13_cis9_1"] = H1_13_cis9_1;
	string_to_enum_["H1_13_CIS9_1"] = H1_13_cis9_1;
	string_to_enum_["H1_14_1"] = H1_14_1;
	string_to_enum_["H1_15_1"] = H1_15_1;
	string_to_enum_["H1_16_1"] = H1_16_1;
	string_to_enum_["H2_8_1"] = H2_8_1;
	string_to_enum_["H2_10_1"] = H2_10_1;
	string_to_enum_["H2_10_2"] = H2_10_2;
	string_to_enum_["H2_10_3"] = H2_10_3;
	string_to_enum_["H2_10_4"] = H2_10_4;
	string_to_enum_["H2_10_5"] = H2_10_5;
	string_to_enum_["H2_10_6"] = H2_10_6;
	string_to_enum_["H2_10_7"] = H2_10_7;
	string_to_enum_["H2_10_8"] = H2_10_8;
	string_to_enum_["H2_10_9"] = H2_10_9;
	string_to_enum_["H2_12_1"] = H2_12_1;
	string_to_enum_["H2_15_1"] = H2_15_1;
	string_to_enum_["H2_9_1"] = H2_9_1;
	string_to_enum_["H2_9_2"] = H2_9_2;
	string_to_enum_["H2_9_3"] = H2_9_3;
	string_to_enum_["H3_13_1"] = H3_13_1;
	string_to_enum_["H3_10_1"] = H3_10_1;
	string_to_enum_["H3_11_1"] = H3_11_1;
	string_to_enum_["H3_12_1"] = H3_12_1;
	string_to_enum_["H3_13_2"] = H3_13_2;
	string_to_enum_["H3_14_1"] = H3_14_1;
	string_to_enum_["H3_16_2"] = H3_16_2;
	string_to_enum_["H3_12_2"] = H3_12_2;
	string_to_enum_["H3_15_2"] = H3_15_2;
	string_to_enum_["H3_14_3"] = H3_14_3;
	string_to_enum_["H3_5_1"] = H3_5_1;
	string_to_enum_["H3_11_2"] = H3_11_2;
	string_to_enum_["H3_12_cis9_1"] = H3_12_cis9_1;
	string_to_enum_["H3_12_CIS9_1"] = H3_12_cis9_1;
	string_to_enum_["H3_9_1"] = H3_9_1;
	string_to_enum_["H3_7_1"] = H3_7_1;
	string_to_enum_["H3_13_3"] = H3_13_2;
	string_to_enum_["H3_19_1"] = H3_19_1;
	string_to_enum_["H3_16_1"] = H3_16_1;
	string_to_enum_["H3_7_cis4_1"] = H3_7_cis4_1;
	string_to_enum_["H3_7_CIS4_1"] = H3_7_cis4_1;
	string_to_enum_["H3_6_1"] = H3_6_1;
	string_to_enum_["H3_18_2"] = H3_18_2;
	string_to_enum_["H3_21_1"] = H3_21_1;
	string_to_enum_["H3_8_1"] = H3_8_1;
	string_to_enum_["H3_17_1"] = H3_17_1;
	string_to_enum_["H3_9_2"] = H3_9_2;
	string_to_enum_["H3_7_3"] = H3_7_3;
	string_to_enum_["H3_5_2"] = H3_5_2;
	string_to_enum_["H3_10_2"] = H3_10_2;
	string_to_enum_["H3_13_cis7_1"] = H3_13_cis7_1;
	string_to_enum_["H3_13_CIS7_1"] = H3_13_cis7_1;
	string_to_enum_["H3_10_cis5_1"] = H3_10_cis5_1;
	string_to_enum_["H3_10_CIS5_1"] = H3_10_cis5_1;
	string_to_enum_["H3_14_2"] = H3_14_2;
	string_to_enum_["H3_7_2"] = H3_7_2;
	string_to_enum_["H3_15_1"] = H3_15_1;
	string_to_enum_["H3_10_3"] = H3_10_3;
	string_to_enum_["H3_20_1"] = H3_20_1;
	string_to_enum_["H3_14_cis7_1"] = H3_14_cis7_1;
	string_to_enum_["H3_14_CIS7_1"] = H3_14_cis7_1;
	string_to_enum_["H3_18_1"] = H3_18_1;
	string_to_enum_["H3_24_1"] = H3_24_1;
	string_to_enum_["H3_26_1"] = H3_26_1;
	string_to_enum_["H3_8_2"] = H3_8_2;
	string_to_enum_["H3_14_4"] = H3_14_4;
	string_to_enum_["H3_13_cis8_1"] = H3_13_cis8_1;
	string_to_enum_["H3_13_CIS8_1"] = H3_13_cis8_1;
	string_to_enum_["H3_22_1"] = H3_22_1;
	string_to_enum_["L1_10_1"] = L1_10_1;
	string_to_enum_["L1_10_2"] = L1_10_2;
	string_to_enum_["L1_11_1"] = L1_11_1;
	string_to_enum_["L1_11_2"] = L1_11_2;
	string_to_enum_["L1_11_3"] = L1_11_3;
	string_to_enum_["L1_12_1"] = L1_12_1;
	string_to_enum_["L1_12_2"] = L1_12_2;
	string_to_enum_["L1_12_3"] = L1_12_3;
	string_to_enum_["L1_13_1"] = L1_13_1;
	string_to_enum_["L1_13_2"] = L1_13_2;
	string_to_enum_["L1_14_1"] = L1_14_1;
	string_to_enum_["L1_14_2"] = L1_14_2;
	string_to_enum_["L1_15_1"] = L1_15_1;
	string_to_enum_["L1_15_2"] = L1_15_2;
	string_to_enum_["L1_16_1"] = L1_16_1;
	string_to_enum_["L1_17_1"] = L1_17_1;
	string_to_enum_["L2_8_1"] = L2_8_1;
	string_to_enum_["L2_8_2"] = L2_8_2;
	string_to_enum_["L2_8_3"] = L2_8_3;
	string_to_enum_["L2_8_4"] = L2_8_4;
	string_to_enum_["L2_8_5"] = L2_8_5;
	string_to_enum_["L2_12_1"] = L2_12_1;
	string_to_enum_["L2_12_2"] = L2_12_2;
	string_to_enum_["L3_10_1"] = L3_10_1;
	string_to_enum_["L3_11_1"] = L3_11_1;
	string_to_enum_["L3_12_1"] = L3_12_1;
	string_to_enum_["L3_13_1"] = L3_13_1;
	string_to_enum_["L3_7_1"] = L3_7_1;
	string_to_enum_["L3_8_1"] = L3_8_1;
	string_to_enum_["L3_8_2"] = L3_8_2;
	string_to_enum_["L3_8_cis6_1"] = L3_8_cis6_1;
	string_to_enum_["L3_8_CIS6_1"] = L3_8_cis6_1;
	string_to_enum_["L3_9_1"] = L3_9_1;
	string_to_enum_["L3_9_2"] = L3_9_2;
	string_to_enum_["L3_9_cis6_1"] = L3_9_cis6_1;
	string_to_enum_["L3_9_cis7_1"] = L3_9_cis7_1;
	string_to_enum_["L3_9_cis7_2"] = L3_9_cis7_2;
	string_to_enum_["L3_10_cis8_1"] = L3_10_cis8_1;
	string_to_enum_["L3_10_cis7_8_1"] = L3_10_cis7_8_1;
	string_to_enum_["L3_11_cis7_1"] = L3_11_cis7_1;
	string_to_enum_["L3_9_CIS6_1"] = L3_9_cis6_1;
	string_to_enum_["L3_9_CIS7_1"] = L3_9_cis7_1;
	string_to_enum_["L3_9_CIS7_2"] = L3_9_cis7_2;
	string_to_enum_["L3_10_CIS8_1"] = L3_10_cis8_1;
	string_to_enum_["L3_10_CIS7_8_1"] = L3_10_cis7_8_1;
	string_to_enum_["L3_11_CIS7_1"] = L3_11_cis7_1;
}

std::string
CDRClusterEnumManager::cdr_cluster_enum_to_string(CDRClusterEnum const cluster) const {

	//Help protect from bad memory access that will give the enum crazy values.
	if ( cluster > CDRClusterEnum_total ) { //|| cluster < 0 ) {
		utility_exit_with_message("Bogus CDRClusterEnum passed to cdr_cluster_enum_to_string" + utility::to_string(int(cluster)));
	}
	return enum_to_string_[cluster];
}

CDRClusterEnum
CDRClusterEnumManager::cdr_cluster_string_to_enum(std::string const & cluster) const {

	//This is here due to const correctness issues with [] operator
	auto iter( string_to_enum_.find( cluster ) );
	if ( iter == string_to_enum_.end() ) {
		utility_exit_with_message("Cluster not found: " + cluster);
	}
	return iter->second;
}

bool
CDRClusterEnumManager::cdr_cluster_is_present(std::string const & cluster) const {
	auto iter( string_to_enum_.find( cluster ) );
	return iter != string_to_enum_.end();
}

utility::vector1<std::string>
CDRClusterEnumManager::get_recognized_cluster_definitions() const {
	utility::vector1<std::string> recognized_clusters;
	for ( auto& cluster : string_to_enum_ ) {
		recognized_clusters.push_back(cluster.first);
	}

	return recognized_clusters;
}

}//clusters
}//Antibody
}//Protocols

