// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/clusters/CDRClusterEnum.hh
/// @brief Enumerators for antibody clusters
/// @author Jared Adolf_Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_clusters_CDRClusterEnum_hh
#define INCLUDED_protocols_antibody_clusters_CDRClusterEnum_hh

namespace protocols {
namespace antibody {
namespace clusters {


/// Identified antibody CDR Clusters
/// See North, Lehmann, Dunbrack. (2011). JMB 406(2): 228-256.
/// More clusters will be added as the PDB grows.
enum CDRClusterEnum {
    H1_10_1 = 1,
    H1_12_1,
    H1_13_1,
    H1_13_2,
    H1_13_3,
    H1_13_4,
    H1_13_5,
    H1_13_6,
    H1_13_7,
    H1_13_8,
    H1_13_9,
    H1_13_10,
    H1_13_11,
    H1_13_cis9_1,
    H1_14_1,
    H1_15_1,
    H1_16_1,
    H2_8_1,
    H2_9_1,
    H2_9_2,
    H2_9_3,
    H2_10_1,
    H2_10_2,
    H2_10_3,
    H2_10_4,
    H2_10_5,
    H2_10_6,
    H2_10_7,
    H2_10_8,
    H2_10_9,
    H2_12_1,
    H2_15_1,
    H3_13_1,
    H3_10_1,
    H3_11_1,
    H3_12_1,
    H3_13_2,
    H3_14_1,
    H3_16_2,
    H3_12_2,
    H3_15_2,
    H3_14_3,
    H3_5_1,
    H3_11_2,
    H3_12_cis9_1,
    H3_9_1,
    H3_7_1,
    H3_13_3,
    H3_19_1,
    H3_16_1,
    H3_7_cis4_1,
    H3_6_1,
    H3_18_2,
    H3_21_1,
    H3_8_1,
    H3_17_1,
    H3_9_2,
    H3_7_3,
    H3_5_2,
    H3_10_2,
    H3_13_cis7_1,
    H3_10_cis5_1,
    H3_14_2,
    H3_7_2,
    H3_15_1,
    H3_10_3,
    H3_20_1,
    H3_14_cis7_1,
    H3_18_1,
    H3_24_1,
    H3_26_1,
    H3_8_2,
    H3_14_4,
    H3_13_cis8_1,
    H3_22_1,
    L1_10_1,
    L1_10_2,
    L1_11_1,
    L1_11_2,
    L1_11_3,
    L1_12_1,
    L1_12_2,
    L1_12_3,
    L1_13_1,
    L1_13_2,
    L1_14_1,
    L1_14_2,
    L1_15_1,
    L1_15_2,
    L1_16_1,
    L1_17_1,
    L2_8_1,
    L2_8_2,
    L2_8_3,
    L2_8_4,
    L2_8_5,
    L2_12_1,
    L2_12_2,
    L3_10_1,
    L3_11_1,
    L3_12_1,
    L3_13_1,
    L3_7_1,
    L3_8_1,
    L3_8_2,
    L3_8_cis6_1,
    L3_9_1,
    L3_9_2,
    L3_9_cis6_1,
    L3_9_cis7_1,
    L3_9_cis7_2,
    L3_9_cis7_3,
    L3_10_cis8_1,
    L3_10_cis7_8_1,
    L3_11_cis7_1,
    NA,
    CDRClusterEnum_start = H1_10_1,
    CDRClusterEnum_stop = L3_11_cis7_1,
    CDRClusterEnum_total = NA
};

} //clusters
} //antibody
} //protocols

#endif	//INCLUDED_protocols_antibody_CDRClusterEnum.hh

