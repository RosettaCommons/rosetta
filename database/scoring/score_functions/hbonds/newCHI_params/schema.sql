-- -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
-- vi: set ts=2 noet:
--
-- (c) Copyright Rosetta Commons Member Institutions.
-- (c) This file is part of the Rosetta software suite and is made available under license.
-- (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
-- (c) For more information, see http://www.rosettacommons.org. Questions about this can be
-- (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
--
-- @file scoring/score_functions/hbonds/newCHI_params/schema.sql
-- @brief Schema for Hydrogen bonds parameters database
-- @author Andrew Leaver-Fay (aleaverfay@gmail.com), Matthew OMeara (mattjomeara@gmail.com)
--



DROP TABLE IF EXISTS HBondWeightType;
CREATE TABLE HBondWeightType (
    id INTEGER PRIMARY KEY,
    name TEXT,
    comment TEXT);

DROP TABLE IF EXISTS HBAccChemType;
CREATE TABLE HBAccChemType (
    id INTEGER PRIMARY KEY,
    name TEXT,
    name_long TEXT,
    comment TEXT);

DROP TABLE IF EXISTS HBDonChemType;
CREATE TABLE HBDonChemType (
    id INTEGER PRIMARY KEY,
    name TEXT,
    name_long TEXT,
    comment TEXT);

DROP TABLE IF EXISTS HBAccHybridization;
CREATE TABLE HBAccHybridization (
    acc_chem_type TEXT PRIMARY KEY,
    hybridization TEXT,
    comment TEXT);

DROP TABLE IF EXISTS HBSeqSep;
CREATE TABLE HBSeqSep (
    id INTEGER PRIMARY KEY,
    name TEXT,
    comment TEXT);

DROP TABLE IF EXISTS HybridizationType;
CREATE TABLE HybridizationType (
    id INTEGER PRIMARY KEY,
    name TEXT,
    comment TEXT);

DROP TABLE IF EXISTS HBPoly1D;
CREATE TABLE HBPoly1D (
    id INTEGER PRIMARY KEY,
    name TEXT,
    classic_name TEXT,
    dimension TEXT,
    xmin REAL,
    xmax REAL,
    min_val REAL,
    max_val REAL,
    root1 REAL,
    root2 REAL,
    degree INTEGER,
    c_a REAL,
    c_b REAL,
    c_c REAL,
    c_d REAL,
    c_e REAL,
    c_f REAL,
    c_g REAL,
    c_h REAL,
    c_i REAL,
    c_j REAL,
    c_k REAL);

DROP TABLE IF EXISTS HBEval;
CREATE TABLE HBEval (
    don_chem_type TEXT,
    acc_chem_type TEXT,
    separation TEXT,
    AHdist_short_fade TEXT,
    AHdist_long_fade TEXT,
    cosBAH_fade TEXT,
    cosAHD_fade TEXT,
    chi_fade TEXT,
    AHdist TEXT,
    cosBAH_short TEXT,
    cosBAH_long TEXT,
    cosAHD_short TEXT,
    cosAHD_long TEXT,
    --chi TEXT,
    weight_type TEXT,
    comment TEXT);

DROP TABLE IF EXISTS HBFadeIntervals;
CREATE TABLE HBFadeIntervals(
		id INTEGER PRIMARY KEY,
    name TEXT,
		junction_type TEXT,
		min0 REAL,
		fmin REAL,
		fmax REAL,
		max0 REAL,
		comment TEXT);


--
-- Load Data into tables
--
-- Note that the formate is very picky: no spaces after the separating ',' are allowed!
--
.separator ','
.import HybridizationType.csv HybridizationType
.import HBAccChemType.csv HBAccChemType
.import HBAccHybridization.csv HBAccHybridization
.import HBDonChemType.csv HBDonChemType
.import HBEval.csv HBEval
.import HBPoly1D.csv HBPoly1D
.import HBSeqSep.csv HBSeqSep
.import HBondWeightType.csv HBondWeightType
.import HBFadeIntervals.csv HBFadeIntervals
