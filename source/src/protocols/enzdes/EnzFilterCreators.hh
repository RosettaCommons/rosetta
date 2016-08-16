// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/enzdes/EnzFilterCreators.hh
/// @brief  FilterCreators for the enzyme design filters.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_enzdes_EnzFilterCreators_hh
#define INCLUDED_protocols_enzdes_EnzFilterCreators_hh

// Package Headers
#include <protocols/filters/FilterCreator.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace protocols {
namespace enzdes {

class DiffAtomSasaFilterCreator : public filters::FilterCreator
{
public:
	virtual filters::FilterOP create_filter() const;
	virtual std::string keyname() const;
};


class EnzScoreFilterCreator : public filters::FilterCreator
{
public:
	virtual filters::FilterOP create_filter() const;
	virtual std::string keyname() const;
};

class LigBurialFilterCreator : public filters::FilterCreator
{
public:
	virtual filters::FilterOP create_filter() const;
	virtual std::string keyname() const;
};

class LigDSasaFilterCreator : public filters::FilterCreator
{
public:
	virtual filters::FilterOP create_filter() const;
	virtual std::string keyname() const;
};

class LigInterfaceEnergyFilterCreator : public filters::FilterCreator
{
public:
	virtual filters::FilterOP create_filter() const;
	virtual std::string keyname() const;
};

class RepackWithoutLigandFilterCreator : public filters::FilterCreator
{
public:
	virtual filters::FilterOP create_filter() const;
	virtual std::string keyname() const;
};

class EnzdesScorefileFilterCreator : public filters::FilterCreator
{
public:
	virtual filters::FilterOP create_filter() const;
	virtual std::string keyname() const;
};

class ResidueConformerFilterCreator : public filters::FilterCreator
{
public:
	virtual filters::FilterOP create_filter() const;
	virtual std::string keyname() const;

};

} //namespace enzdes
} //namespace protocols

#endif
