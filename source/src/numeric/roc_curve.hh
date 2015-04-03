// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/roc_curve.hh
/// @author Sam DeLuca

#ifndef INCLUDED_numeric_roc_curve_HH
#define INCLUDED_numeric_roc_curve_HH

//unit headers
#include <numeric/roc_curve.fwd.hh>

//platform headers
#include <platform/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

//stl headers
#include <string>
namespace numeric {


class RocPoint : public utility::pointer::ReferenceCount {

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~RocPoint();
	RocPoint(bool predicted, bool actual, std::string const & tag,platform::Real const & score);
	RocPoint(RocStatus const & status,std::string const & tag, platform::Real const & score);

	RocStatus status() const;
	void status(RocStatus const & status);

	std::string tag() const;
	void tag(std::string const & tag);

	platform::Real score() const;
	void score(platform::Real const & score);


	bool operator<(RocPoint const & that) const;

private:
	RocStatus status_;
	std::string tag_;
	platform::Real score_;

};

class RocCurve : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~RocCurve();
	RocCurve();

	void insert_point(RocPoint const & roc_point);
	void insert_point(bool predicted, bool actual, std::string const & tag, platform::Real const & score);
	void generate_roc_curve();

	utility::vector1<std::pair<platform::Real, platform::Real> > roc_curve();

	void print_roc_curve();
	platform::Real calculate_auc();

private:

	platform::Size true_positive_count_;
	platform::Size false_positive_count_;
	platform::Size true_negative_count_;
	platform::Size false_negative_count_;

	utility::vector1<RocPoint> roc_point_vector_;
	utility::vector1<std::pair<platform::Real, platform::Real> > roc_curve_;

};

}


#endif /* INCLUDED_numeric_roc_curve_HH */
