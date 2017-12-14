// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/roc_curve.cc
/// @author Sam DeLuca

//stl headers
#include <iostream>

//unit headers
#include <numeric/roc_curve.hh>

//utility headers
#include <utility>
#include <utility/exit.hh>


//stl headers
#include <string>
#include <algorithm>

namespace numeric {

/// @details Auto-generated virtual destructor
RocCurve::~RocCurve() = default;

/// @details Auto-generated virtual destructor
RocPoint::~RocPoint() = default;


RocPoint::RocPoint(
	bool predicted,
	bool actual,
	std::string const & tag,
	platform::Real const & score
) :
	tag_(tag),
	score_(score)
{
	if ( predicted && actual ) {
		status_ = true_positive;
	} else if ( predicted && !actual ) {
		status_ = false_positive;
	} else if ( !predicted && !actual ) {
		status_ = true_negative;
	} else {
		status_ = false_negative;
	}
}

RocPoint::RocPoint(
	RocStatus const & status,
	std::string const & tag,
	platform::Real const & score
) :
	status_(status),
	tag_(tag),
	score_(score)
{
	//
}

RocStatus RocPoint::status() const
{
	return status_;
}

void RocPoint::status(RocStatus const & status)
{
	status_ = status;
}

std::string RocPoint::tag() const
{
	return tag_;
}

void RocPoint::tag(std::string const & tag)
{
	tag_ = tag;
}

platform::Real RocPoint::score() const
{
	return score_;
}

void RocPoint::score(platform::Real const & score)
{
	score_ = score;
}


bool RocPoint::operator <(RocPoint const & that) const
{
	return this->score_ < that.score_;
}

RocCurve::RocCurve() :
	true_positive_count_(0),
	false_positive_count_(0),
	true_negative_count_(0),
	false_negative_count_(0)
{

}

void RocCurve::insert_point(RocPoint const & roc_point)
{
	RocStatus status(roc_point.status());
	roc_point_vector_.push_back(roc_point);

	switch(status)
			{
			case true_positive :
				++true_positive_count_;
				break;
			case false_positive :
				++false_positive_count_;
				break;
			case true_negative :
				++true_negative_count_;
				break;
			case false_negative :
				++false_negative_count_;
				break;
			default :
				utility_exit_with_message("invalid RocPoint status, I have no idea how this happened");
				break;
			}

}

void RocCurve::insert_point(bool predicted, bool actual, std::string const & tag, platform::Real const & score)
{
	RocPoint new_point(predicted,actual,tag,score);
	insert_point(new_point);
}

void RocCurve::generate_roc_curve()
{
	roc_curve_.clear();

	platform::Real true_positives = 0.0;
	platform::Real false_negatives = 0.0;

	platform::Real false_positives = 0.0;
	platform::Real true_negatives = 0.0;


	//push an initial 0,0 point into the array
	roc_curve_.push_back(std::make_pair(0.0,0.0));

	std::sort(roc_point_vector_.begin(),roc_point_vector_.end());

	for ( auto & roc_it : roc_point_vector_ ) {
		RocStatus status(roc_it.status());
		switch(status)
				{
				case true_positive :
					++true_positives;
					break;
				case false_positive :
					++false_positives;
					break;
				case true_negative :
					++true_negatives;
					break;
				case false_negative :
					++false_negatives;
					break;
				default :
					utility_exit_with_message("invalid RocPoint status, I have no idea how this happened");
				}

		platform::Real true_positive_rate = 0.0;
		platform::Real false_positive_rate = 0.0;
		if ( (true_positive_count_+false_negative_count_) != 0 ) {
			true_positive_rate = true_positives / (true_positive_count_+ false_negative_count_);
		}

		if ( (false_positive_count_+ true_negative_count_) != 0 ) {
			false_positive_rate = false_positives / (false_positive_count_ + true_negative_count_);
		}

		roc_curve_.push_back(std::make_pair(true_positive_rate,false_positive_rate));
	}

	roc_curve_.push_back(std::make_pair(1.0,1.0));
}

utility::vector1<std::pair<platform::Real, platform::Real> > RocCurve::roc_curve()
{
	return roc_curve_;
}

void RocCurve::print_roc_curve()
{
	for ( auto & auc_it : roc_curve_ ) {
		std::cout << auc_it.first << " " <<auc_it.second <<std::endl;
	}
}

platform::Real RocCurve::calculate_auc()
{
	platform::Real AUC = 0.0;
	platform::Real last_TPR = 0.0;
	platform::Real last_FPR = 0.0;
	for ( auto & auc_it : roc_curve_ ) {

		platform::Real TPR = auc_it.first;
		platform::Real FPR = auc_it.second;
		AUC += (FPR-last_FPR)*((TPR+last_TPR)/2);
		//AUC += auc_it->first*(auc_it->second - last_FPR);
		last_TPR = TPR;
		last_FPR = FPR;

	}

	return AUC;
}

}


