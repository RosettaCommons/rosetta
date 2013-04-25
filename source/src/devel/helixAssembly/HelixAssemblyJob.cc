// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file HelixAssemblyJob.cc
///
/// @brief

/// @author Tim jacobs

//Unit
#include <devel/helixAssembly/HelixAssemblyJob.hh>

HelixAssemblyJob::HelixAssemblyJob():
id_(0),
name_(""),
remaining_rounds_(0),
query_structure_(""),
search_structure_(""),
search_index_(0),
query_frag_1_index_(0),
query_frag_2_index_(0),
fragments_(),
direction_needed_(true),
first_round_(false),
n_term_growth_(false)
{};
HelixAssemblyJob::~HelixAssemblyJob(){};

core::Size HelixAssemblyJob::get_id() const{
  return id_;
}

std::string HelixAssemblyJob::get_name() const{
  return name_;
}

core::Size HelixAssemblyJob::get_remaining_rounds() const{
  return remaining_rounds_;
}

core::Size HelixAssemblyJob::get_search_index() const{
  return search_index_;
}

std::string HelixAssemblyJob::get_query_structure() const{
  return query_structure_;
}

std::string HelixAssemblyJob::get_search_structure() const{
  return search_structure_;
}

core::Size HelixAssemblyJob::get_query_frag_1_index() const{
  return query_frag_1_index_;
}

core::Size HelixAssemblyJob::get_query_frag_2_index() const{
  return query_frag_2_index_;
}

bool HelixAssemblyJob::get_direction_needed() const{
  return direction_needed_;
}

bool HelixAssemblyJob::get_first_round() const{
  return first_round_;
}

bool HelixAssemblyJob::get_n_term_growth() const{
  return n_term_growth_;
}

protocols::features::helixAssembly::HelicalFragment HelixAssemblyJob::get_query_frag_1() const{
  return fragments_[query_frag_1_index_];
}

protocols::features::helixAssembly::HelicalFragment HelixAssemblyJob::get_query_frag_2() const{
  return fragments_[query_frag_2_index_];
}

std::vector<protocols::features::helixAssembly::HelicalFragment> HelixAssemblyJob::get_fragments() const{
  return fragments_;
}

void HelixAssemblyJob::add_fragment(protocols::features::helixAssembly::HelicalFragment new_fragment){
  fragments_.push_back(new_fragment);
}

void HelixAssemblyJob::set_id(core::Size id){
  this->id_ = id;
}

void HelixAssemblyJob::set_name(std::string name){
  this->name_ = name;
}

void HelixAssemblyJob::set_remaining_rounds(core::Size remaining_rounds){
  this->remaining_rounds_ = remaining_rounds;
}

void HelixAssemblyJob::set_query_structure(std::string query_structure){
  this->query_structure_ = query_structure;
}

void HelixAssemblyJob::set_search_structure(std::string search_structure){
  this->search_structure_ = search_structure;
}

void HelixAssemblyJob::set_search_index(core::Size search_index){
  this->search_index_ = search_index;
}

void HelixAssemblyJob::set_query_frag_1_index(core::Size query_frag_1_index){
  this->query_frag_1_index_ = query_frag_1_index;
}

void HelixAssemblyJob::set_query_frag_2_index(core::Size query_frag_2_index){
  this->query_frag_2_index_ = query_frag_2_index;
}

void HelixAssemblyJob::set_fragments(std::vector<protocols::features::helixAssembly::HelicalFragment> fragments){
  this->fragments_ = fragments;
}

void HelixAssemblyJob::set_direction_needed(bool direction_needed){
  this->direction_needed_ = direction_needed;
}

void HelixAssemblyJob::set_first_round(bool first_round){
  this->first_round_=first_round;
}

void HelixAssemblyJob::set_n_term_growth(bool n_term_growth){
  this->n_term_growth_=n_term_growth;
}

std::string HelixAssemblyJob::printBundleResidues() const{

  std::string output = "BUNDLE: " + name_ + "\n";

  for(core::Size i=0; i<fragments_.size(); ++i){
      output+=fragments_[i].print() + "\n";
  }
  return output;
}
