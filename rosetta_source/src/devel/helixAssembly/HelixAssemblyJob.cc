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

#include <devel/helixAssembly/HelixAssemblyJob.hh>

HelixAssemblyJob::HelixAssemblyJob(){};

//  HelixAssemblyJob clone() const {
//    return new HelixAssemblyJob( *this );
//  }

//copy constructor
HelixAssemblyJob::HelixAssemblyJob(HelixAssemblyJob const & old_job):
        job_id_(old_job.job_id_),
        job_name_(old_job.job_name_),
        remaining_rounds_(old_job.remaining_rounds_),
        query_structure_(old_job.query_structure_),
        search_structure_(old_job.search_structure_),
        search_index_(old_job.search_index_),
        query_frag_1_index_(old_job.query_frag_1_index_),
        query_frag_2_index_(old_job.query_frag_2_index_),
        fragments_(old_job.fragments_),
        direction_(old_job.direction_)
{}

HelixAssemblyJob::HelixAssemblyJob(core::Size job_id, std::string job_name, core::Size remaining_rounds, bool direction,
    std::string query_structure, std::string search_structure, core::Size search_index, core::Size query_frag_1_index,
    core::Size query_frag_2_index, std::vector<HelicalFragment> fragments):
    job_id_(job_id),
    job_name_(job_name),
    remaining_rounds_(remaining_rounds),
    direction_(direction),
    query_structure_(query_structure),
    search_structure_(search_structure),
    search_index_(search_index),
    query_frag_1_index_(query_frag_1_index),
    query_frag_2_index_(query_frag_2_index),
    fragments_(fragments)
{}

core::Size HelixAssemblyJob::get_job_id() const{
  return job_id_;
}

std::string HelixAssemblyJob::get_job_name() const{
  return job_name_;
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

bool HelixAssemblyJob::get_direction() const{
  return direction_;
}

HelicalFragment HelixAssemblyJob::get_query_frag_1() const{
  return fragments_[query_frag_1_index_];
}

HelicalFragment HelixAssemblyJob::get_query_frag_2() const{
  return fragments_[query_frag_2_index_];
}

std::vector<HelicalFragment> HelixAssemblyJob::get_fragments() const{
  return fragments_;
}

void HelixAssemblyJob::set_job_id(core::Size job_id){
  this->job_id_ = job_id;
}

void HelixAssemblyJob::set_job_name(std::string job_name){
  this->job_name_ = job_name;
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

void HelixAssemblyJob::set_fragments(std::vector<HelicalFragment> fragments){
  this->fragments_ = fragments;
}

void HelixAssemblyJob::set_direction(bool direction){
  this->direction_ = direction;
}
