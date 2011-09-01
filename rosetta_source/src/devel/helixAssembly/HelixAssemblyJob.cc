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

HelixAssemblyJob::HelixAssemblyJob(HelixAssemblyJob const & old_job):
        job_id_(old_job.job_id_),
        job_name_(old_job.job_name_),
        round_(old_job.round_),
        query_structure_(old_job.query_structure_),
        search_structure_(old_job.search_structure_),
        frag1_start_(old_job.frag1_start_),
        frag1_end_(old_job.frag1_end_),
        frag2_start_(old_job.frag2_start_),
        frag2_end_(old_job.frag2_end_)
{}

HelixAssemblyJob::HelixAssemblyJob(core::Size job_id, std::string job_name, core::Size round, std::string query_structure, std::string search_structure,
    core::Size frag1_start, core::Size frag1_end, core::Size frag2_start, core::Size frag2_end):
    job_id_(job_id),
    job_name_(job_name),
    round_(round),
    query_structure_(query_structure),
    search_structure_(search_structure),
    frag1_start_(frag1_start),
    frag1_end_(frag1_end),
    frag2_start_(frag2_start),
    frag2_end_(frag2_end)
{}

core::Size HelixAssemblyJob::get_job_id() const{
  return job_id_;
}

std::string HelixAssemblyJob::get_job_name() const{
  return job_name_;
}

core::Size HelixAssemblyJob::get_round() const{
  return round_;
}

core::Size HelixAssemblyJob::get_frag1_end() const{
  return frag1_end_;
}

core::Size HelixAssemblyJob::get_frag1_start() const{
  return frag1_start_;
}

core::Size HelixAssemblyJob::get_frag2_end() const{
  return frag2_end_;
}

core::Size HelixAssemblyJob::get_frag2_start() const{
  return frag2_start_;
}

std::string HelixAssemblyJob::get_query_structure() const{
  return query_structure_;
}

std::string HelixAssemblyJob::get_search_structure() const{
  return search_structure_;
}

void HelixAssemblyJob::set_job_id(core::Size job_id){
  this->job_id_ = job_id;
}

void HelixAssemblyJob::set_job_name(std::string job_name){
  this->job_name_ = job_name;
}

void HelixAssemblyJob::set_round(core::Size round){
  this->round_ = round;
}

void HelixAssemblyJob::set_frag1_end(core::Size frag1_end){
  this->frag1_end_ = frag1_end;
}

void HelixAssemblyJob::set_frag1_start(core::Size frag1_start){
  this->frag1_start_ = frag1_start;
}

void HelixAssemblyJob::set_frag2_end(core::Size frag2_end){
  this->frag2_end_ = frag2_end;
}

void HelixAssemblyJob::set_frag2_start(core::Size frag2_start){
  this->frag2_start_ = frag2_start;
}

void HelixAssemblyJob::set_query_structure(std::string query_structure){
  this->query_structure_ = query_structure;
}

void HelixAssemblyJob::set_search_structure(std::string search_structure){
  this->search_structure_ = search_structure;
}
