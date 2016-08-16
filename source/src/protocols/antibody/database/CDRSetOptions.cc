// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/AntibodyCDRSetInstructions.hh
/// @brief Create and hold AntibodyCDRSetInstructions for the graft design step
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#include <protocols/antibody/database/CDRSetOptions.hh>
#include <protocols/antibody/AntibodyEnum.hh>



namespace protocols {
namespace antibody {

using namespace protocols::antibody::clusters;

CDRSetOptions::CDRSetOptions(bool load):
	utility::pointer::ReferenceCount(),
	load_(load)
{
	set_defaults();
}
CDRSetOptions::CDRSetOptions(CDRNameEnum cdr, bool load):
	utility::pointer::ReferenceCount(),
	cdr_(cdr),
	load_(load)
{
	set_defaults();
}

CDRSetOptions::CDRSetOptions(const CDRSetOptions& src):
	utility::pointer::ReferenceCount(src),
	cdr_(src.cdr_),
	load_(src.load_),
	only_current_cluster_(src.only_current_cluster_),
	only_center_clusters_(src.only_center_clusters_),
	length_types_(src.length_types_),
	exclude_pdb_ids_(src.exclude_pdb_ids_),
	include_only_pdb_ids_(src.include_only_pdb_ids_),
	exclude_clusters_(src.exclude_clusters_),
	include_only_clusters_(src.include_only_clusters_),
	min_length_(src.min_length_),
	max_length_(src.max_length_),
	exclude_species_(src.exclude_species_),
	include_only_species_(src.include_only_species_),
	exclude_germlines_(src.exclude_germlines_),
	include_only_germlines_(src.include_only_germlines_),
	sampling_cutoff_(src.sampling_cutoff_)
{

}

CDRSetOptions::~CDRSetOptions(){}

void
CDRSetOptions::set_defaults() {
	length_types_.clear();
	length_types_.resize(3, true);
	min_length_ = 1;
	max_length_ = 50;
	only_current_cluster_ = false;
	only_center_clusters_ = false;
	exclude_pdb_ids_.clear();
	include_only_pdb_ids_.clear();
	exclude_clusters_.clear();
	include_only_clusters_.clear();
	exclude_species_.clear();
	include_only_species_.clear();
	exclude_germlines_.clear();
	include_only_germlines_.clear();

	sampling_cutoff_ = 0;

}

CDRSetOptionsOP
CDRSetOptions::clone() const {
	return CDRSetOptionsOP( new CDRSetOptions(*this) );
}

void
CDRSetOptions::set_cdr(CDRNameEnum cdr){
	cdr_ = cdr;
}

void
CDRSetOptions::load(bool load) {
	load_ = load;
}

void
CDRSetOptions::include_only_current_cluster(bool only_current_cluster) {
	only_current_cluster_  = only_current_cluster;
}

bool
CDRSetOptions::include_only_current_cluster() const {
	return only_current_cluster_;
}

void
CDRSetOptions::include_only_center_clusters(bool centers_only) {
	only_center_clusters_ = centers_only;
}

bool
CDRSetOptions::include_only_center_clusters() const{
	return only_center_clusters_;
}

void
CDRSetOptions::exclude_clusters(utility::vector1<CDRClusterEnum> exclude) {
	exclude_clusters_ = exclude;
}

utility::vector1<CDRClusterEnum>
CDRSetOptions::exclude_clusters() const {
	return exclude_clusters_;
}

void
CDRSetOptions::exclude_clusters_add(CDRClusterEnum exclude) {
	exclude_clusters_.push_back(exclude);
}

void
CDRSetOptions::exclude_clusters_clear(){
	exclude_clusters_.clear();
}

void
CDRSetOptions::exclude_pdbs(utility::vector1<std::string> exclude_pdbids) {
	exclude_pdb_ids_ = exclude_pdbids;
}

utility::vector1<std::string>
CDRSetOptions::exclude_pdbs() const{
	return exclude_pdb_ids_;
}

void
CDRSetOptions::exclude_pdbs_add(std::string pdbid) {
	exclude_pdb_ids_.push_back(pdbid);
}

void
CDRSetOptions::exclude_pdbs_clear(){
	exclude_pdb_ids_.clear();
}

void
CDRSetOptions::exclude_species(utility::vector1<std::string> exclude) {
	exclude_species_ = exclude;
}

utility::vector1<std::string>
CDRSetOptions::exclude_species() const {
	return exclude_species_;
}

void
CDRSetOptions::exclude_species_add(std::string exclude) {
	exclude_species_.push_back(exclude);
}

void
CDRSetOptions::exclude_species_clear(){
	exclude_species_.clear();
}

void
CDRSetOptions::exclude_germlines(utility::vector1<std::string> exclude) {
	exclude_germlines_ = exclude;
}

utility::vector1<std::string>
CDRSetOptions::exclude_germlines() const {
	return exclude_germlines_;
}

void
CDRSetOptions::exclude_germlines_add(std::string exclude){
	exclude_germlines_.push_back(exclude);
}

void
CDRSetOptions::exclude_germlines_clear(){
	exclude_germlines_.clear();
}

void
CDRSetOptions::include_only_clusters(utility::vector1<CDRClusterEnum> include_only) {
	include_only_clusters_ = include_only;
}

utility::vector1<CDRClusterEnum>
CDRSetOptions::include_only_clusters() const {
	return include_only_clusters_;
}

void
CDRSetOptions::include_only_clusters_add(CDRClusterEnum include_only) {
	include_only_clusters_.push_back(include_only);
}

void
CDRSetOptions::include_only_clusters_clear(){
	include_only_clusters_.clear();
}

void
CDRSetOptions::include_only_germlines(utility::vector1<std::string> include_only) {
	include_only_germlines_ = include_only;
}

utility::vector1<std::string>
CDRSetOptions::include_only_germlines() const {
	return include_only_germlines_;
}

void
CDRSetOptions::include_only_germlines_add(std::string include_only) {
	include_only_germlines_.push_back(include_only);
}

void
CDRSetOptions::include_only_germlines_clear() {
	include_only_germlines_.clear();
}

void
CDRSetOptions::include_only_pdbs(utility::vector1<std::string> include_pdbids){
	include_only_pdb_ids_ = include_pdbids;
}

utility::vector1<std::string>
CDRSetOptions::include_only_pdbs() const {
	return include_only_pdb_ids_;
}

void
CDRSetOptions::include_only_pdbs_add(std::string pdbid){
	include_only_pdb_ids_.push_back(pdbid);
}

void
CDRSetOptions::include_only_pdbs_clear(){
	include_only_pdb_ids_.clear();
}

void
CDRSetOptions::include_only_species(utility::vector1<std::string> include_only) {
	include_only_species_ = include_only;
}

utility::vector1<std::string>
CDRSetOptions::include_only_species() const {
	return include_only_species_;
}

void
CDRSetOptions::include_only_species_add(std::string include_only) {
	include_only_species_.push_back(include_only);
}

void
CDRSetOptions::include_only_species_clear(){
	include_only_species_.clear();
}

void
CDRSetOptions::length_type(const core::Size type, const bool setting) {
	length_types_[type] = setting;
}

void
CDRSetOptions::max_length(core::Size length){
	max_length_ = length;
}

void
CDRSetOptions::min_length(core::Size length) {
	min_length_ =length;
}

void
CDRSetOptions::cluster_sampling_cutoff(core::Size cutoff){

	sampling_cutoff_ = cutoff;

}

core::Size
CDRSetOptions::cluster_sampling_cutoff() const {
	return sampling_cutoff_;
}

}
}
