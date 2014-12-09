// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/HeaderInformation.cc
///
/// @brief Information stored in the HEADER record in the PDB format
/// @author Matthew O'Meara

// Unit headers
#include <core/io/pdb/HeaderInformation.hh>
#include <core/io/pdb/Field.hh>

// Platform headers
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// Boost headers
#include <boost/lexical_cast.hpp>

// C++ Headers
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>


namespace core {
namespace io {
namespace pdb {

using std::string;
using std::list;
using std::endl;
using std::pair;
using ObjexxFCL::rstrip_whitespace;   // by reference
using ObjexxFCL::strip_whitespace;    // by reference
using ObjexxFCL::rstripped_whitespace;// copy
using ObjexxFCL::stripped_whitespace; // copy

static thread_local basic::Tracer TR( "core.io.pdb.HeaderInformation" );

HeaderInformation::HeaderInformation() : utility::pointer::ReferenceCount(),
	classification_(""),
	dep_year_(0),
	dep_month_(0),
	dep_day_(0),
	idCode_(""),
	title_(""),
	keywords_(),
	keyword_in_progress_(false),
	compounds_(),
	compound_in_progress_(false),
	experimental_techniques_(),
	experimental_technique_in_progress_("")
{}

HeaderInformation::HeaderInformation(
	HeaderInformation const & src) : utility::pointer::ReferenceCount(),
	classification_(src.classification_),
	dep_year_(src.dep_year_),
	dep_month_(src.dep_month_),
	dep_day_(src.dep_day_),
	idCode_(src.idCode_),
	title_(""),
	keywords_(src.keywords_),
	keyword_in_progress_(src.keyword_in_progress_),
	compounds_(src.compounds_),
	compound_in_progress_(src.compound_in_progress_),
	experimental_techniques_(src.experimental_techniques_),
	experimental_technique_in_progress_(src.experimental_technique_in_progress_)
{}

HeaderInformation::~HeaderInformation() {}


void
HeaderInformation::store_record(Record & R){
	string const & type = R["type"].value;
	if(type == "HEADER"){
		store_classification(R["classification"].value);
		store_deposition_date(R["depDate"].value);
		store_idCode(R["idCode"].value);
	} else if(type == "TITLE "){
		store_title(R["title"].value);
	} else if(type == "KEYWDS"){
		store_keywords(R["keywords"].value);
	} else if(type == "COMPND"){
		store_compound(R["compound"].value);
	} else if(type == "EXPDTA"){
		store_experimental_techniques(R["technique"].value);
	} else {
		std::stringstream err_msg;
		err_msg
			<< "Attempting to add unrecognized record type '" << type << "' "
			<< "to header information.";
		utility_exit_with_message(err_msg.str());
	}
}

void
HeaderInformation::finalize_parse() {
	finalize_keyword_records();
	finalize_compound_records();
	finalize_experimental_technique_records();
}

bool
HeaderInformation::parse_in_progress() const {
	return
		keyword_in_progress() ||
		compound_in_progress() ||
		experimental_technique_in_progress();
}

void
HeaderInformation::fill_records(
	std::vector<Record> & VR
) const {
	fill_header_record(VR);
	fill_title_records(VR);
	fill_keyword_records(VR);
	fill_compound_records(VR);
	fill_experimental_technique_records(VR);
}

////////////// HEADER ///////////////////

void
HeaderInformation::store_classification(string const & classification){
	classification_ = classification;
	rstrip_whitespace(classification_);

	// TODO: and that the classification is on the list
	// http://www.wwpdb.org/documentation/wwpdb20070104appendices_c.pdf

}

string
HeaderInformation::classification() const {
	return classification_;
}

void
HeaderInformation::store_deposition_date(string const & depDate) {

	dep_day_ = atoi(depDate.substr(0,2).c_str());
	if(dep_day_ > 31 || dep_day_ < 1){
		TR.Warning << "Deposition day not in range [1, 31]: " << dep_day_ << endl;
	}

	string const & mon(depDate.substr(3,3));
	if( mon == "JAN" ) dep_month_ = 1;
	else if( mon == "FEB" ) dep_month_ = 2;
	else if( mon == "MAR" ) dep_month_ = 3;
	else if( mon == "APR" ) dep_month_ = 4;
	else if( mon == "MAY" ) dep_month_ = 5;
	else if( mon == "JUN" ) dep_month_ = 6;
	else if( mon == "JUL" ) dep_month_ = 7;
	else if( mon == "AUG" ) dep_month_ = 8;
	else if( mon == "SEP" ) dep_month_ = 9;
	else if( mon == "OCT" ) dep_month_ = 10;
	else if( mon == "NOV" ) dep_month_ = 11;
	else if( mon == "DEC" ) dep_month_ = 12;
	else {
		TR.Warning << "Unrecognized month in HEADER deposition date " + depDate << mon << std::endl;
	}

	dep_year_ = boost::lexical_cast<Size>(depDate.substr(7,4));
}


void
HeaderInformation::store_deposition_date(
	Size yy,
	Size mm,
	Size dd
) {

	dep_year_ = yy;
	if(dep_month_ > 99 || dep_day_ < 1){
		TR.Warning << "Deposition month not in range [01, 99]: " << dep_month_ << endl;
	}

	dep_month_ = mm;
	if(dep_month_ > 12 || dep_day_ < 1){
		TR.Warning << "Deposition month not in range [1, 12]: " << dep_month_ << endl;
	}

	dep_day_ = dd;
	if(dep_day_ > 31 || dep_day_ < 1){
		TR.Warning << "Deposition day not in range [1, 31]: " << dep_day_ << endl;
	}
}


string
HeaderInformation::deposition_date() const {
	std::stringstream dep_date;

	if(dep_day_ > 31 || dep_day_ < 1){
		utility_exit_with_message("deposition day is outside of range [1,31]: " +
				boost::lexical_cast<std::string>(dep_day_));
	}
	dep_date << dep_day_ << "-";
	switch(dep_month_){
	case 1: dep_date  << "JAN"; break;
	case 2: dep_date  << "FEB"; break;
	case 3: dep_date  << "MAR"; break;
	case 4: dep_date  << "APR"; break;
	case 5: dep_date  << "MAY"; break;
	case 6: dep_date  << "JUN"; break;
	case 7: dep_date  << "JUL"; break;
	case 8: dep_date  << "AUG"; break;
	case 9: dep_date  << "SEP"; break;
	case 10: dep_date << "OCT"; break;
	case 11: dep_date << "NOV"; break;
	case 12: dep_date << "DEC"; break;
	default:
		utility_exit_with_message("Unrecognized deposition month index " +
				boost::lexical_cast<std::string>(dep_month_));
	}
	if( dep_year_ > 99 || dep_year_ < 1){
		utility_exit_with_message("Deposition year is out side of range [01,99]: " +
				boost::lexical_cast<std::string>(dep_year_));
	}
	dep_date << "-" << (dep_year_ < 10 ? "0" : "") << dep_year_;
	return dep_date.str();
}

void
HeaderInformation::deposition_date(
	Size & yy,
	Size & mm,
	Size & dd
) const {
	yy = dep_year_;
	mm = dep_month_;
	dd = dep_day_;
}

std::string
HeaderInformation::idCode() const {
	return idCode_;
}

void
HeaderInformation::store_idCode(string const & idCode) {
	idCode_ = idCode;
}

void
HeaderInformation::fill_header_record(
	std::vector< Record > & VR
) const {
	if(!classification_.empty() &&
		dep_year_ && dep_month_ && dep_day_ && !idCode_.empty()) {
		Record R = Field::getRecordCollection()["HEADER"];
		R["type"].value = "HEADER";
		R["classification"].value = classification();
		R["depDate"].value = deposition_date();
		R["idCode"].value = idCode();
		VR.push_back(R);
	}
}


////////////// TITLE ///////////////////

/// @details Append title, strip off white space on the left for the
/// first record and on the right for all records.
void
HeaderInformation::store_title(string const & title){
	if(title.empty()){
		TR.Warning << "Attempting to store empty title record field." << endl;
		return;
	}

	if(title_.empty()) {
		title_ = title;
		strip_whitespace(title_);
	} else {

		title_.append(rstripped_whitespace(title));
	}
}

void
HeaderInformation::clear_title(){
	title_.clear();
}

std::string const &
HeaderInformation::title() const {
	return title_;
}

void
HeaderInformation::fill_title_records(
	std::vector< Record > & VR
) const {

	if(!title_.empty())	{
		Size line_no(1);
		fill_wrapped_records("TITLE ", "title", title_, line_no, VR);
	}
}

////////////// KEYWDS ///////////////////

void
HeaderInformation::store_keywords(string const & keywords){
	if(keywords.empty()){
		TR.Warning << "Attempting to add empty keywords string." << endl;
		return;
	}

	size_t i(keywords.find_first_not_of(' '));
	size_t j(i);
	while(i != std::string::npos) {
		j = keywords.find(',', i);
		if(keyword_in_progress_){
			keywords_.back().append(
				" " + rstripped_whitespace(keywords.substr(i, j-i)));
			keyword_in_progress_ = false;
		} else {
			keywords_.push_back(rstripped_whitespace(keywords.substr(i, j-i)));
		}
		if(j != std::string::npos){
			i = keywords.find_first_not_of(' ', j+1);
		} else {
			keyword_in_progress_ = true;
			return;
		}
	}
}

list< string > const &
HeaderInformation::keywords() const {
	return keywords_;
}

void
HeaderInformation::finalize_keyword_records() {
	keyword_in_progress_ = false;
}

bool
HeaderInformation::keyword_in_progress() const {
	return keyword_in_progress_;
}

void
HeaderInformation::clear_keywords() {
	keywords_.clear();
}

void
HeaderInformation::fill_keyword_records(
	std::vector< Record > & VR
) const {
	if(keywords_.empty()) return;

	string keywords;
	list< string >::const_iterator k = keywords_.begin(), ke = keywords_.end();
	for(; k!= ke; ++k){
		if(!keywords.empty()) keywords.append(", ");
		keywords.append(*k);
	}
	Size line_no(1);
	fill_wrapped_records("KEYWDS", "keywords", keywords, line_no, VR);
}

///////////// COMPND ///////////////////

std::string
HeaderInformation::compound_token_to_string(CompoundToken token) {
	string token_str;
	switch(token){
	case MOL_ID:        token_str = "MOL_ID";        break;
	case MOLECULE:      token_str = "MOLECULE";      break;
	case CHAIN:         token_str = "CHAIN";         break;
	case FRAGMENT:      token_str = "FRAGMENT";      break;
	case SYNONYM:       token_str = "SYNONYM";       break;
	case EC:            token_str = "EC";            break;
	case ENGINEERED:    token_str = "ENGINEERED";    break;
	case MUTATION:      token_str = "MUTATION";      break;
	case OTHER_DETAILS: token_str = "OTHER_DETAILS"; break;
	default:
		TR.Error << "Unrecognized compound token '" << token << "'" << endl;
		utility_exit();
	}
	return token_str;
}

HeaderInformation::CompoundToken
HeaderInformation::string_to_compound_token(std::string const & token) {
	if(token == "MOL_ID")        return MOL_ID;
	else if(token == "MOLECULE")      return MOLECULE;
	else if(token == "CHAIN")         return CHAIN;
	else if(token == "FRAGMENT")      return FRAGMENT;
	else if(token == "SYNONYM")       return SYNONYM;
	else if(token == "EC")            return EC;
	else if(token == "ENGINEERED")    return ENGINEERED;
	else if(token == "MUTATION")      return MUTATION;
	else if(token == "OTHER_DETAILS") return OTHER_DETAILS;
	else {
		TR.Error << "Unrecognized compound token string '" << token << "'" << endl;
		utility_exit();
	}
	return CompoundToken_max;
}


/// @details Assume each new compound token/value pair begins on a new
/// line but the value can be multiple lines. So, if a compound record
/// is encountered when "in progress" then append the results to the
/// value of the previous pair.
void
HeaderInformation::store_compound(std::string const & compound) {

	if(compound_in_progress_){
		size_t v_end(compound.find(';'));
		compound_in_progress_ = (v_end == std::string::npos);
		compounds_[compounds_.size()].second.append(
			rstripped_whitespace(compound.substr(0,v_end)));
		return;
	}

	size_t t_begin(compound.find_first_not_of(' '));
	size_t t_end(compound.find(':', t_begin));
	if(t_end == std::string::npos) {
		TR.Error
			<< "Attempting to add compound to header information "
			<< "but no compund token was found in '" << compound << "'" << endl;
		utility_exit();
	}
	CompoundToken token(
		string_to_compound_token(compound.substr(t_begin,t_end - t_begin)));

	size_t v_begin(compound.find_first_not_of(' ', t_end + 1));
	if(v_begin == std::string::npos){
		TR.Error
			<< "Attempting to add compound to header information "
			<< "but no compund value was found in '" << compound << "'" << endl;
		utility_exit();
	}
	size_t v_end(compound.find(';', v_begin));
	if(v_end == std::string::npos){
		compound_in_progress_ = true;
	}
	compounds_.push_back(
		make_pair(token, compound.substr(v_begin, v_end - v_begin)));
}

void
HeaderInformation::store_compound(
	HeaderInformation::CompoundToken token,
	string const & value
) {
	compounds_.push_back(make_pair(token, value));
}

utility::vector1< pair< HeaderInformation::CompoundToken, string > > const &
HeaderInformation::compounds() const{
	return compounds_;
}

void
HeaderInformation::finalize_compound_records() {
	compound_in_progress_ = false;
}

bool
HeaderInformation::compound_in_progress() const {
	return compound_in_progress_;
}

void
HeaderInformation::clear_compounds() {
	compounds_.clear();
}

void
HeaderInformation::fill_compound_records(
	std::vector< Record > & VR
) const {

	Size line_no(1);

	for(Size t=1, te = compounds_.size(); t <= te; ++t){

		std::stringstream comp_field;
		comp_field
			// defacto standard in PDB is to add a space after a continuation field
			<< (line_no == 1 ? "" : " ")
			<< compound_token_to_string(
				static_cast<CompoundToken>(compounds_[t].first))
			<< ": "
			<< compounds_[t].second
			// only add ';' to separate compound records
			<< (t < compounds_.size() ? ";" : "");

		fill_wrapped_records("COMPND", "compound", comp_field.str(), line_no, VR);
	}
}

////////////// EXPDTA ///////////////////

string
HeaderInformation::experimental_technique_to_string(
	ExperimentalTechnique technique
) {
	string t;
	switch(technique){
	case X_RAY_DIFFRACTION:        t = "X-RAY DIFFRACTION";        break;
	case FIBER_DIFFRACTION:        t = "FIBER DIFFRACTION";        break;
	case NEUTRON_DIFFRACTION:      t = "NEUTRON DIFFRACTION";      break;
	case ELECTRON_CRYSTALLOGRAPHY: t = "ELECTRON CRYSTALLOGRAPHY"; break;
	case ELECTRON_MICROSCOPY:      t = "ELECTRON MICROSCOPY";      break;
	case SOLID_STATE_NMR:          t = "SOLID-STATE NMR";          break;
	case SOLUTION_NMR:             t = "SOLUTION NMR";             break;
	case SOLUTION_SCATTERING:      t = "SOLUTION SCATTERING";      break;
	case THEORETICAL_MODEL:        t = "THEORETICAL MODEL";        break;

	case ELECTRON_DEFRACTION:
		t = "ELECTRON DEFRACTION";
		TR.Warning
			<< "Encountered obsolete experimental technqiue coding '"
			<< t << "'" << endl;
		break;

	case CRYO_ELECTRON_MICROSCOPY:
		t = "CRYO-ELECTRON MICROSCOPY";
		TR.Warning
			<< "Encountered obsolete experimental technqiue coding '"
			<< t << "'" << endl;
		break;

	case SOLUTION_SCATTERING_THEORETICAL_MODEL:
		t = "SOLUTION SCATTERING, THEORETICAL MODEL";
		TR.Warning
			<< "Encountered obsolete experimental technqiue coding '"
			<< t << "'" << endl;
		break;

	case FLORECENCE_TRANSFER:
		t = "FLORECENCE TRANSFER";
		TR.Warning
			<< "Encountered obsolete experimental technqiue coding '"
			<< t << "'" << endl;
		break;

	case NMR:
		t = "NMR";
		TR.Warning
			<< "Encountered obsolete experimental technqiue coding '"
			<< t << "'" << endl;
		break;

	default:
		TR.Error
			<< "Unrecognized experimental technique value '"
			<< technique << "'" << endl;
		utility_exit();
	}
	return t;
}

HeaderInformation::ExperimentalTechnique
HeaderInformation::string_to_experimental_technique(
	string const & technique
) {
	if(technique == "X-RAY DIFFRACTION")        return X_RAY_DIFFRACTION;
	else if(technique == "FIBER DIFFRACTION")   return FIBER_DIFFRACTION;
	else if(technique == "NEUTRON DIFFRACTION") return NEUTRON_DIFFRACTION;
	else if(technique == "ELECTRON CRYSTALLOGRAPHY")
		return ELECTRON_CRYSTALLOGRAPHY;
	else if(technique == "ELECTRON MICROSCOPY") return ELECTRON_MICROSCOPY;
	else if(technique == "SOLID-STATE NMR")     return SOLID_STATE_NMR;
	else if(technique == "SOLUTION NMR")        return SOLUTION_NMR;
	else if(technique == "SOLUTION SCATTERING") return SOLUTION_SCATTERING;
	else if(technique == "THEORETICAL MODEL")   return THEORETICAL_MODEL;

	// Handle obsolete technique strings
	else if(technique == "ELECTRON DEFRACTION") {
		TR.Warning
			<< "Encountered obsolete experimental technqiue string '"
			<< technique << "'" << endl;
		return ELECTRON_DEFRACTION;
	} else if(technique == "CRYO-ELECTRON MICROSCOPY") {
		TR.Warning
			<< "Encountered obsolete experimental technqiue string '"
			<< technique << "'" << endl;
		return CRYO_ELECTRON_MICROSCOPY;
	} else if(technique == "FLORECENCE TRANSFER") {
		TR.Warning
			<< "Encountered obsolete experimental technqiue string '"
			<< technique << "'" << endl;
		return FLORECENCE_TRANSFER;
	} else if(technique == "NMR") {
		TR.Warning
			<< "Encountered obsolete experimental technqiue string '"
			<< technique << "'" << endl;
		return NMR;
	} else {
		TR.Error
			<< "Unrecognized experimental technique string '"
			<< technique << "'" << endl;
		utility_exit();
	}
	return THEORETICAL_MODEL;
}

void
HeaderInformation::store_experimental_techniques(
	string const & exp) {
	if(exp.empty()){
		TR.Error << "Attempting to add empty experimental technique string." << endl;
		utility_exit();
	}

	size_t t_begin, t_len(0);
	SSize t_end(-1);

	while(true){
		t_begin = exp.find_first_not_of(' ', t_end+1);
		if(t_begin == std::string::npos) return;

		t_end = exp.find(';', t_begin);
		if(t_end == SSize(std::string::npos)){
			experimental_technique_in_progress_ =
				rstripped_whitespace(exp.substr(t_begin, t_len));
			return;
		} else if(exp.length() - t_begin >= 3 && exp.compare(t_begin, 3, "NMR") == 0){
			// The obsolete NMR tag took extra information that is ignored here
			t_len = 3;
		} else {
			t_len = t_end - t_begin;
		}
		if(experimental_technique_in_progress_.empty()){
			experimental_techniques_.push_back(
				string_to_experimental_technique(exp.substr(t_begin, t_len)));
		} else {
			experimental_technique_in_progress_.append(" ");
			experimental_technique_in_progress_.append(exp.substr(t_begin, t_len));
			experimental_techniques_.push_back(
				string_to_experimental_technique(experimental_technique_in_progress_));
			experimental_technique_in_progress_.clear();
		}

	}
	return;
}

void
HeaderInformation::store_experimental_technique(
	HeaderInformation::ExperimentalTechnique technique) {
	experimental_techniques_.push_back(technique);
}

list< HeaderInformation::ExperimentalTechnique > const &
HeaderInformation::experimental_techniques() const {
	return experimental_techniques_;
}

void
HeaderInformation::finalize_experimental_technique_records() {
	if(experimental_technique_in_progress()){
		experimental_techniques_.push_back(
			string_to_experimental_technique(experimental_technique_in_progress_));
		experimental_technique_in_progress_.clear();
	}
}

bool
HeaderInformation::experimental_technique_in_progress() const {
	return !experimental_technique_in_progress_.empty();
}

void
HeaderInformation::clear_experimental_techniques() {
	experimental_techniques_.clear();
}

bool
HeaderInformation::is_experimental_technique(
	HeaderInformation::ExperimentalTechnique technique
) const {
	list< HeaderInformation::ExperimentalTechnique >::const_iterator
		t = find(experimental_techniques_.begin(), experimental_techniques_.end(),
			technique);

	return t != experimental_techniques_.end();
}

void
HeaderInformation::fill_experimental_technique_records(
	std::vector< Record > & VR
) const {
	if(parse_in_progress()){
		TR.Error
			<< "Attempting to fill experimental technique records the "
			<< "HeaderInformation is in the middle of parsing. If you think the "
			<< "parsing is complete and you have reached this recording in error, "
			<< "please call finalize_parse()";
		utility_exit();
	}
	
	if(experimental_techniques_.empty()) return;
	string techniques;
	ExperimentalTechniques::const_iterator
		k = experimental_techniques_.begin(),
		ke= experimental_techniques_.end();
	for(; k != ke; ++k){
		if(!techniques.empty()) techniques.append("; ");
		techniques.append(experimental_technique_to_string(*k));
	}
	Size line_no(1);
	fill_wrapped_records("EXPDTA", "technique", techniques, line_no, VR);
}



////////// Helper Functions /////////////


void
HeaderInformation::fill_wrapped_records(
	string const & record_type,
	string const & field_name,
	string const & contents,
	Size & line_no,
	std::vector< Record > & VR
) const {
	// Assume contents string is stripped of white space
	size_t l_begin(0), l_len(0), l_end(0);
	size_t field_width(60);
	while(l_begin != contents.length()){
		Record R = Field::getRecordCollection()[record_type];
		R["type"].value = record_type;
		set_line_continuation(R, line_no);

		//Will the remainder of the contents fit on this line?
		if(contents.length() - l_begin <= field_width){
			l_len = contents.length() - l_begin;
		} else {
			// Walk back from end where the field would truncate to locate
			// a reasonable place to word wrapping.

			l_end = l_begin + field_width;
			// Note: Since the rest of the contents don't fit in the field,
			// l_end < contents.length()
			while(true){
				if(l_end == l_begin){
					// We have walked all the way to l_begin. The next word is
					// so big it cannot fit in the field
					TR.Error
						<< "The for record type '" << record_type << "', "
						<< "field '" << field_name << "' "
						<< "contains a word that has more than 59 characters and "
						<< "is too long to fit on one line." << endl;
					TR.Error << field_name << ": " << contents << endl;
					utility_exit();
				}
				if (contents[l_end] == ' ' || contents[l_end - 1] == '-'){
					break;
				} else {
					--l_end;
				}
			}
			l_len = l_end - l_begin;
		}
		R[field_name].value = contents.substr(l_begin, l_len);
		VR.push_back(R);
		++line_no;

		// Note this puts l_begin at a ' ' which is how wrapped records
		// are written in the PDB
		l_begin = l_begin + l_len;
	}
}


void
HeaderInformation::set_line_continuation(
	Record & R,
	Size const line_no
) const {
	std::string & con_field = R["continuation"].value;
	if(line_no == 0){
		TR.Error << "Attempting to write a line continuation record for line 0, please begin the line continuation count at 1." << endl;
		utility_exit();
	}
	if(line_no == 1){
		con_field = "  ";
		return;
	} else if(line_no > 99){
		TR.Error << "Attempting to write record that takes more than 99 lines, which overflows the continuation field in the." << endl;
		utility_exit();
	} else {
		con_field.resize(2);
		sprintf(&con_field[0], "%2d", static_cast<int>(line_no));
	}
}

} // namespace pdb
} // namespace io
} // namespace core


