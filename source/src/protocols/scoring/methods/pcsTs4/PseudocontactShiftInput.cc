// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/PseudocontactShiftInput.cc
///
/// @brief Read input .npc input file
///
/// @details The following classes are responsable to read / parse the PCS input file (.npc format)
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references C Schmitz et.al. J Mol Biol. Mar 9, 2012; 416(5): 668–677 ; Yagi H et.al Structure, 2013, 21(6):883-890
///
/// @authorv Christophe Schmitz , Kala Bharath Pilla
///
////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcsTs4/PseudocontactShiftInput.hh>

// Package headers

// Project headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>

// Numeric headers

// Objexx headers

// C++ headers
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace scoring {
namespace methods {
namespace pcsTs4 {

static basic::Tracer TR_PCS_d_i_Ts4( "protocols.scoring.methods.pcsTs4.PCS_data_input_Ts4" );

PCS_line_data_Ts4::PCS_line_data_Ts4(PCS_line_data_Ts4 const & other):
	residue_num_(other.residue_num_),
	atom_name_(other.atom_name_),
	PCS_experimental_(other.PCS_experimental_),
	PCS_tolerance_(other.PCS_tolerance_)
{
}

PCS_line_data_Ts4 &
PCS_line_data_Ts4::operator=( PCS_line_data_Ts4 const & other )
{
	if ( this != &other ) {
		//All data member are const, nothing to copy
	}
	return *this;
}

PCS_line_data_Ts4::~PCS_line_data_Ts4(){
}

PCS_file_data_Ts4::~PCS_file_data_Ts4(){
}

PCS_file_data_Ts4::PCS_file_data_Ts4(PCS_file_data_Ts4 const & other):
	filename_(other.filename_), weight_(other.weight_)
{
	PCS_data_line_all_ = other.PCS_data_line_all_;
}

PCS_file_data_Ts4 &
PCS_file_data_Ts4::operator=( PCS_file_data_Ts4 const & other ){
	if ( this != &other ) {
		PCS_data_line_all_ = other.PCS_data_line_all_;
	}
	return *this;
}


PCS_data_input_Ts4::PCS_data_input_Ts4(){
	utility_exit_with_message( "You shouldn't call the empty constructor for PCS_data_input_Ts4 class" );
}

PCS_data_input_Ts4::~PCS_data_input_Ts4(){
}

PCS_data_input_Ts4::PCS_data_input_Ts4(PCS_data_input_Ts4 const & other){
	PCS_filename_and_data_ = other.PCS_filename_and_data_;
}

PCS_data_input_Ts4 &
PCS_data_input_Ts4::operator=( PCS_data_input_Ts4 const & other ){
	if ( this != &other ) {
		PCS_filename_and_data_ = other.PCS_filename_and_data_;
	}
	return *this;
}

core::Size
PCS_line_data_Ts4::residue_num() const{
	return residue_num_;
}

std::string
PCS_line_data_Ts4::atom_name() const{
	return atom_name_;
}

core::Real
PCS_line_data_Ts4::PCS_experimental() const{
	return PCS_experimental_;
}


core::Real
PCS_line_data_Ts4::PCS_tolerance() const{
	return PCS_tolerance_;
}

PCS_line_data_Ts4::PCS_line_data_Ts4(core::Size residue_num,
	std::string atom_name,
	core::Real PCS_experimental,
	core::Real PCS_tolerance
) :
	residue_num_( residue_num ),
	atom_name_(atom_name),
	PCS_experimental_(PCS_experimental),
	PCS_tolerance_(PCS_tolerance)
{
}

std::map< std::string, PCS_file_data_Ts4 > &
PCS_data_input_Ts4::get_PCS_data_input_reference(){
	return  PCS_filename_and_data_;
}

std::string
PCS_file_data_Ts4::get_filename() const{
	return(filename_);
}

core::Real
PCS_file_data_Ts4::get_weight() const{
	return(weight_);
}

PCS_file_data_Ts4::PCS_file_data_Ts4(std::string const & filename, core::Real const my_weight ):
	filename_(std::string(filename)), weight_(my_weight)
{
	read_PCS_file();
}

PCS_data_input_Ts4::PCS_data_input_Ts4(utility::vector1<std::string> const & filenames, utility::vector1<core::Real> const & weight){
	// utility::vector1<std::string>::const_iterator it;
	core::Real weight_sum;
	core::Size i;

	weight_sum = 0;
	for ( i = 1; i <= filenames.size(); i++ ) {
		weight_sum += weight[i];
	}

	for ( i = 1; i <= filenames.size(); i++ ) {
		//core::Real my_weight(weight[i]/weight_sum);
		//TODO correct the weighting scheme. For the moment it is one automatically
		core::Real my_weight(weight[i]);
		PCS_file_data_Ts4 pcs_f_d_temp(filenames[i], my_weight);
		PCS_filename_and_data_.insert ( std::pair< std::string, PCS_file_data_Ts4 >(filenames[i], pcs_f_d_temp) );
	}
}

void
PCS_file_data_Ts4::read_PCS_file(){
	core::Size residue_num;
	std::string atom_name;
	core::Real PCS_experimental;
	core::Real PCS_tolerance;
	std::ifstream myfile;
	std::string line;
	core::Size line_number(0);

	TR_PCS_d_i_Ts4 << "Opening file '" << get_filename().c_str() << "'" << std::endl;
	myfile.open (get_filename().c_str(), std::ios::in);
	if ( !myfile.is_open () ) {
		std::cerr << "Unable to open the file '" << get_filename().c_str()  <<"'" << std::endl;
		utility_exit();
	}

	while ( getline( myfile, line ) ) {
		line_number++;
		std::istringstream line_stream( line ,std::istringstream::in);
		if ( (line_stream >> residue_num >> atom_name >> PCS_experimental >> PCS_tolerance).fail() ) {
			TR_PCS_d_i_Ts4 << "Ignoring line " <<line_number << ": `" << line << "` from file " << get_filename() <<std::endl;
			continue;
		}

		PCS_data_line_all_.push_back( PCS_line_data_Ts4( residue_num, atom_name, PCS_experimental, PCS_tolerance ) );
	}

	myfile.close();
}

utility::vector1<PCS_line_data_Ts4> &
PCS_file_data_Ts4::get_PCS_data_line_all_reference(){
	return (PCS_data_line_all_);
}

std::ostream &
operator<<(std::ostream& out, const PCS_line_data_Ts4 &PCS_l_d){
	out << "Residue: " << std::setw(4) << PCS_l_d.residue_num();
	out << "   Atom: " << std::setw(4) << PCS_l_d.atom_name();
	out << "   PCS: " << std::setw(7) << PCS_l_d.PCS_experimental();
	out << "   Tolerance: " << std::setw(7) << PCS_l_d.PCS_tolerance()<< std::endl;
	return out;
}

std::ostream &
operator<<(std::ostream& out, const PCS_file_data_Ts4 &PCS_f_d){
	utility::vector1<PCS_line_data_Ts4>::iterator it;
	utility::vector1<PCS_line_data_Ts4> PCS_d_l_a;
	PCS_d_l_a = PCS_f_d.PCS_data_line_all_;

	for ( it = PCS_d_l_a.begin(); it != PCS_d_l_a.end(); ++it ) {
		out << *it;
	}
	out<< PCS_d_l_a.size() << " PCS in total for this file" << std::endl;
	return out;
}

std::ostream &
operator<<(std::ostream & out,  const PCS_data_input_Ts4 &PCS_d_i ){

	std::map< std::string, PCS_file_data_Ts4 >::iterator it;
	std::map< std::string, PCS_file_data_Ts4 > mymap;
	mymap = PCS_d_i.PCS_filename_and_data_;

	for ( it = mymap.begin(); it != mymap.end(); ++it ) {
		out << "For the file '" << it->first <<"' the PCS are:" << std::endl;
		out << it->second;// << std::endl;
		out << "The relative weight is " << (it->second).get_weight() << std::endl;
	}
	return out;
}


PCS_data_input_manager_Ts4::PCS_data_input_manager_Ts4(){

}

PCS_data_input_manager_Ts4 *
PCS_data_input_manager_Ts4::get_instance(){
	if ( instance_ == 0 ) {
		instance_ = new PCS_data_input_manager_Ts4();
	}
	return instance_;
}

PCS_data_input_Ts4
PCS_data_input_manager_Ts4::get_input_data(utility::vector1<std::string> const & filenames, utility::vector1<core::Real> const & weight){
	std::string id;
	core::Size i;

	for ( i = 1; i <= filenames.size(); ++i ) {
		id += filenames[i];
	}

	std::map< std::string, PCS_data_input_Ts4 >::iterator it;

	for ( it = file_2_data_map_.begin(); it != file_2_data_map_.end(); ++it ) {
		if ( it->first == id ) {
			return(it->second);
		}
	}

	PCS_data_input_Ts4 pcs_d_i(filenames, weight);
	it = file_2_data_map_.begin();
	file_2_data_map_.insert(it, std::pair< std::string , PCS_data_input_Ts4 >( id ,pcs_d_i));

	TR_PCS_d_i_Ts4 << pcs_d_i << std::endl;

	return(pcs_d_i);
}

PCS_data_input_manager_Ts4 * PCS_data_input_manager_Ts4::instance_( 0 );

}//namespace pcsTs4
}//namespace methods
}//namespace scoring
}//namespace protocols

#ifdef    SERIALIZATION
/// @brief Bogus default constructor to be used only for deserialization
protocols::scoring::methods::pcsTs4::PCS_line_data_Ts4::PCS_line_data_Ts4() :
	residue_num_( 0 ), atom_name_( "temp" ), PCS_experimental_( 1234 ), PCS_tolerance_( 1234 )
{}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::scoring::methods::pcsTs4::PCS_line_data_Ts4::save( Archive & arc ) const {
	arc( CEREAL_NVP( residue_num_ ) ); // const core::Size
	arc( CEREAL_NVP( atom_name_ ) ); // const std::string
	arc( CEREAL_NVP( PCS_experimental_ ) ); // const core::Real
	arc( CEREAL_NVP( PCS_tolerance_ ) ); // const core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::scoring::methods::pcsTs4::PCS_line_data_Ts4::load( Archive & arc ) {
	arc( const_cast< core::Size & > (residue_num_) ); // const core::Size
	arc( const_cast< std::string & > (atom_name_) ); // const std::string
	arc( const_cast< core::Real & > (PCS_experimental_) ); // const core::Real
	arc( const_cast< core::Real & > (PCS_tolerance_) ); // const core::Real
}
SAVE_AND_LOAD_SERIALIZABLE( protocols::scoring::methods::pcsTs4::PCS_line_data_Ts4 );
#endif // SERIALIZATION
