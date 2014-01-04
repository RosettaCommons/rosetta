// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/chemical_shift/RNA_ChemicalShiftPotential.cc
/// @brief  The real workhorse behind RNA_ChemicalShiftEnergy.cc
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)


// Unit headers
#include <core/scoring/rna/chemical_shift/RNA_ChemicalShiftPotential.hh>


// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>


// Utility headers

//Auto Headers
#include <core/id/AtomID.hh>

////////////////////////////////////////////////////////
#include <basic/Tracer.hh>
////////////////////////////////////////////////////////
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh> 
//#include <core/io/database/open.hh>
#include <utility/file/file_sys_util.hh> 
#include <utility/io/izstream.hh>
#include <ObjexxFCL/format.hh>
#include <core/scoring/rna/chemical_shift/RNA_CS_Util.hh>
#include <core/scoring/rna/chemical_shift/RNA_CS_RingCurrent.hh>
#include <core/scoring/rna/chemical_shift/RNA_CS_MagneticAnisotropy.hh>
#include <math.h>
#include <numeric/xyzVector.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>



// C++

static basic::Tracer TR( "core.scoring.rna.chemical_shift.RNA_ChemicalShiftPotential" );

namespace core {
namespace scoring {
namespace rna {
namespace chemical_shift {



	/// c-tor
	RNA_ChemicalShiftPotential::RNA_ChemicalShiftPotential():
		rna_cs_params_( RNA_CS_parameters() ),
		verbose_( true ),
		include_ring_current_effect_( true ),
		include_magnetic_anisotropy_effect_( true ),
		total_exp_chemical_shift_data_points_( 0 )
	{

		if ( include_ring_current_effect_ ){
			std::cout << "include_ring_current_effect_ = true " << std::endl;
		} else{
			std::cout << "include_ring_current_effect_ = false " << std::endl;
		}

		if ( include_magnetic_anisotropy_effect_ ){
			std::cout << "include_magnetic_anisotropy_effect_ = true " << std::endl;
		} else{
			std::cout << "include_magnetic_anisotropy_effect_ = false " << std::endl;
		}

		if ( basic::options::option[ basic::options::OptionKeys::score::rna_chemical_shift_exp_data].user() == false ){
			utility_exit_with_message( "User need to pass in score:rna_chemical_shift_exp_data" );
		}

		std::string const exp_CS_data_filename = basic::options::option[ basic::options::OptionKeys::score::rna_chemical_shift_exp_data]();

		utility::vector1< core::Size > include_res_list;

		if ( basic::options::option[ basic::options::OptionKeys::score::rna_chemical_shift_include_res].user() ){
			include_res_list = basic::options::option[ basic::options::OptionKeys::score::rna_chemical_shift_include_res]();
		} else{
			TR << "User did not pass in score:rna_chemical_shift_include_res, including all residue!" << std::endl;
			if ( include_res_list.size() != 0 ) utility_exit_with_message( "User did not pass in score:rna_chemical_shift_include_res option but include_res_list.size() != 0!" );
		}

		if ( basic::options::option[ basic::options::OptionKeys::score::rna_chemical_shift_H5_prime_mode].user() ){
			H5_prime_mode_ = basic::options::option[ basic::options::OptionKeys::score::rna_chemical_shift_H5_prime_mode]();
			TR << "Using user - specified H5_prime_mode_ = " << H5_prime_mode_ << std::endl;
		} else{
			H5_prime_mode_ = "LEAST_SQUARE_IGNORE_DUPLICATE"; //DEFAULT!
			TR << "Using default H5_prime_mode_ = " << H5_prime_mode_ << std::endl;
		}


		utility::vector1< utility::vector1< std::string > > proton_entry_list;

		proton_entry_list.push_back( string_list( "H1'" ) );
		proton_entry_list.push_back( string_list( " H2'" ) );
		proton_entry_list.push_back( string_list( "H3'" ) );
		proton_entry_list.push_back( string_list( "H4'" ) );

		if ( H5_prime_mode_ == "UNIQUE" ){
			proton_entry_list.push_back( string_list( " H5'" ) );
			proton_entry_list.push_back( string_list( "H5''" ) );

		} else if ( H5_prime_mode_ == "LEAST_SQUARE" || H5_prime_mode_ == "LEAST_SQUARE_IGNORE_DUPLICATE" ){
			proton_entry_list.push_back( string_list( " H5'", "H5''" ) ); 

		} else{
			utility_exit_with_message( "Invalid H5_prime_mode_ ( " + H5_prime_mode_ + " )!" );
		}

		/////Non-polar base protons/////
		proton_entry_list.push_back( string_list( "H2" ) );
		proton_entry_list.push_back( string_list( "H5" ) );
		proton_entry_list.push_back( string_list( "H6" ) );
		proton_entry_list.push_back( string_list( "H8" ) );


		import_exp_chemical_shift_data( exp_CS_data_filename, include_res_list, proton_entry_list );

	}



	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////
	chemical::AA
	get_res_aa_from_BASE_name( std::string BASE_name, std::string const text_line ){

		chemical::AA res_aa = chemical::aa_unk;

		if ( BASE_name == "G" ){

			res_aa = chemical::na_rgu;

		} else if ( BASE_name == "A" ){

			res_aa = chemical::na_rad;

		} else if ( BASE_name == "C" ){

			res_aa = chemical::na_rcy;

		} else if ( BASE_name == "U" ){

			res_aa = chemical::na_ura;

		} else{
			utility_exit_with_message( "Invalid BASE_name ( " + BASE_name + " ) | line = ( " + text_line + " )" );
		}
		
		return res_aa;

	}

	/////////////////////////////////////////////////////////////////////////////
	std::string
	remove_whitespaces( std::string const in_atom_name ){

		std::string out_atom_name = "";

		for ( Size n = 0; n < in_atom_name.size(); n++ ){

			if ( in_atom_name[n] != ' ' ) out_atom_name += in_atom_name[n];

		}

		return out_atom_name;
	}


	/////////////////////////////////////////////////////////////////////////////
	bool
	is_polar_hydrogen( std::string const input_atom_name ){

		if ( input_atom_name == "HO2'" || input_atom_name == "HO2'" ) return true;

		if ( input_atom_name == "H22"  || input_atom_name == "H21" ) return true;

		if ( input_atom_name == "H21"  || input_atom_name == "H22" ) return true;

		if ( input_atom_name == "H41"  || input_atom_name == "H41" ) return true;

		if ( input_atom_name == "H42"  || input_atom_name == "H42" ) return true;

		if ( input_atom_name == "H61"  || input_atom_name == "H61" ) return true;

		if ( input_atom_name == "H62"  || input_atom_name == "H62" ) return true;

		if ( input_atom_name == "H1" ) return true;

		if ( input_atom_name == "H3" ) return true;

		if ( input_atom_name == "HO3'" ) return true; //Can occur at 3' ends of RNA's chain. This is not part of the standard rosetta atom-set

		if ( input_atom_name == "HO5'" ) return true; //Can occur at 5' ends of RNA's chain. this is not part of the standard rosetta atom-set

		return false;

	}

	/////////////////////////////////////////////////////////////////////////////
	std::string
	get_rosetta_hatom_name( std::string const input_atom_name, std::string const text_line, utility::vector1< std::string > const & flat_proton_entry_list ){ 

		using namespace ObjexxFCL;

		std::string rosetta_atom_name = "";

		//Assume that input_atom_name is a non_polar hydrogen atom! Other atoms should be filtered before reaching this point!
		if ( input_atom_name == "H1'" ){
			rosetta_atom_name = "H1'";

		} else if ( input_atom_name == "H2'" ){
			rosetta_atom_name = " H2'";

		} else if ( input_atom_name == "H3'" ){
			rosetta_atom_name = "H3'";

		} else if ( input_atom_name == "H4'" ){
			rosetta_atom_name = "H4'";

		} else if ( input_atom_name == "H5'" ){
			rosetta_atom_name = " H5'";

		} else if ( input_atom_name == "H5''" ){
			rosetta_atom_name = "H5''";
		} else{
			rosetta_atom_name = input_atom_name;
		}

		Size num_matching_atom_name = 0;

		for ( Size ii = 1; ii <= flat_proton_entry_list.size(); ii++ ){
			if ( rosetta_atom_name == flat_proton_entry_list[ii] ) num_matching_atom_name++;
		}

		if ( num_matching_atom_name != 1 ){
			std::cout << "ERROR: num_matching_atom_name = " << num_matching_atom_name << std::endl;
			std::cout << "ERROR: input_atom_name = " << input_atom_name << std::endl;
			std::cout << "ERROR: rosetta_atom_name = " << rosetta_atom_name << std::endl;
			utility_exit_with_message( "num_matching_atom_name != 1 for input_atom_name ( " + input_atom_name + " ) | text_line ( " + text_line + " )" );
		}

		return rosetta_atom_name;

	}


	/////////////////////////////////////////////////////////////////////////////
	void
	print_chemical_shift_data( std::string prestring, ChemicalShiftData const & CS_data, bool const print_data_line ){
	
		std::cout << prestring  << "seq_num = " << std::setw( 3 ) << CS_data.seq_num;
		std::cout << " | res_aa = " << std::setw( 3 ) << name_from_aa( CS_data.res_aa );
		std::cout << " | atom_name = " << std::setw( 5 ) << CS_data.atom_name;
		std::cout << " | realatomdata_index = " << std::setw( 3 ) << CS_data.realatomdata_index;
		std::cout << " | exp_shift = " << std::setw( 8 ) << CS_data.exp_shift;
		if ( print_data_line ) std::cout << " | data_line = " << CS_data.data_line;
		std::cout << std::endl;

	}

	////////////////////////////////copied from protocols/swa/rna/StepWiseRNA_Util.cc/////////////////////////////
	core::Size
	string_to_int( std::string const input_string ){

		Size int_of_string; //misnomer
		std::stringstream ss ( std::stringstream::in | std::stringstream::out );

		ss << input_string;

		if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss << input_string | string ( " + input_string + " )" );

		ss >> int_of_string;

		if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss >> int_of_string | string ( " + input_string + " )" );

		return int_of_string;
	}

	////////////////////////////////copied from protocols/swa/rna/StepWiseRNA_Util.cc/////////////////////////////
	core::Real
	string_to_real( std::string const input_string ){

		Real real_of_string;
		std::stringstream ss ( std::stringstream::in | std::stringstream::out );

		ss << input_string;

		if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss << input_string | string ( " + input_string + " )" );

		ss >> real_of_string;

		if ( ss.fail() ) utility_exit_with_message( "In string_to_real(): ss.fail() for ss >> real_of_string | string ( " + input_string + " )" );

		return real_of_string;

	}

	////////////////////////////////copied from protocols/swa/rna/StepWiseRNA_Util.cc/////////////////////////////
	bool
	Contain_seq_num( Size const & seq_num, utility::vector1< core::Size > const & residue_list ){
		for ( Size j = 1; j <= residue_list.size(); j++ ){
			if ( seq_num == residue_list[j] ) {
				return true;
			}
		}
		return false;
	}

	/////////////////////////////////////////////////////////////////////////////

	utility::vector1 < ChemicalShiftData >
	filter_chem_shift_data_list( utility::vector1 < ChemicalShiftData > const & flat_EXP_chem_shift_data_list, core::Size const seq_num, utility::vector1 < std::string > const & proton_entry ){

		utility::vector1 < ChemicalShiftData > filtered_CS_data_list;
	
		for ( Size ii = 1; ii <= flat_EXP_chem_shift_data_list.size(); ii++ ){

			ChemicalShiftData const & CS_data = flat_EXP_chem_shift_data_list[ii];

			if ( CS_data.seq_num != seq_num ) continue;

			Size num_matching_atom_name = 0;

			for ( Size jj = 1; jj <= proton_entry.size(); jj++ ){	
				if ( proton_entry[jj] == CS_data.atom_name ) num_matching_atom_name++;
			}

			if ( num_matching_atom_name > 1 ) utility_exit_with_message( "num_matching_atom_name > 1!" );

			if ( num_matching_atom_name == 1 ) filtered_CS_data_list.push_back( CS_data );

		}

		if ( filtered_CS_data_list.size() > 0 ){

			if ( filtered_CS_data_list.size() != proton_entry.size() ){
				std::cout << "ERROR: filtered_CS_data_list.size() > 0 at seq_num = " << seq_num << " proton_entry: ";
				for ( Size jj = 1; jj <= proton_entry.size(); jj++ ){
					std::cout << "# " << jj << " :" << proton_entry[jj] <<std::endl;
				}
				std::cout << "ERROR: filtered_CS_data_list: " << std::endl;
				for ( Size jj = 1; jj <= filtered_CS_data_list.size(); jj++ ){
					print_chemical_shift_data( "#" + ObjexxFCL::lead_zero_string_of(jj, 3) + " :", filtered_CS_data_list[jj], true );
				}

				utility_exit_with_message( "filtered_CS_data_list.size() > 0 BUT filtered_CS_data_list.size() != proton_entry.size()!" );
			}
		}

		return filtered_CS_data_list;

	}

	/////////////////////////////////////////////////////////////////////////////
	Size
	RNA_ChemicalShiftPotential::get_total_exp_chemical_shift_data_points() const{

		if ( total_exp_chemical_shift_data_points_ == 0 ) utility_exit_with_message( "total_exp_chemical_shift_data_points_ == 0!" );

		return total_exp_chemical_shift_data_points_;

	}

	/////////////////////////////////////////////////////////////////////////////

	utility::vector1< std::string > 
	RNA_ChemicalShiftPotential::string_list( std::string const string_one ) const{
	
		utility::vector1< std::string > string_list;
	
		string_list.push_back( string_one );

		return string_list;
	}

	/////////////////////////////////////////////////////////////////////////////

	utility::vector1< std::string > 
	RNA_ChemicalShiftPotential::string_list( std::string const string_one, const std::string string_two ) const{
	
		utility::vector1< std::string > string_list;
	
		string_list.push_back( string_one );
		string_list.push_back( string_two );

		return string_list;
	}


	/////////////////////////////////////////////////////////////////////////////
	Size
	RNA_ChemicalShiftPotential::get_realatomdata_index( std::string const & in_atom_name, chemical::AA const res_aa ) const{

		using namespace ObjexxFCL;

		std::string const atom_name = remove_whitespaces( in_atom_name );

		RNA_CS_residue_parameters const & rna_cs_rsd_params = rna_cs_params_.get_RNA_CS_residue_parameters( res_aa );

		Size const maxatoms = rna_cs_rsd_params.get_atomnames_size();

		Size num_matching_atom_name = 0;

		Size realatomdata_index = 0;

		for ( Size count = 1; count <= maxatoms; count++ ){

			if ( rna_cs_rsd_params.get_atomname( count ) == atom_name ){
				num_matching_atom_name++;
				realatomdata_index = count;
			}

		}

		if ( num_matching_atom_name != 1 ) utility_exit_with_message( "num_matching_atom_name = ( " + string_of( num_matching_atom_name ) + " ) != 1 | atom_name ( " + atom_name + " )!" );

		return realatomdata_index;

	}


	/////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChemicalShiftPotential::assert_is_calc_chem_shift_atom( ChemicalShiftData const & CS_data ) const{

		using namespace ObjexxFCL;

		//OK this is another layer of consistency check////
		if ( CS_data.atom_name == "H8" ){
			if ( CS_data.res_aa != chemical::na_rgu && CS_data.res_aa != chemical::na_rad ){
				utility_exit_with_message( "CS_data.atom_name == \"H8\" but CS_data.res_aa != na_rgu && CS_data.res_aa != na_rad!" );
			}
		}

		if ( CS_data.atom_name == "H5" || CS_data.atom_name == "H6" ){
			if ( CS_data.res_aa != chemical::na_rcy && CS_data.res_aa != chemical::na_ura ){
				utility_exit_with_message( "CS_data.atom_name == \"" + CS_data.atom_name + "\" but CS_data.res_aa != na_rcy && CS_data.res_aa != na_ura!" );
			}
		}

		if ( CS_data.atom_name == "H2" ){
			if ( CS_data.res_aa != chemical::na_rad ) utility_exit_with_message( "CS_data.atom_name == \"H2\" but CS_data.res_aa != na_rad!" );
		}


		////////////////////////////This is imino/amino (polar) proton [not current in use!]////////////////////////////////////////
		if ( CS_data.atom_name == "H1" ){ //Still need to account for the possibility of protonated Adenosine at N1 pos!
			if ( CS_data.res_aa != chemical::na_rgu ) utility_exit_with_message( "CS_data.atom_name == \"H1\" but CS_data.res_aa != na_rgu!" );
		}

		if ( CS_data.atom_name == "H3" ){ //Still need to account for the possibility of protonated Cytosine at N3 pos!
			if ( CS_data.res_aa != chemical::na_ura ) utility_exit_with_message( "CS_data.atom_name == \"H3\" but CS_data.res_aa != na_ura!" );
		}

		if ( CS_data.atom_name == "H22" || CS_data.atom_name == "H21" ){
			if ( CS_data.res_aa != chemical::na_rgu ){
				utility_exit_with_message( "CS_data.atom_name == \"" + CS_data.atom_name + "\" but CS_data.res_aa != na_rgu!" );
			}
		}

		if ( CS_data.atom_name == "H41" || CS_data.atom_name == "H42" ){
			if ( CS_data.res_aa != chemical::na_rcy ){
				utility_exit_with_message( "CS_data.atom_name == \"" + CS_data.atom_name + "\" but CS_data.res_aa != na_rcy!" );
			}
		}

		if ( CS_data.atom_name == "H61" || CS_data.atom_name == "H62" ){
			if ( CS_data.res_aa != chemical::na_rad ){
				utility_exit_with_message( "CS_data.atom_name == \"" + CS_data.atom_name + "\" but CS_data.res_aa != na_rad!" );
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		RNA_CS_residue_parameters const & rna_cs_rsd_params = rna_cs_params_.get_RNA_CS_residue_parameters( CS_data.res_aa );	
	
		if ( rna_cs_rsd_params.atom_data( CS_data.realatomdata_index, csca ) == false ){
			 print_chemical_shift_data( "ERROR CS_data:", CS_data, true );
			 utility_exit_with_message( "CS_data.realatomdata_index ( " + string_of( CS_data.realatomdata_index ) + " ) | CS_data.atom_name ( " + CS_data.atom_name + " ) is not a calc_chem_shift_atom!" );
		}

	}



	/////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ChemicalShiftPotential::Is_magnetic_anisotropy_source_atom( core::conformation::Residue const & rsd, Size const atomno ) const{

		using namespace ObjexxFCL;

		std::string const atom_name = remove_whitespaces( rsd.atom_name( atomno ) );

		RNA_CS_residue_parameters const & rna_cs_rsd_params = rna_cs_params_.get_RNA_CS_residue_parameters( rsd.aa() );

		Size const realatomdata_index = get_realatomdata_index( atom_name, rsd.aa() );		

		bool const Is_MA_source_atom = ( dround( rna_cs_rsd_params.atom_data( realatomdata_index, maca ) ) == 1 );

		return Is_MA_source_atom;

	}

	/////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ChemicalShiftPotential::atom_has_exp_chemical_shift_data( core::conformation::Residue const & rsd, Size const atomno ) const{

		using namespace ObjexxFCL;

		Size num_matching_CS_data = 0;

		Size const seq_num = rsd.seqpos();

		std::string const atom_name = remove_whitespaces( rsd.atom_name( atomno ) );

		for ( Size outer_data_ID = 1; outer_data_ID <= EXP_chem_shift_data_list_.size(); outer_data_ID++ ){

			for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_list_[outer_data_ID].size(); inner_data_ID++ ){

				ChemicalShiftData const & CS_data = EXP_chem_shift_data_list_[outer_data_ID][inner_data_ID];

				if ( CS_data.seq_num != seq_num ) continue;

				if ( CS_data.atom_name != atom_name ) continue;

				num_matching_CS_data++;

			}
		}

		if ( num_matching_CS_data > 1 ) utility_exit_with_message( "num_matching_CS_data > 1 for seq_num = ( " + string_of( seq_num ) + " ) | atom_name = ( " + atom_name + " )" );

		return ( num_matching_CS_data == 1 );

	}


	/////////////////////////////////////////////////////////////////////////////
	utility::vector1 < ChemicalShiftData > const & 
	RNA_ChemicalShiftPotential::get_matching_CS_data_entry( Size const seq_num, std::string const in_atom_name ) const {

		using namespace ObjexxFCL;

		Size num_matching_CS_data = 0;

		Size matching_outer_data_ID = 0;

		std::string const atom_name = remove_whitespaces( in_atom_name );

		for ( Size outer_data_ID = 1; outer_data_ID <= EXP_chem_shift_data_list_.size(); outer_data_ID++ ){

			for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_list_[outer_data_ID].size(); inner_data_ID++ ){

				ChemicalShiftData const & CS_data = EXP_chem_shift_data_list_[outer_data_ID][inner_data_ID];

				if ( CS_data.seq_num != seq_num ) continue;

				if ( CS_data.atom_name != atom_name ) continue;

				num_matching_CS_data++;

				matching_outer_data_ID = outer_data_ID;

			}

		}

		if ( num_matching_CS_data != 1 ) utility_exit_with_message( "num_matching_CS_data != 1 for seq_num = ( " + string_of( seq_num ) + " ) | atom_name = ( " + atom_name + " )" );

		return EXP_chem_shift_data_list_[matching_outer_data_ID];

	}

	/////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChemicalShiftPotential::import_exp_chemical_shift_data( std::string exp_CS_data_filename, 
													utility::vector1 < Size > include_res_list,
													utility::vector1< utility::vector1< std::string > > const & proton_entry_list ){

		using namespace ObjexxFCL;

		//////////////////////////////////////////////////////////////////////

		if ( utility::file::file_exists( exp_CS_data_filename ) == false ){
			utility_exit_with_message( "exp_CS_data_filename " + exp_CS_data_filename + " doesn't exist!" );
		}

		utility::io::izstream data( exp_CS_data_filename );
		if ( !data.good() ){
			utility_exit_with_message( "Unable to open exp_CS_data_file: " + exp_CS_data_filename );
		} else{
			std::cout << "Successfully opened open exp_CS_data_file: \"" << exp_CS_data_filename << "\" !" << std::endl;
		}
		//////////////////////////////////////////////////////////////////////
		utility::vector1 < ChemicalShiftData > flat_EXP_chem_shift_data_list;

		utility::vector1< std::string > flat_proton_entry_list;	
		for ( Size ii = 1; ii <= proton_entry_list.size(); ii++ ){
			for ( Size jj = 1; jj <= proton_entry_list[ii].size(); jj++ ){
				flat_proton_entry_list.push_back( proton_entry_list[ii][jj] );
			}
		}

		//////////////////////////////////////////////////////////////////////

		std::string text_line;

		Size max_seq_num = 0;

		while ( getline( data, text_line ) ){

			std::istringstream text_stream( text_line );

			utility::vector1 < std::string > text_line_list;

			while ( true ){
			
				std::string str_element;

				text_stream >> str_element;

				if ( text_stream.fail() ) break;
		
				text_line_list.push_back( str_element );
			}

			if ( text_line_list.size() != 9 ){
				utility_exit_with_message( "text_line_list.size() != 9 for line ( " + text_line + " )" ); 
			}

			//Consistency_check
			if ( text_line_list[2] != text_line_list[3] ){
				utility_exit_with_message( "text_line_list[2] = ( " + text_line_list[2] + " ) != ( " + text_line_list[3] + " ) = text_line_list[3]" );
			}

			Size const seq_num = string_to_int( text_line_list[3] );


			chemical::AA res_aa = get_res_aa_from_BASE_name( text_line_list[4], text_line );

			// different member versions of find in the same order as above:
			std::string const input_atom_name = text_line_list[5];

			size_t found = input_atom_name.find( "H" );

			if ( found == std::string::npos ){ //Filter out potential carbon, nitrogen and phosphorus chemical_shift line
				if ( false ) std::cout << "Ignoring non - hydrogen   chemical shift line ( " << text_line << " )" << std::endl;
				continue;
			} 

			if ( is_polar_hydrogen( input_atom_name ) ){
				if ( false ) std::cout << "Ignoring polar hydrogen chemical shift line ( " << text_line << " )" << std::endl;
				continue;
			}

			//Assume that input_atom_name is a non_polar hydrogen atom! Other atoms should be filtered before reaching this point!
			std::string const atom_name = get_rosetta_hatom_name( input_atom_name, text_line, flat_proton_entry_list ); 

			Size const realatomdata_index = get_realatomdata_index( atom_name, res_aa );

			Real const exp_shift = string_to_real( text_line_list[7] );

			ChemicalShiftData const chem_shift_data( seq_num, res_aa, atom_name, realatomdata_index, exp_shift, text_line );

			assert_is_calc_chem_shift_atom( chem_shift_data );

			flat_EXP_chem_shift_data_list.push_back( chem_shift_data );

			if ( max_seq_num < seq_num ) max_seq_num = seq_num;

		}

		//////////////////////////////////////////////////////////////////////
		if ( include_res_list.size() == 0 ){
			TR << "User did not pass in score:rna_chemical_shift_include_res, including all residue:";		
			for ( Size seq_num = 1; seq_num <= max_seq_num; seq_num++ ){
				include_res_list.push_back( seq_num );
				TR <<  " " << seq_num;
			} 
			TR <<  std::endl;
		} else{
			TR << "User pass in score:rna_chemical_shift_include_res:";
			for ( Size res_ID = 1; res_ID <= include_res_list.size(); res_ID++ ){
				TR <<  " " << include_res_list[res_ID];			
			}
			TR << std::endl;
		}

		//////////////////////////////////////////////////////////////////////
		EXP_chem_shift_data_list_.clear();

		for ( Size proton_entry_ID = 1; proton_entry_ID <= proton_entry_list.size(); proton_entry_ID++ ){

			for ( Size res_ID = 1; res_ID <= include_res_list.size(); res_ID++ ){
	
				Size const seq_num = include_res_list[res_ID];
				utility::vector1< std::string > const proton_entry = proton_entry_list[proton_entry_ID];

				utility::vector1 < ChemicalShiftData > const & filterd_CS_data_list = filter_chem_shift_data_list( flat_EXP_chem_shift_data_list, seq_num, proton_entry );

				if ( filterd_CS_data_list.size() != 0 ){
					EXP_chem_shift_data_list_.push_back( filterd_CS_data_list );
				}
			}

		}

		//////////////////////////////////////////////////////////////////////
		total_exp_chemical_shift_data_points_ = 0;

		for ( Size outer_data_ID = 1; outer_data_ID <= EXP_chem_shift_data_list_.size(); outer_data_ID++ ){

			utility::vector1 < ChemicalShiftData > const & EXP_chem_shift_data_entry = EXP_chem_shift_data_list_[outer_data_ID];

			utility::vector1 < Real > mock_calc_chem_shift_entry;
			for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_entry.size(); inner_data_ID++ ){
				mock_calc_chem_shift_entry.push_back( 0.0 );
			}

			utility::vector1 < Real > actual_exp_chem_shift_entry;
			utility::vector1 < bool > do_include_CS_data;		 

			get_best_exp_to_calc_chem_shift_mapping( EXP_chem_shift_data_entry, mock_calc_chem_shift_entry, actual_exp_chem_shift_entry, do_include_CS_data );


			for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_entry.size(); inner_data_ID++ ){

				if ( do_include_CS_data[inner_data_ID] == false ) continue;

				total_exp_chemical_shift_data_points_++;

			}
		}
		std::cout << "total_exp_chemical_shift_data_points_ ( possibly ignoring duplicates ) = " << total_exp_chemical_shift_data_points_ << std::endl;

		//////////////////////////////////////////////////////////////////////
		Size data_count = 0;

		if ( verbose_ ){
		
			std::cout << "------------------Imported exp_chem_shift_data_list------------------" << std::endl;
			for ( Size ii = 1; ii <= EXP_chem_shift_data_list_.size(); ii++ ){
				for ( Size jj = 1; jj <= EXP_chem_shift_data_list_[ii].size(); jj++ ){
					data_count++;

					std::string const prefix = "#" + lead_zero_string_of(ii, 3) + "." + lead_zero_string_of(jj, 1) + " | count= " + lead_zero_string_of(data_count, 3) + " :";

					print_chemical_shift_data( prefix, EXP_chem_shift_data_list_[ii][jj], true );
				}
			}
			std::cout << "---------------------------------------------------------------------" << std::endl;
		}
		//////////////////////////////////////////////////////////////////////

	}

	//////////////////////////////////////////////////////////////////////
	Real
	RNA_ChemicalShiftPotential::get_calc_chem_shift_value( ChemicalShiftData const & CS_data, pose::Pose const & pose ) const{

		if ( ( CS_data.seq_num < 1 ) || ( CS_data.seq_num > pose.total_residue() ) ){
			std::cout << "ERROR: CS_data.seq_num = " << CS_data.seq_num << std::endl;
			std::cout << "ERROR: pose.total_residue() = " << pose.total_residue() << std::endl;
			utility_exit_with_message( "( CS_data.seq_num < 1 ) || ( CS_data.seq_num > pose.total_residue() )" );
		}

		if ( CS_data.res_aa != pose.residue( CS_data.seq_num ).aa() ){
			std::cout << "ERROR: CS_data.res_aa = " << CS_data.res_aa << std::endl;
			std::cout << "ERROR: CS_data.seq_num = " << CS_data.seq_num << std::endl;
			std::cout << "ERROR: pose.residue( CS_data.seq_num ).aa() = " << pose.residue( CS_data.seq_num ).aa() << std::endl;
			utility_exit_with_message( "chem_shift_data.res_aa != pose1.residue( CS_data.seq_num ).aa()" );
		}

		Real calc_chem_shift = 0.0;

		core::conformation::Residue const & curr_rsd = pose.residue( CS_data.seq_num );

		RNA_CS_residue_parameters const & rna_cs_curr_rsd_params = rna_cs_params_.get_RNA_CS_residue_parameters( curr_rsd.aa() );

		Size const curr_atom_index = curr_rsd.atom_index( CS_data.atom_name );

		numeric::xyzVector< core::Real > const & curr_atom_xyz = curr_rsd.xyz( curr_atom_index );

		bool curr_atom_is_sugar = ( dround( rna_cs_curr_rsd_params.atom_data( CS_data.realatomdata_index, suga ) ) == 1 ); //NUCHEMIC defines phosphate as part of sugar!

		bool curr_atom_is_base = ( curr_atom_is_sugar == false );

		calc_chem_shift += rna_cs_curr_rsd_params.atom_data( CS_data.realatomdata_index, oshi ); //reference offset

		for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

			core::conformation::Residue const & source_rsd = pose.residue( seq_num );

			RNA_CS_residue_parameters const & source_rsd_CS_params = rna_cs_params_.get_RNA_CS_residue_parameters( source_rsd.aa() );

			if ( ( CS_data.seq_num == seq_num ) && ( curr_atom_is_base ) ) continue;

			if ( include_ring_current_effect_ ) 				calc_chem_shift += ring_current_effect(       curr_atom_xyz, source_rsd, source_rsd_CS_params );
			if ( include_magnetic_anisotropy_effect_ ) calc_chem_shift += magnetic_anisotropy_effect( curr_atom_xyz, source_rsd, source_rsd_CS_params );

		}


		if ( false ){ /*verbose_*/

			std::cout << rna_cs_curr_rsd_params.base_name() << " " << CS_data.seq_num << " " << CS_data.atom_name << "\n";

			std::cout << "duetores  : rc ma\n";
			for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

				core::conformation::Residue const & source_rsd = pose.residue( seq_num );

				RNA_CS_residue_parameters const & source_rsd_CS_params = rna_cs_params_.get_RNA_CS_residue_parameters( source_rsd.aa() );

				std::cout << " " << source_rsd_CS_params.base_name() << " " << seq_num << " : ";

				if ( ( CS_data.seq_num == seq_num ) && ( curr_atom_is_base ) ){
					std::cout <<  "0 0\n";
				} else{
					std::cout << ring_current_effect(       curr_atom_xyz, source_rsd, source_rsd_CS_params ) << " ";
					std::cout << magnetic_anisotropy_effect( curr_atom_xyz, source_rsd, source_rsd_CS_params ) << "\n";
				}

			}
			std::cout << "offset : " << rna_cs_curr_rsd_params.atom_data( CS_data.realatomdata_index, oshi ) << "\n";
			std::cout << "charge effect : " << 0.0 <<  "\n\n\n"; //Currently electric field effect is not yet implemented.
		}	

		return calc_chem_shift;

	}

	//////////////////////////////////////////////////////////////////////
	void
	RNA_ChemicalShiftPotential::update_calc_chem_shift_list( pose::Pose const & pose, utility::vector1 < utility::vector1 < Real > > & calc_chem_shift_list ) const{

		calc_chem_shift_list.clear();

		for ( Size outer_data_ID = 1; outer_data_ID <= EXP_chem_shift_data_list_.size(); outer_data_ID++ ){
			utility::vector1 < Real > calc_chem_shift_entry;
			calc_chem_shift_entry.clear();

			for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_list_[outer_data_ID].size(); inner_data_ID++ ){

				ChemicalShiftData const & CS_data = EXP_chem_shift_data_list_[outer_data_ID][inner_data_ID];

				Real const calc_chem_shift = get_calc_chem_shift_value( CS_data, pose );

				calc_chem_shift_entry.push_back( calc_chem_shift );

			}

			calc_chem_shift_list.push_back( calc_chem_shift_entry );
		}

		//Consistency_check:
		if ( calc_chem_shift_list.size() != EXP_chem_shift_data_list_.size() ){
			utility_exit_with_message( "calc_chem_shift_list.size() != EXP_chem_shift_data_list_.size()" );
		}

		for ( Size ii = 1; ii <= EXP_chem_shift_data_list_.size(); ii++ ){
			if ( EXP_chem_shift_data_list_[ii].size() != calc_chem_shift_list[ii].size() ){
				utility_exit_with_message( "EXP_chem_shift_data_list_[ii].size() != calc_chem_shift_list[ii].size()" );
			}
		}

	}

	/////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChemicalShiftPotential::get_best_exp_to_calc_chem_shift_mapping( utility::vector1 < ChemicalShiftData > const & EXP_chem_shift_data_entry, 
																											 				  utility::vector1 < Real > const & calc_chem_shift_entry, 
															 												 					utility::vector1 < Real > & actual_exp_chem_shift_entry,
															 												 					utility::vector1 < bool > & do_include_CS_data ) const {

		using namespace ObjexxFCL;

		actual_exp_chem_shift_entry.clear();
		do_include_CS_data.clear();

		if ( EXP_chem_shift_data_entry.size() != calc_chem_shift_entry.size() ){
			utility_exit_with_message( "EXP_chem_shift_data_entry.size() != calc_chem_shift_entry.size()" );
		}


		if ( EXP_chem_shift_data_entry.size() == 2 ){ //	Two choice either direct or flipped match.

			//Choose the one that minimize the error (best fit "least-square regression")!
			//"Least squares" means that the overall solution minimizes the sum of the squares of the errors made in solving every single equation
			//Least squares corresponds to the maximum likelihood criterion if the experimental errors have a normal distribution and can also be derived as a method of moments estimator.

			ChemicalShiftData const & CS_data_one = EXP_chem_shift_data_entry[1];
			ChemicalShiftData const & CS_data_two = EXP_chem_shift_data_entry[2];

			Real const calc_shift_one = calc_chem_shift_entry[1];
			Real const calc_shift_two = calc_chem_shift_entry[2];

			/////Consistency checks..comment out if this slows down code!/////////
			if ( CS_data_one.seq_num != CS_data_two.seq_num ) utility_exit_with_message( "CS_data_one.seq_num != CS_data_two.seq_num" );

			Size num_matching_atom_pairs_found = 0;

			if ( CS_data_one.atom_name == " H5'" && CS_data_two.atom_name == "H5''" ) num_matching_atom_pairs_found++;

			if ( CS_data_one.atom_name == "H5''" && CS_data_two.atom_name == " H5'" ) num_matching_atom_pairs_found++;

			if ( num_matching_atom_pairs_found != 1 ) utility_exit_with_message( "num_matching_atom_pairs_found = ( " + string_of( num_matching_atom_pairs_found ) + " ) != 1" );
			///////////////////////////////////////////////////////////////////////


			if ( ( H5_prime_mode_ == "LEAST_SQUARE_IGNORE_DUPLICATE" ) && ( std::fabs( CS_data_one.exp_shift - CS_data_two.exp_shift ) < 0.00001 ) ){ //A duplicated entry!

				if ( std::pow( calc_shift_one - CS_data_one.exp_shift, 2 ) < std::pow( calc_shift_two - CS_data_two.exp_shift, 2 ) ){ //Choose the one that minimize the difference!

					if ( false ){ /*verbose_*/
						std::cout << "---------------------------------------------------" << std::endl;
						std::cout << "Keeping the better_fit of the duplicated H5_prime data | calc_shift_one = " << calc_shift_one; print_chemical_shift_data( "CS_data_one:", CS_data_one, true );
						std::cout << "Ignoring the worst_fit of the duplicated H5_prime data | calc_shift_two = " << calc_shift_two; print_chemical_shift_data( "CS_data_two:", CS_data_two, true );
						std::cout << "---------------------------------------------------" << std::endl;
					}

					do_include_CS_data.push_back( true );
					do_include_CS_data.push_back( false );
					actual_exp_chem_shift_entry.push_back( CS_data_one.exp_shift );
					actual_exp_chem_shift_entry.push_back( 0.0 );


				} else{

					if ( false ){ /*verbose_*/
						std::cout << "---------------------------------------------------" << std::endl;
						std::cout << "Keeping the better_fit of the duplicated H5_prime data | calc_shift_two = " << calc_shift_two; print_chemical_shift_data( "CS_data_two:", CS_data_two, true );
						std::cout << "Ignoring the worst_fit of the duplicated H5_prime data | calc_shift_one = " << calc_shift_one; print_chemical_shift_data( "CS_data_one:", CS_data_one, true );
						std::cout << "---------------------------------------------------" << std::endl;
					}

					//Include the first exp_CS_data in the entry but not the second.
					do_include_CS_data.push_back( false );
					do_include_CS_data.push_back( true );
					actual_exp_chem_shift_entry.push_back( 0.0 );
					actual_exp_chem_shift_entry.push_back( CS_data_one.exp_shift );

				}

			} else{

				do_include_CS_data.push_back( true );
				do_include_CS_data.push_back( true );

				Real const sum_error_square_choice_one = std::pow( calc_shift_one - CS_data_one.exp_shift, 2 ) + std::pow( calc_shift_two - CS_data_two.exp_shift, 2 ); 
				Real const sum_error_square_choice_two = std::pow( calc_shift_one - CS_data_two.exp_shift, 2 ) + std::pow( calc_shift_two - CS_data_one.exp_shift, 2 );

				if ( sum_error_square_choice_one < sum_error_square_choice_two ){
					actual_exp_chem_shift_entry.push_back( CS_data_one.exp_shift );
					actual_exp_chem_shift_entry.push_back( CS_data_two.exp_shift );
				} else{
					actual_exp_chem_shift_entry.push_back( CS_data_two.exp_shift );
					actual_exp_chem_shift_entry.push_back( CS_data_one.exp_shift );	
				}
			}

		} else if ( EXP_chem_shift_data_entry.size() == 1 ){

			do_include_CS_data.push_back( true );
			actual_exp_chem_shift_entry.push_back( EXP_chem_shift_data_entry[1].exp_shift );

		} else{
			std::cout << "ERROR: EXP_chem_shift_data_entry.size() = " << EXP_chem_shift_data_entry.size() << std::endl;
			utility_exit_with_message( "EXP_chem_shift_data_entry.size() != 1 and EXP_chem_shift_data_entry.size() != 2" );
		}


	}

	/////////////////////////////////////////////////////////////////////////////
	core::Real
	RNA_ChemicalShiftPotential::get_chemical_shift_energy( utility::vector1 < utility::vector1 < Real > > const & calc_chem_shift_list ) const{

		using namespace ObjexxFCL;

		//Consistency_check:
		if ( calc_chem_shift_list.size() != EXP_chem_shift_data_list_.size() ){
			utility_exit_with_message( "cal_chem_shift_list.size() != EXP_chem_shift_data_list_.size()" );
		}

		Real chem_shift_energy = 0.0;

		for ( Size outer_data_ID = 1; outer_data_ID <= EXP_chem_shift_data_list_.size(); outer_data_ID++ ){

			utility::vector1 < ChemicalShiftData > const & EXP_chem_shift_data_entry = EXP_chem_shift_data_list_[outer_data_ID];
			utility::vector1 < Real > const & calc_chem_shift_entry = calc_chem_shift_list[outer_data_ID];

			utility::vector1 < Real > actual_exp_chem_shift_entry;
			utility::vector1 < bool > do_include_CS_data;		 

			get_best_exp_to_calc_chem_shift_mapping( EXP_chem_shift_data_entry, calc_chem_shift_entry, actual_exp_chem_shift_entry, do_include_CS_data );

			for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_entry.size(); inner_data_ID++ ){

				if ( do_include_CS_data[inner_data_ID] == false ) continue;

				Real const exp_chem_shift = actual_exp_chem_shift_entry[inner_data_ID];

				Real const calc_chem_shift = calc_chem_shift_entry[inner_data_ID];

				chem_shift_energy += std::pow( calc_chem_shift - exp_chem_shift, 2 );

			}

		}

		return chem_shift_energy;
	}

	///////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChemicalShiftPotential::finalize_total_energy( pose::Pose const & pose, EnergyMap & totals ) const {
		using namespace conformation;

		for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){ //Right now assume pure RNA molecule.
			if ( pose.residue( seq_num ).is_RNA() == false ) return;
		}

		utility::vector1 < utility::vector1 < Real > > cal_chem_shift_list;
		cal_chem_shift_list.clear();

		update_calc_chem_shift_list( pose, cal_chem_shift_list );

		Real const chem_shift_score = get_chemical_shift_energy( cal_chem_shift_list );

		totals[ rna_chem_shift ] = chem_shift_score;


	} // finalize_total_energy

	///////////////////////////////////////////////////////////////////////////////
	///Derivative of the specified CS_data atom (non-polar proton). Both ring current and magnetic_anisotropy effects
	///This function should be called once for each CS_data atom!

	void
	RNA_ChemicalShiftPotential::get_deriv_for_chemical_shift_data_atom(
		pose::Pose const & pose,
		conformation::Residue const & CS_data_rsd,
		Size const CS_data_atomno,
		Vector & f1,
		Vector & f2
	) const {

		using namespace ObjexxFCL;

		std::string const CS_data_atom_name = remove_whitespaces( CS_data_rsd.atom_name( CS_data_atomno ) );

		utility::vector1 < ChemicalShiftData > const & EXP_chem_shift_data_entry = get_matching_CS_data_entry( CS_data_rsd.seqpos(), CS_data_atom_name );

		utility::vector1 < Real > calc_chem_shift_entry;
		for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_entry.size(); inner_data_ID++ ){
			ChemicalShiftData const & CS_data = EXP_chem_shift_data_entry[inner_data_ID];	
			Real const calc_chem_shift = get_calc_chem_shift_value( CS_data, pose );
			calc_chem_shift_entry.push_back( calc_chem_shift );
		}

		utility::vector1 < Real > actual_exp_chem_shift_entry;
		utility::vector1 < bool > do_include_CS_data;		 

		//RIGHT NOW ANALYTICAL DERIV OVER PREDICTS NUMERICAL DERIV AT places where chem shift of H5' ~ chem_shift of H5''. 
		//FIX THIS by using a fade function?
		get_best_exp_to_calc_chem_shift_mapping( EXP_chem_shift_data_entry, calc_chem_shift_entry, actual_exp_chem_shift_entry, do_include_CS_data );

		for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_entry.size(); inner_data_ID++ ){

			ChemicalShiftData const & possible_CS_data = EXP_chem_shift_data_entry[inner_data_ID]; 

			if ( EXP_chem_shift_data_entry.size() != 1 ){ //DEAL WITH SPECIAL CASE (H5'/H5'' and etcs).
				if ( CS_data_atom_name != possible_CS_data.atom_name ) continue;
				//if(CS_data_atomno!=CS_data_rsd.atom_index(CS_data.atom_name)) continue;  
			}

			ChemicalShiftData const & CS_data = possible_CS_data;

			if ( do_include_CS_data[inner_data_ID] == false ) continue;

			Real const calc_chem_shift = calc_chem_shift_entry[inner_data_ID];

			Real const act_exp_chem_shift = actual_exp_chem_shift_entry[inner_data_ID];

			if ( false ){ /*verbose_*/
				std::string const pre_string = "get_deriv_for_chemical_shift_data_atom | jj = " + string_of( inner_data_ID ) + " | EXP_CS_data_entry.size() = " + string_of( EXP_chem_shift_data_entry.size() ) + ": ";
				print_chemical_shift_data( pre_string, CS_data, false /*print_data_line*/ );
			}

			if ( CS_data.res_aa != pose.residue( CS_data.seq_num ).aa() ){
				print_chemical_shift_data( "ERROR CS_data:", CS_data, true );
				std::cout << "ERROR: pose.residue( CS_data.seq_num ).aa() = " << pose.residue( CS_data.seq_num ).aa() << std::endl;
				utility_exit_with_message( "CS_data.res_aa != pose1.residue( CS_data.seq_num ).aa()" );
			}
	

			RNA_CS_residue_parameters const & CS_data_rsd_CS_params = rna_cs_params_.get_RNA_CS_residue_parameters( CS_data.res_aa );

			bool const CS_data_atom_is_sugar = ( dround( CS_data_rsd_CS_params.atom_data( CS_data.realatomdata_index, suga ) ) == 1 ); //NUCHEMIC defines phosphate as part of sugar!

			bool const CS_data_atom_is_base = ( CS_data_atom_is_sugar == false );

			Size const CS_data_atom_index = CS_data_rsd.atom_index( CS_data.atom_name );

			numeric::xyzVector< core::Real > const & CS_data_atom_xyz = CS_data_rsd.xyz( CS_data_atom_index );

			for ( Size source_seq_num = 1; source_seq_num <= pose.total_residue(); source_seq_num++ ){

				core::conformation::Residue const & source_rsd = pose.residue( source_seq_num );

				RNA_CS_residue_parameters const & source_rsd_CS_params = rna_cs_params_.get_RNA_CS_residue_parameters( source_rsd.aa() );

				if ( ( CS_data.seq_num == source_seq_num ) && ( CS_data_atom_is_base ) ) continue;

				if ( include_ring_current_effect_ ){
					///Ring current effects
					for ( Size ring_ID = 1; ring_ID <= source_rsd_CS_params.num_rings(); ring_ID++ ){

						//get_ring_current_deriv() gives gradient of ring_current_effect() wrt to r_vector 
						// +1.0 since r_vector = CS_data_atom_xyz - molecular_ring_center.
						numeric::xyzVector< core::Real > const f2_calc_chem_shift = + 1.0 * get_ring_current_deriv( CS_data_atom_xyz, source_rsd, ring_ID, source_rsd_CS_params ); 
						numeric::xyzVector< core::Real > const f1_calc_chem_shift =  cross( f2_calc_chem_shift, CS_data_atom_xyz ); 

						f1 += 2.0*( calc_chem_shift - act_exp_chem_shift )*( f1_calc_chem_shift ); //Convert deriv_vector of calc_chem_shift to deriv_vector of chemical_shift_energy
						f2 += 2.0*( calc_chem_shift - act_exp_chem_shift )*( f2_calc_chem_shift ); //Convert deriv_vector of calc_chem_shift to deriv_vector of chemical_shift_energy

					}
				}

				if ( include_magnetic_anisotropy_effect_ ){
					//Magnetic_anisotropy effects 
					Size const source_rsd_maxatoms = source_rsd_CS_params.get_atomnames_size();

					numeric::xyzMatrix< core::Real > const source_base_coordinate_matrix = get_rna_base_coordinate_system_from_CS_params( source_rsd, source_rsd_CS_params );

					for ( Size source_realatomdata_index = 1; source_realatomdata_index <= source_rsd_maxatoms; source_realatomdata_index++ ){

						if ( dround( source_rsd_CS_params.atom_data( source_realatomdata_index, maca ) ) != 1 ) continue;

						Size const source_atom_index = source_rsd.atom_index( source_rsd_CS_params.get_atomname( source_realatomdata_index ) );

						numeric::xyzVector< core::Real > const & source_atom_xyz = source_rsd.xyz( source_atom_index );

						//get_delta_magnetic_anisotropy_deriv() gives gradient of delta_magnetic_anisotropy() with respect to r_vector.
						// +1.0 since r_vector = CS_data_atom_xyz - source_atom_xyz.
						numeric::xyzVector< core::Real > const f2_calc_chem_shift = + 1.0 * get_delta_magnetic_anisotropy_deriv( CS_data_atom_xyz, source_atom_xyz, source_base_coordinate_matrix, source_rsd_CS_params, source_realatomdata_index );
											 							
						numeric::xyzVector< core::Real > const f1_calc_chem_shift =  cross( f2_calc_chem_shift, CS_data_atom_xyz );

						f1 += 2.0*( calc_chem_shift - act_exp_chem_shift )*( f1_calc_chem_shift ); //Convert deriv_vector of calc_chem_shift to deriv_vector of chemical_shift_energy
						f2 += 2.0*( calc_chem_shift - act_exp_chem_shift )*( f2_calc_chem_shift ); //Convert deriv_vector of calc_chem_shift to deriv_vector of chemical_shift_energy


					}
				}


			}

		}

	}


	///////////////////////////////////////////////////////////////////////////////
	///Derivative due to the ring_current effect of the source base. Include contribution from all CS_data atoms (non-polar protons). 
	///Note, there are 1 or 2 ring centers per base.
	///This function should be called once for each residue only at the first_base_atomno. 

	void
	RNA_ChemicalShiftPotential::get_ring_current_deriv_for_src_base(	
		pose::Pose const & pose,
		conformation::Residue const & rc_source_rsd,
		Size const chi1_torsion_atomnno,
		Vector & f1,
		Vector & f2
	) const {

		if ( include_ring_current_effect_ == false ) return;

		Size const rc_source_seq_num = rc_source_rsd.seqpos();
	
		RNA_CS_residue_parameters const & rc_source_rsd_CS_params = rna_cs_params_.get_RNA_CS_residue_parameters( rc_source_rsd.aa() );

		std::string const chi1_torsion_atomn_name = rc_source_rsd.atom_name( chi1_torsion_atomnno );

		if ( false ){ /*verbose_*/
			std::cout << "get_deriv_for_ring_current_center for chi1_torsion_atomn_name = " << chi1_torsion_atomn_name;	
			std::cout << " | name_from_aa( rc_source_rsd.aa() ) = " << name_from_aa( rc_source_rsd.aa() );  
			std::cout << " | rc_source_seq_num = " << rc_source_seq_num;	
		}

		//Enumerate through all the chemical_shift data points (right now only non_polar protons).
		for ( Size outer_data_ID = 1; outer_data_ID <= EXP_chem_shift_data_list_.size(); outer_data_ID++ ){

			utility::vector1 < ChemicalShiftData > const & EXP_chem_shift_data_entry = EXP_chem_shift_data_list_[outer_data_ID];

			utility::vector1 < Real > calc_chem_shift_entry;
			for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_entry.size(); inner_data_ID++ ){

					ChemicalShiftData const & CS_data = EXP_chem_shift_data_entry[inner_data_ID];	
					Real const calc_chem_shift = get_calc_chem_shift_value( CS_data, pose );
					calc_chem_shift_entry.push_back( calc_chem_shift );
			}

			utility::vector1 < Real > actual_exp_chem_shift_entry;
			utility::vector1 < bool > do_include_CS_data;		 

			get_best_exp_to_calc_chem_shift_mapping( EXP_chem_shift_data_entry, calc_chem_shift_entry, actual_exp_chem_shift_entry, do_include_CS_data );

			for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_entry.size(); inner_data_ID++ ){

				if ( do_include_CS_data[inner_data_ID] == false ) continue;

				ChemicalShiftData const & CS_data = EXP_chem_shift_data_entry[inner_data_ID];	

				Real const calc_chem_shift = calc_chem_shift_entry[inner_data_ID];

				Real const act_exp_chem_shift = actual_exp_chem_shift_entry[inner_data_ID];

				core::conformation::Residue const & CS_data_rsd = pose.residue( CS_data.seq_num );

				Size const CS_data_atom_index = CS_data_rsd.atom_index( CS_data.atom_name );

				numeric::xyzVector< core::Real > const & CS_data_atom_xyz = CS_data_rsd.xyz( CS_data_atom_index );

				RNA_CS_residue_parameters const & CS_data_rsd_CS_params = rna_cs_params_.get_RNA_CS_residue_parameters( CS_data.res_aa );

				bool const CS_data_atom_is_sugar = ( dround( CS_data_rsd_CS_params.atom_data( CS_data.realatomdata_index, suga ) ) == 1 ); //NUCHEMIC defines phosphate as part of sugar!

				bool const CS_data_atom_is_base = ( CS_data_atom_is_sugar == false );

				if ( ( CS_data.seq_num == rc_source_seq_num ) && ( CS_data_atom_is_base ) ) continue;

				for ( Size rc_source_ring_ID = 1; rc_source_ring_ID <= rc_source_rsd_CS_params.num_rings(); rc_source_ring_ID++ ){

					//get_ring_current_deriv() gives gradient of ring_current_effect() wrt to r_vector
					// -1.0 since r_vector = CS_data_atom_xyz - molecular_ring_center.
					numeric::xyzVector< core::Real > const f2_calc_chem_shift = -1.0 * get_ring_current_deriv( CS_data_atom_xyz, rc_source_rsd, rc_source_ring_ID, rc_source_rsd_CS_params ); 
					numeric::xyzVector< core::Real > const f1_calc_chem_shift =  cross( f2_calc_chem_shift, CS_data_atom_xyz ); //cross( f2_calc_chem_shift, ring_center_xyz); 

					f1 += 2.0*( calc_chem_shift - act_exp_chem_shift )*( f1_calc_chem_shift ); //Convert deriv_vector of calc_chem_shift to deriv_vector of chemical_shift_energy
					f2 += 2.0*( calc_chem_shift - act_exp_chem_shift )*( f2_calc_chem_shift ); //Convert deriv_vector of calc_chem_shift to deriv_vector of chemical_shift_energy
			
				}
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	///Derivative due to the manisotropy effect of the source base. Include contribution from all CS_data atoms (non-polar protons). 
	///This function should be called once for each residue only at the chi1_torsion_atomno. 
	void
	RNA_ChemicalShiftPotential::get_magnetic_anisotropy_deriv_for_src_base(
		pose::Pose const & pose,
		conformation::Residue const & ma_source_rsd,
		Size const chi1_torsion_atomno,
		Vector & f1,
		Vector & f2
	) const {

		if ( include_magnetic_anisotropy_effect_ == false ) return;

		Size const ma_source_seq_num = ma_source_rsd.seqpos();
	
		RNA_CS_residue_parameters const & ma_source_rsd_CS_params = rna_cs_params_.get_RNA_CS_residue_parameters( ma_source_rsd.aa() );

		std::string const chi1_torsion_atom_name = ma_source_rsd.atom_name( chi1_torsion_atomno );

		Size const ma_source_rsd_maxatoms = ma_source_rsd_CS_params.get_atomnames_size();

		numeric::xyzMatrix< core::Real > const ma_source_base_coordinate_matrix = get_rna_base_coordinate_system_from_CS_params( ma_source_rsd, ma_source_rsd_CS_params );

		if ( false ){ /*verbose_*/
			std::cout << "get_deriv_for_magnetic_anisotropy_src_atom for chi1_torsion_atom_name = " << chi1_torsion_atom_name;
			std::cout << " | name_from_aa( ma_source_rsd.aa() ) = " << name_from_aa( ma_source_rsd.aa() );  
			std::cout << " | ma_source_seq_num = " << ma_source_seq_num << std::endl;
		}

		//Enumerate through all the chemical_shift data points (right now only non_polar protons).
		for ( Size outer_data_ID = 1; outer_data_ID <= EXP_chem_shift_data_list_.size(); outer_data_ID++ ){

			utility::vector1 < ChemicalShiftData > const & EXP_chem_shift_data_entry = EXP_chem_shift_data_list_[outer_data_ID];

			utility::vector1 < Real > calc_chem_shift_entry;
			for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_entry.size(); inner_data_ID++ ){

					ChemicalShiftData const & CS_data = EXP_chem_shift_data_entry[inner_data_ID];	
					Real const calc_chem_shift = get_calc_chem_shift_value( CS_data, pose );
					calc_chem_shift_entry.push_back( calc_chem_shift );
			}

			utility::vector1 < Real > actual_exp_chem_shift_entry;
			utility::vector1 < bool > do_include_CS_data;		 

			//RIGHT NOW ANALYTICAL DERIV OVER PREDICTS NUMERICAL DERIV AT places where chem shift of H5' ~ chem_shift of H5''. 
			//FIX THIS by using a fade function?
			get_best_exp_to_calc_chem_shift_mapping( EXP_chem_shift_data_entry, calc_chem_shift_entry, actual_exp_chem_shift_entry, do_include_CS_data );

			for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_entry.size(); inner_data_ID++ ){

				if ( do_include_CS_data[inner_data_ID] == false ) continue;

				ChemicalShiftData const & CS_data = EXP_chem_shift_data_entry[inner_data_ID];	

				Real const calc_chem_shift = calc_chem_shift_entry[inner_data_ID];

				Real const act_exp_chem_shift = actual_exp_chem_shift_entry[inner_data_ID];

				core::conformation::Residue const & CS_data_rsd = pose.residue( CS_data.seq_num );

				Size const CS_data_atom_index = CS_data_rsd.atom_index( CS_data.atom_name );

				numeric::xyzVector< core::Real > const & CS_data_atom_xyz = CS_data_rsd.xyz( CS_data_atom_index );

				RNA_CS_residue_parameters const & CS_data_rsd_CS_params = rna_cs_params_.get_RNA_CS_residue_parameters( CS_data.res_aa );

				bool const CS_data_atom_is_sugar = ( dround( CS_data_rsd_CS_params.atom_data( CS_data.realatomdata_index, suga ) ) == 1 ); //NUCHEMIC defines phosphate as part of sugar!

				bool const CS_data_atom_is_base = ( CS_data_atom_is_sugar == false );

				if ( ( CS_data.seq_num == ma_source_seq_num ) && ( CS_data_atom_is_base ) ) continue;

				for ( Size source_realatomdata_index = 1; source_realatomdata_index <= ma_source_rsd_maxatoms; source_realatomdata_index++ ){

					if ( dround( ma_source_rsd_CS_params.atom_data( source_realatomdata_index, maca ) ) != 1 ) continue;

					Size const ma_source_atom_index = ma_source_rsd.atom_index( ma_source_rsd_CS_params.get_atomname( source_realatomdata_index ) );

					numeric::xyzVector< core::Real > const & ma_source_atom_xyz = ma_source_rsd.xyz( ma_source_atom_index );

					//get_delta_magnetic_anisotropy_deriv() gives gradient of delta_magnetic_anisotropy() with respect to r_vector.
					// -1.0 since r_vector = CS_data_atom_xyz - source_atom_xyz.
					numeric::xyzVector< core::Real > const f2_calc_chem_shift = -1.0 * get_delta_magnetic_anisotropy_deriv( CS_data_atom_xyz, ma_source_atom_xyz, ma_source_base_coordinate_matrix, ma_source_rsd_CS_params, source_realatomdata_index );
										 							
					numeric::xyzVector< core::Real > const f1_calc_chem_shift =  cross( f2_calc_chem_shift, CS_data_atom_xyz );

					f1 += 2.0*( calc_chem_shift - act_exp_chem_shift )*( f1_calc_chem_shift ); //Convert deriv_vector of calc_chem_shift to deriv_vector of chemical_shift_energy
					f2 += 2.0*( calc_chem_shift - act_exp_chem_shift )*( f2_calc_chem_shift ); //Convert deriv_vector of calc_chem_shift to deriv_vector of chemical_shift_energy


				}
			}
		}

	}

	///////////////////////////////////////////////////////////////////////////////
	///Derivative at the magnetic_anisotropy_src_atom (right how only base heavy-atoms) due to the magnetic_anistropy_effect.
	///Include contribution from all CS_data atoms (non-polar protons). 
	void
	RNA_ChemicalShiftPotential::get_magnetic_anisotropy_deriv_for_src_atom(
		pose::Pose const & pose,
		conformation::Residue const & ma_source_rsd,
		Size const ma_source_atomno,
		Vector & f1,
		Vector & f2
	) const {

		if ( include_magnetic_anisotropy_effect_ == false ) return;

		Size const ma_source_seq_num = ma_source_rsd.seqpos();
	
		RNA_CS_residue_parameters const & ma_source_rsd_CS_params = rna_cs_params_.get_RNA_CS_residue_parameters( ma_source_rsd.aa() );

		std::string const ma_source_atom_name = ma_source_rsd.atom_name( ma_source_atomno );

		numeric::xyzMatrix< core::Real > const ma_source_base_coordinate_matrix = get_rna_base_coordinate_system_from_CS_params( ma_source_rsd, ma_source_rsd_CS_params );

		numeric::xyzVector< core::Real > const & ma_source_atom_xyz = ma_source_rsd.xyz( ma_source_atomno );

		Size const ma_source_realatomdata_index = get_realatomdata_index( ma_source_atom_name, ma_source_rsd.aa() );

		if ( dround( ma_source_rsd_CS_params.atom_data( ma_source_realatomdata_index, maca ) ) != 1 ){
			utility_exit_with_message( "dround( ma_source_rsd_CS_params.atom_data( ma_source_realatomdata_index, maca ) ) != 1" );
		}

		if ( false ){ /*verbose_*/
			std::cout << "get_deriv_for_magnetic_anisotropy_src_atom for ma_source_atom_name = " << ma_source_atom_name;	
			std::cout << " | name_from_aa( ma_source_rsd.aa() ) = " << name_from_aa( ma_source_rsd.aa() );  
			std::cout << " | ma_source_seq_num = " << ma_source_seq_num << std::endl;
		}

		//Enumerate through all the chemical_shift data points (right now only non_polar protons).
		for ( Size outer_data_ID = 1; outer_data_ID <= EXP_chem_shift_data_list_.size(); outer_data_ID++ ){

			utility::vector1 < ChemicalShiftData > const & EXP_chem_shift_data_entry = EXP_chem_shift_data_list_[outer_data_ID];

			utility::vector1 < Real > calc_chem_shift_entry;
			for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_entry.size(); inner_data_ID++ ){

					ChemicalShiftData const & CS_data = EXP_chem_shift_data_entry[inner_data_ID];	
					Real const calc_chem_shift = get_calc_chem_shift_value( CS_data, pose );
					calc_chem_shift_entry.push_back( calc_chem_shift );
			}

			utility::vector1 < Real > actual_exp_chem_shift_entry;
			utility::vector1 < bool > do_include_CS_data;		 

			//RIGHT NOW ANALYTICAL DERIV OVER PREDICTS NUMERICAL DERIV AT places where chem shift of H5' ~ chem_shift of H5''. 
			//FIX THIS by using a fade function?
			get_best_exp_to_calc_chem_shift_mapping( EXP_chem_shift_data_entry, calc_chem_shift_entry, actual_exp_chem_shift_entry, do_include_CS_data );

			for ( Size inner_data_ID = 1; inner_data_ID <= EXP_chem_shift_data_entry.size(); inner_data_ID++ ){

				if ( do_include_CS_data[inner_data_ID] == false ) continue;

				ChemicalShiftData const & CS_data = EXP_chem_shift_data_entry[inner_data_ID];	

				Real const calc_chem_shift = calc_chem_shift_entry[inner_data_ID];

				Real const act_exp_chem_shift = actual_exp_chem_shift_entry[inner_data_ID];

				core::conformation::Residue const & CS_data_rsd = pose.residue( CS_data.seq_num );

				Size const CS_data_atom_index = CS_data_rsd.atom_index( CS_data.atom_name );

				numeric::xyzVector< core::Real > const & CS_data_atom_xyz = CS_data_rsd.xyz( CS_data_atom_index );

				RNA_CS_residue_parameters const & CS_data_rsd_CS_params = rna_cs_params_.get_RNA_CS_residue_parameters( CS_data.res_aa );

				bool const CS_data_atom_is_sugar = ( dround( CS_data_rsd_CS_params.atom_data( CS_data.realatomdata_index, suga ) ) == 1 ); //NUCHEMIC defines phosphate as part of sugar!

				bool const CS_data_atom_is_base = ( CS_data_atom_is_sugar == false );

				if ( ( CS_data.seq_num == ma_source_seq_num ) && ( CS_data_atom_is_base ) ) continue;

				//get_delta_magnetic_anisotropy_deriv() gives gradient of delta_magnetic_anisotropy() with respect to r_vector.
				// -1.0 since r_vector = CS_data_atom_xyz - ma_source_atom_xyz.
				numeric::xyzVector< core::Real > const f2_calc_chem_shift = -1.0 * get_delta_magnetic_anisotropy_deriv( CS_data_atom_xyz, ma_source_atom_xyz, ma_source_base_coordinate_matrix, ma_source_rsd_CS_params, ma_source_realatomdata_index );
									 							
				numeric::xyzVector< core::Real > const f1_calc_chem_shift =  cross( f2_calc_chem_shift, CS_data_atom_xyz );

				//RIGHT NOW NUMERICAL AND ANALYTICAL DERIV DOESN'T AGREE NEAR zero due to switch in H5'/H5'' atom pairs.
				//DO TO: Include a fade function to fix this!
				f1 += 2.0*( calc_chem_shift - act_exp_chem_shift )*( f1_calc_chem_shift ); //Convert deriv_vector of calc_chem_shift to deriv_vector of chemical_shift_energy
				f2 += 2.0*( calc_chem_shift - act_exp_chem_shift )*( f2_calc_chem_shift ); //Convert deriv_vector of calc_chem_shift to deriv_vector of chemical_shift_energy
			

			}
		}

	}



	///////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChemicalShiftPotential::eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	 	) const
	{

		using namespace conformation;

		Size const seq_num = atom_id.rsd();
		Size const atomno = atom_id.atomno();

		conformation::Residue const & rsd = pose.residue( seq_num );
		if ( !rsd.is_RNA() ) return;

		//OK two possible force contributions (the two should be mutually exclusive!).
		//i.  Contain chemical_shift data on this atom (non-polar protons)
		//ii. Atom is source of ring_current B-field (most base heavy-atoms) or the magnetic_anisotropy B-field (all base heavy-atoms)

		//Vector f1( 0.0 ),f2( 0.0 );

		numeric::xyzVector< core::Real > f1( 0.0, 0.0, 0.0 );
		numeric::xyzVector< core::Real > f2( 0.0, 0.0, 0.0 );

		if ( false ){
			std::cout << "eval_atom_derivative() for name_from_aa( rsd.aa() ) = " << name_from_aa( rsd.aa() );
			std::cout << " | seq_num = " << seq_num;
			std::cout << " | atomno = " << atomno;
			std::cout << " | atom_name = " << rsd.atom_name( atomno );
			std::cout << std::endl;
		}

		if ( atom_has_exp_chemical_shift_data( rsd, atomno ) ){

			get_deriv_for_chemical_shift_data_atom( pose, rsd, atomno, f1, f2 );

		}


		if ( atomno == chemical::rna::chi1_torsion_atom_index( rsd ) )	{ //first chi1_torsion_atom serves as 'proxy' for the ring_center!

			get_ring_current_deriv_for_src_base( pose, rsd, atomno, f1, f2 );

			get_magnetic_anisotropy_deriv_for_src_base( pose, rsd, atomno, f1, f2 );

		}


		//I don't fully understand this. But if I get_deriv at the actual location of the ma source atoms, then the analytical/numeric deriv ratio at the first chi1_torsion_atom will be slightly by ~5%-35% percent. To get the correct derivative, need get deriv at the chi1_torsion_atom instead (see above).
 	
		//if(Is_magnetic_anisotropy_source_atom(rsd, atomno)){
			//get_deriv_for_magnetic_anisotropy_src_atom(pose, rsd, atomno, f1, f2);
		//}


		F1 += weights[ rna_chem_shift ] * f1;
		F2 += weights[ rna_chem_shift ] * f2;


	} // eval atom derivative



} // chemical_shift
} // rna
} // scoring
} // core
