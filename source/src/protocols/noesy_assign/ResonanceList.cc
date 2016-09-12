// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @details
///
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/Resonance.hh>
#include <protocols/noesy_assign/ResonanceList.hh>
#include <protocols/noesy_assign/FloatingResonance.hh>
#include <protocols/noesy_assign/ProtonResonance.hh>
#include <protocols/noesy_assign/LabelResonance.hh>


// Package Headers
#include <protocols/noesy_assign/Exceptions.hh>
#include <core/id/NamedAtomID.hh>
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>

// Utility headers
#include <ObjexxFCL/format.hh>
#include <utility>
#include <utility/string_util.hh>
// #include <utility/exit.hh>
// #include <utility/excn/Exceptions.hh>
// #include <utility/vector1.fwd.hh>
// #include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
// #include <basic/prof.hh>
#include <basic/Tracer.hh>
// #include <basic/options/option.hh>
// #include <basic/options/keys/abinitio.OptionKeys.gen.hh>
// #include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <deque>
#include <sstream>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.noesy_assign.resonances" );
static THREAD_LOCAL basic::Tracer tr_labels( "protocols.noesy_assign.resonances.labels" );

using core::Real;
using namespace core;
using namespace basic;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

namespace protocols {
namespace noesy_assign {


ResonanceList::ResonanceList( std::string  sequence ) : sequence_(std::move( sequence ))
{}

ResonanceList::~ResonanceList() = default;

/// translate sequence information into AA
core::chemical::AA ResonanceList::aa_from_resid( core::Size resi ) const {
	runtime_assert(  resi <= sequence_.size() );
	return core::chemical::aa_from_oneletter_code( sequence_[ resi-1 ] );
}


/* helper function for 'read_from_stream'
each element is first put into a DEQUE such that atoms like methyl-protons can be combined if they are assigned the same frequency.
-HB1, HB2 -- > QB

example:
input:
HB1  1.00
HB2  1.00
CD  13.00

output:
QB   1.00
CD  13.00
*/
void process_last_resonances( std::deque< ResonanceOP >& last_resonances, bool drain=false ) {
	ResonanceOP first_res = last_resonances.front();
	last_resonances.pop_front();
	first_res->combine( last_resonances, drain );
	last_resonances.push_front( first_res );
}

void ResonanceList::read_from_stream( std::istream& is ) {
	using namespace core::chemical; //for AA
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	std::string line;
	std::deque< ResonanceOP > last_resonances;
	bool has_ambiguity_column( false );
	std::string sequence_from_header( "" );
	while ( getline( is, line ) ) {
		core::Size label;
		core::Real freq;
		core::Real error;
		std::string name;
		core::Size resn;
		AA aa;
		if ( line.find( "##" ) < 10 ) {
			std::istringstream line_stream( line.substr( line.find("##")+2 ) );
			std::string comment_type;
			line_stream >> comment_type;
			if ( comment_type == "VARS" ) {
				has_ambiguity_column = line.find( "AMBIGUITY" );
				tr.Debug << "found ambiguity column" << std::endl;
			}
			if ( comment_type == "DATA" ) {
				std::string data_type;
				line_stream >> data_type;
				if ( data_type == "SEQUENCE" ) {
					while ( line_stream.good() ) {
						std::string sequence_bit;
						line_stream >> sequence_bit;
						sequence_from_header+=sequence_bit;
					}
				}
			}
			continue;
		}
		std::istringstream line_stream( line );
		//read from stream...
		tr.Trace << "read line: " << line << std::endl;
		line_stream >> label >> freq >> error >> name;
		if ( params.ignore_resonancefile_tolerances_ ) error = 0.0;

		//check for stream-error
		if ( !line_stream.good() ) {
			if ( line.size() ) {
				tr.Info << "ignore line : " << line << std::endl;
			}
			continue; //ignore weird lines
		}

		if ( name=="" ) break;//got an empty line at end of file ?

		//read optional fields
		if ( has_ambiguity_column ) {
			core::Size ambiguity_code;
			line_stream >> ambiguity_code;
			//we ignore this code, but should be able to deal if it is present.
			//what counts for us are the float-groups...
		}
		line_stream >> resn;
		std::string aa_name;
		line_stream >> aa_name;

		//process optional fields...
		if ( line_stream.good() ) { //optional field present?
			if ( aa_name.size() == 1 ) {
				aa = aa_from_oneletter_code( aa_name[ 0 ] );
			} else if ( aa_name.size() == 3 ) {
				aa = aa_from_name( aa_name );
			} else {
				throw utility::excn::EXCN_BadInput( "did not recognize aminoacid: " + aa_name);
			}
			if ( sequence_.size() < resn ) {
				while ( sequence_.size() < resn-1 ) sequence_.push_back('X');
				sequence_[ resn-1 ]=oneletter_code_from_aa( aa );
			} else if ( sequence_[ resn-1 ]!=oneletter_code_from_aa( aa ) ) {
				tr.Warning << "sequence mismatch!! check your data: found " << name_from_aa( aa ) << " in line " << line
					<< " which does not match " << sequence_[ resn-1 ] << "\n" << sequence_ << std::endl;
			}
		} else { //optional field was not present... use the pre-supplied fasta-sequence
			aa = aa_from_resid( resn );
		}
		line_stream >> aa_name; //1-letter

		Real intensity( 1.0 );
		// 3/8/13: this is totally buggy, as input files are not curated accordingly
		// remove completely for now.
		// line_stream >> intensity;
		//  if ( !line_stream ) { //option field intensity present ?
		//   intensity=1.0;
		//  }

		std::string fl_tag;
		line_stream >> fl_tag;
		typedef std::set< core::Size > FloatList;
		FloatList floats;
		if ( line_stream ) {
			if ( fl_tag!="None" ) {
				if ( fl_tag[0]!='#' && fl_tag[0]!='[' ) {
					throw utility::excn::EXCN_BadInput( "did not recognize item " +fl_tag + " in line " + line );
				} else { //that means fl_tag[0]=='[' || fl_tag[0]=='#'
					std::stringstream float_info;
					fl_tag[0]=' ';
					bool closed( false );
					if ( fl_tag[ fl_tag.size()-1 ]==']' ) {
						fl_tag.replace( fl_tag.size()-1, 1, " ");
						closed = true;
					}
					float_info << fl_tag;
					while ( line_stream && !closed ) {
						line_stream >> fl_tag;
						if ( fl_tag[ fl_tag.size()-1 ]==']' ) {
							fl_tag.replace( fl_tag.size()-1, 1, " ");
							closed = true;
						}
						float_info << fl_tag;
					}
					if ( !closed ) {
						throw utility::excn::EXCN_BadInput( "expect closing ] in line " + line );
					}
					char cstr[100];
					while ( float_info.getline(cstr, 50, ',' ) ) {
						tr.Debug << "add " << cstr << " to float group for resonance " << label << std::endl;
						floats.insert( utility::string2int( cstr ));
						tr.Debug << floats.size() << std::endl;
					}
				}
			}
		}

		// replace some names...
		if ( name == "HN" ) name ="H";
		if ( aa == aa_leu && name == "HD1" ) name = "QD1";
		if ( aa == aa_leu && name == "HD2" ) name = "QD2";
		if ( aa == aa_val && name == "HG1" ) name = "QG1";
		if ( aa == aa_val && name == "HG2" ) name = "QG2";
		if ( aa == aa_thr && name == "HG2" ) name = "QG2";
		if ( aa == aa_ala && name == "HB"  ) name = "QB";
		if ( aa == aa_gly && name == "HA"  ) name = "QA";
		if ( aa == aa_ile && name == "HG1" ) name = "QG1";
		if ( aa == aa_ile && name == "HG2" ) name = "QG2";
		if ( aa == aa_ile && name == "HD1" ) name = "QD1";

		bool is_proton = ( name[ 0 ]=='Q' || ( name.find("H") != std::string::npos && name[ 0 ] != 'C' /*not CH2 on TRP*/ ) );
		ResonanceOP save_resonance( nullptr );
		if ( is_proton ) {
			save_resonance = ResonanceOP( new ProtonResonance( label, freq, error, core::id::NamedAtomID( name, resn ), aa, intensity ) );
		} else {
			save_resonance = ResonanceOP( new LabelResonance( label, freq, error,  core::id::NamedAtomID( name, resn ), aa, intensity ) );
		}

		if ( !floats.empty() ) { //size() ) {
			save_resonance = ResonanceOP( new FloatingResonance( *save_resonance, floats, this ) );
		}

		// before assigning to Resonance List put it in DEQUE (push_back), maybe it will get combined with next resonance...
		last_resonances.push_back( save_resonance );
		if ( freq != last_resonances.front()->freq() ) { ///if we have just read a new frequency we need to finalize what is in the DEQUE
			if ( last_resonances.size() > 2 ) process_last_resonances( last_resonances ); //just 2 --> 1 old freq and 1 new freq
			map_[ last_resonances.front()->label() ] = last_resonances.front(); ///assign most forward value in our DEQUE
			last_resonances.pop_front(); ///remove the just-assigned value
		}
	}

	//still unprocessed data in DEQUE ?
	if ( last_resonances.size() > 1 ) process_last_resonances( last_resonances, true /*drain*/ );
	map_[ last_resonances.front()->label() ]= last_resonances.front();

	//ASSERT: we have captured everything...
	runtime_assert( last_resonances.size() == 1 );

	///no problems ?
	if ( is.fail() && is.eof() && is.bad() ) {
		tr.Error << "[ERROR WHILE READING]" << std::endl;
	}

	//post-process input...
	update_residue_map();

	//fix intensities...
	if ( params.ignore_resonancefile_intensities_ ) {
		for ( ResonanceIDs::const_iterator it = map_.begin(); it != map_.end(); ++it ) {
			using core::scoring::constraints::parse_NMR_name;
			core::scoring::constraints::NamedAtoms atoms;
			parse_NMR_name( it->second->name(), it->second->resid(), it->second->aa(), atoms );
			it->second->set_intensity( atoms.size() );
		}
	}

	if ( sequence_from_header.size() ) {
		tr.Debug << "sequence information obtained from file-header\n";
		tr.Debug << sequence_from_header << std::endl;
		if ( sequence_.size() ) {
			if ( sequence_from_header.size() >= sequence_.size() ) {
				for ( core::Size i=0; i < sequence_.size(); i++ ) {
					if ( sequence_[ i ] == 'X' ) {
						sequence_[ i ] = sequence_from_header[ i ];
					} else if ( sequence_[ i ] != sequence_from_header[ i ] ) {
						throw utility::excn::EXCN_BadInput( std::string("sequence information in header is not consistent with sequence letter ")+
							"in resonance lines at residue "+ObjexxFCL::string_of( i+1 ) );
					}
				}
			}
		}
		tr.Debug << "";
	}
	update_bond_connections();
}

/// @brief write ResonanceList to stream
void ResonanceList::write_to_stream( std::ostream& os   ) const {
	for ( auto it = map_.begin(); it != map_.end(); ++it ) {
		runtime_assert( it->first == it->second->label() );
		if ( sequence_.size() ) {
			using namespace core::chemical; //for AA
			AA aa( aa_from_resid( it->second->resid() ) );
			it->second->write_to_stream( os, aa );
		} else {
			it->second->write_to_stream( os );
		}
		os << std::endl;
	}
}

/// @brief write ResonanceList in TALOS format
void ResonanceList::write_talos_format( std::ostream& os, bool backbone_only ) const {
	if ( sequence_.size() == 0 ) {
		utility_exit_with_message( "sequence information required to write TALOS format -- use -in:file:fasta" );
	}

	//write sequence section
	Size const TALOS_SEQ_LINE_LEN( 50 );
	Size const TALOS_SEQ_BLCK_SIZE( 10 );
	for ( Size ct=0; ct<sequence_.size(); ct++ ) {
		if ( ct % TALOS_SEQ_LINE_LEN == 0 && ct > 1 ) os <<"\n";
		if ( ct % TALOS_SEQ_LINE_LEN == 0 ) os << "DATA SEQUENCE";
		if ( ct % TALOS_SEQ_BLCK_SIZE == 0 ) os << " ";
		os << sequence_[ ct ];
	}

	///write format section
	os << "\n\n";
	os << "VARS RESID RESNAME ATOMNAME SHIFT\n";
	os << "FORMAT %4d %1s %4s %8.3f\n";
	os << "\n";


	///write resonances
	for ( auto it = map_.begin(); it != map_.end(); ++it ) {
		runtime_assert( it->first == it->second->label() );
		if ( sequence_.size() < it->second->resid() ) {
			tr.Error << " no sequence information for residue " << it->second->resid() << std::endl;
			utility_exit_with_message( "sequence information required for all residues to write TALOS format -- use -in:file:fasta" );
		}
		using namespace core::chemical; //for AA
		std::string atom( it->second->name() );

		//are we backbone
		bool const is_backbone( //backbone im Sinne von SPARTA
			atom=="H" || atom=="HN" ||
			atom=="CA" || atom=="CB" ||
			atom=="N" ||
			atom=="1HA" || atom=="2HA" || atom=="3HA" );
		//write individual atom
		if ( is_backbone || !backbone_only ) {
			if ( atom=="1HA" ) atom="HA3";
			if ( atom=="2HA" ) atom="HA2";
			if ( atom=="3HA" ) atom="HA1";
			AA aa( aa_from_resid( it->second->resid() ) );
			os << ObjexxFCL::format::RJ( 5, it->second->resid() ) << " ";
			os << oneletter_code_from_aa( aa ) << " ";
			os << ObjexxFCL::format::RJ( 3, atom=="H" ? "HN" : atom ) << " ";
			os << ObjexxFCL::format::F( 7, 3, it->second->freq() ) << std::endl;
		}
	}
	os << std::endl;
}


///retrieve Resonance --- throws EXCN_UnknonwResonance if atom not found
Resonance const& ResonanceList::operator[] ( core::id::NamedAtomID const& atom ) const {
	auto it_res( by_resid_.find( atom.rsd() ) );
	if ( it_res != by_resid_.end() ) {
		Resonances const& reso_list( it_res->second );
		for ( auto const & it : reso_list ) {
			if ( it->atom() == atom ) return *it;
		}
	}
	throw EXCN_UnknownResonance( atom, "can't find atom ");
	return *(map_.begin()->second); //to make compiler happy
}

///retrieve Resonance ---  throws EXCN_UnknonwResonance if atom not found
Resonance const& ResonanceList::operator[] ( core::Size key ) const {
	auto iter = map_.find( key );
	if ( iter == map_.end() ) {
		throw EXCN_UnknownResonance( id::BOGUS_NAMED_ATOM_ID, "can't find resonance " + ObjexxFCL::string_of( key ) );
	}
	return *(iter->second);
}

///create map with all resonances sorted by residue number
void ResonanceList::update_residue_map() {
	by_resid_.clear();
	for ( ResonanceIDs::const_iterator it = map_.begin(); it != map_.end(); ++it ) {
		runtime_assert( it->first == it->second->label() );
		by_resid_[ it->second->resid() ].push_back( it->second );
	}
}

///retrieve list of Resonance at certain residue ---  throws EXCN_UnknonwResonance if residue number not found
ResonanceList::Resonances const& ResonanceList::resonances_at_residue( core::Size resid ) const {
	auto it_res( by_resid_.find( resid ) );
	if ( it_res != by_resid_.end() ) {
		return it_res->second;
	}
	throw EXCN_UnknownResonance( id::BOGUS_NAMED_ATOM_ID, "can't find resonance with residue " + ObjexxFCL::string_of( resid ) );
	return by_resid_.begin()->second; //to make compile happy
}


bool ResonanceList::has_residue( core::Size resid ) const {
	return by_resid_.find( resid ) != by_resid_.end();
}

std::string label_atom_name( std::string const& proton_name, core::chemical::AA aa ) {
	using namespace core::chemical; //for AA
	std::string name;
	if ( aa == aa_arg ) {
		if ( proton_name == "HE" ) return "NE";
		if ( proton_name.substr(0,2) == "HH" ) return "N"+proton_name.substr(1,2); //HH11 HH12 HH21 HH22
	}
	if ( aa == aa_lys ) {
		if ( proton_name.substr(0,2) == "HZ" ) return "NZ";
	}
	if ( aa == aa_gln || aa == aa_his ) {
		if ( proton_name.substr(0,3) == "HE2" ) return "NE2"; //HE21, HE22
	}
	if ( aa == aa_asn ) {
		if ( proton_name.substr(0,3) == "HD2" ) return "ND2"; //HD21, HD22
	}
	if ( aa == aa_trp ) {
		if ( proton_name == "HE1" ) return "NE1";
	}
	if ( aa == aa_his ) {
		if ( proton_name == "HD1" ) return "ND1";
	}
	if ( proton_name == "H" ) return "N";
	if ( proton_name[ 0 ] == 'Q' && proton_name[ 1 ] == 'Q' ) { //QQX
		name = "C" + proton_name.substr(2,1);
		return name;
	}
	if ( proton_name[ 0 ] == 'Q' ) { //Qxx
		name = proton_name;
		name[ 0 ] ='C';
		return name;
	} if ( proton_name.substr(1,2) == "HB" ) {
		return "CB";
	} if ( proton_name.substr(1,2) == "HA" ) {
		return "CA";
	}
	if ( aa == aa_trp ) {
		if ( proton_name == "HH2" ) return "CH2";
		if ( proton_name == "HZ2" ) return "CZ2";
		if ( proton_name == "HZ3" ) return "CZ3";
		if ( proton_name == "HE3" ) return "CE3";
		if ( proton_name == "HD1" ) return "CD1";
	}
	if ( aa == aa_phe || aa == aa_tyr || aa == aa_his ) {
		if ( proton_name == "HZ" || proton_name == "HH" ) return "CZ";
		if ( proton_name.substr(0,2) == "HD" || proton_name.substr(0,2) == "HE"  ) return "C"+proton_name.substr(1,2);
	}

	if ( proton_name.substr(0,2) == "HA" ) return "CA";
	if ( proton_name.substr(0,2) == "HB" ) return "CB";
	if ( aa != aa_asn ) { //don't make HD21 -> CD substition for ASN
		Size len=proton_name.size()-2;
		if ( proton_name.substr(0,2) == "HG" || proton_name.substr(0,2)=="HD" || proton_name.substr(0,2)=="HE" ) return "C"+proton_name.substr(1,len < 1 ? 1 : len );
	}
	if ( tr.Trace.visible() ) {
		tr.Trace << "no label for proton_name " + proton_name + " on " + name_from_aa( aa ) << std::endl;
	}
	throw EXCN_UnknownAtomname("");
	return "no_atom";
}

void ResonanceList::update_bond_connections() {
	std::set< core::id::NamedAtomID > unknown_resonances_;
	for ( auto const & it : *this ) {
		it.second->clear_connected_resonances();
	}
	for ( auto it = begin(); it != end(); ++it ) {
		if ( !it->second->is_proton() ) continue; //do this from the proton
		ResonanceOP proton( it->second );
		Size resid( it->second->atom().rsd() );
		try {
			id::NamedAtomID atomID( label_atom_name( proton->atom().atom(), aa_from_resid( resid ) ), resid );
			Resonance const& label_reso( (*this)[ atomID ] ); //use this operator to have full-checking
			if ( tr_labels.Trace.visible() ) {
				tr_labels.Trace << it->second->atom() << " has label " << atomID << std::endl;
			}
			ResonanceOP label_reso_ptr( map_[ label_reso.label() ] );
			proton->add_connected_resonance( ResonanceAP( label_reso_ptr ) );
			label_reso_ptr->add_connected_resonance( proton );
		} catch ( EXCN_UnknownResonance& exception ) {
			if ( !unknown_resonances_.count( exception.atom() ) ) {
				unknown_resonances_.insert( exception.atom() );
				exception.show( tr.Warning );
				tr.Warning << " as label for atom " << it->second->atom().atom() << " " <<  aa_from_resid( resid ) << std::endl;
			}
			continue;
		} catch ( EXCN_UnknownAtomname& exception ) { //this happens if we try to assign a proton that can't have a label: i.e., a H in a HCH spectrum
			utility_exit_with_message( "cannot find label atom for resid: " + it->second->atom().atom() + " " + ObjexxFCL::string_of( resid ) );
		}
	}
}

}
}
