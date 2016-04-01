// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/SS_Killhairpins_Info.hh
/// @brief  Scoring manager class header
/// @author Robert Vernon (rvernon@u.washington.edu)


/// Unit headers
#include <core/scoring/SS_Killhairpins_Info.hh>

/// Package headers

#include <ObjexxFCL/format.hh>

// Numeric Headers

// utility headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <string>
#include <utility/exit.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>


/// C++ Headers
#include <iostream>

#include <utility/vector1.hh>

//Auto Headers
#include <numeric/random/random.fwd.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/utility.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

static THREAD_LOCAL basic::Tracer trKillHairpinsIO( "core.score.SS_Killhairpins_Info" );

//////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief default constructor
Hairpin::Hairpin()
{}

Hairpin::Hairpin( core::Size s1_1, core::Size s1_2, core::Size s2_1, core::Size s2_2)
{
	range_pair_ = std::make_pair( std::make_pair( s1_1, s1_2 ), std::make_pair( s2_1, s2_2 ) );
}

Hairpin::~Hairpin()
{}

core::Size
Hairpin::s1_start() const { return range_pair_.first.first; }

core::Size
Hairpin::s1_end() const { return range_pair_.first.second; }

core::Size
Hairpin::s2_start() const { return range_pair_.second.first; }

core::Size
Hairpin::s2_end() const { return range_pair_.second.second; }

/// @brief copy assignment
Hairpin &
Hairpin::operator =( Hairpin const & s )
{
	if ( this != &s ) {
		range_pair_ = s.range_pair_;
	}
	return *this;
}


std::ostream &
operator<< ( std::ostream & out, Hairpin const & s )
{
	using ObjexxFCL::format::I;
	out << "Hairpin - Strand1: " << I(5,s.s1_start()) << "-" << I(5,s.s1_end()) << "  Strand2: "
		<< I(5,s.s2_start()) << "-" << I(5,s.s2_end()) << "\n";
	return out;
}


/// @brief default constructor
Hairpins::Hairpins()
{}

/// @brief copy constructor
Hairpins::Hairpins(
	Hairpins const & s
) : hairpin_list_( s.hairpin_list_ )
{}

/// @brief default destructor
Hairpins::~Hairpins()
{}

/// @brief copy assignment
Hairpins &
Hairpins::operator =( Hairpins const & s )
{
	if ( this != &s ) {
		hairpin_list_ = s.hairpin_list_;
	}
	return *this;
}

void
Hairpins::append_hairpin( core::Size s1_1, core::Size s1_2, core::Size s2_1, core::Size s2_2)
{
	Hairpin new_hairpin(s1_1, s1_2, s2_1, s2_2);

	hairpin_list_.push_back(new_hairpin);
}

void
Hairpins::clear()
{
	hairpin_list_.clear();
}

utility::vector1< Hairpin >
Hairpins::list() const
{
	return hairpin_list_;
}

core::Size
Hairpins::size() const
{
	return hairpin_list_.size();
}

std::ostream &
operator<< ( std::ostream & out, Hairpins const & s )
{
	using ObjexxFCL::format::I;

	utility::vector1< Hairpin > temp_list( s.list() );
	core::Size count(0);

	for ( utility::vector1< Hairpin >::iterator it= temp_list.begin(),
			ite= temp_list.end(); it != ite; ++it ) {
		count++;
	}
	return out;
}

SS_Killhairpins_Info::SS_Killhairpins_Info() :
	CacheableData(),
	kill_parallel_( false ),
	kill_antiparallel_( true )
{}


SS_Killhairpins_Info::SS_Killhairpins_Info( SS_Killhairpins_Info const & src ) :
	CacheableData(),
	kill_parallel_( src.kill_parallel_ ),
	kill_antiparallel_( src.kill_antiparallel_ ),
	hairpins_ ( src.hairpins_ )
{}

bool
SS_Killhairpins_Info::check_hairpin( core::Size const & strand1_res, core::Size const & strand2_res )
{
	for ( core::Size i = 1; i <= hairpins_.list().size(); ++i ) {
		if ( (strand1_res >= hairpins_.list()[i].s1_start()) && (strand1_res <= hairpins_.list()[i].s1_end())
				&& (strand2_res >= hairpins_.list()[i].s2_start()) && (strand2_res <= hairpins_.list()[i].s2_end()) ) {
			return true;
		}
	}
	return false;
}

void
SS_Killhairpins_Info::setup_from_psipred(utility::io::izstream & input_file)
{
	utility::vector1< Hairpin > all_hairpins;

	std::string line;

	utility::vector1< core::Real > strand_prob;
	utility::vector1< core::Real > helix_prob;
	utility::vector1< core::Real > loop_prob;

	while ( !input_file.eof() ) {
		getline(input_file, line);
		if ( line.length() > 5 ) {
			std::istringstream line_stream(line);

			char aa, ss;
			core::Real E_freq, L_freq, H_freq;
			core::Size res;

			line_stream >> res >> aa >> ss >> L_freq >> H_freq >> E_freq;

			strand_prob.push_back(E_freq);
			helix_prob.push_back(H_freq);
			loop_prob.push_back(L_freq);

			if ( (strand_prob.size() != res) || (strand_prob.size() != helix_prob.size())
					|| (loop_prob.size() != helix_prob.size()) ) {
				utility_exit_with_message("[ERROR] -kill_hairpins can't parse psipred_ss2 file");
			}
		}
	}

	///////////////////
	// The following picks out potential hairpins from the psipred_ss2 file
	// Rules:
	// -Strands start when the E_pred jumps above 0.5, and end when it drops below 0.5
	// -Strands must be at least 3 residues long
	// -Adjacent strands are considered hairpin candidates if the H_pred values between
	//  them sum up to less than 3.0 (eg: < 3 100% helix residues or < 6 50% helix)
	//
	// Also, this nest of if statements made sense at the time.
	// -rmv
	///////////////////

	core::Real helix_prob_sum(0.0);
	bool has_helix(false), on_strand(false), has_strand(false);
	std::pair< core::Size, core::Size > working_strand;
	std::pair< core::Size, core::Size > last_strand;
	for ( core::Size i=1; i <= strand_prob.size(); ++i ) {

		if ( !on_strand ) {
			if ( i <= strand_prob.size() - 2 ) {
				if ( ( strand_prob[i] >= 0.5 ) &&
						( strand_prob[i+1] >= 0.5 ) &&
						( strand_prob[i+2] >= 0.5) ) {

					on_strand = true;
					working_strand.first = i;

					if ( kill_parallel_ ) {
						has_strand = true;
					} else {
						if ( !has_helix ) {
							has_strand = true;
						} else {
							has_strand = false;
						}
					}

					has_helix = false;
					helix_prob_sum = 0.0;

				} else {
					helix_prob_sum += helix_prob[i];

					if ( helix_prob_sum >= 3.0 ) {
						has_helix = true;
					}
				}
			}
		} else {

			if ( ( strand_prob[i] < 0.5 ) || ( i == strand_prob.size()) ) {
				working_strand.second = i-1;

				if ( ( has_strand == true ) && ( last_strand.first != last_strand.second ) ) {
					runtime_assert( last_strand != working_strand );

					runtime_assert(   (last_strand.first<=last_strand.second)
						&& (last_strand.second<=working_strand.first)
						&& (working_strand.first<=working_strand.second));

					Hairpin new_hairpin( last_strand.first, last_strand.second,
						working_strand.first, working_strand.second );
					all_hairpins.push_back( new_hairpin );
				}

				last_strand = working_strand;
				working_strand.first = 0;
				working_strand.second = 0;
				on_strand = false;
				has_strand = true;
			}
		}
	}

	core::Real freq_limit( basic::options::option[ basic::options::OptionKeys::abinitio::kill_hairpins_frequency ] );

	for ( core::Size i=1; i<= all_hairpins.size(); ++i ) {

		core::Real freq_attempt(numeric::random::uniform());

		trKillHairpinsIO.Info << "KHP: Hairpin Detected in Psipred File: " << all_hairpins[i].s1_start() << "-" << all_hairpins[i].s1_end() << " " << all_hairpins[i].s2_start() << "-" << all_hairpins[i].s2_end();
		if ( freq_attempt <= freq_limit ) {
			hairpins_.append_hairpin(all_hairpins[i].s1_start(),
				all_hairpins[i].s1_end(),
				all_hairpins[i].s2_start(),
				all_hairpins[i].s2_end());

			trKillHairpinsIO.Info << " KILLED!";
		}
		trKillHairpinsIO.Info << std::endl;
	}
}

void
SS_Killhairpins_Info::setup_from_kill_hairpins_file(utility::io::izstream & input_file)
{
	std::string line;

	while ( !input_file.eof() ) {
		getline(input_file, line);

		if ( line.length() <= 5 ) continue;
		
		std::istringstream line_stream(line);
		
		core::Real frequency;
		core::Size res1, res2, res3, res4;
		
		line_stream >> frequency >> res1 >> res2 >> res3 >> res4;
		
		runtime_assert((res1<=res2) && (res2<=res3) && (res3<=res4));
		
		core::Real freq_attempt(numeric::random::uniform());
		
		trKillHairpinsIO.Info << "KHP: Hairpin Read from Kill Hairpins File: " << res1 << "-" << res2 << " " << res3 << "-" << res4;
		
		if ( freq_attempt <= frequency ) {
			hairpins_.append_hairpin( res1, res2, res3, res4 );
			trKillHairpinsIO.Info << " KILLED!";
		}
		trKillHairpinsIO.Info << std::endl;
	}
}

void
SS_Killhairpins_Info::setup_killhairpins()
{
	hairpins_.clear();

	if ( ! basic::options::option[ basic::options::OptionKeys::abinitio::kill_hairpins ].user() ) return;
	
	utility::io::izstream input_file( basic::options::option[ basic::options::OptionKeys::abinitio::kill_hairpins ] );
	
	if ( !input_file ) {
		utility_exit_with_message("[ERROR] Unable to open kill_hairpins file");
	}
	
	std::string line, a(""), b("");
	getline(input_file, line);
	utility::vector1< std::string > tokens ( utility::split( line ) );
	if ( tokens.size() == 4 ) {
		a = tokens[ 3 ];
		b = tokens[ 4 ];
	} else if ( tokens.size() == 3 ) {
		a = tokens[ 3 ];
	} else if ( tokens.size() == 2 ) {}
	else {
		utility_exit_with_message("[ERROR] invalid header input for kill_hairpins file. ");
	}
	
	if ( a != "ANTI" && a != "PARA" && a != "" ) {
		utility_exit_with_message("[ERROR] invalid kill type for kill_hairpins. ANTI or PARA, is required. ");
	}
	if ( b != "ANTI" && b != "PARA" && b != "" ) {
		utility_exit_with_message("[ERROR] invalid kill type for kill_hairpins. ANTI or PARA, is required. ");
	}
	
	if ( a == "ANTI" || b == "ANTI" ) {
		kill_antiparallel_ = true; // as default, kill_parallel_ is false
	}
	if ( a == "PARA" || b == "PARA" ) {
		kill_parallel_ = true;  // as default, kill_antiparallel_ is true
	}
	
	std::string keyword = tokens[ 2 ];
	if ( keyword == "PSIPRED" ) {
		setup_from_psipred(input_file);
	} else {
		if ( keyword == "KILL" ) {
			setup_from_kill_hairpins_file(input_file);
		} else {
			utility_exit_with_message("[ERROR] Unknown file type for kill_hairpins file (need psipred_ss2 or kill).  File type is autodetected by file header; check your header.");
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////


} // ns scoring
} // ns core



#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::Hairpin::save( Archive & arc ) const {
	arc( CEREAL_NVP( range_pair_ ) ); // std::pair<std::pair<core::Size, core::Size>, std::pair<core::Size, core::Size> >
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::Hairpin::load( Archive & arc ) {
	arc( range_pair_ ); // std::pair<std::pair<core::Size, core::Size>, std::pair<core::Size, core::Size> >
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::Hairpin );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::Hairpins::save( Archive & arc ) const {
	arc( CEREAL_NVP( hairpin_list_ ) ); // utility::vector1<Hairpin>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::Hairpins::load( Archive & arc ) {
	arc( hairpin_list_ ); // utility::vector1<Hairpin>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::Hairpins );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::SS_Killhairpins_Info::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( kill_parallel_ ) ); // _Bool
	arc( CEREAL_NVP( kill_antiparallel_ ) ); // _Bool
	arc( CEREAL_NVP( hairpins_ ) ); // struct core::scoring::Hairpins
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::SS_Killhairpins_Info::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( kill_parallel_ ); // _Bool
	arc( kill_antiparallel_ ); // _Bool
	arc( hairpins_ ); // struct core::scoring::Hairpins
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::SS_Killhairpins_Info );
CEREAL_REGISTER_TYPE( core::scoring::SS_Killhairpins_Info )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_SS_Killhairpins_Info )
#endif // SERIALIZATION

