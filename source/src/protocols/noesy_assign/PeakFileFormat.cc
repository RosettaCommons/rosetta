// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @details
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/PeakFileFormat.hh>
#include <protocols/noesy_assign/CrossPeak.hh>
#include <protocols/noesy_assign/CrossPeakInfo.hh>
#include <protocols/noesy_assign/ResonanceList.hh>
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>

// Package Headers
#include <protocols/noesy_assign/Exceptions.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/AA.hh>

// Utility headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/string_util.hh>
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

#include <basic/options/option_macros.hh>

#include <utility/vector1.hh>


OPT_1GRP_KEY( RealVector, noesy_weights, tolerances )

bool protocols::noesy_assign::PeakFileFormat::options_registered_( false );

void protocols::noesy_assign::PeakFileFormat::register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  if ( options_registered_ ) return;
	NEW_OPT3( noesy_weights::tolerances, "if no #TOLERANCE is found use this for [indirect_H, direct_H, label ]", 0.04, 0.03, 0.3 );
	options_registered_ = true;
}


static thread_local basic::Tracer tr( "protocols.noesy_assign.io" );

namespace protocols {
namespace noesy_assign {

using namespace core;
PeakFileFormat::PeakFileFormat()
{
	runtime_assert( options_registered_ );
}

PeakFileFormat::~PeakFileFormat() {}

Size const TOL_H_INDIRECT( 1 );
Size const TOL_H_DIRECT( 2 );
Size const TOL_LABEL( 3 );

void PeakFileFormat::write_peak( std::ostream& os, Size ct, CrossPeak const& cp ) const {
	std::ostringstream line_end;
	using namespace ObjexxFCL::format;
	line_end << " | ";
	line_end << F( 5, 2, cp.cumulative_peak_volume() ) << " "
					 << F( 5, 2, cp.probability() ) << " "
					 << I( 3, (int)cp.quality_class() ) << " "
		       << A( 15, cp.quality_class_str() ) << " "
					 << F( 5, 2, cp.smallest_native_violation() ) << " ";

	line_end << " #d " << cp.distance_bound();
	if ( cp.eliminated( false /*recompute*/, true /*do_not_compute*/) ) line_end << " #eliminated: " << cp.elimination_reason();

  os << ObjexxFCL::format::RJ( 6, ct ) << " ";
  // cp.write_to_stream( os );
  write_resonances( os, cp );
	os << " 1 U ";
  write_strength( os, cp );
	os << " e 0 ";
  write_assignments( os, cp, line_end.str() );
}

void PeakFileFormat::write_header( std::ostream& os ) {
  Size const dim( col2proton_.size() );
  utility::vector1< std::string > atom_names;
  utility::vector1< core::Real > tolerances;
  //dimension 1 - label
  if ( info1_->has_label() ) {
    atom_names.push_back( ObjexxFCL::lowercased( info1_->label_atom_type() ) );
    tolerances.push_back( info1_->label_tolerance() );
  }

  if ( info2_->has_label() ) { //dimension 2 - label
    atom_names.push_back( ObjexxFCL::uppercased( info2_->label_atom_type() ) );
    tolerances.push_back( info2_->label_tolerance() );
  }

  //dimension 2
  atom_names.push_back( ObjexxFCL::uppercased( info2_->main_atom() ) );
  tolerances.push_back( info2_->proton_tolerance() );

	//dimension 1
  atom_names.push_back( ObjexxFCL::lowercased( info1_->main_atom() ) );
  tolerances.push_back( info1_->proton_tolerance() );

  os << "# Number of dimensions " << dim << std::endl;
	os << "#FILENAME " << filename() << std::endl;
	std::string format_str = "xeasy" + utility::to_string( dim ) + "D";
  os << "#FORMAT " << format_str << std::endl;
  std::string cyana_str;
  for ( Size ct = 1; ct <= atom_names.size(); ct++ ) {
    os << "#INAME " << ct << " " << atom_names[ ct ] << std::endl;
		if ( atom_names[ ct ].size()>1 ) {
			cyana_str += "["+atom_names[ct]+"]";
		}	else {
			cyana_str += atom_names[ ct ];
		}
  }
	if ( info1_->fold_proton_resonance().is_folded() ) {
		os << "#FOLD "<< dim << " " << info1_->fold_proton_resonance().start() << " " << info1_->fold_proton_resonance().end() << std::endl;
	}
	if ( info1_->fold_label_resonance().is_folded() ) {
		os << "#FOLD "<< 1 << " " << info1_->fold_label_resonance().start() << " " << info1_->fold_label_resonance().end() << std::endl;
	}

	if ( info2_->fold_proton_resonance().is_folded() ) {
		os << "#FOLD "<< dim-1 << " " << info2_->fold_proton_resonance().start() << " " << info2_->fold_proton_resonance().end() << std::endl;
	}
	if ( info2_->fold_label_resonance().is_folded() ) {
		os << "#FOLD "<< 2 << " " << info2_->fold_label_resonance().start() << " " << info2_->fold_label_resonance().end() << std::endl;
	}
	os << "#MAX_NOE_DIST " << info1_->max_noe_distance() << std::endl;
  os << "#CYANAFORMAT " << cyana_str << std::endl;
  os << "#TOLERANCE ";
  for ( Size ct = 1; ct <= tolerances.size(); ct++ ) {
    os << ObjexxFCL::format::RJ( 8, tolerances[ ct ] );
  }
  os << std::endl;
}

bool PeakFileFormat::compatible_with_current_format( CrossPeak const& cp ) const {
  return ( *info1_ == cp.info( 1 ) && *info2_ == cp.info( 2 ) );
}

void PeakFileFormat::set_format_from_peak( CrossPeak const& cp ) {
  info1_ = CrossPeakInfoOP( new CrossPeakInfo( cp.info( 1 ) ) );
  info2_ = CrossPeakInfoOP( new CrossPeakInfo( cp.info( 2 ) ) );
  col2proton_.clear();
  col2islabel_.clear();

	//dimension 1 - label
  if ( info1_->has_label() ) {
    col2proton_.push_back( 1 );
    col2islabel_.push_back( true );
  }

	//dimension 2 - label
  if ( info2_->has_label() ) {
    col2proton_.push_back( 2 );
    col2islabel_.push_back( true );
  }

  //dimension 2
  col2proton_.push_back( 2 );
  col2islabel_.push_back( false );

  //dimension 1
  col2proton_.push_back( 1 );
  col2islabel_.push_back( false );
	set_filename( info1_->filename() );

}


void PeakFileFormat::read_header( std::istream& is, std::string& next_line ) {
	col2proton_.clear();
	col2islabel_.clear();
	column_labels_.clear();
	info1_ = info2_ = NULL;
	using namespace ObjexxFCL;
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
  Size dim( 0 );
  utility::vector1< std::string > atom_names;
  utility::vector1< core::Real > tolerances;
	utility::vector1< core::Real > fold_starts( 4, 0);
	utility::vector1< core::Real > fold_ends( 4, 0);

	Real max_noe_dist( params.calibration_max_noe_dist_ );

	bool HN_column_labels( false ); //true if we find a HC or HN ( instead of h vs H )
  std::string line;
	std::string cyana_string("none");
	bool simnoesy( false );
  while ( true ) {
    std::istringstream line_stream;
		if ( next_line.size() ) {
			line=next_line;
			next_line = "";
		} else {
// 			char c='#';
// 			is.get(c);
// 			is.unget();
// 			//		std::cout << "c: " << c << std::endl;
// 			if ( c!='#' ) break;
			getline( is, line );
		 	//		std::cout << " line: " << line << std::endl;
		}
		tr.Trace << "header line: " << line << std::endl;
		line_stream.str( line );
    std::string tag;
    line_stream >> tag;
		if ( tag[0]!='#' ) { next_line = line; }

    if ( line.find( "Number of dimensions" ) != std::string::npos ) {
      line_stream >> tag >> tag >> tag >> dim; //read Number of dimensions N
      atom_names.resize( dim, "" );
      tolerances.resize( dim, 0.0 );
		} else if ( tag == "#FILENAME" ) {
			std::string filename;
			line_stream >> filename;
			set_filename( filename );
    } else if ( tag == "#FORMAT" ) {
      std::string format;
      line_stream >> format;
      std::string expected_format = "xeasy" + utility::to_string( dim ) + "D";
      if ( format != expected_format ) {
				tr.Warning << "[WARNING] Format inconsistent: " << expected_format << " found in file: " << format << std::endl;
      }
    } else if ( tag == "#INAME" ) {
      Size index;
      std::string name;
      line_stream >> index;
			if ( !line_stream.good() ) {
				throw utility::excn::EXCN_BadInput(" problem reading peak file, column index and atom-name expected after key-word #INAME ");
			}
			line_stream >> name;
			if ( name == "HN" || name =="HC" || name =="H1" || name =="1H" ) HN_column_labels = true;
			if ( name == "H1" ) name = "H";
			if ( name == "1H" ) name = "HN";
			if ( name == "C13" || name == "13C" ) name = "C";
			if ( name == "N15" || name == "15N" ) name = "N";
			if ( name == "SIM" ) name = "NC";
			if ( name == "sim" ) name = "nc";
			if ( name == "CN" ) name = "NC" ;
			if ( name == "cn" ) name = "nc" ;
      if ( atom_names.size() < index ) {
				tr.Error << "only " << dim << "D  format; but " << index << " index found in #INAME line " << line << std::endl;
				continue;
      }
			if ( name == "nc" || name =="NC" ) { simnoesy = true; }
			if ( atom_names[ index ] != "" ) tr.Warning << "found index" << index << "in two different #INAME lines "<< std::endl;
      atom_names[ index ] = name;
    } else if ( tag == "#CYANAFORMAT" ) {
			line_stream >> cyana_string;
	// 	} else if (tag == "#SIMNOESY") {
// 			simnoesy = true;
		}	else if ( tag == "#FOLD" ) {
			Size fold_dim;
			Real start;
			Real end;
			line_stream >> fold_dim >> start >> end;
			fold_starts[ fold_dim ]=start;
			fold_ends[ fold_dim ]=end;
		} else if ( tag == "#IGNORE_NEGATIVE_INTENSITY" ) {
			set_ignore_negative_intensity();
		} else if ( tag == "#MINIMUM_PEAK_INTENSITY" ) {
			core::Real val;
			line_stream >> val;
			set_minimum_peak_intensity( val );
		} else if ( tag == "#MAX_NOE_DIST" ) {
			if ( max_noe_dist < 0.01 ) {
				tr.Warning << "MAX_NOE_DIST flag in peak-file ignored because of 0.0 in -noesy::calibration::max_noe_dist" << std::endl;
			} else {
				line_stream >> max_noe_dist;
			}
		} else if ( tag == "#TOLERANCE" ) {
      for ( Size i = 1; i <= dim; i++ ) {
				core::Real val;
				line_stream >> val;
				tolerances[ i ] = val;
      }
      break; //last line of header ...
    } else if ( tag[0] == '#' ) {
      tr.Warning << "[PeakFileFormat]: ignore header line: " << line << std::endl;
      continue; //ignore unknown comments
    } else {
      tr.Error << "reading line as header: " << line << std::endl
							 << "expect TOLERANCE as last header entry " << std::endl;
      //ERROR
    }
  } //while

	tr.Debug << "finished with header found " << dim << " dimensions" << std::endl;
	using namespace basic::options;
  using namespace basic::options::OptionKeys;
	Real const default_tolerance_h( (dim < 4) ?
		option[ noesy_weights::tolerances ][ TOL_H_INDIRECT ] :
		option[ noesy_weights::tolerances ][ TOL_H_DIRECT ]
	);
	Real const default_tolerance_H( (dim < 3) ?
		option[ noesy_weights::tolerances ][ TOL_H_INDIRECT ] :
		option[ noesy_weights::tolerances ][ TOL_H_DIRECT ]
	);
	if ( dim < 2 ) {
		throw utility::excn::EXCN_BadInput(" problem reading peak file, no or incomplete header ");
	}

  // 	info2_ = new CrossPeakInfo( "N", 0.03, 0.5 );
  col2proton_.resize( dim );
  col2islabel_.resize( dim, false );
	if ( HN_column_labels && cyana_string != "none" ) {
		tr.Info << "use CYANA string to work out column order" << std::endl;
	}
  for ( Size i = 1; i<=dim; i++ ) {
		if ( HN_column_labels && cyana_string != "none" ) {
			if ( simnoesy ) throw utility::excn::EXCN_BadInput("cannot use HN for protons when SimNOESY ( i.e., NC for label atom name");
			atom_names[ i ]=cyana_string[ i-1 ];
		}
    if ( atom_names[ i ] == "h"
			|| ( HN_column_labels && atom_names[ i ]=="H" )
			|| ( atom_names[ i ]=="H" && i==2 && dim==2 && col2proton_[ 1 ] == 1 )
		) {  //indirect in 3D
      col2proton_[ i ] = 2;
			if ( tolerances[ i ]==0.0 ) tolerances[ i ]=default_tolerance_h;
      if ( !info2_ ) {
				info2_ = CrossPeakInfoOP( new CrossPeakInfo( uppercased( atom_names[ i ] ), "", max_noe_dist, tolerances[ i ], 0.0 ) );
      } else {
				info2_->set_proton( uppercased( atom_names[ i ] ), tolerances[ i ] );
      }
    } else if ( atom_names[ i ] == "H" || ( HN_column_labels && atom_names[ i ]=="HN" ) || ( HN_column_labels && atom_names[ i ]=="HC" )) {
      col2proton_[ i ] = 1;  //direct h in 3D
			if ( HN_column_labels ) atom_names[ i ] = "h";
			if ( tolerances[ i ]==0.0 ) tolerances[ i ]=default_tolerance_H;
      if ( !info1_ ) {
				info1_ = CrossPeakInfoOP( new CrossPeakInfo( atom_names[ i ], "", max_noe_dist, tolerances[ i ], 0.0 ) );
      } else {
				info1_->set_proton( atom_names[ i ], tolerances[ i ] );
      }
    } else if ( atom_names[ i ] == "c" || atom_names[ i ] == "n" || atom_names[ i ] == "nc" ) {
      col2proton_[ i ] = 2;
      col2islabel_[ i ] = true;
			if ( tolerances[ i ]==0.0 ) tolerances[ i ]=	option[ noesy_weights::tolerances ][ TOL_LABEL ];
      if ( !info2_ ) {
				info2_ = CrossPeakInfoOP( new CrossPeakInfo( "", uppercased( atom_names[ i ] ), max_noe_dist, 0.0, tolerances[ i ] ) );
      } else {
				info2_->set_label( uppercased( atom_names[ i ] ), tolerances[ i ] );
      }
    } else if ( atom_names[ i ] == "C" || atom_names[ i ] == "N" || atom_names[ i ] == "NC"  ) {
      col2proton_[ i ] = 1;
      col2islabel_[ i ] = true;
			if ( tolerances[ i ]==0.0 ) tolerances[ i ]=	option[ noesy_weights::tolerances ][ TOL_LABEL ];
      if ( !info1_ ) {
				info1_ = CrossPeakInfoOP( new CrossPeakInfo( "", atom_names[ i ], max_noe_dist, 0.0, tolerances[ i ] ) );
      } else {
				info1_->set_label( atom_names[ i ], tolerances[ i ] );
      }
    }
  }
	for ( Size i = 1; i<=dim; i++ ) {
		CrossPeakInfoOP info = col2proton_[ i ]==1 ? info1_ : info2_;
		info->set_folding_window( fold_starts[ i ], fold_ends[ i ], col2islabel_[ i ] );
	}

	if ( !info2_  || !info1_ ) {
		throw utility::excn::EXCN_BadInput(" problem reading peak file, no or errorenous header ");
	}

	tr.Debug << " cross-peak infos: " << *info1_ << " and " << *info2_ << std::endl;

	if ( col2proton_.size() == 3 ) { //only one label.. make sure that it is with proton 1
		output_diagnosis( tr.Debug );
		tr.Debug << " check if we need to swap columns " << std::endl;
		Size ct_col_1( 0 );
		for ( Size i=1; i<=3; i++ ) {
			if ( col2proton_[ i ] == 1 ) ct_col_1++;
		}
		if ( ct_col_1 == 1 ) {//ok we have label on 2 switch numbering around...
			for ( Size i=1; i<=3; i++ ) {
				col2proton_[ i ] = col2proton_[ i ]==1 ? 2 : 1;
			}
			CrossPeakInfoOP info = info2_;
			info2_=info1_;
			info1_=info;
		}
	} //swapping of columns
	info1_->set_filename( filename() );
	info2_->set_filename( filename() );
}

CrossPeakOP PeakFileFormat::read_peak( std::istream& is, std::string& next_line ) const {

  CrossPeakOP cp;
  Size const ncol( col2proton_.size() );
  //CrossPeak factory
  runtime_assert( ncol >=2 && ncol <= 4 );
  if ( ncol == 2 ) {
    cp = CrossPeakOP( new CrossPeak );
  } else if ( ncol == 3 ) {
    cp = CrossPeakOP( new CrossPeak3D );
  } else {
    cp = CrossPeakOP( new CrossPeak4D );
  }
  cp->set_info( 1, info1_ );
  cp->set_info( 2, info2_ );

	if ( !next_line.size() ) {
		getline( is, next_line );
	}
	tr.Trace << " next_line: " << next_line << std::endl;
	std::istringstream line_stream( next_line );

  core::Size id;
  line_stream >> id;
  if ( !line_stream.good() ) {
    return NULL;
  }
	next_line = ""; //now we are consuming this line
  cp->set_peak_id( id );

  read_resonances( line_stream, *cp );

  //eat two weird columns
  std::string tag;
  line_stream >> tag; line_stream >> tag;

	///precautions for different ways to finish the line: in particular with and without assignments
  read_strength( line_stream, *cp );

  //read e 0 -- two columns
  line_stream >> tag; line_stream >> tag;
	read_assignments( is, line_stream, *cp, next_line );

  return cp;
}

void PeakFileFormat::read_resonances( std::istream& is, CrossPeak &cp ) const {
  Size const ncol( col2proton_.size() );
  for ( Size icol=1; icol<=ncol; ++icol ) {
    Real val;
    Size iproton( col2proton_[ icol ] );
    bool is_label( col2islabel_[ icol ] );
    is >> val;
    if ( !is.good() ) {
      throw EXCN_FileFormat( "expected resonance value" );
    }
    if ( !is_label ) cp.proton( iproton ).set_freq( val );
    else {
      runtime_assert( cp.has_label( iproton ) );
      cp.label( iproton ).set_freq( val );
    }
  }
}


void PeakFileFormat::write_resonances( std::ostream& os, CrossPeak const& cp ) const {
  //  cp.write_to_stream( os );
  Size const ncol( col2proton_.size() );
  runtime_assert( col2islabel_.size() == ncol );

  for ( Size icol=1; icol<=ncol; ++icol ) {
    Real val;
    Size iproton( col2proton_[ icol ] );
    bool is_label( col2islabel_[ icol ] );
    if ( !is_label ) val = cp.proton( iproton ).freq();
    else {
      runtime_assert( cp.has_label( iproton ) );
      val = cp.label( iproton ).freq();
    }
    os << ObjexxFCL::format::F( 8, 3, val ) << " ";
  }
}

void PeakFileFormat::read_strength( std::istream& is, CrossPeak& cp ) const {
  core::Real val;
  is >> val;
	//val = val < 0 ? -val : val;
  cp.set_volume( val );
  is >> val; //read 0.00E+00
}

void PeakFileFormat::write_strength( std::ostream& os, CrossPeak const& cp ) const {
  os << ObjexxFCL::format::E( 10, 3, cp.volume() ) << " " << ObjexxFCL::format::E( 10, 3, 0.0 ) << " ";
}


void PeakFileFormat::read_assignments( std::istream& is, std::istream& rest_is, CrossPeak& cp, std::string& new_peak_line ) const {

  Size const ncol( col2proton_.size() );
  runtime_assert( col2islabel_.size() == ncol );
  core::Real weight( 0.0 );
  std::string line;
	if ( !getline( rest_is, line ) ) return;
	tr.Trace << "rest_of_line: --" << line << "--" << std::endl;
	bool rest( true );
	bool first( true );
  while ( true ) {
		if ( !rest ) {
			if ( !getline( is, line ) ) {
				if ( first ) return;
				std::ostringstream errstr;
				errstr << cp;
				throw EXCN_FileFormat( "expected assignment value for " + errstr.str() );
			}
		}
		rest = false;
		std::istringstream line_stream( line );
		new_peak_line = line;
		{ //test if this next cross-peak
			std::istringstream line_stream_test_new_peak( line );
			Size dummyI;
			Real dummyR;
			std::string tag;
			line_stream_test_new_peak >> dummyI >> dummyR >> dummyR >> dummyR >> dummyI >> tag;
			if ( tag == "U" ) return;
		}
		{ //test -- new header ?
			std::istringstream line_stream_test_new_peak( line );
			std::string tag;
			line_stream_test_new_peak >> tag;
			if ( tag[0]=='#' ) return;
		}
		{ //test -- for a decimal point in third value...
			std::istringstream line_stream_test_new_peak( line );
			std::string tag1, tag2, tag3;
			line_stream_test_new_peak >> tag1 >> tag2 >> tag3;
			if ( tag3.find(".") != std::string::npos ) return;
		}
    //only assign if all spins are assigned.
    Size vals[ 5 ];
    for ( Size icol=1; icol<=ncol; ++icol ) {
      Size val;
      line_stream >> val;
      if ( !line_stream ) {
				if ( first ) {
					new_peak_line = ""; //if first then it is maybe just a space at end of line...
					return;
				}
				std::ostringstream errstr;
				errstr << cp;
				throw EXCN_FileFormat( "expected assignment value for " + errstr.str() );
      } else {
				if ( first ) tr.Trace << " read assignments: ";
			}
			first = false;
			tr.Trace << val << " ";
      if ( val == 0 ) {
				line_stream.setstate( std::ios_base::failbit );
				break;
      }
      vals[ icol ]=val;
    }
		tr.Trace << std::endl;
    Size reorder[ 5 ];//spin1 spin2 label1 label2
    if ( line_stream && !ignore_assignments() ) {
      for ( Size icol=1; icol<=ncol; ++icol ) {
				Size val = vals[ icol ];
				Size iproton( col2proton_[ icol ] );
				bool is_label( col2islabel_[ icol ] );
				Size index = iproton + ( is_label ? 2 : 0 );
				runtime_assert( !is_label || cp.has_label( iproton ) );
				reorder[ index ] = val;
      }
			tr.Trace << " add assignment " << std::endl;
			cp.add_full_assignment( reorder );
    }
    std::string tag;
    line_stream >> tag;
    if ( !line_stream.good() ) {
      weight = 1.0; //no VC -- no partial weights.
      break; //finish while loop
    }
    if ( tag == "#VC" ) {
      core::Real w;
      line_stream >> w;
      weight += w;
      line_stream >> tag;
    } else weight += 1.0;
		if ( !line_stream.good() ) break;
  //   if ( tag != "#QU" ) {
//       throw EXCN_FileFormat( "expected #QU " );
//     }
//     core::Real qu, sup;
//     line_stream >> qu;
//     line_stream >> tag;
//     if ( line_stream.good() && tag != "#SUP" ) {
//       throw EXCN_FileFormat( "expected #SUP" );
//     }
		//   line_stream >> sup;
  }
	new_peak_line = "";
}

void PeakFileFormat::write_assignment_indent( std::ostream& os, CrossPeak const& cp ) const {
	os << std::endl << "                                                                  ";
	if ( cp.has_label( 1 ) && cp.has_label( 2 ) ) os << "         ";
}

void PeakFileFormat::write_assignment( std::ostream& os, PeakAssignment const& pa ) const {
	for ( Size icol=1; icol<=ncol(); ++icol ) {
		Size val;
		Size iproton( col2proton_[ icol ] );
		bool is_label( col2islabel_[ icol ] );
		if ( !is_label ) val = pa.resonance_id( iproton );
		else {
			val = pa.label_resonance_id( iproton );
		}
		if ( write_atom_names_ ) {
			if ( !is_label ) os << pa.atom( iproton ) << "   ";
			else os << pa.label_atom( iproton )  << "   ";
		} else {
			os << ObjexxFCL::format::RJ( 6, val ) << " ";
		}
	}
}

void PeakFileFormat::write_assignment_stats( std::ostream& os, PeakAssignment& pa ) const {
	Real val( pa.normalized_peak_volume() );
	os << "#VC " << ObjexxFCL::format::F( 5, 3, val );
	os << " #W ";
	pa.dump_weights( os );
}

void PeakFileFormat::write_assignments( std::ostream& os, CrossPeak const& cp, std::string const& line_end ) const {
#ifndef WIN32
	Size assignments_written( 0 );

	if ( write_only_highest_VC() ) {
		Real bestVC( -1 );
		CrossPeak::PeakAssignments::const_iterator best_VC_it = cp.assignments().end();
		for ( CrossPeak::PeakAssignments::const_iterator it = cp.assignments().begin(); it != cp.assignments().end(); ++it ) {
			Real val( (*it)->normalized_peak_volume() );
			if ( val > bestVC ) {
				bestVC = val;
				best_VC_it = it;
			}
		}
		if ( best_VC_it != cp.assignments().end() ) {
			write_assignment( os, **best_VC_it );
			write_assignment_stats( os, **best_VC_it );
			++assignments_written;
			os << line_end;
		}
	} else {
		for ( CrossPeak::PeakAssignments::const_iterator it = cp.assignments().begin(); it != cp.assignments().end(); ++it ) {
			Real val( (*it)->normalized_peak_volume() );
			if ( val >= min_VC_to_write() ) {
				++assignments_written;
				if ( assignments_written > 1 ) write_assignment_indent( os, cp );
				write_assignment( os, **it );
				write_assignment_stats( os, **it );
			}
			if ( assignments_written == 1 )	os << line_end;
		}
	}
	if ( assignments_written == 0 ) write_nil_assignment( os );
#endif
}

void PeakFileFormat::output_diagnosis( std::ostream& os ) const {
  os << "columns: ";
  for ( Size i = 1; i <= col2proton_.size(); i++ ) {
    os << col2proton_[ i ] << " ";
  }
  os << "\nlabels: ";
  for ( Size i = 1; i <= col2islabel_.size(); i++ ) {
    os << col2islabel_[ i ] << " ";
  }
  os << "\ninfo1: " << info1_->label_atom_type() << " " << info1_->proton_tolerance() << " " << info1_->label_tolerance();
  os << "\ninfo2: " << info2_->label_atom_type() << " " << info2_->proton_tolerance() << " " << info2_->label_tolerance();
  os << std::endl;
}

}
}
