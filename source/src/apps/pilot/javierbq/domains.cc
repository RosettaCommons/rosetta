// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/dssp/Dssp.hh>


// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/pose/PDBInfo.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/util.hh>
#include <protocols/fldsgn/topology/Sheet.hh>
#include <string>
#include <algorithm>
#include <boost/xpressive/xpressive.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <set>


#include <sstream>

#include <utility/excn/Exceptions.hh>

template <typename T>
std::string vec2str(const T& vec){
		std::ostringstream out;
		for(typename T::const_iterator it = vec.begin(); it != vec.end();  ++it) {
			out << *it << "\t";
		}
		return out.str();
};

using namespace core;
using namespace core::pose;
using namespace utility;
using namespace basic::options;
using namespace basic::options::OptionKeys;
namespace boostx = boost::xpressive;
namespace topology = protocols::fldsgn::topology;

static thread_local basic::Tracer TR( "domain_test" );


class Domain : public utility::pointer::ReferenceCount {
		public:
				Domain();
				Domain(const std::string& name, const PoseOP & pose, Size start, Size end) :
						name_(name),
						pose_(pose),
						start_(start),
						end_(end)		{ 		}


		private:
				PoseOP pose_;
				Size start_;
				Size end_;
				std::string name_;
};

typedef utility::pointer::owning_ptr< Domain > DomainOP;

typedef boost::tuple<Size, Size> DomainMatch;
class DomainDescription : public utility::pointer::ReferenceCount {

		// the three vectors on the tuple represent strand pairs of the sheet ordered from
		// one edge of the sheet to the oposite edge.
		// _1 = strand number vector
		// _2 = orientation vector( 'P' for parallel, 'A' for antiparallel)
		// _3 = comparison vector, obtained by applying less_than_sequencial_comparison to
		// 		  the strand number vector.
		typedef boost::tuple< vector1<Size>, vector1<char>, vector1<bool> > SheetDescription;

		public:
			DomainDescription(const std::string& name, const std::string& sheet_descriptions,
							const std::string& regex_pattern, bool barrel = false)
				 : name_(name),
				 regex_pattern_(regex_pattern),
				 barrel_(barrel)
		 	{
				 regex_expression_ = boostx::sregex::compile(regex_pattern_);
				 ss_info_ = new topology::SS_Info2;
				 boostx::sregex sheets_delim = boostx::as_xpr(';');
				 boostx::sregex strand_pair_delim = boostx::as_xpr(',');
				 boostx::sregex strand_pair_regex = boostx::sregex::compile("(\\d+)-(\\d+)\\.(\\w)");
				 boostx::smatch what;
				 //split the sheet descriptions by sheet
				 boostx::sregex_token_iterator end;
				 for( boostx::sregex_token_iterator sheet_token( sheet_descriptions.begin(), sheet_descriptions.end(), sheets_delim, -1 ) ;sheet_token != end  ; ++sheet_token) {
						 std::string sheet = *sheet_token;
						 vector1<Size> sheet_strand_order;
						 vector1<char> sheet_strand_orientation;
				 		 for( boostx::sregex_token_iterator strand_pair_token( sheet.begin(), sheet.end(), strand_pair_delim, -1 ) ;strand_pair_token != end  ; ++strand_pair_token) {
								 std::string token = *strand_pair_token;
								 bool found = boostx::regex_match(token, what, strand_pair_regex);
								 //TR <<  what[1] << "\t" << what[2] << "\t" << what[3] << std::endl;
								 Size s1 = boost::lexical_cast<Size>(what[1]);
								 Size s2 = boost::lexical_cast<Size>(what[2]);
								 char orientation = what.str(3).at(0);
								 if(sheet_strand_order.empty()) {
										 sheet_strand_order.push_back(s1);
										 sheet_strand_order.push_back(s2);
										 sheet_strand_orientation.push_back(orientation);

								 } else {
										 assert(s1 == sheet_strand_order.back());
										 sheet_strand_order.push_back(s2);
										 sheet_strand_orientation.push_back(orientation);
								 }
						 }
						 sheet_descriptions_.push_back(boost::make_tuple(sheet_strand_order, sheet_strand_orientation, less_than_sequencial_comparison(sheet_strand_order)));
				 }
			}
		  // @brief return true if the pose secondary structure is compatible with the  domain.
			bool
		  ss_match(const Pose& p, vector1<DomainMatch>& ss_matches){
					std::string secstruct = p.secstruct();
					boostx::sregex_iterator cur(secstruct.begin(), secstruct.end(), regex_expression_);
					boostx::sregex_iterator end;
					bool found(cur != end);
					if(found) {
							TR << name_ << " compatible secondary structure found in " << p.pdb_info()->name() << std::endl;
							TR << "#\tstart\tlength\tsecstruct" << std::endl;
							for(Size m = 1; cur != end; ++cur, ++m)
							{
								boostx::smatch const &what = *cur;
								TR << m << "\t" << what.position(0) + 1 << "\t" << what.length() << "\t" << what[0] << std::endl;
								ss_matches.push_back(boost::make_tuple(what.position(0) + 1, what.length()));
							}
					}
					return found;
			}

			vector1<DomainOP>
			sheet_match(const PoseOP& p, vector1<DomainMatch>& ss_matches){
					vector1<DomainOP> results;
					Size nmatch = 1;
					for(vector1<DomainMatch>::iterator it = ss_matches.begin(); it != ss_matches.end(); ++it) {
							Size start  = boost::get<0>(*it);
							Size length = boost::get<1>(*it);
							Pose subpose(*p, start, start + length - 1);
							for(Size i = start; i <= start + length - 1; ++i) {
								subpose.set_secstruct(i, p->secstruct().at(i));
							}
							TR << "match "<< nmatch << ": " << subpose.secstruct() << std::endl;;
							ss_info_->initialize(subpose, subpose.secstruct());
							topology::StrandPairingSetOP spairset = new topology::StrandPairingSet( topology::calc_strand_pairing_set( subpose, ss_info_ ) );
						  topology::SheetSet sheetset(ss_info_, spairset);
							topology::Sheets sheets = sheetset.sheets();
							std::set<Size> found_sheets;
						  for(topology::Sheets::iterator it = sheets.begin(); it != sheets.end(); ++it) {
									Size sheet_des_index = compatible_sheet(*it, found_sheets);
									if(sheet_des_index != 0) {
										TR << "Compatible sheet " << sheet_des_index << std::endl;
											found_sheets.insert(sheet_des_index);

									}
							}
							// if all the sheets are present then domain matches the sheet description.
							if(found_sheets.size() == sheet_descriptions_.size())
									results.push_back( new Domain(name_,p, start, start + length -1) );
					}
					return results;
			}

	  //methods
	  private:
			/// @brief return the index of a compatible SheetDescription for a sheet, 0 if non of the
			/// sheets are compatible.
			Size compatible_sheet(const topology::SheetOP& sheet, std::set<Size> exclude_sheets) {
					for(Size i = 1; i <= sheet_descriptions_.size(); ++i){
							if(exclude_sheets.count(i))
									continue;
							SheetDescription& desc = sheet_descriptions_.at(i);
							//TR << "ref:\t\t" << vec2str(boost::get<2>(desc))  << std::endl;
							// check the strands in the default order
							vector1<bool> sheet_seq_comp = less_than_sequencial_comparison(sheet->order_strands());
							//TR << "sheet_order:\t" << vec2str(sheet_seq_comp) << std::endl;
							// ADD STRAND ORITATION
							if( boost::get<2>(desc) == sheet_seq_comp) {
									return i;
							} else { // try inverting the input sheet.
						  		std::reverse(sheet_seq_comp.begin(), sheet_seq_comp.end());
									// also the values of the compared bools have to be inverted since
									// it is a less than comparison.
									for(vector1<bool>::iterator it = sheet_seq_comp.begin(); it != sheet_seq_comp.end(); ++it)
											*it = !(*it);
							//TR << "sheet_reversed:\t" << vec2str(sheet_seq_comp) << std::endl;
							// ADD STRAND ORITATION
									if(boost::get<2>(desc) == sheet_seq_comp)
											return i;
							}

					}
					return 0;
			}

		vector1<bool> less_than_sequencial_comparison(const vector1<Size>& elements) {
  				vector1<bool> results;
  				Size prev = 0;
  				for(vector1<Size>::const_iterator it = elements.begin(); it != elements.end(); ++it) {
  						if(prev == 0) {
  								prev = *it;
  								continue;
  						} else {
  								results.push_back( prev < *it);
  								prev = *it;
  						}
  				}
  				return results;
		}

	  //data
		private:
			std::string name_;
			std::string regex_pattern_;
			boostx::sregex regex_expression_;
			vector1<SheetDescription> sheet_descriptions_;
			topology::SS_Info2_OP ss_info_;
			bool barrel_;
			Size	nstrands_;
			Size  nhelices_;
			Size  nloops_;
};

typedef utility::pointer::owning_ptr< DomainDescription > DomainDescriptionOP;

class DomainFinder : public utility::pointer::ReferenceCount {
	public:
			DomainFinder();
			void add_domain(const DomainDescription& domain);
			utility::vector1<Domain>	find_domains(const Pose& pose);
	private:
			std::set<DomainDescriptionOP> domains_;
};

class ThisApplication {
    public:
        ThisApplication(){};
        static void register_options();
		};

		OPT_KEY( String , pdb)

		void ThisApplication::register_options() {
				NEW_OPT( pdb, "pdb file", "");
		}


int
main ( int argc, char* argv[] ){
	try {

		ThisApplication::register_options();
		devel::init(argc, argv);
		TR << "Reaading File " << option[pdb]() << std::endl;
		core::pose::PoseOP p;
		import_pose::pose_from_pdb(*p, option[pdb]());
		core::scoring::dssp::Dssp dssp(*p);
		dssp.insert_ss_into_pose(*p);
		DomainDescription dd("Ploop2x3", "5-4.P,4-1.P,1-3.P,3-2.P", "E+L+H+L+E+L+H+L+E+L+H+L+E+L+H+L+E+L+H+");
		TR << "dssp: " << p->secstruct() << std::endl;
		vector1<DomainMatch> ss_matches;
		TR << "match: " << dd.ss_match(*p, ss_matches) << std::endl;
		dd.sheet_match(p, ss_matches);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

    return 0;
}
