// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/scs_blast.cc
/// @brief Structural Component Selector (SCS) implementation with NCBI-BLAST+
/// @author Sergey Lyskov

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/scs_blast.hh>
#include <protocols/antibody/grafting/exception.hh>
#include <protocols/antibody/grafting/chothia_numberer.hh>
#include <protocols/antibody/grafting/scs_functor.hh>
#include <protocols/antibody/grafting/util.hh>

#include <basic/Tracer.hh>
#include <basic/execute.hh>

#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <utility/string_util.hh>
#include <utility/stream_util.hh>

#include <algorithm>
#include <fstream>
#include <locale>

namespace protocols {
namespace antibody {
namespace grafting {

using std::string;

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.grafting");



/// @brief Create result set from 'row' element of each SCS_ResultVector vector
///        if strict is true then throw if no resuls found otherwise use empty OP
///
/// @throw std::out_of_range if for some of the row is not present and strict==true
SCS_ResultSet SCS_Results::get_result_set(uint row, bool strict)
{
	if( h1.size()  <= row  and  h2.size()  <= row  and  h3.size() <= row  and
		l1.size()  <= row  and  l2.size()  <= row  and  l3.size() <= row  and
		frh.size() <= row  and  frl.size() <= row  and  orientation.size() <= row  and  strict) throw std::out_of_range("SCS_Results::get_result_set could not find ALL results for row:" + std::to_string(row));

	SCS_ResultSet r;

	r.h1 = h1.size() > row ? h1[row] : SCS_ResultOP();
	r.h2 = h2.size() > row ? h2[row] : SCS_ResultOP();
	r.h3 = h3.size() > row ? h3[row] : SCS_ResultOP();

	r.l1 = l1.size() > row ? l1[row] : SCS_ResultOP();
	r.l2 = l2.size() > row ? l2[row] : SCS_ResultOP();
	r.l3 = l3.size() > row ? l3[row] : SCS_ResultOP();

	r.frh = frh.size() > row ? frh[row] : SCS_ResultOP();
	r.frl = frl.size() > row ? frl[row] : SCS_ResultOP();

	r.orientation = orientation.size() > row ? orientation[row] : SCS_ResultOP();

	return r;
}



void SCS_Base::add_filter(SCS_FunctorCOP filter)
{
	filters_.push_back(filter);
}

void SCS_Base::set_sorter(SCS_FunctorCOP sorter)
{
	sorter_ = sorter;
}

struct ResultItem {
	string name;
	SCS_ResultVector & result;
};


void SCS_Base::report(SCS_ResultsOP r, uint n)
{
	// text report
	for(uint i=0; i<n; ++i) {
		SCS_ResultSet s = r->get_result_set(i);

		*this << "Template №" << i << "\n";

		struct {
			string name;
			SCS_ResultOP SCS_ResultSet::*region;
		} J[] {
			{"h1", &SCS_ResultSet::h1}, {"h2", &SCS_ResultSet::h2}, {"h3", &SCS_ResultSet::h3},
			{"l1", &SCS_ResultSet::l1}, {"l2", &SCS_ResultSet::l2}, {"l3", &SCS_ResultSet::l3},
			{"frh", &SCS_ResultSet::frh}, {"frl", &SCS_ResultSet::frl}, {"orientation", &SCS_ResultSet::orientation},
		};

		for(auto &j : J) {
			if( SCS_ResultOP p = s.*j.region ) {
				SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p.get() );

				*this << "\t" << j.name << "\n";
				*this << "\t\tpdb: " << br->pdb << "\n";
				*this << "\t\tpadded: " << br->padded << "\n";
				*this << "\t\tstructure source: " << br->struct_source << "\n";
				*this << "\t\tlight type iso-type: " << br->light_type << "\n";
				*this << "\t\tbio type: " << br->bio_type << "\n";
				*this << "\n";
			}
			else *this << "\tNO TEMPLATE FOR REGION: " << j.name << " SELECTED!\n";
		}
	}

	// json report
	struct {
		string name;
		SCS_ResultVector SCS_Results::*region;
	} J[] {
		{"h1", &SCS_Results::h1}, {"h2", &SCS_Results::h2}, {"h3", &SCS_Results::h3},
		{"l1", &SCS_Results::l1}, {"l2", &SCS_Results::l2}, {"l3", &SCS_Results::l3},
		{"frh", &SCS_Results::frh}, {"frl", &SCS_Results::frl}, {"orientation", &SCS_Results::orientation},
	};

	utility::json_spirit::Object root;
	for(auto &j : J) {
		utility::json_spirit::Array region;

		for(auto &p : (*r).*j.region) {
			if(p) {
				SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p.get() );

				utility::json_spirit::Object jbr;
				jbr.push_back( utility::json_spirit::Pair("pdb", br->pdb) );
				jbr.push_back( utility::json_spirit::Pair("padded", br->padded) );
				jbr.push_back( utility::json_spirit::Pair("struct_source", br->struct_source) );
				jbr.push_back( utility::json_spirit::Pair("iso-type", br->light_type) );
				jbr.push_back( utility::json_spirit::Pair("bio-type", br->bio_type) );
				jbr.push_back( utility::json_spirit::Pair("resolution", br->resolution) );
				jbr.push_back( utility::json_spirit::Pair("identity", br->identity) );
				jbr.push_back( utility::json_spirit::Pair("bit_score", br->bit_score) );

				region.push_back(jbr);
			}
		}
		root.push_back( utility::json_spirit::Pair(j.name, region) );
	}
	set("scs-results", root);
}


SCS_ResultsOP SCS_Base::select(uint n, AntibodySequence const &A)
{
	SCS_ResultsOP r = raw_select(A);

	for(auto &filter : filters_) filter->apply(A, r);

	if(sorter_) sorter_->apply(A, r);

	pad_results(n, A, *r);

	TR.Debug << "SCS_Base::select:                              Final results count: " << TR.Red << result_sizes(r) << TR.Reset << std::endl;

	report(r, 1);  // by default report only one template and allow multi-template sub-class to report all templates if any

	return r;
}

void SCS_LoopOverSCs::init_from_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	TR << "Setting antibody grafting database path from command line-options to: " << TR.bgWhite << TR.Black << option[OptionKeys::antibody::grafting_database]() << std::endl;
	set_database_path( option[OptionKeys::antibody::grafting_database]() );

}


void populate_results_from_db( SCS_Antibody_Database_ResultOP const & result, std::map< string, std::map<string, string> > const & db )
{
	result->bio_type      = db.at(result->pdb).at("BioType");
	result->light_type    = db.at(result->pdb).at("LightType");
	result->struct_source = db.at(result->pdb).at("StructSource");
	result->resolution    = std::stod( db.at(result->pdb).at("resolution") );

	// include the sequences of the other SCs in the template PDB.
	result->h1 = db.at( result->pdb ).at("h1");
	result->h2 = db.at( result->pdb ).at("h2");
	result->h3 = db.at( result->pdb ).at("h3");
	result->frh = db.at( result->pdb ).at("frh");

	result->l1 = db.at( result->pdb ).at("l1");
	result->l2 = db.at( result->pdb ).at("l2");
	result->l3 = db.at( result->pdb ).at("l3");
	result->frl = db.at( result->pdb ).at("frl");
}

/// @details Parse output of NCBI Blast+
///          Expected input:
///          # BLASTP 2.2.31+
///          # Query: h1.fasta
///          # Database: ./../../tools/antibody/blast_database/database.H1.10
///          # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
///          # 423 hits found
///          h1.fasta	pdb4d9r_chothia.pdb	100.00	10	0	0	1	10	1	10	6e-04	22.7
///          h1.fasta	pdb4d9q_chothia.pdb	100.00	10	0	0	1	10	1	10	6e-04	22.7
///          ...
/// @return utility::vector0< std::map<string field, string value> >

SCS_ResultVector parse_blastp_output(string const & file_name, string const & query, std::map< string, std::map<string, string> > const & db)
{
	SCS_ResultVector results;

	auto lines( parse_plain_text_with_columns(file_name, "# Fields: ", ',', query, '\t') );

	for(auto fields : lines ) {
		// Converting text based results filed to SCS_BlastMetric
		SCS_BlastResultOP r = std::make_shared<SCS_BlastResult>();

		r->pdb = fields["subject-id"].substr(3,4);  // pdb2adf_chothia.pdb → 2adf

		r->alignment_length = std::stoi( fields.at("alignment-length") );

		r->identity  = std::stod( fields.at("%-identity") );
		r->bit_score = std::stod( fields.at("bit-score") );

		populate_results_from_db( r, db );

		/*r->bio_type      = db.at(r->pdb).at("BioType");
		r->light_type    = db.at(r->pdb).at("LightType");
		r->struct_source = db.at(r->pdb).at("StructSource");

		// include the sequences of the other SCs in the template PDB.
		r->h1 = db[r->pdb].at("h1");  r->h2 = db[r->pdb].at("h2");  r->h3 = db[r->pdb].at("h3");  r->frh = db[r->pdb].at("frh");
		r->l1 = db[r->pdb].at("l1");  r->l2 = db[r->pdb].at("l2");  r->l3 = db[r->pdb].at("l3");  r->frl = db[r->pdb].at("frl");

		r->sequence = db[r->pdb][results.sequence];
		//TR << "PDB:" << r->pdb << std::endl; */
		//TR << "pdb:" << r->pdb << " l3: " << r->l3 << "  db.l3:" << db.at( r->pdb ).at("l3") << std::endl;

		results.push_back(r);
	}

	return results;
}


SCS_LoopOverSCs::Antibody_SCS_Database & SCS_LoopOverSCs::antibody_scs_database()
{
	if( antibody_scs_database_.size() == 0 ) antibody_scs_database_ = std::move( parse_plain_text_with_columns(database_path_ + "/info/antibody.info") );
	return antibody_scs_database_;
}


SCS_ResultsOP SCS_LoopOverSCs::raw_select(AntibodySequence const &A) //, AntibodyNumbering const &N)
{
	//AntibodyFramework const frh( A.heavy_framework() );
	//AntibodyFramework const frl( A.light_framework() );

	FRH_FRL fr( calculate_frh_frl(A) );

	string frh = fr.frh1 + fr.frh2 + fr.frh3 + fr.frh4;
	string frl = fr.frl1 + fr.frl2 + fr.frl3 + fr.frl4;

	string heavy_orientation = fr.frh1 + A.h1_sequence() + fr.frh2 + A.h2_sequence() + fr.frh3 + A.h3_sequence() + fr.frh4;
	string light_orientation = fr.frl1 + A.l1_sequence() + fr.frl2 + A.l2_sequence() + fr.frl3 + A.l3_sequence() + fr.frl4;

	string orientation = light_orientation + heavy_orientation; // Note: must be in this order

	TR.Trace << "frh:" << frh << std::endl;
	TR.Trace << "frl:" << frl << std::endl;

	// _framework_names_ = ['FRL', 'FRH', 'light', 'heavy', 'L1', 'L2', 'L3', 'H1', 'H2', 'H3', 'light_heavy']

	SCS_ResultsOP results(new SCS_Results);
	//results->antibody_sequence = A;


	Result J[] {
		{"frh", frh, results->frh}, {"h1", A.h1_sequence(), results->h1}, {"h2", A.h2_sequence(), results->h2}, {"h3", A.h3_sequence(), results->h3},
		{"frl", frl, results->frl}, {"l1", A.l1_sequence(), results->l1}, {"l2", A.l2_sequence(), results->l2}, {"l3", A.l3_sequence(), results->l3},
	    {"orientation", orientation, results->orientation},
	};

	string blast_database( database_path_ + "/blast_database/database" );

	Antibody_SCS_Database & antibody_info_lines( antibody_scs_database() );
	//TR << antibody_info_lines << std::endl;
	std::map< string, std::map<string, string> > antibody_info_db;
	for(auto fields : antibody_info_lines) antibody_info_db[fields["pdb"]] = fields;

	for(auto &j : J) {

		std::stringstream sq;  sq << j.sequence.size();
		string db_suffix = j.name + '.' + sq.str();

		if (j.name.find("fr") < string::npos) { db_suffix = j.name; }
		std::transform(db_suffix.begin(), db_suffix.end(), db_suffix.begin(), toupper);
		if( j.name.find("orientation")  < string::npos ) { db_suffix = "light_heavy"; }

		string db_to_query = blast_database + '.' + db_suffix;

		select_template( j, db_to_query, antibody_info_db );
	}

	return results;
}


void SCS_BlastPlus::init_from_options()
{
	SCS_LoopOverSCs::init_from_options();

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	TR << "Setting output prefix from command line-options to: " << TR.bgWhite << TR.Black << option[OptionKeys::antibody::prefix]() << std::endl;
	set_output_prefix( option[OptionKeys::antibody::prefix]() );

	TR << "Setting path to NCBI-Blast+ executable from command line-options to: " << TR.bgWhite << TR.Black << option[OptionKeys::antibody::blastp]() << std::endl;
	set_blastp_executable( option[OptionKeys::antibody::blastp]() );
}


void SCS_BlastPlus::select_template(
	Result & j,
	string const & db_to_query,
	std::map< string, std::map<string, string> > const & ab_db ) const
{
	string fasta_file_name(prefix_ + j.name + ".fasta");
	string align_file_name(prefix_ + j.name + ".align");

	//TR << "Working on: " << j.name << std::endl;

	std::ofstream f(fasta_file_name);  f << "> " << fasta_file_name << '\n' << j.sequence << '\n';  f.close();

	// wordsize, e_value, matrix = (2, 0.00001, 'BLOSUM62') if k.count('FR') or k.count('heavy') or k.count('light') else (2, 2000, 'PAM30')
	//string extra = (j.name.find("fr") < string::npos  or  j.name.find("heavy") < string::npos  or  j.name.find("light") < string::npos) ? "-evalue 0.00001 -matrix BLOSUM62" : "-evalue 2000 -matrix PAM30";
	string extra = (j.name.find("fr") < string::npos  or  \
					j.name.find("orientation") < string::npos or \
					j.name.find("heavy") < string::npos  or  j.name.find("light") < string::npos) ? "-evalue 0.00001 -matrix BLOSUM62" : "-evalue 2000 -matrix PAM30";

	//string command_line ("cd "+working_dir+" && "+blast+" -db "+blast_database+'.'+db_suffix+" -query "+fasta_file_name+" -out "+align_file_name+" -word_size 2 -outfmt 7 -max_target_seqs 1024 "+extra);
	string command_line (blastp_+" -db "+db_to_query+" -query "+fasta_file_name+" -out "+align_file_name+" -word_size 2 -outfmt 7 -max_target_seqs 1024 "+extra);

	basic::execute("Running blast+ for " + j.name + "...", command_line);

	//j.results.name = j.name;
	//j.results.sequence = j.sequence;
	j.results = parse_blastp_output(align_file_name, fasta_file_name, ab_db);

}

/// Pad results vectors for each region (if possible) by adding arbitraty but compatible templates so at least n templates for each region is avalible
void SCS_BlastPlus::pad_results(uint N, AntibodySequence const &A, SCS_Results &results)
{
	Antibody_SCS_Database & antibody_info_lines( antibody_scs_database() );

	FRH_FRL fr( calculate_frh_frl(A) );

	string frh = fr.frh1 + fr.frh2 + fr.frh3 + fr.frh4;
	string frl = fr.frl1 + fr.frl2 + fr.frl3 + fr.frl4;

	string heavy_orientation = fr.frh1 + A.h1_sequence() + fr.frh2 + A.h2_sequence() + fr.frh3 + A.h3_sequence() + fr.frh4;
	string light_orientation = fr.frl1 + A.l1_sequence() + fr.frl2 + A.l2_sequence() + fr.frl3 + A.l3_sequence() + fr.frl4;

	string orientation = light_orientation + heavy_orientation; // Note: must be in this order

	Result J[] {
		{"frh", frh, results.frh}, {"h1", A.h1_sequence(), results.h1}, {"h2", A.h2_sequence(), results.h2}, {"h3", A.h3_sequence(), results.h3},
		{"frl", frl, results.frl}, {"l1", A.l1_sequence(), results.l1}, {"l2", A.l2_sequence(), results.l2}, {"l3", A.l3_sequence(), results.l3},
		{"orientation", orientation, results.orientation},
	};

	for(auto &j : J) {
		for(auto i = antibody_info_lines.begin(); j.results.size() < N  and  i != antibody_info_lines.end(); ++i) {
			if( j.name == "orientation" || j.sequence.size() == i->at(j.name).size() ) {
				SCS_BlastResultOP r = std::make_shared<SCS_BlastResult>();

				r->padded = true;

				r->pdb = i->at("pdb");

				r->alignment_length = 0;
				r->identity  = 0;
				r->bit_score = 0;

				r->bio_type      = i->at("BioType");
				r->light_type    = i->at("LightType");
				r->struct_source = i->at("StructSource");

				r->h1 = i->at("h1");  r->h2 = i->at("h2");  r->h3 = i->at("h3");  r->frh = i->at("frh");
				r->l1 = i->at("l1");  r->l2 = i->at("l2");  r->l3 = i->at("l3");  r->frl = i->at("frl");

				j.results.push_back(r);

				TR << TR.Red << "Not enought templates for region " << TR.bgRed << TR.Black << j.name << TR.Reset << TR.Red << " was selected... Adding arbitrary template " << TR.bgRed << TR.Black << r->pdb << TR.Reset /*<< " sequence: " << j.sequence << " → " << i->at(j.name)*/ << std::endl;
			}
		}
	}
}

void trim_framework(AntibodySequence const &A, AntibodyFramework &heavy_fr, AntibodyFramework &light_fr)
{
	if( heavy_fr.fr1_end - heavy_fr.fr1_begin > 25 ) heavy_fr.fr1_begin = heavy_fr.fr1_end - heavy_fr.fr1_begin - 25;
	heavy_fr.fr4_begin = A.heavy.cdr3.end;  heavy_fr.fr4_end = A.heavy.cdr3.end + 12; // FR_H4 = heavy_chain[H3_end + 1: H3_end + 1 + 12]

	if( heavy_fr.fr4_end > A.heavy.sequence.size() ) heavy_fr.fr4_end = A.heavy.sequence.size();

	heavy_fr.update_sequences(A.heavy.sequence);

	//TR << "heavy_fr heavy_fr.fr4_begin:" << heavy_fr.fr4_begin << " heavy_fr.fr4_end:" << heavy_fr.fr4_end << " frh4:" << heavy_fr.fr4 << ":" << heavy_fr.fr4.size() << std::endl;

	if( light_fr.fr1_end - light_fr.fr1_begin > 23 ) light_fr.fr1_begin = light_fr.fr1_end - light_fr.fr1_begin - 23;
	light_fr.fr2_begin = A.light.cdr1.end;  light_fr.fr2_end = A.light.cdr1.end + 15;  // FR_L2 = light_chain[ L1_end + 1  :  L1_end + 1 + 15]
	light_fr.fr4_begin = A.light.cdr3.end;  light_fr.fr4_end = A.light.cdr3.end + 12;  // FR_L4 = light_chain[ L3_end + 1  :  L3_end + 1 + 12 ]

	if( light_fr.fr4_end > A.light.sequence.size() ) light_fr.fr4_end = A.light.sequence.size();

	light_fr.update_sequences(A.light.sequence);  //TR << "Light framework: " << light_fr << std::endl;
}



FRH_FRL calculate_frh_frl(AntibodySequence const &A)
{
	AntibodyFramework heavy_fr = A.heavy_framework();
	AntibodyFramework light_fr = A.light_framework();

	trim_framework(A, heavy_fr, light_fr);

	AntibodyNumbering an( Chothia_Numberer().number(A, heavy_fr, light_fr) );

	FRH_FRL r;

	struct {
		string sequence;
		AntibodyChainNumbering::NumberingVector numbering;
		uint begin, end;
		string &dest;
	} tasks[] {
		{heavy_fr.fr1, an.heavy.fr1,  10,  26, r.frh1},
		{heavy_fr.fr2, an.heavy.fr2,  36,  40, r.frh2},  {heavy_fr.fr2, an.heavy.fr2, 46, 50, r.frh2},
		{heavy_fr.fr3, an.heavy.fr3,  66,  95, r.frh3},
		{heavy_fr.fr4, an.heavy.fr4, 103, 110, r.frh4},

		{light_fr.fr1, an.light.fr1,  10,  24, r.frl1},
		{light_fr.fr2, an.light.fr2,  35,  39, r.frl2},  {light_fr.fr2, an.light.fr2, 45, 50, r.frl2},
		{light_fr.fr3, an.light.fr3,  57,  67, r.frl3},  {light_fr.fr3, an.light.fr3, 71, 89, r.frl3},
		{light_fr.fr4, an.light.fr4,  98, 105, r.frl4},
	};

	for(auto &t : tasks ) {
		auto s( t.sequence.begin() );  auto n( t.numbering.begin() );

		for(; s!=t.sequence.end(); ++s, ++n) {
			uint i = utility::string2Size( *n );
			if( i >= t.begin and i < t.end ) t.dest.push_back( *s );
		}
	}

	//r.frh = r.frh1 + r.frh2 + r.frh3 + r.frh4;
	//r.frl = r.frl1 + r.frl2 + r.frl3 + r.frl4;

	return r;
}


} // namespace grafting
} // namespace antibody
} // namespace protocols


#endif // __ANTIBODY_GRAFTING__
