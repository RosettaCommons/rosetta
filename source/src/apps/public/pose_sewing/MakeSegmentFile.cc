// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/pose_sewing/MakeSegmentFile.cc
/// @brief an app to make segment files
/// @author frankdt (frankdt@email.unc.edu)
/// @author Jared Adolf--Bryfogle (jadolfbr@gmail.com)


// MPI headers
#ifdef USEMPI
#include <mpi.h>
#endif

#include <devel/init.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/simple_metrics/SimpleMetric.hh>

#include <algorithm>
#include <regex>

#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/pointer/memory.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <core/conformation/Atom.hh>
#include <numeric/random/random.hh>
#include <numeric/HomogeneousTransform.hh>


#include <protocols/pose_sewing/data_storage/TerminalDSSPSortedPoseVector.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/rosetta_scripts/XmlObjects.hh>
#include <protocols/filters/Filter.hh>

#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

#include <apps/public/pose_sewing/MakeSegmentFile.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "apps.public.pose_sewing.MakeSegmentFile" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

FileOptionKey const motifs( "motif_file" );
FileOptionKey const pdbs( "pdb_file" );
FileOptionKey const segment_path( "segment_path" );
FileOptionKey const seg_script_opt( "seg_xml" );
StringOptionKey const seg_filters("seg_filters");
StringOptionKey const seg_metrics("seg_metrics");
BooleanOptionKey const output_elements_opt("output_elements");

// BooleanOptionKey const store_pdbs( "store_pdbs" );
// BooleanOptionKey const use_pdbs_only( "use_pdbs_only");

class MakeSegmentFileMPI : public protocols::moves::Mover {
public:
	MakeSegmentFileMPI();

	MakeSegmentFileMPI(
		std::string & motifs_list,
		std::string & segment_file_path,
		utility::vector1< protocols::filters::FilterOP > filters,
		utility::vector1< core::simple_metrics::SimpleMetricOP > metrics,
		bool store_pdbs = false ,
		bool use_pdbs_only = false,
		bool output_elements = false);

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override { return "MakeSegmentFileMPI"; }

	protocols::moves::MoverOP clone() const override {
		return protocols::moves::MoverOP( new MakeSegmentFileMPI( *this ) );
	}

	protocols::moves::MoverOP fresh_instance() const override {
		return protocols::moves::MoverOP( new MakeSegmentFileMPI );
	}

private:
	bool use_pdbs_only_;
	bool store_pdbs_;
	std::string segment_file_path_;
	std::string motifs_list_;
	utility::vector1< protocols::filters::FilterOP > filters_;
	utility::vector1< core::simple_metrics::SimpleMetricOP > metrics_;
	bool output_elements_;
};

MakeSegmentFileMPI::MakeSegmentFileMPI( std::string & motifs_list, std::string & segment_file_path, utility::vector1< protocols::filters::FilterOP > filters, utility::vector1< core::simple_metrics::SimpleMetricOP > metrics, bool store_pdbs, bool use_pdbs_only,  bool output_elements):
	protocols::moves::Mover(),
	use_pdbs_only_(use_pdbs_only),
	store_pdbs_(store_pdbs),
	segment_file_path_(segment_file_path),
	motifs_list_(motifs_list),
	filters_(filters),
	metrics_(metrics),
	output_elements_(output_elements)
{
	store_pdbs_ = store_pdbs;
	use_pdbs_only_ = use_pdbs_only;
	output_elements_ = output_elements;
}

MakeSegmentFileMPI::MakeSegmentFileMPI():
	protocols::moves::Mover()
{
	store_pdbs_ = false;
	use_pdbs_only_ = false;
	output_elements_ = false;
}


void
MakeSegmentFileMPI::apply( core::pose::Pose & pose ) {
	protocols::pose_sewing::data_storage::TerminalDSSPSortedPoseVectorOP test_vector = protocols::pose_sewing::data_storage::TerminalDSSPSortedPoseVectorOP(new protocols::pose_sewing::data_storage::TerminalDSSPSortedPoseVector);


	test_vector->simultaneously_populate_and_write_segment_file_pdb(pose, motifs_list_, segment_file_path_, store_pdbs_, use_pdbs_only_, filters_,metrics_,  output_elements_);

}



int
main( int argc, char * argv [] )
{
	try{
		option.add( motifs, "motif file name");
		option.add( segment_path, "Path where the segment file will go.");
		option.add( seg_script_opt, "Path to XML for optional filtering and metrics.");
		option.add(seg_filters,"Comma-separated list of filters to run from a provided parser:protocol.  These will be run on each potential segment matching motif.");
		option.add(seg_metrics,"Comma-separated list of metrics to run from a provided parser:protocol.  These will be added to the output pdb and in the future used for on-the-fly filtering.");
		option.add(output_elements_opt, "Output elements into pdb labels?  Useful for sheets.  Default false.");

		// option.add( store_pdbs, "Store PDBs as a reference");
		// option.add( use_pdbs_only, "Only use PDBs instead of mmtf with optional reference");
		devel::init(argc,argv);
		std::string motifs_list = option[motifs].value();

		if ( ! option[in::file::l].user() ) {
			utility_exit_with_message("PDBLIST must be set with the in:file:l option");
		}
		if ( option[in::file::s].user() ) {
			utility_exit_with_message("option in:file:s not supported!  Please use -l");
		}

		std::string segment_file_path = option[segment_path].value();

		//Make directories if they don't exist. I hate making directories ahead of time.

		if ( ! utility::file::is_directory(segment_file_path) ) {
			utility::file::create_directory_recursive(segment_file_path);
		}

		// bool store_pdbs_ = option[store_pdbs]();
		// bool use_pdbs_only_ = option[use_pdbs_only]();

		//We use our own output.  We don't want score files or echo PDBs.
		option[jd2::no_output].value(true);

		//Cannot declare make_shared in an MPI environment, so leave it to MPI JD to figure it all out.


		//Default an empty filter and metric list.
		utility::vector1< protocols::filters::FilterOP > filters;
		utility::vector1< core::simple_metrics::SimpleMetricOP > metrics;

		bool output_elements = false;
		if ( option[output_elements_opt].user() ) {
			output_elements = option[output_elements_opt].value();
		}
		if ( option[seg_script_opt].user() ) {
			std::string xml_protocol = option[seg_script_opt].value();
			TR << "Parsing Protocol for filters and metrics: " << xml_protocol << std::endl;
			protocols::rosetta_scripts::XmlObjectsCOP xml_ob = protocols::rosetta_scripts::XmlObjects::create_from_file(xml_protocol);
			if ( option[seg_filters].user() ) {
				std::string seg_filter_string = option[seg_filters].value();
				TR << "Filtering segments by: " << seg_filter_string << std::endl;
				utility::vector1< std::string > seg_filtersSP = utility::string_split(seg_filter_string, ',');
				for ( auto const & filter_string : seg_filtersSP ) {
					protocols::filters::FilterOP current_filter = xml_ob->get_filter(filter_string);
					filters.push_back(current_filter);
				}
			}
			if ( option[seg_metrics].user() ) {
				std::string seg_metrics_string = option[seg_metrics].value();
				TR << "Outputting metrics: " << seg_metrics_string << std::endl;
				utility::vector1< std::string > seg_metricsSP = utility::string_split(seg_metrics_string, ',');
				for ( auto const & metric_string : seg_metricsSP ) {
					core::simple_metrics::SimpleMetricOP current_metric = xml_ob->get_simple_metric(metric_string);
					metrics.push_back(current_metric);
				}
			}
		}

#ifdef USEMPI
			//JD2 machinery.
				// protocols::moves::MoverOP seg_wrapper( new MakeSegmentFileMPI (
				// motifs_list,
				// segment_file_path,
				// store_pdbs_,
				// use_pdbs_only_));

				protocols::moves::MoverOP seg_wrapper( new MakeSegmentFileMPI (
				motifs_list,
				segment_file_path,
				filters,
				metrics,
				true,
				true,
				output_elements));

				protocols::jd2::JobDistributor::get_instance()->go(seg_wrapper);



			//Start of Manual implementation
			//MPI_Comm_rank( MPI_COMM_WORLD, &MPI_rank );
			//MPI_Comm_size( MPI_COMM_WORLD, &MPI_n_procs );

			//if (TR.Debug.visible()) {
			//	TR.Debug << "MPI rank " << MPI_rank << " of " << MPI_n_procs << " reporting for duty." << std::endl;
			//}
			//wait_for_proc_zero();


#endif


#ifndef USEMPI
		//std::string pdb_list = option[in::file::l].value().vector()[1]; //Vector1 of FileNames.  If not MPI, we use the first one.
		std::string pdb_list = option[in::file::l].value().vector().at(0).name(); //Vector1 of FileNames.  If not MPI, we use the first one.

		protocols::pose_sewing::data_storage::TerminalDSSPSortedPoseVectorOP test_vector = protocols::pose_sewing::data_storage::TerminalDSSPSortedPoseVectorOP(new protocols::pose_sewing::data_storage::TerminalDSSPSortedPoseVector);


		// test_vector->simultaneously_populate_and_write_segment_file(motifs_list, pdb_list, segment_file_path, store_pdbs_, use_pdbs_only_);
		test_vector->simultaneously_populate_and_write_segment_file(motifs_list, pdb_list, segment_file_path, true, true, filters, metrics, output_elements);
#endif


		return 0;

	} catch (utility::excn::Exception& excn ) {
		excn.display();
		std::exit( 1 );
	}
}


//int
//main( int argc, char * argv [] )
//{
// try{
//  main( argc, argv );
//  return 0;
// } catch (utility::excn::Exception& excn ) {
//  excn.display();
//  std::exit( 1 );
// }
//}




