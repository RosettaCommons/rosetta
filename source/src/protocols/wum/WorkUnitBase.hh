// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/wum/WorkUnitBase.hh
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_protocols_wum_WorkUnitBase_hh
#define INCLUDED_protocols_wum_WorkUnitBase_hh

#include <protocols/wum/WorkUnitBase.fwd.hh>

//#include <protocols/wum/WorkUnitList.hh>
// AUTO-REMOVED #include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
//#include <protocols/loops/BackboneDB.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <string>
#include <vector>
// AUTO-REMOVED #include <iostream>

#include <utility/vector1.hh>


namespace protocols {
namespace wum {

/// @brief  The base class for all work units

class WorkUnitBase  : public utility::pointer::ReferenceCount   {
  public:
    friend class MPI_WorkUnitManager;
    friend class MPI_WorkUnitManager_Slave;
    friend class WorkUnitManager;
		friend class WorkUnitQueue;

		WorkUnitBase();

    virtual ~WorkUnitBase (){}

		virtual protocols::wum::WorkUnitBaseOP clone() const {
			return protocols::wum::WorkUnitBaseOP( new WorkUnitBase( *this ) );
		}

    /// @brief Remove all data, make sure that the memory is also cleared, hence the cals to reserve
    virtual void clear_serial_data(){
      serial_data_.clear();
      serial_data_.reserve(2);
    	runtime_assert( serial_data_.capacity() == 2 );
		};

    /// @brief Run the workunit - overloaded by children of this class
    virtual void run();

    /// @brief Print header information to the stream, single line by default or verbose if verbose is set to true
		void print( std::ostream & out, bool verbose = false ) const ;

    /// @brief Accessor to the ID of the WorkUnit
		core::Size id(){ return header.id_; }

    /// @brief Accessor to the extra_data_1 and 3 field of the header
		//  Needed for blacklisting based on start_ir and ssid
		core::Size extra_data_1(){ return header.extra_data_1_; }
		core::Size extra_data_2(){ return header.extra_data_2_; }
		core::Size extra_data_3(){ return header.extra_data_3_; }

  	void set_extra_data_1( core::Size const value ){ header.extra_data_1_ = value; }
  	void set_extra_data_2( core::Size const value ){ header.extra_data_2_ = value; }
  	void set_extra_data_3( core::Size const value ){ header.extra_data_3_ = value; }

		/// @brief Adds to the blacklist
		void add_blacklist( int mpi_rank );

		/// @brief Erases the blacklist
		void clear_blacklist();
 
		/// @brief Finds in blacklist, true if is, false if it isn't
		bool in_blacklist( int mpi_rank );

    /// @brief Accesor to the "options" text field
		void set_options( const std::string &text );

  protected:
    /// @brief Make ready for sending
    virtual void serialize() {};

    /// @brief Make ready for working - i.o.w. take information in serial_data_ and turn it into whatever real data the derivative class has.
    virtual void deserialize() {} ;

    /// @brief Make a unique number out of Processor Number and unix timestamp ?
    virtual void create_unique_id(){};

protected:

    /// @brief Accessor to the serial data field
    std::string &serial_data() { return serial_data_; }

    /// @brief Accessor to the serial data field
    const std::string &serial_data() const { return serial_data_; }

    /// @brief Set the unixtime of the start of the execution of this WorkUnit
		void set_run_start();

    /// @brief Set the unixtime of the stop of the execution of this WorkUnit
		void set_run_stop();


public:

    /// @brief Returns the differrence between unix start and stop times
		core::Size get_run_time();

    /// @brief Accessor to header structure, return the WorkUnit Type
		std::string get_wu_type() const;

    /// @brief Accessor to header structure, sets the WorkUnit Type
		void set_wu_type( const std::string &text );

    /// @brief Optain the options string from the header
		std::string get_options() const;


private: // functions for serialization and sending of data
		unsigned int raw_data_size() const;

		/// @brief  This allocates and returns a pointer to a block of memory with the
    ///         workunit data totally serialized. Note that it is that *callers*
		///         responsibility to call delete []  on this memory when they're done with it.
		unsigned int raw_data_dump( unsigned char ** raw_data_ptr ) const;

    /// Read in Header and serial data from a raw block of memeory, like one you'd get from a lowlevel network communication
		void raw_data_load( const unsigned char * raw_data_ptr, unsigned int size );



public:

    /// @brief Return the memory usage of this WorkUnit
		virtual core::Size mem_footprint() const {
			return sizeof( WorkUnitBase::WU_Header ) + serial_data().capacity();
		}

		/// @brief this structure can contain any non-dynamicly allocated data.
		///  Any simple data types can be used here, ints, real, floats, char, etc..
		struct WU_Header{
			char  wu_type_[128];
      // Some unique id
      core::Size  id_;

      /// Unixtime of when this WU was created
      core::Size  unixtime_creation_;

      /// Unixtime of when this workunit began execution
      core::Size  unixtime_start_;

      /// Unixtime of when this workunit finished execution
      core::Size  unixtime_stop_;
      core::Size  extra_data_1_;
      core::Size  extra_data_2_;
      core::Size  extra_data_3_;
      core::Size  extra_data_4_;

      /// protocols can put arbitrary small header data here
      char  options_[128];
    };

		core::Size last_received_from(){ return last_received_from_; }
protected:
    /// @brief The header data
		WU_Header header;

    /// @brief Contains the serial number of whatever Rank/Node this WU was last receeived from
    core::Size last_received_from_;

		/// @brief Contains blacklist of nodes.  This data is NOT sent, and is only used on the sending side to determine where not to send a workunit.
		std::vector< int > blacklist_;

private: // data

    /// @brief Contains more data, such as decoys, silent structucts, .. whatever
    ///        really. Must be serialized, i.e. in plain text form
    ///        read for transmission.
    std::string serial_data_;
};




class WorkUnit_Wait: public WorkUnitBase {
  public:
     WorkUnit_Wait( long wait_time = 2  ):
      WorkUnitBase ()
     {
			header.extra_data_1_ = wait_time;
     }

		virtual ~WorkUnit_Wait(){}

		virtual protocols::wum::WorkUnitBaseOP clone() const {
			return protocols::wum::WorkUnitBaseOP( new WorkUnit_Wait( *this ) );
		}

    // @brief Run the workunit - overloaded by children of this class
    virtual void run();
	private:
};





/// @brief This WorkUnit type has structures in it. Most Workunits should derive from this one rather
/// THe the Base class.
class WorkUnit_SilentStructStore : public WorkUnitBase {
  public:
     WorkUnit_SilentStructStore():
      WorkUnitBase ()
     {
     }

		virtual ~WorkUnit_SilentStructStore(){}

		virtual protocols::wum::WorkUnitBaseOP clone() const {
			return protocols::wum::WorkUnitBaseOP( new WorkUnit_SilentStructStore( *this ) );
		}

    /// @brief This Work unit doesnt do *anything* - its just keeps the structures
    virtual void run(){};

    /// @brief write decoys into serial data store overwritinge whatever was there before. It basically syncs the silent struct store with the derial data
    virtual void serialize();

    /// @brief Make ready for working
    virtual void deserialize();

    /// @brief Accessor for decoy store
    const protocols::wum::SilentStructStore& decoys() const { return decoys_; }

    /// @brief Accessor for decoy store
    protocols::wum::SilentStructStore& decoys(){ return decoys_; }

  private:
    /// @brief This WorkUnit type has structures in it.
    protocols::wum::SilentStructStore decoys_;

};





/// @brief This WorkUnit type can encapsulate any MoverOP. When registering this WOrkunit
/// provide it with a MoverOP and then, when executed on the slaves, this workunit will run the mover
/// On every single input structure and return the results.
class WorkUnit_MoverWrapper : public WorkUnit_SilentStructStore {
  public:
    WorkUnit_MoverWrapper( protocols::moves::MoverOP the_mover );

		virtual ~WorkUnit_MoverWrapper(){}

		virtual protocols::wum::WorkUnitBaseOP clone() const {
			return protocols::wum::WorkUnitBaseOP( new WorkUnit_MoverWrapper( *this ) );
		}

    // @brief Run the workunit - overloaded by children of this class
    virtual void run();

  protected:

		void set_defaults();
	private:
		protocols::moves::MoverOP the_mover_;
};







}
}

#endif

