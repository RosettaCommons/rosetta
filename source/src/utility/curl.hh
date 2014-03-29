// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/curl.hh
/// @brief  Implements a OO interface to CURL wrapping and hiding the ugly C callbacks.
/// @author Mike Tyka (mike.tyka@gmail.com)

#ifndef INCLUDED_utility_curl_hh
#define INCLUDED_utility_curl_hh

#ifdef WITHCURL

#include <curl/curl.h>
#include <string>


namespace utility {

// must include -lcurl during linking if you're using -dwithcurl

class CurlGet {
 public:
  CurlGet();
  static int writer(char *data, std::size_t size, std::size_t nmemb, std::string *buffer);
 private:
  char *getErrorBuffer();
  std::string * getbuffer();
 public:
  std::string get( const std::string &url );
 private:
  // Write any errors in here
  char errorBuffer[CURL_ERROR_SIZE];
  // Write all expected data in here
  std::string buffer;
};


class CurlPost {
 public:
  CurlPost();
  // This is the writer call back function used by curl
  static int writer(char *data, std::size_t size, std::size_t nmemb, std::string *buffer);
 private:
  char *getErrorBuffer();
  std::string * getreadbuffer();
  std::string * getwritebuffer();
 public:
  std::string post( const std::string &url, const std::string &data, const std::string &fields );
 private:
  // Write any errors in here
  char errorBuffer[CURL_ERROR_SIZE];
  // Write all expected input and out data in here
  std::string readbuffer;
  std::string writebuffer;
};

}

#endif // WITHCURL

#endif // INCLUDED_utility_curl_hh
