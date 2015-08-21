// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/curl.cc
/// @brief  Implements a OO interface to CURL wrapping and hiding the ugly C callbacks.
/// @author Mike Tyka (mike.tyka@gmail.com)

#ifdef WITHCURL

#include <utility/curl.hh>
#include <utility/excn/Exceptions.hh>


namespace utility {

// must include -lcurl during linking if you're using -dwithcurl

  CurlGet::CurlGet() 
  {
  }

  // This is the writer call back function used by curl
  int 
  CurlGet::writer(char *data, size_t size, size_t nmemb, std::string *buffer)
  {
    // What we will return
    int result = 0;

    // Is there anything in the buffer?
    if (buffer != NULL)
    {
      // Append the data to the buffer
      buffer->append(data, size * nmemb);

      // How much did we write?
      result = size * nmemb;
    }

    return result;
  }

  char *
  CurlGet::getErrorBuffer() 
  { 
    return &errorBuffer[0]; 
  }

  std::string * 
  CurlGet::getbuffer() 
  { 
    return &buffer; 
  }

  std::string 
  CurlGet::get( const std::string &url )
  {
    // Our curl objects
    CURL *curl;
    CURLcode result;

    // Create our curl handle
    curl = curl_easy_init();

    if (curl)
    {
      // Now set up all of the curl options
      curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, getErrorBuffer() );
      curl_easy_setopt(curl, CURLOPT_URL, url.c_str() );
      curl_easy_setopt(curl, CURLOPT_HEADER, 0);
      curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1);
      curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, this->writer);
      curl_easy_setopt(curl, CURLOPT_WRITEDATA, getbuffer() );

      // Attempt to retrieve the remote page
      result = curl_easy_perform(curl);

      // Always cleanup
      curl_easy_cleanup(curl);

      // Did we succeed?
      if (result == CURLE_OK)
      {
        return buffer;
      }
      else
      {
        std::stringstream ss;
        ss << "Error: [" << result << "] - " << errorBuffer;
        throw std::string( ss.str() );
      }
    }
  }


  CurlPost::CurlPost() 
  { 
  }

  // This is the writer call back function used by curl
  int 
  CurlPost::writer(char *data, size_t size, size_t nmemb, std::string *buffer)
  {
    // What we will return
    int result = 0;
    //std::cout << "Buffer: " << *buffer << std::endl;
    // Is there anything in the buffer?
    if (buffer != NULL)
    {
      // Append the data to the buffer
      buffer->append(data, size * nmemb);

      // How much did we write?
      result = size * nmemb;
    }

    return result;
  }

  char *
  CurlPost::getErrorBuffer() 
  { 
    return &errorBuffer[0]; 
  }

  std::string * 
  CurlPost::getreadbuffer() 
  { 
    return &readbuffer; 
  }

  std::string * 
  CurlPost::getwritebuffer() 
  { 
    return &writebuffer; 
  }

  std::string 
  CurlPost::post( const std::string &url, const std::string &data, const std::string &fields )
  {
    // Our curl objects
    CURL *curl;
    CURLcode result;

    // Create our curl handle
    curl = curl_easy_init();

    if (curl)
    {
      // fill read buffer
      readbuffer = data;

      // Now set up all of the curl options
      curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, getErrorBuffer() );
      curl_easy_setopt(curl, CURLOPT_URL, url.c_str() );
      curl_easy_setopt(curl, CURLOPT_POSTFIELDS, fields.c_str() );
      curl_easy_setopt(curl, CURLOPT_HEADER, 0);
      curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1);
      curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, this->writer);
      curl_easy_setopt(curl, CURLOPT_WRITEDATA, getwritebuffer() );

      // Attempt to retrieve the remote page
      result = curl_easy_perform(curl);

      // Always cleanup
      curl_easy_cleanup(curl);

      // Did we succeed?
      if (result == CURLE_OK)
      {
        return writebuffer;
      }
      else
      {
        std::stringstream ss;
        ss << "Error: [" << result << "] - " << errorBuffer;
        throw std::string( ss.str() );
      }
    }
  }


}

#endif

