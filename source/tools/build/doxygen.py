# -*- mode:python;indent-tabs-mode:nil;show-trailing-whitespace:t; -*-
#
# Doxygen SCons tool
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import os


def generate(env):
   """
   Add builders and construction variables for the Doxygen tool.
   """

   # Create a builder for
   doxygen_builder = env.Builder(action = env.Action("$DOXYGEN_COM"))

   # Add the Doxygen() method to the environment
   env.Append(BUILDERS = { "Doxygen" : doxygen_builder })

   # Set environment variables appropriately

   # The Doxygen command
   env["DOXYGEN"] = "doxygen"
   # Doxygen has no command line input: all configuration is done via a file.
   env["DOXYGEN_FLAGS"] = ""
   env["DOXYGEN_COM"] = "$DOXYGEN $DOXYGEN_FLAGS $SOURCE" # $DOXYFILE_FILE


def exists(env):
   """
   Make sure doxygen exists.
   """
   return env.Detect("doxygen")
