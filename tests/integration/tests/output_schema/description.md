# output\_schema integration test

## Author

Vikram K. Mulligan, Centre for Computational Biology, Flatiron Institute (vmulligan@flatironinstitute.org)

## Description

This tests the ability of the rosetta\_scripts application to generate and write out the XML schema definition (XSD) for the RosettaScripts scripting language.  Changes to this test are expected when there are changes to the description or interface of existing scriptable modules (movers, filters, task operations, etc.) or when there are new scriptable modules added.  The test also helps to ensure that all movers, filters, task operations, etc. that are registered with their respective factories can be instantiated without specialized configuration.  Without this, these modules are not usable in RosettaScripts in any case (and it wouldn't make sense to have them registered with the factories).
