import os
import sys
import json
import logging

import pyrosetta.rosetta as rosetta

_logger = logging.getLogger("rosetta")
_logging_tracer = None
_initialized_handler = False


def maybe_initialize_handler(mode):
    """Initialize stdout logging handler if not already initialized."""
    global _initialized_handler

    # Check if rosetta or root logger have had handlers initialized.
    # Preceeding calls to logging.basicConfig will initialize root handlers.

    # If not, add an output handler to logger 'rosetta' targeting stdout and
    # set log level to INFO.
    if _logger.handlers:
        # A rosetta logging handler is already initialized
        return

    if mode is "logging":
        # We're in "quiet mode", handle all logs via the standard logging config.
        return

    if logging.root.handlers:
        if mode is not "interactive":
            # A root logger is already configured and we're in non-interactive mode.
            return

        if mode is "interactive" and logging.root.isEnabledFor(logging.INFO):
            # Default ipython notebook configuration registers a root handler,
            # but does not enable it at the "INFO" level.

            # If info-level handler is enabled we can assume logging is
            # already configured is some way.
            return

    if _initialized_handler:
        # We've already initialized a handler.
        return

    # There no logger configured that's currently displaying log messages
    # Enable a default handler that pushes to stdout for the default case.
    _logger.addHandler(logging.StreamHandler(sys.stdout))
    _logger.setLevel(logging.INFO)
    _logger.propagate = False
    _initialized_handler = True


class PythonLoggingSink(rosetta.basic.PyTracer):
    """Logging sink to dispatch all log messages through the logging module."""

    def __init__(self):
        super(PythonLoggingSink, self).__init__()
        self.set_ios_hook = rosetta.basic.Tracer.set_ios_hook

        # Set the tracer hook, do not enable 'raw' in order to
        # allow rosetta's tracing settings to have priority.
        self.set_ios_hook(
            self, rosetta.basic.Tracer.get_all_channels_string(), False)

    def __del__(self):
        self.set_ios_hook(None, '')

    def output_callback(self, logstring):
        _logger.info(logstring.strip())


def set_logging_sink():
    """Redirect all rosetta trace output through the logging.Logger 'rosetta'."""
    global _logging_tracer

    if _logging_tracer:
        # Tracer hook is already established
        return

    _logging_tracer = PythonLoggingSink()

    # Mute rosetta's standard logging output
    rosetta.basic.Tracer.super_mute(True)
