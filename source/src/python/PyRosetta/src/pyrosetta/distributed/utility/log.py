from __future__ import division

import datetime
import logging
import collections
from decorator import decorator
import pprint


class classproperty(object):
    """Simple decorator for class-level dynamic properies."""

    def __init__(self, getter):
        self.getter = getter

    def __get__(self, instance, owner):
        return self.getter(owner)


class LoggerMixin(object):
    """Mixin class adding 'logger' property initialized logger '<module>.<name>'"""

    @classproperty
    def logger(cls):
        return cls._get_logger()

    @classmethod
    def _get_logger(cls):
        if "_logger" not in cls.__dict__:
            cls._logger = logging.getLogger(
                "%s.%s" % (cls.__module__, cls.__name__))
        return cls._logger


def log_method(level=logging.INFO, log_params=True, log_enter=True, log_exit=False, log_return=None):
    """Logging decorator which logs method calls with class 'logger' property."""

    if log_return is None:
        log_return = log_params

    def log_method(func, self, *args, **kwargs):
        if not self.logger.isEnabledFor(level):
            return func(self, *args, **kwargs)

        if log_enter and log_params:
            self.logger.log(level, "%s(args=%r, kwargs=%r)", func.__name__, args, kwargs)
        elif log_enter:
            self.logger.log(level, "%s(...)", func.__name__)

        f_result = func(self, *args, **kwargs)

        if log_exit and log_return:
            self.logger.log(level, "%s return %r", func.__name__, f_result)
        elif log_exit:
            self.logger.log(level, "%s return ...", func.__name__)

        return f_result

    return decorator(log_method)


class ProgressLogger(object):
    """Issues incremental log messages for a set of tasks."""

    def __init__(self, logger, prefix, level=logging.INFO, seconds=None, iterations=None, fraction=None, total=None):
        """Initialize progress logger.

        logger - target logger
        level - target log level
        prefix - string prefix in log message
        seconds - emit log message after given seconds seconds
        iterations - emit log message after given interations
        fraction - emit log message at given fraction of total run, valid only if total provided
        total - total number of tasks
        """
        self.logger = logger
        self.level = level
        self.prefix = prefix
        self.seconds = seconds
        self.iterations = iterations
        self.fraction = fraction
        self.total = total

        self.iteration = 0
        self.last_iteration = 0

        self.start_time = datetime.datetime.now()
        self.last_time = self.start_time

    def increment(self, by=1):
        """Indicate that 'by' tasks have been completed and log if needed."""

        self.iteration += by
        now = datetime.datetime.now()

        perform_logging = False
        if self.seconds and (now - self.last_time).total_seconds() >= self.seconds:
            perform_logging = True
        elif self.iterations and (self.iteration - self.last_iteration) >= self.iterations:
            perform_logging = True
        elif self.fraction and self.total and (self.iteration - self.last_iteration) / self.total >= self.fraction:
            perform_logging = True
        elif self.total and self.iteration >= self.total and self.last_iteration < self.total:
            perform_logging = True

        if perform_logging:
            self.last_time = now
            self.last_iteration = self.iteration

            self.logger.log(self.level, "%s - %i/%s", self.prefix, self.iteration, self.total if self.total else "?")

    def __iadd__(self, other):
        self.increment(other)
        return self

    def logging_iter(self, iterable):
        """Wrap iterator in logging increment."""
        for i in iterable:
            self += 1
            yield i

    def count(self, start=0, step=1):
        """Generator returning evenly spaced values starting with n, incrementing on iteration."""

        n = start
        while True:
            self += step

            yield n
            n += step

    def iterate(self, sequence):
        """Return a generator that wraps an iterable, performing increment calls on each iteration.

        If sequence is sized and self.total is None, self.total will be updated to length of sequence.
        """

        if self.total is None and isinstance(sequence, collections.Sized):
            self.total = len(sequence)

        for v in sequence:
            self.increment()
            yield v

    def enumerate(self, sequence):
        """Return a generator that enumerates an iterable, performing increment calls on each iteration.

        If sequence is sized and self.total is None, self.total will be updated to length of sequence.
        """

        if self.total is None and isinstance(sequence, collections.Sized):
            self.total = len(sequence)

        return zip(self.count(), sequence)


class Plog(object):
    """Lazy pprint formatter for logging."""

    def __init__(self, obj, *args, **kwargs):
        """Initialize Plog container for given object.

        Args:
            obj: target object
            *kwargs: pprint.pformat kwargs


        """
        self.obj = obj
        self.kwargs = kwargs

    def __repr__(self):
        """pformat obj."""
        return pprint.pformat(self.obj, **self.kwargs)


class LoggingContext(object):
    """Context handler to modify logging configuration with context.

    Example:
        with LoggingContext(rosetta.logger, level=logging.WARNING):
            #...very chatty functions...

        # Back to normal logging
        """
    def __init__(self, logger, level=None, handler=None, close=True):
        self.logger = logger
        self.level = level
        self.handler = handler
        self.close = close

    def __enter__(self):
        if self.level is not None:
            self.old_level = self.logger.level
            self.logger.setLevel(self.level)
        if self.handler:
            self.logger.addHandler(self.handler)

    def __exit__(self, et, ev, tb):
        if self.level is not None:
            self.logger.setLevel(self.old_level)
        if self.handler:
            self.logger.removeHandler(self.handler)
        if self.handler and self.close:
            self.handler.close()
        # implicit return of None => don't swallow exceptions
