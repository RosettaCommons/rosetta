#!/usr/bin/env python

from __future__ import division

from sqlalchemy import *
from sqlalchemy.orm import *
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class Job (Base):
    __tablename__ = 'jobs'

    id = Column(BigInteger, primary_key=True)
    command = Column(Text)
    revision = Column(Text)
    algorithm = Column(Text)
    start_time = Column(Text)
    stop_time = Column(Text)
    input_file = Column(Text)
    loop_begin = Column(Integer)
    loop_end = Column(Integer)
    iterations = Column(Integer)
    frames = Column(Integer)
    description = Column(Text)

    moves = relationship('Move', backref='job')
    temperatures = relationship('Temperature', backref='job')
    trajectories = relationship('Trajectory', backref='job')
    numpy_caches = relationship('NumpyCache', backref='job')

    def __repr__(self):
        return '<job id=%d>' % self.id

    @property
    def age(self):
        from helpers import plural
        from datetime import datetime

        time_string = self.stop_time.strip()
        job_time = datetime.strptime(time_string, '%a %b %d %H:%M:%S %Y')
        job_age = datetime.now() - job_time

        seconds = job_age.seconds
        minutes = seconds // 60
        hours = minutes // 60
        days = job_age.days
        weeks = days // 7
        months = days // 30
        years = days // 365

        if years:       return "%d year%s" % plural(years)
        elif months:    return "%d month%s" % plural(months)
        elif weeks:     return "%d week%s" % plural(weeks)
        elif days:      return "%d day%s" % plural(days)
        elif hours:     return "%d hour%s" % plural(hours)
        elif minutes:   return "%d min%s" % plural(minutes)
        else:           return "%d sec%s" % plural(seconds)

    @property
    def duration(self):
        from datetime import datetime
        
        start_time = datetime.strptime(
                self.start_time.strip(), '%a %b %d %H:%M:%S %Y')
        stop_time = datetime.strptime(
                self.stop_time.strip(), '%a %b %d %H:%M:%S %Y')

        duration = stop_time - start_time
        seconds = duration.seconds % 60
        minutes = (duration.seconds // 60) % 60
        hours = duration.seconds // 3600

        if hours:   return '%d:%02d:%02d' % (hours, minutes, seconds)
        else:       return '%d:%02d' % (minutes, seconds)

    @property
    def protein(self):
        from os.path import basename
        return basename(self.input_file)

    @property
    def name(self):
        from os.path import splitext
        return splitext(self.protein)[0]

    @property
    def loop(self):
        return '%d-%d' % (self.loop_begin, self.loop_end)


class Move (Base):
    __tablename__ = 'moves'
    __table_args__ = (
            ForeignKeyConstraint(
                ['job_id', 'temp_level'],
                ['temperatures.job_id', 'temperatures.level']),
    )

    job_id = Column(BigInteger, ForeignKey('jobs.id'), primary_key=True)
    type = Column(String, primary_key=True)
    temp_level = Column(Integer, primary_key=True)
    num_trials = Column(Integer)
    num_accepted = Column(Integer)

    temperature = relationship("Temperature", uselist=False)

    def __repr__(self):
        return '<move job=%d type=%s temp=%d rate=%f>' % (
                self.job_id, self.type, self.temp_level, self.efficiency)

    @property
    def efficiency(self):
        return self.num_accepted / self.num_trials


class Temperature (Base):
    __tablename__ = 'temperatures'

    job_id = Column(BigInteger, ForeignKey('jobs.id'), primary_key=True)
    level = Column(Integer, primary_key=True)
    temperature = Column(Float, nullable=False)

    def __repr__(self):
        return '<temperature job=%d level=%d temp=%f>' % (
                self.job_id, self.level, self.temperature)

    @property
    def formatted(self):
        return '%.2f' % self.temperature


class Trajectory (Base):
    __tablename__ = 'trajectories'

    job_id = Column(BigInteger, ForeignKey('jobs.id'), primary_key=True)
    iteration = Column(Integer, primary_key=True)
    score = Column(Float)
    silent_pose = Column(LargeBinary)

    def __repr__(self):
        return '<trajectory job=%d iter=%d>' % (self.job_id, self.iteration)


class NumpyCache (Base):
    __tablename__ = 'numpy_cache'
    __table_args__ = {'mysql_engine': 'InnoDB'}

    job_id = Column(BigInteger, ForeignKey('jobs.id'), primary_key=True)
    type = Column(String(50), primary_key=True)
    data = Column(LargeBinary(2**24))   # 16MB

    def __repr__(self):
        return '<numpy_cache job=%d type=%s>' % (self.job_id, self.type)


