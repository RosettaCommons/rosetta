# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import base64
import collections
import pickle


class PoseScoreSerializerBase(object):
    """Base class for `PoseScoreSerializer` methods."""

    @staticmethod
    def to_pickle(value):
        try:
            return pickle.dumps(value)
        except (TypeError, OverflowError, MemoryError, pickle.PicklingError) as ex:
            raise TypeError(
                "Only pickle-serializable object types are allowed to be set "
                + "as score values. Received: %r. %s" % (type(value), ex)
            )

    @staticmethod
    def from_pickle(value):
        try:
            return pickle.loads(value)
        except (TypeError, OverflowError, MemoryError, EOFError, pickle.UnpicklingError) as ex:
            raise TypeError(
                "Could not deserialize score value of type %r. %s" % (type(value), ex)
            )

    @staticmethod
    def to_base64(value):
        return base64.b64encode(value).decode()

    @staticmethod
    def from_base64(value):
        return base64.b64decode(value)

    @staticmethod
    def to_base64_pickle(value):
        return PoseScoreSerializerBase.to_base64(PoseScoreSerializerBase.to_pickle(value))

    @staticmethod
    def from_base64_pickle(value):
        return PoseScoreSerializerBase.from_pickle(PoseScoreSerializerBase.from_base64(value))

    @staticmethod
    def bool_from_str(value):
        if value == "True":
            return True
        elif value == "False":
            return False
        else:
            raise NotImplementedError(value)


class PoseScoreSerializer(PoseScoreSerializerBase):
    """
    Serialize and deserialize score values for CustomStringValueMetric SimpleMetric.

    Examples:
        Automatically serialize an arbitrary score value:
            `pose.cache["foo"] = value`
        Automatically deserialize an arbitrary score value:
            `value = pose.cache["foo"]`
        Manually serialize an arbitrary score value:
            `value = PoseScoreSerializer.maybe_encode(value)`
        Manually deserialize an arbitrary score value:
            `value = PoseScoreSerializer.maybe_decode(value)`
    """
    # Define different data types with human-readable custom prefixes in case anyone
    # accesses the serialized score values outside the scope of the `PoseScoreSerializer`
    _CustomTypeMetric = collections.namedtuple(
        "CustomTypeMetric", ["type", "prefix", "encode_func", "decode_func"],
    )
    _reserved_types = (str, float)
    _custom_type_metrics = {
        "bool": _CustomTypeMetric(
            type=bool,
            prefix="[CustomBooleanValueMetric]",
            encode_func=str,
            decode_func=PoseScoreSerializerBase.bool_from_str,
        ),
        "int": _CustomTypeMetric(
            type=int,
            prefix="[CustomDiscreteValueMetric]",
            encode_func=str,
            decode_func=int,
        ),
        "bytes": _CustomTypeMetric(
            type=bytes,
            prefix="[CustomBinaryValueMetric]",
            encode_func=PoseScoreSerializerBase.to_base64,
            decode_func=PoseScoreSerializerBase.from_base64,
        ),
        "object": _CustomTypeMetric(
            type=object,
            prefix="[CustomArbitraryValueMetric]",
            encode_func=PoseScoreSerializerBase.to_base64_pickle,
            decode_func=PoseScoreSerializerBase.from_base64_pickle,
        ),
    }

    @staticmethod
    def maybe_encode(value):
        """Serialize the input value into a `str` object if it's not a `str` or `float` object."""
        if not isinstance(value, PoseScoreSerializer._reserved_types):
            for _custom_type_metric in PoseScoreSerializer._custom_type_metrics.values():
                if isinstance(value, _custom_type_metric.type):
                    value = "{0}{1}".format(
                        _custom_type_metric.prefix,
                        _custom_type_metric.encode_func(value)
                    )
                    break
            else:
                raise NotImplementedError("Unsupported object type: {0}".format(type(value)))

        return value

    @staticmethod
    def maybe_decode(value):
        """Deserialize the input value if it's serialized."""
        if isinstance(value, str):
            for _custom_type_metric in PoseScoreSerializer._custom_type_metrics.values():
                if value.startswith(_custom_type_metric.prefix):
                    value = _custom_type_metric.decode_func(
                        value[len(_custom_type_metric.prefix):]
                    )
                    break
            else:
                pass # Return without decoding because string doesn't start with a prefix

        return value
