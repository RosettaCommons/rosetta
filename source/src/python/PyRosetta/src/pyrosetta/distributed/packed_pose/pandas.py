import functools
import pandas
import numpy


def _is_dataframe_index_boring(index):
    if not isinstance(index, pandas.Int64Index):
        return False

    if pandas.Index(numpy.arange(0, len(index))).equals(index):
        return True
    elif pandas.Index(numpy.repeat(0, len(index))).equals(index):
        return True
    else:
        return False


def register_pandas_container_traversal(generic_func, dict_func):
    @generic_func.register(pandas.DataFrame)
    def dataframe_traveral(dataframe):
        if _is_dataframe_index_boring(dataframe.index):
            return generic_func(dataframe.to_dict("records"))
        else:
            if dataframe.index.has_duplicates:
                raise ValueError(
                    "Unable to coerce duplicate-indexed dataframe for traversal."
                    " Consider '.reset_index'."
                )
            return {i: generic_func(v) for i, v in dataframe.to_dict("index").items()}

    @generic_func.register(pandas.Series)
    def series_traversal(series):
        return generic_func(dict(series))

    return generic_func
