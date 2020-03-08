class RedshiftWarning(UserWarning):
    # Warn user that redshift is not being used.
    pass

class NotADictError(Exception):
    # Must use a dict.
    pass