from fulmar.fulmar_constants import FULMAR_VERSION_STR
##############################################################################


class FulmarError(Exception):
    """To raise exceptions related to FULMAR
    """
    pass


class FulmarWarning(Warning):
    """Class from warning to be displayed as
    "FulmarWarning"
    """
    pass


def warning_on_one_line(message, category, filename, lineno,
                        file=None, line=None):
    """Function to display warnings on one line, as to avoid displaying
    'warnings.warn('warning message')' under the 'warning message'"""
    return ' %s:%s: %s: %s' % (filename, lineno, category.__name__, message)


def print_version():
    """Prints the version of FULMAR used."""
    print(FULMAR_VERSION_STR)
    return
