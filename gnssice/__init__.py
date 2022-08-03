from gnssice import gnss, pp

try:
    from gnssice.version import version as __version__  # noqa
except ImportError:  # pragma: no cover
    raise ImportError(
        "package is not properly installed. If you are "
        "running from the source directory, please instead "
        "create a new virtual environment (using conda or "
        "virtualenv) and then install it in-place by running: "
        "pip install -e ."
    )
