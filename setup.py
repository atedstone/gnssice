from glob import glob
from os import path, listdir
from typing import Optional

from setuptools import setup

FULLVERSION = "0.2.0"
VERSION = FULLVERSION
    
write_version = True


def write_version_py(filename: Optional[str] = None) -> None:
    cnt = """\
version = '%s'
short_version = '%s'
"""
    if filename is None:
        filename = path.join(path.dirname(__file__), "pygps", "version.py")

    a = open(filename, "w")
    try:
        a.write(cnt % (FULLVERSION, VERSION))
    finally:
        a.close()


if write_version:
    write_version_py()

with open("README.md", encoding="utf-8") as fh:
    long_description = fh.read()

scripts = [os.path.join('./bin', f) for f in listdir('./bin') if path.isfile(os.path.join('./bin', f))]

setup(
    name="PyGPS",
    version=FULLVERSION,
    description="Kinematic GNSS processing for ice flow",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    author="Andrew Tedstone",
    license="BSD-3",
    packages=["pygps"],
    install_requires=["pandas", "matplotlib", "scipy", "statsmodels"],
    scripts=scripts,
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: BSD License",
    ],
)
