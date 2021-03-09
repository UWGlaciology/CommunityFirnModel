import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '1.0.8'
PACKAGE_NAME = 'CommunityFirnModel'
AUTHOR = 'Max Stevens'
AUTHOR_EMAIL = 'maxstev@umd.edu'
URL = 'https://github.com/UWGlaciology/CommunityFirnModel'

LICENSE = 'MIT'
DESCRIPTION = 'An open-source model to simulate firn.'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "markdown"

INSTALL_REQUIRES = [
      'numpy>1.12.0',
      'pandas',
      'scipy>1.0.0',
      'h5py'
]

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"],
      python_requires='>=3.6',
      packages=find_packages()
      )



