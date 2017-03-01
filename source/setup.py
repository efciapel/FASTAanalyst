"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))


setup(
    name='FASTAanalyst',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='1.0.0',

    description='Biblioteka do analizy sekwencji DNA z plików FASTA',



    # Author details
    author='M. Kepska, E. Krol, A. Osina',
    author_email='anna.osina93@gmail.com',

    # Choose your license
    license='open-source',


    classifiers=[
        'Development Status :: 3 - Alpha',

        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.5'
    ],


    keywords='fasta analise',

    packages=find_packages(exclude=['__init__', 'fa']),

    py_modules=['fa'],

    install_requires=['fa'],

    package_data={
        'fastaanalyst': ['__init__', 'fa'],
    },

    entry_points={
        'console_scripts': [
            'sample=sample:main',
        ],
    },
)