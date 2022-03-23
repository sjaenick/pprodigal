
from os import path
from setuptools import setup
import pprodigal


setup(
    name='pprodigal',
    version=pprodigal.__version__,
    description='PProdigal: Parallelized gene prediction based on Prodigal.',
    keywords=['bioinformatics', 'gene prediction'],
    long_description=open('README.txt').read(),
    license='GPLv3',
    author='Sebastian Jaenicke',
    author_email='sebastian.jaenicke@computational.bio.uni-giessen.de',
    url='https://github.com/sjaenick/pprodigal',
    packages=['pprodigal'],
    python_requires='>=3.0',
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'pprodigal=pprodigal.pprodigal:main'
        ]
    },
    classifiers=[
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3 :: Only',
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English'
    ],
    project_urls={
        'Bug Reports': 'https://github.com/sjaenick/pprodigal/issues',
        'Source': 'https://github.com/sjaenick/pprodigal'
    },
)
