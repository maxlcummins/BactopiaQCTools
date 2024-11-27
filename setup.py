# setup.py
from setuptools import setup, find_packages

setup(
    name='bactQC',
    version='0.0.2',
    packages=find_packages(),
    author='Max Cummins',
    author_email='max.l.cummins@gmail.com',
    description='A package for bacterial genome quality control of bactopia outputs',
    url='https://github.com/maxlcummins/bactopiaQCtools',
    install_requires=[
        'pandas',
        'requests',
        'click',
        'emoji'
    ],
    entry_points={
        'console_scripts': [
            'bactQC = bactQC.cli:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
    python_requires='>=3.6',
)
