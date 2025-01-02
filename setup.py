# setup.py
import os
import re
from setuptools import setup, find_packages

# Path to the version.py file
version_file = os.path.join(os.path.dirname(__file__), 'bactQC', 'version.py')

# Read the version string from version.py without importing the package
with open(version_file, 'r') as vf:
    version_content = vf.read()

# Use regex to extract the __version__ variable
version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_content, re.M)
if version_match:
    version = version_match.group(1)
else:
    raise RuntimeError("Unable to find __version__ string in version.py.")


setup(
    name='bactQC',
    version=version,  # Dynamically set the version
    packages=find_packages(),
    author='Max Cummins',
    author_email='max.l.cummins@gmail.com',
    description='A package for bacterial genome quality control of bactopia outputs',
    url='https://github.com/maxlcummins/bactQC',
    install_requires=[
        'pandas',
        'requests',
        'click',
        'emoji',
        'rich'
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
