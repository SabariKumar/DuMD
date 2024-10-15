from setuptools import setup, find_packages

setup(
    name='DuMD',
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    version='0.0.1',
    description='Mindlessly run OpenMM simulations',
    author='Sabari Kumar',
    author_email='sabarik@colostate.edu',
)
