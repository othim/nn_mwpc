from setuptools import setup, Extension

setup(
    name='nn_mwpc',
    version=1.0,
    author='Oliver Thim',
    author_email='oliver.thim@chalmers.se',
    description='Module for acessing the nn_mwpc c++ code.',
    packages=[''],
    package_dir={'':'.'},
    package_data={'':['nn_mwpc.so']},
)
