from setuptools import setup, Extension

setup(
    name='example',
    version=0.1,
    author='Oliver Thim',
    author_email='oliver.thim@chalmers.se',
    description='Just a test',
    packages=[''],
    package_dir={'':'.'},
    package_data={'':['example.so']},
)
