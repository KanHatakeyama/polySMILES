from setuptools import setup, find_packages
import sys

sys.path.append('./polysmiles')

setup(name='polysmiles',
        version='2021.7.03',
        description='PolySMILES',
        long_description="README",
        author='Kan Hatakeyama',
        license="MIT",
        packages = find_packages(),
        #packages = ["PolySMILES"],
        #package_dir = {'': 'polysmiles'}
        #package_dir = {'PolySMILES': 'polysmiles'}
    )