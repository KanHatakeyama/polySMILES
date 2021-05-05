from setuptools import setup, find_packages
import sys

sys.path.append('polysmiles')

setup(name='polysmiles',
        version='2021.5.04',
        description='PolySMILES',
        long_description=README,
        author='Kan Hatakeyama',
        license=license,
        packages = ["PolySMILES"],
        package_dir = {'PolySMILES': 'polysmiles'}
    )