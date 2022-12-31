from setuptools import setup,find_packages
from os.path import splitext
from os.path import basename
from glob import glob


with open('README.md') as f:
    readme = f.read()

setup(name='STREET',
      version='0.0.3',
      description='SpaTial Representation Error EstimaTor',
      long_description=readme,
      long_description_content_type='text/markdown',
      author='Amir Souri',
      author_email='ahsouri@gmail.com',
      license='MIT',
      packages=['street'],
      install_requires=[
          'numpy','matplotlib','scipy','netCDF4','scikit-gstat'
      ],
      zip_safe=False)
