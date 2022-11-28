
##################################################################

from setuptools import setup, find_packages

##################################################################

__authors__ = ["Thomas E.J. Moxham"]
__license__ = "MIT"
__date__ = "28/11/2022"

setup(
    name='atwavpy',
    version='0.0.1',
    description='Simualtion of optical elements.',
    author='Thomas E.J. Moxham',
    author_email='tej.moxham@gmail.com',
    url='https://github.com/tejmoxham/atwavpy/',
    packages=find_packages(),
    install_requires=['numpy', 'scipy'],
    setup_requires=['setuptools']
)

##################################################################