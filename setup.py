# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 09:41:14 2021


"""
from setuptools import setup
from setuptools import find_packages

with open(file='README.md',mode='r') as readme_handle:
    ATopologyOptimisationInPytonDescript = readme_handle.read()

setup(
      
      name='PyTOpt',
      
      author = 'Daniel Pettersson, Erik Säterskog',
      
      author_email ='danielpettersson974@gmail.com , eric.saterskog@gmail.com',
      
      version = '1.0.0a1',
      
      description = 'A python library used in topology optimisation in 2D created at Chalmers University of Technology',
      
      long_description = ATopologyOptimisationInPytonDescript,
      
      long_description_content_type = 'text/markdown',
      
      url = 'https://github.com/ErikSaterskog/Thesis_Repo',
      
      install_requires=[
          'visvis',
          'PyQt5',
          'pyvtk',
          'calfem-python',
          'sphinx==1.6.6',
          'scipy',
          'numpy',
          ],
      
      keywords = 'optimisation, FEM, optimization, topology',
      
      packages = find_packages(
          include=['Pytopt']
          ),
      
          
      classifiers=[

        
        'Development Status :: Alpha',

        
        'Intended Audience :: Engineering Students',

       
        'Natural Language :: English',

        
        'Programming Language :: Python :: 3'
        ],

   
        python_requires='>=3.7'
)