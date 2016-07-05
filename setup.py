from setuptools import setup

setup(name='pyradsim',
      version='0.1',
      description='Simple polarimetric radar simulator',
      url='http://gitlab.epfl.ch/wolfensb/radar_simulator/',
      author='Daniel Wolfensberger - LTE EPFL',
      author_email='daniel.wolfensberger@epfl.ch',
      license='MIT',
      packages=['pyradsim'],
      package_data={'pyradsim': ['.defaults.yml']},
      include_package_data=True,
      install_requires=[
          'pyyaml',
          'numpy',
          'scipy',
	  'pytmatrix',
      ],
      zip_safe=False)

