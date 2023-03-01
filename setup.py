from setuptools import setup

setup(name='ALifeStdDev',
      version='0.1',
      description='Python development tools for working with standardized ALife data.',
      url='https://github.com/alife-data-standards/alife-std-dev-python',
      author='Emily Dolson, Alex Lalejini, Matthew Andres Moreno',
      author_email='dolsonem@msu.edu, lalejini@msu.edu, morenoma@umich.edu',
      license='MIT',
      packages=['ALifeStdDev'],
      install_requires=['networkx', 'pandas'],
      zip_safe=False,
)
