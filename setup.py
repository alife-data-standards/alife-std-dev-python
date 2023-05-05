from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(name='ALifeStdDev',
      version='0.3.3',
      description='Python development tools for working with standardized ALife data.',
      url='https://github.com/alife-data-standards/alife-std-dev-python',
      author='Emily Dolson, Alex Lalejini, Matthew Andres Moreno',
      author_email='dolsonem@msu.edu, lalejini@msu.edu, morenoma@umich.edu',
      license='MIT',
      python_requires='>=3.7',
      classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        # python 3.11 dependencies aren't compatible with 3.7
        # wait to sunset 3.7 before formally testing for/supporting 3.11
        # 'Programming Language :: Python :: 3.11',
      ],
      long_description=readme,
      long_description_content_type='text/markdown',
      include_package_data=True,
      keywords='artificial life',
      test_suite='tests',
      packages=find_packages(include=['ALifeStdDev', 'ALifeStdDev.*']),
      install_requires=['networkx', 'pandas'],
      tests_require=['pytest'],
      zip_safe=False,
)
