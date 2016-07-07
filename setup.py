try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Time-domain FMR analysis tool',
    'author': 'Alex Devonport',
    'url': 'https://github.com/alexdevonport/pypimm',
    'author_email': 'alex.devonport@asu.edu',
    'version': '0.1',
    'install_requires': ['scipy', 'numpy', 'matplotlib', 'nose'],
    'packages': ['pypimm', 'pypimm.fitutil'],
    'scripts': [],
    'name': 'pypimm'
}

setup(**config)
