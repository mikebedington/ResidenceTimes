from setuptools import setup

setup(name='residence_time',
        version='0.1.0',
        description='Calculating different residence times from PyLAG output',
        url='https://gitlab.ecosystem-modelling.pml.ac.uk/mbe/bodc_tide_db',
        author='mbe',
        author_email='mbe@pml.ac.uk',
        packages=['residence_time'],
        #install_requires=['numpy', 'sqlite3', 'datetime', 'subprocess', 'gpxpy.geo'],
        zip_safe=False)

