from setuptools import setup, find_packages
#from Cython.Build import cythonize


setup(name='cLoops',
    version='0.5',
    author='Yaqiang Cao',
    author_email='caoyaqiang@picb.ac.cn',
    url='https://github.com/YaqiangCao/cLoops',
    description='Loops calling for ChIA-PET and HiChIP data.',
    classifiers=[
        'Environment :: Console',
        'Programming Language :: Python :: 2.7',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    packages=find_packages(exclude=['tests','docs']),
    long_description=open('README').read(),
    #setup_requires=["joblib","numpy","seaborn","pandas","scipy","HTSeq"],
    entry_points={
        'console_scripts': [
            'cLoops=cLoops.pipe:main',
                ],
        },
    scripts = ["scripts/deLoops","scripts/jd2juice","scripts/jd2washU","scripts/jd2saturation","scripts/jd2fingerprint"],

    )
