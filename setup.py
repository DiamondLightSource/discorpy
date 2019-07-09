from setuptools import setup, find_packages

VERSION=1.1
setup(
    name='vounwarp',
    packages=find_packages(exclude=['test*']),
    version=VERSION,
    include_package_data=True,
    author='Nghia Vo',
    author_email='nghia.vo@diamond.ac.uk',
    description='Distortion correction in Python.',
    keywords=['Distortion correction','Tomography','Imaging'],
    download_url='https://github.com/DiamondLightSource/vounwarp.git',
    license='Apache 2.0',
    platforms='Any',
    classifiers=[
        'Development Status :: 1.1',
        'License :: OSI Approved :: Apache License 2.0',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
	]
)
