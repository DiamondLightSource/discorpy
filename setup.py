import setuptools

setuptools.setup(
    name="vounwarp",
    version="1.2",
    author="Nghia Vo",
    author_email="nghia.vo@diamond.ac.uk",
    description='Radial lens distortion correction in Python',
    keywords=['Distortion correction', 'Tomography', 'Radial lens distortion'],
    url="https://github.com/nghia-vo/vounwarp",
    download_url='https://github.com/nghia-vo/vounwarp.git',
    license="Apache 2.0",
    platforms="Any",
    packages=setuptools.find_packages(include=["vounwarp", "vounwarp.*"],
                                      exclude=['test*', 'doc*', 'data*',
                                               'example*']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Distortion correction"
    ],
    install_requires=[
        "matplotlib",
        "numpy",
        "scipy",
        "scikit-image",
        "h5py",
        "pillow"
    ],
)