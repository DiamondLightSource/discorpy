import pathlib
import setuptools
import sys

py_ver = sys.version.split(".")[:2]
python_version = py_ver[0] + "." + py_ver[1]

if python_version <= "3.7":
    dependencies = [
        "h5py",
        "pillow",
        "matplotlib<3.6",
        "scikit-image<0.18",
        "scipy<=1.7",
        "numpy<1.22"
    ]
else:
    dependencies = [
        "h5py",
        "pillow",
        "matplotlib",
        "scikit-image",
        "scipy",
        "numpy"
    ]

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

setuptools.setup(
    name="discorpy",
    version="1.4",
    author="Nghia Vo",
    author_email="nghia.vo@diamond.ac.uk",
    description='Correction for radial distortion and perspective distortion '
                'in Python',
    long_description=README,
    long_description_content_type="text/markdown",
    keywords=['Distortion correction', 'Tomography', 'Radial lens distortion',
              'Camera calibration', 'Perspective distortion'],
    url="https://github.com/DiamondLightSource/discorpy",
    download_url='https://github.com/DiamondLightSource/discorpy.git',
    license="Apache 2.0",
    platforms="Any",
    packages=setuptools.find_packages(include=["discorpy", "discorpy.*"],
                                      exclude=['test*', 'doc*', 'data*',
                                               'example*']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Image Processing"
    ],
    install_requires= dependencies,
    python_requires='>=3.7',
)
