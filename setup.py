import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lsfmpy",
    version="0.0.1",
    author="Sylwia Bednarek, Piotr Majka",
    author_email="s.bednarek@nencki.gov.pl",
    description="HDF5-based software for storing and managing voluminous 3D imaging data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Neuroinflab/lsfm_image_server",
    licence="GPLv3",
    packages=setuptools.find_packages(),
    scripts=['bin/lsfmpy', 'bin/dump_metadata'],
    classifiers=[
        "Development Status :: 3 - Alpha"
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Natural Language :: English",
        "Environment :: Console",
    ],
    install_requires=[
        "SimpleITK",
        "enum34",
        "fire",
        "h5py",
        "h5py_cache",
        "imageio",
        "lxml",
        "nibabel",
        "numpy",
        "scipy",
        "tifffile==0.12.1"
        ]
)
