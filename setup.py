import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lsfm_image_server",
    version="0.0.1",
    author="Sylwia Bednarek, Piotr Majka",
    author_email="s.bednarek@nencki.gov.pl",
    description="HDF5-based software for storing and managing voluminous 3D imaging data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research"
    ],
)