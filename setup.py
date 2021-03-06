import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="commmodelpy",
    version="0.0.3",
    author="Paulocracy",
    author_email="bekiaris@mpi-magdeburg.mpg.de",
    description="The commmodelpy package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ARB-Lab/commmodelpy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    install_requires=[
        "cobra",
    ],
)
