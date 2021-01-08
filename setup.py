import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="stochastic_cosmo", 
    version="0.0.1",
    author="Amanda Farah",
    author_email="afarah@uchicago.edu",
    description="Package with some functions for cosmological calculations using the stochastic GW background",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/afarah18/stochastic_cosmo",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
