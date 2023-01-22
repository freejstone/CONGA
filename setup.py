import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CONGA",
    version="1.0.0",
    author="Jack Freestone",
    author_email="jfre0619@uni.sydney.edu.au",
    description="Combined Open and Narrow searches via Group Analysis",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/freejstone/CONGA",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)