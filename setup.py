import setuptools
import versioneer

with open("README.rst", "r") as fh:
    long_description = fh.read()

exec(open('CONGA/version.py').read())

setuptools.setup(
    name="CONGA",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
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
    python_requires='>=3.9.0, <=3.9.16',
    install_requires=[
        'pyteomics',
        'pandas',
        'numpy',
        'scipy',
        'tdqm',
        'matplotlib'
    ],
)