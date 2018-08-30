import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

# Checking the PyGMQL installation
try:
    import gmql
except ImportError:
    sys.exit("Installation requires: 'gmql'.\n"
             "Use 'pip install gmql' first")
	
setuptools.setup(
    name="genereg",
    version="0.0.1",
    author="Kia23",
    author_email="kia34.r18.dev@gmail.com",
    description="A library for analyzing the regulation systems of target genes belonging to genomic pathways relevant for a specific cancer type.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Kia23/genereg",
    packages=setuptools.find_packages(),
	classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
		'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)