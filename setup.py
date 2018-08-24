import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open(os.path.join(".", "genexpreg_cancer", "resources", "version"), "r") as f_ver:
    version = f_ver.read().strip()

with open(os.path.join(".", "genexpreg_cancer", "resources", "github_url"), "r") as f_ver:
	github_url = f_ver.read().strip()

# Checking the PyGMQL installation
try:
    import gmql
except ImportError:
    sys.exit("Installation requires: 'gmql'.\n"
             "Use 'pip install gmql' first")
	
setuptools.setup(
    name="genexpreg_cancer",
    version=version,
    author="Kia23",
    author_email="kia34.r18.dev@gmail.com",
    description="A library for analyzing the regulation systems of target genes belonging to genomic pathways relevant for a specific cancer type.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=github_url,
    packages=setuptools.find_packages(),
	install_requires=['pandas', 'numpy'],
	classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
		'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)