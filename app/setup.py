from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name = 'krepe',
    version = '0.2.0',
    author = 'Aengus McGuinness, and Erika Pedersen',
    author_email = 'erikalaraine2@gmail.com',
    license = 'GPLv3',
    description = 'A package to count k-mers and calculate/plot De Bruijn graphs, occurence bar graphs, venn diagrams, and dendrograms.',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = 'https://github.com/RGBwalnut/Kmer-Counting-Analysis',
    py_modules = ['krepe'],
    packages = find_packages(),
    install_requires = [
        'numpy',
        'toyplot',
        'matplotlib',
        'psutil',
        'matplotlib_venn',
        'sourmash'
    ],
    python_requires='>=3',
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    #scripts=['src/genome-visualization', 'src/organism-comparison', 'src/compare-all']
    entry_points={'console_scripts': ['genome-visualization=krepe:genome_visualization','organism-comparison=krepe:organism_comparison', 'compare-all=krepe:compare_all',]}
     #entry_points={'console_scripts': ['genome-visualization=genomeVisualization:main','organism-comparison=organismComparison.py:main', 'compare-all=compareAll:main',]}
)

