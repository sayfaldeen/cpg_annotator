from setuptools import setup, find_packages

setup(
    name="cpg_annotator",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "polars>=0.20.0",
        "requests>=2.28.0",
    ],
    python_requires=">=3.8",
    author="Your Name",
    author_email="your.email@example.com",
    description="A tool for annotating CpG sites with gene information",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/cpg-annotator",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={
        "console_scripts": [
            "cpg-annotator=cpg_annotator.annotator:main",
        ],
    },
) 