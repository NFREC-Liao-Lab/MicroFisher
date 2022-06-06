import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="MicroFisher",
    version="0.6",
    author=["Steven Wu", "Haihua Wang"],
    # author_email="author@example.com",
    description="MicroFisher pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NFREC-Liao-Lab/MicroFisher-Fungal-Profiling/",
    project_urls={
        "Bug Tracker": "https://github.com/NFREC-Liao-Lab/MicroFisher-Fungal-Profiling/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        # "License :: OSI Approved :: MIT License",
        # "Operating System :: OS Independent",
    ],
    keywords=["bioinformatics", "centrifuge"],
    package_dir={"": "."},
    packages=setuptools.find_packages(where="."),
    # python_requires=">=3.6",
)
