import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="findres",
    version="1.0",
    author="Monica Sugan, Stefano Campanella",
    author_email="msugan@inogs.it, scampanella@inogs.it",
    description="Repeating earthquake discovery using cross-correlation and differential arrival times",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/msugan/pyres",
    packages=setuptools.find_packages(),
    scripts=['bin/findres'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'numpy',
        'tqdm',
        'obspy',
        'pyyaml',
        'multitaper'
    ],
    python_requires='~=3.8',
)
