from setuptools import setup, find_packages

setup(
    name="make_prg",
    version="0.3.0",
    packages=find_packages(),
    url="https://github.com/rmcolq/make_prg",
    license="MIT",
    entry_points={"console_scripts": ["make_prg = make_prg.__main__:main"]},
    test_suite="nose.collector",
    tests_require=["nose >= 1.3", "hypothesis >= 4.0"],
    install_requires=[
        "biopython==1.78",
        "numpy==1.20.0",
        "scikit-learn==0.24.1",
        "intervaltree==3.1.0",
        "loguru~=0.5.3",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)
