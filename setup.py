from setuptools import setup

with open("README.md", encoding="utf-8", errors="ignore") as fp:
    long_description = fp.read()

setup(
    name="Eugene",
    version="1.0.0",
    description="A library models to experiment with breeding programs",
    author="Kelvin Davis",
    author_email="Kelvin.Davis1@monash.edu",
    python_requires=">=3.7",
    packages=["eugene"],
    install_requires=[
        "bitarray"
        "gurobipy=10.0.1",
        "minizinc",
        "eugene_rs",
    ]
)
