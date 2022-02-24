import setuptools

setuptools.setup(
    name="testpack",
    version = "0.0.1",
    packages=['testpack'],
    package_dir = {'testpack':'src'},
    python_requires=">=3.6",
)
