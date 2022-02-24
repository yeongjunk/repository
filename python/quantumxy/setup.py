import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="quantumxy",
    version="0.0.1",
    author="Example Author",
    author_email="yeongjun@ust.ac.kr",
    description="Tools for quantum XY chain with periodic boundary condition",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/xinzhao-jungle/",
    project_urls={
        "Bug Tracker": "https://github.com/xinzhao-jungle/repository/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir = {'' : 'src'},
    packages= [''],
    python_requires=">=3.6",
)
