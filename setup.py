from setuptools import setup, find_packages
import os


def get_version():
    with open(os.path.join("tmad", "VERSION"), "r") as f:
        return f.readline().strip()


def get_long_description():
    with open("README.md", "r") as f:
        description = f.read()
    return description


setup(
    name="tmad",
    packages=find_packages(),
    python_requires=">=3.8",
    description="Xenium TMA separator script",
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    version=get_version(),
    author="Luke Zhang",
    author_email="luke.zhang@stjude.org",
    include_package_data=True,
    install_requires=[
        "pandas",
        "numpy",
        "scipy",
        "click"
    ],
    entry_points={
        "console_scripts": [
            "tmad=tmad.tmad:main"
        ]
    }
)
