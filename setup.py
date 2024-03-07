from setuptools import setup, find_packages
import os

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join(path, filename))
    return paths

with open('README.md') as f:
    long_description = f.read()

setup(
    name="bridgerna2024",
    version='0.0.1_dev',
    description='Code and workflows accompanying Durrant and Perry et al. 2024',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/hsulab-arc/BridgeRNA2024',
    author="Matt Durrant",
    author_email="matthew@arcinstitute.org",
    license="MIT",
    packages=find_packages(),
    package_data={
        "bridgerna2024": package_files("bin/linux/") + package_files("snakemake/") + package_files("scripts/"),
    },
    include_package_data=True,
    install_requires=[
        'numpy',
        'click',
        'biopython',
        'tqdm',
        'forgi',
        'pysam'
    ],
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'brna2024 = bridgerna2024.main:cli'
        ]
}
)