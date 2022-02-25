import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

version = '0'

setuptools.setup(
    name='counting_molecules',
    version=version,
    author='Jack Hellerstedt and Tyler Hennen',
    author_email='hellerstedt.jack@gmail.com',
    description='Library to count molecules in STM images',
    long_description=long_description,
    url='https://github.com/thennen/counting-molecules',
    project_urls = {
        "Bug Tracker": "https://github.com/thennen/counting-molecules/issues"
    },
    license='MIT',
    packages=['counting_molecules', 'pairwise_chirality'],
    install_requires=['numpy', 'matplotlib', 'matplotlib_scalebar', 'nanonispy', 'mahotas', 'scipy', 'sklearn', 'skimage'],
)