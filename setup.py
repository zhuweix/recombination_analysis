from setuptools import setup

with open("README.md", "r") as filep:
    long_description = filep.read()

setup(
    name='Recombination_analysis',
    version='0.1',
    description='Analyze recombination events using kmers',
    url='https://github.com/zhuweix/recombination_analysis',
    author='Zhuwei Xu',
    author_email='zhuweix8@gmail.com',
    license='MIT',
    long_description=long_description,
    long_description_content_type="text/markdown",      
    packages=setuptools.find_packages(),
    zip_safe=False)