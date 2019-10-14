import setuptools 

with open("README.md", "r") as filep:
    long_description = filep.read()

setuptools.setup(
    name='Recombination_analysis',
    version='0.1.4',
    description='Analyze recombination events using kmers',
    url='https://github.com/zhuweix/recombination_analysis',
    author='Zhuwei Xu',
    author_email='zhuweix8@gmail.com',
    license='MIT',
    long_description=long_description,
    long_description_content_type="text/markdown",      
    packages=setuptools.find_packages(),
    install_requires=[
        'matplotlib',
        'bitarray',
        'numpy'],
    zip_safe=False,
    python_requires='>=3.6')