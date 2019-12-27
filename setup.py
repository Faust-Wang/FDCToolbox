from setuptools import setup

def readme():
    with open('README.txt') as f:
        return f.read()

setup(
    name='FDCToolbox',
    version='0.1',
    packages=['FDCToolbox'],
    long_description=readme(),
    long_description_content_type='text',
    url='',
    author='Parth Viradiya',
    author_email='parthviradiya08@gmail.com',
    license='MIT',
)