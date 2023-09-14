"""python ./full_bubble_model/setup.py install --user"""

from setuptools import setup, find_packages

setup(
    name='full_bubble_model',
    summary='',
    author='Kozák Áron',
    author_email='kozi0223@gmail.com',
    url='https://github.com/hihihi2001/Bubble_dynamics_simulation',
    license='MIT',
    python_requires='>=3.8, <3.12',
    packages=find_packages(),
    install_requires=[
        'matplotlib',
        'pandas',
        'numpy',
        'scipy',
        'termcolor',
        'func_timeout',
        'numba',
    ],
)