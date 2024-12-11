import sys
from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext

# Define the extension
ext_modules = [
    Pybind11Extension(
        "molsim._molsim",
        sources=[
            "src/pybindings.cpp",
            "src/hardDisks/hardDisks.cpp",
            "src/molecularDynamics/md.cpp",
            "src/molecularDynamics/sample.cpp",
            "src/molecularDynamics/thermostats.cpp",
            "src/monteCarlo/mc.cpp",
            "src/monteCarlo/moves.cpp",
        ],
        include_dirs=[
            "src/hardDisks",
            "src/molecularDynamics",
            "src/monteCarlo",
            "src/utils",
        ],
        cxx_std=20,  # Requires C++20 as per original CMake
    )
]

setup(
    name="molsim",
    version="1.0.0",
    description="Code accompanying exercises for the Molecular Simulation course.",
    author="Drs. Youri Ran, Drs. Rik Breebaart, Dr. Jocelyne Vreede, Dr. David Dubbeldam",
    author_email="d.dubbeldam@uva.nl",
    python_requires=">=3.10",
    install_requires=[
        "pybind11",
        "numpy",
        "matplotlib",
        "scipy",
    ],
    packages=find_packages(),
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
