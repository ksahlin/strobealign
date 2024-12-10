from skbuild import setup
import nanobind

setup(
    name="strobealign",
    description="Python bindings for strobealign",
    license="MIT",
    version="0.15.0",
    packages=["strobealign"],
    package_dir={"": "src/python"},
    cmake_install_dir="src/python",
    cmake_args=["-DPYTHON_BINDINGS=ON"],
)
