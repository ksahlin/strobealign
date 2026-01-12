from skbuild import setup
import nanobind

setup(
    name="strobealign",
    description="Python bindings for strobealign",
    license="MIT",
    version="0.17.0",
    packages=["strobealign"],
    package_dir={"": "cpp/python"},
    cmake_install_dir="cpp/python",
    cmake_args=["-DPYTHON_BINDINGS=ON"],
)
