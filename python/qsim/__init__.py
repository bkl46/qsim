"""
qsim - pimc c++ implemntation with python interface
 
"""

try:
    from .qsim_py import *  
except ImportError as e:
    raise ImportError("qsim module not found.") from e


__version__ = "0.1.0"

