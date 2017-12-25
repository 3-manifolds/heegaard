"""
heegaard: A Python interface to Berge's Heegaard program, or more
precisely Marc Culler's Unix port of it.  

Provides the single function 'is_realizable' which checks if the given
words (representing relators of a group) can be realized on the
boundary of the appropriate handlebody.

Python interface written by Nathan Dunfield (nathan@dunfield.info).
This code is in the public domain.  
"""

__version__ = '0.9'

def version():
    return __version__

from .heegaard import is_realizable

__all__ = ["is_realizable"]
