"""
heegaard: A Python interface to Berge's Heegaard program, or more
precisely Marc Culler's Unix port of it.  

Provides the function 'is_realizable' which checks if the given
words (representing relators of a group) can be realized on the
boundary of the appropriate handlebody.

You'll need the additional 'pexpect' library.  Probably your system
has a package for this, but if not you can get it here

http://www.noah.org/wiki/Pexpect

or better yet install 'easy_install' from

http://peak.telecommunity.com/DevCenter/EasyInstall

and then do 'easy_install pexpect' at the command line.

Finally, if the 'heegaard' binary is not in your path, you'll need to
edit the variable 'heegaard_program' below.

To test if everything is working fine, do 'python heegaard.py'.  

Written by Nathan Dunfield (nathan@dunfield.info).  This code is in
the public domain.
"""

import pexpect, re

# If the 'heegaard' binary is not in your path, edit the below to
# reflect it's full location, e.g. "/Users/dunfield/work/heegaard/heegaard"

heegaard_program = "heegaard"

def start_heegaard(relations):
    H = pexpect.spawn(heegaard_program)
    H.expect("HIT")
    H.send("k")
    for rel in relations:
        H.send(rel + "\n")
    H.send("\n\n")
    return H

def is_realizable(relations, maxtime = 10, full_answer=False):
    """
    Tests to see if the given words can be realized as curves on the
    boundary of the appropriate handlebody.

    Usage examples:

    >>> is_realizable(['aabbAbbaabbABAAABAbb', 'aabbaabbbaabbABAbbaabbbaabbaabbb'])
    True
    >>> is_realizable(['bbaaccabc'])
    False

    Despite my best efforts, every once and awhile heegaard get
    confused at the input and hang; not sure why. In this case,
    check_realizability throws a RuntimeError; the optional 'maxtime'
    argument is how long it waits before deciding this has happened
    (default is 10 seconds).

    Finally, if the 'full_answer' flag is set, the function returns
    the pair (answer, heegaard output).
    """
    
    H = start_heegaard(relations)
    H.timeout = maxtime
    for command in ["c", "y", "n", "\n", "\n", "y", "q", "\n", "y", "q", "\n", "y", "q", "\n", "y", "q"]:
        H.send(command)
    try:
        raw_data = H.read()
        H.close()
        ans = re.search("is not realizable",  raw_data) == None
        return (ans, raw_data) if full_answer else ans
    except pexpect.TIMEOUT:
        raise RuntimeError

# Code for testing

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

__all__ = ["is_realizable"]
