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

Finally, you need to build the extra program 'heegaard_is_realizable'
by doing 'make all' in the 'heegaard/src' directory.  Moreover, if
'heegaard_is_realizable' is not in your path, you'll need to edit the
variable 'heegaard_program' below.

To test if everything is working fine, do 'python heegaard.py'.  

Written by Nathan Dunfield (nathan@dunfield.info).  This code is in
the public domain.
"""

import pexpect, re, os, tempfile

# If the 'heegaard' binary is not in your path, edit the below to
# reflect it's full location, e.g. "/Users/dunfield/work/heegaard/heegaard_is_realizable"

heegaard_program = "heegaard_is_realizable"

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

    temp_file_name = tempfile.mktemp()
    tempf = open(temp_file_name, "w")
    tempf.write("0\n")
    for rel in relations:
        tempf.write("\t" + rel + "\n")
    tempf.write("\n")
    tempf.close()

    H = pexpect.spawn(heegaard_program + "  " + temp_file_name)
    H.timeout = maxtime
    try:
        raw_data = H.read()
        H.close()
        ans = re.search("HEEGAARD_ANS:\s*\(([A-Z]+)\)", raw_data).groups(1)
    except pexpect.TIMEOUT:
        raise RuntimeError, "Heegaard timedout"

    if not ans in [("NO",), ("YES",)]:
        raise RuntimeError, "Heegaard failed" 
    ans = ans[0] == "YES"
    os.remove(temp_file_name)
    return (ans, raw_data) if full_answer else ans


# Code for testing

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)

__all__ = ["is_realizable"]
