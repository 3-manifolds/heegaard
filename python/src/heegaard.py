import pexpect, re, os, tempfile
from .bin import __path__ as binary_directory
heegaard_program = os.path.join(binary_directory[0],
                                'heegaard_is_realizable')


def is_realizable(relations, maxtime = 10, full_answer=False):
    """
    Test to see if the given words can be realized as curves on the
    boundary of the appropriate handlebody.

    Usage examples:

    >>> is_realizable(['aabbAbbaabbABAAABAbb', 'aabbaabbbaabbABAbbaabbbaabbaabbb'])
    True
    >>> is_realizable(['bbaaccabc'])
    False
    >>> is_realizable(['AABCab','ACCBcb'])
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
    with open(temp_file_name, "w") as tempf:
        tempf.write("0\n")
        for rel in relations:
            tempf.write("\t" + rel + "\n")
        tempf.write("\n")

    H = pexpect.spawn(heegaard_program + "  " + temp_file_name)
    H.timeout = maxtime
    try:
        raw_data = str(H.read())
        H.close()
        ans = re.search("HEEGAARD_ANS:\s*\(([A-Z]+)\)", raw_data).groups(1)
    except pexpect.TIMEOUT:
        raise RuntimeError("Heegaard timedout")

    if ans not in [("NO",), ("YES",)]:
        raise RuntimeError("Heegaard failed")
    ans = ans[0] == "YES"
    os.remove(temp_file_name)
    return (ans, raw_data) if full_answer else ans
