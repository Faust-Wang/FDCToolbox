from numpy import floor, ceil

# from Ref [3] https://github.com/stanleybak/AeroBenchVVPython/blob/master/code/util.py
def fix(ele):
    'round towards zero'

    assert isinstance(ele, float)

    if ele > 0:
        rv = int(floor(ele))
    else:
        rv = int(ceil(ele))

    return rv

# from Ref [3] https://github.com/stanleybak/AeroBenchVVPython/blob/master/code/util.py
def sign(ele):
    'sign of a number'

    if ele < 0:
        rv = -1
    elif ele == 0:
        rv = 1  # sign of 0 is positive 
    else:
        rv = 1

    return rv


def print_matrix(matrix):
    print("[")
    for i in matrix:
        print("[", end="")
        for j in i:
            print(f"{j:2f}\t", end="")
        print("]")
    print("]")


def print_array(arr):
    print("[")
    for i in arr:
        print(f"{i},")
    print("]")

