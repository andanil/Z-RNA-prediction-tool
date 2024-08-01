from collections import defaultdict
import numpy as np


def magnitude(vec):
    return np.sqrt(np.dot(vec, vec))


def vec_distance(vec1, vec2):
    """
    The (euclidean) distance between two points vec1 and vec2.
    This is guaranteed to work for arbitrary but equal dimensions of vectors.
    :param vec1, vec2: A list or np.array of floats
    :returns: A float
    """
    vec1 = np.asarray(vec1)
    vec2 = np.asarray(vec2)
    direction = vec2 - vec1
    return np.sqrt(np.dot(direction, direction))


def inverse_brackets(bracket):
    res = defaultdict(int)
    for i, a in enumerate(bracket):
        res[a] = i
    return res


def dotbracket_to_pairtable(struct):
    if len(struct) == 0:
        raise ValueError("Cannot convert empty structure to pairtable")
    pt = [-1] * len(struct)

    stack = defaultdict(list)
    inverse_bracket_left = inverse_brackets("([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    inverse_bracket_right = inverse_brackets(")]}>abcdefghijklmnopqrstuvwxyz")

    i = 0
    for a in struct:
        if a == ".":
            pt[i] = -1
        else:
            if a in inverse_bracket_left:
                stack[inverse_bracket_left[a]].append(i)
            else:
                assert a in inverse_bracket_right
                if len(stack[inverse_bracket_right[a]]) == 0:
                    raise ValueError('Too many closing brackets!')
                j = stack[inverse_bracket_right[a]].pop()
                pt[i] = j
                pt[j] = i
        i += 1

    if len(stack[inverse_bracket_left[a]]) != 0:
        raise ValueError('Too many opening brackets!')
    return pt


def pairtable_to_tuples(pt):
    '''
    Convert a pairtable to a list of base pair tuples.
    i.e. [4,3,4,1,2] -> [(1,3),(2,4),(3,1),(4,2)]
    :param pt: A pairtable
    :return: A list paired tuples
    '''
    pt = iter(pt)
    tuples = []
    for i, p in enumerate(pt):
        tuples += [(i, p)]
    return tuples


def get_define_len(d):
    val = 0
    for i in range(0, len(d), 2):
        val += d[i + 1] - d[i] + 1
    return val


def any_difference_of_one(stem, bulge):
    """
    See if there's any difference of one between the two
    ends of the stem [(a,b),(c,d)] and a bulge (e,f)
    :param stem: A couple of couples (2 x 2-tuple) indicating the start and end
                 nucleotides of the stem in the form ((s1, e1), (s2, e2))
    :param bulge: A couple (2-tuple) indicating the first and last position
                  of the bulge.
    :return: True if there is an overlap between the stem nucleotides and the
                  bulge nucleotides. False otherwise
    """
    for stem_part in stem:
        for part in stem_part:
            for bulge_part in bulge:
                if abs(bulge_part - part) == 1:
                    return True
    return False


def parens_count(struct):
    return struct.count('(') == struct.count(')')


def dotbracket_to_bp(struct):
    if not parens_count(struct):
        print('Error in input structure.')
        return False

    open_parens = []
    bps = []

    for i, x in enumerate(struct):
        if x == '(':
            open_parens.append(i)
        elif x == ')':
            if len(open_parens) > 0:
                bps.append([open_parens.pop(), i])
            else:
                print('Error in input structure.')
                return False

    return sorted(bps)
