from typing import List
from itertools import groupby


def remove_duplicated_consecutive_elems_from_list(the_list: List) -> List:
    return [elem[0] for elem in groupby(the_list)]