from __future__ import print_function
from __future__ import absolute_import
from __future__ import division


__all__ = [
    'are_items_in_list'
    'list_split',
    'sublist_from_to_items_in_closed_list'
]


def are_items_in_list(items, in_list):
    """Check if items are in a list.

    Parameters
    ----------
    items : list
            A list of items (order does not matter).
    in_list : list
            A list.

    Returns
    -------
    bool
            True if all items are in the list. False otherwise.
    """

    for i in items:
        if i not in in_list:
            return False
    return True


def list_split(in_list, indices):
    """Split list at given indices.
    Closed lists have the same first and last elements.
    If the list is closed, splitting wraps around if the first or last index is not in the indices to split.


    Parameters
    ----------
    in_list : list
            A list.
    indices : list
            A list of indices to split.

    Returns
    -------
    split_lists : list
            Nest lists from splitting the list at the given indices.

    """

    n = len(in_list)

    if in_list[0] == in_list[-1]:
        closed = True
        if n - 1 in indices:
            indices.remove(n - 1)
            if 0 not in indices:
                indices.append(0)
    else:
        closed = False

    indices = list(sorted(set(indices)))

    split_lists = []
    current_list = []
    for index, item in enumerate(in_list):
        current_list.append(item)
        if (index in indices and index != 0) or index == n - 1:
            split_lists.append(current_list)
            current_list = [item]

    if closed:
        if 0 not in indices:
            start = split_lists.pop(0)[1:]
            split_lists[-1] += start

    return split_lists


def sublist_from_to_items_in_closed_list(l, from_item, to_item):
    """Return sublist between oe item to another.

    Parameters
    ----------
    l : list
            A list.
    from_item
            An item to be found in the list. The beginning of the sublist.
    to_item
            An item to be found in the list. The end of the sublist.

    Returns
    -------
    sublist : list
            A sublist from the input list, between from_item and to_item.
    """

    if from_item == to_item:
        return [from_item]
    if l[0] != l[-1]:
        l.append(l[0])
    from_idx = l.index(from_item)
    to_idx = l.index(to_item)
    sublists = list_split(l, [from_idx, to_idx])

    for sublist in sublists:
        if sublist[0] == from_item:
            return sublist