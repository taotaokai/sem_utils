# A python program to sort an array according to the order
# defined by another array
import functools

# function that sorts the first array based on order of
# them in second array
def sortA1ByA2(arr1, arr2):

    # map to store the indices of second array
    # so that we can easily judge the position of two
    # elements in first array
    index = {}

    for i in range(len(arr2)):
        # Consider only first occurrence of element
        if arr2[i] not in index.keys():
            # Assign value of i+1
            index[arr2[i]] = i + 1

    def cmp(a, b):

        # If both a and b are not present in the second array,
        # sort them in ascending order
        if a not in index.keys() and b not in index.keys():
            return a-b

        # If only b is present in the second array, b comes before a
        if a not in index.keys():
            return 1

        # If only a is present in the second array, a comes before b
        if b not in index.keys():
            return -1

        # If both a and b are present in the second array,
        # sort them according to their respective indices
        return index[a]-index[b]

    arr1.sort(key=functools.cmp_to_key(cmp))
