# Christopher Ryba, 2016

# This code calculates the character values \chi_\mu^\lambda of the 
# symmetric groups. Here \lambda and \mu are partitions of the same size,
# represented as either lists or tuples of integers without trailing zeros.
# The function that returns the character value is char_val, so an example
# use might be char_val([4,1],[3,1,1]), which would return 1, indicating
# that the character value of a 3-cycle on the standard representation of
# S_5 is equal to 1. The algorithm is a combination of the 
# Murnaghan-Nakayama rule and the hook-length formula.

# This helper function takes a partition and a skew-hook, and returns
# the partition obtained by removing the skew-hook from the partition.
def rectify(char, skew):
	newchar = sorted(char)[::-1]
	sub = {}
	for i in skew:
		y = i[1]
		if y in sub:
			sub[y] += 1
		else:
			sub[y] = 1
	for i in sub:
		newchar[i] -= sub[i]
	ctr = len(newchar) - 1
	while newchar[ctr] == 0 and ctr >= 0:
		ctr -= 1
	return newchar[:ctr+1]

# For the sake of efficiency, calculated character values are stored in a
# dictionary so they do not need to be recomputed later ("memoisation").
# Initially we only have the base case of the trivial group S_0.
chardict = {}
chardict[( tuple([]) , tuple([]) )] = 1

# This helper function finds the dual (or transpose) partition.
def dual(part):
	out = []
	cursor = len(part)-1
	for i in xrange(1,part[0]+1):
		while part[cursor] < i:
			cursor -= 1
		out.append(cursor+1)
	return out

# Given a partition, this function returns the dimension of the corresponding
# irreducible representation of the symmetric group using the hook-length formula.
from math import factorial
def hook_formula(part):
	prod = 1
	loc = dual(part)
	for i in xrange(len(part)):
		for j in xrange(part[i]):
			prod *= (part[i] - j+loc[j]-i-1)
	return factorial(sum(part))/prod


# This function calculuates character values of the symmetric groups for the
# irreducible representation char at an element of cycle type elt.
def char_val(char, elt):
	# Check whether the calculation has already been done.
	key = (tuple(char), tuple(elt))
	if key in chardict:
		return chardict[key]

	# Check whether we are asked to find the dimension of an irreducible
	# representation, in which case we use the hook-length formula, and
	# record the answer before returning.
	if elt == [1]*len(elt):
		value = hook_formula(char)
		chardict[key] = value
		return value

	# We begin the recursion of the Murnaghan-Nakayama rule. The list
	# borderstrip contains the coordinates of all border cells in the
	# diagram of the partition char. These are candidate locations for
	# rim-hooks. The variable value stores the progressive value of
	# the character.
	borderstrip = []
	value = 0

	# We zig-zag our way from the bottom-left cell of the diagram (in
	# English notation) to the top-right, going up a row whenever we
	# reach the end of a row.
	x, y = (0, len(char)-1)
	for i in xrange(len(char) + char[0] - 1): 
		borderstrip.append((x,y))
		if char[y] == x + 1:
			y -= 1
		else:
			x += 1

	# Each skew-hook of size elt[0] (the largest part of the cycle type)
	# must begin somewhere on the border strup, we consider each possible
	# location and check whether removing elt[0] contiguous blocks results
	# in a valid partition.
	for i in xrange(len(borderstrip)-elt[0]+1):
		start = borderstrip[i]
		end = borderstrip[i+elt[0]-1]
		# Check validity of resulting diagram/partition. If invalid, skip.
		if i > 0:
			prev = borderstrip[i-1]
			if prev[0] == start[0]:
				continue
		if i + elt[0] - 1 < len(borderstrip) - 1:
			next = borderstrip[i+elt[0]]
			if next[1] == end[1]:
				continue
		# Recursive step; calculates character value in a smaller instance.
		augment = char_val(rectify(char, borderstrip[i:i+elt[0]]), elt[1:])
		# Add to the character value with the correct sign.
		if (start[1] - end[1]) % 2 == 1:
			value -= augment
		else:
			value += augment
	# Add to memoisation table before returning.
	chardict[key] = value
	return value