# Christopher Ryba, 2018

# This code calculates the Littlewood-Richardson coefficient
# c_{p1, p2}^{p3} associated to the partitions p1, p2, p3.
# The paritions should be input as lists or tuples of their 
# constituent parts without trailing zeros, e.g. [3,2,1] or (4,3,1).
# The algorithm uses the Littlewood-Richardson rule, which states
# that c_{p1, p2}^{p3} is equal to the number of semistandard
# tableaux of (skew) shape p3/p2 and weight p1 satisfying the
# lattice word property.

def LR_coeff(p1, p2, p3):
	# Check that the partitions have compatible sizes.
	assert sum(p1) + sum(p2) == sum(p3)
	# Check that the diagram of p1 fits inside the diagram of p3.
	for index, row in enumerate(p1):
		if index >= len(p3) or row > p3[index]:
			return 0
	# Check that the diagram of p2 fits inside the diagram of p3.
	for index, row in enumerate(p2):
		if index >= len(p3) or row > p3[index]:
			return 0
	# If |p1| is larger than |p2|, swap p1 and p2, as we'll be
	# constructing lattice words of size |p1|, so we expect it to be 
	# faster this way (the LR coefficient is symmetric in p1, p2).
	if sum(p2) < sum(p1):
		p1, p2 = p2, p1

	# We will construct the lattice words appearing in the 
	# Littlewood-Richardson rule. For that purpose we store which cell
	# (i.e. (row, column) coordinate pair) each entry in the lattice
	# word corresponds to, and vice versa. This allows us to check the
	# semistandard tableau property easily.
	entries_to_coords = {}
	coords_to_entries = {}

	# The variable ctr indicates the position in the lattice word. We
	# iterate over the rows in the skew diagram of p3/p2, and consider
	# the cells in reverse order to (i.e. in the order they appear in
	# the associated lattice word).
	ctr = 0
	for row in xrange(len(p3)):
		# offset is the number of cells we skip in the current row
		# by removing p2 from p3.
		offset = 0
		if row < len(p2):
			offset = p2[row]
		# We consider the squares from right to left, record their row
		# and column, as well as their position in the lattice word.
		for tmp in xrange(p3[row] - offset):
			entries_to_coords[ctr] = (row+1, p3[row] - tmp)
			coords_to_entries[(row+1, p3[row] - tmp)] = ctr
			ctr += 1

	# lattice_word stores the lattice word we've (partially) constructed.
	# We use 0-indexed entries for convenience.
	# weight_count tracks the letters we've used for the purpose of
	# checking the lattice word condition; the i-th entry counts the
	# number of times the entry i has occurred so far.
	lattice_word = [None] * sum(p1)
	weight_count = [0] * len(p1)

	# We recursively consider all possible ways lattice words. We recurse
	# on the location in the lattice word, keeping track of the number of
	# solutions we've found so far.
	def recursive_fill(location, count):
		# If the location is past the end of the word, we have a solution;
		# this is the base case of the recursion; increment the number of
		# solutions and return this.
		if location >= sum(p1):
			# If you wanted the lattice words themselves (from which you
			# could reconstruct the associated tableaux), you could print
			# the variable tableau at this point.
			return count + 1
		# We find all possible values of the next entry in the lattive word.
		# This is subject to three things: (1) the lattice word condition,
		# (2) the semistandard tableau condition, (3) having weight p1.
		viable = []
		# The possible values are bounded by condition (2) if there are
		# adjacent cells to the current one. We calculate these bounds.
		lower_bound = 0
		upper_bound = len(p1) - 1
		row, column = entries_to_coords[location]
		if (row, column + 1) in coords_to_entries:
			upper_bound = lattice_word[coords_to_entries[(row, column + 1)]]
		if (row - 1, column) in coords_to_entries:
			lower_bound = lattice_word[coords_to_entries[(row - 1, column)]] + 1
		# We now check whether adding a possible value is compatible with
		# (firstly) condition (1), and (secondly) condition (2).
		for candidate in xrange(lower_bound, upper_bound + 1):
			if weight_count[candidate] == p1[candidate]:
				continue
			if candidate > 0 and weight_count[candidate] == weight_count[candidate - 1]:
				continue
			# If no condition is voilated, we record the value.
			viable.append(candidate)
		# We now recurse over all possible values by writing the value to
		# lattice_word, and adjusting weight_count accordingly. We proceed
		# to the next letter in the lattice word, keeping track of the number
		# of solutions we found along the way.
		for option in viable:
			lattice_word[location] = option
			weight_count[option] += 1
			count = recursive_fill(location + 1, count)
			weight_count[option] -= 1
		# Finally, return the number of solutions
		return count
	# Use the recursion to count lattice words, starting at the beginning,
	# with an initial solution count of zero.
	return recursive_fill(0,0)

