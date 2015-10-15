# # Splits first argument, vector v into v1 and v2
# # where v1 is all the elements of v that are in a
# # and v2 is all the elements of v that are in b
# # Returns error if:
# # there exists an element in v that is in a and b
# # there exists an element in v that is in a (or b) more than once
# # there exists an element in v that is in neither a nor b
# SplitVecByMembership <- function(v, a, b) {
# 
# 	a.idx <- match(v, a, nomatch = 0)
# 	b.idx <- match(v, b, nomatch = 0)
# 
# 	# three kinds of errors to catch
# 	# i gave up implementing this one.  TODO
# # 	if () {
# # 		stop(paste('An element in the vector to split',
# # 								'occurs more than once in a subvector'))
# # 	}
# 
# 	if (any(a.idx & b.idx)) {
# 		stop(paste('An element of the vector to split occurs',
# 								'in more than one subvector'))
# 	}
# 
# 	if (any(!(a.idx | b.idx))) {
# 		stop(paste('An element of the vector to split',
# 								'occurs in neither subvector'))
# 	}
# 
# 	# if no errors, do the split and return
# 	# these are the v1 and v2 referred to in the function description
# 	return(list(a[a.idx], b[b.idx]))
# }
