def newton_solve(f, x0, tol=1e-7, max_iter=1000, h=1e-5):
	"""
	Newton method for finding roots of a function without knowing its derivative.
	
	Parameters:
	f : The function for which we are trying to find a root.
	x0 : Initial guess for the root.
	tol : Tolerance for convergence.
	max_iter : Maximum number of iterations.
	h : Step size for finite difference approximation of the derivative.
	
	"""
	
	iter_count = 0
	err = abs(f(x0))

	while iter_count < max_iter and err > tol:

		# Approximate the derivative using finite differences
		f_prime_x0 = (f(x0 + h) - f(x0)) / h
		if f_prime_x0 == 0:
			raise ValueError("Derivative is zero. No solution found.")

		# Update the guess
		x = x0 - f(x0) / f_prime_x0
		x0 = x	

		err = abs(f(x0))
		iter_count += 1

	if iter_count == max_iter:
		raise ValueError("Newton's method did not converge.")
	
	return x0