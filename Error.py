""" Python Code for Error Analysis Homework """

import numpy as np 
import math 
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit

def one_normal(n): # n = 5,10,50,100,1000
	measurements = np.random.standard_normal(n)
	mean = measurements.mean(axis = 0)
	std = measurements.std(axis = 0)
	mean_error = std / math.sqrt(n)
	print("Mean: " + str(mean))
	print("Standard Deviation: " + str(std))
	print("Error in the Mean: " + str(mean_error))

def normal(m,n): # m = 1000, n = 5,10,50,100,1000

	" Calculate the stats for the measurements per experiment. "
	experiments = np.random.randn(m,n)
	means = experiments.mean(axis = 1)
	standard_deviations = experiments.std(axis = 1)
	variances = [x*x for x in standard_deviations]
	errors = [x / math.sqrt(n) for x in standard_deviations]

	" Now calculate the stats of all the experiments. " 
	total_mean = sum(means)/m
	total_standard_deviation = math.sqrt(sum(variances))/m
	total_errors = total_standard_deviation/math.sqrt(m)
	print("Mean: " + str(total_mean))
	print("Standard Deviation: " + str(total_standard_deviation))
	print("Error: " + str(total_errors))
	plt.hist(means, 30)
	plt.show()

def one_exponential(n): # n = 100, 1000, 10000
	measurements = np.random.exponential(1,n)
	mean = measurements.mean(axis = 0)
	std = measurements.std(axis = 0)
	mean_error = std / math.sqrt(n)
	print("Mean: " + str(mean))
	print("Standard Deviation: " + str(std))
	print("Error in the Mean: " + str(mean_error))


def exponential(m,n): # m = 1000, n = 100, 1000, 10000

	" Calculate the stats for the measurements per experiment. "
	experiments = np.random.exponential(1,(m,n))
	means = experiments.mean(axis = 1)
	standard_deviations = experiments.std(axis = 1)
	variances = [x*x for x in standard_deviations]
	errors = [x / math.sqrt(n) for x in standard_deviations]

	" Now calculate the stats of all the experiments. " 
	total_mean = sum(means)/m
	total_standard_deviation = math.sqrt(sum(variances))/m
	total_errors = total_standard_deviation/math.sqrt(m)
	print("Mean: " + str(total_mean))
	print("Standard Deviation: " + str(total_standard_deviation))
	print("Error: " + str(total_errors))
	plt.hist(means, 30)
	plt.show()

def uniform(m,n): # m = 1000, n = 100, 1000, 10000

	" Calculate the stats for the measurements per experiment. "
	experiments = np.random.uniform(0,5,(m,n))
	means = experiments.mean(axis = 1)
	standard_deviations = experiments.std(axis = 1)
	variances = [x*x for x in standard_deviations]
	errors = [x / math.sqrt(n) for x in standard_deviations]

	" Now calculate the stats of all the experiments. " 
	total_mean = sum(means)/m
	total_standard_deviation = math.sqrt(sum(variances))/m
	total_errors = total_standard_deviation/math.sqrt(m)
	print("Mean :" + str(total_mean))
	print("Standard Deviation: " + str(total_standard_deviation))
	print("Error :" + str(total_errors))
	plt.hist(means, 30)
	plt.show()

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def problem5(filename):
	gamma_ray = np.loadtxt(filename)
	print(gamma_ray)
	mean = gamma_ray.mean(axis = 0)
	std = gamma_ray.std(axis = 0)
	mean_error = std / math.sqrt(gamma_ray.size)
	counts, bin_edges = np.histogram(gamma_ray,20)
	digit = np.digitize(gamma_ray, bin_edges)
	list_gamma_ray = gamma_ray.tolist()
	bins = {}
	for i in range(gamma_ray.size):
		if digit[i] in bins:
			bins[digit[i]].append(list_gamma_ray[i])
		else:
			bins[digit[i]] = [list_gamma_ray[i]]
	for key in bins:
		bins[key] = (np.asarray(bins[key]).std(axis=0))/(math.sqrt(len(bins[key])))
	plt.hist(gamma_ray, 50)
	plt.show()
