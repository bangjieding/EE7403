#!/usr/bin/python3
# -*- coding:UTF-8 -*-
# AUTHOR: Ding Bangjie
# Contact: bangjieding@outlook.com
# FILE: ~/Documents/Genetic Algoritham/SGA.py
# DATE: 2021/01/12 Tue
# TIME: 10:20:42

# DESCRIPTION: Simple Genetic Algoritham

# Maximizing f(x) = x + 10sin(5x) + 7cos(4x)
# Assume x is an integer between 0 and 9

## Assuming that the accuracy of the solution requires **four decimal** places
### At the first assuming that delta is 10^-1

import numpy as np
import math

def get_choromosome_length(upper_bound, lower_bound, delta):
    return len(bin((upper_bound - lower_bound)/delta)[2:])

def initial_population(popu_size, choro_len):
    chromsomes = np.zeros((popu_size, choro_len), dtype=np.int)
    for i in popu_size:
        chromsomes[i, :] = np.random.randint(0, 2, choro_len)
    return chromsomes
    
def fitness_function():
    pass

def get_genetype():
    pass

def get_phenotype(popu_size, choromosomes, choro_len, upper_bound, lower_bound):
    phenotypes = np.zeros((popu_size), dtype=float)
    for i in popu_size:
        # temp_list = choromosomes[i].tolist()
        phenotypes[i] = int(''.join(str(x) for x in choromosomes[i].tolist()), 2) / (2**choro_len - 1)
