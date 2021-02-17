#!/usr/bin/python3
# -*- coding:UTF-8 -*-
# AUTHOR: Ding Bangjie
# Contact: bangjieding@outlook.com
# FILE: ~/Documents/Genetic Algoritham/test.py
# DATE: 2021/01/12 Tue
# TIME: 10:20:12

# DESCRIPTION: Simple Genetic Algorithams

import numpy as np
import random
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import heapq
import time

# 求染色体长度
def getEncodeLength(decisionvariables, delta):
    # 将每个变量的编码长度放入数组
    # 这个要看有几个自变量，假如自变量的个数是3，那么返回的lengths的元素就有3个
    lengths = []
    for decisionvar in decisionvariables:
        uper = decisionvar[1]
        low = decisionvar[0]
        # res返回一个数组
        res = fsolve(lambda x:((uper-low)/delta-2**x+1), 30)
        # ceil()向上取整
        length = int(np.ceil(res[0]))
        lengths.append(length)
    return lengths


def getinitialPopulation(length, populationSize):
    # 根据DNA长度和种群大小初始化种群
    chromsomes = np.zeros((populationSize, length), dtype=np.int)
    for popusize in range(populationSize):
        # np.random.randint产生[0,2)之间的随机整数，第三个数表示随机数的数量
        chromsomes[popusize, :] = np.random.randint(0, 2, length)
    return chromsomes

# 染色体解码得到表现形的解
def getDecode(population,encodeLength,decisionVariables, delta):
    # 把染色体计算成数字
    populationsize = population.shape[0]
    length = len(encodeLength)
    decodeVariables = np.zeros((populationsize, length), dtype=np.float)
    for i,populationchild in enumerate(population):
        start = 0
        for j, lengthchild in enumerate(encodeLength):
            power = lengthchild - 1
            decimal = 0
            for k in range(start, start+lengthchild):
                decimal += populationchild[k] * (2**power)
                power = power - 1
            start = lengthchild
            lower = decisionVariables[j][0]
            uper = decisionVariables[j][1]
            decodevalue = lower+decimal*(uper-lower)/(2**lengthchild-1)
            decodeVariables[i][j] = decodevalue
    return decodeVariables

# 得到每个个体的适应度值及累计概率
def getFitnessValue(func, decode):
    popusize, decisionvar = decode.shape
    fitnessValue = np.zeros((popusize, 1))
    for popunum in range(popusize):
        fitnessValue[popunum][0] = func(decode[popunum][0])
    probability = fitnessValue/np.sum(fitnessValue)
    cum_probability = np.cumsum(probability)
    return fitnessValue, cum_probability

def selectNewPopulation(decodepopu,cum_probability):
    m,n = decodepopu.shape
    newPopulation = np.zeros((m,n))
    for i in range(m):
        randomnum = np.random.random()
        for j in range(m):
            if(randomnum < cum_probability[j]):
                newPopulation[i] = decodepopu[j]
                break
    return newPopulation


def crossNewPopulation(newpopu, prob):
    m, n = newpopu.shape
    numbers = np.uint8(m*prob)
    if numbers %2 != 0:
        numbers = numbers + 1
    updatepopulation = np.zeros((m,n),dtype=np.uint8)
    index = random.sample(range(m),numbers)
    for i in range(m):
        if not index.__contains__(i):
            updatepopulation[i] = newpopu[i]
    j = 0
    while j < numbers:
        crosspoint = np.random.randint(0,n,1)
        crosspoint = crosspoint[0]
        updatepopulation[index[j]][0:crosspoint] = newpopu[index[j]][0:crosspoint]
        updatepopulation[index[j]][crosspoint:] = newpopu[index[j+1]][crosspoint:]
        updatepopulation[index[j+1]][0:crosspoint] = newpopu[index[j+1]][0:crosspoint]
        updatepopulation[index[j+1]][crosspoint:] = newpopu[index[j]][crosspoint:]
        j = j + 2
    return updatepopulation


def mutation(crosspopulation,mutaprob):
    mutationpopu = np.copy(crosspopulation)
    m,n = crosspopulation.shape
    mutationnums = np.uint8(m*n*mutaprob)
    mutationindex = random.sample(range(m*n),mutationnums)
    for geneindex in mutationindex:
        row = np.uint8(np.floor(geneindex/n))
        colume = geneindex % n
        if mutationpopu[row][colume] == 0:
            mutationpopu[row][colume] =1
        else:
            mutationpopu[row][colume] = 0
    return mutationpopu


def findMaxPopulation(population, maxevaluation, maxSize):
    maxevalue = maxevaluation.flatten()
    maxevaluelist = maxevalue.tolist()
    maxIndex = map(maxevaluelist.index, heapq.nlargest(100, maxevaluelist))
    index = list(maxIndex)
    colume = population.shape[1]
    maxPopulation = np.zeros((maxSize, colume))
    i = 0
    for ind in index:
        maxPopulation[i] = population[ind]
        i = i + 1
    return maxPopulation


def fitnessFunction():
    return lambda x: 10*np.sin(5*x) + 7*np.abs(x-5) + 10


def main():
    optimalvalue = []
    optimalvariables = []

    decisionVariables = [[0.0, 10.0]]
    delta = 0.0001
    EncodeLength = getEncodeLength(decisionVariables, delta)
    print(EncodeLength)

    initialPopuSize = 100
    population = getinitialPopulation(sum(EncodeLength), initialPopuSize)
    print(population)
    # 最大进化代数
    maxgeneration = 30

    # 交叉概率
    prob = 0.8

    # 变异概率
    mutationprob = 0.1

    # 新生成的种群数量
    maxPopuSize = 100

    for generation in range(maxgeneration):
        decode = getDecode(population, EncodeLength, decisionVariables, delta)

        evaluation, cum_proba = getFitnessValue(fitnessFunction(), decode)

        newpopulations = selectNewPopulation(population, cum_proba)

        crossPopulation = crossNewPopulation(newpopulations, prob)

        mutationpopulation = mutation(crossPopulation, mutationprob)

        totalpopalation = np.vstack((population, mutationpopulation))

        final_decode = getDecode(totalpopalation, EncodeLength, decisionVariables,delta)

        final_evaluation, final_cumprob = getFitnessValue(fitnessFunction(), final_decode)

        population = findMaxPopulation(totalpopalation, final_evaluation, maxPopuSize)

        optimalvalue.append(np.max(final_evaluation))

        index = np.where(final_evaluation == max(final_evaluation))

        optimalvariables.append(list(final_decode[index[0][0]]))
        decodeTemp = getDecode(population, EncodeLength, decisionVariables, delta)
        evaluationTemp, cum_probaTemp = getFitnessValue(fitnessFunction(), decodeTemp)
        plt.plot(decodeTemp,evaluationTemp,"*")
        x = np.arange(0.0, 10, 0.001)
        y = 10 * np.sin(5 * x) + 7 * np.abs(x - 5) + 10
        plt.plot(x, y)
        plt.show()
        time.sleep(1)


    x = [i for i in range(maxgeneration)]
    y = [optimalvalue[i] for i in range(maxgeneration)]
    plt.plot(x,y)
    plt.show()
    optimalval = np.max(optimalvalue)
    index = np.where(optimalvalue == max(optimalvalue))

    optimavar = optimalvariables[index[0][0]]
    return optimalval, optimavar


if __name__ == '__main__':
    optval, optvar = main()
    print("x1:",optvar[0])
    print("maxValue:",optval)