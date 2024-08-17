import numpy as np
from sympy.ntheory import factorint, primerange
import math
import time
import ecc

def factorBase(n):
    limit = int(3.38 * np.exp(0.5 * np.sqrt(np.log2(n) * np.log2(np.log2(n)))))
    print("limit: ", limit)
    
    factorBase = list(primerange(2, limit + 1))
    
    return factorBase

def isSystemSolvable(a, b):
    A = np.array(a)
    B = np.array(b).reshape(-1, 1)
    rankA = np.linalg.matrix_rank(A)
    rankAB = np.linalg.matrix_rank(np.hstack((A, B)))
    
    return rankA == rankAB == len(A)

def isSquare(m):
    return all(len(row) == len(m) for row in m)

def solveSystem2(A, B, n):
    sle = np.hstack((A, B.reshape(-1, 1))) 
    numRows, numCols = sle.shape

    processed = []
    for j in range(numCols - 1):
        gcdColMod = [math.gcd(elem, (n - 1)) for elem in sle[:, j]]
        
        for i, gcdVal in enumerate(gcdColMod):
            if i in processed:
                continue
        
            if gcdVal == 1:
                invElem = pow(int(sle[i, j]), -1, (n - 1))
                processed.append(i)

                sle[i] = (sle[i] * invElem) % (n - 1)
                
                mask = np.arange(numRows) != i
                sle[mask] = (sle[mask] - np.outer(sle[mask, j], sle[i])) % (n - 1)
                break
        
    solution = []
    for j in range(numCols - 1):
        nonZeroIndices = np.nonzero(sle[:, j])[0]
        if nonZeroIndices.size > 0:
            solution.append(sle[nonZeroIndices[0], -1])
        else:
            solution.append(0)
    
    return solution

def linearSystem(factorBase, a, n):
    system = [] # Ax = B
    rightPart = []
    
    i = 0
    while True:
        subLogarithm = pow(a, i, n)
        
        logDecomposition = factorint(subLogarithm)
        isSmoothLog = all(primeNumber in factorBase for primeNumber in logDecomposition.keys())
        
        if not isSmoothLog:
            i += 1
            continue
        
        equation = []
        for primeNumber in factorBase:
            if primeNumber in logDecomposition:
                equation.append(logDecomposition[primeNumber])
            else:
                equation.append(0)
                
        system.append(equation)
        rightPart.append(i)
        
        A = np.array(system)
        B = np.array(rightPart)
        if len(system) >= len(factorBase) + 20:
            print("system: ", system)
            print("right part: ", rightPart)
            solution = solveSystem2(A, B, n)
            
            # if solution is not None:
            return solution
                
        i += 1
        
def evaluateLogarithm(factorBase, systemSolution, a, b, n):
    i = 0
    while True:
        number = (b * pow(a, i, n)) % n
        
        numberDecomposition = factorint(number)
        isSmoothNumber = all(primeNumber in factorBase for primeNumber in numberDecomposition.keys())
        
        if not isSmoothNumber:
            i += 1
            continue
        
        result = 0
        for j in range(len(factorBase)):
            if factorBase[j] in numberDecomposition:
                result = (result + (numberDecomposition[factorBase[j]] * systemSolution[j])) % (n-1)
        
        return ((result - i) % (n-1))

def indexCalculus(a, b, n):
    base = factorBase(n)
    solution = linearSystem(base, a, n)
    print("solution: ", solution)
    result = evaluateLogarithm(base, solution, a, b, n)
        
    return result

if __name__ == "__main__":
    a = 5
    b = 10
    n = 17
    
    startTime = time.time()

    print("result of index calculus:", indexCalculus(a, b, n))

    endTime = time.time()
    executionTime = endTime - startTime
    print("Execution time:", executionTime, "seconds")
