import itertools

# Write Includes - Write some standard C++ headers to output file.

def WriteIncludes(file) :

    file.write("#include <vector>")
    file.write("\n#include <map>")
    file.write("\n#include <stdio>")
    file.write("\nusing namespace std;\n")

# GenerateBasis - Given an algebra signature, generate all basis vectors
# for this algebra.

def GenerateBasis(signature, basisVecs) :

    oneVecs         = []
    basisCounter    = 0

    for item in signature :

        oneVecs.append("e" + str(basisCounter))
        basisVecs.append(["e" + str(basisCounter)])
        basisCounter = basisCounter + 1
    
    for l in range(2, len(signature) + 1) :

        for cmb in itertools.combinations(oneVecs, l) :
      
            basisVecs.append(list(cmb))

# Get parity -

def GetParity(vec) :

    swaps = 0

    if vec == None :

        return 1

    for i in range(len(vec), 1, -1) :

        maxVec = max(vec[0 : i])
        idx = vec.index(maxVec)
        swaps = swaps + (i - idx - 1)
        vec.remove(maxVec)
        vec.insert(i - 1, maxVec)

    return pow(-1, swaps)


# BasisVecInnerProduct - Compute the inner product of two basis vectors
# in an algebra.

def BasisVecInnerProduct(vec1, vec2, signature) :

    intsc = set(vec1).intersection(vec2) 
    if len(intsc) == 0 :

        return (0, [])

    tmp1 = list(vec1)
    tmp2 = list(vec2)
    parity = 1

    for item in intsc :

        swaps1 = len(tmp1) - tmp1.index(item) - 1
        swaps2 = tmp2.index(item)
        sig = signature[int(item[1])]

        parity = parity * pow(-1, swaps1 + swaps2) * sig
        
        tmp1.remove(item)
        tmp2.remove(item)

    for oneVec in tmp2 : 

        if oneVec != '' :

            tmp1.append(oneVec)

    parity = parity * GetParity(tmp1)

    if tmp1 == None : 

        tmp1 = ['']

    return (parity, tmp1)

# BasisVecOuterProduct - Compute the outer product of two basis vectors
# in an algebra.

def BasisVecOuterProduct(vec1, vec2, signature) :

    if len(set(vec1).intersection(vec2)) > 0 :

        return (0, [])

    tmp1 = list(vec1)
    tmp2 = list(vec2)

    for oneVec in tmp2 : 

        if oneVec != '' :

            tmp1.append(oneVec)

    parity = GetParity(tmp1)

    return (parity, tmp1)

# GenerateProductTables - Given an algebra and its basis vectors, generate its
# inner and outer multiplication tables.

def GenerateProductTables(basisVecs, signature, innerProductTable, outerProductTable) :

    r = 0

    for vec1 in basisVecs :

        c = 0

        for vec2 in basisVecs :

            innerProd = BasisVecInnerProduct(vec1, vec2, signature)
            outerProd = BasisVecOuterProduct(vec1, vec2, signature)

            print(vec1)
            print(vec2)
            print("Inner:")
            print(innerProd)
            print("Outer:")
            print(outerProd)

            print("\n")

            innerProductTable[r][c] = innerProd
            outerProductTable[r][c] = outerProd
            c = c + 1

        r = r + 1

# WriteBasisEnum - Write a C++ enum from basis vecs to file.

def WriteBasisEnum(file, algebraName, basisVecs) :

    file.write("\nenum " + algebraName + "Basis")
    file.write("\n{")
   
    for vec in basisVecs :

        file.write("\n\t")

        for oneVec in vec :

            file.write(oneVec)

        file.write(",")

    file.write("\n};")

# GenerateAlgebra - Given an algebra signature and input name, generate its geometric
# algebra and generate C++ code for this algebra.

def GenerateAlgebra(algebraName, signature) :

    print("Generating Algebra")
    outFile = open(algebraName + ".h", 'w')
    
    basisVecs = []
    dim = pow(2, len(signature))

    innerProductTable = [[0 for x in range(dim)] for y in range(dim)]
    outerProductTable = [[0 for x in range(dim)] for y in range(dim)]

    GenerateBasis(signature, basisVecs)
    GenerateProductTables(basisVecs, signature, innerProductTable, outerProductTable)

    WriteIncludes(outFile)
    WriteBasisEnum(outFile, algebraName, basisVecs)

# Main

signature = [ 1, -1, -1, -1 ]

GenerateAlgebra("STA", signature)
print("Done")