import itertools

# Write Includes - Write some standard C++ headers to output file.

def WriteIncludes(file) :

    file.write("#include <vector>\n")
    file.write("#include <map>\n")
    file.write("#include <stdio>\n")
    file.write("#include <stdint.h>\n\n")
    file.write("using namespace std;\n\n")

# GenerateBasis - Given an algebra signature, generate all basis vectors
# for this algebra.

def GenerateBasis(signature, basisVecs) :

    oneVecs         = []
    basisCounter    = 0
    basisVecs.append("s")

    for item in signature :

        oneVecs.append("e" + str(basisCounter))
        basisVecs.append(["e" + str(basisCounter)])
        basisCounter = basisCounter + 1
    
    for l in range(2, len(signature) + 1) :

        for cmb in itertools.combinations(oneVecs, l) :
      
            basisVecs.append(list(cmb))

# Get parity - Determine sign of permutation required to transform a basis
# vec into a basis vec with indices in ascending order. E.g., e2e1e3 -> e1e2e3
# has parity -1.

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

    if vec1 == 's' :

        return (1, list(vec2))

    if vec2 == 's' :

        return (1, list(vec1))

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
# in a geometric algebra.

def BasisVecOuterProduct(vec1, vec2, signature) :

    if vec1 == 's' :

        return (0, '')

    if vec2 == 's' :

        return (0, '')

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

    file.write("enum " + algebraName + "Basis\n")
    file.write("{\n")
   
    for vec in basisVecs :

        file.write("\t")

        for oneVec in vec :

            file.write(oneVec)

        file.write(",\n")

    file.write("};\n\n")

# WriteMultiVecStruct - Write a structure definition for a general multivector
# in a given geometric algebra.

def WriteMultivecStruct(file, algebraName, basisVecs) :

    file.write("struct " + algebraName + "MV\n")
    file.write("{\n")
    file.write("\tconst uint32_t numCoeffs = " + str(len(basisVecs)) + ";\n")
    file.write("\tfloat coeffs[numCoeffs];\n\n")

    WriteFullOperators(file, algebraName, basisVecs)

    file.write("};\n\n")

# WriteMultiVecStruct - Write a structure definition for an even multivector
# in a given geometric algebra.

def WriteRotorStruct(file, algebraName, basisVecs) :

    file.write("struct " + algebraName + "Rotor\n")
    file.write("{\n")
    file.write("\tconst uint32_t numCoeffs = " + str(len(basisVecs) / 2) + ";\n")
    file.write("\tfloat coeffs[numCoeffs];\n\n")
    file.write("};\n\n")

# WriteFullOperators - Write arithmetic operators for general multivectors.

def WriteFullOperators(file, algebraName, basisVecs) :

    WriteFullEqualityOperators(file, algebraName, basisVecs)
    WriteFullAddOperators(file, algebraName, basisVecs)
    WriteFullSubtractOperators(file, algebraName, basisVecs)
    WriteFullMultiplyOperators(file, algebraName, basisVecs)
    WriteFullDivideOperators(file, algebraName, basisVecs)

# WriteFullEqualityOperators -

def WriteFullEqualityOperators(file, algebraName, basisVecs) :

    structName = algebraName + "MV"

    file.write("\t" + structName + "& operator=(const " + structName + "& rhs)\n")
    file.write("\t{\n")
    file.write("\t\tif (this != &rhs)\n")
    file.write("\t\t{\n")

    file.write("\t\t\tfor (uint32_t i = 0; i < numCoeffs; i++)\n")
    file.write("\t\t\t{\n")
    file.write("\t\t\t\tcoeffs[i] = rhs.coeffs[i];\n")
    file.write("\t\t\t}\n")

    file.write("\t\t}\n")
    file.write("\t\treturn *this;\n")
    file.write("\t}\n\n")

# WriteFullAddOperators - Write addition operators for general multivectors.

def WriteFullAddOperators(file, algebraName, basisVecs) :

    structName = algebraName + "MV"
    file.write("\t" + structName + "& operator+=(const " + structName + "& rhs)\n")
    file.write("\t{\n")
    file.write("\t\tfor (uint32_t i = 0; i < numCoeffs; i++)\n")
    file.write("\t\t{\n")
    file.write("\t\t\tcoeffs[i] += rhs.coeffs[i];\n")
    file.write("\t\t}\n")
    file.write("\t\treturn *this;\n")
    file.write("\t}\n\n")
    
    return

# WriteFullSubtractOperators - Write subtraction operators for general multivectors.

def WriteFullSubtractOperators(file, algebraName, basisVecs) :

    structName = algebraName + "MV"
    file.write("\t" + structName + "& operator-=(const " + structName + "& rhs)\n")
    file.write("\t{\n")
    file.write("\t\tfor (uint32_t i = 0; i < numCoeffs; i++)\n")
    file.write("\t\t{\n")
    file.write("\t\t\tcoeffs[i] -= rhs.coeffs[i];\n")
    file.write("\t\t}\n")
    file.write("\t\treturn *this;\n")
    file.write("\t}\n\n")

    return

# WriteFullMultiplyOperators - Write multiplication operators for general multivectors.

def WriteFullMultiplyOperators(file, algebraName, basisVecs) :

    structName = algebraName + "MV"
    file.write("\t" + structName + "& operator|=(const " + structName + "& rhs)\n")
    file.write("\t{\n")
    file.write("\t\tfor (uint32_t i = 0; i < numCoeffs; i++)\n")
    file.write("\t\t{\n")
    file.write("\t\t\tcoeffs[i] -= rhs.coeffs[i];\n")
    file.write("\t\t}\n")
    file.write("\t\treturn *this;\n")
    file.write("\t}\n\n")
    return

# WriteFullDivideOperators - Write division operators for general multivectors.

def WriteFullDivideOperators(file, algebraName, basisVecs) :

    return

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
    WriteMultivecStruct(outFile, algebraName, basisVecs)
    WriteRotorStruct(outFile, algebraName, basisVecs)

# Main

signature = [ 1, -1, -1, -1 ]
GenerateAlgebra("STA", signature)
print("Done")