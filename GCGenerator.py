import itertools

# Data structure for holding information about a geometric algebra, e.g., its signature,
# basis vectors, and signs for multiplication between basis vectors.

class Algebra:

    def __init__(self, name, sig, even):

        self.name           = name;
        self.sig            = sig;
        self.basis          = []
        self.even           = even

        GenerateBasis(self.sig, self.basis)

        self.dim                = pow(2, len(self.sig))
        self.innerProductTable  = [[0 for x in range(self.dim)] for y in range(self.dim)]
        self.outerProductTable  = [[0 for x in range(self.dim)] for y in range(self.dim)]
        self.innerProdImg       = []
        self.outerProdImg       = []

        GenerateProductTables(self.basis, self.sig, self.innerProductTable, self.outerProductTable)
        FlipProdTables(self.innerProductTable, self.outerProductTable, self.innerProdImg, self.outerProdImg)

# Write Includes - Write some standard C++ headers to output file.
#
# @param [in] File to write includes to.

def WriteIncludes(file, algebra) :

    file.write("#include <vector>\n")
    file.write("#include <map>\n")
    file.write("#include <stdio>\n")
    file.write("#include <stdint.h>\n\n")
    file.write("using namespace std;\n\n")

# GenerateBasis - Given an algebra signature, generate all basis vectors
# for this algebra, e.g., [s e1 e2 e1e2].
#
# @param signature [in] A list of signs of inner products of basis one-vectors, e.g., [e1^2 e2^2] = [1 1].
# @param basis     [out] Full basis generated from signature.

def GenerateBasis(signature, basisVecs) :

    oneVecs         = []
    basisCounter    = 0
    basisVecs.append("S")

    for item in signature :

        oneVecs.append("E" + str(basisCounter))
        basisVecs.append(["E" + str(basisCounter)])
        basisCounter = basisCounter + 1
    
    for l in range(2, len(signature) + 1) :

        for cmb in itertools.combinations(oneVecs, l) :
      
            basisVecs.append(list(cmb))

# Get parity - Determine sign of permutation required to transform a basis
# vec into a basis vec with indices in ascending order. E.g., e2e1e3 -> e1e2e3
# requires 1 swap and has parity -1. Used when determining the sign of basis vector
# products. 
#
# @param vec The parity of a vector, e.g., e1e2 = 1 and e2e1 = -1.

def GetParity(vec) :

    swaps = 0

    if vec == None :

        return 1

    for i in range(len(vec), 1, -1) :

        maxVec  = max(vec[0 : i])
        idx     = vec.index(maxVec)
        swaps   = swaps + (i - idx - 1)
        
        vec.remove(maxVec)
        vec.insert(i - 1, maxVec)

    return pow(-1, swaps)


# BasisVecInnerProduct - Compute the inner product of two basis vectors
# in an algebra.

def BasisVecInnerProduct(vec1, vec2, signature) :

    # Convention here is the inner product of vector with scalar returns scaled vector.

    if vec1 == 'S' :

        return (1, list(vec2))

    if vec2 == 'S' :

        return (1, list(vec1))

    # Dot product removes common basis one-vecs of operands. Count swaps to move each
    # common basic one-vecs to end of first operand and beginning of second operand
    # and compute running parity.

    intsc = set(vec1).intersection(vec2) 
    
    if len(intsc) == 0 :

        return (0, [])

    tmp1    = list(vec1)
    tmp2    = list(vec2)
    parity  = 1

    for item in intsc :

        swaps1  = len(tmp1) - tmp1.index(item) - 1
        swaps2  = tmp2.index(item)
        sig     = signature[int(item[1])]
        parity  = parity * pow(-1, swaps1 + swaps2) * sig
        
        tmp1.remove(item)
        tmp2.remove(item)

    # Take union of basis one-vecs after common one-vecs have been dotted out.

    for oneVec in tmp2 : 

        if oneVec != '' :

            tmp1.append(oneVec)

    # Rearrange resulting vector into ascending form.

    parity = parity * GetParity(tmp1)

    if tmp1 == [] : 

        tmp1 = ['S']

    return (parity, tmp1)

# BasisVecOuterProduct - Compute the outer product of two basis vectors
# in a geometric algebra.
#
# @param vec1 [in]
# @param vec2 [in]
# @param signature [in]

def BasisVecOuterProduct(vec1, vec2, signature) :

    # Convention here is that outer product with scalar is null.

    if vec1 == 'S' :

        return (0, '')

    if vec2 == 'S' :

        return (0, '')

    # Outher product is empty if operands have common basis one-vecs.

    if len(set(vec1).intersection(vec2)) > 0 :

        return (0, [])

    # Take union of input vecs and rearrange to ascending form. 

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

            innerProductTable[r][c] = BasisVecInnerProduct(vec1, vec2, signature)
            outerProductTable[r][c] = BasisVecOuterProduct(vec1, vec2, signature)
            c = c + 1

        r = r + 1

# FlipProdTables - When doing products, need all pairs of input components that will map to an output component and appropriate sign.
# E.g., for complex numbers out.real = in1.real * in2.real - in1.imag * in2.imag, or out[0] = in1[0] * in2[0] - in1[1] * in2[1].
# Build this table of mappings: {0} : {[0, 0], 1}, {[1, 1], -1}
#
# @param innerProductTable
# @param outerProductTable
# @param innerProdImg
# @param outerProdImg

def FlipProdTables(innerProductTable, outerProductTable, innerProdImg, outerProdImg):


    return

# WriteMultivecStruct - 
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteMultivecStruct(file, algebra):

    WriteOpenMultivecStruct(file, algebra)
    WriteConstructors(file, algebra)
    WriteAssignmentOperators(file, algebra)
    WriteCloseMultivecStruct(file, algebra)


# WriteOpenMultiVecStruct - Write a the opening of a full multivector
# struct in C++.
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteOpenMultivecStruct(file, algebra):

    file.write("static const uint32_t numCoeffs = " + str(algebra.dim) + ";\n\n");
    
    file.write("template<typename T> struct " + algebra.name + "MV\n")
    file.write("{\n")
    file.write("    T coeffs[" + str(len(algebra.basis)) + "];\n\n")

# WriteOpenMultiVecStruct -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteCloseMultivecStruct(file, algebra):

    file.write("};")

# WriteConstructors - Write constructors for a given algebra.
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteConstructors(file, algebra):

    # Default constructor. All basis vector coefficients set to 0.

    file.write("    " + algebra.name + "<T>() { memset(&coeffs[0], 0, numCoeffs * sizeof(T)); }\n\n")
    
    # Constructor with all coefficients specified.

    file.write("    " + algebra.name + "<T>(");

    for i, basisVec in enumerate(algebra.basis):

        vecStr = ""

        for oneVec in basisVec :

            vecStr = vecStr + oneVec

        file.write("T " + vecStr)
        
        if i < len(algebra.basis) - 1:

            file.write(", ")

    file.write(")\n")
    file.write("    {\n")
    
    for i, basisVec in enumerate(algebra.basis):

        vecStr = ""

        for oneVec in basisVec :

            vecStr = vecStr + oneVec

        file.write("        coeffs[" + str(i) + "] = " + vecStr + ";\n")

    file.write("    }\n\n")

# WriteFullOperators - Write arithmetic operators for general multivectors.
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteAssignmentOperators(file, algebra):

    WriteAssignOperator(file, algebra)
    WriteAddAssignOperator(file, algebra)
    WriteSubAssignOperator(file, algebra)

    #WriteDotAssignOperator(file, algebra)
    #WriteWedgeAssignOperator(file, algebra)
    #WriteMultAssignOperator(file, algebra)
    #WriteDivAssignOperator(file, algebra)

# WriteAssignOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteAssignOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator=(const " + structName + "<T>& rhs)\n")
    file.write("    {\n")
    file.write("        if (this != &rhs)\n")
    file.write("            memcpy(&coeffs[0], &rhs.coeffs[0], numCoeffs * sizeof(T));\n");
    file.write("\n");
    file.write("        return *this;\n")
    file.write("    }\n\n")

# WriteAddAssignOperator - Write addition operators for general multivectors.
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteAddAssignOperator(file, algebra):

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator+=(const " + structName + "<T>& rhs)\n")
    file.write("    {\n")
    file.write("        for (uint32_t i = 0; i < numCoeffs; i++)\n")
    file.write("            coeffs[i] += rhs.coeffs[i];\n")
    file.write("\n");
    file.write("        return *this;\n")
    file.write("    }\n\n")
    
    return

# WriteSubAssignOperator - Write subtraction operators for general multivectors.
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteSubAssignOperator(file, algebra):

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator-=(const " + structName + "<T>& rhs)\n")
    file.write("    {\n")
    file.write("        for (uint32_t i = 0; i < numCoeffs; i++)\n")
    file.write("            coeffs[i] -= rhs.coeffs[i];\n")
    file.write("\n");
    file.write("        return *this;\n")
    file.write("    }\n\n")
    
    return

    return

# WriteDotAssignOperator - Write division operators for general multivectors.
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteDotAssignOperators(file, algebra):

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator-=(const " + structName + "<T>& rhs)\n")
    file.write("    {\n")
    file.write("        for (uint32_t i = 0; i < numCoeffs; i++)\n")
    file.write("            coeffs[i] -= rhs.coeffs[i];\n")
    file.write("\n");
    file.write("        return *this;\n")
    file.write("    }\n\n")
    
    return

    return

# WriteWedgeAssignOperator - Write division operators for general multivectors.
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteWedgeAssignOperators(file, algebra):

    structName = algebra.name + "MV<T>"

    prodTbl = algebra.innerProductTable;

    file.write("    " + structName + "& operator|=(const " + structName + "<T>& rhs)\n")
    file.write("    {\n")
    file.write("        for (uint32_t i = 0; i < numCoeffs; i++)\n")
    file.write("        {\n")
    file.write("            coeffs[i] -= rhs.coeffs[i];\n")
    file.write("        }\n")
    file.write("        return *this;\n")
    file.write("    }\n\n")

    return

# WriteMultAssignOperator - Write multiplication operators for general multivectors.
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteMultAssignOperators(file, algebra):

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator*=(const " + structName + "<T>& rhs)\n")
    file.write("    {\n")
    file.write("        for (uint32_t i = 0; i < numCoeffs; i++)\n")
    file.write("        {\n")
    file.write("            coeffs[i] -= rhs.coeffs[i];\n")
    file.write("        }\n")
    file.write("        return *this;\n")
    file.write("    }\n\n")

    return

# WriteDivAssignOperators - Write division operators for general multivectors.
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteDivAssignOperators(file, algebra):

    return

# GenerateAlgebra - Given an algebra signature and input name, generate its geometric
# algebra and generate C++ code for this algebra.

def GenerateAlgebra(algebraName, signature, even) :

    file    = open(algebraName + ".h", 'w')
    algebra = Algebra("E2", signature, even);

    WriteMultivecStruct(file, algebra);

# Main

GenerateAlgebra("E2", [1, 1], False)