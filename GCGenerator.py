import itertools

# Data structure for holding information about a geometric algebra, e.g., its signature,
# basis vectors, and signs for multiplication between basis vectors.

class Algebra:

    def __init__(self, name, sig, even):

        if even == True : name = name + "_E"

        self.name           = name;
        self.sig            = sig;
        self.basis          = []
        self.even           = even

        GenerateBasis(self.sig, self.basis, even)

        if even == True : self.dim = pow(2, len(self.sig) - 1)
        else : self.dim = pow(2, len(self.sig))

        self.innerProductTable  = [[0 for x in range(self.dim)] for y in range(self.dim)]
        self.outerProductTable  = [[0 for x in range(self.dim)] for y in range(self.dim)]

        GenerateProductTables(self.basis, self.sig, self.innerProductTable, self.outerProductTable)

        self.innerProdImg       = [[] for x in range(self.dim)]
        self.outerProdImg       = [[] for x in range(self.dim)]

        FlipProdTables(self)

# WritePreamble - Write some #pragmas and #includes to top of file.
#
# @param [in] File to write includes to.

def WritePreamble(file, algebra) :

    file.write("#pragma once\n")
    file.write("#include <math.h>\n")
    file.write("#include <stdint.h>\n")
    file.write("#include \"gcclib.h\"\n\n")

# GenerateBasis - Given an algebra signature, generate all basis vectors
# for this algebra, e.g., given [e1 e2], generate [s e1 e2 e1e2].
#
# @param signature [in] A list of signs of inner products of basis one-vectors, e.g., [e1^2 e2^2] = [1 1].
# @param basis     [out] Full basis generated from signature.
# @param even      [in] Should this algebra only contain even basis vectors.

def GenerateBasis(signature, basisVecs, even) :

    oneVecs         = []
    basisCounter    = 0
    basisVecs.append(["S"])

    for item in signature :

        oneVecs.append("E" + str(basisCounter))
        if even == False : basisVecs.append(["E" + str(basisCounter)])
        basisCounter = basisCounter + 1
    
    for l in range(2, len(signature) + 1) :

        for cmb in itertools.combinations(oneVecs, l) : 

            if even : 
            
                if len(cmb) % 2 == 0 or cmb[0] == 'S' : basisVecs.append(list(cmb))

            else :  basisVecs.append(list(cmb))


# Get parity - Determine sign of permutation required to transform a basis
# vec into a basis vec with indices in ascending order. E.g., e2e1e3 -> e1e2e3
# requires 1 swap and has parity -1. Used when determining the signs of basis vector
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

    if vec1 == ['S'] : return (1, list(vec2))
    if vec2 == ['S'] : return (1, list(vec1))

    # Dot product removes common basis one-vecs of operands. First, find all common basis
    # one-vecs. For each common one-vec, count swaps to move it to end of first operand and beginning of 
    # second operand. For every swap, toggle the sign of the final product.

    intsc = set(vec1).intersection(vec2) 
    
    if len(intsc) == 0 : return (0, [])

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

        if oneVec != '' : tmp1.append(oneVec)

    # Reorder resulting basis vector into ascending form and multiply by its parity.

    parity = parity * GetParity(tmp1)

    if tmp1 == [] : tmp1 = ['S']

    return (parity, tmp1)

# BasisVecOuterProduct - Compute the outer product of two basis vectors.
#
# @param vec1 [in]
# @param vec2 [in]
# @param signature [in]

def BasisVecOuterProduct(vec1, vec2, signature) :

    # Convention here is that outer product of vector with scalar is null.

    if vec1 == ['S'] : return (0, '')
    if vec2 == ['S'] : return (0, '')

    # Outer product is empty if operands have any common basis one-vecs.

    if len(set(vec1).intersection(vec2)) > 0 : return (0, [])

    # Take union of input vecs and reorder to ascending form. 

    tmp1 = list(vec1)
    tmp2 = list(vec2)

    for oneVec in tmp2 : 

        if oneVec != '' : tmp1.append(oneVec)

    parity = GetParity(tmp1)

    return (parity, tmp1)

# GenerateProductTables - Given an algebra and its basis vectors, generate its
# inner and outer multiplication tables e.g., dot(e2, e1) = {e1e2, - 1), etc.

def GenerateProductTables(basisVecs, signature, innerProductTable, outerProductTable) :

    r = 0

    for vec1 in basisVecs :

        c = 0

        for vec2 in basisVecs :

            innerProductTable[r][c] = BasisVecInnerProduct(vec1, vec2, signature)
            outerProductTable[r][c] = BasisVecOuterProduct(vec1, vec2, signature)
            c = c + 1

        r = r + 1

# FlipProdTables - When doing products, we need all pairs of input components that will map to an output component along with sign.
# E.g., for complex numbers out.real = in1.real * in2.real - in1.imag * in2.imag, or out[0] = in1[0] * in2[0] - in1[1] * in2[1].
# Build this table of mappings: {0} : {[0, 0], 1}, {[1, 1], -1}
#
# @param algebra            [in] Algebra with inner/outer product tables defined and flipped tables to be populated.

def FlipProdTables(algebra):

    innerProductTable   = algebra.innerProductTable
    outerProductTable   = algebra.outerProductTable
    innerProdImg        = algebra.innerProdImg
    outerProdImg        = algebra.outerProdImg

    for i, row in enumerate(innerProductTable):

        for j, col in enumerate(row):

            basisVec = col[1]
            sign     = col[0]

            if sign != 0 :

                outIdx = algebra.basis.index(basisVec)
                innerProdImg[outIdx].append([i, j, sign])

    for i, row in enumerate(outerProductTable):

        for j, col in enumerate(row):

            basisVec = col[1]
            sign     = col[0]

            if sign != 0 :

                outIdx = algebra.basis.index(basisVec)
                outerProdImg[outIdx].append([i, j, sign])

    return

# WriteMultivecStruct - 
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteMultivecStruct(file, algebra):

    WritePreamble(file, algebra)
    WriteOpenMultivecStruct(file, algebra)
    WriteConstructors(file, algebra)
    WriteAssignmentOperators(file, algebra)
    WriteBinaryOperators(file, algebra)
    WriteMiscFunctions(file, algebra)
    WriteCloseMultivecStruct(file, algebra)


# WriteOpenMultiVecStruct - Write a the opening of a full multivector
# struct in C++.
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteOpenMultivecStruct(file, algebra):

    file.write("static const uint32_t numCoeffs = " + str(algebra.dim) + ";\n\n");
    
    file.write("template<typename T> struct GCCLIB_API " + algebra.name + "MV\n")
    file.write("{\n")
    file.write("    T coeffs[numCoeffs];\n\n")

# WriteCloseMultiVecStruct -
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

    structName = algebra.name + "MV<T>"

    # Default constructor. All basis vector coefficients set to 0.

    file.write("    " + structName + "() { memset(&coeffs[0], 0, numCoeffs * sizeof(T)); }\n\n")
    
    # Constructor with all coefficients specified.

    file.write("    " + structName + "(");

    for i, basisVec in enumerate(algebra.basis):

        vecStr = ""

        for oneVec in basisVec : vecStr = vecStr + oneVec
        file.write("T " + vecStr)
        if i < len(algebra.basis) - 1: file.write(", ")

    file.write(")\n")
    file.write("    {\n")
    
    for i, basisVec in enumerate(algebra.basis):

        vecStr = ""

        for oneVec in basisVec : vecStr = vecStr + oneVec

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
    WriteDotAssignOperator(file, algebra)
    WriteWedgeAssignOperator(file, algebra)
    WriteMultAssignOperator(file, algebra)
    WriteScalarMultAssignOperator(file, algebra)
    WriteDivAssignOperator(file, algebra)
    WriteScalarDivAssignOperator(file, algebra)

# WriteBinaryOperators - 
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteBinaryOperators(file, algebra) :

    WriteAddOperator(file, algebra)
    WriteSubOperator(file, algebra)
    WriteDotOperator(file, algebra)
    WriteWedgeOperator(file, algebra)
    WriteMultOperator(file, algebra)
    WriteScalarMultOperator(file, algebra)
    WriteDivOperator(file, algebra)
    WriteScalarDivOperator(file, algebra)

# WriteMiscFunctions - 
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteMiscFunctions(file, algebra) :

    WriteReverseOperator(file, algebra)
    WriteAccessors(file, algebra)
    WriteMagnitudeOperator(file, algebra)

# WriteScalarMultAssignOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteScalarMultAssignOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator*=(const T rhs)\n")
    file.write("    {\n")
    file.write("        for (uint32_t i = 0; i < numCoeffs; i++)\n")
    file.write("            coeffs[i] *= rhs;\n")
    file.write("\n");
    file.write("        return *this;\n")
    file.write("    }\n\n")

    return

# WriteScalarDivAssignOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteScalarDivAssignOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator/=(const T rhs)\n")
    file.write("    {\n")
    file.write("        for (uint32_t i = 0; i < numCoeffs; i++)\n")
    file.write("            coeffs[i] /= rhs;\n")
    file.write("\n");
    file.write("        return *this;\n")
    file.write("    }\n\n")

    return

# WriteMagnitudeOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteMagnitudeOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    T Mag()\n")
    file.write("    {\n")
    file.write("        T mag = 0.0;\n")
    file.write("        for (uint32_t i = 0; i < nCoeffs; i++) mag += coeffs[i] * coeffs[i];\n")
    file.write("        return sqrt(mag);\n")
    file.write("    }\n\n")

    return

# WriteReverseOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteReverseOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + " operator~()\n")
    file.write("    {\n")
    file.write("        " + structName + " tmp = *this;\n")

    for i, basisVec in enumerate(algebra.basis) :

        if ((len(basisVec) * (len(basisVec) - 1)) / 2) % 2 == 1 : 

            file.write("        tmp.coeffs[" + str(i) + "] = -tmp.coeffs[" + str(i) + "];\n")

    file.write("\n")
    file.write("        return tmp;\n")
    file.write("    }\n\n")

    return

# WriteAccessors -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteAccessors(file, algebra) :

    for i, basisVec in enumerate(algebra.basis) :

        accessorName = ""
        for oneVec in basisVec : accessorName += oneVec

        file.write("    T " + accessorName + "(){ return coeffs[" + str(i) + "]; }\n")
        file.write("    void " + accessorName + "(T val){ coeffs[" + str(i) + "] = val; }\n\n")

    return

# WriteAssignOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteAssignOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator=(const " + structName + "& rhs)\n")
    file.write("    {\n")
    file.write("        if (this != &rhs)\n")
    file.write("            memcpy(&coeffs[0], &rhs.coeffs[0], numCoeffs * sizeof(T));\n");
    file.write("\n");
    file.write("        return *this;\n")
    file.write("    }\n\n")

# WriteAddAssignOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteAddAssignOperator(file, algebra):

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator+=(const " + structName + "& rhs)\n")
    file.write("    {\n")
    file.write("        for (uint32_t i = 0; i < numCoeffs; i++)\n")
    file.write("            coeffs[i] += rhs.coeffs[i];\n")
    file.write("\n");
    file.write("        return *this;\n")
    file.write("    }\n\n")

    return

# WriteSubAssignOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteSubAssignOperator(file, algebra):

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator-=(const " + structName + "& rhs)\n")
    file.write("    {\n")
    file.write("        for (uint32_t i = 0; i < numCoeffs; i++)\n")
    file.write("            coeffs[i] -= rhs.coeffs[i];\n")
    file.write("\n");
    file.write("        return *this;\n")
    file.write("    }\n\n")

    return

# WriteDotAssignOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteDotAssignOperator(file, algebra):

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator|=(const " + structName + "& rhs)\n")
    file.write("    {\n")
    file.write("        " + structName + " tmp = *this;\n")
    file.write("\n");

    for r, row in enumerate(algebra.innerProdImg):

        if (len(row) == 0) : continue

        file.write("        coeffs[" + str(r) + "] = ")

        for c, col in enumerate(row):

            if col[2] == 1: 

                if c == 0 : file.write("tmp.coeffs[" + str(col[0]) + "] * rhs.coeffs[" + str(col[1]) + "]")
                else : file.write(" + tmp.coeffs[" + str(col[0]) + "] * rhs.coeffs[" + str(col[1]) + "]")

            else:

                if c == 0 : file.write("-tmp.coeffs[" + str(col[0]) + "] * rhs.coeffs[" + str(col[1]) + "]")
                else : file.write(" - tmp.coeffs[" + str(col[0]) + "] * rhs.coeffs[" + str(col[1]) + "]")

        file.write(";\n")

    file.write("\n")
    file.write("        return *this;\n")
    file.write("    }\n")
    file.write("\n")

    return

# WriteWedgeAssignOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteWedgeAssignOperator(file, algebra):

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator^=(const " + structName + "& rhs)\n")
    file.write("    {\n")
    file.write("        " + structName + " tmp = *this;\n")
    file.write("\n");

    for r, row in enumerate(algebra.outerProdImg):

        if (len(row) == 0) : continue

        file.write("        coeffs[" + str(r) + "] = ")

        for c, col in enumerate(row):

            if col[2] == 1: 

                if c == 0 : file.write("tmp.coeffs[" + str(col[0]) + "] * rhs.coeffs[" + str(col[1]) + "]")
                else : file.write(" + tmp.coeffs[" + str(col[0]) + "] * rhs.coeffs[" + str(col[1]) + "]")

            else:

                if c == 0 : file.write("-tmp.coeffs[" + str(col[0]) + "] * rhs.coeffs[" + str(col[1]) + "]")
                else : file.write(" - tmp.coeffs[" + str(col[0]) + "] * rhs.coeffs[" + str(col[1]) + "]")

        file.write(";\n")

    file.write("\n")
    file.write("        return *this;\n")
    file.write("    }\n")
    file.write("\n")

    return

# WriteMultAssignOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteMultAssignOperator(file, algebra):

    geomProdTable = algebra.innerProdImg

    for r, row in enumerate(algebra.outerProdImg) :

        for c, col in enumerate(row) : geomProdTable[r].append(col)

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator*=(const " + structName + "& rhs)\n")
    file.write("    {\n")
    file.write("        " + structName + " tmp = *this;\n")
    file.write("\n");

    for r, row in enumerate(geomProdTable):

        if (len(row) == 0) : continue

        file.write("        coeffs[" + str(r) + "] = ")

        for c, col in enumerate(row):

            if col[2] == 1: 

                if c == 0 : file.write("tmp.coeffs[" + str(col[0]) + "] * rhs.coeffs[" + str(col[1]) + "]")
                else : file.write(" + tmp.coeffs[" + str(col[0]) + "] * rhs.coeffs[" + str(col[1]) + "]")

            else:

                if c == 0 : file.write("-tmp.coeffs[" + str(col[0]) + "] * rhs.coeffs[" + str(col[1]) + "]")
                else : file.write(" - tmp.coeffs[" + str(col[0]) + "] * rhs.coeffs[" + str(col[1]) + "]")

        file.write(";\n")

    file.write("\n")
    file.write("        return *this;\n")
    file.write("    }\n")
    file.write("\n")

    return

# WriteDivAssignOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion.

def WriteDivAssignOperator(file, algebra):

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + "& operator/=(const " + structName + "& rhs)\n")
    file.write("    {\n")
    file.write("        " + structName + "tmp = ~rhs;\n")
    file.write("        T mag = rhs.Mag();\n")
    file.write("        (*this) *= tmp;\n")
    file.write("        (*this) /= (mag * mag);\n")
    file.write("        return *this;\n")
    file.write("    }\n")
    file.write("\n");

    return

# WriteAddOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion

def WriteAddOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + " operator+(const " + structName + "& rhs)\n")
    file.write("    {\n")
    file.write("        " + structName + " tmp = *this;\n")
    file.write("        tmp += rhs;\n")
    file.write("        return tmp;\n")
    file.write("    }\n")
    file.write("\n");

    return

# WriteSubOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion

def WriteSubOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + " operator-(const " + structName + "& rhs)\n")
    file.write("    {\n")
    file.write("        " + structName + " tmp = *this;\n")
    file.write("        tmp -= rhs;\n")
    file.write("        return tmp;\n")
    file.write("    }\n")
    file.write("\n");

    return

# WriteDotOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion

def WriteDotOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + " operator|(const " + structName + "& rhs)\n")
    file.write("    {\n")
    file.write("        " + structName + " tmp = *this;\n")
    file.write("        tmp |= rhs;\n")
    file.write("        return tmp;\n")
    file.write("    }\n")
    file.write("\n");

    return

# WriteWedgeOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion

def WriteWedgeOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + " operator^(const " + structName + "& rhs)\n")
    file.write("    {\n")
    file.write("        " + structName + " tmp = *this;\n")
    file.write("        tmp ^= rhs;\n")
    file.write("        return tmp;\n")
    file.write("    }\n")
    file.write("\n");

    return


# WriteMultOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion

def WriteMultOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + " operator*(const " + structName + "& rhs)\n")
    file.write("    {\n")
    file.write("        " + structName + " tmp = *this;\n")
    file.write("        tmp *= rhs;\n")
    file.write("        return tmp;\n")
    file.write("    }\n")
    file.write("\n");

    return

# WriteScalarMultOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion

def WriteScalarMultOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + " operator*(const T rhs)\n")
    file.write("    {\n")
    file.write("        " + structName + " tmp = *this;\n")
    file.write("        tmp += rhs;\n")
    file.write("        return tmp;\n")
    file.write("    }\n")
    file.write("\n");

    return

# WriteDivOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion

def WriteDivOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + " operator/(const " + structName + "& rhs)\n")
    file.write("    {\n")
    file.write("        " + structName + " tmp = *this;\n")
    file.write("        tmp /= rhs;\n")
    file.write("        return tmp;\n")
    file.write("    }\n")
    file.write("\n");

    return

# WriteScalarDivOperator -
#
# @param file       [in] File to write.
# @param algebra    [in] Algebra defintion

def WriteScalarDivOperator(file, algebra) :

    structName = algebra.name + "MV<T>"

    file.write("    " + structName + " operator/(const T rhs)\n")
    file.write("    {\n")
    file.write("        " + structName + " tmp = *this;\n")
    file.write("        tmp /= rhs;\n")
    file.write("        return tmp;\n")
    file.write("    }\n")
    file.write("\n");

    return

# GenerateAlgebra - Given an algebra signature and input name, generate its geometric
# algebra and generate C++ code for this algebra.

def GenerateAlgebra(algebraName, signature, even) :

    algebra = Algebra(algebraName, signature, even);
    file    = open(algebra.name + ".h", 'w')
    WriteMultivecStruct(file, algebra);

# Main

GenerateAlgebra("E3", [1, 1, 1], False)
GenerateAlgebra("E3", [1, 1, 1], True)
GenerateAlgebra("E2", [1, 1], False)
GenerateAlgebra("E2", [1, 1], True)