def matmul(mata, matb):
    if len(mata[0]) != len(matb):
        raise Exception("must be same length")
    retm = [[0 for x in range(0, len(matb[0]))] for y in range(0, len(mata))]
    for y in range(0, len(mata)):
        for x in range(0, len(matb[0])):
            for i in range(0, len(matb)):
                retm[y][x] += mata[y][i]*matb[i][x]
    return retm

class Matrix:
    def __init__(self, matrix, vecc = None):
        self.matrix = self.deepcopy(matrix)
        self.matrixaug = None
        self.solvedmat = None
        self.pivotmat = None
        self.lowermat = None
        self.uppermat = None
        self.pivotlist = None
        self.vecc = self.deepcopy(vecc)
        self.size = len(matrix)

    def makeidentityaugmented(self):
        self.matrixaug = [self.matrix[i].copy() for i in range(0, len(self.matrix))]
        width = len(self.matrix)
        if width != len(self.matrix[0]):
            raise Exception("must be square matrix")
        imat = self.makeidentity(width)
        for i in range(0, width):
            self.matrixaug[i].extend(imat[i])

    def makeaugmented(self, augment):
        self.matrixaug = [self.matrix[i].copy() for i in range(0, len(self.matrix))]
        width = len(self.matrix)
        if width != len(self.matrix[0]):
            raise Exception("must be square matrix")
        for i in range(0, width):
            self.matrixaug[i].extend(augment[i])

    def makeidentity(self, size):
        mat = [[ 1 if w == h else 0 for w in range(0, size)] for h in range(0,size)]
        return mat

    def getmultiplier(self, matrix, row, mulrow):
        m = matrix[mulrow][row] / matrix[row][row]
        return m

    def sub(self, matrix, row, line):
        if len(line) != len(matrix[0]):
            raise Exception("must be same length")
        for i in range(0, len(matrix[0])):
            matrix[row][i] -= line[i]

    def mul(self, matrix, row, multiplier):
        m = []
        for i in range(0, len(matrix[0])):
            m.append(matrix[row][i]*multiplier)
        return m
    
    def replaceline(self, matrix, row, line):
        if len(line) != len(matrix[0]):
            raise Exception("must be same length")
        for i in range(0, len(matrix[0])):
            matrix[row][i] = line[i]

    def changeline(self, matrix, dest, src):
        for i in range(0, len(matrix[0])):
            tmp = matrix[src][i]
            matrix[src][i] = matrix[dest][i]
            matrix[dest][i] = tmp

    def pivot(self, matrix, row, pivot=None, lower=None):
        idx = row
        for i in range(row+1, self.size):
            if abs(matrix[i][row]) > abs(matrix[idx][row]):
                idx = i
        if idx != row:
            self.changeline(matrix, idx, row)
            if pivot:
                self.changeline(pivot, idx, row)
            if lower:
                for i in range(0, row):
                    tmp = lower[idx][i]
                    lower[idx][i] = lower[row][i]
                    lower[row][i] = tmp

    def deepcopy(self, matrix):
        mcp =  [ i.copy() for i in matrix ]
        return mcp

    def gauss(self):
        self.solvedmat = self.deepcopy(self.matrixaug)
        try:
            for i in range(0, self.size-1):
                for w in range(i, self.size-1):
                    multiplier = self.getmultiplier(self.solvedmat, i, w+1)
                    mulline = self.mul(self.solvedmat, i, multiplier)
                    self.sub(self.solvedmat, w+1, mulline)
        except ZeroDivisionError as e:
            print(e)

    def gausswithpivot(self):
        self.solvedmat = self.deepcopy(self.matrixaug)
        for i in range(0, self.size-1):
            self.pivot(self.solvedmat, i)
            for w in range(i, self.size-1):
                multiplier = self.getmultiplier(self.solvedmat, i, w+1)
                mulline = self.mul(self.solvedmat, i, multiplier)
                self.sub(self.solvedmat, w+1, mulline)

    def gaussjordan(self):
        self.makeidentityaugmented()
        try:
            self.gauss()
            for i in range(self.size-1, 0, -1):
                for w in range(i-1, -1, -1):
                    multiplier = self.getmultiplier(self.solvedmat, i, w)
                    mulline = self.mul(self.solvedmat, i, multiplier)
                    self.sub(self.solvedmat, w, mulline)
            for i in range(0, self.size):
                multiplier = 1 / self.solvedmat[i][i]
                rline = self.mul(self.solvedmat, i, multiplier)
                self.replaceline(self.solvedmat, i, rline)
            inv = [[self.solvedmat[y][x] for x in range(self.size, len(self.solvedmat[0])) ] for y in range(0,self.size)]
            return matmul(inv, self.vecc)
        except ZeroDivisionError as e:
            print(e)

    def gaussjordanwithpivot(self):
        self.makeidentityaugmented()
        self.gausswithpivot()
        for i in range(self.size-1, 0, -1):
            for w in range(i-1, -1, -1):
                multiplier = self.getmultiplier(self.solvedmat, i, w)
                mulline = self.mul(self.solvedmat, i, multiplier)
                self.sub(self.solvedmat, w, mulline)
        for i in range(0, self.size):
            multiplier = 1 / self.solvedmat[i][i]
            rline = self.mul(self.solvedmat, i, multiplier)
            self.replaceline(self.solvedmat, i, rline)
        inv = [[self.solvedmat[y][x] for x in range(self.size, len(self.solvedmat[0])) ] for y in range(0,self.size)]
        return matmul(inv, self.vecc)

    def ludecomposition(self):
        self.solvedmat = self.deepcopy(self.matrix)
        try:
            lower = self.makeidentity(self.size)
            for i in range(0, self.size-1):
                for w in range(i, self.size-1):
                    multiplier = self.getmultiplier(self.solvedmat, i, w+1)
                    lower[w+1][i] = multiplier
                    mulline = self.mul(self.solvedmat, i, multiplier)
                    self.sub(self.solvedmat, w+1, mulline)
            self.lowermat = lower
            self.uppermat = self.solvedmat
        except ZeroDivisionError as e:
            print(e)

    def ludecompositionwithpivot(self):
        self.solvedmat = self.deepcopy(self.matrix)
        pivot = []
        lower = self.makeidentity(self.size)
        for i in range(0, self.size-1):
            pivot.append(self.makeidentity(self.size))
            self.pivot(self.solvedmat,i ,pivot[i], lower)
            for w in range(i, self.size-1):
                multiplier = self.getmultiplier(self.solvedmat, i, w+1)
                lower[w+1][i] = multiplier
                mulline = self.mul(self.solvedmat, i, multiplier)
                self.sub(self.solvedmat, w+1, mulline)
        self.lowermat = lower
        self.uppermat = self.solvedmat
        self.pivotlist = pivot

        p = self.makeidentity(self.size)
        for i in range(len(pivot)-1, -1, -1):
            p = matmul(p, pivot[i])
        d = matmul(p, self.vecc)
        y = [[] for i in range(0, self.size)]

        for i in range(0, self.size):
            for w in range(0, i+1):
                if y[w]:
                    d[i][0] -= lower[i][w]*y[w][0]
                else:
                    y[i].append(d[i][0]/lower[i][w])
        x = [[] for i in range(0, self.size)]
        for i in range(self.size-1, -1, -1):
            for w in range(self.size-1, i-1, -1):
                if x[w]:
                    y[i][0] -= self.uppermat[i][w]*x[w][0]
                else:
                    x[i].append(y[i][0]/self.uppermat[i][w])
        return x

    def __matrixprint__(self, matrix):
        for h in range(0, len(matrix)):
            for w in range(0, len(matrix[0])):
                print("{:>5.2f}".format(matrix[h][w]), end='')
            print()
        print()

    def matrixprint(self):
        self.__matrixprint__(self.matrix)
    def augmatrixprint(self):
        self.__matrixprint__(self.matrixaug)
    def solvedprint(self):
        self.__matrixprint__(self.solvedmat)
    def luprint(self):
        print("lower matrix")
        self.__matrixprint__(self.lowermat)
        print("upper matrix")
        self.__matrixprint__(self.uppermat)
    def pivotprint(self):
        for i in range(0, len(self.pivotlist)):
            print("pivot matrix {} :".format(i))
            self.__matrixprint__(self.pivotlist[i])

if __name__ == "__main__":
    mm = [[2.0, 4.0, -2.0, -2.0],
    [1.0, 2.0, 4.0, -3.0],
    [-3.0, -3.0, 8.0, -2.0],
    [-1.0, 1.0, 6.0, -3.0]]
    mc = [[-4.0],
    [5.0],
    [7.0],
    [7.0]]
    
    mat1 = Matrix(mm, mc)
    print("gaussjordan :")
    x = mat1.gaussjordan()
    mat1.solvedprint()
    if x:
        print('x :')
        mat1.__matrixprint__(x)
    print("gaussjordan with pivot :")
    x = mat1.gaussjordanwithpivot()
    mat1.solvedprint()
    print('x :')
    mat1.__matrixprint__(x)
    print("ludecompoistion :")
    x = mat1.ludecomposition()
    if x:
        mat1.__matrixprint__(x)
    print("ludecompoistion with pivot :")
    x = mat1.ludecompositionwithpivot()
    mat1.luprint()
    mat1.pivotprint()
    print('x :')
    mat1.__matrixprint__(x)

    mm = [[1.0, 5.0, -1.0, 6.0],
    [2.0, -1.0, 1.0, -2.0],
    [-1.0, 4.0, -1.0, 3.0],
    [3.0, -7.0, -2.0, 1.0]]
    mc = [[19.0],
    [7.0],
    [20.0],
    [-75.0]]

    mat2 = Matrix(mm, mc)
    print("gaussjordan :")
    x = mat2.gaussjordan()
    mat2.solvedprint()
    if x:
        print('x :')
        mat2.__matrixprint__(x)
    print("gaussjordan with pivot :")
    x = mat2.gaussjordanwithpivot()
    mat2.solvedprint()
    print('x :')
    mat2.__matrixprint__(x)
    print("ludecompoistion :")
    x = mat2.ludecomposition()
    if x:
        mat2.__matrixprint__(x)
    print("ludecompoistion with pivot :")
    x = mat2.ludecompositionwithpivot()
    mat2.luprint()
    mat2.pivotprint()
    print('x :')
    mat2.__matrixprint__(x)