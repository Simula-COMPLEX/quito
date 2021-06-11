import math

def run(qc):
    qc.h(2)
    qc.cp(math.pi/2,3,2)#Mutant
    qc.p(math.pi / 3, 2)
    qc.h(2)


    qc.cswap(2,3,5)
    qc.cswap(2,4,6)

    qc.swap(0,3)
    qc.swap(1,4)

    qc.cnot(0,1)
    qc.x(0)

    qc.measure([0,1],[0,1])






