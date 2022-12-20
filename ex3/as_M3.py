import math

def run(qc):
    qc.mct([0, 2, 3, 4], 5)  # M3

    qc.h(2)
    qc.p(math.pi / 4, 2)

    qc.x(6)
    qc.h(7)
    qc.p(math.pi / 2, 7)

    for i in range(5):
        control = []
        control.append(6)
        for j in range(5 - i):
            control.append(j)
        qc.mct(control, 5 - i)
    qc.cnot(6, 0)

    qc.barrier([0, 1, 2, 3, 4, 5])

    for i in range(4):
        control = []
        control.append(6)
        control.append(7)
        for j in range(1, 5 - i):
            control.append(j)
        qc.mct(control, 5 - i)
    qc.ccx(6, 7, 1)

    qc.barrier([0, 1, 2, 3, 4, 5])

    for i in range(4):
        control = []
        control.append(6)
        control.append(7)
        for j in range(1, 5 - i):
            control.append(j)
        qc.mct(control, 5 - i)
    qc.ccx(6, 7, 1)

    qc.barrier([0, 1, 2, 3, 4, 5])

    for i in range(3):
        control = []
        control.append(7)
        for j in range(2, 5 - i):
            control.append(j)
        qc.mct(control, 5 - i)
    qc.cnot(7, 2)

    qc.barrier([0, 1, 2, 3, 4, 5])

    qc.measure([0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5])



