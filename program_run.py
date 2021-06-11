from qiskit import (
    #IBMQ,
    QuantumCircuit,
    QuantumRegister,
    ClassicalRegister,
    execute,
    Aer,
)
from qiskit.tools.visualization import circuit_drawer

def run(qc):
    qc.h(2)
    qc.cswap(2,0,1)
    qc.h(2)
    qc.x(2)
    qc.measure(2,0)

    #circuit_drawer(qc, filename='./swap_circuit')
    





