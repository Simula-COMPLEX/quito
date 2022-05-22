def run(qc):
    qc.h(0)
    qc.cx(0,1)
    qc.draw(output='mpl',filename='mpl')
    qc.measure([0,1],[0,1])
    






