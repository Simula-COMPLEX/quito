#from quito import quito

#quito("/home/xinyi/quito/ex3/as.ini")


from quito.quito_coverage import quito
import time
start = time.time()
i = 0
quito('/home/xinyi/quito/ex3/as.ini')
end = time.time()

print(end - start)

